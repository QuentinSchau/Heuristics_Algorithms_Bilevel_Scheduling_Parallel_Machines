#!/usr/bin/env python
# coding: utf-8

# Installation 
# 
# 1.  We recommend to use conda for its installation cf. [the documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
# 2.  Create and use a virtual environment. `conda create -n TP python==3.9 -y`
# 3.  You need to use the environment with : `conda activate TP`
# 4.  Install the dependencies `pip install -r requirements.txt`
# 
# You should change the version of PyTorch to compute on your GPU, it depends on your system and CUDA driver [cf official site web](https://pytorch.org/get-started/locally/).
# 
# You can also run the following pip installation :
# 1. `pip install torch`
# 2. `pip install tensorboard`
# 3. `pip install tensorboardX`
# 4. `pip install matplotlib`
# 5. `pip install lightning`
# 

# HYPER PARAMS

HIDDEN_UNITS = []
NAME="P-v1-b32-MSE-loss-cst-lr"
VERSION=1
INPUT_SIZE = 96
BATCH_SIZE = 32
NUM_EPOCHS = 6000
LEARNING_RATE = 0.000001
NUM_WORKERS = 4
PREFETCH_FACTOR=2

import numpy as np
import torch    
from torch import nn
import lightning as L
import torchmetrics
import matplotlib
matplotlib.use('Agg')
import time
from lightning.pytorch.callbacks import ModelCheckpoint,LearningRateMonitor
from lightning.pytorch.loggers import TensorBoardLogger
import pathlib
import os

class MLP(nn.Module):
    def __init__(self, input_size, hidden_units):
        super().__init__()
    
        # Initialize MLP layers
        all_layers = []
        for hidden_unit in hidden_units:
            layer = nn.Linear(input_size, hidden_unit, bias=True)
            all_layers.append(layer)
            all_layers.append(nn.ReLU())
            input_size = hidden_unit
    
        output_layer = nn.Linear(
            in_features=input_size,
            out_features=1)
        
        all_layers.append(output_layer)
    
        # add logSoftmax layer
        self.layers = torch.nn.Sequential(*all_layers)

    def forward(self, x):
        logits = self.layers(x)
        return  logits

# LightningModule that receives a PyTorch model as input
class LightningModel(L.LightningModule):
    def __init__(self, model,name, learning_rate):
        super().__init__()
        self.name = name
        self.learning_rate = learning_rate
        # The inherited PyTorch module
        self.model = model
        self.logDir = "Log"
        # Save settings and hyperparameters to the log directory
        # but skip the model parameters
        self.save_hyperparameters(ignore=["model"])

        self.train_mse = torchmetrics.MeanSquaredError()
        self.train_mae = torchmetrics.MeanAbsoluteError()
        self.val_mse = torchmetrics.MeanSquaredError()
        self.val_mae = torchmetrics.MeanAbsoluteError()
        self.test_mse = torchmetrics.MeanSquaredError()
        self.test_mae = torchmetrics.MeanAbsoluteError()

        
    # Defining the forward method is only necessary
    # if you want to use a Trainer's .predict() method (optional)
    def forward(self, x):
        return self.model(x)

    # A common forward step to compute the loss and labels
    # this is used for training, validation, and testing below
    def _shared_step(self, batch):
        features, true_labels = batch
        logits = self(features)
        loss = torch.nn.functional.mse_loss(logits, true_labels)
        return loss, true_labels, logits

    def training_step(self, batch, batch_idx):
        loss, true_labels, logits = self._shared_step(batch)
        self.log("train_loss",loss,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)

        # Do another forward pass in .eval() mode to compute accuracy
        # while accountingfor Dropout, BatchNorm etc. behavior
        # during evaluation (inference)           
        self.train_mae(logits,true_labels)
        self.train_mse(logits,true_labels)
        
        self.log("train_mae",self.train_mae,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        self.log("train_mse",self.train_mse,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        self.model.train()

        return loss  # this is passed to the optimzer for training

    def validation_step(self, batch, batch_idx):
        loss, true_labels, logits = self._shared_step(batch)
        self.log("val_loss",loss,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        
        self.val_mae(logits,true_labels)
        self.val_mse(logits,true_labels)
        
        self.log("val_mae",self.val_mae,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        self.log("val_mse",self.val_mse,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        

    def test_step(self, batch, batch_idx):
        loss, true_labels, logits = self._shared_step(batch)
        self.log("test_loss",loss,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        
        self.test_mae(logits,true_labels)
        self.test_mse(logits,true_labels)
        
        self.log("test_mae",self.test_mae,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)
        self.log("test_mse",self.test_mse,on_epoch=True,on_step=False,prog_bar=True,sync_dist=True)

    #########################
    #       OPTIMIZER       #
    #########################
    
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        # scheduler = torch.optim.lr_scheduler.LinearLR(optimizer, start_factor=1.0,end_factor=0.0001, total_iters=800)
        # lr_scheduler = {
        #     'scheduler':scheduler,
        #     'name': 'learning rate'
        # }
        # return [optimizer], [lr_scheduler]
        return [optimizer]
      
if __name__ == "__main__":
    
    ###########  A completer #############

    # Definie un model pytorch
    pytorch_model = MLP(
        input_size=INPUT_SIZE,
        hidden_units=HIDDEN_UNITS,
    )

    print("Load model")
    lightning_model = LightningModel(pytorch_model,NAME,LEARNING_RATE)
    checkpoint_path = pathlib.Path(str(os.path.dirname(__file__))+ f"/Log/{NAME}/version_{VERSION}/checkpoints/mae-epoch-epoch=00029.ckpt")
    torch.set_float32_matmul_precision('medium')
    model = LightningModel.load_from_checkpoint(checkpoint_path, model=pytorch_model).model.to("cpu").eval()
    # save for use in production environment
    torch.jit.trace(model,torch.ones((1,96)).to("cpu")).save("predictor.pt")
    print(f"Model saved")
