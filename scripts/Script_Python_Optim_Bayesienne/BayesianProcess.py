################################################################################
##### This file contains the code that implement the bayesian process
################################################################################
##### Initial version : Thomas CIRON, 2022-2023
##### Revision : Ronan BOCQUILLON, september 2024
##### Revision : Vincent T'KINDT, october 2024
##### Revision : Quentin SCHAU, january 2026
################################################################################

#Inclusion of mandatory librairies
import json
import bayes_opt
import os
import numpy as np
from pathlib import Path

#Inclusion of other librairies from the project
import Evaluation


class BayesianModel:

    def __init__(self, json_conf_path, function_to_maximize):
        """!
        Constructor of BayesianModel. 
        @param json_conf_path the json configuration file path, we read almost all the parameters from this file.
        @param function_to_maximize the function we want to maximize.
        """
        # parse json file
        with open(json_conf_path, "r") as f:
            data = json.load(f)
        self.pbounds = {}
        for hp in data["hyperparameters"]:
            if 'isInteger' in data["hyperparameters"][hp] and data["hyperparameters"][hp]["isInteger"]:
                self.pbounds[hp] = (data["hyperparameters"][hp]["min"],data["hyperparameters"][hp]["max"],int)
            else:
                self.pbounds[hp] = (data["hyperparameters"][hp]["min"],data["hyperparameters"][hp]["max"])
        print(self.pbounds)

        self.init_points = data["bayesian optimization parameters"]["number of init points"]
        self.nb_iterations = data["bayesian optimization parameters"]["number of iterations"]
        self.acquisition_function_kind = data["bayesian optimization parameters"]["acquisition function"]
        self.xi = data["bayesian optimization parameters"]["xi"]
        self.kappa = data["bayesian optimization parameters"]["kappa"]
        self.kappa_decay = data["bayesian optimization parameters"]["kappa_decay"]
        self.kappa_decay_delay = data["bayesian optimization parameters"]["kappa_decay_delay"]
        self.allowDuplicatePoints = data["bayesian optimization parameters"]["allow_duplicate_points"]
        self.probes=data["probe"]
        if self.acquisition_function_kind == "ucb":
            self.acquisition_function = bayes_opt.acquisition.UpperConfidenceBound(kappa=self.kappa, exploration_decay = self.kappa_decay, exploration_decay_delay=self.kappa_decay_delay)
        elif self.acquisition_function_kind == "ei":
            self.acquisition_function = bayes_opt.acquisition.ExpectedImprovement(xi=self.xi, exploration_decay = self.kappa_decay, exploration_decay_delay=self.kappa_decay_delay)
        elif self.acquisition_function_kind == "pi":
            self.acquisition_function = bayes_opt.acquisition.ProbabilityOfImprovement(xi=self.xi, exploration_decay = self.kappa_decay, exploration_decay_delay=self.kappa_decay_delay)
        # bayesian optimizer
        self.optimizer = bayes_opt.BayesianOptimization(
            f=function_to_maximize,
            pbounds=self.pbounds,
            allow_duplicate_points=self.allowDuplicatePoints,
            acquisition_function=self.acquisition_function,
            verbose=0
        )
    def Train(self):
        """!
        Train our bayesian optimizer model on a base of sets. Saves the result in a .json file.
        """
        print("############################################")
        print("Learning process starting")
        print("############################################")
        my_file = Path("./logs.json")
        if my_file.is_file():
            self.optimizer.load_state("./logs.json")
            # if there is user-defined points
            if len(self.probes["points"])>0:
                for point in self.probes["points"]:
                    self.optimizer.probe(point,lazy=True)
                self.optimizer.maximize(init_points=0, n_iter=0)
            else:
                print("Starts from logs.json with nb iter:",self.nb_iterations)
                self.optimizer.maximize(init_points=0, n_iter=self.nb_iterations)
        else:
            # if there is user-defined points
            if len(self.probes["points"])>0:
                for point in self.probes["points"]:
                    self.optimizer.probe(point,lazy=True)
                self.optimizer.maximize(init_points=0, n_iter=0)
            else:
                print("Starts from nothing with nb init: ",self.init_points," and nb iter:",self.nb_iterations)
                self.optimizer.maximize(init_points=self.init_points, n_iter=self.nb_iterations)

        print("############################################")
        print("Save Model")
        print("############################################")
        self.optimizer.save_state("./logs.json")
        # reset logger to verbose mode and to iteration 0
        self.optimizer.logger._verbose = 2
        self.optimizer.logger._iterations = 0
        self.optimizer.logger.log_optimization_start(self.optimizer._space.keys) 
        for res in self.optimizer._space.res():
            self.optimizer.logger.log_optimization_step(
                self.optimizer._space.keys, res, self.optimizer._space.params_config, self.optimizer.max
            )
        
        self.optimizer.logger.log_optimization_end() 
        
    def Validate(self):
        """!
        Test the model on a validation base.
        """
        print("############################################")
        print("Computing the results on the validation base")
        print("############################################")
        my_file = Path("./logs.json")
        if my_file.is_file():
            self.optimizer.load_state("./logs.json")
        best_parameters = self.optimizer.max["params"]
        with open("best_params.json", "w") as f:
            json.dump(best_parameters, f)
        Evaluation.Compute_deviation(best_parameters)

    def Get_sum_uncertainties(self):
        """!
        Get the sum of uncertainty.
        @return the sum of uncertainty.
        """
        n = 10000
        X = np.array([self.tmp_optimizer._space.random_sample() for i in range(n)])
        mu, sigma = self.tmp_optimizer._gp.predict(X, return_std=True)
        return np.sum(sigma **2)

    # def Get_loss_curve(self):
    #     # parse json file
    #     print("############################################")
    #     print("Computing the loss curve")
    #     print("############################################")
    #     self.tmp_optimizer = bayes_opt.BayesianOptimization(f=lambda x : x,pbounds=self.pbounds)
    #     noises = [self.Get_sum_uncertainties()]
    #     logs = Logs.get_logs()
    #     for log in logs:
    #         with open(log, "r") as f:
    #             data = f.read().split("\n")
    #         for iteration in data[:-1]:
    #             try:
    #                 a = json.loads(iteration)
    #                 self.tmp_optimizer.register(params=a["params"], target=a["target"])
    #                 self.tmp_optimizer._gp.fit(self.tmp_optimizer._space.params, self.tmp_optimizer._space.target)
    #                 noises.append(self.Get_sum_uncertainties())
    #             except KeyError:
    #                 pass
    #     plt.plot(noises[2:]/max(noises[2:]))
    #     plt.title("Courbe de perte de l'optimisation bayésienne")
    #     plt.show()
