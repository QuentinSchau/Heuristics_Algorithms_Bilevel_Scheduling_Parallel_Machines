################################################################################
##### This file contains the main code for optimizing parameters of an algorithm 
################################################################################
##### Initial version : Thomas CIRON, 2022-2023
##### Revision : Ronan BOCQUILLON, september 2024
##### Revision : Vincent T'KINDT, october 2024
##### Revision : Quentin SCHAU, january 2026
################################################################################

#Inclusion of mandatory librairies
import subprocess
import json
import os
import numpy as np

#Inclusion of other librairies from the project
import Algorithm, BayesianProcess, Evaluation

#Local variables
DEBUG=True

################################################################################
##### Function to initialize the optimization process
################################################################################

def InitOptimizationProcess(json_conf_path):
    """!
    Reads the file conf.json, Tag  "bayesian optimization parameters" to grab the values of the parameters to run the learning process
    @param json_conf_path, str, name of the configuration file
    @return AlgName, str, the name of the executable file associated to the algorithm which parameters have to be optimized.
    @return AlgRepetitions, int, the number of execution to perform for the algorithm on any instance
    @return TrainDD, str, the path of the directory containing the training instances
    @return TrainRD, str, the path of the directory containing the baseling results on the training instances
    @return TrainNumInst, int, number of instances in the training database
    @return ValidDD, str, the path of the directory containing the validation instances
    @return ValidRD, str, the path of the directory containing the baseling results on the validation instances
    @return ValidNumInst, int, number of instances in the validation database
    """
    # parse json file
    with open(json_conf_path, "r") as f:
        data = json.load(f)
        if (DEBUG==True):
            print("############################################")
            print("Parameters of the learning process")
            print("############################################")
            print("AlgName=",data["bayesian optimization parameters"]["Algorithm"])
            print("AlgResFile=",data["bayesian optimization parameters"]["AlgorithmResFile"])
            print("AlgRepetitions=",data["bayesian optimization parameters"]["AlgRepetitions"])
            print("MinNumJob=",data["bayesian optimization parameters"]["MinNumJob"])
            print("MaxNumJob=",data["bayesian optimization parameters"]["MaxNumJob"])
            print("StepNumJob=",data["bayesian optimization parameters"]["StepNumJob"])
            print("TrainDD=",data["bayesian optimization parameters"]["TrainDataDirectory"])
            print("TrainRD=",data["bayesian optimization parameters"]["TrainResDirectory"])
            print("TrainNumInst=",data["bayesian optimization parameters"]["TrainNumInst"])
            print("TrainNumInstpJob=",data["bayesian optimization parameters"]["TrainNumInstpJob"])
            print("ValidDD=", data["bayesian optimization parameters"]["ValidDataDirectory"])
            print("ValidRD=",data["bayesian optimization parameters"]["ValidResDirectory"])
            print("ValidNumInst=",data["bayesian optimization parameters"]["ValidNumInst"])
            print("ValidNumInstpJob=",data["bayesian optimization parameters"]["ValidNumInstpJob"])
        return (data["bayesian optimization parameters"]["Algorithm"], data["bayesian optimization parameters"]["AlgorithmResFile"],
                data["bayesian optimization parameters"]["AlgRepetitions"], data["bayesian optimization parameters"]["MinNumJob"],data["bayesian optimization parameters"]["MaxNumJob"], 
                data["bayesian optimization parameters"]["StepNumJob"], data["bayesian optimization parameters"]["TrainDataDirectory"],
                data["bayesian optimization parameters"]["TrainResDirectory"],data["bayesian optimization parameters"]["TrainNumInst"],
                data["bayesian optimization parameters"]["TrainNumInstpJob"],data["bayesian optimization parameters"]["ValidDataDirectory"],
                data["bayesian optimization parameters"]["ValidResDirectory"],data["bayesian optimization parameters"]["ValidNumInst"],
                data["bayesian optimization parameters"]["ValidNumInstpJob"])


################################################################################
##### Main function of the program
################################################################################

if __name__=="__main__":
    #Reading the configuration file
    (Algorithm.ALGName,Algorithm.ALGResFile, Algorithm.ALGRepetitions,
     Evaluation.EVLMinNumJob, Evaluation.EVLMaxNumJob,Evaluation.EVLStepNumJob,
     Evaluation.EVLTrainDD,Evaluation.EVLTrainRD,Evaluation.EVLTrainNumInst,Evaluation.EVLTrainNumInstpJob, Evaluation.EVLValDD,Evaluation.EVLValRD, 
     Evaluation.EVLValNumInst,Evaluation.EVLValidNumInstpJob )=InitOptimizationProcess("conf.json")
    #Initializing the Bayesian process
    model = BayesianProcess.BayesianModel("conf.json", Evaluation.Function)

    #Running the Bayesian process
    model.Train()
        
    #Run the validation tests on the validation database
    # model.Validate()
