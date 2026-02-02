#!/usr/bin/env python3
#  Copyright (C) 2024
#  Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
#
#  DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
#  This file is part of bilevel-scheduling.
#
#  bilevel-scheduling is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  bilevel-scheduling is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
import numpy as np
import os
import pandas as pd
import pathlib
import json
import multiprocessing

# CONSTANT
M = [(1,1),(2,2)]
NList = [40,45,50,60,75,80,100,125]
TF = [0.2,0.4,0.6,0.8,1.0]
RDD = [0.2,0.4,0.6,0.8,1.0]
Vmax = 2.0
V0 = 1.0
NB_INST_GENERATED = 100
NAME="Train"

NB_PROC = multiprocessing.cpu_count()
seed = 1
base_dir = pathlib.Path(str(os.path.dirname(__file__))+ f"/{NAME}")
pathProgramSolver = pathlib.Path(str(pathlib.Path(__file__).parent.parent.parent.resolve())+"/bin/Release/bilevel-scheduling")

if not pathProgramSolver.is_file():
    print(f"Bilevel Scheduler solver does not exist as the path {pathProgramSolver}.")
    exit(1)
# heuristic 
BEAMSIZE = 10

def runConfigFile(configFile):
    commandLine = str(pathProgramSolver) + " "+ str(configFile)
    os.system(commandLine)

def generateInstanceFile():
    

    # create JSON object to generate instance
    
    listNameInstances = []
    listConfigInstance = []
    generate = {
        "seed": seed
    }
    print(f"The database {NAME} will be generated with {len(M)*len(NList)*len(RDD)*len(TF)*NB_INST_GENERATED} instances")
    
    for (mMax, M0) in M:
        for N in NList:
            n = N
            for tf in TF:
                for rdd in RDD:
                    # Keep the list of instances for make the solve config file
                    instancePath = "_n_" + str(n) + "_N_" + str(N) + "_tf_" + str(tf) + "_rdd_" + str(
                        rdd) + "_mMax_" + str(mMax) + "_m0_" + str(M0) + "_Vmax_" + str(Vmax) + "_V0_" + str(V0)
                    listNameInstances.append(instancePath)

                    instance = {
                        "basePath": str(base_dir)+ "/instances/",
                        "numberInstance": NB_INST_GENERATED}
                    param = {
                            "n": n,
                            "N": N,
                            "M_max": mMax,
                            "M_0": M0,
                            "V_max": Vmax,
                            "V_0": V0,
                            "pi": {
                                "inf": 1,
                                "sup": 100
                            },
                            "di": {
                                "TF": tf,
                                "RDD": rdd
                            },
                            "wi": {
                                "inf": 1,
                                "sup": 10
                            }
                            }
                    instance["paramInstance"] = param

                    listConfigInstance.append(instance)

    generate["instances"] = listConfigInstance
    configFile = {}
    configFile["generate"] = generate
    
    configPathGenInstances = pathlib.Path(str(base_dir) + f"/configuration_file_for_generation_database.json").resolve()
    configPathGenInstances.parent.mkdir(parents=True, exist_ok=True)
    with open(configPathGenInstances, 'w') as file:
        file.write(json.dumps(configFile))
    
    print(f"Generate configuration file for database:\n- {configPathGenInstances}")
    runConfigFile(configPathGenInstances)

    # create JSON object to solve all instances with heuristic to get upper bound
    print(f"Generate configuration files for upper bound")
    for x in listNameInstances:
        for i in range(NB_INST_GENERATED):
            configSolve = {"solve": {"verbose": 0, "output": str(base_dir) + f"/", "methods": []}}
            configSolve["solve"]["methods"].append({
                    "name": "BeamSearch",
                    "verbose": 0,
                    "beamSize": BEAMSIZE,
                    "timeLimits": 300,
                    "recovering": True,
                    "reco-strategy": "local-search",
                    "nbSolutionForMSLS": 1000,
                    "version": 3,
                    "alpha": 0.0,
                    "LocalSearch_parameters": {
                        "name": "LocalSearch",
                        "timeLimits": 100,
                        "version": 1,
                        "verbose": 0,
                        "useBestNeighbor": False,
                        "maxIter": 100
                    },                
                    "instances": [{"path": str(base_dir) + f"/instances/instance{i}{x}.txt"}]
                })    
            configPathGenUB = pathlib.Path(str(base_dir) + f"/config_files/config_file_gen_UB_i{i}{x}.json").resolve()
            configPathGenUB.parent.mkdir(parents=True, exist_ok=True)
            with open(configPathGenUB, "w") as file:
                file.write(json.dumps(configSolve))
            
    # create JSON object to compute features vectors
    print(f"Generate configuration file for features:")
    # create header of the csv file 
    with open(str(base_dir)+f"/resultsLocalSearch.csv",'w+') as file:
        file.write("".join(["InstanceName\t"] + [f"feat_{x}\t" for x in range(1, 96)] + ["feat_96"]) + "\n")
    for x in listNameInstances:
        for i in range(NB_INST_GENERATED):
            configSolve = {"solve": {"verbose": 0, "output": str(base_dir) + f"/", "methods": []}}
            configSolve["solve"]["methods"].append({
                    "name": "LocalSearch",
                    "verbose": 0,
                    "version": 6,
                    "timeLimits": 60,
                    "genDatabase": True,
                    "instances": [{"path": str(base_dir) + f"/instances/instance{i}{x}.txt"}]
                })    
            configPathGenUB = pathlib.Path(str(base_dir) + f"/config_files/config_file_gen_UB_i{i}{x}_features.json").resolve()
            configPathGenUB.parent.mkdir(parents=True, exist_ok=True)
            with open(configPathGenUB, "w") as file:
                file.write(json.dumps(configSolve))
    
    print("Generate database...")
    files = [str(f) for f in pathlib.Path(str(base_dir) + f"/config_files/").iterdir() if f.is_file()]
    with multiprocessing.Pool(processes=NB_PROC) as pool:
        pool.map(runConfigFile, files)
    
def createDatabase():
    # merge dataframe
    df_ub = pd.read_csv(str(base_dir)+f"/resultsBeamSearch.csv",sep="\t")
    df_feat = pd.read_csv(str(base_dir)+f"/resultsLocalSearch.csv",sep="\t")
    df = pd.merge(df_ub,df_feat,on="InstanceName")
    headers = ['InstanceName','InstancePath','N','n','m_Max','m_0','V_max', 'V_0','Method', 'Time', 'TimeInLS', 'LimitTime','First_UB','BeamSize', 'Alpha', 'UseRecovering', 'RecoStrategy', 'VersionRBS', 'VersionLS', 'UseBestImprovement', 'NBCallsReco', 'NB_Best_NBH_FindWithReco', 'NB_MSLS_SOLUTIONS', 'NB_MSLS_BEST_SOL', 'NB_MSLS_LEAF_SOL', 'NBNodes', 'NBMinColum', 'GenerateCol', 'MaxNbCallHeuristics', 'NBComputeCost', 'NBCallsHeu', 'NBCallDP', 'NBCallSubProcessCG']
    featuresHeader = [f"feat_{i}" for i in range(1,97)]
    labels_array = df.drop(columns=headers+featuresHeader).values.astype('float32')
    features_array = df.drop(columns=["Objective"]+headers).values.astype('float32')
    np.save(f"{NAME}/features.npy", features_array)
    np.save(f"{NAME}/labels.npy", labels_array)
    

if __name__ == "__main__":
    generateInstanceFile()
    createDatabase()
