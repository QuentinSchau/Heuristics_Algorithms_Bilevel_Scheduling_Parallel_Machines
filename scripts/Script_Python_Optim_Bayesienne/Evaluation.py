################################################################################
##### This file contains the code which evaluates the optimization function f
################################################################################
##### Initial version : Thomas CIRON, 2022-2023
##### Revision : Ronan BOCQUILLON, september 2024
##### Revision : Vincent T'KINDT, october 2024
##### Revision : Quentin SCHAU, january 2026
################################################################################

#Inclusion of mandatory librairies
import numpy as np
from multiprocessing import Pool

#Inclusion of other librairies from the project
import Algorithm
from tqdm import tqdm
#Local variables with their default value
EVLNumCPU=104
EDEBUG=False
EVLMinNumJob=0               #Minimum number of jobs in an instance
EVLMaxNumJob=0               #Maximum number of jobs in an instance
EVLStepNumJob=0              #Increment in the number of jobs
EVLTrainDD=""                #Directory containing all the instance files of the training database
EVLTrainRD=""                #Directory containing all the reference results of the instance files from the training database
EVLTrainNumInst=0            #Number of instances in the training database
EVLTrainNumInstpJob=0        #Number of instances, having the same number of jobs, in the training database
EVLValDD=""                  #Directory containing all the instance files of the validation database
EVLValRD=""                  #Directory containing all the reference results of the instance files from the validation database
EVLNumInst=0                 #Number of instances in the validation database
EVLValidNumInstpJob=0        #Number of instances, having the same number of jobs, in the validation database
################################################################################
### Compute_Improvement: compute the improvement of the algorithm for a current
### Set of parameters wrt to the references results in the training database
################################################################################
### Notice that the algorithm must have been ran on all the instances of the
### training database before running this function
################################################################################
def Compute_Improvement(nb_jobs:int)->float:
    """!
    Get the amelioration score on the train database. Compute this for a given number of jobs.
    Let vd_i be the reference result with the original settings of the instance i and va_i be the result with the tuned settings of the instance i.
    The deviation is computed so : dev = sum (vd_i - va_i)/vd_i for i from 0 to NB_INSTANCES_PER_JOBS.
    @param nb_jobs the number of jobs we want the deviations from.
    @return the deviation
    """
    improvement = 0.0
    for i in (tqdm(range(EVLTrainNumInstpJob),desc="Compute Improvement") if EDEBUG else range(EVLTrainNumInstpJob)):
        ub = Algorithm.Get_UpperBound(EVLTrainDD,nb_jobs,i)
        mhb_objective_function_value = Algorithm.Get_objective_function_value(f"{Algorithm.ALGResFile}/donnees_{nb_jobs}_{i+1}.dat.seq")
        mh_objective_function_value = Algorithm.Get_objective_function_value(f"{EVLTrainRD}/donnees_{nb_jobs}_{i+1}.dat.seq")
        improvement += (mh_objective_function_value-mhb_objective_function_value)/ub
    return improvement

################################################################################
### This is required for multiprocessing to pickle the function correctly.
################################################################################
def run_worker(args):
    """
    Unpacks arguments and calls Algorithm.Run.
    args is a tuple: (EVLDatabase, nb_jobs, alpha1...alpha9, instance_idx)
    """
    # Unpack the tuple. Order must match how we pack it below.
    (evl_dd, nb_jobs, a1, a2, a3, a4, a5, a6,a7,a8,instance_idx) = args
    
    # Call the actual algorithm
    Algorithm.Run(evl_dd, nb_jobs, a1, a2, a3, a4, a5, a6,a7,a8, instance_idx)

################################################################################
### Function: Compute the function that is maximized by the Bayesian process
################################################################################
### Notice that the algorithm must have been ran on all the instances of the
### training database before running this function
### Note: some parts of the code below has to be adapted to your problem
################################################################################
def Function(alpha1:float,alpha2:float,alpha3:float,alpha4:float,alpha5:float,alpha6:float,alpha7:float,alpha8:float)->float:
    """!
    Compute the function value for a given set of learned values.
    """
    print("Evaluate at : ",alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8)
    
    # if we have the lowest value of the term alpha that is greater than 1.0 then return negative value. this allow to implement the Constrained Optimization cf "https://bayesian-optimization.github.io/BayesianOptimization/3.2.0/constraints.html"
    if alpha5 * 40 + alpha6*10 + alpha7*2 + alpha8 > 1.0:
        return -1000.0
    if alpha5 * 80 + alpha6*60 + alpha7*4 + alpha8 < 0.0:
        return -1000.0
    # 0. For each repetition
    # ---------------------------------------------------------
    improvement = 0.0
    # For every repetition and every instance, create a task
    for _ in range(Algorithm.ALGRepetitions):
        # 1. Generate the list of ALL tasks to be run

        # ---------------------------------------------------------
        all_tasks = []

        # We keep track of the nb_jobs steps to calculate improvement later
        nb_jobs_steps = [] 

        current_nb_jobs = EVLMinNumJob
        while current_nb_jobs <= EVLMaxNumJob:
            nb_jobs_steps.append(current_nb_jobs)
            
            for instance_idx in range(EVLTrainNumInstpJob):
                # Pack all necessary arguments into a tuple
                task_args = (
                    EVLTrainDD, 
                    current_nb_jobs, 
                    alpha1, alpha2, alpha3, alpha4, alpha5, alpha6,alpha7,alpha8,
                    instance_idx
                )
                all_tasks.append(task_args)
            
            current_nb_jobs += EVLStepNumJob

        if EDEBUG:
            print(f"  Queueing {len(all_tasks)} total tasks across {len(nb_jobs_steps)} job configurations.")

        # 2. Run ALL tasks in parallel using a single Pool
        # ---------------------------------------------------------
        # This saturates the CPUs because we don't wait between nb_jobs steps
        with Pool(EVLNumCPU) as p:
            if EDEBUG:
                # imap_unordered is often slightly faster if order doesn't matter, 
                # and allows for a smooth progress bar
                list(tqdm(p.imap_unordered(run_worker, all_tasks), total=len(all_tasks), desc="Processing all instances"))
            else:
                p.map(run_worker, all_tasks)

        # 3. Compute Improvement (Post-Processing)
        # ---------------------------------------------------------
        # Now that all files/results are generated, we calculate the score.
        for nb_jobs in nb_jobs_steps:
            improvement += Compute_Improvement(nb_jobs)
                
    return improvement / Algorithm.ALGRepetitions


################################################################################
### Compute_deviation: Compute statistics on the validation database for
### the best learned parameters
################################################################################
### Notice that this can be called only after the learning process stopped
### Note: some parts of the code below has to be adapted to your problem
################################################################################
def Compute_deviation(Best_parameter):
    """!
    Compute the deviation on the validation database of the algorithm with the best learned parameters wrt to the reference results in the validation database
    @param best_parameters, a structure, containing the best learned values
    """

    # initializing arrays for the datas
    dev_res = [[] for _ in range(EVLMinNumJob,EVLMaxNumJob+EVLStepNumJob,EVLStepNumJob)]
    dev_opt = [[] for _ in range(EVLMinNumJob,EVLMaxNumJob+EVLStepNumJob,EVLStepNumJob)]
    times_opt = [[] for _ in range(EVLMinNumJob,EVLMaxNumJob+EVLStepNumJob,EVLStepNumJob)]
    times_res = [[] for _ in range(EVLMinNumJob,EVLMaxNumJob+EVLStepNumJob,EVLStepNumJob)]
    # For every repetition and every instance, create a task
    for _ in range(Algorithm.ALGRepetitions):
        # 1. Generate the list of ALL tasks to be run

        # ---------------------------------------------------------
        all_tasks = []

        # We keep track of the nb_jobs steps to calculate improvement later
        nb_jobs_steps = [] 

        current_nb_jobs = EVLMinNumJob
        while current_nb_jobs <= EVLMaxNumJob:
            nb_jobs_steps.append(current_nb_jobs)
            
            for instance_idx in range(EVLTrainNumInstpJob):
                # Pack all necessary arguments into a tuple
                task_args = (
                    EVLValDD, 
                    current_nb_jobs, 
                    Best_parameter["alpha1"], Best_parameter["alpha2"], Best_parameter["alpha3"], Best_parameter["alpha4"], Best_parameter["alpha5"], Best_parameter["alpha6"],Best_parameter["alpha7"],Best_parameter["alpha8"],
                    instance_idx
                )
                all_tasks.append(task_args)
            
            current_nb_jobs += EVLStepNumJob

        if EDEBUG:
            print(f"  Queueing {len(all_tasks)} total tasks across {len(nb_jobs_steps)} job configurations.")

        # 2. Run ALL tasks in parallel using a single Pool
        # ---------------------------------------------------------
        # This saturates the CPUs because we don't wait between nb_jobs steps
        with Pool(EVLNumCPU) as p:
            if EDEBUG:
                # imap_unordered is often slightly faster if order doesn't matter, 
                # and allows for a smooth progress bar
                list(tqdm(p.imap_unordered(run_worker, all_tasks), total=len(all_tasks), desc="Processing all instances"))
            else:
                p.map(run_worker, all_tasks)

        # 3. Stats (Post-Processing)
        # ---------------------------------------------------------
        # Now that all files/results are generated, we calculate the score.
        indexLoopStatsPerJobs=0
        for nb_jobs in range(EVLMinNumJob,EVLMaxNumJob+EVLStepNumJob,EVLStepNumJob):
            for index_instance in (tqdm(range(EVLValidNumInstpJob), desc="Compute Instances") if EDEBUG else range(EVLValidNumInstpJob)):
                res_opt = Algorithm.Get_objective_function_value(f"{Algorithm.ALGResFile}/donnees_{nb_jobs}_{index_instance+1}.dat.seq")
                res_ref = Algorithm.Get_objective_function_value(f"{EVLValRD}/donnees_{nb_jobs}_{index_instance+1}.dat.seq")
                max_UB = Algorithm.Get_UpperBound(EVLValDD,nb_jobs,index_instance)
                dev_res[indexLoopStatsPerJobs].append((max_UB-res_ref)/max_UB*100)
                dev_opt[indexLoopStatsPerJobs].append((max_UB-res_opt)/max_UB*100)
                times_opt[indexLoopStatsPerJobs].append(Algorithm.Get_execution_time(f"{Algorithm.ALGResFile}/donnees_{nb_jobs}_{index_instance+1}.dat.seq"))
                times_res[indexLoopStatsPerJobs].append(Algorithm.Get_execution_time(f"{EVLValRD}/donnees_{nb_jobs}_{index_instance+1}.dat.seq"))
            indexLoopStatsPerJobs+=1

    # save our results in a csv file
    with open("validation_results.csv", "w") as f:
        f.write(";Baseline;;;;;Learned;;;;\n")
        f.write("n;t_min;t_moy;t_max;dev min (%);dev moy(%);t_min;t_moy;t_max;dev min (%);dev moy(%)\n")
        nb_jobs=EVLMinNumJob
        indexLoopStatsPerJob = 0
        while nb_jobs<=EVLMaxNumJob:
            refTimesPerJobs = times_res[indexLoopStatsPerJob] #list reference times per jobs
            optTimesPerJobs = times_res[indexLoopStatsPerJob] #list optimized times per jobs
            refDevPerJobs = dev_res[indexLoopStatsPerJob] #list reference deviation per jobs
            optDevPerJobs = dev_opt[indexLoopStatsPerJob] #list optimized deviation per jobs
            f.write(f"{nb_jobs};{min(refTimesPerJobs)};{np.mean(refTimesPerJobs)};{max(refTimesPerJobs)};{min(refDevPerJobs)};{np.mean(refDevPerJobs)};{min(optTimesPerJobs)};{np.mean(optTimesPerJobs)};{max(optTimesPerJobs)};{min(optDevPerJobs)};{np.mean(optDevPerJobs)}\n")
            indexLoopStatsPerJob += 1
            nb_jobs+=EVLStepNumJob

    

