import argparse
import json
import pathlib
import math

def main(args):

    """
    :param args: all arguments given to the python script
    :return: Nothing, generate all configuration file for solve or generate instances
    """

    # Initialize parameters from command-line arguments

    """
    GENERATION INSTANCE
    """
    # characteristics of instance
    noIdenticalJobs = False if args.no_identical_jobs == 0 else True
    TF = args.tf
    RDD = args.rdd
    NList = args.n_list
    fracOfN = args.frac_of_n
    M = [tuple(map(int, pair.split(','))) for pair in args.m]
    Vmax = args.vmax
    V0 = args.v0
    # seed and nb instance to generate
    nbInstanceToGenerate = args.nb_instance_to_generate
    seed = args.seed
    # paths
    pathSaveInstance = pathlib.Path(str(pathlib.Path(__file__).parent.resolve()) + "/.."+ args.path_save_instance ).resolve()
    configFileNameGenerate = args.config_file_name_generate

    """
    GENERAL SOLVER PARAM
    """
    verbose = args.verbose
    outputResult = pathlib.Path(str(pathlib.Path(__file__).parent.resolve()) + "/.."+ args.output_result).resolve()
    configFileNameSolve = args.config_file_name_solve
    methods = args.methods
    timeLimit = int(args.timeLimit)

    """
    BRANCH AND BOUND PARAM
    """
    useHeuristic = args.use_heuristic_gen_col
    maxNbCallHeuristic = args.maxNbCallHeuristic
    nbMinStateDP = args.nbMinStateDP

    """
    BEAM SEARCH PARAM
    """
    beamSize = args.beamSize
    autoSetting = args.autoSetting == 1
    recovering = args.recovering
    alpha = args.alpha
    recoStrategy = args.strategy
    nbSolutionForMSLS = args.nbSolutionForMSLS
    versionRBS = int(args.versionRBS)

    """
    LOCAL SEARCH PARAM
    """
    bestImproveLS = args.bestImproveLS
    versionLS = int(args.versionLS)
    timeLimitLS = int(args.timeLimitLS)

    """
    Generate instances config 
    """
    # list of names for each instance
    listNameInstances = []

    # JSON
    configGenerate = {}

    generate = {
        "seed": seed
    }

    # the list of configuration of each instance
    listConfigInstance = []
    # define the config files for generate instance
    for (mMax, M0) in M:
        for N in NList:
            for n in [int(N * x) for x in fracOfN]:
                for tf in TF:
                    for rdd in RDD:
                        # Keep the list of instances for make the solve config file
                        instancePath = "_n_" + str(n) + "_N_" + str(N) + "_tf_" + str(tf) + "_rdd_" + str(
                            rdd) + "_mMax_" + str(mMax) + "_m0_" + str(M0) + "_Vmax_" + str(Vmax) + "_V0_" + str(V0)
                        listNameInstances.append(instancePath)

                        instance = {
                            "basePath": str(pathSaveInstance)+ "/",
                            "numberInstance": nbInstanceToGenerate}
                        param = {
                                "noIdenticalJobs": noIdenticalJobs ,
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
    configGenerate["generate"] = generate
    if configFileNameGenerate is not None:
        configPath = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameGenerate}.json").resolve()
        configPath.parent.mkdir(parents=True, exist_ok=True)
        with open(configPath, 'w') as file:
            file.write(json.dumps(configGenerate))
        print(f"Generated config files:\n- {configPath}")

    def clip(x,a,b):
        return min(max(a,x),b)

    def norm(x,min_x,max_x):
        return (x-min_x)/(max_x-min_x)

    def w(N,m,alpha1,alpha2,alpha3,alpha4):
        return clip(math.floor(math.exp(alpha1*norm(N,40,100)+alpha2)+alpha3*norm(m,2,10))+alpha4,1,100)
    def a(N,n,m,alpha5,alpha6,alpha7,alpha8):
        return clip(alpha5*norm(N,40,100)+alpha6*norm(n,10,75)+alpha7*norm(m,2,10)+alpha8,0.0,1.0)

    """
    Generate solver config file 
    """
    if configFileNameSolve is not None:
        # Add methods based on user selection
        for method in methods:
            for x in listNameInstances:
                params = x.split("_")[1:]
                M0 = int(params[params.index('m0') + 1])
                mMax = int(params[params.index('mMax') + 1])
                N = int(params[params.index('N') + 1])
                n = int(params[params.index('n') + 1])
                for i in range(nbInstanceToGenerate):
                    # Solve instances config
                    configSolve = {"solve": {"verbose": verbose, "output": str(outputResult), "methods": []}}
                    if method == "CG":
                        configSolve["solve"]["methods"].append({
                            "name": "CG",
                            "verbose": verbose,
                            "debug": True,
                            "gen_columns": 0,
                            "nbMinStateDP": nbMinStateDP,
                            "thresholdSetCol": 0.8,
                            "nbTimeNotUsed": 20,
                            "timeLimits": timeLimit,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "BeamSearch":
                        m = (mMax + M0)
                        # Normalization on [0,1]

                        W = w(N, m, -1.94010, 2.55139, -7.58156, 3) if nbSolutionForMSLS > 0 else w(N,m,-2.52765,3.03010,-1.38515,0.0)
                        # If W must be an integer
                        W = int(round(W))
                        alpha = a(N, n,m,1.0,-1.0,-1.0,-1.0) if (autoSetting and nbSolutionForMSLS==0) else alpha
                        TLS = math.ceil( 0.524*timeLimit)
                        K = 989
                        configSolve["solve"]["methods"].append({
                            "name": "BeamSearch",
                            "verbose": verbose,
                            "beamSize": W if autoSetting else beamSize,
                            "timeLimits": timeLimit - TLS if (nbSolutionForMSLS > 0 and autoSetting) else timeLimit,
                            "recovering": recovering==1,
                            "reco-strategy": recoStrategy,
                            "nbSolutionForMSLS": K if (nbSolutionForMSLS > 0 and autoSetting) else nbSolutionForMSLS,
                            "version": versionRBS,
                            "LocalSearch_parameters": {
                                "name": "LocalSearch",
                                "version": versionLS,
                                "useBestNeighbor": bestImproveLS==1,
                                "timeLimits": TLS if autoSetting else timeLimitLS,
                                "maxIter": 100
                            },
                            "alpha": alpha,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })
                    elif method == "LocalSearch":
                        configSolve["solve"]["methods"].append({
                            "name": "LocalSearch",
                            "verbose": verbose,
                            "timeLimits": timeLimit,
                            "maxIter": 100,
                            "predictor": (str(pathlib.Path(__file__).parent.resolve())+"/Predictor/predictor.pt"),
                            "version": versionLS,
                            "useBestNeighbor": bestImproveLS==1,
                            "instances": [
                                {"path": f"{pathSaveInstance}/instance{i}{x}.txt"}
                            ]
                        })


                    solveConfigPath = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameSolve}i{i}{x}.json").resolve()
                    solveConfigPath.parent.mkdir(parents=True, exist_ok=True)
                    with open(solveConfigPath, "w") as file:
                        file.write(json.dumps(configSolve))

        dirOfConfiFiles = pathlib.Path(str(pathlib.Path(str(pathlib.Path(__file__).parent.resolve()))) + f"/../instances/configs/{configFileNameSolve}").resolve()
        if len(methods) > 0:
            print(f"- {dirOfConfiFiles}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and solve scheduling instances.")

    # Add arguments
    parser.add_argument("--no-identical-jobs",type=int, default=0, help="Disable identical jobs, 0 for disable 1 for enable. By default, identical jobs are generated." )
    parser.add_argument("--tf", type=float, nargs='+', default=[0.2, 0.4, 0.6, 0.8, 1.0], help="TF values.")
    parser.add_argument("--rdd", type=float, nargs='+', default=[0.2, 0.4, 0.6, 0.8, 1.0], help="RDD values.")
    parser.add_argument("--n-list", type=int, nargs='+', default=[20], help="List of N values.")
    parser.add_argument("--frac-of-n", type=float, nargs='+', default=[0.5], help="Fraction of N jobs.")
    parser.add_argument("--m", type=str, nargs='+', default=["1,1", "2,2"], help="List of M values as 'mMax,M0'.")
    parser.add_argument("--vmax", type=float, default=2.0, help="Maximum velocity (Vmax).")
    parser.add_argument("--timeLimit", type=float, default=300, help="Maximum Time Limit.")
    parser.add_argument("--timeLimitLS", type=float, default=300, help="Maximum Time Limit for local search in MSLS.")
    parser.add_argument("--v0", type=float, default=1.0, help="Initial velocity (V0).")
    parser.add_argument("--nb-instance-to-generate", type=int, default=5, help="Number of instances to generate.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed.")
    parser.add_argument("--path-save-instance", type=str, default= "/instances/default/", help="Path to save instances. By default it's in directory ./instances/")
    parser.add_argument("--config-file-name-generate", type=str, default=None, help="Name of the generate config file.")
    parser.add_argument("--verbose", type=int, default=0, help="Verbosity level.")
    parser.add_argument("--output-result", type=str, default= "/instances/default/", help="Path to save results. By default it's in directory ./instances/default/")
    parser.add_argument("--config-file-name-solve", type=str, default=None, help="Name of the solve config file. If not set, then we don't generate solver configuration file")
    parser.add_argument("--methods", type=str, nargs='+', choices=["MIP", "LP", "CG", "BaB_CG","BaB_MIP","LocalSearch","BeamSearch"], default=[], help="Methods to use for solving.")
    parser.add_argument('--strategy',type=str,choices=['best-insert','local-search'],default='local-search',help="The strategy to use for the recovering. Options: 'best-insert','local-search'. Defaults to 'local-search'. Only applicable if method is 'BaB'.")
    parser.add_argument('--use-heuristic-gen-col',type=int, default=0, help="Use an heuristic to generate columns in the Column Generation. By default, we don't use it." )
    parser.add_argument('--maxNbCallHeuristic',type=int, default=1, help="The maximum number of calls where the heuristic failed before use only dynamic programming . By default is 1" )
    parser.add_argument('--nbMinStateDP',type=int, default=1, help="The number of column that will be generated. By default is 1" )
    parser.add_argument('--beamSize',type=int, default=1, help="The beam size of the beamSearch algorithm. By default is 1" )
    parser.add_argument('--autoSetting',type=int, default=1, help="Enable the default setting mode, i.e. adapt beam size and other parameter in function of N,n,m and coefficient. When enabled, 'beamSize' parameter are ignored. By default is 1" )
    parser.add_argument('--versionRBS',type=int, default=1, help="The version of the recovering beam search to use" )
    parser.add_argument('--versionLS',type=int, default=1, help="The version of the local search to use" )
    parser.add_argument('--recovering',type=int, default=1, help="Diseable the recovering, 0 to disable 1 to enable. By default is 1" )
    parser.add_argument('--alpha',type=float, default=0.5, help="The coefficient in the evaluation use in the recovering with local search strategy")
    parser.add_argument("--nbSolutionForMSLS", type=int, default=10, help="Enable using local search at the end of the heuristic (as multi start local search) on the best solutions founds, the value correspond to the number of solution we have to keep. By default is 10" )
    parser.add_argument("--bestImproveLS", type=int, default=0, help="Enable using the best improvement strategy in local search. By default it's disabled" )
    args, leftovers = parser.parse_known_args()
    main(args)
