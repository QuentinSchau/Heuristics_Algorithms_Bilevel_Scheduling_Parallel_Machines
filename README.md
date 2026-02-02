# Bilevel-Scheduling-I4.0

## Description
This project implements various algorithms to solve bilevel scheduling problems. All the algorithms are developed as part of Quentin SCHAU's thesis, conducted in collaboration between LIFAT (Tours, France) and Politecnico di Torino (Torino, Italy).
This repo is based on the heuristics algorithms developed by Quentin Schau (https://github.com/QuentinSchau/Heuristics_Algorithms_Bilevel_Scheduling_Parallel_Machines).
It is based on the repo https://github.com/QuentinSchau/Exact_Algorithms_Bilevel_Scheduling_Industry_Futur.

## Installation

### Dependencies

* Download and install **nlohmann/json** version 3.11.3 from GitHub: https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.zip
  + Place the downloaded archive in `./lib/all/BiSchedSolver/lib/all` directory.
  + Unzip the archive using the command `unzip v3.11.3.zip`.
* Download and install **Eigen** version 3.4.0 from GitHub: https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
  + Place the downloaded archive in `./lib/all/BiSchedSolver/lib/all` directory.
  + Unzip the archive using the command `unzip eigen-3.4.0.zip`.
* Download and install **Boost** version 1.88.0 from Boost.org: https://archives.boost.io/release/1.88.0/source/boost_1_88_0.tar.bz2
  + Place the downloaded archive in `./lib/all/BiSchedSolver/lib/all` directory.
  + Decompress the archive using the command `tar -xjf boost_1_88_0.tar.bz2 -C .`.
  + You can use the tool bcp to only keep Dynamic Bitset (reduce the size of boost)
* Download and install **CPLEX Solver** from https://www.ibm.com/fr-fr/products/ilog-cplex-optimization-studio
  + Place the CPLEX library in `./lib/BiSchedSolver/lib/Linux/CPLEX` directory.
  + For Windows users, you will need to set up the CPLEX library in `./lib/BiSchedSolver/lib/Windows/CPLEX`. You may need to modify the `CPLEXlibConfig.cmake` file.
* (if you want to use predictor) Download and install **Libtorch** version 2.9.1 (c++ version of Pytorch) from https://pytorch.org/get-started/locally/:
  + Example download with Linux: `wget https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.9.1%2Bcpu.zip`
  + Place the downloaded archive in `./lib/all/BiSchedSolver/lib/Linux` directory.
  + Decompress the zip using the command `unzip libtorch-shared-with-deps-2.9.1+cpu.zip`.
  + Change the directory for Windows users

The project uses CMake for installation. To generate and build the project, follow these steps:

* Navigate to the project root directory.
* Run the following command to generate the project in release mode: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
* Run the following command to generate the project in debug mode: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
* Change into the `/build` directory and run the command `make`.

The compiled binary program will be located in either `/bin/release` or `/bin/debug`, depending on the build mode chosen.

**Note:** The project was not compiled on Windows, some error may appear. Please feel free to open an issue if you encounter any issues during the compilation process.

### Scripts

The project includes several Python and Bash scripts to facilitate generating configuration files, creating instances, and running experiments. You can use these scripts as-is or modify them to fit your needs.

For python installation:
1. We recommend to use conda for its installation cf. [the documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Create and use a virtual environment. `conda create -n BilevelScheduling -y`
3. You need to use the environment with : `conda activate BilevelScheduling`
4. You need to use the environment with : `conda install pip -y`
5. Install the dependencies `pip install -r requirements.txt`

Scripts:
* The `generateInstances.py` script generates instance config files.
* The `createInstance.sh` Bash script creates instances with customizable settings, allowing you to modify variables such as values.
* The `createConfigFiles.sh` Bash script creates configuration JSON file with customizable settings, allowing you to modify variables such as values.
* The `runInstances.sh` script runs instances locally.


## Specification

### Instances

An instance file must be in "txt" format and follow the specific pattern:
```
name:Instance Name
N:Integer # Number of jobs
M_max:Integer # Number of high speed machines
M_0:Integer # Number of low speed machines
V_max:Float # High-speed machine processing rate
V_0:Float # Low-speed machine processing rate
Jobs:
...
pj \t dj \t wj
...
```
In the job description, the processing time `pj`, due date `dj`, and weighted value `wj` must be separated by a `\t` character.

## Usage

To run the Bilevel-Scheduling-I4.0 project, execute the binary program located in either `./bin/release` or `./bin/debug` directory. You must provide a configuration file. For reference, you can find two example files: `example_config_generate.json` and `example_config_solve.json`.

### Configuration File (JSON)

The project uses JSON files as configuration files. Below, we describe the structure of these files and how to use them.

#### Generate

To generate instances, you'll need to provide a JSON configuration file with the following structure:
```
{
    "generate": {
        // Seed for generating instances
        "seed": <int>,
        
        // Parameters for generating one instance
        "instances": [
            {
                // Base path where to save instances. It creates the directory if it does not exist.
                "basePath": "<string>",
                
                // Number of instances to generate
                "numberInstance": <int>,
                
                "paramInstance": {
                    // Whether to generate no identical jobs (true or false)
                    "noIdenticalJobs": <boolean>,
                    
                    // The number of jobs that the leader has to select.
                    "n": <int>,
                    
                    // Total number of jobs in the instance
                    "N": <int>,
                    
                    // Number of high-speed machines
                    "M_max": <int>,
                    
                    // Number of low-speed machines
                    "M_0": <int>,
                    
                    // High-speed value
                    "V_max": <float>,
                    
                    // Low-speed value
                    "V_0": <float>,
                    
                    // Processing time range for each job
                    "pi": {
                        "inf": <int>, // Lower bound of processing time
                        "sup": <int>  // Upper bound of processing time
                    },
                    
                    // Due date range for each job
                    "di": {
                        "inf": <int>, // Lower bound of due date
                        "sup": <int>  // Upper bound of due date
                    },
                    
                    // Weight range for each job
                    "wi": {
                        "inf": <int>, // Lower bound of weight
                        "sup": <int>  // Upper bound of weight
                    }
                }
            }
        ]
    }
}
```
This configuration file specifies the parameters needed to generate instances, including the seed value, instance base path, and parameters for generating a single instance. The `paramInstance` object contains additional parameters that control the generation process, such as whether to generate no identical jobs, the number of leader's jobs, and the ranges for processing times, due dates, and weights.



Here's an updated version with the names of each method:

#### Solve (Heuristics Algorithms)

To solve instances, there are different parameters described below:

```
"solve": {
    // Level of verbose
    "verbose": <int>,
    // Path where to save the results. If directory does not exist, it will be created.
    "output": "<string>",
    // List of methods used to define parameters for each method. This is described below.
    "methods": [
        {
            // Name of the method
            "name": "<string> (LocalSearch, BeamSearch)",
            // Level of verbose
            "verbose": <int>,
            // Time limits for the method to solve an instance, in seconds
            "timeLimits": <int>,
            // List of instances to solve. Each object is composed of only one attribute:
            "instances": [
                {
                    // Path to the instance file to solve
                    "path": "<string>"
                }
            ]
        }
    ]
}
```

Each method has its own set of parameters, which are described below:

**BS: Beam Search**

* Optional:
  + `LB_parameters` (JSON object): Parameterize the lower bound using a JSON object that defines a lower bound, accepted values: "CG" cf bellow.
  + `LocalSearch_parameters` (JSON object): Parameterize the local search using a JSON object, accepted values: "LocalSearch" cf bellow. 
  + `beamSize` (int): The width of the beam (number of partial solutions kept at each level).
  + `version` (int): Defines the specific variant of the Beam Search algorithm to execute:
    - 1: try to use exact recovering
    - 2: use heuristic recovering, try to improve the complementary schedule found in the evaluation phase.
    - 3: use heuristic recovering, by keeping the complementary schedule computed in evaluation phase. 
  + `alpha` (float): A weight parameter used during node evaluation.
  + `recovering` (boolean): If True, enables the recovering procedure to improve partial solutions.
  + `reco-strategy` (string): The strategy used for recovering. Can you "local-search" or "best-insert". Default: "local-search"
  + `nbSolutionForMSLS` (int): Number of best solutions to keep for Multi-Start Local Search (MSLS). If set to 0, the method skips the local search phase and operates purely as a Recovering Beam Search (RBS)

**LS: Local Search**

* Optional:
  + `maxIter` (int): Maximum number of iterations allowed before the search terminates.
  + `useBestNeighbor` (boolean): If True, the solver will move to the best neighbor found; if False, it may use a different selection strategy (e.g., first improvement). 
  + `version` (int): Define which specific variant of the local search algorithm to execute:
    - 1, 4, 5: Correspond to $LS_s$, $LS_{fa}$, and $LS_a$ as described in the paper.
    - 2, 3: Variants using the leader neighborhood to avoid radical changes to the solution.
    - 6, 7, 8, 9: Variants that integrate a predictor to guide the search.
    - 10: Random baseline version for performance comparison.
  + `genDatabase` (boolean): If True, the solver will generate and save data during the search (typically for training future models).
  + `ratioNeighbor` (float): Define the proportion of the neighborhood to be explored (between 0 and 1). Use only in version 9.
  + `nbPrediction` (int): Number of predictions to be made (used in conjunction with the predictor model). 
  + `batchSize` (int): Size of the batch used for processing neighbors or predictions.
  + `predictor` (string): Path to the predictor model file. Providing this automatically enables the use of a predictor.
  
  
**CG: Column Generation**
This method is used as lower bound in the beam search
* Optional:
  + `Debug` (boolean): Display more information when True.
  + `nbMinStateDP` (int): Define the number of columns generated for each machine schedule in the pricing problem.
  + `gen_columns` (int): Define which method to use to generate columns.
    0 -> Dynamic programming (DP) is used to solve the pricing problem
    1 -> Use heuristic to compute some columns before calling DP.
    ... -> You can define your own heuristic to generate some columns or other methods.
  + `maxNbCallHeuristic` (int): Define the number of call of the heuristic that failed before switch to only use DP.
  + `thresholdSetCol` (float): Define the threshold used when removing some columns, i.e., if the number of unfeasible
    columns over all columns is greater than this threshold then we apply a cleaning
    procedure.
  + `nbTimeNotUsed` (int): Define the number of times that we allow a column to be not used before being removed with the cleaning procedure. (Linked with the threshold defined above)
  + `timeLimitsMS` (int): time limit in milli seconds.

## Advanced

To add more debugging information, follow these steps:

* If you want to add more debugging information in some methods, you can uncomment the flag "DEBUG_X" in file "X.h". This will allow you to see additional debug messages when running the program. Need a recompilation of the program.
* To enable the debug mode, provide the "debug" parameter. This will prompt you with more detailed information about what's happening during execution.

All instances use to compare methods was generated with seed 0. All instance use to generate train database for the predictor or bayesian optimization use seed 1.
## Contributing

The main contributor is Quentin SCHAU. If you want to contribute to this project, you should reach out to Quentin SCHAU at quentin.schau@univ-tours.fr or quentin.schau@polito.it .

## Authors and acknowledgment

This code implements the problem and exact algorithm defined in several communication : 

## License
This project is under GPU Licence. 

