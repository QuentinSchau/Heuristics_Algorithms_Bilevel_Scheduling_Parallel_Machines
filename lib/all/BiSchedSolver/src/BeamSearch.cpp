// Copyright (C) 2024
// Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//
// DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
// This file is part of bilevel-scheduling.
//
// bilevel-scheduling is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// bilevel-scheduling is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

//
// Created by schau on 6/30/25.
//


#include "BeamSearch.h"

BeamSearch::BeamSearch() : globalUB(std::numeric_limits<double>::max()), time_elapsed_in_LS(0), nbNodeLoc(0) {}

BeamSearch::BeamSearch(Instance *instance) : ISolver(instance), columnGeneration(instance), heuristicSolver(instance), localSearchSolver(instance), globalUB(instance->getSumWj() + 1), time_elapsed_in_LS(0), nbNodeLoc(0) {}

BeamSearch::BeamSearch(Instance *instance, nlohmann::json &object) :
        ISolver(instance), columnGeneration(instance), heuristicSolver(instance),localSearchSolver(instance),globalUB(instance->getSumWj()+1), time_elapsed_in_LS(0),nbNodeLoc(0) {
    if (object.contains("name")) {
        if (object["name"] == "BeamSearch") {
            // first set the column generation if we have parameters
            bool CG_haveParameters = false;
            bool LS_haveParameters = false;
            if (object.contains("LB_parameters")) {
                if (object["LB_parameters"].contains("name")) {
                    if (object["LB_parameters"]["name"].is_string()) {
                        if (object["LB_parameters"]["name"] == "CG") {
                            columnGeneration.setParameters(object["LB_parameters"]);
                            CG_haveParameters = true;
                        } else
                            throw std::invalid_argument(
                                    R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have 'name' that not corresponding to current implemented lower bound. We use the Column Generation "CG")");
                    } else throw std::invalid_argument(R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have 'name' attribute that is not a string)");
                } else throw std::invalid_argument(R"(The JSON object for setting lower bound ,i.e. 'LB_parameters', have not a name, check so documentation to set the lower bound)");
            }
            // check if we have configuration for local search solver
            if (object.contains("LocalSearch_parameters")) {
                if (object["LocalSearch_parameters"].contains("name")) {
                    if (object["LocalSearch_parameters"]["name"].is_string()) {
                        if (object["LocalSearch_parameters"]["name"] == "LocalSearch") {
                            localSearchSolver.setParameters(object["LocalSearch_parameters"]);
                            LS_haveParameters = true;
                        } else throw std::invalid_argument(R"(The JSON object for setting local search solver ,i.e. 'LocalSearch_parameters', have 'name' that not corresponding to current implemented local search solver. Check so documentation to set the Local Search JSON object)");
                    } else throw std::invalid_argument(R"(The JSON object for setting local search solver ,i.e. 'LocalSearch_parameters', have 'name' attribute that is not a string)");
                } else throw std::invalid_argument(R"(The JSON object for setting local search solver ,i.e. 'LocalSearch_parameters', have not a name, check so documentation to set the Local Search JSON object)");
            }
            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) {
                    setVerbose(object["verbose"].template get<char>());
                    if (!CG_haveParameters && verbose >= 1)
                        std::cout << "Using default LB, i.e., Column Generation" << std::endl;
                    if (!CG_haveParameters) columnGeneration.setVerbose(object["verbose"].template get<char>());
                    if (!LS_haveParameters && verbose >= 1)
                        std::cout << "Using default LocalSearch solver" << std::endl;
                    if (!LS_haveParameters) localSearchSolver.setVerbose(object["verbose"].template get<char>());
                } else
                    throw std::invalid_argument(
                            R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of BeamSearch solver)");
            }

            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) {
                    setTimeLimit(object["timeLimits"]);
                    columnGeneration.setTimeLimit(object["timeLimits"]);
                    if (!LS_haveParameters) localSearchSolver.setTimeLimit(object["timeLimits"]);
                } else
                    throw std::invalid_argument(
                            R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of BeamSearch solver)");
            }
            if (object.contains("beamSize")) {
                if (object["beamSize"].is_number_unsigned()) setBeamSize(object["beamSize"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "beamSize" of JSON object must be an "unsigned int" value for the constructor of BeamSearch solver)");
            }if (object.contains("version")) {
                if (object["version"].is_number_unsigned()) setVersion(object["version"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "version" of JSON object must be an "unsigned int" value for the constructor of BeamSearch solver)");
            }if (object.contains("alpha")) {
                if (object["alpha"].is_number_float()) setAlpha(object["alpha"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "alpha" of JSON object must be an "double" value for the constructor of BeamSearch solver)");
            }
            if (object.contains("recovering")) {
                if (object["recovering"].is_boolean()) setUseRecovering(object["recovering"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "recovering" of JSON object must be an "boolean" value for the constructor of BeamSearch solver)");
                if (object.contains("reco-strategy")) {
                    if (object["reco-strategy"].is_string()) setRecoStrategy(object["reco-strategy"]);
                    else
                        throw std::invalid_argument(
                                R"(The attribute "reco-strategy" of JSON object must be an "string" value for the constructor of BeamSearch solver)");
                }else setRecoStrategy("local-search");
            }if (object.contains("nbSolutionForMSLS")) {
                if (object["nbSolutionForMSLS"].is_number_unsigned()) {
                    setNbBestSolutionKeep(object["nbSolutionForMSLS"].template get<unsigned int>());
                }else
                    throw std::invalid_argument(
                            R"(The attribute "nbSolutionForMSLS" of JSON object must be an "unsigned int" value for the constructor of BeamSearch solver)");
            }
        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a BeamSearch solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a BeamSearch solver)");
}

BeamSearch::~BeamSearch() = default;

void BeamSearch::initialize() {
    //set timeUp to all other used solver
    columnGeneration.setTimeUp(timeUp);
    heuristicSolver.setTimeUp(timeUp);
    localSearchSolver.setTimeUp(timeUp);
    localSearchSolver.setTimeUpHeuristicSolver(timeUp);

    if (nbBestSolutionKeep>0) {
        // create the priority queue for the leaves solutions
        std::vector<std::pair<double,Solution>> containerLeaves;
        //use resize to silent false positive warning from GCC’s aggressive optimization.
        containerLeaves.resize(nbBestSolutionKeep+1); // add +1 to the desired size, because we will pop it
        containerLeaves.clear();
        listBestFoundLeaf = decltype(listBestFoundSolutions)(Compare<>(), std::move(containerLeaves));
        std::vector<std::pair<double,Solution>> containerSol;
        //use resize to silent false positive warning from GCC’s aggressive optimization.
        containerSol.resize(nbBestSolutionKeep+1); // add +1 to the desired size, because we will pop it
        containerSol.clear();
        listBestFoundSolutions = decltype(listBestFoundSolutions)(Compare<>(), std::move(containerSol));
    }
    listAvailableJobForNode.reserve(instance->getNbJobs());
    jobsInNewBlockStructure.reserve(instance->getNbJobs());

    // creating the root node
    if (verbose >= 2) std::cout << "Creating the root node" << std::endl;
    Node rootNode = Node(instance);
    //init the column generation if alpha > 0.0 or alpha == 0.0, and we use recovering with version 1
    if (isSmaller(0.0,alpha) || (useRecovering  && version == 1 )) {
        columnGeneration.generateStartingColumns(rootNode);
        columnGeneration.clearConstraintOfModel();
        columnGeneration.initializeModel(rootNode);
        columnGeneration.updateValueOfLmax(rootNode, 0);
        columnGeneration.updateValueOfLmax(rootNode, 1);
        columnGeneration.setStartTime(start);
    }

    #ifdef DEBUG_BeamSearch
    std::string name = "rootNode";
    rootNode.stateDebug.emplace_back(name);
    #endif

    isWithinTimeLimit();
    addNode(rootNode);
}

void BeamSearch::solve() {
    #ifdef DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    dot << "graph search {" << std::endl;
    #endif
    try {
        initialize();
        bool haveActiveNode = !heap.empty();
        std::vector<BeamSearchNode> nodeInBeam;
        nodeInBeam.reserve(beamSize);
        while (haveActiveNode && isSmaller(0.0,globalUB)) {
            isWithinTimeLimit();
            for (unsigned int indexLoopBeamSize = 0; indexLoopBeamSize < beamSize && !heap.empty(); indexLoopBeamSize++) {
                nodeInBeam.emplace_back(heap.top());
                heap.pop();
            }
            // clear the heap, and branch on all node in the beam
            heap = std::priority_queue<BeamSearchNode, std::vector<BeamSearchNode>>();
            for (auto & BSNode: nodeInBeam){
                isWithinTimeLimit();
                #if defined DEBUG_BaB && defined DEBUG_DOT
                dot << BSNode.node.id << DOT_BEAM << ";";
                #endif
                if (useRecovering) {
                    recovering(BSNode);
                    #if defined DEBUG_BaB && defined DEBUG_DOT
                    dot << BSNode.node.id << DOT_RECO << ";";
                    #endif
                }
                branchingLocation(BSNode); // generate children from the given node
            }
            nodeInBeam.clear(); //clear the node in the beam
            haveActiveNode = !heap.empty();
            //check the time limit
            isWithinTimeLimit();

        }
    } catch (const BiSchTimeOutException &e) {} //do nothing


    const auto end{std::chrono::steady_clock::now()};
    time_elapsed = std::chrono::duration<double>{end - start};
    solution->evaluate();
    if (verbose >= 1)
        std::cout << "Beam Search is over after " << time_elapsed.count() << " seconds" << std::endl
                  << "The objective value is : " << solution->getSumWjUj() << std::endl;
    if (verbose >= 2) {
        std::cout << "Nb nodes Loc: " << nbNodeLoc << std::endl;
        if (nbBestSolutionKeep > 0)
            std::cout << "Try to improve with local search" << std::endl;
    }
    #ifdef DEBUG_DOT
    dot << "}" << std::endl;
    #endif
    if (nbBestSolutionKeep > 0) {
        //keep the start time, we change it for the local search and will reset it after
        auto tempStartTime = start;
        auto tempTimeLimit = time_limits;
        const auto startLS{std::chrono::steady_clock::now()};
        //restart the start time and the time limit of the method
        resetTimer(localSearchSolver.getTimeLimits().count());
        try {
            static std::mt19937 g(0); // random generator (always the same to get reproducible result)
            if (verbose >= 1) std::cout << "Run local search on " << listBestFoundLeaf.size() << " leaves" << std::endl;
            NB_MULTI_START_SOLUTIONS_FROM_LEAF = listBestFoundLeaf.size();
            std::vector<std::pair<double,Solution>> &solutions = Container(listBestFoundLeaf);
            std::shuffle(solutions.begin(), solutions.end(),g);
            while (not solutions.empty()){
                auto [objValue,sol] = solutions.back();
                solutions.pop_back();
                auto blockStruct = sol.toBlockStruct(instance);
                localSearchSolver.setStartTime();
                localSearchSolver.localSearchBySwapAllJobs(&blockStruct,objValue);
                if (isSmaller(localSearchSolver.getSolution()->getSumWjUj(),objValue) && isSmaller(localSearchSolver.getSolution()->getSumWjUj(),globalUB)) {
                    if (verbose >= 2) std::cout << "Local search found better solution" << std::endl;
                    *solution = *localSearchSolver.getSolution();
                    if (solution->empty()) {
                        throw BiSchException("Erreur Local Search get empty solution");
                    }
                    globalUB = localSearchSolver.getSolution()->getSumWjUj();
                }
            }
            if (verbose >= 1) std::cout << "Run local search on " << std::min(static_cast<unsigned int>(listBestFoundSolutions.size()),nbBestSolutionKeep - NB_MULTI_START_SOLUTIONS_FROM_LEAF) << " solutions" << std::endl;
            NB_MULTI_START_SOLUTIONS_FROM_BEST_SOL = 0;
            std::shuffle(solutions.begin(), solutions.end(),g);
            solutions = Container(listBestFoundSolutions);
            while (not solutions.empty() && ++NB_MULTI_START_SOLUTIONS_FROM_BEST_SOL < nbBestSolutionKeep - NB_MULTI_START_SOLUTIONS_FROM_LEAF ){
                auto [objValue,sol] = solutions.back();
                solutions.pop_back();
                auto blockStruct = sol.toBlockStruct(instance);
                localSearchSolver.setStartTime();
                localSearchSolver.localSearchBySwapAllJobs(&blockStruct,objValue);
                if (isSmaller(localSearchSolver.getSolution()->getSumWjUj(),objValue) && isSmaller(localSearchSolver.getSolution()->getSumWjUj(),globalUB)) {
                    if (verbose >= 2) std::cout << "Local search found better solution" << std::endl;
                    *solution = *localSearchSolver.getSolution();
                    if (solution->empty()) {
                        throw BiSchException("Erreur Local Search get empty solution");
                    }
                    globalUB = localSearchSolver.getSolution()->getSumWjUj();
                }
            }
        } catch (const BiSchTimeOutException &e) {} //do nothing

        //restore time limit and start time of the Beam Search Solver
        start=tempStartTime;
        time_limits=tempTimeLimit;
        const auto endLS{std::chrono::steady_clock::now()};
        time_elapsed_in_LS = std::chrono::duration<double>{endLS - startLS};
        // if (verbose >= 2)
        //     std::cout << "Local Search is over after " << time_elapsed_in_LS.count() << " seconds" << std::endl;
        time_elapsed = std::chrono::duration<double>{endLS - start};
    }
    solution->evaluate();
    if (nbBestSolutionKeep > 0 && verbose >= 2)
        std::cout << "Multi Start Local Search is over after " << time_elapsed.count() << " seconds" << std::endl
        << "The objective value is : " << solution->getSumWjUj() << std::endl;
}

void BeamSearch::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists)
        outputFile << "InstanceName"
                   << "\t" << "InstancePath"
                   << "\t" << "N"
                   << "\t" << "n"
                   << "\t" << "m_Max"
                   << "\t" << "m_0"
                   << "\t" << "V_max"
                   << "\t" << "V_0"
                   << "\t" << "Method"
                   << "\t" << "Time"
                   << "\t" << "TimeInLS"
                   << "\t" << "LimitTime"
                   << "\t" << "First_UB"
                   << "\t" << "BeamSize"
                   << "\t" << "Alpha"
                   << "\t" << "UseRecovering"
                   << "\t" << "RecoStrategy"
                   << "\t" << "VersionRBS"
                   << "\t" << "VersionLS"
                   << "\t" << "UseBestImprovement"
                   << "\t" << "NBCallsReco"
                   << "\t" << "NB_Best_NBH_FindWithReco"
                   << "\t" << "NB_MSLS_SOLUTIONS"
                   << "\t" << "NB_MSLS_BEST_SOL"
                   << "\t" << "NB_MSLS_LEAF_SOL"
                   << "\t" << "NBNodes"
                   << "\t" << "NBMinColum"
                   << "\t" << "GenerateCol"
                   << "\t" << "MaxNbCallHeuristics"
                   << "\t" << "NBComputeCost"
                   << "\t" << "NBCallsHeu"
                   << "\t" << "NBCallDP"
                   << "\t" << "NBCallSubProcessCG"
                   << "\t" << "Objective" << std::endl;

    solution->evaluate();
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "BeamSearch"
               << "\t" << time_elapsed.count()
               << "\t" << (nbBestSolutionKeep >0 ? time_elapsed_in_LS.count() : 0.0)
               << "\t" << time_limits.count()
               << "\t" << firstUB
               << "\t" << beamSize
               << "\t" << alpha
               << "\t" << (useRecovering ? 1u : 0u)
               << "\t" << (useRecovering ? getRecoStrategy() : "null")
               << "\t" << (useRecovering ? std::to_string(version).c_str() : "null")
               << "\t" << (nbBestSolutionKeep >0 ? std::to_string(localSearchSolver.getVersion()).c_str() : "null")
               << "\t" << (nbBestSolutionKeep >0 ? std::to_string(localSearchSolver.isUseBestNeighbor()).c_str() : "null")
               << "\t" << NB_reco
               << "\t" << NB_best_sol_find_reco
               << "\t" << nbBestSolutionKeep
               << "\t" << NB_MULTI_START_SOLUTIONS_FROM_BEST_SOL
               << "\t" << NB_MULTI_START_SOLUTIONS_FROM_LEAF
               << "\t" << nbNodeLoc
               << "\t" << columnGeneration.getNbMinStateDp()
               << "\t" << static_cast<unsigned int>(columnGeneration.getGenerateColumn())
               << "\t" << columnGeneration.getMaxNbCallHeuristic()
               << "\t" << columnGeneration.getNbCallComputeCost()
               << "\t" << columnGeneration.getNbCallsHeu()
               << "\t" << columnGeneration.getNbCallsDp()
               << "\t" << columnGeneration.getNbCallSubProcessCG()
               << "\t" << solution->getSumWjUj() << std::endl;
    outputFile.close();
}
