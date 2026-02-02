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

//
// Created by schau on 1/22/25.
//

#include "LocalSearch.h"

#include "BeamSearch.h"

LocalSearch::LocalSearch() = default;

LocalSearch::LocalSearch(Instance* instance) : ISolver(instance), heuristicSolver(instance) { initializeStructure(); }

LocalSearch::LocalSearch(Instance* instance, nlohmann::json& object) : ISolver(instance), heuristicSolver(instance) {
    setParameters(object);
    initializeStructure();
}

void LocalSearch::initializeStructure() {
    auto& E = instance->getE();
    std::vector<double> completionTimeForEachBlock(E.size(), 0.0);
    // for each type of machines and for each block we compute the min, max and average completion time of the processing time in the block
    minMaxAvgCompletionTime = std::vector<std::vector<double>>(6, completionTimeForEachBlock);
    // for each block, we compute the min, max, avg and std of the wj Uj for jobs that can be scheduled in the block. We use the min,max and average completion time of the block
    minMaxAvgStdForEachBlock = std::vector<std::vector<double>>(4, completionTimeForEachBlock);
}

LocalSearch::~LocalSearch() = default;

void LocalSearch::setParameters(nlohmann::json& object) {
    if (object.contains("name")) {
        if (object["name"] == "LocalSearch") {
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else throw std::invalid_argument(R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of LocalSearch Solver)");
            }
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of LocalSearch Solver)");
            }if (object.contains("version")) {
                if (object["version"].is_number_unsigned()) setVersion(object["version"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "version" of JSON object must be an "unsigned integer" value for the constructor of LocalSearch Solver)");
            }
            if (object.contains("maxIter")) {
                if (object["maxIter"].is_number_unsigned()) setMaxIter(object["maxIter"].template get<unsigned int>());
                else throw std::invalid_argument(R"(The attribute "maxIter" of JSON object must be an "unsigned integer" value for the constructor of LocalSearch Solver)");
            }if (object.contains("nbPrediction")) {
                if (object["nbPrediction"].is_number_unsigned()) setNbPrediction(object["nbPrediction"].template get<unsigned int>());
                else throw std::invalid_argument(R"(The attribute "nbPrediction" of JSON object must be an "unsigned integer" value for the constructor of LocalSearch Solver)");
            }if (object.contains("batchSize")) {
                if (object["batchSize"].is_number_unsigned()) setBatchSize(object["batchSize"].template get<unsigned int>());
                else throw std::invalid_argument(R"(The attribute "batchSize" of JSON object must be an "unsigned integer" value for the constructor of LocalSearch Solver)");
            }
            if (object.contains("predictor")) {
                if (object["predictor"].is_string()) {
                    pathPredictor = std::filesystem::path(object["predictor"].template get<std::string>());
                    setUsePredictor(true);
                }
                else throw std::invalid_argument(R"(The attribute "predictor" of JSON object must be an "string" value to specify the predictor model for the constructor of LocalSearch Solver)");
            }
            if (object.contains("useBestNeighbor")) {
                if (object["useBestNeighbor"].is_boolean()) setUseBestNeighbor(object["useBestNeighbor"].template get<bool>());
                else throw std::invalid_argument(R"(The attribute "useBestNeighbor" of JSON object must be an "boolean" value for the constructor of LocalSearch Solver)");
            }if (object.contains("ratioNeighbor")) {
                if (object["ratioNeighbor"].is_number_float()) setRatioNeighbor(object["ratioNeighbor"].template get<double>());
                else throw std::invalid_argument(R"(The attribute "ratioNeighbor" of JSON object must be an "float" value for the constructor of LocalSearch Solver)");
            }
            if (object.contains("genDatabase")) {
                if (object["genDatabase"].is_boolean()) setGenDatabase(object["genDatabase"].template get<bool>());
                else throw std::invalid_argument(R"(The attribute "genDatabase" of JSON object must be an "boolean" value for the constructor of LocalSearch Solver)");
            }
        }
        else throw std::invalid_argument(R"(The JSON object have not the right name to instance a LocalSearch Solver)");
    }
    else throw std::invalid_argument(R"(The JSON object have not a name to instance a LocalSearch Solver)");
}

void LocalSearch::swapJobsInSelection(std::vector<unsigned int>& listIndexOfJobs, unsigned int indexJobToInsert, unsigned int indexJobToRemove) {
    if (indexJobToInsert == indexJobToRemove)
        throw BiSchException("Error in function that swaps jobs in leader selection, try to swap same job");
    // if we want to insert smaller job
    if (indexJobToInsert < indexJobToRemove) {
        auto itPosToInsert = std::lower_bound(listIndexOfJobs.begin(), listIndexOfJobs.end(), indexJobToInsert);
        auto itPosToRemove = std::lower_bound(listIndexOfJobs.begin(), listIndexOfJobs.end(), indexJobToRemove);
        assert(itPosToRemove != listIndexOfJobs.end());
        auto itRotateJobs = std::rotate(itPosToInsert, itPosToRemove, std::next(itPosToRemove));
        if (itRotateJobs != listIndexOfJobs.begin()) --itRotateJobs;
        *itRotateJobs = indexJobToInsert;
    }
    else {
        // insert greater job
        auto itPosToInsert = std::lower_bound(listIndexOfJobs.begin(), listIndexOfJobs.end(), indexJobToInsert);
        auto itPosToRemove = std::lower_bound(listIndexOfJobs.begin(), listIndexOfJobs.end(), indexJobToRemove);
        assert(itPosToRemove != listIndexOfJobs.end());
        auto itNext = std::next(itPosToRemove);
        if (itNext == listIndexOfJobs.end() || itPosToRemove == itPosToInsert) {
            // the job to delete is already the last one
            *itPosToRemove = indexJobToInsert; // change it
        } else {
            auto itRotateJobs = std::rotate(itPosToRemove, itNext, itPosToInsert);
            *itRotateJobs = indexJobToInsert;
        }
    }
    assert(std::is_sorted(listIndexOfJobs.begin(),listIndexOfJobs.end()));
    assert(listIndexOfJobs.size() == static_cast<size_t>(std::distance(listIndexOfJobs.begin(),std::unique(listIndexOfJobs.begin(),listIndexOfJobs.end()))));
}

std::pair<Solution::BlockStructure, double> LocalSearch::updateBlockStructureWithSwap(std::vector<unsigned int>& listIndexOfJobs,Solution::BlockStructure& blockStructure, unsigned int indexJobToInsert, unsigned int indexJobToRemove) {
    bool blockHaveToChange = false;
    unsigned int indexBlock = 0;
    unsigned int indexLoopJob = 0;
    unsigned int minIndexJob = std::min(indexJobToInsert, indexJobToRemove);
    auto& E = instance->getE();
    auto& listOfJobs = instance->getListJobs();
    unsigned int nbJobToSelectOnBlocks = 0;
    while (indexBlock < E.size() && not blockHaveToChange) {
        unsigned int nbJobToScheduleInBlock = indexBlock == 0
                                                  ? instance->getNbJobsToScheduleOnFirstBlock()
                                                  : E[indexBlock].size();
        unsigned int indexLoopInBlock = 0;
        while (indexLoopInBlock < nbJobToScheduleInBlock && not blockHaveToChange) {
            if (not isSmaller(listOfJobs[listIndexOfJobs[indexLoopJob + indexLoopInBlock]].getPi(), listOfJobs[minIndexJob].getPi())) blockHaveToChange = true;
            ++indexLoopInBlock;
        }
        if (not blockHaveToChange) {
            indexLoopJob += nbJobToScheduleInBlock;
            nbJobToSelectOnBlocks += nbJobToScheduleInBlock;
            indexBlock++;
        }
    }
    Solution::BlockStructure newBlockStructure = blockStructure;
    if (indexJobToInsert == 495 && indexJobToRemove == 494)
        std::cout << "";

    swapJobsInSelection(listIndexOfJobs, indexJobToInsert, indexJobToRemove);
    // fill the rest with the block structure
    for (; indexBlock < E.size(); indexBlock++) {
        //get the number of job to schedule in the block
        nbJobToSelectOnBlocks += indexBlock == 0 ? instance->getNbJobsToScheduleOnFirstBlock() : E[indexBlock].size();
        for (auto& [indexMachine,indexInMachine] : E[indexBlock]) {
            // if there is a job scheduled in the block structure then update it
            if (blockStructure[indexMachine][indexInMachine].first != nullptr) {
                assert(indexLoopJob < listIndexOfJobs.size());
                const Job* job = &listOfJobs[listIndexOfJobs[indexLoopJob]];
                double completionTime = indexInMachine == 0 ? 0.0 : newBlockStructure[indexMachine][indexInMachine - 1].second;
                double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
                completionTime += job->getPi() / speed;
                newBlockStructure[indexMachine][indexInMachine].first = job;
                newBlockStructure[indexMachine][indexInMachine].second = completionTime;
                ++indexLoopJob; // pass next job
                if (nbJobToSelectOnBlocks == indexLoopJob) break; // if we have schedule all job in the block stop
            }
        }
    }

    //manage the identical jobs
    // use a map to store, the group of identical job and their available location thanks to a TreeCj. We can solve optimally the sub case of identical jobs
    std::map<unsigned int, std::pair<std::vector<Job>,TreeCj>> groupIdenticalJobsWithTheirAvailableLocation;
    // get from the list of selected jobs, then identical ones
    for (unsigned int indexJob : listIndexOfJobs) {
        unsigned int indexGroupIdJobs = instance->getIndexIdenticalGroupOfJob(indexJob);
        if (instance->getListGrpJobs()[indexGroupIdJobs].size() >= 2) {
            auto itGroupIdJob = groupIdenticalJobsWithTheirAvailableLocation.find(indexGroupIdJobs);
            if (itGroupIdJob != groupIdenticalJobsWithTheirAvailableLocation.end()) {
                itGroupIdJob->second.first.push_back(listOfJobs[indexJob]);
            }else {
                groupIdenticalJobsWithTheirAvailableLocation.insert({indexGroupIdJobs,{{listOfJobs[indexJob]},TreeCj()}});
            }
        }
    }
    // loop over the block structure and add the location corresponding to a idencal job
    for (indexBlock = 0; indexBlock < E.size(); ++indexBlock) {
        for (auto &[indexMachine,indexInMachine] : E[indexBlock]) {
            // if there is a job
            auto &jobWithCj = newBlockStructure[indexMachine][indexInMachine];
            if (jobWithCj.first != nullptr) {
                unsigned int indexGroupIdJobs = instance->getIndexIdenticalGroupOfJob(jobWithCj.first->getIndex());
                if (instance->getListGrpJobs()[indexGroupIdJobs].size() >= 2) {
                    auto itGroupIdJob = groupIdenticalJobsWithTheirAvailableLocation.find(indexGroupIdJobs);
                    assert(itGroupIdJob != groupIdenticalJobsWithTheirAvailableLocation.end());
                    itGroupIdJob->second.second.insert({jobWithCj.second,{indexMachine,indexInMachine}});
                }
            }
        }
    }
    // solve all sub problem with job and available location
    for (auto &[_,groupIdJobsAndLocation] : groupIdenticalJobsWithTheirAvailableLocation) {
        auto &[listOfIdenticalJobs, listCjAndAvailablePosition] = groupIdJobsAndLocation;
        assert(listOfIdenticalJobs.size() == listCjAndAvailablePosition.size());
        solveProblemWithFixedCompletionTime(nullptr, &newBlockStructure, listCjAndAvailablePosition, listOfIdenticalJobs);
    }

    return {newBlockStructure, Solution::evaluate(newBlockStructure,instance)};
}

bool LocalSearch::improveBlockStructureLocallyBySwap(Solution::BlockStructure& blockStructure, std::vector<unsigned int>& listScheduledJobs,std::vector<unsigned int> * listNotScheduledJobs, std::vector<bool> &alreadySelectedJobs, double& bestObjValue,bool solveIdenticalJobs) {
    // assume block structure is sorted by completion time increasing
    // we keep the best permutation in a block B_k, apply it then continue to neighborhoods on this best permutation
    unsigned int indexBlock = 0;
    bool foundBetterNeighbor = false;
    double UB = std::numeric_limits<double>::infinity();
    while (indexBlock < instance->getE().size()) {
        isWithinTimeLimit();
        // use neighborhood on the follower decisions
        Neighborhoods neighborhoods(instance, &blockStructure, &listScheduledJobs, 2);
        auto itSol = neighborhoods.beginFN();
        auto itEndNeighbor = neighborhoods.endFN();
        itSol.setIndexBlock(indexBlock); // define the current index of block
        itSol.initialize(); // initialize the iterator with the new index block
        for (; itSol != itEndNeighbor; ++itSol) {
            isWithinTimeLimit();
            // check if we have not pass to the next block with the neighborhood
            if (indexBlock < itSol.getIndexBlock()) { break; }
            Solution::BlockStructure blockStructureFromNeighbor = *itSol;
            #ifdef DEBUG_HEURISTIC
            Solution testSol(instance);
            testSol.fromBlockStruct(blockStructureFromNeighbor);
            assert(testSol.feasible(instance));
            #endif
            UB = Solution::evaluate(blockStructureFromNeighbor, instance); // evaluate the neighbor
            // if the permutation reduce the weighted number of tardy jobs
            if (isSmaller(UB, bestObjValue)) {
                blockStructure = std::move(blockStructureFromNeighbor);
                bestObjValue = UB;
                foundBetterNeighbor = true;
                if (not useBestNeighbor) break; // stop if we don't use the best improvement strategy
            }
            else if (isEqual(UB,bestObjValue)) {
                Solution::sortBlockStructurePerMakespan(blockStructureFromNeighbor,instance);
                bool isBetter = true; // the neighbor is better if it has the minimal makespan

                unsigned int nbEgalMakespan = 0;
                for (auto &[indexMachine,indexInMachine] : instance->getE()[indexBlock]){
                    if (isSmaller(blockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second)) {
                        isBetter = false;
                    }else if (isEqual(blockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second))
                        ++nbEgalMakespan;
                }
                if (isBetter and nbEgalMakespan < instance->getE()[indexBlock].size()) {
                    foundBetterNeighbor = true;
                    blockStructure = blockStructureFromNeighbor;
                    if (not useBestNeighbor) break; // stop if we don't use the best improvement strategy
                }
            }
        }
        indexBlock++;
        if (solveIdenticalJobs) solveProblemWithFixedCompletionTime(blockStructure, indexBlock);
        if (UB = Solution::evaluate(blockStructure, instance); isSmaller(UB, bestObjValue)) {
            bestObjValue = UB;
            foundBetterNeighbor = true;
        }
    }
    if (solveIdenticalJobs) {
        listScheduledJobs.clear();
        // update the already selected jobs
        std::fill(alreadySelectedJobs.begin(), alreadySelectedJobs.end(), false);
        // update list of selected jobs because some of them must change with identical sub problem
        for (auto &machine : blockStructure){
            for (auto &jobWithCj : machine) {
                if (jobWithCj.first != nullptr) {
                    listScheduledJobs.push_back(jobWithCj.first->getIndex());
                    alreadySelectedJobs[jobWithCj.first->getIndex()].flip();
                }
            }
        }
        std::sort(listScheduledJobs.begin(),listScheduledJobs.end());
        if (listNotScheduledJobs != nullptr) {
            unsigned int indexLoopJobsNotSelected = 0;unsigned int indexLoopInListJobNotSelected=0;
            for ( bool isSelected : alreadySelectedJobs) {
                if (not isSelected) {
                    (*listNotScheduledJobs)[indexLoopInListJobNotSelected]=indexLoopJobsNotSelected;
                    indexLoopInListJobNotSelected++;
                }
                indexLoopJobsNotSelected++;
            }
        }
    }
    return foundBetterNeighbor;
}

bool LocalSearch::improveBlockStructureLocallyByAssignment(Solution::BlockStructure &blockStructure,std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs,double &bestObjValue, bool solveIdenticalJobs) {
    // assume block structure is sorted by completion time increasing
    // we keep the best permutation in a block B_k, apply it then continue to neighborhoods on this best permutation

    bool foundBetterNeighbor = false;
    double UB;
    bool foundBetterSol = true;
    unsigned int nbIterImprovBlockStruct = 0;
    while (foundBetterSol && nbIterImprovBlockStruct < maxIter) {
        isWithinTimeLimit();
        foundBetterSol = false;
        unsigned int indexBlock = 0;
        nbIterImprovBlockStruct++;
        // solve assigment problem on each block from the left to the right
        while (indexBlock < instance->getE().size()) {
            isWithinTimeLimit();
            Solution::BlockStructure currentBlockStructure = blockStructure;
            heuristicSolver.freeAndAssignmentBlock(currentBlockStructure, indexBlock, nullptr);
            indexBlock++;
            if (solveIdenticalJobs) solveProblemWithFixedCompletionTime(currentBlockStructure, indexBlock);
            if (UB = Solution::evaluate(currentBlockStructure, instance); isSmaller(UB, bestObjValue)) {
                bestObjValue = UB;
                blockStructure = std::move(currentBlockStructure);
                foundBetterNeighbor = true;
                foundBetterSol = true;
            }
        }
        if (foundBetterSol) {
            // foundBetterSol = false;
            // solve assigment problem on each block from the right to the left
            while (indexBlock-->0){
                isWithinTimeLimit();
                assert(indexBlock < instance->getE().size());
                Solution::BlockStructure currentBlockStructure = blockStructure;
                heuristicSolver.freeAndAssignmentBlock(currentBlockStructure, indexBlock, nullptr);
                if (solveIdenticalJobs) solveProblemWithFixedCompletionTime(currentBlockStructure, indexBlock);
                if (UB = Solution::evaluate(currentBlockStructure, instance); isSmaller(UB, bestObjValue)) {
                    bestObjValue = UB;
                    blockStructure = std::move(currentBlockStructure);
                    foundBetterNeighbor = true;
                    // foundBetterSol = true;
                }
            }
        }
    }
    if (solveIdenticalJobs) {
        listScheduledJobs.clear();
        // update the already selected jobs
        std::fill(alreadySelectedJobs.begin(), alreadySelectedJobs.end(), false);
        // update list of selected jobs because some of them must change with identical sub problem
        for (auto &machine : blockStructure){
            for (auto &jobWithCj : machine) {
                if (jobWithCj.first != nullptr) {
                    listScheduledJobs.push_back(jobWithCj.first->getIndex());
                    alreadySelectedJobs[jobWithCj.first->getIndex()].flip();
                }
            }
        }
        std::sort(listScheduledJobs.begin(),listScheduledJobs.end());
    }
    return foundBetterNeighbor;
}

void LocalSearch::localSearchOnlySwapV1() {

    // make local search, by exploring the neighbourhood
    if (verbose >= 2) std::cout << "Begin Local search " << std::endl;

    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    unsigned int nbToSelect = instance->getNbToSelectJob();
    std::vector<unsigned int> listScheduledJobs(nbToSelect); // use a vector of index of job to compute a solution
    std::vector<unsigned int> listJobsNotSelected; // use a vector of index of job that are not selected
    listJobsNotSelected.reserve(instance->getNbJobs() - nbToSelect);
    // First step, found an initial solution
    auto [bestBlockStructure,bestObjValue] = computeInitialSolution(alreadySelectedJobs, listScheduledJobs);

    for (unsigned int indexLoopJobsNotSelected = 0; bool isSelected : alreadySelectedJobs) {
        if (not isSelected) {
            listJobsNotSelected.push_back(indexLoopJobsNotSelected);
        }
        indexLoopJobsNotSelected++;
    }
#ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
#endif
    try {
        nbIter = 0; // the number of iteration of the method
        while (nbIter < maxIter) {
            isWithinTimeLimit();
            bool foundBetterSol = false;
            nbIter++;
            // if it's first iteration try to improve with arrangement
            if (nbIter == 1) {
                improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,&listJobsNotSelected,alreadySelectedJobs, bestObjValue,true);
            }
            if (verbose >= 2) {
                std::cout << "NbIter:" << nbIter << " best know solution:" << bestObjValue << std::endl;
            }
            // loop over leader decisions and find solution in this neighborhood
            bool foundBetterSolLeader = true;
            Solution::BlockStructure bestLeaderBlockStructure = bestBlockStructure;
            while (foundBetterSolLeader) {
                isWithinTimeLimit();
                foundBetterSolLeader = false;
                unsigned int indexJobToAddInSelection = 0;
                unsigned int indexJobToRemoveOfSelection = 0;
                // use neighborhood on the leader decisions
                Neighborhoods neighborhoods(instance, &bestBlockStructure, &listJobsNotSelected,2);
                auto itLastNeighbor = neighborhoods.endLN();
                auto itNeighbor = neighborhoods.beginLN();
                #ifdef DEBUG_HEURISTIC
                auto lastBlockStruct = bestBlockStructure;
                #endif
                for (; itNeighbor != itLastNeighbor ; ++itNeighbor) {
                    isWithinTimeLimit();
                    //check if we don't have exceeded the time limit
                    isWithinTimeLimit();
                    Solution::BlockStructure blockStructFromLeaderN = *itNeighbor;
                    #ifdef DEBUG_HEURISTIC
                    testSol.fromBlockStruct(blockStructFromLeaderN);
                    if (not testSol.feasible(instance))
                        throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                    if (lastBlockStruct == blockStructFromLeaderN) {
                        throw BiSchException("Error in compute leader neighborhood, same neighbor from previous iterator");
                    }
                    lastBlockStruct = blockStructFromLeaderN;
                    #endif
                    if (double UB = Solution::evaluate(blockStructFromLeaderN, instance); isSmaller(UB,bestObjValue)) {
                        foundBetterSolLeader = true;
                        foundBetterSol = true;
                        bestObjValue = UB;
                        bestLeaderBlockStructure = std::move(blockStructFromLeaderN);
                        indexJobToRemoveOfSelection = itNeighbor.getIndexLastSwapJobFromBaseBlockStruct();
                        indexJobToAddInSelection = listJobsNotSelected[itNeighbor.getIndexJobToSwapFromList()];
                        if (not useBestNeighbor) break;
                    }
                }
                #ifdef DEBUG_HEURISTIC
                assert(itNeighbor.nbIter >= neighborhoods.LBNbIterLeaderNeighbor());
                #endif

                // if we have made a swap
                if (indexJobToRemoveOfSelection != indexJobToAddInSelection) {
                    // update the list of selected jobs
                    swapJobsInSelection(listScheduledJobs, indexJobToAddInSelection,indexJobToRemoveOfSelection);
                    // update the list of unselected jobs
                    swapJobsInSelection(listJobsNotSelected, indexJobToRemoveOfSelection,indexJobToAddInSelection);
                    bestBlockStructure = std::move(bestLeaderBlockStructure);
                }
            }
            assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
            // loop over follower decisions and find solution in this neighborhood
            bool foundBetterSolFollower = true;
            #ifdef DEBUG_HEURISTIC
            testSol.fromBlockStruct(bestBlockStructure);
            if (not testSol.feasible(instance))
                throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
            auto lastBlockStruct = bestBlockStructure;
            #endif
            while (foundBetterSolFollower) {
                isWithinTimeLimit();
                //check if we don't have exceeded the time limit
                isWithinTimeLimit();
                // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
                foundBetterSolFollower = improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,&listJobsNotSelected,alreadySelectedJobs, bestObjValue,true);
                if (foundBetterSolFollower) {
                    foundBetterSol = true;
                    if (not useBestNeighbor) break;
                }
            }
            // if we have not found better solution we can stop
            if (not foundBetterSol) break;
        }
    }catch (BiSchTimeOutException &e){} //do nothing

    const auto end{std::chrono::steady_clock::now()};
    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    solution->fromBlockStruct(bestBlockStructure);
}

void LocalSearch::localSearchOnlySwapV2() {
    
    // make local search, by exploring the neighbourhood
    if (verbose >= 2) std::cout << "Begin Local search " << std::endl;

    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    unsigned int nbToSelect = instance->getNbToSelectJob();
    std::vector<unsigned int> listScheduledJobs(nbToSelect); // use a vector of index of job to compute a solution
    std::vector<unsigned int> listJobsNotSelected; // use a vector of index of job that are not selected
    listJobsNotSelected.reserve(instance->getNbJobs() - nbToSelect);
    // First step, found an initial solution
    auto [bestBlockStructure,bestObjValue] = computeInitialSolution(alreadySelectedJobs, listScheduledJobs);

    for (unsigned int indexLoopJobsNotSelected = 0; bool isSelected : alreadySelectedJobs) {
        if (not isSelected) {
            listJobsNotSelected.push_back(indexLoopJobsNotSelected);
        }
        indexLoopJobsNotSelected++;
    }
#ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
#endif

    nbIter = 0; // the number of iteration of the method
    while (nbIter < maxIter) {
        isWithinTimeLimit();
        bool foundBetterSol = false;
        nbIter++;
        // if it's first iteration try to improve with arrangement
        if (nbIter == 1) {
            improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,&listJobsNotSelected,alreadySelectedJobs, bestObjValue,true);
        }
        if (verbose >= 2) {
            std::cout << "NbIter:" << nbIter << " best know solution:" << bestObjValue << std::endl;
        }
        // loop over leader decisions and find solution in this neighborhood
        Solution::BlockStructure bestLeaderBlockStructure = bestBlockStructure;

        unsigned int indexJobToAddInSelection = 0;
        unsigned int indexJobToRemoveOfSelection = 0;
        // use neighborhood on the leader decisions
        Neighborhoods neighborhoods(instance, &bestBlockStructure, &listJobsNotSelected,2);
        auto itLastNeighbor = neighborhoods.endLN();
        auto itNeighbor = neighborhoods.beginLN();
#ifdef DEBUG_HEURISTIC
        auto lastBlockStruct = bestBlockStructure;
#endif
        for (; itNeighbor != itLastNeighbor ; ++itNeighbor) {
            isWithinTimeLimit();
            //check if we don't have exceeded the time limit
            isWithinTimeLimit();
            Solution::BlockStructure blockStructFromLeaderN = *itNeighbor;
#ifdef DEBUG_HEURISTIC
            testSol.fromBlockStruct(blockStructFromLeaderN);
            if (not testSol.feasible(instance))
                throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
            if (lastBlockStruct == blockStructFromLeaderN) {
                throw BiSchException("Error in compute leader neighborhood, same neighbor from previous iterator");
            }
            lastBlockStruct = blockStructFromLeaderN;
#endif
            if (double UB = Solution::evaluate(blockStructFromLeaderN, instance); isSmaller(UB,bestObjValue)) {
                foundBetterSol = true;
                bestObjValue = UB;
                bestLeaderBlockStructure = std::move(blockStructFromLeaderN);
                indexJobToRemoveOfSelection = itNeighbor.getIndexLastSwapJobFromBaseBlockStruct();
                indexJobToAddInSelection = listJobsNotSelected[itNeighbor.getIndexJobToSwapFromList()];
            }
        }

        #ifdef DEBUG_HEURISTIC
        assert(itNeighbor.nbIter >= neighborhoods.LBNbIterLeaderNeighbor());
        #endif

        // if we have made a swap
        if (indexJobToRemoveOfSelection != indexJobToAddInSelection) {
            // update the list of selected jobs
            swapJobsInSelection(listScheduledJobs, indexJobToAddInSelection,indexJobToRemoveOfSelection);
            // update the list of unselected jobs
            swapJobsInSelection(listJobsNotSelected, indexJobToRemoveOfSelection,indexJobToAddInSelection);
            bestBlockStructure = std::move(bestLeaderBlockStructure);
        }

        assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
        // loop over follower decisions and find solution in this neighborhood
#ifdef DEBUG_HEURISTIC
        testSol.fromBlockStruct(bestBlockStructure);
        if (not testSol.feasible(instance))
            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
#endif

        // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
        if (improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,&listJobsNotSelected,alreadySelectedJobs, bestObjValue,true))
            foundBetterSol = true;

#ifdef DEBUG_HEURISTIC
        testSol.fromBlockStruct(bestBlockStructure);
        if (not testSol.feasible(instance))
            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
#endif

        // if we have not found better solution we can stop
        if (not foundBetterSol) break;
    }
    const auto end{std::chrono::steady_clock::now()};
    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    solution->fromBlockStruct(bestBlockStructure);
}

void LocalSearch::localSearchOnlyAssigment() {
    
    // make local search, by exploring the neighbourhood
    if (verbose >= 2) std::cout << "Begin Local search " << std::endl;

    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    unsigned int nbToSelect = instance->getNbToSelectJob();
    std::vector<unsigned int> listScheduledJobs(nbToSelect); // use a vector of index of job to compute a solution
    // First step, found an initial solution
    auto [bestBlockStructure,bestObjValue] = computeInitialSolution(alreadySelectedJobs, listScheduledJobs);

    std::vector<unsigned int> listIndexJobs(instance->getNbJobs());
    std::iota(listIndexJobs.begin(), listIndexJobs.end(), 0);

    #ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
    #endif
    bool foundBetterSol = true;
    nbIter = 0; // the number of iteration of the method
    while (nbIter < maxIter && foundBetterSol) {
        isWithinTimeLimit();
        //check if we don't have exceeded the time limit
        isWithinTimeLimit();
        if (verbose >= 2) {
            std::cout << "NbIter:" << nbIter << " best know solution:" << bestObjValue << std::endl;
        }
        foundBetterSol = false;
        Solution::BlockStructure currentBlockStructure = bestBlockStructure;
        heuristicSolver.upgradeSolutionWithHeuristic(currentBlockStructure, &listIndexJobs);
        double objValue = Solution::evaluate(currentBlockStructure, instance);
        if (isSmaller(objValue,bestObjValue)) {
            bestObjValue = objValue;
            bestBlockStructure = std::move(currentBlockStructure);
            foundBetterSol = true;
        }
        nbIter++;
    }
    const auto end{std::chrono::steady_clock::now()};
    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    solution->fromBlockStruct(bestBlockStructure);
}

void LocalSearch::localSearchPredictor() {
    useBestNeighbor = true;
    // make local search, by exploring the neighbourhood
    if (verbose >= 2) std::cout << "Begin Local search with predictor" << std::endl;

    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    std::vector<unsigned int> listScheduledJobs(instance->getNbToSelectJob()); // use a vector of index of job to compute a solution
    Solution::BlockStructure bestBlockStructure;
    double bestObjValue;
    // First step, found an initial solution
    std::tie(bestBlockStructure,bestObjValue) = computeInitialSolution(alreadySelectedJobs, listScheduledJobs);
    try {
        (version == 6 || version ==7) ? exploreNeighborhoodWithPredictorNotFullTrust(bestBlockStructure,listScheduledJobs,alreadySelectedJobs,bestObjValue)
        : version == 8 ? exploreNeighborhoodWithPredictorFullTrust(bestBlockStructure,listScheduledJobs,alreadySelectedJobs,bestObjValue)
        : exploreRandomNeighborhoodWithPredictorFullTrust(bestBlockStructure,listScheduledJobs,alreadySelectedJobs,bestObjValue);
    }
    catch (BiSchTimeOutException &e){} //do nothing
    auto end{std::chrono::steady_clock::now()};
    end = std::chrono::steady_clock::now();
    // stop time to measure performance
    time_elapsed += std::chrono::duration<double>{end - start};
    solution->fromBlockStruct(bestBlockStructure);
}

void LocalSearch::exploreNeighborhoodWithPredictorNotFullTrust(Solution::BlockStructure& bestBlockStructure, std::vector<unsigned int>& listScheduledJobs, std::vector<bool>& alreadySelectedJobs, double& bestObjValue) {
    torch::NoGradGuard noGradGuard;
    c10::InferenceMode guard;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    if (pathPredictor.empty()) throw BiSchException("No path are provided to load the predictor. See the documentation to use predictor.");
    if (not std::filesystem::path(pathPredictor).has_filename()) throw BiSchException(("Valid path must be used to load the predictor. Path uses is: " + pathPredictor +". See the documentation to use predictor.").c_str());
    auto modelToPredict = torch::jit::load(pathPredictor);
    modelToPredict.eval();
    modelToPredict.to(torch::kCPU);
    std::vector<torch::jit::IValue> inputs(1);
    double prediction;
    if (verbose >=2) {
        std::cout << "Improve using local search only on schedule, initial solution before: wjUj="<< bestObjValue << std::endl;;
    }
    improveBlockStructureLocallyByAssignment(bestBlockStructure, listScheduledJobs, alreadySelectedJobs, bestObjValue, true);
    typedef std::tuple<double,unsigned int,unsigned int> predictionType;
    auto cmp = [](predictionType & left, predictionType & right) {
        return isSmaller(std::get<0>(left),std::get<0>(right));
    };

    // we keep 5 best predictions, if one of them do not lead to improving the solution, then we try the second one and so forth
    std::priority_queue<predictionType,std::vector<predictionType>,decltype(cmp)> listBestPrediction;

    while (nbIter < maxIter) {
        bool foundBetterSol = false;
        nbIter++;
        // try to improve with arrangement and we solve optimally the identical jobs
        if (verbose >= 2) {
            std::cout << "NbIter:" << nbIter << " best know solution:" << bestObjValue << std::endl;
        }
        // if n < N, then we try all swap of selected and not selected jobs
        if (instance->getNbToSelectJob() < instance->getNbJobs()) {
            #ifdef DEBUG_HEURISTIC
            Solution testSol(instance);
            #endif
            if (verbose >= 2) std::cout << "Sample the neighborhood with predictor" << std::endl;
            // we try to swap 2 by 2 all jobs
            for (unsigned int indexJobInSelection = instance->getNbToSelectJob(); indexJobInSelection-- > 0;) {
                // keep the index of the job that we remove in the selection
                unsigned int indexRemovedJobFromSelection = listScheduledJobs[indexJobInSelection];
                alreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
                // now try to add on the leader selection one of the removed jobs
                for (unsigned int indexJobNotSelected = 0; indexJobNotSelected < instance->getNbJobs(); indexJobNotSelected++) {
                    // if the job is not already selected, and jobs we want to swap not belonging to same identical group of jobs
                    if (indexRemovedJobFromSelection != indexJobNotSelected && not alreadySelectedJobs[indexJobNotSelected] && instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                        //check if we don't have exceeded the time limit
                        isWithinTimeLimit();
                        // select the job new job
                        alreadySelectedJobs[indexJobNotSelected] = true;
                        swapJobsInSelection(listScheduledJobs, indexJobNotSelected,indexRemovedJobFromSelection);

                        // make a prediction with the list of scheduled jobs
                        double sumWjWeight = std::accumulate(listScheduledJobs.begin(), listScheduledJobs.end(), 0.0, [&](double sum, unsigned int indexJob) { return sum + instance->getListJobs()[indexJob].getWi(); });

                        auto [blockStruct,sumWjUj] = Solution::blockStructureBySolvingSumCj(listScheduledJobs, instance);
                        computeFeatures(blockStruct);

                        inputs[0] = torch::tensor(std::vector<float>(features.begin(),features.end()),torch::TensorOptions().dtype(torch::kFloat32)).view({1, 96});

                        torch::jit::IValue result = modelToPredict.forward(inputs);
                        prediction = result.toTensor().item<double>();
                        listBestPrediction.emplace(prediction * sumWjWeight,indexRemovedJobFromSelection, indexJobNotSelected);
                        if (listBestPrediction.size()>nbPrediction) listBestPrediction.pop();

                        // unselect the job the new job
                        alreadySelectedJobs[indexJobNotSelected] = false;
                        swapJobsInSelection(listScheduledJobs, indexRemovedJobFromSelection, indexJobNotSelected);
                    }
                }
                alreadySelectedJobs[indexRemovedJobFromSelection] = true;
            }
            if (verbose >= 2) std::cout << "Try to improve each neighbor" << std::endl;
            // get the prediction in the reverse order because we have a max-heap
            std::vector<predictionType> &res = Container(listBestPrediction);
            std::sort(res.begin(), res.end());
            auto currentBlockStructure = bestBlockStructure;
            // we try to swap 2 by 2 all jobs
            std::pair<unsigned int, unsigned int> bestSwap({0, 0}); // pair for the best swap between the removed job and the new one in the selection
            for (auto &[_,indexJobInsert,indexJobRemove] : res){
                // we add to make the best swapped job as a selected one
                alreadySelectedJobs[indexJobInsert] = false;
                alreadySelectedJobs[indexJobRemove] = true;

                auto [newBlockStructure,newObjValue] = updateBlockStructureWithSwap(
                            listScheduledJobs, currentBlockStructure, indexJobRemove, indexJobInsert);
                #ifdef DEBUG_HEURISTIC
                testSol.fromBlockStruct(newBlockStructure);
                if (not testSol.feasible(instance))
                    throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                #endif

                // sort machines by makespan
                Solution::sortBlockStructurePerMakespan(newBlockStructure,instance);
                if (isSmaller(newObjValue, bestObjValue)) {
                    foundBetterSol = true;
                    bestSwap = {indexJobInsert, indexJobRemove};
                    bestBlockStructure = std::move(newBlockStructure);
                    bestObjValue = newObjValue;
                    assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
                    // if we found better solution and use first improvement
                    if (not useBestNeighbor) break;
                    // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
                    version == 6 ? improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr,alreadySelectedJobs, bestObjValue)
                    : improveBlockStructureLocallyByAssignment(bestBlockStructure, listScheduledJobs,alreadySelectedJobs, bestObjValue);
                }
                // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
                else if (
                    version == 6 ?improveBlockStructureLocallyBySwap(newBlockStructure, listScheduledJobs,nullptr,alreadySelectedJobs, bestObjValue)
                        : improveBlockStructureLocallyByAssignment(newBlockStructure, listScheduledJobs,alreadySelectedJobs, bestObjValue)
                    ) {
                    bestSwap = {indexJobInsert, indexJobRemove};
                    bestBlockStructure = std::move(newBlockStructure);
                    assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
                    foundBetterSol = true;
                    if (not useBestNeighbor) break;
                }
                if (foundBetterSol && not useBestNeighbor) break;
                //else reset the original list of jobs
                alreadySelectedJobs[indexJobInsert] = true;
                alreadySelectedJobs[indexJobRemove] = false;
                swapJobsInSelection(listScheduledJobs, indexJobInsert,indexJobRemove);
            }
            //if we don't have found better solution then stop
            if (!foundBetterSol) break;
            if (useBestNeighbor) {
                assert(bestSwap.first != bestSwap.second);
                // we add to make the best swapped job as a selected one
                alreadySelectedJobs[bestSwap.first] = false;
                alreadySelectedJobs[bestSwap.second] = true;
                swapJobsInSelection(listScheduledJobs, bestSwap.second, bestSwap.first);
            }
            if (verbose >= 2) std::cout << res.size() << " solutions tried" << std::endl;
            res.clear();
        }
    }
    if (verbose >=2) {
        std::cout << "Improve using local search only on schedule, final solution before: wjUj="<< bestObjValue;
    }
    // try to improve a last using the other heuristic
    improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr, alreadySelectedJobs, bestObjValue, true);
    if (verbose >=2) {
        std::cout << "final solution after: wjUj="<< bestObjValue << std::endl;
    }
}

inline bool LocalSearch::computePredictions(torch::jit::Module& model,torch::Tensor& batchTensor, std::vector<std::vector<float>>& batchFeatures,std::vector<BatchItem> &batchItems,std::priority_queue<predictionWithSetJob,std::vector<predictionWithSetJob>,ComparePredictionWithSet<>>&listBestPrediction) {
    if (batchFeatures.empty()) {
        return {};
    }

    const size_t batchSize = batchFeatures.size();

    // Get direct pointer to tensor data for efficient copying
    float* tensorData = batchTensor.data_ptr<float>();

    // Copy all features contiguously into tensor memory
    for (size_t i = 0; i < batchSize; ++i) {
        const size_t featureSize = 96;
        // Safety check: ensure feature vector has correct size
        if (batchFeatures[i].size() != featureSize) {
            throw BiSchException("Feature vector must contain exactly 96 elements");
        }
        // Copy directly to the appropriate position in tensor memory
        std::memcpy(tensorData + i * featureSize, batchFeatures[i].data(), featureSize * sizeof(float));
    }

    // Single forward pass for the entire batch
    inputs.emplace_back(batchTensor);
    torch::jit::IValue result = model.forward(inputs);
    inputs.clear();

    // Convert tensor to vector of doubles
    const torch::Tensor& predictionsTensor = result.toTensor();

    // Extract predictions efficiently
    auto predictionsAccessor = predictionsTensor.accessor<float, 2>();
    // Process each prediction in the batch
    bool foundBetterPred = false;
    for (unsigned int i = 0; i < batchSize; ++i) {
        double predictedObjValue = predictionsAccessor[i][0] * batchItems[i].weightSum;
        if (isSmaller(predictedObjValue, std::get<0>(listBestPrediction.top()))) {
            foundBetterPred = true;
            listBestPrediction.emplace(predictedObjValue,std::move(batchItems[i].scheduledJobs),std::move(batchItems[i].selectedJobs));
            if (listBestPrediction.size() > nbPrediction) {
                listBestPrediction.pop();
            }
        }
    }

    // Clear batch data for next iteration
    batchFeatures.clear();
    batchItems.clear();
    return foundBetterPred;
}

void LocalSearch::exploreRandomNeighborhoodWithPredictorFullTrust(Solution::BlockStructure& bestBlockStructure, std::vector<unsigned int>& listScheduledJobs, std::vector<bool>& alreadySelectedJobs, double& bestObjValue) {
    if (verbose >=2) {
        std::cout << "Improve using local search only on schedule, initial solution before: wjUj="<< bestObjValue << std::endl;;
    }
    improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr, alreadySelectedJobs, bestObjValue, true);

    //random part
    std::mt19937 gen(0);
    std::vector<unsigned int> sampledIndexToRemove(instance->getNbToSelectJob()); //vector to get the sample elements
    unsigned int ratioNeighborIndexToRemove = std::floor(instance->getNbToSelectJob()*ratioNeighbor);
    unsigned int ratioNeighborIndexToInsert = std::floor(instance->getNbJobs()*ratioNeighbor);
    std::vector<unsigned int> sampledIndexToInsert(instance->getNbJobs()); //vector to get the sample elements
    unsigned int nbNoImprovement = 0; // count number of no improvement and increase the ratio if it reach 5
    bool incr = true;
    // Pytorch Model
    torch::NoGradGuard noGradGuard;
    c10::InferenceMode guard;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    if (pathPredictor.empty()) throw BiSchException("No path are provided to load the predictor. See the documentation to use predictor.");
    if (not std::filesystem::path(pathPredictor).has_filename()) throw BiSchException(("Valid path must be used to load the predictor. Path uses is: " + pathPredictor +". See the documentation to use predictor.").c_str());
    auto modelToPredict = torch::jit::load(pathPredictor);
    modelToPredict.eval();
    modelToPredict.to(torch::kCPU);

    // Configuration constants
    auto options = torch::TensorOptions()
    .dtype(torch::kFloat32)
    .layout(torch::kStrided)
    .device(torch::kCPU)
    .requires_grad(false);
    // Create an empty tensor with the correct shape (no initialization)
    torch::Tensor batchTensor = torch::zeros({static_cast<long>(BATCH_SIZE),96}, options);
    inputs.emplace_back(batchTensor); // add the tensor to reserve the good size once
    inputs.clear(); // remove the dummy tensor

    // Vectors for batch accumulation
    std::vector<std::vector<float>> batchFeatures;
    std::vector<BatchItem> batchItems;
    batchFeatures.reserve(BATCH_SIZE);
    batchItems.reserve(BATCH_SIZE);

    // we have a list of active prediction
    std::vector<predictionWithSetJob> listActivePrediction;
    listActivePrediction.reserve(nbPrediction);
    listActivePrediction.emplace_back(std::numeric_limits<double>::infinity(), listScheduledJobs, alreadySelectedJobs);
    // we keep 'nbPrediction' best predictions
    std::priority_queue<predictionWithSetJob,std::vector<predictionWithSetJob>,ComparePredictionWithSet<>> listBestPrediction;
    Container(listBestPrediction).reserve(nbPrediction+1);
    auto dummyPred = predictionWithSetJob();
    std::get<0>(dummyPred) = std::numeric_limits<double>::infinity();
    listBestPrediction.emplace(dummyPred);

    // create copy of block structure and sum Wj Uj
    auto blockStruct = bestBlockStructure;
    auto sumWjUj = bestObjValue;
    try{
        // if n < N, then we try all swap of selected and not selected jobs
        if (instance->getNbToSelectJob() < instance->getNbJobs()) {
            while (nbIter < maxIter) {
                nbIter++;
                // try to improve with arrangement and we solve optimally the identical jobs
                if (verbose >= 2) {
                    std::cout << "NbIter:" << nbIter << " with " << listActivePrediction.size() << " predictions where the largest:" << std::get<0>(listBestPrediction.top()) << std::endl;
                }
                bool foundBetterPred = false;

                for (auto &[activePrediction,activeListScheduledJobs,activeAlreadySelectedJobs] : listActivePrediction){
                    // we try to swap 2 by 2 all jobs
                    // but we keep only 'ratioNeighbor' percent
                    assert(ratioNeighborIndexToRemove <= instance->getNbToSelectJob());
                    std::ranges::sample(std::views::iota(0u, instance->getNbToSelectJob()), sampledIndexToRemove.begin(), ratioNeighborIndexToRemove, gen);
                    for (unsigned int indexJobInSelection : std::ranges::reverse_view(sampledIndexToRemove)) {
                        assert(not activeListScheduledJobs.empty());
                        // keep the index of the job that we remove in the selection
                        unsigned int indexRemovedJobFromSelection = activeListScheduledJobs[indexJobInSelection];
                        activeAlreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
                        // now try to add on the leader selection one of the removed jobs
                        //keep only a subset
                        assert(ratioNeighborIndexToInsert <= instance->getNbJobs());
                        std::ranges::sample(std::views::iota(0u, instance->getNbJobs()), sampledIndexToInsert.begin(), ratioNeighborIndexToInsert, gen);
                        for (unsigned int indexJobNotSelected : sampledIndexToInsert) {
                            // if the job is not already selected, and jobs we want to swap not belonging to same identical group of jobs
                            if (indexRemovedJobFromSelection != indexJobNotSelected && not activeAlreadySelectedJobs[indexJobNotSelected] && instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                                //check if we don't have exceeded the time limit
                                isWithinTimeLimit();
                                // select the job new job
                                activeAlreadySelectedJobs[indexJobNotSelected] = true;
                                swapJobsInSelection(activeListScheduledJobs, indexJobNotSelected,indexRemovedJobFromSelection);

                                // make a prediction with the list of scheduled jobs
                                double sumWjWeight = std::accumulate(activeListScheduledJobs.begin(), activeListScheduledJobs.end(), 0.0, [&](double sum, unsigned int indexJob) { return sum + instance->getListJobs()[indexJob].getWi(); });
                                Solution::blockStructureBySolvingSumCj(blockStruct,sumWjUj,activeListScheduledJobs, instance);
                                computeFeatures(blockStruct);

                                // Add to batch
                                batchFeatures.emplace_back(features.begin(), features.end());
                                batchItems.emplace_back(
                                    activeListScheduledJobs,
                                    activeAlreadySelectedJobs,
                                    sumWjWeight
                                );

                                // Process batch if it reaches optimal size
                                if (batchFeatures.size() >= BATCH_SIZE) {
                                    foundBetterPred |= computePredictions(modelToPredict,batchTensor,batchFeatures,batchItems,listBestPrediction);
                                }
                                // unselect the job the new job
                                activeAlreadySelectedJobs[indexJobNotSelected] = false;
                                swapJobsInSelection(activeListScheduledJobs, indexRemovedJobFromSelection, indexJobNotSelected);
                            }
                        }
                        activeAlreadySelectedJobs[indexRemovedJobFromSelection] = true;
                    }
                    // Process any remaining items in the batch processRemainingBatch:
                    foundBetterPred |= computePredictions(modelToPredict,batchTensor,batchFeatures,batchItems,listBestPrediction);
                }
                if (not foundBetterPred && ++nbNoImprovement == 5 ) {
                    if (incr){
                        if (isEqual(0.1,ratioNeighbor)) {
                            incr = false;
                            ratioNeighbor -= 0.01;
                        }else ratioNeighbor += 0.01;
                    }else{
                        if (isEqual(0.01,ratioNeighbor)) {
                            incr = true;
                            ratioNeighbor += 0.01;
                        }else ratioNeighbor -= 0.01;
                    }
                    nbNoImprovement = 0;
                    if (verbose >= 2) std::cout << "No improvement since 5 iterations : increase up ratio to " << ratioNeighbor << std::endl;
                    ratioNeighborIndexToRemove = std::floor(instance->getNbToSelectJob()*ratioNeighbor);
                    ratioNeighborIndexToInsert = std::floor(instance->getNbJobs()*ratioNeighbor);
                    continue;
                }

                std::cout << "nb neighbors " << listBestPrediction.size() << std::endl;
                listActivePrediction.clear();
                //keep the greatest prediction to re add it after to ensure we're always looking for better prediction
                auto greatestPrediction = listBestPrediction.top();

                listActivePrediction = std::move(Container(listBestPrediction));
                Container(listBestPrediction).clear();
                // add the greatest prediction to the list
                listBestPrediction.push(std::move(greatestPrediction));
            }
        }
    }catch (BiSchTimeOutException &e) {
        resetTimer(5);
        std::cout << "Give 5s to improve best found sol" << std::endl;
    }
    while (not listActivePrediction.empty()){
            auto &[pred,activeListScheduledJobs,activeAlreadySelectedJobs] = listActivePrediction.back();
            // compute the corresponding schedule and run local search on it
            auto [newBlockStructure,objValue] = Solution::blockStructureBySolvingSumCj(activeListScheduledJobs, instance);
            if (verbose >=2) {
                std::cout << "Improve using local search only on schedule, solution before: wjUj="<< objValue << " where prediction was: " << pred;
            }

            improveBlockStructureLocallyByAssignment(newBlockStructure, activeListScheduledJobs,activeAlreadySelectedJobs, objValue, true);
            // improveBlockStructureLocallyBySwap(newBlockStructure, activeListScheduledJobs,nullptr,activeAlreadySelectedJobs, objValue, true);
            // unsigned int nbImprove = 0;
            // while (nbImprove++ < maxIter && improveBlockStructureLocallyBySwap(newBlockStructure, activeListScheduledJobs,nullptr,activeAlreadySelectedJobs, objValue, true))

            if (objValue < bestObjValue) {
                bestBlockStructure = newBlockStructure;
                bestObjValue = objValue;
            }
            if (verbose >=2) {
                std::cout << " solution after: wjUj="<< objValue  << std::endl;
            }
            listActivePrediction.pop_back();
    }
}

void LocalSearch::exploreNeighborhoodWithPredictorFullTrust(Solution::BlockStructure& bestBlockStructure, std::vector<unsigned int>& listScheduledJobs, std::vector<bool>& alreadySelectedJobs, double& bestObjValue) {

    if (verbose >=2) {
        std::cout << "Improve using local search only on schedule, initial solution before: wjUj="<< bestObjValue << std::endl;;
    }
    improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr, alreadySelectedJobs, bestObjValue, true);

    // Pytorch Model
    torch::NoGradGuard noGradGuard;
    c10::InferenceMode guard;
    // Deserialize the ScriptModule from a file using torch::jit::load().
    if (pathPredictor.empty()) throw BiSchException("No path are provided to load the predictor. See the documentation to use predictor.");
    if (not std::filesystem::path(pathPredictor).has_filename()) throw BiSchException(("Valid path must be used to load the predictor. Path uses is: " + pathPredictor +". See the documentation to use predictor.").c_str());
    auto modelToPredict = torch::jit::load(pathPredictor);
    modelToPredict.eval();
    modelToPredict.to(torch::kCPU);

    // Configuration constants
    auto options = torch::TensorOptions()
    .dtype(torch::kFloat32)
    .layout(torch::kStrided)
    .device(torch::kCPU)
    .requires_grad(false);
    // Create an empty tensor with the correct shape (no initialization)
    torch::Tensor batchTensor = torch::zeros({static_cast<long>(BATCH_SIZE),96}, options);
    inputs.emplace_back(batchTensor); // add the tensor to reserve the good size once
    inputs.clear(); // remove the dummy tensor

    // Vectors for batch accumulation
    std::vector<std::vector<float>> batchFeatures;
    std::vector<BatchItem> batchItems;
    batchFeatures.reserve(BATCH_SIZE);
    batchItems.reserve(BATCH_SIZE);

    // we have a list of active prediction
    std::vector<predictionWithSetJob> listActivePrediction;
    listActivePrediction.reserve(nbPrediction);
    listActivePrediction.emplace_back(std::numeric_limits<double>::infinity(), listScheduledJobs, alreadySelectedJobs);
    // we keep 'nbPrediction' best predictions
    std::priority_queue<predictionWithSetJob,std::vector<predictionWithSetJob>,ComparePredictionWithSet<>> listBestPrediction;
    Container(listBestPrediction).reserve(nbPrediction+1);
    auto dummyPred = predictionWithSetJob();
    std::get<0>(dummyPred) = std::numeric_limits<double>::infinity();
    listBestPrediction.emplace(dummyPred);
    // std::unordered_set<std::vector<bool>> listBestPredictionInQueue;
    // listBestPredictionInQueue.reserve(nbPrediction*maxIter);

    // create copy of block structure and sum Wj Uj
    auto blockStruct = bestBlockStructure;
    auto sumWjUj = bestObjValue;
    // if n < N, then we try all swap of selected and not selected jobs
    if (instance->getNbToSelectJob() < instance->getNbJobs()) {
        while (nbIter < maxIter) {
            nbIter++;
            // try to improve with arrangement and we solve optimally the identical jobs
            if (verbose >= 2) {
                std::cout << "NbIter:" << nbIter << " with " << listActivePrediction.size() << " predictions where the largest:" << std::get<0>(listBestPrediction.top()) << std::endl;
            }
            bool foundBetterPred = false;
            for (auto &[activePrediction,activeListScheduledJobs,activeAlreadySelectedJobs] : listActivePrediction){
                // we try to swap 2 by 2 all jobs
                for (unsigned int indexJobInSelection = instance->getNbToSelectJob(); indexJobInSelection-- > 0;) {
                    assert(not activeListScheduledJobs.empty());
                    // keep the index of the job that we remove in the selection
                    unsigned int indexRemovedJobFromSelection = activeListScheduledJobs[indexJobInSelection];
                    activeAlreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
                    // now try to add on the leader selection one of the removed jobs
                    for (unsigned int indexJobNotSelected = 0; indexJobNotSelected < instance->getNbJobs(); indexJobNotSelected++) {
                        // if the job is not already selected, and jobs we want to swap not belonging to same identical group of jobs
                        if (indexRemovedJobFromSelection != indexJobNotSelected && not activeAlreadySelectedJobs[indexJobNotSelected] && instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                            //check if we don't have exceeded the time limit
                            isWithinTimeLimit();
                            // select the job new job
                            activeAlreadySelectedJobs[indexJobNotSelected] = true;
                            swapJobsInSelection(activeListScheduledJobs, indexJobNotSelected,indexRemovedJobFromSelection);

                            // make a prediction with the list of scheduled jobs
                            double sumWjWeight = std::accumulate(activeListScheduledJobs.begin(), activeListScheduledJobs.end(), 0.0, [&](double sum, unsigned int indexJob) { return sum + instance->getListJobs()[indexJob].getWi(); });
                            Solution::blockStructureBySolvingSumCj(blockStruct,sumWjUj,activeListScheduledJobs, instance);
                            computeFeatures(blockStruct);

                            // Add to batch
                            batchFeatures.emplace_back(features.begin(), features.end());
                            batchItems.emplace_back(
                                activeListScheduledJobs,
                                activeAlreadySelectedJobs,
                                sumWjWeight
                            );

                            // Process batch if it reaches optimal size
                            if (batchFeatures.size() >= BATCH_SIZE) {
                                foundBetterPred |= computePredictions(modelToPredict,batchTensor,batchFeatures,batchItems,listBestPrediction);
                            }
                            // unselect the job the new job
                            activeAlreadySelectedJobs[indexJobNotSelected] = false;
                            swapJobsInSelection(activeListScheduledJobs, indexRemovedJobFromSelection, indexJobNotSelected);
                        }
                    }
                    activeAlreadySelectedJobs[indexRemovedJobFromSelection] = true;
                }
                // Process any remaining items in the batch processRemainingBatch:
                foundBetterPred |= computePredictions(modelToPredict,batchTensor,batchFeatures,batchItems,listBestPrediction);
            }
            if (not foundBetterPred) break;
            std::cout << "nb neighbors " << listBestPrediction.size() << std::endl;
            listActivePrediction.clear();
            //keep the greatest prediction to re add it after to ensure we're always looking for better prediction
            auto greatestPrediction = listBestPrediction.top();

            listActivePrediction = std::move(Container(listBestPrediction));
            Container(listBestPrediction).clear();
            // add the greatest prediction to the list
            listBestPrediction.push(std::move(greatestPrediction));
        }
    }
    while (not listActivePrediction.empty()){
            auto &[pred,activeListScheduledJobs,activeAlreadySelectedJobs] = listActivePrediction.back();
            // compute the corresponding schedule and run local search on it
            auto [newBlockStructure,objValue] = Solution::blockStructureBySolvingSumCj(activeListScheduledJobs, instance);
            if (verbose >=2) {
                std::cout << "Improve using local search only on schedule, solution before: wjUj="<< objValue << " where prediction was: " << pred;
            }

            improveBlockStructureLocallyByAssignment(newBlockStructure, activeListScheduledJobs,activeAlreadySelectedJobs, objValue, true);
            // improveBlockStructureLocallyBySwap(newBlockStructure, activeListScheduledJobs,nullptr,activeAlreadySelectedJobs, objValue, true);
            // unsigned int nbImprove = 0;
            // while (nbImprove++ < maxIter && improveBlockStructureLocallyBySwap(newBlockStructure, activeListScheduledJobs,nullptr,activeAlreadySelectedJobs, objValue, true))

            if (objValue < bestObjValue) {
                bestBlockStructure = newBlockStructure;
                bestObjValue = objValue;
            }
            if (verbose >=2) {
                std::cout << " solution after: wjUj="<< objValue  << std::endl;
            }
            listActivePrediction.pop_back();
    }
}

void LocalSearch::localSearchRandom() {
    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    std::vector<unsigned int> listScheduledJobs(instance->getNbToSelectJob()); // use a vector of index of job to compute a solution
    auto &listOfJobs = instance->getListJobs();
    // First step, found an initial solution
    std::vector<size_t> indicesInitialSolution(listOfJobs.size());
    std::iota(indicesInitialSolution.begin(), indicesInitialSolution.end(), 0);
    std::mt19937 g(0);
    std::shuffle(indicesInitialSolution.begin(), indicesInitialSolution.end(),g);
    // add this job to the tree structure.
    for (unsigned int indexLoopJobInitialSol = 0; indexLoopJobInitialSol < instance->getNbToSelectJob(); ++indexLoopJobInitialSol) {
        const Job& selectedJob = listOfJobs[indicesInitialSolution[indexLoopJobInitialSol]];
        alreadySelectedJobs[selectedJob.getIndex()].flip();
        listScheduledJobs[indexLoopJobInitialSol] = selectedJob.getIndex();
    }

    // sort the index in listOfJob to get SPT order
    std::sort(listScheduledJobs.begin(),listScheduledJobs.end());
    // get the block structure solution from the list of index
    auto [blockStructure,objValue] = Solution::blockStructureBySolvingSumCj(listScheduledJobs, instance);
    improveBlockStructureLocallyBySwap(blockStructure, listScheduledJobs,nullptr, alreadySelectedJobs, objValue, true);
    solution->fromBlockStruct(blockStructure);
}

void LocalSearch::exploreNeighborhood(Solution::BlockStructure & bestBlockStructure, std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs, double &bestObjValue) {
    //sort the machine according non-decreasing order of makespan
    Solution::sortBlockStructurePerMakespan(bestBlockStructure,instance);

    bool foundBetterSol = false;

    nbIter = 0; // the number of iteration of the method
    while (nbIter < maxIter) {
        isWithinTimeLimit();
        nbIter++;
        // try to improve with arrangement and we solve optimally the identical jobs
        foundBetterSol = version == 1 ?
             improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr,alreadySelectedJobs, bestObjValue,true)
            :improveBlockStructureLocallyByAssignment(bestBlockStructure, listScheduledJobs,alreadySelectedJobs, bestObjValue,true);
        if (verbose >= 2) {
            std::cout << "NbIter:" << nbIter << " best know solution:" << bestObjValue << std::endl;
        }
        // if n < N, then we try all swap of selected and not selected jobs
        if (instance->getNbToSelectJob() < instance->getNbJobs()) {
            foundBetterSol = false; // try to find better solution
            #ifdef DEBUG_HEURISTIC
            Solution testSol(instance);
            #endif
            unsigned int nbSol = 0;
            Solution::BlockStructure currentBlockStructure = bestBlockStructure; // use a current block structure that it change
            // we try to swap 2 by 2 all jobs
            std::pair<unsigned int, unsigned int> bestSwap({0, 0}); // pair for the best swap between the removed job and the new one in the selection
            for (unsigned int indexJobInSelection = instance->getNbToSelectJob(); indexJobInSelection-- > 0;) {
                // keep the index of the job that we remove in the selection
                unsigned int indexRemovedJobFromSelection = listScheduledJobs[indexJobInSelection];
                alreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
                // now try to add on the leader selection one of the removed jobs
                for (unsigned int indexJobNotSelected = 0; indexJobNotSelected < instance->getNbJobs(); indexJobNotSelected++) {
                    // if the job is not already selected, and jobs we want to swap not belonging to same identical group of jobs
                    if (indexRemovedJobFromSelection != indexJobNotSelected && not alreadySelectedJobs[indexJobNotSelected] && instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                        //check if we don't have exceeded the time limit
                        isWithinTimeLimit();
                        // select the job new job
                        alreadySelectedJobs[indexJobNotSelected] = true;
                        if (indexRemovedJobFromSelection == 494 && indexJobNotSelected == 495)
                            std::cout << "";
                        nbSol++;
                        auto [newBlockStructure,newObjValue] = updateBlockStructureWithSwap(
                            listScheduledJobs, currentBlockStructure, indexJobNotSelected, indexRemovedJobFromSelection);
                        #ifdef DEBUG_HEURISTIC
                        testSol.fromBlockStruct(newBlockStructure);
                        if (not testSol.feasible(instance))
                            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                        #endif

                        // sort machines by makespan
                        Solution::sortBlockStructurePerMakespan(newBlockStructure,instance);
                        if (isSmaller(newObjValue, bestObjValue)) {
                            foundBetterSol = true;
                            bestSwap = {indexRemovedJobFromSelection, indexJobNotSelected};
                            bestBlockStructure = std::move(newBlockStructure);
                            bestObjValue = newObjValue;
                            assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
                            // if we found better solution and use first improvement
                            if (not useBestNeighbor) break;
                            // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
                            version == 1 ? improveBlockStructureLocallyBySwap(bestBlockStructure, listScheduledJobs,nullptr,alreadySelectedJobs, bestObjValue)
                                : improveBlockStructureLocallyByAssignment(bestBlockStructure, listScheduledJobs,alreadySelectedJobs, bestObjValue);
                        }
                        // try to improve solution with arrangement, i.e. try follower neighborhood on the leader decisions
                        else if (version == 1 ? improveBlockStructureLocallyBySwap(newBlockStructure, listScheduledJobs,nullptr,alreadySelectedJobs, bestObjValue)
                                : improveBlockStructureLocallyByAssignment(newBlockStructure, listScheduledJobs,alreadySelectedJobs, bestObjValue)) {
                            bestSwap = {indexRemovedJobFromSelection, indexJobNotSelected};
                            bestBlockStructure = std::move(newBlockStructure);
                            assert(bestObjValue == Solution::evaluate(bestBlockStructure,instance));
                            foundBetterSol = true;
                            if (not useBestNeighbor) break;
                        }
                        if (indexRemovedJobFromSelection == 4 && indexJobNotSelected == 0)
                            std::cout << "";
                        // unselect the job the new job
                        alreadySelectedJobs[indexJobNotSelected] = false;
                        swapJobsInSelection(listScheduledJobs, indexRemovedJobFromSelection, indexJobNotSelected);
                    }
                }
                alreadySelectedJobs[indexRemovedJobFromSelection] = true;
                // if we are in first improvement then stop
                if (foundBetterSol && not useBestNeighbor) break;
            }
            if (not foundBetterSol) break;
            if (useBestNeighbor) {
                assert(bestSwap.first != bestSwap.second);
                // we add to make the best swapped job as a selected one
                alreadySelectedJobs[bestSwap.first] = false;
                alreadySelectedJobs[bestSwap.second] = true;
                swapJobsInSelection(listScheduledJobs, bestSwap.second, bestSwap.first);
            }
        }
        if (not foundBetterSol) break;
    }
}

void LocalSearch::localSearchBySwapAllJobs(Solution::BlockStructure * initialBlockStructure, double initialSolutionValue) {
    
    // make local search, by exploring the neighbourhood
    if (verbose >= 2) std::cout << "Begin Local search " << std::endl;

    // vector of bool to know which jobs that have been already selected
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    std::vector<unsigned int> listScheduledJobs(instance->getNbToSelectJob()); // use a vector of index of job to compute a solution
    Solution::BlockStructure bestBlockStructure;
    double bestObjValue;
    if (initialBlockStructure == nullptr) {
        // First step, found an initial solution
        std::tie(bestBlockStructure,bestObjValue) = computeInitialSolution(alreadySelectedJobs, listScheduledJobs);
    }else {
        bestBlockStructure = *initialBlockStructure;
        bestObjValue = initialSolutionValue;
        unsigned int indexLoopListScheduledJobs = 0;
        // update list index and bool
        for (auto &machine: bestBlockStructure) {
            for (auto &job : machine) {
                if (job.first != nullptr) {
                    listScheduledJobs[indexLoopListScheduledJobs] = job.first->getIndex();
                    alreadySelectedJobs[job.first->getIndex()].flip();
                    indexLoopListScheduledJobs++;
                }
            }
        }
        std::sort(listScheduledJobs.begin(),listScheduledJobs.end());
    }
    try {
        exploreNeighborhood(bestBlockStructure,listScheduledJobs,alreadySelectedJobs,bestObjValue);
    }catch (BiSchTimeOutException &e) {}
    auto end{std::chrono::steady_clock::now()};
    end = std::chrono::steady_clock::now();
    // stop time to measure performance
    time_elapsed += std::chrono::duration<double>{end - start};
    solution->fromBlockStruct(bestBlockStructure);
}

std::pair<Solution::BlockStructure, double> LocalSearch::computeInitialSolution(std::vector<bool>& alreadySelectedJobs, std::vector<unsigned int>& listScheduledJobs) {
    const std::vector<Job>& listOfJobs = instance->getListJobs();
    unsigned int nbToSelect = instance->getNbToSelectJob();
    // First step, found an initial solution
    std::vector<size_t> indicesInitialSolution(listOfJobs.size());
    std::iota(indicesInitialSolution.begin(), indicesInitialSolution.end(), 0);
    // sort jobs according to (p_j - d_j) / w_j rules
    auto sortingRule = [&listOfJobs](size_t left, size_t right) {
        return isSmaller((listOfJobs[left].getPi() - listOfJobs[left].getDi()) / listOfJobs[left].getWi(), (listOfJobs[right].getPi() - listOfJobs[right].getDi()) / listOfJobs[right].getWi()
        );
    };
    std::partial_sort(indicesInitialSolution.begin(), indicesInitialSolution.begin() + nbToSelect, indicesInitialSolution.end(), sortingRule);
    // add this job to the tree structure.
    for (unsigned int indexLoopJobInitialSol = 0; indexLoopJobInitialSol < nbToSelect; ++indexLoopJobInitialSol) {
        const Job& selectedJob = listOfJobs[indicesInitialSolution[indexLoopJobInitialSol]];
        alreadySelectedJobs[selectedJob.getIndex()].flip();
        listScheduledJobs[indexLoopJobInitialSol] = selectedJob.getIndex();
    }

    // sort the index in listOfJob to get SPT order
    std::sort(listScheduledJobs.begin(),listScheduledJobs.end());
    // get the block structure solution from the list of index
    auto [blockStructure,objValue] = Solution::blockStructureBySolvingSumCj(listScheduledJobs, instance);
    if (isSmaller(firstUB,0)) {
        firstUB = objValue;
    }
    return {blockStructure,objValue};
}

void LocalSearch::solve() {
    // start time to measure performance

    heuristicSolver.setTimeUp(timeUp);
    if (genDatabase) {
        auto listOfJobs = instance->getListJobs();
        std::sort(listOfJobs.begin(), listOfJobs.end(), std::greater<>());
        if (listOfJobs.size() != instance->getNbToSelectJob())
            throw BiSchException("To generate database, we need to use instance where N=n.");
        auto blockStruct = Solution::solveSumCjCriteria(listOfJobs, instance).toBlockStruct(instance);
        computeFeatures(blockStruct);
    } else {
        try{
            switch (version) {
            case 1: case 5:
                localSearchBySwapAllJobs();
                break;
            case 2:
                localSearchOnlySwapV1();
                break;
            case 3:
                localSearchOnlySwapV2();
                break;
            case 4:
                localSearchOnlyAssigment();
                break;
            case 6: case 7 : case 8: case 9:
                localSearchPredictor();
                break;
            case 10:
                localSearchRandom();
                break;
            default:
                throw BiSchException("The version of the local search is not implemented");
            }
        }catch (const BiSchTimeOutException &e) {} //do nothing

        solution->evaluate();
        const auto end{std::chrono::steady_clock::now()};
        // stop time to measure performance
        time_elapsed = std::chrono::duration<double>{end - start};
        if (verbose >= 1) std::cout << "Heuristic is over after " << time_elapsed.count() << " seconds" << std::endl << "The objective value is : " << solution->getSumWjUj() << std::endl;
    }
}

void LocalSearch::printOutput(std::string& fileOutputName, std::ofstream& outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists) {
        outputFile << "InstanceName" << "\t";
        if (genDatabase) {
            for (unsigned int indexFeatures = 1; indexFeatures <= 95; ++indexFeatures) { outputFile << "feat_" << indexFeatures << "\t"; }
            outputFile << "feat_96" << std::endl;
        }else {
            outputFile << "InstancePath"
                << "\t" << "N"
                << "\t" << "n"
                << "\t" << "m_Max"
                << "\t" << "m_0"
                << "\t" << "V_max"
                << "\t" << "V_0"
                << "\t" << "Method"
                << "\t" << "Predictor"
                << "\t" << "Version"
                << "\t" << "UseBestImprovement"
                << "\t" << "Time"
                << "\t" << "LimitTime"
                << "\t" << "First_UB"
                << "\t" << "NbIter"
                << "\t" << "MaxNbIter"
                << "\t" << "Objective" << std::endl;
        }
    }
    outputFile << instance->getInstanceName() << "\t";
    if (genDatabase) {
        for (unsigned int indexFeatures = 0; indexFeatures < 95; ++indexFeatures) { outputFile << features[indexFeatures] << "\t"; }
        outputFile << features[95] << std::endl;
    }
    else {
        solution->evaluate();
        // write value
        outputFile << instance->getInstancePath().string()
        << "\t" << instance->getNbJobs()
        << "\t" << instance->getNbToSelectJob()
        << "\t" << instance->getNbOfHighSpeedMachines()
        << "\t" << instance->getNbOfLowSpeedMachines()
        << "\t" << instance->getHighSpeed()
        << "\t" << instance->getLowSpeed()
        << "\t" << "LocalSearch"
        << "\t" << usePredictor
        << "\t" << std::to_string(version).c_str()
        << "\t" << useBestNeighbor
        << "\t" << time_elapsed.count()
        << "\t" << time_limits.count()
        << "\t" << firstUB
        << "\t" << nbIter
        << "\t" << maxIter
        << "\t" << solution->getSumWjUj() << std::endl;
    }
    outputFile.close();
}


