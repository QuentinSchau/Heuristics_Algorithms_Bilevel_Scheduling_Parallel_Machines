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
// Created by schau on 7/1/25.
//

#ifndef BILEVEL_SCHEDULING_BEAMSEARCH_IMP_H
#define BILEVEL_SCHEDULING_BEAMSEARCH_IMP_H

inline BeamSearch::BeamSearchNode BeamSearch::evaluation(Node &&node) {
    isWithinTimeLimit();
    // compute an upper bound
    // Get the list of available job, i.e. job that have not been scheduled or removed, to compute upper bound
    listAvailableJobForNode.clear();
    auto viewsSelectNoRemovedAndNoScheduledJob = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex());});
    std::ranges::copy(viewsSelectNoRemovedAndNoScheduledJob,std::back_inserter(listAvailableJobForNode));

    Solution::BlockStructure upperBoundSolution = heuristicSolver.constructSolutionWithHeuristic(node.getBlockStruc(),node.getIndexBlock(),node.getSelectedJobCount(), listAvailableJobForNode);

    #ifdef DEBUG_BaB
    Solution testUB = Solution(instance);
    testUB.fromBlockStruct(upperBoundSolution);
    if (not testUB.feasible(instance))
        throw BiSchException("Error in the compute of upper bound by completing an partial solution with the heuristic");
    #endif
    double upperBound = Solution::evaluate(upperBoundSolution, instance);
    if (isSmaller(firstUB,0)) {
        firstUB = upperBound;
    }
    if (upperBound < globalUB) {
        globalUB = upperBound;
        solution->fromBlockStruct(upperBoundSolution);
    }
    addSolToListBestSolutions(solution,false);
    // compute a lower bound
    double lowerBound = 0.0;
    if (isSmaller(0.0,alpha)) {
        columnGeneration.solve(node);
        lowerBound = columnGeneration.getSumWjUj();
    }
    double eval = lowerBound * alpha + (1-alpha) * upperBound;
    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << eval;
    node.stateDebug.back().append(",eval:").append(stream.str());
    #endif
    // remove the current block structure and keep only the complement
    for (unsigned int indexBlock = 0; indexBlock < node.getIndexBlock(); indexBlock++) {
        for (auto &[indexMachine,indexInMachine]: instance->getE()[indexBlock]) upperBoundSolution[indexMachine][indexInMachine] = {nullptr, 0.0};
    }

    for (auto &[indexMachine,indexInMachine]: instance->getE()[node.getIndexBlock()])
        if (node.getBlockStruc()[indexMachine][indexInMachine].first != nullptr) upperBoundSolution[indexMachine][indexInMachine] = {nullptr, 0.0};

    auto BSnode = BeamSearchNode(lowerBound,upperBound,eval,std::move(node),upperBoundSolution);

    return BSnode;
}

inline bool BeamSearch::checkRules(Node &node, unsigned int indexBlock,unsigned int indexJob, const Solution::BlockStructure &blockStruct, std::vector<Job> &listRemovedJobs, char rule) {
    auto &E = instance->getE();
    bool foundBetterInsertion = false; // boolean to know if we found a better insertion with removed job
    unsigned int bestIndexMachine = 0; // the index of the machine for the better insertion
    unsigned int bestIndexBlockInBlockStruct = 0; // the index in the best machine for the better insertion
    double maxValueSwapped = 0; // if there are several possible better insertion, we take the max value (regarding the case)
    auto &removedJob = listRemovedJobs[indexJob]; // Get the removed job from the list of removed jobs

    // Iterate over each machine and block in the current block structure
    for (auto &[indexMachine, indexBlockInBlockStruct]: E[indexBlock]) {
        auto &jobWithCj = blockStruct[indexMachine][indexBlockInBlockStruct];
        if (jobWithCj.first == nullptr) continue; // Skip this iteration if there is no job on this machine and block
        double previousCompletionTime = indexBlockInBlockStruct == 0 ? 0.0 : blockStruct[indexMachine][indexBlockInBlockStruct - 1].second;
        double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
        bool canBeSwapped = false;
        switch (rule) {
            case 1:
                // Check rule 1 conditions: scheduled job is late, removed job is on-time, and has smaller processing time
                canBeSwapped = isSmaller(jobWithCj.first->getDi(),jobWithCj.second) // scheduled is late
                               &&  isSmallerOrEqual(previousCompletionTime + removedJob.getPi()/speed,removedJob.getDi()) // removed is ontime
                               &&  isSmallerOrEqual(removedJob.getPi(),jobWithCj.first->getPi()); // removed have smaller pj{
                break;
            case 2:
                // Check rule 2 conditions: both scheduled and removed jobs are late, removed job has strict smaller weight and processing time
                canBeSwapped = isSmaller(jobWithCj.first->getDi(),jobWithCj.second)
                &&  isSmaller(removedJob.getDi(),previousCompletionTime + removedJob.getPi()/speed)
                &&  isSmallerOrEqual(removedJob.getPi(),jobWithCj.first->getPi())
                &&  isSmaller(removedJob.getWi(),jobWithCj.first->getWi());
                break;
            case 3:
                // Check rule 3 conditions: both scheduled and removed jobs are on-time, and removed job has strict smaller processing time
                canBeSwapped = isSmaller(jobWithCj.second,jobWithCj.first->getDi())
                               &&  isSmaller(previousCompletionTime + removedJob.getPi()/speed,removedJob.getDi())
                               &&  isSmaller(removedJob.getPi(),jobWithCj.first->getPi());
                break;
            default:
                canBeSwapped = false;
                break;
        }

        if (canBeSwapped){
            foundBetterInsertion = true;
            bool updateMaxSwap = false;
            switch (rule) {
                case 1:
                case 2:
                    // Update max value for rule 1 and 2: swap the job with the greatest weight
                    updateMaxSwap = isSmaller(maxValueSwapped, jobWithCj.first->getWi());
                    if (updateMaxSwap) maxValueSwapped = jobWithCj.first->getWi();
                    break;
                case 3:
                    // Update max value for rule 3: swap the job with the smallest processing time
                    updateMaxSwap = isSmaller(maxValueSwapped,jobWithCj.first->getPi());
                    if (updateMaxSwap) maxValueSwapped = jobWithCj.first->getPi();
                    break;
                default:
                    updateMaxSwap = false;
                    break;
            }

            if (updateMaxSwap) {
                bestIndexMachine = indexMachine;
                bestIndexBlockInBlockStruct = indexBlockInBlockStruct;
            }
        }
    }
    if (foundBetterInsertion) {
        Job swappedJob = *blockStruct[bestIndexMachine][bestIndexBlockInBlockStruct].first; // Get the swapped job
        node.swapRemovedAndScheduledJob(bestIndexMachine, bestIndexBlockInBlockStruct, removedJob.getIndex()); // Swap the jobs in the node
        listRemovedJobs.insert(listRemovedJobs.begin() + indexJob, swappedJob); // Insert the swapped job in the list of removed jobs
        listRemovedJobs.erase(listRemovedJobs.begin() + indexJob + 1); // Erase from the list the initially removed job
    }
    return foundBetterInsertion;
}

inline void BeamSearch::recoveringBestInsert(Node &node) {
    // we add the removed jobs from the node and try to schedule to improve the partial solution
    listAvailableJobForNode.clear();
    auto viewsSelectNoRemovedAndNoScheduledJob = instance->getListJobs() | std::views::filter([node](const Job & job) { return node.isRemoved(job.getIndex());});
    std::ranges::copy(viewsSelectNoRemovedAndNoScheduledJob,std::back_inserter(listAvailableJobForNode));

    auto & blockStruct = node.getBlockStruc();
    // we have the list of removed jobs sorted (by construction) by SPT
    for (unsigned int indexBlock = 0; indexBlock <= node.getIndexBlock(); ++indexBlock) {
        isWithinTimeLimit();
        double maxPj = -1.0, minPj = std::numeric_limits<double>::infinity();
        Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance,blockStruct, indexBlock, maxPj, minPj);
        // removed jobs are sorted by SPT, loop over them
        unsigned int indexJob = 0;
        while( indexJob < listAvailableJobForNode.size()) {
            isWithinTimeLimit();
            // if the job have processing time smaller than maxPj, than mean we have to pass to next job
            if (isSmaller(listAvailableJobForNode[indexJob].getPi(),maxPj)) ++indexJob;
            // else if the job have processing time greater than minPj, than mean all other next job will not fit in the block, then stop
            else if (isSmaller(minPj,listAvailableJobForNode[indexJob].getPi())) break;
            // else check if the job respect one rule to swap it
            else if (checkRules(node, indexBlock,indexJob, blockStruct, listAvailableJobForNode, 1)
                || checkRules(node, indexBlock,indexJob, blockStruct, listAvailableJobForNode, 2)
                || checkRules(node, indexBlock,indexJob, blockStruct, listAvailableJobForNode, 3)){
                continue;
            }else ++indexJob;
        }
    }
    node.updatePartialSumWj();
}

inline void BeamSearch::recovering(BeamSearchNode & BSNode) {
    //check if we don't have exceeded the time limit
    isWithinTimeLimit();
#ifdef DEBUG_BaB
    check_recovering(BSNode.node);
#endif
    switch (recoStrategy) {
        case BEST_INSERT:
            recoveringBestInsert(BSNode.node);
            break;
        case LOCAL_SEARCH:
                switch (version) {
                    case 1:
                        recoveringLocalSearchV1(BSNode);
                        break;
                    case 2:
                       recoveringLocalSearchV2(BSNode);
                        break;
                    case 3:
                        recoveringLocalSearchV3(BSNode);
                        break;
                    default: throw BiSchException("The version of the recovering in beam search is not implemented.");
                }
            break;
    }
#ifdef DEBUG_BaB
    check_recovering(BSNode.node);
#endif
}

inline bool BeamSearch::findBestSolutionReachableFromPartialSolution(Solution::BlockStructure &blockStructure, double & bestObjValue) {
    auto copyListAvailableJobFromNode = listAvailableJobForNode;
    // remove the jobs that are already in the block structure
    Solution::removeExistingJobsFromSolution(copyListAvailableJobFromNode, blockStructure);
    bool foundBetterNeighbor = false;
    double UB;
    bool foundBetterSol = true;
    unsigned int nbIter = 0;
    while (foundBetterSol && nbIter < maxIterLocalSearch) {
        foundBetterSol = false;
        unsigned int indexBlock = 0;
        nbIter++;
        // solve assigment problem on each block from the left to the right
        while (indexBlock < instance->getE().size()) {
            Solution::BlockStructure currentBlockStructure = blockStructure;
            heuristicSolver.freeAndAssignmentBlock(currentBlockStructure, indexBlock, nullptr);
            indexBlock++;
            if (UB = Solution::evaluate(currentBlockStructure, instance); isSmaller(UB, bestObjValue)) {
                bestObjValue = UB;
                blockStructure = std::move(currentBlockStructure);
                foundBetterNeighbor = true;
                foundBetterSol = true;
            }
        }
        if (foundBetterSol) {
            foundBetterSol = false;
            // solve assigment problem on each block from the right to the left
            while (indexBlock-->0){
                assert(indexBlock < instance->getE().size());
                Solution::BlockStructure currentBlockStructure = blockStructure;
                heuristicSolver.freeAndAssignmentBlock(currentBlockStructure, indexBlock, nullptr);
                if (UB = Solution::evaluate(currentBlockStructure, instance); isSmaller(UB, bestObjValue)) {
                    bestObjValue = UB;
                    blockStructure = std::move(currentBlockStructure);
                    foundBetterNeighbor = true;
                    foundBetterSol = true;
                }
            }
        }
    }
    #ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
    testSol.fromBlockStruct(blockStructure);
    assert(testSol.feasible(instance));
    #endif
    return foundBetterNeighbor;
}

inline void BeamSearch::recoveringLocalSearchV1(BeamSearchNode & BSnode){
    Node &node = BSnode.node;
    isWithinTimeLimit();
    // Get the list of available job, i.e. job that have not been scheduled or removed, to compute upper bound
    listAvailableJobForNode.clear();
    auto viewsSelectNoRemovedAndNoScheduledJob = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex());});
    std::ranges::copy(viewsSelectNoRemovedAndNoScheduledJob,std::back_inserter(listAvailableJobForNode));

    // use copy of the list of available job from the node because the heuristic alter it
    std::vector<Job> copyListAvailableJobFromNode;
    copyListAvailableJobFromNode.reserve(listAvailableJobForNode.size());

    // Get the list of job's index that have been scheduled, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listJobsSelected;
    listJobsSelected.reserve(node.getSelectedJobCount());
    auto viewsSelectScheduledJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isScheduled(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectScheduledJobs,std::back_inserter(listJobsSelected));

    // use a vector of boolean, true if a job is selected false otherwise
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    for (unsigned int indexJobSelected : viewsSelectScheduledJobs) {
        alreadySelectedJobs[indexJobSelected] = true;
    }

    // Get the list of job's index that have been removed, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listRemovedJobs;
    listRemovedJobs.reserve(node.getRemovedJobsCount());
    auto viewsSelectRemovedJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isRemoved(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectRemovedJobs,std::back_inserter(listRemovedJobs));

    bool haveFindBestSol = true;
    double LB; // the lower bound of the node
    unsigned int nbBestSolFind = 0; // the number of times we found a best solution in the neighborhood
    unsigned int nbIter = 0; // the number of iteration of the method
    isWithinTimeLimit();
    while (haveFindBestSol) {
        haveFindBestSol = false;

        //update LB if we have found a better solution or if alpha =0.0 (we don't have compute LB in evaluation
        if (nbBestSolFind > 0 || isSmallerOrEqual(alpha,0.0)) {
            columnGeneration.setTimeElapsed(time_elapsed);
            try{
                columnGeneration.solve(node);
            } catch (const BiSchTimeOutException &e) {
                this->setTimeElapsed(columnGeneration.getTimeElapsed());
                throw BiSchTimeOutException(e);
            }
            this->setTimeElapsed(columnGeneration.getTimeElapsed());
            LB = columnGeneration.getSumWjUj();
        }else LB = BSnode.lowerbound;
        // define the best UB
        double bestUB = BSnode.upperbound;

        Solution::BlockStructure bestBlockStructure; // keep the best block structure
        Solution::BlockStructure blockStructFromNode = node.getBlockStruc(); // get the block structure of the node

        /** Start local search **/

        bool foundBetterSol = false;
        // if (node.id == 68)
            // std::cout << "";

        #ifdef DEBUG_HEURISTIC
        Solution testSol(instance);
        testSol.fromBlockStruct(node.getBlockStruc());
        auto testListJob = testSol.extractListOfJobs();
        auto viewTest = std::views::reverse( testListJob | std::views::transform([](Job &job){return job.getIndex();}));
        std::vector<unsigned int> testListJobIndices(viewTest.begin(),viewTest.end());
        #endif
        // keep the best swap to adjust the list of selected jobs, the pair is (index of the job to remove, index of the job to insert)
        std::pair<unsigned int,unsigned int> bestSwap;
        while (nbIter < maxIterLocalSearch) {
            isWithinTimeLimit();
            nbIter++;
            foundBetterSol = false;
            if (verbose >= 3) {
                std::cout << "NbIter:" << nbIter << " best UB:" << bestUB << std::endl;
            }

            // loop over selected jobs in reverse order, i.e. in LPT to try to swap first the greater pj
            for (unsigned int indexRemovedJobFromSelection : std::ranges::reverse_view{listJobsSelected}) {
                isWithinTimeLimit();
                // keep the index of the job that we remove in the selection
                alreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
                // now we try to select one of the removed jobs
                for (unsigned int indexJobNotSelected : listRemovedJobs) {
                    isWithinTimeLimit();
                    assert(indexRemovedJobFromSelection != indexJobNotSelected);
                    assert(not alreadySelectedJobs[indexJobNotSelected]);
                    // if jobs we want to swap not belonging to same identical group of jobs, because we use a exact algorithm to manage this jobs.
                    if (instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                        // select the removed job
                        alreadySelectedJobs[indexJobNotSelected] = true;
                        if (indexRemovedJobFromSelection == 22 && indexJobNotSelected == 20)
                            std::cout << "";
                        // we swap the job and update the block structure, leading to newBlockStructure
                        auto [newBlockStructure,newObjValue] = localSearchSolver.updateBlockStructureWithSwap(listJobsSelected, blockStructFromNode, indexJobNotSelected, indexRemovedJobFromSelection);
                        isWithinTimeLimit();
                        copyListAvailableJobFromNode = listAvailableJobForNode;
                        // we fill the partial schedule with a heuristic to get a feasible solution thus an upper bound
                        newBlockStructure = heuristicSolver.constructSolutionWithHeuristic(newBlockStructure, node.getIndexBlock(),node.getSelectedJobCount(), copyListAvailableJobFromNode);
                        newObjValue = Solution::evaluate(newBlockStructure, instance);
                        #ifdef DEBUG_HEURISTIC
                        testSol.fromBlockStruct(newBlockStructure);
                        if (not testSol.feasible(instance))
                            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                        #endif

                        if (isSmaller(newObjValue, bestUB)) {
                            foundBetterSol = true;
                            bestBlockStructure = std::move(newBlockStructure);
                            bestUB = newObjValue;
                            bestSwap = {indexRemovedJobFromSelection,indexJobNotSelected};
                            assert(bestUB == Solution::evaluate(bestBlockStructure,instance));
                            // try to improve solution with assigment the bestBlockStructure, i.e. try follower neighborhood on the leader decisions
                            findBestSolutionReachableFromPartialSolution(bestBlockStructure,bestUB);
                        }
                        // try to improve solution (newBlockStructure) with assigment, i.e. try follower neighborhood on the leader decisions
                        else if (findBestSolutionReachableFromPartialSolution(newBlockStructure,bestUB)) {
                            bestBlockStructure = std::move(newBlockStructure);
                            assert(bestUB == Solution::evaluate(bestBlockStructure,instance));
                            foundBetterSol = true;
                            bestSwap = {indexRemovedJobFromSelection,indexJobNotSelected};
                        }
                        if (indexRemovedJobFromSelection == 4 && indexJobNotSelected == 0)
                            std::cout << "";
                        // reset the initial state
                        alreadySelectedJobs[indexJobNotSelected] = false;
                        LocalSearch::swapJobsInSelection(listJobsSelected, indexRemovedJobFromSelection, indexJobNotSelected);
                        #ifdef DEBUG_HEURISTIC
                        // check if the listJobsSelected is reset
                        assert(std::ranges::equal(listJobsSelected,testListJobIndices));
                        #endif

                    }
                }
                // re-select the job
                alreadySelectedJobs[indexRemovedJobFromSelection] = true;
            }
            if (not foundBetterSol) break;
        }
        if (isSmaller(bestUB,LB)) {
            // first update the solution if the UB is better than the global UB
            if (isSmaller(bestUB,globalUB)) {
                this->solution->fromBlockStruct(bestBlockStructure);
            }
            #ifdef DEBUG_HEURISTIC
            testSol.fromBlockStruct(bestBlockStructure);
            testListJob = testSol.extractListOfJobs();
            testListJobIndices = std::vector<unsigned int>(viewTest.begin(),viewTest.begin() + node.getSelectedJobCount());
            #endif
            // update the list of removed jobs, the list of selected jobs and the vector of boolean
            LocalSearch::swapJobsInSelection(listJobsSelected, bestSwap.second, bestSwap.first);
            LocalSearch::swapJobsInSelection(listRemovedJobs,  bestSwap.first,bestSwap.second);
            alreadySelectedJobs[bestSwap.first] = false;
            alreadySelectedJobs[bestSwap.second] = true;
            haveFindBestSol = true;
            nbBestSolFind++;
            //we keep the schedule in bestBlockStructure until node.getIndexBlock()
            // then for the last block in the block structure of the node, we loop over it and keep jobs with smaller processing time than node.getCurrentPj()

            for (unsigned int indexBlock = node.getIndexBlock(); indexBlock < instance->getE().size(); indexBlock++) {
                for (auto &[indexMachine, indexInMachine] : instance->getE()[indexBlock]) {
                    if (bestBlockStructure[indexMachine][indexInMachine].first != nullptr && isSmallerOrEqual(node.getCurrentPj(),bestBlockStructure[indexMachine][indexInMachine].first->getPi())) {
                        bestBlockStructure[indexMachine][indexInMachine].first = nullptr;
                        bestBlockStructure[indexMachine][indexInMachine].second = 0.0;
                    }
                }
            }
            node.setBlockStructure(std::move(bestBlockStructure));
#ifdef DEBUG_BaB
            try {
                check_recovering(node);
            }catch (BiSchException &e) {
                Solution::printB(node.getBlockStruc());
                throw BiSchException(e);
            }
#endif
        }
    }
    NB_reco++;
    NB_best_sol_find_reco += nbBestSolFind;
}

inline void BeamSearch::recoveringLocalSearchV2(BeamSearchNode & BSnode){
    Node &node = BSnode.node;
    isWithinTimeLimit();
    // Get the list of available job, i.e. job that have not been scheduled or removed, to compute upper bound
    listAvailableJobForNode.clear();
    auto viewsSelectNoRemovedAndNoScheduledJob = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex());});
    std::ranges::copy(viewsSelectNoRemovedAndNoScheduledJob,std::back_inserter(listAvailableJobForNode));

    // use copy of the list of available job from the node because the heuristic alter it
    std::vector<Job> copyListAvailableJobFromNode;
    copyListAvailableJobFromNode.reserve(listAvailableJobForNode.size());

    // Get the list of job's index that have been scheduled, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listJobsSelected;
    listJobsSelected.reserve(node.getSelectedJobCount());
    auto viewsSelectScheduledJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isScheduled(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectScheduledJobs,std::back_inserter(listJobsSelected));

    // use a vector of boolean, true if a job is selected false otherwise
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    for (unsigned int indexJobSelected : viewsSelectScheduledJobs) {
        alreadySelectedJobs[indexJobSelected] = true;
    }

    // Get the list of job's index that have been removed, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listRemovedJobs;
    listRemovedJobs.reserve(node.getRemovedJobsCount());
    auto viewsSelectRemovedJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isRemoved(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectRemovedJobs,std::back_inserter(listRemovedJobs));
    isWithinTimeLimit();
    // define the best UB
    double bestUB = BSnode.upperbound;

    Solution::BlockStructure bestBlockStructure = node.getBlockStruc(); // keep the best block structure
    Solution::BlockStructure blockStructFromNode = node.getBlockStruc(); // get the block structure of the node
    // sort the block structure by makespan
    Solution::sortBlockStructurePerMakespan(blockStructFromNode, instance);

    // if (node.id == 68)
    //     std::cout << "";

    #ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
    testSol.fromBlockStruct(node.getBlockStruc());
    auto testListJob = testSol.extractListOfJobs();
    auto viewTest = std::views::reverse( testListJob | std::views::transform([](Job &job){return job.getIndex();}));
    std::vector<unsigned int> testListJobIndices(viewTest.begin(),viewTest.end());
    #endif
    // loop over selected jobs in reverse order, i.e. in LPT to try to swap first the greater pj
    for (unsigned int indexRemovedJobFromSelection : std::ranges::reverse_view{listJobsSelected}) {
        isWithinTimeLimit();
        // keep the index of the job that we remove in the selection
        alreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
        // now we try to select one of the removed jobs
        for (unsigned int indexJobNotSelected : listRemovedJobs) {
            isWithinTimeLimit();
            assert(indexRemovedJobFromSelection != indexJobNotSelected);
            assert(not alreadySelectedJobs[indexJobNotSelected]);
            // if jobs we want to swap not belonging to same identical group of jobs, because we use a exact algorithm to manage this jobs.
            if (instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                // select the removed job
                alreadySelectedJobs[indexJobNotSelected] = true;
                if (indexRemovedJobFromSelection == 22 && indexJobNotSelected == 20)
                    std::cout << "";
                // we swap the job and update the block structure, leading to newBlockStructure
                auto [newBlockStructure,newObjValue] = localSearchSolver.updateBlockStructureWithSwap(listJobsSelected, blockStructFromNode, indexJobNotSelected, indexRemovedJobFromSelection);
                unsigned int indexBlock = 0;
                double UB;
                Solution::BlockStructure bestBlockStructureWithPermutation = newBlockStructure;
                bool foundBetterSol = false;
                // use neighborhood on the follower decisions, i.e. try all permutations
                while (indexBlock <= node.getIndexBlock()) {
                    isWithinTimeLimit();
                    Neighborhoods neighborhoods(instance, &bestBlockStructureWithPermutation, &listJobsSelected, 2);
                    auto itSol = neighborhoods.beginFN();
                    auto itEndNeighbor = neighborhoods.endFN();
                    for (; itSol != itEndNeighbor; ++itSol) {
                        isWithinTimeLimit();
                        // check if we have not pass to the next block with the neighborhood
                        if (indexBlock < itSol.getIndexBlock()) { break; }
                        Solution::BlockStructure blockStructureFromNeighbor = *itSol;
                        copyListAvailableJobFromNode = listAvailableJobForNode;
                        // we fill the partial schedule with a heuristic to get a feasible solution thus an upper bound
                        auto newBlockStructureFromNeighbor = heuristicSolver.constructSolutionWithHeuristic(blockStructureFromNeighbor, node.getIndexBlock(),node.getSelectedJobCount(), copyListAvailableJobFromNode);
                        UB = Solution::evaluate(newBlockStructureFromNeighbor, instance);
                        #ifdef DEBUG_HEURISTIC
                        testSol.fromBlockStruct(newBlockStructureFromNeighbor);
                        if (not testSol.feasible(instance))
                            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                        #endif
                        // if the permutation reduce the weighted number of tardy jobs
                        if (isSmaller(UB, bestUB)) {
                            bestBlockStructureWithPermutation = std::move(blockStructureFromNeighbor);
                            bestUB = UB;
                            foundBetterSol =true;
                        }
                        else if (isEqual(UB,bestUB)) {
                            Solution::sortBlockStructurePerMakespan(blockStructureFromNeighbor,instance);
                            bool isBetter = true; // the neighbor is better if it has the minimal makespan
                            unsigned int nbEgalMakespan = 0;
                            for (auto &[indexMachine,indexInMachine] : instance->getE()[node.getIndexBlock()]){
                                if (isSmaller(bestBlockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second)) {
                                    isBetter = false;
                                }else if (isEqual(bestBlockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second))
                                    ++nbEgalMakespan;
                            }
                            if (isBetter and nbEgalMakespan < instance->getE()[node.getIndexBlock()].size()) {
                                bestBlockStructureWithPermutation = std::move(blockStructureFromNeighbor);
                                foundBetterSol = true;
                            }
                        }
                    }
                    indexBlock++;
                }
                if (foundBetterSol) {
                    bestBlockStructure = bestBlockStructureWithPermutation;
                    NB_best_sol_find_reco++;
                }
                if (indexRemovedJobFromSelection == 4 && indexJobNotSelected == 0)
                    std::cout << "";
                // reset the initial state
                alreadySelectedJobs[indexJobNotSelected] = false;
                LocalSearch::swapJobsInSelection(listJobsSelected, indexRemovedJobFromSelection, indexJobNotSelected);
                #ifdef DEBUG_HEURISTIC
                // check if the listJobsSelected is reset
                assert(std::ranges::equal(listJobsSelected,testListJobIndices));
                #endif
            }
        }
        // re-select the job
        alreadySelectedJobs[indexRemovedJobFromSelection] = true;
    }

    node.setBlockStructure(std::move(bestBlockStructure));
#ifdef DEBUG_BaB
    try {
        check_recovering(node);
    }catch (BiSchException &e) {
        Solution::printB(node.getBlockStruc());
        throw BiSchException(e);
    }
#endif
    NB_reco++;
}

inline void BeamSearch::recoveringLocalSearchV3(BeamSearchNode & BSnode){
    Node &node = BSnode.node;

    // Get the list of available job, i.e. job that have not been scheduled or removed, to compute upper bound
    listAvailableJobForNode.clear();
    auto viewsSelectNoRemovedAndNoScheduledJob = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex());});
    std::ranges::copy(viewsSelectNoRemovedAndNoScheduledJob,std::back_inserter(listAvailableJobForNode));

    // use copy of the list of available job from the node because the heuristic alter it
    std::vector<Job> copyListAvailableJobFromNode;
    copyListAvailableJobFromNode.reserve(listAvailableJobForNode.size());

    // Get the list of job's index that have been scheduled, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listJobsSelected;
    listJobsSelected.reserve(node.getSelectedJobCount());
    auto viewsSelectScheduledJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isScheduled(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectScheduledJobs,std::back_inserter(listJobsSelected));

    // use a vector of boolean, true if a job is selected false otherwise
    std::vector<bool> alreadySelectedJobs(instance->getNbJobs(), false);
    for (unsigned int indexJobSelected : viewsSelectScheduledJobs) {
        alreadySelectedJobs[indexJobSelected] = true;
    }

    // Get the list of job's index that have been removed, we use it in the exploration of the neighborhood
    std::vector<unsigned int> listRemovedJobs;
    listRemovedJobs.reserve(node.getRemovedJobsCount());
    auto viewsSelectRemovedJobs = instance->getListJobs()
    | std::views::filter([node](const Job & job) { return node.isRemoved(job.getIndex());})
    | std::views::transform([](const Job & job){return job.getIndex();});
    std::ranges::copy(viewsSelectRemovedJobs,std::back_inserter(listRemovedJobs));

    // define the best UB
    double bestUB = BSnode.upperbound;
    isWithinTimeLimit();
    Solution::BlockStructure bestBlockStructure = node.getBlockStruc(); // keep the best block structure
    Solution::BlockStructure blockStructFromNode = node.getBlockStruc(); // get the block structure of the node
    // sort the block structure by makespan
    Solution::sortBlockStructurePerMakespan(blockStructFromNode, instance);

    // if (node.id == 136)
    //     std::cout << "";

    #ifdef DEBUG_HEURISTIC
    Solution testSol(instance);
    testSol.fromBlockStruct(node.getBlockStruc());
    auto testListJob = testSol.extractListOfJobs();
    auto viewTest = std::views::reverse( testListJob | std::views::transform([](Job &job){return job.getIndex();}));
    std::vector<unsigned int> testListJobIndices(viewTest.begin(),viewTest.end());
    #endif
    // loop over selected jobs in reverse order, i.e. in LPT to try to swap first the greater pj
    for (unsigned int indexRemovedJobFromSelection : std::ranges::reverse_view{listJobsSelected}) {
        isWithinTimeLimit();
        // keep the index of the job that we remove in the selection
        alreadySelectedJobs[indexRemovedJobFromSelection] = false; // make the job unselected
        // now we try to select one of the removed jobs
        for (unsigned int indexJobNotSelected : listRemovedJobs) {
            isWithinTimeLimit();
            assert(indexRemovedJobFromSelection != indexJobNotSelected);
            assert(not alreadySelectedJobs[indexJobNotSelected]);
            // if jobs we want to swap not belonging to same identical group of jobs, because we use an exact algorithm to manage this jobs.
            if (instance->getIndexIdenticalGroupOfJob(indexRemovedJobFromSelection) != instance->getIndexIdenticalGroupOfJob(indexJobNotSelected)) {
                // select the removed job
                alreadySelectedJobs[indexJobNotSelected] = true;
                if (indexRemovedJobFromSelection == 48 && indexJobNotSelected == 4)
                    std::cout << "";
                // we swap the job and update the block structure, leading to newBlockStructure
                auto [newBlockStructure,newObjValue] = localSearchSolver.updateBlockStructureWithSwap(listJobsSelected, blockStructFromNode, indexJobNotSelected, indexRemovedJobFromSelection);
                unsigned int indexBlock = 0;
                double UB;
                Solution::BlockStructure bestBlockStructureWithPermutation = newBlockStructure;
                bool foundBetterSol = false;
                // use neighborhood on the follower decisions, i.e. try all permutations
                while (indexBlock < node.getIndexBlock()) {
                    isWithinTimeLimit();
                    Neighborhoods neighborhoods(instance, &bestBlockStructureWithPermutation, &listJobsSelected, 2);
                    auto itSol = neighborhoods.beginFN();
                    auto itEndNeighbor = neighborhoods.endFN();
                    for (; itSol != itEndNeighbor; ++itSol) {
                        isWithinTimeLimit();
                        // check if we have not pass to the next block with the neighborhood
                        if (indexBlock < itSol.getIndexBlock()) { break; }
                        Solution::BlockStructure blockStructureFromNeighbor = *itSol;
                        Solution::BlockStructure newBlockStructureFromNeighbor = blockStructureFromNeighbor;
                        // we fill the partial schedule with a heuristic to get a feasible solution thus an upper bound
                        Solution::mergeTwoBlockStructure(newBlockStructureFromNeighbor, BSnode.completingSchedule);

                        UB = Solution::evaluate(newBlockStructureFromNeighbor, instance);
                        addSolToListBestSolutions(newBlockStructureFromNeighbor, UB);
                        #ifdef DEBUG_HEURISTIC
                        testSol.fromBlockStruct(newBlockStructureFromNeighbor);
                        if (not testSol.feasible(instance))
                            throw BiSchException("Error in compute leader neighborhood, block structure not feasible");
                        #endif
                        // if the permutation reduce the weighted number of tardy jobs
                        if (isSmaller(UB, bestUB)) {
                            bestBlockStructureWithPermutation = std::move(blockStructureFromNeighbor);
                            bestUB = UB;
                            foundBetterSol =true;
                        }
                        else if (isEqual(UB,bestUB)) {
                            Solution::sortBlockStructurePerMakespan(blockStructureFromNeighbor,instance);
                            bool isBetter = true; // the neighbor is better if it has the minimal makespan
                            unsigned int nbEgalMakespan = 0;
                            for (auto &[indexMachine,indexInMachine] : instance->getE()[node.getIndexBlock()]){
                                if (isSmaller(bestBlockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second)) {
                                    isBetter = false;
                                }else if (isEqual(bestBlockStructure[indexMachine][indexInMachine].second,blockStructureFromNeighbor[indexMachine][indexInMachine].second))
                                    ++nbEgalMakespan;
                            }
                            if (isBetter and nbEgalMakespan < instance->getE()[node.getIndexBlock()].size()) {
                                bestBlockStructureWithPermutation = std::move(blockStructureFromNeighbor);
                                foundBetterSol = true;
                            }
                        }
                    }
                    indexBlock++;
                }
                if (foundBetterSol) {
                    bestBlockStructure = bestBlockStructureWithPermutation;
                    NB_best_sol_find_reco++;
                }
                if (indexRemovedJobFromSelection == 4 && indexJobNotSelected == 0)
                    std::cout << "";
                // reset the initial state
                alreadySelectedJobs[indexJobNotSelected] = false;
                LocalSearch::swapJobsInSelection(listJobsSelected, indexRemovedJobFromSelection, indexJobNotSelected);
                #ifdef DEBUG_HEURISTIC
                // check if the listJobsSelected is reset
                assert(std::ranges::equal(listJobsSelected,testListJobIndices));
                #endif
            }
        }
        // re-select the job
        alreadySelectedJobs[indexRemovedJobFromSelection] = true;
    }
    node.setBlockStructure(std::move(bestBlockStructure));
#ifdef DEBUG_BaB
    try {
        check_recovering(node);
    }catch (BiSchException &e) {
        Solution::printB(node.getBlockStruc());
        throw BiSchException(e);
    }
#endif
    NB_reco++;
}

#ifdef DEBUG_BaB
inline void BeamSearch::check_recovering(Node &node) {
    // check there are no duplicate scheduled jobs
    std::vector<Job> listJobsScheduled;
    auto &blockStructure = node.getBlockStruc();
    for (unsigned int indexMachine = 0; indexMachine < blockStructure.size(); ++indexMachine) {
        auto &machine = node.getBlockStruc()[indexMachine];
        double completionTimeOfMachine = 0.0;
        double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
        char machineSpeed = indexMachine < instance->getNbOfHighSpeedMachines() ? 0 : 1;
        unsigned int nbJobsScheduled = 0;
        for (auto &jobWithCj: machine) {
            if (jobWithCj.first != nullptr) {
                nbJobsScheduled++;
                auto job = *jobWithCj.first;
                listJobsScheduled.push_back(job);
                completionTimeOfMachine += job.getPi() / speed;
                if (not isEqual(completionTimeOfMachine, jobWithCj.second)) {
                    std::string message("Error in recovering: the completion time of jobs J");
                    message.append(std::to_string(job.getIndex())).append(" is ")
                            .append(std::to_string(jobWithCj.second)).append(" instead of ")
                            .append(std::to_string(completionTimeOfMachine));
                    throw BiSchException(message);
                }
                if (not node.isScheduled(job.getIndex())) throw BiSchException("Error in recovering: the encoding of selected jobs is not correct.");
                if (node.isScheduledOnOtherMachines(job.getIndex(), machineSpeed))
                    throw BiSchException("Error in recovering: the encoding of machine is not correct. A job is encoded on several machines");
                if (node.isRemoved(job.getIndex())) throw BiSchException("Error in recovering: the encoding of removed jobs is not correct. A job is encoded as removed whereas it is scheduled");
            }
        }
        if (node.getEncodingSelectedJobOnMachine()[indexMachine].count() != nbJobsScheduled)
            throw BiSchException("The number of job encoded on a machine is different that the one on the block structure");
    }
    std::sort(listJobsScheduled.begin(), listJobsScheduled.end());
    auto ifFindDuplicate = std::adjacent_find(listJobsScheduled.begin(), listJobsScheduled.end());
    if (ifFindDuplicate != listJobsScheduled.end()) throw BiSchException("Error in recovering: there are duplicate scheduled jobs in the block structure of the node.");
    if (listJobsScheduled.size() != node.getSelectedJobCount())
        throw BiSchException("Error in recovering: The number of job in block structure is not the same as the number of jobs counted in the node");
    if (node.getEncodingRemoveDecision().count() != node.getRemovedJobsCount())
        throw BiSchException("Error in recovering: The number of removed job in the encoding is not the same as the number of jobs counted in the node");
}
#endif
inline void BeamSearch::createChildrenNodeWithIdenticalJobs(Node &pNode) {
    // case where we take no jobs:
    Node cNode = Node(pNode);
    #ifdef DEBUG_BaB
    std::stringstream stream;
    stream << std::fixed << std::setprecision(0) << cNode.getCurrentPj();
    std::string name = std::string("p=").append(stream.str()).append(
            ",|s|=").append(std::to_string(cNode.getSelectedJobCount())).append(",L=").append(std::to_string(cNode.getIndexBlock())).append(",X");
    cNode.stateDebug.emplace_back(name);
    #endif
    unsigned int numberIdenticalJobsInGroups = cNode.removeGroupOfIdenticalJobs();
    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    addNode(cNode);
    #if defined DEBUG_BaB && defined DEBUG_DOT
    // if cNode have not the same id as its parents, means we add it to the list
    if (cNode.id != pNode.id)
        dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
            << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << "\"];"
            << std::endl;
    #endif

    // take at least one jobs
    for (unsigned int numberOfJobToSchedule = 1; numberOfJobToSchedule <= numberIdenticalJobsInGroups; ++numberOfJobToSchedule) {
        isWithinTimeLimit();
        // check if we respect the follower problem, i.e., we must select n jobs
        if (numberOfJobToSchedule + pNode.getSelectedJobCount() > instance->getNbToSelectJob()) break;
        // n_L the number of selected machine in the last block
        auto [indexLastBlock, n_L] = computeNumberJobsOnLastBlock(pNode, numberOfJobToSchedule);
        std::set<std::vector<unsigned int>> setOfAssignmentOnFirstBlock;
        // if the first block is the block 0, we need to fill it, that means we need to use different assignment for this block because it's the only one where there
        // the number of available machines a_0 > eta_0 the number of machines we must select. We have exactly a_0 + eta_0 - instance->getE()[0] available location
        if (pNode.getIndexBlock() == 0 and indexLastBlock > 0) {
            unsigned int numberAvailablePosition = pNode.getAvailableLocations() + pNode.getNumJobsToScheduleOnBlock() - instance->getE()[0].size();
            auto assignmentFirstBlock = std::vector<unsigned int>({});
            computeSetCombinationWithoutSymmetry(numberAvailablePosition, 0, pNode, assignmentFirstBlock, setOfAssignmentOnFirstBlock);
        } else setOfAssignmentOnFirstBlock.emplace();

        for (auto &assignmentFirstBlock: setOfAssignmentOnFirstBlock) {
            isWithinTimeLimit();
            std::set<std::vector<unsigned int>> setOfAssignmentOnLastBlock;
            computeSetCombinationWithoutSymmetry(n_L, indexLastBlock, pNode, assignmentFirstBlock, setOfAssignmentOnLastBlock);
            for (auto &assignmentLastBlock: setOfAssignmentOnLastBlock) {
                isWithinTimeLimit();
                createNodeWithAssignment(pNode, indexLastBlock, assignmentLastBlock, assignmentFirstBlock);
            }

        }
    }

}

inline std::pair<unsigned int, unsigned int> BeamSearch::computeNumberJobsOnLastBlock(Node &pNode, unsigned int numberJobsToSchedule) {
    unsigned int startIndexBlock = pNode.getIndexBlock();
    std::vector<unsigned int> availableMachines;
    availableMachines.reserve(numberJobsToSchedule);
    unsigned int numJobsToScheduleOnBlock = pNode.getNumJobsToScheduleOnBlock(); // number of machine that must be selected in the block k
    unsigned int availableLocations = pNode.getAvailableLocations(); // number of available machine in block k
    unsigned int n_L = 0; // the number of job in last block
    while (numberJobsToSchedule > 0) {
        unsigned int numberAvailablePosition = availableLocations + numJobsToScheduleOnBlock - instance->getE()[startIndexBlock].size();
        if (numberJobsToSchedule > numberAvailablePosition) {
            numberJobsToSchedule -= numberAvailablePosition;
            ++startIndexBlock;
            numJobsToScheduleOnBlock = instance->getE()[startIndexBlock].size();
            availableLocations = instance->getE()[startIndexBlock].size();
        } else {
            n_L = numberJobsToSchedule;
            numberJobsToSchedule = 0;
        }
    }
    return {startIndexBlock, n_L};
}

inline void BeamSearch::computeSetCombinationWithoutSymmetry(unsigned int nbSelect, unsigned int indexLastBlock, Node &node, const std::vector<unsigned int> &assigmentOnStartingBlock
                                                             , std::set<std::vector<unsigned int>> &setOfAssignmentOnLastBlock) {
    isWithinTimeLimit();
    // create vector for the assignment on the first block
    std::vector<unsigned int> assignmentOnFirstBlock;
    assignmentOnFirstBlock.reserve(instance->getE()[node.getIndexBlock()].size());
    // get the processing time of identical jobs
    double pj = node.getCurrentPj();
    //if not assigment was defined for the first block, then construct it with available machine
    if (assigmentOnStartingBlock.empty()) {
        for (auto &location: instance->getE()[node.getIndexBlock()]) {
            // if not jobs where schedule, that means we have null value for the pointer
            if (!node.getBlockStruc()[location.first][location.second].first)
                assignmentOnFirstBlock.push_back(location.first);
        }
    } else
        assignmentOnFirstBlock.insert(assignmentOnFirstBlock.end(), assigmentOnStartingBlock.begin(), assigmentOnStartingBlock.end());


    // list of completion time of available machine use to compute the set of combinations
    std::vector<unsigned int> availableCompletionTime;
    // We store the list of elegant pair with completion time and the speed (it's bijection N^2 -> N)
    availableCompletionTime.reserve(instance->getE()[indexLastBlock].size());
    // use a map, where the key is elegant pair with completion time and the speed and the value is a vector of
    // index machine. The idea is to keep with this map the machines with same completion time on same speed.
    std::unordered_map<unsigned int, std::vector<unsigned int>> indexMachinesGroupByCompletionTime;
    for (auto &location: instance->getE()[indexLastBlock]) {
        isWithinTimeLimit();
        // location = (index Machine, index Block)
        // in (Job, Ci) in block structure /!\ Ci in the blockStructure are already divided by the speed
        double speed = (location.first < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
        auto locationAssigned = node.getBlockStruc()[location.first][location.second];
        // if we have null value for the pointer that means no job was schedule, in other words the location is available
        if (!locationAssigned.first) {
            // compute the completion time of the machine.
            unsigned int Ci;
            // Check if we have several blocks.
            if (indexLastBlock > 0) {
                // compute the number of block that was filled,
                unsigned int numberOfFilledBlock = indexLastBlock - node.getIndexBlock() + 1;
                // if the machine was not present in the first block, then we have one less block to fulfill
                if (std::find(assignmentOnFirstBlock.begin(), assignmentOnFirstBlock.end(), location.first) == assignmentOnFirstBlock.end())
                    --numberOfFilledBlock;
                // We need to get the last block index where the machine is present and where there is a scheduled job that is
                unsigned int indexBlockWithSameMachine = indexLastBlock - 1;
                auto predFindIndexMachine = [&location](std::pair<unsigned int, unsigned int> locationInBlock) {
                    return locationInBlock.first == location.first;
                };
                bool foundLastBlockIndex = false;
                // create a location with index last block equal to 0
                auto locationLastScheduledJob = std::pair(location.first, 0);
                while (!foundLastBlockIndex && indexBlockWithSameMachine > 0) {
                    // looking for the machine's index in the last block
                    auto itLocationInBlockWithSameMachine = std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine);
                    // if we have machine's index in the block, and we have a job that is already scheduled
                    if (itLocationInBlockWithSameMachine != instance->getE()[indexBlockWithSameMachine].end() &&
                        node.getBlockStruc()[itLocationInBlockWithSameMachine->first][itLocationInBlockWithSameMachine->second].first) {
                        locationLastScheduledJob = *itLocationInBlockWithSameMachine;
                        foundLastBlockIndex = true;
                    } else {
                        // go to the block before
                        --indexBlockWithSameMachine;
                    }
                }

                Ci = static_cast<unsigned int>(node.getBlockStruc()[locationLastScheduledJob.first][locationLastScheduledJob.second].second * speed +
                                               numberOfFilledBlock * pj); // get the completion time
            } else Ci = static_cast<unsigned int>(pj);
            // the key is the bijection N^2 -> N of (Ci, speed)
            unsigned int keyCi = (location.first < instance->getNbOfHighSpeedMachines()) ? elegantPair(Ci, 0u) : elegantPair(Ci, 1u);
            auto itIndexMachineGroupByCi = indexMachinesGroupByCompletionTime.find(keyCi);
            if (itIndexMachineGroupByCi == indexMachinesGroupByCompletionTime.end())
                indexMachinesGroupByCompletionTime.insert({keyCi, {location.first}});
            else itIndexMachineGroupByCi->second.emplace_back(location.first);
            availableCompletionTime.push_back(keyCi);
        }
    }

    //sort completion time in order to don't have duplicate combination
    std::sort(availableCompletionTime.begin(), availableCompletionTime.end());
    // compute the set of combinations with nbSelect elements
    std::set<std::vector<unsigned int>> setOfCombinationFromCompletionTime;
    findCombinations(availableCompletionTime, nbSelect, setOfCombinationFromCompletionTime);
    // Each combination contains elegant pair with completion time and the speed, we just want a combination of machine's index
    for (auto &combination: setOfCombinationFromCompletionTime) {
        isWithinTimeLimit();
        // we use 'indexListMachineWithCi' that correspond of the index in the list of machine's index
        unsigned int indexListMachineWithCi = 0;
        std::vector<unsigned int> newCombinationWithMachineIndex(combination);
        for (unsigned int indexInCombination = 0; indexInCombination < combination.size(); ++indexInCombination) {
            // if it's not the first element and the completion time is same that the previous element, we increase the index in the list of machine's index
            if (indexInCombination > 0 && combination[indexInCombination - 1] == combination[indexInCombination])
                ++indexListMachineWithCi;
            else
                // else we reset the index to the first element, i.e. 0
                indexListMachineWithCi = 0;
            // we get the real machine's index, by finding the completion time and get the index machine stored in the position 'indexListMachineWithCi'
            unsigned int indexMachine = indexMachinesGroupByCompletionTime.find(combination[indexInCombination])->second[indexListMachineWithCi];
            // change the completion time in combination by the machine's index
            newCombinationWithMachineIndex[indexInCombination] = indexMachine;
        }
        // add the new combination
        setOfAssignmentOnLastBlock.insert(newCombinationWithMachineIndex);
    }
}

inline void
BeamSearch::createNodeWithAssignment(Node &pNode, unsigned int indexLastBlock, const std::vector<unsigned int> &assignmentLastBlock, const std::vector<unsigned int> &assignmentFirstBlock) {
    Node cNode = Node(pNode);

    // get the processing time of identical jobs
    double pj = cNode.getCurrentPj();

    // declare the binary tree of completion times that is not assigned with the corresponding position,
    // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
    std::multimap<double, std::pair<unsigned int, unsigned int>> listCjAndAvailablePosition;

    if (!assignmentFirstBlock.empty()) {
        // compute completion of the first block and from assignmentFirstBlock std::vector
        for (unsigned int indexMachine: assignmentFirstBlock) {
            double machineSpeed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
            // use the real position in the machine with instance->getE()[index block][index machine]
            double newCj = (cNode.getIndexBlock() > 0) ? cNode.getBlockStruc()[indexMachine][instance->getE()[cNode.getIndexBlock() - 1][indexMachine].second].second + pj / machineSpeed : pj /
                                                                                                                                                                                            machineSpeed;
            // update the completion time in the partial scheduling
            unsigned int positionInMachine = instance->getE()[cNode.getIndexBlock()][indexMachine].second;
            cNode.setCompletionTimeOfBlockStructure(indexMachine, positionInMachine, newCj);
            cNode.updateNbJobToSchedule(indexMachine);
            listCjAndAvailablePosition.insert({newCj, {indexMachine, positionInMachine}});
        }
        cNode.updateChangeBlock(instance->getE()[cNode.getIndexBlock() + 1].size());
    }

    while (cNode.getIndexBlock() < indexLastBlock) {
        for (auto &location: instance->getE()[cNode.getIndexBlock()]) {
            isWithinTimeLimit();
            if (!cNode.getBlockStruc()[location.first][location.second].first) {
                double machineSpeed = location.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
                double newCj = (location.second > 0) ? cNode.getBlockStruc()[location.first][location.second - 1].second + pj / machineSpeed : pj / machineSpeed;
                // update the completion time in the partial scheduling
                cNode.setCompletionTimeOfBlockStructure(location.first, location.second, newCj);
                cNode.updateNbJobToSchedule(location.first);
                listCjAndAvailablePosition.insert({newCj, location});
            }
        }
        cNode.updateChangeBlock(instance->getE()[cNode.getIndexBlock() + 1].size());
    }

    // compute completion of the last block and from assignmentLastBlock
    for (unsigned int indexMachine: assignmentLastBlock) {
        isWithinTimeLimit();
        auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> location) {
            return location.first == indexMachine;
        };
        double machineSpeed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
        // use the real position in the machine with instance->getE()[index block][index location] = (position in schedule, index machine), first find the last block (from the
        // left) where we have the machine
        double newCj;
        if (cNode.getIndexBlock() > 0) {
            // go on the left in order to find a block where the indexMachine appear for the first time
            unsigned indexBlockWithSameMachine = cNode.getIndexBlock() - 1;

            while (std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine) ==
                   instance->getE()[indexBlockWithSameMachine].end() && indexBlockWithSameMachine > 0)
                --indexBlockWithSameMachine;
            // if we indexBlockWithSameMachine == 0 and the machine 'indexMachine' is not in the block then the job appear for
            // the first time on this machine
            if (std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine) ==
                instance->getE()[indexBlockWithSameMachine].end() && indexBlockWithSameMachine == 0)
                newCj = pj / machineSpeed;
                //else we can get the completion time of the given machine
            else {
                auto location = std::find_if(instance->getE()[indexBlockWithSameMachine].begin(), instance->getE()[indexBlockWithSameMachine].end(), predFindIndexMachine);
                newCj = cNode.getBlockStruc()[location->first][location->second].second + pj / machineSpeed;
            }
        } else newCj = pj / machineSpeed;
        unsigned int positionInMachine = std::find_if(instance->getE()[cNode.getIndexBlock()].begin(), instance->getE()[cNode.getIndexBlock()].end(), predFindIndexMachine)->second;
        listCjAndAvailablePosition.insert({newCj, {indexMachine, positionInMachine}});
        cNode.updateNbJobToSchedule(indexMachine);
    }

    /********************************************/
    /*      Compute weighted tardy Jobs         */
    /********************************************/

    // use a copy of the list of grouped jobs, because this list is sorted according SPT-EDD
    std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[cNode.getIndexGroup()]);
    solveProblemWithFixedCompletionTime(&cNode, nullptr, listCjAndAvailablePosition, listOfIdenticalJobs);

    // remove other jobs in the group
    for (auto &remainJob: listOfIdenticalJobs) {
        cNode.removeOneJob(remainJob);
    }

    cNode.passNextGroup();

    // remove the number of assignment in the last block to the number of available machines.
    cNode.setAvailableLocations(cNode.getAvailableLocations() - assignmentLastBlock.size());
    changeBlock(cNode);
    #ifdef DEBUG_BaB
    std::string name;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(0) << static_cast<unsigned int>(pj);
    name.append("p=").append(stream.str()).append(",");
    for (auto assignment: assignmentFirstBlock) {
        name.append(std::to_string(assignment));
    }
    name.append("i€|");
    for (auto assignment: assignmentLastBlock) {
        name.append(std::to_string(assignment));
    }
    name.append(",|s|=").append(std::to_string(cNode.getSelectedJobCount()));
    name.append(",L=").append(std::to_string(indexLastBlock));
    cNode.stateDebug.emplace_back(name);
    #endif
    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    addNode(cNode);
    #if defined DEBUG_BaB && defined DEBUG_DOT
    // if cNode have not the same id as its parents, means we add it to the list
    if (cNode.id != pNode.id)
        dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
            << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << "\"];"
            << std::endl;
    #endif
}

inline void BeamSearch::addSolToListBestSolutions(Solution* solution,bool isLeafSolution) {
    auto & priorityQueue = isLeafSolution ? listBestFoundLeaf : listBestFoundSolutions;
    priorityQueue.emplace(solution->getSumWjUj(), *solution);
    // we have more than nbBestSolutionKeep solutions, remove the greater one
    if (priorityQueue.size() > nbBestSolutionKeep) priorityQueue.pop();
}

inline void BeamSearch::addSolToListBestSolutions(Solution::BlockStructure & blockStructure,double objValue){
    if (not listBestFoundSolutions.empty() && isSmaller(objValue, listBestFoundSolutions.top().first)) {
        Solution sol(instance);
        sol.fromBlockStruct(blockStructure);
        listBestFoundSolutions.emplace(objValue, sol);
        if (listBestFoundSolutions.size() > nbBestSolutionKeep) listBestFoundSolutions.pop();
    }
}

inline void BeamSearch::addNode(Node &node) {
    #if defined DEBUG_BaB
    node.id = nbNodeLoc+1;
    #if defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif
    #endif
    // if not enough job or if we are not in the last group then cut
    if (instance->getNbJobs() - node.getRemovedJobsCount() < instance->getNbToSelectJob()) {
        #if defined DEBUG_BaB
        ++nbNodeLoc;
            #if defined DEBUG_DOT
            dot << node.id << DOT_CUT << ";";
            dot << node.id << "[label=\""<<Node::debug_dot_node(node)<<",LB: NE_JOBS\"];" << std::endl;
            #endif
        #endif
        return;
    }
    ++nbNodeLoc; // incr the nb node explored
    // if we have to select and schedule jobs on the last block, we can solve it by assignment algorithm instead of branching, we obtain a leaf node
    if (node.getSelectedJobCount() == instance->getNbToSelectJob() - instance->getE().back().size()) {
        Solution solFromNode(instance);
        auto blockStruct = node.getBlockStruc();
        solFromNode.fromBlockStruct(blockStruct);
        listAvailableJobForNode.clear();
        for (auto job: instance->getListJobs()) {
            if (not node.isScheduled(job.getIndex()) && not node.isRemoved(job.getIndex())) {
                listAvailableJobForNode.push_back(job);
            }
        }
        //here the method freeAndAssignmentBlock, is optimal for the last block since we solve an assignment problem with all available jobs on the last block
        heuristicSolver.freeAndAssignmentBlock(blockStruct, node.getIndexBlock(), &listAvailableJobForNode);
        solFromNode.fromBlockStruct(blockStruct);
        if (solFromNode.feasible(instance)) {
            // if we run local search on leaf node
            addSolToListBestSolutions(&solFromNode, true);
            if (solFromNode.getSumWjUj() < globalUB) {
                #if defined DEBUG_BaB && defined DEBUG_DOT
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("LB:").append(std::to_string(solFromNode.getSumWjUj())));
                #endif
                *solution = solFromNode;
                globalUB = solFromNode.getSumWjUj();
            }
            #if defined DEBUG_BaB && defined DEBUG_DOT
            else {
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("SOL:").append(std::to_string(solFromNode.getSumWjUj())));
            }
            #endif
            return;
        } else throw BiSchException("Select n job, unfeasible where as branching scheme make feasible");
    } else // if we have selected the right nb of job
    if (node.getSelectedJobCount() == instance->getNbToSelectJob()) {
        Solution solFromNode(instance);
        solFromNode.fromBlockStruct(node.getBlockStruc());
        if (solFromNode.feasible(instance)) {
            addSolToListBestSolutions(solution,true);
            if (solFromNode.getSumWjUj() < globalUB) {
                #if defined DEBUG_BaB && defined DEBUG_DOT
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("LB:").append(std::to_string(solFromNode.getSumWjUj())));
                #endif
                *solution = solFromNode;
                globalUB = solFromNode.getSumWjUj();
            }
            #if defined DEBUG_BaB && defined DEBUG_DOT
            else {
                dot << node.id << DOT_OPT << ";";
                node.stateDebug.emplace_back(std::string("SOL:").append(std::to_string(solFromNode.getSumWjUj())));
            }
            #endif
            return;
        } else throw BiSchException("Select n job, unfeasible where as branching scheme make feasible");
    }
    auto beamNode = evaluation(std::move(node));
    // add the node
    #if defined DEBUG_BaB && defined DEBUG_DOT
    // copy the node to get same stateDebug after the call of the method
    heap.push(evaluation(node));
    #else
    if (isSmallerOrEqual(beamNode.lowerbound,beamNode.upperbound))
        heap.push(std::move(beamNode));
    #endif
}

inline void BeamSearch::branchingLocation(BeamSearch::BeamSearchNode &beamSearchNode) {

    #if defined DEBUG_BaB && defined DEBUG_DOT
    std::ofstream dot(DEBUG_DOT, std::ios::app);
    #endif

    double nodeEvaluation = beamSearchNode.evaluation;
    Node &pNode = beamSearchNode.node;

    if (verbose >= 2) {
        unsigned int sizeActive = heap.size();
        std::cout << "Active size:" << sizeActive << " UB: " << globalUB << " evaluation: " << nodeEvaluation << " nbCleaningSetCol: " << columnGeneration.getNbCleaningSetCol() << " setCol size: "
                  << columnGeneration.getXs().getSize() << " nbNode: " << nbNodeLoc << std::endl;
    }

    // get the index of the group, it's the max between the indexGroup of both kind of machine.
    unsigned int indexGroup = pNode.getIndexGroup();
    auto &listGroupedJobs = instance->getListGrpJobs();
    if (listGroupedJobs[indexGroup].size() > 1) {
        createChildrenNodeWithIdenticalJobs(pNode);
    } else {
        // use a boolean to know if the current job is on time whatever its position.
        bool jobIsEarly = true;
        // fixe a location
        for (auto itLocation = instance->getE()[pNode.getIndexBlock()].begin(); itLocation != instance->getE()[pNode.getIndexBlock()].end(); ++itLocation) {
            isWithinTimeLimit();
            // create child node
            Node cNode = Node(pNode);
            // if the location is available
            if (!cNode.getBlockStruc()[itLocation->first][itLocation->second].first) {

                // the machine speed
                char machineSpeed = itLocation->first < instance->getNbOfHighSpeedMachines() ? 0 : 1;
                double speed = (machineSpeed == 0) ? instance->getHighSpeed() : instance->getLowSpeed();

                // we add the jobs to the fixed location
                const Job &scheduledJob = listGroupedJobs[indexGroup].back();
                double C_j = (itLocation->second > 0) ? cNode.getBlockStruc()[itLocation->first][itLocation->second - 1].second + scheduledJob.getPi() / speed : scheduledJob.getPi() / speed;
                // the job is late ?
                if (jobIsEarly)
                    jobIsEarly = isSmallerOrEqual(C_j, scheduledJob.getDi());
                //add the pointer from the list of job (of instance)
                cNode.scheduleOneJob(itLocation->first, itLocation->second, &instance->getListJobs()[scheduledJob.getIndex()], C_j);
                cNode.updateNbJobToSchedule(itLocation->first);

                // update constant that count the number of available location
                cNode.setAvailableLocations(cNode.getAvailableLocations() - 1);
                changeBlock(cNode);
                // add the node to the list of active node
                #if defined DEBUG_BaB
                std::string name = std::string("J").append(std::to_string(scheduledJob.getIndex())).append(
                        ",i").append(
                        std::to_string(itLocation->first)).append(",k=").append(
                        std::to_string(itLocation->second));
                cNode.stateDebug.emplace_back(name);
                #endif
                addNode(cNode);
                #if defined DEBUG_BaB && defined DEBUG_DOT
                // if cNode have not the same id as its parents, means we add it to the list
                if (cNode.id != pNode.id)
                    dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id
                        << "[label=\"" << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str()
                        << "\"];" << std::endl;
                #endif
                // if we are in the first block, then  don't try for the same job to schedule it on identical
                // machine.
                if (pNode.getIndexBlock() == 0) {
                    // So if we are on high speed machine, then go to low speed machine
                    if (itLocation->first < instance->getNbOfHighSpeedMachines()) {
                        // it's the nb of high speed machine - the index of current machine - 1 (-1 because the loop
                        //increase of one)
                        unsigned int nbJump = instance->getNbOfHighSpeedMachines() - itLocation->first - 1;
                        std::advance(itLocation, nbJump);
                    }//else we are on low speed machine then we can break
                    else {
                        break;
                    }
                }
            }
        }
        // if we are in first block and there is no late jobs and the current job is on time, then don't create the node where we remove a job
        if (jobIsEarly && pNode.getIndexBlock() == 0 && isSmallerOrEqual(pNode.getPartialSumWjUj(), 0.0)) {
            return;
        } else {
            // don't select job
            Node cNode = Node(pNode);
            #if defined DEBUG_BaB
            unsigned int indexOfRemovedJob =
            #endif
            cNode.removeCurrentJob();
            #if defined DEBUG_BaB
            std::string name = std::to_string(indexOfRemovedJob).append("X");
            cNode.stateDebug.emplace_back(name);
            #endif
            addNode(cNode);
            #if defined DEBUG_BaB && defined DEBUG_DOT
            // if cNode have not the same id as its parents, means we add it to the list
            if (cNode.id != pNode.id)
                dot << pNode.id << " -- " << cNode.id << std::endl << cNode.id << "[label=\""
                    << "(id:" << cNode.id << "," << cNode.stateDebug.back().append(")").c_str() << ","
                    << indexOfRemovedJob << "X\"];" << std::endl;
            #endif
        }
    }
}

inline void BeamSearch::changeBlock(Node &node) {
    // if we have totally filled the block, then we prepare the next (if is not last, i.e., indexBlock < |E|-1).
    if (node.getAvailableLocations() + node.getNumJobsToScheduleOnBlock() == instance->getE()[node.getIndexBlock()].size() && node.getIndexBlock() < instance->getE().size() - 1) {
        //update other jobs with release date that is the minimum completion time of the block
        double minCj = std::numeric_limits<double>::infinity();
        // for all machine, compute the completion time of each ones and update the release date of other jobs
        for (unsigned int i = 0; i < instance->getNbMachines(); ++i) {
            //if there is no jobs, then go block before until we find a job
            unsigned int indexBlock = node.getIndexBlock();
            while (!node.getBlockStruc()[i][indexBlock].first && indexBlock > 0)
                --indexBlock;
            double C_ik = node.getBlockStruc()[i][indexBlock].second;
            if (C_ik < minCj) {
                minCj = C_ik;
            }
        }
        node.updateChangeBlock(instance->getE()[node.getIndexBlock() + 1].size());
    }
}


#endif //BILEVEL_SCHEDULING_BEAMSEARCH_IMP_H
