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
// Created by schau on 1/22/25.
//

#ifndef BILEVEL_SCHEDULING_HEURISTIC_H
#define BILEVEL_SCHEDULING_HEURISTIC_H


#include "ISolver.h"
#include "HungarianAlgorithm.h"

// uncomment this line to debug the Class
// #define DEBUG_HEURISTIC

class Heuristic : public ISolver {
private:
    // The matrix of cost for each machine and each block. costAtEachBlock[i][k] is the cost of the machine 'i' when the job ended at the block 'k'. The last element of 'costAtEachBlock[i]' is the
    // cost of the machine schedule.
    std::vector<std::vector<double>> costAtEachBlock;
    // The list of jobs that can be schedule on a block that we release
    std::vector<unsigned int> listAvailableIndexJob;
    // The matrix of cost use to solve the assigment problem
    std::vector<std::vector<double>> costMatrix;
    double smallEpsilon = 0.1;

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit Heuristic();

    explicit Heuristic(Instance *instance);

    explicit Heuristic(Instance *instance, nlohmann::json &object);

    void initializeStructure();

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~Heuristic() override;

    void solve() override;

    /**********************/
    /*      Heuristic     */
    /**********************/

    /**
     * Method that clear costAtEachBlock.
     */
    void clearCostBlocks();

    /**
    * Method that clear the cost matrix used to solve assigment problem.
    */
    void clearCostMatrix();

    /**
     * Method that computes the cost of the schedule at each block and for each machines.
     * @param blockStruct The block structure from where we release jobs.
     * @param indexBlock The index of the block to release.
     */
    void computeAllWeightsAtEachBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock);


    /**
    * Method that computes the jobs that can be scheduled in the block we want to release. It will look in the list of jobs given in parameter (if it provides) and add jobs that are already in the
    * block structure.
    * @param blockStruct The block structure from where we release jobs.
    * @param indexBlock The index of the block to release.
    * @param listOfJobsAvailable The list of jobs from which we try to fill the block that we have freed.
    * @param listIndexUnchangedMachines The list of indices of machines from which we don't free the job on it (machines to avoid).
    */
    void computeListJobsCanBeScheduleInBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines = nullptr);

    /**
     * Method that compute small epsilon to apply to the weight w_{i,j}. This epsilon aims to get the SPT order for machines and for a job we prefer first indexed machines both when we solve assignment
     * problem.
     * @param i The index of the machine.
     * @param j The index of the jobs.
     * @return The epsilon to apply to the weight w_{i,j}.
     */
    double changeWeightWithEpsilon(unsigned int i, unsigned int j);

    /**
     * Method that computes the matrix of cost used to solve the assignment problem.
     * The matrix C[i,j] is the cost of scheduling job Jj in machine i of the block that we want to release. We don't compute the weight for unchanged machines.
     * @param blockStruct The block structure from where we release jobs.
     * @param indexBlock The index of the block to release.
     * @param listIndexUnchangedMachines The list of indices of machines from which we don't compute costs
     * @param numJobsToScheduleOnBlock The number of job to schedule in the block, by default is -1 meaning all jobs in the block
     */
    void computeCostMatrixWithUnchangedMachines(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<unsigned int> *listIndexUnchangedMachines,int numJobsToScheduleOnBlock = -1);

    void computeCostMatrixWithUnchangedMachinesV1(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<unsigned int> *listIndexUnchangedMachines,int numJobsToScheduleOnBlock = -1);

    /**
     * Method that releases jobs in a block and tries to assign new jobs using an assignment problem.
     * @param blockStruct The block structure from where we release jobs.
     * @param indexBlock The index of the block to release.
     * @param listJobsAvailable The list of available jobs from which we try to fill the released block. The method change this list at this end, by removing the jobs that have been assigned.
     * @param listIndexUnchangedMachines The list of indices of machines from which we don't free the job on it (machines to avoid).
     * @param numJobsToScheduleOnBlock The number of job to schedule in the block, by default is -1 meaning all jobs in the block
     */
    void freeAndAssignmentBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines = nullptr, int numJobsToScheduleOnBlock = -1);

    void freeAndAssignmentBlockV1(Solution::BlockStructure &blockStruct, unsigned int indexBlock, std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines = nullptr, int numJobsToScheduleOnBlock = -1);

    /**
     * Method that upgrade the solution given in parameter with heuristic.
     * @param sol The solution to improve.
     * @param listJobsAvailable The list of jobs from which we try to fill the block that we have freed.
     * @param leftToRight boolean to know if we must go from the first block 0 to the last block (left to right order)
     * @deprecated Use the version with block structure.
     */
    void upgradeSolutionWithHeuristic(Solution &sol, const std::vector<Job> &listJobsAvailable,bool leftToRight=true);

    /**
     * Method that upgrade the solution given in parameter with heuristic.
     * @param blockStructure The block structure representing a solution to improve.
     * @param listIndexJobsAvailable The list of index of jobs from which we try to fill the block that we have freed.
     * @param leftToRight boolean to know if we must go from the first block 0 to the last block (left to right order)
     */
    void upgradeSolutionWithHeuristic(Solution::BlockStructure &blockStructure, std::vector<unsigned int> * listIndexJobsAvailable,bool leftToRight=true);

    /**
     * Method that constructs a solution from the partial solution given by the block structure.
     * @param blockStruct The block structure representing a partial solution. The parameter is passed as copy to work on it.
     * @param indexBlock The index of the last black that is not fulfilled.
     * @param nbScheduledJob The number of jobs that have been scheduled.
     * @param listJobsAvailable The list of jobs available to fill the partial solution. Must be sorted according SPT
     * @return A complete block structure constructed by heuristic from the partial solution
     */
    Solution::BlockStructure constructSolutionWithHeuristic(Solution::BlockStructure blockStruct,unsigned int indexBlock,unsigned int nbScheduledJob, std::vector<Job> &listJobsAvailable);


    /********************/
    /*      GETTER      */
    /********************/

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);
};

inline void Heuristic::clearCostBlocks() {
    for (auto &costOfBlock: costAtEachBlock) {
        costOfBlock.assign(instance->getMaxNbJobsOnHS(), 0.0);
    }
}

inline void Heuristic::clearCostMatrix() {
    for (auto &costOfMachine: costMatrix) {
        costOfMachine.assign(instance->getNbJobs()+instance->getNbMachines(), instance->getSumWj());
    }
}

inline void Heuristic::computeAllWeightsAtEachBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock) {
    clearCostBlocks();
    auto &E = instance->getE();
    for (unsigned int indexLoopMachine = 0; indexLoopMachine < E[indexBlock].size(); indexLoopMachine++) {
        auto &[indexMachine,_] = E[indexBlock][indexLoopMachine];
        auto &machine = blockStruct[indexMachine];
        double costMachine = 0.0;
        for (unsigned int indexLoopInBlock = 0; indexLoopInBlock < machine.size(); indexLoopInBlock++) {
            isWithinTimeLimit();
            auto &[job, CompletionTime] = machine[indexLoopInBlock];
            if (job != nullptr) {
                costMachine += std::isless(job->getDi(), CompletionTime) ? job->getWi() : 0.0;
            }
            costAtEachBlock[indexLoopMachine][indexLoopInBlock] = costMachine;
        }
    }
}


inline void Heuristic::computeListJobsCanBeScheduleInBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines) {
    auto &E = instance->getE();
    // get the max pj from the block before and the min from the block after
    double maxPj = -1.0, minPj = std::numeric_limits<double>::infinity();
    Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance,blockStruct, indexBlock, maxPj, minPj);

    listAvailableIndexJob.clear();
    if (listOfJobsAvailable) {
        // add jobs that can be scheduled on the block
        for (auto &job: *listOfJobsAvailable) {
            isWithinTimeLimit();
            // if the job is not already present and respect block structure
            if (isSmallerOrEqual(maxPj,job.getPi()) && isSmallerOrEqual(job.getPi(),minPj) && std::find(listAvailableIndexJob.begin(), listAvailableIndexJob.end(), job.getIndex()) == listAvailableIndexJob.end()) {
                listAvailableIndexJob.emplace_back(job.getIndex());
            }
        }
    }
    // add the job from the block
    for (auto &location: E[indexBlock]) {
        isWithinTimeLimit();
        if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(),location.first)!= listIndexUnchangedMachines->end()) continue;
        auto [job, cj] = blockStruct[location.first][location.second];
        if (job != nullptr) {
            // check if we have not added the job
            auto itFindJob = std::find_if(listAvailableIndexJob.begin(), listAvailableIndexJob.end(), [job](const unsigned int indexJob) { return indexJob == job->getIndex(); });
            if (itFindJob == listAvailableIndexJob.end())
                listAvailableIndexJob.push_back(job->getIndex());
        }
    }

}

inline double Heuristic::changeWeightWithEpsilon(unsigned int i, unsigned int j) {
    //since index start from 0 and the formula is true for starting to 1 we increment both
    i++;
    j++;
    double epsilon;
    if (i == j && i == 1) {
        epsilon = smallEpsilon;
    } else if (i <= j) {
        epsilon = (static_cast<double>(j*(j-1))/2.0+static_cast<double>((j-i+1.0)*(j-1)))*smallEpsilon;
    } else {
        auto N = instance->getNbJobs();
        auto m = instance->getNbMachines();
        epsilon = (static_cast<double>(N*(N+1))/2.0+static_cast<double>(N*N + (m - i) * j + (j - 1) * m))*smallEpsilon;
    }
    assert(isSmaller(epsilon,1.0));
    return epsilon;
}

// inline void Heuristic::computeCostMatrixWithUnchangedMachines(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<unsigned int> *listIndexUnchangedMachines,int numJobsToScheduleOnBlock) {
//     clearCostMatrix();
//     auto &E = instance->getE();
//     //compute the matrix of cost
//     for (auto &[indexMachine, indexBlockInStruct]: E[indexBlock]) {
//         if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(), indexMachine) != listIndexUnchangedMachines->end()) continue;
//         double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
//         std::vector<double> &costOfMachine = costMatrix[indexMachine]; // the cost of the machine starting from indexBlock
//
//         // loop over each jobs
//         for (unsigned int indexJob: listAvailableIndexJob) {
//             // get the completion time of the location before
//             double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct - 1].second;
//             // adjust the completion time with the current job
//             auto &job = instance->getListJobs()[indexJob];
//             completionTime += job.getPi() / speed;
//
//             // initiate the cost of schedule with the sum of costs in the previous blocks.
//             double costSchedule = (indexBlock > 0) ? costAtEachBlock[indexMachine][indexBlock - 1] : 0.0;
//
//             // add the cost of the job if is late
//             double costJob = isSmallerOrEqual(job.getDi(),completionTime) ? job.getWi() : 0.0;
//             costSchedule += costJob;
//             // compute the cost of the machine schedule because we change the completion time of other next blocks
//             // So loop over the next block
//             unsigned int indexNextBlock = indexBlock + 1;
//             auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> locationInBlock) { return locationInBlock.first == indexMachine; };
//             while (indexNextBlock < E.size()) {
//                 // if a job is schedule on the machine in the next block
//                 if (std::find_if(E[indexNextBlock].begin(), E[indexNextBlock].end(), predFindIndexMachine) != E[indexNextBlock].end()) {
//                     auto &[_, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
//                     auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
//                     if (nextJob != nullptr) {
//                         completionTime += nextJob->getPi() / speed;
//                         costJob = (nextJob->getDi() < completionTime) ? nextJob->getWi() : 0.0;
//                         costSchedule += costJob;
//                     }
//                 }
//                 ++indexNextBlock;
//             }
//             costOfMachine[indexJob] = costSchedule;
//             // if we are in the first block or if we have to not fulfill a block, we shift the cost by +1 because the cost 0 is reserved for dummy machines
//             if (indexBlock == 0 || numJobsToScheduleOnBlock >= 0) costOfMachine[indexJob] += 1;
//         }
//         // if we are in the first block, we have to select k machines among m, so we create m-k dummy jobs with weight 0.0, such as we will have k assignments with weights non-negative
//         if (indexBlock == 0 || numJobsToScheduleOnBlock >= 0) {
//             // the number of jobs that must be scheduled on the block
//             if (indexBlock == 0 && numJobsToScheduleOnBlock == -1) numJobsToScheduleOnBlock = (int) instance->getNbJobsToScheduleOnFirstBlock();
//             int nbMachineInBlock = static_cast<int>(E[indexBlock].size()); // the number of machines in the block
//             assert(nbMachineInBlock >= numJobsToScheduleOnBlock);
//             for (int indexDummyJob = 0; indexDummyJob < nbMachineInBlock - numJobsToScheduleOnBlock; indexDummyJob++) {
//                 assert(instance->getNbJobs() + indexDummyJob < costOfMachine.size());
//                 costOfMachine[instance->getNbJobs() + indexDummyJob] = 0.0;
//             }
//         }
//         // apply small difference on each cost in order to take the smallest jobs when equal weight appear
//         for (unsigned int indexLoopWeight = 0; indexLoopWeight < costOfMachine.size(); indexLoopWeight++) {
//             costOfMachine[indexLoopWeight] += static_cast<double>(indexLoopWeight)/static_cast<double>(costOfMachine.size());
//         }
//     }
// }

inline void Heuristic::computeCostMatrixWithUnchangedMachines(Solution::BlockStructure &blockStruct, unsigned int indexBlock, const std::vector<unsigned int> *listIndexUnchangedMachines,int numJobsToScheduleOnBlock) {
    clearCostMatrix();
    auto &E = instance->getE();
    //compute the matrix of cost
    for (unsigned int indexLoopMachine =0 ; indexLoopMachine < E[indexBlock].size(); indexLoopMachine++) {
        auto &[indexMachine, indexBlockInStruct] = E[indexBlock][indexLoopMachine];
        if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(), indexMachine) != listIndexUnchangedMachines->end()) continue;
        double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
        std::vector<double> &costOfMachine = costMatrix[indexLoopMachine]; // the cost of the machine starting from indexBlock

        // loop over each jobs
        for (unsigned int indexLoopJob = 0; indexLoopJob < listAvailableIndexJob.size(); indexLoopJob++) {
            isWithinTimeLimit();
            // get the completion time of the location before
            double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct - 1].second;
            // adjust the completion time with the current job
            auto &job = instance->getListJobs()[listAvailableIndexJob[indexLoopJob]];
            completionTime += job.getPi() / speed;

            // initiate the cost of schedule with the sum of costs in the previous blocks.
            double costSchedule = (indexBlock > 0) ? costAtEachBlock[indexLoopMachine][indexBlock - 1] : 0.0;

            // add the cost of the job if is late
            double costJob = isSmallerOrEqual(job.getDi(),completionTime) ? job.getWi() : 0.0;
            costSchedule += costJob;
            // compute the cost of the machine schedule because we change the completion time of other next blocks
            // So loop over the next block
            unsigned int indexNextBlock = indexBlock + 1;
            auto predFindIndexMachine = [&indexMachine](std::pair<unsigned int, unsigned int> locationInBlock) { return locationInBlock.first == indexMachine; };
            while (indexNextBlock < E.size()) {
                isWithinTimeLimit();
                // if a job is schedule on the machine in the next block
                if (std::find_if(E[indexNextBlock].begin(), E[indexNextBlock].end(), predFindIndexMachine) != E[indexNextBlock].end()) {
                    auto &[_, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
                    auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
                    if (nextJob != nullptr) {
                        completionTime += nextJob->getPi() / speed;
                        costJob = (nextJob->getDi() < completionTime) ? nextJob->getWi() : 0.0;
                        costSchedule += costJob;
                    }
                }
                ++indexNextBlock;
            }
            costOfMachine[indexLoopJob] = costSchedule;
            // if we are in the first block or if we have to not fulfill a block, we shift the cost by +1 because the cost 0 is reserved for dummy machines
            if (indexBlock == 0 || numJobsToScheduleOnBlock >= 0) costOfMachine[indexLoopJob] += 1;
        }
        // if we are in the first block, we have to select k machines among m, so we create m-k dummy jobs with weight 0.0, such as we will have k assignments with weights non-negative
        if (indexBlock == 0 || numJobsToScheduleOnBlock >= 0) {
            // the number of jobs that must be scheduled on the block
            if (indexBlock == 0 && numJobsToScheduleOnBlock == -1) numJobsToScheduleOnBlock = (int) instance->getNbJobsToScheduleOnFirstBlock();
            int nbMachineInBlock = static_cast<int>(E[indexBlock].size()); // the number of machines in the block
            assert(nbMachineInBlock >= numJobsToScheduleOnBlock);
            for (int indexDummyJob = 0; indexDummyJob < nbMachineInBlock - numJobsToScheduleOnBlock; indexDummyJob++) {
                assert(listAvailableIndexJob.size() + indexDummyJob < costOfMachine.size());
                costOfMachine[listAvailableIndexJob.size() + indexDummyJob] = 0.0;
            }
        }
        // apply small difference on each cost in order to take the smallest jobs when equal weight appear
        for (unsigned int indexLoopWeight = 0; indexLoopWeight < listAvailableIndexJob.size(); indexLoopWeight++) {
            costOfMachine[indexLoopWeight] += static_cast<double>(listAvailableIndexJob[indexLoopWeight])/static_cast<double>(costOfMachine.size());
        }
    }
}

// inline void Heuristic::freeAndAssignmentBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines, int numJobsToScheduleOnBlock) {
//     auto copyBlock= blockStruct;
//     auto copyListAvailJob = (*listOfJobsAvailable);
//     computeAllWeightsAtEachBlock(blockStruct);
//     auto &E = instance->getE();
//
//     //compute the jobs that can be schedule in the released block
//     computeListJobsCanBeScheduleInBlock(blockStruct, indexBlock, listOfJobsAvailable,listIndexUnchangedMachines);
//
//     //compute the cost matrix
//     computeCostMatrixWithUnchangedMachines(blockStruct, indexBlock, listIndexUnchangedMachines,numJobsToScheduleOnBlock);
//
//     std::vector<int> assignment;
//     auto H = HungarianAlgorithm();
//     H.Solve(costMatrix, assignment);
//
//     // free the block
//     for (auto &[indexMachine, indexBlockInStruct]: E[indexBlock]) {
//         // if the machine is unchanged then continue
//         if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(),indexMachine)!= listIndexUnchangedMachines->end()) continue;
//         // add the job to the list of available jobs. The job will be removed if it is again assigned
//         if (blockStruct[indexMachine][indexBlockInStruct].first != nullptr && listOfJobsAvailable != nullptr)
//             listOfJobsAvailable->push_back(*blockStruct[indexMachine][indexBlockInStruct].first);
//         blockStruct[indexMachine][indexBlockInStruct] = {nullptr, 0.0};
//     }
//     // apply the new assigment
//     for (unsigned int indexMachine = 0; indexMachine < E[indexBlock].size(); indexMachine++) {
//         if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(),indexMachine)!= listIndexUnchangedMachines->end()) continue;
//         auto [_, indexBlockInStruct] = E[indexBlock][indexMachine];
//         unsigned int indexAssignJob = assignment[indexMachine];
//         // if the index of assign jobs is greater or equal to N, that means it's a dummy jobs, so we can pass to the next machines because we have to select fewer machines
//         if (indexAssignJob >= instance->getNbJobs()) {
//             assert(indexBlock == 0 || numJobsToScheduleOnBlock >= 0); // check we are in first block
//             continue;
//         }
//         blockStruct[indexMachine][indexBlockInStruct].first = &instance->getListJobs()[indexAssignJob];
//         // update the completion time of the structure and for the next block
//         double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct - 1].second;
//         double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
//         unsigned int indexNextBlock = indexBlock;
//         while (indexNextBlock < E.size()) {
//             // if a job is schedule on the machine in the next block
//             if (indexMachine < E[indexNextBlock].size()) {
//                 auto &[_d, indexNextBlockInStruct] = E[indexNextBlock][indexMachine];
//                 auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
//                 if (nextJob != nullptr) {
//                     completionTime += nextJob->getPi() / speed;
//                     blockStruct[indexMachine][indexNextBlockInStruct].second = completionTime;
//                 }
//             }
//             ++indexNextBlock;
//         }
//     }
//     // remove assigned jobs from the list of available jobs
//     if (listOfJobsAvailable != nullptr) {
//         for (unsigned int indexMachine = 0; indexMachine < E[indexBlock].size(); indexMachine++) {
//             if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(), indexMachine) != listIndexUnchangedMachines->end()) continue;
//             unsigned int indexAssignJob = assignment[indexMachine];
//             auto itRemovedAssignedJob = std::remove_if(listOfJobsAvailable->begin(), listOfJobsAvailable->end(), [indexAssignJob](Job &job) {return job.getIndex() == indexAssignJob;});
//             listOfJobsAvailable->erase(itRemovedAssignedJob, listOfJobsAvailable->end());
//         }
//     }
//     freeAndAssignmentBlockV1(copyBlock,indexBlock,&copyListAvailJob,listIndexUnchangedMachines,numJobsToScheduleOnBlock);
//     if (copyBlock != blockStruct) {
//         Solution::printB(copyBlock);
//         std::cout << std::endl;
//         Solution::printB(blockStruct);
//     }
//
// }

inline void Heuristic::freeAndAssignmentBlock(Solution::BlockStructure &blockStruct, unsigned int indexBlock, std::vector<Job> *listOfJobsAvailable,const std::vector<unsigned int> *listIndexUnchangedMachines, int numJobsToScheduleOnBlock) {

    computeAllWeightsAtEachBlock(blockStruct,indexBlock);
    auto &E = instance->getE();

    //compute the jobs that can be schedule in the released block
    computeListJobsCanBeScheduleInBlock(blockStruct, indexBlock, listOfJobsAvailable,listIndexUnchangedMachines);

    //compute the cost matrix
    computeCostMatrixWithUnchangedMachines(blockStruct, indexBlock, listIndexUnchangedMachines,numJobsToScheduleOnBlock);
    int nbCols = static_cast<int>(listAvailableIndexJob.size());
    // if we are in the first block, we have to select k machines among m, so we create m-k dummy jobs with weight 0.0, such as we will have k assignments with weights non-negative
    if (indexBlock == 0 || numJobsToScheduleOnBlock >= 0) {
        if (indexBlock == 0 && numJobsToScheduleOnBlock == -1) numJobsToScheduleOnBlock = (int) instance->getNbJobsToScheduleOnFirstBlock();
        int nbMachineInBlock = static_cast<int>(E[indexBlock].size()); // the number of machines in the block
        int nbDummyJobs = nbMachineInBlock - numJobsToScheduleOnBlock;
        nbCols += nbDummyJobs;
    }
    std::vector<int> assignment;
    auto H =  HungarianAlgorithm();
    H.Solve(costMatrix, assignment,nbCols,E[indexBlock].size());

    // free the block
    for (auto &[indexMachine, indexBlockInStruct]: E[indexBlock]) {
        isWithinTimeLimit();
        // if the machine is unchanged then continue
        if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(),indexMachine)!= listIndexUnchangedMachines->end()) continue;
        // add the job to the list of available jobs. The job will be removed if it is again assigned
        if (blockStruct[indexMachine][indexBlockInStruct].first != nullptr && listOfJobsAvailable != nullptr)
            listOfJobsAvailable->push_back(*blockStruct[indexMachine][indexBlockInStruct].first);
        blockStruct[indexMachine][indexBlockInStruct] = {nullptr, 0.0};
    }
    // apply the new assigment
    for (unsigned int indexLoopMachine = 0; indexLoopMachine < E[indexBlock].size(); indexLoopMachine++) {
        isWithinTimeLimit();
        auto [indexMachine, indexBlockInStruct] = E[indexBlock][indexLoopMachine];
        if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(),indexMachine)!= listIndexUnchangedMachines->end()) continue;
        // if the index of assign jobs is greater or equal to N, that means it's a dummy jobs, so we can pass to the next machines because we have to select fewer machines
        if (assignment[indexLoopMachine] >= (int) listAvailableIndexJob.size()) {
            assert(indexBlock == 0 || numJobsToScheduleOnBlock >= 0); // check we are in first block
            continue;
        }
        unsigned int indexAssignJob = listAvailableIndexJob[assignment[indexLoopMachine]];
        blockStruct[indexMachine][indexBlockInStruct].first = &instance->getListJobs()[indexAssignJob];
        // update the completion time of the structure and for the next block
        double completionTime = (indexBlockInStruct == 0) ? 0.0 : blockStruct[indexMachine][indexBlockInStruct - 1].second;
        double speed = (indexMachine < instance->getNbOfHighSpeedMachines()) ? instance->getHighSpeed() : instance->getLowSpeed();
        unsigned int indexNextBlock = indexBlock;
        while (indexNextBlock < E.size()) {
            isWithinTimeLimit();
            // if a job is schedule on the machine in the next block
            if (indexMachine <= E[indexNextBlock].back().first && indexMachine >= E[indexNextBlock].front().first) {
                auto &[_d, indexNextBlockInStruct] = E[indexNextBlock][indexLoopMachine];
                auto nextJob = blockStruct[indexMachine][indexNextBlockInStruct].first;
                if (nextJob != nullptr) {
                    completionTime += nextJob->getPi() / speed;
                    blockStruct[indexMachine][indexNextBlockInStruct].second = completionTime;
                }
            }
            ++indexNextBlock;
        }
    }
    // remove assigned jobs from the list of available jobs
    if (listOfJobsAvailable != nullptr) {
        for (unsigned int indexLoopMachine = 0; indexLoopMachine < E[indexBlock].size(); indexLoopMachine++) {
            isWithinTimeLimit();
            auto [indexMachine, indexBlockInStruct] = E[indexBlock][indexLoopMachine];
            if (listIndexUnchangedMachines != nullptr && std::find(listIndexUnchangedMachines->begin(), listIndexUnchangedMachines->end(), indexMachine) != listIndexUnchangedMachines->end()) continue;
            assert(indexMachine < assignment.size());
            // if we index in assigment is feasible, i.e. if we are not in first block because we can select less jobs and we have introduced some dummy job
            if (assignment[indexLoopMachine] >= (int) listAvailableIndexJob.size()) {
                assert(indexBlock == 0 || numJobsToScheduleOnBlock >= 0); // check we are in first block
                continue;
            }
            unsigned int indexAssignJob = listAvailableIndexJob[assignment[indexLoopMachine]];
            auto itRemovedAssignedJob = std::remove_if(listOfJobsAvailable->begin(), listOfJobsAvailable->end(), [indexAssignJob](Job &job) {return job.getIndex() == indexAssignJob;});
            listOfJobsAvailable->erase(itRemovedAssignedJob, listOfJobsAvailable->end());
        }
    }

}

inline void Heuristic::upgradeSolutionWithHeuristic(Solution &sol, const std::vector<Job> &listJobsAvailable,bool leftToRight) {
    Solution::BlockStructure blockStructure = sol.toBlockStruct(instance);
    std::vector<unsigned int> listIndexJobs(listJobsAvailable.size());
    for (unsigned int indexLoopListJobs =0 ; indexLoopListJobs < listJobsAvailable.size(); indexLoopListJobs ++) {
        listIndexJobs[indexLoopListJobs] = listJobsAvailable[indexLoopListJobs].getIndex();
    }
    upgradeSolutionWithHeuristic(blockStructure, &listIndexJobs,leftToRight);
}

inline void Heuristic::upgradeSolutionWithHeuristic(Solution::BlockStructure &blockStructure, std::vector<unsigned int> * listIndexJobsAvailable, bool leftToRight) {
    std::vector<Job> listOfJobs;
    if (listIndexJobsAvailable != nullptr) {
        listOfJobs.resize(listIndexJobsAvailable->size());
        std::transform(listIndexJobsAvailable->begin(),
                       listIndexJobsAvailable->end(),
                       listOfJobs.begin(),
                       [&](unsigned int indexJob) {
                           return instance->getListJobs()[indexJob];
                       });
        Solution::removeExistingJobsFromSolution(listOfJobs, blockStructure);
    }
    #ifdef DEBUG_HEURISTIC
    double objValBefore = Solution::evaluate(blockStructure,instance); // keep the value of solution before apply the heuristic
    #endif
    // try to free and re-assign jobs from the RIGHT to LEFT
    unsigned int indexBlock = 0;
    // try now from the LEFT to RIGHT
    while (leftToRight ? indexBlock < instance->getE().size() : indexBlock-- > 0) {
        freeAndAssignmentBlock(blockStructure, indexBlock, &listOfJobs);
        // then solve optimally all identical jobs on blocks before
        solveProblemWithFixedCompletionTime(blockStructure,indexBlock+1);
        if (leftToRight) ++indexBlock;
    }
    #ifdef DEBUG_HEURISTIC
    Solution sol(instance);
    sol.fromBlockStruct(blockStructure);
    if (not sol.feasible(instance)) {
        std::cout << "Error Heuristic on:" << instance->getInstancePath() << std::endl;
        throw BiSchException("Error after assignment heuristic");
    }
    if (verbose >= 3) std::cout << "UB before upgrade: " << objValBefore << " after: " << sol.getSumWjUj() << std::endl;
    #endif
}

inline Solution::BlockStructure Heuristic::constructSolutionWithHeuristic(Solution::BlockStructure blockStruct, const unsigned int indexBlock, const unsigned int nbScheduledJob, std::vector<Job> &listJobsAvailable) {
    isWithinTimeLimit();
    // first, identify the block that is not totally filled
    auto &E = instance->getE();
    TreeCj treeCj; // use structure with tree Cj to compute using fixe completion time
    // The minimal number of job that can be schedule on previous blocks
    unsigned int minNbJobToScheduleOnPrevBlock = instance->getNbToSelectJob() - nbScheduledJob;
    // we fill block from the end (right -> left)
    unsigned int indexLoopBlock = E.size()-1;
    //defined predicat to find jobs that are in the block
    auto pred = [&indexLoopBlock,&blockStruct,&E](const Job& job) {
        for (auto &[indexMachine, positionInMachine]: E[indexLoopBlock]) {
            const Job * scheduledJobs = blockStruct[indexMachine][positionInMachine].first;
            if ( scheduledJobs != nullptr && (*scheduledJobs) == job ) return true;
        }
        return false;
    };
    while (indexLoopBlock > indexBlock) {
        isWithinTimeLimit();
        // compute the mean processing time on the list of available jobs
        double meanPj = 0;
        for (auto &job : listJobsAvailable) meanPj += job.getPi();
        meanPj /= static_cast<double>(listJobsAvailable.size());

        // compute the estimation of the completion time on the block for each machines
        treeCj.clear(); // clear the tree
        for (auto &[indexMachine,indexInMachine] : E[indexLoopBlock]) {
            // if there is no already assigned job
            if (blockStruct[indexMachine][indexInMachine].first == nullptr) {
                double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();// speed of the machine
                // loop back until the beginning (left) to find the number of unsigned position and the first scheduled jobs with its completion time
                int indexLoopBackBlock = indexLoopBlock == 0 ? 0 : indexLoopBlock -1;
                double completionTime = 0.0;
                unsigned int indexInMachineBack = indexInMachine;
                while (indexLoopBackBlock >= 0) {
                    // if the block loopBackBlock have the machine, i.e. the index of hte machine is between the first and the last index of machine in the 'indexLoopBackBlock'
                    if (E[indexLoopBackBlock].front().first <= indexMachine && indexMachine <= E[indexLoopBackBlock].back().first) {
                        indexInMachineBack = E[indexLoopBackBlock][indexMachine].second;
                        auto & jobWithCj = blockStruct[indexMachine][indexInMachineBack];
                        if (jobWithCj.first != nullptr) {
                            completionTime = jobWithCj.second;
                            break;
                        }
                    }
                    indexLoopBackBlock--;
                }
                completionTime += static_cast<double>(indexInMachine-indexInMachineBack)*meanPj;
                // if we have indexLoopBackBlock == -1, that mean the is no job on the machine so add meanPj to the completion time
                if (indexLoopBackBlock == -1) completionTime += meanPj;
                completionTime /= speed;
                treeCj.insert({completionTime, {indexMachine, indexInMachine}});
            }
        }
        // update the minNbJobToScheduleOnPrevBlock
        minNbJobToScheduleOnPrevBlock -= E[indexLoopBlock].size();
        // extract from the list of available jobs, the one that can be schedule, i.e., we don't take the first minNbJobToScheduleOnPrevBlock jobs
        assert(minNbJobToScheduleOnPrevBlock < listJobsAvailable.size());
        // create a range based copy of the list of job
        std::vector<Job> listJobThatCanBeSchedule(listJobsAvailable.cbegin() + minNbJobToScheduleOnPrevBlock, listJobsAvailable.cend());
        // sort the list of jobs that can be schedule according EDD rule
        std::sort(listJobThatCanBeSchedule.begin(),listJobThatCanBeSchedule.end(),Job::EDD);
        assigmentOnTimeJobToCompletionTime(nullptr, &blockStruct, treeCj, listJobThatCanBeSchedule);
        if (not listJobThatCanBeSchedule.empty() && not treeCj.empty()) {
            // sort the list of jobs according weight and assign the smallest weights to remaining completion times
            std::sort(listJobThatCanBeSchedule.begin(),listJobThatCanBeSchedule.end(),[](Job &lhs, Job &rhs){return isSmaller(lhs.getWi(),rhs.getWi());});
            auto itLoopRemainingJobs = listJobThatCanBeSchedule.begin();
            auto itLoopRemainingCj  = treeCj.begin();
            while (itLoopRemainingCj != treeCj.end() && itLoopRemainingJobs != listJobThatCanBeSchedule.end()) {
                auto &[indexMachine,indexInMachine] = itLoopRemainingCj->second;
                auto &jobWithCj = blockStruct[indexMachine][indexInMachine];
                jobWithCj.first = &instance->getListJobs()[itLoopRemainingJobs->getIndex()];
                ++itLoopRemainingJobs;
                ++itLoopRemainingCj;
            }
        }
        // find the minimal processing time
        double minPj = std::numeric_limits<double>::infinity();
        for (auto &[indexMachine,indexInMachine]: instance->getE()[indexLoopBlock]) {
            auto &jobWithCj = blockStruct[indexMachine][indexInMachine];
            if (jobWithCj.first != nullptr) minPj = std::min(minPj, jobWithCj.first->getPi());
        }
        // remove from the list of available jobs the ones that are already in the block
        auto itRemove = std::remove_if(listJobsAvailable.begin(), listJobsAvailable.end(), pred);
        listJobsAvailable.erase(itRemove, listJobsAvailable.end());
        // remove from the list of available jobs the jobs that have greater processing time than the min pj
        itRemove = std::remove_if(listJobsAvailable.begin(), listJobsAvailable.end(), [minPj](Job &jobToRem){return isSmaller(minPj,jobToRem.getPi());});
        listJobsAvailable.erase(itRemove, listJobsAvailable.end());
        indexLoopBlock--;
    }
    // know we have the case where indexLoopBlock == indexBlock, some we can compute exactly the optimal assigment because we know all completion time
    std::vector<unsigned int> listUnchangedMachine;
    listUnchangedMachine.reserve(instance->getNbMachines());
    for (auto &[indexMachine,indexInMachine] : E[indexBlock]) {
        // if there is already assigned job, then don't change the machines
        if (blockStruct[indexMachine][indexInMachine].first != nullptr) listUnchangedMachine.emplace_back(indexMachine);
    }
    freeAndAssignmentBlock(blockStruct,indexBlock,&listJobsAvailable,&listUnchangedMachine);

    return blockStruct;
}

#endif //BILEVEL_SCHEDULING_HEURISTIC_H
