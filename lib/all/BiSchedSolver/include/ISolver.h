//  Copyright (C) 2024
//  Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
//
//  DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
//  This file is part of bilevel-scheduling.
//
//  bilevel-scheduling is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published
//  by the Free Software Foundation, either version 3 of the License,
//  or (at your option) any later version.
//
//  bilevel-scheduling is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty
//  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//  See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

#ifndef BILEVEL_SCHEDULING_ISOLVER_H
#define BILEVEL_SCHEDULING_ISOLVER_H

#include <iostream>
#include <fstream>
#include "Solution.h"
#include "Node.h"
#include <chrono>
#include <unordered_set>
#include <nlohmann/json.hpp>
#include <thread>
#include <memory>

// epsilon for equality between two double : a - b < EPSILON <=> a~=b
#define EPSILON 1E-6

/**
 * Interface for the Solver. Each new method or algorithm to solve the bilevel scheduling problem must implement with
 * its own class that use ISolver interface
 */
class ISolver {
protected:
    Solution *solution; // the solution of the problem
    Instance *instance; // a reference on instance
    char verbose;
    std::chrono::steady_clock::time_point start;
    std::chrono::duration<double> time_elapsed;
    std::chrono::duration<double> time_limits = std::chrono::seconds(60);
    std::shared_ptr<std::atomic<bool>> timeUp; //boolean to stop solve when time limit is reached
    std::thread timerThread; //thread only for time
    bool timerActive = false; //boolean to make the timer active

public:
    // declare the binary tree of completion times that is not assigned with the corresponding position,
    // i.e. key is the completion time and the value is (index of machine, index of block). This is sorted in increasing order
    typedef std::multimap<double, std::pair<unsigned int, unsigned int>> TreeCj;

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit ISolver() : solution(new Solution()), instance(nullptr), verbose(0), time_elapsed(0) {}

    /**
     * Constructor of a solver from a instance.
     * @param instance The instance to solve
     */
    explicit ISolver(Instance *instance) : solution(new Solution(instance)), instance(instance), verbose(0), time_elapsed(0) {}

    /**
     * Constructor of a solver from an json object that solve the instance.
     * @param instance The instance to solve
     * @param object The json object with the parameters of the solver
     */
    explicit ISolver(Instance *instance, [[maybe_unused]] nlohmann::json &object) :
            solution(new Solution(instance)), instance(instance), verbose(0){}

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    virtual ~ISolver() {
        // save the solution in a "solutions" directory at same location as the instance file
        std::string strPath = instance->getInstancePath().parent_path().string();
        strPath.append("/solutions/");
        strPath.append(instance->getInstanceName()).append(".sol");

        std::filesystem::path fileSolutionPath = std::filesystem::path(strPath);
        std::filesystem::create_directories(fileSolutionPath.lexically_normal().parent_path());
        std::ofstream outputFileStream;
        outputFileStream.open(fileSolutionPath, std::ios::out | std::ios::trunc);
        outputFileStream << *solution;
        outputFileStream.close();
        delete solution;
    }


    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Runs the solver with the timer thread.
     * This method starts the timer, executes the solver algorithm, and then stops the timer.
     * It ensures all derived solvers automatically respect the time limit.
     */
    void run(std::string nameMethod) {
        setStartTime();       // optional: keep start time for elapsed reporting
        timeUp = std::make_shared<std::atomic<bool>>(false); // the the flag to stop if time limit is reached
        startTimerThread();   // launch timer thread

        try {
            solve();          // call the virtual solve() implemented by the derived class
        } catch (BiSchTimeOutException &e) {
            if (verbose >= 1)
                std::cout << "Solver stopped due to time limit." << std::endl;
        }catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
        }catch (...) {
            std::cerr << "Exception unknown" << std::endl;
        }
        const auto endSolve{std::chrono::steady_clock::now()};
        time_elapsed = std::chrono::duration<double>{endSolve - start};
        stopTimerThread();    // safely join timer thread
        if (verbose >=2) std::cout << nameMethod.c_str() << " is over after " << time_elapsed.count() << " seconds" << std::endl;
    }

    /**
     * Pure Virtual Method that implement a algorithm to solve the current instance.
     */
    virtual void solve() = 0;

    /**
     * Assign the jobs with greatest weights to the completion times when jobs are on time.
     * @param node The current node to solve for.
     * @param blockStruct A pointer to the block structure that will be updated.
     * @param listCjAndAvailablePosition A multimap of completion times and their corresponding positions in the block.
     * @param listOfIdenticalJobs A vector of jobs that are identical.
     */
    void assigmentOnTimeJobToCompletionTime(Node *node, Solution::BlockStructure *blockStruct, std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition
                                              , std::vector<Job> &listOfIdenticalJobs) {
        isWithinTimeLimit();
        unsigned int numberElementToAdd = listCjAndAvailablePosition.size();
        // use iterator to iterate over completion time from the back (in order to get High Cj first)
        auto itCompletionTime = listCjAndAvailablePosition.rbegin();
        auto itJobs = listOfIdenticalJobs.rbegin();

        // declare the queue for max weighted early job
        auto queueEarlyJob = std::priority_queue<std::pair<std::pair<double,double>, std::vector<Job>::reverse_iterator>, std::vector<std::pair<std::pair<double,double>, std::vector<Job>::reverse_iterator >>>();

        // we use this list to delete the jobs after determine on with completion time the job will be scheduled.
        std::vector<unsigned int> indicesOfJobsToDelete;
        indicesOfJobsToDelete.reserve(listOfIdenticalJobs.size());
        std::vector<std::vector<Job>::reverse_iterator> indicesOfJobsTemporaryWaiting; // keep the list of index of job where pj/Vi > Ci
        indicesOfJobsTemporaryWaiting.reserve(listOfIdenticalJobs.size());
        // compute the set of all early jobs for each completion times, we start from the end because we remove element
        while (itCompletionTime != listCjAndAvailablePosition.rend() && itJobs != listOfIdenticalJobs.rend() && queueEarlyJob.size() < numberElementToAdd) {
            double speed = itCompletionTime->second.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
            // if pj/Vmax <= Ci then the job can be schedule (it depends on the machine)
            if (isSmallerOrEqual(itJobs->getPi()/instance->getHighSpeed(),itCompletionTime->first)) {
                // if pj/Vi <= Ci then the job can be schedule
                if (isSmallerOrEqual(itJobs->getPi()/speed,itCompletionTime->first)) {
                    // If the job can be schedule on time, i.e., dj >= Ci, we add it
                    if (isSmallerOrEqual(itCompletionTime->first, itJobs->getDi())) {
                        queueEarlyJob.emplace(std::make_pair(-itJobs->getPi(),itJobs->getWi()), itJobs);
                        ++itJobs;
                    } else if (!queueEarlyJob.empty()) {
                        // else if we have some jobs that early (means there are in the queue), so we affect the one with max weight
                        auto itAssignedJob = queueEarlyJob.top().second;
                        queueEarlyJob.pop();
                        // add the job to the block structure
                        if (node) {
                            node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first);
                        }
                        if (blockStruct) {
                            (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first};
                        }

                        // remove the completion time and the job
                        itCompletionTime = std::reverse_iterator(
                                listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
                        //add the index to the list of indices of job to delete
                        indicesOfJobsToDelete.emplace_back(itAssignedJob->getIndex());
                        // if the list of temporary indexes is not empty, add the job that can be in queue with the new Cj
                        if (itCompletionTime != listCjAndAvailablePosition.rend() && not indicesOfJobsTemporaryWaiting.empty()) {
                            // remove the jobs that can be no more schedule since the processing time is greater than the completion time
                            auto itRemove = std::remove_if(indicesOfJobsTemporaryWaiting.begin(), indicesOfJobsTemporaryWaiting.end(), [&](auto job){return isSmaller(itCompletionTime->first,job->getPi()/instance->getHighSpeed());});
                            indicesOfJobsTemporaryWaiting.erase(itRemove,indicesOfJobsTemporaryWaiting.end());
                            speed = itCompletionTime->second.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
                            // add the job that can be schedule in the queue
                            // use remove_if to swap the job that can be put in queueEarlyJob at the end of the list
                            auto itRemovedJobInWaitingList = std::remove_if(indicesOfJobsTemporaryWaiting.begin(), indicesOfJobsTemporaryWaiting.end(),
                                [&itCompletionTime,&speed](auto& job)
                                { return isSmallerOrEqual(itCompletionTime->first, job->getDi()) && isSmallerOrEqual(job->getPi() / speed, itCompletionTime->first); });
                            // add the early job that are at the end of the list to queueEarlyJob
                            for (auto itLoopWaitingJob = itRemovedJobInWaitingList; itLoopWaitingJob != indicesOfJobsTemporaryWaiting.end(); itLoopWaitingJob++){
                                queueEarlyJob.emplace(std::make_pair(-(*itLoopWaitingJob)->getPi(),(*itLoopWaitingJob)->getWi()), (*itLoopWaitingJob));
                            }
                            // then remove the early job from the list of waiting job
                            indicesOfJobsTemporaryWaiting.erase(itRemovedJobInWaitingList, indicesOfJobsTemporaryWaiting.end());
                        }
                        --numberElementToAdd;
                    } else {
                        // we pass to the next Cj
                        ++itCompletionTime;
                        // if the list of temporary indexes is not empty
                        if (itCompletionTime != listCjAndAvailablePosition.rend() && not indicesOfJobsTemporaryWaiting.empty()) {
                            // remove the jobs that can be no more schedule since the processing time is greater than the completion time
                            auto itRemove = std::remove_if(indicesOfJobsTemporaryWaiting.begin(), indicesOfJobsTemporaryWaiting.end(), [&](auto job){return isSmaller(itCompletionTime->first,job->getPi()/instance->getHighSpeed());});
                            indicesOfJobsTemporaryWaiting.erase(itRemove,indicesOfJobsTemporaryWaiting.end());
                            speed = itCompletionTime->second.first < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
                            // add the job that can be schedule in the queue
                            // use remove_if to swap the job that can be put in queueEarlyJob at the end of the list
                            auto itRemovedJobInWaitingList = std::remove_if(indicesOfJobsTemporaryWaiting.begin(), indicesOfJobsTemporaryWaiting.end(),
                                [&itCompletionTime,&speed](auto& job)
                                { return isSmallerOrEqual(itCompletionTime->first, job->getDi()) && isSmallerOrEqual(job->getPi() / speed, itCompletionTime->first); });
                            // add the early job that are at the end of the list to queueEarlyJob
                            for (auto itLoopWaitingJob = itRemovedJobInWaitingList; itLoopWaitingJob != indicesOfJobsTemporaryWaiting.end(); itLoopWaitingJob++){
                                queueEarlyJob.emplace(std::make_pair(-(*itLoopWaitingJob)->getPi(),(*itLoopWaitingJob)->getWi()), (*itLoopWaitingJob));
                            }
                            // then remove the early job from the list of waiting job
                            indicesOfJobsTemporaryWaiting.erase(itRemovedJobInWaitingList, indicesOfJobsTemporaryWaiting.end());
                        }
                    }
                }else {
                    indicesOfJobsTemporaryWaiting.emplace_back(itJobs);
                    ++itJobs;
                }
            }else {
                ++itJobs;
            }
        }

        // We clear out the queue if it not empty
        while (!queueEarlyJob.empty() && itCompletionTime != listCjAndAvailablePosition.rend()) {
            auto itAssignedJob = queueEarlyJob.top().second;
            queueEarlyJob.pop();
            if (node) {
                node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first);
            }
            if (blockStruct) {
                (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itAssignedJob->getIndex()], itCompletionTime->first};
            }
            // remove the completion time and the job
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
            //add the index to the list of indices of job to delete
            indicesOfJobsToDelete.emplace_back(itAssignedJob->getIndex());
        }
        //remove all jobs that are already scheduled
        auto pred_remove_job = [&indicesOfJobsToDelete](const Job &job) {
            return std::find(indicesOfJobsToDelete.begin(), indicesOfJobsToDelete.end(), job.getIndex()) != indicesOfJobsToDelete.end();
        };
        auto itRemovedJobs = std::remove_if(listOfIdenticalJobs.begin(), listOfIdenticalJobs.end(), pred_remove_job);
        listOfIdenticalJobs.erase(itRemovedJobs, listOfIdenticalJobs.end());
    }

    /**
     * Assign the jobs with smallest weights to the completion times when there are no on-time jobs.
     * @param node The current node to solve for.
     * @param blockStruct A pointer to the block structure that will be updated.
     * @param listCjAndAvailablePosition A multimap of completion times and their corresponding positions in the block.
     * @param listOfIdenticalJobs A vector of jobs that are identical.
     */
    void assignSmallestWeightToCompletionTime(Node *node, Solution::BlockStructure *blockStruct, std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition
                                              , std::vector<Job> &listOfIdenticalJobs) {
        isWithinTimeLimit();
        // now, we schedule not assigned jobs, i.e., the late ones, so iterate through the list of jobs
        auto itCompletionTime = listCjAndAvailablePosition.rbegin();
        // sort the remaining job according non-increasing weight, and loop over it from the end (to get the smallest weight first)
        std::sort(listOfIdenticalJobs.begin(), listOfIdenticalJobs.end(), [](const Job &lhs, const Job &rhs) { return isSmaller(rhs.getWi(),lhs.getWi()); });
        auto itJobs = listOfIdenticalJobs.rbegin();
        while (itCompletionTime != listCjAndAvailablePosition.rend() && itJobs != listOfIdenticalJobs.rend()) {
            if (node) {
                node->scheduleOneJob(itCompletionTime->second.first, itCompletionTime->second.second, &instance->getListJobs()[itJobs->getIndex()], itCompletionTime->first);
            }
            if (blockStruct) {
                (*blockStruct)[itCompletionTime->second.first][itCompletionTime->second.second] = {&instance->getListJobs()[itJobs->getIndex()], itCompletionTime->first};
            }
            itJobs = std::reverse_iterator(
                    listOfIdenticalJobs.erase(std::next(itJobs).base()));
            itCompletionTime = std::reverse_iterator(
                    listCjAndAvailablePosition.erase(std::next(itCompletionTime).base()));
        }
    }

    /**
     * Method that solves all sub-problems with identical jobs, by rearranging them in an optimal way so that the weighted sum of tardy equal-size jobs is minimal
     * @param sol The solution to be improved.
     */
    void solveProblemWithFixedCompletionTime(Solution &sol) {
        auto blockStruct = sol.toBlockStruct(instance);
        // first we determine all identical jobs in the schedule, and we minimize the weight optimally
        // use a vector of (index in list grouped jobs, tree with completion time) we will solve optimally all elements
        // of this vector
        std::vector<std::pair<unsigned int, TreeCj>> listIdenticalJobAndTheirTree;
        unsigned int indexBlock = instance->getE().size(); // start from this end of the schedule
        auto &mapJobToItsGroupIdenticalJobs = instance->getMapListJobToListGroupedJobs();
        while (indexBlock) {
            --indexBlock;
            for (auto [indexMachine, indexBlockInStruct]: instance->getE()[indexBlock]) {
                auto job = blockStruct[indexMachine][indexBlockInStruct].first;
                // check if a job is scheduled
                if (job != nullptr) {
                    // check if there exist identical jobs
                    unsigned int indexInListOfGroupedJobs = mapJobToItsGroupIdenticalJobs[job->getIndex()];
                    assert(indexInListOfGroupedJobs < instance->getListGrpJobs().size());
                    auto &groupJobs = instance->getListGrpJobs()[indexInListOfGroupedJobs];
                    if (groupJobs.size() > 1) {
                        double completionTime = blockStruct[indexMachine][indexBlockInStruct].second; // get the completion time of the job
                        auto predIndexIdenticalGrp = [indexInListOfGroupedJobs](const std::pair<unsigned int, TreeCj> &pair) { return pair.first == indexInListOfGroupedJobs; };
                        // if we have already a tree with Cj and location for this group of jobs
                        auto itGroupAndItsTree = std::find_if(listIdenticalJobAndTheirTree.begin(), listIdenticalJobAndTheirTree.end(), predIndexIdenticalGrp);
                        if (itGroupAndItsTree != listIdenticalJobAndTheirTree.end()) {
                            itGroupAndItsTree->second.insert({completionTime, {indexMachine, indexBlockInStruct}});
                        }//else we add the index and its group
                        else {
                            TreeCj treeCj;
                            treeCj.insert({completionTime, {indexMachine, indexBlockInStruct}});
                            listIdenticalJobAndTheirTree.emplace_back(indexInListOfGroupedJobs, treeCj);
                        }
                    }

                }
            }
        }
        for (unsigned int indexLoopTreeCj = 0; indexLoopTreeCj < listIdenticalJobAndTheirTree.size(); ++indexLoopTreeCj) {
            auto &[indexListOfIdenticalJobs, listCjAndAvailablePosition] = listIdenticalJobAndTheirTree[indexLoopTreeCj];
            // use a copy, because it's modifying by the method, of the list of grouped jobs, that is sorted according SPT-EDD order
            std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[indexListOfIdenticalJobs]);
            solveProblemWithFixedCompletionTime(nullptr, &blockStruct, listCjAndAvailablePosition, listOfIdenticalJobs);
        }
        sol.fromBlockStruct(blockStruct);
    }

    /**
     * Method that solves all sub-problems with identical jobs, by rearranging them in an optimal way so that the weighted sum of tardy equal-size jobs is minimal
     * @param blockStruct The solution to be improved.
     */
    void solveProblemWithFixedCompletionTime(Solution::BlockStructure &blockStruct) {
        solveProblemWithFixedCompletionTime(blockStruct, instance->getE().size());
    }

    /**
     * Method that solves all sub-problems with identical jobs,
     * by rearranging them in an optimal way so that the weighted sum of tardy equal-size jobs is minimal.
     *
     * @param blockStruct The solution to be improved (reference to a BlockStructure object).
     * @param indexBlock The index of the last block from the end, from which the exact algorithm was applied to obtain the identical job schedules.
     * @param indexBlockStop The index of the first block from the beginning where we must stop applying the exact algorithm. Defaults to 0, i.e., the first block.
     */
    void solveProblemWithFixedCompletionTime(Solution::BlockStructure &blockStruct,unsigned int indexBlock, unsigned int indexBlockStop = 0) {
        bool solveOnAllTheSchedule = indexBlock == instance->getE().size(); // define a boolean to know if we solve the problem on all the schedule
        // first we determine all identical jobs in the schedule, and we minimize the weight optimally
        // use a vector of (index in list grouped jobs, tree with completion time) we will solve optimally all elements
        // of this vector
        std::vector<std::pair<unsigned int, TreeCj>> listIdenticalJobAndTheirTree;
        ; // start from this end of the schedule
        auto &mapJobToItsGroupIdenticalJobs = instance->getMapListJobToListGroupedJobs();
        auto loopIndexBlock = indexBlock;
        while (loopIndexBlock-- > indexBlockStop) {
            for (auto [indexMachine, indexBlockInStruct]: instance->getE()[loopIndexBlock]) {
                auto job = blockStruct[indexMachine][indexBlockInStruct].first;
                // check if a job is scheduled
                if (job != nullptr) {
                    // check if there exist identical jobs
                    unsigned int indexInListOfGroupedJobs = mapJobToItsGroupIdenticalJobs[job->getIndex()];
                    assert(indexInListOfGroupedJobs < instance->getListGrpJobs().size());
                    auto &groupJobs = instance->getListGrpJobs()[indexInListOfGroupedJobs];
                    if (groupJobs.size() > 1) {
                        double completionTime = blockStruct[indexMachine][indexBlockInStruct].second; // get the completion time of the job
                        auto predIndexIdenticalGrp = [indexInListOfGroupedJobs](const std::pair<unsigned int, TreeCj> &pair) { return pair.first == indexInListOfGroupedJobs; };
                        // if we have already a tree with Cj and location for this group of jobs
                        auto itGroupAndItsTree = std::find_if(listIdenticalJobAndTheirTree.begin(), listIdenticalJobAndTheirTree.end(), predIndexIdenticalGrp);
                        if (itGroupAndItsTree != listIdenticalJobAndTheirTree.end()) {
                            itGroupAndItsTree->second.insert({completionTime, {indexMachine, indexBlockInStruct}});
                        }//else we add the index and its group
                        else {
                            TreeCj treeCj;
                            treeCj.insert({completionTime, {indexMachine, indexBlockInStruct}});
                            listIdenticalJobAndTheirTree.emplace_back(indexInListOfGroupedJobs, treeCj);
                        }
                    }

                }
            }
        }
        for (unsigned int loopOverListCjAndPos = 0; loopOverListCjAndPos < listIdenticalJobAndTheirTree.size(); ++loopOverListCjAndPos) {
            auto [indexListOfIdenticalJobs, listCjAndAvailablePosition] = listIdenticalJobAndTheirTree[loopOverListCjAndPos];
            // use a copy of the list of grouped jobs, because it's modifying by the method
            std::vector<Job> listOfIdenticalJobs(instance->getListGrpJobs()[indexListOfIdenticalJobs]);
            solveProblemWithFixedCompletionTime(nullptr, &blockStruct, listCjAndAvailablePosition, listOfIdenticalJobs);
            // if we don't solve the identical job problem on all the schedule
            if (not solveOnAllTheSchedule) {
                listCjAndAvailablePosition = listIdenticalJobAndTheirTree[loopOverListCjAndPos].second;
				// a solution respect the forme A B where we have solved the problem with known completion times with identical jobs. So we
                // can have scheduled job in A and their still schedule in B. So, we have to correct schedule B.
                std::unordered_set<unsigned int> setIndexIdenticalJobsNoScheduled; //use a set of index to know which don't appear in schedule B
                for (auto &job: instance->getListGrpJobs()[indexListOfIdenticalJobs]) setIndexIdenticalJobsNoScheduled.insert(job.getIndex());
                std::unordered_set<std::pair<unsigned int,unsigned int>> setPositionInBlockStructure; // set of all position where a job is missing in B
                // remove from the set of index, jobs that have been scheduled in A
                auto itCjWithPosition = listCjAndAvailablePosition.begin();
                while (itCjWithPosition != listCjAndAvailablePosition.end()) {
                    // remove jobs' index from the set of index that come from schedule part A
                    unsigned int indexScheduledJob = blockStruct[itCjWithPosition->second.first][itCjWithPosition->second.second].first->getIndex();
                    setIndexIdenticalJobsNoScheduled.erase(indexScheduledJob);
                    ++itCjWithPosition;
                }
                // loop over the end of the schedule (i.e. the part B), and if there is a scheduled job that belongs to the group of identical jobs, then remove from the set of index
                loopIndexBlock = instance->getE().size();
                while (loopIndexBlock > indexBlock) {
                    --loopIndexBlock;
                    for (auto [indexMachine, indexBlockInStruct]: instance->getE()[loopIndexBlock]) {
                        auto job = blockStruct[indexMachine][indexBlockInStruct].first;
                        // check if a job is scheduled
                        if (job != nullptr && instance->getMapListJobToListGroupedJobs()[job->getIndex()] == indexListOfIdenticalJobs) {
                            if (setIndexIdenticalJobsNoScheduled.contains(job->getIndex())) {
                                setIndexIdenticalJobsNoScheduled.erase(job->getIndex());
                            }else {
                                setPositionInBlockStructure.insert({indexMachine, indexBlockInStruct});
                            }
                        }
                    }
                }
                // loop over the set of positions from B to assign an unscheduled job
                auto itLoopSetIndex= setIndexIdenticalJobsNoScheduled.begin();
                auto itLoopSetPosition = setPositionInBlockStructure.begin();
                // we have in worst case more unscheduled jobs than position
                assert(setPositionInBlockStructure.size() <= setIndexIdenticalJobsNoScheduled.size());
                while (itLoopSetPosition != setPositionInBlockStructure.end()) {
                    blockStruct[itLoopSetPosition->first][itLoopSetPosition->second].first = &instance->getListJobs()[*itLoopSetIndex];
                    ++itLoopSetIndex;
                    ++itLoopSetPosition;
                }
            }
        }
    }

    /**
     * Solves the problem with fixed completion times.
     * @param node The current node to solve for.
     * @param blockStruct A pointer to the block structure that will be updated.
     * @param listCjAndAvailablePosition A multimap of completion times and their corresponding positions in the block.
     * @param listOfIdenticalJobs A vector of jobs that are identical.
     */
    void solveProblemWithFixedCompletionTime(Node *node, Solution::BlockStructure *blockStruct, std::multimap<double, std::pair<unsigned int, unsigned int>> &listCjAndAvailablePosition
                                              , std::vector<Job> &listOfIdenticalJobs) {

        // assign on time job to the completion times
        assigmentOnTimeJobToCompletionTime(node,blockStruct,listCjAndAvailablePosition,listOfIdenticalJobs);
        // assign smallest weights to the completion time when there is no on time jobs
        assignSmallestWeightToCompletionTime(node,blockStruct,listCjAndAvailablePosition,listOfIdenticalJobs);


    }

    /**
     * Reset the timer for the solver.
     * This updates the start time and the time limit atomically, ensuring that
     * any ongoing timer thread sees a consistent state. Should be called
     * before starting a new solving phase or after catching a timeout exception.
     * @param newTimeLimitSec The new time limit use in the solver.
     */
    void resetTimer(double newTimeLimitSec) {
        //stop the old timer
        timerActive = false;
        // join the old thread if it is still joinable
        if (timerThread.joinable()) {
            timerThread.join();
        }
        setStartTime();
        time_limits = std::chrono::duration<double>(newTimeLimitSec);
        timeUp->store(false, std::memory_order_relaxed);
        //rerun the thread timer
        startTimerThread();

    }

    /**
     * Starts a separate thread that monitors the time limit.
     * The timer thread sleeps most of the time and sets an atomic flag once
     * the allowed time has elapsed. This allows the solver to check the time
     * without repeatedly calling chrono::steady_clock::now(), which is expensive.
     */
    void startTimerThread() {
        timeUp->store(false, std::memory_order_relaxed);
        timerActive = true;
        timerThread = std::thread([this]() {
            try {
                while (timerActive) {
                    if (std::chrono::steady_clock::now() >= start + time_limits) {
                        timeUp->store(true, std::memory_order_relaxed);
                        break;
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(5));
                }
            }catch (const std::exception& e) {
                std::cerr << "Thread error: " << e.what() << std::endl;
            }
        });
    }

    /**
     * Checks whether the time limit has been reached.
     * This method is extremely lightweight because it only reads an atomic flag.
     * If raiseException is true, a BiSchTimeOutException is thrown.
     *
     * @param raiseException If true, throws a BiSchTimeOutException when the time limit is exceeded.
     */
    void isWithinTimeLimit(bool raiseException=true) {
        assert(timeUp != nullptr);
        if (timeUp->load(std::memory_order_relaxed)) {
            if (verbose >= 2) std::cout << "Stop: time limit reached\n";
            if (raiseException) throw BiSchTimeOutException();
        }
    }

    /**
     * Stops the timer thread safely by joining it.
     * Should be called in the destructor to avoid dangling threads.
     * The timer thread will exit immediately if the time limit was reached,
     * otherwise it will naturally finish its sleep loop.
     */
    void stopTimerThread() {
        timerActive=false;
        if (timerThread.joinable()) {
            timerThread.join();
        }
        timeUp = nullptr;
    }

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] Solution *getSolution() const { return solution; }

    [[nodiscard]] const Instance *getInstance() const { return instance; }

    [[nodiscard]] char levelVerbose() const { return verbose; }

    [[nodiscard]] const std::chrono::duration<double> &getTimeElapsed() const { return time_elapsed; }

    [[nodiscard]] std::chrono::duration<double> getTimeLimits() const { return time_limits; }

    [[nodiscard]] std::shared_ptr<std::atomic<bool>> getTimeUp() const { return timeUp; }
    /********************/
    /*      SETTER      */
    /********************/

    void setSolution(Solution *solution) { ISolver::solution = solution; }

    void setInstance(Instance *instance) { ISolver::instance = instance; }

    void setVerbose(char verbose) { ISolver::verbose = verbose; }

    void setTimeLimit(unsigned int seconds) { time_limits = std::chrono::seconds(seconds); };

    void setStartTime() { start = std::chrono::steady_clock::now(); }

    void setStartTime(std::chrono::steady_clock::time_point newStart) { start = newStart; }

    void setTimeLimitInMilliSecond(unsigned int newMilli_seconds) { time_limits = std::chrono::milliseconds(newMilli_seconds); };

    void setTimeElapsed(std::chrono::duration<double> newTimeElapsed) { time_elapsed = newTimeElapsed; }

    // set timeUp using a copy
    void setTimeUp(std::shared_ptr<std::atomic<bool>> newTimeUp) {timeUp = std::move(newTimeUp);}
};

template<typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &vector) {

    os << "[";
    if (vector.empty()) os << "]";
    else {
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            const T &element = vector[indexLoopVector];
            os << element << ",";
        }
        os << vector.back() << "]";
    }
    return os;
}


#endif //BILEVEL_SCHEDULING_ISOLVER_H
