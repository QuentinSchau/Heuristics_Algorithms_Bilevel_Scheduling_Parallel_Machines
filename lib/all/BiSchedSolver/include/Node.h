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
// Created by schau on 3/26/24.
//

#ifndef BILEVEL_SCHEDULING_NODE_H
#define BILEVEL_SCHEDULING_NODE_H

#include "Instance.h"
#include "Solution.h"
#include <boost/dynamic_bitset.hpp>
#include "Math.h"

// uncomment this line to debug the method
//#define DEBUG_BaB

#ifdef DEBUG_BaB
// uncomment this line to output the search tree
// k is the position in machine schedule. k is in range [0, max Nb element in a machine schedule]
// L is the index of the block
// 'i€' means we assign on the machine index the last block of identical jobs
// i is the index of the machine
// #define DEBUG_DOT "g.dot"

#define DOT_OPT "[style=filled,fillcolor=blue]"
#define DOT_CUT "[style=filled,fillcolor=red]"
#define DOT_RECO "[style=filled,fillcolor=greenyellow]"
#define DOT_BEAM "[style=filled,fillcolor=orange]"
#endif


class Node {
private:

    /****************************/
    /*      Location scheme     */
    /****************************/

    Instance *instance;
    // partial schedule [i][k] where i is the index of machine and k the index of the block. At position [i,k] we have a job and its completion time
    Solution::BlockStructure blockStruct;
    double partial_sum_wj_Uj; // the partial value of the objective function
    // the list of boolean to know if a group of identical jobs is available
    boost::dynamic_bitset<> listGroupedJobRemoved;
    // The number of jobs that must be fixed in a block. This is equals to the number of location in a block
    // expect perhaps the first one that can have fewer jobs, there are exactly |B1| + n - |E|
    unsigned int indexBlock = 0; // the index of the block where we scheduled jobs
    unsigned int indexGroup = 0; // the index of the group of identical jobs
    unsigned int numJobsToScheduleOnBlock = 0; // the number of jobs that must be scheduled on the block
    unsigned int numJobsToScheduleOnHSMachine = 0; // the number of jobs that must be scheduled on high speed machines
    unsigned int numJobsToScheduleOnLSMachine = 0; // the number of jobs that must be scheduled on low speed machines
    unsigned int selectedJobCount = 0; // number of jobs that have been selected
    unsigned int availableLocations = 0; // number of available locations in the block indexBlock
    unsigned int removedJobsCount = 0; // number of jobs that have been removed
    boost::dynamic_bitset<> encodingRemoveDecision; // encoding vector that represents the removed jobs
    std::vector<boost::dynamic_bitset<>> encodingSelectedJobOnMachines; // encoding vector that represents the selected jobs on each machine

    // constants une in the Generation Column
    unsigned int nbCallHeuristicFailedHighSpeed = 0; // the number of calls of heuristic in columns generation with a failed, i.e. no negative reduced cost was found for high speed machine
    unsigned int nbCallHeuristicFailedLowSpeed = 0; // the number of calls of heuristic in columns generation with a failed, i.e. no negative reduced cost was found for high speed machine

public :

    #ifdef DEBUG_BaB
    std::vector<std::string> stateDebug;
    unsigned long id=0;
    #endif


    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit Node() = default;

    explicit Node(Instance *instance) {
        #if defined DEBUG_BaB && defined DEBUG_DOT
        stateDebug.emplace_back();
        #endif
        this->instance = instance;
        partial_sum_wj_Uj = 0.0;
        encodingSelectedJobOnMachines.resize(instance->getNbMachines());
        for (auto& encoding : encodingSelectedJobOnMachines) encoding = boost::dynamic_bitset<>(instance->getNbJobs());
        encodingRemoveDecision = boost::dynamic_bitset<>(instance->getNbJobs());
        selectedJobCount = 0;
        indexBlock = 0;
        removedJobsCount = 0;
        availableLocations = instance->getE()[0].size();
        numJobsToScheduleOnBlock = instance->getNbJobsToScheduleOnFirstBlock();
        // initialize completion time and block struct
        for (unsigned int i = 0; i < instance->getNbMachines(); ++i) {
            blockStruct.emplace_back(instance->getMaxNbJobsOnHS());
        }
        // compute the rest of number of jobs to schedule on each kind of machines
        for (unsigned int indexLoopBlock = 0; indexLoopBlock < instance->getE().size(); indexLoopBlock++) {
            auto &block = instance->getE()[indexLoopBlock];
            // if we are in the first block
            if (indexLoopBlock == 0) {
                numJobsToScheduleOnHSMachine += numJobsToScheduleOnBlock;
                numJobsToScheduleOnLSMachine += numJobsToScheduleOnBlock;
            } else {
                // if we can have both machines
                if (block.front().first < instance->getNbOfHighSpeedMachines() &&
                    block.back().first >= instance->getNbOfHighSpeedMachines()) {

                    numJobsToScheduleOnHSMachine += instance->getNbOfHighSpeedMachines();
                    numJobsToScheduleOnLSMachine += instance->getNbOfLowSpeedMachines();
                } else if (block.front().first >= instance->getNbOfHighSpeedMachines()) {
                    // else if we have only low speed machine for the first block
                    numJobsToScheduleOnHSMachine += 0;
                    numJobsToScheduleOnLSMachine += instance->getNbOfLowSpeedMachines();
                } else {
                    //we have only high speed machines for the first block
                    numJobsToScheduleOnHSMachine += instance->getNbOfHighSpeedMachines();
                    numJobsToScheduleOnLSMachine += 0;
                }
            }
        }
        indexGroup = 0;
        listGroupedJobRemoved = boost::dynamic_bitset<>(instance->getNbJobs());
    };

    /*******************/
    /*      SETTER     */
    /*******************/

    /**
     * Method that sets the new completion time of a pair (Job *, completion Time) in the block structure.
     * @param indexOfMachine The index of the machine.
     * @param indexOfBlock The index of the block.
     * @param newCompletionTime The completion time associated to a job.
     */
    void setCompletionTimeOfBlockStructure(unsigned int indexOfMachine, unsigned int indexOfBlock, double newCompletionTime) {
        blockStruct[indexOfMachine][indexOfBlock].second = newCompletionTime;
    }

    void setBlockStructure(Solution::BlockStructure && newBlock_structure) {
        /**
         * Let R0 the encoding of removed jobs in the node before updating block structure
         * Let S0 the encoding of selected jobs in the node before updating block structure
         * Let R1 the encoding of removed jobs in the node after updating block structure
         * Let S1 the encoding of selected jobs in the node after updating block structure
         * A = S1 AND R0 represents jobs that have been removed and now they are scheduled.
         * B = S0 XOR S1 represents jobs that change state, i.e., that have been schedule then removed or removed then scheduled.
         * C = B XOR A represents jobs that have been scheduled and now they are removed.
         * D = R XOR A represents jobs that are still removed.
         * R1 is given by: R1 = D OR C, which can be simply in : R1 = ((R0 XOR S1) OR S0) AND (R0 XOR S1 XOR S0)
         */
        assert(newBlock_structure.size() == blockStruct.size());
        partial_sum_wj_Uj = 0;
        boost::dynamic_bitset<> newEncodingSelection(instance->getNbJobs(),false);
        boost::dynamic_bitset<> encodingCurrentSelection(instance->getNbJobs(),false);
        boost::dynamic_bitset<> newEncodingRemove(instance->getNbJobs(),false);
        unsigned int nbScheduledJobs = 0;
        for (unsigned int indexMachine = 0; indexMachine < newBlock_structure.size(); ++indexMachine) {
            auto &machine = newBlock_structure[indexMachine];
            boost::dynamic_bitset<> newEncodingMachine(instance->getNbJobs(), false);
            for (auto pairJobCj: machine) {
                if (pairJobCj.first != nullptr) {
                    nbScheduledJobs++;
                    newEncodingMachine.set(pairJobCj.first->getIndex());
                    if (isSmaller(pairJobCj.first->getDi(), pairJobCj.second))
                        partial_sum_wj_Uj += pairJobCj.first->getWi();
                }
            }
            newEncodingSelection |= newEncodingMachine;
            encodingCurrentSelection |= encodingSelectedJobOnMachines[indexMachine];
            //change the encoding of the machine
            encodingSelectedJobOnMachines[indexMachine] = newEncodingMachine;
        }
        // update the encoding of removed jobs
        newEncodingRemove = ((encodingRemoveDecision ^ newEncodingSelection) | encodingCurrentSelection) & (encodingRemoveDecision ^ encodingCurrentSelection ^ newEncodingSelection);
        assert(nbScheduledJobs == selectedJobCount);
        assert(encodingCurrentSelection.count() == nbScheduledJobs);
        assert(newEncodingRemove.count() == removedJobsCount);
        // change the encoding of removed jobs
        encodingRemoveDecision = newEncodingRemove;
        // change the block structure
        blockStruct = std::move(newBlock_structure);

    }

    /**
     * Swap for the machine at the index 'indexMachine' the job at position 'positionInMachine' with the job given by the index 'indexRemovedJobs'. This method
     * updates only encodings for the case where we swap a scheduled job with a removed job.
     * @param indexMachine The index of the machine where we make the swap.
     * @param positionInMachine The position in the machine where the job must be swapped.
     * @param indexRemovedJobs The index of the removed jobs.
     */
    void swapRemovedAndScheduledJob(unsigned int indexMachine,unsigned int positionInMachine, unsigned int indexRemovedJobs){
        auto &jobWithCj = blockStruct[indexMachine][positionInMachine];
        encodingSelectedJobOnMachines[indexMachine].flip(jobWithCj.first->getIndex()).flip(indexRemovedJobs);
        encodingRemoveDecision.flip(jobWithCj.first->getIndex()).flip(indexRemovedJobs);
        jobWithCj.first = &instance->getListJobs()[indexRemovedJobs];
        // update the completion time of the machine
        double completionTime = 0.0;
        double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
        for (auto &jobInSchedule: blockStruct[indexMachine]) {
            if (jobInSchedule.first != nullptr) {
                completionTime += jobInSchedule.first->getPi()/speed;
                jobInSchedule.second = completionTime;
            }
        }
    }

    void setAvailableLocations(unsigned int availableLocations) {
        Node::availableLocations = availableLocations;
    }

    void setNbCallHeuristicFailed(unsigned int nbCallHeuristicFailed, char machineSpeed) {
        if (machineSpeed == 0)
            Node::nbCallHeuristicFailedHighSpeed = nbCallHeuristicFailed;
        else
            Node::nbCallHeuristicFailedLowSpeed = nbCallHeuristicFailed;
    }


    /*******************/
    /*      GETTER     */
    /*******************/


    [[nodiscard]] unsigned int getIndexBlock() const {
        return indexBlock;
    }

    [[nodiscard]] const Solution::BlockStructure &getBlockStruc() const {
        return blockStruct;
    }

    [[nodiscard]] const Solution::BlockStructure * getPointerOnBlockStruc() const {
        return &blockStruct;
    }

    [[nodiscard]] unsigned int getIndexGroup() const {
        return indexGroup;
    }

    /**
     * Getter to know if a group of identical jobs was removed.
     * @param indexGroupOfIdenticalJobs The index of the group of identical jobs.
     * @return True if the group of identical jobs was removed, false otherwise.
     */
    [[nodiscard]] bool isGroupIdenticalJobRemoved(unsigned int indexGroupOfIdenticalJobs) const {
        return listGroupedJobRemoved[indexGroupOfIdenticalJobs];
    }

    /**
     * Getter that returns the processing time of the current group of identical jobs
     */
    [[nodiscard]] double getCurrentPj() const { return instance->getListGrpJobs()[indexGroup].back().getPi(); }

    [[nodiscard]] unsigned int getNumJobsToScheduleOnBlock() const {
        return numJobsToScheduleOnBlock;
    }

    [[nodiscard]] unsigned int getSelectedJobCount() const {
        return selectedJobCount;
    }

    [[nodiscard]] unsigned int getAvailableLocations() const {
        return availableLocations;
    }

    [[nodiscard]] unsigned int getRemovedJobsCount() const {
        return removedJobsCount;
    }

    [[nodiscard]] const boost::dynamic_bitset<> &getEncodingRemoveDecision() const {
        return encodingRemoveDecision;
    }

    [[nodiscard]] const std::vector<boost::dynamic_bitset<>> &getEncodingSelectedJobOnMachine() const {
        return encodingSelectedJobOnMachines;
    }

    [[nodiscard]] double getPartialSumWjUj() const {
        return partial_sum_wj_Uj;
    }

    [[nodiscard]] unsigned int getNbCallHeuristicFailed(char machineSpeed) const {
        return machineSpeed == 0 ? nbCallHeuristicFailedHighSpeed : nbCallHeuristicFailedLowSpeed;
    }

    /**
     * Method checks if a jobs is removed.
     * @param indexJob The index of the jobs.
     * @return True if the job is removed false otherwise.
     */
    [[nodiscard]] bool isRemoved(unsigned int indexJob) const { return encodingRemoveDecision.test(indexJob); }

    /**
     * Method checks if a jobs is scheduled.
     * @param indexJob The index of the jobs.
     * @return True if the job is scheduled false otherwise.
     */
    [[nodiscard]] bool isScheduled(unsigned int indexJob) const {
        return std::any_of(encodingSelectedJobOnMachines.begin(), encodingSelectedJobOnMachines.end(), [indexJob](auto &encodingMachine) { return encodingMachine.test(indexJob); });
    }

    /**
     * Method that checks if the job is already scheduled on other machines.
     * @param indexJob The index of the job.
     * @param machineSpeed The machine speed of the job where it must be scheduled. Thus 0 for high speed 1 for low speed. If we want to check if the job is scheduled on low speed machine, that mean
     * we assume that it is scheduled on high speed therefor machineSpeed must be equal to 0.
     * @return True if the job is scheduled on the other machines.
     */
    [[nodiscard]] bool isScheduledOnOtherMachines(unsigned int indexJob, char machineSpeed) const {
        // check if the job is not on other machine
        bool isScheduledOnOtherMachine = false;
        unsigned int indexMachine = machineSpeed == 0 ? instance->getNbOfHighSpeedMachines() : 0;
        unsigned int indexLastMachine = machineSpeed == 0 ? instance->getNbMachines() : instance->getNbOfHighSpeedMachines();
        for (; indexMachine < indexLastMachine; indexMachine++) {
            isScheduledOnOtherMachine = encodingSelectedJobOnMachines[indexMachine].test(indexJob);
            if (isScheduledOnOtherMachine) break;
        }
        return isScheduledOnOtherMachine;
    }

    /**
     * Getter to have the number of jobs that must be schedule on high speed and low speed machines.
     * @return A pair of number of job that we have to schedule on high and low speed machines.
     */
    [[nodiscard]] std::pair<unsigned int, unsigned int> getNumberJobsToScheduleOnMachines() const {
        return {numJobsToScheduleOnHSMachine, numJobsToScheduleOnLSMachine};
    }

    /**
 * Method that computes the Cmax for each machine in a schedule respecting the block structure property, and adds
 * these values to the given list of completion times. The list is then sorted in decreasing order for high-speed
 * machines (first elements) and again for low-speed machines (second elements), resulting in a vector with two
 * separate parts sorted independently.
 *
 * @param listCompletionTime The list where we add the completion times
 * @param blockStructure The block structure representing a partial solution
 */
    void getListOfCompletion(std::vector<double> &listCompletionTime, Solution::BlockStructure &blockStructure) {
        listCompletionTime.assign(blockStructure.size(), 0.0);
        for (unsigned int indexOfMachine = 0; indexOfMachine < blockStructure.size(); ++indexOfMachine) {
            auto const &machine = blockStructure[indexOfMachine];
            // find the last insert job
            for (auto itFromTheEndOfMachine = machine.rbegin(); itFromTheEndOfMachine != machine.rend();) {
                if (itFromTheEndOfMachine->first != nullptr) {
                    listCompletionTime[indexOfMachine] = itFromTheEndOfMachine->second;
                    break;
                } else
                    ++itFromTheEndOfMachine;
            }
        }
        // sort high speed machine first
        std::sort(listCompletionTime.begin(), listCompletionTime.begin() + instance->getNbOfHighSpeedMachines(), std::greater<>());
        // then sort low speed machine
        std::sort(listCompletionTime.begin(), listCompletionTime.begin() + instance->getNbOfHighSpeedMachines(), std::greater<>());

    };

    /**
     * Method that computes the Cmax for each machine in a schedule respecting the block structure property, and adds
     * these values to the given list of completion times. The list is then sorted in decreasing order for high-speed
     * machines (first elements) and again for low-speed machines (second elements), resulting in a vector with two
     * separate parts sorted independently.
     *
     * @param listCompletionTime The list where we add the completion times
     */
    void getListOfCompletion(std::vector<double> &listCompletionTime) {
        getListOfCompletion(listCompletionTime, blockStruct);
    };

    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Method that sets all bits corresponding to a given job index in an encoding bitset to true.
     * @param indexOfJob The index of the job whose position is set to true
     * @param encoding The encoding bitset where the bits are set to true
     */
    static void setEncodingForJob(unsigned int indexOfJob, boost::dynamic_bitset<> &encoding) {
        encoding.set(indexOfJob);
    };

    /**
     * Method that removes a group of identical jobs. This method returns the count of elements in the group being removed.
     * @return The number of elements in the group being removed
     */
    unsigned int removeGroupOfIdenticalJobs() {
        auto &listOfGroupedJobs = instance->getListGrpJobs();
        unsigned int numberOfIdenticalJobs = listOfGroupedJobs[indexGroup].size();
        assert(indexGroup < listOfGroupedJobs.size());
        for (auto &job: listOfGroupedJobs[indexGroup]) {
            setEncodingForJob(job.getIndex(), encodingRemoveDecision);
            ++removedJobsCount;
        }
        listGroupedJobRemoved.set(indexGroup, true);
        ++indexGroup;
        return numberOfIdenticalJobs;
    }

    /**
     * Method that removes the only job in the current group of identical jobs. This method throws an exception if the
     * group contains more than one job, as this is not a valid scenario. It increments the number of removed jobs.
     * @return The index of the removed job
     */
    unsigned int removeCurrentJob() {
        assert(indexGroup < instance->getListGrpJobs().size());
        if (instance->getListGrpJobs()[indexGroup].size() > 1)
            throw BiSchException("The current group of identical jobs must contain one element");
        auto &scheduledJob = instance->getListGrpJobs()[indexGroup].back();
        unsigned int indexOfRemovedJob = scheduledJob.getIndex();
        setEncodingForJob(indexOfRemovedJob, encodingRemoveDecision);;
        ++removedJobsCount;
        ++indexGroup;
        return indexOfRemovedJob;
    }

    /**
     * Method that adds a removed job to the list of removed jobs. It increments the number of removed jobs.
     * @param removedJob The job being removed
     */
    void removeOneJob(Job &removedJob) {
        setEncodingForJob(removedJob.getIndex(), encodingRemoveDecision);
        ++removedJobsCount;
    }

    /**
     * Method that adds the pair (*Job, double), i.e., a pointer to a job and the associated completion time.
     * This pair is added to the structure of partial scheduling in a node.
     * @param indexOfMachine The index of the machine.
     * @param indexOfBlock The index of the block.
     * @param pointerToScheduledJob A pointer to a job.
     * @param newCompletionTime The completion time associated to a job.
     */
    void scheduleOneJob(unsigned int indexOfMachine, unsigned int indexOfBlock, const Job *pointerToScheduledJob, double newCompletionTime) {
        // add the job to the block structure
        blockStruct[indexOfMachine][indexOfBlock] = {pointerToScheduledJob, newCompletionTime};
        if (newCompletionTime > pointerToScheduledJob->getDi()) partial_sum_wj_Uj += pointerToScheduledJob->getWi();
        ++selectedJobCount; //incr nb scheduled job
        setEncodingForJob(pointerToScheduledJob->getIndex(), encodingSelectedJobOnMachines[indexOfMachine]);
        // if the current group of jobs is a single job, then pass to the next group otherwise, we have to manage the case of several identical jobs
        if (instance->getListGrpJobs()[indexGroup].size() == 1) ++indexGroup;
    };

    /**
     * Update the constant used in the computing of the number of job that must be scheduled.
     * @param indexOfMachine the index of machine where we schedule a job
     */
    void updateNbJobToSchedule(unsigned int indexOfMachine) {
        if (indexBlock == 0) {
            numJobsToScheduleOnLSMachine--;
            numJobsToScheduleOnHSMachine--;
        } else {
            if (indexOfMachine < instance->getNbOfHighSpeedMachines()) {
                numJobsToScheduleOnHSMachine--;
            } else {
                numJobsToScheduleOnLSMachine--;
            }
        }
    };

    /**
     * Method that updates some constants of this class when the block changes.
     * @param numberOfNewAvailableLocations The number of new available location in the next block.
     */
    void updateChangeBlock(unsigned int numberOfNewAvailableLocations) {
        ++indexBlock;
        availableLocations = numberOfNewAvailableLocations;
        numJobsToScheduleOnBlock = numberOfNewAvailableLocations;
    }

    void updatePartialSumWj(){
        partial_sum_wj_Uj = 0.0;
        for (unsigned int indexMachine = 0; indexMachine < blockStruct.size(); ++indexMachine) {
            double completionTime = 0.0;
            double speed = indexMachine < instance->getNbOfHighSpeedMachines() ? instance->getHighSpeed() : instance->getLowSpeed();
            for (auto &jobInSchedule: blockStruct[indexMachine]) {
                if (jobInSchedule.first != nullptr) {
                    completionTime += jobInSchedule.first->getPi() / speed;
                    partial_sum_wj_Uj += isSmaller(jobInSchedule.first->getDi(), completionTime) ? jobInSchedule.first->getWi() : 0.0;
                }
            }
        }

    }
    void passNextGroup() { indexGroup++; }

    void incrementNbCallHeuristic(char machineSpeed) { machineSpeed == 0 ? ++nbCallHeuristicFailedHighSpeed : ++nbCallHeuristicFailedLowSpeed; }

    #if defined DEBUG_BaB && defined DEBUG_DOT
    static std::string debug_dot_node(Node & node){
        std::stringstream output;
        std::string name = node.stateDebug.back();
        output << "(id:" << node.id<<"," <<name.append(")").c_str();
        return output.str();
    }

    #endif

};
#endif //BILEVEL_SCHEDULING_NODE_H
