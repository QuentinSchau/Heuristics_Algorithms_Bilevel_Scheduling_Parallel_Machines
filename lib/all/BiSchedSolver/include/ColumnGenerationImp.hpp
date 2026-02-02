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


#define BILEVEL_SCHEDULING_COLUMNGENERATION_IMP_H

#include "Heuristic.h"

/*******************/
/*      GETTER     */
/*******************/

inline double ColumnGeneration::getSumWjUj() const {
    return lowerBound;
}

/********************/
/*      METHODS     */
/********************/

inline unsigned int ColumnGeneration::nextIndex(const std::vector<std::vector<Job>> &listOfJobs, unsigned int pairedIndexOfJob) {
    auto [indexGroup, indexInGroup] = getGroupAndJobIndex(pairedIndexOfJob);
    ++indexInGroup;
    if (indexInGroup == listOfJobs[indexGroup].size()) {
        ++indexGroup;
        indexInGroup = 0;
    }
    return elegantPair(indexGroup, indexInGroup);
}

inline std::pair<unsigned int, unsigned int> ColumnGeneration::getGroupAndJobIndex(unsigned int pairedIndexOfJob) {
    std::pair<unsigned int, unsigned int> indexes = elegantUnpair(pairedIndexOfJob);
    return indexes;
}

inline void ColumnGeneration::updateValueOfLmax(Node &node, char machineSpeed) {
    auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
    std::fill(lmax.begin(), lmax.end(), 0);
    unsigned int indexMachine = (machineSpeed == 0) ? 0 : instance->getNbOfHighSpeedMachines();
    auto &listGroupedJobs = instance->getListGrpJobs();
    unsigned int maxNbJobMachine = (machineSpeed == 0) ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
    unsigned int numberJobToSelect = 1;
    auto predFindMachineInBlock = [&indexMachine](std::pair<unsigned int, unsigned int> locationInBlock) {
        return locationInBlock.first == indexMachine;
    };
    // if we have more than 0 jobs at maximum of the machine, we will compute lmax for 0 to maxNbJobMachine - 1 jobs selected
    for (; numberJobToSelect <= maxNbJobMachine; numberJobToSelect++) {
        for (unsigned int indexGroup = 0; indexGroup < listGroupedJobs.size(); indexGroup++) {
            if (node.isGroupIdenticalJobRemoved(indexGroup)) continue;
            auto &listGroupedJob = listGroupedJobs[indexGroup];
            unsigned int numberJobsInGroup = listGroupedJob.size();
            unsigned int nbJobToSelectOnTheMachine = 0; // we count the number of jobs that we have to select before get numberJobToSelect
            unsigned int nbPosWeighted = 0; // the number of positional weighted that we have visited
            // start from the end
            auto itLoopBlockFirstPass = instance->getE().rbegin();
            for (; itLoopBlockFirstPass != instance->getE().rend(); itLoopBlockFirstPass++) {
                if (std::find_if(itLoopBlockFirstPass->begin(), itLoopBlockFirstPass->end(), predFindMachineInBlock) != itLoopBlockFirstPass->end()) {
                    nbJobToSelectOnTheMachine++;
                }
                if (nbJobToSelectOnTheMachine == numberJobToSelect) {
                    nbPosWeighted++;
                    break;
                }
                nbPosWeighted += itLoopBlockFirstPass->size();
            }
            unsigned int smallestPosWeighted = numberJobsInGroup <= nbPosWeighted ? nbPosWeighted - numberJobsInGroup + 1 : 0; // the smallest positional weighted we can reach with our identical jobs
            // reset the number of explored positional weighted
            nbPosWeighted = 0;
            auto itLoopBlockSecondPass = instance->getE().rbegin();
            for (; itLoopBlockSecondPass != instance->getE().rend(); itLoopBlockSecondPass++) {
                nbPosWeighted += itLoopBlockSecondPass->size();
                if (nbPosWeighted >= smallestPosWeighted) {
                    break;
                }
            }
            // count the number of machine that appear between itLoopBlockFirstPass and itLoopBlockSecondPass
            unsigned int newLvalue = 1;
            for (auto itLoopBlock = itLoopBlockSecondPass; itLoopBlock != itLoopBlockFirstPass; itLoopBlock++) {
                if (std::find_if(itLoopBlock->begin(), itLoopBlock->end(), predFindMachineInBlock) != itLoopBlock->end())
                    newLvalue++;
            }
            assert(hashLmax(numberJobToSelect, indexGroup) < lmax.size());
            lmax[hashLmax(numberJobToSelect, indexGroup)] = newLvalue;
            maxLValue = std::max(maxLValue, newLvalue);
        }
    }
}

/********************/
/*      Columns     */
/********************/

inline double ColumnGeneration::computeScheduleCost(Node &node, std::vector<unsigned int> &machineSchedule, char machineSpeed) {
    ++nbCallComputeCost;
    double sumWeightedTardyJobs = 0.0;
    double completionTime = 0.0;
    // get the speed
    double speed = (machineSpeed == 0) ? instance->getHighSpeed() : instance->getLowSpeed();
    // compute the cost
    for (unsigned int &indexInMachineSchedule: machineSchedule) {
        // check if the job can be schedule
        const Job &job = instance->getListJobs()[indexInMachineSchedule];
        // check if the job is not on other machine
        bool isScheduledOnOtherMachine = node.isScheduledOnOtherMachines(job.getIndex(), machineSpeed);
        if (not node.isRemoved(job.getIndex()) && not isScheduledOnOtherMachine) {
            completionTime += job.getPi();
            // check if the job is late, whenever add its weighted
            if (isSmaller(job.getDi(), completionTime / speed)) sumWeightedTardyJobs += job.getWi();
        } else { return instance->getSumWj(); }

    }
    return sumWeightedTardyJobs;
}

inline double ColumnGeneration::subRoutine(unsigned int nbJobsToSelected, unsigned int t, unsigned int g, char machineSpeed) {
    // find the 'nbJobsToSelected' jobs that have the smallest contribution to the reduced cost
    double minContributionToReducedCost = std::numeric_limits<double>::infinity();
    auto &listIdenticalJobs = instance->getListGrpJobs()[g];
    unsigned int nbIdJob = listIdenticalJobs.size();
    unsigned int p = static_cast<unsigned int>(listIdenticalJobs.back().getPi());
    //if nbJobsToSelected == 1, we can achieve it in O(n)
    if (nbJobsToSelected == 1) {
        auto jobToSelect = std::min_element(listIdenticalJobs.begin(), listIdenticalJobs.end(), [&](const Job &J1, const Job &J2) {
            return isSmaller(beta(t, J1, machineSpeed), beta(t, J2, machineSpeed));
        });
        minContributionToReducedCost = beta(t, *jobToSelect, machineSpeed);
    } else if (nbJobsToSelected == 2) {
        for (unsigned int indexLoopFirstPosJob = 0; indexLoopFirstPosJob < nbIdJob; ++indexLoopFirstPosJob) {
            for (unsigned int indexLoopSecondPosJob = 0; indexLoopSecondPosJob < nbIdJob; ++indexLoopSecondPosJob) {
                if (indexLoopFirstPosJob == indexLoopSecondPosJob) continue;
                auto &J1 = listIdenticalJobs[indexLoopFirstPosJob];
                auto &J2 = listIdenticalJobs[indexLoopSecondPosJob];
                double contribution = beta(t, J1, machineSpeed) + beta(t + p, J2, machineSpeed);
                if (isSmaller(contribution, minContributionToReducedCost)) {
                    minContributionToReducedCost = contribution;
                }
            }
        }
    } else if (nbJobsToSelected == 3) {
        for (unsigned int indexLoopFirstPosJob = 0; indexLoopFirstPosJob < nbIdJob; ++indexLoopFirstPosJob) {
            for (unsigned int indexLoopSecondPosJob = 0; indexLoopSecondPosJob < nbIdJob; ++indexLoopSecondPosJob) {
                if (indexLoopFirstPosJob == indexLoopSecondPosJob) continue;
                for (unsigned int indexLoopThirdPosJob = 0; indexLoopThirdPosJob < nbIdJob; ++indexLoopThirdPosJob) {
                    if (indexLoopFirstPosJob == indexLoopThirdPosJob) continue;
                    if (indexLoopSecondPosJob == indexLoopThirdPosJob) continue;
                    auto &J1 = listIdenticalJobs[indexLoopFirstPosJob];
                    auto &J2 = listIdenticalJobs[indexLoopSecondPosJob];
                    auto &J3 = listIdenticalJobs[indexLoopThirdPosJob];
                    double contribution = beta(t, J1, machineSpeed) + beta(t + p, J2, machineSpeed) + beta(t + 2 * p, J3, machineSpeed);
                    if (isSmaller(contribution, minContributionToReducedCost)) {
                        minContributionToReducedCost = contribution;
                    }
                }
            }
        }
    } else {
        nbCallSubProcessCG++;
        // use Jippe algorithm
        auto &constForIdenticalJob = list_const_identical_job[g];
        if (nbJobsToSelected == 4) //if we select the first smallest number of job which is used by Jippe Algorithm
            constForIdenticalJob.initIteratorLoopMValue();
        assert(nbJobsToSelected > 0);
        auto &itLoopMValue = constForIdenticalJob.itLoopMValue;// we loop over M value by decreasing order, thus the number of selected jobs is increasing
        for (; itLoopMValue != constForIdenticalJob.list_M_value.rend(); itLoopMValue++) {
            auto const &M_value = (*itLoopMValue);
            constForIdenticalJob.computeWeights(nbIdJob, M_value, listIdenticalJobs);

            //set all job to be late
            std::fill(constForIdenticalJob.isOnTime.begin(), constForIdenticalJob.isOnTime.end(), true);
            constForIdenticalJob.computeMinReducedCostIdenticalJobs(M_value, t, machineSpeed, instance, listIdenticalJobs);

            std::vector<unsigned int> onTimeJobInEDD;
            constForIdenticalJob.constructSchedule(onTimeJobInEDD, M_value, nbIdJob);
            if (onTimeJobInEDD.size() == nbJobsToSelected) {
                double redCostOfCol = 0;//reset contribution
                unsigned int currentCompletionTime = t;
                // loop over all job in the new column, we compute the contribution to the reduced cost of the column
                // if a job can be schedule, we remove it
                for (unsigned int &indexJobInGroup: onTimeJobInEDD) {
                    assert(indexJobInGroup < listIdenticalJobs.size());
                    auto &job = listIdenticalJobs[indexJobInGroup];
                    indexJobInGroup = job.getIndex();
                    double contributionOfTheJob = beta(currentCompletionTime, job, machineSpeed);
                    // if the job can be scheduled
                    redCostOfCol += contributionOfTheJob;
                    currentCompletionTime += p;
                }
                minContributionToReducedCost = std::min(redCostOfCol, minContributionToReducedCost);
            } else if (onTimeJobInEDD.size() > nbJobsToSelected) {
                // if we don't have found minReduced cost throw error
                if (minContributionToReducedCost == std::numeric_limits<double>::infinity()) {
                    std::string message = "Error in sub routine problem, t=" + std::to_string(t) + " g=" + std::to_string(g) + " l=" + std::to_string(nbJobsToSelected);
                    std::cerr << std::setprecision(15) << "M_values: " << constForIdenticalJob.list_M_value << std::endl;
                    std::vector<double> values;
                    for (unsigned int indexLoopJobs = 0; indexLoopJobs < nbIdJob; ++indexLoopJobs) {
                        auto itVjFirstValue = std::find_if(constForIdenticalJob.list_Vj_Values.begin(), constForIdenticalJob.list_Vj_Values.end(), [&indexLoopJobs](auto const &pairValueIndex) {
                            return pairValueIndex.second == indexLoopJobs;
                        });
                        values.emplace_back(itVjFirstValue->first);
                    }
                    std::cerr << std::setprecision(15) << "Vj_values: " << values << std::endl;
                    values.clear();
                    for (unsigned int indexLoopJobs = 0; indexLoopJobs < nbIdJob; ++indexLoopJobs) {
                        auto itVjFirstValue = std::find_if(constForIdenticalJob.list_Qj_Values.begin(), constForIdenticalJob.list_Qj_Values.end(), [&indexLoopJobs](auto const &pairValueIndex) {
                            return pairValueIndex.second == indexLoopJobs;
                        });
                        values.emplace_back(itVjFirstValue->first);
                    }
                    std::cerr << std::setprecision(15) << "Qj_values: " << values << std::endl;
                    std::cerr << std::setprecision(15) << "Wj_values: " << constForIdenticalJob.list_weights << std::endl;
                    std::cerr << std::setprecision(15) << "on times: " << constForIdenticalJob.isOnTime << std::endl;
                    std::cerr << std::setprecision(15) << "M_value with error: " << M_value << std::endl;
                    std::cerr << std::setprecision(15) << "machine speed: " << int(machineSpeed) << std::endl;
                    std::cerr << std::setprecision(15) << "Jobs: ";
                    for (auto &job: listIdenticalJobs) {
                        std::cerr << job << " / ";
                    }
                    throw BiSchException(message);
                }
                break;
            }
        }
    }
    // if we have not return that means there is an error
    if (minContributionToReducedCost == std::numeric_limits<double>::infinity()) {
        std::string message = "Error in sub routine problem, t=" + std::to_string(t) + " g=" + std::to_string(g) + " l=" + std::to_string(nbJobsToSelected);
        throw BiSchException(message);
    }

    return minContributionToReducedCost;
}

inline void ColumnGeneration::subRoutine(unsigned int nbJobsToSelected, unsigned int t, unsigned int g, char machineSpeed, std::vector<unsigned int> &newColumn) {
    double minContributionToReducedCost = std::numeric_limits<double>::infinity();
    auto &listIdenticalJobs = instance->getListGrpJobs()[g];
    unsigned int nbIdJob = listIdenticalJobs.size();
    std::vector<unsigned int> bestCol(nbJobsToSelected);
    unsigned int p = static_cast<unsigned int>(listIdenticalJobs.back().getPi());
    // find the 'nbJobsToSelected' jobs that have the smallest contribution to the reduced cost
    //if nbJobsToSelected == 1, we can achieve it in O(n)
    if (nbJobsToSelected == 1) {
        for (unsigned int indexLoopFirstPosJob = 0; indexLoopFirstPosJob < nbIdJob; ++indexLoopFirstPosJob) {
            auto &J1 = listIdenticalJobs[indexLoopFirstPosJob];
            double contribution = beta(t, J1, machineSpeed);
            if (isSmaller(contribution, minContributionToReducedCost)) {
                minContributionToReducedCost = contribution;
                bestCol[0] = J1.getIndex();
            }
        }
    } else if (nbJobsToSelected == 2) {
        for (unsigned int indexLoopFirstPosJob = 0; indexLoopFirstPosJob < nbIdJob; ++indexLoopFirstPosJob) {
            for (unsigned int indexLoopSecondPosJob = 0; indexLoopSecondPosJob < nbIdJob; ++indexLoopSecondPosJob) {
                if (indexLoopFirstPosJob == indexLoopSecondPosJob) continue;
                auto &J1 = listIdenticalJobs[indexLoopFirstPosJob];
                auto &J2 = listIdenticalJobs[indexLoopSecondPosJob];
                double contribution = beta(t, J1, machineSpeed) + beta(t + p, J2, machineSpeed);
                if (isSmaller(contribution, minContributionToReducedCost)) {
                    minContributionToReducedCost = contribution;
                    bestCol[0] = J1.getIndex();
                    bestCol[1] = J2.getIndex();
                }
            }
        }
    } else if (nbJobsToSelected == 3) {
        for (unsigned int indexLoopFirstPosJob = 0; indexLoopFirstPosJob < nbIdJob; ++indexLoopFirstPosJob) {
            for (unsigned int indexLoopSecondPosJob = 0; indexLoopSecondPosJob < nbIdJob; ++indexLoopSecondPosJob) {
                if (indexLoopFirstPosJob == indexLoopSecondPosJob) continue;
                for (unsigned int indexLoopThirdPosJob = 0; indexLoopThirdPosJob < nbIdJob; ++indexLoopThirdPosJob) {
                    if (indexLoopFirstPosJob == indexLoopThirdPosJob) continue;
                    if (indexLoopSecondPosJob == indexLoopThirdPosJob) continue;
                    auto &J1 = listIdenticalJobs[indexLoopFirstPosJob];
                    auto &J2 = listIdenticalJobs[indexLoopSecondPosJob];
                    auto &J3 = listIdenticalJobs[indexLoopThirdPosJob];
                    double contribution = beta(t, J1, machineSpeed) + beta(t + p, J2, machineSpeed) + beta(t + 2 * p, J3, machineSpeed);
                    if (isSmaller(contribution, minContributionToReducedCost)) {
                        minContributionToReducedCost = contribution;
                        bestCol[0] = J1.getIndex();
                        bestCol[1] = J2.getIndex();
                        bestCol[2] = J3.getIndex();
                    }
                }
            }
        }
    } else {
        // use Jippe algorithm
        auto &constForIdenticalJob = list_const_identical_job[g];
        if (nbJobsToSelected == 4) //if we select the first smallest number of job which is used by Jippe Algorithm
            constForIdenticalJob.initIteratorLoopMValue();
        assert(nbJobsToSelected > 0);
        auto &itLoopMValue = constForIdenticalJob.itLoopMValue;// we loop over M value by decreasing order, thus the number of selected jobs is increasing
        for (; itLoopMValue != constForIdenticalJob.list_M_value.rend(); itLoopMValue++) {
            auto const &M_value = (*itLoopMValue);
            constForIdenticalJob.computeWeights(nbIdJob, M_value, listIdenticalJobs);

            //set all job to be late
            std::fill(constForIdenticalJob.isOnTime.begin(), constForIdenticalJob.isOnTime.end(), true);
            constForIdenticalJob.computeMinReducedCostIdenticalJobs(M_value, t, machineSpeed, instance, listIdenticalJobs);

            std::vector<unsigned int> onTimeJobInEDD;
            constForIdenticalJob.constructSchedule(onTimeJobInEDD, M_value, nbIdJob);
            if (onTimeJobInEDD.size() == nbJobsToSelected) {
                double redCostOfCol = 0;//reset contribution
                unsigned int currentCompletionTime = t;
                // loop over all job in the new column, we compute the contribution to the reduced cost of the column
                // if a job can be schedule, we remove it
                for (unsigned int &indexJobInGroup: onTimeJobInEDD) {
                    assert(indexJobInGroup < listIdenticalJobs.size());
                    auto &job = listIdenticalJobs[indexJobInGroup];
                    // update the vector onTimeJobInEDD with the index of the jobs
                    indexJobInGroup = job.getIndex();
                    double contributionOfTheJob = beta(currentCompletionTime, job, machineSpeed);
                    // if the job can be scheduled
                    redCostOfCol += contributionOfTheJob;
                    currentCompletionTime += p;
                }
                if (redCostOfCol < minContributionToReducedCost) {
                    minContributionToReducedCost = redCostOfCol;
                    for (unsigned int indexJobInCol = 0; indexJobInCol < onTimeJobInEDD.size(); indexJobInCol++) {
                        bestCol[indexJobInCol] = onTimeJobInEDD[indexJobInCol]; // use the index of the jobs
                    }
                }
            } else if (onTimeJobInEDD.size() > nbJobsToSelected) {
                break;
            }
        }
    }
    // if we have not return that means there is an error
    if (minContributionToReducedCost == std::numeric_limits<double>::infinity())
        throw BiSchException("Error in sub routine problem");
    for (unsigned int indexJob: bestCol) {
        newColumn.emplace_back(indexJob);
    }
}

inline void ColumnGeneration::computeFirstsStateOfDynamicProgramming(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int &indexStartingJobForDP) {
    indexStartingJobForDP = node.getIndexGroup();
    // create state from branching scheme, so we get the right index of machine depending on "machineSpeed"
    unsigned int indexLoopMachine = (machineSpeed == 0) ? 0 : instance->getNbOfHighSpeedMachines();
    unsigned int indexLastMachine = (machineSpeed == 0) ? instance->getNbOfHighSpeedMachines() : instance->getNbMachines();
    // if the current block in the node is not the first one, then we already know that a machine schedule contains max or min jobs.
    // Indeed if index block > 0, if there is no job on the first position in block structure then the machine schedule contains a min of jobs
    for (; indexLoopMachine < indexLastMachine; ++indexLoopMachine) {
        auto &machine = node.getBlockStruc()[indexLoopMachine];
        double redCostOfSchedule = 0.0;
        unsigned int completionTime = 0;
        unsigned int nbSelectedJobs = 0;
        char typeMachineSchedule = 2;
        // if the machine can be in the first block and there is a job schedule in it
        if (indexLoopMachine >= instance->getE()[0].front().first && indexLoopMachine <= instance->getE()[0].back().first && node.getBlockStruc()[indexLoopMachine][0].first != nullptr) {
            typeMachineSchedule = 1;
        } else if (node.getIndexBlock() > 0){
            // if there is the same number of max and min jobs in machine schedule, then we have to create both machine schedule
            unsigned nbMaxJob = machineSpeed == 0 ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
            unsigned nbMinJob = machineSpeed == 0 ? instance->getMinNbJobsOnHS() : instance->getMinNbJobsOnLS();
            if (nbMaxJob==nbMinJob) typeMachineSchedule = 2;
            else typeMachineSchedule = 0;
        }
        for (auto jobInSchedule: machine) {
            if (jobInSchedule.first != nullptr) {
                redCostOfSchedule += beta(completionTime, *jobInSchedule.first, machineSpeed);
                completionTime += static_cast<unsigned int>(jobInSchedule.first->getPi());
                nbSelectedJobs++;
            }
        }
        listStartingStates.emplace_back(nbSelectedJobs, completionTime, redCostOfSchedule, indexLoopMachine, typeMachineSchedule);
    }
}

inline void ColumnGeneration::updateRedCostOfStartingStateOfDP(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates) {
    for (auto &[completionTime, nbSelectedJobs, redCostOfState, indexOfMachine, typeMachineSchedule]: listStartingStates) {
        auto &machine = node.getBlockStruc()[indexOfMachine];
        //check if the machine is empty, i.e. the first and second position are not filled
        if (machine[0].first == nullptr && machine[1].first == nullptr) {
            //we continue because reduced cost is zero
            continue;
        }//the machine is not empty
        else {
            redCostOfState = 0.0;
            unsigned int completionTimeOfSchedule = 0;
            for (auto jobInSchedule: machine) {
                if (jobInSchedule.first != nullptr) {
                    redCostOfState += beta(completionTimeOfSchedule, *jobInSchedule.first, machineSpeed);
                    completionTimeOfSchedule += static_cast<unsigned int>(jobInSchedule.first->getPi());
                }
            }
        }
    }
}

inline void ColumnGeneration::update_M_Vj_Qj_Values(Node &node, char machineSpeed) {
    auto &listGroupsIdenticalJob = instance->getListGrpJobs();

    //update vj and vj-qj list
    for (unsigned int indexLoopGroupIdJobs = 0; indexLoopGroupIdJobs < listGroupsIdenticalJob.size(); ++indexLoopGroupIdJobs) {
        if (node.isGroupIdenticalJobRemoved(indexLoopGroupIdJobs)) continue;
        auto const &groupIdJobs = listGroupsIdenticalJob[indexLoopGroupIdJobs];
        unsigned int nbIdenticalJobs = groupIdJobs.size();
        //if we have only on job pass next group
        if (nbIdenticalJobs > 1) {
            unsigned int maxNbJobMachine = (machineSpeed == 0) ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
            auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
            unsigned int valueLMax = 1;
            // compute constant if the max of the L value for the corresponding group Gg is greater than 2
            for (unsigned int numberJobSelected = 0; numberJobSelected < maxNbJobMachine; numberJobSelected++) {
                assert(hashLmax(maxNbJobMachine - numberJobSelected, indexLoopGroupIdJobs) < lmax.size());
                valueLMax = std::max(valueLMax, lmax[hashLmax(maxNbJobMachine - numberJobSelected, indexLoopGroupIdJobs)]);
            }
            if (valueLMax >= 2 && nbIdenticalJobs >= 4)
                list_const_identical_job[indexLoopGroupIdJobs].computeConstants(groupIdJobs, dualsValues);
        }
    }
}

inline double ColumnGeneration::beta(unsigned int t, const Job &jobJ, char machineSpeed) {
    double result;
    double dualValues = dualsValues[jobJ.getIndex()];
    double speed = machineSpeed == 0 ? instance->getHighSpeed() : instance->getLowSpeed();
    double newT = (static_cast<double>(t) + jobJ.getPi()) / speed;
    if (isSmallerOrEqual(newT, jobJ.getDi())) result = -dualValues;
    else result = jobJ.getWi() - dualValues;
    return result;
}

inline void ColumnGeneration::constructColumnByBacktracking(Node &node, const std::vector<StateBacktracking> &listStartingState, unsigned int startingIndexJob, char machineSpeed) {
    unsigned int maxJobOnCol;
    // Create the 'nbMinStateDP' minimal columns.
    for (unsigned int indexLoopColumn = 0; indexLoopColumn < nbMinStateDP; ++indexLoopColumn) {
        unsigned int indexStartingState = 0; // use index to know on which starting state we are
        for (const auto &startingState: listStartingState) {
            //clear the list and fill it with the current forward
            auto &listState = listStatesDPForBacktracking[indexStartingState];
            listState.clear();
            std::vector<unsigned int> remainingColumn{};
            double redCostFromState;
            double constReducedValue;
            if (startingState.typeMachineSchedule == 0 || startingState.typeMachineSchedule == 2) {
                TYPE_COLUMN typeOfColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MIN : TYPE_COLUMN::LS_MIN;
                maxJobOnCol = (machineSpeed == 0) ? instance->getMinNbJobsOnHS() : instance->getMinNbJobsOnLS();
                if (startingState.alreadySelected < maxJobOnCol) {
                    remainingColumn.reserve(maxJobOnCol);
                    // Forward using the 'memo' to get the minimal reduced cost and create the first column with the minimal number of jobs
                    redCostFromState = createColumnByForward(maxJobOnCol - startingState.alreadySelected, startingState.t, startingIndexJob, remainingColumn, false, typeOfColumn, node, listState);
                    assert(std::adjacent_find(remainingColumn.begin(), remainingColumn.end()) == remainingColumn.end());
                    // add the cost of the starting state
                    redCostFromState += startingState.reducedCost;
                    constReducedValue = machineSpeed == 0 ? dualsValues[instance->getNbJobs() + 1] : dualsValues[instance->getNbJobs() + 3];
                    redCostFromState -= constReducedValue; // add the dual value of the good constraint (max or min nb of jobs)
                    #ifdef DEBUG_CG
                    if (redCostFromState != std::numeric_limits<double>::infinity()) {
                        TYPE_COLUMN typeOfColumnToTest = machineSpeed == 0 ? TYPE_COLUMN::HS_MIN : TYPE_COLUMN::LS_MIN;
                        auto copyColForTestRedCost = remainingColumn;
                        unsigned int indexOfInsertedJob = 0; // use index to insert job to the beginning of the column
                        //add the first part of the column given by the machine of the best state
                        for (auto jobInSchedule: node.getBlockStruc()[startingState.indexMachine]) {
                            if (jobInSchedule.first != nullptr) {
                                copyColForTestRedCost.insert(copyColForTestRedCost.begin() + indexOfInsertedJob, jobInSchedule.first->getIndex());
                                ++indexOfInsertedJob;
                            }
                        }
                        double computeReducedCost = computeScheduleReducedCost(node, copyColForTestRedCost, machineSpeed, typeOfColumnToTest);
                        if (std::abs(computeReducedCost - redCostFromState) > EPSILON) {
                            env.error() << "compute reduced cost failed : diff= " << std::abs(computeReducedCost - redCostFromState) << std::endl;
                            if (std::abs(computeReducedCost - redCostFromState) > 0.01)
                                throw BiSchException(std::string("The value of generated column does not correspond :").append(std::to_string(std::abs(computeReducedCost - redCostFromState))).c_str());
                        }
                        assert(copyColForTestRedCost.size() == maxJobOnCol);
                    }
                    #endif
                    if (redCostFromState < -EPSILON && not remainingColumn.empty()) {
                        unsigned int indexOfInsertedJob = 0; // use index to insert job to the beginning of the column
                        //add the first part of the column given by the machine of the best state
                        for (auto jobInSchedule: node.getBlockStruc()[startingState.indexMachine]) {
                            if (jobInSchedule.first != nullptr) {
                                remainingColumn.insert(remainingColumn.begin() + indexOfInsertedJob, jobInSchedule.first->getIndex());
                                ++indexOfInsertedJob;
                            }
                        }
                        newSetOfColumns.emplace_back(typeOfColumn, std::move(remainingColumn));
                    }
                }
            }
            remainingColumn.clear();
            if (startingState.typeMachineSchedule == 1 || startingState.typeMachineSchedule == 2) {
                maxJobOnCol = (machineSpeed == 0) ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
                TYPE_COLUMN typeOfColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::LS_MAX;
                if (startingState.alreadySelected < maxJobOnCol) {
                    remainingColumn.reserve(maxJobOnCol);
                    // Forward using the 'memo' to get the minimal reduced cost and create the first column with the maximal number of jobs, here we will update the list of state
                    redCostFromState = createColumnByForward(maxJobOnCol - startingState.alreadySelected, startingState.t, startingIndexJob, remainingColumn, true, typeOfColumn, node, listState);
                    assert(std::adjacent_find(remainingColumn.begin(), remainingColumn.end()) == remainingColumn.end());
                    // add the cost of the starting state
                    redCostFromState += startingState.reducedCost;
                    constReducedValue = machineSpeed == 0 ? dualsValues[instance->getNbJobs()] : dualsValues[instance->getNbJobs() + 2];
                    if (instance->isFirstBlockIsOnBothTypeMachine()) {
                        constReducedValue += dualsValues[instance->getNbJobs() + 1 + machineSpeed * 2];
                        constReducedValue += machineSpeed == 0 ? dualsValues[instance->getNbJobs() + 2] : 0.0;
                    }
                    redCostFromState -= constReducedValue; // add the dual value of the good constraint (max or min nb of jobs)
                    #ifdef DEBUG_CG
                    if (redCostFromState != std::numeric_limits<double>::infinity()) {
                        TYPE_COLUMN typeOfColumnToTest = machineSpeed == 0 ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::LS_MAX;
                        auto copyColForTestRedCost = remainingColumn;
                        unsigned int indexOfInsertedJob = 0; // use index to insert job to the beginning of the column
                        //add the first part of the column given by the machine of the best state
                        for (auto jobInSchedule: node.getBlockStruc()[startingState.indexMachine]) {
                            if (jobInSchedule.first != nullptr) {
                                copyColForTestRedCost.insert(copyColForTestRedCost.begin() + indexOfInsertedJob, jobInSchedule.first->getIndex());
                                ++indexOfInsertedJob;
                            }
                        }
                        double computeReducedCost = computeScheduleReducedCost(node, copyColForTestRedCost, machineSpeed, typeOfColumnToTest);
                        if (std::abs(computeReducedCost - redCostFromState) > EPSILON) {
                            env.error() << "compute reduced cost failed : diff= " << std::abs(computeReducedCost - redCostFromState) << std::endl;
                            if (std::abs(computeReducedCost - redCostFromState) > 0.01)
                                throw BiSchException(std::string("The value of generated column does not correspond :").append(std::to_string(std::abs(computeReducedCost - redCostFromState))).c_str());
                        }
                        assert(copyColForTestRedCost.size() == maxJobOnCol);
                    }
                    #endif
                    if (redCostFromState < -EPSILON && not remainingColumn.empty()) {
                        unsigned int indexOfInsertedJob = 0; // use index to insert job to the beginning of the column
                        //add the first part of the column given by the machine of the best state
                        for (auto jobInSchedule: node.getBlockStruc()[startingState.indexMachine]) {
                            if (jobInSchedule.first != nullptr) {
                                remainingColumn.insert(remainingColumn.begin() + indexOfInsertedJob, jobInSchedule.first->getIndex());
                                ++indexOfInsertedJob;
                            }
                        }
                        newSetOfColumns.emplace_back(typeOfColumn, std::move(remainingColumn));
                    }
                }
            }
            ++indexStartingState; // pass to the next state
        }
        // now, using the list of state compute by forward, we change the value of min reduced cost by backward
        indexStartingState = 0;
        for (; indexStartingState < listStartingState.size(); ++indexStartingState) {
            TYPE_COLUMN typeColumn = indexStartingState < instance->getNbOfHighSpeedMachines() ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::LS_MAX;
            // change the state in memory of the dynamic programming only if the list of state is not empty
            if (!listStatesDPForBacktracking[indexStartingState].empty())
                createColumnByBackward(typeColumn, listStatesDPForBacktracking[indexStartingState]);
        }
    }
}

inline double ColumnGeneration::createColumnByForward(unsigned int k, unsigned int t, unsigned int g, std::vector<unsigned int> &newColumn, bool keepTrackState, TYPE_COLUMN typeColumn, Node &node
                                                      , std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &listStates) {
    char machineSpeed = (typeColumn == TYPE_COLUMN::LS_MAX || typeColumn == TYPE_COLUMN::LS_MIN) ? 1 : 0;
    if (keepTrackState) listStates.emplace_back(k, t, g);
    double reducedCostOfColumn = 0.0;// the reduced cost of the column
    //reserve the size of the vector, there is at most N jobs. We add the job at the end of the column, we will reverse it
    // at the end
    auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
    unsigned int sum_n_high_pj = (machineSpeed == 0) ? sum_n0p_high_pj : sum_n1p_high_pj;
    auto &listOfJobs = instance->getListGrpJobs();
    //we stop if g is greater or equal than the number of group of identical jobs
    while (g < listOfJobs.size() && t <= sum_n_high_pj && k > 0) {
        if (not node.isGroupIdenticalJobRemoved(g)) {
            // get the reduced cost of not taking the job
            double redCostNoTakeJob = memo[hashMemo(k, t, g + 1)];
            // if we have one job to schedule
            if (listOfJobs[g].size() == 1) {
                unsigned int newT = t + static_cast<unsigned int>(listOfJobs[g].back().getPi());
                // get the reduced cost of taking the job
                double redCostTakeJob = std::numeric_limits<double>::infinity();
                double contribution = 0.0; // contribution to the reduced cost
                // get its reduced cost if it can be scheduled, i.e., it respects the release date and the deadline
                newT = t + static_cast<unsigned int>(listOfJobs[g].back().getPi());
                contribution = beta(t, listOfJobs[g].back(), machineSpeed);
                redCostTakeJob = memo[hashMemo(k - 1, newT, g + 1)] + contribution;// the cost of taking the job plus its contribution

                // take the job if the resulting reduced cost is smaller or equal to not taking it
                if (isSmallerOrEqual(redCostTakeJob, redCostNoTakeJob)) {
                    newColumn.push_back(listOfJobs[g].back().getIndex());
                    // update the completion time
                    t = newT;
                    // update the reduced cost
                    reducedCostOfColumn += contribution;
                    k--;
                }
                if (keepTrackState) listStates.emplace_back(k, t, g + 1);
            } else {
                // We have several identical jobs to schedule. We need to try all subset with 'l' jobs in the group of identical jobs.
                auto p = static_cast<unsigned int>(listOfJobs[g].back().getPi()); // the processing time of all jobs
                double minReducedCostFromTheGroup = std::numeric_limits<double>::infinity(); // the miniman reduced cost from all possibilities

                // loop over all value possible value of k
                assert(hashLmax(k, g) < lmax.size());
                unsigned int maxValueL = static_cast<int>(lmax[hashLmax(k, g)]);
                unsigned int bestValueL = 0;
                // we have several jobs in Gg, and try to select several (with the maximum value l)
                for (unsigned int l = 1; l <= maxValueL; l++) {
                    if (l > k) break;
                    unsigned int newT = t + l * static_cast<unsigned int>(listOfJobs[g].back().getPi());
                    double contributionSelection = subRoutine(l, t, g, machineSpeed);
                    double reducedCostSelectLJobs = memo[hashMemo(k - l, newT, g + 1)] + contributionSelection; // the value of the reduced cost of select the column
                    if (reducedCostSelectLJobs < minReducedCostFromTheGroup) {
                        minReducedCostFromTheGroup = reducedCostSelectLJobs;
                        bestValueL = l;
                    }
                }
                // we need to reset the iterator of M value since we loop over all M value
                list_const_identical_job[g].initIteratorLoopMValue();
                //add the best k jobs if the reduced cost is smaller or equal than take no jobs
                if (isSmallerOrEqual(minReducedCostFromTheGroup, redCostNoTakeJob)) {
                    subRoutine(bestValueL, t, g, machineSpeed, newColumn);
                    reducedCostOfColumn += (minReducedCostFromTheGroup - memo[hashMemo(k - bestValueL, t + bestValueL * p, g + 1)]); // add its contribution
                    // update the completion time
                    t = t + bestValueL * p;
                    k -= bestValueL;
                }
                if (keepTrackState) listStates.emplace_back(k, t, g + 1);
            }
        }
        // pass to the next group of identical jobs
        ++g;
    }
    if (k == 0) return reducedCostOfColumn;
    else return std::numeric_limits<double>::infinity();
}

inline void ColumnGeneration::createColumnByBackward(TYPE_COLUMN typeColumn, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &listStates) {
    char machineSpeed = (typeColumn == TYPE_COLUMN::LS_MAX || typeColumn == TYPE_COLUMN::LS_MIN) ? 1 : 0;
    // first we need to know if the last state we take the jobs
    auto [k_last, t_last, g_last] = listStates[listStates.size() - 1];
    //set the value of memo for this last state to infinity
    memo[hashMemo(k_last, t_last, g_last)] = std::numeric_limits<double>::infinity();
    auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
    //remove the last state
    listStates.pop_back();
    auto &listOfJobs = instance->getListGrpJobs();
    for (auto itLoopState = listStates.rbegin(); itLoopState != listStates.rend(); ++itLoopState) {
        auto [k, t, g] = *itLoopState;
        // get the reduced cost of not taking the job
        double redCostNoTakeJob = memo[hashMemo(k, t, g + 1)];
        // if we have one job to schedule
        if (listOfJobs[g].size() == 1) {
            unsigned int newT = t + static_cast<unsigned int>(listOfJobs[g].back().getPi());
            // get the reduced cost of taking the job
            double redCostTakeJob = std::numeric_limits<double>::infinity();
            double contribution; // contribution to the reduced cost
            // get its reduced cost if it can be scheduled, i.e., it respects the release date and the deadline
            contribution = beta(t, listOfJobs[g].back(), machineSpeed);
            redCostTakeJob = memo[hashMemo(k - 1, newT, g + 1)] + contribution;// the cost of taking the job plus its contribution
            memo[hashMemo(k, t, g)] = std::min(redCostTakeJob, redCostNoTakeJob);
        } else {
            // We have several identical jobs to schedule. We need to try all subset with 'l' jobs in the group of identical jobs.
            double minReducedCostFromTheGroup = std::numeric_limits<double>::infinity(); // the miniman reduced cost from all possibilities
            assert(hashLmax(k, g) < lmax.size());
            // loop over all value possible value of k
            unsigned int maxValueL = static_cast<int>(lmax[hashLmax(k, g)]);
            // we have several jobs in Gg, and try to select several (with the maximum value l)
            for (unsigned int l = 1; l < maxValueL; l++) {
                if (l > k) break;
                unsigned int newT = t + l * static_cast<unsigned int>(listOfJobs[g].back().getPi());
                double contributionSelection = subRoutine(l, t, g, machineSpeed);
                double reducedCostSelectLJobs = memo[hashMemo(k - l, newT, g + 1)] + contributionSelection; // the value of the reduced cost of select the column
                if (reducedCostSelectLJobs < minReducedCostFromTheGroup) {
                    minReducedCostFromTheGroup = reducedCostSelectLJobs;
                }
            }
            //add the best l jobs if the reduced cost is smaller than take no jobs
            memo[hashMemo(k, t, g)] = std::min(minReducedCostFromTheGroup, redCostNoTakeJob);
        }
    }
}

inline bool ColumnGeneration::generateColumns(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int &indexStartingJobForDP, std::tuple<unsigned int,unsigned int,double,unsigned int,double> &minReducedCosts) {
    #if defined DEBUG_BaB && defined DEBUG_CG
    try {
    #endif
    bool continueColumnGeneration = false;
    auto &listGrpJobs = instance->getListGrpJobs();
    if (!listGrpJobs.empty()) {
        // if we use dynamic programming
        if (generate_Column == 0 || node.getNbCallHeuristicFailed(machineSpeed) >= maxNbCallHeuristic) {
            minReducedCosts = solvePricingProblem(node, machineSpeed, listStartingStates, indexStartingJobForDP);
            newSetOfColumns.clear();
            constructColumnByBacktracking(node, listStartingStates, indexStartingJobForDP, machineSpeed);
        } else if (generate_Column == 1) {
            ++nbCallsHeu;
            newSetOfColumns.clear();
            double minReducedCostFromHeuristic = 0.0;
            for (auto &[completionTime, nbSelectedJobs, redCostOfSchedule, indexLoopMachine, typeMachineSchedule]: listStartingStates) {
                std::vector<unsigned int> machineSchedule;
                unsigned int indexOfInsertedJob = 0; // use index to insert job to the beginning of the column
                //add the first part of the column given by the machine of the best state
                for (auto jobInSchedule: node.getBlockStruc()[indexLoopMachine]) {
                    if (jobInSchedule.first != nullptr) {
                        machineSchedule.insert(machineSchedule.begin() + indexOfInsertedJob, jobInSchedule.first->getIndex());
                        ++indexOfInsertedJob;
                    }
                }

                double reducedCostMinJobs = std::numeric_limits<double>::infinity();
                double reducedCostMaxJobs = std::numeric_limits<double>::infinity();
                if (typeMachineSchedule == 0 || typeMachineSchedule == 2)
                    reducedCostMaxJobs = generateColumnByHeuristic(node, machineSchedule, indexStartingJobForDP, machineSpeed, 0);
                if (typeMachineSchedule == 1 || typeMachineSchedule == 2)
                    reducedCostMinJobs = generateColumnByHeuristic(node, machineSchedule, indexStartingJobForDP, machineSpeed, 1);
                minReducedCostFromHeuristic = std::min(std::min(reducedCostMaxJobs, reducedCostMinJobs), minReducedCostFromHeuristic);

            }
            // if the reduced cost is not negative, then use the dynamic programming
            if (minReducedCostFromHeuristic > -EPSILON) {
                node.incrementNbCallHeuristic(machineSpeed);
                minReducedCosts = solvePricingProblem(node, machineSpeed, listStartingStates, indexStartingJobForDP);
                newSetOfColumns.clear();
                constructColumnByBacktracking(node, listStartingStates, indexStartingJobForDP, machineSpeed);
            }

        } else throw BiSchException("Method to generate columns is not implemented");

        if (newSetOfColumns.empty()) return false;
        // add new columns
        for (auto &[typeColumn, newColumn]: newSetOfColumns) {

            // if we have reduced cost negative and a not empty column, then we can continue the generation of column
            if (!newColumn.empty() && !continueColumnGeneration) {
                continueColumnGeneration = true;
            }
            auto cost = computeScheduleCost(node, newColumn, machineSpeed);
            unsigned indexOnVariable = xs.getSize();
            IloNumArray coefConstraint(env, instance->getNbJobs() + 4);
            for (unsigned int indexJob: newColumn) {
                coefConstraint[indexJob] = 1.0;
            }
            switch (typeColumn) {
                case TYPE_COLUMN::HS_MAX:coefConstraint[instance->getNbJobs()] = 1.0;
                    if (instance->isFirstBlockIsOnBothTypeMachine()) {
                        coefConstraint[instance->getNbJobs() + 1] = 1.0;
                        coefConstraint[instance->getNbJobs() + 2] = 1.0;
                    }
                    break;
                case TYPE_COLUMN::HS_MIN:coefConstraint[instance->getNbJobs() + 1] = 1.0;
                    break;
                case TYPE_COLUMN::LS_MAX:coefConstraint[instance->getNbJobs() + 2] = 1.0;
                    if (instance->isFirstBlockIsOnBothTypeMachine()) coefConstraint[instance->getNbJobs() + 3] = 1.0;
                    break;
                case TYPE_COLUMN::LS_MIN:coefConstraint[instance->getNbJobs() + 3] = 1.0;
                    break;
            }

            xs.add(IloNumVar(objConstr(cost) + UComputeConstr(cost) + constraints(coefConstraint), 0.0, 1.0, ILOFLOAT));
            // we shift the index of machine to fit the range [0..nb speed machine]
            setMachineSchedules.emplace_back(std::move(newColumn), indexOnVariable, 0, cost, typeColumn,instance->getNbJobs());
        }
    }
    return continueColumnGeneration;
    #if defined DEBUG_BaB && defined DEBUG_CG
    } catch (BiSchException &e) {
        std::cerr << e.what() << " node id: " << node.id << " nb Gen:" << nbGen << " machine speed: "<<int(machineSpeed)<< std::endl;
        throw; // if you like
    }
    #endif
}


inline double
ColumnGeneration::generateColumnByHeuristic(Node &node, std::vector<unsigned int> &machineScheduleOnMachine, unsigned int indexStartingJobForDP, char machineSpeed, char maxOrMinJobsOnCol) {

    auto &listGroupedJobs = instance->getListGrpJobs();
    unsigned int nbJobsToSchedule = 0;
    TYPE_COLUMN typeColumn;
    if (maxOrMinJobsOnCol == 1) {
        nbJobsToSchedule = (machineSpeed == 0) ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
        typeColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::LS_MAX;
    } else {
        nbJobsToSchedule = (machineSpeed == 0) ? instance->getMinNbJobsOnHS() : instance->getMinNbJobsOnLS();
        typeColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MIN : TYPE_COLUMN::LS_MIN;
    }
    if (nbJobsToSchedule == 0) return 0.0;
    unsigned int nbJobsToSelect = nbJobsToSchedule - machineScheduleOnMachine.size();
    if (nbJobsToSelect == 0) return std::numeric_limits<double>::infinity();
    double reducedCost = 0.0;
    switch (typeColumn) {
        case TYPE_COLUMN::HS_MAX:reducedCost = -dualsValues[instance->getNbJobs()];
            if (instance->isFirstBlockIsOnBothTypeMachine()) {
                reducedCost -= dualsValues[instance->getNbJobs() + 1];
                reducedCost -= dualsValues[instance->getNbJobs() + 2];
            }
            break;
        case TYPE_COLUMN::HS_MIN:reducedCost = -dualsValues[instance->getNbJobs() + 1];
            break;
        case TYPE_COLUMN::LS_MAX:reducedCost = -dualsValues[instance->getNbJobs() + 2];
            if (instance->isFirstBlockIsOnBothTypeMachine()) reducedCost -= dualsValues[instance->getNbJobs() + 3];
            break;
        case TYPE_COLUMN::LS_MIN:reducedCost = -dualsValues[instance->getNbJobs() + 3];
            break;
    }

    unsigned int startingPairedIndex = elegantPair(indexStartingJobForDP, 0u);
    unsigned int lastPairedIndex = elegantPair(listGroupedJobs.size(), 0ul);
    // compute the range to visit jobs, i.e. (N - indexFirstJobForDP + 1)/ nbJobToSchedule
    unsigned int nbJobToExplore = listGroupedJobs.back().back().getIndex() - listGroupedJobs[indexStartingJobForDP].front().getIndex() + 1;
    unsigned int rangeOfSelection = std::max(1u, static_cast<unsigned int>(std::floor(static_cast<double>(nbJobToExplore) / static_cast<double>(nbJobsToSelect))));
    std::vector<unsigned int> newColumn = machineScheduleOnMachine;
    newColumn.reserve(nbJobsToSchedule);
    unsigned int completionTime = 0;
    for (unsigned int indexJobInCol: newColumn) {
        auto &job = instance->getListJobs()[indexJobInCol];
        reducedCost += beta(completionTime, job, machineSpeed);
        completionTime += static_cast<unsigned int>(job.getPi());
    }
    if (nbJobsToSelect == 0) return reducedCost;
    while (startingPairedIndex != lastPairedIndex) {

        double minEstimatedContribution = std::numeric_limits<double>::infinity();
        // use paired index for indentify the best jobs and for the loop
        unsigned int bestPairedIndex = startingPairedIndex;
        unsigned int currentPairedIndex = startingPairedIndex;
        unsigned int nextPairedIndex = currentPairedIndex;
        // if we have a last job to add, but with the rangeOfSelection, we do not cover all jobs, then we increase it with the good number of jobs
        if (newColumn.size() == nbJobsToSchedule - 1 && nbJobsToSelect * rangeOfSelection < nbJobToExplore) {
            rangeOfSelection += (nbJobToExplore - nbJobsToSelect * rangeOfSelection);
        }
        for (unsigned int itLoopNextPairedIndex = 0; itLoopNextPairedIndex < rangeOfSelection && nextPairedIndex != lastPairedIndex; itLoopNextPairedIndex++) {
            nextPairedIndex = nextIndex(listGroupedJobs, nextPairedIndex);
        }
        for (unsigned int nbExploredJobs = 0; nbExploredJobs < rangeOfSelection; nbExploredJobs++) {
            auto [indexGroup, indexInGroup] = elegantUnpair(currentPairedIndex);
            auto &currentJob = listGroupedJobs[indexGroup][indexInGroup];
            if (not node.isScheduledOnOtherMachines(currentJob.getIndex(), machineSpeed) && not node.isRemoved(currentJob.getIndex())) {
                unsigned int pj = static_cast<unsigned int>(currentJob.getPi());
                // compute the contribution of selection the job using the estimator
                double contributionSelectJob = beta(completionTime, listGroupedJobs[indexGroup][indexInGroup], machineSpeed);
                assert(nbJobsToSchedule >= newColumn.size() + 1);
                double estimationTakeJob = contributionSelectJob + estimationOfReducedCost(node, machineSpeed, nbJobsToSchedule - newColumn.size() - 1, nextPairedIndex, completionTime + pj);
                if (isSmaller(estimationTakeJob, minEstimatedContribution)) {
                    bestPairedIndex = currentPairedIndex;
                    minEstimatedContribution = estimationTakeJob;
                }
            }
            if (currentJob.getIndex() == listGroupedJobs.back().back().getIndex()) {
                currentPairedIndex = lastPairedIndex;
                break;
            } else currentPairedIndex = nextIndex(listGroupedJobs, currentPairedIndex);
        }
        auto [bestIndexGroup, bestIndexInGroup] = elegantUnpair(bestPairedIndex);
        auto &bestJob = listGroupedJobs[bestIndexGroup][bestIndexInGroup];
        newColumn.push_back(bestJob.getIndex());
        reducedCost += beta(completionTime, bestJob, machineSpeed);
        completionTime += static_cast<unsigned int>(bestJob.getPi());
        startingPairedIndex = nextPairedIndex;
    }

    // if the column have reduced cost, and it's not empty
    if (reducedCost < -EPSILON && !newColumn.empty()) {
        newSetOfColumns.emplace_back(typeColumn, std::move(newColumn));
        #ifdef DEBUG_CG
        double computeReducedCost = computeScheduleReducedCost(node, newSetOfColumns.back().second, machineSpeed, newSetOfColumns.back().first);
        if (std::abs(computeReducedCost - reducedCost) > EPSILON) {
            env.error() << "compute reduced cost failed : diff= " << std::abs(computeReducedCost - reducedCost) << std::endl;
            if (std::abs(computeReducedCost - reducedCost) > 0.01)
                throw BiSchException(std::string("The value of generated column does not correspond :").append(std::to_string(std::abs(computeReducedCost - reducedCost))).c_str());
        }
        assert(newSetOfColumns.back().second.size() == nbJobsToSchedule);
        #endif
    }
    return reducedCost;
}


inline double ColumnGeneration::estimationOfReducedCost(Node &node, char machineSpeed, unsigned int sizeOfWindows, unsigned int startingPairingIndex, unsigned int t) {
    double lowerBoundRedCost = 0.0;
    auto &listGroupedJobs = instance->getListGrpJobs();
    // use vector to know if the job are already explored
    std::fill(alreadyExplored.begin(), alreadyExplored.end(), false);
    auto [startingIndexGroup, startingIndexInGroup] = elegantUnpair(startingPairingIndex);
    // estimate the reduced cost in the window
    while (startingIndexGroup < listGroupedJobs.size() && sizeOfWindows > 0) {
        if (not node.isGroupIdenticalJobRemoved(startingIndexGroup)) {

            assert(startingIndexGroup < listGroupedJobs.size());
            assert(startingIndexInGroup < listGroupedJobs[startingIndexGroup].size());
            unsigned int pj = static_cast<unsigned int>(listGroupedJobs[startingIndexGroup][startingIndexInGroup].getPi());
            // find the minimal contribution among all jobs
            // the job with the minimal contribution to the reduced cost
            double minContribution = std::numeric_limits<double>::infinity();
            unsigned int minIndexGroup = 0;
            unsigned int minIndexJob = 0;
            for (unsigned int indexGroup = startingIndexGroup; indexGroup < listGroupedJobs.size(); ++indexGroup) {
                for (unsigned int indexJob = 0; indexJob < listGroupedJobs[indexGroup].size(); ++indexJob) {
                    auto &jobToSchedule = listGroupedJobs[indexGroup][indexJob];
                    // if job respect the constraints
                    double jobContribution = beta(t, jobToSchedule, machineSpeed);
                    if (jobContribution < minContribution) {
                        minContribution = jobContribution;
                        minIndexGroup = indexGroup;
                        minIndexJob = indexJob;
                    }
                }
            }

            alreadyExplored[listGroupedJobs[minIndexGroup][minIndexJob].getIndex()] = true;
            lowerBoundRedCost += minContribution;
            startingIndexGroup = minIndexGroup;
            startingIndexInGroup = minIndexJob;
            //        }
            t += pj;
            --sizeOfWindows;
        } else startingIndexGroup++;
    }
    return lowerBoundRedCost;
}

inline void ColumnGeneration::addColumnFromSolution(Solution &solution, Node &node, char machineSpeed) {
    unsigned int indexLoopMachine = (machineSpeed == 0) ? 0 : instance->getNbOfHighSpeedMachines();
    unsigned int indexMaxLoopOverSetColMachine = (machineSpeed == 0) ? instance->getNbOfHighSpeedMachines() : instance->getNbMachines();
    for (; indexLoopMachine < indexMaxLoopOverSetColMachine; ++indexLoopMachine) {
        //check if the column is not empty
        if (solution[indexLoopMachine].empty()) continue;
        // then add the column to the xs variable of the model
        auto schedule = static_cast<std::vector<unsigned int>>(solution[indexLoopMachine]);
        assert(std::adjacent_find(schedule.begin(), schedule.end()) == schedule.end());
        assert(!schedule.empty());
        TYPE_COLUMN typeCol =
                machineSpeed == 0 ? schedule.size() == instance->getMaxNbJobsOnHS() ? TYPE_COLUMN::HS_MAX : HS_MIN : schedule.size() == instance->getMaxNbJobsOnLS() ? TYPE_COLUMN::LS_MAX : LS_MIN;
        auto cost = computeScheduleCost(node, schedule, machineSpeed);
        unsigned int indexOnVariable = xs.getSize();
        xs.add(IloNumVar(env, 0.0, 1.0, ILOFLOAT));
        setMachineSchedules.emplace_back(std::move(schedule), indexOnVariable, 0, cost, typeCol,instance->getNbJobs());
    }
}

inline void ColumnGeneration::generateStartingColumns(Node &node) {
    if (heuristicSolver.getTimeUp() == nullptr) heuristicSolver.setTimeUp(timeUp);
    // start time to measure performance
    start = std::chrono::steady_clock::now();

    // add column from the upper bound solution
    addColumnFromSolution(*solution, node, 0);
    addColumnFromSolution(*solution, node, 1);
    auto listOfJobs = instance->getListJobs();
    std::vector<Job> selectedJobs;
    selectedJobs.reserve(instance->getNbToSelectJob());
    // use a windows of size n to get all subset of n jobs from the list of jobs and create optimal schedule for sum C_j
    unsigned int startIndexWindow = 0;
    while (startIndexWindow + instance->getNbToSelectJob() <= listOfJobs.size()) {
        isWithinTimeLimit();
        selectedJobs.clear();
        //loop over jobs to assign them
        for (unsigned int indexJob = startIndexWindow; indexJob < startIndexWindow + instance->getNbToSelectJob(); ++indexJob) {
            selectedJobs.push_back(listOfJobs[listOfJobs.size() - 1 - indexJob]);
        }
        Solution newSol = Solution::solveSumCjCriteria(selectedJobs, instance);
        heuristicSolver.upgradeSolutionWithHeuristic(newSol, instance->getListJobs());
        addColumnFromSolution(newSol, node, 0);
        addColumnFromSolution(newSol, node, 1);
        ++startIndexWindow;
    }

    unsigned int nbSolutionToGenerate = 0;
    //generate 'nbSolutionToGenerate' number of solution, we shuffle the set of jobs and take the 'n' first jobs to construct a solution
    std::mt19937 generator(0); // use seed '0' to have same results if we run twice the algorithm
    for (unsigned int loopOverShuffle = 0; loopOverShuffle < nbSolutionToGenerate; ++loopOverShuffle) {
        std::ranges::shuffle(listOfJobs, generator);
        std::copy_n(listOfJobs.begin(), instance->getNbToSelectJob(), selectedJobs.begin());
        std::sort(selectedJobs.begin(), selectedJobs.end(), std::greater<>());
        Solution newSol = Solution::solveSumCjCriteria(selectedJobs, instance);
        heuristicSolver.upgradeSolutionWithHeuristic(newSol, instance->getListJobs());
        addColumnFromSolution(newSol, node, 0);
        addColumnFromSolution(newSol, node, 1);
    }

}

inline void ColumnGeneration::initializeModel(Node &node) {

    /*************************************/
    /*      OBJECTIVE AND VARIABLES      */
    /*************************************/

    boost::dynamic_bitset<> encodingScheduledJob(instance->getNbJobs());
    for (auto &encodingMachine: node.getEncodingSelectedJobOnMachine()) {
        encodingScheduledJob |= encodingMachine;
    }
    IloExpr ExprObj(env);

    for (auto &infoColumn: setMachineSchedules) {
        double cost = infoColumn.costOfSchedule;
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        ExprObj += (xs[indexOnVariable] * cost);
    }

    ExprObj -= obj;
    objConstr = IloRange(env, 0.0, ExprObj, 0.0, "Objective");
    model.add(objConstr);
    ExprObj.end();

    /*************************/
    /*      Constraints      */
    /*************************/

    //constraint to compute the value of the weighted tardy jobs
    IloExpr ExprComputeUValue(env);


    for (auto &infoColumn: setMachineSchedules) {
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        ExprComputeUValue += (xs[indexOnVariable] * infoColumn.costOfSchedule);
    }


    ExprComputeUValue -= U;
    UComputeConstr = IloRange(env, 0.0, ExprComputeUValue, 0.0, "Compute_U_Value");
    model.add(UComputeConstr);
    ExprComputeUValue.end();

    // One job at most on machine schedule
    for (unsigned int j = 0; j < instance->getNbJobs(); ++j) {
        IloExpr ExprC1(env);

        for (auto &infoColumn: setMachineSchedules) {
            unsigned int indexOnVariable = infoColumn.indexOfVariable;
            double ajs = double(infoColumn.encodingColumn.test(j));
            ExprC1 += ajs * xs[indexOnVariable];
        }

        // if the job is deleted then set constraint == 0
        if (node.getEncodingRemoveDecision().test(j)) {
            constraints.add(IloRange(env, 0.0, ExprC1, 0.0, std::string("Job").append(std::to_string(j)).c_str()));
        }
            // else if we have schedule the job set the lower bound to 1
        else if (encodingScheduledJob.test(j)) {
            constraints.add(IloRange(env, 1.0, ExprC1, 1.0, std::string("Job").append(std::to_string(j)).c_str()));
        } else {
            constraints.add(IloRange(env, 0.0, ExprC1, 1.0, std::string("Job").append(std::to_string(j)).c_str()));
        }
        ExprC1.end();
    }

    // High speed machines :
    // set the constraint on number of machine with the max of jobs inside
    IloExpr ExprC2(env);
    for (auto &infoColumn: setMachineSchedules) {
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        double bs0p = infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX ? 1.0 : 0.0;
        ExprC2 += bs0p * xs[indexOnVariable];
    }

    // if the first block is not on both type machine then set the N^th constraint as
    // the one with lower bound on the number of HS machine with a maximum of jobs (last constraint in the paper in annexe)
    if (instance->isFirstBlockIsOnBothTypeMachine()) {
        constraints.add(IloRange(env, instance->getNbMachineWithMaxJobsOnHS(), ExprC2, instance->getNbOfHighSpeedMachines(), "M0+"));
    } else {
        constraints.add(IloRange(env, instance->getNbMachineWithMaxJobsOnHS(), ExprC2, instance->getNbMachineWithMaxJobsOnHS(), "M0+"));
    }

    ExprC2.end();

    // set the constraint on number of machine with the min of jobs inside
    IloExpr ExprC3(env);
    for (auto &infoColumn: setMachineSchedules) {
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        // need to adjust regarding the first block
        double bs0m = instance->isFirstBlockIsOnBothTypeMachine() ? infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN || infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX ? 1.0 : 0.0 :
                      infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN ? 1.0 : 0.0;
        ExprC3 += bs0m * xs[indexOnVariable];
    }
    // if the first block is not on both type machine then set the (N+1)^th constraint as
    // select M0 HS machines
    if (instance->isFirstBlockIsOnBothTypeMachine()) {
        constraints.add(IloRange(env, instance->getNbOfHighSpeedMachines(), ExprC3, instance->getNbOfHighSpeedMachines(), "M0-"));
    } else {
        constraints.add(IloRange(env, instance->getNbMachineWithMinJobsOnHS(), ExprC3, instance->getNbMachineWithMinJobsOnHS(), "M0-"));
    }

    ExprC3.end();

    // Low speed machines :
    // set the constraint on number of machine with the max of jobs inside
    IloExpr ExprC4(env);
    for (auto &infoColumn: setMachineSchedules) {
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        double bs1p = instance->isFirstBlockIsOnBothTypeMachine() ? infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::LS_MAX ? 1.0 : 0.0 :
                      infoColumn.typeOfColumn == TYPE_COLUMN::LS_MAX ? 1.0 : 0.0;
        ExprC4 += bs1p * xs[indexOnVariable];
    }

    // if the first block is not on both type machine then set the (N+2)^th constraint as
    // select Q number of HS and LS machines with a max of jobs
    if (instance->isFirstBlockIsOnBothTypeMachine()) {
        constraints.add(IloRange(env, instance->getNbJobsToScheduleOnFirstBlock(), ExprC4, instance->getNbJobsToScheduleOnFirstBlock(), "M+"));
    } else {
        constraints.add(IloRange(env, instance->getNbMachineWithMaxJobsOnLS(), ExprC4, instance->getNbMachineWithMaxJobsOnLS(), "M1+"));
    }

    ExprC4.end();

    // set the constraint on number of machine with the min of jobs inside
    IloExpr ExprC5(env);
    for (auto &infoColumn: setMachineSchedules) {
        unsigned int indexOnVariable = infoColumn.indexOfVariable;
        // need to adjust regarding the first block
        double bs1m = instance->isFirstBlockIsOnBothTypeMachine() ? infoColumn.typeOfColumn == TYPE_COLUMN::LS_MIN || infoColumn.typeOfColumn == TYPE_COLUMN::LS_MAX ? 1.0 : 0.0 :
                      infoColumn.typeOfColumn == TYPE_COLUMN::LS_MIN ? 1.0 : 0.0;
        ExprC5 += bs1m * xs[indexOnVariable];
    }

    // if the first block is not on both type machine then set the (N+3)^th constraint as
    // select M1 LS machines
    if (instance->isFirstBlockIsOnBothTypeMachine()) {
        // if we can have empty schedule for the low speed machines, then we just impose that nb min machines is between M1+ and M1
        if (instance->getMinNbJobsOnLS() == 0)
            constraints.add(IloRange(env, instance->getNbMachineWithMaxJobsOnLS(), ExprC5, instance->getNbOfLowSpeedMachines(), "M1-"));
        else
            constraints.add(IloRange(env, instance->getNbOfLowSpeedMachines(), ExprC5, instance->getNbOfLowSpeedMachines(), "M1-"));
    } else {
        constraints.add(IloRange(env, instance->getNbMachineWithMinJobsOnLS(), ExprC5, instance->getNbMachineWithMinJobsOnLS(), "M1-"));
    }
    ExprC5.end();
    model.add(constraints);

}

inline void ColumnGeneration::clearConstraintOfModel() {
    if (constraints.getSize() > 0) {
        model.remove(objConstr);
        model.remove(UComputeConstr);
        model.remove(constraints);
        constraints.clear();
    }
}

inline void ColumnGeneration::updateModel(Node &node) {

    if (checkColumns(node)) {
        #ifdef DEBUG_CG
        if (debug && verbose >= 3) {
            std::cout << "rebuild model: "<< std::endl;
        }
        #endif
        nbCleaningSetCol++;
        if (nbCleaningSetCol > 30000) {
            //end cplex env to free memory
            env.end();
            //recreate the new env
            env = IloEnv();
            nbCleaningSetCol = 0;
            model = IloModel(env);
            cplex = IloCplex(model);
            obj = IloNumVar(env, "obj");
            U = IloNumVar(env, "U");
            xs = IloNumVarArray(env);
            constraints = IloRangeArray(env);
            duals = IloNumArray(env);
            parametrize();
            // recreate variables
            for (auto itColumn = setMachineSchedules.begin(); itColumn != setMachineSchedules.end();) {
                unsigned int indexOnVariable = xs.getSize();
                xs.add(IloNumVar(env, 0.0, itColumn->UB_var, ILOFLOAT));
                itColumn->indexOfVariable = indexOnVariable;
                ++itColumn;
            }
            model.add(IloMinimize(env, obj));
            objConstr = IloRange(env, 0.0, 0.0, "Objective");
            model.add(objConstr);
            UComputeConstr = IloRange(env, 0.0, 0.0, "Compute_U_Value");
            model.add(UComputeConstr);

        }
        // clear all constraints of  the model
        clearConstraintOfModel();
        // initialize the master problem
        initializeModel(node);
    } else {
        // check deleted and scheduled jobs to update right-hand side
        for (unsigned int indexJob = 0; indexJob < instance->getNbJobs(); ++indexJob) {
            // if the job is deleted then set constraint == 0
            if (node.getEncodingRemoveDecision().test(indexJob)) {
                constraints[indexJob].setBounds(0.0, 0.0);
            } else {
                // check if a job is schedule regarding the branching scheme
                auto predCheckJobIsScheduled = [indexJob](const boost::dynamic_bitset<> &encodingMachineSchedule) { return encodingMachineSchedule.test(indexJob); };
                if (std::any_of(node.getEncodingSelectedJobOnMachine().begin(), node.getEncodingSelectedJobOnMachine().end(), predCheckJobIsScheduled)) {
                    constraints[indexJob].setBounds(1.0, 1.0);
                } else
                    constraints[indexJob].setBounds(0.0, 1.0);
            }
        }
    }
    // If during the branch and bound, we assign jobs on the first block and the first block is on both machines (low/high), then we can adjust the right-hand side of inequality
    if (instance->isFirstBlockIsOnBothTypeMachine()) {
        unsigned int nbHSMachineWithMaxJob = 0;
        // we can count the number of high speed machine with the max number of jobs, i.e. machine where there is a job on first position
        for (auto &[indexMachine, indexInMachine]: instance->getE()[0]) {
            //if we are in low speed machine then break
            if (indexMachine >= instance->getNbOfHighSpeedMachines()) break;
            if (indexMachine < instance->getNbOfHighSpeedMachines() && node.getBlockStruc()[indexMachine][indexInMachine].first != nullptr) {
                nbHSMachineWithMaxJob++;
            }
        }
        constraints[instance->getNbJobs()].setLB(nbHSMachineWithMaxJob);
    }
}

inline bool ColumnGeneration::checkColumns(Node &node) {

    // remove unused columns
    unsigned int nbUnfeasibleColumns = 0;
    // count nb columns
    unsigned int nbColumns = 0;

    // create a bitset that corresponding to all job that are scheduled regarding the branching scheme decisions
    boost::dynamic_bitset<> encodingSelectedJobOnMachineSpeed(instance->getNbJobs());

    for (const auto &encodingMachine: node.getEncodingSelectedJobOnMachine())
        encodingSelectedJobOnMachineSpeed |= encodingMachine;

    // here use the index for loop over machines that correspond to the given speed
    unsigned int indexLoopMachine = 0;
    unsigned int indexMaxLoopOverMachine = 0;
    bool existEmptyMachineHighSpeed = false;
    bool existEmptyMachineLowSpeed = false;
    auto emptyBitset = boost::dynamic_bitset<>(instance->getNbJobs());
    for (; indexLoopMachine < instance->getNbOfHighSpeedMachines(); indexLoopMachine++) {
        auto &encodingMachine = node.getEncodingSelectedJobOnMachine()[indexLoopMachine];
        if (encodingMachine == emptyBitset) {
            existEmptyMachineHighSpeed = true;
            break;
        }
    }
    for (indexLoopMachine = instance->getNbOfHighSpeedMachines(); indexLoopMachine < instance->getNbMachines(); indexLoopMachine++) {
        auto &encodingMachine = node.getEncodingSelectedJobOnMachine()[indexLoopMachine];
        if (encodingMachine == boost::dynamic_bitset<>(instance->getNbJobs())) {
            existEmptyMachineLowSpeed = true;
            break;
        }
    }

    #ifdef DEBUG_CG
    if (debug && verbose >= 3)  {
        Solution::printB(node.getBlockStruc());
        std::cout << "remove job: " << node.getEncodingRemoveDecision() << std::endl;
    }
    #endif

    IloNumArray newUBValues(env,xs.getSize());
    IloNumArray newLBValues(env,xs.getSize());
    for (auto &infoColumn: setMachineSchedules) {
        #ifdef DEBUG_CG
        if (debug && verbose >= 3) std::cout << "check col: " << infoColumn.column << " is on " << infoColumn.typeOfColumn << std::flush;
        #endif
        nbColumns++;
        indexLoopMachine = (infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN) ? 0 : instance->getNbOfHighSpeedMachines();
        indexMaxLoopOverMachine = (infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN) ? instance->getNbOfHighSpeedMachines() : instance->getNbMachines();
        bool existEmptyMachine = (infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN) ? existEmptyMachineHighSpeed : existEmptyMachineLowSpeed;
        // check if the column is unfeasible
        // first does it contain removed jobs ?
        bool isFeasible = true;
        if ((infoColumn.encodingColumn & node.getEncodingRemoveDecision()).any()) isFeasible = false;

        if (isFeasible) {
            // first we check if the column contain a partial schedule that come from the branching scheme
            bool containPartialSchedule = false;
            for (; indexLoopMachine < indexMaxLoopOverMachine; indexLoopMachine++) {
                // the column need to include a machine schedule given by the branching scheme decision
                auto &encodingMachine = node.getEncodingSelectedJobOnMachine()[indexLoopMachine];
                // if the column contains the same set of job that is present in the machine, and the machine is not empty, then compare it
                // I.e., check if the encoding of the machine is included in the encoding of the column
                if (encodingMachine != emptyBitset && (encodingMachine & infoColumn.encodingColumn) == encodingMachine) {
                    unsigned int indexInColumn = 0;
                    for (auto &jobInMachine: node.getBlockStruc()[indexLoopMachine]) {
                        if (indexInColumn < infoColumn.column.size()) {
                            if (jobInMachine.first != nullptr && infoColumn.column[indexInColumn] != jobInMachine.first->getIndex()) {
                                containPartialSchedule = false;
                                break;
                            }
                                // if we have the indexInColumn greater than 0, and the remaining job on machine are nullptr that me we have already compared all jobs on machine
                            else if (indexInColumn > 0 && jobInMachine.first == nullptr) {
                                containPartialSchedule = true;
                                break;
                            } // if indexInColumn == 0 and jobInMachine == null, that means we are in the first block and no jobs was assign to the machine
                            else if (indexInColumn == 0 && jobInMachine.first == nullptr) {
                                continue;
                            } else {
                                // else we have the same job in column and machine
                                ++indexInColumn;
                            }
                        } else {
                            // indexInColumn is equal to the size of the column that me we have already compared all jobs on machine, and it's fully included
                            containPartialSchedule = true;
                            break;
                        }
                    }
                    if (indexInColumn == infoColumn.column.size()) {
                        // indexInColumn is equal to the size of the column that me we have already compared all jobs on machine, and it's fully included
                        containPartialSchedule = true;
                        break;
                    }
                }
                // if the column contain the partial schedule, we stop the loop
                if (containPartialSchedule) break;
            }

            // if the column don't contain any partial schedule and there exist a list a machine with no jobs (on branching scheme)
            // then the column must contain jobs that are not coming from branching scheme
            if (!containPartialSchedule && existEmptyMachine && (encodingSelectedJobOnMachineSpeed & infoColumn.encodingColumn) == emptyBitset) {
                isFeasible = true;
            } else if (containPartialSchedule) {
                // we found that the column contains the partial schedule from 'indexLoopMachine', we now check that it don't have any jobs from other machine
                for (unsigned int indexLoopOtherMachine = 0; indexLoopOtherMachine < instance->getNbMachines() && isFeasible; ++indexLoopOtherMachine) {
                    if (containPartialSchedule && indexLoopOtherMachine == indexLoopMachine) continue;
                    auto &encodingMachine = node.getEncodingSelectedJobOnMachine()[indexLoopOtherMachine];
                    if ((infoColumn.encodingColumn & encodingMachine).any()) {
                        isFeasible = false;
                        break;
                    }
                }
            } else {
                // the column don't have a partial machine schedule, but there is no empty machine schedule so it's unfeasible
                isFeasible = false;
            }
            if (isFeasible) {
                #ifdef DEBUG_CG
                if (debug && verbose >= 3) std::cout << " feasible" << std::endl;
                #endif
                newUBValues[infoColumn.indexOfVariable] = 1.0;
                newLBValues[infoColumn.indexOfVariable] = 0.0;
                infoColumn.UB_var = 1.0;
                infoColumn.nbUnfeasible = 0;
            } else {
                #ifdef DEBUG_CG
                if (debug && verbose >= 3) std::cout << " NOT feasible" << std::endl;
                #endif
                ++infoColumn.nbUnfeasible;
                newUBValues[infoColumn.indexOfVariable] = 0.0;
                newLBValues[infoColumn.indexOfVariable] = 0.0;
                infoColumn.UB_var = 0.0;
                if (infoColumn.nbUnfeasible >= nbTimeNotUsed) {
                    ++nbUnfeasibleColumns;
                }
            }
        } else {
            #ifdef DEBUG_CG
            if (debug && verbose >= 3) std::cout << " NOT feasible" << std::endl;
            #endif
            ++infoColumn.nbUnfeasible;
            newUBValues[infoColumn.indexOfVariable] = 0.0;
            newLBValues[infoColumn.indexOfVariable] = 0.0;
            infoColumn.UB_var = 0.0;
            if (infoColumn.nbUnfeasible >= nbTimeNotUsed) {
                ++nbUnfeasibleColumns;
            }
        }
    }

    xs.setBounds(newLBValues,newUBValues);
    if (double(nbUnfeasibleColumns) / double(nbColumns) > thresholdSetCol) {

        // remove all variables
        xs.clear();
        model.remove(xs);
        auto pred_remove_not_used_col = [&](const ColumnGeneration::InfoColumn &infoColumn) {
            return infoColumn.nbUnfeasible >= nbTimeNotUsed;
        };

        auto itNotUsedCol = std::remove_if(setMachineSchedules.begin(), setMachineSchedules.end(), pred_remove_not_used_col);
        setMachineSchedules.erase(itNotUsedCol, setMachineSchedules.end());

        // create new variable

        for (auto itColumn = setMachineSchedules.begin(); itColumn != setMachineSchedules.end();) {
            unsigned int indexOnVariable = xs.getSize();
            xs.add(IloNumVar(env, 0.0, itColumn->UB_var, ILOFLOAT));
            itColumn->indexOfVariable = indexOnVariable;
            ++itColumn;
        }
        return true;
    }
    return false;
}

inline void ColumnGeneration::initRestrictedMasterProblem(Node &node) {
    // create column from the block Structure
    Solution::BlockStructure blockStructure = node.getBlockStruc();
    auto &E = instance->getE();
    unsigned int indexBlock = 0;
    double maxPjInBlock = -1.0;
    std::vector<unsigned int> alreadyUsedJobs;
    alreadyUsedJobs.reserve(instance->getNbToSelectJob());
    while (indexBlock < E.size()) {
        // the number of jobs that must be scheduled on the block k
        unsigned int numJobsToScheduleOnBlock = (indexBlock == 0) ? instance->getNbJobsToScheduleOnFirstBlock() : E[indexBlock].size();
        // get the maximum of the pj in the block and reduced numJobsToScheduleOnBlock if there is a job scheduled on machine
        for (auto &location: E[indexBlock]) {
            auto [job, cj] = blockStructure[location.first][location.second];
            if (job != nullptr) {
                maxPjInBlock = std::max(maxPjInBlock, job->getPi());
                alreadyUsedJobs.push_back(job->getIndex());
                --numJobsToScheduleOnBlock;
            }
        }
        // if we have not scheduled all jobs in the block, we add jobs with greater pj than the max
        for (unsigned int indexLoopLocation = 0; indexLoopLocation < E[indexBlock].size() && numJobsToScheduleOnBlock != 0; indexLoopLocation++) {
            auto &location = E[indexBlock][indexLoopLocation];
            auto [job, cj] = blockStructure[location.first][location.second];
            if (job == nullptr) {
                //find a job that can be schedule on the block
                double machineSpeed = (location.first < instance->getNbOfHighSpeedMachines()) ? 0 : 1;
                double speed = (machineSpeed == 0) ? instance->getHighSpeed() : instance->getLowSpeed();
                auto &listGroupedJobs = instance->getListGrpJobs();
                for (unsigned int indexLoopGrpJobs = 0; indexLoopGrpJobs < listGroupedJobs.size(); indexLoopGrpJobs++) {
                    if (node.isGroupIdenticalJobRemoved(indexLoopGrpJobs)) continue;
                    auto &groupOfIdenticalJob = listGroupedJobs[indexLoopGrpJobs];
                    // if the processing time is not enough higher
                    if (maxPjInBlock > groupOfIdenticalJob.back().getPi()) continue;
                    //else check if there exist at least a job that are not already scheduled
                    for (auto &jobInGroup: groupOfIdenticalJob) {
                        auto itFindJob = std::find(alreadyUsedJobs.begin(), alreadyUsedJobs.end(), jobInGroup.getIndex());
                        // if the job is not already scheduled then add it to the block structure
                        if (itFindJob == alreadyUsedJobs.end() && not node.isRemoved(jobInGroup.getIndex())) {
                            auto *newJob = &instance->getListJobs()[jobInGroup.getIndex()];
                            double newCj = (location.second > 0) ? blockStructure[location.first][location.second - 1].second + newJob->getPi() / speed : newJob->getPi() / speed;
                            blockStructure[location.first][location.second] = {newJob, newCj};
                            alreadyUsedJobs.push_back(jobInGroup.getIndex());
                            maxPjInBlock = std::max(maxPjInBlock, jobInGroup.getPi());
                            numJobsToScheduleOnBlock--; //decrease the number of jobs to schedule
                            break;
                        }
                    }
                    // if we add a job then stop
                    if (blockStructure[location.first][location.second].first != nullptr) break;
                }
            }
        }
        ++indexBlock;
    }
    #ifdef DEBUG_CG
    if (verbose >= 3 && debug) {
        std::cout << "Add column to make master feasible, solution found:" << std::endl;
        Solution::printB(blockStructure);
    }
    #endif

    // add each high-speed machine of this solution
    for (unsigned int indexLoopHSMachine = 0; indexLoopHSMachine < instance->getNbOfHighSpeedMachines(); indexLoopHSMachine++) {
        std::vector<unsigned int> newColumn;
        newColumn.reserve(instance->getMaxNbJobsOnHS());
        for (auto &jobAndCj: blockStructure[indexLoopHSMachine]) {
            // if there is a job, then add it to the column
            if (jobAndCj.first) {
                newColumn.push_back(jobAndCj.first->getIndex());
            }
        }
        // add the column on high speed machines
        auto cost = computeScheduleCost(node, newColumn, 0);
        unsigned int indexOnVariable = xs.getSize();
        IloNumArray coefConstraint(env, instance->getNbJobs() + 4);
        for (unsigned int indexJob: newColumn) {
            coefConstraint[indexJob] = 1.0;
        }
        TYPE_COLUMN typeColumn = newColumn.size() == instance->getMaxNbJobsOnHS() ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::HS_MIN;
        switch (typeColumn) {
            case TYPE_COLUMN::HS_MAX:coefConstraint[instance->getNbJobs()] = 1.0;
                if (instance->isFirstBlockIsOnBothTypeMachine()) {
                    coefConstraint[instance->getNbJobs() + 1] = 1.0;
                    coefConstraint[instance->getNbJobs() + 2] = 1.0;
                }
                break;
            case TYPE_COLUMN::HS_MIN:coefConstraint[instance->getNbJobs() + 1] = 1.0;
                break;
            case LS_MAX:break;
            case LS_MIN:break;
        }

        xs.add(IloNumVar(objConstr(cost) + UComputeConstr(cost) + constraints(coefConstraint), 0.0, 1.0, ILOFLOAT, "newVmax"));
        setMachineSchedules.emplace_back(std::move(newColumn), indexOnVariable, 0, cost, typeColumn,instance->getNbJobs());
    }

    // add each low-speed machine of this solution
    for (unsigned int indexLoopLSMachine = 0; indexLoopLSMachine < instance->getNbOfLowSpeedMachines(); indexLoopLSMachine++) {
        std::vector<unsigned int> newColumn;
        newColumn.reserve(instance->getMaxNbJobsOnLS());
        for (auto &jobAndCj: blockStructure[indexLoopLSMachine + instance->getNbOfHighSpeedMachines()]) {
            // if there is a job, then add it to the column
            if (jobAndCj.first) {
                newColumn.push_back(jobAndCj.first->getIndex());
            }
        }
        if (newColumn.empty()) continue;
        // add the column on high speed machines
        auto cost = computeScheduleCost(node, newColumn, 1);
        unsigned int indexOnVariable = xs.getSize();
        IloNumArray coefConstraint(env, instance->getNbJobs() + 4);
        for (unsigned int indexJob: newColumn) {
            coefConstraint[indexJob] = 1.0;
        }
        TYPE_COLUMN typeColumn = newColumn.size() == instance->getMaxNbJobsOnLS() ? TYPE_COLUMN::LS_MAX : TYPE_COLUMN::LS_MIN;
        switch (typeColumn) {
            case TYPE_COLUMN::HS_MAX:break;
            case TYPE_COLUMN::HS_MIN:break;
            case TYPE_COLUMN::LS_MAX:coefConstraint[instance->getNbJobs() + 2] = 1.0;
                if (instance->isFirstBlockIsOnBothTypeMachine()) coefConstraint[instance->getNbJobs() + 3] = 1.0;
                break;
            case TYPE_COLUMN::LS_MIN:coefConstraint[instance->getNbJobs() + 3] = 1.0;
                break;
        }
        xs.add(IloNumVar(objConstr(cost) + UComputeConstr(cost) + constraints(coefConstraint), 0.0, 1.0, ILOFLOAT, "newV0"));
        setMachineSchedules.emplace_back(std::move(newColumn), indexOnVariable, 0, cost, typeColumn,instance->getNbJobs());
    }
}

inline void ColumnGeneration::solve(Node &node) {
    try {
        bool haveGenerateColumn = false; // boolean to know if we need to generate a column to start the column generation
        lowerBound = node.getPartialSumWjUj();

        // if there is no columns, then generate them
        updateModel(node);

        bool continueGenerationHighSpeedColumn = true;
        bool continueGenerationLowSpeedColumn = true;
        // create the first state of the Dynamic Programming using the Branching Scheme decision
        std::vector<StateBacktracking> listStartingStateForHighSpeedDP;
        std::vector<StateBacktracking> listStartingStateForLowSpeedDP;
        unsigned int indexStartingJobForHighSpeedDP;
        unsigned int indexStartingJobForLowSpeedDP;
        computeFirstsStateOfDynamicProgramming(node, 0, listStartingStateForHighSpeedDP, indexStartingJobForHighSpeedDP);
        computeFirstsStateOfDynamicProgramming(node, 1, listStartingStateForLowSpeedDP, indexStartingJobForLowSpeedDP);
        // count the number of generation
        unsigned int nbGeneration = 1;
        do {
            isWithinTimeLimit();
            if (!cplex.solve()) {
                if (!haveGenerateColumn) {
                    #ifdef DEBUG_CG
                    if (debug) {
                        std::string modelName = instance->getInstancePath().parent_path().string();
                        modelName.append("/model_debug/model_unfeasible.lp");
                        std::filesystem::path fileSolutionPath = std::filesystem::path(modelName);
                        std::filesystem::create_directories(fileSolutionPath.lexically_normal().parent_path());
                        cplex.exportModel(modelName.c_str());
                        IloNumArray prefs(env, constraints.getSize());
                        for (unsigned int i = 0; i < constraints.getSize(); i++) {
                            prefs[i] = 1.0;
                        }
                        if (cplex.refineConflict(constraints,prefs)) {
                            modelName = instance->getInstancePath().parent_path().string();
                            modelName.append("/model_debug/model_unfeasible.clp");
                            cplex.writeConflict(modelName.c_str());
                        } else {
                            std::cerr << "Échec de l'affinage du conflit." << std::endl;
                        }
                    }
                    #endif
                    initRestrictedMasterProblem(node);
                    haveGenerateColumn = true;
                    continue;
                } else {
                    failedSolveMasterProblem = true;
                    #ifdef DEBUG_CG
                    if (debug) {
                        std::string modelName = instance->getInstancePath().parent_path().string();
                        modelName.append("/model_debug/model_unfeasible.lp");
                        std::filesystem::path fileSolutionPath = std::filesystem::path(modelName);
                        std::filesystem::create_directories(fileSolutionPath.lexically_normal().parent_path());
                        cplex.exportModel(modelName.c_str());
                        IloNumArray prefs(env, constraints.getSize());
                        for (unsigned int i = 0; i < constraints.getSize(); i++) {
                            prefs[i] = 1.0;
                        }
                        if (cplex.refineConflict(constraints,prefs)) {
                            modelName = instance->getInstancePath().parent_path().string();
                            modelName.append("/model_debug/model_unfeasible.clp");
                            cplex.writeConflict(modelName.c_str());
                        } else {
                            std::cerr << "Échec de l'affinage du conflit." << std::endl;
                        }
                    }
                    #endif
                    throw BiSchException("Erreur failed master");
                    break;
                }
            } else {
                #ifdef DEBUG_CG
                nbGen = nbGeneration;
                if (nbGeneration > 10000) break;
                #endif

                failedSolveMasterProblem = false;
                if (verbose >= 2)
                    std::cout << "Generation " << std::to_string(nbGeneration) << " minimum master problem : " << std::to_string(cplex.getObjValue()) << " LB: " << std::to_string(lowerBound)
                              << std::endl;
                double optValueRMP = cplex.getObjValue(); // optimal value of the restriced master problem
                // get dual values
                cplex.getDuals(duals, constraints);
                for (unsigned int indexDual = 0; indexDual < duals.getSize(); ++indexDual) {
                    double lambda_j = duals[indexDual];
                    dualsValues[indexDual] = lambda_j;
                }

                if (isSmaller(cplex.getObjValue(), lowerBound)) {
                    throw BiSchException(
                            std::string("Error LB > UB i.e., ").append(std::to_string(cplex.getObjValue())).append("<").append(std::to_string(lowerBound)).append("in column generation process"));
                }
                // if the diff UB - LB < 1 then we can stop the column generation
                if (isSmaller(optValueRMP - lowerBound, 1)) {
                    // since we are interested in the wj Uj value, which is an integer, then if UB - LB < 1 then we take LB = floor(UB)
                    // where U is the value of weighted sum of tardy jobs of the restricted master problem. Here, we add epsilon to avoid floating representation error
                    lowerBound = std::floor(cplex.getValue(U)+EPSILON);
                    break;
                }


                //***********************
                //      High Speed      *
                //***********************

                std::tuple<unsigned int,unsigned int,double,unsigned int,double> minRedCostsForHSMachines(0,0,std::numeric_limits<double>::infinity(),0,std::numeric_limits<double>::infinity());
                update_M_Vj_Qj_Values(node, 0);
                updateRedCostOfStartingStateOfDP(node, 0, listStartingStateForHighSpeedDP);
                continueGenerationHighSpeedColumn = generateColumns(node, 0, listStartingStateForHighSpeedDP, indexStartingJobForHighSpeedDP, minRedCostsForHSMachines);

                //***********************
                //      Low Speed      *
                //***********************
                std::tuple<unsigned int,unsigned int,double,unsigned int,double> minRedCostsForLSMachines(0,0,std::numeric_limits<double>::infinity(),0,std::numeric_limits<double>::infinity());
                update_M_Vj_Qj_Values(node, 1);
                updateRedCostOfStartingStateOfDP(node, 1, listStartingStateForLowSpeedDP);
                continueGenerationLowSpeedColumn = generateColumns(node, 1, listStartingStateForLowSpeedDP, indexStartingJobForLowSpeedDP, minRedCostsForLSMachines);
                // we can compute the LB, if we have solved the pricing problem for both kind of machines by the DP
                // we know that we have solve pricing pb for all machines when the sum of m_{g+/-} are equals to number of machines
                unsigned int checkNbMachines = std::get<0>(minRedCostsForHSMachines) + std::get<0>(minRedCostsForLSMachines) + std::get<1>(minRedCostsForHSMachines) + std::get<1>(minRedCostsForLSMachines) + std::get<3>(minRedCostsForHSMachines) + std::get<3>(minRedCostsForLSMachines);
                if (checkNbMachines > instance->getNbMachines()) throw BiSchException("Error in compute nb machines with solving pricing by DP");
                // can't have machine schedule generate for both max and min of jobs and in same time machine schedule generate for only a min of jobs
                if (std::get<0>(minRedCostsForHSMachines) > 0 && std::get<3>(minRedCostsForHSMachines) > 0) throw BiSchException("Error in compute nb high speed machine with m_{0+/-} machines with solving pricing by DP");
                if (std::get<0>(minRedCostsForLSMachines) > 0 && std::get<3>(minRedCostsForLSMachines) > 0) throw BiSchException("Error in compute nb high speed machine with m_{0+/-} machines with solving pricing by DP");
                if (checkNbMachines > instance->getNbMachines()) throw BiSchException("Error in compute nb machines with solving pricing by DP");
                if (checkNbMachines == instance->getNbMachines()) {
                    // if the first element of the tuple is not 0 => we have machines that we don't know if there is max or min jobs so take min between both
                    double lbRMP = optValueRMP;
                    // high speed machines
                    double minReducedCostHS = std::min(std::get<2>(minRedCostsForHSMachines), std::get<4>(minRedCostsForHSMachines));
                    double minReducedCostLS = std::min(std::get<2>(minRedCostsForLSMachines), std::get<4>(minRedCostsForLSMachines));
                    if (isSmaller(minReducedCostHS, 0.0)){
                        lbRMP += (static_cast<double>(std::get<0>(minRedCostsForHSMachines)) * minReducedCostHS);
                        if (isSmaller(std::get<2>(minRedCostsForHSMachines),0.0)) {
                            lbRMP += (static_cast<double>(std::get<1>(minRedCostsForHSMachines)) * std::get<2>(minRedCostsForHSMachines));
                        }if (isSmaller(std::get<4>(minRedCostsForHSMachines),0.0)) {
                            lbRMP += (static_cast<double>(std::get<3>(minRedCostsForHSMachines)) * std::get<4>(minRedCostsForHSMachines));
                        }
                    }if (isSmaller(minReducedCostLS, 0.0)){
                        lbRMP += (static_cast<double>(std::get<0>(minRedCostsForLSMachines)) * minReducedCostLS);
                        if (isSmaller(std::get<2>(minRedCostsForLSMachines),0.0)) {
                            lbRMP += (static_cast<double>(std::get<1>(minRedCostsForLSMachines)) * std::get<2>(minRedCostsForLSMachines));
                        }if (isSmaller(std::get<4>(minRedCostsForLSMachines),0.0)) {
                            lbRMP += (static_cast<double>(std::get<3>(minRedCostsForLSMachines)) * std::get<4>(minRedCostsForLSMachines));
                        }
                    }
                    lowerBound = std::max(lbRMP, lowerBound);
                }
                nbGeneration++;
            }
        } while (continueGenerationHighSpeedColumn || continueGenerationLowSpeedColumn);

        // if we have finished the column generation, set the lower bound as the value of the restricted master problem
        if (not continueGenerationHighSpeedColumn && not continueGenerationLowSpeedColumn) {
            lowerBound = cplex.getValue(U);
            solution->setSumWjUj(lowerBound);
        }
        isWithinTimeLimit();
        if (verbose >= 2) {
            if (failedSolveMasterProblem)
                std::cout << "No solution for master problem" << std::endl;
            else
                std::cout << "The objective value is : " << cplex.getObjValue() << std::endl << "The sum wj Uj is : " << getSumWjUj() << std::endl;
        }
        // reset the number of call of the node to zero
        node.setNbCallHeuristicFailed(0ul, 0);
        node.setNbCallHeuristicFailed(0ul, 1);
        #ifdef DEBUG_CG
        if (debug)  {
            double ownCost = 0.0;
            if (verbose >= 3) std::cout << "All high-speed machine schedules : " << std::endl;
            for (auto &infoColumn: setMachineSchedules) {
                if (infoColumn.typeOfColumn == TYPE_COLUMN::LS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::LS_MIN) continue;
                auto &schedule = infoColumn.column;
                unsigned int indexVariables = infoColumn.indexOfVariable;
                if (verbose >= 3) {
                    std::cout << "X_" << indexVariables << " : ";
                    if (!failedSolveMasterProblem) std::cout << cplex.getValue(xs[indexVariables]);
                }
                double cost = infoColumn.costOfSchedule;
                if (verbose >= 3){
                    std::cout << " | C:" << cost << " - > ";
                    for (unsigned int job: schedule) std::cout << std::to_string(job) << " ";
                    std::cout << " UB: " << xs[infoColumn.indexOfVariable].getUB();
                    std::cout << std::endl;
                }
                if (xs[infoColumn.indexOfVariable].getUB() != infoColumn.UB_var) throw BiSchException("Erreur UB variables");
                if (infoColumn.UB_var == 0) continue;
                ownCost += (cplex.getValue(xs[indexVariables]) * cost);
                auto redFromCplex = cplex.getReducedCost(xs[infoColumn.indexOfVariable]);
                auto ownRed = computeScheduleReducedCost(node,infoColumn.column,0,infoColumn.typeOfColumn);
                if (std::fabs(redFromCplex - ownRed) > EPSILON) {
                    throw BiSchException(std::string("Erreur red cost").append(std::to_string(redFromCplex - ownRed)).c_str());
                }
            }

            if (verbose >= 3) std::cout << std::endl << "All low-speed machine schedules : " << std::endl;
            for (auto &infoColumn: setMachineSchedules) {
                if (infoColumn.typeOfColumn == TYPE_COLUMN::HS_MAX || infoColumn.typeOfColumn == TYPE_COLUMN::HS_MIN) continue;
                auto &schedule = infoColumn.column;
                unsigned int indexVariables = infoColumn.indexOfVariable;
                if (verbose >= 3) {
                    std::cout << "X_" << (indexVariables + xs.getSize()) << " : ";
                    if (!failedSolveMasterProblem) std::cout << cplex.getValue(xs[indexVariables]);
                }
                double cost = infoColumn.costOfSchedule;
                if (verbose >= 3) {
                    std::cout << " | C:" << cost << " - > ";
                    for (unsigned int job: schedule) std::cout << std::to_string(job) << " ";
                    std::cout << " UB: " << xs[infoColumn.indexOfVariable].getUB();
                    std::cout << std::endl;
                }
                if (xs[infoColumn.indexOfVariable].getUB() != infoColumn.UB_var) throw BiSchException("Erreur UB variables");
                if (infoColumn.UB_var == 0) continue;
                ownCost += (cplex.getValue(xs[indexVariables]) * cost);
                auto redFromCplex = cplex.getReducedCost(xs[infoColumn.indexOfVariable]);
                auto ownRed = computeScheduleReducedCost(node,infoColumn.column,1,infoColumn.typeOfColumn);
                if (std::fabs(redFromCplex - ownRed) > EPSILON) {
                    throw BiSchException(std::string("Erreur red cost").append(std::to_string(redFromCplex - ownRed)).c_str());
                }
            }
            if (not isEqual(ownCost, cplex.getObjValue())) {
                throw BiSchException("Error with cplex");
            }
            if (nbGeneration > 100000) {
                throw BiSchException("Stop cause too many iteration");
            }
        }
        #endif
    } catch (const IloException &e) {
        env.error() << e;
        throw; // if you like
    }
}

inline std::tuple<unsigned int,unsigned int,double,unsigned int,double> ColumnGeneration::solvePricingProblem(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int indexStartingJobForDP) {
    std::tuple<unsigned int,unsigned int,double,unsigned int,double> results(0,0,std::numeric_limits<double>::infinity(),0,std::numeric_limits<double>::infinity());
    // clear memoization
    clearMemo();
    ++nbCallsDP;
    for (auto &[completionTime, nbSelectedJobs, redCostOfSchedule, indexLoopMachine, typeMachineSchedule]: listStartingStates) {
        double redCostCol = std::numeric_limits<double>::infinity();
        unsigned int maxNbJobs;
        if (typeMachineSchedule == 2) std::get<0>(results) += 1;
        if (typeMachineSchedule == 1 || typeMachineSchedule == 2) {
            if (typeMachineSchedule == 1) std::get<1>(results) += 1;
            // try with max nb of jobs on it
            maxNbJobs = machineSpeed == 0 ? instance->getMaxNbJobsOnHS() : instance->getMaxNbJobsOnLS();
            if (nbSelectedJobs <= maxNbJobs) {
                TYPE_COLUMN typeColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MAX : TYPE_COLUMN::LS_MAX;
                redCostCol = computeMinReducedCost(maxNbJobs - nbSelectedJobs, typeColumn, completionTime, indexStartingJobForDP, node) + redCostOfSchedule;
                redCostCol -= (machineSpeed == 0) ? dualsValues[instance->getNbJobs()] : dualsValues[instance->getNbJobs() + 2];
                if (instance->isFirstBlockIsOnBothTypeMachine()) {
                    redCostCol -= dualsValues[instance->getNbJobs() + 1 + machineSpeed * 2];
                    redCostCol -= machineSpeed == 0 ? dualsValues[instance->getNbJobs() + 2] : 0.0;
                }
                #ifdef DEBUG_MIP_PRICING
                std::vector<unsigned int> bestCol;
                double minFromMIPPricing = computeMinReducedCost(node, typeColumn, indexLoopMachine, bestCol);
                assert(bestCol.size() == maxNbJobs);
                if (not isEqual(minFromMIPPricing, redCostCol, EPSILON)) {
                    throw BiSchException("Found better col with MIP Pricing");
                }
                #endif
                std::get<2>(results) = std::min(redCostCol, std::get<2>(results));
            }
        }
        if (typeMachineSchedule == 0 || typeMachineSchedule == 2) {
            if (typeMachineSchedule == 0) std::get<3>(results) += 1;
            TYPE_COLUMN typeColumn = machineSpeed == 0 ? TYPE_COLUMN::HS_MIN : TYPE_COLUMN::LS_MIN;
            // try with min nb of jobs on it
            maxNbJobs = machineSpeed == 0 ? instance->getMinNbJobsOnHS() : instance->getMinNbJobsOnLS();
            if (nbSelectedJobs <= maxNbJobs) {
                redCostCol = computeMinReducedCost(maxNbJobs - nbSelectedJobs, typeColumn, completionTime, indexStartingJobForDP, node) + redCostOfSchedule;
                redCostCol -= (machineSpeed == 0) ? dualsValues[instance->getNbJobs() + 1] : dualsValues[instance->getNbJobs() + 3];
            }
            #ifdef DEBUG_MIP_PRICING
            std::vector<unsigned int> bestCol;
            double minFromMIPPricing = computeMinReducedCost(node, typeColumn, indexLoopMachine, bestCol);
            assert(bestCol.size() == maxNbJobs);
            if (not isEqual(minFromMIPPricing, redCostCol, EPSILON)) {
                throw BiSchException("Found better col with MIP Pricing");
            }
            #endif
            std::get<4>(results) = std::min(redCostCol, std::get<4>(results));
        }
    }
    return results;
}
