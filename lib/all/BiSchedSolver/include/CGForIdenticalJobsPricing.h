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
// Created by schau on 10/22/24.
//

#ifndef BILEVEL_SCHEDULING_CGFORIDENTICALJOBSPRINCING_H
#define BILEVEL_SCHEDULING_CGFORIDENTICALJOBSPRINCING_H

#include <algorithm>
#include "ColumnGeneration.h"

struct ConstForIdenticalJobs {
    std::vector<double>::reverse_iterator itLoopMValue;
    std::vector<double>::reverse_iterator itLoopBeginMValue;
    std::vector<std::pair<double, unsigned int>> list_Vj_Values;
    std::vector<std::pair<double, unsigned int>> list_Qj_Values;
    std::vector<std::pair<double, unsigned int>> list_Vj_minus_Qj_Values;
    std::vector<std::pair<double, unsigned int>> list_weights;
    // a vector of boolean to know if the job in non-increasing weight order is on time or not
    std::vector<bool> isOnTime;
    // a vector that maps EDD order to vj vector order
    std::vector<unsigned int> matchingEDD_orderToVj_Desc_order;
    // a vector that maps EDD order to Qj vector order
    std::vector<unsigned int> matchingEDD_orderToQj_Desc_order;
    // a vector that maps Weight non-increasing order to EDD order
    std::vector<unsigned int> matchingWeightToEDDOrder;
    // a vector that maps EDD order to Weight non-increasing order
    std::vector<unsigned int> matchingEDDToWeightOrder;
    // use vector to know how many times we visit a job, when we visit it twice then we update the corresponding weight
    std::vector<unsigned int> nbTimeVisitedJob;
    std::vector<double> list_M_value;
    DSU *dsu;
    std::vector<JumpPoint> *listJumpPoint;

    inline ConstForIdenticalJobs() : dsu(nullptr), listJumpPoint(nullptr) {}

    inline ConstForIdenticalJobs(unsigned int nbIdenticalJobs, DSU *dsu, std::vector<JumpPoint> *listJumpPoint) : dsu(dsu), listJumpPoint(listJumpPoint) {
        list_Vj_Values = std::vector<std::pair<double, unsigned int>>();
        list_Vj_minus_Qj_Values = std::vector<std::pair<double, unsigned int>>();
        list_Qj_Values = std::vector<std::pair<double, unsigned int>>();
        list_weights = std::vector<std::pair<double, unsigned int>>(nbIdenticalJobs);
        list_M_value = std::vector<double>();
        isOnTime = std::vector<bool>(nbIdenticalJobs);
        matchingEDD_orderToVj_Desc_order = std::vector<unsigned int>(nbIdenticalJobs);
        matchingEDD_orderToQj_Desc_order = std::vector<unsigned int>(nbIdenticalJobs);
        matchingEDDToWeightOrder = std::vector<unsigned int>(nbIdenticalJobs);
        matchingWeightToEDDOrder = std::vector<unsigned int>(nbIdenticalJobs);
        nbTimeVisitedJob = std::vector<unsigned int>(nbIdenticalJobs);
        // reserve the size
        list_Vj_Values.reserve(nbIdenticalJobs);
        list_Vj_minus_Qj_Values.reserve(nbIdenticalJobs);
        list_Qj_Values.reserve(nbIdenticalJobs);
        list_M_value.reserve(4 * nbIdenticalJobs + nbIdenticalJobs * nbIdenticalJobs * 2 + 1); // there are at most 2*n + 1 values
        itLoopMValue = list_M_value.rbegin();
        itLoopBeginMValue = list_M_value.rbegin();
    }

    inline void clear(unsigned int nbIdenticalJobs) {
        list_M_value.assign(4 * nbIdenticalJobs + nbIdenticalJobs * nbIdenticalJobs * 2 + 1, 0.0);
        list_Vj_Values.assign(nbIdenticalJobs, {0, 0});
        list_Vj_minus_Qj_Values.assign(nbIdenticalJobs, {0, 0});
        list_Qj_Values.assign(nbIdenticalJobs, {0, 0});
        std::fill(matchingEDD_orderToQj_Desc_order.begin(), matchingEDD_orderToQj_Desc_order.end(), 0);
        std::fill(matchingEDD_orderToVj_Desc_order.begin(), matchingEDD_orderToVj_Desc_order.end(), 0);
    }

    inline void initIteratorLoopMValue() { itLoopMValue = itLoopBeginMValue; };

    static double increase_safely(double value, double delta = 1e-10) {
        return std::max(value + delta, std::nextafter(value, std::numeric_limits<double>::infinity()));
    }

    static double decrease_safely(double value, double delta = 1e-10) {
        return std::min(value - delta, std::nextafter(value, -std::numeric_limits<double>::infinity()));
    }

    /**
     * Compute constants use in Jippe algorithm.
     * @param groupIdJobs The group of identical jobs.
     * @param dualsValues The dual values.
     * @param instance A pointer to the current instance.
     */
    inline void
    computeConstants(const std::vector<Job> &groupIdJobs, const std::vector<double> &dualsValues) {
        unsigned int nbIdenticalJobs = groupIdJobs.size();
        clear(nbIdenticalJobs);

        assert(list_Vj_Values.size() == nbIdenticalJobs);
        assert(list_Vj_minus_Qj_Values.size() == nbIdenticalJobs);
        assert(list_Qj_Values.size() == nbIdenticalJobs);

        // compute each vj and qj
        for (unsigned int indexLoopJobInGrp = 0; indexLoopJobInGrp < nbIdenticalJobs; ++indexLoopJobInGrp) {
            auto const &job = groupIdJobs[indexLoopJobInGrp];
            list_Vj_Values[indexLoopJobInGrp] = {dualsValues[job.getIndex()], indexLoopJobInGrp};
            list_Vj_minus_Qj_Values[indexLoopJobInGrp] = {job.getWi(), indexLoopJobInGrp};
            // for same reasons, we change with little - epsilon the q_j value
            list_Qj_Values[indexLoopJobInGrp] = {dualsValues[job.getIndex()] - job.getWi(), indexLoopJobInGrp};
        }

        // sort vj by non-increasing value if there are equality sort by non-increasing wj value if there are equality sort by EDD
        std::sort(list_Vj_Values.begin(), list_Vj_Values.end(), [&](const std::pair<double, unsigned int> &left, const std::pair<double, unsigned int> &right) {
            return isEqual(left.first, right.first) ?
                   isEqual(groupIdJobs[left.second].getWi(), groupIdJobs[right.second].getWi()) ?
                   isSmaller(groupIdJobs[left.second].getDi(), groupIdJobs[right.second].getDi())
                                                                                                : isSmaller(groupIdJobs[right.second].getWi(), groupIdJobs[left.second].getWi())
                                                    : isSmaller(right.first, left.first);
        });
        // sort qj by non-increasing value if there are equality sort by non-decreasing wj value if there are equality sort by inv EDD
        std::sort(list_Qj_Values.begin(), list_Qj_Values.end(), [&](const std::pair<double, unsigned int> &left, const std::pair<double, unsigned int> &right) {
            return isEqual(left.first, right.first) ?
                   isEqual(groupIdJobs[left.second].getWi(), groupIdJobs[right.second].getWi()) ?
                   isSmaller(groupIdJobs[left.second].getDi(), groupIdJobs[right.second].getDi()) :
                   isSmaller(groupIdJobs[left.second].getWi(), groupIdJobs[right.second].getWi())
                                                    : isSmaller(right.first, left.first);
        });
        // sort vj-qj by non-increasing value and by non-decreasing index (i.e. EDD order) value if there are equality
        std::sort(list_Vj_minus_Qj_Values.begin(), list_Vj_minus_Qj_Values.end(), [&](const std::pair<double, unsigned int> &left, const std::pair<double, unsigned int> &right) {
            return isEqual(left.first, right.first) ?
                   left.second < right.second
                                                    : left.first > right.first;
        });

        // for stability reasons, we change all by adding an epsilon >0 (conserving the order), in order to get all pairwise distinct values for vj and qj,
        // to avoid some case where there exist a value M such as v_j - M = w_k' the weight of other job k
        // that can lead to instability because the algorithm can keep to many jobs or remove too many jobs
        unsigned int indexLoopVj = 0;
        unsigned int indexLoopQj = 0;
        while (indexLoopQj + indexLoopVj < 2 * nbIdenticalJobs) {
            double vjValue = indexLoopVj < list_Vj_Values.size() ? list_Vj_Values[indexLoopVj].first : std::numeric_limits<double>::infinity();
            double qjValue = indexLoopQj < list_Qj_Values.size() ? list_Qj_Values[indexLoopQj].first : std::numeric_limits<double>::infinity();
            // we reduce the value by the number of value unmodified times 3 (because we want to add +eps -eps for the M_value)
            unsigned int nbValueUnmodified = (nbIdenticalJobs - indexLoopVj) + (nbIdenticalJobs - indexLoopQj);
            // if the smallest value is vj
            if (isSmaller(vjValue, qjValue)) {
                //reduce the vj value nbValueUnmodified*3 times
                for (unsigned int iLoopDecrease = 0; iLoopDecrease < nbValueUnmodified * 3; ++iLoopDecrease) {
                    vjValue = increase_safely(vjValue, 1e-6);
                }
                list_Vj_Values[indexLoopVj].first = vjValue;
                //pass to the next vj
                indexLoopVj++;
            }
                //else if the smallest value is qj
            else if (isSmaller(qjValue, vjValue)) {
                //reduce the qj value nbValueUnmodified*3 times
                for (unsigned int iLoopDecrease = 0; iLoopDecrease < nbValueUnmodified * 3; ++iLoopDecrease) {
                    qjValue = increase_safely(qjValue, 1e-6);
                }
                list_Qj_Values[indexLoopQj].first = qjValue;
                //pass to the next vj
                indexLoopQj++;
            } else {
                //the both values are equals
                //reduce the vj value nbValueUnmodified*3 times
                for (unsigned int iLoopDecrease = 0; iLoopDecrease < nbValueUnmodified * 3; ++iLoopDecrease) {
                    vjValue = increase_safely(vjValue, 1e-6);
                }
                list_Vj_Values[indexLoopVj].first = vjValue;
                //pass to the next vj
                indexLoopVj++;
                //update the number of unmodified values
                nbValueUnmodified = (nbIdenticalJobs - indexLoopVj) + (nbIdenticalJobs - indexLoopQj);
                //reduce the qj value nbValueUnmodified*3 times
                for (unsigned int iLoopDecrease = 0; iLoopDecrease < nbValueUnmodified * 3; ++iLoopDecrease) {
                    qjValue = increase_safely(qjValue, 1e-6);
                }
                list_Qj_Values[indexLoopQj].first = qjValue;
                //pass to the next vj
                indexLoopQj++;
            }
        }

        // create the map between EDD order and vj and qj vector order
        for (unsigned int indexLoopJobsInVector = 0; indexLoopJobsInVector < nbIdenticalJobs; ++indexLoopJobsInVector) {
            unsigned int indexJobInQjVectorOrder = list_Qj_Values[indexLoopJobsInVector].second;
            matchingEDD_orderToQj_Desc_order[indexJobInQjVectorOrder] = indexLoopJobsInVector;
            unsigned int indexJobInVjVectorOrder = list_Vj_Values[indexLoopJobsInVector].second;
            matchingEDD_orderToVj_Desc_order[indexJobInVjVectorOrder] = indexLoopJobsInVector;
        }

        // add the first value of M, i.e. the min of qj
        auto minQj = std::min_element(list_Qj_Values.begin(), list_Qj_Values.end());
        unsigned int indexLoopMValue = 0;
        assert(indexLoopMValue < list_M_value.size());
        list_M_value[indexLoopMValue] = minQj->first - 1;
        indexLoopMValue++;
        // add M values that come from the cross of w_j and w_i, used a high modification of 1e-6
        for (unsigned int indexLoopJobsInVector = 0; indexLoopJobsInVector < nbIdenticalJobs; ++indexLoopJobsInVector) {
            auto vjValueWithIndex = list_Vj_Values[indexLoopJobsInVector];
            for (unsigned int indexLoopOtherJobsInVector = 0; indexLoopOtherJobsInVector < nbIdenticalJobs; indexLoopOtherJobsInVector++) {
                auto vjMinusQjValueWithIndex = list_Vj_minus_Qj_Values[indexLoopOtherJobsInVector];
                if (vjValueWithIndex.second == vjMinusQjValueWithIndex.second) continue;

                assert(indexLoopMValue < list_M_value.size());
                // add + epsilon
                list_M_value[indexLoopMValue] = increase_safely(vjValueWithIndex.first - vjMinusQjValueWithIndex.first, 1e-6);
                indexLoopMValue++;
                assert(indexLoopMValue < list_M_value.size());
                // add - epsilon
                list_M_value[indexLoopMValue] = decrease_safely(vjValueWithIndex.first - vjMinusQjValueWithIndex.first, 1e-6);
                indexLoopMValue++;

            }
            auto qjValueWithIndex = list_Qj_Values[indexLoopJobsInVector];
            assert(indexLoopMValue < list_M_value.size());
            // add q_j + epsilon
            list_M_value[indexLoopMValue] = increase_safely(qjValueWithIndex.first);
            indexLoopMValue++;
            assert(indexLoopMValue < list_M_value.size());
            // add q_j - epsilon
            list_M_value[indexLoopMValue] = decrease_safely(qjValueWithIndex.first);
            indexLoopMValue++;
            assert(indexLoopMValue < list_M_value.size());
            // add v_j + epsilon
            list_M_value[indexLoopMValue] = increase_safely(vjValueWithIndex.first);
            indexLoopMValue++;
            assert(indexLoopMValue < list_M_value.size());
            // add v_j - epsilon
            list_M_value[indexLoopMValue] = decrease_safely(vjValueWithIndex.first);
            indexLoopMValue++;
        }
        std::sort(list_M_value.begin(), list_M_value.end(), [](const double &left, const double &right) { return left < right; });
        itLoopBeginMValue = std::make_reverse_iterator(std::unique(list_M_value.begin(), list_M_value.end()));

    }

    /**
     * Compute weight for given 'nbIdJob' identical jobs
     * @param nbIdJob The number of identical jobs
     * @param M_value The value considered for the procedure.
     * @param groupIdJobs The group of identical jobs.
     */
    inline void computeWeights(unsigned int nbIdJob, double M_value, const std::vector<Job> &groupIdJobs) {

        //clear all list because we use the value for M
        std::fill(nbTimeVisitedJob.begin(), nbTimeVisitedJob.end(), 0);

        // create the list of weight from the value of M
        // we loop over to list simultaneously
        unsigned int indexLoopOnVj = 0;
        unsigned int indexLoopOnVj_m_Qj = 0;

        // use index on list of weight to update it if is not same weight
        unsigned int indexOnWeight = 0;
        while (indexLoopOnVj < nbIdJob || indexLoopOnVj_m_Qj < nbIdJob) {
            //if both indexes are not at the end of the list of jobs
            if (indexLoopOnVj < nbIdJob && indexLoopOnVj_m_Qj < nbIdJob) {
                auto [valueFromVj, indexFromVj] = list_Vj_Values[indexLoopOnVj];
                auto [valueFromVjMinusQj, indexFromVjMinusQj] = list_Vj_minus_Qj_Values[indexLoopOnVj_m_Qj];
                // if both value are minimals
                if (isEqual(valueFromVj - M_value, valueFromVjMinusQj)) {
                    ++nbTimeVisitedJob[indexFromVj];
                    ++nbTimeVisitedJob[indexFromVjMinusQj];
                    // if we see for the second time a job then we add it
                    if (nbTimeVisitedJob[indexFromVj] == 2 && nbTimeVisitedJob[indexFromVjMinusQj] == 2) {
                        // jobs have same weight, so in case of equality, we sort weight by EDD order
                        unsigned int firstMinIndex = isSmaller(groupIdJobs[indexFromVj].getDi(), groupIdJobs[indexFromVjMinusQj].getDi()) ? indexFromVj : indexFromVjMinusQj;
                        // then the second min jobs
                        unsigned int secondMinIndex = isSmaller(groupIdJobs[indexFromVj].getDi(), groupIdJobs[indexFromVjMinusQj].getDi()) ? indexFromVjMinusQj : indexFromVj;
                        // the first min value
                        double firstMinValue = isSmaller(groupIdJobs[indexFromVj].getDi(), groupIdJobs[indexFromVjMinusQj].getDi()) ? valueFromVj - M_value : valueFromVjMinusQj;
                        double secondMinValue = isSmaller(groupIdJobs[indexFromVj].getDi(), groupIdJobs[indexFromVjMinusQj].getDi()) ? valueFromVjMinusQj : valueFromVj - M_value;
                        changeWeight(indexOnWeight, secondMinIndex, secondMinValue);
                        // if is not the same job, then add it
                        if (indexFromVj != indexFromVjMinusQj)
                            changeWeight(indexOnWeight, firstMinIndex, firstMinValue);
                    } else if (nbTimeVisitedJob[indexFromVj] == 2 || nbTimeVisitedJob[indexFromVjMinusQj] == 2) {
                        // if at least one is seen twice
                        unsigned int minIndex = nbTimeVisitedJob[indexFromVj] == 2 ? indexFromVj : indexFromVjMinusQj;
                        double minValue =
                                nbTimeVisitedJob[indexFromVj] == 2 ? valueFromVj - M_value : valueFromVjMinusQj;
                        changeWeight(indexOnWeight, minIndex, minValue);
                    }
                    ++indexLoopOnVj;
                    ++indexLoopOnVj_m_Qj;
                } else if (valueFromVjMinusQj < valueFromVj - M_value) {
                    // if the minimal come from the v_j - q_j
                    ++nbTimeVisitedJob[indexFromVj];
                    // if at least one is seen twice
                    if (nbTimeVisitedJob[indexFromVj] == 2)
                        changeWeight(indexOnWeight, indexFromVj, valueFromVj - M_value);
                    ++indexLoopOnVj;

                } else if (valueFromVj - M_value < valueFromVjMinusQj) {
                    // if the min come from v_j - M
                    ++nbTimeVisitedJob[indexFromVjMinusQj];
                    // if at least one is seen twice
                    if (nbTimeVisitedJob[indexFromVjMinusQj] == 2)
                        changeWeight(indexOnWeight, indexFromVjMinusQj, valueFromVjMinusQj);
                    ++indexLoopOnVj_m_Qj;
                }
            } else {
                //else at least one index is at the end of its list of jobs. So we just
                auto [minValue, minIndex] = indexLoopOnVj < nbIdJob ? list_Vj_Values[indexLoopOnVj]
                                                                    : list_Vj_minus_Qj_Values[indexLoopOnVj_m_Qj];
                // minus M_value if we use vj
                if (indexLoopOnVj < nbIdJob) minValue -= M_value;
                changeWeight(indexOnWeight, minIndex, minValue);
                ++nbTimeVisitedJob[minIndex];
                if (indexLoopOnVj < nbIdJob) ++indexLoopOnVj;
                else ++indexLoopOnVj_m_Qj;
            }
        }
        for (unsigned int indexOfJobInWeightOrder = 0;
             indexOfJobInWeightOrder < nbIdJob; ++indexOfJobInWeightOrder) {
            assert(indexOfJobInWeightOrder < list_weights.size());
            auto &[wj, indexJ] = list_weights[indexOfJobInWeightOrder];
            matchingWeightToEDDOrder[indexOfJobInWeightOrder] = indexJ;
            matchingEDDToWeightOrder[indexJ] = indexOfJobInWeightOrder;
        }
    }

    /**
     * Method that changes the value of the weight (if it's not the same) in the list of weights.
     *
     * @param indexInTheListOfWeight The index of the weight to change in the list.
     * @param indexOfJobOfTheWeight The index of the job used to compute the weight.
     * @param newWeightValue The new value if a change is needed.
     */
    inline void changeWeight(unsigned int &indexInTheListOfWeight, unsigned int indexOfJobOfTheWeight, double newWeightValue) {
        list_weights[indexInTheListOfWeight].first = newWeightValue;
        list_weights[indexInTheListOfWeight].second = indexOfJobOfTheWeight;
        ++indexInTheListOfWeight;
    };

    /**
     * Appy the procedure for a given M value.
     * @param M_value The value considered for the procedure.
     * @param t The completion time of the last scheduled job on the machine schedule. Must be divided by the speed to
     * get the real completion time
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listOfGroupedJobs A list containing the jobs that can be added.
     */
    inline void computeMinReducedCostIdenticalJobs(double M_value, unsigned int t, char machineSpeed, const Instance *instance, const std::vector<Job> &listOfGroupedJobs) {
        double speed = (machineSpeed == 0) ? instance->getHighSpeed() : instance->getLowSpeed();
        double p = listOfGroupedJobs.back().getPi() / speed;
        unsigned int nbIdJobs = listOfGroupedJobs.size();
        assert(dsu != nullptr);
        dsu->reset(nbIdJobs + 1);
        // we set the first jumpPoint with the previous equal to -1 to says it's the origin point
        (*listJumpPoint)[0] = JumpPoint(1, -1, 1);
        assert(nbIdJobs < listJumpPoint->size());
        //create the other list of jump points
        for (unsigned int indexJumpPoint = 1; indexJumpPoint < nbIdJobs + 1; ++indexJumpPoint) {
            (*listJumpPoint)[indexJumpPoint] = JumpPoint(indexJumpPoint + 1, indexJumpPoint - 1, 1);
        }

        for (unsigned int indexLoopOnWeight = 0; indexLoopOnWeight < nbIdJobs; ++indexLoopOnWeight) {
            // get the job in EDD order corresponding to the index of weight
            auto &job = listOfGroupedJobs[matchingWeightToEDDOrder[indexLoopOnWeight]];
            // if the value of the job is < 0 then we can pass to next job
            if (list_Vj_Values[matchingWeightToEDDOrder[indexLoopOnWeight]].first - M_value < 0.0) continue;
            double d = job.getDi() - t / speed;
            // if the job is late then pass next job
            if (d <= -EPSILON) {
                isOnTime[indexLoopOnWeight] = false;
                continue;
            }

            // if dj is greater than all t then we can add the job add mt doesn't change
            unsigned int indexInDsu = std::floor(d / p);
            auto t1 = dsu->find_parent(indexInDsu);
            auto mt1 = dsu->minSet[t1];
            assert(mt1 < listJumpPoint->size());
            // get the corresponding jump point
            auto jumpPoint = (*listJumpPoint)[mt1];
            // if the current jump point is the origin, i.e. the previous index of the jump point is -1
            if (jumpPoint.previousT == -1) {
                isOnTime[indexLoopOnWeight] = false;
                continue;
            } else {
                --jumpPoint.difference;
                if (jumpPoint.difference == 0) {
                    //merge both segment in DSU, i.e. the one from mt1 and the previous segment given by mt1 - 1
                    auto t0 = dsu->find_parent(mt1 - 1);
                    dsu->union_sets(t1, t0);
                    (*listJumpPoint)[jumpPoint.previousT].nextT = jumpPoint.nextT;
                    if (jumpPoint.nextT < listJumpPoint->size())
                        (*listJumpPoint)[jumpPoint.nextT].previousT = jumpPoint.previousT;
                }
            }
        }
    }

    /**
     * Constructs the schedule by executing the jobs in the on-time set in order of the earliest due date (EDD).
     * If a job is tardy and has a positive weight, it will be scheduled at the end of the column.
     * @param onTimeJobInEDD The set of on-time jobs in EDD order
     * @param M_value The value considered for the procedure
     * @param nbIdJob The number of identical jobs
     */
    inline void constructSchedule(std::vector<unsigned int> &onTimeJobInEDD, double M_value, unsigned int nbIdJob) {
        onTimeJobInEDD.reserve(nbIdJob);
        //use index to know where the on time jobs are in the current schedule
        unsigned int indexOnTimeJob = 0;
        // iterate through the jobs in order of EDD while knowing if they are on time or not thanks to the matching between EDD order to Weight Order
        for (unsigned int indexLoopJobInEDD = 0; indexLoopJobInEDD < nbIdJob; ++indexLoopJobInEDD) {
            assert(indexLoopJobInEDD < matchingWeightToEDDOrder.size());
            assert(indexLoopJobInEDD < matchingEDDToWeightOrder.size());
            assert(indexLoopJobInEDD < matchingEDD_orderToVj_Desc_order.size());
            assert(indexLoopJobInEDD < matchingEDD_orderToQj_Desc_order.size());
            // get the index of the job when sorted on weight
            unsigned int indexLoopJobInWeight = matchingEDDToWeightOrder[indexLoopJobInEDD];
            unsigned int indexLoopJobInVjOrder = matchingEDD_orderToVj_Desc_order[indexLoopJobInEDD];
            // if the job has a negative value when is on time then pass to next one
            if (list_Vj_Values[indexLoopJobInVjOrder].first - M_value < 0.0) continue;

            // if the job is on time, then schedule it in EDD order, here isOnTime is sorted according the weight,
            //so we need to use the indexLoopJobInWeight and indexLoopJobInQjOrder to get the value of Qj
            unsigned indexLoopJobInQjOrder = matchingEDD_orderToQj_Desc_order[indexLoopJobInEDD];
            if (isOnTime[indexLoopJobInWeight]) {
                onTimeJobInEDD.insert(onTimeJobInEDD.begin() + indexOnTimeJob, indexLoopJobInEDD);
                ++indexOnTimeJob;
            }
                // else if the job is tardy, and it has positive value when tardy then schedule it
            else if (0.0 <= list_Qj_Values[indexLoopJobInQjOrder].first - M_value) {
                onTimeJobInEDD.emplace_back(indexLoopJobInEDD); // else schedule it as tardy job
            }
        }
    }
};

// add ww operator to display list v_j and q_j value
template<typename T1, typename T2>
inline std::ostream &operator<<(std::ostream &os, const std::vector<std::pair<T1, T2>> &vector) {

    os << "[";
    if (vector.empty()) os << "]";
    else {
        for (size_t indexLoopVector = 0; indexLoopVector < vector.size() - 1; ++indexLoopVector) {
            const T1 &firstElement = vector[indexLoopVector].first;
            const T2 &secondElement = vector[indexLoopVector].second;
            os << "(" << firstElement << "," << secondElement << "),";
        }
        const T1 &firstElement = vector.back().first;
        const T2 &secondElement = vector.back().second;
        os << "(" << firstElement << "," << secondElement << ")]";
    }
    return os;
}


#endif //BILEVEL_SCHEDULING_CGFORIDENTICALJOBSPRINCING_H
