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
// Created by schau on 9/11/25.
//
#ifndef BILEVEL_SCHEDULING_NEIGHBORHOOD_H
#define BILEVEL_SCHEDULING_NEIGHBORHOOD_H
#include <iterator>
#include "Solution.h"

// uncomment this line to debug the Class
// #define DEBUG_NEIGHBORHOOD

/**
 * Class that implement an interface for make a neighborhood iterator
 */
class INeighborhoodIterator {
protected:
    const Solution::BlockStructure* baseBlockStructure;
    Instance* instance;
    unsigned int indexBlock = 0;

public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Solution::BlockStructure;
    using difference_type = std::ptrdiff_t;
    using pointer = const Solution::BlockStructure*;
    using reference = Solution::BlockStructure;
#ifdef DEBUG_NEIGHBORHOOD
    unsigned int nbIter = 0;
#endif
    explicit INeighborhoodIterator() : baseBlockStructure(nullptr), instance(nullptr) {}

    explicit INeighborhoodIterator(const Solution::BlockStructure* baseBlockStructure,
                                   Instance* instance) : baseBlockStructure(baseBlockStructure), instance(instance) {}

    virtual ~INeighborhoodIterator() = default;

    [[nodiscard]] virtual bool equalsBase(const INeighborhoodIterator& other) const {
        return other.indexBlock == indexBlock
                && other.instance == instance
                && other.baseBlockStructure == baseBlockStructure;
    }

    /**
     * Method that compute and update the completion on the indexMachine from the indexInMachine until the end
     * @param indexInMachine The index of the modified job in the machine.
     * @param indexMachine The index of the machine.
     * @param newBlockStructure The newBlockStructure to update.
     */
    void updateCompletionTimeOnMachineAndBlock(unsigned int indexInMachine, unsigned int indexMachine,
                                               Solution::BlockStructure& newBlockStructure) const {
        // update the completion time
        double completionTime = indexInMachine == 0 ? 0.0 : newBlockStructure[indexMachine][indexInMachine - 1].second;
        double speed = indexMachine < instance->getNbOfHighSpeedMachines()
                           ? instance->getHighSpeed()
                           : instance->getLowSpeed();
        for (unsigned int indexLoopInMachine = indexInMachine; indexLoopInMachine < newBlockStructure[indexMachine].
             size(); indexLoopInMachine++) {
            auto& jobWithCj = newBlockStructure[indexMachine][indexLoopInMachine];
            if (jobWithCj.first != nullptr) {
                completionTime += jobWithCj.first->getPi() / speed;
                jobWithCj.second = completionTime;
            }else if (indexLoopInMachine == 0) jobWithCj.second = 0.0;
                // we can have the first jobs that is NULL pointer because we don't have fulfilled the first block
            else break;
        }
    }

    /***************
     *  Getters
     ***************/

    [[nodiscard]] unsigned int getIndexBlock() const { return indexBlock; }

    /***************
     *  Setters
     ***************/

    void setIndexBlock(unsigned int indexBlock) { this->indexBlock = indexBlock; }

};

class LeaderNeighborhoodIterator : public INeighborhoodIterator {
public:
    explicit LeaderNeighborhoodIterator() : listJobToSwap(nullptr) {}

    explicit LeaderNeighborhoodIterator(Instance* instance, const Solution::BlockStructure* block_structure,
                                        const std::vector<unsigned int>*
                                        listJobToSwapInSolution) : INeighborhoodIterator(block_structure, instance),
                                                                   listJobToSwap(listJobToSwapInSolution) {
        initialize();
    }

    void initialize() {
        Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance, *baseBlockStructure, indexBlock, maxPj, minPj);
    }

    void updateFirstIndexBlock() {
        // if there is list of job to swap, then find the first block where we can insert it
        if (not listJobToSwap->empty()) {
            Job jobToSwap = instance->getListJobs()[(*listJobToSwap)[indexJobToSwapFromList]];
            // pass to the block
            while (not(isSmallerOrEqual(maxPj, jobToSwap.getPi()) && isSmallerOrEqual(jobToSwap.getPi(), minPj)) &&
                indexBlock < instance->getE().size()) {
                ++indexBlock;
                minPj = std::numeric_limits<double>::infinity();
                Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance, *baseBlockStructure, indexBlock, maxPj,
                                                                   minPj);
            }
            // now we find the first position in the block where there is a job
            positionInBlock = 0;
            while (positionInBlock < static_cast<int>(instance->getE()[indexBlock].size())) {
                auto& [indexMachine,indexInMachine] = instance->getE()[indexBlock][positionInBlock];
                if ((*baseBlockStructure)[indexMachine][indexInMachine].first != nullptr) break;
                ++positionInBlock;
            }
            // if there is no job, make the begin iterator as end iterator
            if (positionInBlock >= static_cast<int>(instance->getE()[indexBlock].size())) {
                indexJobToSwapFromList = listJobToSwap->size();
                positionInBlock = 0;
            }
        }
    }

    void nextNeighbor() {
#ifdef DEBUG_NEIGHBORHOOD
        ++nbIter;
#endif
        bool haveFoundJob = false;
        while (not haveFoundJob) {
            if (indexJobToSwapFromList < listJobToSwap->size()) {
                Job jobToSwap = instance->getListJobs()[(*listJobToSwap)[indexJobToSwapFromList]];
                // if the job to swap can go in the current block
                if (isSmallerOrEqual(maxPj, jobToSwap.getPi()) && isSmallerOrEqual(jobToSwap.getPi(), minPj)) {
                    // if there is no scheduled job, pass to the next position
                    bool thereIsScheduledJobToSwap = false;
                    while (not thereIsScheduledJobToSwap && positionInBlock < static_cast<int>(instance->getE()[indexBlock].size() - 1)) {
                        auto& [indexMachine,indexInMachine] = instance->getE()[indexBlock][positionInBlock+1];
                        thereIsScheduledJobToSwap = (*baseBlockStructure)[indexMachine][indexInMachine].first !=
                                nullptr;
                        if (not thereIsScheduledJobToSwap) { ++positionInBlock; }
                    }
                    // if we are at the end of a block, then we need to pass to the next job
                    if (positionInBlock >= static_cast<int>(instance->getE()[indexBlock].size() - 1)) {
                        positionInBlock = -1;
                        ++indexJobToSwapFromList;
                    }else {
                        // pass to the position in the block
                        ++positionInBlock;
                        haveFoundJob = true;
                    }
                }else {
                    // pass to the block
                    while (not(isSmallerOrEqual(maxPj, jobToSwap.getPi()) && isSmallerOrEqual(jobToSwap.getPi(), minPj))
                        && indexBlock < instance->getE().size()) {
                        ++indexBlock;
                        minPj = std::numeric_limits<double>::infinity();
                        Solution::computeMaxMinProcessingTimeOnNearbyBlock(
                            instance, *baseBlockStructure, indexBlock, maxPj, minPj);
                    }
                }
                // if we are at the end of the block structure
                if (indexBlock >= instance->getE().size()) {
                    indexBlock = 0;
                    indexJobToSwapFromList = listJobToSwap->size();
                }
            }else {
                indexBlock = 0; // set indexBlock = 0 to match the .end() iterator
                break;
            }
        }
    }

    /***************
     *  Setters
     ***************/

    void setIndexJobToSwapFromList(unsigned int indexJobToSwapFromList) { this->indexJobToSwapFromList = indexJobToSwapFromList; }

    /**************
     *  Getters
     **************/

    [[nodiscard]] int getIndexLastSwapJobFromBaseBlockStruct() const { return indexLastSwapJobFromBaseBlockStruct; }

    [[nodiscard]] unsigned int getIndexJobToSwapFromList() const { return indexJobToSwapFromList; }

    /***************
     *  Operators
     ***************/

    Solution::BlockStructure operator*() {
        if (indexJobToSwapFromList >= listJobToSwap->size()) {
            throw std::out_of_range(
                "Leader's NeighborHood out of range, indexJobToSwapFromList >= listJobToSwap.size()");
        }
        Solution::BlockStructure newBlockStructure = *baseBlockStructure;
        auto& [indexMachine,indexInMachine] = instance->getE()[indexBlock][positionInBlock];
        //swap the jobs
        assert(newBlockStructure[indexMachine][indexInMachine].first != nullptr);
        indexLastSwapJobFromBaseBlockStruct = newBlockStructure[indexMachine][indexInMachine].first->getIndex();
        newBlockStructure[indexMachine][indexInMachine].first = &instance->getListJobs()[(*listJobToSwap)[indexJobToSwapFromList]];
        // update the completion time
        updateCompletionTimeOnMachineAndBlock(indexInMachine, indexMachine, newBlockStructure);
        return newBlockStructure;
    }

    LeaderNeighborhoodIterator& operator++() {
        nextNeighbor();
        return *this;
    }

    LeaderNeighborhoodIterator operator++(int) {
        LeaderNeighborhoodIterator tmp = *this;
        nextNeighbor();
        return tmp;
    }

    bool operator==(const LeaderNeighborhoodIterator& other) const {
        return other.indexJobToSwapFromList == indexJobToSwapFromList
                && other.listJobToSwap == listJobToSwap
                && other.positionInBlock == positionInBlock
                && equalsBase(other);
    }

    bool operator!=(const LeaderNeighborhoodIterator& other) const { return !(*this == other); }

private:
    unsigned int indexLastSwapJobFromBaseBlockStruct = 0;
    unsigned int indexJobToSwapFromList = 0;
    const std::vector<unsigned int>* listJobToSwap;
    int positionInBlock = -1;
    // get the max pj from the block before and the min from the block after
    double maxPj = -1.0, minPj = std::numeric_limits<double>::infinity();
};

class FollowerNeighborhoodIterator : public INeighborhoodIterator {
public:
    FollowerNeighborhoodIterator() = default;

    FollowerNeighborhoodIterator(Instance* instance,
                                 const Solution::BlockStructure* block_structure) : INeighborhoodIterator(
        block_structure, instance) {}

    FollowerNeighborhoodIterator(Instance* instance, const Solution::BlockStructure* block_structure,
                                 unsigned int nbElementsToArrange) : INeighborhoodIterator(block_structure, instance),
                                                                     nbElementsToArrange(nbElementsToArrange),
                                                                     currentPermutation() { initialize(); }

    void initialize() {

        bool isEmpty = true;
        assert(baseBlockStructure != nullptr);
        assert(indexBlock < instance->getE().size());
        // check if the first block is not empty
        for (auto& [indexMachine,indexInMachine] : instance->getE()[indexBlock]) {
            assert(indexMachine < (*baseBlockStructure).size());
            assert(indexInMachine < (*baseBlockStructure)[indexMachine].size());
            if ((*baseBlockStructure)[indexMachine][indexInMachine].first != nullptr) {
                isEmpty = false;
                break;
            }
        }
        currentPermutation.clear();
        allTranspositions.clear();
        allTranspositions.reserve(instance->getNbMachines());
        haveBeenPermuted = false;
        if (not isEmpty) nextNeighbor();
        else indexBlock = instance->getE().size();
    }

    void nextNeighbor() {
#ifdef DEBUG_NEIGHBORHOOD
        ++nbIter;
#endif

        // If we are at this end then stop
        if (indexBlock >= instance->getE().size()) throw std::out_of_range(
            "Follower's NeighborHood out of range, indexBlock >= the number of block");
        unsigned int blockSize = static_cast<unsigned int>(instance->getE()[indexBlock].size());
        // Generate the next permutation
        unsigned int sizeArrangement = std::min(nbElementsToArrange, blockSize);
        if (currentPermutation.empty()) currentPermutation.resize(blockSize, 1);
        // Initialize for the current block
        if (not haveBeenPermuted) {
            // if the sizeArrangement is smaller than the blockSize, then compute real arrangement
            if (sizeArrangement != blockSize) std::fill_n(currentPermutation.begin(), blockSize - sizeArrangement, 0);
            else std::iota(currentPermutation.begin(), currentPermutation.end(), 0);
            // fill with index the permutation to compute a real permutation
        }
        do {
            bool hasNext = std::next_permutation(currentPermutation.begin(), currentPermutation.end());
            if (!hasNext) {
                // no more permutation
                haveBeenPermuted = false;
                indexBlock++; // pass to the next block
                if (indexBlock < instance->getE().size()) {
                    blockSize = static_cast<unsigned int>(instance->getE()[indexBlock].size());
                    // Generate the next permutation
                    sizeArrangement = std::min(nbElementsToArrange, blockSize);
                    currentPermutation.resize(blockSize, 1);
                    // if the sizeArrangement is different from the blockSize, then compute real arrangement i.e. sizeArrangement choose among blockSize
                    if (sizeArrangement != blockSize) std::fill_n(currentPermutation.begin(), sizeArrangement, 0);
                    else std::iota(currentPermutation.begin(), currentPermutation.end(), 0);
                    // fill with index in order to compute a real permutation of index

                    bool thereIsAJob = false;
                    // check if there is at least one job in the block
                    for (auto& [indexMachine,indexInMachine] : instance->getE()[indexBlock]) {
                        if ((*baseBlockStructure)[indexMachine][indexInMachine].first != nullptr) {
                            thereIsAJob = true;
                        }
                    }
                    // if there is no jobs, then stop the iterator
                    if (not thereIsAJob) {
                        indexBlock = instance->getE().size(); // set indexBlock to max value to be equal to end()
                        break;
                    }
                    std::next_permutation(currentPermutation.begin(), currentPermutation.end()); // apply at least one permutation
                }else break; // we have no more block
            }
            haveBeenPermuted = false;
            // loop over currentPermutation, if we have a permutation, that mean sizeArrangement == blockSize, then the permutation is valid
            if (sizeArrangement == blockSize) haveBeenPermuted = true;
            else {
                // we have an arrangement, i.e. we need to select sizeArrangement machines among block size, check if
                // the arrangement is valid, i.e. if there is at least one job on each machines
                for (unsigned int indexLoopPermutation = 0; indexLoopPermutation < currentPermutation.size();
                     indexLoopPermutation++) {
                    // by construction, we have currentPermutation.size() == blockSize
                    assert(indexLoopPermutation < blockSize);
                    auto& [indexMachine,indexInMachine] = instance->getE()[indexBlock][indexLoopPermutation];
                    if (currentPermutation[indexLoopPermutation] > 0 && (*baseBlockStructure)[indexMachine][
                        indexInMachine].first != nullptr) {
                        haveBeenPermuted = true;
                        break;
                    }
                }
            }
        }while (not haveBeenPermuted);
    }

    /***************
     *  Operators
     ***************/

    void swapJobs(unsigned int indexMachine, unsigned int indexMachineToSwap,
                  Solution::BlockStructure& newBlockStructure) const {
        assert(instance->getE()[indexBlock].size() > indexMachine);
        unsigned int originalPosition = instance->getE()[indexBlock][indexMachine].second;
        assert(instance->getE()[indexBlock].size() > indexMachineToSwap);
        unsigned int swappedPosition = instance->getE()[indexBlock][indexMachineToSwap].second;
        std::swap(newBlockStructure[indexMachine][originalPosition].first, newBlockStructure[indexMachineToSwap][swappedPosition].first);
        updateCompletionTimeOnMachineAndBlock(originalPosition, indexMachine, newBlockStructure);
        updateCompletionTimeOnMachineAndBlock(swappedPosition, indexMachineToSwap, newBlockStructure);
    }

    Solution::BlockStructure operator*() {
        if (indexBlock >= instance->getE().size()) throw std::out_of_range(
            "Follower's NeighborHood out of range, indexBlock >= the number of block");
        Solution::BlockStructure newBlockStructure = *baseBlockStructure;
        bool isArrangement = static_cast<unsigned int>(instance->getE()[indexBlock].size()) != nbElementsToArrange;
        allTranspositions.clear();
        // apply the permutation
        //find the first index
        unsigned int indexLoopPermutation = 0;
        int firstIndex = -1;
        for (; indexLoopPermutation < currentPermutation.size(); indexLoopPermutation++) {
            if (not isArrangement || currentPermutation[indexLoopPermutation] > 0) {
                if (firstIndex == -1) { firstIndex = static_cast<int>(indexLoopPermutation); }
                else {
                    allTranspositions.emplace_back(firstIndex, indexLoopPermutation);
                    firstIndex = static_cast<int>(indexLoopPermutation);
                }
            }
        }
        for (auto& [indexMachine,indexMachineToSwap] : allTranspositions) swapJobs(indexMachine, indexMachineToSwap, newBlockStructure);

        return newBlockStructure;
    }

    FollowerNeighborhoodIterator& operator++() {
        nextNeighbor();
        return *this;
    }

    FollowerNeighborhoodIterator operator++(int) {
        FollowerNeighborhoodIterator tmp = *this;
        nextNeighbor();
        return tmp;
    }

    /***************
     *  Setters
     ***************/

    void setHaveBeenPermuted(bool haveBeenPermuted) { this->haveBeenPermuted = haveBeenPermuted; }

    bool operator==(const FollowerNeighborhoodIterator& other) const {
        return other.haveBeenPermuted == haveBeenPermuted
                && equalsBase(other);
    }

    bool operator!=(const FollowerNeighborhoodIterator& other) const { return !(*this == other); }

private:
    unsigned int nbElementsToArrange = 2;
    // the nb of element to arrange, if greater than the nb of machine then we compute all permutation
    std::vector<unsigned int> currentPermutation;
    bool haveBeenPermuted = false;
    std::vector<std::pair<unsigned int, unsigned int>> allTranspositions;
    // the list of transposition, i.e. swap, that represent a permutation
};

class Neighborhoods {
public:
    Neighborhoods(Instance* const instance, const Solution::BlockStructure* blockStructure,
                  const std::vector<unsigned int>* listJobsToSwapInSolution,
                  unsigned int nbElementsToArrange) : blockStructure(blockStructure),
                                                      listJobsToSwapInSolution(listJobsToSwapInSolution),
                                                      instance(instance),
                                                      nbElementsToArrange(nbElementsToArrange) {}

    [[nodiscard]] LeaderNeighborhoodIterator beginLN() const {
        LeaderNeighborhoodIterator it(instance, blockStructure, listJobsToSwapInSolution);
        it.updateFirstIndexBlock();
        return it;
    }

    [[nodiscard]] LeaderNeighborhoodIterator endLN() const {
        LeaderNeighborhoodIterator it(instance, blockStructure, listJobsToSwapInSolution);
        it.setIndexJobToSwapFromList(listJobsToSwapInSolution->size());
        return it;
    }

    [[nodiscard]] FollowerNeighborhoodIterator beginFN() const {
        return FollowerNeighborhoodIterator(instance, blockStructure, nbElementsToArrange);
    }

    [[nodiscard]] FollowerNeighborhoodIterator endFN() const {
        FollowerNeighborhoodIterator it(instance, blockStructure);
        it.setIndexBlock(instance->getE().size());
        it.setHaveBeenPermuted(false);
        return it;
    }

#ifdef DEBUG_NEIGHBORHOOD

    unsigned int LBNbIterLeaderNeighbor() {
        unsigned int nbMaxIter = 0;
        for (unsigned int indexJob : (*listJobsToSwapInSolution)) {
            auto& jobToSwap = instance->getListJobs()[indexJob];
            // found the block of the job where it can be schedule
            unsigned int indexBlock = 0;
            double minPjBlock = std::numeric_limits<double>::infinity();
            double maxPjBlock = 0;
            Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance, *blockStructure, indexBlock, maxPjBlock, minPjBlock);
            // pass to the block
            while (not(isSmallerOrEqual(maxPjBlock, jobToSwap.getPi()) && isSmallerOrEqual(jobToSwap.getPi(), minPjBlock)) && indexBlock < instance->getE().size()) {
                ++indexBlock;
                minPjBlock = std::numeric_limits<double>::infinity();
                Solution::computeMaxMinProcessingTimeOnNearbyBlock(instance, *blockStructure, indexBlock, maxPjBlock, minPjBlock);
            }
            nbMaxIter += indexBlock == 0 ? instance->getNbJobsToScheduleOnFirstBlock() : instance->getE()[indexBlock].size();
        }
        return nbMaxIter;
    }
#endif

private:
    const Solution::BlockStructure* blockStructure = nullptr;
    const std::vector<unsigned int>* listJobsToSwapInSolution;
    Instance* instance = nullptr;
    unsigned int nbElementsToArrange = 2;
};

#endif //BILEVEL_SCHEDULING_NEIGHBORHOOD_H
