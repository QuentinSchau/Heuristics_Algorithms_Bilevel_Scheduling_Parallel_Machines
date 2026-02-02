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

#ifndef BILEVEL_SCHEDULING_LOCALSEARCH_H
#define BILEVEL_SCHEDULING_LOCALSEARCH_H


#include "ISolver.h"
#include "Heuristic.h"
#include "Neighborhood.h"
#include <torch/script.h>
#include <torch/torch.h>

class LocalSearch : public ISolver {
private:
    Heuristic heuristicSolver;

    /**
     * Predictor
     */
    // Batch processing structure to hold all necessary data for each prediction
    struct BatchItem {
        std::vector<unsigned int> scheduledJobs;
        std::vector<bool> selectedJobs;
        double weightSum;

        BatchItem(
            std::vector<unsigned int> jobs,
            std::vector<bool> selected,
            double weight
        ) : scheduledJobs(std::move(jobs)),
            selectedJobs(std::move(selected)),
            weightSum(weight) {}
    };
    unsigned int BATCH_SIZE = 128;  // Optimal for CPU inference
    typedef std::tuple<double,std::vector<unsigned int>,std::vector<bool>> predictionWithSetJob;
    //define our own comparator for prediction with set
    template<class L=std::less<>>
    struct ComparePredictionWithSet {
        bool operator()(predictionWithSetJob & left, predictionWithSetJob & right) {
            return isSmaller(std::get<0>(left),std::get<0>(right));
        }
    };
    std::vector<torch::jit::IValue> inputs; // inputs use to compute prediction

    std::string pathPredictor;
    bool usePredictor = false;
    bool genDatabase = false;
    std::array<double,96> features;
    std::vector<std::vector<double>> minMaxAvgCompletionTime;
    std::vector<std::vector<double>> minMaxAvgStdForEachBlock;

    double ratioNeighbor = 0.1;
    char version = 1;
    double firstUB = -1; // set the first UB equals to -1, we change it when we compute for the first time an upper bound
    bool useBestNeighbor = false;
    unsigned int maxIter = 20;
    unsigned int nbIter = 0;
    unsigned int nbPrediction = 10;

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit LocalSearch();

    explicit LocalSearch(Instance *instance);

    explicit LocalSearch(Instance *instance, nlohmann::json &object);

    void initializeStructure();

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~LocalSearch() override;

    void solve() override;

    /**********************/
    /*      Heuristic     */
    /**********************/

    /**
     * Method that compute an initial solution for the local search.
     * @param alreadySelectedJobs A vector of bool to know which jobs that have been already selected.
     * @param listScheduledJobs We use a vector of index of job to compute a solution.
     * @return A pair with the initial solution through a block structure object and its value.
     */
    std::pair<Solution::BlockStructure, double> computeInitialSolution(std::vector<bool> &alreadySelectedJobs,std::vector<unsigned int> &listScheduledJobs);

    /**
     * Methods that insert and remove a job's index. Doing it in O(n) operations.
     *
     * @param listIndexOfJobs The list of indices of jobs, indexed by SPT.
     * @param indexJobToInsert The index of the job that will be inserted into the list.
     * @param indexJobToRemove The index of the job to remove from the list.
     */
    void static swapJobsInSelection(std::vector<unsigned int> &listIndexOfJobs, unsigned int indexJobToInsert, unsigned int indexJobToRemove);

    /**
     * Method that updates the block structure by removing and inserting the given job's index.
     * This method reschedules jobs in the impacted block according to the SPT-FAM rule.
     *
     * @param listIndexOfJobs The list of indices of jobs, indexed by SPT.
     * @param blockStructure The block structure to be altered with the modification.
     * @param indexJobToInsert The index of the job that will be inserted into the list.
     * @param indexJobToRemove The index of the job to remove from the list.
     *
     * @return A pair containing the updated changed block structure and its corresponding evaluation of the weighted number of tardy jobs.
     */
    std::pair<Solution::BlockStructure,double> updateBlockStructureWithSwap(std::vector<unsigned int> &listIndexOfJobs, Solution::BlockStructure &blockStructure, unsigned int indexJobToInsert, unsigned int indexJobToRemove);

    /**
     * Method that tries to find a better swap of 2 jobs for each block. Two configurations:
     * find the best improving swap or the first one found; then move to the next block.
     *
     * This method can also try to solve the assignment of all identical jobs optimally.
     * It updates some data structures used in the local search, such as lists of jobs (which can be swapped) and a vector of booleans indicating which jobs are
     * already scheduled.
     *
     * @param blockStructure The block structure of a given solution that must be improved.
     * @param listScheduledJobs The list of indices of jobs that are selected/scheduled and can be swapped in the solution.
     * @param listNotScheduledJobs The list of indices of jobs that are not selected. (Optional, may be null)
     * @param alreadySelectedJobs The vector of boolean indicating which jobs are already scheduled in the solution.
     * @param bestObjValue The best objective value, updated if a better solution is found.
     * @param solveIdenticalJobs Boolean to know whether to try to reassign identical jobs optimally. Defaults to false.
     *
     * @return True if a better solution was found; False otherwise.
     */
    bool improveBlockStructureLocallyBySwap(Solution::BlockStructure &blockStructure,std::vector<unsigned int> &listScheduledJobs,std::vector<unsigned int> * listNotScheduledJobs, std::vector<bool> &alreadySelectedJobs,double &bestObjValue, bool solveIdenticalJobs=false);

    /**
     * Method that tries to find a better assignment for each block. The method makes two passes:
     * 1) From left to right (first block to last block).
     * 2) If an improvement is found, it tries to do the same in reverse order: from right to left.
     * This process continues until the maximum number of iterations is reached or no further improvements are made.
     *
     * This method can also try to solve the assignment of all identical jobs optimally.
     * It updates some data structures used in the local search, such as lists of jobs (which can be swapped) and a vector of booleans indicating which jobs are
     * already scheduled.
     *
     * @param blockStructure The block structure of a given solution that must be improved.
     * @param listScheduledJobs The list of indices of jobs that are selected/scheduled and can be swapped in the solution.
     * @param alreadySelectedJobs The vector of boolean indicating which jobs are already scheduled in the solution.
     * @param bestObjValue The best objective value, updated if a better solution is found.
     * @param solveIdenticalJobs Boolean to know whether to try to reassign identical jobs optimally. Defaults to false.
     *
     * @return True if a better solution was found; False otherwise.
     */
    bool improveBlockStructureLocallyByAssignment(Solution::BlockStructure &blockStructure,std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs,double &bestObjValue, bool solveIdenticalJobs=false);

    /**
     * Method that explores the neighborhood of a solution.
     *
     * @param bestBlockStructure The block structure of a given partial solution that must be improved.
     * @param listScheduledJobs The list of indices of jobs that are selected/scheduled and can be swapped in the solution.
     * @param alreadySelectedJobs The vector of boolean indicating which jobs are already scheduled in the solution.
     * @param bestObjValue The best objective value, updated if a better solution is found.
     */
    void exploreNeighborhood(Solution::BlockStructure & bestBlockStructure,std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs, double &bestObjValue);

    /**
     * Local search that starts from an initial solution,
     * the method makes a maximum iteration 'maxIter', for each iteration we have:
     * 1) improve the best solution locally, i.e. only the schedule, by considering swap for version = 1 or assignment for version = 5
     * 2) try to swap a non-selected job with a selected one, then rebuild a feasible schedule by keeping as much as possible the original schedule given by best solution,
     *    then try to improve locally by considering swap for version = 1 or assignment for version = 5.
     *
     * @param initialBlockStructure Pointer to an initial block structure, i.e. an initial solution, that we try to improve.
     *                               By default, there is no initial solution, then we compute one.
     * @param initialSolutionValue The value of the given initial solution.
     *                             This value must be equal to the value of the given initial solution.
     */
    void localSearchBySwapAllJobs(Solution::BlockStructure* initialBlockStructure = nullptr, double initialSolutionValue = 0.0);

    void localSearchOnlySwapV1();

    void localSearchOnlySwapV2();

    // Version Only Assigment:

    void localSearchOnlyAssigment();

    // Version Predictor:
    void localSearchPredictor();

    void exploreRandomNeighborhoodWithPredictorFullTrust(Solution::BlockStructure & bestBlockStructure, std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs, double &bestObjValue);

    void exploreNeighborhoodWithPredictorFullTrust(Solution::BlockStructure & bestBlockStructure, std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs, double &bestObjValue);

    void exploreNeighborhoodWithPredictorNotFullTrust(Solution::BlockStructure & bestBlockStructure, std::vector<unsigned int> &listScheduledJobs, std::vector<bool> &alreadySelectedJobs, double &bestObjValue);

    bool computePredictions(torch::jit::Module& model,torch::Tensor& batchTensor, std::vector<std::vector<float>>& batchFeatures,std::vector<BatchItem> &batchItems,std::priority_queue<predictionWithSetJob,std::vector<predictionWithSetJob>,ComparePredictionWithSet<>> &listBestPrediction);

    // Version Ramdom
    /**
     * Method that create take an random selection and try to improve the block structure
     */
    void localSearchRandom();

    /*********************/
    /*      Features     */
    /*********************/

    /**
     * Compute features from a block structure.
     * @param blockStructure The block structure used to compute features.
     */
    void computeFeatures(Solution::BlockStructure &blockStructure);

    /**
     * Method that computes features, i.e., min, max, average, and standard deviation from the list of completion times computed with 'computeFeatures'.
     * @param blockStructure The block structure used to compute features.
     */
    void computeFeaturesFromCompletionTime(Solution::BlockStructure &blockStructure);

    /**
     * Computes the min, max, average, and standard deviation of a list of values.
     * @param listValues The list of values to compute features from.
     * @return An array containing the minimum, maximum, average, and standard deviation of the input list.
     */
    std::array<double,4> computeMinMaxAvgStdForListValues(std::vector<double> &listValues);

    /********************/
    /*      SETTER      */
    /********************/

    void setParameters(nlohmann::json &object);

    void setUsePredictor(bool newUsePredictor) {this->usePredictor = newUsePredictor;}

    void setGenDatabase(bool newGen_data_base) {genDatabase = newGen_data_base;}

    void setUseBestNeighbor(bool newUseBestNeighbor){this->useBestNeighbor = newUseBestNeighbor;}

    void setVersion(char newVersion) { this->version = newVersion; }

    void setMaxIter(unsigned int newMaxIter) {this->maxIter = newMaxIter;}

    void setNbPrediction(unsigned int nbPrediction) { this->nbPrediction = nbPrediction; }

    void setBatchSize(unsigned int batchSize) {this->BATCH_SIZE = batchSize;}

    void setRatioNeighbor(double ratioNeighbor) { this->ratioNeighbor = ratioNeighbor; }

    void setTimeUpHeuristicSolver(std::shared_ptr<std::atomic<bool>> &newTimeUp){heuristicSolver.setTimeUp(newTimeUp);}
    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] bool isUsePredictor() const {return usePredictor;}

    bool isGenDatabase() const {return genDatabase;}

    [[nodiscard]] char getVersion() const { return version; }
    [[nodiscard]] bool isUseBestNeighbor() const { return useBestNeighbor; }

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);
};

void inline LocalSearch::computeFeatures(Solution::BlockStructure &blockStructure) {
    auto &E = instance->getE();
    // compute the inverse of speed to only compute multiplication instead of division
    const double invHS = 1.0 / instance->getHighSpeed();
    const double invLS = 1.0 / instance->getLowSpeed();
    const unsigned int nbHS = instance->getNbOfHighSpeedMachines();

    // clear the structure uses to compute min, max and average completion time
    for (auto &stat : minMaxAvgCompletionTime) std::fill(stat.begin(), stat.end(), 0.0);

    // compute the min, max and avg for each block
    for (unsigned int indexBlock = 0; indexBlock < E.size(); indexBlock++) {
        auto &block = E[indexBlock];

        // Détermination du type une seule fois
        int typeBlock = 2;
        if (block.front().first < nbHS && block.back().first < nbHS) typeBlock = 0;
        else if (block.front().first >= nbHS && block.back().first >= nbHS) typeBlock = 1;

        double minPjOnHS = 1e18, maxPjOnHS = 0, avgPjOnHS = 0;
        double minPjOnLS = 1e18, maxPjOnLS = 0, avgPjOnLS = 0;

        for (auto &[mIdx, jIdx] : block) {
            auto job = blockStructure[mIdx][jIdx].first;
            if (!job) continue;

            double pi = job->getPi();
            if (typeBlock == 0 || typeBlock == 2) {
                double p = pi * invHS;
                if (p < minPjOnHS) minPjOnHS = p;
                if (p > maxPjOnHS) maxPjOnHS = p;
                avgPjOnHS += p;
            }
            if (typeBlock == 1 || typeBlock == 2) {
                double p = pi * invLS;
                if (p < minPjOnLS) minPjOnLS = p;
                if (p > maxPjOnLS) maxPjOnLS = p;
                avgPjOnLS += p;
            }
        }
        // set the min of pj that was not explored to 0.0
        if (minPjOnHS > instance->getMaxPj()) minPjOnHS = 0.0;
        if (minPjOnLS > instance->getMaxPj()) minPjOnLS = 0.0;

        avgPjOnHS /= static_cast<double>(indexBlock == 0 ? instance->getNbJobsToScheduleOnFirstBlock() : block.size());
        avgPjOnLS /= static_cast<double>(indexBlock == 0 ? instance->getNbJobsToScheduleOnFirstBlock() : block.size());
        minMaxAvgCompletionTime[0][indexBlock] = indexBlock == 0 ? minPjOnHS : minMaxAvgCompletionTime[0][indexBlock-1]+minPjOnHS;
        minMaxAvgCompletionTime[1][indexBlock] = indexBlock == 0 ? maxPjOnHS : minMaxAvgCompletionTime[1][indexBlock-1]+maxPjOnHS;
        minMaxAvgCompletionTime[2][indexBlock] = indexBlock == 0 ? avgPjOnHS : minMaxAvgCompletionTime[2][indexBlock-1]+avgPjOnHS;
        minMaxAvgCompletionTime[3][indexBlock] = indexBlock == 0 ? minPjOnLS : minMaxAvgCompletionTime[3][indexBlock-1]+minPjOnLS;
        minMaxAvgCompletionTime[4][indexBlock] = indexBlock == 0 ? maxPjOnLS : minMaxAvgCompletionTime[4][indexBlock-1]+maxPjOnLS;
        minMaxAvgCompletionTime[5][indexBlock] = indexBlock == 0 ? avgPjOnLS : minMaxAvgCompletionTime[5][indexBlock-1]+avgPjOnLS;
    }
    computeFeaturesFromCompletionTime(blockStructure);
    //normalize the vector
    std::transform(features.begin(), features.end(), features.begin(), [&](double &feature){return feature / instance->getSumWj();});
}

void inline LocalSearch::computeFeaturesFromCompletionTime(Solution::BlockStructure &blockStructure){
    auto &E = instance->getE();
    const double sumWjLimit = instance->getSumWj();

    // Process each of the 6 completion time variants
    for (unsigned int k = 0; k < 6; ++k) {
        // Clear working structure for current metric
        for (auto &stat : minMaxAvgStdForEachBlock) {
            std::fill(stat.begin(), stat.end(), 0.0);
        }

        for (unsigned int indexBlock = 0; indexBlock < E.size(); indexBlock++) {
            auto &block = E[indexBlock];
            double minValue = 1e18, maxValue = 0.0, sumValue = 0.0, squareSumValue = 0.0;
            double currentThreshold = minMaxAvgCompletionTime[k][indexBlock];

            for (auto &[mIdx, jIdx] : block) {
                auto job = blockStructure[mIdx][jIdx].first;
                // Condition: Job is considered delayed based on the estimated completion time
                if (job && isSmaller(job->getDi(), currentThreshold)) {
                    double wi = job->getWi();
                    if (wi < minValue) minValue = wi;
                    if (wi > maxValue) maxValue = wi;
                    sumValue += wi;
                    squareSumValue += wi * wi;
                }
            }

            double nbElements = static_cast<double>(indexBlock == 0 ? instance->getNbJobsToScheduleOnFirstBlock() : block.size());
            double avgValue = sumValue / nbElements;
            // Use std::max to prevent negative values due to floating point precision errors before sqrt
            double stdValue = std::sqrt(std::max(0.0, (squareSumValue / nbElements) - (avgValue * avgValue)));

            if (minValue > sumWjLimit) minValue = 0.0;

            // Compute cumulative values for min, max, and avg
            minMaxAvgStdForEachBlock[0][indexBlock] = (indexBlock == 0) ? minValue : minMaxAvgStdForEachBlock[0][indexBlock - 1] + minValue;
            minMaxAvgStdForEachBlock[1][indexBlock] = (indexBlock == 0) ? maxValue : minMaxAvgStdForEachBlock[1][indexBlock - 1] + maxValue;
            minMaxAvgStdForEachBlock[2][indexBlock] = (indexBlock == 0) ? avgValue : minMaxAvgStdForEachBlock[2][indexBlock - 1] + avgValue;
            minMaxAvgStdForEachBlock[3][indexBlock] = stdValue;
        }

        // Map the results back to the global features vector
        unsigned int featIdx = k * 16;
        for (unsigned int s = 0; s < 4; ++s) {
            auto stats = computeMinMaxAvgStdForListValues(minMaxAvgStdForEachBlock[s]);
            features[featIdx++] = stats[0];
            features[featIdx++] = stats[1];
            features[featIdx++] = stats[2];
            features[featIdx++] = stats[3];
        }
    }
}

std::array<double,4> inline LocalSearch::computeMinMaxAvgStdForListValues(std::vector<double> &listValues) {
    double minValue = std::numeric_limits<double>::max(), maxValue = 0.0, sumValue = 0.0, squareSumValue = 0.0;

    for (auto value : listValues){
        minValue = std::min(minValue, value);
        maxValue = std::max(maxValue, value);
        sumValue += value;
        squareSumValue += value*value;
    }
    double nbElement = static_cast<double>(listValues.size());
    double avgValue = sumValue / nbElement;
    double stdValue = std::sqrt(squareSumValue/nbElement - avgValue*avgValue);
    return {minValue, maxValue, avgValue, stdValue};
}

#endif //BILEVEL_SCHEDULING_LOCALSEARCH_H
