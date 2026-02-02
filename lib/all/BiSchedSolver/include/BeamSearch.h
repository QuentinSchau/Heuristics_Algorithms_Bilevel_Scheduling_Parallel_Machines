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

#ifndef BILEVEL_SCHEDULING_BEAMSEARCH_H
#define BILEVEL_SCHEDULING_BEAMSEARCH_H

#include "ISolver.h"
#include "ColumnGeneration.h"
#include "LocalSearch.h"
#include "Neighborhood.h"
#include <unordered_set>

/**
 * Hack to get access to the container of priority queue, found in
 * https://stackoverflow.com/questions/1185252/is-there-a-way-to-access-the-underlying-container-of-stl-container-adaptors
 */
template <class T, class S, class C>
    S& Container(std::priority_queue<T, S, C>& q) {
    struct HackedQueue : private std::priority_queue<T, S, C> {
        static S& Container(std::priority_queue<T, S, C>& q) {
            return q.*&HackedQueue::c;
        }
    };
    return HackedQueue::Container(q);
}

class BeamSearch : public ISolver {

public:
    /**
     * Enum to define the recovery strategy.
     */
    enum RecoStrategy { BEST_INSERT , LOCAL_SEARCH };

    struct BeamSearchNode {
        double lowerbound{0};
        double upperbound{std::numeric_limits<double>::infinity()};
        double evaluation{std::numeric_limits<double>::infinity()};
        Node node;
        Solution::BlockStructure completingSchedule;

        inline friend bool operator<(BeamSearchNode const &lhs, BeamSearchNode const &rhs) {
            return lhs.evaluation > rhs.evaluation;
        }

        inline friend std::ostream &operator<<(std::ostream &os, BeamSearchNode const &e) {
            Solution solOfNode;
            solOfNode.fromBlockStruct(e.node.getBlockStruc());
            return os << '{' << e.lowerbound << "," << e.upperbound << "," << e.evaluation << ", " << solOfNode << "'}";
        }
    };

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit BeamSearch();

    explicit BeamSearch(Instance *instance);

    explicit BeamSearch(Instance *instance, nlohmann::json &object);


    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~BeamSearch() override;

    /*******************/
    /*      GETTER     */
    /*******************/

    unsigned int getBeamSize() const { return beamSize;}

    std::string getRecoStrategy() const {
        switch (recoStrategy) {
            case BEST_INSERT:
                return "best-insert";
            case LOCAL_SEARCH:
                return "localSearch";
            default:
                return "Not Implemented";
        }
    }

    [[nodiscard]] double getGlobalUb() const { return globalUB; }

    [[nodiscard]] double getFirstUb() const { return firstUB; }

    /*******************/
    /*      SETTER     */
    /*******************/

    void setBeamSize(unsigned int newBeamSize) {beamSize = newBeamSize;}

    void setUseRecovering(bool newUseRecovering) {useRecovering = newUseRecovering;}

    void setNbBestSolutionKeep(unsigned int nbBestSolutionKeep) { this->nbBestSolutionKeep = nbBestSolutionKeep; }

    void setAlpha(double newAlpha) {
        if (isSmallerOrEqual(0.0,alpha) || isSmallerOrEqual(alpha,1.0)) alpha = newAlpha;
        else {
            std::string errorMessage("The alpha value in BeamSearch heuristic must be between 0.0 and 1.0. Current value is: ");
            errorMessage.append(std::to_string(alpha));
            throw BiSchException(errorMessage);
        }
    }

    void setRecoStrategy(std::string strat){
        if (strat == "best-insert") recoStrategy = BEST_INSERT;
        else if (strat == "local-search") recoStrategy = LOCAL_SEARCH;
        else throw BiSchException("Strategy is not known for the beam search recovering, read \"README\" file for more details on which strategy to use.");
    }

    void setVersion(char version) { this->version = version; }

    /********************/
    /*      METHODS     */
    /********************/

    BeamSearchNode evaluation(Node && node);

    /**
     * Method that checks the condition given by the number of rule 'rule'. If this condition is true,
     * it means we make the best insertion regarding this rule in the recovering method.
     *
     * @param node The node used in the recovering.
     * @param indexBlock The index of the block being considered.
     * @param indexJob The index of the removed jobs that we try to insert.
     * @param blockStruct The block structure of the partial solution.
     * @param listRemovedJobs The list of removed jobs from where we try to insert the job given by 'indexJob'.
     * @param rule The number index of the rule, cf. the paper to identify which rule it's used.
     *
     * @return True if there was a swap, false otherwise.
     */
    bool checkRules(Node &node, unsigned int indexBlock,unsigned int indexJob,const Solution::BlockStructure &blockStruct, std::vector<Job> &listRemovedJobs, char rule);

    void recoveringBestInsert(Node & node);

    // Recovering local search use LB evaluation, LB in local search

    bool findBestSolutionReachableFromPartialSolution(Solution::BlockStructure &blockStructure, double & bestObjValue);

    /**
     * Use Local Search to find a solution with UB < LB, i.e. get real domination
     * @param BSnode The BeamSearch Node from which we apply the recovering
     */
    void recoveringLocalSearchV1(BeamSearchNode & BSnode);

    /**
     * Use Local Search to find a solution with a better UB than the one given by the current node.
     * @param BSnode The BeamSearch Node from which we apply the recovering
     */
    void recoveringLocalSearchV2(BeamSearchNode & BSnode);

    /**
     * Use Local Search to find a solution with a better UB than the one given by the current node.
     * @param BSnode The BeamSearch Node from which we apply the recovering
     */
    void recoveringLocalSearchV3(BeamSearchNode & BSnode);

    void recovering(BeamSearchNode & BSnode);

    #ifdef DEBUG_BaB
    /**
     * Method that check if the state of the node are correct after recovering. Compiling with debug flag
     * @param node The node to check.
     */
    void check_recovering(Node & node);
    #endif

    /***************************************/
    /*      Sub-Problem Identical Jobs     */
    /***************************************/

    /**
     * Method that creates children nodes.
     * @param pNode The parent node.
     * @return True if a node have been added False otherwise
     */
    void createChildrenNodeWithIdenticalJobs(Node &pNode);

    /**
     * Method that computes the number of jobs on last block if we need to schedule 'numberJobsToSchedule' of identical jobs.
     * @param pNode The parent node
     * @param numberJobsToSchedule The number of jobs that we must scheduled
     * @return (listAvailableMachines,k',n_L) the list of available machines from where we can assign n_L jobs on last last block k'
     */
    std::pair<unsigned int, unsigned int> computeNumberJobsOnLastBlock(Node &pNode, unsigned int numberJobsToSchedule);

    /**
     * Computes a set combination without symmetry.
     *
     * @param nbSelect The number of elements to select from the set.
     * @param indexLastBlock The index of the last block.
     * @param node A reference to a Node object.
     * @param assigmentOnStartingBlock A vector of unsigned integers representing the assignment on the starting block.
     * @param setOfAssignmentOnLastBlock A set of combination, represented as vectors of machine indexes, representing the set of assignments on the last block.
     */
    void computeSetCombinationWithoutSymmetry(unsigned int nbSelect, unsigned int indexLastBlock, Node &node, const std::vector<unsigned int> &assigmentOnStartingBlock
                                              , std::set<std::vector<unsigned int>> &setOfAssignmentOnLastBlock);

    /**
    * Method that creates a child node from a given parent node `pNode` with identical jobs using an assignment.
    * @param pNode The parent node from which we create child nodes.
    * @param indexLastBlock The index of the last block where scheduled jobs are located.
    * @param assignmentLastBlock A vector of machine indices where jobs must be scheduled on the last block.
    * @param assignmentFirstBlock A vector of machine indices where jobs must be scheduled on the first block. This vector is not empty when the first block is the leftmost block (block 0).
    *
    */
    void createNodeWithAssignment(Node &pNode, unsigned int indexLastBlock, const std::vector<unsigned int> &assignmentLastBlock, const std::vector<unsigned int> &assignmentFirstBlock);

    /********************/
    /*      General     */
    /********************/

    /**
     * Method that add the given solution to the right priority queue, regarding 'isLeafSolution'.
     * @param solution The solution to add to 'listBestFoundSolutions' or 'listBestFoundLeaf'
     * @param isLeafSolution Boolean, true add the solution to 'listBestFoundLeaf', false add the solution to 'listBestFoundSolutions'.
     */
    void addSolToListBestSolutions(Solution *solution,bool isLeafSolution);

    /**
     * Method that add the solution corresponding to the block structure to the list of founded solutions if and only if it's better upper bound.
     * @param blockStructure The block structure that can be added.
     * @param objValue The objective value of the block structure.
     */
    void addSolToListBestSolutions(Solution::BlockStructure & blockStructure,double objValue);

    /**
     * Method that adds the node to the right structure depending on the walk strategy.
     * @param node The node to add.
     */
    void addNode(Node &node);

    /**
     * Method that generates new nodes by branching from a given node on the all types of machine. The branching scheme
     * use the location. The new nodes are then added to the active list.
     * @param nodeWithLb The parent node from where all children nodes are created with its lower bound.
     */
    void branchingLocation(BeamSearch::BeamSearchNode &beamSearchNode);

    /**
     * Checks if a block has been filled. If so, it updates all relevant constants.
     * @param node The current node in the tree of branch and bound.
     */
    void changeBlock(Node &node);

    /**
     * Method that initialize the tree branch and bound algorithm
     */
    void initialize();

    /**
     * Method that solves the model by using a branch and bound algorithm
     */
    void solve() override;

    /**
     * Method that saves the result of the instance in a file
     * @param fileOutputName The name of the file
     * @param outputFile The stream of the file
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);

private:

    std::vector<Job> listAvailableJobForNode; // use global list of available jobs for a node that we fill whenever we want to add this node
    std::unordered_set<unsigned int> jobsInNewBlockStructure; // use global set to have the job in the new block structure in the recovering assignment

    unsigned int beamSize = 1;
    // The heap to store active node
    std::priority_queue<BeamSearchNode, std::vector<BeamSearchNode>> heap;
    // column generation model to solve the relaxation of a node
    ColumnGeneration columnGeneration;
    Heuristic heuristicSolver;
    LocalSearch localSearchSolver;
    // best know solution
    double globalUB;
    double firstUB = -1; // set the first UB equals to -1, we change it when we compute for the first time an upper bound
    double alpha = 0.0; //use in evaluation

    // local search at the end
    unsigned int nbBestSolutionKeep = 0;
    std::chrono::duration<double> time_elapsed_in_LS; //time elapsed in local search
    //define our own comparator
    template<class L=std::less<>>
    struct Compare {
        bool operator()(std::pair<double,Solution> const& lhs,
                        std::pair<double,Solution> const& rhs) const {
            return isSmallerOrEqual(lhs.first, rhs.first);
        }
    };

    // Use the comparator explicitly here:
    std::priority_queue<
        std::pair<double,Solution>,
        std::vector<std::pair<double,Solution>>,
        Compare<>
    > listBestFoundSolutions;

    // Use the comparator explicitly here:
    std::priority_queue<
        std::pair<double,Solution>,
        std::vector<std::pair<double,Solution>>,
        Compare<>
    > listBestFoundLeaf;

    // recovering
    bool useRecovering = true;
    RecoStrategy recoStrategy = LOCAL_SEARCH;
    char version = 1;
    // reco v1
    unsigned int maxIterLocalSearch = 10;
    // metrics
    unsigned long nbNodeLoc; // number of nodes that was created by the method
    unsigned int NB_reco = 0;
    unsigned int NB_best_sol_find_reco = 0;
    unsigned int NB_MULTI_START_SOLUTIONS_FROM_LEAF = 0; // the number of solution corresponding at a leaf for the multi start local search
    unsigned int NB_MULTI_START_SOLUTIONS_FROM_BEST_SOL = 0; // the number of solution corresponding at a best solution for the multi start local search
};

#ifndef BILEVEL_SCHEDULING_BEAMSEARCH_IMP_H
// DO NOT INCLUDE IN OTHER FILE ! IT'S IMPLEMENTATION OF BeamSearch.h METHODS INLINED
#include "BeamSearchImp.hpp"

#endif

#endif //BILEVEL_SCHEDULING_BEAMSEARCH_H
