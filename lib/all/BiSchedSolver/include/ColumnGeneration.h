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

#ifndef BILEVEL_SCHEDULING_COLUMNGENERATION_H
#define BILEVEL_SCHEDULING_COLUMNGENERATION_H

// uncomment this line to debug the method
//#define DEBUG_CG
// uncomment this line to use a MIP to solve the pricing problem the method
//#define DEBUG_MIP_PRICING

#include <ranges>
#include "ISolver.h"
#include "Heuristic.h"
#include "ilconcert/ilomodel.h"
#include "ilcplex/ilocplex.h"
#include "Node.h"
#include "HungarianAlgorithm.h"
#include "DSU.h"
#include "CGForIdenticalJobsPricing.h"

class ColumnGeneration : public ISolver {
private:
    #ifdef DEBUG_CG
    unsigned int nbGen =0;
    #endif
    double lowerBound; // Lower bound of the objective function

    /*      CPLEX      */

    IloEnv env; // Environment for CPLEX solver
    IloModel model; // Model for optimization problem
    IloCplex cplex; // CPLEX solver instance
    IloNumVar obj; // Objective function variable
    IloRange objConstr; // Constraint for objective function
    IloRange UComputeConstr; // Constraint for U (weighted tardy jobs)
    IloNumVar U; // Value of weighted tardy jobs
    IloNumVarArray xs; // Array of variables in the model
    IloRangeArray constraints; // Array of constraints in the model
    IloNumArray duals; // Dual values for constraints
    std::vector<double> dualsValues;

    /*      Dynamic Programming      */

    // Memorization data structure to store already explored states
    std::vector<double> memo;
    std::vector<bool> memoExplored;
    std::vector<unsigned int> l_max_HS; // Maximum identical jobs on high-speed machines
    std::vector<unsigned int> l_max_LS; // Maximum identical jobs on low-speed machines
    unsigned int maxLValue = 0;

    /**
     * Hash function to get the value of lmax (l_max_HS or l_max_LS).
     * @param k Number of jobs already selected on the machine.
     * @param g Index of the group of identical jobs.
     * @return Lmax value.
     */
    size_t hashLmax(unsigned int k, unsigned int g) {
        return g + instance->getListGrpJobs().size() * k;
    }

    /**
     * Hash function to get the index in memorization data structure.
     * @param k Number of jobs already selected on the machine.
     * @param t Completion time.
     * @param j Job index.
     * @return Index in memorization data structure.
     */
    size_t hashMemo(unsigned int k, unsigned int t, unsigned int j) {
        size_t hash = j + (instance->getListGrpJobs().size() + 1) * k + (instance->getListGrpJobs().size() + 1) * (instance->getMaxNbJobsOnHS() + 1) * t;
        assert(hash < memo.size());
        assert(hash < memoExplored.size());
        return hash;
    }

    // Sums of max processing time that can be scheduled on high-speed and low-speed machines
    unsigned int sum_n0p_high_pj = 0; // High-speed machines
    unsigned int sum_n1p_high_pj = 0; // Low-speed machines

    // Vector to store constant values for identical jobs algorithm
    std::vector<ConstForIdenticalJobs> list_const_identical_job;

    // DSU (Disjoint Set Union) structure with number of identical jobs
    DSU *dsu = nullptr;
    std::vector<JumpPoint> listJumpPoint;

    /*      Heuristic      */

    Heuristic heuristicSolver;

    // Vector to keep track of already explored jobs in estimation function
    std::vector<bool> alreadyExplored;

    // Maximum number of heuristic calls before stopping
    unsigned int maxNbCallHeuristic = 1; // Default value

    // Constants for removing unused columns

    // Threshold for cleaning policy (ratio of unfeasible columns)
    double thresholdSetCol = 0.8;
    // Number of times a column is not used before being considered unfeasible
    unsigned int nbTimeNotUsed = 20;

    /*      General      */

    // Debug flag
    bool debug;

    // Flag to define which method generates columns (DP or heuristic)
    char generate_Column; // '0' for DP, '1' for Heuristic

    // Flag to indicate if master problem solving failed
    bool failedSolveMasterProblem;

    // Number of minimum reduced costs computed
    unsigned int nbMinStateDP = 5;

    // Counters for various operations
    unsigned int nbCallsDP = 0;
    unsigned int nbCallsHeu = 0;
    unsigned long nbCallComputeCost = 0;
    unsigned long nbCallSubProcessCG = 0;
    unsigned int nbCleaningSetCol = 0;

    // Enum for column types (high/low speed, max/min)
    enum TYPE_COLUMN { HS_MAX, HS_MIN, LS_MAX, LS_MIN };

    // Vector to store new set of columns generated by backtracking
    std::vector<std::pair<TYPE_COLUMN, std::vector<unsigned int>>> newSetOfColumns; // Column type and corresponding values

    // Lists to store states for backtracking (forward computation)
    std::vector<std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>> listStatesDPForBacktracking;


    struct InfoColumn {
        std::vector<unsigned int> column; // The column vector representing a machine schedule

        boost::dynamic_bitset<> encodingColumn; // Bitset representation of the column for efficient storage and comparison

        unsigned int indexOfVariable = 0; // Index of the corresponding variable in the model
        unsigned int nbUnfeasible = 0; // Number of times this column has been deemed unfeasible

        double costOfSchedule = 0.0; // Associated cost of the schedule (two values for \sum wj U_J and \alpha * \sum Cj)

        TYPE_COLUMN typeOfColumn; // Type of column (high/low speed, max/min)

        double UB_var; // Upper bound value for this variable

        inline InfoColumn(std::vector<unsigned int> &&col, unsigned int indexOfVariable, unsigned int nbUnfeasible, double costOfSchedule, TYPE_COLUMN typeColumn, unsigned int nbJobs) :
                column(std::move(col)), // Move the vector of column into the struct
                encodingColumn(nbJobs), // Initialize bitset with default values (all set to false)
                indexOfVariable(indexOfVariable), // Assign index of variable
                nbUnfeasible(nbUnfeasible), // Set number of unfeasibilities
                costOfSchedule(costOfSchedule), // Move the pair of costs into the struct
                typeOfColumn(typeColumn), // Set type of column
                UB_var(1.0) { // Initialize upper bound value
            for (unsigned int &indexJobInSchedule: column) { // Iterate over the vector of column indices
                encodingColumn.set(indexJobInSchedule, true); // Set corresponding bit in the bitset to true
            }
        }
    };

    struct StateBacktracking {
        unsigned int t; // Completion time
        unsigned int alreadySelected; // Number of jobs already selected on the machine
        double reducedCost; // Reduced cost of this state
        unsigned int indexMachine; // Index of the machine in question
        char typeMachineSchedule; // Type of schedule for this machine (0: min jobs, 1: max jobs, 2: unknown)
        /**
         * Constructor for StateBacktracking.
         * @param alreadySelected Number of jobs already selected on the machine.
         * @param t Completion time.
         * @param reducedCost Reduced cost of this state.
         * @param indexMachine Index of the machine in question.
         * @param typeMachineSchedule Type of schedule for this machine (0: min jobs, 1: max jobs, 2: unknown).
         */
        inline StateBacktracking(unsigned int alreadySelected, unsigned int t, double reducedCost, unsigned int indexMachine, char typeMachineSchedule) :
                t(t), // Initialize Completion time
                alreadySelected(alreadySelected), // Initialize number of selected jobs
                reducedCost(reducedCost), // Initialize reduced cost
                indexMachine(indexMachine), // Initialize machine index
                typeMachineSchedule(typeMachineSchedule) {} // Initialize schedule type
    };


    // The set of columns is represented using a vector. For each entry in the vector, we store:
    //  - The machine schedule (column)
    //  - Its bitset encoding for efficient storage and comparison (encodingColumn)
    //  - The index of the corresponding variable in the model (indexOfVariable)
    //  - The number of times this column has been deemed unfeasible (nbUnfeasible)
    //  - The associated cost of the schedule (costOfSchedule)
    //  - The type of column (high/low speed machine with max/min jobs) (typeOfColumn)
    //  - The upper bound value for this variable (UB_var)
    std::vector<InfoColumn> setMachineSchedules;

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    explicit ColumnGeneration();

    explicit ColumnGeneration(Instance *instance);

    explicit ColumnGeneration(Instance *instance, nlohmann::json &object);

    /***********************/
    /*      DESTRUCTOR     */
    /***********************/

    ~ColumnGeneration() override;

    /*******************/
    /*      GETTER     */
    /*******************/


    /**
     * This method returns the lower bound given by the column generation for the sum of wj * Uj.
     * @return The lower bound of the sum of wj * Uj given by the column generation.
    */
    [[nodiscard]] double getSumWjUj() const;

    [[nodiscard]] const IloNumVarArray &getXs() const { return xs; }

    [[nodiscard]] unsigned int getNbCallsDp() const { return nbCallsDP; }

    [[nodiscard]] unsigned int getNbCallsHeu() const { return nbCallsHeu; }

    [[nodiscard]] unsigned int getNbCallSubProcessCG() const { return nbCallSubProcessCG; }

    [[nodiscard]] unsigned int getNbMinStateDp() const { return nbMinStateDP; }

    [[nodiscard]] unsigned long getNbCallComputeCost() const { return nbCallComputeCost; }

    [[nodiscard]] char getGenerateColumn() const { return generate_Column; }

    [[nodiscard]] unsigned int getMaxNbCallHeuristic() const { return maxNbCallHeuristic; }

    [[nodiscard]] unsigned int getNbCleaningSetCol() const { return nbCleaningSetCol; }

    /*******************/
    /*      SETTER     */
    /*******************/

    void setParameters(nlohmann::json &object);

    void setDebug(bool isDebug) { debug = isDebug; }

    void setGenColumns(char methodToGenerateColumn) { generate_Column = methodToGenerateColumn; }

    void setMaxNbCallHeuristic(unsigned int newMaxNbCallHeuristic) { maxNbCallHeuristic = newMaxNbCallHeuristic; }

    void setThresholdSetCol(double newThresholdSetCol) { thresholdSetCol = newThresholdSetCol; }

    void setNbTimeNotUsed(unsigned int newNbTimeNotUsed) { nbTimeNotUsed = newNbTimeNotUsed; }

    void setNbMinStateDp(unsigned int nbMinStateDp) { nbMinStateDP = nbMinStateDp; }



    /********************/
    /*      METHODS     */
    /********************/

    void initialize();

    /**
     * Method that saves the current MIP model to a file in LP format.
     * @param modelPath The path where the model will be saved.
     */
    void saveModel(const std::string &modelPath) {
        cplex.exportModel(modelPath.c_str()); // Export CPLEX model to LP file at specified path
    };

    /**
     * Computes the index of the next job in a list of grouped jobs using a bijection from N^2 to N.
     * The bijection maps (group index, job index) to a unique index.
     * @param listOfJobs A list containing the jobs, sorted by Shortest Processing Time (SPT), and grouped by same processing time.
     * @param indexOfJob The index of the current job in its group.
     * @return The new index of the next job.
     */
    static unsigned int nextIndex(const std::vector<std::vector<Job>> &listOfJobs, unsigned int indexOfJob);

    /**
     * Inverse of the bijection used to compute indices. Given an index computed by the bijection, returns the group
     * index and job index that correspond to it.
     * @param indexOfJob The index computed using the bijection.
     * @return A pair containing the group index and job index corresponding to the given index.
     */
    static std::pair<unsigned int, unsigned int> getGroupAndJobIndex(unsigned int indexOfJob);

    void updateValueOfLmax(Node &node, char machineSpeed);


    /********************/
    /*      Columns     */
    /********************/

    /**
     * Calculates the cost associated with a given machine scheduling.
     * @param node The node in the tree of branch and bound.
     * @param machineSchedule A list of job indices, arranged in Shortest Processing Time (SPT) order, representing the assigned jobs for each time slot on the machine.
     *  @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @return The cost associated with this machine scheduling, consisting of two values: the sum of weighted tardy job costs and alpha times the Total completion time.
    */
    double computeScheduleCost(Node &node, std::vector<unsigned int> &machineSchedule, char machineSpeed);

    /**
     * Calculates the reduced cost associated with a given machine scheduling.
     * @param node The node in the tree of branch and bound.
     * @param machineSchedule A list of job indices, arranged in Shortest Processing Time (SPT) order, representing the assigned jobs for each time slot on the machine.
     *  @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param typeColumn The type of column (high/low speed machine with max/min jobs).
     * @return The reduced cost corresponding to the provided machine scheduling.
    */
    double computeScheduleReducedCost(Node &node, std::vector<unsigned int> &machineSchedule, char machineSpeed, TYPE_COLUMN typeColumn);

    /**
    * Computes the minimum reduced cost for scheduling a job on a machine using dynamic programming. The corresponding schedule can be constructed by backtracking.
    * All jobs in this machine schedule are sorted according to the Shortest Processing Time (SPT) rule.
    * It is required that the dual values are obtained before using this method with "cplex.getDuals(duals, constraints)".
    * @param k The number of jobs that need to be scheduled on the machine schedule..
    * @param typeColumn The type of column (high/low speed machine with max/min jobs).
    * @param t The completion time of the last scheduled job on the machine schedule. Must be divided by the speed to get the real completion time.
    * @param g The index of the current group of jobs that are being scheduled.
    * @param node The node in the tree of branch and bound.
    * @return The minimum reduced cost.
    */
    double computeMinReducedCost(unsigned int k, TYPE_COLUMN typeColumn, unsigned int t, unsigned int g, const Node &node);

    #ifdef DEBUG_MIP_PRICING
    double computeMinReducedCost(Node &node, TYPE_COLUMN typeColumn, unsigned int indexMachine,std::vector<unsigned int> &bestCol);
    #endif

    /**
     * Method that solves the pricing problem.
     * @param node The node in the tree of branch and bound.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listStartingStates A vector to store all starting states.
     * @param indexStartingJobForDP Reference to the index of the starting job in dynamic programming.
     * @return a tuple (m_{g+-},m_{g+},r_{g+},m_{g+},r_{g+}) where m_{g+-} is the number of machine schedule where we don't know if there is max or min of jobs on it.
     * m_{g+} is the number of machine schedules with max of jobs and r_{g+} is the corresponding minimal reduced cost for this kind of machine schedule.
     * Similarly m_{g-} is the number of machine schedules with min of jobs and r_{g-} is the corresponding minimal reduced cost for this kind of machine schedule
     */
    std::tuple<unsigned int,unsigned int,double,unsigned int,double> solvePricingProblem(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int indexStartingJobForDP);

    /**
     * Compute the minimal reduced cost by selecting 'nbJobsToSelected' job among all identical jobs.
     * @param nbJobsToSelected The number of jobs to select among all identical jobs.
     * @param t The completion time where to start.
     * @param g The index of the current group of jobs that are being scheduled.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listOfJobs The list of identical jobs.
     * @return The minimal reduced cost by selecting the given number of jobs among all identical jobs.
     */
    double subRoutine(unsigned int nbJobsToSelected, unsigned int t, unsigned int g, char machineSpeed);

    /**
     * Compute the minimal reduced cost by selecting 'nbJobsToSelected' job among all identical jobs. Adding the best selection to the current 'newColumn' given in parameters
     * @param nbJobsToSelected The number of jobs to select among all identical jobs.
     * @param t The completion time where to start.
     * @param g The index of the current group of jobs that are being scheduled.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listOfJobs The list of identical jobs.
     * @param newColumn The column where the best selection will be added.
     */
    void subRoutine(unsigned int nbJobsToSelected, unsigned int t, unsigned int g, char machineSpeed, std::vector<unsigned int> &newColumn);

    /**
     * Method to compute all starting states for dynamic programming based on a given node, machine speed, and index of starting job.
     *
     * @param node The node in the tree of branch and bound.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listStartingStates A vector to store all starting states.
     * @param indexStartingJobForDP Reference to the index of the starting job in dynamic programming.
     */
    void computeFirstsStateOfDynamicProgramming(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int &indexStartingJobForDP);

    /**
     * Updates the reduced cost of starting states in dynamic programming based on a node, machine speed, and list of starting states.
     * @param node The node in the tree of branch and bound.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listStartingStates A vector to store all starting states.
     */
    void updateRedCostOfStartingStateOfDP(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates);

    /**
     * Method that updates the constants V_j and Q_j before solving the dynamic programming. These constants are used in the sub-procedure of Jippe.
     * @param node The node in the tree of branch and bound.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     */
    void update_M_Vj_Qj_Values(Node &node, char machineSpeed);


    /**
     * Computes the contribution of the job to the reduced cost depending on whether job j is late or not on
     * the machine with the given speed.
     * @param t The completion time of the last scheduled job on the machine schedule.
     * @param jobJ The job j.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @return The contribution of this job to the reduced cost.
     */
    double beta(unsigned int t, const Job &jobJ, char machineSpeed);


    /**
     * Constructs the column associated to the minimum reduced cost using backtracking.
     * Before using this method, computeScheduleReducedCost must have been run, and both memoization maps cannot be empty.
     * @param node The node in the tree of branch and bound.
     * @param listStartingState The list of starting states of the dynamic programming from where we compute the backtracking.
     * @param startingIndexJob The index of the first job that needs to be considered in the dynamic programming.
     *  @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     */
    void
    constructColumnByBacktracking(Node &node, const std::vector<StateBacktracking> &listStartingState, unsigned int startingIndexJob, char machineSpeed);

    /**
     * Creates a column by forward computation.
     * @param k The number of jobs that need to be scheduled on the machine schedule.
     * @param t The completion time of the last scheduled job on the machine schedule. Must be divided by the speed to get the real completion time.
     * @param g The index of the current group of jobs that are being scheduled.
     * @param newColumn The column created by the forward computation.
     * @param keepTrackState Boolean flag indicating whether to keep track of states leading to the optimal value.
     * @param typeColumn The type of column (high/low speed machine with max/min jobs).
     * @param node The node in the tree of branch and bound.
     * @param listStates The list of states that lead to the minimal reduced cost.
     * @return Minimal reduced cost of the column (without the constant part, cf. formula).
    */
    double createColumnByForward(unsigned int k, unsigned int t, unsigned int g, std::vector<unsigned int> &newColumn, bool keepTrackState, TYPE_COLUMN typeColumn, Node &node
                                 , std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &listStates);

    /**
     * This method changes the value in the memoization of the last state from 'listStates'. After changing this value,
     * the method updates by backward, all states in the list 'listStates' of the memorization with a new value.
     * The idea is when we recall the forward method, we will get a new minimal reduced column.
     * @param typeColumn The type of column (high/low speed machine with max/min jobs).
     * @param listStates The list of states that lead to the minimal reduced cost.
     */
    void createColumnByBackward(TYPE_COLUMN typeColumn, std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> &listStates);

    /**
     * Generates new columns with a reduced cost from the values of dual variables.
     * @param node The node in the tree used to generate a set of columns.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param listStartingStates The list of starting states for the dynamic programming.
     * @param indexStartingJobForDP Reference to the index of job from which we need to start the dynamic programming.
     * @param minReducedCost Minimal reduced costs obtained by solving the pricing problem.
     * @return True if there exist column with negative reduced cost, false otherwise.
     */
    bool generateColumns(Node &node, char machineSpeed, std::vector<StateBacktracking> &listStartingStates, unsigned int &indexStartingJobForDP, std::tuple<unsigned int,unsigned int,double,unsigned int,double> &minRedCosts);

    /**
     * Method that clear memorisation.
     */
    void clearMemo() {
        memoExplored.assign((instance->getListGrpJobs().size() + 1) * static_cast<size_t>(sum_n0p_high_pj + 1) * (instance->getMaxNbJobsOnHS() + 1), false);
    };

    /*********************************/
    /*      Column - Heuristic       */
    /*********************************/


    /**
     * Generates new column with a negative reduced cost. The main idea, is to compute a lower bound of the reduced cost
     * in a window and generate the min colum.
     * @param node The node in the tree use to generate a set of columns with the heuristic.
     * @param machineScheduleOnMachine The machine schedule already assign during the branching scheme.
     * @param indexStartingJobForDP The index of job from which we need to start the dynamic programming.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param maxOrMinJobsOnCol Flag to define which kind of column to generate: 0 for min jobs on a high-speed machine or 1 for max jobs on a low-speed machine.
     * @return The most negative reduced cost compute by the heuristic.
     */
    double generateColumnByHeuristic(Node &node, std::vector<unsigned int> &machineScheduleOnMachine, unsigned int indexStartingJobForDP, char machineSpeed, char maxOrMinJobsOnCol);

    /**
     * Method that calculates a lower bound on the reduced cost by utilizing a window of size 'sizeOfWindows'.
     * The idea is to utilize the 'sizeOfWindows' smallest processing times and calculate the minimal contribution for each one.
     * Sum all of them, and you get your lower bound.
     * @param node The node in the tree used to generate a set of columns with the heuristic.
     * @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     * @param sizeOfWindows The number of jobs using for the size of the windows.
     * @param startingPairingIndex The index from where we start in list of jobs. Here we use the bijection index to know the group and which jobs in the group.
     * @param t The completion time from where we start the calculation of the contribution.
     */
    double estimationOfReducedCost(Node &node, char machineSpeed, unsigned int sizeOfWindows, unsigned int startingPairingIndex, unsigned int t);

    /**
     * Method that add column, on machine speed, to the master problem using a given solution.
     * @param solution The solution from where we create columns.
     * @param node The node in the tree use to generate a set of columns.
     *  @param machineSpeed The speed consideration, 0 for high speed and 1 for low speed.
     */
    void addColumnFromSolution(Solution &solution, Node &node, char machineSpeed);

    /**
     * Method that generates an initial set of columns for the master problem using a given node.
     * If a job index is not in provided lists from the node, it will not be added to the generated columns.
     * @param node The node in the tree use to generate a set of columns.
    */
    void generateStartingColumns(Node &node);

    /******************/
    /*      Model     */
    /******************/

    /**
     * Method that initializes the basic MIP formulation.
     * @param node The node in the tree of branch and bound. 
     */
    void initializeModel(Node &node);

    /**
     * Method that parametrizes the solver.
     */
    void parametrize();

    /**
     * Method that removes all constraints of the model.
     */
    void clearConstraintOfModel();

    /**
     * Method that updates variables that have a unfeasible column.
     * @param node The node in the tree of branch and bound. 
     */
    void updateModel(Node &node);

    /**
     * Method that check if a column is feasible or note.
     * @param node The node in the tree of branch and bound.
     * @return True if columns are feasible otherwise False.
     */
    bool checkColumns(Node &node);

    /**
     * Method that create columns in order to make the master restricted problem feasible,
     * with the previous decision from the node.
     * @param node The current node from which we are branching.
     */
    void initRestrictedMasterProblem(Node &node);

    /**
     * Method that solves the model that is active in memory.
     */
    void solve() override;

    /**
     * Method that solves the model using a given node, that holds lists of job indices and jobs.
     * @param node The node that contains all lists of indexes and jobs.
     */
    void solve(Node &node);

    /**
     * Method that saves the result of the instance in a file.
     * @param fileOutputName The name of the file.
     * @param outputFile The stream of the file.
     */
    void printOutput(std::string &fileOutputName, std::ofstream &outputFile);

};

#ifndef BILEVEL_SCHEDULING_COLUMNGENERATION_IMP_H
// DO NOT INCLUDE IN OTHER FILE ! IT'S IMPLEMENTATION OF ColumnGeneration.h METHODS INLINED
#include "ColumnGenerationImp.hpp"

#endif


#endif //BILEVEL_SCHEDULING_COLUMNGENERATION_H
