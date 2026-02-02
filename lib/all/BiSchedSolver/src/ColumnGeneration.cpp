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
// Created by schau on 3/15/24.
//


#include "ColumnGeneration.h"

ColumnGeneration::ColumnGeneration() :
        lowerBound(0.0), model(env), cplex(model), obj(env, "obj"), U(env, "U"), xs(env),constraints(env), duals(env), heuristicSolver(instance), debug(false), generate_Column(0), failedSolveMasterProblem(false)
        , nbCallsDP(0), nbCallsHeu(0) {
    model.add(IloMinimize(env, obj));
    parametrize();
}

ColumnGeneration::ColumnGeneration(Instance *instance) :
        ISolver(instance), lowerBound(0.0), model(env), cplex(model), obj(env, "obj"), U(env, "U"), xs(env),constraints(env), duals(env), heuristicSolver(instance), debug(false)
        , generate_Column(0), failedSolveMasterProblem(false), nbCallsDP(0), nbCallsHeu(0) {
    initialize();
}

ColumnGeneration::ColumnGeneration(Instance *instance, nlohmann::json &object) :
        ISolver(instance, object), lowerBound(0.0), model(env), cplex(model), obj(env, "obj"), U(env, "U"), xs(env), constraints(env), duals(env),heuristicSolver(instance), debug(false)
        , generate_Column(0), failedSolveMasterProblem(false), nbCallsDP(0), nbCallsHeu(0) {
    setParameters(object);
    initialize();
}

ColumnGeneration::~ColumnGeneration() {
    env.end();
    delete dsu;
}

/*******************/
/*      SETTER     */
/*******************/

void ColumnGeneration::setParameters(nlohmann::json &object) {
    if (object.contains("name")) {
        if (object["name"] == "CG") {
            // check if we have verbose mode
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of CG solver)");
            }

            if (object.contains("debug")) {
                if (object["debug"].is_boolean()) setDebug(object["debug"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "debug" of JSON object must be an "boolean" value for the constructor of CG solver)");
            }
            if (object.contains("gen_columns")) {
                if (object["gen_columns"].is_number_integer()) setGenColumns(object["gen_columns"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "gen_columns" of JSON object must be an "positive integer" value for the constructor of CG solver)");
            } else
                throw std::invalid_argument(
                        R"(The JSON object have not the attribute method to define which method is use to generate new columns)");
            if (object.contains("maxNbCallHeuristic")) {
                if (object["maxNbCallHeuristic"].is_number_integer())setMaxNbCallHeuristic(object["maxNbCallHeuristic"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "maxNbCallHeuristic" of JSON object must be an "positive integer" value for the constructor of CG solver)");
            }
            if (object.contains("nbStateDP")) {
                if (object["nbStateDP"].is_number_integer())setMaxNbCallHeuristic(object["nbStateDP"].template get<char>());
                else
                    throw std::invalid_argument(
                            R"(The attribute "maxNbCallHeuristic" of JSON object must be an "positive integer" value for the constructor of CG solver)");
            }
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of CG solver)");
            }if (object.contains("timeLimitsMS")) {
                if (object["timeLimitsMS"].is_number_unsigned()) setTimeLimitInMilliSecond(object["timeLimitsMS"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "timeLimitsMS" of JSON object must be an "unsigned int" value for the constructor of CG solver)");
            }
            if (object.contains("nbMinStateDP")) {
                if (object["nbMinStateDP"].is_number_unsigned()) setNbMinStateDp(object["nbMinStateDP"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "nbMinStateDP" of JSON object must be an "unsigned int" value for the constructor of CG solver)");
            }
            if (object.contains("thresholdSetCol")) {
                if (object["thresholdSetCol"].is_number_float()) setThresholdSetCol(object["thresholdSetCol"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "thresholdSetCol" of JSON object must be an "float" value for the constructor of CG solver)");
            }
            if (object.contains("nbTimeNotUsed")) {
                if (object["nbTimeNotUsed"].is_number_unsigned()) setNbTimeNotUsed(object["nbTimeNotUsed"]);
                else
                    throw std::invalid_argument(
                            R"(The attribute "nbTimeNotUsed" of JSON object must be an "unsigned int" value for the constructor of CG solver)");
            }
        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a CG solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a CG solver)");
}

/********************/
/*      METHODS     */
/********************/

void ColumnGeneration::initialize() {
    model.add(IloMinimize(env, obj));
    objConstr = IloRange(env, 0.0, 0.0, "Objective");
    model.add(objConstr);
    UComputeConstr = IloRange(env, 0.0, 0.0, "Compute_U_Value");
    model.add(UComputeConstr);
    unsigned int indexConstraint = 0;
    for (; indexConstraint < instance->getNbJobs(); indexConstraint++) {
        constraints.add(IloRange(env, 0.0, 1.0, std::string("constraint_job_J").append(std::to_string(indexConstraint)).c_str()));
    }
    constraints.add(IloRange(env, 0.0, instance->getNbOfHighSpeedMachines(), "constraint_select_V1_machines"));
    constraints.add(IloRange(env, 0.0, instance->getNbOfLowSpeedMachines(), "constraint_select_V0_machines"));
    constraints.add(IloRange(env, instance->getNbToSelectJob(), instance->getNbToSelectJob(), "constraint_select_n_jobs"));
    model.add(constraints);
    parametrize();

    // create the matrix of cost for identical jobs
    //find the max number of identical jobs
    size_t maxNumberIdenticalJobs = 0;
    for (auto &groupIdenticalJobs: instance->getListGrpJobs()) {
        maxNumberIdenticalJobs = std::max(maxNumberIdenticalJobs, groupIdenticalJobs.size());
    }
    auto listJob = instance->getListJobs();
    std::sort(listJob.begin(), listJob.end(), std::greater<>());
    unsigned int i = 0;
    for (; i < instance->getMaxNbJobsOnLS(); i++) {
        sum_n1p_high_pj += static_cast<unsigned int>(listJob[i].getPi());
    }
    sum_n0p_high_pj = sum_n1p_high_pj;
    for (; i < instance->getMaxNbJobsOnHS(); i++) {
        sum_n0p_high_pj += static_cast<unsigned int>(listJob[i].getPi());
    }
    //initialize structure for identical jobs
    dsu = new DSU(maxNumberIdenticalJobs + 1);
    listJumpPoint.assign(maxNumberIdenticalJobs + 1, JumpPoint(0, 0, 0));
    // create enough list of constants for identical jobs with the maximal identical jobs that we can found
    for (unsigned int indexLoopGroupIdJobs = 0; indexLoopGroupIdJobs <= instance->getNbJobs(); ++indexLoopGroupIdJobs)
        list_const_identical_job.emplace_back(maxNumberIdenticalJobs, dsu, &listJumpPoint);

    dualsValues = std::vector<double>(instance->getNbJobs() + 4);

    //reserve size for l_max
    l_max_HS.resize((instance->getMaxNbJobsOnHS() + 1) * instance->getListGrpJobs().size(), 0);
    l_max_LS.resize((instance->getMaxNbJobsOnLS() + 1) * instance->getListGrpJobs().size(), 0);

    // reserve the size for memoization
    memo.resize((instance->getListGrpJobs().size() + 1) * static_cast<size_t>(sum_n0p_high_pj + 1) * (instance->getMaxNbJobsOnHS() + 1), std::numeric_limits<double>::infinity());

    newSetOfColumns.reserve(nbMinStateDP);
    listStatesDPForBacktracking.resize(instance->getNbMachines(), std::vector<std::tuple<unsigned int, unsigned int, unsigned int>>(instance->getNbJobs()));
    alreadyExplored.resize(instance->getNbJobs(), false);
}

double
ColumnGeneration::computeScheduleReducedCost(Node &node, std::vector<unsigned int> &machineSchedule, char machineSpeed, TYPE_COLUMN typeColumn) {
    auto cost = computeScheduleCost(node, machineSchedule, machineSpeed);
    double machineScheduleReducedCost = cost;
    // the dual value of type machine constraint
    double typeMachineCoef = 0.0;
    bool firstBlockIsOnBothTypeMachine = instance->isFirstBlockIsOnBothTypeMachine();
    switch (typeColumn) {
        case TYPE_COLUMN::HS_MAX:typeMachineCoef = dualsValues[instance->getNbJobs()];
            //add the second multiplier from the other constraint (present twice)
            if (instance->isFirstBlockIsOnBothTypeMachine()) {
                typeMachineCoef += dualsValues[instance->getNbJobs() + 1];
                typeMachineCoef += dualsValues[instance->getNbJobs() + 2];
            }
            break;
        case TYPE_COLUMN::HS_MIN:typeMachineCoef = dualsValues[instance->getNbJobs() + 1];
            break;
        case TYPE_COLUMN::LS_MAX:typeMachineCoef = dualsValues[instance->getNbJobs() + 2];
            if (firstBlockIsOnBothTypeMachine) typeMachineCoef += dualsValues[instance->getNbJobs() + 3];
            break;
        case TYPE_COLUMN::LS_MIN:typeMachineCoef = dualsValues[instance->getNbJobs() + 3];
            break;
    }

    // compute the cost
    for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
        double ajs = std::find(machineSchedule.begin(), machineSchedule.end(), j) != machineSchedule.end();
        machineScheduleReducedCost -= ajs * dualsValues[j];
    }
    machineScheduleReducedCost -= typeMachineCoef;
    return machineScheduleReducedCost;
}


/*
 * NOT INLINED
 */
double
ColumnGeneration::computeMinReducedCost(unsigned int k, TYPE_COLUMN typeColumn, unsigned int t, unsigned int g, const Node &node) {
    auto &listOfJobs = instance->getListGrpJobs();
    double minReducedCost = std::numeric_limits<double>::infinity();
    char machineSpeed = (typeColumn == TYPE_COLUMN::LS_MAX || typeColumn == TYPE_COLUMN::LS_MIN) ? 1 : 0;
    unsigned int sum_n_high_pj = (machineSpeed == 0) ? sum_n0p_high_pj : sum_n1p_high_pj;
    auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
    // TERMINAL STATES
    // if k = max count or min count
    if (k == 0) {
        memo[hashMemo(k, t, g)] = 0.0;
        memoExplored[hashMemo(k, t, g)] = true;
        return 0.0;
    } else {
        // j > N
        if (g >= instance->getListGrpJobs().size()) {
            memo[hashMemo(k, t, g)] = minReducedCost;
            memoExplored[hashMemo(k, t, g)] = true;
            return minReducedCost;
        }
    }
    if (t > sum_n_high_pj) {
        memo[hashMemo(k, t, g)] = minReducedCost;
        memoExplored[hashMemo(k, t, g)] = true;
        return minReducedCost;
    }
        // if we have already solved this state
    else if (memoExplored[hashMemo(k, t, g)]) {
        minReducedCost = memo[hashMemo(k, t, g)];
        return minReducedCost;
    } else {
        // loop over all value possible value of k
        assert(hashLmax(k, g) < lmax.size());
        unsigned int maxValueL = lmax[hashLmax(k, g)];
        auto passNextJob = computeMinReducedCost(k, typeColumn, t, g + 1, node);
        double selectLJobs = std::numeric_limits<double>::infinity();
        if (not node.isGroupIdenticalJobRemoved(g)) {
            unsigned int p = static_cast<unsigned int>(listOfJobs[g].back().getPi());
            // if we have one jobs in Gg
            if (listOfJobs[g].size() == 1 && not node.isScheduledOnOtherMachines(listOfJobs[g].back().getIndex(), machineSpeed)) {
                unsigned int newT = t + p;
                // if the job j respects the release date and deadline
                // select and assign the job
                auto assignJobs = computeMinReducedCost(k - 1, typeColumn, newT, g + 1, node);
                selectLJobs = std::min(selectLJobs, assignJobs + beta(t, listOfJobs[g].back(), machineSpeed));

            } else {
                // we have several jobs in Gg, and try to select several (with the maximum value l)
                for (unsigned int l = 1; l <= maxValueL; l++) {
                    if (l > k) break;
                    unsigned int newT = t + l * p;
                    double contributionSelection = subRoutine(l, t, g, machineSpeed);
                    auto assignJobs = computeMinReducedCost(k - l, typeColumn, newT, g + 1, node);
                    selectLJobs = std::min(selectLJobs, assignJobs + contributionSelection);
                }
            }
        }
        minReducedCost = std::min(passNextJob, selectLJobs);
        memo[hashMemo(k, t, g)] = minReducedCost;
        memoExplored[hashMemo(k, t, g)] = true;
        return minReducedCost;
    }
}

#ifdef DEBUG_MIP_PRICING
double ColumnGeneration::computeMinReducedCost(Node &node, TYPE_COLUMN typeColumn, unsigned int indexMachine,std::vector<unsigned int> &bestCol){
    IloEnv envPricing = IloEnv();
    IloModel modelPricing = IloModel(envPricing);
    IloCplex cplexPricing = IloCplex(modelPricing);
    IloNumVar objPricing = IloNumVar(envPricing,-IloInfinity, IloInfinity, ILOFLOAT,"obj");
    std::vector<IloNumVarArray> e_jk;
    std::vector<IloNumVarArray> t_jk;
    IloNumVarArray c_k = IloNumVarArray(envPricing);
    IloRangeArray constraintsPricing = IloRangeArray(envPricing);

    // create data
    unsigned int maxPj = static_cast<unsigned int>(instance->getMaxPj());
    unsigned int maxJobOnMachine = 0;
    char machineSpeed = 0;
    double constReducedCost=0.0;
    switch (typeColumn) {
        case TYPE_COLUMN::HS_MAX:
            maxJobOnMachine = instance->getMaxNbJobsOnHS();
            constReducedCost = -dualsValues[instance->getNbJobs()];
            if (instance->isFirstBlockIsOnBothTypeMachine()) {
                constReducedCost -= dualsValues[instance->getNbJobs() + 1];
                constReducedCost -= dualsValues[instance->getNbJobs() + 2];
            }
            break;
        case TYPE_COLUMN::HS_MIN:
            maxJobOnMachine = instance->getMinNbJobsOnHS();
            constReducedCost = -dualsValues[instance->getNbJobs() + 1];
            break;
        case TYPE_COLUMN::LS_MAX:
            maxJobOnMachine = instance->getMaxNbJobsOnLS();
            machineSpeed = 1;
            constReducedCost = -dualsValues[instance->getNbJobs() + 2];
            if (instance->isFirstBlockIsOnBothTypeMachine()) constReducedCost -= dualsValues[instance->getNbJobs() + 3];
            break;
        case TYPE_COLUMN::LS_MIN:
            maxJobOnMachine = instance->getMinNbJobsOnLS();
            machineSpeed = 1;
            constReducedCost = -dualsValues[instance->getNbJobs() + 3];
            break;
    }
    double speed = machineSpeed == 0 ? instance->getHighSpeed() : instance->getLowSpeed();
    double Cmax = (machineSpeed == 0) ? sum_n0p_high_pj : sum_n1p_high_pj;
    auto &lmax = (machineSpeed == 0) ? l_max_HS : l_max_LS;
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        c_k.add(IloNumVar(envPricing, 0.0, IloInfinity, ILOFLOAT, std::string("c").append(std::to_string(k)).c_str()));
    }

    //create variables
    for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
        e_jk.emplace_back(envPricing);
        t_jk.emplace_back(envPricing);
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            e_jk[j].add(IloNumVar(envPricing, 0.0, 1.0, ILOBOOL, std::string("e(").append(std::to_string(j)).append(",").append(std::to_string(k)).append(")").c_str()));
            t_jk[j].add(IloNumVar(envPricing, 0.0, 1.0, ILOBOOL, std::string("t(").append(std::to_string(j)).append(",").append(std::to_string(k)).append(")").c_str()));
        }
    }

    // set variable regarding the current node
    auto &machine = node.getBlockStruc()[indexMachine];
    //check if the machine is not empty, i.e. the first or second position are filled
    if (machine[0].first != nullptr || machine[1].first != nullptr) {
        unsigned int position_k = 0;
        for (auto jobInSchedule: machine) {
            // if we have a job we set the corresponding variables
            if (jobInSchedule.first != nullptr) {
                auto job = jobInSchedule.first;
                //set the completion time
                c_k[position_k].setBounds(jobInSchedule.second, jobInSchedule.second);
                //set variables if is late or not
                if (isSmallerOrEqual(jobInSchedule.second, job->getDi())) {
                    e_jk[job->getIndex()][position_k].setBounds(1.0, 1.0);
                    t_jk[job->getIndex()][position_k].setBounds(0.0, 0.0);
                }else{
                    t_jk[job->getIndex()][position_k].setBounds(1.0, 1.0);
                    e_jk[job->getIndex()][position_k].setBounds(0.0, 0.0);
                }
                ++position_k;
            }
        }
    }
    // if job is removed or scheduled on other machines, then make it unavailable
    for (auto &job: instance->getListJobs()) {
        if (node.isRemoved(job.getIndex())) {
            for (unsigned int k = 0; k < maxJobOnMachine; k++) {
                t_jk[job.getIndex()][k].setBounds(0.0, 0.0);
                e_jk[job.getIndex()][k].setBounds(0.0, 0.0);
            }
            continue;
        }
        for (unsigned int indexOtherMachine = 0; indexOtherMachine < instance->getNbMachines(); indexOtherMachine++) {
            if (indexOtherMachine == indexMachine) continue;
            if (node.getEncodingSelectedJobOnMachine()[indexOtherMachine].test(job.getIndex())) {
                for (unsigned int k = 0; k < maxJobOnMachine; k++) {
                    t_jk[job.getIndex()][k].setBounds(0.0, 0.0);
                    e_jk[job.getIndex()][k].setBounds(0.0, 0.0);
                }
                break;
            }
        }
    }

    // define constraint

    // Declaring the objective function
    // the product (t_jk+e_jk)*c_k must be linearised
    std::vector<IloNumVarArray> l_jk;
    //create variables
    for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
        l_jk.emplace_back(envPricing);
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            l_jk[j].add(IloNumVar(envPricing, 0.0,Cmax, ILOFLOAT, std::string("l(").append(std::to_string(j)).append(",").append(std::to_string(k)).append(")").c_str()));
        }
    }
    for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            IloExpr ExprL1(envPricing);
            ExprL1 = l_jk[j][k] - (e_jk[j][k]+t_jk[j][k])*Cmax;
            auto name = "linUB_BIGCmaxJ" + std::to_string(j) + "K" + std::to_string(k);
            constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprL1,0.0,name.c_str()));
            ExprL1.end();

            IloExpr ExprL2(envPricing);
            ExprL2 = l_jk[j][k] - c_k[k];
            name = "linUBCmaxJ" + std::to_string(j) + "K" + std::to_string(k);
            constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprL2,0.0,name.c_str()));
            ExprL2.end();

            IloExpr ExprL3(envPricing);
            name = "linLB_BIGCmaxJ" + std::to_string(j) + "K" + std::to_string(k);
            ExprL3 = l_jk[j][k] - c_k[k] + (1 - (e_jk[j][k]+t_jk[j][k]))*Cmax;
            constraintsPricing.add(IloRange(envPricing,0.0,ExprL3,IloInfinity,name.c_str()));
            ExprL3.end();

            IloExpr ExprL4(envPricing);
            ExprL4 = l_jk[j][k];
            name = "linLB_ZERO_J" + std::to_string(j) + "K" + std::to_string(k);
            constraintsPricing.add(IloRange(envPricing,0.0,ExprL4,IloInfinity,name.c_str()));
            ExprL4.end();
        }
    }

    modelPricing.add(IloMinimize(envPricing, objPricing));
    IloExpr ExprObj(envPricing);
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
            ExprObj += ((t_jk[j][k]+e_jk[j][k])*(-dualsValues[j]) + instance->getListJobs()[j].getWi()*t_jk[j][k]);
        }
    }
    ExprObj -= objPricing;
    modelPricing.add(ExprObj == 0);
    ExprObj.end();

    // One job by position
    for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
        IloExpr ExprC1(envPricing);
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            ExprC1 += (t_jk[j][k] + e_jk[j][k]);
        }
        ExprC1 -= 1;
        auto name = "OneJobJ" + std::to_string(j) + "ByPos";
        constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprC1,0.0,name.c_str()));
        ExprC1.end();
    }

    // One position by job
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        IloExpr ExprC2(envPricing);
        for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
            ExprC2 += (t_jk[j][k] + e_jk[j][k]);
        }
        ExprC2 -= 1;
        auto name = "OnePosK" + std::to_string(k) + "ByJob";
        constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprC2,0.0,name.c_str()));
        ExprC2.end();
    }

    // compute completion time
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        IloExpr ExprC3(envPricing);
        // If were are not at first batch
        if (k > 0) {
            ExprC3 = c_k[k-1];
        }
        for (unsigned int j = 0; j < instance->getNbJobs(); ++j) {
            auto pj = instance->getListJobs()[j].getPi();
            ExprC3 += (t_jk[j][k] + e_jk[j][k]) * pj / speed;
        }
        ExprC3 -= c_k[k];
        auto name = "ComputeC" + std::to_string(k);
        constraintsPricing.add(IloRange(envPricing,0.0,ExprC3,0.0,name.c_str()));
        ExprC3.end();
    }

    // compute tardy jobs
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        IloExpr ExprC4(envPricing);
        for (unsigned int j = 0; j < instance->getNbJobs(); ++j)
            ExprC4 += t_jk[j][k]* instance->getListJobs()[j].getDi();
        ExprC4 -= c_k[k];
        auto name = "ComputeTardyLBC" +std::to_string(k);
        constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprC4,0.0,name.c_str()));
        ExprC4.end();
        double HV = double(k+1) * double(maxPj) / speed;
        IloExpr ExprC5(envPricing);
        for (unsigned int j = 0; j < instance->getNbJobs(); ++j)
            ExprC5 -= (e_jk[j][k]* instance->getListJobs()[j].getDi() + t_jk[j][k]* HV);
        ExprC5 += c_k[k];
        name = "ComputeTardyUBC" + std::to_string(k);
        constraintsPricing.add(IloRange(envPricing,-IloInfinity,ExprC5,0.0,name.c_str()));
        ExprC5.end();
    }

    // select the right nb of jobs
    IloExpr ExprC6(envPricing);
    for (unsigned int k = 0; k < maxJobOnMachine; k++) {
        for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
            ExprC6 += (t_jk[j][k] + e_jk[j][k]);
        }
    }
    ExprC6 -= maxJobOnMachine;
    constraintsPricing.add(IloRange(envPricing,0.0,ExprC6,0.0,std::string("SelectionJobs").c_str()));
    ExprC6.end();

    // make SPT order on the machine
    if (maxJobOnMachine>0) {
        for (unsigned int k = 0; k < maxJobOnMachine - 1; k++) {
            IloExpr ExprC7(envPricing);
            for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
                auto pj = instance->getListJobs()[j].getPi();
                ExprC7 += (t_jk[j][k] + e_jk[j][k]) * pj - (t_jk[j][k + 1] + e_jk[j][k + 1]) * pj;
            }
            auto name = "SPTPosK" + std::to_string(k);
            constraintsPricing.add(IloRange(envPricing, -IloInfinity, ExprC7, 0.0, name.c_str()));
            ExprC7.end();
        }
        // set lmax job for identical ones
        for (unsigned int k = 0; k < maxJobOnMachine - 1; k++) {
            for (unsigned int g = 0; g < instance->getListGrpJobs().size(); g++) {
                if (instance->getListGrpJobs()[g].size() == 1) continue;
                IloExpr ExprC8(envPricing);
                // lmax with the number of jobs to schedule given by maxJobOnMachine - k
                assert(hashLmax(maxJobOnMachine-k, g) < lmax.size());
                unsigned int maxIdJobCanBeSchedule = lmax[hashLmax(maxJobOnMachine-k, g)];
                unsigned int lastPosition = std::min(k + maxIdJobCanBeSchedule + 1, maxJobOnMachine);//+1 because we loop until < lastPosition
                for (unsigned int nextK = k; nextK < lastPosition; nextK++) {
                    for (auto &job: instance->getListGrpJobs()[g]) {
                        ExprC8 += (t_jk[job.getIndex()][nextK] + e_jk[job.getIndex()][nextK]);
                    }
                }

                ExprC8 -= maxIdJobCanBeSchedule; // we can have at most maxIdJobCanBeSchedule
                auto name = "LmaxPos" + std::to_string(k) + "grp" + std::to_string(g);
                constraintsPricing.add(IloRange(envPricing, -IloInfinity, ExprC8, 0.0, name.c_str()));
                ExprC8.end();
            }
        }
    }

    modelPricing.add(constraintsPricing);
    cplexPricing.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
    cplexPricing.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0);
    cplexPricing.setParam(IloCplex::Param::Emphasis::Numerical, true);
    cplexPricing.setOut(env.getNullStream());
    if (!cplexPricing.solve()) {
        Solution::printB(node.getBlockStruc());
        std::string modelName = instance->getInstancePath().parent_path().string();
        modelName.append("/model_debug/model_unfeasible_Pricing.lp");
        std::filesystem::path fileSolutionPath = std::filesystem::path(modelName);
        std::filesystem::create_directories(fileSolutionPath.lexically_normal().parent_path());
        cplexPricing.exportModel(modelName.c_str());
        IloNumArray prefs(env, constraintsPricing.getSize());
        for (unsigned int i = 0; i < constraintsPricing.getSize(); i++) {
            prefs[i] = 1.0;
        }
        if (cplexPricing.refineConflict(constraintsPricing,prefs)) {
            modelName = instance->getInstancePath().parent_path().string();
            modelName.append("/model_debug/model_unfeasible_Pricing.clp");
            cplexPricing.writeConflict(modelName.c_str());
        } else {
            std::cerr << "Échec de l'affinage du conflit." << std::endl;
        }
        throw BiSchException("Error solve MIP Pricing");
    }
    else {
        std::vector<std::pair<unsigned int,unsigned int>> selectedVar;
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
                double x_jk = (cplexPricing.getValue(t_jk[j][k]) + cplexPricing.getValue(e_jk[j][k]));
                // check if x_j,k equals 1.0 with 10^-6 precision
                if (isEqual(1.0,x_jk, 1e-5)) {
                    selectedVar.emplace_back(k, j);
                }
            }
        }
        if (debug && verbose > 3) {
            for (unsigned int k = 0; k < maxJobOnMachine; k++) {
                std::cout << "C(" << k << "]=" << cplexPricing.getValue(c_k[k]) << std::endl;
                for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
                    std::cout << e_jk[j][k].getLB() << " <= e(" << j << "," << k << ")=" << cplexPricing.getValue(e_jk[j][k]) << " <=" << e_jk[j][k].getUB() << std::endl;
                    std::cout << t_jk[j][k].getLB() << " <= t(" << j << "," << k << ")=" << cplexPricing.getValue(t_jk[j][k]) << " <=" << t_jk[j][k].getUB() << std::endl;
                    std::cout << l_jk[j][k].getLB() << " <= l(" << j << "," << k << ")=" << cplexPricing.getValue(l_jk[j][k]) << " <=" << l_jk[j][k].getUB() << std::endl;
                }
            }
        }
        std::sort(selectedVar.begin(), selectedVar.end());
        for (auto &[k, j]: selectedVar) {
            bestCol.emplace_back(j);
        }
    }

    double minRedCost = cplexPricing.getObjValue();
    minRedCost += constReducedCost;
    double redCostOfCol = computeScheduleReducedCost(node, bestCol, machineSpeed, typeColumn);
    if (std::abs(redCostOfCol - minRedCost) > EPSILON) {
        std::cout << "Solver status: " << cplex.getStatus() << std::endl;
        double gap = cplex.getMIPRelativeGap();
        std::cout << "Gap: " << gap << "\n";
        for (unsigned int k = 0; k < maxJobOnMachine; k++) {
            std::cout << "C(" << k << "]=" << cplexPricing.getValue(c_k[k]) << std::endl;
            for (unsigned int j = 0; j < instance->getNbJobs(); j++) {
                std::cout << e_jk[j][k].getLB() << " <= e(" << j << "," << k << ")=" << cplexPricing.getValue(e_jk[j][k]) << " <=" << e_jk[j][k].getUB() << std::endl;
                std::cout << t_jk[j][k].getLB() << " <= t(" << j << "," << k << ")=" << cplexPricing.getValue(t_jk[j][k]) << " <=" << t_jk[j][k].getUB() << std::endl;
                std::cout << l_jk[j][k].getLB() << " <= l(" << j << "," << k << ")=" << cplexPricing.getValue(l_jk[j][k]) << " <=" << l_jk[j][k].getUB() << std::endl;
            }
        }
        throw BiSchException(("compute reduced cost failed : diff= " + std::to_string(std::abs(redCostOfCol - minRedCost))).c_str());
    }
    envPricing.end();
    return minRedCost;
}
#endif

/******************/
/*      Model     */
/******************/

void ColumnGeneration::parametrize() {
    cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0.0);
    cplex.setParam(IloCplex::Param::Simplex::Tolerances::Optimality, 1e-9);
    cplex.setParam(IloCplex::Param::Simplex::Tolerances::Feasibility, 1e-9);
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::Param::TimeLimit, time_limits.count());
    cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Dual);
    cplex.setParam(IloCplex::Param::NodeAlgorithm, IloCplex::Algorithm::Dual);
    if (!debug) cplex.setParam(IloCplex::Param::MIP::Display, 0);
    if (verbose <= 3) cplex.setOut(env.getNullStream());
}

void ColumnGeneration::solve() {
    instance->sort_by_SPT();
    // create a root node of branch and bound
    Node node = Node(instance);
    generateStartingColumns(node);
    clearConstraintOfModel();
    initializeModel(node);
    solve(node);
}

void ColumnGeneration::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists)
        outputFile << "InstanceName"
                   << "\t" << "InstancePath"
                   << "\t" << "N"
                   << "\t" << "n"
                   << "\t" << "m_Max"
                   << "\t" << "m_0"
                   << "\t" << "V_max"
                   << "\t" << "V_0"
                   << "\t" << "Method"
                   << "\t" << "Time"
                   << "\t" << "LimitTime"
                   << "\t" << "nbCallsDP"
                   << "\t" << "nbCallsHeu"
                   << "\t" << "Objective" << std::endl;
    // write value
    outputFile << instance->getInstanceName()
               << "\t" << instance->getInstancePath().string()
               << "\t" << instance->getNbJobs()
               << "\t" << instance->getNbToSelectJob()
               << "\t" << instance->getNbOfHighSpeedMachines()
               << "\t" << instance->getNbOfLowSpeedMachines()
               << "\t" << instance->getHighSpeed()
               << "\t" << instance->getLowSpeed()
               << "\t" << "CG"
               << "\t" << time_elapsed.count()
               << "\t" << time_limits.count()
               << "\t" << nbCallsDP
               << "\t" << nbCallsHeu
               << "\t" << getSumWjUj() << std::endl;
    outputFile.close();
}



