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

#include "Instance.h"


size_t Instance::computeListAvailableLocation(
        std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &listOfLocationForSumCj) const {
    size_t nbLocations = 1;

    unsigned int nbMachines = getNbMachines(); // the nb of machines
    //the weights array
    Eigen::ArrayXXd weights = Eigen::ArrayXXd::Zero(nbMachines, getNbToSelectJob());

    // Compute weights of high speed machine
    for (unsigned int i = 0; i < getNbToSelectJob(); ++i)
        weights.block(0, i, getNbOfHighSpeedMachines(), 1).setConstant(
                double(getNbToSelectJob() - i) / getHighSpeed());

    // Compute weights of low speed machine
    for (unsigned int i = 0; i < getNbToSelectJob(); ++i)
        weights.block(getNbOfHighSpeedMachines(), i, getNbOfLowSpeedMachines(), 1).setConstant(
                double(getNbToSelectJob() - i) / getLowSpeed());

    double maxWeight = weights.maxCoeff(); // the maxWeight
    unsigned nbCols = weights.cols();
    unsigned int row, col;
    double minWeight = weights.minCoeff(&row, &col);
    weights(row, col) = maxWeight;

    //It's the list of location available for one coefficient
    std::vector<std::pair<unsigned int, unsigned int>> locationsForSumCj;
    // HERE col is the column index of a position. It's indexed by the LEFT to the RIGHT
    locationsForSumCj.emplace_back(row, nbCols-col);

    unsigned int maxNbElements = weights.rows() * weights.cols();
    // loop over the number of available locations
    while (nbLocations < maxNbElements) {
        double minValue = weights.minCoeff(&row, &col);
        if (minValue != minWeight) {
            listOfLocationForSumCj.emplace_back(locationsForSumCj);
            locationsForSumCj.clear();
            minWeight = minValue;
            if (nbLocations >= getNbToSelectJob()) {
                break;
            }
        }
        weights(row, col) = maxWeight;
        locationsForSumCj.emplace_back(row, nbCols-col);
        ++nbLocations;
    }

    return nbLocations;
}

std::array<unsigned int, 4> Instance::computeMinMaxNumberJobsOnMachines() {
    std::array<unsigned int,4> minMaxNumberOfJobs{0,0,0,0};

    // the number of position that we have selected
    unsigned int numberPositionSelected = 0;
    for (auto & itBlock : E) {
        numberPositionSelected += itBlock.size();
        // check which machines there are present in the block
        bool haveLowSpeed = false;
        bool haveHighSpeed = false;
        for (auto &[indexMachine, indexBlock]: itBlock) {
            if (indexMachine < getNbOfHighSpeedMachines()) haveHighSpeed = true;
            else if (indexMachine >= getNbOfHighSpeedMachines()) haveLowSpeed = true;
        }
        // if we have not already selected more than necessary jobs then increase both min and max
        if (numberPositionSelected <= getNbToSelectJob() ) {
            if (haveHighSpeed) {
                ++minMaxNumberOfJobs[0];
                ++minMaxNumberOfJobs[1];
            }
            if (haveLowSpeed) {
                ++minMaxNumberOfJobs[2];
                ++minMaxNumberOfJobs[3];
            }
        }else {
            if (haveHighSpeed)
                ++minMaxNumberOfJobs[1];
            if (haveLowSpeed)
                ++minMaxNumberOfJobs[3];
        }
    }
    return minMaxNumberOfJobs;
}

/**
 * Method that compute the nb of machine with the max or the min number of jobs on it
 * @return the array with the number of machine [M1-,M1+,M0-,M0+]
 */
std::array<int, 4> Instance::computeNumberMachineWithMinMaxJobs() {
    std::array<int,4> nbMachineMaxMin{0,0,0,0};
    // to know the number of machine with min or max, we should check first if the
    // first block is on high/low speed machine
    nbJobsToScheduleOnFirstBlock = E[0].size() + nbToSelectJob - maxNBLocation;
    // if we are on high-speed machines
    if (E[0].front().first < nbOfHighSpeedMachines && E[0].back().first < nbOfHighSpeedMachines) {
        // the number of job to schedule on this block will lead to a machine schedule
        nbMachineMaxMin[1] = static_cast<int>(nbJobsToScheduleOnFirstBlock);
        nbMachineMaxMin[0] = static_cast<int>(nbOfHighSpeedMachines) - static_cast<int>(nbJobsToScheduleOnFirstBlock);
        // then on low speed machine there is only one possibility, i.e. a maximal number of jobs, which is egal to the number of low-speed machines
        nbMachineMaxMin[2] = 0;
        nbMachineMaxMin[3] = static_cast<int>(nbOfLowSpeedMachines);
    }
    // we are on low-speed machines
    else if (E[0].front().first >= nbOfHighSpeedMachines && E[0].back().first >= nbOfHighSpeedMachines){
        nbMachineMaxMin[3] = static_cast<int>(nbJobsToScheduleOnFirstBlock);
        nbMachineMaxMin[2] = static_cast<int>(nbOfLowSpeedMachines) - static_cast<int>(nbJobsToScheduleOnFirstBlock);
        // then on high speed machine there is only one possibility, i.e. a maximal number of jobs, which is egal to the number of high-speed machines
        nbMachineMaxMin[0] = 0;
        nbMachineMaxMin[1] = static_cast<int>(nbOfHighSpeedMachines);
    }else{
        // the first block is on both low/speed machines
        nbMachineMaxMin[0] = -1;
        nbMachineMaxMin[2] = -1;
        // if there is less low-speed machine that the nb of jobs to schedule on the first block
        if (nbJobsToScheduleOnFirstBlock > nbOfLowSpeedMachines) {
            // then we have at least "nb jobs on first block - nb low-speed machine" high-speed machines with the highest job count
            nbMachineMaxMin[1] = static_cast<int>(nbJobsToScheduleOnFirstBlock) - static_cast<int>(nbOfLowSpeedMachines);
        }else{
            nbMachineMaxMin[1] = 0;
        }

        // if there is less high-speed machine that the nb of jobs to schedule on the first block
        if (nbJobsToScheduleOnFirstBlock > nbOfHighSpeedMachines) {
            // then we have at least "nb jobs on first block - nb low-speed machine" high-speed machines with the highest job count
            nbMachineMaxMin[3] = static_cast<int>(nbJobsToScheduleOnFirstBlock) - static_cast<int>(nbOfHighSpeedMachines);
        }else{
            nbMachineMaxMin[3] = 0;
        }
    }

    return nbMachineMaxMin;
}

size_t Instance::computeChronologicalLocations(std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &E) const {
    size_t nbLocation = maxNBLocation;
    // reverse the list in order to have batch from the LEFT to the RIGHT
    std::reverse(E.begin(), E.end());

    // keep the last anti-chronologicalLocation
    std::vector<unsigned int> lastAntiChronological(getNbMachines(),0);

    // get the first ones
    for (auto &location: E[0]) {
        location.second = 0;
        lastAntiChronological[location.first] = location.second + 1;
    }

    // change the position of other batches
    for (auto it = E.begin() + 1; it != E.end(); ++it) {
        // change the location
        for (auto &location: *it) {
            location.second = lastAntiChronological[location.first];
            lastAntiChronological[location.first] = location.second +1 ;
        }
    }
    return nbLocation;

}

Job Instance::generateJob(unsigned int infPi, unsigned int supPi, unsigned int infDi, unsigned int supDi,
                          unsigned int infWi, unsigned int supWi) {

    std::uniform_int_distribution<> piDistribution(infPi, supPi);
    std::uniform_int_distribution<> diDistribution(infDi, supDi);
    std::uniform_int_distribution<> wiDistribution(infWi, supWi);
    double pi = double(piDistribution(numGenerator));
    double di = double(diDistribution(numGenerator));
    double wi = double(wiDistribution(numGenerator));
    return Job(pi, di, wi);
}

Job Instance::generateJobWithoutDueDate(unsigned int infPi, unsigned int supPi, unsigned int infWi,
                                        unsigned int supWi) {
    std::uniform_int_distribution<> piDistribution(infPi, supPi);
    std::uniform_int_distribution<> wiDistribution(infWi, supWi);
    double pi = double(piDistribution(numGenerator));
    double wi = double(wiDistribution(numGenerator));
    Job newJob = Job();
    newJob.setPi(pi);
    newJob.setWi(wi);
    return newJob;
}
void Instance::generateInstance(unsigned int nbJobs, unsigned int nbToSelectJob, unsigned int nbOfHighSpeedMachines,
                                unsigned int nbOfLowSpeedMachines, double highSpeed, double lowSpeed,
                                unsigned int seed) {

    // set attributes
    setNbJobs(nbJobs);
    setNbToSelectJob(nbToSelectJob);
    setNbOfHighSpeedMachines(nbOfHighSpeedMachines);
    setNbOfLowSpeedMachines(nbOfLowSpeedMachines);
    setHighSpeed(highSpeed);
    setLowSpeed(lowSpeed);
    setSeed(seed);

    // generate Jobs
    listJobs.reserve(nbJobs);
    for (unsigned int i = 0; i < nbJobs; ++i) {
        Job newJob = generateJob();
        newJob.setIndex(listJobs.size());
        listJobs.push_back(newJob);
        if (newJob.getPi() > maxPj) maxPj = newJob.getPi();
        sumWj += newJob.getWi();
        sumPj += newJob.getPi();
    }
}

void Instance::generateInstance(nlohmann::json &paramInstance) {
    // set the number of N jobs
    if (paramInstance.contains("N")) {
        if (paramInstance["N"].is_number_unsigned()) setNbJobs(paramInstance["N"]);
        else throw std::invalid_argument(R"(The "N" must be an unsigned integer)");
    }
    // set the number of n jobs
    if (paramInstance.contains("n")) {
        if (paramInstance["n"].is_number_unsigned()) setNbToSelectJob(paramInstance["n"]);
        else throw std::invalid_argument(R"(The "n" must be an unsigned integer)");
    }

    if (nbJobs < nbToSelectJob)
        throw BiSchException("The number of job 'N' must be greater or equal to the number of job to select 'n'");
    // set the number of high speed machines
    if (paramInstance.contains("M_max")) {
        if (paramInstance["M_max"].is_number_unsigned()) setNbOfHighSpeedMachines(paramInstance["M_max"]);
        else throw std::invalid_argument(R"(The "M_max" must be an unsigned integer)");
    }
    // set the number of low  speed machines
    if (paramInstance.contains("M_0")) {
        if (paramInstance["M_0"].is_number_unsigned()) setNbOfLowSpeedMachines(paramInstance["M_0"]);
        else throw std::invalid_argument(R"(The "M_0" must be an unsigned integer)");
    }
    // set the high speed
    if (paramInstance.contains("V_max")) {
        if (paramInstance["V_max"].is_number_float()) setHighSpeed(paramInstance["V_max"]);
        else throw std::invalid_argument(R"(The "V_max" must be a float : )");
    }
    // set the low speed
    if (paramInstance.contains("V_0")) {
        if (paramInstance["V_0"].is_number_float()) setLowSpeed(paramInstance["V_0"]);
        else throw std::invalid_argument(R"(The "V_0" must be a float)");
    }

    //set pi distribution
    unsigned int infPi = 1;
    unsigned int supPi = 10;

    if (paramInstance.contains("pi")) {
        if (paramInstance["pi"].contains("inf")) {
            if (paramInstance["pi"]["inf"].is_number_unsigned()) infPi = paramInstance["pi"]["inf"];
            else throw std::invalid_argument(R"(The "inf" must be an unsigned integer in the "pi" object)");
        }
        if (paramInstance["pi"].contains("sup")) {
            if (paramInstance["pi"]["sup"].is_number_unsigned()) supPi = paramInstance["pi"]["sup"];
            else throw std::invalid_argument(R"(The "sup" must be an unsigned integer in the "pi" object)");
        }
    }

    bool noIdenticalJobs = false;
    if (paramInstance.contains("noIdenticalJobs")) {
        if (paramInstance["noIdenticalJobs"].is_boolean()) noIdenticalJobs = paramInstance["noIdenticalJobs"];
        else throw std::invalid_argument(R"(The "noIdenticalJobs" must be a boolean)");
    }

    std::vector<bool> identicalJobs;
    if (noIdenticalJobs) {
        if (supPi < nbJobs)
            throw std::invalid_argument(
                    R"(The max of p_j value must be greater than the number of jobs if you want no identical jobs)");
        else identicalJobs.resize(supPi, false);
    }

    //set di distribution
    unsigned int infDi = 1;
    unsigned int supDi = 10;
    double TF = -1;
    double RDD = -1;
    if (paramInstance.contains("di")) {
        if (paramInstance["di"].contains("inf")) {
            if (paramInstance["di"]["inf"].is_number_unsigned()) infDi = paramInstance["di"]["inf"];
            else throw std::invalid_argument(R"(The "inf" must be an unsigned integer in the "di" object)");
        }
        if (paramInstance["di"].contains("sup")) {
            if (paramInstance["di"]["sup"].is_number_unsigned()) supDi = paramInstance["di"]["sup"];
            else throw std::invalid_argument(R"(The "sup" must be an unsigned integer in the "di" object)");
        }
        // if we use a protocol
        if (paramInstance["di"].contains("TF")) {
            if (paramInstance["di"]["TF"].is_number_float()) TF = paramInstance["di"]["TF"];
            else throw std::invalid_argument(R"(The "TF" must be a float in the "di" object)");
        }
        if (paramInstance["di"].contains("RDD")) {
            if (paramInstance["di"]["RDD"].is_number_float()) RDD = paramInstance["di"]["RDD"];
            else throw std::invalid_argument(R"(The "RDD" must be a float in the "di" object)");
        }
    }

    //set wi distribution
    unsigned int infWi = 1;
    unsigned int supWi = 10;

    if (paramInstance.contains("wi")) {
        if (paramInstance["wi"].contains("inf")) {
            if (paramInstance["wi"]["inf"].is_number_unsigned()) infWi = paramInstance["wi"]["inf"];
            else throw std::invalid_argument(R"(The "inf" must be an unsigned integer in the "wi" object)");
        }
        if (paramInstance["wi"].contains("sup")) {
            if (paramInstance["wi"]["sup"].is_number_unsigned()) supWi = paramInstance["wi"]["sup"];
            else throw std::invalid_argument(R"(The "sup" must be an unsigned integer in the "wi" object)");
        }
    }

    // generate Jobs
    listJobs.reserve(nbJobs);
    sumPj = 0.0;
    sumWj = 0.0;
    for (unsigned int i = 0; i < nbJobs; ++i) {
        Job newJob;
        // if TF and RDD are defined
        do {
            if (TF != -1 && RDD != -1) newJob = generateJobWithoutDueDate(infPi, supPi,infWi, supWi);
            else newJob = generateJob(infPi, supPi, infDi, supDi, infWi, supWi);
        } while (noIdenticalJobs && identicalJobs[static_cast<unsigned int>(newJob.getPi())]);
        if (noIdenticalJobs) identicalJobs[static_cast<unsigned int>(newJob.getPi())] = true;
        newJob.setIndex(listJobs.size());
        listJobs.push_back(newJob);
    if (newJob.getPi() > maxPj) maxPj = newJob.getPi();
        sumWj += newJob.getWi();
        sumPj += newJob.getPi();
    }
    double P = (sumPj/(nbOfHighSpeedMachines*highSpeed+nbOfLowSpeedMachines*lowSpeed))*(double(nbToSelectJob)/double(nbJobs));
    // set the due date according the TF, RDD and sum Pj
    for (auto &job: listJobs) {
        job.setDi(TF, RDD, P, numGenerator);
    }
    setConstant();
}










