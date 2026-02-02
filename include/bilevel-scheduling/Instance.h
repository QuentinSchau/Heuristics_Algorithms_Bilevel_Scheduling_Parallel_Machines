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
//
//

#ifndef BILEVEL_SCHEDULING_INSTANCE_H
#define BILEVEL_SCHEDULING_INSTANCE_H

#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include "Job.h"
#include "BiSchException.h"
#include "Math.h"
#include <random>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

class Instance {

private:

    std::string instanceName;
    std::filesystem::path instancePath; // the path to the instance
    unsigned int nbJobs; // the nb of job
    unsigned int nbToSelectJob; // the number of jobs to select
    unsigned int nbOfHighSpeedMachines; // the nb of high speed machine
    unsigned int nbOfLowSpeedMachines; // the nb of high low machine
    double highSpeed;
    double lowSpeed;
    double maxPj; // the maximum value of processing times
    double sumWj; // the sum of all weights
    double sumPj; // the sum of all processing times
    std::vector<Job> listJobs; // the list of jobs
    std::vector<std::vector<Job>> listGrpJobs; // the list of grouped jobs
    std::vector<unsigned int> mapListJobToListGroupedJobs; // use a vector that map the index of a job to its index in the list of grouped jobs
    // the seed use for generate instance
    std::mt19937 numGenerator;
    // An array containing the minimum and maximum number of jobs for high-speed machines (first and second
    // elements), and the minimum and maximum number of jobs for low-speed machines (third and fourth elements).
    std::array<unsigned int, 4> minMaxNumberOfJobsByMachines;
    // An array containing the number of machines with minimum (first element) and maximum (second element)
    // on high-speed machine, and the number of machines with minimum (third element) and maximum (fourth element)
    // on low-speed machine
    std::array<int, 4> nbMachinesWithMinMaxJobs;
    bool firstBlockIsOnBothTypeMachine;
    unsigned int nbJobsToScheduleOnFirstBlock;
    // max number of location
    unsigned int maxNBLocation;

    //The set of all location where schedule jobs to be optimal for Sum Cj. A location is (indexMachine,Position in the machine)
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> E;

public:


    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    /**
     * Default Constructor
     * Set seed at 0 by default
     */
    explicit Instance() : maxPj(0), sumWj(0), listJobs(std::vector<Job>()) {
        std::random_device rd;//random number engine
        setSeed(rd());
    }

    /**
     * Constructor by instance's path. It construct an instance with setting the attribute path. Moreover,
     * the method check if the path have a corresponding file, otherwise an exception is throw
     * @param newInstancePath The path to set
     */
    explicit Instance(std::string &newInstancePath) : maxPj(0), sumWj(0), listJobs(std::vector<Job>()) {
        setInstancePath(newInstancePath);
        std::random_device rd;//random number engine
        setSeed(rd());
    }

    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Method that computes the list of locations using block structure. Blocks are created from RIGHT TO LEFT.
     * An element in a block represents an available position where a job can be scheduled on a machine.
     * @param listOfLocationForSumCj A list to add the locations to. A location is a pair (indexMachine,indexBlock)
     * @return The number of positions added.
    */
    size_t computeListAvailableLocation(std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &listOfLocationForSumCj) const;


    /**
     * Method that calculates the list of blocks using structure characterization.
     * The blocks are ordered chronologically, meaning all blocks contain positions where jobs can be scheduled from left to right (LEFT TO RIGHT).
     * An element in a block represents an available position where a job can be scheduled on a machine.
     * @param E The set of all locations where schedule jobs to be optimal for Sum Cj
     * @return The number of positions added.
     */
    size_t computeChronologicalLocations(std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &E) const;

    /**
     * This method computes the minimum and maximum number of jobs that can be scheduled on each machine type based on
     * the available positions. The set of location 'E' must be computed before (not checking it)
     * @return An array containing the minimum and maximum number of jobs for high-speed machines (first and second
     * elements), and the minimum and maximum number of jobs for low-speed machines (third and fourth elements).
    */
    std::array<unsigned int, 4> computeMinMaxNumberJobsOnMachines();

    /**
     * This method computes the minimum and maximum number of jobs that can be scheduled on each machine type based on
     * the available positions. If we do not know the min/max because the first block is on both type machine, then we
     * set the min nb of machine to -1 and the max nb of machine as the minimal number of machine with the maximal nb
     * of job inside. The set of location 'E' must be computed before (not checking it)
     * @return An array containing the minimum and maximum number of jobs for high-speed machines (first and second
     * elements), and the minimum and maximum number of jobs for low-speed machines (third and fourth elements).
    */
    std::array<int, 4> computeNumberMachineWithMinMaxJobs();


    /**
     * Method that generate a instance from the given parameters 
     * @param nbJobs The number of jobs to generate
     * @param nbOfHighSpeedMachines The number of high speed machines to create
     * @param nbOfHighLowMachines The number of low speed machines to create
     * @param highSpeed The high speed for machine
     * @param lowSpeed The low speed for machine
     * @param seed The seed to use for generate instance
     */
    void generateInstance(unsigned int nbJobs, unsigned int nbToSelectJob, unsigned int nbOfHighSpeedMachines, unsigned int nbOfLowSpeedMachines, double highSpeed, double lowSpeed, unsigned int seed);

    /**
     * Method that generate a instance from the given parameters with JSON Object
     * @param paramInstance The JSON Object that contains the parameters to generate the instance
     */
    void generateInstance(nlohmann::json &paramInstance);

    /**
     * Method that generates a Job with random values for pi, di, and wi.
     * These values come from a uniform distribution where the
     * lower bound and upper bound are defined by infPi and supPi for pi,
     * infDi and supDi for di, and infWi and supWi for wi.
     * @param infPi Lower bound of the uniform distribution to generate pi
     * @param supPi Upper bound of the uniform distribution to generate pi
     * @param infDi Lower bound of the uniform distribution to generate di
     * @param supDi Upper bound of the uniform distribution to generate di
     * @param infWi Lower bound of the uniform distribution to generate wi
     * @param supWi Upper bound of the uniform distribution to generate wi
     */
    Job generateJob(unsigned int infPi = 1, unsigned int supPi = 10, unsigned int infDi = 1, unsigned int supDi = 10, unsigned int infWi = 1, unsigned int supWi = 20);

    /**
     * Method that generates a Job with random values for pi, di, and wi.
     * These values come from a uniform distribution where the
     * lower bound and upper bound are defined by infPi and supPi for pi,
     *  and infWi and supWi for wi. For di, is generate throw a protocol
     * @param infPi Lower bound of the uniform distribution to generate pi
     * @param supPi Upper bound of the uniform distribution to generate pi
     * @param TF Lower bound of the uniform distribution to generate di
     * @param RDD Upper bound of the uniform distribution to generate di
     * @param infWi Lower bound of the uniform distribution to generate wi
     * @param supWi Upper bound of the uniform distribution to generate wi
     */
    Job generateJobWithoutDueDate(unsigned int infPi, unsigned int supPi, unsigned int infWi, unsigned int supWi);

    /**
     * Method that add a new job to the list of jobs of this instance
     * @param newJob The job to add
     */
    void add_job(Job &newJob) {
        newJob.setIndex(listJobs.size());
        listJobs.push_back(newJob);
    }

    /**
     * Method that clear the list of jobs
     */
    void clearListJobs() { listJobs.clear(); }

    /**
     * Sort job according the Longest Processing Time rule. If processing times are equals then sort according due date increasing.
     */
    void sort_by_LPT() {
        std::sort(listJobs.begin(), listJobs.end(), std::greater<>());
        reindexJobs();
    }

    /**
     * Sort job according the Shortest Processing Time rule. If processing times are equals then sort according due date increasing.
     */
    void sort_by_SPT() {
        std::sort(listJobs.begin(), listJobs.end());
        reindexJobs();
    }

    /**
     * Sort job according the Longest Processing Time rule. If processing times are equals then sort according due date decreasing.
     */
    void sort_by_LPT_inv_EDD() {
        std::sort(listJobs.begin(), listJobs.end(), Job::LPT_inv_EDD);
        reindexJobs();
    }

    /**
     * Sort job according the Weighted Minimum Slack Time rule, i.e., (d_1 - p_1)/w_1 <= (d_2 - p_2)/w_2 <= ... <= (d_n - p_n)/w_n.
     * This rule measure the "urgency" of a job with its weight.
     */
    void sort_by_WMST() {
        std::sort(listJobs.begin(), listJobs.end(), Job::WMST);
        reindexJobs();
    }

    /**
     * Method that re-index job after a sorting
     */
    void reindexJobs() {
        for (unsigned int j = 0; j < listJobs.size(); j++) {
            listJobs[j].setIndex(j);
        }
        setListGroupedJobs();
    }

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] const std::filesystem::path &getInstancePath() const { return instancePath; }

    [[nodiscard]] const std::string &getInstanceName() const { return instanceName; }

    [[nodiscard]] unsigned int getNbJobs() const { return nbJobs; }

    [[nodiscard]] unsigned int getNbToSelectJob() const { return nbToSelectJob; }

    [[nodiscard]] unsigned int getNbOfHighSpeedMachines() const { return nbOfHighSpeedMachines; }

    [[nodiscard]] unsigned int getNbOfLowSpeedMachines() const { return nbOfLowSpeedMachines; }

    [[nodiscard]] unsigned int getNbMachines() const { return (nbOfHighSpeedMachines + nbOfLowSpeedMachines); }

    [[nodiscard]] double getHighSpeed() const { return highSpeed; }

    [[nodiscard]] double getLowSpeed() const { return lowSpeed; }

    [[nodiscard]] double getMaxPj() const { return maxPj; }

    [[nodiscard]] double getSumWj() const { return sumWj; }

    [[nodiscard]] double getSumPj() const { return sumPj; }

    [[nodiscard]] const std::vector<Job> &getListJobs() { return listJobs; }

    [[nodiscard]] const std::vector<std::vector<Job>> &getListGrpJobs() { return listGrpJobs; }

    [[nodiscard]] const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> &getE() const {
        return E;
    }

    [[nodiscard]] bool isFirstBlockIsOnBothTypeMachine() const { return firstBlockIsOnBothTypeMachine; }

    [[nodiscard]] unsigned int getNbJobsToScheduleOnFirstBlock() const { return nbJobsToScheduleOnFirstBlock; }

    /**
     * @return The minimum number of jobs for high-speed machines
     */
    [[nodiscard]] unsigned int getMinNbJobsOnHS() const { return minMaxNumberOfJobsByMachines[0]; }

    /**
     * @return The maximum number of jobs for high-speed machines
     */
    [[nodiscard]] unsigned int getMaxNbJobsOnHS() const { return minMaxNumberOfJobsByMachines[1]; }

    /**
     * @return The minimum number of jobs for low-speed machines
     */
    [[nodiscard]] unsigned int getMinNbJobsOnLS() const { return minMaxNumberOfJobsByMachines[2]; }

    /**
     * @return The maximum number of jobs for low-speed machines
     */
    [[nodiscard]] unsigned int getMaxNbJobsOnLS() const { return minMaxNumberOfJobsByMachines[3]; }

    /**
     * @return The number of machine with the minimal number of jobs for high-speed machines
     */
    [[nodiscard]] int getNbMachineWithMinJobsOnHS() const { return nbMachinesWithMinMaxJobs[0]; }

    /**
     * @return The number of machine with the maximum number of jobs for high-speed machines
     */
    [[nodiscard]] int getNbMachineWithMaxJobsOnHS() const { return nbMachinesWithMinMaxJobs[1]; }

    /**
     * @return The number of machine with the minimal number of jobs for low-speed machines
     */
    [[nodiscard]] int getNbMachineWithMinJobsOnLS() const { return nbMachinesWithMinMaxJobs[2]; }

    /**
     * @return The number of machine with the maximum number of jobs for low-speed machines
     */
    [[nodiscard]] int getNbMachineWithMaxJobsOnLS() const { return nbMachinesWithMinMaxJobs[3]; }


    [[nodiscard]] unsigned int getMaxNbLocation() const { return maxNBLocation; }

    [[nodiscard]] const std::vector<unsigned int> &getMapListJobToListGroupedJobs() const {return mapListJobToListGroupedJobs;}

    [[nodiscard]] bool isIdenticalJob(unsigned int indexJob) const {return listGrpJobs[mapListJobToListGroupedJobs[indexJob]].size() >= 2;}

    [[nodiscard]] unsigned int getIndexIdenticalGroupOfJob(unsigned int indexJob) const {return mapListJobToListGroupedJobs[indexJob];}


    /********************/
    /*      SETTER      */
    /********************/


    void setInstanceName(const std::string &newInstanceName) { Instance::instanceName = newInstanceName; }

    void setInstancePath(const std::string &newInstancePath) {
        std::filesystem::path newPath = std::filesystem::path(newInstancePath);
        instancePath = newPath;
        std::filesystem::directory_entry parentDir{newPath.lexically_normal().parent_path()};
        if (std::filesystem::exists(newPath)) instancePath = newPath;
        else {
            // if the dir doesn't exist
            if (!parentDir.exists()) {
                std::filesystem::create_directories(newPath.lexically_normal().parent_path());
            }
        }
        instanceName = instancePath.stem();
    }

    void setNbJobs(unsigned int newNb_jobs) {
        nbJobs = newNb_jobs;
    }

    void setNbToSelectJob(unsigned int newNb_to_select_job) {
        nbToSelectJob = newNb_to_select_job;
    }

    void setNbOfHighSpeedMachines(unsigned int newNb_of_high_speed_machines) {
        nbOfHighSpeedMachines = newNb_of_high_speed_machines;
    }

    void setNbOfLowSpeedMachines(unsigned int newNb_of_low_speed_machines) {
        nbOfLowSpeedMachines = newNb_of_low_speed_machines;
    }

    void setHighSpeed(double newHigh_speed) {
        highSpeed = newHigh_speed;
    }

    void setLowSpeed(double newLow_speed) {
        lowSpeed = newLow_speed;
    }

    void setMaxPj(double newMax_pj) {
        maxPj = newMax_pj;
    }

    void setSumWj(double newSum_wj) {
        sumWj = newSum_wj;
    }

    void setSumPj(double newSum_pj) {
        sumPj = newSum_pj;
    }

    void setListGroupedJobs() {
        listGrpJobs.clear();
        if (!listJobs.empty()) {
            listGrpJobs.emplace_back(std::vector<Job>({listJobs[0]}));
            for (unsigned int j = 1; j < listJobs.size(); ++j) {
                if (listGrpJobs.back().back().getPi() == listJobs[j].getPi()) listGrpJobs.back().push_back(listJobs[j]);
                else listGrpJobs.emplace_back(std::vector<Job>({listJobs[j]}));
            }
        }
    }

    void setConstant() {
        setListGroupedJobs();
        // create the mapping between the list of jobs and the list of grouped jobs
        mapListJobToListGroupedJobs.assign(nbJobs, 0);
        unsigned int indexInListGrpJobs = 0;
        auto itJob = listJobs.cbegin();
        auto predFindJobInGrp = [&itJob](const Job &jobInGroup) { return isEqual(jobInGroup.getPi(),itJob->getPi()); };
        for (; itJob != listJobs.cend(); itJob++) {
            assert(indexInListGrpJobs < listGrpJobs.size());
            auto& groupIdenticalJobs = listGrpJobs[indexInListGrpJobs];
            // if the job is not in the group then is necessary in the next group
            if (std::find_if(groupIdenticalJobs.cbegin(), groupIdenticalJobs.cend(), predFindJobInGrp) == groupIdenticalJobs.cend()) {
                indexInListGrpJobs++;
            }
            mapListJobToListGroupedJobs[itJob->getIndex()] = indexInListGrpJobs;
        }
        //The set of all location where schedule jobs to be optimal for Sum Cj. A location is (indexMachine,IndexBlock)
        maxNBLocation = computeListAvailableLocation(E);
        minMaxNumberOfJobsByMachines = computeMinMaxNumberJobsOnMachines();
        computeChronologicalLocations(E);
        nbMachinesWithMinMaxJobs = computeNumberMachineWithMinMaxJobs();
        firstBlockIsOnBothTypeMachine = getNbMachineWithMinJobsOnHS() == -1;

    };

    void setSeed(unsigned int seed) { Instance::numGenerator = std::mt19937(seed); }

};

/************************/
/*      OPERATORS       */
/************************/

inline std::ostream &operator<<(std::ostream &os, Instance &instance) {

    os << "Instance Path : " << instance.getInstancePath() << std::endl
       << "Instance name : " << instance.getInstanceName() << std::endl
       << "Number of jobs : " << instance.getNbJobs() << std::endl
       << "Number of high speed machines : " << instance.getNbOfHighSpeedMachines() << " with the speed of : "
       << instance.getHighSpeed() << std::endl
       << "Number of low speed machines : " << instance.getNbOfLowSpeedMachines() << " with the speed of : "
       << instance.getLowSpeed() << std::endl
       << "List of jobs : [" << std::endl;

    for (const Job &job: instance.getListJobs()) {
        os << job << std::endl;
    }
    os << "]" << std::endl;
    return os;
}

#endif //BILEVEL_SCHEDULING_INSTANCE_H
