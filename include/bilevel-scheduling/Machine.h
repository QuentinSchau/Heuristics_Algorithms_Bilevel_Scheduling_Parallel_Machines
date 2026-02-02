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

#ifndef BILEVEL_SCHEDULING_MACHINE_H
#define BILEVEL_SCHEDULING_MACHINE_H

#include <vector>
#include <ostream>
#include <algorithm>
#include "Job.h"
#include "Math.h"
#include "BiSchException.h"

class Machine {

private:
    std::vector<Job> listAffectedJobs; // the list of jobs that be processed on the machine
    double speed; // the speed of the machine
    double sum_wj_Uj; // the sum of the number weighted tardy jobs on the machine
    double sum_Cj; // the sum of the completion time on the machine

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    /**
     * Default constructor
     */
    Machine() : listAffectedJobs(std::vector<Job>()), speed(0), sum_wj_Uj(0.0), sum_Cj(0.0) {}

    /**
     * Constructor a machine with the speed.
     * @param speed The speed of the current machine
     */
    explicit Machine(double speed) : listAffectedJobs(std::vector<Job>()), speed(speed), sum_wj_Uj(0.0), sum_Cj(0.0) {}

    explicit operator std::vector<unsigned int>() const {
        std::vector<unsigned int> machineSchedule(listAffectedJobs.size());
        std::transform(listAffectedJobs.begin(), listAffectedJobs.end(), machineSchedule.begin(), [](const Job &job) { return job.getIndex(); });
        return machineSchedule;
    }
    /********************/
    /*      METHODS     */
    /********************/

    /**
     * Method that resets the machine
     */
    void reset() {
        listAffectedJobs.clear();
        sum_wj_Uj = 0.0;
        sum_Cj = 0.0;
    }

    /**
     * Method that returns the number of scheduled jobs
     * @return The number of scheduled jobs
     */
    [[nodiscard]] size_t size() const { return listAffectedJobs.size(); }

    /**
     * Method that returns if the machine is empty
     * @return Machine is empty
     */
    bool empty() { return listAffectedJobs.empty(); }

    /**
     * Method that evaluate the scheduling on the machine. It compute the sum of completion times and the sum of weighted tardy jobs.
     */
    void evaluate() {
        double cumulativeWeightedTardyJob = 0.0;
        double cumulativeCompletionTimes = 0.0;
        double completionTimes = 0.0;
        size_t nbJobs = listAffectedJobs.size();
        for (size_t i = 0; i < nbJobs; ++i) {
            // compute cumulative completion times, a job at position i count (nbJobs-i) times
            cumulativeCompletionTimes += (double(nbJobs - i) * listAffectedJobs[i].getPi()) / speed;
            completionTimes += listAffectedJobs[i].getPi() / speed;
            cumulativeWeightedTardyJob +=
                    (listAffectedJobs[i].getDi() < completionTimes) ? listAffectedJobs[i].getWi() : 0.0;
        }
        sum_Cj = cumulativeCompletionTimes;
        sum_wj_Uj = cumulativeWeightedTardyJob;
    }


    /**
     * Method that adds a job to be scheduled on machine. This job is added at the given position.
     * @param pos Iterator before which the content will be inserted
     * @param affectedJob Job to insert
     */
    void add_job(unsigned int position, const Job &affectedJob) {
        listAffectedJobs.insert(listAffectedJobs.begin() + position, affectedJob);
    }

    /**
     * Method that adds a job to be scheduled on machine. This job is added at the end of the machine.
     * @param affectedJob The job to insert
     */
    void add_job(const Job &affectedJob) { listAffectedJobs.push_back(affectedJob); }


    /**
     * Method that reverses the job on machine
     */
    void reverse() { std::reverse(listAffectedJobs.begin(), listAffectedJobs.end()); }


    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] const std::vector<Job> &getAffectedJob() const { return listAffectedJobs; }

    [[nodiscard]] double getSpeed() const { return speed; }

    [[nodiscard]] double getSumWjUj() const { return sum_wj_Uj; }

    [[nodiscard]] double getSumCj() const { return sum_Cj; }

    /********************/
    /*      SETTER      */
    /********************/


    void setAffectedJob(const std::vector<Job> &affectedJob) { Machine::listAffectedJobs = affectedJob; }

    void setSpeed(double speed) { Machine::speed = speed; }

    void setSumWjUj(double sumWjUj) { sum_wj_Uj = sumWjUj; }

    void setSumCj(double sumCj) { sum_Cj = sumCj; }

    /************************/
    /*      OPERATORS       */
    /************************/

    Job &operator[](size_t pos) { return listAffectedJobs[pos]; }

    bool operator==(const Machine &M) const { return M.getAffectedJob() == listAffectedJobs && isEqual(M.speed,speed); }

};


inline std::ostream &operator<<(std::ostream &os, const Machine &machine) {

    double completionTimes = 0.0;
    size_t nbJobs = machine.size();
    os << "[";
    for (size_t j = 0; j < nbJobs; ++j) {
        Job currentJob = machine.getAffectedJob()[j];
        completionTimes += currentJob.getPi() / machine.getSpeed();
        os << currentJob;
        if (completionTimes > currentJob.getDi()) os << " **";
    }
    os << "]";
    return os;
}

#endif //BILEVEL_SCHEDULING_MACHINE_H
