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

#ifndef BILEVEL_SCHEDULING_JOB_H
#define BILEVEL_SCHEDULING_JOB_H

#include <iostream>
#include <iomanip>
#include <limits>
#include <random>
#include "BiSchException.h"

class Job {

private:
    double pi;
    double di;
    double wi;
    unsigned int index;
public:


    /**
     * Default constructor
     */
    Job() : pi(0.0), di(0.0), wi(0.0), index(0) {};

    /**
     * Constructor by processing time, due date and weight
     * @param pi The processing times
     * @param di The due date
     * @param wi The weight
     */
    Job(double pi, double di, double wi) : pi(pi), di(di), wi(wi), index(0) {}

    /********************/
    /*      GETTER      */
    /********************/

    [[nodiscard]] double getPi() const { return pi; }

    [[nodiscard]] double getDi() const { return di; }

    [[nodiscard]] double getWi() const { return wi; }

    [[nodiscard]] unsigned int getIndex() const { return index; }

    /*******************/
    /*      SETTER     */
    /*******************/

    void setIndex(unsigned int index) { Job::index = index; }

    void setPi(double pi) { Job::pi = pi; }

    void setDi(double tf, double rdd, double P, std::mt19937 &numGenerator) {
        auto infDi = P * (1 - tf - (rdd / 2));
        auto supDi = P * (1 - tf + (rdd / 2));
        std::uniform_int_distribution<> diDistribution(infDi, supDi);
        double newDi = double(diDistribution(numGenerator));
        unsigned int nbIteration = 0;
        while (newDi < 0 && nbIteration < 1000) {
            newDi = double(diDistribution(numGenerator));
            ++nbIteration;
        }
        // if we make 1000 iteration and the due date is still negative, then we set is value to 1.0
        if (nbIteration >= 1000) {
            newDi = 1.0;
            std::cerr << "Set default value for d_i=1.0 after 1000 random draw with negative value" << std::endl;
        }
        Job::di = newDi;
    }

    void setWi(double wi) { Job::wi = wi; }

    /************************/
    /*      OPERATORS       */
    /************************/

    static bool LPT_inv_EDD(const Job &lhs, const Job &rhs) {
        return (lhs.pi == rhs.pi) ? (lhs.di < rhs.di) : (lhs.pi > rhs.pi);
    }

    static bool WMST(const Job &lhs, const Job &rhs) {
        return double(lhs.di - lhs.pi) / double(lhs.wi) > double(rhs.di - rhs.pi) / double(rhs.wi);
    }

    static bool EDD(const Job &lhs, const Job &rhs) { return lhs.di < rhs.di; }

    bool operator==(const Job &J) const { return pi == J.pi && di == J.di && wi == J.wi && index == J.index; }

    bool operator!=(const Job &J) const { return !(J == *this); }

    /* Inequalities Operators compare processing times and due date in case of equalities */
    bool operator<(const Job &J) const { return (pi == J.pi) ? (di == J.di) ? (wi == J.wi) ? index < J.index : (wi < J.wi) : (di < J.di) : (pi < J.pi); }

    bool operator>(const Job &J) const { return (pi == J.pi) ? (di == J.di) ? (wi == J.wi) ? index > J.index : (wi > J.wi) : (di > J.di) : (pi > J.pi); }

    bool operator<=(const Job &J) const { return (pi == J.pi) ? (di == J.di) ? (wi == J.wi) ? index <= J.index : (wi <= J.wi) : (di <= J.di) : (pi <= J.pi); }

    bool operator>=(const Job &J) const { return (pi == J.pi) ? (di == J.di) ? (wi == J.wi) ? index >= J.index : (wi >= J.wi) : (di >= J.di) : (pi >= J.pi); }

};

inline std::ostream &operator<<(std::ostream &os, const Job &job) {
    os << std::setprecision(8) << "|id: " << job.getIndex() << " pi: " << job.getPi() << " di: " << job.getDi() << " wi: " << job.getWi() << " |";
    return os;
}

#endif //BILEVEL_SCHEDULING_JOB_H
