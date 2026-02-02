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
// Created by schau on 1/15/24.
//

#include "Parser.h"
#include "BiSchException.h"
#include <iostream>
#include <fstream>
#include <stdexcept>


Parser::Parser() {}

Instance Parser::readFromFile(std::string &filePath) const {
    Instance newInstance = Instance(filePath);
    std::fstream fileStream(newInstance.getInstancePath().lexically_normal(), std::fstream::in);
    std::string line; // new line
    // compute the sum of weights
    double sumWj = 0.0;
    double sumPj = 0.0;
    // open and read the file
    if (fileStream.is_open()) {
        while (std::getline(fileStream, line)) {
            auto pos = line.find(":");
            // if we don't find the ':' then we read jobs
            if (pos == std::string::npos) {
                // create the job, the value should be separate with \t
                std::istringstream stream(line);
                double pj, dj, wj;
                stream >> pj;
                stream >> dj;
                stream >> wj;
                Job newJob = Job(pj, dj, wj);
                newInstance.add_job(newJob);
                if (pj > newInstance.getMaxPj()) newInstance.setMaxPj(pj);
                sumWj += wj;
                sumPj += pj;
            } else {
                // we read so attributes
                std::string attribute = line.substr(0, pos);
                if (attribute == "name") newInstance.setInstanceName(line.substr(pos + 1));
                else if (attribute == "N") newInstance.setNbJobs(std::stoul(line.substr(pos + 1)));
                else if (attribute == "n") newInstance.setNbToSelectJob(std::stoul(line.substr(pos + 1)));
                else if (attribute == "M_max") newInstance.setNbOfHighSpeedMachines(std::stoul(line.substr(pos + 1)));
                else if (attribute == "M_0") newInstance.setNbOfLowSpeedMachines(std::stoul(line.substr(pos + 1)));
                else if (attribute == "V_max") newInstance.setHighSpeed(std::stof(line.substr(pos + 1)));
                else if (attribute == "V_0") newInstance.setLowSpeed(std::stof(line.substr(pos + 1)));
            }
        }
    } else
        throw BiSchException(std::string("Can't open the file ").append(newInstance.getInstancePath().lexically_normal().string()).c_str());
    fileStream.close();
    // set the sum wj
    newInstance.setSumWj(sumWj);
    // set the sum pj
    newInstance.setSumPj(sumPj);

    // check if we have the right number of created job
    if (newInstance.getNbJobs() != newInstance.getListJobs().size())
        throw std::invalid_argument("The number of jobs is not equals to N");
    if (newInstance.getNbToSelectJob() == 0) throw std::invalid_argument("The number of jobs to select is not defined");
    // group all identical jobs
    newInstance.setListGroupedJobs();
    newInstance.setConstant();
    return newInstance;
}

void Parser::serializeInstance(Instance &instance) {
    instance.sort_by_SPT();
    std::fstream fileStream(instance.getInstancePath().lexically_normal().string(), std::fstream::out );
    if (fileStream.is_open()) {
        fileStream << std::setprecision(5) << "name:" << instance.getInstanceName() << std::endl
                   << "N:" << instance.getNbJobs() << std::endl
                   << "n:" << instance.getNbToSelectJob() << std::endl
                   << "M_max:" << instance.getNbOfHighSpeedMachines() << std::endl
                   << "M_0:" << instance.getNbOfLowSpeedMachines() << std::endl
                   << "V_max:" << instance.getHighSpeed() << std::endl
                   << "V_0:" << instance.getLowSpeed() << std::endl
                   << "Jobs:" << std::endl;

        for (const Job &job: instance.getListJobs()) {
            fileStream << job.getPi() << "\t" << job.getDi() << "\t" << job.getWi() << "\t" << std::endl;
        }
    } else throw BiSchException(std::string("Can't open the file ").append(instance.getInstancePath().lexically_normal().string()).c_str());
    fileStream.close();
}

void Parser::generateInstance(nlohmann::json &object) {

    Instance newInstance;

    // set the seed for generate
    if (object.contains("seed")) {
        if (object["seed"].is_number_unsigned()) newInstance.setSeed(object["seed"]);
        else throw std::invalid_argument(R"(The seed is not a unsigned int)");
    }

    // loop over each instances
    if (object.contains("instances")) {
        // create each kind of instance
        for (auto &paramInstance: object["instances"]) {
            // check if there is a base Path
            std::string basePath;
            if (paramInstance.contains("basePath")) {
                if (paramInstance["basePath"].is_string()) basePath = paramInstance["basePath"];
                else throw std::invalid_argument(R"(The base Path is not a string)");
            } else{
                basePath = std::filesystem::current_path().lexically_normal().string();
            }
            std::cout << "All instances will be generated at : " << basePath << std::endl;
            unsigned int nbGeneratedInstance = 0;
            unsigned int nbInstanceToGenerate = 1;
            // get the number of instance that we have to generate
            if (paramInstance.contains("numberInstance")) {
                if (paramInstance["numberInstance"].is_number_unsigned()) nbInstanceToGenerate = paramInstance["numberInstance"];
                else throw std::invalid_argument(R"(The "numberInstance" must be an unsigned integer)");
            }
            // we need to have the object "paramInstance"
            if (!paramInstance.contains("paramInstance")) throw std::invalid_argument(R"(The "paramInstance" is not defined)");

            for (unsigned int newInstanceLoop = 0; newInstanceLoop < nbInstanceToGenerate; ++newInstanceLoop) {
                // convert tf and rdd with 1 decimal
                char buffer[4];
                std::snprintf(buffer, 4, "%.1f", paramInstance["paramInstance"]["di"]["TF"].template get<double>());
                std::string tf(buffer);

                std::snprintf(buffer, 4, "%.1f", paramInstance["paramInstance"]["di"]["RDD"].template get<double>());
                std::string rdd(buffer);
                // convert Vmax and V0 with 1 decimal
                std::snprintf(buffer, 4, "%.1f", paramInstance["paramInstance"]["V_max"].template get<double>());
                std::string vmax(buffer);

                std::snprintf(buffer, 4, "%.1f", paramInstance["paramInstance"]["V_0"].template get<double>());
                std::string v0(buffer);
                // set the path
                std::string path = basePath;
                path.append("instance")
                    .append(std::to_string(nbGeneratedInstance))
                    .append("_n_").append(std::to_string(paramInstance["paramInstance"]["n"].template get<unsigned int>()))
                    .append("_N_").append(std::to_string(paramInstance["paramInstance"]["N"].template get<unsigned int>()))
                    .append("_tf_").append(tf)
                    .append("_rdd_").append(rdd)
                    .append("_mMax_").append(std::to_string(paramInstance["paramInstance"]["M_max"].template get<unsigned int>()))
                    .append("_m0_").append(std::to_string(paramInstance["paramInstance"]["M_0"].template get<unsigned int>()))
                    .append("_Vmax_").append(vmax)
                    .append("_V0_").append(v0)
                    .append(".txt");
                newInstance.setInstancePath(path);
                newInstance.generateInstance(paramInstance["paramInstance"]);
                serializeInstance(newInstance);
                ++nbGeneratedInstance;
                newInstance.clearListJobs();
            }
        }
    }else throw std::invalid_argument(R"(The generate config JSON must have an instance object)");

}


