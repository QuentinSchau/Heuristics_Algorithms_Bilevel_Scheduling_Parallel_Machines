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
// Created by schau on 1/22/25.
//

#include "Heuristic.h"

Heuristic::Heuristic() = default;

Heuristic::Heuristic(Instance *instance) : ISolver(instance) {
    initializeStructure();
}

Heuristic::Heuristic(Instance *instance, nlohmann::json &object) : ISolver(instance) {
    if (object.contains("name")) {
        if (object["name"] == "Heuristic") {
            if (object.contains("timeLimits")) {
                if (object["timeLimits"].is_number_unsigned()) setTimeLimit(object["timeLimits"]);
                else throw std::invalid_argument(R"(The attribute "timeLimits" of JSON object must be an "unsigned int" value for the constructor of Heuristic solver)");
            }
            if (object.contains("verbose")) {
                if (object["verbose"].is_number_unsigned()) setVerbose(object["verbose"].template get<char>());
                else throw std::invalid_argument(R"(The attribute "verbose" of JSON object must be an "unsigned integer" value for the constructor of Heuristic solver)");
            }
        } else throw std::invalid_argument(R"(The JSON object have not the right name to instance a Heuristic solver)");
    } else throw std::invalid_argument(R"(The JSON object have not a name to instance a Heuristic solver)");
    initializeStructure();
}

void Heuristic::initializeStructure() {
    for (unsigned int indexLoopMachine = 0; indexLoopMachine < instance->getNbMachines(); indexLoopMachine++) {
        std::vector<double> costForBlocks(instance->getMaxNbJobsOnHS(), 0.0);
        costAtEachBlock.emplace_back(costForBlocks);
        // create matrix of size m * (N+m-1) // m-1 because in worst case we have to select 1 machines among m, so we create m-1 dummy jobs with weight 0
        std::vector<double> costOfJobs(instance->getNbJobs()+instance->getNbMachines(), instance->getSumWj());
        costMatrix.emplace_back(costOfJobs);
    }
    listAvailableIndexJob.reserve(instance->getNbJobs());
    auto N = instance->getNbJobs();
    auto m = instance->getNbMachines();
    smallEpsilon = 1.0 / (static_cast<double>(N * (N + 1) / 2 + N * N + 2 * m * N));
}

Heuristic::~Heuristic() = default;

void Heuristic::solve() {

    // start time to measure performance
    const auto start = std::chrono::steady_clock::now();

    auto listOfJobs = instance->getListJobs();
    std::sort(listOfJobs.begin(), listOfJobs.end(), std::greater<>());
    Solution newSol = Solution::solveSumCjCriteria(listOfJobs, instance);
    upgradeSolutionWithHeuristic(newSol, listOfJobs);
    *solution = newSol;
    solution->evaluate();
    const auto end{std::chrono::steady_clock::now()};

    // stop time to measure performance
    time_elapsed = std::chrono::duration<double>{end - start};
    if (verbose >= 1)
        std::cout << "Heuristic is over after " << time_elapsed.count() << " seconds" << std::endl << "The objective value is : " << solution->getSumWjUj() << std::endl;
}

void Heuristic::printOutput(std::string &fileOutputName, std::ofstream &outputFile) {
    bool fileExists = std::filesystem::exists(fileOutputName);
    auto filePath = std::filesystem::path(fileOutputName);
    std::filesystem::create_directories(filePath.lexically_normal().parent_path());
    outputFile.open(fileOutputName, std::ios::out | std::ios::app | std::ios::ate);
    // print header
    if (!fileExists)
        outputFile << "InstanceName" << "\t" << "InstancePath" << "\t" << "N" << "\t" << "n" << "\t" << "m_Max" << "\t" << "m_0" << "\t" << "V_max" << "\t" << "V_0" << "\t" << "Method" << "\t" <<
                   "Predictor" << "\t" << "Time" << "\t" << "LimitTime" << "\t" << "Objective" << std::endl;

    // write value
    outputFile << instance->getInstanceName() << "\t" << instance->getInstancePath().string() << "\t" << instance->getNbJobs() << "\t" << instance->getNbToSelectJob() << "\t"
               << instance->getNbOfHighSpeedMachines() << "\t" << instance->getNbOfLowSpeedMachines() << "\t" << instance->getHighSpeed() << "\t" << instance->getLowSpeed() << "\t" << "Heuristic"
               << "\t" << time_elapsed.count() << "\t" << time_limits.count() << "\t" << solution->getSumWjUj() << std::endl;
    outputFile.close();
}


