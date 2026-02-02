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

#include "Instance.h"
#include "Parser.h"
#include "Solution.h"
#include "ColumnGeneration.h"
#include "BeamSearch.h"
#include "LocalSearch.h"
#include <iostream>
#include <nlohmann/json.hpp>


int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "ERROR: You need at least one argument." << std::endl;
        return -1;
    }
    char **pargv = argv + 1;
    try {
        for (; *pargv != argv[argc]; pargv++) {
            auto pathFile = std::filesystem::path(std::string(*pargv));
            if (! std::filesystem::exists(pathFile)) {
                throw BiSchException(std::string("The file configuration use do not exist. Path: ").append(pathFile.string()));
            }
            std::ifstream f(*pargv);
            // parse json input config file
            nlohmann::json config = nlohmann::json::parse(f);

            // create instance parser
            Parser parser = Parser();

            /********************************/
            /*      GENERATE INSTANCES      */
            /********************************/

            if (config.contains("generate")) {
                // generate them with Parser
                parser.generateInstance(config["generate"]);
            }

            /*****************************/
            /*      SOLVE INSTANCES      */
            /*****************************/

            if (config.contains("solve")) {
                // set verbose mode
                char verbose = config["solve"].contains("verbose") && config["solve"]["verbose"].is_number_unsigned()
                               ? config["solve"]["verbose"].template get<char>() : 0;

                if (config["solve"].contains("methods")) {
                    // for each method
                    for (auto &method: config["solve"]["methods"]) {
                        // set output path
                        std::string outputPath;
                        if (config["solve"].contains("output") && config["solve"]["output"].is_string())
                            outputPath = std::filesystem::path(config["solve"]["output"].template get<std::string>());
                        else {
                            outputPath = std::filesystem::current_path().parent_path().parent_path().string() +
                                         "/instances/";
                        }

                        // keep the path without the extension and add the name method;
                        outputPath.append("results" + method["name"].template get<std::string>() + ".csv");
                        if (verbose >= 2) std::cout << "Save results in the path : " << outputPath << std::endl;
                        std::ofstream outputFileStream;

                        if (method.contains("instances")) {
                            // loop over each instances
                            for (auto &instance: method["instances"]) {
                                Instance newInstance;
                                if (instance.contains("path")) {
                                    if (instance["path"].is_string()) {
                                        std::string path = instance["path"];
                                        if (verbose >= 2) std::cout << "Parsing instance : " << path << std::endl;
                                        newInstance = parser.readFromFile(path);
                                        newInstance.sort_by_SPT();
                                        if (method["name"] == "CG") {
                                            // sort instance by SPT
                                            newInstance.sort_by_SPT();

                                            if (verbose >= 1)
                                                std::cout << "-------- COLUMN GENERATION ---------" << std::endl
                                                          << std::endl;
                                            ColumnGeneration optSolCG = ColumnGeneration(&newInstance, method);
                                            Node node = Node(&newInstance);
                                            try{
                                            optSolCG.solve(node);
                                            }catch (const BiSchTimeOutException &e) {}//do nothing

                                            optSolCG.printOutput(outputPath, outputFileStream);
                                            double lowerBound = optSolCG.getSumWjUj();
                                            if (verbose >= 2) {
                                                if (!optSolCG.getSolution()->empty()) {
                                                    optSolCG.getSolution()->evaluate();
                                                    if (verbose >= 3)
                                                        std::cout << (*optSolCG.getSolution()) << std::endl;
                                                    std::string isFeasible = optSolCG.getSolution()->feasible(
                                                            &newInstance) ? "true" : "false";
                                                    std::cout << "Is feasible : " << isFeasible << std::endl;
                                                } else {
                                                    std::cout << "Not integer solution " << std::endl;
                                                }
                                            }
                                            std::cout << "Lower bound value given by CG: "<< lowerBound<< std::endl;
                                        }else if (method["name"] == "LocalSearch") {
                                            // sort instance by SPT
                                            newInstance.sort_by_SPT();
                                            if (verbose >= 1)
                                                std::cout << "-------- LOCAL SEARCH ---------" << std::endl << std::endl;
                                            LocalSearch solveLS = LocalSearch(&newInstance, method);
                                            solveLS.run("Local Search");
                                            solveLS.printOutput(outputPath, outputFileStream);
                                            if (not solveLS.isGenDatabase()) {
                                                bool isFeasible = solveLS.getSolution()->feasible(&newInstance);
                                                if (verbose >= 2) {
                                                    if (verbose >= 3)
                                                        std::cout << (*solveLS.getSolution()) << std::endl;
                                                    std::string feasibility = isFeasible ? "true" : "false";
                                                    std::cout << "Is feasible : " << feasibility << std::endl;
                                                }
                                                if (not isFeasible){
                                                    std::cerr << *solveLS.getSolution() << std::endl;
                                                    throw BiSchException("Heuristic don't construct a feasible solution");
                                                 }
                                            }
                                        } else if (method["name"] == "BeamSearch") {
                                            // sort instance by SPT
                                            newInstance.sort_by_SPT();
                                            if (verbose >= 1)
                                                std::cout << "-------- RECOVERING BEAM SEARCH ---------" << std::endl
                                                          << std::endl;
                                            BeamSearch solveBeamSearch = BeamSearch(&newInstance, method);
                                            solveBeamSearch.run("Recovering Beam Search");
                                            bool isFeasible = solveBeamSearch.getSolution()->feasible(&newInstance);
                                            if (verbose >= 2) {
                                                if (verbose >= 3)
                                                    std::cout << *solveBeamSearch.getSolution() << std::endl;
                                                std::string feasibility = isFeasible ? "true" : "false";
                                                std::cout << "Is feasible : " << feasibility << std::endl;
                                            }
                                            if (not isFeasible) {
                                                std::cerr << *solveBeamSearch.getSolution() << " firstUB:" << solveBeamSearch.getFirstUb() << " globalUB:" << solveBeamSearch.getGlobalUb() << " sumPj:" << newInstance.getSumPj() << std::endl;
                                                solveBeamSearch.getSolution()->explainInfeasibility(&newInstance);
                                                throw BiSchException("BeamSearch don't construct a feasible solution");
                                            }
                                            solveBeamSearch.printOutput(outputPath, outputFileStream);
                                        } else std::cout << method["name"] << " IS NOT IMPLEMENTED " << std::endl;
                                    } else throw std::invalid_argument(R"(The instance path is not a string)");
                                } else throw std::invalid_argument(R"(The instance don't have attribute "path")");
                            }
                        }
                        outputFileStream.close();
                    }
                } else throw std::invalid_argument(R"(The config don't have attribute "methods")");
            }
        }
    }catch (const IloException &e) {
        std::cerr << "Error with "<< *pargv << std::endl << "Error: " << e;
        return -1;
    }catch (const BiSchException &e) {
        std::cerr << "Error with "<< *pargv << std::endl << "Error: " << e.what();
        return -1;
    }catch (const std::exception &e) {
        std::cerr << "Error with "<< *pargv << std::endl << "Error: " << e.what();
        return -1;
    }
    catch (...) {
        std::cerr << "Error with "<< *pargv << std::endl << "Error: unknown" ;
    }

    return 0;
}



