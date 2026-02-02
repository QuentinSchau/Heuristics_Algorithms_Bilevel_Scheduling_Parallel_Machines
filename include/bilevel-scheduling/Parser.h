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

#ifndef BILEVEL_SCHEDULING_PARSER_H
#define BILEVEL_SCHEDULING_PARSER_H

#include "Instance.h"
#include <string>
#include <nlohmann/json.hpp>

class Parser {

public:

    /************************/
    /*      CONSTRUCTOR     */
    /************************/

    Parser();

    /********************/
    /*      METHODS     */
    /********************/


    /**
     * Constructor method that parses a file and generates an instance. It constructs an Instance object.
     * @param filePath The path of the file to parse
     * @return A new instance constructed from the file
     */
    Instance readFromFile(std::string &filePath) const;

    /**
     * Method that serializes an instance into a file specified by the attribute Instance::instancePath.
     * @param instance The instance to be serialized
     */
    void serializeInstance(Instance &instance);


    /**
     * Method that constructs an instance from a JSON object.
     * @param jsonObject The JSON object to be parsed and used to construct the instance
     */
    void generateInstance(nlohmann::json &object);
};


#endif //BILEVEL_SCHEDULING_PARSER_H
