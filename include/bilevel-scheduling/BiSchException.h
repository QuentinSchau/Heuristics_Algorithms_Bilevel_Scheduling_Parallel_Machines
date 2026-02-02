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
// Created by schau on 3/27/24.
//

#ifndef BILEVEL_SCHEDULING_BISCHEXCEPTION_H
#define BILEVEL_SCHEDULING_BISCHEXCEPTION_H

#include <exception>
#include <string>

// Define a new exception class that inherits from
// std::exception
class BiSchException : public std::exception {
private:
    std::string message;

public:
    explicit BiSchException(const char *msg) : message(msg) {}

    explicit BiSchException(std::string &msg) : message(msg) {}

    // Override the what() method to return our message
    virtual const char *what() const throw() {
        auto error = message.c_str();
        return error;
    }
};

class BiSchTimeOutException : public std::exception {
private:
    std::string message;

public:
    explicit BiSchTimeOutException() : message() {}

    explicit BiSchTimeOutException(const char *msg) : message(msg) {}

    explicit BiSchTimeOutException(std::string &msg) : message(msg) {}

    // Override the what() method to return our message
    virtual const char *what() const throw() {
        auto error = message.c_str();
        return error;
    }
};

#endif //BILEVEL_SCHEDULING_BISCHEXCEPTION_H
