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
// Created by schau on 9/23/24.
//

#ifndef BILEVEL_SCHEDULING_DSU_H
#define BILEVEL_SCHEDULING_DSU_H

#include <algorithm>
#include <iostream>

struct JumpPoint {
    unsigned int nextT;
    int previousT;
    unsigned int difference;

    inline explicit JumpPoint(unsigned int nextT, unsigned int previousT, unsigned int difference) : nextT(nextT), previousT(previousT), difference(difference) {};
};

// DSU (Disjoint Set Union) datastructure
struct DSU {
    // Array of parent nodes
    unsigned int *parents;
    // Array of sizes for each set
    unsigned int *sizes;
    // Array of minimum values for each set
    unsigned int *minSet;
    // Number of sets in the forest
    unsigned int n;

    inline explicit DSU() {}

    // Constructor, initializes the DSU with a specified number of sets
    inline explicit DSU(unsigned int n) : n(n) {
        parents = new unsigned int[n];
        sizes = new unsigned int[n];
        minSet = new unsigned int[n];
        for (unsigned int i = 0; i < n; ++i) {
            // Initialize parent, size, and minimum value for each set
            parents[i] = i;
            minSet[i] = i;
            sizes[i] = 1;
        }
    }

    // Resets the DSU to a new number of sets
    inline void reset(unsigned int new_n) {
        n = new_n;
        for (unsigned int i = 0; i < n; ++i) {
            parents[i] = i;
            minSet[i] = i;
            sizes[i] = 1;
        }
    }

    // Unites two sets in the DSU
    inline void union_sets(unsigned int u, unsigned int v) {
        u = find_parent(u);
        v = find_parent(v);
        if (sizes[u] > sizes[v]) std::swap(u, v);
        parents[u] = v;
        minSet[v] = std::min(minSet[u], minSet[v]);
        sizes[v] += sizes[u];
    }

    // Finds the parent of a set in the DSU
    inline int find_parent(unsigned int u) {
        if (u >= n) u = n - 1;
        while (parents[u] != u) {
            parents[u] = parents[parents[u]];
            u = parents[u];
        }
        return u;
    }

    // Prints the minimum value of each set in the DSU
    inline void printM_t() {
        std::cout << "mt: [";
        for (unsigned int i = 0; i < n; ++i) {
            auto t = find_parent(i);
            std::cout << minSet[t] << " ";
        }
        std::cout << "]" << std::endl;
    }

    // Destructor, frees memory allocated by the DSU
    inline ~DSU() {
        delete[] parents;
        delete[] sizes;
        delete[] minSet;

    }
};


#endif //BILEVEL_SCHEDULING_DSU_H
