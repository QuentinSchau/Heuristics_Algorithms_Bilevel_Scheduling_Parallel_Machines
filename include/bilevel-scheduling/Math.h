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

#ifndef BILEVEL_SCHEDULING_MATH_H
#define BILEVEL_SCHEDULING_MATH_H

#include <numeric>
#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <stack>
#include <deque>
#include "Job.h"


template<typename type1, typename type2>
static constexpr auto isSmallerOrEqual(const type1 a, const type2 b, const double epsilon = 1e-9) {
    double da = static_cast<double>(a);
    double db = static_cast<double>(b);
    return (da - db <= -epsilon) || (std::fabs(da - db) <= epsilon);
}

template<typename type1, typename type2>
static constexpr auto isEqual(const type1 a, const type2 b, const double epsilon = 1e-9) {
    double da = static_cast<double>(a);
    double db = static_cast<double>(b);
    return (std::fabs(da - db) <= epsilon);
}

template<typename type1, typename type2>
static constexpr auto isSmaller(const type1 a, const type2 b, const double epsilon = 1e-9) {
    double da = static_cast<double>(a);
    double db = static_cast<double>(b);
    return da - db <= -epsilon;
}


template<typename type>
auto computeCartProdRec(std::vector<std::vector<type>> &outputs, std::vector<type> &currentProduct, typename std::vector<std::vector<type>>::const_iterator me
                        , typename std::vector<std::vector<type>>::const_iterator end) {
    if (me == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. Add the current result (rvi)
        // to the total set of results (rvvvi).
        outputs.push_back(currentProduct);
        return;
    }

    // need an easy name for my vector-of-ints
    const std::vector<type> &mevi = *me;
    for (typename std::vector<type>::const_iterator it = mevi.begin(); it != mevi.end(); ++it) {
        // final rvi will look like "a, b, c, ME, d, e, f"
        // At the moment, rvi already has "a, b, c"
        currentProduct.push_back(*it);  // add ME
        computeCartProdRec(outputs, currentProduct, me + 1, end); //add "d, e, f"
        currentProduct.pop_back(); // clean ME off for next round
    }
}

template<typename type>
auto cartesianProduct(const std::vector<std::vector<type>> &lists) {

    //compute the number of product in order to reserve memory space
    size_t numberOfProduct = std::accumulate(lists.begin(), lists.end(), 0, [&](auto nbProd, auto newSet) -> size_t {
        return nbProd + newSet.size();
    });
    std::vector<std::vector<type>> result;
    result.reserve(numberOfProduct);
    std::vector<type> temp;
    computeCartProdRec(result, temp, lists.cbegin(), lists.cend());
    return result;
}

/**
 * Method that compute the factorial
 * @param n An integer
 * @return n!
 */
unsigned int factorial(unsigned int n);

/**
 * GenerateConcatenatedPermutations - Method that generates concatenated permutations for a list of lists.
 * This method takes two parameters: the first list is for the output, and the second list is a parameter.
 * It computes all permutations for each element in the given list of lists and concatenates them.
 * For example, in the first step, we take the first element of 'lists' named 'lists[1]'. We compute all
 * permutations of 'lists[1]' denoted as Plists[1][1], Plists[1][2], ..., Plists[1][p], where p = (lists[1].size())!
 * - the number of permutations of 'lists[1]'. In the second step, we take the second element of 'lists' named
 * 'lists[2]'. We compute all permutations of 'lists[2]' denoted as Plists[2][1], Plists[2][2], ..., Plists[2][p']
 * where p' = (lists[2].size())! - the number of permutations of 'lists[2]'. Then, we concatenate each permutation
 * i.e., Plists[1][1] | Plists[2][1], Plists[1][1] | Plists[2][2], ..., Plists[1][1] | Plists[2][p'], Plists[1][2] | Plists[2][1], ...
 * We repeat this process for all elements in 'lists', and in the end, we get a new list that contains all permutations
 * of each element from 'lists' concatenated.
 * @tparam type The type of each list of lists
 * @param listsOutput The lists where the outputs will be added
 * @param lists The lists from which we compute all permutations
 */
template<typename type>
void
generateConcatenatedPermutations(std::vector<std::vector<type>> &listsOutput, std::vector<std::vector<type>> &lists) {

    //compute the number of permutation for reserve memory for the list of permutations
    size_t nbPermutation = 1;
    for (auto &listPermutation: lists) {
        nbPermutation *= factorial(listPermutation.size());
    }
    listsOutput.reserve(nbPermutation);

    //create a dequeu for create all permutation of list of location
    std::deque<std::vector<type>> dequeListLocation(1);
    // loop over each list of locations
    for (auto &listLocation: lists) {
        // loop over the permutation from the dequeue. We take an element from the dequeue.
        // We compute each permutation from the list of location, and append this permutation to the end of element from dequeue.
        size_t currentSizeOfDequeue = dequeListLocation.size();
        for (size_t i = 0; i < currentSizeOfDequeue; ++i) {
            std::vector<type> currentPermutation = dequeListLocation.back();
            dequeListLocation.pop_back();
            do {
                std::vector<type> temp;
                temp.reserve(currentPermutation.size());
                temp.insert(temp.end(), currentPermutation.begin(), currentPermutation.end());
                temp.insert(temp.end(), listLocation.begin(), listLocation.end());
                dequeListLocation.push_front(std::move(temp));
            } while (std::next_permutation(listLocation.begin(), listLocation.end()));
        }

    }

    for (auto &listLocation: dequeListLocation) {
        listsOutput.push_back(std::move(listLocation));
    }

}

struct HashFunction {
    std::size_t operator()(std::vector<int> const &vec) const {
        std::size_t seed = vec.size();
        for (auto x: vec) {
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = ((x >> 16) ^ x) * 0x45d9f3b;
            x = (x >> 16) ^ x;
            seed ^= x + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

unsigned int nChoosek(unsigned int n, unsigned int k);

/**
 * Compute the elegantPair bijection between N^2 and N
 * @tparam T Must be integer, i.e., unsigned int, unsigned long long, etc.
 * @param x first element
 * @param y second element
 * @return the bijection of the couple in N
 * @throws std::invalid_argument if T is not an integer type
 */
template<typename T>
T elegantPair(T x, T y) {
    // Check if T is an integer type
    static_assert(std::is_integral_v<T>, "T must be an integer type");
    if (x < y) return y * y + x;
    else return x * x + x + y;
};

/**
 * Compute the elegantPair bijection between N^k and N
 * @tparam T Must be integer, i.e., unsigned int, unsigned long long, etc.
 * @tparam Args Must be integer, i.e., unsigned int, unsigned long long, etc. Must be same type as T
 * @param x1 first element
 * @param x2 second element
 * @param args the k-2 element of the tuple
 * @return the bijection of the tuple in N
 * @throws std::invalid_argument if T and Args is not an integer type
 */
template<typename T, typename... Args>
T elegantPair(T x1, T x2, Args... args) {
    // Check if T is an integer type
    static_assert(std::is_integral_v<T>, "T must be an integer type");
    // Base case: If there are no more arguments, compute the bijection for x1 and x2
    if constexpr (sizeof...(args) == 0) {
        return elegantPair(x1, x2);
    } else {
        // Recursive case: Compute the bijection for each group of two elements
        auto result = elegantPair(elegantPair(x1, x2), args...);
        // Return the final result
        return result;
    }
};

/**
 * Compute the list of integer that correspond to the bijection elegantPair
 * @tparam T Must be integer, i.e., unsigned int, unsigned long long, etc.
 * @param z the value in N
 * @param k the dimension of N^k from where the bijection was computed
 * @return the vector of size k corresponding to z
 * @throws std::invalid_argument if T and Args is not an integer type
 */
template<typename T>
std::vector<T> elegantUnpair(T z, unsigned int k) {
    // Check if T is an integer type
    static_assert(std::is_integral_v<T>, "T must be an integer type");
    std::vector<T> output;
    output.reserve(k);
    for (unsigned i = 1; i < k - 1; ++i) {
        T a = static_cast<T>(std::floor(std::sqrt(z)));
        T b = a * a;
        if (z - b < a) {
            z = z - b;
            output.insert(output.begin(), a);
        } else {
            output.insert(output.begin(), z - b - a);
            z = a;
        }
    }
    T a = static_cast<T>(std::floor(std::sqrt(z)));
    T b = a * a;
    if (z - b < a) {
        output.insert(output.begin(), a);
        output.insert(output.begin(), z - b);
    } else {
        output.insert(output.begin(), z - b - a);
        output.insert(output.begin(), a);
    }
    return output;
}

/**
 * Compute the pair of integer that correspond to the bijection elegantPair
 * @tparam T Must be integer, i.e., unsigned int, unsigned long long, etc.
 * @param z the value in N
 * @return the pair corresponding to z
 * @throws std::invalid_argument if T and Args is not an integer type
 */
template<typename T>
std::pair<T, T> elegantUnpair(T z) {
    // Check if T is an integer type
    static_assert(std::is_integral_v<T>, "T must be an integer type");
    std::pair<T, T> output;
    T a = static_cast<T>(std::floor(std::sqrt(z)));
    T b = a * a;
    if (z - b < a) {
        output.first = z - b;
        output.second = a;
    } else {
        output.first = a;
        output.second = z - b - a;
    }
    return output;
}

// Structure to hold state information for the stack
template<typename T>
struct CombinationState {
    unsigned int currentIndex;
    unsigned int remainingElements;
    std::vector<T> currentCombination;
};

/**
 * Function to find all distinct combinations of length `combinationLength`. For instance, if arr = [1,1,2,3,3] and
 * combinationLength = 2, then we will get {[1,1],[1,2],[1,3],[2,3],[3,3]}
 * @param arr the array of element from where we will select subsets of length combinationLength
 * @param combinationLength the length of all subsets
 * @param combinations the set of all combination.
 */
template<typename T>
void findCombinations(const std::vector<T> &arr, unsigned int combinationLength, std::set<std::vector<T>> &combinations) {
    if (arr.empty() || combinationLength == 0) return;

    unsigned int nbDistinctElements = arr.size();
    std::string bitmask(combinationLength, 1); // K leading 1's
    bitmask.resize(nbDistinctElements, 0); // N-K trailing 0's

    do {
        // create the combination
        std::vector<T> combination;
        combination.reserve(combinationLength);
        for (unsigned int i = 0; i < nbDistinctElements; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) combination.push_back(arr[i]);
        }
        combinations.emplace(std::move(combination));
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

}

template<typename T>
struct std::hash<std::pair<T,T>>{
    static_assert(std::is_integral_v<T>, "Hash is not implemented for non-integer types. Please use an integer type (e.g. int, long, etc.) or consider std::hash or boost::hash");
    std::size_t operator()(const std::pair<T,T>& p) const noexcept{
        return elegantPair(p.first, p.second);
    }
};


#endif //BILEVEL_SCHEDULING_MATH_H
