cmake_minimum_required(VERSION 3.25.1)

# -----------------------------------------
# Name of the project
set(PROJECT bilevel-scheduling)
# Version of the project
set(VERSION_PROJECT 1.0.0 )
# Language of the project
set(PROJECT_LANGUAGES CXX)
project(${PROJECT} VERSION ${VERSION_PROJECT} LANGUAGES ${PROJECT_LANGUAGES})

# Generation of the list of sources files.

# -----------------------------------------
# files .cpp
# -----------------------------------------
# source for librairy (create new source for new lib)

set(SRCS_BiSchedSolver
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ISolver.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/ColumnGeneration.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Heuristic.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BeamSearch.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/LocalSearch.cpp
)
# ---------------------------------------

# ---------------------------------------
# files .h
# -----------------------------------------
set(HEADERS_BiSchedSolver
        ${CMAKE_CURRENT_SOURCE_DIR}/include/
)
# ---------------------------------------

# files dll under windows
if(WIN32)
   ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
   set(LIB_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/Windows/BiSchedSolver) # define our path for the lib
endif()

# Flags compilation
# compilation parameters
if (UNIX)
    set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-deprecated-copy -Wuninitialized -Wmaybe-uninitialized -pthread ")
    set(SANITIZER_FLAGS
    -fsanitize=address
    -fsanitize=leak
    -fsanitize=undefined
    -fsanitize=signed-integer-overflow
    -fsanitize=integer-divide-by-zero
    -fsanitize=float-divide-by-zero
    -fsanitize=alignment
    -fsanitize=bounds
    -fsanitize=return
    -fsanitize=shift
    -fsanitize=unreachable
    -fsanitize=vla-bound
    -fsanitize=null
    -fsanitize=pointer-compare
    -fsanitize=pointer-subtract
    -fsanitize=bool
    -fsanitize=enum
    -fsanitize=nonnull-attribute
    -fsanitize=float-cast-overflow
    -fno-sanitize-recover=all
    -fno-omit-frame-pointer
    -fsanitize-address-use-after-scope)

    # Convert the list of sanitizer in string to use as flag
    string(JOIN " " SANITIZER_FLAGS_STR ${SANITIZER_FLAGS})

    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -MMD -DNDEBUG -march=native")
    set(CMAKE_CXX_FLAGS_DEBUG "-g ${SANITIZER_FLAGS_STR} -DDEBUG_HEURISTIC -DDEBUG_NEIGHBORHOOD -DDEBUG_CG -DDEBUG_BaB")
    set(LINKER_OPTIONS -flto -Wl,--no-as-needed )
    set(LIB_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/Linux/BiSchedSolver) # define our path for the lib
endif ()

# ---------------------------------------------------------
# Our Bilevel Solver Lib LIBRAIRY
# ---------------------------------------------------------
include_directories(
        ${HEADERS_BiSchedSolver}
)

add_library(BiSchedSolver STATIC ${SRCS_BiSchedSolver}) # add the library to the sources


set_property(TARGET BiSchedSolver PROPERTY CXX_STANDARD 20)

#  Specify the output directory for the library file (.a)
set_target_properties(BiSchedSolver PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${LIB_OUTPUT_DIR}
)

# --------------------------------------------------------
#                   ADD LIBRARIES
# --------------------------------------------------------

set(LIB_PATH "lib")

# ---------------------------------------------------------
# ------- Inclusion of CPLEX --------
# CPLEX Lib for solver.
if(WIN32)
    set(CPLEXlib_DIR ${LIB_PATH}/Windows/CPLEX)
else()
    set(CPLEXlib_DIR ${LIB_PATH}/Linux/CPLEX)
endif()

find_package(CPLEXlib REQUIRED)

if(CPLEXlib_FOUND)
    message(STATUS "lib CPLEX found")
else()
    message(STATUS "lib CPLEX not found")
endif()
# ---------------------------------------
include_directories(${CPLEX_INCLUDE_DIRS})

target_link_libraries(BiSchedSolver
        ${CPLEX_LIBRARIES}
        ${CMAKE_DL_LIBS}
)

target_include_directories(BiSchedSolver PUBLIC ${CPLEX_INCLUDE_DIRS}) # lib header
message(STATUS "CPLEX_LIBRARIES (debug) = ${CPLEX_LIBRARIES}")
message(STATUS "CPLEX_LIBRARIES_RELEASE (optimized) = ${CPLEX_LIBRARIES_RELEASE}")
# ---------------------------------------
# ------- Inclusion of Eigen --------
# ------ Matrix library --------

message(STATUS "******************************")
message(STATUS "**** ADDING Eigen library ****")
message(STATUS "******************************")

add_subdirectory(${LIB_PATH}/all/eigen-3.4.0)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)

if(Eigen3_FOUND)
    message(STATUS "lib Eigen3.4 found")
else()
    message(STATUS "lib Eigen3.4 not found")
endif()
# ---------------------------------------

# ---------------------------------------
# ------- Inclusion of Boost Dynamic bitset --------
# ------ Dynamic bitset --------

message(STATUS "*************************************")
message(STATUS "**** ADDING Boost Dynamic bitset ****")
message(STATUS "*************************************")

include_directories(
        ${LIB_PATH}/all/boost_1_88_0
        ${LIB_PATH}/all/boost_1_88_0/boost
)

target_include_directories(BiSchedSolver PUBLIC
        ${LIB_PATH}/all/boost_1_88_0
        ${LIB_PATH}/all/boost_1_88_0/boost
)

# ---------------------------------------

# -----------------------------------------------------------
# ------- Inclusion of JSON library --------

message(STATUS "**************************************")
message(STATUS "**** ADDING external JSON library ****")
message(STATUS "**************************************")
add_subdirectory(${LIB_PATH}/all/json-3.11.3)


# specify target link of BiSchedSolver

target_link_libraries(BiSchedSolver
        Eigen3::Eigen
        nlohmann_json::nlohmann_json
)

# ---------------------------------------

# -----------------------------------------------------------
# ------- Inclusion of Libtorch library --------

message(STATUS "**************************************")
message(STATUS "**** ADDING external Libtorch library ****")
message(STATUS "**************************************")
# Torch Lib for predictor.
if(WIN32)
    set(Torch_DIR ${LIB_PATH}/Windows/libtorch/share/cmake/Torch)
else()
    set(Torch_DIR ${LIB_PATH}/Linux/libtorch/share/cmake/Torch)
endif()

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

# specify target link of BiSchedSolver
target_link_libraries(BiSchedSolver
        "${TORCH_LIBRARIES}"
)

# ---------------------------------------
