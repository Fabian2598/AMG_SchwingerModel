cmake_minimum_required(VERSION 3.26.5)
set(CMAKE_CXX_STANDARD 20) 
set(CMAKE_CXX_STANDARD_REQUIRED True)
project(Multigrid VERSION 1.0)

#Source files
set(SRC_FILES 
    src/variables.cpp           #Relevant parameters for the method
    src/gauge_conf.cpp          #Gauge configuration class
    src/dirac_operator.cpp      #Implementation of the Dirac operator
    src/bi_cgstab.cpp           #Bi-CGStab solver for comparison
    src/conjugate_gradient.cpp  #Conjugate gradient for solving normal equations
    src/sap.cpp                 #Schwarz alternating method
    src/amg.cpp                 #Two-grid implementation 
    src/main.cpp                #Main program
    src/level.cpp
    src/fgmres.cpp              #Flexible GMRES solver
)

#--------------------------Lattice dimensions---------------------------#
set(NS "64") #Space
set(NT "64") #Time
#-----------------------------------------------------------------------#

#--------------------------------Levels---------------------------------#
set(LEVELS "3") #Number of AMG levels (Including the finest level)
#-----------------------------------------------------------------------#

#------------------Lattice blocks used for aggregation------------------#
# A_j_0 = L_j x {0}, A_j_1 = L_j x {1} (Volume times Spin component)
set(BLOCK_X "4")
set(BLOCK_T "4")
#N_i / BLOCK_i has to be an integer number
#The number of aggregates is 2 * BlockX * BlockT
#-----------------------------------------------------------------------#

#-------------Lattice blocks for Schwarz alternating method-------------# 
#The blocking does not have to the same as for the aggregation
set(SAP_BLOCK_X "4")
set(SAP_BLOCK_T "4")
#N_i / SAP_BLOCK_i has to be an integer number
#One block contains 2 * N_x / SAP_BLOCK_x * N_t / SAP_BLOCK_t variables 
#For parallelization we fix one block per core i.e. 
#SAP_BLOCK_X * SAP_BLOCK_T / 2 cores
#-----------------------------------------------------------------------#

#Number of test vectors
set(NTEST "10") 
#The coarse grid matrix dimension is NTEST * 2 * BLOCK_X * BLOCK_T

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)

add_executable("SM_${NS}x${NT}" ${SRC_FILES})

target_include_directories("SM_${NS}x${NT}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_options("SM_${NS}x${NT}" PRIVATE
            #-Wall
			-O3
            #-g
			#-fopenmp
)

find_package(MPI REQUIRED)

target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX)

#find_package(BLAS)
#include_directories(" /usr/include/cblas/")
#target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX ${BLAS_LIBRARIES})

#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)

#For windows with MinGW, use the following command to generate the makefile:
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../


