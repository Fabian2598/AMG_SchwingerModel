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
    src/amg.cpp                 #Multilevel method
    src/main.cpp                #Main program
    src/level.cpp               #Class implementation for one level
    src/fgmres.cpp              #Flexible GMRES solver
    src/tests.cpp               #Class implementation for testing solvers
)

#--------------------------Lattice dimensions---------------------------#
set(NS "128") #Space
set(NT "128") #Time
#-----------------------------------------------------------------------#

#--------------------------------Levels---------------------------------#
set(LEVELS "4") #Number of AMG levels (Including the finest and coarsest level)
#-----------------------------------------------------------------------#

configure_file(${CMAKE_SOURCE_DIR}/include/config.h.in ${CMAKE_SOURCE_DIR}/include/config.h)

add_executable("AMG_${NS}x${NT}_l${LEVELS}" ${SRC_FILES})

target_include_directories("AMG_${NS}x${NT}_l${LEVELS}"
        PUBLIC 
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_options("AMG_${NS}x${NT}_l${LEVELS}" PRIVATE
            #-Wall
			-O3
            #-g
			#-fopenmp
)

find_package(MPI REQUIRED)

target_link_libraries("AMG_${NS}x${NT}_l${LEVELS}" PUBLIC MPI::MPI_CXX)

#find_package(BLAS)
#include_directories(" /usr/include/cblas/")
#target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX ${BLAS_LIBRARIES})

#find_package(OpenMP REQUIRED)
#target_link_libraries("SM_${NS}x${NT}" PUBLIC MPI::MPI_CXX OpenMP::OpenMP_CXX)

#For windows with MinGW, use the following command to generate the makefile:
#cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER="C:/msys64/ucrt64/bin/g++.exe" -DCMAKE_C_COMPILER="C:/msys64/ucrt64/bin/gcc.exe" ../
