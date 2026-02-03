#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <random>

typedef std::complex<double> c_double;
constexpr c_double I_number(0, 1);
constexpr double pi=3.14159265359;

extern double coarse_time; //Time spent in the coarse grid solver
extern double smooth_time; //Time spent in the smoother
extern double total_time; //Total time spent in the multigrid solver

namespace mass{
    extern double m0;
}

//------------Lattice parameters--------------//
namespace LV {
    //Finest level
    constexpr int Nx=NS; 
    constexpr int Nt = NT; 
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}

//------------Schwarz alternating procedure parameters--------------//
namespace SAPV {
    using namespace LV; 
    constexpr int sap_block_x = 4; //Default values for SAP as a preconditioner for FGMRES (not multigrid)
    constexpr int sap_block_t = 4; 
    //Parameters for GMRES in SAP
    extern int sap_gmres_restart_length; //GMRES restart length for the Schwarz blocks.
    extern int sap_gmres_restarts; //GMRES iterations for the Schwarz blocks
    extern double sap_gmres_tolerance; //GMRES tolerance for the Schwarz blocks 
    extern double sap_tolerance; //tolerance for the SAP method
    extern int sap_blocks_per_proc; //Number of blocks per process for the parallel SAP method
}

//------------Parameters for AMG--------------//
namespace AMGV{
    constexpr int levels = LEVELS; 
    extern int SAP_test_vectors_iterations; //Number of SAP iterations to smooth test vectors
    //Parameters for the coarse level solver
    extern int gmres_restarts_coarse_level; //restart length for GMRES at the coarse level
    extern int gmres_restart_length_coarse_level; //GMRES restart length for the coarse level
    extern double gmres_tol_coarse_level; //GMRES tolerance for the coarse level
    //Parameters for GMRES as a smoother (the default AMG version uses SAP)
    extern int gmres_restarts_smoother; //GMRES iterations for the smoother

    extern int nu1; //Pre-smoothing iterations
    extern int nu2; //Post-smoothing iterations
    extern int Nit; //Number of iterations for improving the interpolator 

    //Outer fgmres solver
    extern int fgmres_k_cycle_restart_length;
    extern int fgmres_k_cycle_restarts;
    extern double fgmres_k_cycle_tol;
    extern int cycle; //Cycling stratey. Cycle = 0 -> V-cycle, = 1 --> K-cycle
}

namespace LevelV{
    extern int BlocksX[AMGV::levels-1];
    extern int BlocksT[AMGV::levels-1];
    extern int Ntest[AMGV::levels-1]; 
    extern int Nagg[AMGV::levels-1]; 
    extern int NBlocks[AMGV::levels-1];
    
    extern int Nsites[AMGV::levels]; //Number of lattice sites at level l
    extern int NxSites[AMGV::levels];
    extern int NtSites[AMGV::levels];
    extern int DOF[AMGV::levels]; //Number of degrees of freedom at each lattice site.
    //On the finest level, DOF = 2 (only spin), on the coarse levels DOF = 2*Ntest
    extern int Colors[AMGV::levels]; //Number of "colors" at each level 
    constexpr int maxLevel = AMGV::levels - 1; //Maximum level id is levels - 1

    extern int SAP_Block_x[AMGV::levels]; //Number of SAP blocks on the x direction 
    extern int SAP_Block_t[AMGV::levels]; //Number of SAP blocks on the t direction
    extern int SAP_elements_x[AMGV::levels]; //Number of SAP blocks on the x direction 
    extern int SAP_elements_t[AMGV::levels]; //Number of SAP blocks on the t direction
    extern int SAP_variables_per_block[AMGV::levels]; //Number of variables in each SAP block

    extern int GMRES_restart_len[AMGV::levels];
    extern int GMRES_restarts[AMGV::levels];
    extern double GMRES_tol[AMGV::levels];
}

//--------------Parameters for outer FGMRES--------------//
namespace FGMRESV {
    extern double fgmres_tolerance; //Tolerance for FGMRES
    extern int fgmres_restart_length; //Restart length for FGMRES
    extern int fgmres_restarts; //Number of restarts for FGMRES
}

namespace CG{
    extern int max_iter;
    extern double tol;
}


//Coordinates vectorization for the lattice points
extern std::vector<std::vector<int>>Coords; 
void Coordinates();

/*
For the boundary conditions
*/
extern std::vector<std::vector<std::vector<int>>>RightPB_l; //Right periodic boundary
extern std::vector<std::vector<std::vector<int>>>LeftPB_l; //Left periodic boundary
extern std::vector<std::vector<std::vector<c_double>>>SignR_l; //Right fermionic boundary
extern std::vector<std::vector<std::vector<c_double>>>SignL_l; //Left fermionic boundary


void save_vec(const std::vector<double>& vec,const std::string& name); //save vector to .txt file 
void read_rhs(std::vector<std::vector<c_double>>& vec,const std::string& name);
void readBinaryRhs(std::vector<std::vector<c_double>>& vec, const std::string& name);
void save_rhs(std::vector<std::vector<c_double>>& vec,const std::string& name);
void random_rhs(std::vector<std::vector<c_double>>& vec,const int seed);


/*
    Print relevant parameters
*/
void printParameters();

void saveParameters(double *Iter, double *dIter, double *exTime, double *dexTime, const int nMeas,const int nconf);


#endif 