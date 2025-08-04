#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <complex>
#include <vector>

typedef std::complex<double> c_double;
constexpr double pi=3.14159265359;

extern double coarse_time; //Time spent in the coarse grid solver
extern double smooth_time; //Time spent in the smoother
extern double SAP_time;


//------------Lattice parameters--------------//
namespace LV {
    //Finest level
    constexpr int Nx=NS; 
    constexpr int Nt = NT; 
    constexpr int block_x = BLOCK_X; //Read this from a parameters file ... 
    constexpr int block_t = BLOCK_T; //read from parameters file ...
    constexpr int Nblocks = block_x * block_t;
    constexpr int x_elements = Nx/block_x; //Number of elements in the x direction
    constexpr int t_elements = Nt/block_t; //Number of elements in the t direction
    constexpr int lattice_sites_per_block = x_elements * t_elements;
    constexpr int variables_per_agg = lattice_sites_per_block;
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
}

//------------Schwarz alternating procedure parameters--------------//
namespace SAPV {
    using namespace LV; 
    constexpr int sap_block_x = SAP_BLOCK_X; 
    constexpr int sap_block_t = SAP_BLOCK_T; 
    constexpr int sap_x_elements = Nx/sap_block_x; //Number of lattice points in the x direction (without the spin index)
    constexpr int sap_t_elements = Nt/sap_block_t; //Number of lattice points in the t direction (without the spin index)
    constexpr int N_sap_blocks = sap_block_x * sap_block_t; //Number of Schwarz blocks
    constexpr int sap_lattice_sites_per_block = sap_x_elements * sap_t_elements; //Number of lattice points in the block
    constexpr int sap_variables_per_block = 2 * sap_lattice_sites_per_block; //Number of variables in the block
    constexpr int sap_coloring_blocks = N_sap_blocks/2; //Number of red or black blocks 
    extern bool schwarz_blocks; //True if the Schwarz blocks are initialized
    extern int sap_gmres_restart_length; //GMRES restart length for the Schwarz blocks.
    extern int sap_gmres_restarts; //GMRES iterations for the Schwarz blocks
    extern double sap_gmres_tolerance; //GMRES tolerance for the Schwarz blocks 
    extern double sap_tolerance; //tolerance for the SAP method
    extern int sap_blocks_per_proc; //Number of blocks per process for the parallel SAP method
}

//------------Parameters for AMG--------------//
namespace AMGV{
    constexpr int levels = LEVELS; 
    constexpr int Ntest = NTEST; //Number of test vectors
    constexpr int Nagg = 2*LV::block_t*LV::block_x; //Number of aggregates (one spin per aggregate)
    extern bool aggregates_initialized; //True if the aggregates are initialized
    extern int SAP_test_vectors_iterations; //Number of SAP iterations to smooth test vectors
    //Parameters for the coarse level solver
    extern int gmres_restarts_coarse_level; //restart length for GMRES at the coarse level
    extern int gmres_restart_length_coarse_level; //GMRES restart length for the coarse level
    extern double gmres_tol_coarse_level; //GMRES tolerance for the coarse level
    //Parameters for GMRES as a smoother (the default AMG version uses SAP)
    extern int gmres_restarts_smoother; //GMRES iterations for the smoother
    //Paramaters for bi-cgstab as a coarse solver
    extern int bi_cgstab_Dc_iterations;
    extern double bi_cgstab_Dc_iterations_tol; //Tolerance for the bi-cgstab method
    extern int nu1; //Pre-smoothing iterations
    extern int nu2; //Post-smoothing iterations
    extern int Nit; //Number of iterations for improving the interpolator 
    extern bool SetUpDone; 
}

namespace LBlocks{
    extern int BlocksX[AMGV::levels-1];
    extern int BlocksT[AMGV::levels-1];
    extern int Ntest[AMGV::levels-1]; 
    extern int Nagg[AMGV::levels-1]; 
}

//--------------Parameters for FGMRES--------------//
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

//Coordinates vectorization for the aggregates 
//Intialized in Aggregates()
extern std::vector<int> XCoord; 
extern std::vector<int> TCoord;
extern std::vector<int> SCoord;
extern std::vector<std::vector<int>> Agg; 
//Aggregates[number_of_aggregate][vectorized_coordinate of the lattice point]
//This vectorization also takes into account the spin index


//--Coordinates of the neighbors to avoid recomputing them each time the operator D is called--//
//Check dirac_operator.h for the definition of RightPB and LeftPB
extern std::vector<std::vector<int>>RightPB; //Right periodic boundary
extern std::vector<std::vector<int>>LeftPB; //Left periodic boundary
extern std::vector<std::vector<c_double>>SignR; //Right fermionic boundary
extern std::vector<std::vector<c_double>>SignL; //Left fermionic boundary

//Indexes for the SAP blocks
extern std::vector<std::vector<int>> SAP_Blocks; //SAP_Blocks[number_of_block][vectorized_coordinate of the lattice point]
//The vectorization does not take into account the spin index, since both spin indices are in the same block.
extern std::vector<int> SAP_RedBlocks; //Block index for the red blocks
extern std::vector<int> SAP_BlackBlocks; //Block index for the black blocks

extern std::vector<std::vector<c_double>>D_TEMP;

extern std::vector<std::vector<int>> LatticeBlocks;
extern int LeftPB_blocks[LV::Nblocks][2];
extern int RightPB_blocks[LV::Nblocks][2];

//Build lattice blocks 
void MakeBlocks();

void CheckBlocks(); //Check that Nx/block_x and Nt/block_t are integers, the same for Schwarz blocks
void CheckAggregates(); //Check that the aggregates are initialized and have the correct size

#endif 