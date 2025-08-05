#ifndef LEVEL_H
#define LEVEL_H

#include "variables.h"
#include "operator_overloads.h"

/*
    One level of the AMG method
*/
class Level {
public:
    Level(const int& level) : level(level) {
        test_vectors = std::vector<spinor>(LevelV::Ntest[level],
        spinor( LevelV::Nsites[level], c_vector (LevelV::DOF[level],0))); 
	    interpolator_columns = std::vector<spinor>(LevelV::Ntest[level],
        spinor( LevelV::Nsites[level], c_vector (LevelV::DOF[level],0))); 
	    v_chopped = std::vector<spinor>(LevelV::Ntest[level],
        spinor( LevelV::Nsites[level], c_vector (LevelV::DOF[level],0))); 

        x_elements = LevelV::NxSites[level] / LevelV::BlocksX[level]; 
        t_elements = LevelV::NtSites[level] / LevelV::BlocksT[level]; 
        sites_per_block = x_elements * t_elements;
        Agg = new int[2*LevelV::NBlocks[level] * sites_per_block];
        LatticeBlocks = new int[LevelV::NBlocks[level] * sites_per_block];
    };

    ~Level() {
        delete[] Agg;
        delete[] LatticeBlocks;
    }

    std::vector<spinor> test_vectors; //[Ntest][Nsites][degrees of freedom per site]
    std::vector<spinor> interpolator_columns;
    std::vector<spinor> v_chopped;
private:
    int level; 
    int x_elements = LevelV::NxSites[level] / LevelV::BlocksX[level], t_elements = LevelV::NtSites[level] / LevelV::BlocksT[level]; //x and t elements of each lattice block
    int sites_per_block = x_elements * t_elements;

    
    //Coarse gauge links --> Used to assemble the operator for the level l+1
	c_double A_coeff[LV::Nblocks][2][2][AMGV::Ntest][AMGV::Ntest];    //[A(x)]^{alf,bet}_{p,s} --> A_coeff[x][alf][bet][p][s] 
	c_double B_coeff[LV::Nblocks][2][2][AMGV::Ntest][AMGV::Ntest][2]; //[B_mu(x)]^{alf,bet}_{p,s} --> B_coeff[x][alf][bet][p][s][mu]
	c_double C_coeff[LV::Nblocks][2][2][AMGV::Ntest][AMGV::Ntest][2]; //[C_mu(x)]^{alf,bet}_{p,s}  --> C_coeff[x][alf][bet][p][s][mu]
    //Level next_level;

    //Agg[i][j] is accessed as Agg[i * sites_per_block + j]
    //i runs from 0 to 2 NBlocks - 1, j runs from 0 to sites_per_block - 1 
    int* Agg;
    //LatticeBlocks[i][j] is accessed as LatticeBlocks[i * sites_per_block + j]
    //i runs from 0 to NBlocks - 1, j runs from 0 to sites_per_block - 1 
    int *LatticeBlocks; 

    void makeAggregates();

    //---The coarsest level does not need these operators---//
    /*
	Prolongation operator times a spinor x = P v
	x_i = P_ij v_j. dim(P) = DOF Nsites x Ntest Nagg, 
	dim(v) = [NBlocks][2*Ntest], dim(x) = [Nsites][DOF]
    */
    void P_v(const spinor& v,spinor& out);

    /*
	Restriction operator times a spinor on the coarse grid, x = P^H v
	x_i = P^dagg_ij v_j. dim(P^dagg) =  Ntest Nagg x DOF Nsites,
	dim(v) = [Nsites][DOF], dim(x) = [NBlocks][2*Ntest] 
    */
    void Pt_v(const spinor& v,spinor& out);

    void D_operator(const spinor& v,spinor& out); //--> This operator has to be constructed with the gauge links from the 
    //previous level

    
};


#endif