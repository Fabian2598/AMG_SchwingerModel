#ifndef LEVEL_H
#define LEVEL_H

#include "variables.h"
#include "operator_overloads.h"
#include <algorithm>
#include "dirac_operator.h"
#include "gauge_conf.h"

/*
    One level of the AMG method
*/
class Level {
public:
    Level(const int& level, const c_matrix& U) : level(level), U(U) {
        test_vectors = std::vector<spinor>(Ntest,
        spinor( Nsites, c_vector (DOF,0))); 
	    interpolator_columns = std::vector<spinor>(Ntest,
        spinor( Nsites, c_vector (DOF,0))); 
	    
        LatticeBlocks = std::vector<std::vector<int>> (NBlocks, std::vector<int>(sites_per_block,0));

        v_chopped = std::vector<spinor>(Ntest*Nagg, spinor(Nsites, c_vector(DOF,0))); //For orthonormalization
        
        //For level = 0 DOF[level] = 2
        //For level = 1 DOF[level] = 2 * LevelV::Ntest[level-1] = 2 * LevelV::Colors[level]
        Agg = new int[NBlocks * DOF * sites_per_block]; 
        //LatticeBlocks = new int[NBlocks * sites_per_block];
        nCoords = new int[Nsites * 2 * colors];
        sCoords = new int[Nsites * 2 * colors];
        cCoords = new int[Nsites * 2 * colors];

        //Gauge links to define D_operator (matrix problem at this level)
        G1 = c_vector(Nsites*2*2*colors*colors,0);
        G2 = c_vector(Nsites*2*2*colors*colors*2,0);
        G3 = c_vector(Nsites*2*2*colors*colors*2,0);

        if (level == 0){
            makeDirac();
        }
        
        std::cout << "Level " << level << " initialized with the following parameters " << "colors " << colors
        << " NBlocks " << NBlocks << " Nsites " << Nsites
        << " Ntest " << Ntest << " Nagg " << Nagg << std::endl;
        std::cout << "Nxsites " << Nxsites << " NtSites " << Ntsites << " DOF " << DOF << std::endl;
        std::cout << "x_elements " << x_elements << " t_elements " << t_elements << std::endl;
        std::cout << "sites_per_block " << sites_per_block << "\n" << std::endl;
        

    };

    ~Level() {
        delete[] Agg;
        delete[] nCoords;
        delete[] sCoords;
        delete[] cCoords;
    }

    std::vector<spinor> test_vectors; //[Ntest][Nsites][degrees of freedom per site]
    std::vector<spinor> interpolator_columns;
//private:
    const int level; 
    const int x_elements = LevelV::NxSites[level] / LevelV::BlocksX[level], t_elements = LevelV::NtSites[level] / LevelV::BlocksT[level]; //x and t elements of each lattice block
    const int sites_per_block = x_elements * t_elements;
    const int NBlocks = LevelV::NBlocks[level]; //Number of lattice blocks 
    const int colors = LevelV::Colors[level];   //Number of colors at this level
    const int Nsites = LevelV::Nsites[level];   //Number of lattice sites at this level
    const int Ntest = LevelV::Ntest[level];     //Number of test vectors to go to the next level
    const int Nagg = LevelV::Nagg[level];       //Number of aggregates to go to the next level
    const int DOF = LevelV::DOF[level];         //Degrees of freedom at each lattice site at this level
    int Ntsites = LevelV::NtSites[level];       //Number of time sites at this level
    int Nxsites = LevelV::NxSites[level];       //Number of space sites at this level
    const c_matrix U; //gauge configuration

    //At level = 0 these vectors represent the gauge links.
    //At level > 1 they are the coarse gauge links generated in the previous level
    c_vector G1; 
    c_vector G2; 
    c_vector G3; 

    //Index functions for gauge links. These correspond to the current level
	//get index for A_coeff 1D array
    //[A(x)]^{alf,bet}_{c,b} --> A_coeff[x][alf][bet][c][b]
	inline int getG1index(const int& x, const int& alf, const int& bet, const int& c, const int& b){
		return x * 2 * 2 * colors * colors 
        + alf * 2 * colors * colors 
        + bet * colors * colors
        + c * colors 
        + b;
	}
	//[B_mu(x)]^{alf,bet}_{c,b}  --> B_coeff[x][alf][bet][c][b][mu]
    //[C_mu(x)]^{alf,bet}_{c,b}  --> C_coeff[x][alf][bet][c][b][mu]
	inline int getG2G3index(const int& x, const int& alf, const int& bet, const int& c, const int& b, const int& mu){
        return x * 2 * 2 * colors * colors * 2 
        + alf * 2 * colors * colors * 2 
        + bet * colors * colors * 2
        + c * colors * 2 
        + b * 2 
        + mu;
    }
    	
    //Index functions for coarse gauge links. These correspond to the next level, but are generated here (not stored)
	//get index for A_coeff 1D array
    //[A(x)]^{alf,bet}_{c,b} --> A_coeff[x][alf][bet][c][b]
	inline int getAindex(const int& block, const int& alf, const int& bet, const int& c, const int& b){
		return block * 2 * 2 * Ntest * Ntest 
        + alf * 2 * Ntest * Ntest 
        + bet * Ntest * Ntest
        + c * Ntest 
        + b;
	}
	//[B_mu(x)]^{alf,bet}_{c,b}  --> B_coeff[x][alf][bet][c][b][mu]
    //[C_mu(x)]^{alf,bet}_{c,b}  --> C_coeff[x][alf][bet][c][b][mu]
	inline int getBCindex(const int& block, const int& alf, const int& bet, const int& c, const int& b, const int& mu){
        return block * 2 * 2 * Ntest * Ntest * 2 
        + alf * 2 * Ntest * Ntest * 2 
        + bet * Ntest * Ntest * 2
        + c * Ntest * 2 
        + b * 2 
        + mu;
    }

    /*
    For level = 0
    Agg[i][j] is accessed as Agg[i * sites_per_block + j]
    i: 0 to 2 NBlocks - 1, j: 0 to sites_per_block - 1 
    For 0 < level < maxLevel 
    Agg[i][j][k] is accessed as Agg[i * sites_per_block * Colors + j * Colors + k]
    i: 0 to 2 NBlocks - 1, j: 0 to sites_per_block - 1, k: 0 to Colors - 1
    Colors is just the number of test vectors at the previous level
    For the coarsest level we don't need need any aggregation
    */
    int* Agg;

    int* nCoords; int* sCoords; int* cCoords;
    
    std::vector<std::vector<int>> LatticeBlocks;

    //LatticeBlocks[i][j] is accessed as LatticeBlocks[i * sites_per_block + j]
    //i runs from 0 to NBlocks - 1, j runs from 0 to sites_per_block - 1 
    //int *LatticeBlocks; 

    void makeAggregates();
    void printAggregates();

    void makeBlocks();
    void printBlocks();

    //Creates G1, G2 and G3
    void makeDirac();

    //Make coarse gauge links. They will be used in the next level as G1, G2 and G3.
    void makeCoarseLinks(Level& next_level);//& A_coeff,c_vector& B_coeff, c_vector& C_coeff);

    void setUp(); //This is just for testing

    void orthonormalize(); //Local orthonormalization of the test vectors
    void checkOrthogonality(); //Check orthogonality of the test vectors
    std::vector<spinor> v_chopped;


    /*
    Matrix-vector operation that defines the level l.
    For instance, at level = 0, D_operator is just the Dirac operator
    at level = 1 D_operator is Dc
    at level = 2 D_operator is (Dc)_c ...
    */
    void D_operator(const spinor& v, spinor& out);
    
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

    inline void getLatticeBlock(const int& n, int &block) {
        int x = n / Ntsites; //x coordinate of the lattice point 
        int t = n % Ntsites; //t coordinate of the lattice point
        //Reconstructing the block and m index from x and t
        int block_x = x / x_elements; //Block index in the x direction
        int block_t = t / t_elements; //Block index in the t direction
        block = block_x * LevelV::BlocksT[level] + block_t; //Block index in the SAP method

        //int mx = x % x_elements; //x coordinate in the block
        //int mt = t % t_elements; //t coordinate in the block
        //m = mx * t_elements + mt; //Index in the block
    }
};




#endif