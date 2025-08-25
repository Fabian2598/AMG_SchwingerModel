#ifndef LEVEL_H
#define LEVEL_H

#include "variables.h"
#include "operator_overloads.h"
#include <algorithm>
#include "dirac_operator.h"
#include "gauge_conf.h"
#include "sap.h"

/*
    One level of the AMG method
*/
class Level {
public:   
    //SAP for smoothing D_operator
    //-------------------------------Nested class-------------------------------//
    class SAP_level_l : public SAP_C {
    public:
        SAP_level_l(const int& dim1, const int& dim2, const double& tol,const int& Nt, const int& Nx,const int& block_x,const int& block_t,
        const int& spins, const int& colors,Level* parent) :
        SAP_C(dim1, dim2, tol, Nt, Nx, block_x, block_t,spins,colors), parent(parent) {
        }

 

    private: 
        Level* parent; //Parent class
        /*
        Global D operation
        */
        void funcGlobal(const spinor& in, spinor& out) override { 
            parent->D_operator(in, out); //Dirac operator at the current level
        }

        /*
        Local D operations
        */
        void D_local(const spinor& in, spinor& out, const int& block);

        void funcLocal(const spinor& in, spinor& out) override { 
            //std::cout << "funcLocal called for block " << blockMPI << std::endl;
            D_local( in, out,blockMPI);
        }

        /*
            Given a lattice point index n, it returns the corresponding 
            SAP block index and the local index m within that block.
        */
        inline void getMandBlock(const int& n, int &m, int &block) {
            int x = n / Nx; //x coordinate of the lattice point 
            int t = n % Nt; //t coordinate of the lattice point
            //Reconstructing the block and m index from x and t
            int block_x = x / x_elements; //Block index in the x direction
            int block_t = t / t_elements; //Block index in the t direction
            block = block_x * Block_t + block_t; //Block index in the SAP method

            int mx = x % x_elements; //x coordinate in the block
            int mt = t % t_elements; //t coordinate in the block
            m = mx * t_elements + mt; //Index in the block
        }

    };

    SAP_level_l sap_l; 
    //----------------------------------------------------------------------------//
    //GMRES for the current level. We use it for solving the coarsest system.
    class GMRES_level_l : public FGMRES {
	public:
    	GMRES_level_l(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol, Level* parent) : 
		FGMRES(dim1, dim2, m, restarts, tol), parent(parent) {}
    
    	~GMRES_level_l() { };
    
	private:
		Level* parent; //Pointer to the enclosing AMG instance
    	/*
    	Implementation of the function that computes the matrix-vector product for the fine level
    	*/
    	void func(const spinor& in, spinor& out) override {
        	parent->D_operator(in,out);
    	}
		//No preconditioning for the coarsest level
		void preconditioner(const spinor& in, spinor& out) override {
        out = std::move(in); //Identity operation
		}
	};

	GMRES_level_l gmres_l;
    //----------------------------------------------------------------------------//
    
    //Level Constructor
    Level(const int& level, const c_matrix& U) : level(level), U(U),
        sap_l(LevelV::Nsites[level], 
            LevelV::DOF[level], 
            SAPV::sap_tolerance,
            LevelV::NtSites[level], 
            LevelV::NxSites[level],
            LevelV::SAP_Block_x[level],
            LevelV::SAP_Block_t[level],
            2, //two spins
            LevelV::Colors[level],
            this),
         gmres_l(LevelV::Nsites[level], LevelV::DOF[level],
            LevelV::GMRES_restart_len[level],
            LevelV::GMRES_restarts[level],
            LevelV::GMRES_tol[level],
            this) 
    {
        /*
        std::cout << "Level " << level << " initialized with the following parameters " << "colors " << colors
        << " NBlocks " << NBlocks << " Nsites " << Nsites
        << " Ntest " << Ntest << " Nagg " << Nagg << std::endl;
        std::cout << "Nxsites " << Nxsites << " NtSites " << Ntsites << " DOF " << DOF << std::endl;
        std::cout << "x_elements " << x_elements << " t_elements " << t_elements << std::endl;
        std::cout << "sites_per_block " << sites_per_block << "\n" << std::endl;
        */

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
    const int x_elements = (level != LevelV::maxLevel) ?  LevelV::NxSites[level] / LevelV::BlocksX[level]: 1;
    const int t_elements = (level != LevelV::maxLevel) ?  LevelV::NtSites[level] / LevelV::BlocksT[level]: 1; //x and t elements of each lattice block
    const int sites_per_block = x_elements * t_elements;
    const int NBlocks = (level != LevelV::maxLevel) ? LevelV::NBlocks[level]: 1; //Number of lattice blocks 
    const int colors = LevelV::Colors[level];   //Number of colors at this level
    const int Nsites = LevelV::Nsites[level];   //Number of lattice sites at this level
    const int Ntest = (level != LevelV::maxLevel) ? LevelV::Ntest[level]: 1;     //Number of test vectors to go to the next level
    const int Nagg = (level != LevelV::maxLevel) ? LevelV::Nagg[level]: 1;       //Number of aggregates to go to the next level
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


inline void TestCoarseGaugeFieldsV1(const GaugeConf& GConf){
    /*
    This test was done for four levels 
    paramaters.dat 
    //Use this data for parameters.dat
    //0 4 4 10 4 4
    //1 2 2 5 2 2

    //Lattice size 64x64
    
    It could be reused for other set of parameters, but I have to check them carefully
    */
    int level0 = 0;
    int level1 = 1;
    int level2 = 2;
    Level Level0(level0,GConf.Conf);
    Level Level1(level1,GConf.Conf);
    Level Level2(level2,GConf.Conf);


    
    Level0.makeAggregates();
    std::cout << "Aggregates Level 0 built " << std::endl;
    Level1.makeAggregates();
    std::cout << "Aggregates Level 1 built " << std::endl;

    Level0.makeBlocks();
    Level1.makeBlocks();
    std::cout << "Blocks built " << std::endl;

   
    //---------------------------------------//
    //Checking operator at level = 1
    std::cout << "Set up and orthonormalization ... " << std::endl;
    Level0.setUp(); //Build test vectors 
    Level0.checkOrthogonality();
    Level0.makeCoarseLinks(Level1); //D_operator for level 1

    Level1.setUp(); //Build test vectors for level 1
    Level1.checkOrthogonality();
    Level1.makeCoarseLinks(Level2); //D_operator for level 2

    
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Testing that P^dag D P = D_c level by level 

    //---Testing level 1
    spinor in1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1)); //in
    spinor temp1(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor Dphi1(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor out1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //out
    spinor out1_v2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //D_c
    //P^H D P
    Level0.P_v(in1,temp1);
    Level0.D_operator(temp1,Dphi1);
    Level0.Pt_v(Dphi1,out1);

    Level1.D_operator(in1,out1_v2);
    std::cout << "Testing level 1" << std::endl;
    for(int x = 0; x<LevelV::Nsites[level1]; x++){
        for(int dof = 0; dof<LevelV::DOF[level1]; dof++){
            if (std::abs(out1[x][dof]-out1_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 1  different" << std::endl; 
            std::cout << out1[x][dof] << "   /=    " << out1_v2[x][dof] << std::endl;
            return;
            }
        }
    }
    
    std::cout << "P^dag D P coincides with Dc for level 1" << std::endl;
    std::cout << out1[0][0] << "   =    " << out1_v2[0][0] << std::endl;

  
    //----Testing level 2
    std::cout << "Testing level 2" << std::endl;
    spinor in2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],1));
    spinor temp2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor Dphi2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor out2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    spinor out2_v2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    Level1.P_v(in2,temp2);
    Level1.D_operator(temp2,Dphi2);
    Level1.Pt_v(Dphi2,out2);
    
    

    Level2.D_operator(in2,out2_v2);
    for(int x = 0; x<LevelV::Nsites[level2]; x++){
        for(int dof = 0; dof<LevelV::DOF[level2]; dof++){
            if (std::abs(out2[x][dof]-out2_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 2  different" << std::endl; 
            std::cout << out2[x][dof] << "   /=    " << out2_v2[x][dof] << std::endl;
            return;
            }
        }
    }

    std::cout << "P^dag D P coincides with Dc for level 2" << std::endl;
    std::cout << out2[0][0] << "   =    " << out2_v2[0][0] << std::endl;

}




inline void TestCoarseGaugeFieldsV2(const GaugeConf& GConf){
    /*
    This test was done for four levels 
    paramaters.dat 
    0 8 8 10
    1 4 4 5
    2 2 2 5

    and Nt = Nx = 32
    
    It could be reused for other set of parameters, but I have to check them carefully
    */
    int level0 = 0;
    int level1 = 1;
    int level2 = 2;
    int level3 = 3;
    Level Level0(level0,GConf.Conf);
    Level Level1(level1,GConf.Conf);
    Level Level2(level2,GConf.Conf);
    Level Level3(level3,GConf.Conf);

    
    Level0.makeAggregates();
    std::cout << "Aggregates Level 0 built " << std::endl;
    Level1.makeAggregates();
    std::cout << "Aggregates Level 1 built " << std::endl;
    Level2.makeAggregates(); //-->For the coarsest level this is not necessary
    std::cout << "Aggregates Level 2 built " << std::endl;

    Level0.makeBlocks();
    Level1.makeBlocks();
    Level2.makeBlocks(); //--> Not necessary for the coarsest level
    std::cout << "Blocks built " << std::endl;

   
    //---------------------------------------//
    //Checking operator at level = 1
    std::cout << "Set up and orthonormalization ... " << std::endl;
    Level0.setUp(); //Build test vectors 
    Level0.checkOrthogonality();
    Level0.makeCoarseLinks(Level1); //D_operator for level 1

    Level1.setUp(); //Build test vectors for level 1
    Level1.checkOrthogonality();
    Level1.makeCoarseLinks(Level2); //D_operator for level 2

    Level2.setUp(); //Build test vectors for level 2
    Level2.checkOrthogonality();
    Level2.makeCoarseLinks(Level3); //D_operator for level 3
    
    

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Testing that P^dag D P = D_c level by level 

    
    //---Testing level 1
    spinor in1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1)); //in
    spinor temp1(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor Dphi1(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor out1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //out
    spinor out1_v2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //D_c
    //P^H D P
    Level0.P_v(in1,temp1);
    Level0.D_operator(temp1,Dphi1);
    Level0.Pt_v(Dphi1,out1);

    Level1.D_operator(in1,out1_v2);
    std::cout << "Testing level 1" << std::endl;
    for(int x = 0; x<LevelV::Nsites[level1]; x++){
        for(int dof = 0; dof<LevelV::DOF[level1]; dof++){
            if (std::abs(out1[x][dof]-out1_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 1  different" << std::endl; 
            std::cout << out1[x][dof] << "   /=    " << out1_v2[x][dof] << std::endl;
            return;
            }
        }
    }
    
    std::cout << "P^dag D P coincides with Dc for level 1" << std::endl;
    std::cout << out1[0][0] << "   =    " << out1_v2[0][0] << std::endl;

  
    //----Testing level 2
    std::cout << "Testing level 2" << std::endl;
    spinor in2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],1));
    spinor temp2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor Dphi2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor out2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    spinor out2_v2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    Level1.P_v(in2,temp2);
    Level1.D_operator(temp2,Dphi2);
    Level1.Pt_v(Dphi2,out2);
    
    

    Level2.D_operator(in2,out2_v2);
    for(int x = 0; x<LevelV::Nsites[level2]; x++){
        for(int dof = 0; dof<LevelV::DOF[level2]; dof++){
            if (std::abs(out2[x][dof]-out2_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 2  different" << std::endl; 
            std::cout << out2[x][dof] << "   /=    " << out2_v2[x][dof] << std::endl;
            return;
            }
        }
    }

    std::cout << "P^dag D P coincides with Dc for level 2" << std::endl;
    std::cout << out2[0][0] << "   =    " << out2_v2[0][0] << std::endl;

    //----Testing level 3
    std::cout << "Testing level 2" << std::endl;
    spinor in3(LevelV::Nsites[level3],c_vector(LevelV::DOF[level3],1));
    spinor temp3(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    spinor Dphi3(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    spinor out3(LevelV::Nsites[level3],c_vector(LevelV::DOF[level3],0));
    spinor out3_v2(LevelV::Nsites[level3],c_vector(LevelV::DOF[level3],0));
    Level2.P_v(in3,temp3);
    Level2.D_operator(temp3,Dphi3);
    Level2.Pt_v(Dphi3,out3);
    
    

    Level3.D_operator(in3,out3_v2);
    for(int x = 0; x<LevelV::Nsites[level3]; x++){
        for(int dof = 0; dof<LevelV::DOF[level3]; dof++){
            if (std::abs(out3[x][dof]-out3_v2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 3  different" << std::endl; 
            std::cout << out3[x][dof] << "   /=    " << out3_v2[x][dof] << std::endl;
            return;
            }
        }
    }

    std::cout << "P^dag D P coincides with Dc for level 3" << std::endl;
    std::cout << out3[0][0] << "   =    " << out3_v2[0][0] << std::endl;


}


inline void checkSAPV1(const GaugeConf& GConf){
    //Use this data for parameters.dat
    //0 4 4 10 4 4
    //1 2 2 5 2 2

    //Lattice size 64x64
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int level0 = 0;
    int level1 = 1;
    int level2 = 2;
    Level Level0(level0, GConf.Conf);
    
    Level Level1(level1,GConf.Conf);

    Level Level2(level2,GConf.Conf);
    
    Level0.makeAggregates();
    std::cout << "Aggregates Level 0 built " << std::endl;
    Level1.makeAggregates();
    std::cout << "Aggregates Level 1 built " << std::endl;

    Level0.makeBlocks();
    Level1.makeBlocks();

    std::cout << "Set up and orthonormalization ... " << std::endl;
    Level0.setUp(); //Build test vectors 
    Level0.checkOrthogonality();
    Level0.makeCoarseLinks(Level1); //D_operator for level 1

    Level1.setUp(); //Build test vectors for level 1
    Level1.checkOrthogonality();
    Level1.makeCoarseLinks(Level2); //D_operator for level 2
    

    //--------------------Level 0-------------------//
    
    spinor rhs(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],1)); //in
    //for(int n = 0; n< LevelV::Nsites[level0]; n++){
    //for(int dof = 0; dof < LevelV::DOF[level0]; dof++){
    //    rhs[n][dof] = RandomU1();
    //}
    //}

    spinor x(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in
    int iter = 100;
    MPI_Barrier(MPI_COMM_WORLD);
    sap.SAP(rhs,x,iter,SAPV::sap_blocks_per_proc,true);

    
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Testing Level0 " << std::endl;
    //Check D_local implementation
    spinor xtest(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in 
    Level0.sap_l.SAP(rhs,xtest,iter,SAPV::sap_blocks_per_proc,true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level0]; n++){
        for(int alf=0; alf<LevelV::DOF[level0]; alf++){
            if(std::abs(xtest[n][alf] - x[n][alf]) > 1e-8){
                std::cout << "Error in SAP implementation at level " << level0 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest[n][alf] << ", x = " << x[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "Previous SAP implementation coincides with the new one " << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    spinor xgmres(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in 
    Level0.gmres_l.fgmres(rhs,xgmres,xgmres,true);

    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level0]; n++){
        for(int alf=0; alf<LevelV::DOF[level0]; alf++){
            if(std::abs(xtest[n][alf] - xgmres[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << level0 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest[n][alf] << ", x = " << xgmres[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << level0 << std::endl;
    }

    //--------------------Level 1-------------------//

    std::cout << "Testing Level1 " << std::endl;
    spinor rhs1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1)); //in
    //for(int n = 0; n< LevelV::Nsites[level0]; n++){
    //for(int dof = 0; dof < LevelV::DOF[level0]; dof++){
    //    rhs[n][dof] = RandomU1();
    //}
    //}

    
    MPI_Barrier(MPI_COMM_WORLD);
    //Check D_local implementation
    spinor xtest1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //in 
    Level1.sap_l.SAP(rhs1,xtest1,iter,SAPV::sap_blocks_per_proc,true);

    MPI_Barrier(MPI_COMM_WORLD);
    spinor xgmres1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //in 
    Level1.gmres_l.fgmres(rhs1,xgmres1,xgmres1,true);

    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level1]; n++){
        for(int alf=0; alf<LevelV::DOF[level1]; alf++){
            if(std::abs(xtest1[n][alf] - xgmres1[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << level1 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest1[n][alf] << ", x = " << xgmres1[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << level1 << std::endl;
        std::cout << "xtest = " << xtest1[0][0] << ", x = " << xgmres1[0][0] << std::endl;
    }
}

inline void checkSAPV2(const GaugeConf& GConf){
    //lattice size 64 x 64
    //parameters.dat
    //level, xblock, tblock, ntest, sapXblock, sapTblock
    //0 8 8 10 4 4
    //1 4 4 5 4 4
    //2 2 2 7 2 2  


    //We also tried
    //0 16 16 10 4 4
    //1 8 8 5 4 4
    //2 4 4 7 4 4  

    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int level0 = 0;
    int level1 = 1;
    int level2 = 2;
    int level3 = 3;
    Level Level0(level0, GConf.Conf);
    Level Level1(level1,GConf.Conf);
    Level Level2(level2,GConf.Conf);
    Level Level3(level3,GConf.Conf);
    
    Level0.makeAggregates();
    std::cout << "Aggregates Level 0 built " << std::endl;
    Level1.makeAggregates();
    std::cout << "Aggregates Level 1 built " << std::endl;
    Level2.makeAggregates();
    std::cout << "Aggregates Level 2 built " << std::endl;

    Level0.makeBlocks();
    Level1.makeBlocks();
    Level2.makeBlocks();

    std::cout << "Set up and orthonormalization ... " << std::endl;
    Level0.setUp(); //Build test vectors 
    Level0.checkOrthogonality();
    Level0.makeCoarseLinks(Level1); //D_operator for level 1

    Level1.setUp(); //Build test vectors for level 1
    Level1.checkOrthogonality();
    Level1.makeCoarseLinks(Level2); //D_operator for level 2

    Level2.setUp(); //Build test vectors for level 1
    Level2.checkOrthogonality();
    Level2.makeCoarseLinks(Level3); //D_operator for level 2
    

    //--------------------Level 0-------------------//
    
    spinor rhs(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],1)); //in
    //for(int n = 0; n< LevelV::Nsites[level0]; n++){
    //for(int dof = 0; dof < LevelV::DOF[level0]; dof++){
    //    rhs[n][dof] = RandomU1();
    //}
    //}

    spinor x(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in
    int iter = 100;
    MPI_Barrier(MPI_COMM_WORLD);
    sap.SAP(rhs,x,iter,SAPV::sap_blocks_per_proc,true);

    
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Testing Level0 " << std::endl;
    //Check D_local implementation
    spinor xtest(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in 
    Level0.sap_l.SAP(rhs,xtest,iter,SAPV::sap_blocks_per_proc,true);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level0]; n++){
        for(int alf=0; alf<LevelV::DOF[level0]; alf++){
            if(std::abs(xtest[n][alf] - x[n][alf]) > 1e-8){
                std::cout << "Error in SAP implementation at level " << level0 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest[n][alf] << ", x = " << x[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "Previous SAP implementation coincides with the new one " << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    spinor xgmres(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0)); //in 
    Level0.gmres_l.fgmres(rhs,xgmres,xgmres,true);

    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level0]; n++){
        for(int alf=0; alf<LevelV::DOF[level0]; alf++){
            if(std::abs(xtest[n][alf] - xgmres[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << level0 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest[n][alf] << ", x = " << xgmres[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << level0 << std::endl;
    }

    //--------------------Level 1-------------------//

    std::cout << "Testing Level1 " << std::endl;
    spinor rhs1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1)); //in
    //for(int n = 0; n< LevelV::Nsites[level0]; n++){
    //for(int dof = 0; dof < LevelV::DOF[level0]; dof++){
    //    rhs[n][dof] = RandomU1();
    //}
    //}

    
    MPI_Barrier(MPI_COMM_WORLD);
    //Check D_local implementation
    spinor xtest1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //in 
    Level1.sap_l.SAP(rhs1,xtest1,iter,SAPV::sap_blocks_per_proc,true);

    MPI_Barrier(MPI_COMM_WORLD);
    spinor xgmres1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0)); //in 
    Level1.gmres_l.fgmres(rhs1,xgmres1,xgmres1,true);

    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level1]; n++){
        for(int alf=0; alf<LevelV::DOF[level1]; alf++){
            if(std::abs(xtest1[n][alf] - xgmres1[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << level1 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest1[n][alf] << ", x = " << xgmres1[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << level1 << std::endl;
        std::cout << "xtest = " << xtest1[0][0] << ", x = " << xgmres1[0][0] << std::endl;
    }

     //--------------------Level 2-------------------//

    std::cout << "Testing Level2 " << std::endl;
    spinor rhs2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],1)); //in
    //for(int n = 0; n< LevelV::Nsites[level0]; n++){
    //for(int dof = 0; dof < LevelV::DOF[level0]; dof++){
    //    rhs[n][dof] = RandomU1();
    //}
    //}

    
    MPI_Barrier(MPI_COMM_WORLD);
    //Check D_local implementation
    spinor xtest2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0)); //in 
    Level2.sap_l.SAP(rhs2,xtest2,iter,SAPV::sap_blocks_per_proc,true);

    MPI_Barrier(MPI_COMM_WORLD);
    spinor xgmres2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0)); //in 
    Level2.gmres_l.fgmres(rhs2,xgmres2,xgmres2,true);

    if (rank==0){
        //Checking that the results are the same
        for(int n=0; n<LevelV::Nsites[level2]; n++){
        for(int alf=0; alf<LevelV::DOF[level2]; alf++){
            if(std::abs(xtest2[n][alf] - xgmres2[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << level2 << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xtest = " << xtest2[n][alf] << ", x = " << xgmres2[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << level2 << std::endl;
        std::cout << "xtest = " << xtest2[0][0] << ", x = " << xgmres2[0][0] << std::endl;
    }
    


}

#endif