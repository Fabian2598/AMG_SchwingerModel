#include "variables.h"

typedef std::complex<double> c_double;
double coarse_time = 0.0; //Time spent in the coarse grid solver
double smooth_time = 0.0; //Time spent in the smoother
double SAP_time = 0.0; //Time spent in the SAP method


std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
void Coordinates() {
	for (int x = 0; x < LV::Nx; x++) {
		for (int t = 0; t < LV::Nt; t++) {
			Coords[x][t] = x * LV::Nt+ t;
		}
	}
}


//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1}
std::vector<std::vector<int>>Agg = std::vector<std::vector<int>>(2*LV::block_x*LV::block_t, std::vector<int>(LV::x_elements*LV::t_elements, 0));
std::vector<int> XCoord = std::vector<int>(2*LV::Ntot, 0);
std::vector<int> TCoord = std::vector<int>(2*LV::Ntot, 0);
std::vector<int> SCoord = std::vector<int>(2*LV::Ntot, 0);

//--Coordinates of the neighbors to avoid recomputing them each time the operator D is called--//
std::vector<std::vector<int>>LeftPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); //LeftPB[x][t][mu]
std::vector<std::vector<int>>RightPB = std::vector<std::vector<int>>(LV::Ntot, std::vector<int>(2,0)); //RightPB[x][t][mu]


std::vector<std::vector<c_double>>SignL =std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); //SignL[x][t][mu]
std::vector<std::vector<c_double>>SignR = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); ////SignR[x][t][mu]

namespace mass{
    double m0 = 0; //Default mass
}

namespace LevelV{
    int BlocksX[AMGV::levels-1];
    int BlocksT[AMGV::levels-1];
    int Ntest[AMGV::levels-1]; 
    int Nagg[AMGV::levels-1]; 
    int NBlocks[AMGV::levels-1];
    int Nsites[AMGV::levels];
    int NxSites[AMGV::levels];
    int NtSites[AMGV::levels];
    int DOF[AMGV::levels];
    int Colors[AMGV::levels];

    int SAP_Block_x[AMGV::levels];
    int SAP_Block_t[AMGV::levels];
    int SAP_elements_x[AMGV::levels]; 
    int SAP_elements_t[AMGV::levels]; 
    int SAP_variables_per_block[AMGV::levels]; 

    int GMRES_restart_len[AMGV::levels];
    int GMRES_restarts[AMGV::levels];
    double GMRES_tol[AMGV::levels];
};


namespace SAPV {
    bool schwarz_blocks = false; //Schwarz blocks are not initialized by default
    int sap_gmres_restart_length = 5; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    int sap_gmres_restarts = 5; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    double sap_gmres_tolerance = 1e-3; //GMRES tolerance for the Schwarz blocks
    double sap_tolerance = 1e-10; //Tolerance for the SAP method
    int sap_blocks_per_proc = 1; //Number of blocks per process for the parallel SAP method
}

namespace AMGV {
    int SAP_test_vectors_iterations = 1; //Number of SAP iterations to smooth test vectors
    bool aggregates_initialized = false;  //Aggregates are not initialized by default
    //Parameters for the coarse level solver. They can be changed in the main function
    int gmres_restarts_coarse_level = 10; 
    int gmres_restart_length_coarse_level = 20; //GMRES restart length for the coarse level
    double gmres_tol_coarse_level = 0.1; //GMRES tolerance for the coarse level

    int gmres_restarts_smoother = 20; //Iterations for GMRES as a smoother (SAP is the default)

    int bi_cgstab_Dc_iterations= 1000; //Number of iterations for the bi-cgstab method
    double bi_cgstab_Dc_iterations_tol = 1e-10; //Tolerance for the bi-cgstab method
    int nu1 = 0; //Pre-smoothing iterations
    int nu2 = 2; //Post-smoothing iterations
    int Nit = 1; //Number of iterations for improving the interpolator

    bool SetUpDone = false; //Set to true when the setup is done

    //Parameters for FGMRES used in the k-cycle
    int fgmres_k_cycle_restart_length = 5;
    int fgmres_k_cycle_restarts = 2;
    double fgmres_k_cycle_tol = 0.1;

    int cycle = 0; //Cycling stratey. Cycle = 0 -> V-cycle, = 1 --> K-cycle
}

//--------------Parameters for outer FGMRES solver--------------//
namespace FGMRESV {
    double fgmres_tolerance = 1e-10; //Tolerance for FGMRES
    int fgmres_restart_length = 20; //Restart length for FGMRES
    int fgmres_restarts = 50; //Number of restarts for FGMRES
}

namespace CG{
    int max_iter = 100000;
    double tol = 1e-10;
}

std::vector<std::vector<c_double>>D_TEMP = std::vector<std::vector<c_double>>(LV::Ntot, std::vector<c_double>(2,0)); 


void CheckBlocks(){
    bool check = true;
    using namespace SAPV;
    if (Nx % block_x != 0) {
        std::cout << "Error: Ns/block_x is not an integer" << std::endl;
        check = false;
    }
    if (Nt % block_t != 0) {
        std::cout << "Error: Nt/block_t is not an integer" << std::endl;
        check = false;
    }
    if (Nx % sap_block_x != 0) {
        std::cout << "Error: Ns/sap_block_x is not an integer" << std::endl;
        check = false;
    }
    if (Nt % sap_block_t != 0) {
        std::cout << "Error: Nt/sap_block_t is not an integer" << std::endl;
        check = false;
    }
    if (sap_block_t*sap_block_x % 2 != 0 ){
        std::cout << "Error: sap_block_t*sap_block_x is not even" << std::endl;
        std::cout << "Expected an even number of SAP blocks" << std::endl;
        check = false;
    }
    if (check == false){
        exit(1);
    }

}

void CheckAggregates(){
    bool check = true;
    using namespace AMGV;
    if (aggregates_initialized == false) {
        std::cout << "Error: Aggregates are not initialized" << std::endl;
        check = false;
    }
    if (Nagg*Ntest > 2*LV::Ntot) {
        std::cout << "Error: Nagg*Ntest > 2*Ntot" << std::endl;
        check = false;
    }
    if (check == false) {
        exit(1);
    }
}

std::vector<std::vector<int>> LatticeBlocks = std::vector<std::vector<int>> (LV::Nblocks, std::vector<int>(LV::lattice_sites_per_block));
int RightPB_blocks[LV::Nblocks][2];
int LeftPB_blocks[LV::Nblocks][2];

void MakeBlocks(){
	using namespace LV; //Lattice parameters namespace
	int count; 
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
				int block = x * block_t + t; //block index
				int x0 = x * x_elements, t0 = t * t_elements;
				int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
            	count = 0;  
            	for(int x = x0; x < x1; x++) {
                	for (int t = t0; t < t1; t++) {
                    	LatticeBlocks[block][count++] = x * Nt+ t; 
                    	
					}
            	}
		}
	}
}


   std::vector<std::vector<std::vector<int>>> RightPB_l = std::vector<std::vector<std::vector<int>>>
    (LEVELS,std::vector<std::vector<int>>(LV::Ntot,std::vector<int>(2,0)) );
	std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
   std::vector<std::vector<std::vector<int>>> LeftPB_l = std::vector<std::vector<std::vector<int>>>
    (LEVELS,std::vector<std::vector<int>>(LV::Ntot,std::vector<int>(2,0)) );
   std::vector<std::vector<std::vector<c_double>>> SignR_l = std::vector<std::vector<std::vector<c_double>>>
    (LEVELS,std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)) );
   std::vector<std::vector<std::vector<c_double>>> SignL_l = std::vector<std::vector<std::vector<c_double>>>
    (LEVELS,std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)) );


void printParameters(){
    using namespace LV; //Lattice parameters namespace
        std::cout << "******************* AMG for the Dirac matrix in the Schwinger model *******************" << std::endl;
        std::cout << "| Nx = " << Nx << " Nt = " << Nt << std::endl;
        std::cout << "| Lattice dimension = " << (Nx * Nt) << std::endl;
        std::cout << "| Number of entries of the Dirac matrix = (" << (2 * Nx * Nt) << ")^2" << std::endl;
        std::cout << "| Bare mass parameter m0 = " << mass::m0 << std::endl;
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;  
        std::cout << "* Blocks, aggregates and test vectors at each level" << std::endl;
        using namespace LevelV; //Lattice parameters namespace
        for(int l=0; l< AMGV::levels-1; l++){
            std::cout << "| Level " << l << " Block X " << BlocksX[l] 
            << " Block T " << BlocksT[l] << " Ntest " << Ntest[l] << " Nagg " << Nagg[l]
            << " Number of lattice blocks " << NBlocks[l] 
            << " Schwarz Block T " << SAP_Block_t[l] << " Schwarz Block X " << SAP_Block_x[l] << std::endl;
        }
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cout << "* SAP blocks" << std::endl;
        for(int l=0; l< AMGV::levels-1; l++){
            std::cout << "| Level " << l << " Schwarz Block T " << SAP_Block_t[l] << " Schwarz Block X " << SAP_Block_x[l]  
            << " Number of blocks " << SAP_Block_t[l]*SAP_Block_x[l] << 
            " Each block has " << SAP_variables_per_block[l] << " variables" << std::endl;
        }
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cout << "* Sites and degrees of freedom at each level" << std::endl;
        for(int l=0; l< AMGV::levels; l++){
            std::cout << "| Level " << l << " Nsites " << Nsites[l] 
            << " Nxsites " << NxSites[l] << " NtSites " << NtSites[l] << " DOF " << DOF[l]
            << " Colors " << Colors[l] << std::endl;
        }
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cout << "| GMRES restart length for SAP blocks = " << SAPV::sap_gmres_restart_length << std::endl;
        std::cout << "| GMRES iterations for SAP blocks = " << SAPV::sap_gmres_restarts << std::endl;
        std::cout << "| GMRES tolerance for SAP blocks = " << SAPV::sap_gmres_tolerance << std::endl;
        std::cout << "---------------------------------------------------------------------------------------" << std::endl;
        std::cout << "* AMG parameters" << std::endl;
        if (AMGV::cycle == 0)
            std::cout << "| Cycle = " << "V-cycle" << std::endl;
        else if (AMGV::cycle == 1)
            std::cout << "| Cycle = " << "K-cycle" << std::endl;
        std::cout << "| Number of levels = " << AMGV::levels << std::endl;
        std::cout << "| nu1 (pre-smoothing) = " << AMGV::nu1 << " nu2 (post-smoothing) = " << AMGV::nu2 << std::endl;
        std::cout << "| Number of iterations for improving the interpolator = " << AMGV::Nit << std::endl;
        std::cout << "| Restart length of GMRES at the coarse level = " << LevelV::GMRES_restart_len[LevelV::maxLevel] << std::endl;
        std::cout << "| Restarts of GMRES at the coarse level = " << LevelV::GMRES_restarts[LevelV::maxLevel] << std::endl;
        std::cout << "| GMRES tolerance for the coarse level solution = " << LevelV::GMRES_tol[LevelV::maxLevel] << std::endl;
        std::cout << "* FGMRES with AMG preconditioning parameters" << std::endl;
        std::cout << "| FGMRES restart length = " << FGMRESV::fgmres_restart_length << std::endl;
        std::cout << "| FGMRES restarts = " << FGMRESV::fgmres_restarts << std::endl;
        std::cout << "| FGMRES tolerance = " << FGMRESV::fgmres_tolerance << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;

}
