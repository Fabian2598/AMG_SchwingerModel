#include "variables.h"

typedef std::complex<double> c_double;
double coarse_time = 0.0; //Time spent in the coarse grid solver
double smooth_time = 0.0; //Time spent in the smoother


std::vector<std::vector<int>>Coords = std::vector<std::vector<int>>(LV::Nx, std::vector<int>(LV::Nt, 0));
void Coordinates() {
	for (int x = 0; x < LV::Nx; x++) {
		for (int t = 0; t < LV::Nt; t++) {
			Coords[x][t] = x * LV::Nt+ t;
		}
	}
}

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
    int sap_gmres_restart_length = 5; //GMRES restart length for the Schwarz blocks. Set to 20 by default
    int sap_gmres_restarts = 5; //GMRES iterations for the Schwarz blocks. Set to 10 by default.
    double sap_gmres_tolerance = 1e-3; //GMRES tolerance for the Schwarz blocks
    double sap_tolerance = 1e-10; //Tolerance for the SAP method
    int sap_blocks_per_proc = 1; //Number of blocks per process for the parallel SAP method
}

namespace AMGV {
    int SAP_test_vectors_iterations = 2; //Number of SAP iterations to smooth test vectors
    //Parameters for the coarse level solver. They can be changed in the main function
    int gmres_restarts_coarse_level = 10; 
    int gmres_restart_length_coarse_level = 20; //GMRES restart length for the coarse level
    double gmres_tol_coarse_level = 0.1; //GMRES tolerance for the coarse level

    int gmres_restarts_smoother = 20; //Iterations for GMRES as a smoother (SAP is the default)

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

std::vector<std::vector<std::vector<int>>> RightPB_l = std::vector<std::vector<std::vector<int>>>
    (LEVELS,std::vector<std::vector<int>>(LV::Ntot,std::vector<int>(2,0)) );
std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
std::vector<std::vector<std::vector<int>>> LeftPB_l = std::vector<std::vector<std::vector<int>>>
    (LEVELS,std::vector<std::vector<int>>(LV::Ntot,std::vector<int>(2,0)) );
std::vector<std::vector<std::vector<c_double>>> SignR_l = std::vector<std::vector<std::vector<c_double>>>
    (LEVELS,std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)) );
std::vector<std::vector<std::vector<c_double>>> SignL_l = std::vector<std::vector<std::vector<c_double>>>
    (LEVELS,std::vector<std::vector<c_double>>(LV::Ntot,std::vector<c_double>(2,0)) );



void save_vec(const std::vector<double>& vec,const std::string& Name){
    std::ofstream Datfile(Name);
    if (!Datfile.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    int size = vec.size();
    for (int n = 0; n < size; n++) {
        Datfile << n
                << std::setw(30) << std::setprecision(17) << std::scientific << vec[n]
                << "\n";
    }
        
    
}

void read_rhs(std::vector<std::vector<c_double>>& vec,const std::string& name){
    std::ifstream infile(name);
    if (!infile) {
        std::cerr << "File " << name << " not found" << std::endl;
    }
    int x, t, mu;
    double re, im;
    //x, t, mu, real part, imaginary part
    while (infile >> x >> t >> mu >> re >> im) {
        vec[Coords[x][t]][mu] = c_double(re, im); 
    }
    infile.close();
  
}

void save_rhs(std::vector<std::vector<c_double>>& rhs,const std::string& name){
    std::ofstream rhsfile(name);
    if (!rhsfile.is_open()) {
        std::cerr << "Error opening rhs.txt for writing." << std::endl;
    } 
    else {
        int x,t;
        //x, t, mu, real part, imaginary part
        for (int n = 0; n < LV::Ntot; ++n) {
            x = n/LV::Nt;
            t = n%LV::Nt;
            rhsfile << x << std::setw(30) << t << std::setw(30) << 0 << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::real(rhs[n][0]) << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::imag(rhs[n][0]) << "\n";

            rhsfile << x << std::setw(30) << t << std::setw(30) << 1 << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::real(rhs[n][1]) << std::setw(30)
                    << std::setprecision(17) << std::scientific << std::imag(rhs[n][1]) << "\n";
          
        }
        rhsfile.close();
    }

}


void random_rhs(std::vector<std::vector<c_double>>& vec,const int seed){
    c_double I_number(0,1);
    static std::mt19937 randomInt(seed);
	std::uniform_real_distribution<double> distribution(-1.0, 1.0); //mu, standard deviation
    for(int i = 0; i < LV::Ntot; i++) {
        vec[i][0] = distribution(randomInt) + I_number * distribution(randomInt); //RandomU1();
        vec[i][1] = distribution(randomInt) + I_number * distribution(randomInt);
    }

}


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
        std::cout << "| GMRES restarts for SAP blocks = " << SAPV::sap_gmres_restarts << std::endl;
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