#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "params.h"
#include "bi_cgstab.h"
#include "conjugate_gradient.h"
#include "boundary.h"
#include "amg.h"
#include "AMG.h"
#include "mpi.h"



//Formats decimal numbers
//For opening file with confs 
static std::string format(const double& number) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << number;
    std::string str = oss.str();
    str.erase(str.find('.'), 1); //Removes decimal dot 
    return str;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    GaugeConf GConf = GaugeConf(LV::Nx, LV::Nt);
    GConf.initialize();
    readParameters("../parameters.dat");
    srand(19);

    //srand(time(0));
    
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t
    MakeBlocks(); //Makes lattice blocks 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    boundary();

    //double m0 = -0.5;
    mass::m0 = -0.5;//-0.18840579710144945;
    double m0 = mass::m0; 

    MPI_Barrier(MPI_COMM_WORLD);
    //Parameters in variables.cpp
    if (rank == 0){
        using namespace LV; //Lattice parameters namespace
        std::cout << "******************* AMG for the Dirac matrix in the Schwinger model *******************" << std::endl;
        std::cout << " Nx = " << Nx << " Nt = " << Nt << std::endl;
        std::cout << " Lattice dimension = " << (Nx * Nt) << std::endl;
        std::cout << " Number of entries of the Dirac matrix = (" << (2 * Nx * Nt) << ")^2 = " << (2 * Nx * Nt) * (2 * Nx * Nt) << std::endl;
        std::cout << " Bare mass parameter m0 = " << m0 << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        
        std::cout << "Blocks, aggregates and test vectors for each level" << std::endl;
        using namespace LevelV; //Lattice parameters namespace
        for(int l=0; l< AMGV::levels-1; l++){
            std::cout << "Level " << l << " Block X " << BlocksX[l] 
            << " Block T " << BlocksT[l] << " Ntest " << Ntest[l] << " Nagg " << Nagg[l]
            << " Number of lattice blocks " << NBlocks[l] 
            << " Schwarz Block T " << SAP_Block_t[l] << " Schwarz Block X " << SAP_Block_x[l] << std::endl;
        }
        for(int l=0; l< AMGV::levels; l++){
            std::cout << "Level " << l << " Nsites " << Nsites[l] 
            << " Nxsites " << NxSites[l] << " NtSites " << NtSites[l] << " DOF " << DOF[l]
            << " Colors " << Colors[l] << std::endl;
        }
        std::cout << "\n";
    }
    Aggregates();
    sap.set_params(GConf.Conf, m0); 
    //AMG amg(GConf, m0,AMGV::nu1, AMGV::nu2); 
	//amg.setUpPhase(1,AMGV::Nit);
    //amg.initializeCoarseLinks();

    //Three levels tests
    //TestCoarseGaugeFieldsV1(GConf);
    //checkSAPV1(GConf);
    
    //Four levels tests
    //TestCoarseGaugeFieldsV2(GConf);
    //checkSAPV2(GConf);

    spinor rhs(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    spinor x(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));

    AlgebraicMG AMG(GConf, m0,AMGV::nu1, AMGV::nu2);
    AMG.setUpPhase(1,3);
    MPI_Barrier(MPI_COMM_WORLD);

    //AMG.testSetUp();
    AMG.applyMultilevel(50, rhs,x,1e-10,true);
    spinor Dphi(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    D_phi(GConf.Conf,x,Dphi,m0);
    //Check if D_phi and rhs are equal
    for(int n=0; n<LevelV::Nsites[0]; n++){
        for(int c=0; c<LevelV::DOF[0]; c++){
            if (std::abs(Dphi[n][c]-rhs[n][c])>1e-8){
                std::cout << "Error: D_phi and rhs are not equal at site " << n << " component " << c << std::endl;
                std::cout << "D_phi: " << Dphi[n][c] << " rhs: " << rhs[n][c] << std::endl;
                exit(1);
            }
        }
    }
    

    MPI_Finalize();

    return 0;
}


