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
#include "level.h"
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
    mass::m0 = -0.18840579710144945;
    double m0 = mass::m0; 

    //Parameters in variables.cpp
    int rank = 0;
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
            << " Number of lattice blocks" << NBlocks[l] << std::endl;
        }
        for(int l=0; l< AMGV::levels; l++){
            std::cout << "Level " << l << " Nsites " << Nsites[l] 
            << " Nxsites " << NxSites[l] << " NtSites " << NtSites[l] << " DOF " << DOF[l]
            << " Colors " << Colors[l] << std::endl;
        }
    }
    Aggregates();

    int level0 = 0;
    int level1 = 1;
    Level Level0(level0,GConf.Conf);
    Level Level1(level1,GConf.Conf);


    Level0.makeAggregates();
    Level1.makeAggregates();

    Level0.makeBlocks();
    Level1.makeBlocks();

    Level0.makeDirac(); //Builds Dirac operator at level = l

    
    spinor v(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    for(int n = 0; n < LevelV::Nsites[0]; n++){
        v[n][0] = RandomU1(); v[n][1] = RandomU1();
    }
    spinor out(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    spinor outv2(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    D_phi(GConf.Conf,v,out,mass::m0);
    Level0.D_operator(v,outv2);

    for(int n = 0; n < LevelV::Nsites[0]; n++){
    for(int mu : {0,1}){
        if (std::abs(out[n][mu]-outv2[n][mu]) > 1e-8 ){
            std::cout << "[" << n << "][" << mu << "] " << "is different" << std::endl; 
            std::cout << out[n][mu] << "   /=    " << outv2[n][mu] << std::endl;
            return 1;
        }
        std::cout << out[n][mu] << "   =    " << outv2[n][mu] << std::endl;
    }
    }
 

    
    
   
    
    return 0;
}