#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "params.h"
#include "bi_cgstab.h"
#include "conjugate_gradient.h"
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

    Level0.setUp();


    spinor v(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    spinor vNew(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    spinor w(LevelV::Ntest[0],c_vector(LevelV::Nagg[0],0));
    spinor wNew(LevelV::Nsites[1],c_vector(LevelV::DOF[1],0));
    for(int i = 0; i<LevelV::Nsites[0]; i++){
    for(int j = 0; j<LevelV::DOF[0]; j++){
        v[i][j] = RandomU1();
    }
    }


 
    AMG testAMG(GConf,m0, 0, 2);
    testAMG.setUpPhase(1,1);
    testAMG.initializeCoarseLinks();
    testAMG.Pt_v(v,w);
    Level0.Pt_v(v,wNew);
    
    int a, c, n, s;
    for(int c = 0; c<LevelV::Ntest[0]; c++){
        for(int a = 0; a<LevelV::Nagg[0]; a++){
            n = a/2; //Lattice block
            s = a%2; //spin
            std::cout << "w[" << c << "][" << a << "] " << w[c][a] << "  " <<  wNew[n][2*c+s] << std::endl;
            w[c][a] = RandomU1();
            wNew[n][2*c+s] = w[c][a];
        }
    }

        std::cout << " ------------------ " << std::endl;
    /*
    Level0.P_v(wNew,vNew);
    testAMG.P_v(w,v);
    for(int i = 0; i<LevelV::Nsites[0]; i++){
    for(int j = 0; j<LevelV::DOF[0]; j++){
        std::cout << "v[" << i << "][" << j << "] " << v[i][j] << "  " << vNew[i][j] << std::endl;
    }
    }
    */

    //Testing implementation for the two levels case
    //For generating the coarse gauge links we do 
    Level0.makeCoarseLinks(Level1); //Generates the coarse gauge links of level 1
    //Now Level1.D_operator is defined
    spinor in(LevelV::Ntest[0],c_vector(LevelV::Nagg[0],1));
    spinor out(LevelV::Ntest[0],c_vector(LevelV::Nagg[0],0));
    spinor IN(LevelV::Nsites[1],c_vector(LevelV::DOF[1],1));
    spinor OUT(LevelV::Nsites[1],c_vector(LevelV::DOF[1],0));
    
    testAMG.Pt_D_P(in,out);
    Level1.D_operator(IN,OUT);

    
    for(int c = 0; c<LevelV::Ntest[0]; c++){
        for(int a = 0; a<LevelV::Nagg[0]; a++){
            n = a/2; //Lattice block
            s = a%2; //spin
            std::cout << "out[" << c << "][" << a << "] " << out[c][a] << "  " <<  OUT[n][2*c+s] << std::endl;
        }
    }





    //Level1.makeCoarseLinks(Level2); //Uses the coarse gauge links generated on the previous level to generate 
    //the coarse gauge links of the next level
    //Now Level2.D_operator is defined
    
    
    
   
    
    return 0;
}