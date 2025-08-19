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
            << " Number of lattice blocks " << NBlocks[l] << std::endl;
        }
        for(int l=0; l< AMGV::levels; l++){
            std::cout << "Level " << l << " Nsites " << Nsites[l] 
            << " Nxsites " << NxSites[l] << " NtSites " << NtSites[l] << " DOF " << DOF[l]
            << " Colors " << Colors[l] << std::endl;
        }
        std::cout << "\n";
    }
    Aggregates();

    //AMG amg(GConf, m0,AMGV::nu1, AMGV::nu2); 
	//amg.setUpPhase(1,AMGV::Nit);
    //amg.initializeCoarseLinks();
    
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
            return 1;
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
            return 1;
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
            return 1;
            }
        }
    }

    std::cout << "P^dag D P coincides with Dc for level 3" << std::endl;
    std::cout << out3[0][0] << "   =    " << out3_v2[0][0] << std::endl;
    
    
    return 0;
}