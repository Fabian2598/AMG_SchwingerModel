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

    AMG amg(GConf, m0,AMGV::nu1, AMGV::nu2); 
	amg.setUpPhase(1,AMGV::Nit);
    amg.initializeCoarseLinks();
    
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
    //Level2.makeAggregates(); -->For the coarsest level this is not necessary
    //std::cout << "Aggregates Level 2 built " << std::endl;

    Level0.makeBlocks();
    Level1.makeBlocks();
    //Level2.makeBlocks(); //--> Not necessary for the coarsest level
    std::cout << "Blocks built " << std::endl;

    
    //Checking operator at level = 0 
    spinor v(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    for(int n = 0; n < LevelV::Nsites[0]; n++){
        v[n][0] = RandomU1(); v[n][1] = RandomU1();
    }
    spinor out(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    spinor outv2(LevelV::Nsites[0],c_vector(LevelV::DOF[0],1));
    std::cout << "computing Dirac " << std::endl;
    D_phi(GConf.Conf,v,out,mass::m0);
    Level0.D_operator(v,outv2);

    for(int n = 0; n < LevelV::Nsites[0]; n++){
    for(int mu : {0,1}){
        if (std::abs(out[n][mu]-outv2[n][mu]) > 1e-8 ){
            std::cout << "[" << n << "][" << mu << "] " << "for level 0  different" << std::endl; 
            std::cout << out[n][mu] << "   /=    " << outv2[n][mu] << std::endl;
            return 1;
        }
    }
    }
    std::cout << "test passed for level 0 " << out[0][0] << "   =    " << outv2[0][0] << "\n" << std::endl;
    //---------------------------------------//
    //Checking operator at level = 1
    std::cout << "Set up and orthonormalization ... " << std::endl;
    Level0.setUp(); //Build test vectors 
    std::cout << "done " << std::endl;
    Level0.makeCoarseLinks(Level1);


    Level1.setUp(); //Build test vectors for level 1
    Level1.makeCoarseLinks(Level2);
    
    spinor vL1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1));
    spinor v1(AMGV::Ntest,c_vector(AMGV::Nagg,1));
    for(int n = 0; n < LevelV::Nsites[level1]; n++){
        for(int m; m<LevelV::DOF[level1]; m++){
            vL1[n][m] = 1;//RandomU1();
        }  
    }

    spinor outL1(AMGV::Ntest,c_vector(AMGV::Nagg,1));
    spinor outv2L1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1));
    amg.Pt_D_P(v1,outL1);
    
    Level1.D_operator(vL1,outv2L1);

    
    for(int ntest = 0; ntest<AMGV::Ntest; ntest++){
        for(int nagg = 0; nagg<AMGV::Nagg; nagg++){
            int x = nagg / 2; //Lattice block
            int alf = nagg % 2; //spin
            if (std::abs(outL1[ntest][nagg]-outv2L1[x][2*ntest+alf]) > 1e-8 ){
            std::cout << "[" << ntest << "][" << nagg << "] " << "for level 1  different" << std::endl; 
            std::cout << outL1[ntest][nagg] << "   /=    " << outv2L1[x][2*ntest+alf] << std::endl;
            return 1;
            }
        }
    }
    


    std::cout << "Test passed for level 1 " << outL1[0][0] << "   ==    " << outv2L1[0][0] << std::endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Testing that P^dag D P = D_c level by level 

    

    //---Testing level 1
    spinor V1(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],1));
    spinor temp(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor Dphi(LevelV::Nsites[level0],c_vector(LevelV::DOF[level0],0));
    spinor outv3(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    Level0.P_v(V1,temp);
    //std::cout << "All good" << std::endl;
    Level0.D_operator(temp,Dphi);
    Level0.Pt_v(Dphi,outv3);
    std::cout << "Testing level 1" << std::endl;
    for(int x = 0; x<LevelV::Nsites[level1]; x++){
        for(int dof = 0; dof<LevelV::DOF[level1]; dof++){
            if (std::abs(outv3[x][dof]-outv2L1[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 1  different" << std::endl; 
            std::cout << outv3[x][dof] << "   /=    " << outv2L1[x][dof] << std::endl;
            return 1;
            }
        }
    }
    
    std::cout << "P^dag D P coincides with Dc for level 1" << std::endl;
    std::cout << outv3[0][0] << "   =    " << outv2L1[0][0] << std::endl;

    //----This one does not coincide----//
    //----Testing level 2
    std::cout << "Testing level 2" << std::endl;
    spinor V2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],1));
    spinor temp2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor Dphi2(LevelV::Nsites[level1],c_vector(LevelV::DOF[level1],0));
    spinor outv4(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    spinor outv2L2(LevelV::Nsites[level2],c_vector(LevelV::DOF[level2],0));
    Level1.P_v(V2,temp2);
    //std::cout << "All good" << std::endl;
    Level1.D_operator(temp2,Dphi2);
    Level1.Pt_v(Dphi2,outv4);
    
    Level2.D_operator(V2,outv2L2);
    for(int x = 0; x<LevelV::Nsites[level2]; x++){
        for(int dof = 0; dof<LevelV::DOF[level2]; dof++){
            if (std::abs(outv4[x][dof]-outv2L2[x][dof]) > 1e-8 ){
            std::cout << "[" << x << "][" << dof << "] " << "for level 2  different" << std::endl; 
            std::cout << outv4[x][dof] << "   /=    " << outv2L2[x][dof] << std::endl;
            return 1;
            }
        }
    }
   
    
    return 0;
}


 /*
    int x=0, alf=0, bet=0, p=0 ,s=0;
    for(x = 0; x < Level1.Nsites; x++){
    for(alf = 0; alf < 2; alf++){
    for(bet = 0; bet < 2; bet++){
    for(p = 0; p < Level1.colors; p++){
    for(s = 0; s < Level1.colors; s++){
        int indx = Level1.getG1index(x,alf,bet,p,s);
        int indx2 = Level1.getG2G3index(x,alf,bet,p,s,0);
        if (std::abs(amg.A_coeff[x][alf][bet][p][s] -  Level1.G1[indx]) > 1e-8){
            std::cout << "x " << x << " alf " << alf << " bet " << bet << " p " << p << " s " << s << std::endl;
            std::cout << "A Coefficients different " << std::endl;
            exit(1);
        }

    }
    }
    }
    }
    }

    std::cout << "A coefficients checked for level 1 " << std::endl;
    for(x = 0; x < Level1.Nsites; x++){
    for(alf = 0; alf < 2; alf++){
    for(bet = 0; bet < 2; bet++){
    for(p = 0; p < Level1.colors; p++){
    for(s = 0; s < Level1.colors; s++){
        for(int mu : {0,1}){
        int indx2 = Level1.getG2G3index(x,alf,bet,p,s,mu);
        if (std::abs(amg.B_coeff[x][alf][bet][p][s][mu] -  Level1.G2[indx2]) > 1e-8){
            std::cout << "x " << x << " alf " << alf << " bet " << bet << " p " << p << " s " << s << " mu " << mu << std::endl;
            std::cout << "B Coefficients different " << std::endl;
            exit(1);
        }
        }

    }
    }
    }
    }
    }

    std::cout << "B coefficients checked for level 1 " << std::endl;

    for(x = 0; x < Level1.Nsites; x++){
    for(alf = 0; alf < 2; alf++){
    for(bet = 0; bet < 2; bet++){
    for(p = 0; p < Level1.colors; p++){
    for(s = 0; s < Level1.colors; s++){
        for(int mu : {0,1}){
        int indx2 = Level1.getG2G3index(x,alf,bet,p,s,mu);
        if (std::abs(amg.C_coeff[x][alf][bet][p][s][mu] -  Level1.G3[indx2]) > 1e-8){
            std::cout << "x " << x << " alf " << alf << " bet " << bet << " p " << p << " s " << s << " mu " << mu << std::endl;
            std::cout << "C Coefficients different " << std::endl;
            exit(1);
        }
        }

    }
    }
    }
    }
    }

    std::cout << "C coefficients checked for level 1 " << std::endl;
    */
    /*
    std::cout << "amg 1" << amg.A_coeff[x][alf][bet][p][s] << std::endl;
    std::cout << "g1 " << Level1.G1[indx] << std::endl;
    std::cout << "amg 2" << amg.B_coeff[x][alf][bet][p][s][0] << std::endl;
    std::cout << "g2 " << Level1.G2[indx2] << std::endl;
    std::cout << "amg 3" << amg.C_coeff[x][alf][bet][p][s][0] << std::endl;
    std::cout << "g3 " << Level1.G3[indx2] << std::endl;
    */