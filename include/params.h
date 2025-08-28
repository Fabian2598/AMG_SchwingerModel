#include <string>
#include <fstream>
#include <sstream>
#include "variables.h"

/*
    Function for reading the AMG blocks and test vectors for each level
    The parameters file has the following information on each row
    level, block_x, block_t, ntest, sap_block_x, sap_block_t
*/
void readParameters(const std::string& inputFile){
    std::ostringstream NameData;
    NameData << inputFile;
    std::ifstream infile(NameData.str());
    if (!infile) {
        std::cerr << "File " << NameData.str() <<  " not found" << std::endl;
        //MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int block_x, block_t, ntest;
    int level; 
    int maxLevel = LevelV::maxLevel;
    int sap_block_x, sap_block_t;
    //read parameters for level < max_level
    while (infile >> level >> block_x >> block_t >> ntest >> sap_block_x >> sap_block_t) {
        LevelV::BlocksX[level] = block_x;
        LevelV::BlocksT[level] = block_t;
        LevelV::Ntest[level] = ntest;
        LevelV::Nagg[level] = 2*block_x*block_t;
        LevelV::NBlocks[level] = block_x * block_t;
        LevelV::Nsites[level] = (level == 0 ) ? LV::Ntot : LevelV::BlocksX[level-1] * LevelV::BlocksT[level-1];
        LevelV::NxSites[level] = (level == 0 ) ? LV::Nx : LevelV::BlocksX[level-1];
        LevelV::NtSites[level] = (level == 0 ) ? LV::Nt : LevelV::BlocksT[level-1];
		LevelV::DOF[level] = (level == 0 ) ? 2 : 2 * LevelV::Ntest[level-1];
        LevelV::Colors[level] = (level == 0 ) ? 1 : LevelV::Ntest[level-1]; 
        
        LevelV::SAP_Block_x[level] = sap_block_x;
        LevelV::SAP_Block_t[level] = sap_block_t;
        LevelV::SAP_elements_x[level] = LevelV::NxSites[level]/sap_block_x; 
        LevelV::SAP_elements_t[level] = LevelV::NtSites[level]/sap_block_t; 
        LevelV::SAP_variables_per_block[level] = LevelV::DOF[level] * LevelV::SAP_elements_x[level] * LevelV::SAP_elements_t[level] ; 

        LevelV::GMRES_restart_len[level] = 20;
        LevelV::GMRES_restarts[level] = 20;
        LevelV::GMRES_tol[level] = 0.1;

    }
    //Store the number of sites and degrees of freedom for the coarsest lattice as well
    LevelV::Nsites[maxLevel] =   LevelV::BlocksX[maxLevel-1] * LevelV::BlocksT[maxLevel-1];
	LevelV::DOF[maxLevel] =  2 * LevelV::Ntest[maxLevel-1];
    LevelV::NxSites[maxLevel] = LevelV::BlocksX[maxLevel-1];
    LevelV::NtSites[maxLevel] = LevelV::BlocksT[maxLevel-1];
    LevelV::Colors[maxLevel] = LevelV::Ntest[maxLevel-1];
    LevelV::SAP_Block_x[maxLevel] = 1; //For the coarsest level we don't need a blocking, but it is necessary to define some numbers here.
    LevelV::SAP_Block_t[maxLevel] = 1;

    LevelV::GMRES_restart_len[maxLevel] = 20;
    LevelV::GMRES_restarts[maxLevel] = 20;
    LevelV::GMRES_tol[maxLevel] = 0.1;

    infile.close();
    std::cout << "Parameters read from " << NameData.str() << std::endl;
}