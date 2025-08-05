#include <string>
#include <fstream>
#include <sstream>
#include "variables.h"

/*
    Function for reading the AMG blocks and test vectors for each level
    The parameters file has the following information on each row
    level, block_x, block_t, ntest
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
        //read parameters for level < max_level
        while (infile >> level >> block_x >> block_t >> ntest) {
            LevelV::BlocksX[level] = block_x;
            LevelV::BlocksT[level] = block_t;
            LevelV::Ntest[level] = ntest;
            LevelV::Nagg[level] = 2*block_x*block_t;
            LevelV::NBlocks[level] = block_x * block_t;
            LevelV::Nsites[level] = (level == 0 ) ? LV::Ntot : LevelV::BlocksX[level-1] * LevelV::BlocksT[level-1];
            LevelV::NxSites[level] = (level == 0 ) ? LV::Nx : LevelV::BlocksX[level-1];
            LevelV::NtSites[level] = (level == 0 ) ? LV::Nt : LevelV::BlocksT[level-1];
		    LevelV::DOF[level] = (level == 0 ) ? 2 : 2 * LevelV::Ntest[level-1];
        }
        //Store the number of sites and degrees of freedom for the coarsest lattice as well
        LevelV::Nsites[level] = (level == 0 ) ? LV::Ntot : LevelV::NxSites[level-1] * LevelV::NtSites[level-1];
		LevelV::DOF[level] = (level == 0 ) ? 2 : 2 * LevelV::Ntest[level-1];
        LevelV::NxSites[level] = (level == 0 ) ? LV::Nx : LevelV::BlocksX[level-1];
        LevelV::NtSites[level] = (level == 0 ) ? LV::Nt : LevelV::BlocksT[level-1];

        infile.close();
        std::cout << "Parameters read from " << NameData.str() << std::endl;
}