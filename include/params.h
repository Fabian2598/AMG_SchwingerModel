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
        //std::cout << "Reading conf from file: " << NameData.str() << std::endl;
        std::ifstream infile(NameData.str());
        if (!infile) {
            std::cerr << "File not found" << std::endl;
            //MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int block_x, block_t, ntest;
        int level; 
        while (infile >> level >> block_x >> block_t >> ntest) {
            LBlocks::BlocksX[level] = block_x;
            LBlocks::BlocksT[level] = block_t;
            LBlocks::Ntest[level] = ntest;
            LBlocks::Nagg[level] = 2*block_x*block_t;
        }
        infile.close();
        std::cout << "Parameters read from " << NameData.str() << std::endl;
}