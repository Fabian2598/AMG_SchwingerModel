#include <time.h> 
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include "params.h" //Read parameters for lattice blocks, test vectors and SAP blocks
#include "bi_cgstab.h" //BiCGstab for comparison
#include "conjugate_gradient.h" //Conjugate gradient for inverting the normal equations
#include "boundary.h" //Build boundary conditions at every grid level
#include "amg.h" //Algebraic Multigrid Method
#include "mpi.h" //MPI
#include "tests.h" //Class for testing

#include <cstdint>
#include <cstring>

//mean of a vector
template <typename T>
double mean(std::vector<T> x){ 
    double prom = 0;
    for (T i : x) {
        prom += i*1.0;
    }   
    prom = prom / x.size();
    return prom;
}

template <typename T>
double standard_deviation(const std::vector<T>& data) {
    if (data.empty()) return 0.0;
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double sq_sum = 0.0;
    for (const auto& val : data) {
        sq_sum += (static_cast<double>(val) - mean) * (static_cast<double>(val) - mean);
    }
    return std::sqrt(sq_sum / data.size());
}

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

    
    readParameters("../parameters.dat");
    srand(19);
    //srand(time(0));
    
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t
    boundary(); //Boundaries for every level

    AMGV::cycle = 1; //K-cycle = 1, V-cycle = 0
    AMGV::Nit = 0;
    AMGV::SAP_test_vectors_iterations = 4;
    //-0.1023;//-0.0933;//-0.18840579710144945; //0.0709
    double m0 = mass::m0; 


    //Open conf from file//
    GaugeConf GConf = GaugeConf(LV::Nx, LV::Nt);
    GConf.initialize();

    double beta;
    int nconf;
    std::string confFile;
    std::string rhsFile;
    if (rank == 0){
         //---Input data---//
        std::cout << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cout << "beta : ";
        std::cin >> beta;
        std::cout << "m0: ";
        std::cin >> m0;
        std::cout << "Configuration id: ";
        std::cin >> nconf;
        std::cout << "Configuration file path: ";
        std::cin >> confFile;
        std::cout << "RHS file path: ";
        std::cin >> rhsFile;
        std::cout << " " << std::endl;
    }
   
    MPI_Bcast(&beta, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    MPI_Bcast(&nconf, 1, MPI_INT,  0, MPI_COMM_WORLD);
    mass::m0 = m0;

    int filename_len = 0;
    if (rank == 0) {
        filename_len = static_cast<int>(confFile.size()) + 1; // include null terminator
    }
    MPI_Bcast(&filename_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<char> filename_buf(filename_len);
    if (rank == 0) {
        std::memcpy(filename_buf.data(), confFile.c_str(), filename_len);
    }
    MPI_Bcast(filename_buf.data(), filename_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        confFile.assign(filename_buf.data());
    }

    filename_len = 0;
    if (rank == 0) {
        filename_len = static_cast<int>(rhsFile.size()) + 1; // include null terminator
    }
    MPI_Bcast(&filename_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<char> rhsname_buf(filename_len);
    if (rank == 0) {
        std::memcpy(rhsname_buf.data(), rhsFile.c_str(), filename_len);
    }
    MPI_Bcast(rhsname_buf.data(), filename_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        rhsFile.assign(rhsname_buf.data());
    }
           
    
    MPI_Barrier(MPI_COMM_WORLD);
    //Parameters in variables.cpp
    if (rank == 0){
        printParameters();
        std::cout << "Conf read from " << confFile << std::endl;
        std::cout << "rhs read from " << rhsFile << std::endl;
    }
    
    const spinor x0(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0)); //Intial guesss
    spinor rhs(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));

    GConf.readBinary(confFile);
    readBinaryRhs(rhs,rhsFile);
    



    //Solution buffers
    //spinor x_bi(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    //spinor x_cg(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    spinor xFAMG(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    //spinor xAMG(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));

    
    Tests test(GConf, rhs, x0 ,m0);
    //if (rank == 0){
        //test.BiCG(x_bi, 10000,true); //BiCGstab for comparison  
        //test.CG(x_cg); //Conjugate Gradient for inverting the normal equations
    //}

    MPI_Barrier(MPI_COMM_WORLD);
    test.fgmresAMG(xFAMG, true);



    MPI_Finalize();

    return 0;
}

//Four levels multigrid
//0 8 8 10 4 4
//1 4 4 10 4 4
//2 2 2 10 2 2

//Three levels multigrid
//0 8 8 10 4 4
//1 4 4 10 4 4

//or

//0 4 4 10 4 4
//1 2 2 10 2 2



/*
 //Reading Conf
    
    {
        std::ostringstream NameData;
        NameData << "../../SchwingerModel/fermions/SchwingerModel/confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-01023/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
               
        if (rank == 0)
            std::cout << "Reading conf from file: " << NameData.str() << std::endl;
        //GConf.read_conf(NameData.str());
        GConf.readBinary(NameData.str());
    }
*/

/*
    double Iter[3]; double exTime[3];
    double dIter[3]; double dexTime[3];
    const int Meas = 10;
    for(int i = 0; i < 3; i++){

    AMGV::Nit = 2*i;
    if (rank == 0) std::cout << "Number of iterations for improving the interpolator: " << AMGV::Nit << std::endl;

    
    std::vector<double> iterations(Meas,0);
    std::vector<double> times(Meas,0);
    if (rank == 0) std::cout << "--------------Flexible GMRES with AMG preconditioning--------------" << std::endl;

    for(int i = 0; i < Meas; i++){
        if (rank == 0) std::cout << "Meas " << i << std::endl;
        iterations[i] = test.fgmresAMG(xFAMG, true);
        times[i] = total_time;
    }
    if (rank == 0){
        std::cout << "Average iteration number over " << Meas << " runs: " << mean(iterations) << " +- " 
        << standard_deviation(iterations)/sqrt(1.0*Meas) << std::endl;
    
    }

    Iter[i] = mean(iterations);
    dIter[i] = standard_deviation(iterations)/sqrt(1.0*Meas);
    exTime[i] = mean(times);
    dexTime[i] = standard_deviation(times)/sqrt(1.0*Meas);

    }

    for(int i = 0; i < 3; i++){
        if (rank == 0) std::cout << "Nit: " << 2*i << " Iter: " << Iter[i] << " +- " << dIter[i] << std::endl;

    }
    if (rank == 0)
        saveParameters(Iter, dIter, exTime, dexTime, 3,nconf);
*/

    
    //test.multigrid(xAMG,true); //Multigrid as stand-alone solver
    //test.check_solution(xFAMG); //Check that the solution is correct


    /*
        std::ostringstream FileName;
    FileName << "../../SchwingerModel/fermions/SchwingerModel/confs/rhs/rhs_conf" << 
    nconf << "_" << LV::Nx << "x" << LV::Nt << "_b" << format(beta) << "_m" << format(m0) << ".rhs";
    
    //read_rhs(rhs,FileName.str());
    readBinaryRhs(rhs,FileName.str());
    //random_rhs(rhs,10);

    // Save rhs to a .txt file
    if (rank == 0){
        std::ostringstream FileName;
        FileName << "rhs_conf" << nconf << "_" << LV::Nx << "x" << LV::Nt
                 << "_b" << format(beta) << "_m" << format(m0)
                 << ".rhs";
      //  save_rhs(rhs,FileName.str());
    }
    */