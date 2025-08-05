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
   
    using namespace SAPV;
    MPI_Init(&argc, &argv);
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    readParameters("../parameters.dat");
    srand(19);

    //srand(time(0));
    
    Coordinates(); //Builds array with coordinates of the lattice points x * Nt + t
    MakeBlocks(); //Makes lattice blocks 
    periodic_boundary(); //Builds LeftPB and RightPB (periodic boundary for U_mu(n))
    
    //double m0 = -0.5;
    double m0 = -0.18840579710144945;

    //Parameters in variables.cpp
    if (rank == 0){
        std::cout << "******************* Two-grid method for the Dirac matrix in the Schwinger model *******************" << std::endl;
        std::cout << " Nx = " << Nx << " Nt = " << Nt << std::endl;
        std::cout << " Lattice dimension = " << (Nx * Nt) << std::endl;
        std::cout << " Number of entries of the Dirac matrix = (" << (2 * Nx * Nt) << ")^2 = " << (2 * Nx * Nt) * (2 * Nx * Nt) << std::endl;
        std::cout << " Bare mass parameter m0 = " << m0 << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << " Lattice blocking for the aggregates" << std::endl;
        std::cout << "| block_x = " << block_x << " block_t = " << block_t << std::endl;
        std::cout << "| x_elements = " << x_elements << " t_elements = " << t_elements << std::endl;
        std::cout << "| Each aggregate has x_elements * t_elements = " <<  x_elements * t_elements << " elements" << std::endl;
        std::cout << "| Number of aggregates = " << AMGV::Nagg << std::endl;
    } 
    Aggregates(); //build aggregates
    CheckAggregates();
    if (rank == 0){
        std::cout << "----------------------------------" << std::endl;
        std::cout << " Variable blocking for SAP" << std::endl;
        std::cout << "| sap_block_x = " << sap_block_x << " sap_block_t = " << sap_block_t << std::endl;
        std::cout << "| Lattice sites in the x direction = " << sap_x_elements<< " and in the t direction = " << sap_t_elements<< std::endl;
        std::cout << "| Each Schwarz block has " <<  sap_lattice_sites_per_block << " lattice points and " << sap_variables_per_block << " variables" << std::endl;
        std::cout << "| D restricted to each block has (" << 2 * sap_lattice_sites_per_block << ")^2 = " << sap_variables_per_block*sap_variables_per_block << " entries" << std::endl;
        std::cout << "| Number of Schwarz blocks = " << N_sap_blocks << std::endl;
        std::cout << "| Red/Black blocks = " << sap_coloring_blocks << std::endl;
        std::cout << "| Number of processes = " << size << std::endl;
        std::cout << "| Number of blocks per process = " << sap_blocks_per_proc << std::endl;
        std::cout << "| GMRES restart length for SAP blocks = " << sap_gmres_restart_length << std::endl;
        std::cout << "| GMRES iterations for SAP blocks = " << sap_gmres_restarts << std::endl;
        std::cout << "| GMRES tolerance for SAP blocks = " << sap_gmres_tolerance << std::endl;
    }
    SchwarzBlocks(); //Builds the blocks for the Schwarz alternating method
    CheckBlocks(); //Check blocks dimensions

    if (rank == 0){
        std::cout << "----------------------------------" << std::endl;
        std::cout << " Two-grid parameters" << std::endl;
        std::cout << "| Number of test vectors = " << AMGV::Ntest << std::endl;
        std::cout << "| nu1 (pre-smoothing) = " << AMGV::nu1 << " nu2 (post-smoothing) = " << AMGV::nu2 << std::endl;
        std::cout << "| Number of iterations for improving the interpolator = " << AMGV::Nit << std::endl;
        std::cout << "| Coarse grid matrix (Dc) dimension = " << AMGV::Ntest * AMGV::Nagg << std::endl;
        std::cout << "| Number of SAP iterations to smooth test vectors = " << AMGV::SAP_test_vectors_iterations << std::endl; 
        std::cout << "| Restart length of GMRES at the coarse level = " << AMGV::gmres_restart_length_coarse_level << std::endl;
        std::cout << "| Restarts of GMRES at the coarse level = " << AMGV::gmres_restarts_coarse_level << std::endl;
        std::cout << "| GMRES tolerance for the coarse level solution = " << AMGV::gmres_tol_coarse_level << std::endl;
        std::cout << "----------------------------------" << std::endl;
        std::cout << " FGMRES with AMG preconditioning parameters" << std::endl;
        std::cout << "| FGMRES restart length = " << FGMRESV::fgmres_restart_length << std::endl;
        std::cout << "| FGMRES restarts = " << FGMRESV::fgmres_restarts << std::endl;
        std::cout << "| FGMRES tolerance = " << FGMRESV::fgmres_tolerance << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
    }
    
    GaugeConf GConf = GaugeConf(Nx, Nt);
    GConf.initialize(); //Initialize a random gauge configuration
    for(int l; l< AMGV::levels-1; l++){
        std::cout << "Level " << l << " Block X " << LevelV::BlocksX[l] 
        << " Block T " << LevelV::BlocksT[l] << " Ntest " << LevelV::Ntest[l] << " Nagg " << LevelV::Nagg[l]
        << std::endl;
    }
    //Open conf from file//
    
    /*
    {
        double beta = 2;
        int nconf = 20;
        std::ostringstream NameData;
        NameData << "../../confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
        //std::cout << "Reading conf from file: " << NameData.str() << std::endl;
        std::ifstream infile(NameData.str());
        if (!infile) {
            std::cerr << "File not found on rank " << rank << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int x, t, mu;
        double re, im;
        
        c_matrix CONF(Ntot,c_vector(2,0)); 
        while (infile >> x >> t >> mu >> re >> im) {
            CONF[Coords[x][t]][mu] = c_double(re, im); 
        }
        GConf.setGconf(CONF);
        infile.close();
        if (rank == 0){
            std::cout << "Conf read from " << NameData.str() << std::endl;
        }
    }
    */

    gmres_DB.set_params(GConf.Conf,m0); //Setting gauge conf and m0 for GMRES used in the Schwarz blocks

    
    spinor rhs(Ntot, c_vector(2, 0)); //random right hand side 
    spinor x(Ntot, c_vector(2, 0)); //solution vector 
    //Random right hand side
    //rhs[0][0] = 1.0;
    for(int i = 0; i < Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }
    
    clock_t start, end;
    double elapsed_time;
    double startT, endT;

   
   
    if (rank == 0){
        //Bi-cgstab inversion for comparison
        std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
        start = clock();
        spinor x0(Ntot, c_vector(2, 0)); //Initial guess
        int max_iter = 10000;//100000; //Maximum number of iterations
        spinor x_bi = bi_cgstab(&D_phi,Ntot,2,GConf.Conf, rhs, x0, m0, max_iter, 1e-10, true);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;  
        
        /*
        int len = AMGV::gmres_restart_length_coarse_level;
        int restarts = 1000; //If the restart length is too large this could be problematic ...

        spinor xgmres(Ntot,c_vector(2));
        FGMRES_fine_level fgmres_fine_level(Ntot, 2, len, restarts,1e-10,GConf.Conf, m0);
        start = clock();
        fgmres_fine_level.fgmres(rhs,x0,xgmres,true);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for GMRES = " << elapsed_time << " seconds" << std::endl; 
        */

        std::cout << "Inverting the normal equations with CG" << std::endl; 
        spinor xCG(Ntot,c_vector(2,0));
        start = clock();
        conjugate_gradient(GConf.Conf, rhs, xCG, m0);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for CG = " << elapsed_time << " seconds" << std::endl;  
        

    }
  
    

/*
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){std::cout << "--------------Flexible GMRES with SAP preconditioning version --------------" << std::endl;}   
    spinor xSAP(Ntot, c_vector(2, 0)); //Solution vector for SAP
    FGMRES_SAP fgmres_sap(Ntot, 2, FGMRESV::fgmres_restart_length, FGMRESV::fgmres_restarts,FGMRESV::fgmres_tolerance,GConf.Conf, m0);
    startT = MPI_Wtime();
    fgmres_sap.fgmres(rhs,x,xSAP,true);
    endT = MPI_Wtime();
    printf("[rank %d] time elapsed during the job NEW implementation: %.4fs.\n", rank, endT - startT);
*/

    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0){std::cout << "--------------Flexible GMRES with AMG preconditioning--------------" << std::endl;}
    spinor xAMG(Ntot, c_vector(2, 0)); //Solution vector for SAP
    startT = MPI_Wtime();
    FGMRES_two_grid fgmres_two_grid(Ntot, 2, FGMRESV::fgmres_restart_length, FGMRESV::fgmres_restarts,FGMRESV::fgmres_tolerance,GConf, m0);
    fgmres_two_grid.fgmres(rhs,x,xAMG,true);
    endT = MPI_Wtime();
    printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);
    printf("[MPI process %d] coarse time: %.4fs.\n", rank, coarse_time);
    printf("[MPI process %d] smooth time: %.4fs.\n", rank, smooth_time);
    printf("[MPI process %d] SAP time: %.4fs.\n", rank, SAP_time);
    

    MPI_Finalize();

    return 0;
}