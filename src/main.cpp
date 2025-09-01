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

    AMGV::cycle = 0; //K-cycle = 1, V-cycle = 0
    AMGV::Nit = 1;
    mass::m0 = -0.18840579710144945;
    double m0 = mass::m0; 


    //Open conf from file//
    GaugeConf GConf = GaugeConf(LV::Nx, LV::Nt);
    GConf.initialize();
    double beta = 2;
    int nconf = 3;
    if (LV::Nx == 128)
        nconf = 3;
    else if (LV::Nx == 256)
        nconf = 20;
    else if (LV::Nx == 64)
        nconf = 0;
    {
        std::ostringstream NameData;
        //
        NameData << "../../SchwingerModel/fermions/SchwingerModel/confs/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
        //NameData << "C:/Users/jafan/Downloads/Conf/b" << beta << "_" << LV::Nx << "x" << LV::Nt << "/m-018/2D_U1_Ns" << LV::Nx << "_Nt" << LV::Nt << "_b" << 
        //format(beta).c_str() << "_m" << format(m0).c_str() << "_" << nconf << ".ctxt";
        //std::cout << "Reading conf from file: " << NameData.str() << std::endl;
        std::ifstream infile(NameData.str());
        if (!infile) {
            std::cerr << "File " << NameData.str() <<" not found on rank " << rank << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int x, t, mu;
        double re, im;
        
        c_matrix CONF(LV::Ntot,c_vector(2,0)); 
        while (infile >> x >> t >> mu >> re >> im) {
            CONF[Coords[x][t]][mu] = c_double(re, im); 
        }
        GConf.setGconf(CONF);
        infile.close();
        if (rank == 0){
            std::cout << "Conf read from " << NameData.str() << std::endl;
        }
    }
    

    MPI_Barrier(MPI_COMM_WORLD);
    //Parameters in variables.cpp
    if (rank == 0)
        printParameters();
     
    spinor rhs(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    //rhs[0][0] = 1;
    for(int i = 0; i < LV::Ntot; i++) {
        rhs[i][0] = RandomU1();
        rhs[i][1] = RandomU1();
    }

    spinor x0(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0)); //Intial guesss

      
    if (rank == 0){
        double elapsed_time;
        double start, end;
        std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
        start = clock();
        int max_iter = 10000;//100000; //Maximum number of iterations
        spinor x_bi = bi_cgstab(&D_phi,LV::Ntot,2,GConf.Conf, rhs, x0, m0, max_iter, 1e-10, true);
        end = clock();
        elapsed_time = double(end - start) / CLOCKS_PER_SEC;
        std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;  
    }
    
    
    spinor xAMG(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    

    MPI_Barrier(MPI_COMM_WORLD);

    double elapsed_time;
    double startT, endT;
   
    
    startT = MPI_Wtime();
    FGMRES_amg fgmres_amg(LV::Ntot, 2,  FGMRESV::fgmres_restart_length, FGMRESV::fgmres_restarts,FGMRESV::fgmres_tolerance,GConf, m0);
    fgmres_amg.fgmres(rhs,x0,xAMG,true);
    endT = MPI_Wtime();
    printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);
    
    
    
    //Stand alone solver
    /*
    {
        startT = MPI_Wtime();
        spinor x(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
        AlgebraicMG AMG(GConf, m0,AMGV::nu1, AMGV::nu2);
        AMG.setUpPhase(1,AMGV::Nit);
        MPI_Barrier(MPI_COMM_WORLD);
        //AMG.testSetUp();
        AMG.applyMultilevel(100, rhs,x,1e-10,true);
        endT = MPI_Wtime();
        printf("[MPI process %d] time elapsed during the job: %.4fs.\n", rank, endT - startT);
    }
    */
    
    //Check if D_phi and rhs are equal
    /*
    spinor Dphi(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
    D_phi(GConf.Conf,xAMG,Dphi,m0);
    for(int n=0; n<LevelV::Nsites[0]; n++){
        for(int dof=0; dof<LevelV::DOF[0]; dof++){
            if (std::abs(Dphi[n][dof]-rhs[n][dof])>1e-8){
                std::cout << "Error: D_phi and rhs are not equal at site " << n << " component " << dof << std::endl;
                std::cout << "D_phi: " << Dphi[n][dof] << " rhs: " << rhs[n][dof] << std::endl;
                exit(1);
            }
        }
    }
    */
    

    MPI_Finalize();

    return 0;
}


