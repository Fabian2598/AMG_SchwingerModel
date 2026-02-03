#include "tests.h"

void Tests::BiCG(spinor& x,const int max_it, const bool print){   
    const bool save = false;
    std::cout << "--------------Bi-CGstab inversion--------------" << std::endl;
    start = clock();
    FLOPS = 0;
    x = bi_cgstab(&D_phi,LV::Ntot,2,GConf.Conf, rhs, x0, m0, max_it, 1e-10, print,save);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for Bi-CGstab = " << elapsed_time << " seconds" << std::endl;
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
}

void Tests::GMRES(spinor& x, const int len, const int restarts,const bool print){
    const bool save = false;
    std::cout << "--------------GMRES without preconditioning--------------" << std::endl;
    FGMRES_fine_level fgmres_fine_level(LV::Ntot, 2, len, restarts,1e-10,GConf.Conf, m0);
    start = clock();
    FLOPS = 0;
    fgmres_fine_level.fgmres(rhs,x0,x,print,save);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for GMRES = " << elapsed_time << " seconds" << std::endl; 
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
}

void Tests::CG(spinor& x){
    std::cout << "--------Inverting the normal equations with CG----------" << std::endl; 
    start = clock();
    FLOPS = 0;
    conjugate_gradient(GConf.Conf, rhs, x, m0);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for CG = " << elapsed_time << " seconds" << std::endl;  
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
}

void Tests::FGMRES_sap(spinor& x, const bool print){
    const bool save = false;

    std::cout << "--------------Flexible GMRES with SAP preconditioning version --------------" << std::endl;
    start = clock();
    FLOPS = 0;
    FGMRES_SAP fgmres_sap(LV::Ntot, 2, FGMRESV::fgmres_restart_length, FGMRESV::fgmres_restarts,FGMRESV::fgmres_tolerance,GConf.Conf, m0);
    fgmres_sap.fgmres(rhs,x0,x,print,save);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    printf("time elapsed during FGMRES_SAP implementation: %.4fs.\n", elapsed_time);
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
    fflush(stdout);

}

void Tests::SAP(spinor& x,const int iterations, const bool print){
    std::cout << "--------------SAP as stand-alone solver --------------" << std::endl;
    start = clock();
    FLOPS = 0;
    sap.SAP(rhs,x,iterations, SAPV::sap_blocks_per_proc,print);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    printf("time elapsed during SAP implementation: %.4fs.\n", elapsed_time);
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
    fflush(stdout);

}

int Tests::fgmresAMG(spinor& x, const bool print){
    const bool save = false;
    std::cout << "--------------FGMRES with AMG --------------" << std::endl;
    int iter;
    start = clock();
    FLOPS = 0;
    FGMRES_AMG f_amg(LV::Ntot, 2,  FGMRESV::fgmres_restart_length, FGMRESV::fgmres_restarts,FGMRESV::fgmres_tolerance,GConf, m0);
    iter = f_amg.fgmres(rhs,x0,x,print,save);
    end = clock();
    elapsed_time = double(end - start) / CLOCKS_PER_SEC;
    printf("Time elapsed during the job: %.4fs.\n", elapsed_time);
    std::cout << "FLOPS = " << FLOPS << std::endl;
    printFLOPS(FLOPS);
    std::cout << "-----------------------------------------------------------\n" << std::endl;
    return iter;
}

void Tests::multigrid(spinor& x, const bool print){
    std::cout << "--------------Stand-alone AMG --------------" << std::endl;
    AlgebraicMG AMG(GConf, m0,AMGV::nu1, AMGV::nu2);
    AMG.setUpPhase(AMGV::Nit);
    //AMG.testSetUp();
    AMG.applyMultilevel(100, rhs,x,1e-10,true);
}

void Tests::check_solution(const spinor& x_sol){
    spinor xini(LV::Ntot, c_vector(2, 0)); //Initial guess
    D_phi(GConf.Conf, x_sol, xini, m0); //D_phi U x
    for(int i = 0; i< LV::Ntot; i++){
        if (std::abs(xini[i][0] - rhs[i][0]) > 1e-8 || std::abs(xini[i][1] - rhs[i][1]) > 1e-8) {
            std::cout << "Solution not correct at index " << i << ": " << xini[i][0] << " != " << rhs[i][0] << " or " << xini[i][1] << " != " << rhs[i][1] << std::endl;
            return ;
        }
    }
    std::cout << "Solution is correct" << std::endl;
}