#ifndef TESTS_H_INCLUDED
#define TESTS_H_INCLUDED

#include <time.h> 
#include <ctime>
#include "bi_cgstab.h"
#include "conjugate_gradient.h"
#include "amg.h"
#include "mpi.h"

/*
    Class for testing the different methods
*/
class Tests{
public:

    Tests(const GaugeConf& GConf, const spinor& rhs, const spinor& x0 ,const double m0): 
    GConf(GConf), rhs(rhs), x0(x0), m0(m0){

    }

    void BiCG(spinor& x, const int max_it,const bool print);
    void GMRES(spinor& x, const int len, const int restarts,const bool print);
    void CG(spinor& x);
    void FGMRES_sap(spinor &x,const bool print);
    void SAP(spinor& x,const int iterations,const bool print);
    int fgmresAMG(spinor& x, const bool print);
    void multigrid(spinor& x, const bool print);
    
    void check_solution(const spinor& x_sol);

private:

    const GaugeConf GConf;
    const spinor rhs;
    const spinor x0;
    const double m0;
    clock_t start, end;
    double elapsed_time;
    double startT, endT;



};

#endif