#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "dirac_operator.h"
#include <cmath>
#include <iostream>

/*
    Conjugate gradient method for computing (DD^dagger)^-1 phi 
    U: gauge configuration
    phi: right-hand side vector
    x: output solution
    m0: mass parameter for Dirac matrix 

    The right-hand side is used as initial solution
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/
int conjugate_gradient(const c_matrix& U, const spinor& phi, spinor &x, const double& m0); 


#endif
