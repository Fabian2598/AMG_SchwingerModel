#ifndef BI_CGSTAB_H
#define BI_CGSTAB_H

#include "dirac_operator.h"
#include "operator_overloads.h"
#include "variables.h"
#include <cmath>
#include <iostream>

/*
    Bi-cgstab method for comparing the two-grid method
    func: matrix-vector operation 
    dim1: dimension of the first index of the spinor 
  	dim2: dimension of the second index of the spinor 	
        ->The solution is of the form x[dim1][dim2]
    U: gauge configuration
    phi: right-hand side vector
    x0: initial guess vector
    m0: mass parameter for Dirac matrix 
    max_iter: maximum number of iterations
    tol: tolerance for convergence
    print_message: flag for printing convergence messages
        
    The convergence criterion is ||r|| < ||phi|| * tol
*/
spinor bi_cgstab(void (*func)(const c_matrix&, const spinor&, spinor&, const double&), const int& dim1, const int& dim2,
    const c_matrix& U, const c_matrix& phi, const c_matrix& x0, const double& m0, const int& max_iter, const double& tol, const bool& print_message);

#endif
