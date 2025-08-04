#include "conjugate_gradient.h"

//Conjugate gradient for computing (DD^dagger)^-1 phi, where phi is a vector represented by a matrix
//phi[Ntot][2]
int conjugate_gradient(const c_matrix& U, const spinor& phi, spinor& x,const double& m0) {
    int k = 0; //Iteration number
    double err;
    double err_sqr;

  
    spinor r(LV::Ntot, c_vector(2, 0));  //r[coordinate][spin] residual
    spinor d(LV::Ntot, c_vector(2, 0)); //search direction
    spinor Ad(LV::Ntot, c_vector(2, 0)); //DD^dagger*d

    c_double alpha, beta;

	x = phi;
    D_D_dagger_phi(U, x, Ad, m0); //DD^dagger*x
    for(int n = 0; n<LV::Ntot; n++){
        r[n][0] = phi[n][0] - Ad[n][0];
        r[n][1] = phi[n][1] - Ad[n][1];
    }

    d = r; //initial search direction
    c_double r_norm2 = dot(r, r);
    double phi_norm2 = sqrt(std::real(dot(phi, phi)));

    while (k<CG::max_iter) {
        D_D_dagger_phi(U, d,Ad, m0); //DD^dagger*d 
        alpha = r_norm2 / dot(d, Ad); //alpha = (r_i,r_i)/(d_i,Ad_i)

        //x = x + alpha * d; //x_{i+1} = x_i + alpha*d_i
        for(int n = 0; n<LV::Ntot; n++){ 
            x[n][0] += alpha*d[n][0];
            x[n][1] += alpha*d[n][1];
        }
        
        //r = r - alpha * Ad; //r_{i+1} = r_i - alpha*Ad_i
        for(int n = 0; n<LV::Ntot; n++){
            r[n][0] -= alpha*Ad[n][0];
            r[n][1] -= alpha*Ad[n][1];
        }
        
        err_sqr = std::real(dot(r, r)); //err_sqr = (r_{i+1},r_{i+1})
		err = sqrt(err_sqr); // err = sqrt(err_sqr)
        if (err < CG::tol*phi_norm2) {
            std::cout << "Converged in " << k << " iterations" << " Error " << err << std::endl;
            return 1;
        }

        beta = err_sqr / r_norm2; //beta = (r_{i+1},r_{i+1})/(r_i,r_i)

        //d = r + beta * d; //d_{i+1} = r_{i+1} + beta*d_i 
        for(int n = 0; n<LV::Ntot; n++){
            d[n][0] *= beta; 
            d[n][1] *= beta;
            d[n][0] += r[n][0];
            d[n][1] += r[n][1];
        }
               
        r_norm2 = err_sqr;
        k++;
    }
    std::cout << "CG did not converge in " << CG::max_iter << " iterations" << " Error " << err << std::endl;
    return 0;
}