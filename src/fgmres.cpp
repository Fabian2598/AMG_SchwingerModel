#include "fgmres.h"
#include "iomanip"

//------Class FGMRES implementation------//
int FGMRES::fgmres(const spinor& phi, const spinor& x0, spinor& x,const bool& print_message) { 
    setZeros();
    int k = 0; //Iteration number (restart cycle)
    double err;
    x = x0; //initial solution. Perhaps it is better to give a reference to avoid a copy
    func(x, Dx); //Matrix-vector operation
    axpy(phi,Dx, -1.0, r); //r = b - A*x
	double norm_phi = sqrt(std::real(dot(phi, phi))); //norm of the right hand side
    err = sqrt(std::real(dot(r, r))); //Initial error
    
    while (k < restarts) {
        beta = err + 0.0 * I_number;
        scal(1.0/beta, r,VmT[0]); //VmT[0] = r / ||r||
        gm[0] = beta; //gm[0] = ||r||
        //-----Arnoldi process to build the Krylov basis and the Hessenberg matrix-----//
        for (int j = 0; j < m; j++) {
            preconditioner(VmT[j], ZmT[j]); //ZmT[j] = M^-1 VmT[j]
        
            func(ZmT[j],w); 
            //Gram-Schmidt process to orthogonalize the vectors
            for (int i = 0; i <= j; i++) {
                Hm[i][j] = dot(w, VmT[i]); //  (v_i^dagger, w)
                //w = w -  Hm[i][j] * VmT[i];
                for(int n=0; n<dim1; n++){
					for(int l=0; l<dim2; l++){
						w[n][l] -= Hm[i][j] * VmT[i][n][l];
					}
				}
            }

            Hm[j + 1][j] = sqrt(std::real(dot(w, w))); //H[j+1][j] = ||A v_j||
            if (std::real(Hm[j + 1][j]) > 0) {
                scal(1.0 / Hm[j + 1][j], w, VmT[j + 1]); //VmT[j + 1] = w / ||A v_j||
            }
            //----Rotate the matrix----//
            rotation(j);

            //Rotate gm
            gm[j + 1] = -sn[j] * gm[j];
            gm[j] = std::conj(cn[j]) * gm[j];
        }        
        //Solve the upper triangular system//
		solve_upper_triangular(Hm, gm,m,eta);
        
        for (int i = 0; i < dim1 * dim2; i++) {
            int n = i / dim2; int mu = i % dim2;
            for (int j = 0; j < m; j++) {
                x[n][mu] = x[n][mu] + eta[j] * ZmT[j][n][mu]; 
            }
        }
        //Compute the residual
        func(x, Dx);
        axpy(phi,Dx, -1.0, r); //r = b - A*x

        err = sqrt(std::real(dot(r, r)));
        if (err < tol* norm_phi) {
            if (print_message == true) {
             std::cout << "FGMRES converged in " << k + 1 << " iterations" << " Error " << err << std::endl;
            }
            return 1;
        }
        k++;
    }
    if (print_message == true) {
        std::cout << "FGMRES did not converge in " << restarts << " restarts" << " Error " << err << std::endl;
    }
    return 0;
}

void FGMRES::rotation(const int& j) {
    //Rotation of the column elements that are <j
    c_double temp;
    for (int i = 0; i < j; i++) {
		temp = std::conj(cn[i]) * Hm[i][j] + std::conj(sn[i]) * Hm[i + 1][j];
		Hm[i + 1][j] = -sn[i] * Hm[i][j] + cn[i] * Hm[i + 1][j];
		Hm[i][j] = temp;
    }
    //Rotation of the diagonal and element right below the diagonal
    c_double den = sqrt(std::conj(Hm[j][j] ) * Hm[j][j] + std::conj(Hm[j + 1][j]) * Hm[j + 1][j]);
	sn[j] = Hm[j + 1][j] / den; cn[j] = Hm[j][j] / den;
	Hm[j][j] = std::conj(cn[j]) * Hm[j][j] + std::conj(sn[j]) * Hm[j + 1][j];
    Hm[j + 1][j] = 0.0;

}

//x = A^-1 b, A an upper triangular matrix of dimension n
void FGMRES::solve_upper_triangular(const c_matrix& A, const c_vector& b, const int& n, c_vector& out) {
	for (int i = n - 1; i >= 0; i--) {
		out[i] = b[i];
		for (int j = i + 1; j < n; j++) {
			out[i] -= A[i][j] * out[j];
		}
		out[i] /= A[i][i];
	}
}