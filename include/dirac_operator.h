#ifndef DIRAC_OPERATOR_H
#define DIRAC_OPERATOR_H
#include "variables.h"
#include "operator_overloads.h"
#include "omp.h"

extern c_double I_number; //imaginary number

/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

/*
	Periodic boundary conditions used for the link variables U_mu(n).
	This function builds the arrays x_1_t1, x1_t_1, RightPB and LeftPB, which
	store the neighbor coordinates for the periodic boundary conditions.
	This prevents recalculation every time we call the operator D.
	The function is only called once at the beginning of the program.

	right periodic boundary x+hat{mu}
	left periodic boundary x-hat{mu}
	hat_mu[0] = { 1, 0 } --> hat_t
	hat_mu[1] = { 0, 1 } --> hat_x
*/
inline void periodic_boundary() {
	using namespace LV; //Lattice parameters namespace
	//unit vectors in the "mu" direction
	//mu = 0 -> time, mu = 1 -> space
	std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
	for (int x = 0; x < Nx; x++) {
		for (int t = 0; t < Nt; t++) {
			for (int mu = 0; mu < 2; mu++) {
				RightPB[Coords[x][t]][mu] = Coords[mod(x + hat_mu[mu][1], Nx)][mod(t + hat_mu[mu][0], Nt)]; 
				LeftPB[Coords[x][t]][mu] = Coords[mod(x - hat_mu[mu][1], Nx)][mod(t - hat_mu[mu][0], Nt)];
				SignR[Coords[x][t]][mu] = (mu == 0 && t == Nt - 1) ? -1 : 1; //sign for the right boundary in time
				SignL[Coords[x][t]][mu] = (mu == 0 && t == 0) ? -1 : 1; //sign for the left boundary in time
			}
		}
	}

	//Lattice neighbouring for the blocks
	int block;
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
			for (int mu = 0; mu < 2; mu++) {
				block = x * block_t + t;
				RightPB_blocks[block][mu] = mod(x + hat_mu[mu][1], block_x) * block_t + mod(t + hat_mu[mu][0], block_t);
				LeftPB_blocks[block][mu] = mod(x - hat_mu[mu][1], block_x) * block_t + mod(t - hat_mu[mu][0], block_t);
			}
		}
	}
}

/*
	Dirac operator application D phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
	out: spinor with the result of the operator application
*/
void D_phi(const c_matrix& U, const spinor& phi, spinor &out,const double& m0);


/*
	Dirac dagger operator application D^+ phi
	U: gauge configuration
	phi: spinor to apply the operator to
	m0: mass parameter
*/
void D_dagger_phi(const c_matrix& U, const spinor& phi, spinor &Dphi, const double& m0);


/*
	Application of D D^+
	It just calls the previous functions
*/
void D_D_dagger_phi(const c_matrix& U, const spinor& phi, spinor &Dphi,const double& m0);


#endif