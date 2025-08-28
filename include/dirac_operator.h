#ifndef DIRAC_OPERATOR_H
#define DIRAC_OPERATOR_H
#include "variables.h"
#include "lin_alg_op.h"
#include "omp.h"

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