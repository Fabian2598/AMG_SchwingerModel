#include "dirac_operator.h"

c_double I_number(0, 1); //imaginary number

void D_phi(const c_matrix& U, const spinor& phi, spinor &out,const double& m0) {
	using namespace LV;
	for (int n = 0; n < Ntot; n++) {
		//n = x * Nt + t
		out[n][0] = (m0 + 2) * phi[n][0] -0.5 * ( 
			U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] - phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (phi[RightPB[n][1]][0] + I_number * phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi[LeftPB[n][1]][0] - I_number * phi[LeftPB[n][1]][1])
		);

		out[n][1] = (m0 + 2) * phi[n][1] -0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi[RightPB[n][1]][0] + phi[RightPB[n][1]][1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi[LeftPB[n][1]][0] + phi[LeftPB[n][1]][1])
		);
			
	}
	
}


void D_dagger_phi(const c_matrix& U, const spinor& phi, spinor &Dphi,const double& m0) {
	using namespace LV;
	for (int n = 0; n < Ntot; n++) {
		//n = x * Nt + t
		Dphi[n][0] = (m0 + 2) * phi[n][0] -0.5 * ( 
			std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi[LeftPB[n][0]][0] - phi[LeftPB[n][0]][1]) 
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi[LeftPB[n][1]][0] + I_number * phi[LeftPB[n][1]][1])
		+   U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1])
		+   U[n][1] * SignR[n][1] * (phi[RightPB[n][1]][0] - I_number * phi[RightPB[n][1]][1])
		);

		Dphi[n][1] = (m0 + 2) * phi[n][1] -0.5 * ( 
			std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (-phi[LeftPB[n][0]][0] + phi[LeftPB[n][0]][1]) 
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (-I_number*phi[LeftPB[n][1]][0] + phi[LeftPB[n][1]][1])
		+   U[n][0] * SignR[n][0] * (phi[RightPB[n][0]][0] + phi[RightPB[n][0]][1])
		+   U[n][1] * SignR[n][1] * (I_number*phi[RightPB[n][1]][0] + phi[RightPB[n][1]][1])
		);
			
	}

}


void D_D_dagger_phi(const c_matrix& U, const spinor& phi, spinor &Dphi,const double& m0) {
	spinor temporary(LV::Ntot,c_vector(2,0));
	D_dagger_phi(U, phi, temporary, m0);
	D_phi(U,  temporary, Dphi, m0);
}


	