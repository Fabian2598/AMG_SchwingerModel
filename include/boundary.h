#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "variables.h"

/*
	Modulo operation
*/
inline int mod(int a, int b) {
	int r = a % b;
	return r < 0 ? r + b : r;
}

/*
	Boundary conditions for all the levels
	The function is only called once at the beginning of the program.

	right periodic boundary x+hat{mu}
	left periodic boundary x-hat{mu}
	hat_mu[0] = { 1, 0 } --> hat_t
	hat_mu[1] = { 0, 1 } --> hat_x

    CALL ONCE THE BLOCK PARAMETERS ARE READ
*/
inline void boundary() {
	//unit vectors in the "mu" direction
	//mu = 0 -> time, mu = 1 -> space
    using namespace LevelV;
    std::vector<std::vector<int>>hat_mu(2, std::vector<int>(2, 0));
	hat_mu[0] = { 1, 0 }; //hat_t
	hat_mu[1] = { 0, 1 }; //hat_x
    int blockID;
    for(int l = 0; l < AMGV::levels; l++){
	    for (int x = 0; x < NxSites[l]; x++) {
		for (int t = 0; t < NtSites[l]; t++) {
		for (int mu = 0; mu < 2; mu++) {
            blockID = x * NtSites[l] + t;
			RightPB_l[l][blockID][mu] = mod(x + hat_mu[mu][1], NxSites[l]) * NtSites[l] + mod(t + hat_mu[mu][0], NtSites[l]);
			LeftPB_l[l][blockID][mu] = mod(x - hat_mu[mu][1], NxSites[l]) * NtSites[l] + mod(t - hat_mu[mu][0], NtSites[l]);
			SignR_l[l][blockID][mu] = (mu == 0 && t  == NtSites[l] - 1) ? -1 : 1; //sign for the right boundary in time
			SignL_l[l][blockID][mu] = (mu == 0 && t == 0) ? -1 : 1; //sign for the left boundary in time
		}
		}
	    }
    }

}


#endif