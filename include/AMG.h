#ifndef ALGEBRAICMG_H_INCLUDED
#define ALGEBRAICMG_H_INCLUDED

#include "level.h"

class AlgebraicMG{
    /*
	GaugeConf GConf: Gauge configuration
    m0: Mass parameter for the Dirac matrix
    nu1: Number of pre-smoothing iterations
    nu2: Number of post-smoothing iterations
    nlevels: number of levels
	*/
public:
    AlgebraicMG(const GaugeConf & GConf, const double& m0, const int& nu1, const int& nu2) 
	: GConf(GConf), m0(m0), nu1(nu1), nu2(nu2){

    for(int l = 0; l<AMGV::levels; l++){
        Level* level = new Level(l,GConf.Conf);
        levels.push_back(level);
    }

    for(int l = 0; l<AMGV::levels-1; l++){
        levels[l]->makeBlocks();
        levels[l]->makeAggregates();
    }    
    //std::cout << "Lattice blocks and aggregates initialized" << std::endl;
    	
    }


    

    //Pages 84 and 85 of Rottmann's thesis explain how to implement this ...
    void setUpPhase(const double& eps, const int& Nit);
//private:    
    GaugeConf GConf;
	double m0; 
	int nu1, nu2; 
    std::vector<Level*> levels; //If i try to use a vector of objects I will run out of memory

    //Checks orthonormalization 
    void testSetUp();

    // psi_l = V_cycle(l,eta_l)
    void v_cycle(const int& l, const spinor& eta_l, spinor& psi_l);

    //Calls any cycle (for the moment only v)
    void applyMultilevel(const int& it, const spinor&rhs, spinor& out,const double tol,const bool print_message);
    

};

#endif