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

	//FGMRES for the k-cycle
    class FGMRES_k_cycle : public FGMRES {
	public:
    	FGMRES_k_cycle(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol, AlgebraicMG* parent, int l) : 
		FGMRES(dim1, dim2, m, restarts, tol), parent(parent), l(l) {}
    
    	~FGMRES_k_cycle() { };
    
	private:
		AlgebraicMG* parent; //Pointer to the enclosing AMG instance
		int l; //Level
    	/*
    	Implementation of the function that computes the matrix-vector product for the fine level
    	*/
    	void func(const spinor& in, spinor& out) override {
        	parent->levels[l]->D_operator(in,out);
    	}
		//Preconditioning with the k-cycle
		void preconditioner(const spinor& in, spinor& out) override {
            parent->k_cycle(l,in,out); 
		}
	};

    AlgebraicMG(const GaugeConf & GConf, const double& m0, const int& nu1, const int& nu2) 
	: GConf(GConf), m0(m0), nu1(nu1), nu2(nu2){

    	for(int l = 0; l<AMGV::levels; l++){
        	Level* level = new Level(l,GConf.Conf);
        	levels.push_back(level);
			//We don't really need this FGMRES for the coarsest level and the finest level
			FGMRES_k_cycle* fgmres = new FGMRES_k_cycle(LevelV::Nsites[l], 
				LevelV::DOF[l], 
				AMGV::fgmres_k_cycle_restart_length, 
				AMGV::fgmres_k_cycle_restarts, 
				AMGV::fgmres_k_cycle_tol, 
				this,
				l);
			fgmres_k_cycle_l.push_back(fgmres);
    	}

    	for(int l = 0; l<AMGV::levels-1; l++){
        	levels[l]->makeBlocks();
        	levels[l]->makeAggregates();
		}
    }    
    	

    //Pages 84 and 85 of Rottmann's thesis explain how to implement this ...
    void setUpPhase(const double& eps, const int& Nit);
//private:    
    GaugeConf GConf;
	double m0; 
	int nu1, nu2; 
    std::vector<Level*> levels; //If I try to use a vector of objects I will run out of memory
	std::vector<FGMRES_k_cycle*> fgmres_k_cycle_l;

    //Checks orthonormalization and that P^H D P = Dc
    void testSetUp();

    // psi_l = V_cycle(l,eta_l)
    void v_cycle(const int& l, const spinor& eta_l, spinor& psi_l);

	// psi_l = K_cycle(l,eta_l)
	void k_cycle(const int& l, const spinor& eta_l, spinor& psi_l);

    //Calls K or V-cycle depending on the value of AMGV::cycle
    void applyMultilevel(const int& it, const spinor&rhs, spinor& out,const double tol,const bool print_message);
    

};

#endif