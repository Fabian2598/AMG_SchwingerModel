#include "AMG.h"

void AlgebraicMG::setUpPhase(const double& eps, const int& Nit){
    
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	//Test vectors random initialization for each level (except the coarsest level)
	for(int l=0; l<AMGV::levels-1; l++){
		for (int i = 0; i < LevelV::Ntest[l]; i++) {
		for (int n = 0; n < LevelV::Nsites[l]; n++) {
		for (int dof = 0; dof < LevelV::DOF[l]; dof++) {
			levels[l]->test_vectors[i][n][dof] = eps * RandomU1();
		}
		}
		}
        levels[l]->interpolator_columns = levels[l]->test_vectors;
        levels[l]->orthonormalize(); 
	    levels[l]->makeCoarseLinks(*levels[l+1]); //Make coarse gauge links which define the operator D for the next level
	}
	
   
    
    //Smoothing the test vectors
    for(int l=0; l<AMGV::levels-1; l++){
        spinor rhs(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
		for (int i = 0; i < LevelV::Ntest[l]; i++) {
			//Right hand side of the linear system 
			rhs = levels[l]->test_vectors[i];  //I could also leave it as zero
            levels[l]->sap_l.SAP(rhs,levels[l]->test_vectors[i],AMGV::SAP_test_vectors_iterations,SAPV::sap_blocks_per_proc,false);
		}
		levels[l]->interpolator_columns = levels[l]->test_vectors; 
		levels[l]->orthonormalize(); 
		levels[l]->makeCoarseLinks(*levels[l+1]); 
	}

	//Adaptivity part
	
    if (rank == 0)std::cout << "Improving interpolator" << std::endl;
    
	for (int it = 0; it < Nit; it++) {
		if (rank == 0)std::cout << "****** Bootstrap iteration " << it << " ******" << std::endl;
		//for (int l = LevelV::maxLevel - 1; l=0; l--){
		for (int l = 0; l<AMGV::levels-1; l++){
			spinor rhs(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
			for (int i = 0; i < LevelV::Ntest[l]; i++) {
				rhs = levels[l]->test_vectors[i];
				v_cycle(l, rhs, levels[l]->test_vectors[i]);
			}
			//Build the interpolator between level l and l+1
			levels[l]->interpolator_columns = levels[l]->test_vectors; 
			levels[l]->orthonormalize(); 
			levels[l]->makeCoarseLinks(*levels[l+1]); //Make coarse gauge links which define the operator D for the next level
		}
	}
	
    if (rank == 0)std::cout << "Set-up phase finished" << std::endl;
	
	
    
	

}


void AlgebraicMG::testSetUp(){
	//Checking orthogonality
    for(int l = 0; l<AMGV::levels-1;l++){
        levels[l]->checkOrthogonality();
    }

     // Testing that P^dag D P = D_c level by level 
    for(int l = 0; l<AMGV::levels-1;l++){
        spinor in(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],1)); //in
        spinor temp(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0));
        spinor Dphi(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0));
        spinor out(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //out
        spinor out_v2(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //D_c
        //P^H D P
        levels[l]->P_v(in,temp);
        levels[l]->D_operator(temp,Dphi);
        levels[l]->Pt_v(Dphi,out);

        levels[l+1]->D_operator(in,out_v2);
        std::cout << "Testing level " << l+1 << std::endl;
        for(int x = 0; x<LevelV::Nsites[l+1]; x++){
            for(int dof = 0; dof<LevelV::DOF[l+1]; dof++){
                if (std::abs(out[x][dof]-out_v2[x][dof]) > 1e-8 ){
                std::cout << "[" << x << "][" << dof << "] " << "for level " << l+1 << " different" << std::endl; 
                std::cout << out[x][dof] << "   /=    " << out_v2[x][dof] << std::endl;
                return;
                }
            }
        }
        std::cout << "P^dag D P coincides with Dc for level " << l+1 << std::endl;
        std::cout << out[0][0] << "   =    " << out_v2[0][0] << std::endl;
    }
     
}


void AlgebraicMG::v_cycle(const int& l, const spinor& eta_l, spinor& psi_l){
	if (l == LevelV::maxLevel){
		//For the coarsest level we just use GMRES to find a solution
		levels[l]->gmres_l.fgmres(eta_l, eta_l, psi_l, false); //psi_l = D_l^-1 eta_l 
	}
	else{
		//Buffers
		spinor Dpsi(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); //D_l psi_l
		spinor r_l(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); //r_l = eta_l - D_l psi_l
		spinor eta_l_1(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //eta_{l+1}
		spinor psi_l_1(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //eta_{l+1}
		spinor P_psi(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0));  //P_l psi_{l+1}

		//Pre - smoothing
		if (AMGV::nu1 > 0)
			levels[l]->sap_l.SAP(eta_l,psi_l,AMGV::nu1,SAPV::sap_blocks_per_proc,false); 
		

		//Coarse grid correction 
		levels[l]->D_operator(psi_l,Dpsi); 
		for(int n = 0;n < LevelV::Nsites[l]; n++){
		for(int dof = 0; dof < LevelV::DOF[l]; dof++){
			r_l[n][dof] = eta_l[n][dof] - Dpsi[n][dof]; //r_l = eta_l - D_l psi_l
		}
		}
		levels[l]->Pt_v(r_l,eta_l_1); //eta_{l+1} = P^H (eta_l - D_l psi_l)
		this->v_cycle(l+1,eta_l_1,psi_l_1); //psi_{l+1} = V-Cycle(l+1,eta_{l+1})

		levels[l]->P_v(psi_l_1,P_psi); //P_psi = P_l psi_{l+1}
		for(int n = 0;n < LevelV::Nsites[l]; n++){
		for(int dof = 0; dof < LevelV::DOF[l]; dof++){
			psi_l[n][dof] += P_psi[n][dof]; //psi_l = psi_l + P_l psi_{l+1}
		}
		}

		//Post - smoothing
		if (AMGV::nu2 > 0)
			levels[l]->sap_l.SAP(eta_l,psi_l,AMGV::nu2,SAPV::sap_blocks_per_proc,false); 
		

	}	

}


void AlgebraicMG::k_cycle(const int& l, const spinor& eta_l, spinor& psi_l){
	if (l == LevelV::maxLevel){
		//For the coarsest level we just use GMRES to find a solution
		levels[l]->gmres_l.fgmres(eta_l, eta_l, psi_l, false); //psi_l = D_l^-1 eta_l 
	}
	else{
		//Buffers
		spinor Dpsi(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); //D_l psi_l
		spinor r_l(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); //r_l = eta_l - D_l psi_l
		spinor eta_l_1(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //eta_{l+1}
		spinor psi_l_1(LevelV::Nsites[l+1],c_vector(LevelV::DOF[l+1],0)); //eta_{l+1}
		spinor P_psi(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0));  //P_l psi_{l+1}

		//Pre - smoothing
		if (AMGV::nu1 > 0)
			levels[l]->sap_l.SAP(eta_l,psi_l,AMGV::nu1,SAPV::sap_blocks_per_proc,false); 
		

		//Coarse grid correction 
		levels[l]->D_operator(psi_l,Dpsi); 
		for(int n = 0;n < LevelV::Nsites[l]; n++){
		for(int dof = 0; dof < LevelV::DOF[l]; dof++){
			r_l[n][dof] = eta_l[n][dof] - Dpsi[n][dof]; //r_l = eta_l - D_l psi_l
		}
		}
		levels[l]->Pt_v(r_l,eta_l_1); //eta_{l+1} = P^H (eta_l - D_l psi_l)
		this->fgmres_k_cycle_l[l+1]->fgmres(eta_l_1,eta_l_1,psi_l_1,false);
		//psi_{l+1} = fgmres(l+1,eta_{l+1}) with K-cycle(l+1,rhs) as preconditioner

		levels[l]->P_v(psi_l_1,P_psi); //P_psi = P_l psi_{l+1}
		for(int n = 0;n < LevelV::Nsites[l]; n++){
		for(int dof = 0; dof < LevelV::DOF[l]; dof++){
			psi_l[n][dof] += P_psi[n][dof]; //psi_l = psi_l + P_l psi_{l+1}
		}
		}

		//Post - smoothing
		if (AMGV::nu2 > 0)
			levels[l]->sap_l.SAP(eta_l,psi_l,AMGV::nu2,SAPV::sap_blocks_per_proc,false); 
		

	}	

}




void AlgebraicMG::applyMultilevel(const int& it, const spinor&rhs, spinor& out,const double tol,const bool print_message){
	spinor r(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
	spinor Dx(LevelV::Nsites[0],c_vector(LevelV::DOF[0],0));
	double err;
	double norm = sqrt(std::real(dot(rhs, rhs)));

	//If cycle = 0 --> V-cycle
	if (AMGV::cycle == 0){
		for(int i = 0; i<it; i++){
			v_cycle(0, rhs, out);
			levels[0]->D_operator(out,Dx);
			for(int n = 0;n < LevelV::Nsites[0]; n++){
			for(int dof = 0; dof < LevelV::DOF[0]; dof++){
				r[n][dof] = rhs[n][dof] - Dx[n][dof];
			}
			}
		
			err = sqrt(std::real(dot(r, r)));
        	if (err < tol* norm) {
            	if (print_message == true) {
            		std::cout << "V-cycle converged in " << i+1 << " iterations" << " Error " << err << std::endl;
            	}
            	return ;
        	} 
		}
		if (print_message == true) 
        	std::cout << "V-cycle did not converge in " << it << " iterations" << " Error " << err << std::endl;
	}

	else if (AMGV::cycle == 1){
		for(int i = 0; i<it; i++){
			k_cycle(0, rhs, out);
			levels[0]->D_operator(out,Dx);
			for(int n = 0;n < LevelV::Nsites[0]; n++){
			for(int dof = 0; dof < LevelV::DOF[0]; dof++){
				r[n][dof] = rhs[n][dof] - Dx[n][dof];
			}
			}
		
			err = sqrt(std::real(dot(r, r)));
        	if (err < tol* norm) {
            	if (print_message == true) {
            		std::cout << "K-cycle converged in " << i+1 << " iterations" << " Error " << err << std::endl;
            	}
            	return ;
        	} 
		}
		if (print_message == true) 
        	std::cout << "K-cycle did not converge in " << it << " iterations" << " Error " << err << std::endl;
	}

	//If cycle = 1 --> K-cycle

}