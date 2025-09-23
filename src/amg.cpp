#include "amg.h"

void AlgebraicMG::setUpPhase(const int& Nit){
    
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	static std::mt19937 randomInt(50); //Same seed for all the MPI copies
	std::uniform_real_distribution<double> distribution(-1.0, 1.0); //mu, standard deviation
	

	//Generate test vectors at the fine level
	for (int i = 0; i < LevelV::Ntest[0]; i++) {
		for (int n = 0; n < LevelV::Nsites[0]; n++) {
		for (int dof = 0; dof < LevelV::DOF[0]; dof++) {
			levels[0]->interpolator_columns[i][n][dof] = distribution(randomInt) + I_number * distribution(randomInt);
		}
		}
	}
	
	//v_l = P^dagger v_{l-1}
	for(int l=1; l<AMGV::levels-1; l++){
		for(int i = 0; i<LevelV::Ntest[l];i++){
			if (i<LevelV::Ntest[l-1]){
				levels[l-1]->Pt_v(levels[l-1]->interpolator_columns[i],levels[l]->interpolator_columns[i]);
			}
			else{
				for (int n = 0; n < LevelV::Nsites[l]; n++) {
				for (int dof = 0; dof < LevelV::DOF[l]; dof++) {
					levels[l]->interpolator_columns[i][n][dof] =  distribution(randomInt) + I_number * distribution(randomInt);
				}	
				}
			}
		}
	}

	//Smoothing the test vectors
    for(int l=0; l<AMGV::levels-1; l++){
        spinor rhs(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
		for (int i = 0; i < LevelV::Ntest[l]; i++) {
			//Approximately solving D x = 0
            levels[l]->sap_l.SAP(rhs,levels[l]->interpolator_columns[i],AMGV::SAP_test_vectors_iterations,SAPV::sap_blocks_per_proc,false);
		}
		levels[l]->orthonormalize(); 
		levels[l]->makeCoarseLinks(*levels[l+1]); 
	}

	//Adaptivity part
	
    if (rank == 0)std::cout << "Improving interpolator" << std::endl;
    
	for (int it = 0; it < Nit; it++) {
		if (rank == 0)std::cout << "****** Bootstrap iteration " << it << " ******" << std::endl;
		for (int l = 0; l<AMGV::levels-1; l++){
			spinor rhs(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
			spinor Dv(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
			spinor zero(LevelV::Nsites[l], c_vector(LevelV::DOF[l],0));
			for (int i = 0; i < LevelV::Ntest[l]; i++) {
				levels[l]->test_vectors[i] = zero; 

				levels[l]->D_operator(levels[l]->interpolator_columns[i], Dv); //Dv = D v
				for(int n = 0; n < LevelV::Nsites[l]; n++) {
				for(int dof = 0; dof < LevelV::DOF[l]; dof++) {
					rhs[n][dof] = levels[l]->interpolator_columns[i][n][dof] - Dv[n][dof]; //rhs = v - D v
				}
				}

				if (AMGV::cycle == 0)
					v_cycle(l, rhs, levels[l]->test_vectors[i]);
				else if (AMGV::cycle == 1)
					k_cycle(l, rhs, levels[l]->test_vectors[i]);

				for(int n = 0; n < LevelV::Nsites[l]; n++) {
				for(int dof = 0; dof < LevelV::DOF[l]; dof++) {
					levels[l]->test_vectors[i][n][dof] += levels[l]->interpolator_columns[i][n][dof]; //v = v + Cycle(v-Dv)
				}
				}
			}
			//Build the interpolator between level l and l+1
			levels[l]->interpolator_columns = levels[l]->test_vectors; 
			levels[l]->orthonormalize(); 
			levels[l]->makeCoarseLinks(*levels[l+1]); //Make coarse gauge links which define the operator D for the next level
		}
	}
	
    if (rank == 0)std::cout << "Set-up phase finished" << std::endl;
	
}

void AlgebraicMG::v_cycle(const int& l, const spinor& eta_l, spinor& psi_l){
	if (l == LevelV::maxLevel){
		//For the coarsest level we just use GMRES to find a solution
		levels[l]->gmres_l.fgmres(eta_l, eta_l, psi_l, false,false); //psi_l = D_l^-1 eta_l 
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
		for(int n = 0; n < LevelV::Nsites[l]; n++){
		for(int dof = 0; dof < LevelV::DOF[l]; dof++){
			r_l[n][dof] = eta_l[n][dof] - Dpsi[n][dof]; //r_l = eta_l - D_l psi_l
		}
		}
		levels[l]->Pt_v(r_l,eta_l_1); //eta_{l+1} = P^H (eta_l - D_l psi_l)
		v_cycle(l+1,eta_l_1,psi_l_1); //psi_{l+1} = V-Cycle(l+1,eta_{l+1})

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
		levels[l]->gmres_l.fgmres(eta_l, eta_l, psi_l, false,false); //psi_l = D_l^-1 eta_l 
	}
	else{
		//Buffers might be useful to define them somewhere else and just set them to zero here
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
		fgmres_k_cycle_l[l+1]->fgmres(eta_l_1,eta_l_1,psi_l_1,false,false);
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
            		std::cout << "V-cycle converged in " << i+1 << " cycles" << " Error " << err << std::endl;
            	}
            	return ;
        	} 
		}
		if (print_message == true) 
        	std::cout << "V-cycle did not converge in " << it << " cycles" << " Error " << err << std::endl;
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
            		std::cout << "K-cycle converged in " << i+1 << " cycles" << " Error " << err << std::endl;
            	}
            	return ;
        	} 
		}
		if (print_message == true) 
        	std::cout << "K-cycle did not converge in " << it << " cycles" << " Error " << err << std::endl;
	}
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

void AlgebraicMG::testSAP(){
    int rank, size; 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(int l = 0; l<AMGV::levels-1;l++){
    	spinor rhs(LevelV::Nsites[l],c_vector(LevelV::DOF[l],1)); 
    	spinor x(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); 
		spinor xgmres(LevelV::Nsites[l],c_vector(LevelV::DOF[l],0)); 
    	int iter = 100; //SAP iterations
    	MPI_Barrier(MPI_COMM_WORLD);
		levels[l]->sap_l.SAP(rhs,x,iter,SAPV::sap_blocks_per_proc,true); 
    	MPI_Barrier(MPI_COMM_WORLD);
    	levels[l]->gmres_l.fgmres(rhs,xgmres,xgmres,true);

        for(int n=0; n<LevelV::Nsites[l]; n++){
        for(int alf=0; alf<LevelV::DOF[l]; alf++){
            if(std::abs(x[n][alf] - xgmres[n][alf]) > 1e-8){
                std::cout << "GMRES and SAP give something different at level " << l << " at site " << n << " and DOF " << alf << std::endl;
                std::cout << "xSAP = " << x[n][alf] << ", x = " << xgmres[n][alf] << std::endl;
                exit(1);
            }
        }
        }
        std::cout << "GMRES and SAP solution coincide at level " << l << std::endl;
    }

}