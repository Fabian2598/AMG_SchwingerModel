#include "amg.h"

//Aggregates A_j_0 = L_j x {0}, A_j_1 = L_j x {1} (Volume times Spin component)
void Aggregates() {
	using namespace LV; //Lattice parameters namespace
	for (int x = 0; x < block_x; x++) {
		for (int t = 0; t < block_t; t++) {
			for (int s = 0; s < 2; s++) {
				//vectorize the coordinates of these three loops
				int x0 = x * x_elements, t0 = t * t_elements;
				int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
				int aggregate = x * block_t * 2 + t * 2 + s;
				int count = 0;
				//x and t are redefined in the following loop
				for (int x = x0; x < x1; x++) {
					for (int t = t0; t < t1; t++) {
						int i = x * Nt * 2  + t * 2 + s;
						XCoord[i] = x; TCoord[i] = t; SCoord[i] = s;
						Agg[aggregate][count] = i;
						count++;
					}
				}
				if (count != x_elements * t_elements) {
					std::cout << "Aggregate " << aggregate << " has " << count << " elements" << std::endl;
				}
				//Once the loops are finished count should be x_elements*t_elements
			}	
		}
	}
	AMGV::aggregates_initialized = true; //Set aggregates as initialized
}


void PrintAggregates() {
	for (int i = 0; i < AMGV::Nagg; i++) {
		std::cout << "-------Aggregate-----" << i << std::endl;
		for (int j = 0; j < LV::x_elements * LV::t_elements; j++) {
			std::cout << Agg[i][j] << " ";
		}
		std::cout << std::endl;
	}
}




spinor canonical_vector(const int& i, const int& N1, const int& N2) {
	spinor e_i(N1, c_vector (N2,0.0));
	int j = i / N2;
	int k = i % N2;
	e_i[j][k] = 1.0;
	return e_i;
}


void AMG::orthonormalize(){
	/*
	Local orthonormalization of the test vectors
	
	Each test vector is chopped into the Nagg aggregates, which yields Ntest*Nagg columns for the interpolator.
	Each column is orthonormalized with respect to the others that belong to the same aggregate.
	This follows the steps from Section 3.1 of A. Frommer et al "Adaptive Aggregation-Based Domain Decomposition 
	Multigrid for the Lattice Wilson-Dirac Operator", SIAM, 36 (2014).
	*/

	using namespace AMGV; //AMG parameters

	//Getting the columns of the interpolator for the orthonormalization
	spinor e_i(Ntest, c_vector(Nagg,0));
	int index1, index2;
	for(int i = 0; i < Ntest*Nagg; i++){
		//e_i = canonical_vector(i, Ntest, Nagg);
		index1 = i / Nagg;
		index2 = i % Nagg;
		e_i[index1][index2] = 1.0;
		P_v(e_i,v_chopped[i]); //Columns of the interpolator
		e_i[index1][index2] = 0.0;
	}

	//Orthonormalization by applying Gram-Schmidt
	c_double proj; 
	for (int i = 0; i < Nagg; i++) {
		for (int nt = 0; nt < Ntest; nt++) {
			for (int j = 0; j < nt; j++) {
				proj = 0;//dot(v_chopped[nt*Nagg+i], v_chopped[j*Nagg+i]);
				for (int n = 0; n < LV::Ntot; n++) {
        			for (int alf = 0; alf < 2; alf++) {
            		proj += v_chopped[nt*Nagg+i][n][alf] * std::conj(v_chopped[j*Nagg+i][n][alf]);
        			}
    			}

				for(int n=0; n<LV::Ntot; n++){
					for(int alf=0; alf<2; alf++){
						v_chopped[nt*Nagg+i][n][alf] = v_chopped[nt*Nagg+i][n][alf] - proj * v_chopped[j*Nagg+i][n][alf];
					}
				}
				
			}
			normalize(v_chopped[nt*Nagg+i]);
		}
	}

 	//We sum all the columns of the interpolator that belong to the same aggregate and store the result
	//in a single vector. We do that for each aggregate. This enables us to have all the information of 
	//the locally orthonormalized test vectors in a single vector.

	for(int i = 0; i < Ntest; i++){
		for(int n = 0; n < LV::Ntot; n++){
			for(int alf=0; alf<2; alf++){
				interpolator_columns[i][n][alf] = 0;
			}
		}
	}

	for(int i = 0; i < Ntest; i++){
		for(int j = 0; j < Nagg; j++){
			for(int n = 0; n < LV::Ntot; n++){
				for(int alf=0; alf<2; alf++){
					interpolator_columns[i][n][alf] += v_chopped[i*Nagg + j][n][alf];
				}
			}
			
		}
	}
}; 

void AMG::setUpPhase(const double& eps,const int& Nit) {
	for (int i = 0; i < LevelV::Ntest[0]; i++) {
		for (int n = 0; n < LevelV::Nsites[0]; n++) {
			for (int dof = 0; dof < LevelV::DOF[0]; dof++) {
				test_vectors[i][n][dof] = 1;//i*LevelV::Nsites[0]*LevelV::DOF[0] + n*LevelV::DOF[0] + dof;
			}
		}
	}
	interpolator_columns = test_vectors; 
	
	//Call MPI for SAP parallelization
	/*
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	spinor rhs(LV::Ntot, c_vector(2, 0));
	using namespace AMGV; //AMG parameters namespace
	*/

	//Test vectors random initialization for each level (except the coarsest level)
	/*
	for(int l=0; l<AMGV::levels-1; l++){
		for (int i = 0; i < LevelV::Ntest[l]; i++) {
		for (int n = 0; n < LevelV::Nsites[l]; n++) {
		for (int dof = 0; dof < LevelV::DOF[l]; dof++) {
			level_l[l].test_vectors[i][n][dof] = eps * RandomU1();
		}
		}
		}
	}
	*/
	//level_l[0].makeDirac(); //Build Dirac operator for the finest level
	//Improving the test vectors by approximately solving the linear system D test_vectors[i] = rhs 
	/*
	for(int l=0; l<AMGV::levels-1; l++){
		for (int i = 0; i < LevelV::Ntest[l]; i++) {
			//Right hand side of the linear system 
			//rhs = test_vectors[i];  
			//rhs(LV::Ntot, c_vector(2, 0));
			//The result will be stored in test_vectors[i]
			//double startT, endT;
			//startT = MPI_Wtime();
			//I need to rewrite this method as a general class that I can use for every level
			SAP(GConf.Conf, rhs, level_l[l].test_vectors[i], m0, AMGV::SAP_test_vectors_iterations,SAPV::sap_blocks_per_proc);  
			//endT = MPI_Wtime();
			//SAP_time += endT - startT; 	
		}
		//Build the interpolator between level l and l+1
		level_l[l].interpolator_columns = test_vectors; 
		level_l[l].orthonormalize(); 
		level_l[l].makeCoarseLinks(); //Make coarse gauge links which define the operator D for the next level
	}
	*/

	/*
	//Improving the interpolator quality by iterating over the two-grid method defined by the current test vectors
	if (rank == 0){std::cout << "Improving interpolator" << std::endl;}
	for (int it = 0; it < Nit; it++) {
		if (rank == 0){std::cout << "****** Bootstrap iteration " << it << " ******" << std::endl;}
		for(int l = 0; l<AMGV::levels-1; l++){

		for (int i = 0; i < Ntest; i++) {
			//Number of two-grid iterations for each test vector 
			int two_grid_iter = 1; 
			//Tolerance for the two-grid method, but in this case it is not relevant. Still, it is needed to call
			//the function
			double tolerance = 1e-10; 
			rhs = test_vectors[i];
			TwoGrid(two_grid_iter,tolerance, rhs,rhs, test_vectors[i], false); 
			//I have to call a function that applies the method for the other levels.
			//For instance, for l=0 it will apply the full method.
			//For l=1 it will only apply the method for the levels l>1 
			//For l=2 ""                                        "" l>2  
		}
		//"Assemble" the new interpolator
		interpolator_columns = test_vectors; 
		orthonormalize();
		initializeCoarseLinks();
		}

	}

	if (rank == 0){std::cout << "Set-up phase finished" << std::endl;}
	*/	
}

/*
	Prolongation operator times a spinor x = P v
	x_i = P_ij v_j. dim(P) = 2 Ntot x Ntest Nagg, 
	dim(v) = Ntest Nagg, dim(x) = 2 Ntot
*/
void AMG::P_v(const spinor& v,spinor& out) {
	//Loop over columns
	for(int n = 0; n < LV::Ntot; n++){
		for(int alf = 0; alf < 2; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	using namespace AMGV; //AMG parameters namespace
	int x_coord, t_coord, s_coord; //Coordinates of the lattice point
	int k, a;
	int i, j; //Loop indices
	for (j = 0; j < Ntest * Nagg; j++) {
		k = j / Nagg; //Number of test vector
		a = j % Nagg; //Number of aggregate
		for (i = 0; i < Agg[a].size(); i++) {
			x_coord = XCoord[Agg[a][i]], t_coord = TCoord[Agg[a][i]], s_coord = SCoord[Agg[a][i]];
			out[Coords[x_coord][t_coord]][s_coord] += interpolator_columns[k][Coords[x_coord][t_coord]][s_coord] * v[k][a];			
		}
	}
}

/*
	Restriction operator times a spinor on the coarse grid, x = P^H v
	x_i = P^H_ij v_j. dim(P^H) =  Ntest Nagg x 2 Ntot, Nagg = block_x * block_t
	dim(v) = 2 Ntot, dim(x) = Ntest Nagg
*/

void AMG::Pt_v(const spinor& v,spinor& out) {
	//Restriction operator times a spinor
	using namespace AMGV;
	for(int n = 0; n < Ntest; n++){
		for(int alf = 0; alf < Nagg; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	int k, a;
	int x_coord, t_coord, s_coord; //Coordinates of the lattice point
	int i, j; //Loop indices
	for (i = 0; i < Ntest*Nagg; i++) {
		k = i / Nagg; //number of test vector
		a = i % Nagg; //number of aggregate
		for (j = 0; j < Agg[a].size(); j++) {
			x_coord = XCoord[Agg[a][j]], t_coord = TCoord[Agg[a][j]], s_coord = SCoord[Agg[a][j]];
			out[k][a] += std::conj(interpolator_columns[k][Coords[x_coord][t_coord]][s_coord]) * v[Coords[x_coord][t_coord]][s_coord];
		}
	}
}

/*
	Intialize the coarse gauge links for Dc
*/
void AMG::initializeCoarseLinks(){
	c_double P[2][2][2], M[2][2][2]; 
	
	P[0][0][0] = 1.0; P[0][0][1] = 1.0;
	P[0][1][0] = 1.0; P[0][1][1] = 1.0; 

	P[1][0][0] = 1.0; P[1][0][1] = -I_number;
	P[1][1][0] = I_number; P[1][1][1] = 1.0; 

	M[0][0][0] = 1.0; M[0][0][1] = -1.0;
	M[0][1][0] = -1.0; M[0][1][1] = 1.0; 

	M[1][0][0] = 1.0; M[1][0][1] = I_number;
	M[1][1][0] = -I_number; M[1][1][1] = 1.0; 

	c_matrix &U = GConf.Conf;
	std::vector<spinor> &w = interpolator_columns;
	c_double Lm, Lp, R;
	int block_r, block_l;
	for(int x=0; x<LV::Nblocks; x++){
	for(int alf=0; alf<2;alf++){
	for(int bet=0; bet<2;bet++){
	for(int p = 0; p<AMGV::Ntest; p++){
	for(int s = 0; s<AMGV::Ntest; s++){

		A_coeff[x][alf][bet][p][s] = 0;
		B_coeff[x][alf][bet][p][s][0] = 0; B_coeff[x][alf][bet][p][s][1] = 0;
		C_coeff[x][alf][bet][p][s][0] = 0; C_coeff[x][alf][bet][p][s][1] = 0;
		for(int n : LatticeBlocks[x]){
		for(int mu : {0,1}){
			getLatticeBlock(RightPB[n][mu], block_r); //Get the block index for the right periodic boundary
			getLatticeBlock(LeftPB[n][mu], block_l); //Get the block index for the right periodic boundary

			Lm = 0.5 * M[mu][alf][bet] * std::conj(w[p][n][alf]) * U[n][mu];
			Lp = 0.5 * P[mu][alf][bet] * std::conj(w[p][n][alf]) * std::conj(U[LeftPB[n][mu]][mu]);	
			//std::cout << "Lm : " << Lm <<  std::endl;
			//if n+\hat{mu} in Block(x)
			if (block_r == x)
				//[A(x)]^{alf,bet}_{p,s} --> A_coeff[x][alf][bet][p][s] 
				A_coeff[x][alf][bet][p][s] += Lm * w[s][RightPB[n][mu]][bet] * SignR[n][mu];
			//if n+\hat{mu} in Block(x+hat{mu})
			else if (block_r == RightPB_blocks[x][mu])
				//[B_mu(x)]^{alf,bet}_{p,s} --> B_coeff[x][alf][bet][p][s][mu]
				B_coeff[x][alf][bet][p][s][mu] += Lm * w[s][RightPB[n][mu]][bet] * SignR[n][mu];
			//if n-\hat{mu} in Block(x)	
			if (block_l == x)
				A_coeff[x][alf][bet][p][s] += Lp * w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
			//if n-\hat{mu} in Block(x-hat{mu})
			else if (block_l == LeftPB_blocks[x][mu])
				//[C_mu(x)]^{alf,bet}_{p,s} --> C_coeff[x][alf][bet][p][s][mu]
				C_coeff[x][alf][bet][p][s][mu] += Lp * w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
		}//mu 
		}//n 


	//---------Close loops---------//
	} //s
	} //p
	} //bet
	} //alf
	} //x 

}

/*
	Coarse grid matrix operator Dc = P^H D P times a spinor v
	dim(Dc) = Ntest Nagg x Ntest Nagg, dim(v) = Ntest Nagg,
	Dc_{xy} = A(x) d_{xy} + \sum_{mu} B_mu(x) d_{x+mu,y} + \sum_{mu} C_mu(x)  d_{x-mu,y} 
	With x, y lattice blocks.
*/
void AMG::Pt_D_P(const spinor& v,spinor& out){
	for(int x = 0; x<LV::Nblocks; x++){
	for(int alf = 0; alf<2; alf++){
	for(int p = 0; p<AMGV::Ntest; p++){
		out[p][2*x+alf] = (m0+2)*v[p][2*x+alf]; //Mass term
	for(int bet = 0; bet<2; bet++){
	for(int s = 0; s<AMGV::Ntest; s++){
		//out[test_vec][aggregate]
		out[p][2*x+alf] -= A_coeff[x][alf][bet][p][s] * v[s][2*x+bet];

		for(int mu:{0,1}){
			out[p][2*x+alf] -=  (B_coeff[x][alf][bet][p][s][mu] * v[s][2*RightPB_blocks[x][mu]+bet]
							    +C_coeff[x][alf][bet][p][s][mu] * v[s][2*LeftPB_blocks[x][mu]+bet]);
		}
		

	}
	}
	}
	}
	}
}


/*
	Two-grid method for solving the linear system D x = phi
*/
int AMG::TwoGrid(const int& max_iter, const double& tol, const spinor& x0, 
	const spinor& phi, spinor& x, const bool& print_message) {

	x = x0; //Solution spinor
	spinor r(LV::Ntot,c_vector(2,0)); //Residual 
	double err; //Error = sqrt(dot(r,r))
	int k = 0; //Iteration number
	double norm = sqrt(std::real(dot(phi,phi)));

	//The convergence criterion is ||r|| < ||phi|| * tol

	while(k < max_iter){
		//Pre-smoothing
		if (nu1>0){
			//gmres smoothing
			SAP(GConf.Conf, phi, x, m0, nu1,SAPV::sap_blocks_per_proc); 
		} 

		//*************Coarse grid correction*************//
		double startT, endT;
		startT = MPI_Wtime();
		//x = x + P*Dc^-1 * P^H * (phi-D*x)  
		spinor temp(LV::Ntot,c_vector(2,0));
		spinor Dphi(LV::Ntot,c_vector(2,0));
		D_phi(GConf.Conf, x, Dphi,m0); //D x
		for(int n = 0; n<LV::Ntot; n++){
			for(int alf=0; alf<2; alf++){
				temp[n][alf] = phi[n][alf] - Dphi[n][alf]; //temp = phi - D x
			}
		}
		//spinor Pt_r = Pt_v(temp); //Pt_r = P^H (phi - D x)
		spinor Pt_r(AMGV::Ntest, c_vector(AMGV::Nagg, 0));
		Pt_v(temp,Pt_r);
		
		//Using GMRES for the coarse grid solver 
		spinor gmresResult(AMGV::Ntest, c_vector(AMGV::Nagg, 0));
		gmres_c_level.fgmres(Pt_r,Pt_r,gmresResult,false);
		P_v(gmresResult,temp);

		for(int n = 0; n<LV::Ntot; n++){
			for(int alf=0; alf<2; alf++){
				x[n][alf] += temp[n][alf]; 
			}
		}
	
		endT = MPI_Wtime();
		coarse_time += endT - startT; //Measuring time spent for solving the coarse level 
		//************************************************//
		
		//Post-smoothing
		if (nu2>0){
			//Measure time spent smoothing
			double startT, endT;
			startT = MPI_Wtime();
			SAP(GConf.Conf, phi, x, m0, nu2, SAPV::sap_blocks_per_proc); 
			endT = MPI_Wtime();
			smooth_time += endT - startT; //Add post-smoothing time
			SAP_time += endT - startT; //Add post-smoothing time
		}
		D_phi(GConf.Conf, x,Dphi, m0); 

		//r = phi - D_phi(GConf.Conf, x, m0);		
		for(int n = 0; n < LV::Ntot; n++){
			for(int alf=0; alf<2; alf++){
				r[n][alf] = phi[n][alf] - Dphi[n][alf];
			}
		}
		
		
		err = sqrt(std::real(dot(r,r)));

		if (err < tol*norm){
			if (print_message == true){
				std::cout << "Two-grid method converged in " << k+1 << " iterations" << " Error " << err << std::endl;
			}
		return 1;
		}
		k++;
	}
	if (print_message == true){
		std::cout << "Two-grid did not converge in " << max_iter << " iterations" << " Error " << err << std::endl;
	}
	return 0;
}
