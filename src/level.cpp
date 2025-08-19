#include "level.h"

void Level::makeAggregates(){
    for (int x = 0; x < LevelV::BlocksX[level]; x++) {
	for (int t = 0; t < LevelV::BlocksT[level]; t++) {
        int x0 = x * x_elements, t0 = t * t_elements;
		int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
		
	for (int s = 0; s < 2; s++) {    
		int aggregate = x * LevelV::BlocksT[level] * 2 + t * 2 + s;
		int count = 0;
		//x and t are redefined in the following loop
		for (int x = x0; x < x1; x++) {
		for (int t = t0; t < t1; t++) {
		for(int c = 0; c< LevelV::Colors[level]; c++){
			//Vectorization of coordinates, test vector dof and spin
			int i = x * LevelV::NtSites[level] * LevelV::Colors[level] * 2 + t * LevelV::Colors[level] * 2 + c * 2 + s;
			sCoords[i] = s; cCoords[i] = c;  //spin and color
			nCoords[i] = x * LevelV::NtSites[level] + t; //Vectorization of coordinates
			Agg[aggregate * sites_per_block * LevelV::Colors[level] + count] = i;
			count++;
			
        }
		}
		}
		if (count != LevelV::Colors[level] * x_elements * t_elements) {
			std::cout << "Aggregate " << aggregate << " has " << count << " elements" << std::endl;
		}
  
	}	
	}
	}

}

void Level::printAggregates() {
	std::cout << "Aggregates for level " << level << ": " << std::endl;
	int s;
	int var;
	for (int i = 0; i < Nagg; i++) {
		s =  i % 2; //spin
		std::cout << "---Aggregate " << i << "---\n";
		for (int n = 0; n < sites_per_block; n++) {
		for (int c = 0; c < colors; c++){
			var = Agg[i * sites_per_block * LevelV::Colors[level] + n * LevelV::Colors[level] + c ];
			std::cout << "n=" << nCoords[var] << " c=" << cCoords[var]  << " s="  << sCoords[var] << " " << Agg[i * sites_per_block * LevelV::Colors[level] + n * LevelV::Colors[level] + c ] << " \n";
		}	
		}
		std::cout << std::endl;
	}
}

void Level::makeBlocks(){
    for (int x = 0; x < LevelV::BlocksX[level]; x++) {
	for (int t = 0; t < LevelV::BlocksT[level]; t++) {
        int x0 = x * x_elements, t0 = t * t_elements;
		int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
		int block = x * LevelV::BlocksT[level] + t;
		int count = 0;
		//x and t are redefined in the following loop
		for (int x = x0; x < x1; x++) {
		for (int t = t0; t < t1; t++) {
			//LatticeBlocks[block * sites_per_block + count] = i;
			LatticeBlocks[block][count] = x * LevelV::NtSites[level] + t;
			count++;
        }
		}
		if (count != x_elements * t_elements) {
			std::cout << "Block " << block << " has " << count << " elements" << std::endl;
		}
	}	
	}
}

void Level::printBlocks() {
	std::cout << "Blocks for level " << level << ": " << std::endl;
	for (int i = 0; i < LevelV::NBlocks[level]; i++) {
		std::cout << "Block " << i << ":";
		for (int n = 0; n < sites_per_block; n++) {	
		//std::cout << LatticeBlocks[i * sites_per_block + n] << " ";
		std::cout << LatticeBlocks[i][n] << " ";
		}	
		std::cout << std::endl;
	}
}

void Level::makeDirac(){
	c_double P[2][2][2], M[2][2][2]; 
	//P = 1 + sigma
	P[0][0][0] = 1.0; P[0][0][1] = 1.0;
	P[0][1][0] = 1.0; P[0][1][1] = 1.0; 

	P[1][0][0] = 1.0; P[1][0][1] = -I_number;
	P[1][1][0] = I_number; P[1][1][1] = 1.0; 

	//M = 1- sigma
	M[0][0][0] = 1.0; M[0][0][1] = -1.0;
	M[0][1][0] = -1.0; M[0][1][1] = 1.0; 

	M[1][0][0] = 1.0; M[1][0][1] = I_number;
	M[1][1][0] = -I_number; M[1][1][1] = 1.0; 


	for(int x=0; x<Nsites; x++){
	for(int alf=0; alf<2;alf++){
	for(int bet=0; bet<2;bet++){
	for(int c = 0; c<colors; c++){
	for(int b = 0; b<colors; b++){
		G1[getG1index(x,alf,bet,c,b)] = 0;//This coefficient is not used at level 0
		G2[getG2G3index(x,alf,bet,c,b,0)] = 0; G2[getG2G3index(x,alf,bet,c,b,1)] = 0;
		G3[getG2G3index(x,alf,bet,c,b,0)] = 0; G3[getG2G3index(x,alf,bet,c,b,1)] = 0;
		//For level = 0 
		for(int mu : {0,1}){
			G2[getG2G3index(x,alf,bet,c,b,mu)] = 0.5 * M[mu][alf][bet] * U[x][mu];
			G3[getG2G3index(x,alf,bet,c,b,mu)] = 0.5 * P[mu][alf][bet] * std::conj(U[LeftPB_l[level][x][mu]][mu]);
		}
		
	
	}
	}
	}
	}
	}
		
}

void Level::orthonormalize(){
	/*
	Local orthonormalization of the test vectors
	
	Each test vector is chopped into the Nagg aggregates, which yields Ntest*Nagg columns for the interpolator.
	Each column is orthonormalized with respect to the others that belong to the same aggregate.
	This follows the steps from Section 3.1 of A. Frommer et al "Adaptive Aggregation-Based Domain Decomposition 
	Multigrid for the Lattice Wilson-Dirac Operator", SIAM, 36 (2014).
	*/

	//Getting the columns of the interpolator for the orthonormalization
	
	spinor e_i(NBlocks, c_vector(2*Ntest,0));
	int x, s, c, a; //block, spin and color
	for(int i = 0; i < Ntest*Nagg; i++){
		//e_i = canonical_vector(i, Ntest, Nagg);
		c = i / Nagg; //test vector
		a = i % Nagg; //Aggregate
		x = a / 2; //block
		s = a % 2; //spin
		e_i[x][2*c+s] = 1.0;
		P_v(e_i,v_chopped[i]); //Columns of the interpolator
		e_i[x][2*c+s] = 0.0;
	}
	
	//Orthonormalization by applying Gram-Schmidt
	c_double proj; 
	for (int i = 0; i < Nagg; i++) {
		s = i % Nagg; x = s / 2; //spin, lattice block
		for (int nt = 0; nt < Ntest; nt++) {
			for (int j = 0; j < nt; j++) {
				proj = 0;//dot(v_chopped[nt*Nagg+i], v_chopped[j*Nagg+i]);			
				for (int n = 0; n < Nsites; n++) {
        			for (int alf = 0; alf < DOF; alf++) {
            			proj += v_chopped[nt*Nagg+i][n][alf] * std::conj(v_chopped[j*Nagg+i][n][alf]);
        			}
    			}

				for(int n=0; n<Nsites; n++){
					for(int alf=0; alf<DOF; alf++){
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

	for(int nt = 0; nt < Ntest; nt++){
		for(int n = 0; n < Nsites; n++){
			for(int alf=0; alf<DOF; alf++){
				interpolator_columns[nt][n][alf] = 0;
			}
		}
	}

	for(int nt = 0; nt < Ntest; nt++){
		for(int a = 0; a < Nagg; a++){
			for(int n = 0; n < Nsites; n++){
				for(int alf=0; alf<DOF; alf++){
					interpolator_columns[nt][n][alf] += v_chopped[nt*Nagg + a][n][alf];
				}
			}
		}
	}



}; 

void Level::checkOrthogonality() {
	//Check orthogonality of the test vectors
	//aggregate 
	for(int block = 0; block < NBlocks; block++){
	for (int alf = 0; alf < 2; alf++) {
		//checking orthogonality 
		for (int i = 0; i < Ntest; i++) {
		for (int j = 0; j < Ntest; j++) {
			c_double dot_product = 0.0;
			for (int n: LatticeBlocks[block]) {
			for (int c = 0; c<colors; c++){
				dot_product += std::conj(interpolator_columns[i][n][2*c+alf]) * interpolator_columns[j][n][2*c+alf];
			}
			}
			if (std::abs(dot_product) > 1e-8 && i!=j) {
				std::cout << "Block " << block << " spin " << alf << std::endl;
				std::cout << "Level " << level << std::endl;
				std::cout << "Test vectors " << i << " and " << j << " are not orthogonal: " << dot_product << std::endl;
				exit(1);
			}
			else if(std::abs(dot_product-1.0) > 1e-8 && i==j){
				std::cout << "Level " << level << std::endl;
				std::cout << "Test vector " << i << " not orthonormalized " << dot_product << std::endl;
				exit(1);
			}

		}
		}
	}
	}
	std::cout << "Test vectors on level " << level << " are orthonormalized " << std::endl;
}

//Dirac operator at the current level
void Level::D_operator(const spinor& v, spinor& out){	

	for(int x = 0; x<Nsites;x++){
	for(int alf = 0; alf<2; alf++){
	for(int c = 0; c<colors; c++){
		out[x][2*c+alf] = (mass::m0+2)*v[x][2*c+alf];
	for(int bet = 0; bet<2; bet++){
	for(int b = 0; b<colors; b++){
		//For this term the antiperiodic sign is considered in G1
		out[x][2*c+alf] -= G1[getG1index(x,alf,bet,c,b)] * v[x][2*b+bet];

		//For the other terms we write the antiperiodic sign here, not in the coefficients
		for(int mu:{0,1}){
			out[x][2*c+alf] -= ( G2[getG2G3index(x,alf,bet,c,b,mu)] * SignR_l[level][x][mu] * v[RightPB_l[level][x][mu]][2*b+bet]
							+ G3[getG2G3index(x,alf,bet,c,b,mu)] * SignL_l[level][x][mu] * v[LeftPB_l[level][x][mu]][2*b+bet] );
		}
	}
	}
	}
	}
	}

}

void Level::P_v(const spinor& v,spinor& out){
	//Loop over columns
	for(int n = 0; n < Nsites; n++){
		for(int alf = 0; alf < DOF; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	int n, s, c; //Coordinates of the lattice point
	int nc, sc,cc; //ncoarse (block), sc (spin coarse), cc (coarse color)
	int i, j; //Loop indices
	int a;

	for (j = 0; j < Ntest*Nagg; j++) {
		cc = j / Nagg; //Number of test vector
		a = j % Nagg; //Number of aggregate
		nc = a/2; //Number of lattice block
		//Each aggregate has x_elements * t_elements * Colors elements
		for (i = 0; i < colors * x_elements * t_elements; i++) {
			//Agg[a * sites_per_block * Colors + j * Colors + k]
			//Agg[a * sites_per_block * Colors + i], i from LevelV::Colors[level] * x_elements * t_elements - 1
			n = nCoords[Agg[a * sites_per_block * colors + i]];
			s = sCoords[Agg[a * sites_per_block * colors + i]];
			c = cCoords[Agg[a * sites_per_block * colors + i]];
			out[n][2*c + s] += interpolator_columns[cc][n][2*c+s] * v[nc][2*cc+s];//v[k][a];		
		}
	}
    
}


void Level::Pt_v(const spinor& v,spinor& out) {
	//Restriction operator times a spinor
	for(int n = 0; n < NBlocks; n++){
		for(int alf = 0; alf < 2*Ntest; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	int n, s, c; //Fine variables
	int nc, sc,cc; //Coarse variables s = sc
	int i, j; //Loop indices
	int a; //Aggregate
	int var; 

	for (i = 0; i < Ntest*Nagg; i++) {	
		cc = i / Nagg; //Number of test vector
		a = i % Nagg; //Number of aggregate
		nc = a/2; //Number of lattice block
		for (j = 0; j < colors * x_elements * t_elements; j++) {
			//Agg[a * sites_per_block * Colors + j], j from 0 to LevelV::Colors[level] * x_elements * t_elements - 1
			var = Agg[a * sites_per_block * colors + j]; 
			n = nCoords[var]; s = sCoords[var]; c = cCoords[var];
			out[nc][2*cc+s] += std::conj(interpolator_columns[cc][n][2*c+s]) * v[n][2*c+s];
		}
	}

}


void Level::setUp(){
	srand(19);
	for (int i = 0; i < Ntest; i++) {
		for (int n = 0; n < Nsites; n++) {
			for (int dof = 0; dof < DOF; dof++) {
				interpolator_columns[i][n][dof] = RandomU1();
			}
		}
	}

	orthonormalize();
}


/*
	Make coarse gauge links. They will be used in the next level as G1, G2 and G3.
*/
void Level::makeCoarseLinks(Level& next_level){
	//Make gauge links for level l
	std::vector<spinor> &w = interpolator_columns;
	c_double wG2, wG3;
	c_vector &A_coeff = next_level.G1; 
	c_vector &B_coeff = next_level.G2;
	c_vector &C_coeff = next_level.G3;
	int indxA; int indxBC[2]; //Indices for A, B and C coefficients
	int block_r;
	int block_l;
	//p and s are the coarse colors
	//c and b are the colors at the current level
	for(int x=0; x<NBlocks; x++){
	for(int alf=0; alf<2;alf++){
	for(int bet=0; bet<2;bet++){
	for(int p = 0; p<Ntest; p++){
	for(int s = 0; s<Ntest; s++){
		indxA = getAindex(x,alf,bet,p,s); //Indices for the next level
		indxBC[0] = getBCindex(x,alf,bet,p,s,0);
		indxBC[1] = getBCindex(x,alf,bet,p,s,1);
		A_coeff[indxA] = 0;
		B_coeff[indxBC[0]] = 0; B_coeff[indxBC[1]] = 0;
		C_coeff[indxBC[0]] = 0; C_coeff[indxBC[1]] = 0;
		for(int n : LatticeBlocks[x]){
			for(int c = 0; c<colors; c++){
			for(int b = 0; b<colors; b++){
			
				//[w*_p^(block,alf)]_{c,alf}(x) [A(x)]^{alf,bet}_{c,b} [w_s^{block,bet}]_{b,bet}(x)
			A_coeff[indxA] += std::conj(w[p][n][2*c+alf]) * G1[getG1index(n,alf,bet,c,b)] * w[s][n][2*b+bet];
			for(int mu : {0,1}){
				getLatticeBlock(RightPB_l[level][n][mu], block_r); //block_r: block where RightPB_l[n][mu] lives
				getLatticeBlock(LeftPB_l[level][n][mu], block_l); //block_l: block where LeftPB_l[n][mu] lives
				wG2 = std::conj(w[p][n][2*c+alf]) * G2[getG2G3index(n,alf,bet,c,b,mu)]; 
				wG3 = std::conj(w[p][n][2*c+alf]) * G3[getG2G3index(n,alf,bet,c,b,mu)];
				
				//Only diff from zero when n+hat{mu} in Block(x)
				if (block_r == x)
					A_coeff[indxA] += wG2 * w[s][RightPB_l[level][n][mu]][2*b+bet];// * SignR_l[level][n][mu];
				//Only diff from zero when n+hat{mu} in Block(x+hat{mu})
				else if (block_r == RightPB_l[level+1][x][mu])
					B_coeff[indxBC[mu]] += wG2 * w[s][RightPB_l[level][n][mu]][2*b+bet]; //Sign considered in the operator
				//Only diff from zero when n-hat{mu} in Block(x)
				if (block_l == x)
					A_coeff[indxA] += wG3 * w[s][LeftPB_l[level][n][mu]][2*b+bet];// *  SignL_l[level][n][mu];
				//Only diff from zero when n-hat{mu} in Block(x-hat{mu})
				else if (block_l == LeftPB_l[level+1][x][mu])
					C_coeff[indxBC[mu]] += wG3 * w[s][LeftPB_l[level][n][mu]][2*b+bet];
	
			}
			}	
			}
		}
	//---------Close loops---------//
	} //s
	} //p
	} //bet
	} //alf
	} //x 
	
}
