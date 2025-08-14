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
	std::cout << "Pv called for level " << level << std::endl;
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
	for (int i = 0; i < Ntest; i++) {
		for (int n = 0; n < Nsites; n++) {
			for (int dof = 0; dof < DOF; dof++) {
				interpolator_columns[i][n][dof] = 1;//i*Nsites*DOF + n*DOF + dof;
			}
		}
	}
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
	std::cout << "G1 size " << G1.size() << std::endl;
	std::cout << "G2 size " << G2.size() << std::endl;
	std::cout << "G3 size " << G3.size() << std::endl;
	std::cout << "Next level " << next_level.level << std::endl;
	std::cout << "A size " << A_coeff.size() << std::endl;
	std::cout << "B size " << B_coeff.size() << std::endl;
	std::cout << "C size " << C_coeff.size() << std::endl;
	int indxA; int indxBC[2]; //Indices for A, B and C coefficients
	int block_r;
	int block_l;
	/*
	//Printing the test vectors
	for(int ntest = 0; ntest<Ntest; ntest++){
		for(int dof = 0; dof<DOF; dof++){
			std::cout << "w[" << ntest << "][" << dof << "] = ";
			for(int n = 0; n < Nsites; n++){
				std::cout << w[ntest][n][dof] << " ";
			}
			std::cout << std::endl;
		}
	}
	*/
	//p and s are the coarse colors
	//c and b are the colors at the current level
	std::cout << "NBlocks " << NBlocks << " Ntest " << Ntest << " Nagg " << Nagg << std::endl;
	std::cout << "colors " << colors << " DOF " << DOF << std::endl;
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
			A_coeff[indxA] -= std::conj(w[p][n][2*c+alf]) * G1[getG1index(n,alf,bet,c,b)] * std::conj(w[s][n][2*b+bet]);
			for(int mu : {0,1}){
				getLatticeBlock(RightPB_l[level][n][mu], block_r); //block_r: block where RightPB_l[n][mu] lives
				getLatticeBlock(LeftPB_l[level][n][mu], block_l); //block_l: block where LeftPB_l[n][mu] lives
				wG2 = std::conj(w[p][n][2*c+alf]) * G2[getG2G3index(n,alf,bet,c,b,mu)]; 
				wG3 = std::conj(w[p][n][2*c+alf]) * G3[getG2G3index(n,alf,bet,c,b,mu)];
				
				//Only diff from zero when n+hat{mu} in Block(x)
				if (block_r == x)
					A_coeff[indxA] += wG2 * w[s][RightPB_l[level][n][mu]][2*b+bet] * SignR_l[level][n][mu];
				//Only diff from zero when n+hat{mu} in Block(x+hat{mu})
				else if (block_r == RightPB_l[level+1][x][mu])
					B_coeff[indxBC[mu]] += wG2 * w[s][RightPB_l[level][n][mu]][2*b+bet];// * SignR_l[level][n][mu]; //Sign considered in the operator
				//Only diff from zero when n-hat{mu} in Block(x)
				if (block_l == x)
					A_coeff[indxA] += wG3 * w[s][LeftPB_l[level][n][mu]][2*b+bet] *  SignL_l[level][n][mu];
				//Only diff from zero when n-hat{mu} in Block(x-hat{mu})
				else if (block_l == LeftPB_l[level+1][x][mu])
					C_coeff[indxBC[mu]] += wG3 * w[s][LeftPB_l[level][n][mu]][2*b+bet];// * SignL_l[level][n][mu];
	
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
