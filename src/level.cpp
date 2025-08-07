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
		//std::cout << Agg[i * sites_per_block * LevelV::Colors[level] + n * LevelV::Colors[level] + c ] << " ";
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
			int i = x * LevelV::NtSites[level] + t;
			//LatticeBlocks[block * sites_per_block + count] = i;
			LatticeBlocks[block][count] = i;
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


void Level::D_operator(const spinor& v, spinor& out){
	//For the finest level the D_operator is defined as the Dirac operator
	if (level == 0){D_phi(U, v, out,mass::m0);}
	//I could also define "coarse gauge links" for the finest level and then use the same implementation as for level>0,
	//that way I can get rid of the if/else ...
	// In such a case B_coeff[x][alf][bet][c][b][mu] = 0.5 * (1-sigma_mu)^{alf,bet} U_mu(x). c and b would be spurious 
	//for the finest level, since we don't have colors. The same idea applies to C_coeff and A_coeff = (m0+2).
	//Maybe I try it later for improving the code in the future, for the moment I will just leave the if/else. 

	//For level > 0 we use the coarse gauge links generated in level l-1
	else{
		for(int x = 0; x<Nsites;x++){
		for(int alf = 0; alf<2; alf++){
		for(int c = 0; c<colors; c++){
			out[x][2*c+alf] = (mass::m0+2)*v[x][2*c+alf];
			//out[x][dof]

			//out[c][2*x+alf] = (mass::m0+2)*v[c][2*x+alf]; //Mass term
		for(int bet = 0; bet<2; bet++){
		for(int b = 0; b<colors; b++){
			out[x][2*c+alf] -= G1[getG1index(x,alf,bet,c,b)] * v[x][2*b+bet];
			
			//out[x][dof]
			
			//out[test_vec][aggregate]
			//out[c][2*x+alf] -= A_coeff[x][alf][bet][c][b] * v[b][2*x+bet];

			for(int mu:{0,1}){
				out[x][2*c+alf] -= ( G2[getG2G3index(x,alf,bet,c,b,mu)] * v[RightPB_blocks[x][mu]][2*b+bet]
							    + G3[getG2G3index(x,alf,bet,c,b,mu)] * v[LeftPB_blocks[x][mu]][2*b+bet] );			
				//out[x][dof]
				
				//out[c][2*x+alf] -=  (B_coeff[x][alf][bet][c][b][mu] * v[b][2*RightPB_blocks[x][mu]+bet]
				//			    	+C_coeff[x][alf][bet][c][b][mu] * v[b][2*LeftPB_blocks[x][mu]+bet]);
			}
		

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
	for (int i = 0; i < Ntest; i++) {
		for (int n = 0; n < Nsites; n++) {
			for (int dof = 0; dof < DOF; dof++) {
				interpolator_columns[i][n][dof] = 1;//i*Nsites*DOF + n*DOF + dof;
			}
		}
	}
}


void Level::makeCoarseLinks(Level& next_level){
	//Generate links for the next level
	c_double P[2][2][2], M[2][2][2]; 
	
	P[0][0][0] = 1.0; P[0][0][1] = 1.0;
	P[0][1][0] = 1.0; P[0][1][1] = 1.0; 

	P[1][0][0] = 1.0; P[1][0][1] = -I_number;
	P[1][1][0] = I_number; P[1][1][1] = 1.0; 

	M[0][0][0] = 1.0; M[0][0][1] = -1.0;
	M[0][1][0] = -1.0; M[0][1][1] = 1.0; 

	M[1][0][0] = 1.0; M[1][0][1] = I_number;
	M[1][1][0] = -I_number; M[1][1][1] = 1.0; 

	//Make gauge links for level
	std::vector<spinor> &w = interpolator_columns;
	c_double Lm, Lp, R;
	c_vector &A_coeff = next_level.G1; 
	c_vector &B_coeff = next_level.G2;
	c_vector &C_coeff = next_level.G3;
	std::cout << "A size " << A_coeff.size() << std::endl;
	std::cout << "B size " << B_coeff.size() << std::endl;
	std::cout << "C size " << C_coeff.size() << std::endl;

	
	//Here instead of using the gauge links I have to write everything in terms of the G's
	
	for(int x=0; x<LV::Nblocks; x++){
	for(int alf=0; alf<2;alf++){
	for(int bet=0; bet<2;bet++){
	for(int p = 0; p<AMGV::Ntest; p++){
	for(int s = 0; s<AMGV::Ntest; s++){
		A_coeff[getAindex(x,alf,bet,p,s)] = 0;
		B_coeff[getBCindex(x,alf,bet,p,s,0)] = 0; B_coeff[getBCindex(x,alf,bet,p,s,1)] = 0;
		C_coeff[getBCindex(x,alf,bet,p,s,0)] = 0; C_coeff[getBCindex(x,alf,bet,p,s,1)] = 0;
		for(int n : LatticeBlocks[x]){
		for(int mu : {0,1}){
			Lm = 0.5 * M[mu][alf][bet] * std::conj(w[p][n][alf]) * U[n][mu];
			Lp = 0.5 * P[mu][alf][bet] * std::conj(w[p][n][alf]) * std::conj(U[LeftPB[n][mu]][mu]);
			//           [A(x)]^{alf,bet}_{p,s} --> A_coeff[x][alf][bet][p][s] 
			//--------------- 1 - sigma_mu---------------//
			R = 0.0;
			//if n+\hat{mu} in Block(x)
			if (std::find(LatticeBlocks[x].begin(), LatticeBlocks[x].end(), RightPB[n][mu]) != LatticeBlocks[x].end()){
				R = w[s][RightPB[n][mu]][bet] * SignR[n][mu];
			}
			A_coeff[getAindex(x,alf,bet,p,s)] += Lm * R;
			//-------------- 1 + sigma_mu --------------//
			R = 0.0;
			//if n-\hat{mu} in Block(x)
			if (std::find(LatticeBlocks[x].begin(), LatticeBlocks[x].end(), LeftPB[n][mu]) != LatticeBlocks[x].end()){
				R = w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
			}
			A_coeff[getAindex(x,alf,bet,p,s)] += Lp * R;

			//			[B_mu(x)]^{alf,bet}_{p,s} --> B_coeff[x][alf][bet][p][s][mu]
			R = 0.0;
			//if n+\hat{mu} in Block(x+hat{mu})
			if (std::find(LatticeBlocks[RightPB_blocks[x][mu]].begin(), LatticeBlocks[RightPB_blocks[x][mu]].end(), RightPB[n][mu]) != LatticeBlocks[RightPB_blocks[x][mu]].end()){
				R = w[s][RightPB[n][mu]][bet] * SignR[n][mu];
			}
			B_coeff[getBCindex(x,alf,bet,p,s,mu)] += Lm * R;
			
			//			[C_mu(x)]^{alf,bet}_{p,s} --> C_coeff[x][alf][bet][p][s][mu]
			R = 0.0;
			//if n-\hat{mu} in Block(x-hat{mu})
			if (std::find(LatticeBlocks[LeftPB_blocks[x][mu]].begin(), LatticeBlocks[LeftPB_blocks[x][mu]].end(), LeftPB[n][mu]) != LatticeBlocks[LeftPB_blocks[x][mu]].end()){
				R = w[s][LeftPB[n][mu]][bet] * SignL[n][mu];
			}
			C_coeff[getBCindex(x,alf,bet,p,s,mu)] += Lp * R;

		}//mu 
		}//n 


	//---------Close loops---------//
	} //s
	} //p
	} //bet
	} //alf
	} //x 
	
}

