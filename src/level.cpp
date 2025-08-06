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
			//XCoord[i] = x; TCoord[i] = t; SCoord[i] = s;
			LatticeBlocks[block * sites_per_block + count] = i;
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
		std::cout << LatticeBlocks[i * sites_per_block + n] << " ";
		}	
		std::cout << std::endl;
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
				interpolator_columns[i][n][dof] = i*Nsites*DOF + n*DOF + dof;
			}
		}
	}
}