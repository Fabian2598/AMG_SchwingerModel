#include "level.h"

void Level::makeAggregates(){
    if (level == 0){

	for (int x = 0; x < LevelV::BlocksX[level]; x++) {
	for (int t = 0; t < LevelV::BlocksT[level]; t++) {
	for (int s = 0; s < 2; s++) {
		int x0 = x * x_elements, t0 = t * t_elements;
		int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
		int aggregate = x * LevelV::BlocksT[level] * 2 + t * 2 + s;
		int count = 0;
		//x and t are redefined in the following loop
		for (int x = x0; x < x1; x++) {
		for (int t = t0; t < t1; t++) {
			int i = x * LevelV::NtSites[level] * 2  + t * 2 + s;
			XCoord[i] = x; TCoord[i] = t; SCoord[i] = s;
			Agg[aggregate * sites_per_block + count] = i;
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

    }
    else{

    for (int x = 0; x < LevelV::BlocksX[level]; x++) {
	for (int t = 0; t < LevelV::BlocksT[level]; t++) {
        int x0 = x * x_elements, t0 = t * t_elements;
		int x1 = (x + 1) * x_elements, t1 = (t + 1) * t_elements;
	for (int s = 0; s < 2; s++) {
    for(int ntest = 0; ntest< LevelV::Ntest[level]; ntest++){
        //FINISH THIS PART ... WRITTEN ON 05/08/2025
		int aggregate = x * LevelV::BlocksT[level] * 2 + t * 2 + s;
		int count = 0;
		//x and t are redefined in the following loop
		for (int x = x0; x < x1; x++) {
		for (int t = t0; t < t1; t++) {
			int i = x * LevelV::NtSites[level] * 2  + t * 2 + s;
			XCoord[i] = x; TCoord[i] = t; SCoord[i] = s;
			Agg[aggregate * sites_per_block + count] = i;
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
	}

    }

}

void Level::P_v(const spinor& v,spinor& out){
	//Loop over columns
	for(int n = 0; n < LevelV::Nsites[level]; n++){
		for(int alf = 0; alf < LevelV::DOF[level]; alf++){
			out[n][alf] = 0.0; //Initialize the output spinor
		}
	}
	int x_coord, t_coord, s_coord; //Coordinates of the lattice point
	int k, a;
	int i, j; //Loop indices


	for (j = 0; j < LevelV::Ntest[level] * LevelV::Nagg[level]; j++) {
		k = j / Nagg; //Number of test vector
		a = j % Nagg; //Number of aggregate
		for (i = 0; i < Agg[a].size(); i++) {
			x_coord = XCoord[Agg[a][i]], t_coord = TCoord[Agg[a][i]], s_coord = SCoord[Agg[a][i]];
			out[Coords[x_coord][t_coord]][s_coord] += interpolator_columns[k][Coords[x_coord][t_coord]][s_coord] * v[k][a];			
		}
	}
    
}
