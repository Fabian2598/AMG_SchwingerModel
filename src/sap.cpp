#include "sap.h"

void SAP_C::SchwarzBlocks(){
    int count, block;
    int x0, t0, x1, t1;
    for (int x = 0; x < Block_x; x++) {
        for (int t = 0; t < Block_t; t++) {
            x0 = x * x_elements; t0 = t * t_elements;
            x1 = (x + 1) * x_elements; t1 = (t + 1) * t_elements;
            block = x * Block_t + t;
            count = 0;  
            //Filling the block with the coordinates of the lattice points
            for(int x = x0; x < x1; x++) {
                for (int t = t0; t < t1; t++) {
                    Blocks[block][count++] = x * Nt+ t; 
                    //Each block also considers both spin components, 
                    //so we only reference the lattice coordinates here.
                }
            }
            if (count != lattice_sites_per_block) {
                std::cout << "Block " << block << " has " << count << " lattice points" << std::endl;
                std::cout << "Expected " << lattice_sites_per_block << std::endl;
                exit(1);
            }
            //Red-black decomposition for the blocks.
            if (Block_t % 2 == 0) {
                if  (x%2 ==0){
                    (block % 2 == 0) ? RedBlocks[block / 2] = block:BlackBlocks[block / 2] = block; 
                }
                else{
                    (block % 2 == 0) ? BlackBlocks[block / 2] = block:RedBlocks[block / 2] = block; 
                }
            } 
            else {
                (block % 2 == 0) ? RedBlocks[block / 2] = block:BlackBlocks[block / 2] = block;                
            }
        }
    }
}

//A_B = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
void SAP_C::I_D_B_1_It(const spinor& v, spinor& x,const int& block){
    bool print_message = false; //good for testing GMRES   

    spinor temp(lattice_sites_per_block, c_vector(spins*colors, 0)); 

    //temp = I_B^T v
    for (int j = 0; j < lattice_sites_per_block; j++){
        //Writing result to x 
        for(int alf = 0; alf<spins; alf++){
        for(int c = 0; c<colors; c++){
            temp[j][spins*c+alf] = v[Blocks[block][j]][spins*c+alf];
        }
        }
    }
     
    set_zeros(x,lattice_sites_per_block,spins*colors); //Initialize x to zero
    gmres_DB.set_block(block); //Set the block index for the GMRES_D_B operator
    gmres_DB.fgmres(temp,temp,x, print_message,false); //Call the GMRES solver 
}

int SAP_C::SAP(const spinor& v,spinor &x, const int& nu, const int& blocks_per_proc,const bool& print){
    /*
    Solves D x = v using the SAP method
    */
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    //Divide SAP_RedBlocks among processes
    int start = rank * blocks_per_proc;
    int end = std::min(start + blocks_per_proc, coloring_blocks);

    spinor temp(lattice_sites_per_block, c_vector(spins*colors, 0)); 
    spinor r(Ntot, c_vector(spins*colors, 0)); //residual
    spinor Dphi(Ntot, c_vector(spins*colors, 0)); //Temporary spinor for D x
    funcGlobal(x,Dphi);
    axpy(v,Dphi,-1.0,r); //r = v - D x

    //Prepare buffers for MPI communication
    c_vector local_buffer(Ntot * spins * colors, 0);
    c_vector global_buffer(Ntot * spins * colors, 0);
    
    for (int i = 0; i< nu; i++){  
        for(int n = 0; n < Ntot * spins * colors; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }
        //For the coarser levels, when the number of SAP blocks is smaller than for the finest level, notice that
        //if we consider all the ranks from the finest level, we will have cases where start>end, i.e. they don't enter
        //in the for loop. This is convenient, because even though every rank will access this function, only the ranks that
        //satisfy start < end will play a role here. i.e. we ignore those ranks where rank>NBlocks 
        for (int b = start; b < end; b++) {
            int block = RedBlocks[b];
            I_D_B_1_It(r, temp, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < lattice_sites_per_block; n++) {
            for(int alf = 0; alf<spins; alf++){
            for(int c = 0; c<colors; c++){
                local_buffer[colors*spins*Blocks[block][n]+spins*c+alf] = temp[n][spins*c+alf];
            }
            }
            }
        }

        //------MPI communication for red blocks------//
        // Perform single allreduce
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * spins * colors, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        
        //---------------------------------------------//
        //x = x + global_x;
        for(int n = 0; n < Ntot; n++) {
        for(int alf = 0; alf<spins; alf++){
        for(int c = 0; c<colors; c++){
             x[n][spins*c+alf] += global_buffer[colors * spins * n + spins * c + alf]; //global_x[n][2c+alf];
        }
        }
        }

        funcGlobal(x,Dphi);
        //r = v - D x
        axpy(v,Dphi,-1.0,r);
        //r = v - D_phi(U, x, m0); //r = v - D x

        for(int n = 0; n < Ntot * spins * colors; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }

        for (int b = start; b < end; b++) {
            int block = BlackBlocks[b];
            I_D_B_1_It(r, temp, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < lattice_sites_per_block; n++) {
            for(int alf = 0; alf<spins; alf++){
            for(int c = 0; c<colors; c++){
                local_buffer[colors*spins*Blocks[block][n]+spins*c+alf] = temp[n][spins*c+alf];
            }
            }
            }
        }

        //------MPI communication for black blocks------//

        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * spins * colors, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        for(int n = 0; n < Ntot; n++) {
        for(int alf = 0; alf<spins; alf++){
        for(int c = 0; c<colors; c++){
             x[n][spins*c+alf] += global_buffer[colors * spins * n + spins * c + alf]; //global_x[n][2c+alf];
        }
        }
        }

        funcGlobal(x,Dphi);
        //r = v - D x
        axpy(v,Dphi,-1.0,r);

        err = sqrt(std::real(dot(r, r))); 
        if (err < tol * v_norm) {
            if (print == true)
                std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
    if (print == true)
        std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    
    return 0; 
}

void SAP_fine_level::D_B(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    int RightPB_0, blockRPB_0; //Right periodic boundary in the 0-direction
    int RightPB_1, blockRPB_1; //Right periodic boundary in the 1-direction
    int LeftPB_0, blockLPB_0; //Left periodic boundary in the 0-direction
    int LeftPB_1, blockLPB_1; //Left periodic boundary in the 1-direction

    c_vector phi_RPB_0 = c_vector(2, 0);
    c_vector phi_RPB_1 = c_vector(2, 0);
    c_vector phi_LPB_0 = c_vector(2, 0);
    c_vector phi_LPB_1 = c_vector(2, 0);

    std::vector<std::vector<int>> &LeftPB = LeftPB_l[0];
	std::vector<std::vector<int>> &RightPB = RightPB_l[0];
	c_matrix &SignR = SignR_l[0];
	c_matrix &SignL = SignL_l[0];

    for (int m = 0; m < lattice_sites_per_block; m++) {
		//n = x * Nt + t
        int n = Blocks[block][m]; //n is the index of the lattice point in the block
        
        //Get m and block for the neighbors
        getMandBlock(RightPB[n][0], RightPB_0, blockRPB_0);
        getMandBlock(RightPB[n][1], RightPB_1, blockRPB_1); 
        getMandBlock(LeftPB[n][0], LeftPB_0, blockLPB_0); 
        getMandBlock(LeftPB[n][1], LeftPB_1, blockLPB_1); 
        
        if(blockRPB_0 == block){phi_RPB_0 = v[RightPB_0];}
        else {phi_RPB_0[0] = 0; phi_RPB_0[1] = 0; }

        if(blockRPB_1 == block){phi_RPB_1 = v[RightPB_1];}
        else {phi_RPB_1[0] = 0; phi_RPB_1[1] = 0; }

        if(blockLPB_0 == block){phi_LPB_0 = v[LeftPB_0];}
        else {phi_LPB_0[0] = 0; phi_LPB_0[1] = 0; }

        if(blockLPB_1 == block) {phi_LPB_1 = v[LeftPB_1];}
        else {phi_LPB_1[0] = 0; phi_LPB_1[1] = 0; }
 
		x[m][0] = (m0 + 2) * v[m][0] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (phi_RPB_0[0] - phi_RPB_0[1]) 
		+   U[n][1] * SignR[n][1] * (phi_RPB_1[0] + I_number * phi_RPB_1[1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi_LPB_0[0] + phi_LPB_0[1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (phi_LPB_1[0] - I_number * phi_LPB_1[1])
		);

		x[m][1] = (m0 + 2) * v[m][1] - 0.5 * ( 
			U[n][0] * SignR[n][0] * (-phi_RPB_0[0] + phi_RPB_0[1]) 
		+   U[n][1] * SignR[n][1] * (-I_number*phi_RPB_1[0] + phi_RPB_1[1])
		+   std::conj(U[LeftPB[n][0]][0]) * SignL[n][0] * (phi_LPB_0[0] + phi_LPB_0[1])
		+   std::conj(U[LeftPB[n][1]][1]) * SignL[n][1] * (I_number*phi_LPB_1[0] + phi_LPB_1[1])
		);
			
	}   
}
//One color and two spins per lattice, i.e. two degrees of freedom at level 0 
SAP_fine_level sap(LV::Ntot,  2, SAPV::sap_tolerance, LV::Nt, LV::Nx,SAPV::sap_block_x,SAPV::sap_block_t,2,1);

/*
This is the explicit definition of the interpolation and restriction operators for the SAP method.
I leave it here for reference, but it is not used in the code.

//x = I_B^T v --> Restriction of the vector v to the block B
//dim(v) = 2 Ntot, dim(x) = 2 * sap_lattice_sites_per_block
void It_B_v(const spinor& v, spinor& x, const int& block){
    using namespace SAPV;
    set_zeros(x,sap_lattice_sites_per_block,2); //Initialize the output vector to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        x[j][0] = v[SAP_Blocks[block][j]][0];
        x[j][1] = v[SAP_Blocks[block][j]][1];
    }
}

// x = I_B v --> Interpolation of the vector v to the original lattice
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 Ntot
void I_B_v(const spinor& v, spinor& x,const int& block){
    using namespace SAPV;
    set_zeros(x,Ntot,2); //Initialize x to zero
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        x[SAP_Blocks[block][j]][0] += v[j][0];
        x[SAP_Blocks[block][j]][1] += v[j][1];
    }

}

*/
