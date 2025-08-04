#include "sap.h"

void SchwarzBlocks(){
    using namespace SAPV; //SAP parameters namespace
    schwarz_blocks = true; //Schwarz blocks are initialized after calling this function
    int count, block;
    int x0, t0, x1, t1;
    for (int x = 0; x < sap_block_x; x++) {
        for (int t = 0; t < sap_block_t; t++) {
            x0 = x * sap_x_elements; t0 = t * sap_t_elements;
            x1 = (x + 1) * sap_x_elements; t1 = (t + 1) * sap_t_elements;
            block = x * sap_block_t + t;
            count = 0;  
            //Filling the block with the coordinates of the lattice points
            for(int x = x0; x < x1; x++) {
                for (int t = t0; t < t1; t++) {
                    SAP_Blocks[block][count++] = x * Nt+ t; 
                    //Each block also considers both spin components, 
                    //so we only reference the lattice coordinates here.
                }
            }
            if (count != sap_lattice_sites_per_block) {
                std::cout << "Block " << block << " has " << count << " lattice points" << std::endl;
                std::cout << "Expected " << sap_lattice_sites_per_block << std::endl;
                exit(1);
            }
            //Red-black decomposition for the blocks.
            if (sap_block_t % 2 == 0) {
                if  (x%2 ==0){
                    (block % 2 == 0) ? SAP_RedBlocks[block / 2] = block:SAP_BlackBlocks[block / 2] = block; 
                }
                else{
                    (block % 2 == 0) ? SAP_BlackBlocks[block / 2] = block:SAP_RedBlocks[block / 2] = block; 
                }
            } 
            else {
                (block % 2 == 0) ? SAP_RedBlocks[block / 2] = block:SAP_BlackBlocks[block / 2] = block;                
            }
        }
    }
}


//Restriction of the Dirac operator to the block B
//D_B v = I_B^T D I_B v  
//dim(v) = 2 * sap_lattice_sites_per_block, dim(x) = 2 * sap_lattice_sites_per_block
void D_B(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    using namespace SAPV;
    int RightPB_0, blockRPB_0; //Right periodic boundary in the 0-direction
    int RightPB_1, blockRPB_1; //Right periodic boundary in the 1-direction
    int LeftPB_0, blockLPB_0; //Left periodic boundary in the 0-direction
    int LeftPB_1, blockLPB_1; //Left periodic boundary in the 1-direction

    c_vector phi_RPB_0 = c_vector(2, 0);
    c_vector phi_RPB_1 = c_vector(2, 0);
    c_vector phi_LPB_0 = c_vector(2, 0);
    c_vector phi_LPB_1 = c_vector(2, 0);

    for (int m = 0; m < sap_lattice_sites_per_block; m++) {
		//n = x * Nt + t
        int n = SAP_Blocks[block][m]; //n is the index of the lattice point in the block
        
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


GMRES_D_B gmres_DB(SAPV::sap_lattice_sites_per_block, 2, SAPV::sap_gmres_restart_length, SAPV::sap_gmres_restarts, 
    SAPV::sap_gmres_tolerance); //Define the GMRES solver for D_B

//A_B = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
//dim(v) = 2 * Ntot, dim(x) = 2 Ntot
//v: input, x: output
void I_D_B_1_It(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block){
    using namespace SAPV;
    bool print_message = false; //good for testing GMRES   

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 

    //temp = I_B^T v
    for (int j = 0; j < sap_lattice_sites_per_block; j++){
        //Writing result to x 
        temp[j][0] = v[SAP_Blocks[block][j]][0];
        temp[j][1] = v[SAP_Blocks[block][j]][1];
    }
     
    set_zeros(x,sap_lattice_sites_per_block,2); //Initialize x to zero

    gmres_DB.set_block(block); //Set the block index for the GMRES_D_B operator
    gmres_DB.fgmres(temp,temp,x, print_message); //Call the GMRES solver 
}

int SAP(const c_matrix& U, const spinor& v,spinor &x, const double& m0,const int& nu,const int& blocks_per_proc){
    /*
    Solves D x = v using the SAP method
    */
   using namespace SAPV;
   int size, rank;
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    double err;
    double v_norm = sqrt(std::real(dot(v, v))); //norm of the right hand side

    //Divide SAP_RedBlocks among processes
    int start = rank * blocks_per_proc;
    int end = std::min(start + blocks_per_proc, sap_coloring_blocks);

    spinor temp(sap_lattice_sites_per_block, c_vector(2, 0)); 
    spinor r(Ntot, c_vector(2, 0)); //residual

    spinor Dphi(Ntot, c_vector(2, 0)); //Temporary spinor for D x
    D_phi(U, x, Dphi,m0);
    axpy(v,Dphi,-1.0,r); //r = v - D x


    //Prepare buffers for MPI communication
    c_vector local_buffer(Ntot * 2, 0);
    c_vector global_buffer(Ntot * 2, 0);
    
    for (int i = 0; i< nu; i++){  
        for(int n = 0; n < Ntot * 2; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }
        for (int b = start; b < end; b++) {
            int block = SAP_RedBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                local_buffer[2*SAP_Blocks[block][n]]     = temp[n][0];
                local_buffer[2*SAP_Blocks[block][n] + 1] = temp[n][1];
            }
        }

        //------MPI communication for red blocks------//
        // Perform single allreduce
        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        
        //---------------------------------------------//
        //x = x + global_x;

        for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_buffer[2*n]; //global_x[n][0];
            x[n][1] += global_buffer[2*n+1]; //global_x[n][1];
        }

        D_phi(U, x, Dphi,m0);
        //r = v - D x
        axpy(v,Dphi,-1.0,r);
        //r = v - D_phi(U, x, m0); //r = v - D x


        //set_zeros(local_x,Ntot,2); //Initialize local_x to zero
        for(int n = 0; n < Ntot * 2; n++) {
            local_buffer[n] = 0.0; //Initialize local_buffer to zero
        }

        for (int b = start; b < end; b++) {
            int block = SAP_BlackBlocks[b];
            I_D_B_1_It(U, r, temp, m0, block);
            //local_x = local_x + temp; // Local computation
            for(int n = 0; n < sap_lattice_sites_per_block; n++) {
                local_buffer[2*SAP_Blocks[block][n]]     = temp[n][0];
                local_buffer[2*SAP_Blocks[block][n] + 1] =  temp[n][1];
            }
        }

        //------MPI communication for black blocks------//

        MPI_Allreduce(local_buffer.data(), global_buffer.data(), Ntot * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        for(int n = 0; n < Ntot; n++) {
            x[n][0] += global_buffer[2*n]; //global_x[n][0];
            x[n][1] += global_buffer[2*n+1]; //global_x[n][1];
        }

        D_phi(U, x, Dphi,m0);
        //r = v - D x
        axpy(v,Dphi,-1.0,r);

        err = sqrt(std::real(dot(r, r))); 
        if (err < sap_tolerance * v_norm) {
            //std::cout << "SAP converged in " << i << " iterations, error: " << err << std::endl;
            return 1;
        }
    }
   //std::cout << "SAP did not converge in " << nu << " iterations, error: " << err << std::endl;
    
    return 0; 
}

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
