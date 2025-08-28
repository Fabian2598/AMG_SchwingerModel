#ifndef SAP_H
#define SAP_H
#include "fgmres.h"
#include "mpi.h"


/*
    Set all elements of a spinor to zero.
    v: input spinor, N1: first index dimension, N2: second index dimension
*/
inline void set_zeros(spinor& v, const int& N1, const int& N2) {
    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
            v[i][j] = 0;
        }
    }
}

/*
    SAP class 
    The local matrix-vector operation and the global matrix-vector operation have to be defined
    in a subclass.
*/
class SAP_C {
public: 
    //--------- Nested class for GMRES_D_B operator ---------//
    class GMRES_D_B : public FGMRES {
        public:GMRES_D_B(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol,SAP_C* parent) :
        FGMRES(dim1, dim2, m, restarts, tol), parent(parent) {

        };
        ~GMRES_D_B() { };

        //Set block for GMRES_D_B
        void set_block(const int& block_index) { 
            parent->blockMPI = block_index;
        }
        private:
        SAP_C* parent; //Pointer to the parent SAP_C object
        
        void func(const spinor& in, spinor& out) override { 
            parent->funcLocal(in, out);
        }
        void preconditioner(const spinor& in, spinor& out) override { 
            out = std::move(in); //No preconditioning
        }
    };
    GMRES_D_B gmres_DB;
    //--------------------------------------------------------//

    //------------------------------------//
    //Constructor
    SAP_C(const int& dim1, const int& dim2, const double& tol,const int& Nt, const int& Nx,const int& block_x,const int& block_t, const int& spins,
        const int& colors) :
    dim1(dim1), dim2(dim2), tol(tol), Nt(Nt), Nx(Nx), Block_x(block_x), Block_t(block_t), spins(spins), colors(colors),
    gmres_DB(Nx/block_x* Nt/block_t, 
        dim2, 
        SAPV::sap_gmres_restart_length, 
        SAPV::sap_gmres_restarts, 
        SAPV::sap_gmres_tolerance,
        this)  
    {
        
        x_elements = Nx/block_x, t_elements = Nt/block_t; //Number of elements in the x and t direction
        NBlocks = Block_x * Block_t; //Number of Schwarz blocks
        lattice_sites_per_block = x_elements * t_elements; //Number of lattice points in the block
        variables_per_block = spins * colors * lattice_sites_per_block; //Number of variables in the block 
        coloring_blocks = ( NBlocks != 1) ? NBlocks/2: 1; //Number of red or black blocks
        Ntot = Nt * Nx; //Total number of lattice points This is equal to dim1 ...

        //std::cout << "SAP_C initialized with " << NBlocks << " blocks, each with " << lattice_sites_per_block << " lattice points." << std::endl;
        //std::cout << "x_elements: " << x_elements << ", t_elements: " << t_elements << std::endl;
        //std::cout << "variables_per_block: " << variables_per_block << ", coloring_blocks: " << coloring_blocks << std::endl;

        Blocks = std::vector<std::vector<int>>(block_x*block_t, std::vector<int>(x_elements*t_elements, 0));
        RedBlocks = std::vector<int>(coloring_blocks, 0); //Red blocks
        BlackBlocks = std::vector<int>(coloring_blocks, 0); //Black blocks       
        SchwarzBlocks(); //Initialize the Schwarz blocks
     };

    /*
        Parallel version of the SAP method.
        Solves D x = v using the SAP method.
        v: right-hand side,
        x: output --> The initial given value to the method is considered the initial guess.
        nu: number of iterations,
        blocks_per_proc: number of blocks per process

        The convergence criterion is ||r|| < ||phi|| * tol
    */
    int SAP(const spinor& v,spinor &x,const int& nu, const int& blocks_per_proc,const bool& print);

    int dim1, dim2, nu;
    int Nt, Nx; //Dimensions of the lattice
    int Ntot;
    double tol;
    int Block_x, Block_t; //Block dimensions for the SAP method
    int x_elements, t_elements; //Number of elements in the x and t direction
    int NBlocks, lattice_sites_per_block, variables_per_block, coloring_blocks; 
    int blockMPI; //Current block index for the GMRES_D_B operator
    int colors, spins;
    //Number of blocks, lattice sites per block, variables per block and coloring blocks

    std::vector<std::vector<int>> Blocks; //SAP_Blocks[number_of_block][vectorized_coordinate of the lattice point]
    //The vectorization does not take into account the spin index, since both spin indices are in the same block.
    std::vector<int> RedBlocks; //Block index for the red blocks
    std::vector<int> BlackBlocks; //Block index for the black blocks

private: 
    /*
        Build the Schwarz blocks
        Function only has to be called once before using the SAP method.
    */
    void SchwarzBlocks();

    /*
        A_B v = I_B * D_B^-1 * I_B^T v --> Extrapolation of D_B^-1 to the original lattice.
        dim(v) = 2 * Ntot, dim(x) = 2 Ntot
        v: input, x: output 
    */
    void I_D_B_1_It(const spinor& v, spinor& x, const int& block);

    /*
        Matrix-vector global operation
        This is defined in the derived classes
    */
    virtual void funcGlobal(const spinor& in, spinor& out) = 0; 

    /*
        Matrix-vector operation restricted to the blocks
    */
    virtual void funcLocal(const spinor& in, spinor& out) = 0; 
        
};

//SAP for the finest level
class SAP_fine_level : public SAP_C {
public:
    SAP_fine_level(const int& dim1, const int& dim2, const double& tol,const int& Nt, const int& Nx,const int& block_x,const int& block_t,
    const int& spins, const int& colors) :
    SAP_C(dim1, dim2, tol, Nt, Nx, block_x, block_t,spins,colors) {
    }


    void set_params(const c_matrix& conf, const double& bare_mass){
        U = &conf;
        m0 = bare_mass;
    }

private: 
    const c_matrix* U; 
    double m0; 

    void funcGlobal(const spinor& in, spinor& out) override { 
        D_phi(*U, in, out,m0);
    }

    void D_B(const c_matrix& U, const spinor& v, spinor& x, const double& m0,const int& block);

    void funcLocal(const spinor& in, spinor& out) override { 
        D_B(*U, in, out, m0,blockMPI);
    }

    /*
        Given a lattice point index n, it returns the corresponding 
        SAP block index and the local index m within that block.
    */
    inline void getMandBlock(const int& n, int &m, int &block) {
        int x = n / LV::Nt; //x coordinate of the lattice point 
        int t = n % LV::Nt; //t coordinate of the lattice point
        //Reconstructing the block and m index from x and t
        int block_x = x / x_elements; //Block index in the x direction
        int block_t = t / t_elements; //Block index in the t direction
        block = block_x * Block_t + block_t; //Block index in the SAP method

        int mx = x % x_elements; //x coordinate in the block
        int mt = t % t_elements; //t coordinate in the block
        m = mx * t_elements + mt; //Index in the block
    }

};

/*
    SAP solver for D_B
    Defined on sap.cpp
*/
extern SAP_fine_level sap;


/*  
    FGMRES with SAP preconditioner
    This method solves the original Dirac equation
*/
class FGMRES_SAP : public FGMRES {
    public:
    FGMRES_SAP(const int& dim1, const int& dim2, const int& m, const int& restarts, const double& tol,
    const c_matrix& U, const double& m0) : FGMRES(dim1, dim2, m, restarts, tol), U(U), m0(m0), dim1(dim1), dim2(dim2)
    {
    };
    ~FGMRES_SAP() { };
    
private:
    const c_matrix& U; //reference to Gauge configuration. This is to avoid copying the matrix
    const double& m0; //reference to mass parameter
    const int &dim1;
    const int &dim2;
    /*
    Implementation of the function that computes the matrix-vector product for the fine level
    */
    void func(const spinor& in, spinor& out) override {
        D_phi(U, in, out, m0); 
    }


    void preconditioner(const spinor& in, spinor& out) override {
        //Initialize ZmT[j] to zero
        for(int i = 0; i<dim1; i++){
            for(int j = 0; j<dim2; j++){
                out[i][j] = 0;
            }
        }
        //Defined outside of the class
        sap.SAP(in,out,1, SAPV::sap_blocks_per_proc,false);
    }
};


#endif