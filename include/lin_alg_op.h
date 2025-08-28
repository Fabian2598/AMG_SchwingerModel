#ifndef LIN_ALG_OP_H
#define LIN_ALG_OP_H
#include <vector>
#include <complex>
#include <iostream>
//#include <cblas.h>
//#include <omp.h>

//Linear algebra operations

typedef std::complex<double> c_double;
typedef std::vector<c_double> c_vector;
typedef std::vector<c_vector> c_matrix; 
typedef std::vector<c_vector> spinor;


/*
    dot product between two complex vectors 
*/
inline c_double dot(const c_vector& x, const c_vector& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        z += x[i] * std::conj(y[i]);
    }
    return z;
}

/*
    dot product between two spinors of the form psi[ntot][2]
    A.B = sum_i A_i conj(B_i) 
*/
inline c_double dot(const spinor& x, const spinor& y) {
    c_double z = 0;
    for (int i = 0; i < x.size(); i++) {
        for (int j = 0; j < x[i].size(); j++) {
            z += x[i][j] * std::conj(y[i][j]);
        }
    }
    return z;
}

/*
    Scalar times a complex vector
*/
template <typename T>
inline void scal(const T& lambda, const c_vector& X, c_vector& Y) {
    // Y = lambda X
    int size = X.size();
    for (int i = 0; i < size; i++) {
        Y[i] = lambda * X[i];
    }
}

/*
    Complex vector addition
*/
template <typename T>
inline void axpy(const c_vector& X, const c_vector& Y, const T&lambda,  c_vector& out) {
    int size = X.size();
    for (int i = 0; i < size; i++) {
        out[i] = X[i] + lambda * Y[i];
    }
}

/*
    Matrix-vector multiplication 
*/
inline void AtimesV(const c_matrix& A, const c_vector& v, c_vector& w) {
    //Matrix-vector multiplication
    // A * v = w, 
    int size1 = A.size();
    int size2 = A[0].size();
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            w[i] += A[i][j] * v[j];
        }
    }
}

/*
    Scalar times a complex matrix.
    Also works for spinors, since they are just matrices with 2 columns.
*/

template <typename T>
inline void scal(const T& lambda, const c_matrix& X, c_matrix& Y) {
    // alpha times a complex matrix
    int size1 = X.size();
    int size2 = X[0].size();
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            Y[i][j] = lambda * X[i][j];
        }
    }
}

/*
    Complex matrix addition
    Also works for spinors, since they are just matrices with 2 columns.
*/
template <typename T>
inline void axpy(const c_matrix& X, const c_matrix& Y, const T& lambda, c_matrix& out) {
    // A + B = C
    int size1 = X.size();
    int size2 = X[0].size();
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            out[i][j] = X[i][j] + lambda * Y[i][j];
        }
    }
}


/*
    Print a complex matrix
    Also works for spinors
*/
inline void PrintComplexMatrix(const c_matrix& v ){
    for(int i = 0; i < v.size(); i++){
        for(int j = 0; j < v[i].size(); j++){
            std::cout << v[i][j] << " ";
        }
		std::cout << std::endl;
    }
    std::cout << std::endl;
}

/*
    Print a complex vector
*/
inline void PrintComplexVector(const c_vector& v ){
    for(int i = 0; i < v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}


/*
	Normalize a spinor.
*/
inline void normalize(spinor& v){
	c_double norm = sqrt(std::real(dot(v,v))) + 0.0*c_double(0,1); 
	scal(1.0/norm, v, v); //v = v / norm
}


#endif 