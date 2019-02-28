#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>
#include <math.h>
#include <omp.h>

using namespace std;

// ------------------------------------------------------------------------- //
//                                MATRIX CLASS
// ------------------------------------------------------------------------- //
class Matrix
{
public:
    Matrix(const int& n, const int& m);
    ~Matrix(){};

    // make matrix with vector
    explicit Matrix(vector<double> m);

    // returns the size of the matrix
    vector<int> size();

    // setting two matrices equal to eachother
    void operator= (Matrix& m);

    // operators for getting elements
    double operator() (const int& i, const int& j) const;
    double& operator() (const int& i, const int& j);

    // transpose of matrix
    Matrix T();

    // determinatn of 2x2 matrix
    double det();

    // inverse of 2 x 2 matrix
    Matrix inv();

    // multiply two matrices
    Matrix operator* (Matrix m);

protected:
    vector<double> M;
    int dimx;
    int dimy;
};

// ------------------------------------------------------------------------- //
//                                MATRIX CLASS
// ------------------------------------------------------------------------- //
class RotMatrix : public Matrix
{
public:
    RotMatrix(const double& alpha);

};

#endif
