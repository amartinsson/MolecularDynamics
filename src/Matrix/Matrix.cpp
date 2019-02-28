#include "Matrix.hpp"

using namespace std;

// ------------------------------------------------------------------------- //
//                                MATRIX CLASS
// ------------------------------------------------------------------------- //
Matrix::Matrix(const int& n, const int& m) : dimx(n), dimy(m)
{
    M.resize(dimx * dimy, 0.0);
}

// Make matrix from vector
Matrix::Matrix(vector<double> m)
{
    if(m.size() % 2 != 0)
    {
        printf("ERROR: Vector must be divisible by 2!\n");
        exit(-1);
    }

    // set the matrix
    M.resize(m.size(), 0.0);
    M = m;

    // set the dimensions
    dimx = M.size() / 2;
    dimy = M.size() / 2;
}

// get  the size of the matrix
vector<int> Matrix::size()
{
    vector<int> size(2, 0.0);
    size[0] = dimx;
    size[1] = dimy;
    return size;
}

// operator for getting elements
double Matrix::operator() (const int& i, const int& j) const
{
    return M[j + i * dimy];
}

// return the adress of positon for assigning to operator
double& Matrix::operator() (const int& i, const int& j)
{
    return M[j + i * dimy];
}

// Setting two matrices equal to eachother
void Matrix::operator= (Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error cannot set matrices equal!\n");
        exit(-1);
    }

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            M[j + i * dimy] = m(i, j);
}

// Matrix Matrix multiplication
Matrix Matrix::operator* (Matrix m)
{
    if(m.size()[0] != dimy)
    {
        printf("Matrix dimensions must agree!\n");
        exit(-1);
    }

    // make return matrix
    Matrix Mret(dimx, m.size()[1]);

#pragma omp simd collapse(3)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<m.size()[1]; j++)
            for(unsigned k=0; k<dimy; k++)
                Mret(i,j) += M[k + i * dimy] * m(k, j);

    // return the new matrix
    return Mret;
}

// calculate the transpose of a matrix
Matrix Matrix::T()
{
    // make matrix to return
    Matrix Mret(dimy, dimx);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(j, i) = M[j + i * dimy];

    return Mret;
}

double Matrix::det()
{
    if(dimx > 2 or dimy > 2)
    {
        printf("Error: Inverse only works for 2x2 matrix!\n");
        exit(-1);
    }

    if(dimx != dimy)
    {
        printf("Error: Can only do inverse for Square Matrix\n");
        exit(-1);
    }

    return M[0] *  M[1 + dimy] -  M[1] * M[dimy];
}

// calculate the invers of a 2x2 matrix
Matrix Matrix::inv()
{
    if(dimx > 2 or dimy > 2)
    {
        printf("Error: Inverse only works for 2x2 matrix!\n");
        exit(-1);
    }

    if(dimx != dimy)
    {
        printf("Error: Can only do inverse for Square Matrix\n");
        exit(-1);
    }

    Matrix Mret(dimx, dimy);

    double det = M[0] *  M[1 + dimy] -  M[1] * M[dimy];

    Mret(0, 0) = M[1 + dimy] / det;
    Mret(0, 1) = -M[1] / det;

    Mret(1, 0) = -M[dimy] / det;
    Mret(1, 1) = M[0] / det;

    return Mret;
}

// ------------------------------------------------------------------------- //
//                        ROTATION MATRIX CLASS
// ------------------------------------------------------------------------- //
RotMatrix::RotMatrix(const double& alpha) : Matrix(2, 2)
{
    // set the rotation matrix
    M[0] = cos(alpha);
    M[1] = sin(alpha);
    M[dimy] = -sin(alpha);
    M[1 + dimy] = cos(alpha);
}
