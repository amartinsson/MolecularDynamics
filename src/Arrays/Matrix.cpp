#include "Array.hpp"

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

Matrix::Matrix(const int& n, const int& m, NormalGenerator& norm) :
        dimx(n), dimy(m)
{
    M.resize(dimx * dimy, 0.0);

    //#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=i; j<dimy; j++)
            this->get(i,j) = norm();
}

// Matrix::Matrix(const Matrix& m)
// {
//     dimx = m.size()[0];
//     dimy = m.size()[1];
//
//     // resize the vector
//     M.resize(dimx * dimy, 0.0);
//
// #pragma omp simd collapse(2)
//     for(unsigned i=0; i<dimx; i++)
//         for(unsigned j=0; j<dimy; j++)
//             this->get(i, j) = m(i, j);
// }

// get  the size of the matrix
vector<int> Matrix::size() const
{
    vector<int> size(2, 0.0);
    size[0] = dimx;
    size[1] = dimy;
    return size;
}

// operator for getting elements
double Matrix::operator() (const int& i, const int& j) const
{
    return this->get(i,j);
}

// return the adress of positon for assigning to operator
double& Matrix::operator() (const int& i, const int& j)
{
    return this->get(i, j);
}

// zero entire matrix
void Matrix::zero()
{
    std::fill(M.begin(), M.end(), 0.0);
}

// Setting two matrices equal to eachother
void Matrix::operator= (const Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error: Dimension do not agree for operator=!\n");
        exit(-1);
    }

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            this->get(i,j) = m(i, j);
}

// Setting two matrices equal to eachother
void Matrix::operator+= (const Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error: Dimension do not agree for operator+=!\n");
        exit(-1);
    }

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            this->get(i,j) += m(i, j);
}

// Setting two matrices equal to eachother
void Matrix::operator-= (const Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error: Dimension do not agree for operator-=!\n");
        exit(-1);
    }

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            this->get(i,j) -= m(i, j);
}

// subtracting two matrices
Matrix Matrix::operator- (const Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error: Dimension do not agree for operator-!\n");
        exit(-1);
    }

    Matrix Mret(dimx, dimy);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(i,j) = this->get(i,j) - m(i, j);

    return Mret;
}

// adding two matrices
Matrix Matrix::operator+ (const Matrix& m)
{
    if(m.size()[0] != dimx or m.size()[1] != dimy)
    {
        printf("Error cannot set matrices equal!\n");
        exit(-1);
    }

    Matrix Mret(dimx, dimy);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(i,j) = this->get(i,j) + m(i, j);

    return Mret;
}


// subtracting two matrices
Matrix Matrix::sqrt() const
{
    if(dimy != dimx)
    {
        printf("Error cannot calculate sqrt diagonal!\n");
        exit(-1);
    }

    Matrix Mret(dimx, dimy);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
            Mret(i,i) = std::sqrt(this->get(i,i));

    return Mret;
}

// Symmetrix matrix square root
Matrix Matrix::symsqrt() const
{
    // assumes that the matrix is of the form
    //   a b
    //   b c
    ////////////

    if(dimy != dimx)
    {
        printf("Error cannot calculate sym sqrt of non square Matrix!\n");
        exit(-1);
    }

    if(dimx != 2 or dimy != 2)
    {
        printf("Error: sym sqrt only works for 2x2 matrix!\n");
        exit(-1);
    }

    // cache the values
    double a = this->get(0,0);
    double b = this->get(0,1);
    double c = this->get(1,1);

    //  make the two matrices
    Matrix D(dimx, dimy);
    Matrix V(dimx, dimy);
    Matrix Mret(dimx, dimy);

    if(b * b > 0.0)
    {
        // calculate the root
        double root = std::sqrt( pow((a - c), 2.0) + 4.0 * pow(b, 2.0));

        // assign the diagonal values of D
        D(0,0) = std::sqrt(0.5 * (a + c - root));
        D(1,1) = std::sqrt(0.5 * (a + c + root));

        // assign the values of V
        V(0,0) = 0.5 * (a - c - root) / b;
        V(0,1) = 0.5 * (a - c + root) / b;
        V(1,0) = 1.0;
        V(1,1) = 1.0;

        Mret = V * D * V.inv();
    }
    else
    {
        Mret(0,0) = std::sqrt(a);
        Mret(1,1) = std::sqrt(c);
    }

    return Mret;
}


// set the diagonal of matrix
void Matrix::diag(const double& value)
{
    for(unsigned i=0; i<dimx; i++)
        this->get(i,i) = value;
}

// Matrix Matrix multiplication
Matrix Matrix::operator* (const Matrix& m)
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
                Mret(i,j) += this->get(i, k) * m(k, j);

    // return the new matrix
    return Mret;
}

// void Matrix::operator*= (const Matrix& m)
// {
//     if(m.size()[0] != dimy)
//     {
//         printf("Matrix dimensions must agree!\n");
//         exit(-1);
//     }
//
// #pragma omp simd collapse(3)
//     for(unsigned i=0; i<dimx; i++)
//         for(unsigned j=0; j<m.size()[1]; j++)
//             for(unsigned k=0; k<dimy; k++)
//                 this->get(i,j) += this->get(i, k) * m(k, j);
//                 this->get(i, k) *= m(k, j);
// }

// matrix vector multiplication
Vector Matrix::operator* (const Vector& v)
{
    if(v.size() != dimy)
    {
        printf("Dimensions must agree!\n");
        exit(-1);
    }

    if(v._is_transpose())
    {
        printf("Error Cannot right multiply transposed vector with Matrix!\n");
        exit(-1);
    }

    // make return matrix
    Vector Vret(dimy);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
                Vret(i) += this->get(i, j) * v(j);

    // return the new matrix
    return Vret;
}

void Matrix::operator*= (const double& value)
{
#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            this->get(i, j) *= value;
}

Matrix Matrix::operator* (const double& value)
{
    Matrix Mret(dimx, dimy);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(i, j) = value * this->get(i, j);

    return Mret;
}

Matrix Matrix::operator/ (const double& value)
{
    Matrix Mret(dimx, dimy);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(i, j) = this->get(i, j) / value;

    return Mret;
}

// calculate the transpose of a matrix
Matrix Matrix::T() const
{
    // make matrix to return
    Matrix Mret(dimy, dimx);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
            Mret(j, i) = this->get(i, j);

    return Mret;
}

double Matrix::det() const
{
    if(dimx != 2 or dimy != 2)
    {
        printf("Error: determinatn only works for 2x2 matrix!\n");
        exit(-1);
    }

    if(dimx != dimy)
    {
        printf("Error: Can only do Determinant for Square Matrix\n");
        exit(-1);
    }

    return M[0] *  M[1 + dimy] -  M[1] * M[dimy];
}

// calculate the invers of a 2x2 matrix
Matrix Matrix::inv() const
{
    if(dimx != 2 or dimy != 2)
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
