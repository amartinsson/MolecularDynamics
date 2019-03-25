#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <vector>
#include <iostream>
#include <math.h>
#include <omp.h>

#include "Generator.hpp"

using namespace std;

// ------------------------------------------------------------------------- //
//                                MATRIX CLASS
// ------------------------------------------------------------------------- //
class Matrix
{

public:
    Matrix() : dimx(0), dimy(0), M(0) {};
    // make matrix with vector
    explicit Matrix(const int& n, const int& m);
    explicit Matrix(vector<double> m);
    // initialises an upper triangular matrix of random numbers
    explicit Matrix(const int& n, const int& m, NormalGenerator& gen);
    //explicit Matrix(const Matrix& m);
    ~Matrix()
    {
        // for(unsigned i=0; i < M.size(); i++)
        // {
        //     delete &M[i];
        // }
        // M.clear();
        // delete &dimx;
        // delete &dimy;
    };

    // returns the size of the matrix
    vector<int> size() const;

    // setting two matrices equal to eachother
    void operator= (const Matrix& m);
    void operator*= (const double& value);
    Matrix operator* (const double& value);
    Matrix operator/ (const double& value);

    // operators for getting elements
    double operator() (const int& i, const int& j) const;
    double& operator() (const int& i, const int& j);

    // transpose of matrix
    Matrix T() const;

    // determinatn of 2x2 matrix
    double det() const;

    // inverse of 2 x 2 matrix
    Matrix inv() const;

    // multiply two matrices
    Matrix operator* (const Matrix& m);

    // adding matrices
    void operator+= (const Matrix& m);
    void operator-= (const Matrix& m);

    // subtracting Matrix
    Matrix operator+ (const Matrix& m);
    Matrix operator- (const Matrix& m);

    // vector operations
    class Vector operator* (const Vector& v);

    // zero matrix
    Matrix sqrt() const;
    // Symmetrix matrix square root
    Matrix symsqrt() const;

    void zero();
    void diag(const double& value);

protected:
    vector<double> M;
    vector<double> Minv;
    int dimx;
    int dimy;

    // simple get function to simplify getting element
    double get(const int& i, const int& j) const {return M[j + i * dimy];}
    double& get(const int& i, const int& j) {return M[j + i * dimy];}
};

// ------------------------------------------------------------------------- //
//                                VECTOR CLASS
// ------------------------------------------------------------------------- //
class Vector
{
public:
    Vector() : dimx(0), vec(0) {};
    // constructor with values
    explicit Vector(const int& dimx, const double& val);
    explicit Vector(const int& dim);
    explicit Vector(const int& dim, NormalGenerator& gen);

    // explicit Vector(const Vector& V);
    ~Vector(){};

    // return the size
    int size() const;

    // setting two matrices equal to eachother
    //Vector* operator*();
    // void operator= (const Vector& other) {*this = other;}
    void operator= (const Vector& other);
    void operator+= (const Vector& other);
    void operator-= (const Vector& other);

    void operator*= (const double& value);
    Vector operator* (const double& value);

    // addition and subtraction
    Vector operator+ (const Vector& other) const;
    Vector operator- (const Vector& other) const;

    // getting vector
    double operator() (const int& i) const;
    double& operator() (const int& i);

    // multiply two vectors
    double dot(const Vector& other);
    Matrix out(const Vector& other);

    // multiply vector with matrix
    Vector operator* (const Matrix& Mat);

    // transpose vector
    Vector T() const;
    Vector& _T();

    // zero entire vector
    void zero();

    // check if vector is transposed
    bool _is_transpose() const {return transpose;}

private:
    vector<double> vec;
    int dimx;
    bool transpose;
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
