#include "Array.hpp"

// default constructor
Vector::Vector(const int& dim) : dimx(dim), transpose(false)
{
    vec.resize(dim, 0.0);
}

// second constructor
Vector::Vector(const int& dim, const double& val) : dimx(dim), transpose(false)
{
    vec.resize(dim, val);
}

// third constructor
Vector::Vector(const int& dim, NormalGenerator& norm) : dimx(dim), transpose(false)
{
    vec.resize(dimx);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
        vec[i] = norm();
}

// Vector::Vector(const Vector& V)
// {
//     dimx = V.size();
//     vec.resize(dimx, 0.0);
//
// #pragma omp simd collapse
//     for(unsigned i=0; i<dimx; i++)
//         vec[i] = V(i);
// }

// return the size
int Vector::size() const
{
    int size = vec.size();
    return size;
}

// getting vector
double Vector::operator() (const int& i) const
{
    return vec[i];
}

// getting vector
double& Vector::operator() (const int& i)
{
    return vec[i];
}

// zero entire vector
void Vector::zero()
{
    std::fill(vec.begin(), vec.end(), 0.0);
}

// return ptr to the
// Vector* Vector::operator*() {return this};

// multiply two vectors together
double Vector::dot(const Vector& other)
{
    if(other.size() != dimx)
    {
        printf("Error Vector dimensions don't agree in dot product!\n");
        exit(-1);
    }

    if(other._is_transpose())
    {
        printf("Error Second vector is transposed in dot!\n");
        exit(-1);
    }

    double dot = 0.0;

    #pragma omp simd collapse(1)
        for(unsigned i=0; i<dimx; i++)
            dot += vec[i] * other(i);

    return dot;
}

Matrix Vector::out(const Vector& other)
{
    if(other.size() != dimx)
    {
        printf("Error Vector dimensions don't agree in outer product!\n");
        exit(-1);
    }

    if(other._is_transpose())
    {
        printf("Error Second vector is transposed in outer!\n");
        exit(-1);
    }

    Matrix outer(dimx, dimx);

    #pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimx; j++)
            outer(i, j) = vec[i] * other(j);

    return outer;
}

// multiply left multiply vector with matrix
Vector Vector::operator* (const Matrix& Mat)
{
    if(!transpose)
    {
        printf("Error Cannot left multiply untransposed Vector by Matrix!\n");
        exit(-1);
    }

    if(dimx != Mat.size()[0])
    {
        printf("Error Matrix and Vector Dimnesions dont agree!\n");
        exit(-1);
    }

    // dereference dimension
    int dimy = Mat.size()[1];

    // make return matrix
    Vector Vret(dimx);

#pragma omp simd collapse(2)
    for(unsigned i=0; i<dimx; i++)
        for(unsigned j=0; j<dimy; j++)
                Vret(j) += vec[i] * Mat(i, j);

    // return the new matrix
    return Vret.T();
}

// set vectors to be equal
void Vector::operator= (const Vector& other)
{
    if(this != &other)
    {
        if(other.size() != dimx && dimx != 0)
        {
            printf("Error cannot set vectors equal!\n");
            exit(-1);
        }
        else
        {
            this->dimx = other.size();
            vec.resize(dimx);
        }

#pragma omp simd collapse(1)
        for(unsigned i=0; i<dimx; i++)
            vec[i] = other(i);
    }
}

// addition and subtraction
Vector Vector::operator+ (const Vector& other) const
{
    if(other.size() != dimx)
    {
        printf("Error cannot set vectors equal!\n");
        exit(-1);
    }

    Vector retVec(dimx);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
        retVec(i) = vec[i] + other(i);

    return retVec;
}

// addition and subtraction
Vector Vector::operator- (const Vector& other) const
{
    if(other.size() != dimx)
    {
        printf("Error cannot subtract unequal vectors!\n");
        exit(-1);
    }

    Vector retVec(dimx);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
        retVec(i) = vec[i] - other(i);

    return retVec;
}

// set vectors to be equal
void Vector::operator+= (const Vector& other)
{
    if(other.size() != dimx)
    {
        printf("Error cannot set vectors equal!\n");
        exit(-1);
    }

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
            vec[i] += other(i);
}

// set vectors to be equal
void Vector::operator*= (const double& value)
{
#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
            vec[i] *= value;
}

// set vectors to be equal
Vector Vector::operator* (const double& value)
{
    Vector retVec(dimx);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
        retVec(i) = vec[i] * value;

    return retVec;
}

// set vectors to be equal
void Vector::operator-= (const Vector& other)
{
    if(other.size() != dimx)
    {
        printf("Error cannot set vectors equal!\n");
        exit(-1);
    }

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
            vec[i] -= other(i);
}

// return copy of transposed Matrix
Vector Vector::T() const
{
    Vector retVec(dimx);

#pragma omp simd collapse(1)
    for(unsigned i=0; i<dimx; i++)
        retVec(i) = vec[i];

    return retVec._T();
}

// set transpose to true for vector
Vector& Vector::_T()
{
    if(!transpose)
        transpose = true;
    else
        transpose = false;

    return *this;
}
