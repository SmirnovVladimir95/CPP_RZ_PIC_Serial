#ifndef CPP_RZ_PIC_MATRIX_H
#define CPP_RZ_PIC_MATRIX_H

#include <iostream>
#include <vector>
#include "ProjectTypes.h"
using namespace std;

//template <class Type>
class Matrix {
private:
    size_t _rows;
    size_t _columns;
    //unique_ptr<scalar[]> data; // c++14 and higher
    vector<scalar> data;
public:
    Matrix() : _rows(0), _columns(0) {};
    Matrix(size_t rows, size_t columns);
    Matrix(size_t rows, size_t columns, scalar value);
    size_t rows() const;
    size_t columns() const;
    scalar* data_ptr();
    scalar& operator()(size_t row, size_t column);
    scalar operator()(size_t row, size_t column) const;
    Matrix operator+(const Matrix& other) const;
    Matrix& operator+=(const Matrix& other);
    Matrix& operator-=(const Matrix& other);
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(scalar value) const;
    Matrix operator/(scalar value) const;
    Matrix& operator/=(scalar value);
    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix);
    void print() const;
    void fill(scalar value);
    void copy(Matrix& matrix);
    void resize(size_t rows, size_t columns);
};


#endif //CPP_RZ_PIC_MATRIX_H
