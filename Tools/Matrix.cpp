#include "Matrix.h"

Matrix::Matrix(size_t rows, size_t columns) : _rows(rows), _columns(columns) {
    data.resize(rows*columns);
}

Matrix::Matrix(size_t rows, size_t columns, scalar value) : _rows(rows), _columns(columns) {
    data.resize(rows*columns, value);
}

size_t Matrix::rows() const { return _rows; }

size_t Matrix::columns() const { return _columns; }

scalar* Matrix::data_ptr() { return data.data(); }

scalar& Matrix::operator()(size_t row, size_t column) {
    if (row >= _rows or column >= _columns) {
        cout << "Error, index (" << row << " " << column << ") out of range (return 1st element)" << endl;
        throw;
    }
    return data[row * _columns + column];
}

scalar Matrix::operator()(size_t row, size_t column) const {
    if (row >= _rows or column >= _columns) {
        cout << "Error, index (" << row << " " << column << ") out of range (return 1st element)" << endl;
        throw;
    }
    return data[row * _columns + column];
}

void Matrix::print() const {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            cout << this->operator()(i, j) << " ";
        }
        cout << endl;
    }
}

void Matrix::fill(const scalar value) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            this->operator()(i, j) = value;
        }
    }
}

void Matrix::copy(Matrix& matrix) {
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            this->operator()(i, j) = matrix(i, j);
        }
    }
}

void Matrix::resize(size_t rows, size_t columns) {
    data.resize(rows*columns);
    _rows = rows;
    _columns = columns;
}

Matrix Matrix::operator+(const Matrix &other) const {
    Matrix tmp(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            tmp(i, j) = this->operator()(i, j) + other(i, j);
        }
    }
    return tmp;
}

Matrix Matrix::operator-(const Matrix &other) const {
    Matrix tmp(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            tmp(i, j) = this->operator()(i, j) - other(i, j);
        }
    }
    return tmp;
}

std::ostream& operator<<(std::ostream &out, const Matrix &matrix) {
    for (int i = 0; i < matrix.rows(); i++) {
        for (int j = 0; j < matrix.columns(); j++) {
            out << matrix(i, j) << " ";
        }
        out << endl;
    }
    return out;
}

Matrix Matrix::operator/(scalar value) const {
    Matrix tmp(_rows, _columns);
    assert(value != 0);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            tmp(i, j) = this->operator()(i, j) / value;
        }
    }
    return tmp;
}

Matrix Matrix::operator*(scalar value) const {
    Matrix tmp(_rows, _columns);
    for (int i = 0; i < _rows; i++) {
        for (int j = 0; j < _columns; j++) {
            tmp(i, j) = this->operator()(i, j) * value;
        }
    }
    return tmp;
}
