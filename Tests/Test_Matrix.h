#ifndef CPP_RZ_PIC_TEST_MATRIX_H
#define CPP_RZ_PIC_TEST_MATRIX_H

#include "../Tools/Matrix.h"
#include "../Tools/Logger.h"
#include <fstream>


void test_Matrix() {
    int N = 2;
    const Matrix a(N, N, 1), b(N, N);
    Matrix c(N, N);
    a.print();
    b.print();
    c = a + b + a;
    c.print();
    element_logging(a, "matrix_test.txt", " ", ios::app);
    //out << a << endl;
}

#endif //CPP_RZ_PIC_TEST_MATRIX_H
