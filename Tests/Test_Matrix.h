#ifndef CPP_RZ_PIC_TEST_MATRIX_H
#define CPP_RZ_PIC_TEST_MATRIX_H

#include "../Tools/Matrix.h"
#include "../Tools/Logger.h"
#include <fstream>


void test_Matrix() {
    cout << "test_Matrix:" << endl;
    int N = 2;
    scalar value = 1;
    const Matrix a(N, N, value), b(N, N);
    Matrix c(N, N);
    a.print();
    b.print();
    c = a + b + a;
    c.print();
    //c = c * 4;
    c = c / 4 * 4;
    c.print();
    c += a;
    c.print();
    c -= a;
    c.print();
    c /= 4;
    c.print();
    //element_logging(a, "matrix_test.txt", "\n", ios::app);
    //out << a << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_MATRIX_H
