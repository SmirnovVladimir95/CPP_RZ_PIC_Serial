#ifndef CPP_RZ_PIC_TEST_NEUTRALGAS_H
#define CPP_RZ_PIC_TEST_NEUTRALGAS_H

#include <iostream>
#include "../ElementaryProcesses/NeutralGas.h"
using namespace std;

void test_NeutralGas() {
    double n = 1e-18, mAr = 6.6335209e-26, T = 500;
    NeutralGas gas(n, mAr, T);
    array<double, 3> vel = gas.generate_velocity();
    cout << vel[0] << " " << vel[1] << " " << vel[2] << endl;
    cout << sqrt(1.38e-23*T/mAr) << endl;
}

#endif //CPP_RZ_PIC_TEST_NEUTRALGAS_H
