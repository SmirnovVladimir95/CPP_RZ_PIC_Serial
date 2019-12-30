#ifndef CPP_RZ_PIC_TEST_NEUTRALGAS_H
#define CPP_RZ_PIC_TEST_NEUTRALGAS_H

#include <iostream>
#include "../ElementaryProcesses/NeutralGas.h"
using namespace std;

void test_NeutralGas() {
    size_t Ntot = 1e6;
    double n = 1e-18, mAr = 6.6335209e-26, T = 500;
    NeutralGas gas(n, mAr, T);
    array<double, 3> vel = gas.generate_velocity();
    cout << vel[0] << " " << vel[1] << " " << vel[2] << endl;
    cout << sqrt(1.38e-23*T/mAr) << endl;
    clock_t start = clock();
    vector<array<double, 3>> velocity;
    for (int i = 0; i < Ntot; i++)
        vel = gas.generate_velocity();
        velocity.push_back(vel);
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/Ntot << endl;
    cout << velocity[0][0] << endl;
}

#endif //CPP_RZ_PIC_TEST_NEUTRALGAS_H
