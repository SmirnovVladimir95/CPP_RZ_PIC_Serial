#ifndef CPP_RZ_PIC_TEST_FIELDINTERPOLATION_H
#define CPP_RZ_PIC_TEST_FIELDINTERPOLATION_H

#include <iostream>
#include <ctime>
#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
using namespace std;

void test_Field_interpolation() {
    size_t Nz = 100, Nr = 50;
    double dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(1, 1.6e-19, 1e1, grid);
    Matrix Ez(Nz, Nr);
    Matrix Er(Nz, Nr);
    Ez.fill(-100);
    Er.fill(-100);
    double start, end, summ_time=0, num = 10;
    for (int i = 0; i < num; i++) {
        Ez.fill(-100);
        Er.fill(-100);
        //start = omp_get_wtime();
        clock_t start = clock();
        ptcls.electric_field_interpolation(Ez, Er);
        //end = omp_get_wtime();
        clock_t end = clock();
        //summ_time += end - start;
        summ_time += (double)(end - start) / CLOCKS_PER_SEC;
        Ez.fill(0);
        Er.fill(0);
    }
    cout << summ_time/num << endl;
    for(int i = 0; i < ptcls.Ntot; i++)
        cout << ptcls.efz[i] << " ";
}


#endif //CPP_RZ_PIC_TEST_FIELDINTERPOLATION_H
