#ifndef CPP_RZ_PIC_TEST_PUSHER_H
#define CPP_RZ_PIC_TEST_PUSHER_H

#include <iostream>
#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
using namespace std;

void test_Pusher() {
    int seed = 0;
    size_t Nz = 100, Nr = 50;
    double dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, 1e6, grid);
    cout << ptcls.get_Ntot() << endl;
    ptcls.generate_velocities(1*1.6e-19, seed);
    array<double, 2> z_bounds = {dz, (Nz-2)*dz};
    array<double, 2> r_bounds = {dz, (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);

    Matrix Ez(Nz, Nr);
    Matrix Er(Nz, Nr);
    Ez.fill(-100);
    Er.fill(-100);
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> distribution_z(0, 100);
    std::uniform_real_distribution<double> distribution_r(0, 100);
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            Ez(i, j) = distribution_z(generator);
            Er(i, j) = distribution_r(generator);
        }
    }
    ptcls.electric_field_interpolation(Ez, Er);
    ptcls.charge_interpolation();
    vector<double> Bz, Br;
    Bz.assign(ptcls.get_Ntot(), 10);
    Br.assign(ptcls.get_Ntot(), 10);
    ptcls.set_const_magnetic_field(Bz, Br);


    double dt = 1e-11;
    int num = 10;
    clock_t start = clock();
    //double start = omp_get_wtime();
    ptcls.vel_pusher(-0.5*dt);
    ptcls.vel_pusher(-0.5*dt);
    for (int i = 0; i < num; i++) {
        ptcls.pusher(dt);
    }
    //double end = omp_get_wtime();
    //double seconds = end - start;
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/num << endl;
}


#endif //CPP_RZ_PIC_TEST_PUSHER_H
