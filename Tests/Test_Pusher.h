#ifndef CPP_RZ_PIC_TEST_PUSHER_H
#define CPP_RZ_PIC_TEST_PUSHER_H

#include <iostream>
#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
using namespace std;

void test_Pusher() {
    cout << "test_Pusher: ";

    int seed = 0;
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, 1e6, grid);
    //cout << ptcls.get_Ntot() << endl;
    ptcls.generate_velocities(1*1.6e-19, seed);
    array<scalar, 2> z_bounds = {dz, (Nz-2)*dz};
    array<scalar, 2> r_bounds = {0, (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);
    //cout << *min_element(ptcls.r.data(), ptcls.r.data()+ptcls.r.size()) << endl;
    //cout << *min_element(ptcls.node_volume, ptcls.node_volume+Nr) << endl;

    Matrix Ez(Nz, Nr);
    Matrix Er(Nz, Nr);
    Ez.fill(-100);
    Er.fill(-100);
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<scalar> distribution_z(0, 100);
    std::uniform_real_distribution<scalar> distribution_r(0, 100);
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            Ez(i, j) = distribution_z(generator);
            Er(i, j) = distribution_r(generator);
        }
    }
    ptcls.electric_field_interpolation(Ez, Er);
    ptcls.charge_interpolation();
    //ptcls.rho.print();
    scalar Bz = 10, Br = 10;
    ptcls.set_const_magnetic_field(Bz, Br);

    scalar dt = 1e-11;
    int num = 10;
    clock_t start = clock();
    //scalar start = omp_get_wtime();
    ptcls.vel_pusher(-0.5*dt);
    for (int i = 0; i < num; i++) {
        ptcls.pusher(dt);
    }
    //scalar end = omp_get_wtime();
    //scalar seconds = end - start;
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/num << endl;
    cout << "OK" << endl;
}

void test_single_particle_motion() {
    cout << "test_single_particle_motion" << endl;

    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-5, dr = 2e-5;
    scalar dt = 1e-1;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(1, 1, 1, grid);
    ptcls.z[0] = ptcls.r[0] = ptcls.vz[0] = ptcls.vr[0] = ptcls.vy[0] =  1;
    ptcls.efz[0] = ptcls.efr[0] = ptcls.mfz[0] = 1;
    ptcls.mfr[0] = 0;
    ptcls.vel_pusher(-0.5*dt);
    ptcls.pusher(dt);
    cout << ptcls.z[0] << " " << ptcls.r[0] << " " << ptcls.vz[0] << " " << ptcls.vr[0] << " " << ptcls.vy[0] << endl;

    cout << "OK" << endl;
}


#endif //CPP_RZ_PIC_TEST_PUSHER_H
