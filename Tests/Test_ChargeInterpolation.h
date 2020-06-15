#ifndef CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H
#define CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"

void test_Charge_interpolation() {
    // Particles Init
    cout << "test_Charge_interpolation: ";
    int seed = 3;
    int Ntot = 1e6;
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-4, dr = 2e-4;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, Ntot, grid);
    ptcls.generate_velocities(1, seed);
    array<scalar, 2> z_bounds = {dz, (Nz-2)*dz};
    array<scalar, 2> r_bounds = {dr*0., (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);
    ptcls.mfz.resize(Ntot, 0.1);

    ptcls.rho.print();
    cout << endl;
    ptcls.charge_interpolation();
    ptcls.rho.print();
    cout << ptcls.rho(10, 10) << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H
