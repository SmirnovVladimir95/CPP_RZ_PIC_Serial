#ifndef CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H
#define CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"

void test_Charge_interpolation() {
    // Particles Init
    int seed = 2;
    int Ntot = 1e1;
    size_t Nz = 10, Nr = 5;
    double dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, Ntot, grid);
    cout << ptcls.Ntot << endl;
    ptcls.generate_velocities(1, seed);
    array<double, 2> z_bounds = {dr, (Nz-1)*dz};
    array<double, 2> r_bounds = {dr, (Nr-1)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);
    ptcls.mfz.resize(Ntot, 0.1);

    ptcls.rho.print();
    ptcls.charge_interpolation();
    ptcls.rho.print();
}

#endif //CPP_RZ_PIC_TEST_CHARGEINTERPOLATION_H
