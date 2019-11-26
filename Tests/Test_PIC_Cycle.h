#ifndef CPP_RZ_PIC_TEST_PIC_CYCLE_H
#define CPP_RZ_PIC_TEST_PIC_CYCLE_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../Field/PoissonSolver.h"

void test_PIC_Cycle() {
    // Field Init
    size_t Nz = 100;
    size_t Nr = 50;
    Matrix phi_Jacobi(Nz, Nr), phi(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    double dz = 2e-5, dr = 2e-5, CathodeV = -100., AnodeV = 0., tolerance=1e-3, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi, CathodeR, CathodeV, AnodeV, -50);
    Matrix Ez(Nz, Nr), Er(Nz, Nr);

    // Particles Init
    int seed = 0, Ntot=1e6;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, Ntot, grid);
    cout << ptcls.Ntot << endl;
    ptcls.generate_velocities(1*1.6e-19, seed);
    array<double, 2> z_bounds = {dz, (Nz-2)*dz};
    array<double, 2> r_bounds = {dz, (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);
    ptcls.mfz.resize(Ntot, 0.1);

    // PIC Cycle
    int it_num = 10;
    double dt = 1e-12;
    clock_t start = clock();
    ptcls.vel_pusher(-0.5*dt);
    for (int it = 0; it < it_num; it++) {
        //cout << "iter: " << it << endl;
        ptcls.charge_interpolation();
        PoissonSolverSOR(phi, ptcls.rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
        //PoissonSolverJacobi(phi, ptcls.rho, radii, dz, dr, CathodeR, tolerance, max_iter);
        compute_E(Ez, Er, phi, dz, dr);
        ptcls.electric_field_interpolation(Ez, Er);
        ptcls.pusher(dt);
    }
    clock_t end = clock();
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/it_num << endl;
}

#endif //CPP_RZ_PIC_TEST_PIC_CYCLE_H
