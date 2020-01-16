#ifndef CPP_RZ_PIC_TEST_PIC_CYCLE_H
#define CPP_RZ_PIC_TEST_PIC_CYCLE_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../Field/PoissonSolver.h"

void test_PIC_Cycle() {
    // Field Init
    cout << "test_PIC_Cycle: ";
    size_t Nz = 100;
    size_t Nr = 50;
    Matrix phi_Jacobi(Nz, Nr), phi(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    scalar dz = 2e-4, dr = 2e-4, CathodeV = -100., AnodeV = 0., tolerance=1e-4, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi, CathodeR, CathodeV, AnodeV, -50);
    Matrix Ez(Nz, Nr), Er(Nz, Nr);

    // Particles Init
    int seed = 1, Ntot=1e6;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, Ntot, grid);
    cout << ptcls.get_Ntot() << endl;
    ptcls.generate_velocities(1*1.6e-19, seed);
    array<scalar, 2> z_bounds = {dz, (Nz-1)*dz};
    array<scalar, 2> r_bounds = {dz, (Nr-1)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);
    ptcls.mfz.resize(Ntot, 0.1);

    //PoissonSolverSOR(phi, ptcls.rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
    //compute_E(Ez, Er, phi, dz, dr);
    // PIC Cycle
    int it_num = 100;
    scalar dt = 1e-12;
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
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/it_num << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_PIC_CYCLE_H
