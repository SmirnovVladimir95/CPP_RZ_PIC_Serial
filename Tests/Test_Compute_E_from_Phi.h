#ifndef CPP_RZ_PIC_TEST_COMPUTE_E_FROM_PHI_H
#define CPP_RZ_PIC_TEST_COMPUTE_E_FROM_PHI_H

#include "../Tools/Matrix.h"
#include "../Field/PoissonSolver.h"

void test_compute_E() {
    // Field Init
    size_t Nz = 10;
    size_t Nr = 5;
    Matrix phi_Jacobi(Nz, Nr), phi_SOR(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    double dz = 2e-2, dr = 2e-2, CathodeV = -100., AnodeV = 0., tolerance=1e-5, betta=1.5 ;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi_SOR, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverSOR(phi_SOR, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
    phi_SOR.print();

    Matrix Ez(Nz, Nr), Er(Nz, Nr);
    compute_E(Ez, Er, phi_SOR, dz, dr);
    Ez.print();
    Er.print();
}

#endif //CPP_RZ_PIC_TEST_COMPUTE_E_FROM_PHI_H
