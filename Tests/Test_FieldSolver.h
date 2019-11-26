#ifndef CPP_RZ_PIC_TEST_FIELDSOLVER_H
#define CPP_RZ_PIC_TEST_FIELDSOLVER_H

#include "../Tools/Matrix.h"
#include "../Field/PoissonSolver.h"

void test_Field_solver() {
    size_t Nz = 100;
    size_t Nr = 50;
    Matrix phi_Jacobi(Nz, Nr), phi_SOR(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    double dz = 2e-5, dr = 2e-5, CathodeV = -100., AnodeV = 0., tolerance=1e-5, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);

    clock_t start = clock();
    //double start = omp_get_wtime();
    phi_init(phi_Jacobi, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverJacobi(phi_Jacobi, rho, radii, dz, dr, CathodeR, tolerance, max_iter);
    clock_t end = clock();
    //double end = omp_get_wtime();
    //double seconds = end - start;
    double seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout << seconds << endl;

    start =  clock();
    //start = omp_get_wtime();
    phi_init(phi_SOR, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverSOR(phi_SOR, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
    end = clock();
    //end = omp_get_wtime();
    //seconds = end - start;
    seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout << seconds << endl;

    double res;
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            res = phi_Jacobi(i, j) - phi_SOR(i, j);
            if (abs(res) > 1e-1) {
                cout << "wrong" << endl;
                return;
            }
        }
    }
}

#endif //CPP_RZ_PIC_TEST_FIELDSOLVER_H
