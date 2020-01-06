#ifndef CPP_RZ_PIC_TEST_FIELDSOLVER_H
#define CPP_RZ_PIC_TEST_FIELDSOLVER_H

#include "../Tools/Matrix.h"
#include "../Field/PoissonSolver.h"

void test_Field_solver() {
    cout << "test_Field_solver:" << endl;
    size_t Nz = 100;
    size_t Nr = 50;
    Matrix phi_Jacobi(Nz, Nr), phi_SOR(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    scalar dz = 2e-4, dr = 2e-4, CathodeV = -100., AnodeV = 0., tolerance=1e-4, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);

    cout << "Jacobi method:" << endl;
    clock_t start = clock();
    //scalar start = omp_get_wtime();
    phi_init(phi_Jacobi, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverJacobi(phi_Jacobi, rho, radii, dz, dr, CathodeR, tolerance, max_iter);
    clock_t end = clock();
    //scalar end = omp_get_wtime();
    //scalar seconds = end - start;
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds << endl;

    cout << "SOR method:" << endl;
    start =  clock();
    //start = omp_get_wtime();
    phi_init(phi_SOR, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverSOR(phi_SOR, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
    end = clock();
    //end = omp_get_wtime();
    //seconds = end - start;
    seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds << endl;
    //phi_SOR.print();
    //phi_Jacobi.print();
    scalar max = 0, res;
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            res = phi_Jacobi(i, j) - phi_SOR(i, j);
            if (abs(res) > max) {
                max = abs(res);
            }
        }
    }
    cout << "Max difference between Jacobi and SOR methods: " << max << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_FIELDSOLVER_H
