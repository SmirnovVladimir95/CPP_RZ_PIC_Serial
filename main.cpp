#include <iostream>
#include <cmath>
#include "Tools/Matrix.h"
#include "Field/PoissonSolver.h"
#include <ctime>
#include "Particles/Particles.h"
#include <array>
#include "Grid/Grid.h"
using namespace std;

void test_Pusher() {
    int seed = 0;
    size_t Nz = 100, Nr = 50;
    double dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(1, 1.6e-19, 1e6, grid);
    cout << ptcls.Ntot << endl;
    ptcls.generate_velocities(1, seed);
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
    Bz.assign(ptcls.Ntot, 10);
    Br.assign(ptcls.Ntot, 10);
    ptcls.set_const_magnetic_field(Bz, Br);


    double dt = 1e-13;
    int num = 100;
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

void test_Field_interpolation() {
    size_t Nz = 100, Nr = 50;
    double dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(1, 1.6e-19, 1e6, grid);
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
}

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

void test_PIC_Cycle() {
    // Field Init
    size_t Nz = 100;
    size_t Nr = 50;
    Matrix phi_Jacobi(Nz, Nr), phi_SOR(Nz, Nr), rho(Nz, Nr), radii(Nz, Nr);
    double dz = 2e-5, dr = 2e-5, CathodeV = -100., AnodeV = 0., tolerance=1e-5, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi_SOR, CathodeR, CathodeV, AnodeV, -50);
    PoissonSolverSOR(phi_SOR, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
    phi_SOR.print();
    rho.print();

    // Particles Init
    int seed = 0;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, 1e6, grid);
    cout << ptcls.Ntot << endl;
    ptcls.generate_velocities(1, seed);
    array<double, 2> z_bounds = {dz, (Nz-2)*dz};
    array<double, 2> r_bounds = {dz, (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds, seed);

    // PIC Cycle
    int it_num = 100;
    for (int it = 0; it < it_num; it++) {
        //ptcls.electric_field_interpolation(Ez, Er);
    }
}


int main() {
    //test_Field_interpolation();
    //test_Pusher();
    //test_PIC_Cycle();
    test_compute_E();
    return 0;
}