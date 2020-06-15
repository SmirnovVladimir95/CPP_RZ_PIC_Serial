#include "_PoissonSolver.h"
#include <cmath>
#include <iostream>
#define EPSILON_0 8.854187817620389e-12


void _compute_b(scalar b[], const scalar rho[], const int Nz, const int Nr) {
    //#pragma omp parallel for num_threads(100)
    for (int i = 0; i < Nz*Nr; i++) {
        b[i] = rho[i] / EPSILON_0;
    }
}

bool _convergence(const scalar phi[], const scalar b[], const scalar r[], const int Nz, const int Nr, const scalar dz,
                  const scalar dr, const scalar tolerance) {
    scalar res;
    scalar dz2 = dz*dz;
    scalar dr2 = dr*dr;
    //#pragma omp parallel for private(res) num_threads(100)
    for (int i = 1; i < Nz-1; i++) {
        for (int j = 1; j < Nr-1; j++) {
            res = phi[i*Nr+j] - (b[i*Nr+j] + (phi[i*Nr+j-1] + phi[i*Nr+j+1])/dr2 +
                                 (phi[i*Nr+j+1] - phi[i*Nr+j-1])/(2*dr*r[i*Nr+j]) +
                                 (phi[(i-1)*Nr+j] + phi[(i+1)*Nr+j])/dz2) / (2/dr2 + 2/dz2);
            //#pragma omp single
            if (abs(res) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

void _PoissonSolverSOR(scalar phi[], const scalar rho[], const scalar radii[], const int Nz, const int Nr,
                       const scalar dz, const scalar dr, const int CathodeR, const scalar tolerance, const int max_iter,
                       const scalar betta, const int convergence_check) {
    auto* g = new scalar[Nz*Nr];
    auto* b = new scalar[Nz*Nr];
    scalar dz2 = dz*dz;
    scalar dr2 = dr*dr;
    _compute_b(b, rho, Nz, Nr);
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            g[i*Nr+j] = phi[i*Nr+j];
        }
    }
    for (int it = 0; it < max_iter; it++) {
        for (int i = 1; i < Nz-1; i++) {
            for (int j = 1; j < Nr-1; j++) {
                g[i*Nr+j] = betta*((b[i*Nr+j] +
                            (phi[i*Nr+j+1] + g[i*Nr+j-1]) / dr2 +
                            (phi[(i+1)*Nr+j] + g[(i-1)*Nr+j]) / dz2 +
                            (phi[i*Nr+j+1] - g[i*Nr+j-1])/(radii[i*Nr+j]*2*dr)) / (2/dr2 + 2/dz2)) +
                            (1-betta)*phi[i*Nr+j];
            }
        }
        // Neumann boundaries on the axis r = 0
        for(int i = 0; i < Nz; i++) {
            g[i*Nr] = g[i*Nr+1];
        }
        // Neumann boundaries on the left/right wall
        for (int i = CathodeR; i < Nr; i++) {
            g[i] = g[Nr+i];
            g[(Nz-1)*Nr+i] = g[(Nz-2)*Nr+i];
        }
        // Dirichlet boundaries
        for (int i = 0; i < Nz; i++) {
            for (int j = 0; j < Nr; j++) {
                phi[i*Nr+j] = g[i*Nr+j];
            }
        }
        // Convergence check
        if (it % convergence_check == 0 and _convergence(phi, b, radii, Nz, Nr, dz, dr, tolerance)) {
            //std::cout << "Convergence achieved at iteration: " << it << std::endl;
            break;
        }
    }

    delete [] g;
    delete [] b;
}

void _PoissonSolverSOR_Improved(scalar phi[], const scalar rho[], const scalar radii[], int Nz, int Nr, scalar dz,
                                scalar dr, int CathodeR, scalar tolerance, int max_iter, scalar betta,
                                int convergence_check) {
    auto* g = new scalar[Nz*Nr];
    auto* b = new scalar[Nz*Nr];
    scalar dz2 = dz*dz;
    scalar dr2 = dr*dr;
    _compute_b(b, rho, Nz, Nr);
    for (int i = 0; i < Nz*Nr; i++) {
        g[i] = phi[i];
    }
    for (int it = 0; it < max_iter; it++) {
        //#pragma omp parallel num_threads(100)
        //{
        // Red Ordering
        //#pragma omp for
        for (int i = 1; i < Nz-1; i++) {
            for (int j = 1; j < Nr-1; j++) {
                if ((i*Nr+j) % 2 == 0) {
                    g[i * Nr + j] = betta * ((b[i * Nr + j] +
                                              (phi[i * Nr + j + 1] + g[i * Nr + j - 1]) / dr2 +
                                              (phi[(i + 1) * Nr + j] + g[(i - 1) * Nr + j]) / dz2 +
                                              (phi[i * Nr + j + 1] - g[i * Nr + j - 1]) /
                                              (radii[i * Nr + j] * 2 * dr)) / (2 / dr2 + 2 / dz2)) +
                                    (1 - betta) * phi[i * Nr + j];
                }
            }
        }
        // Black Ordering
        //#pragma omp for
        for (int i = 1; i < Nz-1; i++) {
            for (int j = 1; j < Nr-1; j++) {
                if ((i*Nr+j) % 2 == 1) {
                    g[i * Nr + j] = betta * ((b[i * Nr + j] +
                                              (phi[i * Nr + j + 1] + g[i * Nr + j - 1]) / dr2 +
                                              (phi[(i + 1) * Nr + j] + g[(i - 1) * Nr + j]) / dz2 +
                                              (phi[i * Nr + j + 1] - g[i * Nr + j - 1]) /
                                              (radii[i * Nr + j] * 2 * dr)) / (2 / dr2 + 2 / dz2)) +
                                    (1 - betta) * phi[i * Nr + j];
                }
            }
        }
        // Neumann boundaries on the axis r = 0
        //#pragma omp for
        for(int i = 0; i < Nz; i++) {
            g[i*Nr] = g[i*Nr+1];
        }
        // Neumann boundaries on the left/right wall
        //#pragma omp for
        for (int i = CathodeR; i < Nr; i++) {
            g[i] = g[Nr+i];
            g[(Nz-1)*Nr+i] = g[(Nz-2)*Nr+i];
        }
        // Dirichlet boundaries
        //#pragma omp for
        for (int i = 0; i < Nz*Nr; i++) {
            phi[i] = g[i];
        }
        //}
        // Convergence check
        if (it % convergence_check == 0 and _convergence_improved(phi, b, radii, Nz, Nr, dz, dr, tolerance)) {
            //std::cout << "Convergence achieved at iteration: " << it << std::endl;
            break;
        }
    }

    delete [] g;
    delete [] b;
}

void _compute_E(scalar Ez[], scalar Er[], const scalar phi[], const int Nz, const int Nr, const scalar dz,
                const scalar dr) {
    scalar dz2 = dz*2;
    scalar dr2 = dr*2;
    // central difference, not right on walls
    for (int i = 0; i < Nz; i++) {
        for (int j = 0; j < Nr; j++) {
            if (i != 0 and i != Nz-1)
                Ez[i*Nr+j] = (phi[(i-1)*Nr+j] - phi[(i+1)*Nr+j]) / dz2;
            else if (i == 0)
                Ez[i*Nr+j] = (phi[i*Nr+j] - phi[(i+1)*Nr+j]) / dz;
            else
                Ez[i*Nr+j] = (phi[(i-1)*Nr+j] - phi[i*Nr+j]) / dz;
            if (j != 0 and j != Nr-1)
                Er[i*Nr+j] = (phi[i*Nr+j-1] - phi[i*Nr+j+1]) / dr2;
            else if (j == 0)
                Er[i*Nr+j] = (phi[i*Nr+j] - phi[i*Nr+j+1]) / dr;
            else
                Er[i*Nr+j] = (phi[i*Nr+j-1] - phi[i*Nr+j]) / dr;
        }
    }
}

bool _convergence_improved(const scalar phi[], const scalar b[], const scalar r[], const int Nz, const int Nr, const scalar dz,
                           const scalar dr, const scalar tolerance) {
    scalar res;
    scalar dz2 = dz*dz;
    scalar dr2 = dr*dr;
    bool flag = true;
    //#pragma omp parallel for private(res) num_threads(100)
    for (int i = 1; i < Nz-1; i++) {
        for (int j = 1; j < Nr-1; j++) {
            res = phi[i*Nr+j] - (b[i*Nr+j] + (phi[i*Nr+j-1] + phi[i*Nr+j+1])/dr2 +
                                 (phi[i*Nr+j+1] - phi[i*Nr+j-1])/(2*dr*r[i*Nr+j]) +
                                 (phi[(i-1)*Nr+j] + phi[(i+1)*Nr+j])/dz2) / (2/dr2 + 2/dz2);
            //#pragma omp single
            if (abs(res) > tolerance) {
                flag = false;
            }
        }
    }
    return flag;
}