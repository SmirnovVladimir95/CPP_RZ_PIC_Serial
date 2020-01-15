#ifndef CPP_RZ_PIC_POISSONSOLVER_H
#define CPP_RZ_PIC_POISSONSOLVER_H

#include "../Tools/Matrix.h"

void set_radii(Matrix& radii, int Nz, int Nr, scalar dr);
void compute_b(Matrix& b, Matrix& rho);
void phi_init(Matrix& phi, const int CathodeR, const scalar CathodeV, const scalar AnodeV, scalar value=0);
bool convergence(Matrix& phi, Matrix& b, Matrix& r, const scalar dz, const scalar dr, const scalar tolerance);
void PoissonSolverJacobi(Matrix& phi, Matrix& rho, Matrix& radii, const scalar dz, const scalar dr, const int CathodeR,
                         const scalar tolerance, const int max_iter=1e6, int convergence_check=10);
void PoissonSolverSOR(Matrix& phi, Matrix& rho, Matrix& radii, const scalar dz, const scalar dr, const int CathodeR,
                      const scalar tolerance, const int max_iter=1e6, const scalar betta=1.5, int convergence_check=1);
void compute_E(Matrix& Ez, Matrix& Er, Matrix& phi, const scalar dz, const scalar dr);

#endif //CPP_RZ_PIC_POISSONSOLVER_H
