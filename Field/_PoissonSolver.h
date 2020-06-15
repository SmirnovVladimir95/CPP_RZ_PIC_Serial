#ifndef CPP_RZ_PIC__POISSONSOLVER_H
#define CPP_RZ_PIC__POISSONSOLVER_H

#include "../Tools/ProjectTypes.h"

void _compute_b(scalar b[], const scalar rho[], int Nz, int Nr);
bool _convergence(const scalar phi[], const scalar b[], const scalar r[], int Nz, int Nr, scalar dz, scalar dr,
                  scalar tolerance);
void _PoissonSolverSOR(scalar phi[], const scalar rho[], const scalar radii[], int Nz, int Nr, scalar dz, scalar dr,
                       int CathodeR, scalar tolerance, int max_iter, scalar betta, int convergence_check);
void _PoissonSolverSOR_Improved(scalar phi[], const scalar rho[], const scalar radii[], int Nz, int Nr, scalar dz,
                                scalar dr, int CathodeR, scalar tolerance, int max_iter, scalar betta,
                                int convergence_check);
void _compute_E(scalar Ez[], scalar Er[], const scalar phi[], int Nz, int Nr, scalar dz, scalar dr);
bool _convergence_improved(const scalar phi[], const scalar b[], const scalar r[], const int Nz, const int Nr, const scalar dz,
                           const scalar dr, const scalar tolerance);

#endif //CPP_RZ_PIC__POISSONSOLVER_H
