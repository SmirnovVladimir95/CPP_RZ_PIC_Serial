#include "../Grid/Grid.h"
#include <cmath>
#include "Interpolation.h"
//#include <omp.h>
//#define NUM_THREADS 100

void LinearFieldInterpolation(scalar efz[], scalar efr[], const scalar z[], const scalar r[],
                              const scalar Ez[], const scalar Er[], const Grid& grid, const size_t Ntot) {
    int cell_z, cell_r, Nr=grid.Nr;
    scalar hz, hr;
    //#pragma omp parallel for private(hz, hr, cell_z, cell_r) num_threads(NUM_THREADS)
    for (int i = 0; i < Ntot; i++) {
        cell_z = floor(z[i]/grid.dz);
        cell_r = floor(r[i]/grid.dr);
        hz = (z[i] - cell_z*grid.dz) / grid.dz;
        hr = (r[i] - cell_r*grid.dr) / grid.dr;

        efz[i] = Ez[cell_z*Nr + cell_r] * (1 - hz) * (1 - hr);
        efz[i] += Ez[(cell_z+1)*Nr + cell_r] * hz * (1 - hr);
        efz[i] += Ez[(cell_z+1)*Nr + cell_r+1] * hz * hr;
        efz[i] += Ez[cell_z*Nr + cell_r+1] * (1 - hz) * hr;

        efr[i] = Er[cell_z*Nr + cell_r] * (1 - hz) * (1 - hr);
        efr[i] += Er[(cell_z+1)*Nr + cell_r] * hz * (1 - hr);
        efr[i] += Er[(cell_z+1)*Nr + cell_r+1] * hz * hr;
        efr[i] += Er[cell_z*Nr + cell_r+1] * (1 - hz) * hr;
    }
}

void LinearChargeInterpolation(scalar rho[], const scalar z[], const scalar r[], const Grid& grid,
                               const scalar charge, const size_t Ntot, const scalar node_volume[]) {
    int cell_z, cell_r, Nr = grid.Nr;
    scalar hz, hr;
    //#pragma omp parallel private(hz, hr, cell_z, cell_r) num_threads(NUM_THREADS)
    //{
    //#pragma omp for
    for (int i = 0; i < grid.Nz*grid.Nr; i++) {
        rho[i] = 0;
    }
    //#pragma omp for
    for (int i = 0; i < Ntot; i++) {
        cell_z = floor(z[i] / grid.dz);
        cell_r = floor(r[i] / grid.dr);
        hz = (z[i] - cell_z * grid.dz) / grid.dz;
        hr = (r[i] - cell_r * grid.dr) / grid.dr;
        rho[cell_z * Nr + cell_r] += charge * (1 - hz) * (1 - hr) / node_volume[cell_z * Nr + cell_r];
        rho[(cell_z + 1) * Nr + cell_r] += charge * hz * (1 - hr) / node_volume[(cell_z + 1) * Nr + cell_r];
        rho[(cell_z + 1) * Nr + cell_r + 1] += charge * hz * hr / node_volume[(cell_z + 1) * Nr + cell_r + 1];
        rho[cell_z * Nr + cell_r + 1] += charge * (1 - hz) * hr / node_volume[cell_z * Nr + cell_r + 1];
    }
    //}
}

void InitVolume(Matrix& node_volume, const Grid& grid) {
    scalar r_min, r_max, r_middle;
    for (int i = 0; i < grid.Nz; i++) {
        for (int j = 0; j < grid.Nr; j++) {
            if (j > 0 and j < grid.Nr - 1) {
                r_min = (j - 1) * grid.dr;
                r_middle = j * grid.dr;
                r_max = (j + 1) * grid.dr;
                node_volume(i, j) = grid.dz * (M_PI / 3) * (r_max * (r_middle + r_max) - r_min * (r_min + r_middle));
            } else if (j == 0) {
                r_min = j * grid.dr;
                r_middle = grid.dr;
                node_volume(i, j) = grid.dz * (M_PI / 3) * (r_middle - r_min) * (2 * r_min + r_middle);
            } else {
                r_middle = (j - 1) * grid.dr;
                r_max = j * grid.dr;
                node_volume(i, j) = grid.dz * (M_PI / 3) * (r_max - r_middle) * (r_middle + 2 * r_max);
            }
        }
    }
}