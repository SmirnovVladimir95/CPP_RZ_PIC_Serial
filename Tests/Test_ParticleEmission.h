#ifndef CPP_RZ_PIC_TEST_PARTICLEEMISSION_H
#define CPP_RZ_PIC_TEST_PARTICLEEMISSION_H

#include "../Particles/Particles.h"
#include "../InteractionWithMaterial/ParticleEmission.h"


void test_ParticleEmission() {
    cout << "test_ParticleEmission:" << endl;

    int seed = 0;
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles electrons(9.1e-31, 1.6e-19, 1e6, grid);
    Particles ions(6.6e-26, 1.6e-19, 1e6, grid);
    array<scalar, 2> z_bounds = {dz, (Nz-2)*dz};
    array<scalar, 2> r_bounds = {dz, (Nr-2)*dr};
    electrons.generate_positions(z_bounds, r_bounds, seed);
    ions.generate_positions(z_bounds, r_bounds, seed);

    cout << "electrons/ions num before emission: " << electrons.get_Ntot() << "/" << ions.get_Ntot() << endl;

    Matrix domain_condition_left(Nz, Nr);
    for (int row = 0; row < domain_condition_left.rows(); row++) {
        for (int col = 0; col < domain_condition_left.columns(); col++) {
            if (row == 1 and col < Nr/2)
                domain_condition_left(row, col) = 1;
        }
    }
    Matrix domain_condition_right(Nz, Nr);
    for (int row = 0; row < domain_condition_right.rows(); row++) {
        for (int col = 0; col < domain_condition_right.columns(); col++) {
            if (row == Nz-2 and col < Nr/2)
                domain_condition_right(row, col) = 1;
        }
    }

    scalar gamma = 0.1, emission_energy = 10*1.6e-19;
    array<scalar, 3> emission_direction_left = {1, 0, 0};
    ParticleEmission emission_left(ions, electrons, grid, domain_condition_left, emission_direction_left,
                                   gamma, emission_energy);
    array<scalar, 3> emission_direction_right = {-1, 0, 0};
    ParticleEmission emission_right(ions, electrons, grid, domain_condition_right, emission_direction_right,
                                   gamma, emission_energy);

    clock_t start = clock();
    emission_left.emission();
    emission_right.emission();
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;

    cout << "time for emission: " << seconds << endl;
    cout << "electrons/ions num after emission: " << electrons.get_Ntot() << "/" << ions.get_Ntot() << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_PARTICLEEMISSION_H
