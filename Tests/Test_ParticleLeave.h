#ifndef CPP_RZ_PIC_TEST_PARTICLELEAVE_H
#define CPP_RZ_PIC_TEST_PARTICLELEAVE_H

#include "../Particles/Particles.h"
#include "../InteractionWithMaterial/ParticleLeave.h"

void test_ParticleLeave() {
    cout << "test_ParticleLeave:" << endl;

    int seed = 0;
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-5, dr = 2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(9.1e-31, 1.6e-19, 1e6, grid);
    array<scalar, 2> z_bounds = {dz, (Nz-2)*dz};
    array<scalar, 2> r_bounds = {dz, (Nr-2)*dr};
    ptcls.generate_positions(z_bounds, r_bounds);

    cout << "ptcls num before elimination: " << ptcls.get_Ntot() << endl;
    cout << "rho before elimination:" << endl;
    ptcls.charge_interpolation();
    ptcls.rho.print();

    Matrix domain_condition(Nz, Nr);
    for (int row = 0; row < domain_condition.rows(); row++) {
        for (int col = 0; col < domain_condition.columns(); col++) {
            if (((row < 2 or row >= Nz - 3) or (col >= Nr - 3)))
                domain_condition(row, col) = 1;
        }
    }
    cout << "domain_condition:" << endl;
    domain_condition.print();

    ParticleLeave leave(ptcls, grid, domain_condition);
    clock_t start = clock();
    leave.leave();
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    ptcls.charge_interpolation();

    cout << "rho after elimination:" << endl;
    ptcls.rho.print();

    cout << "time for particle elimination: " << seconds << endl;
    cout << "ptcls num after elimination: " << ptcls.get_Ntot() << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_PARTICLELEAVE_H
