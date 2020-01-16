#ifndef CPP_RZ_PIC_TEST_SIMULATION_H
#define CPP_RZ_PIC_TEST_SIMULATION_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../Field/PoissonSolver.h"
#include "../InteractionWithMaterial/ParticleLeave.h"
#define K_B 1.380649e-23

void test_Simulation() {
    cout << "test_Simulation: ";

    // Grid init
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-4, dr = 2e-4;
    Grid grid(Nz, Nr, dz, dr);

    // Overal particles information
    int seed;
    int Ntot = 5e5;
    scalar init_dens = 1e14;
    scalar volume = Nz * dz * M_PI * Nr * dr * Nr * dr;
    scalar ptcls_per_macro = init_dens * volume / Ntot;
    cout << "ptcls_per_macro: " << ptcls_per_macro << endl;

    // Electron Init
    Ntot = 5e5;
    scalar m_e = 9.1e-31;
    Particles electrons(9.1e-31, -1*1.6e-19, Ntot, grid, ptcls_per_macro);
    electrons.generate_velocities(1*1.6e-19, seed);
    array<scalar, 2> z_bounds = {2*dz, (Nz-2)*dz};
    array<scalar, 2> r_bounds = {2*dr, (Nr-2)*dr};
    electrons.generate_positions(z_bounds, r_bounds, seed);
    electrons.set_const_magnetic_field(0.1, 0);

    // Ion init
    seed = 2;
    Ntot = 5e5;
    scalar m_ion = 1e4*m_e;
    Particles ions(m_ion, 1.6e-19, Ntot, grid, ptcls_per_macro);
    scalar T_ion = 600;
    ions.generate_velocities((3 / 2) * K_B * T_ion, 2);
    ions.generate_positions(z_bounds, r_bounds, 2);
    ions.set_const_magnetic_field(0.1, 0);

    // Neutral gas init
    scalar n = 1e20, m_gas = m_ion, T_gas = 500;
    NeutralGas gas(n, m_gas, T_gas);

    // Collisions init
    scalar dt_collision = 1e-10;
    ElectronNeutralElasticCollision electron_elastic(1e-19, dt_collision, gas, electrons);
    scalar ion_threshold = 10*1.6e-19;
    Ionization argon_ionization(1e-20, ion_threshold, dt_collision, gas, electrons, ions);
    IonNeutralElasticCollision ion_elastic(1e-19, dt_collision, gas, ions);

    // Particle leave
    Matrix domain_condition(Nz, Nr);
    for (int row = 0; row < domain_condition.rows(); row++) {
        for (int col = 0; col < domain_condition.columns(); col++) {
            if ((row < 2 or row > Nz - 3) or (col > Nr - 3))
                domain_condition(row, col) = 1;
        }
    }
    ParticleLeave electrons_leave(electrons, grid, domain_condition);
    ParticleLeave ions_leave(ions, grid, domain_condition);

    // Electron emission
    Matrix domain_condition_left(Nz, Nr);
    for (int row = 0; row < domain_condition_left.rows(); row++) {
        for (int col = 0; col < domain_condition_left.columns(); col++) {
            if (row == 2 and col < Nr/2)
                domain_condition_left(row, col) = 1;
        }
    }
    Matrix domain_condition_right(Nz, Nr);
    for (int row = 0; row < domain_condition_right.rows(); row++) {
        for (int col = 0; col < domain_condition_right.columns(); col++) {
            if (row == Nz-3 and col < Nr/2)
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

    // Init for sum of electron and ion rho
    Matrix rho(Nz, Nr);

    // Field Init
    Matrix phi(Nz, Nr), radii(Nz, Nr);
    scalar CathodeV = -100., AnodeV = 0., tolerance=1e-4, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi, CathodeR, CathodeV, AnodeV, -50);
    Matrix Ez(Nz, Nr), Er(Nz, Nr);

    // PIC cycle
    int it_num = 100;
    scalar dt = 1e-12;
    int collision_step = dt_collision / dt;
    int ion_step = 1;
    clock_t start = clock();
    electrons.vel_pusher(-0.5*dt);
    ions.vel_pusher(-0.5*dt);
    for (int it = 0; it < it_num; it++) {
        //cout << "iter: " << it << endl;
        cout << "Ntot ions/electrons: " << ions.get_Ntot() << " " << electrons.get_Ntot() << endl;
        electrons.charge_interpolation();
        if (it % ion_step == 0)
            ions.charge_interpolation();
        rho = electrons.rho + ions.rho;
        PoissonSolverSOR(phi, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
        //PoissonSolverJacobi(phi, rho, radii, dz, dr, CathodeR, tolerance, max_iter);
        compute_E(Ez, Er, phi, dz, dr);
        electrons.electric_field_interpolation(Ez, Er);
        electrons.pusher(dt);
        if (it % ion_step == 0)
            ions.electric_field_interpolation(Ez, Er);
            ions.pusher(dt);
        if (it % collision_step == 0) {
            NanbuElectronCollisionProcess(electron_elastic, argon_ionization, 2);
            NanbuIonCollisionProcess(ion_elastic, 1);
        }
        electrons_leave.leave();
        if (it % ion_step == 0) {
            ions_leave.leave();
            emission_left.emission();
            emission_right.emission();
        }
    }
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds/it_num << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_SIMULATION_H
