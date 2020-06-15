#ifndef CPP_RZ_PIC_TEST_COLLISIONENERGYCROSSSECTION_H
#define CPP_RZ_PIC_TEST_COLLISIONENERGYCROSSSECTION_H


#include "../ElementaryProcesses/Collision.h"
#include <iostream>
#define K_B 1.380649e-23
#define EV 1.6021766208e-19
#define E_M 9.10938356e-31


scalar compute_vel_module(array<scalar, 3> vel) {
    return sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
}


void test_CollisionEnergyCrossSection() {
    cout << "test_CollisionEnergyCrossSection:" << endl;

    // Grid init
    int fine = 2;
    size_t Nz = 100*fine, Nr = 50*fine;
    scalar dz = 2e-4/fine, dr = 2e-4/fine;
    Grid grid(Nz, Nr, dz, dr);

    // Overal particles information
    int seed = 1;
    int ptcls_per_cell = 10;
    scalar init_dens = 1e14;

    int inject_region = 10*fine;
    scalar z_init_min = inject_region * dz, z_init_max = (Nz - inject_region) * dz;
    scalar r_init_min = inject_region * dr, r_init_max = (Nr - inject_region) * dr;
    scalar volume = (z_init_max - z_init_min) * M_PI * (r_init_max*r_init_max - r_init_min*r_init_min);
    int Ntot = (Nz-2*inject_region) * (Nr-2*inject_region) * ptcls_per_cell;
    scalar ptcls_per_macro = init_dens * volume / Ntot;
    cout << "ptcls_per_macro: " << ptcls_per_macro << endl;
    cout << "num_of_macro_ptcls: " << Ntot << endl;

    // Electron Init
    scalar m_e = E_M;
    Particles electron(m_e, -1*EV, Ntot, grid, ptcls_per_macro);
    electron.generate_velocities(50*EV, seed);
    array<scalar, 2> z_bounds = {z_init_min, z_init_max};
    array<scalar, 2> r_bounds = {r_init_min, r_init_max};
    electron.generate_positions(z_bounds, r_bounds, seed);
    electron.set_const_magnetic_field(0.1, 0);

    // Helium init
    scalar m_He = 7296.2965674 * E_M;
    Particles helium(m_He, EV, Ntot, grid, ptcls_per_macro);
    //scalar T_ion = 500;
    scalar T_ion = 50*11604;
    helium.generate_velocities((3. / 2) * K_B * T_ion, seed);
    helium.generate_positions(z_bounds, r_bounds, seed);
    helium.set_const_magnetic_field(0.1, 0);

    // Neutral gas init
    scalar n = 1e21, m_gas = m_He, T_gas = 500;
    NeutralGas gas(n, m_gas, T_gas);

    scalar dt_electron_collision = 5e-11;
    EnergyCrossSection electron_He_elastic_sigma("../ElementaryProcesses/CrossSectionData/e-He_elastic.txt");
    ElectronNeutralElasticCollision electron_elastic(electron_He_elastic_sigma, dt_electron_collision, gas, electron);
    EnergyCrossSection electron_He_ionization_sigma("../ElementaryProcesses/CrossSectionData/e-He_ionization.txt");
    Ionization electron_He_ionization(electron_He_ionization_sigma, dt_electron_collision, gas, electron, helium);

    cout << "electron_vel_old: " << compute_vel_module(electron.get_velocity(0)) << endl;
    electron_elastic.collision(0);
    cout << "electron_vel_new: " << compute_vel_module(electron.get_velocity(0)) << endl;
    cout << "electron_elastic_probability: " << electron_elastic.probability(0) << endl;

    cout << "--------------------------------" << endl;

    int ptcl_idx = 1;
    cout << "electron_Ntot before ionization: " << electron.get_Ntot() << endl;
    cout << "helium_Ntot before ionization: " << helium.get_Ntot() << endl;
    scalar vel = compute_vel_module(electron.get_velocity(ptcl_idx));
    cout << "electron_energy_before_ionization (EV): " << 0.5 * m_e * vel * vel / EV << endl;
    cout << "ionization_probability: " << electron_He_ionization.probability(ptcl_idx) << endl;

    electron_He_ionization.collision(ptcl_idx);

    cout << "electron_Ntot after ionization: " << electron.get_Ntot() << endl;
    cout << "helium_Ntot after ionization: " << helium.get_Ntot() << endl;
    vel = compute_vel_module(electron.get_velocity(ptcl_idx));
    cout << "electron_energy_after_ionization (EV): " << 0.5 * m_e * vel * vel / EV << endl;
    vel = compute_vel_module(electron.get_velocity(Ntot));
    cout << "new_electron_energy_after_ionization (EV): " << 0.5 * m_e * vel * vel / EV << endl;
    vel = compute_vel_module(helium.get_velocity(Ntot));
    cout << "new_helium_energy_after_ionization (EV): " << 0.5 * m_He * vel * vel / EV << endl;

    cout << "OK" << endl;
}


#endif //CPP_RZ_PIC_TEST_COLLISIONENERGYCROSSSECTION_H
