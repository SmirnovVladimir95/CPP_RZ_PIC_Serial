#ifndef CPP_RZ_PIC_TEST_ELEMENTARYPROCESSES_H
#define CPP_RZ_PIC_TEST_ELEMENTARYPROCESSES_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../ElementaryProcesses/Collision.h"
#define K_B 1.380649e-23

void test_ElementaryProcesses() {
    cout << "test_ElementaryProcesses:" << endl;

    // Grid init
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-4, dr = 2e-4;
    Grid grid(Nz, Nr, dz, dr);

    // Overal particles information
    int seed;
    int ptcls_per_cell = 50;
    int Ntot = Nz * Nr * ptcls_per_cell;
    scalar init_dens = 1e14;
    scalar z_init_min = 10 * dz, z_init_max = (Nz - 10) * dz;
    scalar r_init_min = 10 * dr, r_init_max = (Nr - 10) * dr;
    scalar volume = Nz * dz * M_PI * Nr * dr * Nr * dr;
    scalar ptcls_per_macro = init_dens * volume / Ntot;
    cout << "ptcls_per_macro: " << ptcls_per_macro << endl;
    cout << "num_of_macro_ptcls: " << Ntot << endl;

    // Electron Init
    seed = 1;
    scalar m_e = 9.1e-31;
    Particles electrons(9.1e-31, -1*1.6e-19, Ntot, grid, ptcls_per_macro);
    electrons.generate_velocities(1*1.6e-19, seed);
    cout << "v_sigma: " << sqrt(2*1.6e-19/(3*9.1e-31)) << endl;
    scalar v_max = sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]);
    scalar v_min = sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]);
    scalar v_temp;
    for (int i = 0; i < electrons.get_Ntot(); i++) {
        v_temp = sqrt(electrons.vz[i]*electrons.vz[i] + electrons.vr[i]*electrons.vr[i] + electrons.vy[i]*electrons.vy[i]);
        if (v_temp > v_max)
            v_max = v_temp;
        if (v_temp < v_min)
            v_min = v_temp;
    }
    cout << "energy_max: " << m_e*v_max*v_max/2/1.6e-19 << " eV" << endl;
    cout << "energy_min: " << m_e*v_min*v_min/2/1.6e-19 << " eV" << endl;
    array<scalar, 2> z_bounds = {z_init_min, z_init_max};
    array<scalar, 2> r_bounds = {r_init_min, r_init_max};
    electrons.generate_positions(z_bounds, r_bounds, seed);
    electrons.set_const_magnetic_field(0.1, 0);

    // Ion init
    seed = 2;
    scalar m_ion = 1e2*m_e;
    Particles ions(m_ion, 1.6e-19, Ntot, grid, ptcls_per_macro);
    scalar T_ion = 600;
    ions.generate_velocities((3 / 2.) * K_B * T_ion, 2);
    ions.generate_positions(z_bounds, r_bounds, 2);
    ions.set_const_magnetic_field(0.1, 0);

    // Neutral gas init
    scalar n = 1e20, m_gas = 100*m_e, T_gas = 500;
    NeutralGas gas(n, m_gas, T_gas);

    // Collisions init
    scalar dt_collision = 1e-10;
    //ElectronNeutralElasticCollision electron_elastic(1e-19, dt_collision, gas, electrons);
    IonNeutralElasticCollision electron_elastic(1e-19, dt_collision, gas, electrons, false);
    scalar ion_threshold = 0*1.6e-19;
    Ionization argon_ionization(1e-20, ion_threshold, dt_collision, gas, electrons, ions);
    IonNeutralElasticCollision ion_elastic(1e-19, dt_collision, gas, ions);

    cout << "electron_elastic_prob: " << electron_elastic.probability(1) << endl;
    scalar summ1 = 0, summ2 = 0;
    for (int i = 0; i < electrons.get_Ntot(); i++) {
        summ1 += electron_elastic.probability(i);
        summ2 += argon_ionization.probability(i);
    }
    cout << "number of elastic collisions: " << summ1 << endl;
    cout << "number of ionizations: " << summ2 << endl;
}

#endif //CPP_RZ_PIC_TEST_ELEMENTARYPROCESSES_H
