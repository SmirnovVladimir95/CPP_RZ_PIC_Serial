#ifndef CPP_RZ_PIC_TEST_PARTICLESLOAD_H
#define CPP_RZ_PIC_TEST_PARTICLESLOAD_H

#include "../Grid/Grid.h"
#include <iostream>
#include <cmath>
#include "../Particles/Particles.h"
#include "../Tools/ParticlesLoad.h"
using namespace std;


void test_ParticleLoad() {
    cout << "test_ParticleLoad:" << endl;

    // Grid init
    size_t Nz = 100, Nr = 50;
    scalar dz = 2e-4, dr = 2e-4;
    Grid grid(Nz, Nr, dz, dr);

    // Overal particles information
    int seed;
    int ptcls_per_cell = 50;
    scalar init_dens = 1e14;
    int inject_region = 10;
    scalar z_init_min = inject_region * dz, z_init_max = (Nz - inject_region) * dz;
    scalar r_init_min = inject_region * dr, r_init_max = (Nr - inject_region) * dr;
    scalar volume = (z_init_max - z_init_min) * M_PI * (r_init_max - r_init_min) * (r_init_max - r_init_min);
    //int Ntot =  (Nz-2*inject_region) * (Nr-2*inject_region) * ptcls_per_cell;
    int Ntot = 143348;
    scalar ptcls_per_macro = init_dens * volume / Ntot;
    cout << "ptcls_per_macro: " << ptcls_per_macro << endl;
    cout << "num_of_macro_ptcls: " << Ntot << endl;

    // Electron Init
    seed = 1;
    scalar m_e = 9.1e-31;
    Particles electrons(9.1e-31, -1*1.6e-19, Ntot, grid, ptcls_per_macro);

    //Load velocity and position from file
    string pos_path = "../Particles/Data/electron_positions.txt";
    string vel_path = "../Particles/Data/electron_velocities.txt";
    ParticlesLoad load(electrons, pos_path, vel_path);
    cout << "ptcls number: " << electrons.get_Ntot() << endl;
    cout << electrons.get_position(Ntot-1)[0] << " " << electrons.get_position(Ntot-1)[1] << endl;
    cout << electrons.get_velocity(0)[0] << " " << electrons.get_velocity(0)[1] << " " << electrons.get_velocity(0)[2] << endl;

    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_PARTICLESLOAD_H


