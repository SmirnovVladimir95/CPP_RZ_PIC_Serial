#ifndef CPP_RZ_PIC_TEST_SIMULATION_H
#define CPP_RZ_PIC_TEST_SIMULATION_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../Field/PoissonSolver.h"
#include "../InteractionWithMaterial/ParticleLeave.h"
#include "../Tools/Logger.h"
#include "../Tools/ParticlesLogger.h"
#include <numeric>
#include <algorithm>
#define K_B 1.380649e-23
#define EV 1.6021766208e-19
#define E_M 9.10938356e-31

void test_Simulation() {
    cout << "test_Simulation: " << endl;

    // Grid init
    int fine = 2;
    size_t Nz = 100*fine, Nr = 50*fine;
    scalar dz = 2e-4/fine, dr = 2e-4/fine;
    Grid grid(Nz, Nr, dz, dr);

    // Overal particles information
    int seed;
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
    seed = 1;
    scalar m_e = E_M;
    Particles electrons(m_e, -1*EV, Ntot, grid, ptcls_per_macro);
    electrons.generate_velocities(1*EV, seed);
    scalar v_max = sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]);
    scalar v_min = sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]);
    array<scalar, 2> z_bounds = {z_init_min, z_init_max};
    array<scalar, 2> r_bounds = {r_init_min, r_init_max};
    electrons.generate_positions(z_bounds, r_bounds, seed);
    electrons.set_const_magnetic_field(0.1, 0);

    // Ion init
    seed = 3;
    scalar m_ion = 1e2*m_e;
    Particles ions(m_ion, EV, Ntot, grid, ptcls_per_macro);
    //scalar T_ion = 500;
    scalar T_ion = 1*11604;
    ions.generate_velocities((3 / 2) * K_B * T_ion, 2);
    ions.generate_positions(z_bounds, r_bounds, 2);
    ions.set_const_magnetic_field(0.1, 0);

    // Neutral gas init
    scalar n = 1e21, m_gas = m_ion, T_gas = 500;
    NeutralGas gas(n, m_gas, T_gas);

    // Collisions init
    scalar dt_collision_electron = 5e-10, dt_collision_ion = 5e-10;
    //ElectronNeutralElasticCollision electron_elastic(1e-19, dt_collision_electron, gas, electrons);
    IonNeutralElasticCollision electron_elastic(1e-19, dt_collision_electron, gas, electrons, false);
    scalar ion_threshold = 10*EV;
    Ionization argon_ionization(1e-20, ion_threshold, dt_collision_electron, gas, electrons, ions);
    IonNeutralElasticCollision ion_elastic(1e-19, dt_collision_ion, gas, ions);

    // Particle leave
    Matrix domain_condition_electron(Nz, Nr);
    for (int row = 0; row < domain_condition_electron.rows(); row++) {
        for (int col = 0; col < domain_condition_electron.columns(); col++) {
            if ((row == 0 or row == Nz - 1) or (col == Nr - 1))
                domain_condition_electron(row, col) = 1;
            if ((row == 1 or row == Nz - 2) and col >= Nr / 2)
                domain_condition_electron(row, col) = 1;
        }
    }
    Matrix domain_condition_ion(Nz, Nr);
    for (int row = 0; row < domain_condition_ion.rows(); row++) {
        for (int col = 0; col < domain_condition_ion.columns(); col++) {
            if ((row == 0 or row == Nz - 1) or (col == Nr - 1))
                domain_condition_ion(row, col) = 1;
            if ((row == 1 or row == Nz - 2) and col >= Nr / 2)
                domain_condition_ion(row, col) = 1;
        }
    }
    ParticleLeave electrons_leave(electrons, grid, domain_condition_electron);
    ParticleLeave ions_leave(ions, grid, domain_condition_ion);

    // Electron emission
    Matrix domain_condition_left(Nz, Nr);
    for (int row = 0; row < domain_condition_left.rows(); row++) {
        for (int col = 0; col < domain_condition_left.columns(); col++) {
            if (row == 1 and col < Nr / 2)
                domain_condition_left(row, col) = 1;
        }
    }
    Matrix domain_condition_right(Nz, Nr);
    for (int row = 0; row < domain_condition_right.rows(); row++) {
        for (int col = 0; col < domain_condition_right.columns(); col++) {
            if (row == Nz - 2 and col < Nr / 2)
                domain_condition_right(row, col) = 1;
        }
    }
    scalar gamma = 0.1, emission_energy = 10*EV;
    array<scalar, 3> emission_direction_left = {1, 0, 0};
    ParticleEmission emission_left(ions, electrons, grid, domain_condition_left, emission_direction_left,
                                   gamma, emission_energy);
    array<scalar, 3> emission_direction_right = {-1, 0, 0};
    ParticleEmission emission_right(ions, electrons, grid, domain_condition_right, emission_direction_right,
                                    gamma, emission_energy);

    cout << "domain_condition:"<<endl;
    Matrix overal_domain_condition(Nz, Nr);
    overal_domain_condition = domain_condition_ion + domain_condition_left + domain_condition_right;
    //overal_domain_condition.print();

    // Init for sum of electron and ion rho
    Matrix rho(Nz, Nr);

    // Field Init
    Matrix phi(Nz, Nr), radii(Nz, Nr);
    scalar CathodeV = -100., AnodeV = 0., tolerance=1e-4, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi, CathodeR, CathodeV, AnodeV, -50);
    Matrix Ez(Nz, Nr), Er(Nz, Nr);

    // PIC cycle parameters
    int it_num = 5e3*2;
    scalar dt = 5e-12/2;
    int collision_step_electron = dt_collision_electron / dt;
    int collision_step_ion = dt_collision_ion / dt;
    cout << "collision_step_electron: " << collision_step_electron << endl;
    cout << "collision_step_ion: " << collision_step_ion << endl;
    int ion_step = sqrt(m_ion/m_e);
    //int ion_step = 1;
    cout << "ion step: " << ion_step << endl;

    // Logger
    string phi_file = "phi.txt";
    int average_phi_it_num = 100;
    Matrix phi_average(Nz, Nr);
    clear_file(phi_file);

    string emission_file = "emission.txt";
    clear_file(emission_file);

    string ionization_file = "ionization.txt";
    clear_file(ionization_file);

    string electron_leave_file = "electron_leave_file.txt";
    clear_file(electron_leave_file);

    string ion_leave_file = "ion_leave_file.txt";
    clear_file(ion_leave_file);

    string Ntot_file = "Ntot_(iter).txt";
    int Ntot_step = 200*2, pos_step = it_num - 1, vel_step = it_num - 1, energy_step = 100;
    ParticlesLogger electrons_logger(electrons, "electrons");
    ParticlesLogger ions_logger(ions, "ions");

    ParticlesLogger electrons_traj(electrons, "electrons_traj");
    ParticlesLogger ions_traj(ions, "ions_traj");
    int pos_traj_step_electron = 10, pos_traj_step_ion = 100, vel_traj_step = 100;
    vector<int> track_ptcls = {0, 1};

    // PIC cycle
    clock_t start = clock();
    electrons.vel_pusher(-0.5*dt);
    ions.vel_pusher(-0.5*dt*ion_step);
    for (int it = 0; it < it_num; it++) {
        electrons.charge_interpolation();
        if (it % ion_step == 0) {
            ions.charge_interpolation();
        }
        rho = electrons.rho + ions.rho;
        PoissonSolverSOR(phi, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
        //PoissonSolverJacobi(phi, rho, radii, dz, dr, CathodeR, tolerance, max_iter);
        compute_E(Ez, Er, phi, dz, dr);
        electrons.electric_field_interpolation(Ez, Er);
        electrons.pusher(dt);
        if (it % ion_step == 0) {
            ions.electric_field_interpolation(Ez, Er);
            ions.pusher(dt*ion_step);
        }
        if (it % collision_step_electron == 0) {
            //NanbuElectronCollisionProcess(electron_elastic, argon_ionization, 2);
            //int electron_n_total_before_ionization = electrons.get_Ntot();
            //int ion_n_total_before_ionization = ions.get_Ntot();
            NanbuIonCollisionProcess(electron_elastic, argon_ionization, 2);
            //element_logging(it, ionization_file, " ");
            //element_logging(electrons.get_Ntot() - electron_n_total_before_ionization, ionization_file, " ");
            //element_logging(ions.get_Ntot() - ion_n_total_before_ionization, ionization_file);
        }
        if (it % collision_step_ion == 0) {
            NanbuIonCollisionProcess(ion_elastic, 1);
        }
        electrons_leave.leave();
        if (it % ion_step == 0) {
            ions_leave.leave();
            /*
            int electron_n_total_before = electrons.get_Ntot();
            int ion_n_total_before = ions.get_Ntot();
            */
            emission_left.emission();
            emission_right.emission();
            /*
            int electron_n_total_after = electrons.get_Ntot();
            int ion_n_total_after = ions.get_Ntot();
            element_logging(it, emission_file, " ");
            element_logging(electron_n_total_after - electron_n_total_before, emission_file, " ");
            element_logging(ion_n_total_after - ion_n_total_before, emission_file);
            */
        }
        /*
        if (it > it_num - average_phi_it_num) {
            phi_average += phi;
            if (it == it_num - 1) {
                phi_average = phi_average / average_phi_it_num;
                element_logging(phi, phi_file);
            }
        }
        */
        electrons_logger.n_total_log(it, Ntot_step);
        ions_logger.n_total_log(it, Ntot_step);
        electrons_logger.mean_energy_log(it, energy_step);
        ions_logger.mean_energy_log(it, energy_step);
        /*
        electrons_logger.position_log(it, pos_step);
        ions_logger.position_log(it, pos_step);
        electrons_logger.velocity_log(it, vel_step);
        ions_logger.velocity_log(it, vel_step);
        */
        //electrons_traj.position_log(it, pos_traj_step_electron, track_ptcls);
        //ions_traj.position_log(it, pos_traj_step_ion, track_ptcls);
    }
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << "time: " << seconds/it_num << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_SIMULATION_H
