#ifndef CPP_RZ_PIC_TEST_SINGLEPTCLSIMULATION_H
#define CPP_RZ_PIC_TEST_SINGLEPTCLSIMULATION_H

#include "../Tools/Matrix.h"
#include "../Particles/Particles.h"
#include "../Field/PoissonSolver.h"
#include "../InteractionWithMaterial/ParticleLeave.h"
#include "../Tools/Logger.h"
#include "../Tools/ParticlesLogger.h"
#define K_B 1.380649e-23
#define EV 1.6021766208e-19
#define E_M 9.10938356e-31


void test_SinglePtclSimulation() {
    cout << "test_SinglePtclSimulation: " << endl;

    // Grid init
    int fine = 2;
    size_t Nz = 100*fine, Nr = 50*fine;
    scalar dz = 2e-4/fine, dr = 2e-4/fine;
    Grid grid(Nz, Nr, dz, dr);

    // Electron Init
    scalar electrons_per_macro = 1;
    int electron_Ntot = 1;
    scalar m_e = E_M;
    Particles electrons(m_e, -1*EV, electron_Ntot, grid, electrons_per_macro);
    scalar T_e = 1*11604;
    scalar vel_e = sqrt(2 * (3 / 2.) * K_B * T_e / m_e);
    array<scalar, 3> velocity_e = {vel_e, 0, 0};
    array<scalar, 2> position_e = {0.01, 0.002};
    electrons.set_position(0, position_e);
    electrons.set_velocity(0, velocity_e);
    electrons.set_const_magnetic_field(0.1, 0);

    // Ion init
    int ion_Ntot = 1;
    double ions_per_macro = 1;
    scalar m_ion = 1e2*m_e;
    Particles ions(m_ion, EV, ion_Ntot, grid, ions_per_macro);
    scalar T_ion = 10000*11604;
    scalar vel = sqrt(2 * (3 / 2.) * K_B * T_ion / m_ion);
    array<scalar, 3> velocity = {-1*vel, 0, 0};
    array<scalar, 2> position = {0.001, 0.002};
    ions.set_position(0, position);
    ions.set_velocity(0, velocity);
    ions.set_const_magnetic_field(0.1, 0);

    // Neutral gas init
    scalar n = 1e21, m_gas = m_ion, T_gas = 500;
    NeutralGas gas(n, m_gas, T_gas);

    // Collisions init
    scalar dt_collision_electron = 1e-9, dt_collision_ion = 1e-9;
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
        }
    }
    Matrix domain_condition_ion(Nz, Nr);
    for (int row = 0; row < domain_condition_ion.rows(); row++) {
        for (int col = 0; col < domain_condition_ion.columns(); col++) {
            if ((row == 0 or row == Nz - 1) or (col == Nr - 1))
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

    // Field Init
    Matrix phi(Nz, Nr), radii(Nz, Nr);
    scalar CathodeV = -100., AnodeV = 0., tolerance=1e-4, betta=1.93;
    int CathodeR = Nr/2, max_iter=1e6;
    set_radii(radii, Nz, Nr, dr);
    phi_init(phi, CathodeR, CathodeV, AnodeV, -50);
    Matrix Ez(Nz, Nr), Er(Nz, Nr);

    // PIC cycle parameters
    int it_num = 1e4*2;
    scalar dt = 4e-12/2;
    int collision_step_electron = dt_collision_electron / dt;
    int collision_step_ion = dt_collision_ion / dt;
    int ion_step = sqrt(m_ion/m_e);

    // Logger
    ParticlesLogger electrons_traj(electrons, "electrons_traj");
    ParticlesLogger ions_traj(ions, "ions_traj");
    int pos_traj_step_electron = 1, pos_traj_step_ion = 100, vel_traj_step = 100;
    vector<int> track_ptcls = {0, 1};

    // PIC cycle
    Matrix rho(Nz, Nr);
    electrons.vel_pusher(-0.5*dt);
    ions.vel_pusher(-0.5*dt*ion_step);
    for (int it = 0; it < it_num; it++) {
        if (it % 100 == 0) {
            cout << "electrons_Ntot: " << electrons.get_Ntot() << endl;
            cout << "ions_Ntot: " << ions.get_Ntot() << endl;
        }
        //electrons.charge_interpolation();
        if (it % ion_step == 0) {
            //ions.charge_interpolation();
        }
        rho = electrons.rho + ions.rho;
        PoissonSolverSOR(phi, rho, radii, dz, dr, CathodeR, tolerance, max_iter, betta);
        compute_E(Ez, Er, phi, dz, dr);

        electrons.electric_field_interpolation(Ez, Er);
        electrons.pusher(dt);
        if (it % ion_step == 0) {
            ions.electric_field_interpolation(Ez, Er);
            ions.pusher(dt*ion_step);
        }
        if (it % collision_step_electron == 0) {
            //NanbuIonCollisionProcess(electron_elastic, argon_ionization, 2);
        }
        if (it % collision_step_ion == 0) {
            //NanbuIonCollisionProcess(ion_elastic, 1);
        }
        electrons_leave.leave();
        if (it % ion_step == 0) {
            int n = ions.get_Ntot();
            array<scalar, 2> ion_pos = ions.get_position(0);
            emission_left.emission();
            emission_right.emission();
            if (n - ions.get_Ntot() != 0) {
                cout << "new electron pos: " << electrons.get_position(1)[0] << " " << electrons.get_position(1)[1] << endl;
                cout << "ion pos: " << ion_pos[0] << " " << ion_pos[1] << endl;
                cout << "electrons_Ntot: " << electrons.get_Ntot() << endl;
            }
            ions_leave.leave();
        }
        electrons_traj.position_log(it, pos_traj_step_electron, track_ptcls);
        ions_traj.position_log(it, pos_traj_step_ion, track_ptcls);
        ions_traj.mean_energy_log(it, pos_traj_step_ion);
    }


    cout << "OK" << endl;
}


#endif //CPP_RZ_PIC_TEST_SINGLEPTCLSIMULATION_H
