#include "ParticlesLogger.h"
#include "Logger.h"
#define EV 1.60217662e-19


ParticlesLogger::ParticlesLogger(Particles& particles, const string& observerID) : particles(&particles),
                                                                                   observerID(observerID) {
    mean_energy_file = observerID + "_mean_energy.txt";
    velocity_file = observerID + "_velocities.txt";
    position_file = observerID + "_particles.txt";
    n_total_file = observerID + "_Ntot.txt";
    // file check
    string add = "new_";
    while (check_file(mean_energy_file))
        mean_energy_file = add + mean_energy_file;
    while (check_file(velocity_file))
        velocity_file = add + velocity_file;
    while (check_file(position_file))
        position_file = add + position_file;
    while (check_file(n_total_file))
        n_total_file = add + n_total_file;
    //open file
    mean_energy_output.open(mean_energy_file);
    velocity_output.open(velocity_file);
    position_output.open(position_file);
    n_total_output.open(n_total_file);
}

void ParticlesLogger::n_total_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    n_total_output << iter << " ";
    n_total_output << particles->get_Ntot() << endl;
}

void ParticlesLogger::mean_energy_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    scalar mean_energy = 0, vel_2;
    int Ntot = particles->get_Ntot();
    for (int ptcl_idx = 0; ptcl_idx < Ntot; ptcl_idx++) {
        vel_2 = 0;
        for (auto v : particles->get_velocity(ptcl_idx)) {
            vel_2 += v*v;
        }
        mean_energy += particles->get_mass() * vel_2 / (2 * EV);
    }
    mean_energy /= Ntot;
    mean_energy_output << iter << " ";
    mean_energy_output << mean_energy << endl;
}

void ParticlesLogger::velocity_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step) or iter == 0)
        return;
    velocity_output << iter << endl;
    for (int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
        for (auto v: particles->get_velocity(ptcl_idx)) {
            velocity_output << v << " ";
        }
        velocity_output << endl;
    }
}

void ParticlesLogger::position_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step) or iter == 0)
        return;
    position_output << iter << endl;
    for (int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
        for (auto pos: particles->get_position(ptcl_idx)) {
            position_output << pos << " ";
        }
        position_output << endl;
    }
}

void ParticlesLogger::log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    n_total_log(iter, step, mode);
    mean_energy_log(iter, step, mode);
    velocity_log(iter, step, mode);
    position_log(iter, step, mode);
}

bool ParticlesLogger::iter_check(int iter, int step) {
    return iter % step != 0;
}

ParticlesLogger::~ParticlesLogger() {
    mean_energy_output.close();
    velocity_output.close();
    position_output.close();
    n_total_output.close();
}
