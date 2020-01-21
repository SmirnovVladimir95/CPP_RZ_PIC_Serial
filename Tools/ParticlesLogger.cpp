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
    if (check_file(mean_energy_file))
        mean_energy_file = add + mean_energy_file;
    if (check_file(velocity_file))
        velocity_file = add + velocity_file;
    if (check_file(position_file))
        position_file = add + position_file;
    if (check_file(n_total_file))
        n_total_file = add + n_total_file;
}

void ParticlesLogger::n_total_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    element_logging(iter, n_total_file, " ", mode);
    element_logging(particles->get_Ntot(), n_total_file, "\n");
}

void ParticlesLogger::mean_energy_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    scalar average_energy = 0;
    scalar vel_2 = 0;
    int Ntot = particles->get_Ntot();
    for (int ptcl_idx = 0; ptcl_idx < Ntot; ptcl_idx++) {
        for (auto v : particles->get_velocity(ptcl_idx)) {
            vel_2 += v*v;
        }
        average_energy += particles->get_mass() * vel_2 / (2 * EV);
    }
    element_logging(iter, mean_energy_file, " ", mode);
    element_logging(average_energy, mean_energy_file, "\n");
}

void ParticlesLogger::velocity_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    element_logging(iter, velocity_file, "\n", mode);
    for (int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
        for (auto v: particles->get_velocity(ptcl_idx)) {
            element_logging(v, velocity_file, " ");
        }
        element_logging("", velocity_file, "\n");
    }
}

void ParticlesLogger::position_log(int iter, int step, unsigned int mode) {
    if (iter_check(iter, step))
        return;
    element_logging(iter, position_file, "\n", mode);
    for (int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
        for (auto pos: particles->get_position(ptcl_idx)) {
            element_logging(pos, position_file, " ");
        }
        element_logging("", position_file, "\n");
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
