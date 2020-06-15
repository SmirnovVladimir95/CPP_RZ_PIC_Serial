#include "ParticlesLoad.h"
#include <fstream>
#include "../Tools/Helpers.h"


ParticlesLoad::ParticlesLoad(Particles &particles, const string& pos_path, const string& vel_path) : particles(&particles), pos_path(pos_path), vel_path(vel_path) {
    if (not pos_path.empty()) {
        position_load();
    } else {
        cout << "Warning no particle position filename is specified in constructor" << endl;
    }
    if (not vel_path.empty()) {
        velocity_load();
    } else {
        cout << "Warning no particle velocity filename is specified in constructor" << endl;
    }
}

void ParticlesLoad::velocity_load() {
    ifstream input(vel_path);
    vector<scalar> tmp;
    array<scalar, 3> vel{};
    string line;
    if (input) {
        for(int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
            for(int i = 0; i < vel.size(); i++) {
                if (input) {
                    input >> vel[i];
                } else {
                    cout << "velocities in file is not enough to initialize particles!" << endl;
                    throw;
                }
            }
            particles->set_velocity(ptcl_idx, vel);
        }
    } else {
        cout << "can't load velocities from file";
        throw;
    }
}

void ParticlesLoad::position_load() {
    ifstream input(pos_path);
    vector<scalar> tmp;
    array<scalar, 2> pos{};
    string line;
    if (input) {
        for(int ptcl_idx = 0; ptcl_idx < particles->get_Ntot(); ptcl_idx++) {
            for(int i = 0; i < pos.size(); i++) {
                if (input) {
                    input >> pos[i];
                } else {
                    cout << "positions in file is not enough to initialize particles!" << endl;
                    throw;
                }
            }
            particles->set_position(ptcl_idx, pos);
        }
    } else {
        cout << "can't load positions from file";
        throw;
    }
}