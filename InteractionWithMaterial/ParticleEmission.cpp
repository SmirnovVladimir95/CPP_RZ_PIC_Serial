//
// Created by Vladimir Smirnov on 30.12.2019.
//

#include "ParticleEmission.h"


ParticleEmission::ParticleEmission(Particles &incident_particles, Particles &emitted_particles, Grid &grid,
                                   Matrix &domain_condition, array<scalar, 3> emission_direction, scalar gamma,
                                   scalar emission_energy) :
                                   incident_particles(&incident_particles), emitted_particles(&emitted_particles),
                                   grid(&grid), domain_condition(&domain_condition),
                                   emission_direction(emission_direction), gamma(gamma),
                                   emission_energy(emission_energy) {}

void ParticleEmission::emission(int seed) {
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<scalar> distribution(0.0,1.0);
    scalar vel_module = sqrt(2*emission_energy/emitted_particles->get_mass());
    int ptcl_idx = 0;
    while (ptcl_idx < incident_particles->get_Ntot()) {
        if (emission_condition(ptcl_idx)) {
            if (distribution(generator) < gamma) {
                array<scalar, 2> pos = incident_particles->get_position(ptcl_idx);
                array<scalar, 3> vel = {emission_direction[0]*vel_module,
                                        emission_direction[1]*vel_module,
                                        emission_direction[2]*vel_module};
                emitted_particles->append(pos, vel);
            }
            incident_particles->pop(ptcl_idx);
        } else {
            ptcl_idx++;
        }
    }
}

bool ParticleEmission::emission_condition(int ptcl_idx) {
    int cell_z, cell_r;
    const size_t dim = 2;
    array<scalar, dim> pos = incident_particles->get_position(ptcl_idx);
    cell_z = floor(pos[0] / grid->dz);
    cell_r = floor(pos[1] / grid->dr);
    return domain_condition->operator()(cell_z, cell_r);
}
