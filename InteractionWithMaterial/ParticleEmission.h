#ifndef CPP_RZ_PIC_PARTICLEEMISSION_H
#define CPP_RZ_PIC_PARTICLEEMISSION_H

#include "../Grid/Grid.h"
#include "../Particles/Particles.h"

class ParticleEmission {
private:
    Grid* grid;
    Particles* incident_particles;
    Particles* emitted_particles;
    array<scalar, 3> emission_direction;
    Matrix* domain_condition;
    scalar gamma, emission_energy;
    bool emission_condition(int ptcl_idx);
public:
    ParticleEmission(Particles& incident_particles, Particles& emitted_particles, Grid& grid, Matrix& domain_condition,
                     array<scalar, 3> emission_direction, scalar gamma, scalar emission_energy);
    void emission(int seed = 0);
};


#endif //CPP_RZ_PIC_PARTICLEEMISSION_H
