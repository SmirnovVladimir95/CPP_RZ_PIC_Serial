#include "ParticleLeave.h"


ParticleLeave::ParticleLeave(Particles &particles, Grid& grid, Matrix& domain_condition) :
particles(&particles), grid(&grid), domain_condition(&domain_condition) {}

void ParticleLeave::leave() {
    int ptcl_idx = 0;
    while (ptcl_idx < particles->get_Ntot()) {
        if (not out_of_domain(ptcl_idx))
            ptcl_idx++;
        else
            particles->pop(ptcl_idx);
    }
}

bool ParticleLeave::out_of_domain(int ptcl_idx) {
    int cell_z, cell_r;
    const size_t dim = 2;
    array<scalar, dim> pos = particles->get_position(ptcl_idx);
    cell_z = floor(pos[0] / grid->dz);
    cell_r = floor(pos[1] / grid->dr);
    return domain_condition->operator()(cell_z, cell_r);
}
