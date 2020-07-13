#include "Dielectric.h"

Dielectric::Dielectric(Particles &incident_particles, Grid &grid, Matrix &domain_condition, scalar dielectric_const) :
incident_particles(&incident_particles), grid(&grid), domain_condition(&domain_condition),
dielectric_const(dielectric_const) {
    rho.resize(grid.Nz, grid.Nr);
}

void Dielectric::update_rho() {
    int cell_z, cell_r;
    scalar hz, hr;
    scalar charge = incident_particles->get_charge() * incident_particles->get_ptcl_per_macro() / dielectric_const;
    const size_t dim = 2;
    array<scalar, dim> pos{};
    int ptcl_idx = 0;
    while (ptcl_idx < incident_particles->get_Ntot()) {
        pos = incident_particles->get_position(ptcl_idx);
        cell_z = floor(pos[0] / grid->dz);
        cell_r = floor(pos[1] / grid->dr);
        if (domain_condition->operator()(cell_z, cell_r) == RIGHT) {
            hr = (pos[1] - cell_r * grid->dr) / grid->dr;
            rho(cell_z+1, cell_r) += charge * (1 - hr) / incident_particles->node_volume(cell_z+1, cell_r);
            rho(cell_z+1, cell_r+1) += charge * hr / incident_particles->node_volume(cell_z+1, cell_r+1);
            incident_particles->pop(ptcl_idx);
        } else if (domain_condition->operator()(cell_z, cell_r)  == LEFT) {
            hr = (pos[1] - cell_r * grid->dr) / grid->dr;
            rho(cell_z, cell_r) += charge * (1 - hr) / incident_particles->node_volume(cell_z, cell_r);
            rho(cell_z, cell_r+1) += charge * hr / incident_particles->node_volume(cell_z, cell_r+1);
            incident_particles->pop(ptcl_idx);
        } else if (domain_condition->operator()(cell_z, cell_r)  == UP) {
            hz = (pos[0] - cell_z * grid->dz) / grid->dz;
            rho(cell_z, cell_r+1) += charge * (1 - hz) / incident_particles->node_volume(cell_z, cell_r+1);
            rho(cell_z+1, cell_r+1) += charge * hz / incident_particles->node_volume(cell_z+1, cell_r+1);;
            incident_particles->pop(ptcl_idx);
        } else if (domain_condition->operator()(cell_z, cell_r)  == DOWN) {
            hz = (pos[0] - cell_z * grid->dz) / grid->dz;
            rho(cell_z, cell_r) += charge * (1 - hz) / incident_particles->node_volume(cell_z, cell_r);
            rho(cell_z+1, cell_r) += charge * hz / incident_particles->node_volume(cell_z+1, cell_r);
            incident_particles->pop(ptcl_idx);
        }
        else {
            ptcl_idx++;
        }
    }
}
