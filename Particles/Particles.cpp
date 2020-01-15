#include "Particles.h"
#include <random>
#include <cmath>
#include "../Tools/Matrix.h"
#include "Pusher.h"
#include "Interpolation.h"

Particles::Particles(scalar m, scalar q, size_t N, const Grid& init_grid, bool volume_linear_correction) {
    mass = m;
    charge = q;
    z.resize(N, 0);
    r.resize(N, 0);
    vz.resize(N, 0);
    vr.resize(N, 0);
    vy.resize(N, 0);
    efz.resize(N, 0);
    efr.resize(N, 0);
    mfz.resize(N, 0);
    mfr.resize(N, 0);
    Ntot = N;
    grid = init_grid;
    node_volume.resize(init_grid.Nz, init_grid.Nr);
    init_node_volume(node_volume);
    rho.resize(init_grid.Nz, init_grid.Nr);
}

void Particles::generate_velocities(scalar energy, int seed) {
    std::default_random_engine generator(seed);
    std::normal_distribution<scalar> distribution(0.0, sqrt(2*energy/(3*mass)));
    for (int i = 0; i < Ntot; i++) {
        vz[i] = distribution(generator);
        vr[i] = distribution(generator);
        vy[i] = distribution(generator);
    }
}

void Particles::generate_positions(const array<scalar, 2> &z_bounds, const array<scalar, 2> &r_bounds,
                                   const int seed) {
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<scalar> distribution_z(z_bounds[0],z_bounds[1]);
    std::uniform_real_distribution<scalar> distribution_r(r_bounds[0],r_bounds[1]);
    for (int i = 0; i < Ntot; i++) {
        z[i] = distribution_z(generator);
        r[i] = distribution_r(generator);
    }
}

void Particles::vel_pusher(scalar  dt) {
    auto vel_z = vz.data();
    auto vel_r = vr.data();
    auto vel_y = vy.data();
    auto Ez = efz.data();
    auto Er = efr.data();
    auto Bz = mfz.data();
    auto Br = mfr.data();
    UpdateVelocity(vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, charge, mass, Ntot);
}

void Particles::pusher(scalar dt) {
    auto pos_z = z.data();
    auto pos_r = r.data();
    auto vel_z = vz.data();
    auto vel_r = vr.data();
    auto vel_y = vy.data();
    auto Ez = efz.data();
    auto Er = efr.data();
    auto Bz = mfz.data();
    auto Br = mfr.data();
    ParticlePush(pos_z, pos_r, vel_z, vel_r, vel_y, Ez, Er, Bz, Br, dt, charge, mass,
                 Ntot, grid.dr);
}

void Particles::init_node_volume(Matrix& node_volume) {
    InitVolume(node_volume, grid);
}

void Particles::electric_field_interpolation(Matrix& Ez, Matrix& Er) {
    LinearFieldInterpolation(efz.data(), efr.data(), z.data(), r.data(), Ez.data_ptr(), Er.data_ptr(),
                             grid, Ntot);
}

void Particles::charge_interpolation() {
    LinearChargeInterpolation(rho.data_ptr(), z.data(), r.data(), grid, charge, Ntot, node_volume.data_ptr());
}

void Particles::set_const_magnetic_field(scalar Bz, scalar Br) {
    mfz.assign(Ntot, Bz);
    mfr.assign(Ntot, Br);
    mfz_const = Bz;
    mfr_const = Br;
}

void Particles::append(const array<scalar, 2> &position,const array<scalar, 3> &velocity) {
    z.push_back(position[0]);
    r.push_back(position[1]);
    vz.push_back(velocity[0]);
    vr.push_back(velocity[1]);
    vy.push_back(velocity[2]);
    efz.push_back(0);
    efr.push_back(0);
    mfz.push_back(mfz_const);
    mfr.push_back(mfr_const);
    Ntot++;
}

void swap(scalar& a, scalar& b) {
    scalar temp = a;
    a = b;
    b = temp;
}

void Particles::pop(int ptcl_idx) {
    swap(z[ptcl_idx], z[Ntot-1]);
    swap(r[ptcl_idx], r[Ntot-1]);
    swap(vz[ptcl_idx], vz[Ntot-1]);
    swap(vr[ptcl_idx], vr[Ntot-1]);
    swap(vy[ptcl_idx], vy[Ntot-1]);
    swap(efz[ptcl_idx], efz[Ntot-1]);
    swap(efr[ptcl_idx], efr[Ntot-1]);
    swap(mfz[ptcl_idx], mfz[Ntot-1]);
    swap(mfr[ptcl_idx], mfr[Ntot-1]);
    z.pop_back();
    r.pop_back();
    vz.pop_back();
    vr.pop_back();
    vy.pop_back();
    efz.pop_back();
    efr.pop_back();
    mfz.pop_back();
    mfr.pop_back();
    Ntot--;
}

array<scalar, 2> Particles::get_position(int ptcl_idx) const {
    array<scalar, 2> pos = {z[ptcl_idx], r[ptcl_idx]};
    return pos;
}

array<scalar, 3> Particles::get_velocity(int ptcl_idx) const {
    array<scalar, 3> vel = {vz[ptcl_idx], vr[ptcl_idx], vy[ptcl_idx]};
    return vel;
}

size_t Particles::get_Ntot() const {
    return Ntot;
}

void Particles::set_velocity(int ptcl_idx, array<scalar, 3> velocity) {
    vz[ptcl_idx] = velocity[0];
    vr[ptcl_idx] = velocity[1];
    vy[ptcl_idx] = velocity[2];
}
