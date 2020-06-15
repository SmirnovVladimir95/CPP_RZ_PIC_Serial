#ifndef CPP_RZ_PIC_PARTICLES_H
#define CPP_RZ_PIC_PARTICLES_H

#include <cstddef>
#include <vector>
#include <array>
#include <random>
#include <ctime>
#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"
#include "../Tools/ProjectTypes.h"
using namespace std;

class Particles {
private:
    Grid grid;
    size_t Ntot; // number of macro particles
    scalar ptcls_per_macro; // number of particles in one macro particle
    scalar mfz_const = 0, mfr_const = 0;
    Matrix node_volume;
    void init_node_volume(Matrix& node_volume);
    scalar mass;
    scalar charge;
    static constexpr int pos_dim = 2;
    static constexpr int vel_dim = 3;
public:
    vector<scalar> z;
    vector<scalar> r;
    vector<scalar> vz;
    vector<scalar> vr;
    vector<scalar> vy;
    vector<scalar> efz;
    vector<scalar> efr;
    vector<scalar> mfz;
    vector<scalar> mfr;
    Matrix rho;
    Particles(scalar m, scalar q, size_t N, const Grid& grid, scalar N_per_macro = 1);
    void generate_velocities(scalar energy, int seed=time(nullptr));
    void generate_positions(const array<scalar, pos_dim>& z_bounds, const array<scalar, pos_dim>& r_bounds, int seed=time(nullptr));
    array<scalar, 2> get_position(int ptcl_idx) const;
    array<scalar, 3> get_velocity(int ptcl_idx) const;
    void append(const array<scalar, pos_dim>& position, const array<scalar, vel_dim>& velocity);
    void pop(int ptcl_idx);
    void pusher(scalar  dt);
    void vel_pusher(scalar  dt);
    void electric_field_interpolation(Matrix& Ez, Matrix& Er);
    void magnetic_field_interpolation(Matrix& Bz, Matrix& Br);
    void set_const_magnetic_field(scalar Bz, scalar Br);
    void charge_interpolation();
    int get_Ntot() const;
    void set_position(int ptcl_idx, array<scalar, pos_dim> position);
    void set_velocity(int ptcl_idx, array<scalar, vel_dim> velocity);
    scalar get_mass() const;
    scalar get_charge() const;
    scalar get_ptcl_per_macro() const;
};

scalar uniform_cylindrical(scalar value);

#endif //CPP_RZ_PIC_PARTICLES_H