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
    size_t Ntot;
    scalar mfz_const = 0, mfr_const = 0;
    scalar* node_volume;
    //Matrix node_volume;
    void init_node_volume(scalar node_volume[]);
public:
    scalar mass;
    scalar charge;
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
    Particles(scalar m, scalar q, size_t N, const Grid& grid, bool volume_linear_correction = true);
    ~Particles();
    void generate_velocities(scalar energy, int seed=time(nullptr));
    void generate_positions(const array<scalar, 2>& z_bounds, const array<scalar, 2>& r_bounds, int seed=time(nullptr));
    array<scalar, 2> get_position(int ptcl_idx) const;
    array<scalar, 3> get_velocity(int ptcl_idx) const;
    void append(const array<scalar, 2>& position, const array<scalar, 3>& velocity);
    void pop(int ptcl_idx);
    void pusher(scalar  dt);
    void vel_pusher(scalar  dt);
    void electric_field_interpolation(Matrix& Ez, Matrix& Er);
    void magnetic_field_interpolation(Matrix& Bz, Matrix& Br);
    void set_const_magnetic_field(scalar Bz, scalar Br);
    void charge_interpolation();
    size_t get_Ntot() const;
    void set_velocity(int ptcl_idx, array<scalar, 3> velocity);
};


#endif //CPP_RZ_PIC_PARTICLES_H