#ifndef CPP_RZ_PIC_PARTICLES_H
#define CPP_RZ_PIC_PARTICLES_H

#include <cstddef>
#include <vector>
#include <array>
#include <random>
#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"
using namespace std;

typedef double type_double; // c++0x

class Particles {
private:
    Grid grid;
    size_t Ntot;
    void init_node_volume(Matrix& node_volume);
public:
    type_double mass;
    type_double charge;
    vector<type_double> z;
    vector<type_double> r;
    vector<type_double> vz;
    vector<type_double> vr;
    vector<type_double> vy;
    vector<type_double> efz;
    vector<type_double> efr;
    vector<type_double> mfz;
    vector<type_double> mfr;
    Matrix rho;
    Matrix node_volume;
    Particles(type_double m, type_double q, size_t N, const Grid& grid);
    void generate_velocities(type_double energy, int seed=time(nullptr));
    void generate_positions(const array<type_double, 2>& z_bounds, const array<type_double, 2>& r_bounds, int seed=time(nullptr));
    //vector<vector<type_double>> get_positions() const;
    //vector<vector<type_double>> get_velocities() const;
    void append(const array<type_double, 2>& position, const array<type_double, 3>& velocity);
    void pusher(type_double  dt);
    void vel_pusher(type_double  dt);
    void electric_field_interpolation(Matrix& Ez, Matrix& Er);
    void magnetic_field_interpolation(Matrix& Bz, Matrix& Br);
    void set_const_magnetic_field(const vector<type_double>& Bz, const vector<type_double>& Br);
    void charge_interpolation();
    size_t get_Ntot() const;
};


#endif //CPP_RZ_PIC_PARTICLES_H