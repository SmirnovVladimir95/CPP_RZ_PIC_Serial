#ifndef CPP_RZ_PIC_COLLISION_H
#define CPP_RZ_PIC_COLLISION_H


#include <unordered_map>
#include "../Particles/Particles.h"
#include "NeutralGas.h"
using type_double = double;

class Collision {
protected:
    type_double sigma=0;
    std::unordered_map<type_double, type_double> sigma_energy;
    type_double dt;
    NeutralGas* gas;
    type_double velocity_module(array<type_double, 3> vel) const;
    array<type_double, 3> sum(array<type_double, 3> vel1, array<type_double, 3> vel2) const;
    array<type_double, 3> subtraction(array<type_double, 3> vel1, array<type_double, 3> vel2) const;
    array<type_double, 3> multiplication_by_constant(array<type_double, 3> vel, type_double value) const;
    array<type_double, 3> isotropic_velocity(type_double vel_module);
public:
    Particles* particles;
    Collision(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    Collision(const std::unordered_map<type_double, type_double>& sigma, type_double dt, NeutralGas& gas,
              Particles& particles);
};

class ElectronNeutralElasticCollision : public Collision {
public:
    ElectronNeutralElasticCollision(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    ElectronNeutralElasticCollision(std::unordered_map<type_double, type_double> sigma, type_double dt,
            NeutralGas& gas, Particles& particles);
    void collision(int ptcl_idx);
    type_double probability(int ptcl_idx) const;
};

class IonNeutralElasticCollision : public Collision {
    // Hard-Sphere Interaction
public:
    IonNeutralElasticCollision(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    IonNeutralElasticCollision(std::unordered_map<type_double, type_double> sigma, type_double dt,
            NeutralGas& gas, Particles& particles);
    void collision(int ptcl_idx);
    type_double probability(int ptcl_idx) const;
};

class Ionization : public Collision {
private:
    type_double ion_threshold = -1;
    Particles* ionized_particles;
public:
    Ionization(type_double sigma, type_double ion_threshold, type_double dt, NeutralGas& gas, Particles& electrons,
            Particles& ions);
    Ionization(const std::unordered_map<type_double, type_double>& sigma, type_double dt, NeutralGas& gas,
            Particles& electrons, Particles& ions);
    void collision(int ptcl_idx);
    type_double probability(int ptcl_idx) const;
};

#endif //CPP_RZ_PIC_COLLISION_H
