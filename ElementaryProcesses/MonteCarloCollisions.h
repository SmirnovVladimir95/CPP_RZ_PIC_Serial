#ifndef CPP_RZ_PIC_MONTECARLOCOLLISIONS_H
#define CPP_RZ_PIC_MONTECARLOCOLLISIONS_H


#include <unordered_map>
#include "../Particles/Particles.h"
#include "NeutralGas.h"
using type_double = double;

class MonteCarloCollisions {
private:
    type_double sigma;
    std::unordered_map<type_double, type_double> sigma_energy;
    type_double dt;
    NeutralGas* gas;
    Particles* particles;
public:
    MonteCarloCollisions(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    MonteCarloCollisions(std::unordered_map<type_double, type_double> sigma, type_double dt, NeutralGas& gas,
                         Particles& particles);
};

class ElectronNeutralElasticCollisions : public MonteCarloCollisions {
public:
    ElectronNeutralElasticCollisions(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    ElectronNeutralElasticCollisions(std::unordered_map<type_double, type_double> sigma, type_double dt,
            NeutralGas& gas, Particles& particles);
    void velocity_update();
};

class IonNeutralElasticCollisions : public MonteCarloCollisions {
public:
    IonNeutralElasticCollisions(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles);
    IonNeutralElasticCollisions(std::unordered_map<type_double, type_double> sigma, type_double dt,
            NeutralGas& gas, Particles& particles);
    void velocity_update();
};

#endif //CPP_RZ_PIC_MONTECARLOCOLLISIONS_H
