#include "MonteCarloCollisions.h"

MonteCarloCollisions::MonteCarloCollisions(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles) :
        sigma(sigma), particles(&particles), dt(dt), gas(&gas) {}

MonteCarloCollisions::MonteCarloCollisions(std::unordered_map<type_double, type_double> sigma, type_double dt,
        NeutralGas &gas, Particles &particles) : sigma_energy(sigma), dt(dt), gas(&gas), particles(&particles) {}

ElectronNeutralElasticCollisions::ElectronNeutralElasticCollisions(type_double sigma, type_double dt, NeutralGas& gas,
        Particles& particles) : MonteCarloCollisions(sigma, dt, gas, particles) {}

ElectronNeutralElasticCollisions::ElectronNeutralElasticCollisions(std::unordered_map<type_double, type_double> sigma,
        type_double dt, NeutralGas& gas, Particles& particles) : MonteCarloCollisions(sigma, dt, gas, particles)
        {}

void ElectronNeutralElasticCollisions::velocity_update() {}

IonNeutralElasticCollisions::IonNeutralElasticCollisions(type_double sigma, type_double dt, NeutralGas& gas,
        Particles& particles) : MonteCarloCollisions(sigma, dt, gas, particles) {}

IonNeutralElasticCollisions::IonNeutralElasticCollisions(std::unordered_map<type_double, type_double> sigma,
        type_double dt, NeutralGas& gas, Particles& particles) : MonteCarloCollisions(sigma, dt, gas, particles)
        {}

void IonNeutralElasticCollisions::velocity_update() {}