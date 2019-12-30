#include "NeutralGas.h"
#include <random>
#include <cmath>
#include <array>
#define K_b 1.380649e-23

using namespace std;

NeutralGas::NeutralGas(type_double n, type_double mass, type_double T) : n(n), mass(mass), T(T) {}

array<type_double, 3> NeutralGas::generate_velocity() const {
    std::random_device seed;
    type_double sigma = sqrt(K_b*T/mass);
    std::default_random_engine generator(seed());
    std::normal_distribution<type_double> distribution(0.0, sigma);
    array<type_double, 3> velocity = {distribution(generator), distribution(generator), distribution(generator)};
    return velocity;
}