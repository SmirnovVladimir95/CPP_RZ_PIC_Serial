#include "NeutralGas.h"
#include <random>
#include <cmath>
#define K_b 1.380649e-23

using namespace std;

NeutralGas::NeutralGas(scalar n, scalar mass, scalar T) : n(n), mass(mass), T(T) {}

array<scalar, 3> NeutralGas::generate_velocity() const {
    std::random_device seed;
    scalar sigma = sqrt(K_b*T/mass);
    std::default_random_engine generator(seed());
    std::normal_distribution<scalar> distribution(0.0, sigma);
    array<scalar, 3> velocity = {distribution(generator), distribution(generator), distribution(generator)};
    return velocity;
}