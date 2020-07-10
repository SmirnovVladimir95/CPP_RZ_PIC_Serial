#include "Collision.h"
#include <cmath>

#define K_b 1.380649e-23
#define EV 1.6021766208e-19

Collision::Collision(scalar sigma, scalar dt, NeutralGas& gas, Particles& particles) :
        sigma(sigma), particles(&particles), dt(dt), gas(&gas) {
    EnergyCrossSection default_energy_sigma;
    energy_sigma = &default_energy_sigma;
}

Collision::Collision(EnergyCrossSection& energy_sigma, scalar dt,
                     NeutralGas &gas, Particles &particles) : energy_sigma(&energy_sigma), dt(dt), gas(&gas), particles(&particles) {
    sigma = -1;
}

scalar Collision::velocity_module(array<scalar, 3> vel) const {
    return sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
}

array<scalar, 3> Collision::subtraction(array<scalar, 3> vel1, array<scalar, 3> vel2) const {
    return {vel1[0] - vel2[0], vel1[1] - vel2[1], vel1[2] - vel2[2]};
}

array<scalar, 3> Collision::sum(array<scalar, 3> vel1, array<scalar, 3> vel2) const {
    return {vel1[0] + vel2[0], vel1[1] + vel2[1], vel1[2] + vel2[2]};
}

array<scalar, 3> Collision::multiplication_by_constant(array<scalar, 3> vel, scalar value) const {
    return {vel[0]*value, vel[1]*value, vel[2]*value};
}

array<scalar, 3> Collision::isotropic_velocity(scalar vel_module) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<scalar> distribution(0.0, 1.0);
    scalar theta = M_PI*distribution(generator);
    scalar phi = 2*M_PI*distribution(generator);
    return {vel_module*sin(theta)*cos(phi), vel_module*sin(theta)*sin(phi), vel_module*cos(theta)};
}




ElectronNeutralElasticCollision::ElectronNeutralElasticCollision(scalar sigma, scalar dt, NeutralGas& gas,
        Particles& particles) : Collision(sigma, dt, gas, particles) {}

ElectronNeutralElasticCollision::ElectronNeutralElasticCollision(EnergyCrossSection& energy_sigma,
        scalar dt, NeutralGas& gas, Particles& particles) : Collision(energy_sigma, dt, gas, particles) {}

void ElectronNeutralElasticCollision::collision(int ptcl_idx) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<scalar> distribution(0.0, 1.0);
    scalar theta = M_PI*distribution(generator);
    scalar phi = 2*M_PI*distribution(generator);
    scalar vel_module = velocity_module(particles->get_velocity(ptcl_idx));
    if (1 - 2*particles->get_mass()/gas->mass*(1-cos(theta)) < 0) {
        cout << "ion_mass/gas_mass is too small to use electron neutral collision" << endl;
        throw;
    }
    scalar new_vel_module = vel_module*sqrt(1 - 2*particles->get_mass()/gas->mass*(1-cos(theta)));
    array<scalar, 3> new_vel = {new_vel_module*sin(theta)*cos(phi), new_vel_module*sin(theta)*sin(phi),
                                new_vel_module*cos(theta)};
    particles->set_velocity(ptcl_idx, new_vel);
}

scalar ElectronNeutralElasticCollision::probability(int ptcl_idx) const {
    array<scalar, 3> vel = particles->get_velocity(ptcl_idx);
    scalar vel_module = velocity_module(vel);
    scalar energy = 0.5 * particles->get_mass() * vel_module * vel_module;
    if (sigma != -1) {
        return sigma * gas->n * vel_module * dt;
    } else if (not energy_sigma->empty()) {
        return energy_sigma->get_cross_section(energy) * gas->n * vel_module * dt;
    }
    else {
        cout << "cross section is not set!" << endl;
        throw;
    }
}




IonNeutralElasticCollision::IonNeutralElasticCollision(scalar sigma, scalar dt, NeutralGas& gas,
        Particles& particles, bool charge_exchange) : Collision(sigma, dt, gas, particles),
        charge_exchange(charge_exchange) {}

IonNeutralElasticCollision::IonNeutralElasticCollision(EnergyCrossSection& energy_sigma,
        scalar dt, NeutralGas& gas, Particles& particles, bool charge_exchange) :
        Collision(energy_sigma, dt, gas, particles), charge_exchange(charge_exchange) {}

void IonNeutralElasticCollision::collision(int ptcl_idx) {
    //hard sphere approximation
    //cout << "IonNeutralElasticCollision" << endl;
    array<scalar, 3> new_ion_vel{}, gas_vel{};
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0, 1.0);
    float random_num;
    if (charge_exchange)
        random_num = distribution(generator);
    else
        random_num = 0;
    if (random_num <= 0.5) {
        //elastic collision
        array<scalar, 3> ion_vel{}, delta_vel{}, delta_p{};
        scalar ion_mass, gas_mass;
        ion_vel = particles->get_velocity(ptcl_idx);
        gas_vel = gas->generate_velocity();
        ion_mass = particles->get_mass();
        gas_mass = gas->mass;
        for (int i = 0; i < 3; i++) {
            new_ion_vel[i] = (ion_vel[i] * ion_mass + gas_vel[i] * gas_mass + gas_mass * (gas_vel[i] - ion_vel[i])) /
                    (ion_mass + gas_mass);
        }
        particles->set_velocity(ptcl_idx, new_ion_vel);
    } else {
        //charge exchange
        new_ion_vel = gas->generate_velocity();
        particles->set_velocity(ptcl_idx, new_ion_vel);
    }
}

scalar IonNeutralElasticCollision::probability(int ptcl_idx) const {
    scalar vel_module = velocity_module(particles->get_velocity(ptcl_idx));
    if (sigma != -1) {
        scalar relative_vel, R_B, nu_A;
        R_B = K_b/gas->mass;
        nu_A = vel_module/(sqrt(2*R_B*gas->T));
        relative_vel = sqrt(2*R_B*gas->T)*((nu_A+1/(2*nu_A))*erf(nu_A) + 1/(sqrt(M_PI))*exp(-1*nu_A*nu_A));
        return gas->n*relative_vel*sigma*dt;
    } else if (not energy_sigma->empty()) {
        scalar energy = 0.5 * particles->get_mass() * vel_module * vel_module;
        scalar g_sigma = energy_sigma->get_cross_section(energy);
        return g_sigma * gas->n * dt;
    } else {
        cout << "cross section is not set!" << endl;
        throw;
    }
}




Ionization::Ionization(scalar sigma, scalar ion_threshold, scalar dt, NeutralGas& gas,
        Particles& electrons, Particles& ions) : Collision(sigma, dt, gas, electrons),
        ion_threshold(ion_threshold), ionized_particles(&ions) {}

Ionization::Ionization(EnergyCrossSection& energy_sigma, scalar dt, NeutralGas &gas,
        Particles& electrons, Particles& ions) : Collision(energy_sigma, dt, gas, electrons),
        ionized_particles(&ions) {
    ion_threshold = energy_sigma.get_energy_cross_section(0)[0];
}

void Ionization::collision(int ptcl_idx) {
    //cout << "Ionization process" << endl;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    scalar delta_E, mu;
    array<scalar, 3> electron_vel{}, gas_vel{}, electron_vel_new{}, delta_vel{}, emitted_electron_vel{}, ion_vel{},
            center_mass_vel{};

    mu = particles->get_mass() * gas->mass / (particles->get_mass() + gas->mass);
    electron_vel = particles->get_velocity(ptcl_idx);
    gas_vel = gas->generate_velocity();
    delta_vel = {electron_vel[0] - gas_vel[0], electron_vel[1] - gas_vel[1], electron_vel[2] - gas_vel[2]};
    scalar delta_vel_module = velocity_module(delta_vel);
    delta_E = 0.5 * mu * delta_vel_module * delta_vel_module - ion_threshold;

    if (delta_E <= 0) {
        return;
    }

    scalar random_num = distribution(generator);
    electron_vel_new = isotropic_velocity(sqrt(random_num * delta_E * 2 / particles->get_mass()));
    emitted_electron_vel = isotropic_velocity(sqrt((1 - random_num) * delta_E * 2 / particles->get_mass()));
    scalar m_div_M = particles->get_mass() / ionized_particles->get_mass();
    ion_vel = multiplication_by_constant(sum(electron_vel_new, emitted_electron_vel), -1 * m_div_M);

    array<scalar, 3> p_electron = multiplication_by_constant(electron_vel, particles->get_mass());
    array<scalar, 3> p_gas = multiplication_by_constant(gas_vel, gas->mass);
    center_mass_vel = multiplication_by_constant(sum(p_electron, p_gas), 1 / (particles->get_mass() + gas->mass));

    electron_vel_new = sum(electron_vel_new, center_mass_vel);
    emitted_electron_vel = sum(emitted_electron_vel, center_mass_vel);
    ion_vel = sum(ion_vel, center_mass_vel);

    particles->set_velocity(ptcl_idx, electron_vel_new);
    particles->append(particles->get_position(ptcl_idx), emitted_electron_vel);
    ionized_particles->append(particles->get_position(ptcl_idx), ion_vel);

}

scalar Ionization::probability(int ptcl_idx) const {
    array<scalar, 3> vel = particles->get_velocity(ptcl_idx);
    scalar vel_module = velocity_module(vel);
    scalar energy = 0.5*vel_module*vel_module*particles->get_mass();
    if (sigma != -1) {
        if (energy >= ion_threshold) {
            return sigma * gas->n * vel_module * dt;
        }
        return 0;
    } else if (not energy_sigma->empty()) {
        return energy_sigma->get_cross_section(energy) * gas->n * vel_module * dt;
    } else {
        cout << "cross section is not set!" << endl;
        throw;
    }
}
