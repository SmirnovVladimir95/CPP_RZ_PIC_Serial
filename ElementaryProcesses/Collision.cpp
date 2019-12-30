#include "Collision.h"
#include <cmath>
#define K_b 1.380649e-23

Collision::Collision(type_double sigma, type_double dt, NeutralGas& gas, Particles& particles) :
        sigma(sigma), particles(&particles), dt(dt), gas(&gas) {}

Collision::Collision(const std::unordered_map<type_double, type_double>& sigma, type_double dt,
                     NeutralGas &gas, Particles &particles) : sigma_energy(sigma), dt(dt), gas(&gas), particles(&particles) {}

type_double Collision::velocity_module(array<type_double, 3> vel) const {
    return sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
}

array<type_double, 3> Collision::subtraction(array<type_double, 3> vel1, array<type_double, 3> vel2) const {
    return {vel1[0] - vel2[0], vel1[1] - vel2[1], vel1[2] - vel2[2]};
}

array<type_double, 3> Collision::sum(array<type_double, 3> vel1, array<type_double, 3> vel2) const {
    return {vel1[0] + vel2[0], vel1[1] + vel2[1], vel1[2] + vel2[2]};
}

array<type_double, 3> Collision::multiplication_by_constant(array<type_double, 3> vel, type_double value) const {
    return {vel[0]*value, vel[1]*value, vel[2]*value};
}

array<type_double, 3> Collision::isotropic_velocity(type_double vel_module) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    float theta = M_PI*distribution(generator);
    float phi = M_PI*distribution(generator);
    return {vel_module*cos(theta)*cos(phi), vel_module*cos(theta)*sin(phi), vel_module*sin(theta)};
}




ElectronNeutralElasticCollision::ElectronNeutralElasticCollision(type_double sigma, type_double dt, NeutralGas& gas,
        Particles& particles) : Collision(sigma, dt, gas, particles) {}

ElectronNeutralElasticCollision::ElectronNeutralElasticCollision(std::unordered_map<type_double, type_double> sigma,
        type_double dt, NeutralGas& gas, Particles& particles) : Collision(sigma, dt, gas, particles) {}

void ElectronNeutralElasticCollision::collision(int ptcl_idx) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    float theta = M_PI*distribution(generator);
    float phi = M_PI*distribution(generator);

    type_double vel_module = velocity_module(particles->get_velocity(ptcl_idx));
    type_double new_vel_module = vel_module*sqrt(1 - 2*particles->mass/gas->mass*(1-cos(theta)));
    particles->vz[ptcl_idx] = new_vel_module*cos(theta)*cos(phi);
    particles->vr[ptcl_idx] = new_vel_module*cos(theta)*sin(phi);
    particles->vy[ptcl_idx] = new_vel_module*sin(theta);
}

type_double ElectronNeutralElasticCollision::probability(int ptcl_idx) const {
    if (sigma != 0) {
        array<type_double, 3> vel = particles->get_velocity(ptcl_idx);
        type_double vel_module = velocity_module(vel);
        return sigma * gas->n * vel_module * dt;
    } else if (sigma_energy.size() != 0) {
        // do smth
        return 0;
    }
    else {
        return -1;
    }
}




IonNeutralElasticCollision::IonNeutralElasticCollision(type_double sigma, type_double dt, NeutralGas& gas,
        Particles& particles) : Collision(sigma, dt, gas, particles) {}

IonNeutralElasticCollision::IonNeutralElasticCollision(std::unordered_map<type_double, type_double> sigma,
        type_double dt, NeutralGas& gas, Particles& particles) : Collision(sigma, dt, gas, particles)
        {}

void IonNeutralElasticCollision::collision(int ptcl_idx) {
    //hard sphere approximation
    array<type_double, 3> new_ion_vel, gas_vel;
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0, 1.0);
    float random_num = distribution(generator);

    if (random_num <= 0.5) {
        //elastic collision
        array<type_double, 3> ion_vel, delta_vel, delta_p;
        type_double delta_p_module;
        ion_vel = particles->get_velocity(ptcl_idx);
        gas_vel = gas->generate_velocity();
        delta_vel = {ion_vel[0] - gas_vel[0], ion_vel[1] - gas_vel[1], ion_vel[2] - gas_vel[2]};
        delta_p_module = velocity_module(delta_vel) * gas->mass;
        delta_p = isotropic_velocity(delta_p_module);
        new_ion_vel[0] =
                (ion_vel[0] * particles->mass + gas_vel[0] * gas->mass + delta_p[0]) / (particles->mass + gas->mass);
        new_ion_vel[1] =
                (ion_vel[1] * particles->mass + gas_vel[1] * gas->mass + delta_p[1]) / (particles->mass + gas->mass);
        new_ion_vel[2] =
                (ion_vel[2] * particles->mass + gas_vel[2] * gas->mass + delta_p[2]) / (particles->mass + gas->mass);
        particles->set_velocity(ptcl_idx, new_ion_vel);
    } else {
        //charge exchange
        new_ion_vel = gas->generate_velocity();
        particles->set_velocity(ptcl_idx, new_ion_vel);
    }
}

type_double IonNeutralElasticCollision::probability(int ptcl_idx) const {
    type_double relative_vel, R_B, nu_A, vel_module;
    R_B = K_b/gas->mass;
    vel_module = velocity_module(particles->get_velocity(ptcl_idx));
    nu_A = vel_module/(sqrt(2*R_B*gas->T));
    relative_vel = sqrt(2*R_B*gas->T)*((nu_A+1/(2*nu_A))*erf(nu_A) + 1/(sqrt(M_PI))*exp(-1*nu_A*nu_A));
    return gas->n*relative_vel*sigma*dt;
}




Ionization::Ionization(type_double sigma, type_double ion_threshold, type_double dt, NeutralGas& gas,
        Particles& electrons, Particles& ions) : Collision(sigma, dt, gas, electrons),
        ion_threshold(ion_threshold), ionized_particles(&ions) {}

Ionization::Ionization(const std::unordered_map<type_double, type_double> &sigma, type_double dt, NeutralGas &gas,
        Particles& electrons, Particles& ions) : Collision(sigma, dt, gas, electrons),
        ionized_particles(&ions) {}

void Ionization::collision(int ptcl_idx) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0, 1.0);
    type_double delta_E, mu;
    array<type_double, 3> electron_vel, gas_vel, electron_vel_new, delta_vel, emitted_electron_vel, ion_vel,
    center_mass_vel;

    mu = particles->mass * gas->mass/(particles->mass + gas->mass);
    electron_vel = particles->get_velocity(ptcl_idx);
    gas_vel = gas->generate_velocity();
    delta_vel = {electron_vel[0]-gas_vel[0], electron_vel[1]-gas_vel[1], electron_vel[2]-gas_vel[2]};
    type_double delta_vel_module = velocity_module(delta_vel);
    delta_E = 0.5 * mu * delta_vel_module * delta_vel_module - ion_threshold;

    type_double random_num = distribution(generator);
    electron_vel_new = isotropic_velocity(sqrt(random_num * delta_E * 2 / particles->mass));
    emitted_electron_vel = isotropic_velocity(sqrt((1-random_num) * delta_E * 2 / particles->mass));
    type_double m_div_M = particles->mass/ionized_particles->mass;
    ion_vel = multiplication_by_constant(sum(electron_vel_new, emitted_electron_vel), -1*m_div_M);

    array<type_double, 3> p_electron = multiplication_by_constant(electron_vel, particles->mass);
    array<type_double, 3> p_gas = multiplication_by_constant(gas_vel, gas->mass);
    center_mass_vel = multiplication_by_constant(sum(p_electron, p_gas), 1/(particles->mass+gas->mass));

    electron_vel_new = sum(electron_vel_new, center_mass_vel);
    emitted_electron_vel = sum(emitted_electron_vel, center_mass_vel);
    ion_vel = sum(ion_vel, center_mass_vel);

    particles->set_velocity(ptcl_idx, electron_vel_new);
    particles->append(particles->get_position(ptcl_idx), emitted_electron_vel);
    ionized_particles->append(particles->get_position(ptcl_idx), ion_vel);
}

type_double Ionization::probability(int ptcl_idx) const {
    if (sigma != 0) {
        array<type_double, 3> vel = particles->get_velocity(ptcl_idx);
        type_double vel_module = velocity_module(vel);
        if (0.5*vel_module*vel_module*particles->mass >= ion_threshold)
            return sigma * gas->n * vel_module * dt;
        return 0;
    } else if (sigma_energy.size() != 0) {
        //do smth
    }
    return 0;
}
