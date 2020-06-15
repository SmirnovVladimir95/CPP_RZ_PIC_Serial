#ifndef CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H
#define CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H

#include "../ElementaryProcesses/Collision.h"
#include "../ElementaryProcesses/NeutralGas.h"
#include "../Grid/Grid.h"
#include "../ElementaryProcesses/NanbuCollisions.h"
#include <ctime>
#define K_B 1.380649e-23

void test_NanbuCollisionChoice() {
    cout << "test_NanbuCollisionChoice: " << endl;
    vector<scalar> prob = {1e-2, 1e-3};
    int prob_num = 2;
    int num = 1e6;
    int count_null = 0, count_one = 0, count_second = 0, coll_type;
    for (int i = 0; i < num; i++) {
        coll_type = NanbuCollisionChoice(prob, prob_num);
        //cout << coll_type << " ";
        if (coll_type == 0)
            count_null++;
        if (coll_type == 1)
            count_one++;
        if (coll_type == 2)
            count_second++;
    }
    cout << "no collision: " << count_null << endl;
    cout << "collision 1: " << count_one << endl;
    cout << "collision 2: " << count_second << endl;
    cout << "no collision + collision 1 + collision 2: " << count_null+count_one+count_second << endl;
    cout << "OK" << endl;
}

void test_MonteCarloCollision() {
    cout << "test_MonteCarloCollision: ";

    // Electron Elastic Collision
    scalar sigma=1e-19, dt=1e-9, mass=9.1e-31, charge=1.6e-19;
    scalar n = 1e20, m_gas = 100*mass, T = 500;
    NeutralGas gas(n, m_gas, T);
    size_t Ntot=1e3, Nz=100, Nr=50;
    scalar dz=2e-5, dr=2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles electrons(mass, -1*charge, Ntot, grid);
    electrons.generate_velocities(1*1.6e-19);
    ElectronNeutralElasticCollision electron_elastic(sigma, dt, gas, electrons);
    cout << electron_elastic.probability(0) << endl;
    cout << "vel before collision: " <<
    sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]) << endl;
    cout << electrons.vz[0] << " " << electrons.vr[0] << " " << electrons.vy[0] << endl;
    electron_elastic.collision(0);
    cout << "vel after collision: " <<
    sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]) << endl;
    cout << electrons.vz[0] << " " << electrons.vr[0] << " " << electrons.vy[0] << endl;

    // Ion Elastic Collision
    Particles ions(m_gas, charge, Ntot, grid);
    scalar energy = (3 / 2) * K_B * 500;
    ions.generate_velocities(energy);
    IonNeutralElasticCollision ion_elastic(1e-19, dt, gas, ions, true);
    cout << ion_elastic.probability(0) << endl;
    cout << "vel before collision: " <<
         sqrt(ions.vz[0]*ions.vz[0] + ions.vr[0]*ions.vr[0] + ions.vy[0]*ions.vy[0]) << endl;
    cout << ions.vz[0] << " " << ions.vr[0] << " " << ions.vy[0] << endl;
    ion_elastic.collision(0);
    cout << "vel after collision: " <<
        sqrt(ions.vz[0]*ions.vz[0] + ions.vr[0]*ions.vr[0] + ions.vy[0]*ions.vy[0]) << endl;
    cout << ions.vz[0] << " " << ions.vr[0] << " " << ions.vy[0] << endl;

    // Electron Ionization
    scalar ion_threshold = 0.1*1.6e-19, dt_collision = 1e-9;
    Ionization gas_ionization(1e-20, ion_threshold, dt_collision, gas, electrons, ions);
    cout << "ionization prob: " << gas_ionization.probability(0) << endl;
    cout << "vel of incident electron before ionization: " <<
         sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]) << endl;
    cout << electrons.vz[0] << " " << electrons.vr[0] << " " << electrons.vy[0] << endl;
    gas_ionization.collision(0);
    cout << "Ntot after ion ion/electron: " << ions.get_Ntot() << " " << electrons.get_Ntot() << endl;
    cout << "vel of new ion after ionization: " <<
         sqrt(ions.vz[Ntot-1]*ions.vz[Ntot-1] + ions.vr[Ntot-1]*ions.vr[Ntot-1] + ions.vy[Ntot-1]*ions.vy[Ntot-1]) << endl;
    cout << ions.vz[Ntot-1] << " " << ions.vr[Ntot-1] << " " << ions.vy[Ntot-1] << endl;
    cout << "vel of incident electron after ionization: " <<
        sqrt(electrons.vz[0]*electrons.vz[0] + electrons.vr[0]*electrons.vr[0] + electrons.vy[0]*electrons.vy[0]) << endl;
    cout << electrons.vz[0] << " " << electrons.vr[0] << " " << electrons.vy[0] << endl;
    cout << "vel of new electron after ionization: " <<
         sqrt(electrons.vz[Ntot-1]*electrons.vz[Ntot-1] + electrons.vr[Ntot-1]*electrons.vr[Ntot-1] +
         electrons.vy[Ntot-1]*electrons.vy[Ntot-1]) << endl;
    cout << electrons.vz[Ntot-1] << " " << electrons.vr[Ntot-1] << " " << electrons.vy[Ntot-1] << endl;

    cout << "OK" << endl;
}

void test_MonteCarloCollisions() {
    cout << "test_MonteCarloCollisions: ";
    scalar n = 1e20, mAr = 6.6335209e-26, T = 500;
    NeutralGas gas(n, mAr, T);

    scalar dz=2e-5, dr=2e-5;
    size_t Nz=100, Nr=50;
    Grid grid(Nz, Nr, dz, dr);

    Particles electrons(9.1e-31, 1e-19, 1e6, grid);
    electrons.generate_velocities(10*1.6e-19);
    array<scalar, 2> z_bounds = {0, Nz*dz}, r_bounds = {0, Nr*dr};
    electrons.generate_positions(z_bounds, r_bounds);

    Particles ions(mAr, 1.6e-19, 1e6, grid);
    ions.generate_velocities(T*1.38e-23);
    ions.generate_positions(z_bounds, r_bounds);

    ElectronNeutralElasticCollision electron_elastic(1e-19, 1e-10, gas, electrons);
    scalar ion_threshold = 10*1.6e-19;
    Ionization argon_ionization(1e-20, ion_threshold, 1e-10, gas, electrons, ions);
    IonNeutralElasticCollision ion_elastic(1e-19, 1-10, gas, ions);

    cout << "Number of ptcls before collisions: " << electrons.get_Ntot() << " " << ions.get_Ntot() << endl;
    NanbuElectronCollisionProcess(electron_elastic, argon_ionization, 2);
    NanbuIonCollisionProcess(ion_elastic, 1);
    cout << "Number of ptcls after collisions: " << electrons.get_Ntot() << " " << ions.get_Ntot() << endl;
    cout << "OK" << endl;
}


void test_Collisions() {
    cout << "test_Collisions: ";
    clock_t start = clock();
    test_NanbuCollisionChoice();
    clock_t end = clock();
    scalar seconds = (scalar)(end - start) / CLOCKS_PER_SEC;
    cout << seconds << endl;
    test_MonteCarloCollision();
    clock_t start1 = clock();
    test_MonteCarloCollisions();
    clock_t end1 = clock();
    scalar seconds1 = (scalar)(end1 - start1) / CLOCKS_PER_SEC;
    cout << seconds1 << endl;
    cout << "OK" << endl;
}

#endif //CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H
