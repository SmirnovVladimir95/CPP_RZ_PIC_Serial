#ifndef CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H
#define CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H

#include "../ElementaryProcesses/MonteCarloCollisions.h"
#include "../ElementaryProcesses/NeutralGas.h"
#include "../Grid/Grid.h"

void test_MonteCarloCollision() {
    double n = 1e-18, mAr = 6.6335209e-26, T = 500;
    NeutralGas gas(n, mAr, T);
    double sigma=1e-20, dt=1e-12, mass=9.1e-31, charge=1.6e-19;
    size_t Ntot=1e3, Nz=100, Nr=50;
    double dz=2e-5, dr=2e-5;
    Grid grid(Nz, Nr, dz, dr);
    Particles ptcls(mass, charge, Ntot, grid);
    ElectronNeutralElasticCollisions collisions(sigma, dt, gas, ptcls);
    //collisions.particles.rho.print();
}

#endif //CPP_RZ_PIC_TEST_MONTECARLOCOLLISION_H
