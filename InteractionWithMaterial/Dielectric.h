#ifndef CPP_RZ_PIC_DIELECTRIC_H
#define CPP_RZ_PIC_DIELECTRIC_H


#include "../Grid/Grid.h"
#include "../Particles/Particles.h"

class Dielectric {
private:
    //       UP
    // LEFT      RIGHT
    //      DOWN
    Grid* grid;
    Particles* incident_particles;
    Matrix* domain_condition;
    scalar dielectric_const;
    const scalar RIGHT = 1, DOWN = 2, LEFT = 3, UP = 4, None = 0;
public:
    Matrix rho;
    Dielectric(Particles& incident_particles, Grid& grid, Matrix& domain_condition, scalar dielectric_const);
    void update_rho();
};


#endif //CPP_RZ_PIC_DIELECTRIC_H
