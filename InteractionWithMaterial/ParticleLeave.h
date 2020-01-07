#ifndef CPP_RZ_PIC_PARTICLELEAVE_H
#define CPP_RZ_PIC_PARTICLELEAVE_H

#include "../Grid/Grid.h"
#include "../Particles/Particles.h"

class ParticleLeave {
private:
    Grid* grid;
    Particles* particles;
    Matrix* domain_condition;
    bool out_of_domain(int ptcl_idx);
public:
    ParticleLeave(Particles& particles, Grid& grid, Matrix& domain_condition);
    void leave();
};


#endif //CPP_RZ_PIC_PARTICLELEAVE_H
