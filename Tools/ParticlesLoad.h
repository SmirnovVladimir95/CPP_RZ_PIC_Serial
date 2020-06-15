#ifndef CPP_RZ_PIC_PARTICLESLOAD_H
#define CPP_RZ_PIC_PARTICLESLOAD_H

#include "../Particles/Particles.h"
#include <string>


class ParticlesLoad {
private:
    Particles* particles;
    string vel_path, pos_path;
    void position_load();
    void velocity_load();
public:
    ParticlesLoad(Particles& particles, const string& pos_path, const string& vel_path);
};

#endif //CPP_RZ_PIC_PARTICLESLOAD_H
