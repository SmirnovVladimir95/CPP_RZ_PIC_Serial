#ifndef CPP_RZ_PIC_PARTICLESLOGGER_H
#define CPP_RZ_PIC_PARTICLESLOGGER_H


#include "ProjectTypes.h"
#include "../Particles/Particles.h"
#include "Matrix.h"

class ParticlesLogger {
private:
    Particles* particles;
    string observerID;
    string mean_energy_file, velocity_file, position_file, n_total_file;
public:
    ParticlesLogger(Particles& particles, const string& observerID);
    void mean_energy_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void velocity_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void position_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void n_total_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    bool iter_check(int iter, int step);
    void log(int iter = 0, int step = 1, unsigned int mode = ios::app);
};


#endif //CPP_RZ_PIC_PARTICLESLOGGER_H
