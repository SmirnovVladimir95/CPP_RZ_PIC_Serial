#ifndef CPP_RZ_PIC_PARTICLESLOGGER_H
#define CPP_RZ_PIC_PARTICLESLOGGER_H


#include <fstream>
#include "ProjectTypes.h"
#include "../Particles/Particles.h"
#include "Matrix.h"

class ParticlesLogger {
private:
    Particles* particles;
    string observerID;
    string mean_energy_file, velocity_file, position_file, n_total_file;
    ofstream mean_energy_output, velocity_output, position_output, n_total_output;
    bool iter_check(int iter, int step);
public:
    ParticlesLogger(Particles& particles, const string& observerID);
    ~ParticlesLogger();
    void mean_energy_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void velocity_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void position_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void n_total_log(int iter = 0, int step = 1, unsigned int mode = ios::app);
    void log(int iter = 0, int step = 1, unsigned int mode = ios::app);
};


#endif //CPP_RZ_PIC_PARTICLESLOGGER_H
