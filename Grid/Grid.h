#ifndef CPP_RZ_PIC_GRID_H
#define CPP_RZ_PIC_GRID_H

#include <iostream>
#include "../Tools/ProjectTypes.h"


struct Grid {
    size_t Nz, Nr;
    scalar dz, dr;
    Grid();
    Grid(size_t Nz, size_t Nr, scalar dz, scalar  dr);
};


#endif //CPP_RZ_PIC_GRID_H
