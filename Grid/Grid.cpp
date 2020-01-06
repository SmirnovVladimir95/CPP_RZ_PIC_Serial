#include "Grid.h"


Grid::Grid() : Nz(0), Nr(0), dz(0), dr(0) {}
Grid::Grid(size_t Nz, size_t Nr, scalar dz, scalar dr) : Nz(Nz), Nr(Nr), dz(dz), dr(dr) {}