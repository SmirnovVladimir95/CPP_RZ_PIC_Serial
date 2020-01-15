#ifndef CPP_RZ_PIC_INTERPOLATION_H
#define CPP_RZ_PIC_INTERPOLATION_H

#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"

void LinearFieldInterpolation(scalar efz[], scalar efr[], const scalar z[], const scalar r[],
                              const scalar Ez[], const scalar Er[], const Grid& grid, size_t Ntot);

void LinearChargeInterpolation(scalar rho[], const scalar z[], const scalar r[], const Grid& grid,
                               scalar charge, size_t Ntot, const scalar node_volume[]);

void InitVolume(Matrix& node_volume, const Grid& grid);


#endif //CPP_RZ_PIC_INTERPOLATION_H
