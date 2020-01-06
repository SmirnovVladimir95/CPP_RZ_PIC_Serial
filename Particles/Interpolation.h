//
// Created by Vladimir Smirnov on 06.11.2019.
//

#ifndef CPP_RZ_PIC_INTERPOLATION_H
#define CPP_RZ_PIC_INTERPOLATION_H

#include "../Tools/Matrix.h"
#include "../Grid/Grid.h"

void LinearFieldInterpolation(scalar efz[], scalar efr[], const scalar z[], const scalar r[],
                              const scalar Ez[], const scalar Er[], const Grid& grid, const size_t Ntot);

void LinearChargeInterpolation(scalar rho[], const scalar z[], const scalar r[], const Grid& grid,
                               const scalar charge, const size_t Ntot, const scalar node_volume[]);


#endif //CPP_RZ_PIC_INTERPOLATION_H
