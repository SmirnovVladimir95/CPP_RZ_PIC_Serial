#ifndef CPP_RZ_PIC_NEUTRALGAS_H
#define CPP_RZ_PIC_NEUTRALGAS_H

#include <array>
#include "../Tools/ProjectTypes.h"


class NeutralGas {
public:
    const scalar n;
    const scalar mass;
    const scalar T;
    NeutralGas(scalar n, scalar mass, scalar T);
    std::array<scalar, 3> generate_velocity() const;
};


#endif //CPP_RZ_PIC_NEUTRALGAS_H
