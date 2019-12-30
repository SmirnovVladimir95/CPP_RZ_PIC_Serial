#ifndef CPP_RZ_PIC_NEUTRALGAS_H
#define CPP_RZ_PIC_NEUTRALGAS_H

#include <array>

using type_double = double;

class NeutralGas {
public:
    const type_double n;
    const type_double mass;
    const type_double T;
    NeutralGas(type_double n, type_double mass, type_double T);
    std::array<type_double, 3> generate_velocity() const;
};


#endif //CPP_RZ_PIC_NEUTRALGAS_H
