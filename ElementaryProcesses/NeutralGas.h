#ifndef CPP_RZ_PIC_NEUTRALGAS_H
#define CPP_RZ_PIC_NEUTRALGAS_H

#include <array>

using type_double = double;

class NeutralGas {
private:
    const type_double n;
    const type_double mass;
    const type_double T;
public:
    NeutralGas(const type_double n, const type_double mass, const type_double T);
    std::array<type_double, 3> generate_velocity(const int seed=time(nullptr)) const;
};


#endif //CPP_RZ_PIC_NEUTRALGAS_H
