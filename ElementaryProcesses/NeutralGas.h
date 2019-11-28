#ifndef CPP_RZ_PIC_NEUTRALGAS_H
#define CPP_RZ_PIC_NEUTRALGAS_H

using type_double = double;

struct NeutralGas {
    type_double n;
    type_double mass;
    type_double T;
    void generate_velocity();
};


#endif //CPP_RZ_PIC_NEUTRALGAS_H
