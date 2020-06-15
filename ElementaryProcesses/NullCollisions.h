#ifndef CPP_RZ_PIC_NULLCOLLISIONS_H
#define CPP_RZ_PIC_NULLCOLLISIONS_H

#include "Collision.h"


void NullElectronCollisions(ElectronNeutralElasticCollision& collision1, Ionization& collision2,
                            int collisions_num, scalar max_probability);

void NullIonCollisions(IonNeutralElasticCollision& collision1, int collisions_num, scalar max_probability);

#endif //CPP_RZ_PIC_NULLCOLLISIONS_H
