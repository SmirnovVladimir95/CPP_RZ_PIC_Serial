#ifndef CPP_RZ_PIC_NANBUCOLLISIONS_H
#define CPP_RZ_PIC_NANBUCOLLISIONS_H

#include "Collision.h"
#include <ctime>
#include <random>
#include <iostream>

void NanbuElectronCollisionProcess(ElectronNeutralElasticCollision& collision1, Ionization& collision2,
        int collisions_num);

void NanbuIonCollisionProcess(IonNeutralElasticCollision& collision1, int collisions_num);

void NanbuIonCollisionProcess(IonNeutralElasticCollision& collision1, Ionization& collision2, int collisions_num);

int NanbuCollisionChoice(vector<scalar> collision_prob, int collisions_num);

#endif //CPP_RZ_PIC_NANBUCOLLISIONS_H
