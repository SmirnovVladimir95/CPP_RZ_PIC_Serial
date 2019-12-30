#ifndef CPP_RZ_PIC_MONTECARLOCOLLISION_H
#define CPP_RZ_PIC_MONTECARLOCOLLISION_H

#include <vector>
#include "Collision.h"

class MonteCarloCollision {
public:
    vector<Collision> collision_list;
    MonteCarloCollision(vector<Collision> collision_list);
    int nanbu_method();
    void collisions();
};


#endif //CPP_RZ_PIC_MONTECARLOCOLLISION_H
