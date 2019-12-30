#include "NanbuCollisions.h"
#include <random>


void NanbuElectronCollisionProcess(ElectronNeutralElasticCollision& collision1, Ionization& collision2,
                                   int collisions_num) {
    int collision_type, Ntot;
    vector<type_double> prob;
    Ntot = collision1.particles->get_Ntot();
    for (int ptcl_idx = 0; ptcl_idx < Ntot; ptcl_idx++) {
        prob.push_back(collision1.probability(ptcl_idx));
        prob.push_back(collision2.probability(ptcl_idx));
        collision_type = NanbuCollisionChoice(prob, collisions_num);
        if (collision_type == 1)
            collision1.collision(ptcl_idx);
        else if (collision_type == 2)
            collision2.collision(ptcl_idx);
        prob.pop_back();
        prob.pop_back();
    }
}

void NanbuIonCollisionProcess(IonNeutralElasticCollision& collision1, int collisions_num) {
    int collision_type, Ntot;
    //vector<type_double> prob = {0};
    vector<type_double> prob;
    Ntot = collision1.particles->get_Ntot();
    for (int ptcl_idx = 0; ptcl_idx < Ntot; ptcl_idx++) {
        //prob[ptcl_idx] = collision1.probability(ptcl_idx);
        prob.push_back(collision1.probability(ptcl_idx));
        collision_type = NanbuCollisionChoice(prob, collisions_num);
        if (collision_type == 1)
            collision1.collision(ptcl_idx);
        prob.pop_back();
    }
}

int NanbuCollisionChoice(vector<type_double> collision_prob, int collisions_num) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    float U = distribution(generator);
    float i = floor(1+U*collisions_num);
    if (U > i/collisions_num - collision_prob[i-1])
        return i;
    return 0;
}