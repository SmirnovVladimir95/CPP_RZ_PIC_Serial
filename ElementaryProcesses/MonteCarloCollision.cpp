#include "MonteCarloCollision.h"

MonteCarloCollision::MonteCarloCollision(vector<Collision> collision_list) : collision_list(collision_list) {}

int MonteCarloCollision::nanbu_method() {
    /*U = np.random.uniform(0, 1)
    i = int(np.floor(U*len(self.collision_list)))
    if (U > float(i)/len(self.collision_list) - self.get_collision_prob(ptcl_vel, i)) {
        return i;
    }
    return None*/
    return 0;
}