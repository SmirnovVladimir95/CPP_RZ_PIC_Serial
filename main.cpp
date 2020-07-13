#include <fstream>
#include "Tests/Test_ChargeInterpolation.h"
#include "Tests/Test_Compute_E_from_Phi.h"
#include "Tests/Test_FieldInterpolation.h"
#include "Tests/Test_FieldSolver.h"
#include "Tests/Test_PIC_Cycle.h"
#include "Tests/Test_Pusher.h"
#include "Tests/Test_MonteCarloCollision.h"
#include "Tests/Test_NeutralGas.h"
#include "Tests/Test_ParticleLeave.h"
#include "Tests/Test_ParticleEmission.h"
#include "Tests/Test_Simulation.h"
#include "Tests/Test_ElementaryProcesses.h"
#include "Tests/Test_Logger.h"
#include "Tests/Test_Matrix.h"
#include "Tests/Test_SinglePtclSimulation.h"
#include "Tests/Test_EnergyCrossSection.h"
#include "Tests/Test_CollisionEnergyCrossSection.h"
#include "Tests/Test_ParticlesLoad.h"
#include "Tests/Test_SimulationEnergyCrossSection.h"
#include "Tests/Test_ReloadSimulation.h"
#include "Tests/Test_SimulationWithDielectric.h"
using namespace std;


int main() {
    //test_Field_interpolation();
    //test_Pusher();
    //test_single_particle_motion();
    //test_PIC_Cycle();
    //test_compute_E();
    //test_PIC_Cycle();
    //test_Field_solver();
    //test_Charge_interpolation();
    //test_Collisions();
    //test_NeutralGas();
    //test_ParticleLeave();
    //test_ParticleEmission();
    //test_Simulation();
    //test_SinglePtclSimulation();
    //test_MonteCarloCollision();
    //test_NanbuCollisionChoice();
    //test_MonteCarloCollisions();
    //test_ElementaryProcesses();
    //test_Logger();
    //test_Matrix();
    //test_EnergyCrossSection();
    //test_ParticleLoad();
    //test_CollisionEnergyCrossSection();
    //test_SimulationEnergyCrossSection();
    //test_ReloadSimulation();
    test_SimulationWithDielectric();
    return 0;
}