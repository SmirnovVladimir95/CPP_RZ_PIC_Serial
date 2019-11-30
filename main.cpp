#include "Tests/Test_ChargeInterpolation.h"
#include "Tests/Test_Compute_E_from_Phi.h"
#include "Tests/Test_FieldInterpolation.h"
#include "Tests/Test_FieldSolver.h"
#include "Tests/Test_PIC_Cycle.h"
#include "Tests/Test_Pusher.h"
#include "Tests/Test_MonteCarloCollision.h"
#include "Tests/Test_NeutralGas.h"
using namespace std;

int main() {
    //test_Field_interpolation();
    //test_Pusher();
    //test_PIC_Cycle();
    //test_compute_E();
    //test_PIC_Cycle();
    //test_Field_solver();
    //test_Charge_interpolation();
    //test_MonteCarloCollision();
    test_NeutralGas();
    return 0;
}