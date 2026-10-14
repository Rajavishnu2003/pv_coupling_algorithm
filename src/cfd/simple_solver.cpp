#include "cfd/coupling_solver.h"
namespace cfd {
    void CouplingSolver::simpleIteration(const FlowFields& previous) {
        momentumPredictor(previous,0.0,c_.alpha_u);
        pressure_.assembleMatrix(momentum_);
        PressureLinearSolver solver(pressure_.A,c_,g_.nx,g_.ny);
        summary_.pressure_setup_seconds+=solver.setupSeconds();
        pressureCorrect(solver,c_.alpha_p);
    }
}
