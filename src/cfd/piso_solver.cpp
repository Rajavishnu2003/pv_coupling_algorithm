#include "cfd/coupling_solver.h"
namespace cfd {
    void CouplingSolver::transientOuterIteration(const FlowFields& old_time,double dt,double alpha) {
        const auto assembly_pressure=f_.p;
        momentumPredictor(old_time,dt,alpha);
        pressure_.assembleMatrix(momentum_);
        PressureLinearSolver solver(pressure_.A,c_,g_.nx,g_.ny);
        summary_.pressure_setup_seconds+=solver.setupSeconds();
        for(int k=0;k<c_.pressure_correctors;++k) {
            if(k>0)updateMomentumH(g_,momentum_,assembly_pressure,f_);
            pressureCorrect(solver,1.0);
        }
    }
    void CouplingSolver::pisoTimeStep(const FlowFields& old_time,double dt) {
        // One momentum assembly/solve per step; >=2 neighbour/pressure corrections.
        // PISO ignores SIMPLE relaxation: alpha_u=alpha_p=1.
        transientOuterIteration(old_time,dt,1.0);
    }
}
