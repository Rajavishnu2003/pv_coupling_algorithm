#include "cfd/coupling_solver.h"
namespace cfd {
    void CouplingSolver::pimpleTimeStep(const FlowFields& old_time,double dt) {
        // Picard outer loops update convection and pressure while preserving time-n data.
        for(int outer=0;outer<c_.outer_correctors;++outer)
        transientOuterIteration(old_time,dt,c_.alpha_u);
    }
}
