#pragma once
#include "cfd/convergence_monitor.h"
#include "cfd/pressure_assembler.h"
#include "cfd/linear_solver_adapter.h"
#include <fstream>
#include <cstdint>
namespace cfd {
    struct RunSummary {
        bool converged=false,completed=false;
        int steps=0;
        double time=0,wall_seconds=0,momentum_solve_seconds=0,pressure_solve_seconds=0,pressure_setup_seconds=0;
        std::int64_t momentum_iterations=0,pressure_iterations=0,pressure_solves=0,momentum_solves=0;
        FlowResiduals residuals;
    };
    // Shared discretization and linear solvers; separate coupling controllers below.
    class CouplingSolver {
        public:
        explicit CouplingSolver(Config config,bool write_output=true);
        RunSummary run();
        const FlowFields& fields() const {
            return f_;
        }
        const StaggeredGrid& grid() const {
            return g_;
        }
        private:
        Config c_;
        StaggeredGrid g_;
        CavityBoundary bc_;
        FlowFields f_;
        MomentumAssembler momentum_,residual_workspace_;
        PressureAssembler pressure_;
        RunSummary summary_;
        bool output_;
        void momentumPredictor(const FlowFields& old_time,double dt,double alpha);
        void pressureCorrect(PressureLinearSolver& solver,double alpha_p);
        void simpleIteration(const FlowFields& old_iter);
        void pisoTimeStep(const FlowFields& old_time,double dt);
        void pimpleTimeStep(const FlowFields& old_time,double dt);
        void transientOuterIteration(const FlowFields& old_time,double dt,double alpha);
    };
}
