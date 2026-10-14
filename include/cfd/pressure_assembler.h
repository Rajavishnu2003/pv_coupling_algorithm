#pragma once
#include "cfd/momentum_assembler.h"
namespace cfd {
    class PressureAssembler {
        public:
        CSRMatrix A;
        std::vector<double> b;
        explicit PressureAssembler(const StaggeredGrid& g);
        void assembleMatrix(const MomentumAssembler& momentum);
        void assembleRHS(const FlowFields& predicted);
        private: const StaggeredGrid& g_;
    };
    // Net volume flux outward from each pressure cell (rho=1).
    std::vector<double> massImbalance(const StaggeredGrid& g,const FlowFields& f);
    void correctPressureVelocity(const StaggeredGrid& g,const MomentumAssembler& m,
    const std::vector<double>& pp,double alpha_p,FlowFields& f);
}
