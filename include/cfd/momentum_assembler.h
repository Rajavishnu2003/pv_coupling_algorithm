#pragma once
#include "cfd/flow_fields.h"
#include "cfd/boundary_conditions.h"
#include "cfd/config.h"
#include "cfd/convection_scheme.h"
#include "csr_matrix.h"
namespace cfd {
    struct MomentumSystem {
        CSRMatrix A;
        std::vector<double> b,d;
        SchemeStatistics stats;
        MomentumSystem(int nx,int ny);
    };
    class MomentumAssembler {
        public:
        explicit MomentumAssembler(const StaggeredGrid& g);
        MomentumSystem u,v;
        // convecting: Picard-iteration state; old_time: fixed time-n state.
        // dt=0 -> steady; alpha=1 -> unrelaxed physical matrix.
        void assemble(const FlowFields& convecting,const FlowFields& old_time,
        const CavityBoundary& bc,const Config& c,double dt,double alpha);
        private:
        const StaggeredGrid& g_;
    };
    // Frozen-matrix neighbour correction used after the first PISO pressure correction.
    // This recalculates H(U), not a repetition of an unchanged continuity solve.
    void updateMomentumH(const StaggeredGrid& g,const MomentumAssembler& m,
    const std::vector<double>& assembly_pressure,FlowFields& f);
}
