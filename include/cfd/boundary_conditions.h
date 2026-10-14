#pragma once
#include "cfd/flow_fields.h"
namespace cfd {
    struct CavityBoundary {
        double lid_speed=1.0;
        void enforceNormalVelocities(const StaggeredGrid& g,FlowFields& f) const;
        // Tangential Dirichlet values are eliminated at half-cell distance by the assembler.
    };
}
