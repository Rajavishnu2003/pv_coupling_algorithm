#pragma once
#include "cfd/coupling_solver.h"
namespace cfd {
    void writeResults(const StaggeredGrid& g,const FlowFields& f,const Config& c,const RunSummary& s,
    const MomentumAssembler& m,const PressureAssembler& p);
}
