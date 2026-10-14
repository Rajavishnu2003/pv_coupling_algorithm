#include "cfd/convergence_monitor.h"
#include "cfd/pressure_assembler.h"
#include "cfd/csr_tools.h"
#include <cmath>
#include <algorithm>
namespace cfd {
    FlowResiduals evaluateResiduals(const StaggeredGrid& g,const FlowFields& f,const FlowFields& old,const CavityBoundary& bc,const Config& c,MomentumAssembler& m,double dt) {
        // Reassemble the UNRELAXED STEADY operator, not the just-solved relaxed predictor.
        m.assemble(f,f,bc,c,0.0,1.0);
        FlowResiduals r;
        r.u_momentum=normInf(computeResidual(m.u.A,packU(g,f),m.u.b))/g.volume();
        r.v_momentum=normInf(computeResidual(m.v.A,packV(g,f),m.v.b))/g.volume();
        r.divergence=normInf(massImbalance(g,f))/g.volume();
        for(std::size_t i=0;i<f.u.size();++i)r.velocity_change=std::max(r.velocity_change,std::abs(f.u[i]-old.u[i]));
        for(std::size_t i=0;i<f.v.size();++i)r.velocity_change=std::max(r.velocity_change,std::abs(f.v[i]-old.v[i]));
        if(dt>0)for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            const double co=dt*(0.5*(std::abs(f.u[g.uf(i,j)])+std::abs(f.u[g.uf(i+1,j)]))/g.dx+
            0.5*(std::abs(f.v[g.vf(i,j)])+std::abs(f.v[g.vf(i,j+1)]))/g.dy);
            r.courant=std::max(r.courant,co);
        }
        r.schemes=m.u.stats;
        r.schemes.merge(m.v.stats);
        return r;
    }
    bool steadyConverged(const FlowResiduals& r,const Config& c,double dt) {
        const double change=dt>0?r.velocity_change/dt:r.velocity_change;
        return std::max(r.u_momentum,r.v_momentum)<=c.momentum_tolerance&&r.divergence<=c.continuity_tolerance&&change<=c.change_tolerance;
    }
}
