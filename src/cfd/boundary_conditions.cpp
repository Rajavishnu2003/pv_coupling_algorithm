#include "cfd/boundary_conditions.h"
namespace cfd {
    void CavityBoundary::enforceNormalVelocities(const StaggeredGrid& g,FlowFields& f) const {
        for(int j=0;j<g.ny;++j) {
            f.u[g.uf(0,j)]=0;
            f.u[g.uf(g.nx,j)]=0;
        }
        for(int i=0;i<g.nx;++i) {
            f.v[g.vf(i,0)]=0;
            f.v[g.vf(i,g.ny)]=0;
        }
    }
}
