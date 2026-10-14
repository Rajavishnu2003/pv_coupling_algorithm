#include "cfd/pressure_assembler.h"
#include "cfd/csr_tools.h"
#include <algorithm>
#include <stdexcept>
namespace cfd {
    PressureAssembler::PressureAssembler(const StaggeredGrid& g):A(fivePointPattern(g.nx,g.ny)),b(g.np(),0),g_(g){}
    void PressureAssembler::assembleMatrix(const MomentumAssembler& m) {
        zeroValues(A);
        const auto& g=g_;
        // Each pressure edge is added ONCE to both rows, guaranteeing symmetry.
        auto edge=[&](int p,int q,double a) {
            if(!(a>0))throw std::runtime_error("Pressure conductance must be positive.");
            A.values[A.diagonalIndex[p]]+=a;
            A.values[A.diagonalIndex[q]]+=a;
            addCoefficient(A,p,q,-a);
            addCoefficient(A,q,p,-a);
        };
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i)edge(g.p(i-1,j),g.p(i,j),g.dy*m.u.d[g.u(i,j)]);
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i)edge(g.p(i,j-1),g.p(i,j),g.dx*m.v.d[g.v(i,j)]);
        std::fill(b.begin(),b.end(),0);
        pinZeroReference(A,b,0);
    }
    std::vector<double> massImbalance(const StaggeredGrid& g,const FlowFields& f) {
        std::vector<double> b(g.np());
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i)
        b[g.p(i,j)]=(f.u[g.uf(i+1,j)]-f.u[g.uf(i,j)])*g.dy+(f.v[g.vf(i,j+1)]-f.v[g.vf(i,j)])*g.dx;
        return b;
    }
    void PressureAssembler::assembleRHS(const FlowFields& f) {
        b=massImbalance(g_,f);
        for(double& a:b)a=-a;
        b[0]=0;
    }
    void correctPressureVelocity(const StaggeredGrid& g,const MomentumAssembler& m,const std::vector<double>& pp,double alpha_p,FlowFields& f) {
        if(static_cast<int>(pp.size())!=g.np())throw std::invalid_argument("Wrong pressure correction size.");
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i)
        f.u[g.uf(i,j)]+=m.u.d[g.u(i,j)]*(pp[g.p(i-1,j)]-pp[g.p(i,j)]);
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i)
        f.v[g.vf(i,j)]+=m.v.d[g.v(i,j)]*(pp[g.p(i,j-1)]-pp[g.p(i,j)]);
        for(int k=0;k<g.np();++k)f.p[k]+=alpha_p*pp[k];
        const double ref=f.p[0];
        for(double& p:f.p)p-=ref;
    }
}
