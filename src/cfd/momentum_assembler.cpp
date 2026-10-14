#include "cfd/momentum_assembler.h"
#include "cfd/csr_tools.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd {
    MomentumSystem::MomentumSystem(int nx,int ny):A(fivePointPattern(nx,ny)),b(nx*ny,0),d(nx*ny,0){}
    MomentumAssembler::MomentumAssembler(const StaggeredGrid& g):u(g.nx-1,g.ny),v(g.nx,g.ny-1),g_(g){}
    namespace {
        void reset(MomentumSystem& s) {
            zeroValues(s.A);
            std::fill(s.b.begin(),s.b.end(),0);
            s.stats={};
        }
        void face(MomentumSystem& s,int r,int neighbour,double F,double D,double known,const std::string& scheme) {
            const auto a=faceCoefficients(F,D,scheme);
            s.A.values[s.A.diagonalIndex[r]]+=a.diagonal;
            if(neighbour>=0)addCoefficient(s.A,r,neighbour,a.neighbour);
            else s.b[r]-=a.neighbour*known;
            if(a.cds)++s.stats.cds_faces;
            else ++s.stats.upwind_faces;
            s.stats.max_abs_peclet=std::max(s.stats.max_abs_peclet,std::abs(a.peclet));
        }
        void finish(MomentumSystem& s,int r,double area,double mass_dt,double old_time,double old_iter,double alpha) {
            auto& diag=s.A.values[s.A.diagonalIndex[r]];
            diag+=mass_dt;
            s.b[r]+=mass_dt*old_time;
            const double unrelaxed=diag;
            if(!(unrelaxed>0))throw std::runtime_error("Nonpositive momentum diagonal. Check CDS/cell Peclet and relaxation.");
            diag=unrelaxed/alpha;
            s.b[r]+=(1.0-alpha)*diag*old_iter;
            s.d[r]=area/diag;
            // Includes alpha ONCE via the relaxed diagonal.
        }
    }
    void MomentumAssembler::assemble(const FlowFields& f,const FlowFields& old,
    const CavityBoundary& bc,const Config& c,double dt,double alpha) {
        if(alpha<=0||alpha>1||dt<0)throw std::invalid_argument("Invalid assembly controls.");
        const auto& g=g_;
        const double mu=1.0/c.reynolds;
        const double Dx=mu*g.dy/g.dx,Dy=mu*g.dx/g.dy;
        const double mass_dt=dt>0?g.volume()/dt:0;
        reset(u);
        reset(v);
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i) {
            int r=g.u(i,j);
            double E=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)])*g.dy;
            double W=-0.5*(f.u[g.uf(i-1,j)]+f.u[g.uf(i,j)])*g.dy;
            double N=0.5*(f.v[g.vf(i-1,j+1)]+f.v[g.vf(i,j+1)])*g.dx;
            double S=-0.5*(f.v[g.vf(i-1,j)]+f.v[g.vf(i,j)])*g.dx;
            face(u,r,i<g.nx-1?g.u(i+1,j):-1,E,Dx,0,c.convection);
            face(u,r,i>1?g.u(i-1,j):-1,W,Dx,0,c.convection);
            face(u,r,j<g.ny-1?g.u(i,j+1):-1,N,j==g.ny-1?2*Dy:Dy,j==g.ny-1?bc.lid_speed:0,c.convection);
            face(u,r,j>0?g.u(i,j-1):-1,S,j==0?2*Dy:Dy,0,c.convection);
            u.b[r]+=(f.p[g.p(i-1,j)]-f.p[g.p(i,j)])*g.dy;
            finish(u,r,g.dy,mass_dt,old.u[g.uf(i,j)],f.u[g.uf(i,j)],alpha);
        }
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            int r=g.v(i,j);
            double E=0.5*(f.u[g.uf(i+1,j-1)]+f.u[g.uf(i+1,j)])*g.dy;
            double W=-0.5*(f.u[g.uf(i,j-1)]+f.u[g.uf(i,j)])*g.dy;
            double N=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)])*g.dx;
            double S=-0.5*(f.v[g.vf(i,j-1)]+f.v[g.vf(i,j)])*g.dx;
            face(v,r,i<g.nx-1?g.v(i+1,j):-1,E,i==g.nx-1?2*Dx:Dx,0,c.convection);
            face(v,r,i>0?g.v(i-1,j):-1,W,i==0?2*Dx:Dx,0,c.convection);
            face(v,r,j<g.ny-1?g.v(i,j+1):-1,N,Dy,0,c.convection);
            face(v,r,j>1?g.v(i,j-1):-1,S,Dy,0,c.convection);
            v.b[r]+=(f.p[g.p(i,j-1)]-f.p[g.p(i,j)])*g.dx;
            finish(v,r,g.dx,mass_dt,old.v[g.vf(i,j)],f.v[g.vf(i,j)],alpha);
        }
    }
    void updateMomentumH(const StaggeredGrid& g,const MomentumAssembler& m,const std::vector<double>& p0,FlowFields& f) {
        auto ux=packU(g,f),vx=packV(g,f),newu=ux,newv=vx;
        auto off=[](const CSRMatrix& A,const std::vector<double>& x,int r) {
            double s=0;
            for(int k=A.rowPtr[r];k<A.rowPtr[r+1];++k)if(A.colIndex[k]!=r)s+=A.values[k]*x[A.colIndex[k]];
            return s;
        };
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i) {
            int r=g.u(i,j),l=g.p(i-1,j),e=g.p(i,j);
            double dp=(f.p[l]-f.p[e])-(p0[l]-p0[e]);
            newu[r]=(m.u.b[r]+g.dy*dp-off(m.u.A,ux,r))/m.u.A.values[m.u.A.diagonalIndex[r]];
        }
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            int r=g.v(i,j),s=g.p(i,j-1),n=g.p(i,j);
            double dp=(f.p[s]-f.p[n])-(p0[s]-p0[n]);
            newv[r]=(m.v.b[r]+g.dx*dp-off(m.v.A,vx,r))/m.v.A.values[m.v.A.diagonalIndex[r]];
        }
        unpackU(g,newu,f);
        unpackV(g,newv,f);
    }
}
