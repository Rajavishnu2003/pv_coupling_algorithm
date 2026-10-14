#include "cfd/output_writer.h"
#include "cfd/csr_tools.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
namespace cfd {
    namespace {
        std::ofstream open(const std::string& p) {
            std::ofstream f(p);
            if(!f)throw std::runtime_error("Cannot write "+p);
            f<<std::setprecision(17);
            return f;
        }
    }
    void writeResults(const StaggeredGrid& g,const FlowFields& f,const Config& c,const RunSummary& s,const MomentumAssembler& m,const PressureAssembler& p) {
        auto field=open(c.output+"/fields.csv");
        field<<"x,y,u,v,p\n";
        auto vtk=open(c.output+"/cavity.vtk");
        vtk<<"# vtk DataFile Version 3.0\nCSR staggered cavity (cell-centred output)\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<g.nx<<' '<<g.ny<<" 1\nORIGIN "<<0.5*g.dx<<' '<<0.5*g.dy<<" 0\nSPACING "<<g.dx<<' '<<g.dy<<" 1\nPOINT_DATA "<<g.np()<<"\nVECTORS velocity double\n";
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            double u=0.5*(f.u[g.uf(i,j)]+f.u[g.uf(i+1,j)]),v=0.5*(f.v[g.vf(i,j)]+f.v[g.vf(i,j+1)]);
            field<<(i+0.5)*g.dx<<','<<(j+0.5)*g.dy<<','<<u<<','<<v<<','<<f.p[g.p(i,j)]<<'\n';
            vtk<<u<<' '<<v<<" 0\n";
        }
        vtk<<"SCALARS pressure double 1\nLOOKUP_TABLE default\n";
        for(double a:f.p)vtk<<a<<'\n';
        auto uf=open(c.output+"/u_faces.csv");
        uf<<"i,j,x,y,u\n";
        for(int j=0;j<g.ny;++j)for(int i=0;i<=g.nx;++i)uf<<i<<','<<j<<','<<i*g.dx<<','<<(j+.5)*g.dy<<','<<f.u[g.uf(i,j)]<<'\n';
        auto vf=open(c.output+"/v_faces.csv");
        vf<<"i,j,x,y,v\n";
        for(int j=0;j<=g.ny;++j)for(int i=0;i<g.nx;++i)vf<<i<<','<<j<<','<<(i+.5)*g.dx<<','<<j*g.dy<<','<<f.v[g.vf(i,j)]<<'\n';
        auto uc=open(c.output+"/centerline_u.csv");
        uc<<"y,u\n0,0\n";
        int il=g.nx/2;
        double ax=0.5*g.nx-il;
        for(int j=0;j<g.ny;++j)uc<<(j+.5)*g.dy<<','<<(1-ax)*f.u[g.uf(il,j)]+ax*f.u[g.uf(il+1,j)]<<'\n';
        uc<<"1,1\n";
        auto vc=open(c.output+"/centerline_v.csv");
        vc<<"x,v\n0,0\n";
        int jl=g.ny/2;
        double ay=0.5*g.ny-jl;
        for(int i=0;i<g.nx;++i)vc<<(i+.5)*g.dx<<','<<(1-ay)*f.v[g.vf(i,jl)]+ay*f.v[g.vf(i,jl+1)]<<'\n';
        vc<<"1,0\n";
        auto json=open(c.output+"/summary.json");
        json<<std::boolalpha<<"{\n";
        json<<"  \"algorithm\": \""<<c.algorithm<<"\",\n  \"pressure_solver\": \""<<c.pressure_solver<<"\",\n  \"convection\": \""<<c.convection<<"\",\n";
#define J(key,val) json<<"  \""<<key<<"\": "<<(val)<<",\n";
        J("nx",g.nx) J("ny",g.ny) J("Re",c.reynolds) J("dt",c.algorithm=="simple"?0:c.dt)
        J("converged",s.converged) J("completed",s.completed) J("steps",s.steps) J("time",s.time)
        J("wall_seconds",s.wall_seconds) J("momentum_solve_seconds",s.momentum_solve_seconds) J("pressure_solve_seconds",s.pressure_solve_seconds) J("pressure_setup_seconds",s.pressure_setup_seconds)
        J("momentum_iterations",s.momentum_iterations) J("pressure_iterations",s.pressure_iterations) J("pressure_solves",s.pressure_solves) J("momentum_solves",s.momentum_solves)
        J("u_momentum_inf",s.residuals.u_momentum) J("v_momentum_inf",s.residuals.v_momentum) J("divergence_inf",s.residuals.divergence)
        J("velocity_change",s.residuals.velocity_change) J("courant",s.residuals.courant)
        J("cds_faces",s.residuals.schemes.cds_faces) J("upwind_faces",s.residuals.schemes.upwind_faces) J("max_abs_peclet",s.residuals.schemes.max_abs_peclet)
        J("u_unknowns",g.nu()) J("v_unknowns",g.nv()) J("p_unknowns",g.np()) J("u_nnz",m.u.A.nnz()) J("v_nnz",m.v.A.nnz())
#undef J
        json<<"  \"p_nnz\": "<<p.A.nnz()<<"\n}\n";
        if(c.dump_matrices) {
            dumpCSR(m.u.A,m.u.b,c.output+"/u_predictor");
            dumpCSR(m.v.A,m.v.b,c.output+"/v_predictor");
            dumpCSR(p.A,p.b,c.output+"/pressure_correction");
        }
    }
}
