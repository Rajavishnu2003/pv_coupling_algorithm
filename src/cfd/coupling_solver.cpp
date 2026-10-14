#include "cfd/coupling_solver.h"
#include "cfd/output_writer.h"
#include <filesystem>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
namespace cfd {
    namespace {
        Config validated(Config c) {
            c.validate();
            return c;
        }
    }
    CouplingSolver::CouplingSolver(Config c,bool output):c_(validated(std::move(c))),g_(c_.nx,c_.ny),f_(g_),momentum_(g_),residual_workspace_(g_),pressure_(g_),output_(output) {
        c_.validate();
        bc_.enforceNormalVelocities(g_,f_);
    }
    void CouplingSolver::momentumPredictor(const FlowFields& old_time,double dt,double alpha) {
        momentum_.assemble(f_,old_time,bc_,c_,dt,alpha);
        auto u=packU(g_,f_),v=packV(g_,f_);
        auto ru=solveMomentum(momentum_.u.A,momentum_.u.b,u,c_);
        auto rv=solveMomentum(momentum_.v.A,momentum_.v.b,v,c_);
        unpackU(g_,u,f_);
        unpackV(g_,v,f_);
        bc_.enforceNormalVelocities(g_,f_);
        summary_.momentum_iterations+=ru.iterations+rv.iterations;
        summary_.momentum_solves+=2;
        summary_.momentum_solve_seconds+=ru.seconds+rv.seconds;
    }
    void CouplingSolver::pressureCorrect(PressureLinearSolver& solver,double alpha_p) {
        pressure_.assembleRHS(f_);
        std::vector<double> pp(g_.np(),0);
        const auto r=solver.solve(pressure_.b,pp);
        summary_.pressure_iterations+=r.iterations;
        ++summary_.pressure_solves;
        summary_.pressure_solve_seconds+=r.seconds;
        correctPressureVelocity(g_,momentum_,pp,alpha_p,f_);
        bc_.enforceNormalVelocities(g_,f_);
    }
    RunSummary CouplingSolver::run() {
        const auto start=std::chrono::steady_clock::now();
        std::ofstream history;
        if(output_) {
            std::filesystem::create_directories(c_.output);
            writeConfig(c_,c_.output+"/input.cfg");
            history.open(c_.output+"/residuals.csv");
            if(!history)throw std::runtime_error("Cannot open convergence history.");
            history<<"step,time,u_momentum_inf,v_momentum_inf,divergence_inf,velocity_change,courant,cds_faces,upwind_faces,max_abs_peclet,momentum_iterations,pressure_iterations\n"<<std::setprecision(17);
        }
        for(int step=1;step<=c_.max_steps;++step) {
            const FlowFields previous=f_;
            double dt=c_.algorithm=="simple"?0.0:c_.dt;
            if(!c_.stop_at_steady)dt=std::min(dt,c_.end_time-summary_.time);
            if(c_.algorithm=="simple")simpleIteration(previous);
            else if(c_.algorithm=="piso")pisoTimeStep(previous,dt);
            else pimpleTimeStep(previous,dt);
            summary_.time+=dt;
            summary_.steps=step;
            summary_.residuals=evaluateResiduals(g_,f_,previous,bc_,c_,residual_workspace_,dt);
            const auto& r=summary_.residuals;
            if(history)history<<step<<','<<summary_.time<<','<<r.u_momentum<<','<<r.v_momentum<<','<<r.divergence<<','<<r.velocity_change<<','<<r.courant<<','<<r.schemes.cds_faces<<','<<r.schemes.upwind_faces<<','<<r.schemes.max_abs_peclet<<','<<summary_.momentum_iterations<<','<<summary_.pressure_iterations<<'\n';
            if(step>2&&steadyConverged(r,c_,dt))summary_.converged=true;
            if(output_&&(step==1||step%c_.log_every==0||summary_.converged))
            std::cout<<"step="<<step<<" t="<<summary_.time<<" mom="<<std::max(r.u_momentum,r.v_momentum)<<" div="<<r.divergence<<" change="<<r.velocity_change<<" upwind="<<r.schemes.upwind_faces<<'\n';
            if((c_.stop_at_steady&&summary_.converged)||(!c_.stop_at_steady&&summary_.time>=c_.end_time-1e-12)) {
                summary_.completed=true;
                break;
            }
            if(!c_.stop_at_steady)summary_.converged=false;
        }
        summary_.wall_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(output_)writeResults(g_,f_,c_,summary_,momentum_,pressure_);
        return summary_;
    }
}
