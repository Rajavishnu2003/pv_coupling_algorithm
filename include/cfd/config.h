#pragma once
#include <string>
namespace cfd {
    // Nondimensional unit cavity: rho=1, U_lid=1, L=1, nu=1/Re.
    struct Config {
        int nx=16, ny=16;
        double reynolds=100.0;
        std::string algorithm="simple";
        // simple, piso, pimple
        std::string pressure_solver="cg_ic0";
        // cg_ic0, cg, sd_ic0, bicgstab, multigrid
        std::string convection="peclet";
        // peclet, cds, upwind (latter two for comparison)
        double alpha_u=0.7, alpha_p=0.3;
        double dt=0.01;
        int pressure_correctors=2, outer_correctors=3;
        int max_steps=20000, max_linear_iterations=3000;
        double momentum_tolerance=1e-6, continuity_tolerance=1e-8, change_tolerance=1e-8;
        double linear_rtol=1e-10, linear_atol=1e-13;
        bool stop_at_steady=true;
        double end_time=1.0;
        std::string mg_cycle="V";
        int mg_pre=2, mg_post=2;
        bool mg_fmg_start=false;
        int log_every=100;
        std::string output="results/cavity";
        bool dump_matrices=false;
        void validate() const;
    };
    Config readConfig(const std::string& path);
    void writeConfig(const Config& c, const std::string& path);
}
// namespace cfd
