#include "cfd/coupling_solver.h"
#include <iostream>
#include <exception>
int main(int argc,char** argv) {
    if(argc!=2) {
        std::cerr<<"Usage: cavity path/to/case.cfg\n";
        return 1;
    }
    try {
        auto c=cfd::readConfig(argv[1]);
        std::cout<<"CSR cavity: "<<c.algorithm<<", "<<c.nx<<"x"<<c.ny<<", Re="<<c.reynolds<<", convection="<<c.convection<<'\n';
        cfd::CouplingSolver solver(c);
        auto s=solver.run();
        std::cout<<"completed="<<s.completed<<" steady_converged="<<s.converged<<" steps="<<s.steps<<" wall_seconds="<<s.wall_seconds<<'\n';
        return s.completed?0:2;
        // An iteration/time-step cap is NOT a convergence claim.
    }
    catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n';
        return 1;
    }
}
