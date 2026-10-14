#include "cfd/convection_scheme.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace cfd {
    void SchemeStatistics::merge(const SchemeStatistics& s) {
        cds_faces+=s.cds_faces;
        upwind_faces+=s.upwind_faces;
        max_abs_peclet=std::max(max_abs_peclet,s.max_abs_peclet);
    }
    FaceCoefficients faceCoefficients(double F,double D,const std::string& scheme) {
        if(!std::isfinite(F)||!std::isfinite(D)||D<=0)throw std::invalid_argument("Invalid face flux/conductance.");
        if(scheme!="peclet"&&scheme!="cds"&&scheme!="upwind")throw std::invalid_argument("Invalid scheme.");
        const double Pe=F/D;
        const bool cds=scheme=="cds"||(scheme=="peclet"&&std::abs(Pe)<=2.0);
        if(cds)return {
            D+0.5*F,-D+0.5*F,true,Pe
        };
        return {
            D+std::max(F,0.0),-D+std::min(F,0.0),false,Pe
        };
    }
}
