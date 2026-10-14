#pragma once
#include <string>
#include <cstdint>
namespace cfd {
    struct SchemeStatistics {
        std::uint64_t cds_faces=0,upwind_faces=0;
        double max_abs_peclet=0;
        void merge(const SchemeStatistics& s);
    };
    struct FaceCoefficients {
        double diagonal,neighbour;
        bool cds;
        double peclet;
    };
    // F is SIGNED OUTWARD mass flux of this control volume. D=mu*S/distance.
    // Diffusion is retained in full for BOTH convection branches.
    FaceCoefficients faceCoefficients(double F,double D,const std::string& scheme="peclet");
}
