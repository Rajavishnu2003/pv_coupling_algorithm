#pragma once
namespace cfd {
    struct StaggeredGrid {
        int nx,ny;
        double dx,dy;
        explicit StaggeredGrid(int nx_,int ny_);
        int np() const {
            return nx*ny;
        }
        int nu() const {
            return (nx-1)*ny;
        }
        int nv() const {
            return nx*(ny-1);
        }
        // Pressure: i=0..nx-1, j=0..ny-1; physical coordinates ((i+.5)dx,(j+.5)dy).
        int p(int i,int j) const {
            return j*nx+i;
        }
        // Full staggered arrays include prescribed normal velocities on walls.
        int uf(int i,int j) const {
            return j*(nx+1)+i;
        }
        int vf(int i,int j) const {
            return j*nx+i;
        }
        // Linear-system unknowns exclude prescribed normal wall velocities.
        int u(int i,int j) const {
            return j*(nx-1)+(i-1);
        }
        // i=1..nx-1, j=0..ny-1
        int v(int i,int j) const {
            return (j-1)*nx+i;
        }
        // i=0..nx-1, j=1..ny-1
        double volume() const {
            return dx*dy;
        }
    };
}
