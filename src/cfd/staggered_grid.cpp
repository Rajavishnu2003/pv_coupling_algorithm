#include "cfd/staggered_grid.h"
#include <stdexcept>
namespace cfd {
    StaggeredGrid::StaggeredGrid(int x,int y):nx(x),ny(y),dx(1.0/x),dy(1.0/y) {
        if(x<3||y<3||x>2048||y>2048) throw std::invalid_argument("Unsupported grid size.");
    }
}
