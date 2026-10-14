#pragma once

#include "csr_matrix.h"

#include <vector>

// Perform a specified number of Gauss-Seidel smoothing sweeps
// on the linear system:
//
//      A * x = b
//
// x is updated in place.
void gaussSeidel(
    const CSRMatrix& A,
    std::vector<double>& x,
    const std::vector<double>& b,
    int sweeps
);
