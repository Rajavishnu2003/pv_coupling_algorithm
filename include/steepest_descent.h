#pragma once

#include "csr_matrix.h"
#include "linear_algebra.h"
#include "preconditioner.h"
#include <vector>

class steepest_descent_solver{

    public : 

    // Convergence criterion: ||r|| / ||b|| <= tolerance
    double tolerance = 1e-8;

    int maxIterations = 1000000;

    // Keep the original console output unless a comparison driver disables it.
    bool verbose = true;


    double computeAlpha(const CSRMatrix& A,
        const std::vector<double>& r);

    // PSD: alpha = (r^T z) / (z^T A z), where M z = r.
    double computeAlpha(const CSRMatrix& A,
        const std::vector<double>& r,
        const std::vector<double>& z);

    void updatePhi( std::vector<double>& x, 
         double alpha, const std::vector<double>& r);    

    void solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x);

    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x,
        const Preconditioner& preconditioner);

    int iterations() const { return iters_; }


private:
    int iters_ = 0;

};

