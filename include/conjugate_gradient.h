#pragma once

#include "csr_matrix.h"
#include "linear_algebra.h"
#include "preconditioner.h"

#include <vector>

// Conjugate Gradient solver for symmetric positive definite systems A x = b.
// Convergence criterion: ||r|| <= tol * ||b||
class conjugate_gradient_solver
{
public:

    explicit conjugate_gradient_solver(double tol = 1e-8, int maxIter = 10000)
        : tol_(tol), maxIter_(maxIter)
    {
    }

    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x);

    // Fixed SPD preconditioner: applyInverse(r) solves M z = r.
    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x,
        const Preconditioner& preconditioner);

    // Number of iterations performed by the last solve()
    int iterations() const { return iters_; }

private:

    double computeBeta(double rzNew, double rzOld) const;

    void updateDirection(
        const std::vector<double>& z,
        double beta,
        std::vector<double>& d) const;

    double computeAlpha(
        double rzOld,
        const std::vector<double>& d,
        const std::vector<double>& Ad) const;

    void updatePhi(
        std::vector<double>& x,
        double alpha,
        const std::vector<double>& d) const;

    void updateResidual(
        std::vector<double>& r,
        double alpha,
        const std::vector<double>& Ad) const;

    double tol_;
    int maxIter_;
    int iters_ = 0;
};
