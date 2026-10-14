#pragma once

#include "csr_matrix.h"
#include "linear_algebra.h"
#include "preconditioner.h"

#include <vector>

// BiCGSTAB solver for general (non-symmetric) systems A x = b.
// Convergence criterion: ||r|| <= tol * ||b||
class bicgstab_solver
{
public:

    explicit bicgstab_solver(double tol = 1e-8, int maxIter = 10000)
        : tol_(tol), maxIter_(maxIter)
    {
    }

    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x);

    // Fixed right preconditioning: matrix products use A * M^-1.
    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x,
        const Preconditioner& preconditioner);

    // Number of iterations performed by the last solve()
    int iterations() const { return iters_; }

private:

    double computeBeta(
        double rho,
        double rhoOld,
        double alpha,
        double omega) const;

    // p = r + beta * (p - omega * v)
    void updateDirection(
        const std::vector<double>& r,
        double beta,
        double omega,
        const std::vector<double>& v,
        std::vector<double>& p) const;

    double computeAlpha(
        double rho,
        const std::vector<double>& rHat,
        const std::vector<double>& v) const;

    // s = r - alpha * v
    void computeIntermediateResidual(
        const std::vector<double>& r,
        double alpha,
        const std::vector<double>& v,
        std::vector<double>& s) const;

    double computeOmega(
        const std::vector<double>& t,
        const std::vector<double>& s) const;

    // x = x + alpha * pPre + omega * sPre
    void updatePhi(
        std::vector<double>& x,
        double alpha,
        const std::vector<double>& pPre,
        double omega,
        const std::vector<double>& sPre) const;

    // r = s - omega * t
    void updateResidual(
        std::vector<double>& r,
        const std::vector<double>& s,
        double omega,
        const std::vector<double>& t) const;

    double tol_;
    int maxIter_;
    int iters_ = 0;
};
