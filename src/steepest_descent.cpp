#include "steepest_descent.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace
{
void validateSolve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double tolerance,
    int maxIterations)
{
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || maxIterations < 0)
    {
        throw std::invalid_argument(
            "steepest_descent_solver: tolerance must be finite and positive; max iterations must be nonnegative");
    }
    const std::size_t n = b.size();
    if (A.rows < 0 || A.rows != A.cols ||
        static_cast<std::size_t>(A.rows) != n || x.size() != n)
    {
        throw std::invalid_argument(
            "steepest_descent_solver: A must be square and match the sizes of b and x");
    }
    if (A.rowPtr.size() != n + 1 || A.values.size() != A.colIndex.size() ||
        A.rowPtr.front() != 0 || A.rowPtr.back() < 0 ||
        static_cast<std::size_t>(A.rowPtr.back()) != A.values.size())
    {
        throw std::invalid_argument("steepest_descent_solver: invalid CSR storage dimensions");
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        if (A.rowPtr[i] < 0 || A.rowPtr[i] > A.rowPtr[i + 1])
        {
            throw std::invalid_argument("steepest_descent_solver: invalid CSR row pointers");
        }
        if (!std::isfinite(b[i]) || !std::isfinite(x[i]))
        {
            throw std::invalid_argument("steepest_descent_solver: b and x must be finite");
        }
    }
    for (std::size_t k = 0; k < A.values.size(); ++k)
    {
        if (A.colIndex[k] < 0 || A.colIndex[k] >= A.cols ||
            !std::isfinite(A.values[k]))
        {
            throw std::invalid_argument("steepest_descent_solver: invalid CSR column or matrix value");
        }
    }
}

double checkedNorm(const std::vector<double>& v)
{
    const double result = vectorNorm(v);
    if (!std::isfinite(result))
    {
        throw std::runtime_error("steepest_descent_solver: nonfinite vector norm");
    }
    return result;
}

std::vector<double> applyChecked(
    const Preconditioner& preconditioner,
    const std::vector<double>& rhs)
{
    std::vector<double> result = preconditioner.applyInverse(rhs);
    if (result.size() != rhs.size())
    {
        throw std::invalid_argument("steepest_descent_solver: preconditioner returned the wrong vector size");
    }
    for (double value : result)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("steepest_descent_solver: preconditioner returned a nonfinite value");
        }
    }
    return result;
}

} // namespace

// Preserve the original public helper: identity preconditioning means z = r.
double steepest_descent_solver::computeAlpha(
    const CSRMatrix& A,
    const std::vector<double>& r)
{
    return computeAlpha(A, r, r);
}

double steepest_descent_solver::computeAlpha(
    const CSRMatrix& A,
    const std::vector<double>& r,
    const std::vector<double>& z)
{
    const double rz = dotProduct(r, z);
    if (!std::isfinite(rz) || rz <= 0.0)
    {
        throw std::runtime_error(
            "steepest_descent_solver: r.z must be positive; use an SPD preconditioner");
    }
    const std::vector<double> Az = multiply(A, z);
    const double zAz = dotProduct(z, Az);
    if (!std::isfinite(zAz) || zAz <= 0.0)
    {
        throw std::runtime_error(
            "steepest_descent_solver: z.Az must be positive; A must be SPD");
    }
    // PSD uses r^T z, not z^T z, in the numerator.
    const double alpha = rz / zAz;
    if (!std::isfinite(alpha) || alpha <= 0.0)
    {
        throw std::runtime_error("steepest_descent_solver: invalid step size");
    }
    return alpha;
}

void steepest_descent_solver::updatePhi(
    std::vector<double>& x,
    double alpha,
    const std::vector<double>& direction)
{
    axpy(x, alpha, direction);
}

void steepest_descent_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x)
{
    const IdentityPreconditioner identity;
    solve(A, b, x, identity);
}

void steepest_descent_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const Preconditioner& preconditioner)
{
    iters_ = 0;
    validateSolve(A, b, x, tolerance, maxIterations);

    double bNorm = checkedNorm(b);
    if (bNorm == 0.0) bNorm = 1.0; // Absolute tolerance for a zero RHS.

    std::vector<double> r = computeResidual(A, x, b);
    double residualNorm = checkedNorm(r);

    while (residualNorm / bNorm > tolerance && iters_ < maxIterations)
    {
        // Solve M z = r once; construction/factorization happened before solve().
        const std::vector<double> z = applyChecked(preconditioner, r);
        const double alpha = computeAlpha(A, r, z);
        updatePhi(x, alpha, z); // x_(k+1) = x_k + alpha z_k.
        ++iters_;

        // Retain the same true-residual SpMV for both identity and IC(0) runs.
        r = computeResidual(A, x, b);
        residualNorm = checkedNorm(r);
    }

    if (verbose)
    {
        if (residualNorm / bNorm <= tolerance)
        {
            std::cout << "Steepest Descent converged in "
                      << iters_ << " iterations." << std::endl;
        }
        else
        {
            std::cout << "Steepest Descent did not converge within "
                      << maxIterations << " iterations." << std::endl;
        }
        std::cout << "Relative residual ||r||/||b|| = "
                  << residualNorm / bNorm << std::endl;
    }
}
