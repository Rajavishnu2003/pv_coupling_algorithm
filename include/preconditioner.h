#pragma once

#include <vector>

// A preconditioner supplies the action of M^-1 on a vector.
// Its setup is separate from the outer iterative solver.
class Preconditioner
{
public:
    virtual ~Preconditioner() = default;

    virtual std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const = 0;
};


// M = I recovers each original, unpreconditioned algorithm.
class IdentityPreconditioner final : public Preconditioner
{
public:
    std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const override
    {
        return rhs;
    }
};
