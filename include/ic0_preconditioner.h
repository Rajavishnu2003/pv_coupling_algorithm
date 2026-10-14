#pragma once

#include "csr_matrix.h"
#include "preconditioner.h"

#include <vector>


// Incomplete Cholesky with zero additional fill: M = L * L^T.
//
// The input must be a nonempty square CSR matrix with finite entries,
// strictly increasing, unique columns in each row, explicit positive
// diagonals, and numerical symmetry. The constructor builds the factors
// once and leaves A unchanged. Input diagonalIndex entries are not used.
//
// Only positions present in the original lower triangle are retained,
// including explicitly stored zero entries. No shifts or pivoting are
// applied. An incomplete-pivot breakdown throws an exception; it does
// not by itself prove that the original matrix is indefinite.
class IC0Preconditioner final : public Preconditioner
{
public:
    explicit IC0Preconditioner(const CSRMatrix& A);

    std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const override;

    int size() const;

    // Number of stored positions in L, including its diagonal and any
    // retained structural zeros. L^T is not counted a second time.
    int factorNonzeros() const;

private:
    CSRMatrix lowerFactor;
    CSRMatrix transposedLowerFactor;

    void buildTranspose();
};
