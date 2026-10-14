#pragma once

#include <vector>

struct CSRMatrix
{
    int rows = 0;
    int cols = 0;

    // All non-zero matrix values
    std::vector<double> values;

    // Column number corresponding to each value
    std::vector<int> colIndex;

    // rowPtr[i] = starting location of row i
    // Size = rows + 1
    std::vector<int> rowPtr;

    // Location of A[i][i] inside values.
    // Useful later for Gauss-Seidel.
    std::vector<int> diagonalIndex;


    // Default constructor
    CSRMatrix() = default;


    // Constructor
    CSRMatrix(int numberOfRows, int numberOfCols)
        : rows(numberOfRows),
          cols(numberOfCols),
          rowPtr(numberOfRows + 1, 0),
          diagonalIndex(numberOfRows, -1)
    {
    }


    // Number of non-zero elements
    int nnz() const
    {
        return static_cast<int>(values.size());
    }
};


// Matrix-vector multiplication
// result = A * x
std::vector<double> multiply(
    const CSRMatrix& A,
    const std::vector<double>& x
);


// Residual calculation
// r = b - A*x
std::vector<double> computeResidual(
    const CSRMatrix& A,
    const std::vector<double>& x,
    const std::vector<double>& b
);