#pragma once
#include "csr_matrix.h"
#include <string>
namespace cfd {
    CSRMatrix fivePointPattern(int nx,int ny);
    int coefficientPosition(const CSRMatrix& A,int row,int col);
    void addCoefficient(CSRMatrix& A,int row,int col,double value);
    void zeroValues(CSRMatrix& A);
    void validateCSR(const CSRMatrix& A);
    // Eliminate BOTH reference row and column. A retains its existing structural zeros.
    void pinZeroReference(CSRMatrix& A,std::vector<double>& b,int reference=0);
    CSRMatrix transposeCSR(const CSRMatrix& A);
    CSRMatrix productCSR(const CSRMatrix& A,const CSRMatrix& B);
    void dumpCSR(const CSRMatrix& A,const std::vector<double>& b,const std::string& prefix);
    double normInf(const std::vector<double>& x);
}
