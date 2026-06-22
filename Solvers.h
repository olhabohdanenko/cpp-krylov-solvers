#pragma once
#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "Preconditioner.h"
#include "SolverResult.h"

namespace Solvers
{
    SolverResult solveCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);
    SolverResult solveBiCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);
    SolverResult solveBiCGSTAB (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);
    SolverResult solveGMRES (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, int m, double tol);
}