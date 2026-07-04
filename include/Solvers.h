#pragma once
#include <iostream>
#include <vector>

#include "Vector.h"
#include "Matrix.h"
#include "Preconditioners/Preconditioner.h"
#include "SolverResult.h"

namespace Solvers
{
	SolverResult solveCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);
	SolverResult solveBiCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);
	SolverResult solveBiCGSTAB (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol);

	SolverResult solveGMRES (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, int m, double tol);
	static Vector SOL(const std::vector<Vector>& V, const std::vector<std::vector<double>>& R, const Vector& b_, const Vector& x, int nr);
}
