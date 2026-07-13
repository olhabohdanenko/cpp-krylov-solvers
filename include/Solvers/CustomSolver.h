// CustomSolver.h
#pragma once
#include "Solvers/Solver.h"
#include "Matrix.h"
#include "Vector.h"
#include "SolverResult.h"
#include "Preconditioners/Preconditioner.h"

class CustomSolver : public Solver
{
public:
	SolverResult CG(const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter = 1000, double tol = 1e-7);
	SolverResult BiCG(const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter = 1000, double tol = 1e-7);
	SolverResult BiCGSTAB(const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter = 1000, double tol = 1e-7);

	static Vector SOL(const std::vector<Vector>& V, const std::vector<std::vector<double>>& R, const Vector& b_, const Vector& x, int nr);
	SolverResult GMRES(const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter = 1000, int m = 30, double tol = 1e-7);
};
