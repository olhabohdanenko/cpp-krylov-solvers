// TrilinosSolver.h
#pragma once
#include "Solvers/Solver.h"
#include "Matrix.h"
#include "Vector.h"
#include "SolverResult.h"
#include "Preconditioners/Preconditioner.h"
#include <string>
#include <map>

class TrilinosSolver : public Solver
{
public:
	SolverResult solve(const std::string& method_name, const Matrix& A, const Vector& b, const Preconditioner& M, const std::map<std::string, std::string>& params = {}, int max_iter = 1000, double tol = 1e-7);
};
