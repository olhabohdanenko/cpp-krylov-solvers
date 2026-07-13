// SolverResult.h
#pragma once
#include <iostream>
#include <vector>
#include <string.h>

#include "Vector.h"

enum class SolverStatus
{
	Success,
	MaxIterReached,
	Diverged,
	Failed
};

struct SolverResult
{
	Vector x;
	std::vector<double> norm;

	SolverStatus status;
	int iterations = 0;

	double time_solve = 0.0;
	double time_precond = 0.0;

	double final_norm = 0.0;
	double rel_error = 0.0;

	long long total_flops = 0;
	double memory_solve_mb = 0.0;
	double memory_precond_mb = 0.0;

	std::string method_name;
	std::string preconditioner_name;
};
