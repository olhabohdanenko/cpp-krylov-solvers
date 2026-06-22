#pragma once
#include <iostream>
#include <string.h>
#include "Vector.h"

using namespace std;

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

    SolverStatus status;
    int iterations;

    double time_solve;
    double time_precond;

    double final_norm;
    double rel_error;

    long long total_flops;
    double memory_solve_mb;
    double memory_precond_mb;

    string method_name;
    string preconditioner_name;
};