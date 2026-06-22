#include <iostream>
#include <vector>
#include <string.h>
#include <iomanip>
#include "SolverResult.h"
#include "Logger.h"

using namespace std;

void Logger::printHeader ()
{
    cout << left
        << setw(10) << "Method"
        << setw(15) << "Precond"
        << setw(10) << "Status"
        << setw(10) << "Iter"
        << setw(12) << "Time solve"
        << setw(12) << "Time precond"
        << setw(14) << "Memory solve"
        << setw(14) << "Memory precond"
        << setw(10) << "Norm"
        << setw(10) << "Rel_error" << "\n";
}

void Logger::log (const SolverResult& res)
{
    string status_str;
    switch (res.status)
    {
        case SolverStatus::Success:
            status_str = "OK";
            break;
        case SolverStatus::MaxIterReached:
            status_str = "MaxIter";
            break;
        case SolverStatus::Diverged:
            status_str = "Diverged";
            break;
        default:
            status_str = "Fail";
            break;
    }
    cout << left
        << setw(10) << res.method_name
        << setw(15) << res.preconditioner_name
        << setw(10) << status_str
        << setw(10) << res.iterations
        << setw(12) << res.time_solve
        << setw(12) << res.time_precond
        << setw(14) << res.memory_solve_mb
        << setw(14) << res.memory_precond_mb
        << setw(10) << res.final_norm
        << setw(10) << res.rel_error << "\n";
}