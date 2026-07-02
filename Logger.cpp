#include <iostream>
#include <vector>
#include <string.h>
#include <iomanip>

#include "SolverResult.h"
#include "Logger.h"

void Logger::printHeader ()
{
	std::cout << std::left 
          << std::setw(10) << "Method"
          << std::setw(16) << "Precond"
          << std::setw(12) << "Status"
          << std::setw(8)  << "Iter"
          << std::setw(14) << "Time solve"
          << std::setw(14) << "Time prec"
          << std::setw(15) << "Memory solve"
          << std::setw(15) << "Memory prec"
          << std::setw(14) << "Norm"
          << std::setw(14) << "Rel_error" << std::endl;
}

void Logger::log (const SolverResult& res)
{
    std::string status_str;
    switch (res.status)
    {
        case SolverStatus::Success:
            status_str = "Success";
            break;
        case SolverStatus::MaxIterReached:
            status_str = "MaxIter";
            break;
        case SolverStatus::Diverged:
            status_str = "Diverged";
            break;
        default:
            status_str = "Failed";
            break;
    }
    std::cout << std::left << std::scientific << std::setprecision(5)
          << std::setw(10) << res.method_name
          << std::setw(16) << res.preconditioner_name
          << std::setw(12) << status_str
          << std::setw(8)  << res.iterations
          << std::setw(14) << res.time_solve
          << std::setw(14) << res.time_precond
          << std::setw(15) << res.memory_solve_mb
          << std::setw(15) << res.memory_precond_mb
          << std::setw(14) << res.final_norm
          << std::setw(14) << res.rel_error << std::endl;
}
