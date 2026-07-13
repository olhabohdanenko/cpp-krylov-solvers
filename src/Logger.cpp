// Logger.cpp
#include <iostream>
#include <vector>
#include <string.h>
#include <iomanip>
#include <fstream>

#include "SolverResult.h"
#include "Logger.h"

void Logger::printHeader ()
{
	std::cout << std::left 
		<< std::setw(15) << "Method"
		<< std::setw(16) << "Precond"
		<< std::setw(12) << "Status"
		<< std::setw(8)  << "Iter"
		<< std::setw(14) << "Time solve"
		<< std::setw(14) << "Time prec"
		<< std::setw(15) << "Memory solve b"
		<< std::setw(15) << "Memory prec b"
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
		<< std::setw(15) << res.method_name
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

void Logger::saveToJson(const SolverResult& res, const std::string& filename)
{
	bool is_empty = true;
	{
		std::ifstream check_file(filename);
		if (check_file && check_file.peek() != std::ifstream::traits_type::eof())
			is_empty = false;
	}

	if (is_empty)
	{
		std::ofstream file(filename);
		if (!file.is_open()) return;

		file << "[\n";
		
		auto write_content = [&]()
		{
			file << "  {\n";
			file << "    \"method\": \"" << res.method_name << "\",\n";
			file << "    \"preconditioner\": \"" << res.preconditioner_name << "\",\n";
			file << "    \"iterations\": " << res.iterations << ",\n";
			file << "    \"time_solve\": " << res.time_solve << ",\n";
			file << "    \"time_precond\": " << res.time_precond << ",\n";
			file << "    \"memory_solve_b\": " << res.memory_solve_mb << ",\n";
			file << "    \"memory_precond_b\": " << res.memory_precond_mb << ",\n";
			file << "    \"rel_error\": " << res.rel_error << ",\n";
			file << "    \"convergence_history\": [";
			for (size_t i = 0; i < res.norm.size(); i++)
			{
				file << res.norm[i];
				if (i != res.norm.size() - 1) file << ", ";
			}
			file << "]\n";
			file << "  }";
		};

		write_content();
		file << "\n]";
	}
	else
	{
		std::fstream file(filename, std::ios::in | std::ios::out);
		if (!file.is_open()) return;

		file.seekp(0, std::ios::end);
		std::streamoff pos = file.tellp();
		std::streamoff offset = 1;
		char ch = 0;

		while (pos >= offset)
		{
			file.seekg(-offset, std::ios::end);
			file.get(ch);
			if (ch == ']')
			{
				file.seekp(-offset, std::ios::end);
				break;
			}
			offset++;
		}

		file << ",\n";
		file << "  {\n";
		file << "    \"method\": \"" << res.method_name << "\",\n";
		file << "    \"preconditioner\": \"" << res.preconditioner_name << "\",\n";
		file << "    \"iterations\": " << res.iterations << ",\n";
		file << "    \"time_solve\": " << res.time_solve << ",\n";
		file << "    \"time_precond\": " << res.time_precond << ",\n";
		file << "    \"memory_solve_b\": " << res.memory_solve_mb << ",\n";
		file << "    \"memory_precond_b\": " << res.memory_precond_mb << ",\n";
		file << "    \"rel_error\": " << res.rel_error << ",\n";
		file << "    \"convergence_history\": [";
		for (size_t i = 0; i < res.norm.size(); i++)
		{
			file << res.norm[i];
			if (i != res.norm.size() - 1) file << ", ";
		}
		file << "]\n";
		file << "  }\n]";
	}
}
