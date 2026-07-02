#include <iostream>
#include <memory>
#include <random>
#include <fstream>
#include <functional>

#include <petscsys.h>
#include <Tpetra_Core.hpp>

#include "Vector.h"
#include "Matrix.h"

#include "Preconditioner.h"
#include "Solvers.h"

#include "Logger.h"

using namespace std;

int main(int argc, char* argv[])
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	Tpetra::ScopeGuard tpetraScope(&argc, &argv);
	{
		random_device rd;
		mt19937 gen(rd());
		
		Matrix A (0, 0);
		
		string path = "../matrices/bcsstk08.mtx";
		
		if (A.loadFromMtx(path))
			cout << "rows: " << A.rows() << " cols: " << A.cols() << endl;
		else
		{
			cout << "no open file(((((" << endl;
			return 1;
		}
		
		int n = A.rows();
		double eps = 1e-12;
		
		int overlap = 1;
		
		int m = 30;
		
		double omega = 1.0;
		
		int p = 2;
		double tau = 1e-3;
		
		uniform_real_distribution<double> distr(-50.0, 50.0);
		Vector x_exact(n);
		for(int i = 0; i < n; i++)
			x_exact[i] = distr(gen);
		
		Vector b = A * x_exact;
		Vector rel_er(n);
		
		struct PrecondDef
		{
			string name;
			function<unique_ptr<Preconditioner>(const Matrix&)> builder;
		};
		
		//~ vector<PrecondDef> precond_definitions =
		//~ {
		//~ {"None",          [](const Matrix& mat) { return std::make_unique<NoPreconditioner>(); }},
		//~ {"Custom Jacobi", [](const Matrix& mat) { return std::make_unique<JacobiPreconditioner>(mat); }},
		//~ {"PETSc Jacobi",  [](const Matrix& mat) { return std::make_unique<PetscJacobiPreconditioner>(mat); }},
		//~ {"PETSc ASM",     [](const Matrix& mat) { return std::make_unique<PetscASMPreconditioner>(mat, overlap); }},
		//~ {"Trilinos Jac",  [](const Matrix& mat) { return std::make_unique<TrilinosPreconditioner>(mat, TrilinosType::Jacobi); }},
		//~ {"Trilinos ILUT", [](const Matrix& mat) { return std::make_unique<TrilinosPreconditioner>(mat, TrilinosType::ILUT); }}
		//~ };
		
		vector<PrecondDef> precond_definitions;
        precond_definitions.push_back({"None", [](const Matrix& mat) {
            return make_unique<NoPreconditioner>();
        }});
        precond_definitions.push_back({"Custom Jacobi", [](const Matrix& mat) {
            return make_unique<JacobiPreconditioner>(mat);
        }});
        precond_definitions.push_back({"PETSc Jacobi", [](const Matrix& mat) {
            return make_unique<PetscJacobiPreconditioner>(mat);
        }});
        // Тут захоплюємо overlap за значенням
        precond_definitions.push_back({"PETSc ASM", [overlap](const Matrix& mat) {
            return make_unique<PetscASMPreconditioner>(mat, overlap);
        }});
        precond_definitions.push_back({"Trilinos Jac", [](const Matrix& mat) {
            return make_unique<TrilinosPreconditioner>(mat, TrilinosType::Jacobi);
        }});
        precond_definitions.push_back({"Trilinos ILUT", [](const Matrix& mat) {
            return make_unique<TrilinosPreconditioner>(mat, TrilinosType::ILUT);
        }});
		
		cout << "compute matrices patch " << path << endl;
		Logger::printHeader();
		
		for (const auto& item : precond_definitions)
		{
			auto start_prec = std::chrono::high_resolution_clock::now();
			unique_ptr<Preconditioner> M = item.builder(A);
			auto end_prec = std::chrono::high_resolution_clock::now();
			double time_prec = std::chrono::duration<double>(end_prec - start_prec).count();
			
			SolverResult res = Solvers::solveCG(A, b, *M, n, eps);
			
			if (res.x.size() != n)
			{
				cout << "Error: Solver returned vector with wrong size: " << res.x.size() << " (expected " << n << ")" << endl;
				return 1;
			}

			for (int i = 0; i < n; i++)
				rel_er[i] = res.x[i] - x_exact[i];

			res.rel_error = rel_er.norm() / x_exact.norm();
				
			res.preconditioner_name = item.name;
			
			res.time_precond = time_prec;

			Logger::log(res);
		}
		
		Logger::printHeader();
		
		for (const auto& item : precond_definitions)
		{
			auto start_prec = std::chrono::high_resolution_clock::now();
			unique_ptr<Preconditioner> M = item.builder(A);
			auto end_prec = std::chrono::high_resolution_clock::now();
			double time_prec = std::chrono::duration<double>(end_prec - start_prec).count();
			
			SolverResult res = Solvers::solveBiCG(A, b, *M, n, eps);
			
			if (res.x.size() != n)
			{
				cout << "Error: Solver returned vector with wrong size: " << res.x.size() << " (expected " << n << ")" << endl;
				return 1;
			}

			for (int i = 0; i < n; i++)
				rel_er[i] = res.x[i] - x_exact[i];

			res.rel_error = rel_er.norm() / x_exact.norm();
				
			res.preconditioner_name = item.name;
			
			res.time_precond = time_prec;

			Logger::log(res);
		}
		
		Logger::printHeader();
		
		for (const auto& item : precond_definitions)
		{
			auto start_prec = std::chrono::high_resolution_clock::now();
			unique_ptr<Preconditioner> M = item.builder(A);
			auto end_prec = std::chrono::high_resolution_clock::now();
			double time_prec = std::chrono::duration<double>(end_prec - start_prec).count();
			
			SolverResult res = Solvers::solveBiCGSTAB(A, b, *M, n, eps);
			
			if (res.x.size() != n)
			{
				cout << "Error: Solver returned vector with wrong size: " << res.x.size() << " (expected " << n << ")" << endl;
				return 1;
			}

			for (int i = 0; i < n; i++)
				rel_er[i] = res.x[i] - x_exact[i];

			res.rel_error = rel_er.norm() / x_exact.norm();
				
			res.preconditioner_name = item.name;
			
			res.time_precond = time_prec;

			Logger::log(res);
		}
	}

    PetscFinalize();
    
	return 0;
}

