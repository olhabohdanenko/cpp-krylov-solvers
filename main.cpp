// main.cpp
#include <iostream>
#include <random>
#include <functional>

#include <petscsys.h>
#include <Tpetra_Core.hpp>

#include "Vector.h"
#include "Matrix.h"

#include "Preconditioners/Preconditioner.h"
#include "Preconditioners/CustomPreconditioner.h"
#include "Preconditioners/PetscPreconditioner.h"
#include "Preconditioners/TrilinosPreconditioner.h"

#include "Solvers/CustomSolver.h"
#include "Solvers/PetscSolver.h"
#include "Solvers/TrilinosSolver.h"

#include "Logger.h"
#include "MemoryTracker.h"

using namespace std;

struct SolverDef
{
	string name;
	function<SolverResult(const Matrix&, const Vector&, const Preconditioner&, int, double)> solve;
};

struct PrecondDef
{
	string name;
	function<unique_ptr<Preconditioner>(const Matrix&)> builder;
};

int main(int argc, char* argv[])
{
	PetscInitialize(&argc, &argv, nullptr, nullptr);
	Tpetra::ScopeGuard tpetraScope(&argc, &argv);
	{
		random_device rd;
		mt19937 gen(rd());
		
		Matrix A (0, 0);
		
		string path = "../matrices/shallow_water1.mtx";
		string json_file = "benchmark_results.json";
		
		if (A.loadFromMtx(path))
			cout << "rows: " << A.rows() << " cols: " << A.cols() << endl;
		else
		{
			cout << "no open file(((((" << endl;
			return 1;
		}
		
		int n = A.rows();
		double eps = 1e-14;
		
		int overlap = 1;
		int degree = 2;
		
		int m = 50;
		
		double omega = 1.5;
		
		int p = 2;
		double k = 2.0;
		double tau = 1e-3;
		
		uniform_real_distribution<double> distr(-50.0, 50.0);
		Vector x_exact(n);
		for(int i = 0; i < n; i++)
			x_exact[i] = distr(gen);
		
		Vector b = A * x_exact;
		Vector rel_er(n);
		
		vector<PrecondDef> precond_definitions =
		{
			{"P ILU(0)", [](const Matrix& mat) { return make_unique<PetscILUPreconditioner>(mat); }},
			{"T ILU(0)", [](const Matrix& mat) { return make_unique<TrilinosRILUKPreconditioner>(mat); }},
			
			{"P ILU(p)", [p](const Matrix& mat) { return make_unique<PetscILUPPreconditioner>(mat, p); }},
			{"T ILU(p)", [p](const Matrix& mat) { return make_unique<TrilinosRILUKPreconditioner>(mat, p); }},
			
			{"P ILU(p,tau)", [p, tau](const Matrix& mat) { return make_unique<PetscILUPTauPreconditioner>(mat, p, tau); }},
			{"T ILU(p, tau)", [p, tau](const Matrix& mat) { return make_unique<TrilinosILUPTauPreconditioner>(mat, p, tau); }},
			
			{"P ASM", [overlap](const Matrix& mat) { return make_unique<PetscASMPreconditioner>(mat, overlap); }},
			{"P ICC", [](const Matrix& mat) { return make_unique<PetscICCreconditioner>(mat); }},
			{"T ILU(k)", [k](const Matrix& mat) { return make_unique<TrilinosILUTPreconditioner>(mat, k); }},
			{"T Chebyshev", [degree](const Matrix& mat) { return make_unique<TrilinosChebyshevPreconditioner>(mat, degree); }},
			{"T Schwarz(ilut)", [](const Matrix& mat) { return make_unique<TrilinosSchwarzPreconditioner>(mat); }}
		};
		
		CustomSolver custom_backend;
		PetscSolver petsc_backend;
		TrilinosSolver trilinos_backend;
		
		vector<SolverDef> solver_definitions =
		{
			{"Custom_CG",       [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return custom_backend.CG(mat, vec_b, prec, dim, tol); }},
			{"Custom_BiCG",     [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return custom_backend.BiCG(mat, vec_b, prec, dim, tol); }},
			{"Custom_BiCGSTAB", [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return custom_backend.BiCGSTAB(mat, vec_b, prec, dim, tol); }},
			{"Custom_GMRES", [&, m](auto& mat, auto& vec_b, auto& prec, int max_iter, double tol) { return custom_backend.GMRES(mat, vec_b, prec, max_iter, m, tol); }},

			{"Petsc_CG",        [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return petsc_backend.solve("cg", mat, vec_b, prec, {}, dim, tol); }},
			{"Petsc_BiCG",      [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return petsc_backend.solve("bicg", mat, vec_b, prec, {}, dim, tol); }},
			{"Petsc_BiCGSTAB",  [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return petsc_backend.solve("bcgs", mat, vec_b, prec, {}, dim, tol); }},
			{"Petsc_GMRES",     [&, m](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return petsc_backend.solve("gmres", mat, vec_b, prec, {{"-ksp_gmres_restart", to_string(m)}}, dim, tol); }},
			{"Petsc_PipeFGMRES", [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return petsc_backend.solve("pipefgmres", mat, vec_b, prec, {}, dim, tol); }},

			{"Trilinos_CG",       [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return trilinos_backend.solve("Block CG", mat, vec_b, prec, {}, dim, tol); }},
			{"Trilinos_BiCGSTAB", [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return trilinos_backend.solve("BiCGStab", mat, vec_b, prec, {}, dim, tol); }},
			{"Trilinos_TFQMR",    [&](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return trilinos_backend.solve("tfqmr", mat, vec_b, prec, {}, dim, tol); }},
			{"Trilinos_GMRES",    [&, m](auto& mat, auto& vec_b, auto& prec, int dim, double tol) { return trilinos_backend.solve("Block Gmres", mat, vec_b, prec, {{"restart", to_string(m)}}, dim, tol); }}
		};
		
		cout << "compute matrices patch " << path << endl;
		
		for (const auto& solver : solver_definitions)
		{
			cout << endl;
			Logger::printHeader();
			
			for (const auto& precond : precond_definitions)
			{
				auto start_prec = chrono::high_resolution_clock::now();
				size_t mem_before = MemoryTracker::getCurrentMemory();
				
				unique_ptr<Preconditioner> M = precond.builder(A);
				
				size_t mem_after = MemoryTracker::getCurrentMemory();
				auto end_prec = chrono::high_resolution_clock::now();
				
				double memory_precond = (mem_after > mem_before) ? (mem_after - mem_before) : 0;
				double time_prec = chrono::duration<double>(end_prec - start_prec).count();
				
				
				mem_before = MemoryTracker::getCurrentMemory();
				
				SolverResult res = solver.solve(A, b, *M, n, eps);
				
				mem_after = MemoryTracker::getCurrentMemory();
				double memory_solve = (mem_after > mem_before) ? (mem_after - mem_before) : 0;
				
				if (res.x.size() != n)
				{
					cout << "Error: Solver " << solver.name << " returned vector with wrong size!" << endl;
					return 1;
				}
				
				for (int i = 0; i < n; i++)
					rel_er[i] = res.x[i] - x_exact[i];
				res.rel_error = rel_er.norm() / x_exact.norm();
				
				res.preconditioner_name = precond.name;
				res.time_precond = time_prec;
				res.memory_precond_mb = memory_precond;
				res.memory_solve_mb = memory_solve;
				
				Logger::log(res);
				Logger::saveToJson(res, json_file);
			}
		}
	}
	PetscFinalize();

	return 0;
}
