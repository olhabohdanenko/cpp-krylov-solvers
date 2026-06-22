#include <iostream>
#include <chrono>
#include "Vector.h"
#include "Matrix.h"
#include "Preconditioner.h"
#include "SolverResult.h"
#include "Solvers.h"

using namespace std;

SolverResult Solvers::solveCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
    SolverResult res;
    res.method_name = "CG";
    res.preconditioner_name = "Custom";

    auto start_time = chrono::high_resolution_clock::now();

    int n = b.size();
    Vector x(n);
    Vector r = b - (A * x);

    Vector z = M.apply(r);
    Vector p = z;

    double rho_old = r.dot(z);

    for (int i = 0; i < max_iter; i++)
    {
        res.iterations = i + 1;

        Vector Ap = A * p;
        double alpha = rho_old / p.dot(Ap);

        x += p * alpha;
        r -= p * alpha;

        if (r.norm() < tol)
        {
            res.status = SolverStatus::Success;
            break;
        }

        z = M.apply(r);
        double rho_new = r.dot(z);
        double beta = rho_new / rho_old;

        p = z + (p * beta);
        rho_old = rho_new;
    }

    auto end_time = chrono::high_resolution_clock::now();
    res.time_solve = chrono::duration<double>(end_time - start_time).count();
    res.final_norm = r.norm();

    return res;
}