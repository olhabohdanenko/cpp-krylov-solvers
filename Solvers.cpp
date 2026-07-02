#include <iostream>
#include <chrono>
#include <memory>

#include "Vector.h"
#include "Matrix.h"
#include "Preconditioner.h"
#include "SolverResult.h"
#include "Solvers.h"

SolverResult Solvers::solveCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
    SolverResult res;
    res.method_name = "CG";
    res.preconditioner_name = "Custom";
    res.status = SolverStatus::MaxIterReached;

    auto start_time = std::chrono::high_resolution_clock::now();

    int n = b.size();
    Vector x(n);
    Vector r = b - (A * x);

    Vector z(n);
    Vector p(n);

    double rho_old = 0.0;
    double r_norm = r.norm();
	res.norm.push_back(r_norm);
	
	int i;

    for (i = 0; i < max_iter; i++)
    {
		z = M.matvec(r);
		
		double rho_new = r.dot(z);
		
		if (i == 0)
			p = z;
		else
		{
			double beta = rho_new / rho_old;
			p = z + (beta * p);
		}
		
        Vector q = A * p;
        
        double alpha = rho_new / p.dot(q);

        x += alpha * p;
        r -= alpha * q;
		
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = i + 1;
		rho_old = rho_new;

        if (r_norm < tol)
        {
            res.status = SolverStatus::Success;
            break;
        }

		if (r.NotFinite() || x.NotFinite() || p.NotFinite())
		{
			res.status = SolverStatus::Diverged;
			break;
		}
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    res.time_solve = std::chrono::duration<double>(end_time - start_time).count();
    res.final_norm = r.norm();
    
    res.x = x;

    return res;
}

SolverResult Solvers::solveBiCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
    SolverResult res;
    res.method_name = "BiCG";
    res.preconditioner_name = "Custom";
    res.status = SolverStatus::MaxIterReached;

    auto start_time = std::chrono::high_resolution_clock::now();

    int n = b.size();
    Vector x(n);
    Vector r = b - (A * x);
    Vector r_ = r;

    Vector z(n);
    Vector z_(n);
    Vector p(n);
    Vector p_(n);

    double rho_old = 0.0;
    double r_norm = r.norm();
	res.norm.push_back(r_norm);
	
	int i;

    for (i = 0; i < max_iter; i++)
    {
		z = M.matvec(r);
		z_ = M.rmatvec(r_);
		
		double rho_new = z.dot(r_);
		
		if (std::abs(rho_new) < 5e-16)
		{
			res.status = SolverStatus::Failed;
			break;
		}
		
		if (i == 0)
		{
			p = z;
			p_ = z_;
		}
		else
		{
			double beta = rho_new / rho_old;
			p = z + (beta * p);
			p_ = z_ + (beta * p_);
		}
		
        Vector q = A * p;
        Vector q_ = A.transpose() * p_;
        
        double alpha = rho_new / p_.dot(q);

        x += alpha * p;
        r -= alpha * q;
        r_ -= alpha * q_;
		
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = i + 1;
		rho_old = rho_new;

        if (r_norm < tol)
        {
            res.status = SolverStatus::Success;
            break;
        }

		if (r.NotFinite() || x.NotFinite() || p.NotFinite() || r_.NotFinite() || p_.NotFinite())
		{
			res.status = SolverStatus::Diverged;
			break;
		}
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    res.time_solve = std::chrono::duration<double>(end_time - start_time).count();
    res.final_norm = r.norm();
    
    res.x = x;

    return res;
}

SolverResult Solvers::solveBiCGSTAB (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
    SolverResult res;
    res.method_name = "BiCGSTAB";
    res.preconditioner_name = "Custom";
    res.status = SolverStatus::MaxIterReached;

    auto start_time = std::chrono::high_resolution_clock::now();

    int n = b.size();
    Vector x(n);
    Vector r = b - (A * x);
    Vector r_ = r;

    Vector z(n);
    Vector p(n);
    Vector p_(n);
    Vector v(n);

    double rho_old = 0.0;
    double alpha_new = 0.0;
    double omega_new = 0.0;
    
    double r_norm = r.norm();
	res.norm.push_back(r_norm);
	
	int i;

    for (i = 0; i < max_iter; i++)
    {
		double rho_new = r.dot(r_);
		
		if (std::abs(rho_new) < 5E-16)
		{
			res.status = SolverStatus::Failed;
            break;
		}
		
		if (i == 0)
			p = r;
		else
		{
			double beta = (rho_new / rho_old) * (alpha_new / omega_new);
			p = r + beta * (p - omega_new * v);
		}
		
		p_ = M.matvec(p);
		
		v = A * p_;
        
        double alpha = rho_new / r_.dot(v);
        Vector s = r - alpha * v;
		
		double s_norm = s.norm();

        if (s_norm < tol)
        {
			x += alpha * p_;
            res.status = SolverStatus::Success;
            break;
        }

		if (r.NotFinite() || x.NotFinite() || p.NotFinite())
		{
			res.status = SolverStatus::Diverged;
			break;
		}
		
		Vector s_ = M.matvec(s);
		
		Vector t = A * s_;
		double omega = t.dot(s) / t.dot(t);
		
		if (std::abs(omega) < 5e-16)
		{
			break;
		}
		
		x += (alpha * p_) + (omega * s_);
		r = s - omega * t;
		
		rho_old = rho_new;
		alpha_new = alpha;
		omega_new = omega;
		
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = i + 1;
		
		if (r_norm < tol)
        {
            res.status = SolverStatus::Success;
            break;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    res.time_solve = std::chrono::duration<double>(end_time - start_time).count();
    res.final_norm = r.norm();
    
    res.x = x;

    return res;
}
