// CustomSolver.cpp
#include <iostream>
#include <chrono>
#include <memory>
#include <vector>
#include <cmath>

#include "Vector.h"
#include "Matrix.h"
#include "Preconditioners/Preconditioner.h"
#include "SolverResult.h"
#include "Solvers/CustomSolver.h"

SolverResult CustomSolver::CG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
	SolverResult res;
	res.method_name = "C CG";
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

	Vector q(n);

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
		
		q = A * p;
		
		double alpha = rho_new / p.dot(q);

		x += alpha * p;
		r -= alpha * q;
		
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = i + 1;
		rho_old = rho_new;

		if (r_norm / b.norm() < tol)
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

SolverResult CustomSolver::BiCG (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
	SolverResult res;
	res.method_name = "C BiCG";
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

	Vector q(n);
	Vector q_(n);

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
		
		q = A * p;
		q_ = A.transpose() * p_;
		
		double alpha = rho_new / p_.dot(q);

		x += alpha * p;
		r -= alpha * q;
		r_ -= alpha * q_;
		
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = i + 1;
		rho_old = rho_new;

		if (r_norm / b.norm() < tol)
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

SolverResult CustomSolver::BiCGSTAB (const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, double tol)
{
	SolverResult res;
	res.method_name = "C BiCGSTAB";
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

	Vector s(n);
	Vector s_(n);
	Vector t(n);

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
		s = r - alpha * v;
		
		double s_norm = s.norm();

		if (s_norm / b.norm() < tol)
		{
			x += alpha * p_;
			res.norm.push_back(s_norm);
			res.status = SolverStatus::Success;
			break;
		}

		if (r.NotFinite() || x.NotFinite() || p.NotFinite())
		{
			res.status = SolverStatus::Diverged;
			break;
		}
		
		s_ = M.matvec(s);
		
		t = A * s_;
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
		
		if (r_norm / b.norm() < tol)
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

Vector CustomSolver::SOL(const std::vector<Vector>& V, const std::vector<std::vector<double>>& R, const Vector& b_, const Vector& x, int nr)
{
	Vector x_new = x;
	Vector y(nr);

	y[nr - 1] = b_[nr - 1] / R[nr - 1][nr - 1];

	for (int k = nr - 2; k >= 0; k--)
	{
		double sum = 0.0;
		for (int i = k + 1; i < nr; i++)
			sum += R[k][i] * y[i];
		
		y[k] = (b_[k] - sum) / R[k][k];
	}

	for (int i = 0; i < nr; i++)
		x_new += y[i] * V[i];
		
	return x_new;
}

SolverResult CustomSolver::GMRES(const Matrix& A, const Vector& b, const Preconditioner& M, int max_iter, int m, double tol)
{
	SolverResult res;
	res.method_name = "C GMRES(m)";
	res.preconditioner_name = "Custom";
	res.status = SolverStatus::MaxIterReached;

	auto start_time = std::chrono::high_resolution_clock::now();

	bool inner_converged = false;
	int n = b.size();

	Vector x(n);
	Vector delta_x(n);
	Vector init_x(n);

	std::vector<std::vector<double>> H(m + 1, std::vector<double>(m, 0.0));
	std::vector<std::vector<double>> R(m, std::vector<double>(m, 0.0));

	Vector c(m);
	Vector s(m);

	std::vector<Vector> V(m + 1, Vector(n));

	Vector e1(m + 1);
	e1[0] = 1.0;

	Vector r = b - (A * x);

	double r_norm = r.norm();
	double b_norm = b.norm();
	res.norm.push_back(r_norm);

	Vector w(n);
	Vector b_(m + 1);

	int i;
	int nr = m;
	int iter = 0;

	for (i = 0; i < max_iter; i++)
	{
		inner_converged = false;
		double beta = r.norm();
		V[0] = r / beta;
		
		b_ = beta * e1;
		
		for (int j = 0; j < m; j++)
		{
			w = A * M.matvec(V[j]);
			
			for (int k = 0; k < j + 1; k++)
			{
				H[k][j] = V[k].dot(w);
				w -= H[k][j] * V[k];
			}
			
			H[j + 1][j] = w.norm();
			V[j + 1] = w / H[j + 1][j];

			for (int k = 0; k <= j; k++)
				R[k][j] = H[k][j];
			
			for (int k = 1; k < j + 1; k++)
			{
				double gamma = c[k - 1] * R[k - 1][j] + s[k - 1] * R[k][j];
				R[k][j] = -s[k - 1] * R[k - 1][j] + c[k - 1] * R[k][j];
				R[k - 1][j] = gamma;
			}
			
			double delta = std::sqrt(R[j][j] * R[j][j] + H[j + 1][j] * H[j + 1][j]);
			
			c[j] = R[j][j] / delta;
			s[j] = H[j + 1][j] / delta;
			
			R[j][j] = c[j] * R[j][j] + s[j] * H[j + 1][j];
			
			b_[j + 1] = -s[j] * b_[j];
			b_[j] = c[j] * b_[j];
			
			double rho = std::abs(b_[j + 1]);
			iter++;
			
			if (rho / b_norm < tol)
			{
				nr = j + 1;
				delta_x = SOL(V, R, b_, init_x, nr);
				inner_converged = true;
				break;
			}
		}
		
		if (!inner_converged)
		{
			nr = m;
			delta_x = SOL(V, R, b_, init_x, nr);
		}
		
		x += M.matvec(delta_x);
		
		r = b - (A * x);
		r_norm = r.norm();
		res.norm.push_back(r_norm);
		res.iterations = iter;
		
		if (r_norm / b_norm < tol)
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

