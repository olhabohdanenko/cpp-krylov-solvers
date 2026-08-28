// CustomPreconditioner.cpp
#include <iostream>
#include <Eigen/Sparse>

#include "Preconditioners/CustomPreconditioner.h"
#include "Matrix.h"
#include "Vector.h"

// Jacobi
JacobiPreconditioner::JacobiPreconditioner(const Matrix& A)
{
	Eigen::VectorXd diag = A.get_sparse_mat().diagonal();

	for (int i = 0; i < diag.size(); i++)
		if (std::abs(diag(i)) < 1e-14)
			diag(i) = 1e-14;

	inv_diag = Vector(diag.array().inverse());
}

Vector JacobiPreconditioner::apply(const Vector& x) const
{
	return Vector(x.get_Eigen_vec().array() * inv_diag.get_Eigen_vec().array());
}

Vector JacobiPreconditioner::matvec(const Vector& v) const
{
	return Vector(inv_diag.get_Eigen_vec().array() * v.get_Eigen_vec().array());
}

Vector JacobiPreconditioner::rmatvec(const Vector& v) const
{
	return Vector(inv_diag.get_Eigen_vec().array() * v.get_Eigen_vec().array());
}


// GaussSeidel
GaussSeidelPreconditioner::GaussSeidelPreconditioner(const Matrix& A) : mat_ref(A) {}

Vector GaussSeidelPreconditioner::apply(const Vector& x) const
{
	Eigen::VectorXd res = mat_ref.get_sparse_mat().triangularView<Eigen::Lower>().solve(x.get_Eigen_vec());
	return Vector(res);
}

Vector GaussSeidelPreconditioner::matvec(const Vector& v) const
{
	Eigen::VectorXd res = mat_ref.get_sparse_mat().triangularView<Eigen::Lower>().solve(v.get_Eigen_vec());
	return Vector(res);
}

Vector GaussSeidelPreconditioner::rmatvec(const Vector& v) const
{
	Eigen::VectorXd res = mat_ref.get_sparse_mat().triangularView<Eigen::Lower>().transpose().solve(v.get_Eigen_vec());
	return Vector(res);
}


// SOR
SORPreconditioner::SORPreconditioner(const Matrix& A, double omega) : omega(omega)
{
	const auto& A_sparse = A.get_sparse_mat();
	int rows = A.rows();
	int cols = A.cols();

	M_.resize(rows, cols);

	M_ = A_sparse.triangularView<Eigen::Lower>();

	for (int i = 0; i < rows; i++)
	{
		double d_ii = A_sparse.coeff(i, i);
		if (std::abs(d_ii) < 1e-14)
			d_ii = 1e-14;

		M_.coeffRef(i, i) = d_ii / omega;
	}

	M_.makeCompressed();
}

Vector SORPreconditioner::apply(const Vector& x) const
{
	Eigen::VectorXd res = M_.triangularView<Eigen::Lower>().solve(x.get_Eigen_vec());
	return Vector(res);
}

Vector SORPreconditioner::matvec(const Vector& v) const
{
	Eigen::VectorXd res = M_.triangularView<Eigen::Lower>().solve(v.get_Eigen_vec());
	return Vector(res);
}

Vector SORPreconditioner::rmatvec(const Vector& v) const
{
	Eigen::VectorXd res = M_.triangularView<Eigen::Lower>().transpose().solve(v.get_Eigen_vec());
	return Vector(res);
}


// SSOR
SSORPreconditioner::SSORPreconditioner(const Matrix& A, double omega) : omega(omega)
{
	const auto& A_sparse = A.get_sparse_mat();
	int rows = A.rows();
	int cols = A.cols();

	M_lower.resize(rows, cols);
	M_upper.resize(rows, cols);

	M_lower = A_sparse.triangularView<Eigen::Lower>();
	M_upper = A_sparse.triangularView<Eigen::Upper>();

	for (int i = 0; i < rows; i++)
	{
		double d_ii = A_sparse.coeff(i, i) + 1e-4;
		if (std::abs(d_ii) < 1e-14)
			d_ii = 1e-14;

		M_lower.coeffRef(i, i) = d_ii / omega;
		M_upper.coeffRef(i, i) = d_ii / omega;
	}

	M_lower.makeCompressed();
	M_upper.makeCompressed();
}

Vector SSORPreconditioner::apply(const Vector& x) const
{
	Eigen::VectorXd y = M_lower.triangularView<Eigen::Lower>().solve(x.get_Eigen_vec());
	Eigen::VectorXd z = M_upper.triangularView<Eigen::Upper>().solve(y);
	return Vector(z);
}

Vector SSORPreconditioner::matvec(const Vector& v) const
{
	Eigen::VectorXd y = M_lower.triangularView<Eigen::Lower>().solve(v.get_Eigen_vec());
	Eigen::VectorXd z = M_upper.triangularView<Eigen::Upper>().solve(y);
	return Vector(z);
}

Vector SSORPreconditioner::rmatvec(const Vector& v) const
{
	Eigen::VectorXd y = M_upper.triangularView<Eigen::Upper>().transpose().solve(v.get_Eigen_vec());
	Eigen::VectorXd z = M_lower.triangularView<Eigen::Lower>().transpose().solve(y);
	return Vector(z);
}
