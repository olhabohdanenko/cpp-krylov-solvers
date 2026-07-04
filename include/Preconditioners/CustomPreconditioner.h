#pragma once
#include <iostream>

#include <Eigen/Sparse>

#include "Vector.h"
#include "Matrix.h"

#include "Preconditioner.h"

class CustomPreconditioner : public Preconditioner
{
public:
	virtual ~CustomPreconditioner () = default;
};


// Jacobi
class JacobiPreconditioner : public CustomPreconditioner
{
private:
	Vector inv_diag;

public:
	explicit JacobiPreconditioner (const Matrix& A);

	Vector apply (const Vector& x) const override;

	Vector matvec (const Vector& v) const override;
	Vector rmatvec (const Vector& v) const override;
};


// GaussSeidel
class GaussSeidelPreconditioner : public CustomPreconditioner
{
private:
	const Matrix& mat_ref;

public:
	explicit GaussSeidelPreconditioner (const Matrix& A);

	Vector apply (const Vector& x) const override;

	Vector matvec (const Vector& v) const override;
	Vector rmatvec (const Vector& v) const override;
};


// SOR
class SORPreconditioner : public CustomPreconditioner
{
private:
	double omega;
	Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt> M_;

public:
	explicit SORPreconditioner (const Matrix& A, double omega = 1.0);

	Vector apply (const Vector& x) const override;

	Vector matvec (const Vector& v) const override;
	Vector rmatvec (const Vector& v) const override;
};


// SSOR
class SSORPreconditioner : public CustomPreconditioner
{
private:
	double omega;
	Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt> M_lower;
	Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt> M_upper;

public:
	explicit SSORPreconditioner (const Matrix& A, double omega = 1.0);

	Vector apply (const Vector& x) const override;

	Vector matvec (const Vector& v) const override;
	Vector rmatvec (const Vector& v) const override;
};
