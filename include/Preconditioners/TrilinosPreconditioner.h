// TrilinosPreconditioner.h
#pragma once
#include <iostream>

#include "TrilinosTypes.h"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Teuchos_RCP.hpp>

#include "Vector.h"
#include "Matrix.h"

#include "Preconditioner.h"

class TrilinosPreconditioner : public Preconditioner
{
protected:
	Teuchos::RCP<TpetraMatrix> A_trilinos;
	Teuchos::RCP<Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>> prec;
	Teuchos::ParameterList params;
	
	mutable Teuchos::RCP<TpetraVector> tpetra_work_out;

public:
	explicit TrilinosPreconditioner(const Matrix& A);
	virtual ~TrilinosPreconditioner() = default;
	Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getTpetraOp() const { return prec; }

	Vector apply(const Vector& x) const override;
	
	Vector matvec(const Vector& v) const override;
	Vector rmatvec(const Vector& v) const override;
};


// Jacobi
class TrilinosJacobiPreconditioner : public TrilinosPreconditioner
{
public:
	explicit TrilinosJacobiPreconditioner(const Matrix& A);
};


// GaussSeidel
class TrilinosGaussSeidelPreconditioner : public TrilinosPreconditioner
{
public:
	explicit TrilinosGaussSeidelPreconditioner(const Matrix& A);
};


// SOR
class TrilinosSORPreconditioner : public TrilinosPreconditioner
{
public:
	TrilinosSORPreconditioner(const Matrix& A, double omega = 1.0);
};


// SGS
class TrilinosSGSPreconditioner : public TrilinosPreconditioner
{
public:
	explicit TrilinosSGSPreconditioner(const Matrix& A);
};


// Chebyshev
class TrilinosChebyshevPreconditioner : public TrilinosPreconditioner
{
public:
	TrilinosChebyshevPreconditioner(const Matrix& A, int degree = 2);
};


// Schwarz
class TrilinosSchwarzPreconditioner : public TrilinosPreconditioner
{
public:
	explicit TrilinosSchwarzPreconditioner(const Matrix& A, const std::string& inner_prec_type = "ILUT");
};


// ILU(0), ILU(p)
class TrilinosRILUKPreconditioner : public TrilinosPreconditioner
{
public:
	TrilinosRILUKPreconditioner(const Matrix& A, int p = 0);
};



// ILU(tau)
class TrilinosILUTPreconditioner : public TrilinosPreconditioner
{
public:
	TrilinosILUTPreconditioner(const Matrix& A, double tau = 1.0);
};


// ILU(p, tau)
class TrilinosILUPTauPreconditioner : public TrilinosPreconditioner
{
public:
	TrilinosILUPTauPreconditioner(const Matrix& A, int p = 2, double tau = 1.0);
};
