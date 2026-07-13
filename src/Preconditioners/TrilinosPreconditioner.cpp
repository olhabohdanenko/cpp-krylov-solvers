// TrilinosPreconditioner.cpp
#include <iostream>

#include "TrilinosTypes.h"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>
#include <Ifpack2_Details_Factory.hpp>

#include "Preconditioners/TrilinosPreconditioner.h"
#include "Matrix.h"
#include "Vector.h"

TrilinosPreconditioner::TrilinosPreconditioner(const Matrix& A)
{
	A_trilinos = A.createTpetraMat();
	tpetra_work_out = Teuchos::rcp(new TpetraVector(A_trilinos->getMap()));
}

Vector TrilinosPreconditioner::apply(const Vector& x) const
{
	Teuchos::ArrayView<const Scalar> x_view(x.data(), x.size());
	TpetraVector tpetra_x(A_trilinos->getMap(), x_view);
	
	prec->apply(tpetra_x, *tpetra_work_out);
	return Vector::fromTpetraVec(tpetra_work_out);
}

Vector TrilinosPreconditioner::matvec(const Vector& v) const
{
	return apply(v);
}

Vector TrilinosPreconditioner::rmatvec(const Vector& v) const
{
	auto tpetra_v = v.createTpetraVec();
	auto tpetra_res = Teuchos::rcp(new TpetraVector(tpetra_v->getMap()));

	prec->apply(*tpetra_v, *tpetra_res, Teuchos::TRANS, 1.0, 0.0);

	return Vector::fromTpetraVec(tpetra_res);
}


// Jacobi
TrilinosJacobiPreconditioner::TrilinosJacobiPreconditioner(const Matrix& A) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("RELAXATION", A_const);

	params.set("relaxation: type", "Jacobi");
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// GaussSeidel
TrilinosGaussSeidelPreconditioner::TrilinosGaussSeidelPreconditioner(const Matrix& A) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("RELAXATION", A_const);

	params.set("relaxation: type", "Gauss-Seidel");
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// SOR
TrilinosSORPreconditioner::TrilinosSORPreconditioner(const Matrix& A, double omega) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("RELAXATION", A_const);

	params.set("relaxation: type", "Gauss-Seidel");
	params.set("relaxation: damping factor", omega);
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// SGS
TrilinosSGSPreconditioner::TrilinosSGSPreconditioner(const Matrix& A) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("RELAXATION", A_const);

	params.set("relaxation: type", "Symmetric Gauss-Seidel");
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// Chebyshev
TrilinosChebyshevPreconditioner::TrilinosChebyshevPreconditioner(const Matrix& A, int degree) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("CHEBYSHEV", A_const);

	params.set("chebyshev: degree", degree);
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// Schwarz
TrilinosSchwarzPreconditioner::TrilinosSchwarzPreconditioner(const Matrix& A, const std::string& inner_prec_type) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("SCHWARZ", A_const);

	params.set("schwarz: inner preconditioner name", inner_prec_type);

	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// ILU(0), ILU(p)
TrilinosRILUKPreconditioner::TrilinosRILUKPreconditioner(const Matrix& A, int p) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("RILUK", A_const);

	params.set("fact: iluk level-of-fill", p);
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// ILU(p)
TrilinosILUTPreconditioner::TrilinosILUTPreconditioner(const Matrix& A, double p) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("ILUT", A_const);

	params.set("fact: ilut level-of-fill", static_cast<double>(p));
	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}


// ILU(p, tau)
TrilinosILUPTauPreconditioner::TrilinosILUPTauPreconditioner(const Matrix& A, int p, double tau) : TrilinosPreconditioner(A)
{
	auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
	prec = Ifpack2::Factory::create<TpetraMatrix>("ILUT", A_const);

	params.set("fact: ilut level-of-fill", static_cast<double>(p));
	params.set("fact: drop tolerance", tau);

	prec->setParameters(params);
	prec->initialize();
	prec->compute();
}
