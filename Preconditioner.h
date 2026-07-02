#pragma once
#include <iostream>

#include "TrilinosTypes.h"

#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Teuchos_RCP.hpp>

#include "Vector.h"
#include "Matrix.h"

//~ enum class PreconditionerType
//~ {
	//~ None,
	//~ Custom,
	//~ PETSc,
	//~ Trilinos
//~ };
class Preconditioner
{
    public:
    virtual ~Preconditioner () = default;

    virtual Vector apply (const Vector& x) const = 0;

    virtual Vector matvec (const Vector& v) const = 0;
    virtual Vector rmatvec (const Vector& v) const = 0;
};

class NoPreconditioner : public Preconditioner
{
    public:
    Vector apply (const Vector& x) const override { return x; }
    Vector matvec (const Vector& v) const override { return v; }
    Vector rmatvec (const Vector& v) const override { return v; }
};

enum class CustomType
{
	Jacobi
};
class CustomPreconditioner : public Preconditioner
{
	public:
	virtual ~CustomPreconditioner () = default;
};

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

enum class PetscType
{
	Jacobi,
	ASM
};
class PetscPreconditioner : public Preconditioner
{
	protected:
	PC pc;
	Mat A_petsc;
	Vec work_in, work_out;
	
	Vector extractResult(int size) const;
	
	public:
	PetscPreconditioner (const Matrix& A);
	
	virtual ~PetscPreconditioner();
	
	Vector apply (const Vector& x) const override;
	
	Vector matvec (const Vector& v) const override;
    Vector rmatvec (const Vector& v) const override;
};

class PetscJacobiPreconditioner : public PetscPreconditioner
{
	public:
	PetscJacobiPreconditioner(const Matrix& A);
};

class PetscASMPreconditioner : public PetscPreconditioner
{
	public:
	PetscASMPreconditioner(const Matrix& A, int overlap = 1);
};

enum class TrilinosType
{
	Jacobi,
	ILUT
};

class TrilinosPreconditioner : public Preconditioner
{
	protected:
	Teuchos::RCP<TpetraMatrix> A_trilinos;
	Teuchos::RCP<Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>> prec;
	TrilinosType type; /////
	
	void compute();
	
	public:
	TrilinosPreconditioner (const Matrix& A, TrilinosType type = TrilinosType::ILUT);
	
	virtual ~TrilinosPreconditioner() = default;
	
	Vector apply (const Vector& x) const override;
	
	Vector matvec (const Vector& v) const override;
    Vector rmatvec (const Vector& v) const override;
};
