#pragma once
#include <iostream>

#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>

#include "Vector.h"
#include "Matrix.h"

#include "Preconditioner.h"

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


// Jacobi
class PetscJacobiPreconditioner : public PetscPreconditioner
{
public:
	PetscJacobiPreconditioner(const Matrix& A);
};


// GaussSeidel
class PetscGaussSeidelPreconditioner : public PetscPreconditioner
{
public:
	PetscGaussSeidelPreconditioner(const Matrix& A);
};


// SOR
class PetscSORPreconditioner : public PetscPreconditioner
{
public:
	PetscSORPreconditioner(const Matrix& A, double omega = 1.0);
};


// SSOR
class PetscSSORPreconditioner : public PetscPreconditioner
{
public:
	PetscSSORPreconditioner(const Matrix& A, double omega = 1.0);
};


// ASM (Additiw Schwarz Method)
class PetscASMPreconditioner : public PetscPreconditioner
{
public:
	PetscASMPreconditioner(const Matrix& A, int overlap = 1);
};


// ILU(0)
class PetscILUPreconditioner : public PetscPreconditioner
{
public:
	PetscILUPreconditioner(const Matrix& A);
};


// ILU(p)
class PetscILUPPreconditioner : public PetscPreconditioner
{
public:
	PetscILUPPreconditioner(const Matrix& A, int p = 2);
};


// ILU(p, tau)
class PetscILUPTauPreconditioner : public PetscPreconditioner
{
public:
	PetscILUPTauPreconditioner(const Matrix& A, int p = 2, double tau = 1E-3);
};


// icc (incomplete cholecky)
class PetscICCreconditioner : public PetscPreconditioner
{
public:
	PetscICCreconditioner(const Matrix& A);
};
