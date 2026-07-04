#include <iostream>

#include <petscksp.h>
#include <petscvec.h>
#include <petscmat.h>

#include "Preconditioners/PetscPreconditioner.h"
#include "Matrix.h"
#include "Vector.h"

PetscPreconditioner::PetscPreconditioner(const Matrix& A)
{
	A_petsc = A.createPetscMat();
	PCCreate(PETSC_COMM_SELF, &pc);
	PCSetOperators(pc, A_petsc, A_petsc);

	VecCreateSeq(PETSC_COMM_SELF, A.rows(), &work_in);
	VecCreateSeq(PETSC_COMM_SELF, A.rows(), &work_out);
}

PetscPreconditioner::~PetscPreconditioner()
{
	PCDestroy(&pc);
	MatDestroy(&A_petsc);
	VecDestroy(&work_in);
	VecDestroy(&work_out);
}

Vector PetscPreconditioner::extractResult(int size) const
{
	const double* arr_out;
	VecGetArrayRead(work_out, &arr_out);
	Vector result(Eigen::Map<const Eigen::VectorXd>(arr_out, size));
	VecRestoreArrayRead(work_out, &arr_out);
	return result;
}

Vector PetscPreconditioner::apply(const Vector& x) const
{
	VecPlaceArray(work_in, const_cast<double*>(x.data()));
	PCApply(pc, work_in, work_out);
	Vector result = extractResult(x.size());
	VecResetArray(work_in);
	return result;
}

Vector PetscPreconditioner::matvec(const Vector& v) const
{
	VecPlaceArray(work_in, const_cast<double*>(v.data()));
	PCApply(pc, work_in, work_out);
	Vector result = extractResult(v.size());
	VecResetArray(work_in);
	return result;
}

Vector PetscPreconditioner::rmatvec(const Vector& v) const
{
	VecPlaceArray(work_in, const_cast<double*>(v.data()));
	PCApplyTranspose(pc, work_in, work_out);
	Vector result = extractResult(v.size());
	VecResetArray(work_in);
	return result;
}


// Jacobi
PetscJacobiPreconditioner::PetscJacobiPreconditioner(const Matrix& A) : PetscPreconditioner(A)
{
	PCSetType(pc, PCJACOBI);
	PCSetUp(pc);
}


// GaussSeidel
PetscGaussSeidelPreconditioner::PetscGaussSeidelPreconditioner(const Matrix& A) : PetscPreconditioner(A)
{
	PCSetType(pc, PCSOR);
	PCSORSetOmega(pc, 1.0);
	PCSORSetSymmetric(pc, SOR_FORWARD_SWEEP);
	PCSetUp(pc);
}


// SOR
PetscSORPreconditioner::PetscSORPreconditioner(const Matrix& A, double omega) : PetscPreconditioner(A)
{
	PCSetType(pc, PCSOR);
	PCSORSetOmega(pc, omega);
	PCSORSetSymmetric(pc, SOR_FORWARD_SWEEP);
	PCSetUp(pc);
}


// SSOR
PetscSSORPreconditioner::PetscSSORPreconditioner(const Matrix& A, double omega) : PetscPreconditioner(A)
{
	PCSetType(pc, PCSOR);
	PCSORSetOmega(pc, omega);
	PCSORSetSymmetric(pc, SOR_SYMMETRIC_SWEEP); 
	PCSetUp(pc);
}


// ASM (Additiw Schwarz Method)
PetscASMPreconditioner::PetscASMPreconditioner(const Matrix& A, int overlap) : PetscPreconditioner(A)
{
	PCSetType(pc, PCASM);
	PCASMSetOverlap(pc, overlap);
	PCASMSetType(pc, PC_ASM_BASIC);
	PCSetUp(pc);
}


// ILU(0)
PetscILUPreconditioner::PetscILUPreconditioner(const Matrix& A) : PetscPreconditioner(A)
{
	PCSetType(pc, PCILU);
	PCFactorSetLevels(pc, 0);
	PCSetUp(pc);
}


// ILU(p)
PetscILUPPreconditioner::PetscILUPPreconditioner(const Matrix& A, int p) : PetscPreconditioner(A)
{
	PCSetType(pc, PCILU);
	PCFactorSetLevels(pc, p);
	PCSetUp(pc);
}


// ILU(p, tau)
PetscILUPTauPreconditioner::PetscILUPTauPreconditioner(const Matrix& A, int p, double tau) : PetscPreconditioner(A)
{
	PCSetType(pc, PCILU);
	PCFactorSetLevels(pc, p);
	PCFactorSetDropTolerance(pc, tau, tau, PETSC_DEFAULT);
	PCSetUp(pc);
}


// icc (incomplete cholecky)
PetscICCreconditioner::PetscICCreconditioner(const Matrix& A) : PetscPreconditioner(A)
{
	PCSetType(pc, PCICC);
	PCSetUp(pc);
}
