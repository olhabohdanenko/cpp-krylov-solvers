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
#include <Tpetra_Core.hpp>
#include <Ifpack2_Details_Factory.hpp>

#include "Preconditioner.h"
#include "Matrix.h"

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
    return Vector(inv_diag.get_Eigen_vec().array().inverse() * v.get_Eigen_vec().array());
}

PetscPreconditioner::PetscPreconditioner(const Matrix& A)
{
	A_petsc = A.createPetscMat();
	PCCreate(PETSC_COMM_SELF, &pc);
	PCSetOperators(pc, A_petsc, A_petsc);

	VecCreateSeq(PETSC_COMM_SELF, A.rows(), &work_in);
	VecCreateSeq(PETSC_COMM_SELF, A.rows(), &work_out);
}

PetscPreconditioner::~PetscPreconditioner ()
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

PetscJacobiPreconditioner::PetscJacobiPreconditioner( const Matrix& A) : PetscPreconditioner(A)
{
	PCSetType(pc, PCJACOBI);
	PCSetUp(pc);
}

PetscASMPreconditioner::PetscASMPreconditioner(const Matrix& A, int overlap) : PetscPreconditioner(A)
{
	PCSetType(pc, PCASM);
	PCASMSetOverlap(pc, overlap);
	PCASMSetType(pc, PC_ASM_BASIC);
	PCSetUp(pc);
}

TrilinosPreconditioner::TrilinosPreconditioner(const Matrix& A, TrilinosType type) : type(type)
{
    A_trilinos = A.createTpetraMat();
    
    compute();
}

void TrilinosPreconditioner::compute()
{
    auto A_const = Teuchos::rcp_const_cast<const TpetraMatrix>(A_trilinos);
    Ifpack2::Factory factory;
    Teuchos::ParameterList params;
    
    if (type == TrilinosType::ILUT)
    {
        prec = factory.create<TpetraMatrix>("ILUT", A_const);
        params.set("fact: ilut level-of-fill", 1.0);
    }
    else if (type == TrilinosType::Jacobi)
    {
        prec = factory.create<TpetraMatrix>("RELAXATION", A_const);
        params.set("relaxation: type", "Jacobi");
    }
    
    prec->setParameters(params);
    prec->initialize();
    prec->compute();
}

Vector TrilinosPreconditioner::apply(const Vector& x) const
{
    auto tpetra_x = x.createTpetraVec(); 
    auto tpetra_y = Teuchos::rcp(new TpetraVector(tpetra_x->getMap()));
    
    prec->apply(*tpetra_x, *tpetra_y);
    
    return Vector::fromTpetraVec(tpetra_y); 
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
