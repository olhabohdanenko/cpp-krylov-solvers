// PetscSolver.cpp
#include "Solvers/PetscSolver.h"
#include "Preconditioners/PetscPreconditioner.h"
#include <petscksp.h>
#include <chrono>
#include <iostream>
#include <algorithm>

#include "Matrix.h"
#include "Vector.h"

SolverResult PetscSolver::solve(const std::string& method_name, const Matrix& A, const Vector& b, const Preconditioner& M, const std::map<std::string, std::string>& params, int max_iter, double tol)
{
	SolverResult res;
	res.method_name = "P " + method_name;
	res.status = SolverStatus::MaxIterReached;

	auto start_time = std::chrono::high_resolution_clock::now();
	int size = b.size();

	Mat A_petsc = A.createPetscMat();
	Vec b_petsc = const_cast<Vector&>(b).createPetscVecFromVector();

	Vector x_res(size);
	Vec x_petsc = x_res.createPetscVecFromVector();
	
	KSP ksp;
	KSPCreate(PETSC_COMM_SELF, &ksp);
	KSPSetOperators(ksp, A_petsc, A_petsc);
	KSPSetType(ksp, method_name.c_str());

	if (auto* petsc_prec = dynamic_cast<const PetscPreconditioner*>(&M))
		KSPSetPC(ksp, petsc_prec->getPetscPC());
	else
	{
		PC pc;
		KSPGetPC(ksp, &pc);
		PCSetType(pc, PCSHELL);
		
		PCShellSetContext(pc, const_cast<void*>(static_cast<const void*>(&M)));
		
		PCShellSetApply(pc, PetscCustomPCApply);
		PCShellSetApplyTranspose(pc, PetscCustomPCApplyTranspose);
		
		PCShellSetName(pc, "FrameworkGenericPreconditioner");
	}

	KSPSetFromOptions(ksp);
	KSPSetTolerances(ksp, tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);

	KSPSetUp(ksp);
	KSPSetResidualHistory(ksp, NULL, PETSC_DECIDE, PETSC_TRUE);
	KSPSolve(ksp, b_petsc, x_petsc);

	auto end_time = std::chrono::high_resolution_clock::now();
	res.time_solve = std::chrono::duration<double>(end_time - start_time).count();

	PetscInt its;
	KSPGetIterationNumber(ksp, &its);
	res.iterations = its;

	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp, &reason);
	res.status = (reason > 0) ? SolverStatus::Success : SolverStatus::Failed;

	res.x = x_res;
	PetscReal final_norm_;
	KSPGetResidualNorm(ksp, &final_norm_);
	res.final_norm = final_norm_;
	
	const PetscReal *hist;
	PetscInt count = 0;
	KSPGetResidualHistory(ksp, &hist, &count);

	for (PetscInt i = 0; i < count; i++)
		res.norm.push_back(hist[i]);

	VecDestroy(&b_petsc);
	VecDestroy(&x_petsc);
	KSPDestroy(&ksp);
	MatDestroy(&A_petsc);

	return res;
}

PetscErrorCode PetscSolver::PetscCustomPCApply(PC pc, Vec x, Vec y)
{
	void* ctx;
	PCShellGetContext(pc, &ctx);
	const Preconditioner* prec = static_cast<const Preconditioner*>(ctx);

	PetscInt size;
	VecGetSize(x, &size);

	const PetscScalar* x_arr;
	VecGetArrayRead(x, &x_arr);
	Vector framework_x(Eigen::Map<const Eigen::VectorXd>(x_arr, size));
	VecRestoreArrayRead(x, &x_arr);

	Vector framework_y = prec->matvec(framework_x);

	PetscScalar* y_arr;
	VecGetArray(y, &y_arr);
	std::copy(framework_y.data(), framework_y.data() + size, y_arr);
	VecRestoreArray(y, &y_arr);

	return PETSC_SUCCESS;
}

PetscErrorCode PetscSolver::PetscCustomPCApplyTranspose(PC pc, Vec x, Vec y)
{
	void* ctx;
	PCShellGetContext(pc, &ctx);
	const Preconditioner* prec = static_cast<const Preconditioner*>(ctx);

	PetscInt size;
	VecGetSize(x, &size);

	const PetscScalar* x_arr;
	VecGetArrayRead(x, &x_arr);
	Vector framework_x(Eigen::Map<const Eigen::VectorXd>(x_arr, size));
	VecRestoreArrayRead(x, &x_arr);

	Vector framework_y = prec->rmatvec(framework_x);

	PetscScalar* y_arr;
	VecGetArray(y, &y_arr);
	std::copy(framework_y.data(), framework_y.data() + size, y_arr);
	VecRestoreArray(y, &y_arr);

	return PETSC_SUCCESS;
}
