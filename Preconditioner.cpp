#include "Preconditioner.h"
#include "Matrix.h"
#include <iostream>

using namespace std;

JacobiPreconditioner::JacobiPreconditioner(const Matrix& A)
{
    Eigen::VectorXd diag = A.get_sparse_mat().diagonal();
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