#include <cmath>

#include "TrilinosTypes.h"

#include "Vector.h"

#include <petscvec.h>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Core.hpp>
#include <Ifpack2_Details_Factory.hpp>

Vector::Vector () : vec(Eigen::VectorXd::Zero(0)) {}
Vector::Vector (int n) : vec(Eigen::VectorXd::Zero(n)) {}
Vector::Vector (const std::vector<double>& vvec) : vec(Eigen::Map<const Eigen::VectorXd>(vvec.data(), vvec.size())) {}
Vector::Vector (const Eigen::VectorXd& vvec) : vec(vvec) {}

double& Vector::operator[] (int i)
{
    return vec(i);
}
const double& Vector::operator[] (int i) const
{
    return vec(i);
}
int Vector::size () const
{
    return vec.size();
}

Vector Vector::fromTpetraVec(const Teuchos::RCP<TpetraVector>& tpetra_vec)
{
	int n = static_cast<int>(tpetra_vec->getLocalLength());
	Vector v(n);

	auto vec_data = tpetra_vec->getData();
	for (int i = 0; i < n; i++)
		v[i] = vec_data[i];
		
	return v;
}

Vector Vector::operator+ (const Vector& other) const
{
    return Vector(vec + other.vec);
}
Vector Vector::operator- (const Vector& other) const
{
    return Vector(vec - other.vec);
}

Vector Vector::operator* (double scalar) const
{
    return Vector (vec * scalar);
}
Vector operator* (double scalar, const Vector& vec_)
{
    return vec_ * scalar;
}

Vector& Vector::operator+= (const Vector& other)
{
    vec += other.vec; return *this;
}
Vector& Vector::operator-= (const Vector& other)
{
    vec -= other.vec; return *this;
}

double Vector::dot (const Vector& other) const
{
    return vec.dot(other.vec);
}

double Vector::norm () const
{
    return vec.norm();
}

bool Vector::IsFinite() const
{
	for (int i = 0; i < this->size(); i++)
		if (isnan((*this)[i]) || isinf((*this)[i]))
				return false;
			
	return true;
}

bool Vector::NotFinite() const
{
	for (int i = 0; i < this->size(); i++)
		if (isnan((*this)[i]) || isinf((*this)[i]))
				return true;
			
	return false;
}

Vec Vector::createPetscVecFromVector()
{
	Vec petsc_v;
	VecCreateSeqWithArray(PETSC_COMM_SELF, 1, this->size(), this->data(), &petsc_v);
	return petsc_v;
}

Teuchos::RCP<TpetraVector> Vector::createTpetraVec() const
{
	//auto comm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
	auto comm = Tpetra::getDefaultComm();

	//auto map = Tpetra::createContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(static_cast<GlobalOrdinal>(this->size()), 0, comm);
	auto map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(static_cast<Tpetra::global_size_t>(this->size()), static_cast<size_t>(this->size()), 0, comm));

	auto tpetra_vec = Teuchos::rcp(new TpetraVector(map));

	auto vec_data = tpetra_vec->getDataNonConst();
	for (int i = 0; i < this->size(); i++)
		vec_data[i] = (*this)[i];
	return tpetra_vec;
}
