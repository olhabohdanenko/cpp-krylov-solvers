// Vector.h
#pragma once
#include <iostream>
#include <vector>

#include "TrilinosTypes.h"

#include <Eigen/Dense>

#include <petscvec.h>

#include <Tpetra_Vector.hpp>
#include <Teuchos_RCP.hpp>

class Vector
{
	private:
    Eigen::VectorXd vec;
    public:

    Vector ();
    Vector (int n);
    Vector (const std::vector<double>& vvec);
    Vector (const Eigen::VectorXd& vvec);

    double& operator[] (int i);
    const double& operator[] (int i) const;
    int size () const;
    
    double* data () { return vec.data(); }
    const double* data () const {return vec.data(); }
    
    const Eigen::VectorXd& get_Eigen_vec () const { return vec; }
    Eigen::VectorXd& get_Eigen_vec_nonconst () { return vec; }
    static Vector fromTpetraVec(const Teuchos::RCP<TpetraVector>& tpetra_vec);

    Vector operator+ (const Vector& other) const;
    Vector operator- (const Vector& other) const;

    Vector operator* (double scalar) const;
    friend Vector operator* (double scalar, const Vector& vec_);

    Vector& operator+= (const Vector& other);
    Vector& operator-= (const Vector& other);

    double dot (const Vector& other) const;
    double norm() const;
    
    bool IsFinite() const;
    bool NotFinite() const;
    
    Vec createPetscVecFromVector();
    Teuchos::RCP<TpetraVector> createTpetraVec() const;
};
