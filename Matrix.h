// Matrix.h
#pragma once
#include <iostream>
#include <vector>

#include "TrilinosTypes.h"

#include <Eigen/Sparse>

#include <petscmat.h>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_RCP.hpp>

class Vector;

class Matrix
{
    private:
    Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt> mat;
    public:

    Matrix (int rows, int cols);

    void setFromTriplets (const std::vector<Eigen::Triplet<double>>& triplets);
    
    bool loadFromMtx(const std::string& filepath);

    Vector operator* (const Vector& v) const;

    int rows () const { return static_cast<int>(mat.rows()); }
    int cols () const { return static_cast<int>(mat.cols()); }

    const Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt>& get_sparse_mat () const { return mat; }

	Mat createPetscMat () const;
	Teuchos::RCP<TpetraMatrix> createTpetraMat() const;
	Matrix transpose () const;

    void print () const;
};
