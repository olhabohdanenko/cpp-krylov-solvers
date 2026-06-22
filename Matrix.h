// Matrix.h
#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Sparse>

using namespace std;

class Vector;

class Matrix
{
    private:
    Eigen::SparseMatrix<double> mat;
    public:

    Matrix (int rows, int cols);

    void setFromTriplets (const vector<Eigen::Triplet<double>>& triplets);

    Vector operator* (const Vector& v) const;

    int rows () const { return static_cast<int>(mat.rows()); }
    int cols () const { return static_cast<int>(mat.cols()); }

    const Eigen::SparseMatrix<double>& get_sparse_mat () const { return mat; }

    void print () const;
};