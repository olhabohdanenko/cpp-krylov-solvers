// Vector.h
#pragma once
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std;

class Vector
{
	private:
    Eigen::VectorXd vec;
    public:

    Vector ();
    Vector (int n);
    Vector (const vector<double>& vvec);
    Vector (const Eigen::VectorXd& vvec);

    double& operator[] (int i);
    const double& operator[] (int i) const;
    int size () const;
    const Eigen::VectorXd& get_Eigen_vec () const { return vec; }

    Vector operator+ (const Vector& other) const;
    Vector operator- (const Vector& other) const;

    Vector operator* (double scalar) const;
    friend Vector operator* (double scalar, const Vector& vec_);

    Vector& operator+= (const Vector& other);
    Vector& operator-= (const Vector& other);

    double dot (const Vector& other) const;
    double norm() const;
};
