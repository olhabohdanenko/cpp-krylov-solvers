#include "Vector.h"
using namespace std;

Vector::Vector () : vec(Eigen::VectorXd::Zero(0)) {}
Vector::Vector (int n) : vec(Eigen::VectorXd::Zero(n)) {}
Vector::Vector (const vector<double>& vvec) : vec(Eigen::Map<const Eigen::VectorXd>(vvec.data(), vvec.size())) {}
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
