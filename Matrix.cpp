#include "Matrix.h"
#include "Vector.h"

using namespace std;

Matrix::Matrix ( int rows, int cols) : mat(rows, cols) {}

void Matrix::setFromTriplets (const vector<Eigen::Triplet<double>>& triplets)
{
    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.makeCompressed();
}

Vector Matrix::operator* (const Vector& v) const
{
    return Vector(mat * v.get_Eigen_vec());
}

void Matrix::print () const
{
    cout << Eigen::MatrixXd(mat) << endl;
}