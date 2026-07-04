#include "TrilinosTypes.h"

#include "Matrix.h"
#include "Vector.h"

#include <unsupported/Eigen/SparseExtra>

#include <petscmat.h>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>

Matrix::Matrix ( int rows, int cols) : mat(rows, cols) {}

void Matrix::setFromTriplets (const std::vector<Eigen::Triplet<double>>& triplets)
{
	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();
}

bool Matrix::loadFromMtx(const std::string& filepath)
{
	std::ifstream file(filepath);
	if (!file.is_open())
		return false;

	std::string line;
	bool isSymmetric = false;
	int n = 0, m = 0, nnz = 0;

	while (std::getline(file, line))
	{
		if (line.empty())
			continue;
		if (line[0] == '%')
		{
			if (line.find("symmetric") != std::string::npos || line.find("Symmetric") != std::string::npos)
				isSymmetric = true;
			continue;
		}

		std::istringstream iss(line);
		iss >> n >> m >> nnz;
		break;
	}

	if (n != m)
	{
		std::cerr << "Error: matrix is not square\n";
		return false;
	}

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(isSymmetric ? 2 * nnz : nnz);

	int row, col;
	double val;
	for (int i = 0; i < nnz; i++)
	{
		file >> row >> col >> val;
		--row; --col;
		triplets.emplace_back(row, col, val);
		if (isSymmetric && row != col)
			triplets.emplace_back(col, row, val);
	}

	mat = Eigen::SparseMatrix<double, Eigen::RowMajor, PetscInt>(n, m);
	mat.setFromTriplets(triplets.begin(), triplets.end());
	mat.makeCompressed();

	return true;
}

Vector Matrix::operator* (const Vector& v) const
{
	return Vector(mat * v.get_Eigen_vec());
}

Mat Matrix::createPetscMat() const
{
	Mat petsc_mat;
	MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,  this->rows(), this->cols(), const_cast<PetscInt*>(mat.outerIndexPtr()), const_cast<PetscInt*>(mat.innerIndexPtr()), const_cast<PetscScalar*>(mat.valuePtr()), &petsc_mat);
	return petsc_mat;
}

Teuchos::RCP<TpetraMatrix> Matrix::createTpetraMat() const
{
	auto comm = Tpetra::getDefaultComm();

	auto map = Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(
		static_cast<Tpetra::global_size_t>(mat.rows()),
		static_cast<size_t>(mat.rows()), 0, comm));

	Teuchos::Array<size_t> numEntriesPerRow(mat.rows(), 0);
	for (int k = 0; k < mat.outerSize(); ++k)
		for (decltype(mat)::InnerIterator it(mat, k); it; ++it)
			numEntriesPerRow[it.row()]++;

	auto tpetra_mat = Teuchos::rcp(new TpetraMatrix(map, numEntriesPerRow));

	for (int k = 0; k < mat.outerSize(); ++k)
	{
		for (decltype(mat)::InnerIterator it(mat, k); it; ++it)
		{
			GlobalOrdinal row = static_cast<GlobalOrdinal>(it.row());
			GlobalOrdinal col = static_cast<GlobalOrdinal>(it.col());
			double val = it.value();

			Teuchos::Array<GlobalOrdinal> cols = {col};
			Teuchos::Array<double> vals = {val};

			tpetra_mat->insertGlobalValues(row, cols, vals);
		}
	}
	tpetra_mat->fillComplete();
	return tpetra_mat;
}

Matrix Matrix::transpose() const
{
	Matrix T(mat.rows(), mat.cols());

	T.mat = mat.transpose();

	return T;
}

void Matrix::print () const
{
	std::cout << Eigen::MatrixXd(mat) << std::endl;
}
