// TrilinosPreconditionerAdapter.h
#pragma once

#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_RCP.hpp>

#include "Preconditioners/Preconditioner.h"
#include "Vector.h"
#include "TrilinosTypes.h"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TrilinosPreconditionerAdapter : public Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
private:
	const Preconditioner& m_prec;
	Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> m_map;

public:
	TrilinosPreconditionerAdapter(const Preconditioner& prec, Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map) : m_prec(prec), m_map(map) {}

	void apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS,
			   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(), Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override 
	{
		auto X_vec = X.getVector(0);
		
		auto X_non_const = Teuchos::rcp_const_cast<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(X_vec);
		Vector framework_x = Vector::fromTpetraVec(X_non_const);
		
		Vector framework_y;
		if (mode == Teuchos::NO_TRANS)
			framework_y = m_prec.matvec(framework_x);
		else
			framework_y = m_prec.rmatvec(framework_x);

		auto Y_vec = Y.getVectorNonConst(0);
		auto y_view = Y_vec->get1dViewNonConst();
		
		for (size_t i = 0; i < framework_y.size(); ++i)
			y_view[i] = alpha * framework_y[i] + beta * y_view[i];
	}

	Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const override { return m_map; }
	Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const override { return m_map; }
};
