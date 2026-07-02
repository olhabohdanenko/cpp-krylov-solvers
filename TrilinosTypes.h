#pragma once

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Tpetra_ConfigDefs.hpp>

using Scalar = double;
using LocalOrdinal = Tpetra::Map<>::local_ordinal_type;
using GlobalOrdinal = Tpetra::Map<>::global_ordinal_type;
using Node = Tpetra::Map<>::node_type;

using Node = Tpetra::Details::DefaultTypes::node_type;

using TpetraMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using TpetraVector = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
using Ifpack2Prec = Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
