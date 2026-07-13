// TrilinosSolver.cpp
#include "Solvers/TrilinosSolver.h"
#include "Preconditioners/TrilinosPreconditioner.h"
#include "Solvers/TrilinosPreconditionerAdapter.h"
#include "TrilinosTypes.h"
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <regex>

#include <Teuchos_ParameterList.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

#include "Matrix.h"
#include "Vector.h"

SolverResult TrilinosSolver::solve(const std::string& method_name, const Matrix& A, const Vector& b, const Preconditioner& M, const std::map<std::string, std::string>& params, int max_iter, double tol)
{
	SolverResult res;
	res.method_name = "T " + method_name;

	auto start_time = std::chrono::high_resolution_clock::now();

	Teuchos::RCP<TpetraMatrix> A_tpetra = A.createTpetraMat();
	Teuchos::RCP<TpetraVector> b_tpetra = b.createTpetraVec();

	Vector x_initial(b.size());
	Teuchos::RCP<TpetraVector> x_tpetra = x_initial.createTpetraVec();
	
	auto ss_output = Teuchos::rcp(new std::stringstream());
	Teuchos::RCP<std::ostream> belos_output = ss_output;

	Teuchos::ParameterList belosList;
	belosList.set("Maximum Iterations", max_iter);
	belosList.set("Convergence Tolerance", tol);
	belosList.set("Output Stream", belos_output);
	belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::StatusTestDetails);

	if (params.count("restart"))
		belosList.set("Num Blocks", std::stoi(params.at("restart")));
	if (params.count("recycle_blocks"))
		belosList.set("Num Recycled Blocks", std::stoi(params.at("recycle_blocks")));
	if (params.count("orthogonalization"))
		belosList.set("Orthogonalization", params.at("orthogonalization"));

	using ST = double;
	using MV = Tpetra::MultiVector<ST, LocalOrdinal, GlobalOrdinal, Node>;
	using OP = Tpetra::Operator<ST, LocalOrdinal, GlobalOrdinal, Node>;

	Teuchos::RCP<MV> b_mv = b_tpetra;
	Teuchos::RCP<MV> x_mv = x_tpetra;
	Teuchos::RCP<OP> A_op = A_tpetra;

	auto problem = Teuchos::rcp(new Belos::LinearProblem<ST, MV, OP>(A_op, x_mv, b_mv));

	Teuchos::RCP<OP> prec_op;

	if (auto* trilinos_prec = dynamic_cast<const TrilinosPreconditioner*>(&M))
		prec_op = trilinos_prec->getTpetraOp();
	else
		prec_op = Teuchos::rcp(new TrilinosPreconditionerAdapter<ST, LocalOrdinal, GlobalOrdinal, Node>(M, A_tpetra->getDomainMap()));

	std::string lower_method = method_name;
	std::transform(lower_method.begin(), lower_method.end(), lower_method.begin(), ::tolower);

	if (lower_method == "bicgstab")
		problem->setRightPrec(prec_op);
	else
		problem->setLeftPrec(prec_op);

	problem->setProblem();

	Belos::SolverFactory<ST, MV, OP> factory;
	auto solver = factory.create(method_name, Teuchos::rcpFromRef(belosList));
	solver->setProblem(problem);

	Belos::ReturnType ret = solver->solve();

	auto end_time = std::chrono::high_resolution_clock::now();
	res.time_solve = std::chrono::duration<double>(end_time - start_time).count();

	res.status = (ret == Belos::Converged) ? SolverStatus::Success : SolverStatus::Failed;
	res.iterations = solver->getNumIters();
	res.x = Vector::fromTpetraVec(x_tpetra);
	double final_res_norm = solver->achievedTol();
	res.final_norm = final_res_norm;
	
	std::string line;
    std::regex norm_regex(R"((?:Rel Res Norm|Implicit Residual|Direct Residual|Residual)\s*(?:=\s*\[scaled\]|=\s*)?([0-9eE\.\-\+]+))", std::regex_constants::icase);
    std::smatch match;

    ss_output->clear();
    ss_output->seekg(0, std::ios::beg);

    while (std::getline(*ss_output, line))
    {
        if (std::regex_search(line, match, norm_regex))
        {
            if (match.size() > 1)
            {
                try 
                {
                    double val = std::stod(match[1].str());
                    res.norm.push_back(val);
                }
                catch (...) 
                {
                    continue;
                }
            }
        }
    }
    
	return res;
}
