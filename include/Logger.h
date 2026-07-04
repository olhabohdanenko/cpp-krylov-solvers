#pragma once
#include <iostream>
#include <fstream>

#include "SolverResult.h"

class Logger
{
public:
	static void printHeader ();

	static void log (const SolverResult& res);
	
	static void saveToJson(const SolverResult& res, const std::string& filename);
	
	static void finalizeJson(const std::string& filename);
};
