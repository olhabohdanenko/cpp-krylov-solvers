#pragma once
#include <iostream>

#include "SolverResult.h"

class Logger
{
    public:
    static void printHeader ();

    static void log (const SolverResult& res);
};
