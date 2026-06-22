#pragma once
#include <iostream>
#include <vector>
#include <string.h>
#include "SolverResult.h"

using namespace std;

class Logger
{
    public:
    static void printHeader ();

    static void log (const SolverResult& res);
};