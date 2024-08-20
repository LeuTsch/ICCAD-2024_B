#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <list>

#include "inst.h"
#include "parser.h"
#include "solver.h"

namespace Solver
{
    class evaluator
    {
    public:
        evaluator() : _ptrSolver(nullptr) {};
        double evaluate();

    private:
        // data part
        Solver *_ptrSolver;
        pair<double, double> _DieBasePoint; // the Lower left of die
        double _binWidth;
        double _binHeight;
        double _binMaxUtil;
        double _totalTNS;
        double _totalPower;
        double _totalArea;
        vector<vector<double>> _areaInBin;
        vector<vector<double>> _areaOfBin;

        // auxilary function
        };
}
#endif