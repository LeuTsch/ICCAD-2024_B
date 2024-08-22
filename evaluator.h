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
        void setSolverPtr(Solver *ptr) { _ptrSolver = ptr; };
        double evaluate();             // retrun the value of cost metric
        void evaluate(const string &); // write the information of cost metric into the file named by the string
        double evaluateTNS();

    private:
        // data part
        Solver *_ptrSolver;
        pair<double, double> _DieBasePoint; // the Lower left of die
        double _dieWidth;
        double _dieHeight;
        double _binWidth;
        double _binHeight;
        double _binMaxUtil;
        double _totalTNS;
        double _totalPower;
        double _totalArea;
        double _totalUtil;
        double _TNSbyDisplace;
        vector<vector<double>> _areaInBin;
        vector<vector<double>> _areaOfBin;

        // catch
        vector<Inst::FF_D *> _FF_DList;
        vector<Inst::Gate *> _GateList;

        // auxilary function
        void calculateTNS();
        void calculatePower();
        void calculateUtil();
        void calculateArea();
        void calculateTNSbyDisplace();
        pair<double, double> getFFLF(Inst::FF_D *);
        pair<double, double> getFFUR(Inst::FF_D *);
        double calSlack4FF(Inst::FF_D *);
        double calDisplace4FF(Inst::FF_D *);
    };
}
#endif