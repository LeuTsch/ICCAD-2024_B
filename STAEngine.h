#ifndef STAENGINE_H
#define STAENGINE_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <list>

#include "inst.h"
#include "parser.h"
#include "solver.h"

using std::list;
using std::pair;
using std::string;
using std::vector;

namespace Solver
{
    class STAEngine
    {
    public:
        STAEngine() : _ptrSolver(nullptr){};
        double getDistance(size_t &, size_t &) const;                                  // we would expect the first ID is gate's in, the other is ID of out
        void initEngine(const vector<vector<size_t>> &, const vector<Inst::Inst *> &); // the input should be a set of Net and an ID to pin mapping

        void setSolverPtr(Solver *ptr) { _ptrSolver = ptr; }; // change to the class you may need for auxilury data

    private:
        // data part
        vector<pair<size_t, vector<pair<size_t, double>>>> _distanceList;
        vector<vector<size_t>> _InPinList;  // record the In Pin in the i-th Net
        vector<vector<size_t>> _OutPinList; // record the Out Pin in the i-th Net
        Solver *_ptrSolver;

        // auxilury function
        void propagateForward(const size_t &, const size_t &, double);

        // function need to be rewrited with the data structures used
        bool isInPin(const size_t &) const;                  // see whether it is an input pin of gate
        bool isOutPin(const size_t &) const;                 // see whether it is an output pin of gate
        list<size_t> getInPinRelated(const size_t &) const;  // get the input pin in the same gate of an output pin
        list<size_t> getOutPinRelated(const size_t &) const; // get the output pin in the same gate
        size_t getRelatedNet(const size_t &) const;          // get the NetID relate to the pin
        pair<double, double> getPinPosition(const size_t &) const;
    };
}

#endif