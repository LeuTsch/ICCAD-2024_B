#ifndef STAENGINE_H
#define STAENGINE_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <list>
#include <map>

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
        STAEngine() : _ptrSolver(nullptr) {};
        double getDistance(const size_t &, const size_t &);                            // we would expect the first ID is gate's in, the other is ID of out
        void initEngine(const vector<vector<size_t>> &, const vector<Inst::Inst *> &); // the input should be a set of Net and an ID to pin mapping

        void setSolverPtr(Solver *ptr) { _ptrSolver = ptr; }; // change to the class you may need for auxilury data

        ////// should not be modified //////
        vector<pair<size_t, std::map<size_t, double>>> _distanceList;
        ///////////////////////////////////

    private:
        // data part
        vector<vector<size_t>> _InPinList;        // record the In Pin in the i-th Net
        vector<vector<size_t>> _OutPinList;       // record the Out Pin in the i-th Net
        std::map<size_t, int> _InPin2PositionMap; // record the in pin ID to the position in _distanceList
        std::map<size_t, bool> _isOutPinCritical; // record whether an out pin is connected to FF
        Solver *_ptrSolver;
        vector<bool> _checkList; // check for whether the pin got explored

        // auxilury function
        void propagateForward(list<size_t>, const size_t &, double, list<double>);

        // function need to be rewrited with the data structures used
        bool isInPin(const size_t &) const;                  // see whether it is an input pin of gate
        bool isOutPin(const size_t &) const;                 // see whether it is an output pin of gate
        bool isFF_D(const size_t &) const;                   // see whether it is an D pin of FF
        list<size_t> getInPinRelated(const size_t &) const;  // get the input pin in the same gate (may include itself)
        list<size_t> getOutPinRelated(const size_t &) const; // get the output pin in the same gate (may include itself)
        size_t getRelatedNet(const size_t &) const;          // get the NetID relate to the pin
        pair<double, double> getPinPosition(const size_t &) const;
    };
}

#endif