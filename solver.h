#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>

#include "inst.h"
#include "parser.h"

typedef size_t FF_D_ID;
typedef size_t FF_Q_ID;
typedef size_t Gate_ID;
typedef size_t PIO_ID;

using std::pair;
using std::string;
using std::vector;

namespace Solver
{
    class Solver
    {
    public:
        Solver() : _ptr_Parser(nullptr) {}
        void setParserPtr(Parse::Parser *ptr) { _ptr_Parser = ptr; };
        void initSolver();
        void solve();
        void printOutput();

    private:
        // data part
        Parse::Parser *_ptr_Parser;
        size_t _GATE_OFFSET;
        size_t _FF_D_OFFSET;
        size_t _FF_Q_OFFSET;
        vector<Inst::Inst *> _ID_to_instance;
        vector<Inst::PIO> _PIO_arr;
        vector<Inst::Gate> _Gate_arr;
        vector<Inst::Gate *> _GID_to_ptrGate_map;
        vector<Inst::FF_D> _FF_D_arr;
        vector<Inst::FF_Q> _FF_Q_arr;
        vector<string> _Name_to_ID;
        vector<vector<size_t>> _NetList;
        vector<struct PlacementRow> _PlaceRow;
        vector<string> _ClkList;
        vector<vector<vector<size_t>>> _Gate_in_Bin; // inside is the position in _Gate_arr, rather than ID of gate
        vector<vector<vector<size_t>>> _FF_in_Bin;   // only record FF_D
        vector<vector<vector<size_t>>> _PlaceRow_in_Bin;

        // function part
        bool mbffCluster();
        bool prePlace(vector<FF_D_ID>, size_t, pair<double, double>);
        void legalize();

        // interface
        pair<double, double> findPinPosition(const size_t &) const;
        size_t name2ID(string &) const;

        // the function should only be called in initialization
        void findFanin(const FF_D_ID &);
        void findFaninRecur(const FF_D_ID &, const Gate_ID &);
        void findFanout(const FF_Q_ID &);
        void findFanoutRecur(const FF_Q_ID &, const Gate_ID &);
    };

    struct PlacementRow
    {
        double x;
        double y;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };
}
#endif