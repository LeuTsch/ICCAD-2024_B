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
    struct PlacementRow
    {
        double x;
        double y;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };

    class Solver
    {
    public:
        Solver() : _ptr_Parser(nullptr) {}
        void setParserPtr(Parse::Parser *ptr) { _ptr_Parser = ptr; };
        void initSolver();
        void solve();
        void printOutput();
        void test(); // do whatever test you want here

        // friend class declaration
        friend class STAEngine;

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
        vector<size_t> prePlace(vector<FF_D_ID>, pair<double, double>); // return the FF_D_ID be grouped
        void slackDistribute(const double);                             // should be called before mbff clustering
        void legalize();

        // interface
        pair<double, double> findPinPosition(const size_t &) const;
        size_t name2ID(string &) const;
        vector<pair<double, double>> getAdjacentPinPosition(size_t &) const; // the input should be ID of Q, or D pin. Return the position of related pin

        // auxilary function
        vector<struct PlacementRow> getPlacementRow() const { return _PlaceRow; };
        bool placementRowIsUniform(vector<struct PlacementRow> &) const;
        size_t getPlaceRowNum() const { return _PlaceRow.size(); };
        double getSiteHeight() const { return _PlaceRow[0].siteHight; };
        double getSiteWidth() const { return _PlaceRow[0].siteWidth; };
        double getTotalSiteNum() const { return _PlaceRow[0].totalNumOfSites; };
        double getLowerLeftX() const;
        double getLowerLeftY() const;
        pair<double, double> getGateLF(Inst::Gate *) const;
        pair<double, double> getGateUR(Inst::Gate *) const;
        pair<double, double> getFFPosition(Inst::FF_D *) const;
        double getFFWidth(Inst::FF_D *) const;
        double getFFHeight(Inst::FF_D *) const;
        void setFFPosition(Inst::FF_D *, pair<double, double> &);
        vector<size_t> getGroupMem(Inst::FF_D *) const;

        // the function should only be called in initialization
        void findFanin(const FF_D_ID &);
        void findFaninRecur(const FF_D_ID &, const Gate_ID &);
        void findFanout(const FF_Q_ID &);
        void findFanoutRecur(const FF_Q_ID &, const Gate_ID &);
    };
}
#endif