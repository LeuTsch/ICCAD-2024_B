#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <map>

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
    class STAEngine;    // forward declaration
    class legalizer;    // forward declaration
    class GlobalPlacer; ///////////////////////////////

    class Solver
    {
    public:
        Solver() : _ptr_Parser(nullptr), _ptr_STAEngine(nullptr) {}
        void setParserPtr(Parse::Parser *ptr) { _ptr_Parser = ptr; };
        void setSTAEnginePtr(STAEngine *ptr) { _ptr_STAEngine = ptr; };
        void setLegalizerPtr(legalizer *ptr) { _ptr_legalizer = ptr; };
        ///////////////////////////////////////////////
        void setGlobalPlacerPtr(GlobalPlacer *ptr) { _ptr_GlobalPlacer = ptr; };
        ///////////////////////////////////////////////
        void initSolver();
        void solve();
        void printOutput(const string &);
        void test(); // do whatever test you want here

        /////////////////////// victor's part
        void solve_initbuild();
        void solve_findfeasible();
        vector<size_t> solve_findmaximal(const vector<size_t> &, size_t, pair<double, double> &, pair<double, double> &);
        void feasible_cal(const vector<size_t> &);
        void drawpic(const string &); ////////////////////////////
        ///////////////////////

        // friend class declaration
        friend class GlobalPlacer; ////////////////////////////////////////////////

        friend class STAEngine;
        friend class legalizer;

    private:
        // data part
        Parse::Parser *_ptr_Parser;
        STAEngine *_ptr_STAEngine;
        legalizer *_ptr_legalizer;
        ///////////////////////////////////////
        GlobalPlacer *_ptr_GlobalPlacer;
        ///////////////////////////////////////
        size_t _GATE_OFFSET;
        size_t _FF_D_OFFSET;
        size_t _FF_Q_OFFSET;
        vector<Inst::Inst *> _ID_to_instance;
        vector<Inst::PIO> _PIO_arr;
        vector<Inst::Gate> _Gate_arr;
        vector<Inst::Gate *> _GID_to_ptrGate_map;
        vector<Inst::FF_D> _FF_D_arr;
        vector<Inst::FF_Q> _FF_Q_arr;
        std::map<string, size_t> _Name_to_ID;                      // map pin's name(ex: C8763/D) to its ID
        std::map<string, pair<size_t, size_t>> _initFFName2arrPos; // map instace(ex: C87/ ) name to its pin position in FF_D_arr <start, end>
        vector<vector<size_t>> _NetList;
        vector<struct PlacementRow> _PlaceRow;
        vector<string> _ClkList;
        vector<vector<vector<size_t>>> _Gate_in_Bin; // inside is the position in _Gate_arr, rather than ID of gate
        vector<vector<vector<size_t>>> _FF_in_Bin;   // only record FF_D
        vector<vector<vector<size_t>>> _PlaceRow_in_Bin;

        // function part
        vector<size_t> prePlace(const vector<FF_D_ID> &, size_t, pair<double, double>); // return the FF_D_ID be grouped
        vector<pair<size_t, double>> getSlack2ConnectedFF(const size_t &);
        // if the the connected FF has grouped, return the 100% of remaining slack. Otherwise, it would be 50%
        void legalize();

        ////////////// no need to use it more(if you would use the slack in Inst::FF_D or Inst::FF_Q)
        void slackDistribute(const double); // no need to use it more
        /////////////

        // interface
        pair<double, double> findPinPosition(const size_t &) const;
        size_t name2ID(string &) const;
        vector<pair<double, double>> getAdjacentPinPosition(size_t &) const; // the input should be ID of Q, or D pin. Return the position of related pin

        // auxilary function
        void findMaxSlack(); // calculate the max capacity(the critical path + slack) for all FF_D
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
        vector<size_t> getGroupMem(Inst::FF_D *) const; // return the index in FF_D_arr related to this ff

        // the function should only be called in initialization
        void findFaninout4all(const vector<pair<size_t, std::map<size_t, double>>> &);
        vector<pair<size_t, vector<size_t>>> inPin2Out;
        vector<pair<size_t, vector<size_t>>> outPin2In;
        vector<bool> isInPinCritical;
        vector<bool> isOutPinCritical;
        std::map<size_t, vector<size_t>> FFConnect2CriticalPin;
        std::map<size_t, int> _OutPin2PositionMap; // record the in pin ID to the position in outPin2In
        std::map<size_t, int> _InPin2PositionMap;  // record the in pin ID to the position in inPin2Out
        bool isInPin(const size_t &) const;        // see whether it is an input pin of gate
        bool isOutPin(const size_t &) const;       // see whether it is an output pin of gate
        vector<vector<size_t>> _InPinList;         // record the In Pin in the i-th Net
        vector<vector<size_t>> _OutPinList;        // record the Out Pin in the i-th Net
    };
}
#endif