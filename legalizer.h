#ifndef LEGALIZER_H
#define LEGALIZER_H

#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <list>

#include "inst.h"
#include "parser.h"
#include "solver.h"

using std::list;

namespace Solver
{
    struct PlaceRegion // record the information of connected union placement row
    {
        double xStart;
        double xEnd;
        double yStart;
        double yEnd;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };

    class legalizer
    {
    public:
        legalizer() : _ptrSolver(nullptr), isInit(false) {};
        bool legalize(); // return true if succussfully legalize

        void setSolverPtr(Solver *ptr) { _ptrSolver = ptr; }; // change to the class you may need for auxilury data

    private:
        // data part (need to be adjusted to your own data structure)
        Solver *_ptrSolver;
        vector<struct PlacementRow *> _Pid2PlaceRowPtr; // the placement row ID to the placement row's pointer
        vector<Inst::FF_D *> _id2ffPtr;                 // please change the ptr to the type you use for the ff
        vector<Inst::Gate *> _id2gate;                  // please change the ptr to the type you use for gate

        // data part (don't need to change)
        bool isInit;
        vector<struct PlaceRegion> _PlaceRegionSet;
        vector<list<size_t>> _FFinPlaceRegion; // record the id of ff in placeRegion
        list<size_t> _legalizedFFID;           // record the legalized ff id
        // catch for legalization
        vector<vector<bool>> _availPosTable;                 // if the grid point is available, the value would be true
        vector<vector<pair<size_t, int>>> _listWait4Legal;   // catch for the FF wait for legalize (should not be modified during the legalize)
        vector<vector<pair<size_t, int>>> _legalist;         // double copy for _listWait4Legal, can be modified
        vector<vector<pair<int, int>>> _finalSolution;       // catch for record the solution | pair<y index, x index>
        vector<vector<bool>> DPtable;                        // DPtable[ i-th object in legalist ][ position ]
        vector<int> heightConstraint;                        // record the available height for x position
        vector<vector<vector<pair<int, int>>>> solutionList; // pair<int y, int x>
        vector<vector<unsigned int>> totalDisplace;          // totalDisplace[# of object][last position]

        // auxilury function (need to be modified)
        void initLegalizer();                                                                       // initialize the basic data structure of legalizer
        pair<double, double> getLowerLeftPlaceRow(struct PlacementRow *) const;                     // return the lower left of the placement row
        double getSiteWidth(struct PlacementRow *ptr) const { return ptr->siteWidth; };             // return the siteWidth of placement row
        double getSiteHeight(struct PlacementRow *ptr) const { return ptr->siteHight; };            // return the siteHeight of placement row
        double getTotalNumOfSites(struct PlacementRow *ptr) const { return ptr->totalNumOfSites; }; // return the TotalNumOfSites of placement row
        pair<double, double> getFFPosition(Inst::FF_D *) const;                                     // return the lower left position of FF
        vector<size_t> getGroupMem(Inst::FF_D *) const;                                             // return the ID that group with the FF
        double getFFWidth(Inst::FF_D *) const;                                                      // return the Width of FF
        double getFFHeight(Inst::FF_D *) const;                                                     // return the Height of FF
        pair<double, double> getGateLF(Inst::Gate *) const;                                         // get the lower left coordinate of the gate
        pair<double, double> getGateUR(Inst::Gate *) const;                                         // get the upper right coordinate of the gate
        void setFFPosition(Inst::FF_D *, pair<double, double> &);                                   // As the name

        // auxilury function (don't need to be modified)
        void releaseCatch(); // releas the catch occupied during the process of legalization
        void categorizePlaceRegion();
        void categorizeFF();                                          // categorize FF into nearest placeRegion
        bool legalRegion(const list<size_t> &, struct PlaceRegion *); // return true if succussfully legalize
    };
}

#endif