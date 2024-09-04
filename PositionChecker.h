#ifndef POSITIONCHECKER_H
#define POSITIONCHECKER_H

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
    struct PlaceRegion_checker // record the information of connected union placement row
    {
        double xStart;
        double xEnd;
        double yStart;
        double yEnd;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };

    struct PlacementRow_checker
    {
        double x;
        double y;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };

    class PositionChecker
    {
    public:
        PositionChecker(Solver *parser) : _ptrSolver(parser), isInit(false) {};
        vector<pair<double, double>> findAvaiablePosition(size_t fflibID, pair<double, double> position, double radius); // return the available position for the ff lib in the region centered by position with radius "radius"

    private:
        // data part (need to be adjusted to your own data structure)
        Solver *_ptrSolver;
        vector<struct PlacementRow_checker> _Pid2PlaceRow; // the placement row ID to the placement row's pointer
        vector<Inst::Gate *> _id2gate;                     // please change the ptr to the type you use for gate

        // data part (don't need to change)
        bool isInit;
        vector<struct PlaceRegion_checker> _PlaceRegionSet;

        // catch for legalization
        vector<vector<vector<bool>>> _availPosTable; // if the grid point is available, the value would be true _availPosTable[place region][y index][x index]

        // auxilury function (need to be modified)
        pair<double, double> getGateLF(Inst::Gate *) const; // get the lower left coordinate of the gate
        pair<double, double> getGateUR(Inst::Gate *) const; // get the upper right coordinate of the gate
        void initPositionChecker();

        // auxilury function (don't need to be modified)
        void categorizePlaceRegion();
        void constructAvailableTable();
    };
}

#endif