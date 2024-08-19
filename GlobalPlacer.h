#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <string>

#include "solver.h"
#include "parser.h"

namespace Solver
{
    class GlobalPlacer
    {
    public:
        GlobalPlacer(Parse::Parser *parser, Solver *solver_Ptr) : _parser(parser), _solver_Ptr(solver_Ptr) {};
        // void plotPlacementResult(const std::string outfilename, bool isPrompt = false);
        // void plotSpecificNet(const std::string outfilename, bool isPrompt = false);
        // void plotCLKNet(const std::string outfilename, int clkNetID, bool isPrompt = false);
        void placementVisualization(const std::string outfilename, vector<size_t> FFIDs, vector<size_t> GateIDs, vector<size_t> NetIDs, bool isPrompt = false);

    private:
        Parse::Parser *_parser;
        Solver *_solver_Ptr;

        void plotBoxPLT(std::ofstream &stream, double x1, double y1, double x2, double y2);
        // void plotNet(ofstream &stream, int netID, bool connectFFPartOnly);
    };
}

#endif // GLOBALPLACER_H