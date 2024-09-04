#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "parser.h"
#include "solver.h"
#include "STAEngine.h"
#include "legalizer.h"
#include "evaluator.h"
#include "GlobalPlacer.h"
#include "PositionChecker.h"

using namespace std;
int main(int argc, char *argv[])
{
    clock_t a, b;
    a = clock();
    static Parse::Parser parser;
    for (int i = 1; i < 2; ++i)
    {

        // read the input file
        // cout << "start to read file: " << argv[i] << endl;
        parser.readInput(argv[i]);
    }
    b = clock();
    cout << "Time for read file: " << (b - a) << " sec" << endl;

    Solver::Solver solver;
    Solver::STAEngine STAEngine;
    Solver::legalizer legalizer;
    Solver::evaluator evaluator;
    Solver::GlobalPlacer GlobalPlacer(&parser, &solver);
    Solver::PositionChecker PositionChecker(&solver);
    solver.setParserPtr(&parser);
    solver.setSTAEnginePtr(&STAEngine);
    solver.setLegalizerPtr(&legalizer);
    solver.setEvaluatorPtr(&evaluator);
    solver.setGlobalPlacerPtr(&GlobalPlacer);
    solver.setPositionCheckerPtr(&PositionChecker);

    solver.initSolver();
    a = clock();
    cout << "Time for initialize solver: " << (a - b) << " sec" << endl;
    solver.solve_by_window();
    // solver.test();
    b = clock();
    cout << "Time for solve: " << (b - a) << " sec" << endl;
    solver.printOutput(argv[argc - 1]);
    a = clock();
    cout << "Time for output: " << (a - b) << " sec" << endl;

    return 0;
}
