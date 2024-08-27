#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "parser.h"
#include "solver.h"
#include "STAEngine.h"
#include "legalizer.h"
#include "GlobalPlacer.h"
#include "evaluator.h"

using namespace std;
int main(int argc, char *argv[])
{
    clock_t a, b;
    a = clock();
    static Parse::Parser parser;
    for (int i = 1; i < 2; ++i)
    {
        parser.readInput(argv[i]);
    }
    b = clock();
    cout << "Time for read file: " << (b - a) << " sec" << endl;
    for (int slice = 50; slice <= 300; slice += 10)
    {

        for (double low_weight = 0; low_weight <= 1; low_weight += 0.1)
        {
            std::cout << slice << " ," << low_weight << " start--------------------------" << std::endl;

            Solver::Solver solver;
            Solver::STAEngine STAEngine;
            Solver::legalizer legalizer;
            Solver::GlobalPlacer GlobalPlacer(&parser, &solver);
            Solver::evaluator evaluator;
            solver.setParserPtr(&parser);
            solver.setSTAEnginePtr(&STAEngine);
            solver.setLegalizerPtr(&legalizer);
            solver.setGlobalPlacerPtr(&GlobalPlacer);
            solver.setEvaluatorPtr(&evaluator);

            solver.initSolver();
            a = clock();
            cout << "Time for initialize solver: " << (a - b) << " sec" << endl;

            string filename = "tc2_0812_s" + std::to_string(slice) + "_low_" + std::to_string(low_weight);

            solver.solve_by_window(slice, low_weight, filename);

            std::cout << slice << " ," << low_weight << " is done--------------------------" << std::endl;
            // b = clock();
            // cout << "Time for solve: " << (b - a) << " sec" << endl;
            // solver.printOutput(argv[argc - 1]);
            // a = clock();
            // cout << "Time for output: " << (a - b) << " sec" << endl;
        }
    }

    return 0;
}
