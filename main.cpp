#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "parser.h"
#include "solver.h"
#include "STAEngine.h"

using namespace std;
int main(int argc, char *argv[])
{
    clock_t a, b;
    a = clock();
    static Parse::Parser parser;
    for (int i = 1; i < argc; ++i)
    {

        // read the input file
        // cout << "start to read file: " << argv[i] << endl;
        parser.readInput(argv[i]);
    }
    b = clock();
    cout << "Time for read file: " << (b - a) << " sec" << endl;

    Solver::Solver solver;
    solver.setParserPtr(&parser);
    solver.initSolver();
    a = clock();
    cout << "Time for initialize solver: " << (a - b) << " sec" << endl;
    solver.test();
    b = clock();
    cout << "Time for legal and output: " << (b - a) << " sec" << endl;

    return 0;
}
