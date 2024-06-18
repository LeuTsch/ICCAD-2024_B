#include "parser.h"
#include <iostream>
#include <string>
#include <vector>
#include <ctime>

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

    vector<struct Parse::FLIPFLOP> FlipflopLib(parser.getFlipflopLib());
    a = clock();
    cout << "Time for copy the lib: " << (a - b) << " sec" << endl;
    cout << "alpha = " << parser.getAlpha() << endl;
    cout << "beta = " << parser.getBeta() << endl;
    cout << "flipflop size: " << parser.getInputName().size() << endl;
    for (size_t i = 0; i < FlipflopLib.size(); i++)
    {
        cout << "Name:" << FlipflopLib[i].Name << endl;
        cout << "Bit:" << (FlipflopLib)[i].Bit << endl;
        cout << "Width:" << (FlipflopLib)[i].Width << endl;
        cout << "Hight:" << (FlipflopLib)[i].Hight << endl;
        cout << "Power:" << (FlipflopLib)[i].Power << endl;
        cout << "PinDelay:" << (FlipflopLib)[i].PinDelay << endl;
        for (size_t j = 0; j < (FlipflopLib)[i].PinName.size(); j++)
        {
            cout << (FlipflopLib)[i].PinName[j] << ": " << (FlipflopLib)[i].PinCrdnate[j].first << " " << (FlipflopLib)[i].PinCrdnate[j].second << endl;
        }
    }
    b = clock();
    cout << "Time for output the lib: " << (b - a) << " sec" << endl;

    return 0;
}
