#include "parser.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;
void Parse::Parser::readInput(string file)
{
    fstream input;
    cout << "start to read file: " << file << endl;
    input.open(file, ios::in);
    if (!input)
    {
        cerr << "Can't open file with the path: " << file << "\n";
        exit(1);
    }

    // start to read input
    string keyw;
    while (!input.eof())
    {
        input >> keyw;
        if (keyw == "Alpha")
        {
            input >> _alpha;
        }
        else if (keyw == "Beta")
        {
            input >> _beta;
        }
        else if (keyw == "Gamma")
        {
            input >> _gamma;
        }
        else if (keyw == "Delta")
        {
            input >> _delta;
        }
        else if (keyw == "DieSize")
        {
            input >> _dieLx;
            input >> _dieLy;
            input >> _dieRx;
            input >> _dieRy;
        }
        else if (keyw == "NumInput")
        {
            int count;
            input >> count;
            _inputName.reserve(count);
            _inputCrdnate.reserve(count);
            for (size_t i = 0; i < count; i++)
            {
                input >> keyw;
                input >> keyw;
                _inputName.push_back(keyw);
                pair<coord, coord> xyCoord;
                input >> xyCoord.first;
                input >> xyCoord.second;
                _inputCrdnate.push_back(xyCoord);
            }
        }
        else if (keyw == "NumOutput")
        {
            int count;
            input >> count;
            _outputName.reserve(count);
            _outputCrdnate.reserve(count);
            for (int i = 0; i < count; i++)
            {
                input >> keyw;
                input >> keyw;
                _outputName.push_back(keyw);
                pair<coord, coord> xyCoord;
                input >> xyCoord.first;
                input >> xyCoord.second;
                _outputCrdnate.push_back(xyCoord);
            }
        }
        else if (keyw == "FlipFlop")
        {
            struct FLIPFLOP ff;
            input >> ff.Bit;
            input >> ff.Name;
            input >> ff.Width;
            input >> ff.Hight;
            int count;
            input >> count;
            for (int i = 0; i < count; i++)
            {
                input >> keyw;
                input >> keyw;
                ff.PinName.push_back(keyw);
                pair<coord, coord> xyCoord;
                input >> xyCoord.first;
                input >> xyCoord.second;
                ff.PinCrdnate.push_back(xyCoord);
            }
            _flipflopLib.push_back(ff);
        }
        else if (keyw == "Gate")
        {
            struct FLIPFLOP ff;
            ff.Bit = 0;
            input >> ff.Name;
            input >> ff.Width;
            input >> ff.Hight;
            int count;
            input >> count;
            for (int i = 0; i < count; i++)
            {
                input >> keyw;
                input >> keyw;
                ff.PinName.push_back(keyw);
                pair<coord, coord> xyCoord;
                input >> xyCoord.first;
                input >> xyCoord.second;
                ff.PinCrdnate.push_back(xyCoord);
            }
            _gateLib.push_back(ff);
        }
        else if (keyw == "NumInstances")
        {
            int count;
            input >> count;
            _instList.reserve(count);
            for (int i = 0; i < count; i++)
            {
                struct INST ff;
                input >> keyw;
                input >> ff.Name;
                input >> ff.Type;
                input >> ff.x;
                input >> ff.y;
                _instList.push_back(ff);
            }
        }
        else if (keyw == "NumNets")
        {
            int count;
            input >> count;
            _netName.reserve(count);
            _netPin.reserve(count);
            for (int i = 0; i < count; i++)
            {
                input >> keyw;
                input >> keyw;
                _netName.push_back(keyw);
                int n;
                input >> n;
                vector<string> netPin;
                netPin.reserve(n);
                for (size_t j = 0; j < n; j++)
                {
                    input >> keyw;
                    input >> keyw;
                    netPin.push_back(keyw);
                }
                _netPin.push_back(netPin);
            }
        }
        else if (keyw == "BinWidth")
        {
            input >> _binWidth;
        }
        else if (keyw == "BinHeight")
        {
            input >> _binHeight;
        }
        else if (keyw == "BinMaxUtil")
        {
            input >> _binMaxUtil;
        }
        else if (keyw == "PlacementRows")
        {
            struct PLACEROW placeRow;
            input >> placeRow.x;
            input >> placeRow.y;
            input >> placeRow.siteWidth;
            input >> placeRow.siteHight;
            input >> placeRow.totalNumOfSites;
            _placeRow.push_back(placeRow);
        }
        else if (keyw == "DisplacementDelay")
        {
            input >> _displaceDelay;
        }
        else if (keyw == "QpinDelay")
        {
            input >> keyw;
            for (size_t i = 0; i < _flipflopLib.size(); i++)
            {
                if (_flipflopLib[i].Name == keyw)
                {
                    input >> _flipflopLib[i].PinDelay;
                    break;
                }
            }
        }
        else if (keyw == "TimingSlack")
        {
            pair<string, int> PinSlack;
            input >> keyw;
            PinSlack.first = keyw + "/";
            input >> keyw;
            PinSlack.first += keyw;
            input >> PinSlack.second;
            _timeSlack.push_back(PinSlack);
        }
        else if (keyw == "GatePower")
        {
            input >> keyw;
            for (size_t i = 0; i < _flipflopLib.size(); i++)
            {
                if (_flipflopLib[i].Name == keyw)
                {
                    input >> _flipflopLib[i].Power;
                    break;
                }
            }
        }
    }
    input.close();
};