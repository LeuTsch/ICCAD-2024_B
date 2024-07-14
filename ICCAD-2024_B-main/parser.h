#ifndef PARSER_H
#define PARSER_H

#include <string>
#include <vector>
#include <utility>

typedef double coord;

using std::pair;
using std::string;
using std::vector;

namespace Parse
{
    class Parser
    {
    public:
        Parser() {}
        void readInput(string file);
        double getAlpha() const { return _alpha; };
        double getBeta() const { return _beta; };
        vector<string> getInputName() const { return _inputName; };
        vector<string> getOutputName() const { return _outputName; };
        vector<pair<coord, coord> > getInputCrdnate() const { return _inputCrdnate; };
        vector<pair<coord, coord> > getOutputCrdnate() const { return _outputCrdnate; };
        vector<struct FLIPFLOP> getFlipflopLib() const { return _flipflopLib; };
        vector<struct INST> getInstList() const { return _instList; };
        vector<string> getNetName() const { return _netName; };
        vector<vector<string> > getNetPin() const { return _netPin; };
        vector<pair<string, double> > getTimeSlack() const { return _timeSlack; };

        // data part
        double _alpha;
        double _beta;
        double _gamma;
        double _delta;
        double _lambda;
        long _dieLx;
        long _dieLy;
        long _dieRx;
        long _dieRy;
        vector<string> _inputName;
        vector<string> _outputName;
        vector<pair<coord, coord> > _inputCrdnate;
        vector<pair<coord, coord> > _outputCrdnate;
        vector<struct FLIPFLOP> _flipflopLib;
        vector<struct FLIPFLOP> _gateLib;
        vector<struct INST> _instList;
        vector<string> _netName;
        vector<vector<string> > _netPin;
        int _binWidth;
        int _binHeight;
        double _binMaxUtil;
        vector<struct PLACEROW> _placeRow;
        double _displaceDelay;
        vector<pair<string, double> > _timeSlack;
    };

    struct FLIPFLOP
    {
        string Name;
        int Bit;
        double Width;
        double Hight;
        double Power;
        double PinDelay;
        vector<string> PinName;
        vector<pair<coord, coord> > PinCrdnate;
    };

    struct INST
    {
        string Name;
        string Type;
        coord x;
        coord y;
    };

    struct PLACEROW
    {
        double x;
        double y;
        double siteWidth;
        double siteHight;
        double totalNumOfSites;
    };
}
#endif