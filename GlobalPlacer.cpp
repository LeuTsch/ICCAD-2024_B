#define _USE_MATH_DEFINES
#include "GlobalPlacer.h"
#include "inst.h"

#include <cstdio>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <sstream>

// input id should be the local id (remember minus the OFFSET)
void Solver::GlobalPlacer::placementVisualization(const std::string outfilename, vector<size_t> FFIDs, vector<size_t> GateIDs, vector<size_t> NetIDs, bool isPrompt)
{

    // Draw the die
    std::ofstream outfile(outfilename.c_str(), std::ios::out);
    outfile << " " << std::endl;
    outfile << "set size ratio 1" << std::endl;
    outfile << "set nokey" << std::endl;
    outfile << "set term qt size 2560,2560" << std::endl // qt for macos, it might use x11 in windows
            << std::endl;

    // lt 1 is type 1, lw 1 is width 1, the width of lw 5 is bigger than lw 1.
    // black for die bounding box, red for FF, green for gate.
    outfile << "plot[:][:] '-' w l lt 1 lw 1 lc rgb \"black\", '-' w l lt 1 lw 1 lc rgb \"red\", '-' w l lt 1 lw 1 lc rgb \"green\"" << std::endl
            << std::endl;
    outfile << "# bounding box" << std::endl;
    plotBoxPLT(outfile, _parser->_dieLx, _parser->_dieLy, _parser->_dieRx, _parser->_dieRy);
    outfile << "EOF" << std::endl;

    // Draw FFs
    // TODOs: Please use your FF instances.
    outfile << "# FFs" << std::endl
            << "0.00, 0.00" << std::endl
            << std::endl;
    for (size_t i = 0; i < FFIDs.size(); ++i)
    {

        Inst::FF_D &ff = _solver_Ptr->_FF_D_arr[FFIDs[i]];
        pair<double, double> pos = _solver_Ptr->getFFPosition(&ff);
        double x = pos.first;
        double y = pos.second;
        double w = _parser->_flipflopLib[ff.FF_type].Width, h = _parser->_flipflopLib[ff.FF_type].Hight;
        plotBoxPLT(outfile, x, y, (x + w), (y + h));
    }
    outfile << "EOF" << std::endl
            << std::endl;

    // Draw Gates
    // TODOs: Please use your Gate instances.
    if (GateIDs.size() > 0)
    {
        outfile << "# Gates" << std::endl
                << "0.00, 0.00" << std::endl
                << std::endl;
        for (size_t i = 0; i < GateIDs.size(); ++i)
        {
            Inst::Gate &gate = _solver_Ptr->_Gate_arr[GateIDs[i]];
            size_t Gatetype_id = gate.gate_type;
            double x = gate.getPosition().first, y = gate.getPosition().second;
            double w = _parser->_gateLib[Gatetype_id].Width, h = _parser->_gateLib[Gatetype_id].Hight;
            plotBoxPLT(outfile, x, y, (x + w), (y + h));
        }
        outfile << "EOF" << std::endl
                << std::endl;
    }

    // Draw Nets, the input of plotNet function is NetID.
    // TODOs: Please choose the kind of displaying of a net, by giving input

    /*
    if (NetIDs.size() > 0)
    {
        outfile << "# Net " << endl;
        for (size_t i = 0; i < NetIDs.size(); ++i)
        {
            plotNet(outfile, NetIDs[i], true);
        }
        outfile << "EOF" << endl;
    }
    */

    // Don't change here.
    outfile << "pause -1 'Press any key to close.'" << std::endl;
    outfile.close();

    if (isPrompt)
    {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd))
        {
            std::cout << "Fail to execute: \"" << cmd << "\"." << std::endl;
        }
    }
}

void Solver::GlobalPlacer::plotBoxPLT(std::ofstream &stream, double x1, double y1, double x2, double y2)
{
    stream << x1 << ", " << y1 << std::endl
           << x2 << ", " << y1 << std::endl
           << x2 << ", " << y2 << std::endl
           << x1 << ", " << y2 << std::endl
           << x1 << ", " << y1 << std::endl
           << std::endl;
}

/*
void GlobalPlacer::plotNet(ofstream &stream, int netID, bool connectFFPartOnly)
{
    // from point
    bool zeroQ = (_parser._Nets[netID]->_QFFID.size() == 0) ? true : false;
    bool zeroGOUT = (_parser._Nets[netID]->_OUTGATEID.size() == 0) ? true : false;
    bool zeroINPUT = (_parser._Nets[netID]->_INPUTID.size() == 0) ? true : false;
    // to point
    bool zeroD = (_parser._Nets[netID]->_DFFID.size() == 0) ? true : false;
    bool zeroGIN = (_parser._Nets[netID]->_INGATEID.size() == 0) ? true : false;
    bool zeroOUTPUT = (_parser._Nets[netID]->_OUTPUTID.size() == 0) ? true : false;

    // Only draw the part of net that connects FF D & Q pins.
    if (connectFFPartOnly)
    {
        if (!zeroQ && !zeroD)
        {
            for (int i = 0; i < _parser._Nets[netID]->_QFFID.size(); ++i)
            {
                for (int j = 0; j < _parser._Nets[netID]->_DFFID.size(); ++j)
                {
                    stream << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getY() << endl;
                    stream << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getY() << endl;
                    stream << endl;
                }
            }
        }
        return;
    }

    if (!zeroQ && !zeroD)
    {
        for (int i = 0; i < _parser._Nets[netID]->_QFFID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_DFFID.size(); ++j)
            {
                stream << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getY() << endl;
                stream << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroQ && !zeroGIN)
    {
        for (int i = 0; i < _parser._Nets[netID]->_QFFID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_INGATEID.size(); ++j)
            {
                stream << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getY() << endl;
                stream << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroQ && !zeroOUTPUT)
    {
        for (int i = 0; i < _parser._Nets[netID]->_QFFID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_OUTPUTID.size(); ++j)
            {
                stream << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_QFFID[i]]->getY() << endl;
                stream << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_X << " " << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_Y << endl;
                stream << endl;
            }
        }
    }

    if (!zeroGOUT && !zeroD)
    {
        for (int i = 0; i < _parser._Nets[netID]->_OUTGATEID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_DFFID.size(); ++j)
            {
                stream << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getY() << endl;
                stream << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroGOUT && !zeroGIN)
    {
        for (int i = 0; i < _parser._Nets[netID]->_OUTGATEID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_INGATEID.size(); ++j)
            {
                stream << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getY() << endl;
                stream << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroGOUT && !zeroOUTPUT)
    {
        for (int i = 0; i < _parser._Nets[netID]->_OUTGATEID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_OUTPUTID.size(); ++j)
            {
                stream << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_OUTGATEID[i]]->getY() << endl;
                stream << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_X << " " << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_Y << endl;
                stream << endl;
            }
        }
    }

    if (!zeroINPUT && !zeroD)
    {
        for (int i = 0; i < _parser._Nets[netID]->_INPUTID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_DFFID.size(); ++j)
            {
                stream << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_X << " " << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_Y << endl;
                stream << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getX() << " " << _parser._FFs[_parser._Nets[netID]->_DFFID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroINPUT && !zeroGIN)
    {
        for (int i = 0; i < _parser._Nets[netID]->_INPUTID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_INGATEID.size(); ++j)
            {
                stream << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_X << " " << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_Y << endl;
                stream << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getX() << " " << _parser._Gates[_parser._Nets[netID]->_INGATEID[j]]->getY() << endl;
                stream << endl;
            }
        }
    }
    if (!zeroINPUT && !zeroOUTPUT)
    {
        for (int i = 0; i < _parser._Nets[netID]->_INPUTID.size(); ++i)
        {
            for (int j = 0; j < _parser._Nets[netID]->_OUTPUTID.size(); ++j)
            {
                stream << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_X << " " << _parser._Inputs[_parser._Nets[netID]->_INPUTID[i]]->_Y << endl;
                stream << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_X << " " << _parser._Outputs[_parser._Nets[netID]->_OUTPUTID[j]]->_Y << endl;
                stream << endl;
            }
        }
    }
}
*/
