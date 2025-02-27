#include <string>
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <map>
#include <list>
#include <fstream>
#include <iomanip>
#include <ctime>

#include "inst.h"
#include "solver.h"
#include "STAEngine.h"
#include "legalizer.h"
#include "GlobalPlacer.h"
#include "evaluator.h"
#include "PositionChecker.h"

using std::list;
using std::pair;
using std::string;
using std::vector;

void Solver::Solver::initSolver()
{
    /*
    TODO: use the data in Parser to construct the data structure
    */
    assert(_ptr_Parser != nullptr);

    std::cout << "start to initialize solver" << std::endl;
    // step1: realize the instance for PIO, ff, and Gate,and assign ID to them
    // realize PIO
    std::cout << "step1: realize the instance for PIO, ff, and Gate,and assign ID" << std::endl;
    _PIO_arr.reserve(_ptr_Parser->_inputName.size() + _ptr_Parser->_outputName.size());
    for (size_t i = 0; i < _ptr_Parser->_inputName.size(); i++)
    {
        Inst::PIO pi(_ptr_Parser->_inputName[i], _ptr_Parser->_inputCrdnate[i]);
        _PIO_arr.push_back(pi);
    }
    for (size_t i = 0; i < _ptr_Parser->_outputName.size(); i++)
    {
        Inst::PIO po(_ptr_Parser->_outputName[i], _ptr_Parser->_outputCrdnate[i]);
        _PIO_arr.push_back(po);
    }
    // realize Gate and FF
    _FF_D_arr.reserve(_ptr_Parser->_instList.size());
    _FF_Q_arr.reserve(_ptr_Parser->_instList.size());
    _Gate_arr.reserve(_ptr_Parser->_instList.size());
    _initFFName2arrPos.clear();
    for (const auto &instance : _ptr_Parser->_instList)
    {
        bool isFF = false;
        // check the instance is FF or Gate
        for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
        {
            // is FF
            if (_ptr_Parser->_flipflopLib[i].Name == instance.Type)
            {
                string nameFFbit = instance.Name + "/";
                pair<size_t, size_t> rangeInArr;
                rangeInArr.first = _FF_D_arr.size();
                rangeInArr.second = rangeInArr.first + _ptr_Parser->_flipflopLib[i].Bit - 1;
                _initFFName2arrPos[nameFFbit] = rangeInArr;
                // for every pin (D/Q), we create an instance for it
                for (size_t k = 0; k < _ptr_Parser->_flipflopLib[i].PinName.size(); k++)
                {
                    if (_ptr_Parser->_flipflopLib[i].PinName[k][0] == 'D')
                    {
                        string QpinName = _ptr_Parser->_flipflopLib[i].PinName[k];
                        QpinName[0] = 'Q';
                        pair<double, double> position;
                        position.first = instance.x + _ptr_Parser->_flipflopLib[i].PinCrdnate[k].first;
                        position.second = instance.y + _ptr_Parser->_flipflopLib[i].PinCrdnate[k].second;
                        Inst::FF_D FF_D(nameFFbit + _ptr_Parser->_flipflopLib[i].PinName[k], position);
                        FF_D.FF_type = i;
                        FF_D.OriFF_type = i;
                        _FF_D_arr.push_back(FF_D);
                        // search for correspomding Q
                        for (size_t j = 0; j < _ptr_Parser->_flipflopLib[i].PinName.size(); j++)
                        {
                            if (_ptr_Parser->_flipflopLib[i].PinName[j] == QpinName)
                            {
                                position.first = instance.x + _ptr_Parser->_flipflopLib[i].PinCrdnate[j].first;
                                position.second = instance.y + _ptr_Parser->_flipflopLib[i].PinCrdnate[j].second;
                                Inst::FF_Q FF_Q(nameFFbit + QpinName, position);
                                FF_Q.FF_type = i;
                                _FF_Q_arr.push_back(FF_Q);
                            }
                        }
                    }
                }
                isFF = true;
                break;
            }
        }
        if (isFF)
            continue;
        for (size_t i = 0; i < _ptr_Parser->_gateLib.size(); i++)
        {
            // is Gate
            if (_ptr_Parser->_gateLib[i].Name == instance.Type)
            {
                pair<double, double> position;
                position.first = instance.x;
                position.second = instance.y;
                Inst::Gate gate(instance.Name, position);
                gate.gate_type = i;
                for (size_t k = 0; k < _ptr_Parser->_gateLib[i].PinName.size(); k++)
                {
                    position.first = instance.x + _ptr_Parser->_gateLib[i].PinCrdnate[k].first;
                    position.second = instance.y + _ptr_Parser->_gateLib[i].PinCrdnate[k].second;
                    gate.pinPosition.push_back(position);
                    gate.pinName.push_back(_ptr_Parser->_gateLib[i].PinName[k]);
                }
                _Gate_arr.push_back(gate);
                isFF = true;
                break;
            }
        }
        if (!isFF)
            std::cout << "the instance " << instance.Name << "'s type is not defined in either FF or Gate library" << std::endl;
    }

    // start to assign Id
    _ID_to_instance.clear();
    _ID_to_instance.reserve(_PIO_arr.size() + 5 * _Gate_arr.size() + _FF_D_arr.size() + _FF_Q_arr.size());
    _Name_to_ID.clear();
    _GID_to_ptrGate_map.clear();
    _GID_to_ptrGate_map.reserve(10 * _Gate_arr.size());
    // assign Id for PIO
    for (size_t i = 0; i < _PIO_arr.size(); i++)
    {
        Inst::Inst *ptr = &_PIO_arr[i];
        _Name_to_ID[ptr->getName()] = i;
        _ID_to_instance.push_back(ptr);
    }
    _GATE_OFFSET = _ID_to_instance.size();
    // assign Id for gate
    for (size_t i = 0; i < _Gate_arr.size(); i++)
    {
        Inst::Inst *ptr = &_Gate_arr[i];
        Inst::Gate *ptr_gate = &_Gate_arr[i];
        ptr_gate->PIN_OFFSET = _ID_to_instance.size();
        for (size_t j = 0; j < _Gate_arr[i].pinName.size(); j++)
        {
            _Name_to_ID[ptr->getName() + "/" + _Gate_arr[i].pinName[j]] = j + ptr_gate->PIN_OFFSET;
            _ID_to_instance.push_back(ptr);
            _GID_to_ptrGate_map.push_back(ptr_gate);
        }
    }
    _FF_D_OFFSET = _ID_to_instance.size();
    // assign Id for FF_D
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        Inst::Inst *ptr = &_FF_D_arr[i];
        _FF_D_arr[i].grouped_member.push_back(_FF_D_OFFSET + i);
        _Name_to_ID[_FF_D_arr[i].getOriName()] = i + _FF_D_OFFSET;
        _ID_to_instance.push_back(ptr);
    }
    _FF_Q_OFFSET = _ID_to_instance.size();
    // assign Id for FF_Q
    for (size_t i = 0; i < _FF_Q_arr.size(); i++)
    {
        Inst::Inst *ptr = &_FF_Q_arr[i];
        _Name_to_ID[_FF_Q_arr[i].getOriName()] = i + _FF_Q_OFFSET;
        _ID_to_instance.push_back(ptr);
    }

    // step2: establish Net
    std::cout << "step2: establish Net" << std::endl;
    for (const auto &netPinList : _ptr_Parser->_netPin)
    {
        bool isCLKNet = false;
        for (const auto &PinName : netPinList)
        {
            if (PinName.find("CLK") != std::string::npos || PinName.find("clk") != std::string::npos)
            {
                isCLKNet = true;
                break;
            }
        }
        if (isCLKNet)
        {
            bool hasSource = false;
            for (const auto &PinName : netPinList)
            {
                if (PinName.find("OUT") != std::string::npos || PinName.find("/") == std::string::npos || PinName.find("out") != std::string::npos)
                {
                    if (PinName.find("OUTPUT") != std::string::npos)
                    {
                        continue;
                    }
                    hasSource = true;
                    _ClkList.push_back(PinName);
                    break;
                }
            }
            if (!hasSource)
            {
                _ClkList.push_back("_NoSource");
            }
            for (const auto &PinName : netPinList)
            {
                if (PinName.find("/CLK") != std::string::npos || PinName.find("/clk") != std::string::npos)
                {
                    string insName = PinName.substr(0, PinName.find("/") + 1);
                    pair<size_t, size_t> posArr = _initFFName2arrPos.at(insName);
                    for (size_t i = posArr.first; i <= posArr.second; i++)
                    {
                        _FF_D_arr[i].setClk(_ClkList.size() - 1);
                        _FF_Q_arr[i].setClk(_ClkList.size() - 1);
                    }
                }
            }
        }
        else // is not a clk net
        {
            size_t netID = _NetList.size();
            vector<size_t> net;
            net.clear();
            net.reserve(netPinList.size());
            for (const auto &PinName : netPinList)
            {
                size_t pinID = _Name_to_ID[PinName];
                _ID_to_instance[pinID]->addRelatedNet(netID);
                net.push_back(pinID);
            }
            _NetList.push_back(net);
        }
    }

    // step3: initialize the STA engine
    std::cout << "step3: initialize the STA engine" << std::endl;
    assert(_ptr_STAEngine != nullptr);
    _ptr_STAEngine->setSolverPtr(this);
    _ptr_STAEngine->initEngine(_NetList, _ID_to_instance);

    // step4: update the slack for the pin, and distribute slack
    std::cout << "step4: update the slack for the pin, and distribute slack" << std::endl;
    for (const auto &N2Slack : _ptr_Parser->_timeSlack)
    {
        _FF_D_arr[_Name_to_ID.at(N2Slack.first) - _FF_D_OFFSET].setOriSlack(N2Slack.second);
    }

    // step5: init the placement row
    std::cout << "step5: init the placement row" << std::endl;
    for (const auto &plRow : _ptr_Parser->_placeRow)
    {
        struct PlacementRow placeRow;
        placeRow.x = plRow.x;
        placeRow.y = plRow.y;
        placeRow.siteHight = plRow.siteHight;
        placeRow.siteWidth = plRow.siteWidth;
        placeRow.totalNumOfSites = plRow.totalNumOfSites;
        _PlaceRow.push_back(placeRow);
    }

    // step6: categorize instance to the bin
    std::cout << "step6: categorize instance to the bin" << std::endl;
    //  compute the number of bin in x and y direction
    int numOfBin_x = (_ptr_Parser->_dieRx - (_ptr_Parser->_dieRx % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
    int numOfBin_y = (_ptr_Parser->_dieRy - (_ptr_Parser->_dieRy % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
    if (_ptr_Parser->_dieRx % _ptr_Parser->_binWidth)
        numOfBin_x++;
    if (_ptr_Parser->_dieRy % _ptr_Parser->_binHeight)
        numOfBin_y++;
    // initialize _Gate_in_Bin, _FF_in_Bin, _PlaceRow_in_Bin
    vector<size_t> blankList;
    blankList.reserve(512);
    vector<vector<size_t>> binList_y(numOfBin_y, blankList);
    vector<vector<vector<size_t>>> initBinList(numOfBin_x, binList_y);
    _Gate_in_Bin = initBinList;
    _FF_in_Bin = initBinList;
    _PlaceRow_in_Bin = initBinList;
    // sort gate into bin
    for (size_t i = 0; i < _Gate_arr.size(); i++)
    {
        pair<double, double> lowerLeft = _Gate_arr[i].getPosition();
        int x_index = (int(lowerLeft.first) - (int(lowerLeft.first) % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int y_index = (int(lowerLeft.second) - (int(lowerLeft.second) % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
        _Gate_in_Bin[x_index][y_index].push_back(i);
        bool overX = false;
        bool overY = false;
        if ((int(lowerLeft.first) % _ptr_Parser->_binWidth) + _ptr_Parser->_gateLib[_Gate_arr[i].gate_type].Width >= _ptr_Parser->_binWidth)
        {
            overX = true;
            _Gate_in_Bin[x_index + 1][y_index].push_back(i);
        }
        if ((int(lowerLeft.second) % _ptr_Parser->_binHeight) + _ptr_Parser->_gateLib[_Gate_arr[i].gate_type].Hight >= _ptr_Parser->_binHeight)
        {
            overY = true;
            _Gate_in_Bin[x_index][y_index + 1].push_back(i);
        }
        if (overX && overY)
        {
            _Gate_in_Bin[x_index + 1][y_index + 1].push_back(i);
        }
    }
    // sort FF into bin
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        size_t id = i + _FF_D_OFFSET;
        pair<double, double> pos = _FF_D_arr[i].getPosition();
        int x_index = (int(pos.first) - (int(pos.first) % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int y_index = (int(pos.second) - (int(pos.second) % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
        _FF_in_Bin[x_index][y_index].push_back(id);
    }
    // sort Placement Row into bin
    for (size_t i = 0; i < _PlaceRow.size(); i++)
    {
        int Lx_index = (int(_PlaceRow[i].x) - (int(_PlaceRow[i].x) % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int Ly_index = (int(_PlaceRow[i].y) - (int(_PlaceRow[i].y) % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
        int Rx = int(_PlaceRow[i].x + _PlaceRow[i].siteWidth * _PlaceRow[i].totalNumOfSites);
        int Ry = int(_PlaceRow[i].y + _PlaceRow[i].siteHight);
        int Rx_index = (Rx - (Rx % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int Ry_index = (Ry - (Ry % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
        for (int x = Lx_index; x < Rx_index; x++)
        {
            for (int y = Ly_index; y < Ry_index; y++)
            {
                _PlaceRow_in_Bin[x][y].push_back(i);
            }
        }
    }

    // step7: establish the fanin and fanout cone for the FF
    std::cout << "step7: establish the fanin and fanout cone for the FF" << std::endl;
    findFaninout4all(_ptr_STAEngine->_distanceList);

    // step8: initialize the legalizer
    std::cout << "step8: initialize the legalizer" << std::endl;
    assert(_ptr_legalizer != nullptr);
    _ptr_legalizer->setSolverPtr(this);

    // step9: calculate the Max slack for every FF
    std::cout << "step9: calculate the Max slack for every FF" << std::endl;
    findMaxSlack();

    // step10: evaluate the "grade" for ff library
    std::cout << "step10: evaluate the \"grade\" for ff library" << std::endl;
    _libGrade.clear();
    _libGrade.reserve(_ptr_Parser->_flipflopLib.size());
    double QpinOffset = 0;
    double oneBitNum = 0;
    int maxBit = 0;
    double maxBitNum = 0;
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib.at(i).Bit == 1)
        {
            oneBitNum += 1;
            QpinOffset += _ptr_Parser->_flipflopLib.at(i).PinDelay;
        }
        if (_ptr_Parser->_flipflopLib.at(i).Bit > maxBit)
        {
            maxBit = _ptr_Parser->_flipflopLib.at(i).Bit;
        }
    }
    QpinOffset = QpinOffset / oneBitNum;
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        struct libGrade init;
        init.aspectRatio = _ptr_Parser->_flipflopLib.at(i).Hight / _ptr_Parser->_flipflopLib.at(i).Width;
        init.QpinEffect = _ptr_Parser->_flipflopLib.at(i).PinDelay - QpinOffset;
        double area = _ptr_Parser->_flipflopLib.at(i).Hight * _ptr_Parser->_flipflopLib.at(i).Width;
        init.avgCost = (_ptr_Parser->_beta * _ptr_Parser->_flipflopLib.at(i).Power) + (_ptr_Parser->_gamma * area);
        init.avgCost = (init.avgCost / double(_ptr_Parser->_flipflopLib.at(i).Bit)) + (_ptr_Parser->_alpha * init.QpinEffect);
        _libGrade.push_back(init);
    }
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib.at(i).Bit == maxBit)
        {
            maxBitNum += 1;
        }
    }
    _avgQpinDelayDiff = 0;
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib.at(i).Bit == maxBit)
        {
            _avgQpinDelayDiff += _libGrade.at(i).QpinEffect / maxBitNum;
        }
    }

    std::cout << "End initialization!!! \n"
              << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THE NEW PART

void Solver::Solver::solve_initbuild()
{

    // build the coordinate vector
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        double Dx = findPinPosition(i + _FF_D_OFFSET).first;
        _FF_D_arr.at(i).Dx_pos = Dx;

        double Dy = findPinPosition(i + _FF_D_OFFSET).second;
        _FF_D_arr.at(i).Dy_pos = Dy;

        double Qx = findPinPosition(i + _FF_Q_OFFSET).first;
        _FF_Q_arr.at(i).Qx_pos = Qx;

        double Qy = findPinPosition(i + _FF_Q_OFFSET).second;
        _FF_Q_arr.at(i).Qy_pos = Qy;
    }

    std::cout << "initbuild is completed" << std::endl;
}

struct coor_w_se
{
    bool type; // s = 0; e = 1;
    double pos_val;
    int or_ind; // the index to represent each diamond, Dfanin diamond's index should be 0, 1~n is Qfanout diamond.
};
bool compareByPosVal(const coor_w_se &a, const coor_w_se &b)
{
    return a.pos_val < b.pos_val;
}

bool compareByPosVal_2(const Inst::feasible_coor &a, const Inst::feasible_coor &b)
{
    return a.pos_val < b.pos_val;
}

void Solver::Solver::solve_findfeasible()
{
    clock_t a, b;
    a = clock();
    // draw the feasible region, (use fanin fanout pos and the DQ slack ??)
    // rotate the coor first, or draw the feasible region first ?

    // 畫某個FF的feasible
    // feasible 邊長是 wire length + slack
    // 跑過所有fanin_arr(fanout)的開頭，遇到下一個fanin_arr(out)的結尾，就是x的解
    // 如果在所有開頭跑完之前，就遇到某個人的結尾，那無解，放原地

    // rotate the coordinate and find the feasible

    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        double dist = (_FF_Q_arr.at(i).Qx_pos - _FF_D_arr.at(i).Dx_pos) * 0.5;
        vector<coor_w_se> dia_x_arr;
        vector<coor_w_se> dia_y_arr;

        vector<pair<size_t, double>> fi_pinID_slack_Dvec = getSlack2ConnectedFF(i + _FF_D_OFFSET, false, 0);
        if (fi_pinID_slack_Dvec.size() != 0) // avoid no fanin
        {
            _FF_D_arr.at(i).D_fanin_pos = findPinPosition(fi_pinID_slack_Dvec.at(0).first);

            double fi_x_r = _FF_D_arr.at(i).D_fanin_pos.second + _FF_D_arr.at(i).D_fanin_pos.first + 2 * dist; // x'=y+x //D shift right
            double fi_y_r = _FF_D_arr.at(i).D_fanin_pos.second - _FF_D_arr.at(i).D_fanin_pos.first - 2 * dist; // y'=y-x

            double s = fi_pinID_slack_Dvec.at(0).second;

            double D_dia_len = fabs(_FF_D_arr.at(i).Dx_pos - _FF_D_arr.at(i).D_fanin_pos.first) + fabs(_FF_D_arr.at(i).Dy_pos - _FF_D_arr.at(i).D_fanin_pos.second) + s / (_ptr_Parser->_displaceDelay);

            if (D_dia_len < 0)
            {
                _FF_D_arr[i].hasfeasible = 0;
                _FF_D_arr.at(i).grouped = 0;
                continue;
            }

            coor_w_se Dx_s, Dx_e;
            dia_x_arr.reserve(256);
            Dx_s.pos_val = fi_x_r - (D_dia_len);
            Dx_s.type = 0;
            Dx_s.or_ind = 0;
            Dx_e.pos_val = fi_x_r + (D_dia_len);
            Dx_e.type = 1;
            Dx_e.or_ind = 0;

            dia_x_arr.push_back(Dx_s);
            dia_x_arr.push_back(Dx_e);

            coor_w_se Dy_s, Dy_e;
            dia_y_arr.reserve(256);
            Dy_s.pos_val = fi_y_r - (D_dia_len);
            Dy_s.type = 0;
            Dy_s.or_ind = 0;
            Dy_e.pos_val = fi_y_r + (D_dia_len);
            Dy_e.type = 1;
            Dy_e.or_ind = 0;

            dia_y_arr.push_back(Dy_s);
            dia_y_arr.push_back(Dy_e);
        }

        // Q part, run through all fanout.
        vector<pair<size_t, double>> fo_pinID_slack_Qvec = getSlack2ConnectedFF(i + _FF_Q_OFFSET);

        if (fo_pinID_slack_Qvec.size() != 0)
        {
            bool Q_break = 0;

            for (size_t j = 0; j < fo_pinID_slack_Qvec.size(); j++)
            {
                _FF_Q_arr.at(i).Q_fanout_pos.push_back(findPinPosition(fo_pinID_slack_Qvec.at(j).first));

                double fo_x_r = _FF_Q_arr.at(i).Q_fanout_pos.at(j).second + _FF_Q_arr.at(i).Q_fanout_pos.at(j).first - 2 * dist; // x'=y+x //Q shift left
                double fo_y_r = _FF_Q_arr.at(i).Q_fanout_pos.at(j).second - _FF_Q_arr.at(i).Q_fanout_pos.at(j).first + 2 * dist; // y'=y-x

                double s = fo_pinID_slack_Qvec.at(j).second;

                double Q_dia_len = fabs(_FF_Q_arr.at(i).Qx_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(i).Qy_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).second) + s / (_ptr_Parser->_displaceDelay);
                if (Q_dia_len < 0)
                {
                    _FF_D_arr[i].hasfeasible = 0;
                    _FF_D_arr.at(i).grouped = 0;
                    Q_break = 1;
                    break;
                }

                coor_w_se Qx_s, Qx_e;
                Qx_s.pos_val = fo_x_r - (Q_dia_len);
                Qx_s.type = 0;
                Qx_s.or_ind = j + 1;
                Qx_e.pos_val = fo_x_r + (Q_dia_len);
                Qx_e.type = 1;
                Qx_e.or_ind = j + 1;

                dia_x_arr.push_back(Qx_s);
                dia_x_arr.push_back(Qx_e);

                coor_w_se Qy_s, Qy_e;
                Qy_s.pos_val = fo_y_r - (Q_dia_len);
                Qy_s.type = 0;
                Qy_s.or_ind = j + 1;
                Qy_e.pos_val = fo_y_r + (Q_dia_len);
                Qy_e.type = 1;
                Qy_e.or_ind = j + 1;

                dia_y_arr.push_back(Qy_s);
                dia_y_arr.push_back(Qy_e);
            }
            if (Q_break == 1)
            {
                continue;
            }
        }

        // draw feasible by using dia_x_arr, dia_y_arr
        // 對 dia_x_arr 排序
        std::sort(dia_x_arr.begin(), dia_x_arr.end(), compareByPosVal);

        // 對 dia_y_arr 排序
        std::sort(dia_y_arr.begin(), dia_y_arr.end(), compareByPosVal);

        int target_x;
        int bound_ctr_x = 0;
        for (size_t k = 1; k < dia_x_arr.size(); k++)
        {
            if ((dia_x_arr.at(k).type == 1) && (dia_x_arr.at(k - 1).type == 0))
            {
                target_x = k;
                bound_ctr_x++;
            }
        }

        int target_y;
        int bound_ctr_y = 0;
        for (size_t k = 1; k < dia_y_arr.size(); k++)
        {
            if ((dia_y_arr.at(k).type == 1) && (dia_y_arr.at(k - 1).type == 0))
            {
                target_y = k;
                bound_ctr_y++;
            }
        }

        if ((bound_ctr_x != 1) || (bound_ctr_y != 1))
        {
            // no feasible
            _FF_D_arr.at(i).grouped = 0;
            _FF_D_arr.at(i).hasfeasible = 0;
        }
        else
        {
            // feasible
            _FF_D_arr.at(i).grouped = 0;
            _FF_D_arr.at(i).hasfeasible = 1;

            _FF_D_arr.at(i).fea_x_s.type = 0;
            _FF_D_arr.at(i).fea_x_s.pos_val = dia_x_arr.at(target_x - 1).pos_val;
            _FF_D_arr.at(i).fea_x_s.FF_id = i;

            _FF_D_arr.at(i).fea_x_e.type = 1;
            _FF_D_arr.at(i).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
            _FF_D_arr.at(i).fea_x_e.FF_id = i;

            _FF_D_arr.at(i).fea_y_s.type = 0;
            _FF_D_arr.at(i).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
            _FF_D_arr.at(i).fea_y_s.FF_id = i;

            _FF_D_arr.at(i).fea_y_e.type = 1;
            _FF_D_arr.at(i).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
            _FF_D_arr.at(i).fea_y_e.FF_id = i;
        }
        if (_FF_D_arr.at(i).hasfeasible == 1)
        {
            // std::cout << _FF_D_arr.at(i).fea_x_s.pos_val << " " << _FF_D_arr.at(i).fea_x_e.pos_val << std::endl;
        }
    }

    std::cout << "findfeasible is completed" << std::endl;
    b = clock();
    std::cout << "time for finding feasible: " << (b - a) / CLOCKS_PER_SEC << std::endl;
}

void Solver::Solver::solve()
{
    evaluate("init Metrix.txt");
    solve_initbuild();

    solve_findfeasible();

    // drawpic("pic_1.plt"); /////////////////////////////////////////////

    // clock一樣才能綁，外面掛一層迴圈跑過所有clock
    // FFD.getclock 吐 id, id == 0 ... for(int i = 0; i < _ClkList.size(); i++)

    for (size_t k = 0; k < _ClkList.size(); k++) // testing
    {                                            /*run through clk*/

        vector<Inst::feasible_coor> feas_x_clk;
        feas_x_clk.reserve(_FF_D_arr.size() * 2);
        // vector<Inst::feasible_coor> feas_y_clk;
        // feas_y_clk.reserve(_FF_D_arr.size() * 2);

        for (int i = 0; i < _FF_D_arr.size(); i++)
        {

            if ((_FF_D_arr[i].getClk() == k) && (_FF_D_arr[i].hasfeasible == 1))
            {
                // collect
                // may collect the empty value of no-feasible FF
                feas_x_clk.push_back(_FF_D_arr.at(i).fea_x_s);
                feas_x_clk.push_back(_FF_D_arr.at(i).fea_x_e);

                // feas_y_clk.push_back(_FF_D_arr.at(i).fea_y_s);
                // feas_y_clk.push_back(_FF_D_arr.at(i).fea_y_e);
            }
        }

        std::sort(feas_x_clk.begin(), feas_x_clk.end(), compareByPosVal_2);
        // std::sort(feas_y_clk.begin(), feas_y_clk.end(), compareByPosVal_2);

        size_t esssential_ff;
        // find maximal clique
        while (!feas_x_clk.empty())
        {
            vector<size_t> ff_group;

            // std::cout << "Before clustering, the size: " << feas_x_clk.size() << std::endl;
            for (int i = 1; i < feas_x_clk.size(); i++)
            {

                ff_group.push_back(feas_x_clk.at(i - 1).FF_id + _FF_D_OFFSET);

                if ((feas_x_clk.at(i - 1).type) == 0 && (feas_x_clk.at(i).type == 1)) // 0 is start, 1 is end
                {
                    // ff_group 有一個 essential, result_group有，essential已經在裡面了，不用pushback
                    // record the essential is i

                    esssential_ff = feas_x_clk.at(i).FF_id + _FF_D_OFFSET;

                    // consider the y part
                    pair<double, double> x_pos_r, y_pos_r; // final region
                    vector<size_t> final_group;
                    vector<size_t> result_group;
                    // find the position, 還原成正常座標

                    // x_pos_r and y_pos_r is pass by ref and get the pos value by this function
                    result_group = solve_findmaximal(ff_group, esssential_ff, x_pos_r, y_pos_r);

                    /// another way: change the essential ff
                    /// new part to get a more flexible feasible region
                    vector<size_t> new_result_group;
                    new_result_group.reserve(4);
                    if (result_group.size() > 4)
                    {
                        esssential_ff = result_group[0];

                        for (int id = 0; id < 4; id++)
                        {
                            new_result_group.push_back(result_group[id]);
                        }
                        result_group = solve_findmaximal(new_result_group, esssential_ff, x_pos_r, y_pos_r);
                        // std::cout << "group size: " << result_group.size() << std::endl;
                        assert(result_group.size() == 4);
                    }
                    ///////////////////////////////////////////////////

                    // preplace and slack release
                    pair<double, double> pos;
                    pos.first = (x_pos_r.first + x_pos_r.second) * 0.5; // choose the middle point
                    pos.second = (y_pos_r.first + y_pos_r.second) * 0.5;
                    double x_ori = pos.first;
                    double y_ori = pos.second;
                    // change to original coordinate (x'=y+x, y'=y-x)
                    pos.first = (x_ori - y_ori) * 0.5;  // x=(x'-y')/2
                    pos.second = (x_ori + y_ori) * 0.5; // y=(x'+y')/2

                    //   std::cout << "preplace begin" << std::endl;
                    final_group = prePlace(result_group, esssential_ff, pos);
                    // std::cout << "finalgroup size: " << final_group.size() <<std::endl;

                    // final_group有essential

                    //  calculate the feasible

                    feasible_cal(final_group, k);

                    // solve_findfeasible();
                    break;
                }
            }

            // delete the grouped member in feas_x_clk
            for (size_t i = 0; i < feas_x_clk.size(); i++)
            {
                // std::cout << "FF" << feas_x_clk.at(i).FF_id << " group or not: " << _FF_D_arr[feas_x_clk.at(i).FF_id].grouped << std::endl;
                if (_FF_D_arr[feas_x_clk.at(i).FF_id].grouped == 1)
                {
                    feas_x_clk.erase(feas_x_clk.begin() + i);
                    i = i - 1;
                }
            }

            // std::cout << "Clear is completed" << std::endl;
            // std::cout << "After clustering, the size: " << feas_x_clk.size() << std::endl;
            // std::cout << _FF_D_arr[feas_x_clk[0].FF_id].fea_x_s.pos_val << std::endl;
            // std::cout << _FF_D_arr[feas_x_clk[0].FF_id].fea_x_e.pos_val << std::endl;

            // std::cout << _FF_D_arr[feas_x_clk[0].FF_id].getOriName() << std::endl;
            // std::cout << _FF_D_arr[feas_x_clk[1].FF_id].getOriName() << std::endl;

            // std::cout << "FFID is " << feas_x_clk[0].FF_id << std::endl;
        }
        std::cout << "CLKlist " << k << " is finished" << std::endl;
        ////////////////////////////////////////////////////////////////////////////////////////
        if (k % (_ClkList.size() / 5) == 0)
        {
            string s = std::to_string(k / (_ClkList.size() / 5) + 1);
            // drawpic("pic_" + s + ".plt");
        }
        ////////////////////////////////////////////////////////////////////////////////////////
    }

    // 決定後來的DQ 的 pos, 記得Q的正方形要往左移（或D往右, Q往左）
    // preplace完 slack release, 更新DQ正方形，重新畫table
    drawpic("pic_z.plt"); ///////////////////////////////////////////////////////////////////////
    evaluate("before legal Metrix.txt");
    legalize();
    evaluate("after legal Metrix.txt");
    // drawpic("pic_final.plt"); ///////////////////////////////////////////////////////////////////
    std::cout << "Solver is completed !" << std::endl;
}

void Solver::Solver::windowbinfeasible(vector<size_t> &ffbingroup)
{

    if (ffbingroup.size() == 0)
    {
        return;
    }

    for (size_t bin = 0; bin < ffbingroup.size(); bin++)
    {
        size_t i = ffbingroup[bin];
        double dist = (_FF_Q_arr.at(i).Qx_pos - _FF_D_arr.at(i).Dx_pos) * 0.5;
        vector<coor_w_se> dia_x_arr;
        vector<coor_w_se> dia_y_arr;

        vector<pair<size_t, double>> fi_pinID_slack_Dvec = getSlack2ConnectedFF(i + _FF_D_OFFSET, true, 1);
        if (fi_pinID_slack_Dvec.size() != 0) // avoid no fanin
        {
            _FF_D_arr.at(i).D_fanin_pos = findPinPosition(fi_pinID_slack_Dvec.at(0).first);

            double fi_x_r = _FF_D_arr.at(i).D_fanin_pos.second + _FF_D_arr.at(i).D_fanin_pos.first + 2 * dist; // x'=y+x //D shift right
            double fi_y_r = _FF_D_arr.at(i).D_fanin_pos.second - _FF_D_arr.at(i).D_fanin_pos.first - 2 * dist; // y'=y-x

            double s = fi_pinID_slack_Dvec.at(0).second;

            double D_dia_len = fabs(_FF_D_arr.at(i).Dx_pos - _FF_D_arr.at(i).D_fanin_pos.first) + fabs(_FF_D_arr.at(i).Dy_pos - _FF_D_arr.at(i).D_fanin_pos.second) + s / (_ptr_Parser->_displaceDelay);

            if (D_dia_len < 0)
            {
                _FF_D_arr[i].hasfeasible = 0;
                continue;
            }

            coor_w_se Dx_s, Dx_e;
            dia_x_arr.reserve(256);
            Dx_s.pos_val = fi_x_r - (D_dia_len);
            Dx_s.type = 0;
            Dx_s.or_ind = 0;
            Dx_e.pos_val = fi_x_r + (D_dia_len);
            Dx_e.type = 1;
            Dx_e.or_ind = 0;

            dia_x_arr.push_back(Dx_s);
            dia_x_arr.push_back(Dx_e);

            coor_w_se Dy_s, Dy_e;
            dia_y_arr.reserve(256);
            Dy_s.pos_val = fi_y_r - (D_dia_len);
            Dy_s.type = 0;
            Dy_s.or_ind = 0;
            Dy_e.pos_val = fi_y_r + (D_dia_len);
            Dy_e.type = 1;
            Dy_e.or_ind = 0;

            dia_y_arr.push_back(Dy_s);
            dia_y_arr.push_back(Dy_e);
        }

        // Q part, run through all fanout.
        vector<pair<size_t, double>> fo_pinID_slack_Qvec = getSlack2ConnectedFF(i + _FF_Q_OFFSET, true, 1);

        if (fo_pinID_slack_Qvec.size() != 0)
        {
            bool Q_break = 0;

            for (size_t j = 0; j < fo_pinID_slack_Qvec.size(); j++)
            {
                _FF_Q_arr.at(i).Q_fanout_pos.push_back(findPinPosition(fo_pinID_slack_Qvec.at(j).first));

                double fo_x_r = _FF_Q_arr.at(i).Q_fanout_pos.at(j).second + _FF_Q_arr.at(i).Q_fanout_pos.at(j).first - 2 * dist; // x'=y+x //Q shift left
                double fo_y_r = _FF_Q_arr.at(i).Q_fanout_pos.at(j).second - _FF_Q_arr.at(i).Q_fanout_pos.at(j).first + 2 * dist; // y'=y-x

                double s = fo_pinID_slack_Qvec.at(j).second;

                double Q_dia_len = fabs(_FF_Q_arr.at(i).Qx_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(i).Qy_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).second) + s / (_ptr_Parser->_displaceDelay);
                if (Q_dia_len < 0)
                {
                    _FF_D_arr[i].hasfeasible = 0;
                    Q_break = 1;
                    break;
                }

                coor_w_se Qx_s, Qx_e;
                Qx_s.pos_val = fo_x_r - (Q_dia_len);
                Qx_s.type = 0;
                Qx_s.or_ind = j + 1;
                Qx_e.pos_val = fo_x_r + (Q_dia_len);
                Qx_e.type = 1;
                Qx_e.or_ind = j + 1;

                dia_x_arr.push_back(Qx_s);
                dia_x_arr.push_back(Qx_e);

                coor_w_se Qy_s, Qy_e;
                Qy_s.pos_val = fo_y_r - (Q_dia_len);
                Qy_s.type = 0;
                Qy_s.or_ind = j + 1;
                Qy_e.pos_val = fo_y_r + (Q_dia_len);
                Qy_e.type = 1;
                Qy_e.or_ind = j + 1;

                dia_y_arr.push_back(Qy_s);
                dia_y_arr.push_back(Qy_e);
            }
            if (Q_break == 1)
            {
                continue;
            }
        }

        // draw feasible by using dia_x_arr, dia_y_arr
        // 對 dia_x_arr 排序
        std::sort(dia_x_arr.begin(), dia_x_arr.end(), compareByPosVal);

        // 對 dia_y_arr 排序
        std::sort(dia_y_arr.begin(), dia_y_arr.end(), compareByPosVal);

        int target_x;
        int bound_ctr_x = 0;
        for (size_t k = 1; k < dia_x_arr.size(); k++)
        {
            if ((dia_x_arr.at(k).type == 1) && (dia_x_arr.at(k - 1).type == 0))
            {
                target_x = k;
                bound_ctr_x++;
            }
        }

        int target_y;
        int bound_ctr_y = 0;
        for (size_t k = 1; k < dia_y_arr.size(); k++)
        {
            if ((dia_y_arr.at(k).type == 1) && (dia_y_arr.at(k - 1).type == 0))
            {
                target_y = k;
                bound_ctr_y++;
            }
        }

        if ((bound_ctr_x != 1) || (bound_ctr_y != 1))
        {
            // no feasible

            _FF_D_arr.at(i).hasfeasible = 0;
        }
        else
        {
            // feasible

            _FF_D_arr.at(i).hasfeasible = 1;

            _FF_D_arr.at(i).fea_x_s.type = 0;
            _FF_D_arr.at(i).fea_x_s.pos_val = dia_x_arr.at(target_x - 1).pos_val;
            _FF_D_arr.at(i).fea_x_s.FF_id = i;

            _FF_D_arr.at(i).fea_x_e.type = 1;
            _FF_D_arr.at(i).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
            _FF_D_arr.at(i).fea_x_e.FF_id = i;

            _FF_D_arr.at(i).fea_y_s.type = 0;
            _FF_D_arr.at(i).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
            _FF_D_arr.at(i).fea_y_s.FF_id = i;

            _FF_D_arr.at(i).fea_y_e.type = 1;
            _FF_D_arr.at(i).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
            _FF_D_arr.at(i).fea_y_e.FF_id = i;
        }
    }

    // std::cout << "findfeasible is completed" << std::endl;
}

void Solver::Solver::feasible_cal(const vector<size_t> &final_group, const size_t &k)
{

    ////faninCone
    for (size_t i = 0; i < final_group.size(); i++)
    {
        size_t id = final_group.at(i) - _FF_D_OFFSET;
        Inst::FF_D *ptr_FF_D = &_FF_D_arr[id];

        if (ptr_FF_D->faninCone.size() == 0)
        {
            continue;
        }

        for (size_t kk = 0; kk < ptr_FF_D->faninCone.size(); kk++)
        {
            size_t FFid = ptr_FF_D->faninCone[kk] - _FF_Q_OFFSET;
            // bool renew = FFinwindow(&_FF_D_arr[FFid], w, w_r, h, h_u);

            if (_FF_D_arr[FFid].grouped == false && (_FF_D_arr[FFid].getClk() == k || _FF_D_arr[FFid].getClk() == k + 1)) // if the slack release to the grouped FF, what should i do ?
            {

                double dist = (_FF_Q_arr.at(FFid).Qx_pos - _FF_D_arr.at(FFid).Dx_pos) * 0.5;
                vector<coor_w_se> dia_x_arr;
                vector<coor_w_se> dia_y_arr;

                vector<pair<size_t, double>> fi_pinID_slack_Dvec = getSlack2ConnectedFF(FFid + _FF_D_OFFSET, false, 0);
                if (fi_pinID_slack_Dvec.size() != 0) // avoid no fanin
                {
                    _FF_D_arr.at(FFid).D_fanin_pos = findPinPosition(fi_pinID_slack_Dvec.at(0).first);

                    double fi_x_r = _FF_D_arr.at(FFid).D_fanin_pos.second + _FF_D_arr.at(FFid).D_fanin_pos.first + 2 * dist; // x'=y+x //D shift right
                    double fi_y_r = _FF_D_arr.at(FFid).D_fanin_pos.second - _FF_D_arr.at(FFid).D_fanin_pos.first - 2 * dist; // y'=y-x

                    double s = fi_pinID_slack_Dvec.at(0).second;

                    double D_dia_len = fabs(_FF_D_arr.at(FFid).Dx_pos - _FF_D_arr.at(FFid).D_fanin_pos.first) + fabs(_FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).D_fanin_pos.second) + s / (_ptr_Parser->_displaceDelay);
                    if (D_dia_len < 0)
                    {
                        _FF_D_arr[FFid].grouped = 0;
                        _FF_D_arr[FFid].hasfeasible = 0;
                        break;
                    }

                    coor_w_se Dx_s, Dx_e;
                    dia_x_arr.reserve(256);
                    Dx_s.pos_val = fi_x_r - (D_dia_len);
                    Dx_s.type = 0;
                    Dx_s.or_ind = 0;
                    Dx_e.pos_val = fi_x_r + (D_dia_len);
                    Dx_e.type = 1;
                    Dx_e.or_ind = 0;

                    dia_x_arr.push_back(Dx_s);
                    dia_x_arr.push_back(Dx_e);

                    coor_w_se Dy_s, Dy_e;
                    dia_y_arr.reserve(256);
                    Dy_s.pos_val = fi_y_r - (D_dia_len);
                    Dy_s.type = 0;
                    Dy_s.or_ind = 0;
                    Dy_e.pos_val = fi_y_r + (D_dia_len);
                    Dy_e.type = 1;
                    Dy_e.or_ind = 0;

                    dia_y_arr.push_back(Dy_s);
                    dia_y_arr.push_back(Dy_e);
                }

                // Q part, run through all fanout.
                vector<pair<size_t, double>> fo_pinID_slack_Qvec = getSlack2ConnectedFF(FFid + _FF_Q_OFFSET, false, 0);

                if (fo_pinID_slack_Qvec.size() != 0) // avoid PI and no fanout
                {
                    bool Q_break = 0;

                    for (size_t j = 0; j < fo_pinID_slack_Qvec.size(); j++)
                    {

                        _FF_Q_arr.at(FFid).Q_fanout_pos.push_back(findPinPosition(fo_pinID_slack_Qvec.at(j).first));

                        double fo_x_r = _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second + _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first - 2 * dist; // x'=y+x //Q shift left
                        double fo_y_r = _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first + 2 * dist; // y'=y-x

                        double s = fo_pinID_slack_Qvec.at(j).second;

                        double Q_dia_len = fabs(_FF_Q_arr.at(FFid).Qx_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second) + s / (_ptr_Parser->_displaceDelay);

                        if (Q_dia_len < 0)
                        {
                            _FF_D_arr[FFid].hasfeasible = 0;
                            _FF_D_arr.at(FFid).grouped = 0;
                            Q_break = 1;
                            break;
                        }
                        coor_w_se Qx_s, Qx_e;
                        Qx_s.pos_val = fo_x_r - (Q_dia_len);
                        Qx_s.type = 0;
                        Qx_s.or_ind = j + 1;
                        Qx_e.pos_val = fo_x_r + (Q_dia_len);
                        Qx_e.type = 1;
                        Qx_e.or_ind = j + 1;

                        dia_x_arr.push_back(Qx_s);
                        dia_x_arr.push_back(Qx_e);

                        coor_w_se Qy_s, Qy_e;
                        Qy_s.pos_val = fo_y_r - (Q_dia_len);
                        Qy_s.type = 0;
                        Qy_s.or_ind = j + 1;
                        Qy_e.pos_val = fo_y_r + (Q_dia_len);
                        Qy_e.type = 1;
                        Qy_e.or_ind = j + 1;

                        dia_y_arr.push_back(Qy_s);
                        dia_y_arr.push_back(Qy_e);
                    }
                    if (Q_break == 1)
                    {
                        continue;
                    }
                }

                // draw feasible by using dia_x_arr, dia_y_arr
                // 對 dia_x_arr 排序
                std::sort(dia_x_arr.begin(), dia_x_arr.end(), compareByPosVal);

                // 對 dia_y_arr 排序
                std::sort(dia_y_arr.begin(), dia_y_arr.end(), compareByPosVal);

                int target_x;
                int bound_ctr_x = 0;
                for (size_t k = 1; k < dia_x_arr.size(); k++)
                {
                    if ((dia_x_arr.at(k).type == 1) && (dia_x_arr.at(k - 1).type == 0))
                    {
                        target_x = k;
                        bound_ctr_x++;
                    }
                }

                int target_y;
                int bound_ctr_y = 0;
                for (size_t k = 1; k < dia_y_arr.size(); k++)
                {
                    if ((dia_y_arr.at(k).type == 1) && (dia_y_arr.at(k - 1).type == 0))
                    {
                        target_y = k;
                        bound_ctr_y++;
                    }
                }

                if ((bound_ctr_x != 1) || (bound_ctr_y != 1))
                {
                    // no feasible

                    _FF_D_arr.at(FFid).hasfeasible = 0;
                }
                else
                {
                    // feasible
                    _FF_D_arr.at(FFid).grouped = 0;
                    _FF_D_arr.at(FFid).hasfeasible = 1;

                    _FF_D_arr.at(FFid).fea_x_s.type = 0;
                    _FF_D_arr.at(FFid).fea_x_s.pos_val = dia_x_arr.at(target_x - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_x_s.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_x_e.type = 1;
                    _FF_D_arr.at(FFid).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
                    _FF_D_arr.at(FFid).fea_x_e.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_y_s.type = 0;
                    _FF_D_arr.at(FFid).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_y_s.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_y_e.type = 1;
                    _FF_D_arr.at(FFid).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
                    _FF_D_arr.at(FFid).fea_y_e.FF_id = FFid;
                }
            }
        }
    }

    ////fanoutCone
    for (size_t i = 0; i < final_group.size(); i++)
    {
        size_t id = final_group.at(i) - _FF_D_OFFSET;
        Inst::FF_Q *ptr_FF_Q = &_FF_Q_arr[id];

        if (ptr_FF_Q->fanoutCone.size() == 0)
        {
            continue;
        }

        for (size_t kk = 0; kk < ptr_FF_Q->fanoutCone.size(); kk++)
        {
            size_t FFid = ptr_FF_Q->fanoutCone[kk] - _FF_D_OFFSET;
            // bool renew = FFinwindow(&_FF_D_arr[FFid], w, w_r, h, h_u);

            if (_FF_D_arr[FFid].grouped == false && (_FF_D_arr[FFid].getClk() == k || _FF_D_arr[FFid].getClk() == k + 1)) // if the slack release to the grouped FF, what should i do ?
            {

                double dist = (_FF_Q_arr.at(FFid).Qx_pos - _FF_D_arr.at(FFid).Dx_pos) * 0.5;
                vector<coor_w_se> dia_x_arr;
                vector<coor_w_se> dia_y_arr;

                vector<pair<size_t, double>> fo_pinID_slack_Dvec = getSlack2ConnectedFF(FFid + _FF_D_OFFSET, false, 0);

                if (fo_pinID_slack_Dvec.size() != 0) // avoid no fanin
                {
                    _FF_D_arr.at(FFid).D_fanin_pos = findPinPosition(fo_pinID_slack_Dvec.at(0).first);

                    double fi_x_r = _FF_D_arr.at(FFid).D_fanin_pos.second + _FF_D_arr.at(FFid).D_fanin_pos.first + 2 * dist; // x'=y+x //D shift right
                    double fi_y_r = _FF_D_arr.at(FFid).D_fanin_pos.second - _FF_D_arr.at(FFid).D_fanin_pos.first - 2 * dist; // y'=y-x

                    double s = fo_pinID_slack_Dvec.at(0).second;

                    double D_dia_len = fabs(_FF_D_arr.at(FFid).Dx_pos - _FF_D_arr.at(FFid).D_fanin_pos.first) + fabs(_FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).D_fanin_pos.second) + s / (_ptr_Parser->_displaceDelay);
                    if (D_dia_len < 0)
                    {
                        _FF_D_arr[FFid].grouped = 0;
                        _FF_D_arr[FFid].hasfeasible = 0;
                        break;
                    }

                    coor_w_se Dx_s, Dx_e;
                    dia_x_arr.reserve(256);
                    Dx_s.pos_val = fi_x_r - (D_dia_len);
                    Dx_s.type = 0;
                    Dx_s.or_ind = 0;
                    Dx_e.pos_val = fi_x_r + (D_dia_len);
                    Dx_e.type = 1;
                    Dx_e.or_ind = 0;

                    dia_x_arr.push_back(Dx_s);
                    dia_x_arr.push_back(Dx_e);

                    coor_w_se Dy_s, Dy_e;
                    dia_y_arr.reserve(256);
                    Dy_s.pos_val = fi_y_r - (D_dia_len);
                    Dy_s.type = 0;
                    Dy_s.or_ind = 0;
                    Dy_e.pos_val = fi_y_r + (D_dia_len);
                    Dy_e.type = 1;
                    Dy_e.or_ind = 0;

                    dia_y_arr.push_back(Dy_s);
                    dia_y_arr.push_back(Dy_e);
                }

                // Q part, run through all fanout.

                vector<pair<size_t, double>> fo_pinID_slack_Qvec = getSlack2ConnectedFF(FFid + _FF_Q_OFFSET, false, 0);
                if (fo_pinID_slack_Qvec.size() != 0)
                {
                    bool Q_break = 0;

                    for (size_t j = 0; j < fo_pinID_slack_Qvec.size(); j++)
                    {

                        _FF_Q_arr.at(FFid).Q_fanout_pos.push_back(findPinPosition(fo_pinID_slack_Qvec.at(j).first));

                        double fo_x_r = _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second + _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first - 2 * dist; // x'=y+x //Q shift left
                        double fo_y_r = _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first + 2 * dist; // y'=y-x

                        double s = fo_pinID_slack_Qvec.at(j).second;

                        double Q_dia_len = fabs(_FF_Q_arr.at(FFid).Qx_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second) + s / (_ptr_Parser->_displaceDelay);
                        if (Q_dia_len < 0)
                        {
                            _FF_D_arr[FFid].hasfeasible = 0;
                            _FF_D_arr.at(FFid).grouped = 0;
                            Q_break = 1;
                            break;
                        }
                        coor_w_se Qx_s, Qx_e;
                        Qx_s.pos_val = fo_x_r - (Q_dia_len);
                        Qx_s.type = 0;
                        Qx_s.or_ind = j + 1;
                        Qx_e.pos_val = fo_x_r + (Q_dia_len);
                        Qx_e.type = 1;
                        Qx_e.or_ind = j + 1;

                        dia_x_arr.push_back(Qx_s);
                        dia_x_arr.push_back(Qx_e);

                        coor_w_se Qy_s, Qy_e;
                        Qy_s.pos_val = fo_y_r - (Q_dia_len);
                        Qy_s.type = 0;
                        Qy_s.or_ind = j + 1;
                        Qy_e.pos_val = fo_y_r + (Q_dia_len);
                        Qy_e.type = 1;
                        Qy_e.or_ind = j + 1;

                        dia_y_arr.push_back(Qy_s);
                        dia_y_arr.push_back(Qy_e);
                    }
                    if (Q_break == 1)
                    {
                        continue;
                    }
                }

                // draw feasible by using dia_x_arr, dia_y_arr
                // 對 dia_x_arr 排序
                std::sort(dia_x_arr.begin(), dia_x_arr.end(), compareByPosVal);

                // 對 dia_y_arr 排序
                std::sort(dia_y_arr.begin(), dia_y_arr.end(), compareByPosVal);

                int target_x;
                int bound_ctr_x = 0;
                for (size_t k = 1; k < dia_x_arr.size(); k++)
                {
                    if ((dia_x_arr.at(k).type == 1) && (dia_x_arr.at(k - 1).type == 0))
                    {
                        target_x = k;
                        bound_ctr_x++;
                    }
                }

                int target_y;
                int bound_ctr_y = 0;
                for (size_t k = 1; k < dia_y_arr.size(); k++)
                {
                    if ((dia_y_arr.at(k).type == 1) && (dia_y_arr.at(k - 1).type == 0))
                    {
                        target_y = k;
                        bound_ctr_y++;
                    }
                }

                if ((bound_ctr_x != 1) || (bound_ctr_y != 1))
                {
                    // no feasible

                    _FF_D_arr.at(FFid).hasfeasible = 0;
                }
                else
                {
                    // feasible
                    _FF_D_arr.at(FFid).grouped = 0;
                    _FF_D_arr.at(FFid).hasfeasible = 1;

                    _FF_D_arr.at(FFid).fea_x_s.type = 0;
                    _FF_D_arr.at(FFid).fea_x_s.pos_val = dia_x_arr.at(target_x - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_x_s.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_x_e.type = 1;
                    _FF_D_arr.at(FFid).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
                    _FF_D_arr.at(FFid).fea_x_e.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_y_s.type = 0;
                    _FF_D_arr.at(FFid).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_y_s.FF_id = FFid;

                    _FF_D_arr.at(FFid).fea_y_e.type = 1;
                    _FF_D_arr.at(FFid).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
                    _FF_D_arr.at(FFid).fea_y_e.FF_id = FFid;
                }
            }
        }
    }
}

vector<size_t> Solver::Solver::solve_findmaximal(const vector<size_t> &ff_group, size_t essential_ID, pair<double, double> &x_pos_r, pair<double, double> &y_pos_r)
{
    vector<size_t> result_group; // bigger id
    double ess_x_s = _FF_D_arr.at(essential_ID - _FF_D_OFFSET).fea_x_s.pos_val;
    double ess_x_e = _FF_D_arr.at(essential_ID - _FF_D_OFFSET).fea_x_e.pos_val;
    double ess_y_s = _FF_D_arr.at(essential_ID - _FF_D_OFFSET).fea_y_s.pos_val;
    double ess_y_e = _FF_D_arr.at(essential_ID - _FF_D_OFFSET).fea_y_e.pos_val;
    pair<double, double> x_pos_r_temp, y_pos_r_temp;
    x_pos_r_temp.first = ess_x_s;
    x_pos_r_temp.second = ess_x_e;
    y_pos_r_temp.first = ess_y_s;
    y_pos_r_temp.second = ess_y_e;

    for (size_t i = 0; i < ff_group.size(); i++)
    {
        double ff_x_s = _FF_D_arr.at(ff_group.at(i) - _FF_D_OFFSET).fea_x_s.pos_val;
        double ff_x_e = _FF_D_arr.at(ff_group.at(i) - _FF_D_OFFSET).fea_x_e.pos_val;
        double ff_y_s = _FF_D_arr.at(ff_group.at(i) - _FF_D_OFFSET).fea_y_s.pos_val;
        double ff_y_e = _FF_D_arr.at(ff_group.at(i) - _FF_D_OFFSET).fea_y_e.pos_val;
        // std::cout << ff_x_s << " ," << ff_x_e << std::endl;
        size_t id = _FF_D_arr.at(ff_group.at(i) - _FF_D_OFFSET).fea_x_s.FF_id;

        bool valid_x = 0;
        bool valid_y = 0;
        pair<double, double> x_ori, y_ori;
        x_ori = x_pos_r_temp; // record the original best area (pos_temp)
        y_ori = y_pos_r_temp;

        if (x_pos_r_temp.first < ff_x_s)
        {
            if (x_pos_r_temp.second >= ff_x_s && x_pos_r_temp.second < ff_x_e)
            {
                x_pos_r_temp.first = ff_x_s;
                valid_x = 1;
                // std::cout << x_pos_r.first << " ," << x_pos_r.second << std::endl;
            }
            else if (x_pos_r_temp.second >= ff_x_e)
            {
                x_pos_r_temp.first = ff_x_s;
                x_pos_r_temp.second = ff_x_e;
                valid_x = 1;
            }
            else
            {
                valid_x = 0;
            }
        }
        else
        { // temp_x_s >= ff_x_s
            if (ff_x_e >= x_pos_r_temp.first && ff_x_e < x_pos_r_temp.second)
            {
                x_pos_r_temp.second = ff_x_e;
                valid_x = 1;
            }
            else if (x_pos_r_temp.second <= ff_x_e)
            {
                valid_x = 1;
            }
            else
            {
                valid_x = 0;
            }
        }

        if (y_pos_r_temp.first < ff_y_s)
        {
            if (y_pos_r_temp.second >= ff_y_s && y_pos_r_temp.second < ff_y_e)
            {
                y_pos_r_temp.first = ff_y_s;
                valid_y = 1;
            }
            else if (y_pos_r_temp.second >= ff_y_e)
            {
                y_pos_r_temp.first = ff_y_s;
                y_pos_r_temp.second = ff_y_e;
                valid_y = 1;
            }
            else
            {
                valid_y = 0;
            }
        }
        else
        { // temp_y_s >= ff_y_s
            if (ff_y_e >= y_pos_r_temp.first && ff_y_e < y_pos_r_temp.second)
            {
                y_pos_r_temp.second = ff_y_e;
                valid_y = 1;
            }
            else if (y_pos_r_temp.second <= ff_y_e)
            {
                valid_y = 1;
            }
            else
            {
                valid_y = 0;
            }
        }

        if (valid_x == 1 && valid_y == 1)
        {
            result_group.push_back(id + _FF_D_OFFSET);
        }
        else
        {
            x_pos_r_temp = x_ori;
            y_pos_r_temp = y_ori;
        }
    }
    // std::cout << "findmaximal is completed" << std::endl;
    x_pos_r = x_pos_r_temp;
    y_pos_r = y_pos_r_temp;

    return result_group;
}

void Solver::Solver::drawpic(const string &s)
{
    vector<size_t> FFIDs;
    FFIDs.reserve(_FF_D_arr.size());
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        FFIDs.push_back(i);
    }

    vector<size_t> GateIDs;
    FFIDs.reserve(_Gate_arr.size());
    for (size_t i = 0; i < _Gate_arr.size(); i++)
    {
        GateIDs.push_back(i);
    }

    vector<size_t> NetIDs;

    _ptr_GlobalPlacer->placementVisualization(s, FFIDs, GateIDs, NetIDs);
}

void Solver::Solver::solve_by_window()
{
    evaluate("init Metrix.txt");
    // drawpic("pic_init.plt");

    solve_initbuild();
    // solve_findfeasible();

    std::map<size_t, size_t> MBFFSfinal;

    double width = _ptr_Parser->_dieRx - _ptr_Parser->_dieLx;
    double height = _ptr_Parser->_dieRy - _ptr_Parser->_dieLy;
    /////////////////////////////////////////////////////////////change the parameters
    int slice = 160;
    double low_weight = 0.2;
    ////////////////////////////////////////////////////////////

    double window_w = width / slice;
    double window_h = height / slice;
    vector<vector<vector<size_t>>> FF_inwindow(slice, vector<vector<size_t>>(slice)); //[0][0] stores the id of FF in (0,0)
    // construct FF_inwindow
    for (int id = 0; id < _FF_D_arr.size(); id++)
    {
        // test
        // std::cout << "CLK.size(): " << _ClkList.size() << std::endl;
        // std::cout << "FFID: " << id << " ,CLK: " << _FF_D_arr[id].getClk() << std::endl;
        int x = std::floor(getFFPosition(&_FF_D_arr[id]).first / window_w);
        int y = std::floor(getFFPosition(&_FF_D_arr[id]).second / window_h);
        if (x == slice)
        {
            x = slice - 1;
        }
        if (y == slice)
        {
            y = slice - 1;
        }
        FF_inwindow[x][y].push_back(id);
    }

    /*
    vector<size_t> FF_remain; // store the remain ff index
    FF_remain.reserve(_FF_D_arr.size() + 1);
    std::unordered_map<size_t, size_t> val_idxmap; // {FF_D_index: FF_remain_idx}
    */

    for (int m = 0; m < slice; m++)
    {
        for (int n = 0; n < slice; n++)
        {
            vector<vector<Inst::feasible_coor>> CLKtofeas_x(_ClkList.size());
            // calculate the feasible
            windowbinfeasible(FF_inwindow[m][n]);

            for (int c = 0; c < FF_inwindow[m][n].size(); c++)
            {
                size_t ffid = FF_inwindow[m][n][c];
                // std::cout << "FFD.size: " << _FF_D_arr.size() << std::endl;
                // std::cout << "ffid: " << ffid << std::endl;
                // std::cout << "CLK: " << _FF_D_arr[ffid].getClk()<< std::endl;
                if (_FF_D_arr[ffid].hasfeasible == 1) //_FF_D_arr[ffid].getClk() < _ClkList.size() &&
                {
                    CLKtofeas_x.at(_FF_D_arr[ffid].getClk()).push_back(_FF_D_arr[ffid].fea_x_s);
                    CLKtofeas_x[_FF_D_arr[ffid].getClk()].push_back(_FF_D_arr[ffid].fea_x_e);
                }
            }

            for (size_t k = 0; k < _ClkList.size(); k++)
            {
                vector<Inst::feasible_coor> feas_x_clk;
                feas_x_clk = CLKtofeas_x[k];

                std::sort(feas_x_clk.begin(), feas_x_clk.end(), compareByPosVal_2);
                // std::sort(feas_y_clk.begin(), feas_y_clk.end(), compareByPosVal_2);

                size_t esssential_ff;
                // find maximal clique
                while (!feas_x_clk.empty())
                {
                    vector<size_t> ff_group;

                    for (int i = 1; i < feas_x_clk.size(); i++)
                    {

                        ff_group.push_back(feas_x_clk.at(i - 1).FF_id + _FF_D_OFFSET);

                        if ((feas_x_clk.at(i - 1).type) == 0 && (feas_x_clk.at(i).type == 1)) // 0 is start, 1 is end
                        {
                            // ff_group 有一個 essential，不知道位置會在vector的哪, result_group有，essential已經在裡面了，不用pushback
                            // record the essential is i

                            esssential_ff = feas_x_clk.at(i).FF_id + _FF_D_OFFSET;

                            // consider the y part
                            pair<double, double> x_pos_r, y_pos_r; // final region
                            vector<size_t> final_group;
                            vector<size_t> result_group;
                            // find the position, 還原成正常座標

                            // x_pos_r and y_pos_r is pass by ref and get the pos value by this function
                            result_group = solve_findmaximal(ff_group, esssential_ff, x_pos_r, y_pos_r);

                            /// another way: change the essential ff
                            /// new part to get a more flexible feasible region
                            vector<size_t> new_result_group;
                            new_result_group.reserve(4);
                            if (result_group.size() > 4)
                            {
                                esssential_ff = result_group[0];

                                for (int id = 0; id < 4; id++)
                                {
                                    new_result_group.push_back(result_group[id]);
                                }
                                result_group = solve_findmaximal(new_result_group, esssential_ff, x_pos_r, y_pos_r);
                                // std::cout << "group size: " << result_group.size() << std::endl;
                                assert(result_group.size() == 4);
                            }
                            ///////////////////////////////////////////////////

                            // preplace and slack release
                            pair<double, double> pos;
                            double x_ori = (x_pos_r.first + x_pos_r.second) * 0.5; // choose the middle point
                            double y_ori = (y_pos_r.first + y_pos_r.second) * 0.5;
                            // change to original coordinate (x'=y+x, y'=y-x)
                            pos.first = (x_ori - y_ori) * 0.5;  // x=(x'-y')/2
                            pos.second = (x_ori + y_ori) * 0.5; // y=(x'+y')/2
                            // std::cout << "real feasible center" << pos.first << " ," << pos.second << std::endl;

                            pair<double, double> lowest_pos; // left lower of rotate coor
                            x_ori = x_pos_r.first;
                            y_ori = y_pos_r.first;
                            lowest_pos.first = (x_ori - y_ori) * 0.5;  // x=(x'-y')/2
                            lowest_pos.second = (x_ori + y_ori) * 0.5; // y=(x'+y')/2

                            pair<double, double> leftmost_pos; // left upper of rotate coor
                            x_ori = x_pos_r.first;
                            y_ori = y_pos_r.second;
                            leftmost_pos.first = (x_ori - y_ori) * 0.5;  // x=(x'-y')/2
                            leftmost_pos.second = (x_ori + y_ori) * 0.5; // y=(x'+y')/2

                            double w = m * window_w;
                            double w_r = (m + 1) * window_w;
                            double h = n * window_h;
                            double h_u = (n + 1) * window_h;

                            pair<double, double> decision_pos;
                            decision_pos.first = (low_weight) * (lowest_pos.first) + (1 - low_weight) * (pos.first);
                            decision_pos.second = (low_weight) * (lowest_pos.second) + (1 - low_weight) * (pos.second);

                            final_group = prePlace(result_group, esssential_ff, pos);
                            // std::cout << "finalgroupsize: " << final_group.size() << std::endl;

                            auto it2 = MBFFSfinal.find(final_group.size());
                            if (it2 == MBFFSfinal.end())
                            {
                                MBFFSfinal.insert({final_group.size(), 1});
                            }
                            else
                            {
                                MBFFSfinal.at(final_group.size()) += 1;
                            }
                            // std::cout << "before feasible" << std::endl;
                            //  feasible_cal(final_group, k, w, w_r, h, h_u);
                            windowbinfeasible(FF_inwindow[m][n]);
                            // std::cout << "after feasiblel" << std::endl;

                            break;
                        }
                    }

                    // delete the grouped member in feas_x_clk

                    for (size_t ii = 0; ii < feas_x_clk.size(); ii++)
                    {
                        // std::cout << "A: " << feas_x_clk.at(ii).FF_id << " ," << _FF_D_arr[feas_x_clk.at(ii).FF_id].grouped << std::endl;

                        // std::cout << "FF" << feas_x_clk.at(i).FF_id << " group or not: " << _FF_D_arr[feas_x_clk.at(i).FF_id].grouped << std::endl;
                        if (_FF_D_arr[feas_x_clk.at(ii).FF_id].grouped == 1)
                        {

                            feas_x_clk.erase(feas_x_clk.begin() + ii);
                            ii = ii - 1;
                        }
                    }
                }
                // std::cout << "CLKlist " << k << " is finished" << std::endl;
            }

            // string s = std::to_string(10 * w / window_w + h / window_h);
            // drawpic("pic_" + s + ".plt");
            std::cout << "M:" << m << " N:" << n << std::endl;
        }
    }

    MBFFSfinal.insert({0, 0});
    for (int idx = 0; idx < _FF_D_arr.size(); idx++)
    {
        if (_FF_D_arr[idx].grouped == 0)
        {
            MBFFSfinal.at(0) += 1;
        }
    }
    std::cout << "After preplacing" << std::endl;
    for (const auto &pair : MBFFSfinal)
    {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }

    // drawpic("pic_e.plt");
    evaluate("before legal Metrix.txt");
    legalize();
    evaluate("after legal Metrix.txt");
    // drawpic("pic_final.plt");
    std::cout << "Solver is completed !" << std::endl;
}

bool Solver::Solver::FFinwindow(Inst::FF_D *ptr_FFD, double &LLx, double &LLy, double &URx, double &URy)
{
    pair<double, double> pos = ptr_FFD->getPosition();
    if (pos.first >= LLx && pos.first < URx && pos.second >= LLy && pos.second < URy)
    {
        return true;
    }

    return false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Solver::Solver::printOutput(const string &outFileName)
{
    // step1: calculate how many instance need to be output
    vector<bool> isChecked(_FF_D_arr.size(), false);
    int instNum = 0;
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        if (isChecked[i])
        {
            continue;
        }
        instNum++;
        for (const auto &memberID : _FF_D_arr[i].grouped_member)
        {
            size_t id = memberID - _FF_D_OFFSET;
            if (id != i)
            {
                isChecked[id] = true;
            }
        }
    }

    // step2: output the file
    // reset the checked list
    for (size_t i = 0; i < isChecked.size(); i++)
    {
        isChecked[i] = false;
    }
    // start to output
    std::fstream outputFile;
    outputFile.open(outFileName, std::fstream::out);
    outputFile << std::setprecision(15);
    outputFile << "CellInst " << instNum << "\n";
    // output the Inst list
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        if (isChecked[i])
        {
            continue;
        }
        string instName = _FF_D_arr[i].getName().substr(0, _FF_D_arr[i].getName().find('/'));
        pair<double, double> position = getFFPosition(&_FF_D_arr[i]);
        outputFile << "Inst " << instName << " " << _ptr_Parser->_flipflopLib[_FF_D_arr[i].FF_type].Name
                   << " " << position.first << " " << position.second << "\n";
        for (const auto &memberID : _FF_D_arr[i].grouped_member)
        {
            size_t id = memberID - _FF_D_OFFSET;
            if (id != i)
            {
                isChecked[id] = true;
            }
        }
    }
    // output the pin mapping
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        outputFile << _FF_D_arr[i].getOriName() << " map " << _FF_D_arr[i].getName() << "\n";
        outputFile << _FF_Q_arr[i].getOriName() << " map " << _FF_Q_arr[i].getName() << "\n";
        // find the clk name for origin and new pin
        string newInstName = _FF_D_arr[i].getName().substr(0, _FF_D_arr[i].getName().find('/'));
        string oriInstName = _FF_D_arr[i].getOriName().substr(0, _FF_D_arr[i].getOriName().find('/'));
        string oriCLKName;
        string newCLKName;
        for (const auto &name : _ptr_Parser->_flipflopLib[_FF_D_arr[i].OriFF_type].PinName)
        {
            if (name.find("clk") != std::string::npos || name.find("CLK") != std::string::npos)
            {
                oriCLKName = name;
                break;
            }
        }
        for (const auto &name : _ptr_Parser->_flipflopLib[_FF_D_arr[i].FF_type].PinName)
        {
            if (name.find("clk") != std::string::npos || name.find("CLK") != std::string::npos)
            {
                newCLKName = name;
                break;
            }
        }
        outputFile << oriInstName + '/' + oriCLKName << " map " << newInstName + '/' + newCLKName << "\n";
    }

    std::cout << "Output successfully!" << std::endl;
    outputFile.close();
}

void Solver::Solver::test()
{
    // initSolver();
    // legal test
    evaluate("init Metrix.txt");
    solve_initbuild();
    solve_findfeasible();
    _ptr_evaluator->setSolverPtr(this);

    std::cout << "test size: " << _ID_to_instance[1011197]->getRelatedNet().size() << std::endl;
    std::cout << "FFDOFFSET: " << _FF_D_OFFSET << std::endl;
    std::cout << "FFDsize: " << _FF_D_arr.size() << std::endl;
    assert(_ID_to_instance[1011197]->getType() == Inst::INST_FF_D);

    for (int i = 0; i < _FF_D_arr.size(); i++)
    {
        if (_ID_to_instance[i + _FF_D_OFFSET]->getRelatedNet().size() == 0)
        {
            std::cout << "FF_" << i << ", id: " << i + _FF_D_OFFSET << " size: " << _ID_to_instance[i + _FF_D_OFFSET]->getRelatedNet().size() << std::endl;
        }
    }
}

vector<size_t> Solver::Solver::prePlace(const vector<size_t> &ff_group, size_t essential_ID, pair<double, double> pos)
{
    // std::cout << "start preplace" << std::endl;
    //  return the id get grouped
    vector<size_t> ID_group;
    ID_group.reserve(ff_group.size());
    assert(essential_ID >= _FF_D_OFFSET);
    ID_group.push_back(essential_ID);

    // std::cout << "step1: choose the ff size from ff lib" << std::endl;
    //  step1: choose the ff size from ff lib
    int maxSize = 0;
    size_t fftype = 0;
    // here we choose the first FF with largest available size in the ff lib
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib[i].Bit > maxSize && _ptr_Parser->_flipflopLib[i].Bit <= int(ff_group.size()))
        {
            maxSize = _ptr_Parser->_flipflopLib[i].Bit;
            fftype = i;
        }
        if (_ptr_Parser->_flipflopLib[i].Bit == maxSize && _libGrade.at(fftype).avgCost >= _libGrade.at(i).avgCost)
        {
            fftype = i;
        }
    }
    assert(maxSize > 0);

    // std::cout << "step2: choose which ff be grouped" << std::endl;
    //  step2: choose which ff be grouped
    //  here we randomly choose other maxsize - 1 ff to be group
    // std::cout << "ff_group.size = " << ff_group.size() << std::endl;
    for (size_t i = 0; i < ff_group.size(); i++)
    {
        assert(ff_group[i] >= _FF_D_OFFSET);
        if (ff_group.at(i) != essential_ID)
        {
            ID_group.push_back(ff_group[i]);
            // std::cout << "ff_group[" << i << "] = " << ff_group[i] << std::endl;
        }
        if (int(ID_group.size()) == maxSize)
        {
            break;
        }
    }
    assert(int(ID_group.size()) == maxSize);

    // std::cout << "step3: place it and release slack" << std::endl;
    //  step3: place it and release slack
    pair<double, double> center(0, 0);
    vector<pair<size_t, double>> DpinHeight;
    vector<pair<size_t, double>> ID2HeightMap;
    ID2HeightMap.reserve(maxSize);
    DpinHeight.reserve(maxSize);
    // std::cout << "A" << std::endl;
    //  calculate the center of pin
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib[fftype].PinCrdnate.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib[fftype].PinName[i].find("CLK") != std::string::npos || _ptr_Parser->_flipflopLib[fftype].PinName[i].find("clk") != std::string::npos)
        {
            continue;
        }
        if (_ptr_Parser->_flipflopLib[fftype].PinName[i][0] == 'D' || _ptr_Parser->_flipflopLib[fftype].PinName[i][0] == 'd')
        {
            pair<size_t, double> id2yMap;
            id2yMap.first = i;
            id2yMap.second = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].second;
            DpinHeight.push_back(id2yMap);
        }
        center.first += _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].first / (2 * maxSize);
        center.second += _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].second / (2 * maxSize);
    }
    // std::cout << "AA" << std::endl;
    // initialize ID2HeightMap
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        pair<size_t, double> id2yMap;
        id2yMap.first = ID_group.at(i);
        // std::cout << "ID_group[i] = " << ID_group[i] << std::endl;
        id2yMap.second = findPinPosition(ID_group[i]).second;
        ID2HeightMap.push_back(id2yMap);
    }
    // std::cout << "AAA" << std::endl;
    //  sort the grouped ffs and D pins with their y pos
    std::sort(ID2HeightMap.begin(), ID2HeightMap.end(), [](pair<size_t, double> &a, pair<size_t, double> &b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    std::sort(DpinHeight.begin(), DpinHeight.end(), [](pair<size_t, double> &a, pair<size_t, double> &b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    // assign the position to each pin according to their order of height
    // FFName = instanceName of the first pin/
    // std::cout << "AAAA" << std::endl;
    string FFName = _ID_to_instance[ID2HeightMap[0].first]->getName().substr(0, _ID_to_instance[ID2HeightMap[0].first]->getName().find('/') + 1);
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        size_t FFid = ID2HeightMap[i].first;
        size_t pinID = DpinHeight[i].first;
        // std::cout << "FFid = " << FFid << " " << _FF_D_OFFSET << std::endl;
        Inst::FF_D *ptr_FF_D = &_FF_D_arr.at(FFid - _FF_D_OFFSET);
        string pinName = _ptr_Parser->_flipflopLib[fftype].PinName[pinID];
        pair<double, double> Pin2CenterDis = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[pinID];
        Pin2CenterDis.first -= center.first;
        Pin2CenterDis.second -= center.second;
        // set the information for FF_D
        ptr_FF_D->setName(FFName + pinName);
        ptr_FF_D->setPosition(pos.first + Pin2CenterDis.first, pos.second + Pin2CenterDis.second);
        ptr_FF_D->grouped_member = ID_group;
        ptr_FF_D->grouped = true;
        ptr_FF_D->FF_type = fftype;
        // set the information for FF_Q
        Inst::FF_Q *ptr_FF_Q = &_FF_Q_arr[FFid - _FF_D_OFFSET];
        if (pinName[0] == 'D')
        {
            pinName[0] = 'Q';
        }
        else
        {
            pinName[0] = 'q';
        }
        for (size_t j = 0; j < _ptr_Parser->_flipflopLib[fftype].PinName.size(); j++)
        {
            if (_ptr_Parser->_flipflopLib[fftype].PinName[j] == pinName)
            {
                pinID = j;
                break;
            }
        }
        Pin2CenterDis = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[pinID];
        Pin2CenterDis.first -= center.first;
        Pin2CenterDis.second -= center.second;
        ptr_FF_Q->setName(FFName + pinName);
        ptr_FF_Q->setPosition(pos.first + Pin2CenterDis.first, pos.second + Pin2CenterDis.second);
        ptr_FF_Q->grouped = true;
        ptr_FF_Q->FF_type = fftype;
    }

    // std::cout << "step4: release slack to connected ff" << std::endl;
    //  release slack to connected ff
    for (const auto &id : ID_group)
    {
        Inst::FF_D *ptr_FF_D = &_FF_D_arr[id - _FF_D_OFFSET];
        for (size_t i = 0; i < ptr_FF_D->faninCone.size(); i++)
        {
            pair<double, double> OriPos = ptr_FF_D->getOriPosition();
            pair<double, double> NowPos = ptr_FF_D->getPosition();
            pair<double, double> GatePinPos;
            if (ptr_FF_D->outGate2Fanin[i] != _ID_to_instance.size())
            {
                GatePinPos = findPinPosition(ptr_FF_D->outGate2Fanin[i]);
            }
            else
            {
                GatePinPos = findPinPosition(ptr_FF_D->faninCone[i]);
            }
            double Oridis = std::fabs(OriPos.first - GatePinPos.first) + std::fabs(OriPos.second - GatePinPos.second);
            double Nowdis = std::fabs(NowPos.first - GatePinPos.first) + std::fabs(NowPos.second - GatePinPos.second);
            double deltaSlack = (Nowdis - Oridis) * _ptr_Parser->_displaceDelay;
            if (_FF_Q_arr[ptr_FF_D->faninCone[i] - _FF_Q_OFFSET].slack > ptr_FF_D->slack - deltaSlack)
            {
                _FF_Q_arr[ptr_FF_D->faninCone[i] - _FF_Q_OFFSET].slack = ptr_FF_D->slack - deltaSlack;
            }
        }

        Inst::FF_Q *ptr_FF_Q = &_FF_Q_arr[id - _FF_D_OFFSET];
        for (size_t i = 0; i < ptr_FF_Q->fanoutCone.size(); i++)
        {
            pair<double, double> OriPos = ptr_FF_Q->getOriPosition();
            pair<double, double> NowPos = ptr_FF_Q->getPosition();
            pair<double, double> GatePinPos;
            if (ptr_FF_Q->inGate2Fanout[i] != _ID_to_instance.size())
            {
                GatePinPos = findPinPosition(ptr_FF_Q->inGate2Fanout[i]);
            }
            else
            {
                GatePinPos = findPinPosition(ptr_FF_Q->fanoutCone[i]);
            }
            double Oridis = std::fabs(OriPos.first - GatePinPos.first) + std::fabs(OriPos.second - GatePinPos.second);
            double Nowdis = std::fabs(NowPos.first - GatePinPos.first) + std::fabs(NowPos.second - GatePinPos.second);
            double deltaSlack = (Nowdis - Oridis) * _ptr_Parser->_displaceDelay;
            if (_FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].slack > _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].getOriSlack() - deltaSlack)
            {
                _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].slack = _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].getOriSlack() - deltaSlack;
            }
        }
    }
    // std::cout << "end preplace" << std::endl;
    return ID_group;
}

vector<size_t> Solver::Solver::prePlace(const vector<size_t> &ff_group, pair<double, double> pos, size_t fflib_ID)
{
    // std::cout << "start preplace" << std::endl;
    //  return the id get grouped
    vector<size_t> ID_group;
    ID_group.reserve(ff_group.size());
    assert(ff_group.size() < _ptr_Parser->_flipflopLib[fflib_ID].Bit);

    // std::cout << "step1: choose the ff size from ff lib" << std::endl;
    //  step1: choose the ff size from ff lib
    int maxSize = _ptr_Parser->_flipflopLib[fflib_ID].Bit;
    size_t fftype = fflib_ID;
    assert(maxSize > 0);

    // std::cout << "step2: choose which ff be grouped" << std::endl;
    //  step2: choose which ff be grouped
    //  here we randomly choose other maxsize - 1 ff to be group
    // std::cout << "ff_group.size = " << ff_group.size() << std::endl;
    for (size_t i = 0; i < ff_group.size(); i++)
    {
        assert(ff_group[i] >= _FF_D_OFFSET);
        ID_group.push_back(ff_group[i]);
        if (int(ID_group.size()) == maxSize)
        {
            break;
        }
    }
    assert(int(ID_group.size()) == maxSize);

    // std::cout << "step3: place it and release slack" << std::endl;
    //  step3: place it and release slack
    pair<double, double> center(0, 0);
    vector<pair<size_t, double>> DpinHeight;
    vector<pair<size_t, double>> ID2HeightMap;
    ID2HeightMap.reserve(maxSize);
    DpinHeight.reserve(maxSize);
    // std::cout << "A" << std::endl;
    //  calculate the center of pin
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib[fftype].PinCrdnate.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib[fftype].PinName[i].find("CLK") != std::string::npos || _ptr_Parser->_flipflopLib[fftype].PinName[i].find("clk") != std::string::npos)
        {
            continue;
        }
        if (_ptr_Parser->_flipflopLib[fftype].PinName[i][0] == 'D' || _ptr_Parser->_flipflopLib[fftype].PinName[i][0] == 'd')
        {
            pair<size_t, double> id2yMap;
            id2yMap.first = i;
            id2yMap.second = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].second;
            DpinHeight.push_back(id2yMap);
        }
        center.first += _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].first / (2 * maxSize);
        center.second += _ptr_Parser->_flipflopLib[fftype].PinCrdnate[i].second / (2 * maxSize);
    }
    // std::cout << "AA" << std::endl;
    // initialize ID2HeightMap
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        pair<size_t, double> id2yMap;
        id2yMap.first = ID_group.at(i);
        // std::cout << "ID_group[i] = " << ID_group[i] << std::endl;
        id2yMap.second = findPinPosition(ID_group[i]).second;
        ID2HeightMap.push_back(id2yMap);
    }
    // std::cout << "AAA" << std::endl;
    //  sort the grouped ffs and D pins with their y pos
    std::sort(ID2HeightMap.begin(), ID2HeightMap.end(), [](pair<size_t, double> &a, pair<size_t, double> &b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    std::sort(DpinHeight.begin(), DpinHeight.end(), [](pair<size_t, double> &a, pair<size_t, double> &b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    // assign the position to each pin according to their order of height
    // FFName = instanceName of the first pin/
    // std::cout << "AAAA" << std::endl;
    string FFName = _ID_to_instance[ID2HeightMap[0].first]->getName().substr(0, _ID_to_instance[ID2HeightMap[0].first]->getName().find('/') + 1);
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        size_t FFid = ID2HeightMap[i].first;
        size_t pinID = DpinHeight[i].first;
        // std::cout << "FFid = " << FFid << " " << _FF_D_OFFSET << std::endl;
        Inst::FF_D *ptr_FF_D = &_FF_D_arr.at(FFid - _FF_D_OFFSET);
        string pinName = _ptr_Parser->_flipflopLib[fftype].PinName[pinID];
        pair<double, double> Pin2CenterDis = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[pinID];
        Pin2CenterDis.first -= center.first;
        Pin2CenterDis.second -= center.second;
        // set the information for FF_D
        ptr_FF_D->setName(FFName + pinName);
        ptr_FF_D->setPosition(pos.first + Pin2CenterDis.first, pos.second + Pin2CenterDis.second);
        ptr_FF_D->grouped_member = ID_group;
        ptr_FF_D->grouped = true;
        ptr_FF_D->FF_type = fftype;
        // set the information for FF_Q
        Inst::FF_Q *ptr_FF_Q = &_FF_Q_arr[FFid - _FF_D_OFFSET];
        if (pinName[0] == 'D')
        {
            pinName[0] = 'Q';
        }
        else
        {
            pinName[0] = 'q';
        }
        for (size_t j = 0; j < _ptr_Parser->_flipflopLib[fftype].PinName.size(); j++)
        {
            if (_ptr_Parser->_flipflopLib[fftype].PinName[j] == pinName)
            {
                pinID = j;
                break;
            }
        }
        Pin2CenterDis = _ptr_Parser->_flipflopLib[fftype].PinCrdnate[pinID];
        Pin2CenterDis.first -= center.first;
        Pin2CenterDis.second -= center.second;
        ptr_FF_Q->setName(FFName + pinName);
        ptr_FF_Q->setPosition(pos.first + Pin2CenterDis.first, pos.second + Pin2CenterDis.second);
        ptr_FF_Q->grouped = true;
        ptr_FF_Q->FF_type = fftype;
    }

    // std::cout << "end preplace" << std::endl;
    return ID_group;
}

void Solver::Solver::legalize()
{
    if (!_ptr_legalizer->legalize())
    {
        std::cout << "legalize fail!!!!!!!!!!!!!" << std::endl;
    }
}

vector<size_t> Solver::Solver::getGroupMem(Inst::FF_D *ptr) const
{
    // this function should return the ID(position in _id2ffPtr) related to this ff
    vector<size_t> output = ptr->grouped_member;
    for (size_t i = 0; i < output.size(); i++)
    {
        output[i] -= _FF_D_OFFSET;
    }
    return output;
}

void Solver::Solver::setFFPosition(Inst::FF_D *ptr, pair<double, double> &pos)
{

    for (const auto &id : ptr->grouped_member)
    {
        size_t id_arr = id - _FF_D_OFFSET;
        Inst::FF_D *ptrD = &_FF_D_arr[id_arr];
        Inst::FF_Q *ptrQ = &_FF_Q_arr[id_arr];
        string dPinName = (ptrD->getName()).substr(ptrD->getName().find('/') + 1);
        string qPinName = (ptrQ->getName()).substr(ptrQ->getName().find('/') + 1);
        vector<string> PinName = _ptr_Parser->_flipflopLib[ptrD->FF_type].PinName;
        vector<pair<coord, coord>> PinCrdnate = _ptr_Parser->_flipflopLib[ptrD->FF_type].PinCrdnate;
        for (size_t i = 0; i < PinName.size(); i++)
        {
            pair<coord, coord> pinPos = pos;
            if (dPinName == PinName[i])
            {
                ptrD->setPosition(pinPos.first + PinCrdnate[i].first, pinPos.second + PinCrdnate[i].second);
            }
            else if (qPinName == PinName[i])
            {
                ptrQ->setPosition(pos.first + PinCrdnate[i].first, pos.second + PinCrdnate[i].second);
            }
        }
    }
}

double Solver::Solver::getFFWidth(Inst::FF_D *ptr) const
{
    return _ptr_Parser->_flipflopLib[ptr->FF_type].Width;
}

double Solver::Solver::getFFHeight(Inst::FF_D *ptr) const
{
    return _ptr_Parser->_flipflopLib[ptr->FF_type].Hight;
}

pair<double, double> Solver::Solver::getFFPosition(Inst::FF_D *ptr) const
{
    string Name = ptr->getName();
    Name = Name.substr(Name.find('/') + 1);
    pair<double, double> pos = ptr->getPosition();
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib[ptr->FF_type].PinName.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib[ptr->FF_type].PinName[i] == Name)
        {
            pos.first -= _ptr_Parser->_flipflopLib[ptr->FF_type].PinCrdnate[i].first;
            pos.second -= _ptr_Parser->_flipflopLib[ptr->FF_type].PinCrdnate[i].second;
            break;
        }
    }
    return pos;
}

pair<double, double> Solver::Solver::getGateLF(Inst::Gate *ptr) const
{
    return ptr->getPosition();
}

pair<double, double> Solver::Solver::getGateUR(Inst::Gate *ptr) const
{
    pair<double, double> UR = ptr->getPosition();
    UR.first += (_ptr_Parser->_gateLib[ptr->gate_type]).Width;
    UR.second += (_ptr_Parser->_gateLib[ptr->gate_type]).Hight;
    return UR;
}

double Solver::Solver::getLowerLeftX() const
{
    double LFx = _PlaceRow[0].x;
    for (const auto &plR : _PlaceRow)
    {
        if (plR.x < LFx)
        {
            LFx = plR.x;
        }
    }
    return LFx;
}

double Solver::Solver::getLowerLeftY() const
{
    double LFy = _PlaceRow[0].y;
    for (const auto &plR : _PlaceRow)
    {
        if (plR.y < LFy)
        {
            LFy = plR.y;
        }
    }
    return LFy;
}

bool Solver::Solver::placementRowIsUniform(vector<struct PlacementRow> &_PLR) const
{
    double x = _PLR[0].x;
    double height = _PLR[0].siteHight;
    double width = _PLR[0].siteWidth;
    double numOfSite = _PLR[0].totalNumOfSites;
    for (const auto &plR : _PLR)
    {
        if (x != plR.x)
        {
            std::cout << "the start point of placement rows aren't aligned in x direction." << std::endl;
            return false;
        }
        if (height != plR.siteHight)
        {
            std::cout << "the height of placement rows aren't same as each others." << std::endl;
            return false;
        }
        if (width != plR.siteWidth)
        {
            std::cout << "the width of placement rows aren't same as each others." << std::endl;
            return false;
        }
        if (numOfSite != plR.totalNumOfSites)
        {
            std::cout << "the # of site of placement rows aren't same as each others." << std::endl;
            return false;
        }
    }
    return true;
}

void Solver::Solver::slackDistribute(const double k)
{
    // k shold be a number between 0 and 1;
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        bool isNeg = false;
        double oriSlack = _FF_D_arr[i].getOriSlack();
        if (oriSlack < 0)
        {
            isNeg = true;
        }
        // if slack is negative, do not spread it
        if (isNeg)
        {
            _FF_D_arr[i].slack = oriSlack;
        }
        else
        {
            _FF_D_arr[i].slack = k * oriSlack;
            for (const auto &ffq_ID : _FF_D_arr[i].faninCone)
            {
                size_t id = ffq_ID - _FF_Q_OFFSET;
                if (_FF_Q_arr[id].slack == 0)
                {
                    _FF_Q_arr[id].slack = k * oriSlack;
                }
                else
                {
                    if (_FF_Q_arr[id].slack > k * oriSlack)
                    {
                        _FF_Q_arr[id].slack = k * oriSlack;
                    }
                }
            }
        }
    }
}

pair<double, double> Solver::Solver::findPinPosition(const size_t &id) const
{
    if (_ID_to_instance[id]->getType() == Inst::InstType::INST_GATE)
    {
        Inst::Gate *ptr_gate = _GID_to_ptrGate_map[id - _GATE_OFFSET];
        pair<double, double> pos = ptr_gate->pinPosition[id - ptr_gate->PIN_OFFSET];
        return pos;
    }
    else
    {
        return _ID_to_instance[id]->getPosition();
    }
}

size_t Solver::Solver::name2ID(string &name) const
{
    return _Name_to_ID.at(name);
}

vector<pair<double, double>> Solver::Solver::getAdjacentPinPosition(size_t &id) const
{
    vector<pair<double, double>> pinPosition;
    if (_ID_to_instance[id]->getType() == Inst::INST_FF_Q)
    {
        for (const auto &netID : _ID_to_instance[id]->getRelatedNet())
        {
            pinPosition.reserve(_NetList[netID].size());
            for (const auto &relatePinID : _NetList[netID])
            {
                if (relatePinID != id)
                {
                    pinPosition.push_back(findPinPosition(relatePinID));
                }
            }
        }
    }
    else if (_ID_to_instance[id]->getType() == Inst::INST_FF_D)
    {
        for (const auto &netID : _ID_to_instance[id]->getRelatedNet())
        {
            pinPosition.reserve(_NetList[netID].size());
            for (const auto &relatePinID : _NetList[netID])
            {
                if (_ID_to_instance[relatePinID]->getType() != Inst::INST_FF_D)
                {
                    if (_ID_to_instance[relatePinID]->getName().find("OUTPUT") != std::string::npos) // ignore the output pin
                    {
                        continue;
                    }
                    pinPosition.push_back(findPinPosition(relatePinID));
                }
            }
        }
    }
    return pinPosition;
}

void Solver::Solver::findMaxSlack()
{
    // calculate the distance of critical path and convert to slack for every FF_D
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        size_t numOfConnectedFF = _FF_D_arr[i].faninCone.size();
        if (numOfConnectedFF == 0) // for the case _FF_D is directly connect to PI
        {
            double distance = 0;
            if (_FF_D_arr[i].getRelatedNet().size() != 0)
            {
                for (const auto &pinID : _NetList[_FF_D_arr[i].getRelatedNet().at(0)])
                {
                    if (_ID_to_instance[pinID]->getType() == Inst::INST_PIO)
                    {
                        pair<double, double> PIPos = findPinPosition(pinID);
                        pair<double, double> FFPos = _FF_D_arr[i].getPosition();
                        distance += (std::fabs(PIPos.first - FFPos.first) + std::fabs(PIPos.second - FFPos.second));
                        break;
                    }
                }
            }
            _FF_D_arr[i].maxSlack = (distance * _ptr_Parser->_displaceDelay) + _FF_D_arr[i].getOriSlack();
        }
        else
        {
            double longestPath = 0;
            for (size_t j = 0; j < numOfConnectedFF; j++)
            {
                double distance = 0;
                if (_FF_D_arr[i].inGate2Fanin[j] != _ID_to_instance.size()) // if the former FF is not directly connect to FF_D
                {
                    distance += _ptr_STAEngine->getDistance(_FF_D_arr[i].inGate2Fanin[j], _FF_D_arr[i].outGate2Fanin[j]);
                    // calculate the distance between outPin of gate and FF_D
                    pair<double, double> gatePos = findPinPosition(_FF_D_arr[i].outGate2Fanin[j]);
                    pair<double, double> FFPos = _FF_D_arr[i].getPosition();
                    distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                    // calculate the distance between inPin of gate and FF_D
                    gatePos = findPinPosition(_FF_D_arr[i].inGate2Fanin[j]);
                    FFPos = findPinPosition(_FF_D_arr[i].faninCone[j]);
                    distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                    // add the effect of Qpin delay
                    distance += (_ptr_Parser->_flipflopLib[_FF_Q_arr[_FF_D_arr[i].faninCone[j] - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay;

                    // if the path is longer than the former longest path, record it
                    if (distance > longestPath)
                    {
                        longestPath = distance;
                    }
                }
                else // the former FF is directly connect to FF_D
                {
                    pair<double, double> FFPos = _FF_D_arr[i].getPosition();
                    pair<double, double> FFBeforePos = findPinPosition(_FF_D_arr[i].faninCone[j]);
                    distance += (std::fabs(FFBeforePos.first - FFPos.first) + std::fabs(FFBeforePos.second - FFPos.second));
                    // add the effect of Qpin delay
                    distance += (_ptr_Parser->_flipflopLib[_FF_Q_arr[_FF_D_arr[i].faninCone[j] - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay;

                    // if the path is longer than the former longest path, record it
                    if (distance > longestPath)
                    {
                        longestPath = distance;
                    }
                }
            }

            // convert the longestPath into slack
            _FF_D_arr[i].maxSlack = (longestPath * _ptr_Parser->_displaceDelay) + _FF_D_arr[i].getOriSlack();
        }
    }
}

vector<pair<size_t, double>> Solver::Solver::getSlack2ConnectedFF(const size_t &id, bool considerQpinDelay, double weight) // the input should be global ID
{
    vector<pair<size_t, double>> output;
    output.reserve(128);
    // If the input is FF_D, we only return one slack, which is the slack of critical path
    if (_ID_to_instance[id]->getType() == Inst::INST_FF_D)
    {
        size_t numOfConnectedFF = _FF_D_arr[id - _FF_D_OFFSET].faninCone.size();
        // case1: _FF_D is directly connect to PI
        if (numOfConnectedFF == 0) // for the case _FF_D is directly connect to PI
        {
            pair<size_t, double> initPair;
            double distance = 0;
            for (const auto &pinID : _NetList[_ID_to_instance[id]->getRelatedNet()[0]])
            {
                if (_ID_to_instance[pinID]->getType() == Inst::INST_PIO)
                {
                    initPair.first = pinID;
                    pair<double, double> PIPos = findPinPosition(pinID);
                    pair<double, double> FFPos = _ID_to_instance[id]->getPosition();
                    distance += (std::fabs(PIPos.first - FFPos.first) + std::fabs(PIPos.second - FFPos.second));
                    break;
                }
            }
            initPair.second = _FF_D_arr[id - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay);
            output.push_back(initPair);
            return output;
        }

        // case2: _FF_D is not directly connect to PI
        double slackRemain = _FF_D_arr[id - _FF_D_OFFSET].maxSlack;
        bool isdirect2FF = false;
        for (size_t i = 0; i < numOfConnectedFF; i++)
        {
            double distance = 0;
            if (_FF_D_arr[id - _FF_D_OFFSET].inGate2Fanin[i] != _ID_to_instance.size()) // if the former FF is not directly connect to FF_D
            {
                distance += _ptr_STAEngine->getDistance(_FF_D_arr[id - _FF_D_OFFSET].inGate2Fanin[i], _FF_D_arr[id - _FF_D_OFFSET].outGate2Fanin[i]);
                // calculate the distance between outPin of gate and FF_D
                pair<double, double> gatePos = findPinPosition(_FF_D_arr[id - _FF_D_OFFSET].outGate2Fanin[i]);
                pair<double, double> FFPos = _FF_D_arr[id - _FF_D_OFFSET].getPosition();
                distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                // calculate the distance between inPin of gate and FF_D
                gatePos = findPinPosition(_FF_D_arr[id - _FF_D_OFFSET].inGate2Fanin[i]);
                FFPos = findPinPosition(_FF_D_arr[id - _FF_D_OFFSET].faninCone[i]);
                distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                // add the effect of Qpin delay
                distance += (_ptr_Parser->_flipflopLib[_FF_Q_arr[_FF_D_arr[id - _FF_D_OFFSET].faninCone[i] - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay;

                double s;
                // check whether the connected FF has been grouped
                if (_FF_Q_arr[_FF_D_arr[id - _FF_D_OFFSET].faninCone[i] - _FF_Q_OFFSET].grouped)
                {
                    // if is grouped give all slack to FF_D
                    s = _FF_D_arr[id - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay);
                }
                else
                {
                    // if is not grouped give half slack to FF_D
                    s = (_FF_D_arr[id - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay)) / 2;
                    if (considerQpinDelay)
                    {
                        s -= _avgQpinDelayDiff * weight;
                    }
                }

                if (s < slackRemain)
                {
                    slackRemain = s;
                    isdirect2FF = false;
                }
            }
            else // the former FF is directly connect to FF_D
            {
                pair<double, double> FFPos = _FF_D_arr[id - _FF_D_OFFSET].getPosition();
                pair<double, double> FFBeforePos = findPinPosition(_FF_D_arr[id - _FF_D_OFFSET].faninCone[i]);
                distance += (std::fabs(FFBeforePos.first - FFPos.first) + std::fabs(FFBeforePos.second - FFPos.second));
                // add the effect of Qpin delay
                distance += (_ptr_Parser->_flipflopLib[_FF_Q_arr[_FF_D_arr[id - _FF_D_OFFSET].faninCone[i] - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay;

                double s;
                // check whether the connected FF has been grouped
                if (_FF_Q_arr[_FF_D_arr[id - _FF_D_OFFSET].faninCone[i] - _FF_Q_OFFSET].grouped)
                {
                    // if is grouped give all slack to FF_D
                    s = _FF_D_arr[id - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay);
                }
                else
                {
                    // if is not grouped give half slack to FF_D
                    s = (_FF_D_arr[id - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay)) / 2;
                    s -= _avgQpinDelayDiff * weight;
                }

                if (s < slackRemain)
                {
                    slackRemain = s;
                    isdirect2FF = true;
                }
            }
        }
        if (!isdirect2FF)
        {
            pair<size_t, double> initPair(_FF_D_arr[id - _FF_D_OFFSET].outGate2Fanin[0], slackRemain);
            output.push_back(initPair);
        }
        else
        {
            pair<size_t, double> initPair(_FF_D_arr[id - _FF_D_OFFSET].faninCone[0], slackRemain);
            output.push_back(initPair);
        }
        return output;
    }
    else if (_ID_to_instance[id]->getType() == Inst::INST_FF_Q)
    {
        size_t numOfConnectedFF = _FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone.size();
        size_t previousInGate;
        double minSlack;
        for (size_t i = 0; i < numOfConnectedFF; i++)
        {
            if (i == 0)
            {
                previousInGate = _FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i];
                minSlack = _FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack;
            }
            // if previousInGate != next inGate add result to output
            if (previousInGate != _FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i])
            {
                if (previousInGate != _ID_to_instance.size())
                {
                    pair<size_t, double> initPair(previousInGate, minSlack);
                    output.push_back(initPair);
                }
                else
                {
                    pair<size_t, double> initPair(_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i - 1], minSlack);
                    output.push_back(initPair);
                }
                previousInGate = _FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i];
                minSlack = _FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack;
            }

            double distance = 0;
            if (_FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i] != _ID_to_instance.size()) // if the later FF is not directly connect to FF_Q
            {
                distance += _ptr_STAEngine->getDistance(_FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i], _FF_Q_arr[id - _FF_Q_OFFSET].outGate2Fanout[i]);
                // calculate the distance between outPin of gate and FF_D
                pair<double, double> gatePos = findPinPosition(_FF_Q_arr[id - _FF_Q_OFFSET].outGate2Fanout[i]);
                pair<double, double> FFPos = findPinPosition(_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i]);
                distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                // calculate the distance between inPin of gate and FF_Q
                gatePos = findPinPosition(_FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout[i]);
                FFPos = _FF_Q_arr[id - _FF_Q_OFFSET].getPosition();
                distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
                // add the effect of Qpin delay
                distance += ((_ptr_Parser->_flipflopLib[_FF_Q_arr[id - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay);

                // check whether the connected FF has been grouped
                if (_FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].grouped)
                {
                    // if is grouped give all slack to FF_D
                    double s = _FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay);
                    if (s < minSlack)
                    {
                        minSlack = s;
                    }
                }
                else
                {
                    // if is not grouped give half slack to FF_D
                    double s = (_FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay)) / 2;
                    if (s < minSlack)
                    {
                        minSlack = s;
                    }
                }
            }
            else // the later FF is directly connect to FF_D
            {
                pair<double, double> FFPos = _FF_Q_arr[id - _FF_Q_OFFSET].getPosition();
                pair<double, double> FFBeforePos = findPinPosition(_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i]);
                distance += (std::fabs(FFBeforePos.first - FFPos.first) + std::fabs(FFBeforePos.second - FFPos.second));
                // add the effect of Qpin delay
                distance += (_ptr_Parser->_flipflopLib[_FF_Q_arr[id - _FF_Q_OFFSET].FF_type].PinDelay) / _ptr_Parser->_displaceDelay;

                // check whether the connected FF has been grouped
                if (_FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].grouped)
                {
                    // if is grouped give all slack to FF_D
                    double s = (_FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay));
                    if (s < minSlack)
                    {
                        minSlack = s;
                    }
                }
                else
                {
                    // if is not grouped give half slack to FF_D
                    double s = ((_FF_D_arr[_FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone[i] - _FF_D_OFFSET].maxSlack - (distance * _ptr_Parser->_displaceDelay)) / 2);
                    if (s < minSlack)
                    {
                        minSlack = s;
                    }
                }
            }

            // if i == numOfConnectedFF - 1, record to output
            if (i == numOfConnectedFF - 1)
            {
                if (previousInGate != _ID_to_instance.size())
                {
                    pair<size_t, double> initPair(previousInGate, minSlack);
                    output.push_back(initPair);
                }
                else
                {
                    pair<size_t, double> initPair(_FF_Q_arr.at(id - _FF_Q_OFFSET).fanoutCone.at(i), minSlack);
                    output.push_back(initPair);
                }
            }
        }
    }
    return output;
}

void Solver::Solver::findFaninout4all(const vector<pair<size_t, std::map<size_t, double>>> &_distanceList)
{
    // step1: categorize which pin is critial
    // initialize isInPinCritical isOutPinCritical
    isInPinCritical.clear();
    isInPinCritical.reserve(_FF_D_OFFSET);
    isOutPinCritical.clear();
    isOutPinCritical.reserve(_FF_D_OFFSET);
    for (size_t i = 0; i < _FF_D_OFFSET; i++)
    {

        isOutPinCritical.push_back(false);
        isInPinCritical.push_back(false);
    }
    _InPinList.reserve(_NetList.size());
    _OutPinList.reserve(_NetList.size());
    for (size_t i = 0; i < _NetList.size(); i++)
    {
        vector<size_t> initVec;
        initVec.reserve(_NetList[i].size());
        _InPinList.push_back(initVec);
        _OutPinList.push_back(initVec);
        for (const auto &id : _NetList[i])
        {
            if (isInPin(id))
            {
                _InPinList[i].push_back(id);
            }
            else if (isOutPin(id))
            {
                _OutPinList[i].push_back(id);
            }
        }
    }
    for (size_t i = 0; i < _NetList.size(); i++)
    {
        if (_OutPinList[i].size() > 0)
        {
            bool haveFF_D = false;
            vector<size_t> ffConnect;
            ffConnect.reserve(_NetList[i].size());
            for (const auto &pinID : _NetList[i])
            {
                if (_FF_D_OFFSET <= pinID && pinID < _FF_Q_OFFSET)
                {
                    haveFF_D = true;
                    ffConnect.push_back(pinID);
                }
            }
            if (haveFF_D)
            {
                for (const auto &outID : _OutPinList[i])
                {
                    isOutPinCritical[outID] = true;
                    FFConnect2CriticalPin[outID] = ffConnect;
                }
            }
        }
        else if (_OutPinList[i].size() == 0)
        {
            bool haveFF_Q = false;
            vector<size_t> ffConnect;
            ffConnect.reserve(_NetList[i].size());
            for (const auto &pinID : _NetList[i])
            {
                if (_FF_Q_OFFSET <= pinID)
                {
                    haveFF_Q = true;
                    ffConnect.push_back(pinID);
                }
            }
            if (haveFF_Q)
            {
                for (const auto &inID : _InPinList[i])
                {
                    isInPinCritical[inID] = true;
                    FFConnect2CriticalPin[inID] = ffConnect;
                }
            }
        }
    }

    // step2: constuct inPin2Out and outPin2In
    // initialize inPin2Out and outPin2In
    for (size_t i = 0; i < _FF_D_OFFSET; i++)
    {
        if (isInPinCritical[i])
        {
            vector<size_t> initVec;
            initVec.reserve(512);
            pair<size_t, vector<size_t>> initPair(i, initVec);
            inPin2Out.push_back(initPair);
        }
        else if (isOutPinCritical[i])
        {
            vector<size_t> initVec;
            initVec.reserve(512);
            pair<size_t, vector<size_t>> initPair(i, initVec);
            outPin2In.push_back(initPair);
        }
    }
    // update _OutPin2PositionMap and _InPin2PositionMap
    for (int i = 0; i < int(inPin2Out.size()); i++)
    {
        _InPin2PositionMap[inPin2Out.at(i).first] = i;
    }
    for (int i = 0; i < int(outPin2In.size()); i++)
    {
        _OutPin2PositionMap[outPin2In.at(i).first] = i;
    }
    // constuct constuct inPin2Out and outPin2In
    for (const auto &p : _distanceList)
    {
        if (isInPinCritical.at(p.first))
        {
            for (const auto &outlist : p.second)
            {
                if (outlist.second != -1)
                {
                    inPin2Out.at(_InPin2PositionMap.at(p.first)).second.push_back(outlist.first);
                    outPin2In.at(_OutPin2PositionMap.at(outlist.first)).second.push_back(p.first);
                }
            }
        }
    }

    // step3: use the table to find fanin and fanout cone for FF_D and FF_Q
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        for (const auto &netID : _FF_D_arr[i].getRelatedNet())
        {
            // bool hasOutPin = false;
            // bool hasFF_Q = false;
            for (const auto &pinID : _NetList[netID])
            {
                if (isOutPin(pinID))
                {
                    // hasOutPin = true;
                    for (const auto &inPinID : outPin2In.at(_OutPin2PositionMap.at(pinID)).second)
                    {
                        for (const auto &ffConnectID : FFConnect2CriticalPin.at(inPinID))
                        {
                            _FF_D_arr[i].faninCone.push_back(ffConnectID);
                            _FF_D_arr[i].inGate2Fanin.push_back(inPinID);
                            _FF_D_arr[i].outGate2Fanin.push_back(pinID);
                        }
                    }
                    break;
                }
                else if (pinID >= _FF_Q_OFFSET)
                {
                    // hasFF_Q = true;
                    _FF_D_arr[i].faninCone.push_back(pinID);
                    _FF_D_arr[i].inGate2Fanin.push_back(_ID_to_instance.size());
                    _FF_D_arr[i].outGate2Fanin.push_back(_ID_to_instance.size());
                    break;
                }
            }
        }
    }
    for (size_t i = 0; i < _FF_Q_arr.size(); i++)
    {
        // vector<bool> isChecked(_ID_to_instance.size(), false);
        for (const auto &netID : _FF_Q_arr[i].getRelatedNet())
        {
            for (const auto &pinID : _NetList[netID])
            {
                if (isInPin(pinID))
                {
                    for (const auto &outPinID : inPin2Out.at(_InPin2PositionMap.at(pinID)).second)
                    {
                        for (const auto &ffConnectID : FFConnect2CriticalPin.at(outPinID))
                        {

                            _FF_Q_arr[i].fanoutCone.push_back(ffConnectID);
                            _FF_Q_arr[i].inGate2Fanout.push_back(pinID);
                            _FF_Q_arr[i].outGate2Fanout.push_back(outPinID);
                        }
                    }
                }
                else if (pinID >= _FF_D_OFFSET && pinID < _FF_Q_OFFSET)
                {
                    _FF_Q_arr[i].fanoutCone.push_back(pinID);
                    _FF_Q_arr[i].inGate2Fanout.push_back(_ID_to_instance.size());
                    _FF_Q_arr[i].outGate2Fanout.push_back(_ID_to_instance.size());
                }
            }
        }

        // sort the fanoutCone by inGate
        vector<pair<size_t, size_t>> in2PosInFanout; // pair<ingateID, position in fanout>
        in2PosInFanout.reserve(_FF_Q_arr[i].inGate2Fanout.size());
        for (size_t k = 0; k < _FF_Q_arr[i].inGate2Fanout.size(); k++)
        {
            pair<size_t, size_t> initPair(_FF_Q_arr[i].inGate2Fanout[k], k);
            in2PosInFanout.push_back(initPair);
        }
        std::sort(in2PosInFanout.begin(), in2PosInFanout.end(), [](pair<size_t, size_t> &a, pair<size_t, size_t> &b)
                  { return a.first < b.first; });
        vector<size_t> c = _FF_Q_arr[i].fanoutCone;
        for (size_t k = 0; k < _FF_Q_arr[i].inGate2Fanout.size(); k++)
        {
            _FF_Q_arr[i].fanoutCone.at(k) = c.at(in2PosInFanout.at(k).second);
        }
        c = _FF_Q_arr[i].inGate2Fanout;
        for (size_t k = 0; k < _FF_Q_arr[i].inGate2Fanout.size(); k++)
        {
            _FF_Q_arr[i].inGate2Fanout.at(k) = c.at(in2PosInFanout.at(k).second);
        }
        c = _FF_Q_arr[i].outGate2Fanout;
        for (size_t k = 0; k < _FF_Q_arr[i].inGate2Fanout.size(); k++)
        {
            _FF_Q_arr[i].outGate2Fanout.at(k) = c.at(in2PosInFanout.at(k).second);
        }
    }

    // step4: release the memory
    inPin2Out.clear();
    inPin2Out.shrink_to_fit();
    outPin2In.clear();
    outPin2In.shrink_to_fit();
    isInPinCritical.clear();
    isInPinCritical.shrink_to_fit();
    isOutPinCritical.clear();
    isOutPinCritical.shrink_to_fit();
    FFConnect2CriticalPin.clear();
    _OutPin2PositionMap.clear();
    _InPin2PositionMap.clear();
    _InPinList.clear();
    _InPinList.shrink_to_fit();
    _OutPinList.clear();
    _OutPinList.shrink_to_fit();
}

bool Solver::Solver::isInPin(const size_t &id) const
{
    if (_ID_to_instance[id]->getType() == Inst::INST_GATE)
    {
        Inst::Gate *ptrGate = _GID_to_ptrGate_map[id - _GATE_OFFSET];
        if (ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("IN") != std::string::npos || ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("in") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

bool Solver::Solver::isOutPin(const size_t &id) const
{
    if (_ID_to_instance[id]->getType() == Inst::INST_GATE)
    {
        Inst::Gate *ptrGate = _GID_to_ptrGate_map[id - _GATE_OFFSET];
        if (ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("OUT") != std::string::npos || ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("out") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

void Solver::Solver::evaluate(const string &fileName)
{
    _ptr_evaluator->setSolverPtr(this);
    _ptr_evaluator->evaluate(fileName);
}
