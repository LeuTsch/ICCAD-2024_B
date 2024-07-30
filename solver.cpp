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

#include "inst.h"
#include "solver.h"
#include "STAEngine.h"
#include "legalizer.h"

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
            for (const auto &PinName : netPinList)
            {
                if (PinName.find("OUT") != std::string::npos || PinName.find("/") == std::string::npos || PinName.find("out") != std::string::npos)
                {
                    if (PinName.find("OUTPUT") != std::string::npos)
                    {
                        continue;
                    }
                    _ClkList.push_back(PinName);
                    break;
                }
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

    // step3: establish the fanin and fanout cone for the FF
    std::cout << "step3: establish the fanin and fanout cone for the FF" << std::endl;
    for (size_t ffd_id = _FF_D_OFFSET; ffd_id < _FF_D_arr.size(); ffd_id++)
    {
        findFanin(ffd_id);
    }
    for (size_t ffq_id = _FF_Q_OFFSET; ffq_id < _FF_Q_arr.size(); ffq_id++)
    {
        findFanout(ffq_id);
    }

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

    // step7: initialize the STA engine
    std::cout << "step7: initialize the STA engine" << std::endl;
    assert(_ptr_STAEngine != nullptr);
    _ptr_STAEngine->setSolverPtr(this);
    _ptr_STAEngine->initEngine(_NetList, _ID_to_instance);

    // step8: initialize the legalizer
    std::cout << "step8: initialize the legalizer" << std::endl;
    assert(_ptr_legalizer != nullptr);
    _ptr_legalizer->setSolverPtr(this);

    std::cout << "End initialization!!! \n"
              << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// THE NEW PART

void Solver::Solver::solve_initbuild()
{
    slackDistribute(0.6);
    std::cout << "slackdistribute is completed" << std::endl;

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
    // std::cout << "check1" << std::endl;
    // std::cout << "FFDsize: " << _FF_D_arr.size() << std::endl;
    //  build fanin fanout position
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        size_t D_id = i + _FF_D_OFFSET;
        vector<pair<double, double>> pos = getAdjacentPinPosition(D_id);
        // std::cout << "check2" << std::endl;
        //  check the size of pos is 1.
        _FF_D_arr.at(i).D_fanin_pos = pos[0];
    }

    for (size_t i = 0; i < _FF_Q_arr.size(); i++)
    {
        size_t Q_id = i + _FF_Q_OFFSET;
        _FF_Q_arr.at(i).Q_fanout_pos.reserve(_FF_Q_arr.size());
        _FF_Q_arr.at(i).Q_fanout_pos = getAdjacentPinPosition(Q_id);
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
        double D_dia_len = fabs(_FF_D_arr.at(i).Dx_pos - _FF_D_arr.at(i).D_fanin_pos.first) + fabs(_FF_D_arr.at(i).Dy_pos - _FF_D_arr.at(i).D_fanin_pos.second) + _FF_D_arr.at(i).slack / (_ptr_Parser->_displaceDelay);
        double Dx_pos_r = _FF_D_arr.at(i).Dy_pos + _FF_D_arr.at(i).Dx_pos + dist; // x'=y+x //D shift right
        double Dy_pos_r = _FF_D_arr.at(i).Dy_pos - _FF_D_arr.at(i).Dx_pos - dist; // y'=y-x
        coor_w_se Dx_s, Dx_e;
        vector<coor_w_se> dia_x_arr;
        dia_x_arr.reserve(256);
        Dx_s.pos_val = Dx_pos_r - (2 * D_dia_len);
        Dx_s.type = 0;
        Dx_s.or_ind = 0;
        Dx_e.pos_val = Dx_pos_r + (2 * D_dia_len);
        Dx_e.type = 1;
        Dx_e.or_ind = 0;

        dia_x_arr.push_back(Dx_s);
        dia_x_arr.push_back(Dx_e);

        coor_w_se Dy_s, Dy_e;
        vector<coor_w_se> dia_y_arr;
        dia_y_arr.reserve(256);
        Dy_s.pos_val = Dy_pos_r - (2 * D_dia_len);
        Dy_s.type = 0;
        Dy_s.or_ind = 0;
        Dy_e.pos_val = Dy_pos_r + (2 * D_dia_len);
        Dy_e.type = 1;
        Dy_e.or_ind = 0;

        dia_y_arr.push_back(Dy_s);
        dia_y_arr.push_back(Dy_e);

        // Q part, run through all fanout.

        double Qx_pos_r = _FF_Q_arr.at(i).Qy_pos + _FF_Q_arr.at(i).Qx_pos - dist; // Q shift left
        double Qy_pos_r = _FF_Q_arr.at(i).Qy_pos - _FF_Q_arr.at(i).Qx_pos + dist;
        for (size_t j = 0; j < _FF_Q_arr.at(i).Q_fanout_pos.size(); j++)
        {
            double Q_dia_len = fabs(_FF_Q_arr.at(i).Qx_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(i).Qy_pos - _FF_Q_arr.at(i).Q_fanout_pos.at(j).second) + _FF_Q_arr.at(i).slack / (_ptr_Parser->_displaceDelay);
            coor_w_se Qx_s, Qx_e;
            Qx_s.pos_val = Qx_pos_r - (2 * Q_dia_len);
            Qx_s.type = 0;
            Qx_s.or_ind = j + 1;
            Qx_e.pos_val = Qx_pos_r + (2 * Q_dia_len);
            Qx_e.type = 1;
            Qx_e.or_ind = j + 1;

            dia_x_arr.push_back(Qx_s);
            dia_x_arr.push_back(Qx_e);

            coor_w_se Qy_s, Qy_e;
            Qy_s.pos_val = Qy_pos_r - (2 * Q_dia_len);
            Qy_s.type = 0;
            Qy_s.or_ind = j + 1;
            Qy_e.pos_val = Qy_pos_r + (2 * Q_dia_len);
            Qy_e.type = 1;
            Qy_e.or_ind = j + 1;

            dia_y_arr.push_back(Qy_s);
            dia_y_arr.push_back(Qy_e);
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
    }

    std::cout << "findfeasible is completed" << std::endl;
}

void Solver::Solver::solve()
{
    solve_initbuild();

    solve_findfeasible();

    // clock一樣才能綁，外面掛一層迴圈跑過所有clock
    // FFD.getclock 吐 id, id == 0 ... for(int i = 0; i < _ClkList.size(); i++)

    for (size_t k = 0; k < _ClkList.size(); k++)
    { /*run through clk*/

        vector<Inst::feasible_coor> feas_x_clk;
        feas_x_clk.reserve(_FF_D_arr.size() * 2);
        vector<Inst::feasible_coor> feas_y_clk;
        feas_y_clk.reserve(_FF_D_arr.size() * 2);

        for (int i = 0; i < _FF_D_arr.size(); i++)
        {

            if ((_FF_D_arr[i].getClk() == k) && (_FF_D_arr[i].hasfeasible == 1))
            {
                // collect
                // may collect the empty value of no-feasible FF
                feas_x_clk.push_back(_FF_D_arr.at(i).fea_x_s);
                feas_x_clk.push_back(_FF_D_arr.at(i).fea_x_e);

                feas_y_clk.push_back(_FF_D_arr.at(i).fea_y_s);
                feas_y_clk.push_back(_FF_D_arr.at(i).fea_y_e);
            }
        }

        std::sort(feas_x_clk.begin(), feas_x_clk.end(), compareByPosVal_2);

        std::sort(feas_y_clk.begin(), feas_y_clk.end(), compareByPosVal_2);

        size_t esssential_ff;
        // find maximal clique
        while (!feas_x_clk.empty())
        {
            vector<size_t> ff_group;

            // std::cout << "Before clustering, the size: " << feas_x_clk.size() << std::endl;
            for (int i = 1; i < feas_x_clk.size(); i++)
            {

                ff_group.push_back(feas_x_clk.at(i - 1).FF_id + _FF_D_OFFSET);
                if ((feas_x_clk.at(i - 1).type) == 0 && (feas_x_clk.at(i).type == 1))
                {
                    // ff_group 沒有 essential, result_group有
                    // record the essential is i

                    esssential_ff = feas_x_clk.at(i).FF_id + _FF_D_OFFSET;

                    // consider the y part
                    pair<double, double> x_pos_r, y_pos_r; // final region
                    vector<size_t> final_group;
                    vector<size_t> result_group;
                    // find the position, 還原成正常座標

                    result_group = solve_findmaximal(ff_group, esssential_ff, x_pos_r, y_pos_r);

                    // std::cout << x_pos_r.first << " ," << x_pos_r.second << std::endl;

                    result_group.push_back(esssential_ff);

                    // preplace and slack release
                    pair<double, double> pos;
                    pos.first = (x_pos_r.first + x_pos_r.second) * 0.5;
                    pos.second = (y_pos_r.first + y_pos_r.second) * 0.5;
                    double x_ori = pos.first;
                    double y_ori = pos.second;
                    pos.first = (x_ori - y_ori) * 0.5; // change to original coordinate
                    // std::cout << "X: " << pos.first <<std::endl;
                    pos.second = (x_ori + y_ori) * 0.5;
                    // std::cout << "Y: " << pos.second <<std::endl;
                    //  std::cout << pos.first << " ," << pos.second << std::endl;
                    //   std::cout << "preplace begin" << std::endl;
                    final_group = prePlace(result_group, esssential_ff, pos);
                    // std::cout << "finalgroup size: " << final_group.size() <<std::endl;

                    // final_group有essential

                    // std::cout << "preplace is completed" << std::endl;
                    //  calculate the feasible
                    feasible_cal(final_group);
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
    }

    // 決定後來的DQ 的 pos, 記得Q的正方形要往左移（或D往右, Q往左）
    // preplace完 slack release, 更新DQ正方形，重新畫table

    legalize();
    std::cout << "Solver is completed !" << std::endl;
}

void Solver::Solver::feasible_cal(const vector<size_t> &final_group)
{
    // std::cout << "Hello" << std::endl;

    for (size_t i = 0; i < final_group.size(); i++)
    {
        size_t id = final_group.at(i) - _FF_D_OFFSET;
        Inst::FF_D *ptr_FF_D = &_FF_D_arr[id];

        for (size_t kk = 0; kk < ptr_FF_D->faninCone.size(); kk++)
        {
            size_t FFid = ptr_FF_D->faninCone[kk] - _FF_Q_OFFSET;

            if (_FF_D_arr[FFid].grouped == false) // if the slack release to the grouped FF, what should i do ?
            {

                double dist = (_FF_Q_arr.at(FFid).Qx_pos - _FF_D_arr.at(FFid).Dx_pos) * 0.5;
                double D_dia_len = fabs(_FF_D_arr.at(FFid).Dx_pos - _FF_D_arr.at(FFid).D_fanin_pos.first) + fabs(_FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).D_fanin_pos.second) + _FF_D_arr.at(FFid).slack / (_ptr_Parser->_displaceDelay);
                double Dx_pos_r = _FF_D_arr.at(FFid).Dy_pos + _FF_D_arr.at(FFid).Dx_pos + dist; // x'=y+x //D shift right
                double Dy_pos_r = _FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).Dx_pos - dist; // y'=y-x
                coor_w_se Dx_s, Dx_e;
                vector<coor_w_se> dia_x_arr;
                dia_x_arr.reserve(256);
                Dx_s.pos_val = Dx_pos_r - (2 * D_dia_len);
                Dx_s.type = 0;
                Dx_s.or_ind = 0;
                Dx_e.pos_val = Dx_pos_r + (2 * D_dia_len);
                Dx_e.type = 1;
                Dx_e.or_ind = 0;

                dia_x_arr.push_back(Dx_s);
                dia_x_arr.push_back(Dx_e);

                coor_w_se Dy_s, Dy_e;
                vector<coor_w_se> dia_y_arr;
                dia_y_arr.reserve(256);
                Dy_s.pos_val = Dy_pos_r - (2 * D_dia_len);
                Dy_s.type = 0;
                Dy_s.or_ind = 0;
                Dy_e.pos_val = Dy_pos_r + (2 * D_dia_len);
                Dy_e.type = 1;
                Dy_e.or_ind = 0;

                dia_y_arr.push_back(Dy_s);
                dia_y_arr.push_back(Dy_e);

                // Q part, run through all fanout.

                double Qx_pos_r = _FF_Q_arr.at(FFid).Qy_pos + _FF_Q_arr.at(FFid).Qx_pos - dist; // Q shift left
                double Qy_pos_r = _FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Qx_pos + dist;

                for (size_t j = 0; j < _FF_Q_arr.at(FFid).Q_fanout_pos.size(); j++)
                {
                    double Q_dia_len = fabs(_FF_Q_arr.at(FFid).Qx_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second) + _FF_Q_arr.at(FFid).slack / (_ptr_Parser->_displaceDelay);
                    coor_w_se Qx_s, Qx_e;
                    Qx_s.pos_val = Qx_pos_r - (2 * Q_dia_len);
                    Qx_s.type = 0;
                    Qx_s.or_ind = j + 1;
                    Qx_e.pos_val = Qx_pos_r + (2 * Q_dia_len);
                    Qx_e.type = 1;
                    Qx_e.or_ind = j + 1;

                    dia_x_arr.push_back(Qx_s);
                    dia_x_arr.push_back(Qx_e);

                    coor_w_se Qy_s, Qy_e;
                    Qy_s.pos_val = Qy_pos_r - (2 * Q_dia_len);
                    Qy_s.type = 0;
                    Qy_s.or_ind = j + 1;
                    Qy_e.pos_val = Qy_pos_r + (2 * Q_dia_len);
                    Qy_e.type = 1;
                    Qy_e.or_ind = j + 1;

                    dia_y_arr.push_back(Qy_s);
                    dia_y_arr.push_back(Qy_e);
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
                    _FF_D_arr.at(FFid).fea_x_s.FF_id = i;

                    _FF_D_arr.at(FFid).fea_x_e.type = 1;
                    _FF_D_arr.at(FFid).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
                    _FF_D_arr.at(FFid).fea_x_e.FF_id = i;

                    _FF_D_arr.at(FFid).fea_y_s.type = 0;
                    _FF_D_arr.at(FFid).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_y_s.FF_id = i;

                    _FF_D_arr.at(FFid).fea_y_e.type = 1;
                    _FF_D_arr.at(FFid).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
                    _FF_D_arr.at(FFid).fea_y_e.FF_id = i;
                }
            }
        }
    }

    for (size_t i = 0; i < final_group.size(); i++)
    {
        size_t id = final_group.at(i) - _FF_D_OFFSET;
        Inst::FF_Q *ptr_FF_Q = &_FF_Q_arr[id];

        for (size_t kk = 0; kk < ptr_FF_Q->fanoutCone.size(); kk++)
        {
            size_t FFid = ptr_FF_Q->fanoutCone[kk] - _FF_Q_OFFSET;

            if (_FF_D_arr[FFid].grouped == false) // if the slack release to the grouped FF, what should i do ?
            {

                double dist = (_FF_Q_arr.at(FFid).Qx_pos - _FF_D_arr.at(FFid).Dx_pos) * 0.5;
                double D_dia_len = fabs(_FF_D_arr.at(FFid).Dx_pos - _FF_D_arr.at(FFid).D_fanin_pos.first) + fabs(_FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).D_fanin_pos.second) + _FF_D_arr.at(FFid).slack / (_ptr_Parser->_displaceDelay);
                double Dx_pos_r = _FF_D_arr.at(FFid).Dy_pos + _FF_D_arr.at(FFid).Dx_pos + dist; // x'=y+x //D shift right
                double Dy_pos_r = _FF_D_arr.at(FFid).Dy_pos - _FF_D_arr.at(FFid).Dx_pos - dist; // y'=y-x
                coor_w_se Dx_s, Dx_e;
                vector<coor_w_se> dia_x_arr;
                dia_x_arr.reserve(256);
                Dx_s.pos_val = Dx_pos_r - (2 * D_dia_len);
                Dx_s.type = 0;
                Dx_s.or_ind = 0;
                Dx_e.pos_val = Dx_pos_r + (2 * D_dia_len);
                Dx_e.type = 1;
                Dx_e.or_ind = 0;

                dia_x_arr.push_back(Dx_s);
                dia_x_arr.push_back(Dx_e);

                coor_w_se Dy_s, Dy_e;
                vector<coor_w_se> dia_y_arr;
                dia_y_arr.reserve(256);
                Dy_s.pos_val = Dy_pos_r - (2 * D_dia_len);
                Dy_s.type = 0;
                Dy_s.or_ind = 0;
                Dy_e.pos_val = Dy_pos_r + (2 * D_dia_len);
                Dy_e.type = 1;
                Dy_e.or_ind = 0;

                dia_y_arr.push_back(Dy_s);
                dia_y_arr.push_back(Dy_e);

                // Q part, run through all fanout.

                double Qx_pos_r = _FF_Q_arr.at(FFid).Qy_pos + _FF_Q_arr.at(FFid).Qx_pos - dist; // Q shift left
                double Qy_pos_r = _FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Qx_pos + dist;

                for (size_t j = 0; j < _FF_Q_arr.at(FFid).Q_fanout_pos.size(); j++)
                {
                    double Q_dia_len = fabs(_FF_Q_arr.at(FFid).Qx_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).first) + fabs(_FF_Q_arr.at(FFid).Qy_pos - _FF_Q_arr.at(FFid).Q_fanout_pos.at(j).second) + _FF_Q_arr.at(FFid).slack / (_ptr_Parser->_displaceDelay);
                    coor_w_se Qx_s, Qx_e;
                    Qx_s.pos_val = Qx_pos_r - (2 * Q_dia_len);
                    Qx_s.type = 0;
                    Qx_s.or_ind = j + 1;
                    Qx_e.pos_val = Qx_pos_r + (2 * Q_dia_len);
                    Qx_e.type = 1;
                    Qx_e.or_ind = j + 1;

                    dia_x_arr.push_back(Qx_s);
                    dia_x_arr.push_back(Qx_e);

                    coor_w_se Qy_s, Qy_e;
                    Qy_s.pos_val = Qy_pos_r - (2 * Q_dia_len);
                    Qy_s.type = 0;
                    Qy_s.or_ind = j + 1;
                    Qy_e.pos_val = Qy_pos_r + (2 * Q_dia_len);
                    Qy_e.type = 1;
                    Qy_e.or_ind = j + 1;

                    dia_y_arr.push_back(Qy_s);
                    dia_y_arr.push_back(Qy_e);
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
                    _FF_D_arr.at(FFid).fea_x_s.FF_id = i;

                    _FF_D_arr.at(FFid).fea_x_e.type = 1;
                    _FF_D_arr.at(FFid).fea_x_e.pos_val = dia_x_arr.at(target_x).pos_val;
                    _FF_D_arr.at(FFid).fea_x_e.FF_id = i;

                    _FF_D_arr.at(FFid).fea_y_s.type = 0;
                    _FF_D_arr.at(FFid).fea_y_s.pos_val = dia_y_arr.at(target_y - 1).pos_val;
                    _FF_D_arr.at(FFid).fea_y_s.FF_id = i;

                    _FF_D_arr.at(FFid).fea_y_e.type = 1;
                    _FF_D_arr.at(FFid).fea_y_e.pos_val = dia_y_arr.at(target_y).pos_val;
                    _FF_D_arr.at(FFid).fea_y_e.FF_id = i;
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

        if (ess_x_s < ff_x_s)
        {
            if (ess_x_e >= ff_x_s && ess_x_e < ff_x_e)
            {
                x_pos_r.first = ff_x_s;
                x_pos_r.second = ess_x_e;
                valid_x = 1;
                // std::cout << x_pos_r.first << " ," << x_pos_r.second << std::endl;
            }
            else if (ess_x_e >= ff_x_e)
            {
                x_pos_r.first = ff_x_s;
                x_pos_r.second = ff_x_e;
                valid_x = 1;
            }
            else
            {
                x_pos_r.first = ess_x_s;
                x_pos_r.second = ess_x_e;
                valid_x = 0;
            }
        }
        else
        { // ess_x_s >= ff_x_s
            if (ff_x_e >= ess_x_s && ff_x_e < ess_x_e)
            {
                x_pos_r.first = ess_x_s;
                x_pos_r.second = ff_x_e;
                valid_x = 1;
            }
            else if (ess_x_e < ff_x_e)
            {
                x_pos_r.first = ess_x_s;
                x_pos_r.second = ess_x_e;
                valid_x = 1;
            }
            else
            {
                x_pos_r.first = ess_x_s;
                x_pos_r.second = ess_x_e;
                valid_x = 0;
            }
        }

        if (ess_y_s < ff_y_s)
        {
            if (ess_y_e >= ff_y_s && ess_y_e < ff_y_e)
            {
                y_pos_r.first = ff_y_s;
                y_pos_r.second = ess_y_e;
                valid_y = 1;
            }
            else if (ess_y_e >= ff_y_e)
            {
                y_pos_r.first = ff_y_s;
                y_pos_r.second = ff_y_e;
                valid_y = 1;
            }
            else
            {
                y_pos_r.first = ess_y_s;
                y_pos_r.second = ess_y_e;
                valid_y = 0;
            }
        }
        else
        { // ess_x_s >= ff_x_s
            if (ff_y_e >= ess_y_s && ff_y_e < ess_y_e)
            {
                y_pos_r.first = ess_y_s;
                y_pos_r.second = ff_y_e;
                valid_y = 1;
            }
            else if (ess_y_e < ff_y_e)
            {
                y_pos_r.first = ess_y_s;
                y_pos_r.second = ess_y_e;
                valid_y = 1;
            }
            else
            {
                y_pos_r.first = ess_y_s;
                y_pos_r.second = ess_y_e;
                valid_y = 0;
            }
        }

        if (valid_x == 1 && valid_y == 1)
        {
            result_group.push_back(id + _FF_D_OFFSET);
        }
    }
    // std::cout << "findmaximal is completed" << std::endl;

    return result_group;
}

void Solver::Solver::solve_test()
{
    solve_initbuild();

    solve_findfeasible();

    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        vector<size_t> ff_group;
        size_t id = i + _FF_D_OFFSET;
        ff_group.push_back(id);
        prePlace(ff_group, id, _FF_D_arr[i].getPosition());
    }
    legalize();
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
        string oriInstName = _FF_D_arr[i].getOriName().substr(0, _FF_D_arr[i].getName().find('/'));
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
    legalize();

    /*// STA test
    size_t in = _Name_to_ID.at("C9/IN");
    size_t out = _Name_to_ID.at("C12/OUT");
    std::cout << _ptr_STAEngine->getDistance(in, out) << std::endl;*/
}

vector<size_t> Solver::Solver::prePlace(const vector<size_t> &ff_group, size_t essential_ID, pair<double, double> pos)
{
    // return the id get grouped
    vector<size_t> ID_group;
    ID_group.reserve(ff_group.size());
    ID_group.push_back(essential_ID);

    // step1: choose the ff size from ff lib
    int maxSize = 0;
    size_t fftype = 0;
    double area = 0;
    // here we choose the first FF with largest available size in the ff lib
    for (size_t i = 0; i < _ptr_Parser->_flipflopLib.size(); i++)
    {
        if (_ptr_Parser->_flipflopLib[i].Bit > maxSize && _ptr_Parser->_flipflopLib[i].Bit <= int(ff_group.size()))
        {
            maxSize = _ptr_Parser->_flipflopLib[i].Bit;
            fftype = i;
            area = _ptr_Parser->_flipflopLib[i].Hight * _ptr_Parser->_flipflopLib[i].Width;
        }
        if (_ptr_Parser->_flipflopLib[i].Bit == maxSize && area > _ptr_Parser->_flipflopLib[i].Hight * _ptr_Parser->_flipflopLib[i].Width)
        {
            fftype = i;
            area = _ptr_Parser->_flipflopLib[i].Hight * _ptr_Parser->_flipflopLib[i].Width;
        }
    }

    // step2: choose which ff be grouped
    // here we randomly choose other maxsize - 1 ff to be group
    for (size_t i = 0; int(i) < maxSize - 1;)
    {
        if (ff_group[i] != essential_ID)
        {
            ID_group.push_back(ff_group[i]);
            i++;
        }
    }

    // step3: place it and release slack
    pair<double, double> center(0, 0);
    vector<pair<size_t, double>> DpinHeight;
    vector<pair<size_t, double>> ID2HeightMap;
    ID2HeightMap.reserve(maxSize);
    DpinHeight.reserve(maxSize);
    // calculate the center of pin
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
    // initialize ID2HeightMap
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        pair<size_t, double> id2yMap;
        id2yMap.first = ID_group[i];
        id2yMap.second = findPinPosition(ID_group[i]).second;
        ID2HeightMap.push_back(id2yMap);
    }
    // sort the grouped ffs and D pins with their y pos
    std::sort(ID2HeightMap.begin(), ID2HeightMap.end(), [](pair<size_t, double> a, pair<size_t, double> b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    std::sort(DpinHeight.begin(), DpinHeight.end(), [](pair<size_t, double> a, pair<size_t, double> b)
              {
                  return a.second < b.second; // sort in acsending order
              });
    // assign the position to each pin according to their order of height
    // FFName = instanceName of the first pin/
    string FFName = _ID_to_instance[ID2HeightMap[0].first]->getName().substr(0, _ID_to_instance[ID2HeightMap[0].first]->getName().find('/') + 1);
    for (size_t i = 0; int(i) < maxSize; i++)
    {
        size_t FFid = ID2HeightMap[i].first;
        size_t pinID = DpinHeight[i].first;
        Inst::FF_D *ptr_FF_D = &_FF_D_arr[FFid - _FF_D_OFFSET];
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

    // release slack to connected ff
    for (const auto &id : ID_group)
    {
        Inst::FF_D *ptr_FF_D = &_FF_D_arr[id - _FF_D_OFFSET];
        for (size_t i = 0; i < ptr_FF_D->faninCone.size(); i++)
        {
            pair<double, double> OriPos = ptr_FF_D->getOriPosition();
            pair<double, double> NowPos = ptr_FF_D->getPosition();
            pair<double, double> GatePinPos = findPinPosition(ptr_FF_D->outGate2Fanin[i]);
            double Oridis = std::fabs(OriPos.first - GatePinPos.first) + std::fabs(OriPos.second - GatePinPos.second);
            double Nowdis = std::fabs(NowPos.first - GatePinPos.first) + std::fabs(NowPos.second - GatePinPos.second);
            double deltaSlack = (Nowdis - Oridis) / _ptr_Parser->_alpha;
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
            pair<double, double> GatePinPos = findPinPosition(ptr_FF_Q->inGate2Fanout[i]);
            double Oridis = std::fabs(OriPos.first - GatePinPos.first) + std::fabs(OriPos.second - GatePinPos.second);
            double Nowdis = std::fabs(NowPos.first - GatePinPos.first) + std::fabs(NowPos.second - GatePinPos.second);
            double deltaSlack = (Nowdis - Oridis) / _ptr_Parser->_alpha;
            if (_FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].slack > _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].getOriSlack() - deltaSlack)
            {
                _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].slack = _FF_D_arr[ptr_FF_Q->fanoutCone[i] - _FF_D_OFFSET].getOriSlack() - deltaSlack;
            }
        }
    }

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

void Solver::Solver::findFanin(const FF_D_ID &id)
{
    /*
    TODO: run DFS to search the fanin ff for this ff
    */
    assert(_ID_to_instance[id]->getType() == Inst::InstType::INST_FF_D);
    vector<size_t> relateNet = _ID_to_instance[id]->getRelatedNet();
    for (const auto &NetID : relateNet)
    {
        for (const auto &pinID : _NetList[NetID])
        {
            if (pinID == id)
            {
                continue;
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_Q)
            {
                _FF_D_arr[id - _FF_D_OFFSET].faninCone.push_back(pinID);
                _FF_D_arr[id - _FF_D_OFFSET].inGate2Fanin.push_back(_ID_to_instance.size());
                _FF_D_arr[id - _FF_D_OFFSET].outGate2Fanin.push_back(_ID_to_instance.size());
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
            {
                Inst::Gate *ptr_gate = _GID_to_ptrGate_map[pinID - _GATE_OFFSET];
                if (ptr_gate->pinName[pinID - ptr_gate->PIN_OFFSET].find("OUT") != std::string::npos || ptr_gate->pinName[pinID - ptr_gate->PIN_OFFSET].find("out") != std::string::npos)
                {
                    findFaninRecur(id, pinID, pinID);
                }
            }
        }
    }
}

void Solver::Solver::findFaninRecur(const FF_D_ID &ID_ffd, const Gate_ID &gateID, const Gate_ID &baseGID)
{
    Inst::Gate *ptr_gate = _GID_to_ptrGate_map[gateID - _GATE_OFFSET];
    vector<size_t> inputPinID;
    inputPinID.reserve(ptr_gate->pinName.size());
    for (size_t i = 0; i < ptr_gate->pinName.size(); i++)
    {
        if (ptr_gate->pinName[i].find("IN") != std::string::npos || ptr_gate->pinName[i].find("in") != std::string::npos)
        {
            inputPinID.push_back(i + ptr_gate->PIN_OFFSET);
        }
    }
    for (const auto &NetID : ptr_gate->getRelatedNet())
    {
        bool isInputNet = false;
        size_t inpin;
        for (const auto &pinID : _NetList[NetID])
        {
            // check whether _NetList[NetID] is outputPin's Net
            for (const auto &inPinID : inputPinID)
            {
                if (inPinID == pinID)
                {
                    isInputNet = true;
                    inpin = pinID;
                    break;
                }
            }
            // if is outputNet ignore it
            if (isInputNet)
            {
                break;
            }
        }
        if (!isInputNet)
        {
            continue;
        }
        else
        {
            for (const auto &pinID : _NetList[NetID])
            {
                if ((pinID >= ptr_gate->PIN_OFFSET) && (pinID < (ptr_gate->PIN_OFFSET + ptr_gate->pinName.size())))
                {
                    // if pinID is the input of the gate
                    continue;
                }
                if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_Q)
                {
                    _FF_D_arr[ID_ffd - _FF_D_OFFSET].faninCone.push_back(pinID);
                    _FF_D_arr[ID_ffd - _FF_D_OFFSET].inGate2Fanin.push_back(inpin);
                    _FF_D_arr[ID_ffd - _FF_D_OFFSET].outGate2Fanin.push_back(baseGID);
                }
                else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
                {
                    Inst::Gate *ptr = _GID_to_ptrGate_map[pinID - _GATE_OFFSET];
                    if (ptr->pinName[pinID - ptr->PIN_OFFSET].find("OUT") != std::string::npos || ptr->pinName[pinID - ptr->PIN_OFFSET].find("out") != std::string::npos)
                    {
                        findFaninRecur(ID_ffd, pinID, baseGID);
                    }
                }
            }
        }
    }
}

void Solver::Solver::findFanout(const FF_Q_ID &id)
{
    /*
    TODO: run DFS to search the fanout ff for this ff
    */
    assert(_ID_to_instance[id]->getType() == Inst::InstType::INST_FF_Q);
    vector<size_t> relateNet = _ID_to_instance[id]->getRelatedNet();
    for (const auto &NetID : relateNet)
    {
        for (const auto &pinID : _NetList[NetID])
        {
            if (pinID == id)
            {
                continue;
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_D)
            {
                _FF_Q_arr[id - _FF_Q_OFFSET].fanoutCone.push_back(pinID);
                _FF_Q_arr[id - _FF_Q_OFFSET].inGate2Fanout.push_back(_ID_to_instance.size());
                _FF_Q_arr[id - _FF_Q_OFFSET].outGate2Fanout.push_back(_ID_to_instance.size());
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
            {
                findFanoutRecur(id, pinID, pinID);
            }
        }
    }
}

void Solver::Solver::findFanoutRecur(const FF_Q_ID &ID_ffq, const Gate_ID &gateID, const Gate_ID &baseGID)
{
    Inst::Gate *ptr_gate = _GID_to_ptrGate_map[gateID - _GATE_OFFSET];
    vector<size_t> outputPinID;
    outputPinID.reserve(ptr_gate->pinName.size());
    for (size_t i = 0; i < ptr_gate->pinName.size(); i++)
    {
        if (ptr_gate->pinName[i].find("IN") == std::string::npos || ptr_gate->pinName[i].find("in") == std::string::npos)
        {
            outputPinID.push_back(i + ptr_gate->PIN_OFFSET);
        }
    }
    for (const auto &NetID : ptr_gate->getRelatedNet())
    {
        bool isOutputNet = false;
        size_t outpin;
        for (const auto &pinID : _NetList[NetID])
        {
            // check whether _NetList[NetID] is outputPin's Net
            for (const auto &outPinID : outputPinID)
            {
                if (outPinID == pinID)
                {
                    isOutputNet = true;
                    outpin = pinID;
                    break;
                }
            }
            // if is outputNet start to search the net
            if (isOutputNet)
            {
                break;
            }
        }
        if (!isOutputNet)
        {
            continue;
        }
        else
        {
            for (const auto &pinID : _NetList[NetID])
            {
                if ((pinID >= ptr_gate->PIN_OFFSET) && (pinID < (ptr_gate->PIN_OFFSET + ptr_gate->pinName.size())))
                {
                    // if pinID is the pin of the gate
                    continue;
                }
                if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_D)
                {
                    _FF_Q_arr[ID_ffq - _FF_Q_OFFSET].fanoutCone.push_back(pinID);
                    _FF_Q_arr[ID_ffq - _FF_Q_OFFSET].inGate2Fanout.push_back(baseGID);
                    _FF_Q_arr[ID_ffq - _FF_Q_OFFSET].outGate2Fanout.push_back(outpin);
                }
                else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
                {
                    findFaninRecur(ID_ffq, pinID, baseGID);
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