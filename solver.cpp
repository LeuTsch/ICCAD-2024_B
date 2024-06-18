#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <assert.h>

#include "inst.h"
#include "solver.h"

using std::pair;
using std::string;
using std::vector;

void Solver::Solver::initSolver()
{
    /*
    TODO: use the data in Parser to construct the data structure
    */
    assert(_ptr_Parser != nullptr);

    // step1: realize the instance for PIO, ff, and Gate,and assign ID to them
    // realize PIO
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
    _Name_to_ID.clear();
    _GID_to_ptrGate_map.clear();
    // assign Id for PIO
    for (size_t i = 0; i < _PIO_arr.size(); i++)
    {
        Inst::Inst *ptr = &_PIO_arr[i];
        _Name_to_ID.push_back(ptr->getName());
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
            _Name_to_ID.push_back(ptr->getName() + "/" + _Gate_arr[i].pinName[j]);
            _ID_to_instance.push_back(ptr);
            _GID_to_ptrGate_map.push_back(ptr_gate);
        }
    }
    _FF_D_OFFSET = _ID_to_instance.size();
    // assign Id for FF_D
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        Inst::Inst *ptr = &_FF_D_arr[i];
        _Name_to_ID.push_back(_FF_D_arr[i].getOriName());
        _ID_to_instance.push_back(ptr);
    }
    _FF_Q_OFFSET = _ID_to_instance.size();
    // assign Id for FF_Q
    for (size_t i = 0; i < _FF_Q_arr.size(); i++)
    {
        Inst::Inst *ptr = &_FF_Q_arr[i];
        _Name_to_ID.push_back(_FF_Q_arr[i].getOriName());
        _ID_to_instance.push_back(ptr);
    }

    // step2: establish Net
    for (const auto &netPinList : _ptr_Parser->_netPin)
    {
        bool isCLKNet = false;
        for (const auto &PinName : netPinList)
        {
            if (PinName.find("CLK") != std::string::npos)
            {
                isCLKNet = true;
                break;
            }
        }
        if (isCLKNet)
        {
            for (const auto &PinName : netPinList)
            {
                if (PinName.find("OUT") != std::string::npos)
                {
                    _ClkList.push_back(PinName);
                    break;
                }
            }
            for (const auto &PinName : netPinList)
            {
                if (PinName.find("CLK") != std::string::npos)
                {
                    string insName = PinName.substr(0, PinName.find("/") + 1);
                    for (size_t i = 0; i < _FF_D_arr.size(); i++)
                    {
                        if (_FF_D_arr[i].getOriName().find(insName) == 0)
                        {
                            _FF_D_arr[i].setClk(_ClkList.size() - 1);
                            _FF_Q_arr[i].setClk(_ClkList.size() - 1);
                        }
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
                for (size_t i = 0; i < _Name_to_ID.size(); i++)
                {
                    if (PinName == _Name_to_ID[i])
                    {
                        _ID_to_instance[i]->addRelatedNet(netID);
                        net.push_back(i);
                        break;
                    }
                }
            }
            _NetList.push_back(net);
        }
    }

    // step3: establish the fanin and fanout cone for the FF
    for (size_t ffd_id = _FF_D_OFFSET; ffd_id < _FF_D_arr.size(); ffd_id++)
    {
        findFanin(ffd_id);
    }
    for (size_t ffq_id = _FF_Q_OFFSET; ffq_id < _FF_Q_arr.size(); ffq_id++)
    {
        findFanout(ffq_id);
    }

    // step4: update the slack for the pin
    for (const auto &N2Slack : _ptr_Parser->_timeSlack)
    {
        for (size_t i = 0; i < _FF_D_arr.size(); i++)
        {
            if (_FF_D_arr[i].getOriName() == N2Slack.first)
            {
                _FF_D_arr[i].setOriSlack(N2Slack.second);
                break;
            }
        }
    }

    // step5: init the placement row
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
        _Gate_in_Bin[x_index][y_index].push_back(id);
    }
    // sort Placement Row into bin
    for (size_t i = 0; i < _PlaceRow.size(); i++)
    {
        int Lx_index = (int(_PlaceRow[i].x) - (int(_PlaceRow[i].x) % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int Ly_index = (int(_PlaceRow[i].y) - (int(_PlaceRow[i].y) % _ptr_Parser->_binHeight)) / _ptr_Parser->_binHeight;
        int Rx = int(_PlaceRow[i].x + _PlaceRow[i].siteWidth * _PlaceRow[i].totalNumOfSites);
        int Ry = int(_PlaceRow[i].y + _PlaceRow[i].siteHight);
        int Rx_index = (Rx - (Rx % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        int Ry_index = (Ry - (Ry % _ptr_Parser->_binWidth)) / _ptr_Parser->_binWidth;
        for (int x = Lx_index; x <= Rx_index; x++)
        {
            for (int y = Ly_index; y <= Ry_index; y++)
            {
                _PlaceRow_in_Bin[x][y].push_back(i);
            }
        }
    }
}

void Solver::Solver::solve()
{
    /*
    TODO: implement the MBFF clustering algorithm
    */
}

void Solver::Solver::printOutput()
{
    /*
    TODO: output
    */
}

bool Solver::Solver::mbffCluster() // can add parameter to implement the Window-based sequence generation
{
    /*
    TODO: group the possible Multibit FF and use prePlace() function to place it and release slack
    */
    return true;
}

bool Solver::Solver::prePlace(vector<size_t> ff_group, size_t placementRowID, pair<double, double> pos)
{
    // let me know what function do you want it to be
    return true;
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
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
            {
                findFaninRecur(id, pinID);
            }
        }
    }
}

void Solver::Solver::findFaninRecur(const FF_D_ID &ID_ffd, const Gate_ID &gateID)
{
    Inst::Gate *ptr_gate = _GID_to_ptrGate_map[gateID - _GATE_OFFSET];
    vector<size_t> outputPinID;
    outputPinID.reserve(ptr_gate->pinName.size());
    for (size_t i = 0; i < ptr_gate->pinName.size(); i++)
    {
        if (ptr_gate->pinName[i].find("IN") == std::string::npos)
        {
            outputPinID.push_back(i + ptr_gate->PIN_OFFSET);
        }
    }
    for (const auto &NetID : ptr_gate->getRelatedNet())
    {
        bool isOutputNet = false;
        for (const auto &pinID : _NetList[NetID])
        {
            // check whether _NetList[NetID] is outputPin's Net
            for (const auto &outPinID : outputPinID)
            {
                if (outPinID == pinID)
                {
                    isOutputNet = true;
                    break;
                }
            }
            // if is outputNet ignore it
            if (isOutputNet)
            {
                break;
            }
        }
        if (isOutputNet)
        {
            continue;
        }
        else
        {
            for (const auto &pinID : _NetList[NetID])
            {
                if ((pinID >= _GATE_OFFSET) && (pinID < (_GATE_OFFSET + ptr_gate->pinName.size())))
                {
                    // if pinID is the input of the gate
                    continue;
                }
                if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_Q)
                {
                    _FF_D_arr[ID_ffd - _FF_D_OFFSET].faninCone.push_back(pinID);
                }
                else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
                {
                    findFaninRecur(ID_ffd, pinID);
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
            }
            else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
            {
                findFanoutRecur(id, pinID);
            }
        }
    }
}

void Solver::Solver::findFanoutRecur(const FF_Q_ID &ID_ffq, const Gate_ID &gateID)
{
    Inst::Gate *ptr_gate = _GID_to_ptrGate_map[gateID - _GATE_OFFSET];
    vector<size_t> outputPinID;
    outputPinID.reserve(ptr_gate->pinName.size());
    for (size_t i = 0; i < ptr_gate->pinName.size(); i++)
    {
        if (ptr_gate->pinName[i].find("IN") == std::string::npos)
        {
            outputPinID.push_back(i + ptr_gate->PIN_OFFSET);
        }
    }
    for (const auto &NetID : ptr_gate->getRelatedNet())
    {
        bool isOutputNet = false;
        for (const auto &pinID : _NetList[NetID])
        {
            // check whether _NetList[NetID] is outputPin's Net
            for (const auto &outPinID : outputPinID)
            {
                if (outPinID == pinID)
                {
                    isOutputNet = true;
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
                if ((pinID >= _GATE_OFFSET) && (pinID < (_GATE_OFFSET + ptr_gate->pinName.size())))
                {
                    // if pinID is the pin of the gate
                    continue;
                }
                if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_FF_D)
                {
                    _FF_Q_arr[ID_ffq - _FF_Q_OFFSET].fanoutCone.push_back(pinID);
                }
                else if (_ID_to_instance[pinID]->getType() == Inst::InstType::INST_GATE)
                {
                    findFaninRecur(ID_ffq, pinID);
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
        pair<double, double> pos = ptr_gate->getPosition();
        pos.first += ptr_gate->pinPosition[id - ptr_gate->PIN_OFFSET].first;
        pos.second += ptr_gate->pinPosition[id - ptr_gate->PIN_OFFSET].second;
        return pos;
    }
    else
    {
        return _ID_to_instance[id]->getPosition();
    }
}

size_t Solver::Solver::name2ID(string &name) const
{
    for (size_t i = 0; i < _Name_to_ID.size(); i++)
    {
        if (name == _Name_to_ID[i])
        {
            return i;
        }
    }
    return -1;
}
