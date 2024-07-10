#include <string>
#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <map>
#include <list>

#include "inst.h"
#include "solver.h"
#include "STAEngine.h"

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

    std::cout << "End initialization!!! \n"
              << std::endl;
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
    /*for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        pair<double, double> pos = getFFPosition(&_FF_D_arr[i]);
        string name = _FF_D_arr[i].getName();
        std::cout << name << " " << pos.first << " " << pos.second << std::endl;
    }*/
}

void Solver::Solver::test()
{
    // initSolver();
    string a = "C9/IN";
    size_t in = name2ID(a);
    a = "C12/OUT";
    size_t out = name2ID(a);
    std::cout << _ptr_STAEngine->getDistance(in, out) << std::endl;
}

bool Solver::Solver::mbffCluster() // can add parameter to implement the Window-based sequence generation
{
    /*
    TODO: group the possible Multibit FF and use prePlace() function to place it and release slack
    */
    return true;
}

vector<size_t> Solver::Solver::prePlace(vector<size_t> ff_group, pair<double, double> pos)
{
    // let me know what function do you want it to be
    vector<size_t> a;
    return a;
}

void Solver::Solver::legalize()
{
    std::cout << "Start the legalization!!" << std::endl;
    int _MaxIteration = 30;
    std::cout << "start building data structure for legalization." << std::endl;
    // init the data structure for FF, we assume every FF would have 1 instance and we will assign id for it
    vector<Inst::FF_D *> _id2ffPtr; // please change the ptr to the type you use and init it as following
    _id2ffPtr.reserve(_FF_D_arr.size());
    for (size_t i = 0; i < _FF_D_arr.size(); i++)
    {
        _id2ffPtr.push_back(&_FF_D_arr[i]);
    }

    // init the set of Gate
    vector<Inst::Gate *> _id2gate; // please change the ptr to the type you use and init it as following
    _id2gate.reserve(_Gate_arr.size());
    for (size_t i = 0; i < _Gate_arr.size(); i++)
    {
        _id2gate.push_back(&_Gate_arr[i]);
    }

    std::cout << "start to get placement row information." << std::endl;
    // init the set of placement row
    vector<struct PlacementRow> _Pid2PLR = getPlacementRow();
    /* in these legalization, we assume the placement row are uniform.
    That is, the start x position, siteHeight, siteWidth, totalNumOfSites should be same.
    Also, their start y position should be the lower placeRow's y + siteHeight.
    Base on the assumption, we will first check whether the placement rows meet these conditions
    */
    if (!placementRowIsUniform(_Pid2PLR))
    {
        std::cout << "placement rows aren't uniform, the legaliztion might fail." << std::endl;
    }

    // start the algorithm part
    // define some constant for the algorithm
    size_t numPlaceRow = getPlaceRowNum();
    double siteHeight = getSiteHeight();
    double siteWidth = getSiteWidth();
    double siteNum = getTotalSiteNum();
    double LFx_pos = getLowerLeftX();
    double LFy_pos = getLowerLeftY();

    std::cout << "step 1: construct the matrix for placement grid." << std::endl;
    // step 1: construct the matrix for placement grid
    vector<vector<bool>> _availPosTable; // if the grid point is available, the value would be true
    vector<vector<pair<size_t, int>>> _listWait4Legal;
    _availPosTable.reserve(numPlaceRow);
    for (size_t i = 0; i < numPlaceRow; i++)
    {
        vector<bool> initVec(int(siteNum), true);
        _availPosTable.push_back(initVec);
    }
    for (const auto &ptr_gate : _id2gate)
    {
        // get the lower left and upper right's position for the gate
        pair<double, double> lowerLeft = getGateLF(ptr_gate);
        pair<double, double> upperRight = getGateUR(ptr_gate);
        int LLx = std::floor((lowerLeft.first - LFx_pos) / siteWidth);
        if (LLx < 0)
        {
            LLx = 0;
        }
        int LLy = std::floor((lowerLeft.second - LFy_pos) / siteHeight);
        if (LLy < 0)
        {
            LLy = 0;
        }
        int URx = std::ceil((upperRight.first - LFx_pos) / siteWidth);
        if (URx > siteNum - 1)
        {
            URx = int(siteNum - 1);
        }
        int URy = std::ceil((upperRight.second - LFy_pos) / siteHeight);
        if (URy > int(numPlaceRow - 1))
        {
            URy = int(numPlaceRow - 1);
        }
        // set grid related in _availPosTable to false
        for (int y = LLy; y < URy; y++)
        {
            for (int x = LLx; x < URx; x++)
            {
                _availPosTable[y][x] = false;
            }
        }
    }
    _listWait4Legal.reserve(numPlaceRow);
    for (size_t i = 0; i < numPlaceRow; i++)
    {
        vector<pair<size_t, int>> initVec;
        initVec.reserve(512);
        _listWait4Legal.push_back(initVec);
    }
    std::cout << "step 2: throw all FF to the nearests placement grid and sort it to placement row." << std::endl;
    // step 2: throw all FF to the nearests placement grid and sort it to placement row
    vector<bool> check_list4ffInit(_id2ffPtr.size(), false);
    for (size_t i = 0; i < _id2ffPtr.size(); i++)
    {
        if (check_list4ffInit[i])
        {
            continue;
        }
        pair<double, double> ffPos = getFFPosition(_id2ffPtr[i]);
        int LLx = std::floor((ffPos.first - LFx_pos) / siteWidth);
        int ffWidth = std::ceil(getFFWidth(_id2ffPtr[i]) / siteWidth);
        // boundry checking
        if (LLx < 0)
        {
            LLx = 0;
        }
        else if (LLx >= siteNum - ffWidth)
        {
            LLx = int(siteNum - ffWidth);
        }
        int LLy = std::floor((ffPos.second - LFy_pos) / siteHeight);
        int ffHeight = std::floor(getFFHeight(_id2ffPtr[i]) / siteHeight);
        // boundry checking
        if (LLy < 0)
        {
            LLy = 0;
        }
        else if (LLy >= int(numPlaceRow - ffHeight))
        {
            LLy = int(numPlaceRow - ffHeight);
        }
        ffPos.first = LFx_pos + (double(LLx) * siteWidth);
        ffPos.second = LFy_pos + (double(LLy) * siteHeight);
        setFFPosition(_id2ffPtr[i], ffPos);
        // add to waiting list
        pair<size_t, int> pr(i, LLx);
        _listWait4Legal[LLy].push_back(pr);
        // if group with other, mark on them in the check list
        vector<size_t> groupMemberID = getGroupMem(_id2ffPtr[i]); // this function should return the ID(position in _id2ffPtr) related to this ff
        for (const auto &id : groupMemberID)
        {
            check_list4ffInit[id] = true;
        }
    }
    // sort the wait list in the order of x index
    for (size_t i = 0; i < _listWait4Legal.size(); i++)
    {
        std::sort(_listWait4Legal[i].begin(), _listWait4Legal[i].end(), [](pair<size_t, int> a, pair<size_t, int> b)
                  {
                      return a.second < b.second; // sort in acsending order
                  });
    }

    /*// test!!!!!!!!!!!!!!!!!!!
    for (const auto &aaa : _listWait4Legal[25])
    {
        std::cout << "id in _listWait4Legal[25] " << aaa.first << ", xPos: " << aaa.second << std::endl;
    }

    // end test*/

    std::cout << "step 3: DP from lowest row to the highest." << std::endl;
    // step 3: DP from lowest row to the highest
    int iterCount = 1; // use for iteration counter
    while (iterCount <= _MaxIteration)
    {
        // std::cout << "Start the iteration: " << iterCount << std::endl;
        iterCount++;
        // init some info for algorithm
        int _MaxDisplacement = 100 * iterCount;
        vector<vector<pair<size_t, int>>> _legalist = _listWait4Legal;
        vector<vector<pair<int, int>>> _finalSolution(numPlaceRow); // pair<y index, x index>
        vector<vector<bool>> _ffOccupation;                         // record the position occupied by ff
        _ffOccupation.reserve(numPlaceRow);
        for (size_t i = 0; i < numPlaceRow; i++) // init for _ffOccupation
        {
            vector<bool> initVec(int(siteNum), true);
            _ffOccupation.push_back(initVec);
        }
        for (size_t i = 0; i < numPlaceRow; i++) // init _finalSolution
        {
            vector<pair<int, int>> initVec;
            initVec.reserve(2 * _legalist[i].size());
            _finalSolution.push_back(initVec);
        }
        // start the DP
        bool Solvable = true;
        // std::cout << "Start DP process" << std::endl;
        for (size_t i = 0; i < numPlaceRow; i++)
        {
            if (_legalist[i].empty())
            {
                continue;
            }
            // std::cout << "Start to legalize row " << i << std::endl;
            //   sort the legalist before start to legal it(some object may be added to it from lower row)
            std::sort(_legalist[i].begin(), _legalist[i].end(), [](pair<size_t, int> a, pair<size_t, int> b)
                      {
                          return a.second < b.second; // sort in acsending order
                      });
            vector<vector<bool>> DPtable;                                // DPtable[ i-th object in legalist ][ position ]
            vector<int> heightConstraint(int(siteNum), numPlaceRow - i); // record the available height for x position
            vector<vector<vector<pair<int, int>>>> solutionList;         // pair<int y, int x>
            vector<vector<unsigned int>> totalDisplace;                  // totalDisplace[# of object][last position]
            // init totalDisplace
            totalDisplace.reserve(_legalist[i].size());
            for (size_t k = 0; k < _legalist[i].size(); k++)
            {
                vector<unsigned int> initVec(int(siteNum), -1);
                totalDisplace.push_back(initVec);
            }
            // init solutionList
            solutionList.reserve(_legalist[i].size());
            for (size_t k = 0; k < _legalist[i].size(); k++)
            {
                vector<pair<int, int>> initList;
                vector<vector<pair<int, int>>> initVec(int(siteNum), initList);
                solutionList.push_back(initVec);
            }
            // init DP table
            DPtable.reserve(_legalist[i].size());
            for (size_t k = 0; k < _legalist[i].size(); k++)
            {
                vector<bool> initVec(int(siteNum), false);
                DPtable.push_back(initVec);
            }
            // update the heightConstraint for  this placement row
            for (int k = 0; k < int(siteNum); k++)
            {
                for (size_t j = i; j < numPlaceRow; j++)
                {
                    if (_availPosTable[j][k] && _ffOccupation[j][k])
                    {
                        continue;
                    }
                    else
                    {
                        heightConstraint[k] = j - i;
                        break;
                    }
                }
            }
            // std::cout << "Start to construct table for the first item in legalist " << std::endl;
            //   start to construct table for the first item in legalist
            while (true)
            {
                // construct some local variables for iteration
                int ffHeight = std::ceil(getFFHeight(_id2ffPtr[_legalist[i][0].first]) / siteHeight);
                int ffWidth = std::ceil(getFFWidth(_id2ffPtr[_legalist[i][0].first]) / siteWidth);
                int leftLimit = std::max(0, _legalist[i][0].second - _MaxDisplacement);
                int rightLimit = std::min(_legalist[i][0].second + _MaxDisplacement, int(siteNum) - 1);
                bool keepLoop = true;
                if (_legalist[i].size() < 1)
                {
                    keepLoop = false;
                    break;
                }
                for (int k = 0; k < int(siteNum); k++) // k means you can use 0~k-th grid
                {
                    if (k < ffWidth - 1) // space is not enough for put it in
                    {
                        DPtable[0][k] = false;
                        continue;
                    }
                    if (k - ffWidth + 1 < leftLimit) // the placement is not exceeding left limit
                    {
                        DPtable[0][k] = false;
                        continue;
                    }
                    if ((k > rightLimit) && !DPtable[0][k - 1]) // 0-th element has no solution, throw to higher row
                    {
                        /*
                            can implement putting it to the lower row if it is space
                        */
                        if (i + 2 > numPlaceRow) // do not have upper row
                        {
                            Solvable = false;
                            keepLoop = false;
                            break;
                        }
                        else
                        {
                            _legalist[i + 1].push_back(_legalist[i][0]);
                        }
                        _legalist[i].erase(_legalist[i].begin());
                        break;
                    }
                    bool isSafe = true;
                    for (int j = k; j > std::max(k - ffWidth, -1); j--)
                    {
                        if (!_availPosTable[i][j] || !_ffOccupation[i][j])
                        {
                            isSafe = false;
                            break;
                        }
                        if (heightConstraint[j] < ffHeight)
                        {
                            isSafe = false;
                            break;
                        }
                    }
                    if (isSafe)
                    {
                        DPtable[0][k] = true;
                        unsigned int displace = abs(_legalist[i][0].second + ffWidth - 1 - k);
                        if (totalDisplace[0][k - 1] >= displace)
                        {
                            totalDisplace[0][k] = displace;
                            pair<int, int> position(int(i), k - ffWidth + 1);
                            solutionList[0][k].push_back(position);
                        }
                        else
                        {
                            totalDisplace[0][k] = totalDisplace[0][k - 1];
                            solutionList[0][k] = solutionList[0][k - 1];
                        }
                    }
                    else
                    {
                        if (DPtable[0][k - 1])
                        {
                            DPtable[0][k] = true;
                            totalDisplace[0][k] = totalDisplace[0][k - 1];
                            solutionList[0][k] = solutionList[0][k - 1];
                        }
                    }
                }
                if (DPtable[0][int(siteNum) - 1])
                {
                    keepLoop = false;
                }
                else // item 0 cannot be legalize in row i
                {
                    if (i + 2 > siteNum) // do not have upper row
                    {
                        Solvable = false;
                        keepLoop = false;
                    }
                    else
                    {
                        _legalist[i + 1].push_back(_legalist[i][0]);
                    }
                    _legalist[i].erase(_legalist[i].begin());
                }
                if (!keepLoop)
                {
                    break;
                }
            }
            if (!Solvable)
            {
                break;
            }
            // std::cout << "DP for remaining object" << std::endl;
            //   DP for remaining object
            for (size_t j = 1; j < _legalist[i].size(); j++)
            {
                // optimize the memory usage
                if (j >= 2)
                {
                    solutionList[j - 2].clear();
                    totalDisplace[j - 2].clear();
                    DPtable[j - 2].clear();
                }

                while (true)
                {
                    // construct some local variables for iteration
                    int ffHeight = std::ceil(getFFHeight(_id2ffPtr[_legalist[i][j].first]) / siteHeight);
                    int ffWidth = std::ceil(getFFWidth(_id2ffPtr[_legalist[i][j].first]) / siteWidth);
                    int leftLimit = std::max(0, _legalist[i][j].second - _MaxDisplacement);
                    int rightLimit = std::min(_legalist[i][j].second + _MaxDisplacement, int(siteNum) - 1);
                    bool keepLoop = true;
                    if (j >= _legalist[i].size())
                    {
                        break;
                    }
                    for (int k = 0; k < int(siteNum); k++) // k means you can use 0~k-th grid
                    {
                        solutionList[j][k].reserve(_legalist[i].size());
                        if (k < ffWidth - 1) // space is not enough for put it in
                        {
                            DPtable[j][k] = false;
                            continue;
                        }
                        if (k - ffWidth + 1 < leftLimit) // the placement is not exceeding left limit
                        {
                            DPtable[j][k] = false;
                            continue;
                        }
                        if ((k > rightLimit) && !DPtable[j][k - 1]) // j-th element has no solution, throw to higher row
                        {
                            /*
                                can implement putting it to the lower row if it is space
                            */
                            if (i + 2 > numPlaceRow) // do not have upper row
                            {
                                Solvable = false;
                                keepLoop = false;
                                break;
                            }
                            else
                            {
                                _legalist[i + 1].push_back(_legalist[i][j]);
                            }
                            _legalist[i].erase(_legalist[i].begin() + j);
                            break;
                        }
                        bool isSafe = true; // see whether k ~ k-ffWidth have any ostacle
                        for (int l = k; l > std::max(k - ffWidth, -1); l--)
                        {
                            if (!_availPosTable[i][l] || !_ffOccupation[i][l])
                            {
                                isSafe = false;
                                break;
                            }
                            if (heightConstraint[l] < ffHeight)
                            {
                                isSafe = false;
                                break;
                            }
                        }
                        if (isSafe && DPtable[j - 1][k - ffWidth])
                        {
                            DPtable[j][k] = true;
                            unsigned int displace = abs(_legalist[i][j].second + ffWidth - 1 - k);
                            if (totalDisplace[j][k - 1] >= displace + totalDisplace[j - 1][k - ffWidth])
                            {
                                totalDisplace[j][k] = displace + totalDisplace[j - 1][k - ffWidth];
                                pair<int, int> position(int(i), k - ffWidth + 1);
                                solutionList[j][k] = solutionList[j - 1][k - ffWidth];
                                solutionList[j][k].push_back(position);
                            }
                            else
                            {
                                totalDisplace[j][k] = totalDisplace[j][k - 1];
                                solutionList[j][k] = solutionList[j][k - 1];
                            }
                        }
                        else
                        {
                            if (DPtable[j][k - 1])
                            {
                                DPtable[j][k] = true;
                                totalDisplace[j][k] = totalDisplace[j][k - 1];
                                solutionList[j][k] = solutionList[j][k - 1];
                            }
                        }
                    }
                    if (DPtable[j][int(siteNum) - 1])
                    {
                        keepLoop = false;
                    }
                    else // item j cannot be legalize in row i
                    {
                        if (i + 2 > numPlaceRow) // do not have upper row
                        {
                            Solvable = false;
                            keepLoop = false;
                            break;
                        }
                        else
                        {
                            _legalist[i + 1].push_back(_legalist[i][j]);
                        }
                        _legalist[i].erase(_legalist[i].begin() + j);
                    }
                    if (!keepLoop)
                    {
                        break;
                    }
                }
                // if unsolvable, there's no need to DP for remaining item
                if (!Solvable)
                {
                    break;
                }
            }
            if (!Solvable)
            {
                break;
            }
            // std::cout << "successfully solve row " << i << ", record the solution" << std::endl;
            //   as the assign for one row complete,
            //   we start to record the solution and start the legalize for next row
            int counter = 0;

            // test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            /*std::cout << "solutionList.size = " << solutionList.size() << std::endl;
            std::cout << "_legalist[i].size = " << _legalist[i].size() << std::endl;
            std::cout << "solutionList[_legalist[i].size() - 1].size = " << solutionList[_legalist[i].size() - 1].size() << std::endl;
            std::cout << "int(siteNum) = " << int(siteNum) << std::endl;
            int ccc = 0;
            for (const auto &a : solutionList[_legalist[i].size() - 1][int(siteNum) - 1])
            {
                std::cout << "AAA" << std::endl;
                if (a.first != int(i))
                {
                    std::cout << "i = " << i << std::endl;
                    std::cout << "id = " << _legalist[i][ccc].first << std::endl;
                    std::cout << "FFName = " << _id2ffPtr[_legalist[i][ccc].first]->getName() << std::endl;
                    std::cout << "(*it).first = " << a.first << std::endl;
                    std::cout << "(*it).second = " << a.second << std::endl;
                    ccc++;
                }
            }
            // end test*/
            if (_legalist[i].empty())
            {
                continue;
            }
            for (auto it = solutionList[_legalist[i].size() - 1][int(siteNum) - 1].begin();
                 it != solutionList[_legalist[i].size() - 1][int(siteNum) - 1].end(); it++)
            {
                // std::cout << "counter = " << counter << std::endl;
                int ffHeight = std::ceil(getFFHeight(_id2ffPtr[_legalist[i][counter].first]) / siteHeight);
                int ffWidth = std::ceil(getFFWidth(_id2ffPtr[_legalist[i][counter].first]) / siteWidth);
                // std::cout << "FFName = " << _id2ffPtr[_legalist[i][counter].first]->getName() << std::endl;
                // std::cout << "ffHeight = " << ffHeight << std::endl;
                // std::cout << "ffWidth = " << ffWidth << std::endl;
                // std::cout << "(*it).first = " << (*it).first << std::endl;
                // std::cout << "(*it).second = " << (*it).second << std::endl;
                _finalSolution[i].push_back(*it);
                // record the space occupied by the placed ff
                for (int y = 0; y < ffHeight; y++)
                {
                    for (int x = 0; x < ffWidth; x++)
                    {
                        _ffOccupation[(*it).first + y][(*it).second + x] = false;
                    }
                }
                // std::cout << "counter = " << counter << " end" << std::endl;
                counter++;
            }
        }

        std::cout << "step 4: assign the result." << std::endl;
        // step 4: assign the result. If can't be solved, relax the constraint and goto step 3.
        if (!Solvable)
        {
            continue;
        }
        else
        {
            for (size_t i = 0; i < numPlaceRow; i++)
            {
                for (size_t j = 0; j < _legalist[i].size(); j++)
                {
                    pair<double, double> position;
                    position.second = LFy_pos + double(_finalSolution[i][j].first) * siteHeight;
                    position.first = LFx_pos + double(_finalSolution[i][j].second) * siteWidth;
                    setFFPosition(_id2ffPtr[_legalist[i][j].first], position);
                }
            }
            std::cout << "legalization successfully. \n"
                      << std::endl;
            return;
        }
    }
    std::cout << "legalization fail !!!!!!!!! \n"
              << std::endl;
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