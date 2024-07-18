#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <map>

#include "inst.h"
#include "parser.h"
#include "STAEngine.h"

using std::pair;
using std::string;
using std::vector;

template <typename T, typename Container>
// we expect T is an ID, and container would be vector of pair<ID, whatever>
int binarySearch(const Container &container, const T &key)
{
    int left = 0;
    int right = container.size() - 1;

    while (left <= right)
    {
        int mid = left + (right - left) / 2;

        if (container[mid].first == key)
        {
            return mid; // Key found
        }
        else if (container[mid].first < key)
        {
            left = mid + 1; // Search right half
        }
        else
        {
            right = mid - 1; // Search left half
        }
    }

    return -1; // Key not found
}

double Solver::STAEngine::getDistance(const size_t &inID, const size_t &outID) const
{
    assert(isInPin(inID));
    int inIDPos = binarySearch(_distanceList, inID);
    if (inIDPos != -1)
    {
        int outIDPos = binarySearch(_distanceList[inIDPos].second, outID);
        if (outIDPos == -1)
        {
            return -1;
        }
        return _distanceList[inIDPos].second[outIDPos].second;
    }
    return -1;
}

void Solver::STAEngine::initEngine(const vector<vector<size_t>> &netList, const vector<Inst::Inst *> &IDList)
{
    assert(_ptrSolver != nullptr);
    // step1: initialize the variable
    // initilize the InPinList and OutPinList
    _InPinList.reserve(netList.size());
    _OutPinList.reserve(netList.size());
    int numOfIn = 0;
    for (size_t i = 0; i < netList.size(); i++)
    {
        vector<size_t> initVec;
        initVec.reserve(netList[i].size());
        _InPinList.push_back(initVec);
        _OutPinList.push_back(initVec);
        for (const auto &id : netList[i])
        {
            if (isInPin(id))
            {
                _InPinList[i].push_back(id);
                numOfIn++;
            }
            else if (isOutPin(id))
            {
                _OutPinList[i].push_back(id);
            }
        }
    }
    // initilize the _distanceList
    _distanceList.reserve(numOfIn);
    for (const auto &inNet : _InPinList)
    {
        for (const auto &inPinID : inNet)
        {
            pair<size_t, vector<pair<size_t, double>>> initPair;
            vector<pair<size_t, double>> initVec;
            initVec.reserve(512);
            initPair.first = inPinID;
            initPair.second = initVec;
            _distanceList.push_back(initPair);
        }
    }
    // sort _distanceList
    std::sort(_distanceList.begin(), _distanceList.end(), [](pair<size_t, vector<pair<size_t, double>>> a, pair<size_t, vector<pair<size_t, double>>> b)
              {
                  return a.first < b.first; // sort in acsending order
              });
    // initialize _InPin2PositionMap and _isOutPinCritical
    for (size_t i = 0; i < IDList.size(); i++)
    {
        if (isInPin(i))
        {
            pair<size_t, int> initPair;
            initPair.first = i;
            initPair.second = -1;
            _InPin2PositionMap.insert(initPair);
        }
        else if (isOutPin(i))
        {
            pair<size_t, bool> initPair;
            initPair.first = i;
            initPair.second = false;
            _isOutPinCritical.insert(initPair);
        }
    }
    for (size_t i = 0; i < _distanceList.size(); i++)
    {
        _InPin2PositionMap[_distanceList[i].first] = int(i);
    }
    for (size_t i = 0; i < netList.size(); i++)
    {
        if (_OutPinList[i].size() > 0)
        {
            bool haveFF_D = false;
            for (const auto &pinID : netList[i])
            {
                if (isFF_D(pinID))
                {
                    haveFF_D = true;
                    break;
                }
            }
            if (haveFF_D)
            {
                for (const auto &outID : _OutPinList[i])
                {
                    _isOutPinCritical[outID] = true;
                }
            }
        }
    }

    // initialize _checkList
    _checkList.clear();
    _checkList.reserve(IDList.size());
    for (size_t i = 0; i < IDList.size(); i++)
    {
        _checkList.push_back(false);
    }

    // step2: propagate the distance from the in of gate which connect to Q pin
    for (size_t i = 0; i < netList.size(); i++)
    {
        // std::cout << "start to propagate net " << i << std::endl;
        if (_OutPinList[i].size() == 0)
        {
            // std::cout << "start to propagate net " << i << std::endl;
            for (const auto &inPinID : _InPinList[i])
            {
                list<size_t> l;
                l.push_back(inPinID);
                list<double> emptyList;
                emptyList.clear();
                propagateForward(l, inPinID, 0, emptyList);
            }
        }
    }
}

void Solver::STAEngine::propagateForward(list<size_t> InPinIdList, const size_t &InPinIdNext, double d, list<double> disList)
{
    //  if the InPinIdNext has been explored, directly add its fanout to the Inpin passed
    if (_checkList[InPinIdNext])
    {
        size_t InPinIdBase;
        list<double>::iterator iteDis = disList.begin();
        int posNext = binarySearch(_distanceList, InPinIdNext);
        for (const auto &inID : InPinIdList)
        {
            list<size_t> inPinList = getInPinRelated(inID);
            for (const auto &relatedInID : inPinList)
            {
                InPinIdBase = relatedInID;
                int posBase = binarySearch(_distanceList, InPinIdBase);
                if (posBase == -1) // if the pin doesn't include in any net
                {
                    continue;
                }
                for (const auto &DisPair : _distanceList[posNext].second)
                {
                    int outIndex_base = binarySearch(_distanceList[posBase].second, DisPair.first);
                    if (outIndex_base != -1) // if this path already explored
                    {
                        // and new path is longer than former one, replace it
                        if (_distanceList[posBase].second[outIndex_base].second < d - *iteDis)
                        {
                            _distanceList[posBase].second[outIndex_base].second = d - *iteDis;
                        }
                        continue;
                    }
                    // else add out pin to the base pin
                    pair<size_t, double> disPair;
                    disPair.first = DisPair.first;
                    disPair.second = d;
                    _distanceList[posBase].second.push_back(disPair);
                    // sort it
                    std::sort(_distanceList[posBase].second.begin(), _distanceList[posBase].second.end(), [](pair<size_t, double> a, pair<size_t, double> b)
                              {
                                  return a.first < b.first; // sort in acsending order
                              });
                }
            }
            iteDis++;
        }
        return;
    }

    // if the InPinIdNext hasn't been explored
    if (InPinIdList.front() != InPinIdNext)
    {
        InPinIdList.push_back(InPinIdNext);
    }
    disList.push_back(d);
    list<size_t> outPinList = getOutPinRelated(InPinIdNext);
    int posNext = binarySearch(_distanceList, InPinIdNext);
    assert(posNext != -1);
    for (const auto &outID : outPinList)
    {
        // step1: record this out Pin to the distance list of Base in pin
        // if the out pin is connected to FF, record it
        if (_isOutPinCritical.at(outID))
        {
            size_t InPinIdBase;
            list<double>::iterator iteDis = disList.begin();
            for (const auto &inID : InPinIdList)
            {
                list<size_t> inPinList = getInPinRelated(inID);
                for (const auto &relatedInID : inPinList)
                {
                    InPinIdBase = relatedInID;
                    int posBase = binarySearch(_distanceList, InPinIdBase);
                    if (posBase == -1) // if the pin is not include in any Net
                    {
                        continue;
                    }
                    int outIndex_base = binarySearch(_distanceList[posBase].second, outID);
                    if (outIndex_base != -1) // if this path already explored
                    {
                        // and new path is longer than former one, replace it
                        if (_distanceList[posBase].second[outIndex_base].second < d - *iteDis)
                        {
                            _distanceList[posBase].second[outIndex_base].second = d - *iteDis;
                        }
                        continue;
                    }
                    // add out pin to the base pin
                    if (outIndex_base == -1)
                    {
                        pair<size_t, double> disPair;
                        disPair.first = outID;
                        disPair.second = d;
                        _distanceList[posBase].second.push_back(disPair);
                        // sort it
                        std::sort(_distanceList[posBase].second.begin(), _distanceList[posBase].second.end(), [](pair<size_t, double> a, pair<size_t, double> b)
                                  {
                                      return a.first < b.first; // sort in acsending order
                                  });
                    }
                }
                iteDis++;
            }
        }

        // step2: recurrsively explore the in pin connect to this out pin
        vector<size_t> nextInPin = _InPinList[getRelatedNet(outID)];
        for (const auto &ID_next : nextInPin)
        {
            pair<double, double> outPos = getPinPosition(outID);
            pair<double, double> inPos = getPinPosition(ID_next);
            double dis = std::fabs(outPos.first - inPos.first) + std::fabs(outPos.second - inPos.second);
            propagateForward(InPinIdList, ID_next, d + dis, disList);
        }
    }
    //  step3: mark the input pin of gate to be explored
    list<size_t> inPinList = getInPinRelated(InPinIdNext);
    for (const auto &id : inPinList)
    {
        _checkList[id] = true;
    }
}

bool Solver::STAEngine::isInPin(const size_t &id) const
{
    if (_ptrSolver->_ID_to_instance[id]->getType() == Inst::INST_GATE)
    {
        Inst::Gate *ptrGate = _ptrSolver->_GID_to_ptrGate_map[id - _ptrSolver->_GATE_OFFSET];
        if (ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("IN") != std::string::npos || ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("in") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

bool Solver::STAEngine::isOutPin(const size_t &id) const
{
    if (_ptrSolver->_ID_to_instance[id]->getType() == Inst::INST_GATE)
    {
        Inst::Gate *ptrGate = _ptrSolver->_GID_to_ptrGate_map[id - _ptrSolver->_GATE_OFFSET];
        if (ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("OUT") != std::string::npos || ptrGate->pinName[id - ptrGate->PIN_OFFSET].find("out") != std::string::npos)
        {
            return true;
        }
    }
    return false;
}

bool Solver::STAEngine::isFF_D(const size_t &id) const
{
    if (_ptrSolver->_ID_to_instance[id]->getType() == Inst::INST_FF_D)
    {
        return true;
    }
    return false;
}

list<size_t> Solver::STAEngine::getInPinRelated(const size_t &id) const
{
    // assert(isOutPin(id));
    Inst::Gate *ptr = _ptrSolver->_GID_to_ptrGate_map[id - _ptrSolver->_GATE_OFFSET];
    list<size_t> a;
    for (size_t i = 0; i < ptr->pinName.size(); i++)
    {
        if (ptr->pinName[i].find("IN") != std::string::npos || ptr->pinName[i].find("in") != std::string::npos)
        {
            a.push_back(i + ptr->PIN_OFFSET);
        }
    }
    return a;
}

list<size_t> Solver::STAEngine::getOutPinRelated(const size_t &id) const
{
    // assert(isInPin(id));
    Inst::Gate *ptr = _ptrSolver->_GID_to_ptrGate_map[id - _ptrSolver->_GATE_OFFSET];
    list<size_t> a;
    for (size_t i = 0; i < ptr->pinName.size(); i++)
    {
        if (ptr->pinName[i].find("OUT") != std::string::npos || ptr->pinName[i].find("out") != std::string::npos)
        {
            a.push_back(i + ptr->PIN_OFFSET);
        }
    }
    return a;
}

size_t Solver::STAEngine::getRelatedNet(const size_t &id) const
{
    for (const auto &netID : _ptrSolver->_ID_to_instance[id]->getRelatedNet())
    {
        for (const auto &pinID : _ptrSolver->_NetList[netID])
        {
            if (pinID == id)
            {
                return netID;
            }
        }
    }
    return -1;
}

pair<double, double> Solver::STAEngine::getPinPosition(const size_t &id) const
{
    return _ptrSolver->findPinPosition(id);
}
