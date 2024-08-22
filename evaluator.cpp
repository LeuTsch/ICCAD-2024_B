#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <assert.h>

#include "evaluator.h"
#include "STAEngine.h"

double Solver::evaluator::evaluate()
{
    assert(_ptrSolver != nullptr);
    calculateTNS();
    calculatePower();
    calculateUtil();
    calculateArea();
    double output = (_ptrSolver->_ptr_Parser->_alpha * _totalTNS) + (_ptrSolver->_ptr_Parser->_beta * _totalPower) + (_ptrSolver->_ptr_Parser->_gamma * _totalArea) + (_ptrSolver->_ptr_Parser->_lambda * _totalUtil);
    return output;
}

void Solver::evaluator::evaluate(const string &fileName)
{
    assert(_ptrSolver != nullptr);

    calculateTNS();
    calculatePower();
    calculateUtil();
    calculateArea();

    std::fstream outputFile;
    outputFile.open(fileName, std::fstream::out);
    outputFile << std::setprecision(15);
    double output = (_ptrSolver->_ptr_Parser->_alpha * _totalTNS) + (_ptrSolver->_ptr_Parser->_beta * _totalPower) + (_ptrSolver->_ptr_Parser->_gamma * _totalArea) + (_ptrSolver->_ptr_Parser->_lambda * _totalUtil);
    outputFile << "cost metrix: " << output << "\n";
    outputFile << "TNS: " << _totalTNS << "\n";
    outputFile << "Power: " << _totalPower << "\n";
    outputFile << "Area: " << _totalArea << "\n";
    outputFile << "Util: " << _totalUtil << "\n";
    outputFile.close();
}

double Solver::evaluator::evaluateTNS()
{
    assert(_ptrSolver != nullptr);
    calculateTNS();
    return _totalTNS;
}

void Solver::evaluator::calculatePower()
{
    _totalPower = 0;
    // step1: initialize the _FF_DList
    _FF_DList.clear();
    _FF_DList.reserve(_ptrSolver->_FF_D_arr.size());
    for (size_t i = 0; i < _ptrSolver->_FF_D_arr.size(); i++)
    {
        _FF_DList.push_back(&_ptrSolver->_FF_D_arr[i]);
    }
    vector<bool> isChecked(_ptrSolver->_FF_D_arr.size(), false);

    // step2: calculate the power for all FF
    for (size_t i = 0; i < _FF_DList.size(); i++)
    {
        if (isChecked.at(i))
        {
            continue;
        }

        // add the power of FF to _totalPower
        _totalPower += _ptrSolver->_ptr_Parser->_flipflopLib.at(_FF_DList.at(i)->FF_type).Power;
        // mark the group member
        for (const auto &FFID : _FF_DList.at(i)->grouped_member)
        {
            size_t id = FFID - _ptrSolver->_FF_D_OFFSET;
            isChecked.at(id) = true;
        }
    }
}

void Solver::evaluator::calculateArea()
{
    _totalArea = 0;
    // step1: initialize the _FF_DList
    _FF_DList.clear();
    _FF_DList.reserve(_ptrSolver->_FF_D_arr.size());
    for (size_t i = 0; i < _ptrSolver->_FF_D_arr.size(); i++)
    {
        _FF_DList.push_back(&_ptrSolver->_FF_D_arr[i]);
    }
    vector<bool> isChecked(_ptrSolver->_FF_D_arr.size(), false);

    // step2: calculate the area for all FF
    for (size_t i = 0; i < _FF_DList.size(); i++)
    {
        if (isChecked.at(i))
        {
            continue;
        }

        // add the area of FF to _totalArea
        _totalArea += (_ptrSolver->_ptr_Parser->_flipflopLib.at(_FF_DList.at(i)->FF_type).Width) * (_ptrSolver->_ptr_Parser->_flipflopLib.at(_FF_DList.at(i)->FF_type).Hight);
        // mark the group member
        for (const auto &FFID : _FF_DList.at(i)->grouped_member)
        {
            size_t id = FFID - _ptrSolver->_FF_D_OFFSET;
            isChecked.at(id) = true;
        }
    }
}

void Solver::evaluator::calculateUtil()
{
    _totalUtil = 0;
    // step1: initialize all parameter
    _DieBasePoint.first = double(_ptrSolver->_ptr_Parser->_dieLx);
    _DieBasePoint.second = double(_ptrSolver->_ptr_Parser->_dieLy);
    _binWidth = _ptrSolver->_ptr_Parser->_binWidth;
    _binHeight = _ptrSolver->_ptr_Parser->_binHeight;
    _binMaxUtil = _ptrSolver->_ptr_Parser->_binMaxUtil;
    _dieWidth = double(_ptrSolver->_ptr_Parser->_dieRx - _ptrSolver->_ptr_Parser->_dieLx);
    _dieHeight = double(_ptrSolver->_ptr_Parser->_dieRy - _ptrSolver->_ptr_Parser->_dieLy);
    bool isWidthDivisible = (std::floor(_dieWidth / _binWidth) == std::ceil(_dieWidth / _binWidth));
    bool isHeightDivisible = (std::floor(_dieHeight / _binHeight) == std::ceil(_dieHeight / _binHeight));
    int numberOfHeight = int(std::ceil(_dieHeight / _binHeight));
    int numberOfWidth = int(std::ceil(_dieWidth / _binWidth));
    _FF_DList.clear();
    _FF_DList.reserve(_ptrSolver->_FF_D_arr.size());
    for (size_t i = 0; i < _ptrSolver->_FF_D_arr.size(); i++)
    {
        _FF_DList.push_back(&_ptrSolver->_FF_D_arr[i]);
    }
    _GateList.clear();
    _GateList.reserve(_ptrSolver->_Gate_arr.size());
    for (size_t i = 0; i < _ptrSolver->_Gate_arr.size(); i++)
    {
        _GateList.push_back(&_ptrSolver->_Gate_arr[i]);
    }
    for (int x = 0; x < numberOfWidth; x++)
    {
        vector<double> initVec(numberOfHeight, 0);
        _areaInBin.push_back(initVec);
        _areaOfBin.push_back(initVec);
    }
    for (int x = 0; x < numberOfWidth; x++)
    {
        double Width = _binWidth;
        if ((!isWidthDivisible) && (x == numberOfWidth - 1))
        {
            Width = _dieWidth - (std::floor(_dieWidth / _binWidth) * _binWidth);
        }

        for (int y = 0; y < numberOfHeight; y++)
        {
            double Height = _binHeight;
            if ((!isHeightDivisible) && (y == numberOfHeight - 1))
            {
                Height = _dieHeight - (std::floor(_dieHeight / _binHeight) * _binHeight);
            }
            _areaOfBin.at(x).at(y) = Width * Height;
        }
    }

    // step2: calculate the area in each bin for all gate
    for (const auto &ptr_gate : _GateList)
    {
        // get the lower left and upper right's position for the gate
        pair<double, double> lowerLeft = _ptrSolver->getGateLF(ptr_gate);
        pair<double, double> upperRight = _ptrSolver->getGateUR(ptr_gate);
        int LLx = std::floor((lowerLeft.first - _DieBasePoint.first) / _binWidth);
        if (LLx < 0)
        {
            LLx = 0;
        }
        int LLy = std::floor((lowerLeft.second - _DieBasePoint.second) / _binHeight);
        if (LLy < 0)
        {
            LLy = 0;
        }
        int URx = std::ceil((upperRight.first - _DieBasePoint.first) / _binWidth);
        if (URx > numberOfWidth)
        {
            URx = int(numberOfWidth);
        }
        int URy = std::ceil((upperRight.second - _DieBasePoint.second) / _binHeight);
        if (URy > int(numberOfHeight))
        {
            URy = int(numberOfHeight);
        }

        for (int x = LLx; x < URx; x++)
        {
            double Width = std::min((double(x + 1) * _binWidth + _DieBasePoint.first), upperRight.first) - std::max((double(x) * _binWidth + _DieBasePoint.first), lowerLeft.first);
            for (int y = LLy; y < URy; y++)
            {
                double Height = std::min((double(y + 1) * _binHeight + _DieBasePoint.second), upperRight.second) - std::max((double(y) * _binHeight + _DieBasePoint.second), lowerLeft.second);
                _areaInBin.at(x).at(y) += Height * Width;
            }
        }
    }

    // step3: calculate the area in each bin for all FF
    vector<bool> isChecked(_FF_DList.size(), false);
    for (size_t i = 0; i < _FF_DList.size(); i++)
    {
        if (isChecked.at(i))
        {
            continue;
        }

        pair<double, double> lowerLeft = getFFLF(_FF_DList.at(i));
        pair<double, double> upperRight = getFFUR(_FF_DList.at(i));
        int LLx = std::floor((lowerLeft.first - _DieBasePoint.first) / _binWidth);
        if (LLx < 0)
        {
            LLx = 0;
        }
        int LLy = std::floor((lowerLeft.second - _DieBasePoint.second) / _binHeight);
        if (LLy < 0)
        {
            LLy = 0;
        }
        int URx = std::ceil((upperRight.first - _DieBasePoint.first) / _binWidth);
        if (URx > numberOfWidth)
        {
            URx = int(numberOfWidth);
        }
        int URy = std::ceil((upperRight.second - _DieBasePoint.second) / _binHeight);
        if (URy > int(numberOfHeight))
        {
            URy = int(numberOfHeight);
        }
        for (int x = LLx; x < URx; x++)
        {
            double Width = std::min((double(x + 1) * _binWidth + _DieBasePoint.first), upperRight.first) - std::max((double(x) * _binWidth + _DieBasePoint.first), lowerLeft.first);
            for (int y = LLy; y < URy; y++)
            {
                double Height = std::min((double(y + 1) * _binHeight + _DieBasePoint.second), upperRight.second) - std::max((double(y) * _binHeight + _DieBasePoint.second), lowerLeft.second);
                _areaInBin.at(x).at(y) += Height * Width;
            }
        }
        // mark the group member
        for (const auto &FFID : _FF_DList.at(i)->grouped_member)
        {
            size_t id = FFID - _ptrSolver->_FF_D_OFFSET;
            isChecked.at(id) = true;
        }
    }

    // check how many bin exceed the MaxUtil
    for (int x = 0; x < numberOfWidth; x++)
    {
        for (int y = 0; y < numberOfHeight; y++)
        {
            if (_areaInBin.at(x).at(y) / _areaOfBin.at(x).at(y) > _binMaxUtil)
            {
                _totalUtil++;
            }
        }
    }
}

void Solver::evaluator::calculateTNS()
{
    _totalTNS = 0;
    // step1: initialize the _FF_DList
    _FF_DList.clear();
    _FF_DList.reserve(_ptrSolver->_FF_D_arr.size());
    for (size_t i = 0; i < _ptrSolver->_FF_D_arr.size(); i++)
    {
        _FF_DList.push_back(&_ptrSolver->_FF_D_arr[i]);
    }

    // step2: calculate slack for all FF
    for (size_t i = 0; i < _FF_DList.size(); i++)
    {
        double slack = calSlack4FF(_FF_DList.at(i));
        // std::cout << _FF_DList.at(i)->getName() << ": " << _FF_DList.at(i)->maxSlack << " " << slack << std::endl;
        if (slack < 0)
        {
            _totalTNS += -slack;
        }
        // std::cout << _totalTNS << std::endl;
    }
}

pair<double, double> Solver::evaluator::getFFLF(Inst::FF_D *ptr)
{
    return _ptrSolver->getFFPosition(ptr);
}

pair<double, double> Solver::evaluator::getFFUR(Inst::FF_D *ptr)
{
    pair<double, double> pos = _ptrSolver->getFFPosition(ptr);
    pos.first += _ptrSolver->getFFWidth(ptr);
    pos.second += _ptrSolver->getFFHeight(ptr);
    return pos;
}

double Solver::evaluator::calSlack4FF(Inst::FF_D *ptr)
{
    size_t numOfConnectedFF = ptr->faninCone.size();
    // case1: _FF_D is directly connect to PI
    if (numOfConnectedFF == 0) // for the case _FF_D is directly connect to PI
    {
        double output = 0;
        double distance = 0;
        for (const auto &pinID : _ptrSolver->_NetList[ptr->getRelatedNet().at(0)])
        {
            if (_ptrSolver->_ID_to_instance[pinID]->getType() == Inst::INST_PIO)
            {
                pair<double, double> PIPos = _ptrSolver->findPinPosition(pinID);
                pair<double, double> FFPos = ptr->getPosition();
                distance += (std::fabs(PIPos.first - FFPos.first) + std::fabs(PIPos.second - FFPos.second));
                break;
            }
        }
        output = ptr->maxSlack - (distance * _ptrSolver->_ptr_Parser->_displaceDelay);
        return output;
    }

    // case2: _FF_D is not directly connect to PI
    double slackRemain = ptr->maxSlack;
    for (size_t i = 0; i < numOfConnectedFF; i++)
    {
        double distance = 0;
        if (ptr->inGate2Fanin[i] != _ptrSolver->_ID_to_instance.size()) // if the former FF is not directly connect to FF_D
        {
            distance += _ptrSolver->_ptr_STAEngine->getDistance(ptr->inGate2Fanin[i], ptr->outGate2Fanin[i]);
            // calculate the distance between outPin of gate and FF_D
            pair<double, double> gatePos = _ptrSolver->findPinPosition(ptr->outGate2Fanin[i]);
            pair<double, double> FFPos = ptr->getPosition();
            distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
            // calculate the distance between inPin of gate and FF_D
            gatePos = _ptrSolver->findPinPosition(ptr->inGate2Fanin[i]);
            FFPos = _ptrSolver->findPinPosition(ptr->faninCone[i]);
            distance += (std::fabs(gatePos.first - FFPos.first) + std::fabs(gatePos.second - FFPos.second));
            // add the effect of Qpin delay
            distance += (_ptrSolver->_ptr_Parser->_flipflopLib[_ptrSolver->_FF_Q_arr[ptr->faninCone[i] - _ptrSolver->_FF_Q_OFFSET].FF_type].PinDelay) / _ptrSolver->_ptr_Parser->_displaceDelay;

            double s;
            // check whether the connected FF has been grouped
            // if is grouped give all slack to FF_D
            s = ptr->maxSlack - (distance * _ptrSolver->_ptr_Parser->_displaceDelay);

            if (s < slackRemain)
            {
                slackRemain = s;
            }
        }
        else // the former FF is directly connect to FF_D
        {
            pair<double, double> FFPos = ptr->getPosition();
            pair<double, double> FFBeforePos = _ptrSolver->findPinPosition(ptr->faninCone[i]);
            distance += (std::fabs(FFBeforePos.first - FFPos.first) + std::fabs(FFBeforePos.second - FFPos.second));
            // add the effect of Qpin delay
            distance += (_ptrSolver->_ptr_Parser->_flipflopLib[_ptrSolver->_FF_Q_arr[ptr->faninCone[i] - _ptrSolver->_FF_Q_OFFSET].FF_type].PinDelay) / _ptrSolver->_ptr_Parser->_displaceDelay;

            double s;
            // check whether the connected FF has been grouped

            // if is grouped give all slack to FF_D
            s = ptr->maxSlack - (distance * _ptrSolver->_ptr_Parser->_displaceDelay);

            if (s < slackRemain)
            {
                slackRemain = s;
            }
        }
    }

    double output = slackRemain;

    return output;
}
