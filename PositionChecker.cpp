#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <assert.h>

#include "inst.h"
#include "parser.h"
#include "PositionChecker.h"

vector<pair<double, double>> Solver::PositionChecker::findAvaiablePosition(size_t fflibID, pair<double, double> position, double radius)
{
    initPositionChecker();
    vector<pair<double, double>> output;
    output.reserve(512);

    // step1: find the position is in which legal Region
    bool isInCertainRegion = false;
    int regionIndex = -1;
    // check if the ff is in any of placeRegion
    for (size_t i = 0; i < _PlaceRegionSet.size(); i++)
    {
        if (position.first >= _PlaceRegionSet.at(i).xStart && position.first < _PlaceRegionSet.at(i).xEnd)
        {
            if (position.second >= _PlaceRegionSet.at(i).yStart && position.second < _PlaceRegionSet.at(i).yEnd)
            {
                isInCertainRegion = true;
                regionIndex = i;
                break;
            }
        }
    }
    // if ff is not in any of placeRegion
    if (!isInCertainRegion)
    {
        vector<pair<double, size_t>> dis2Region;
        dis2Region.reserve(_PlaceRegionSet.size());
        // (1): calculate the distance to all Region
        for (size_t i = 0; i < _PlaceRegionSet.size(); i++)
        {
            double distance = 0;
            // calculate the y distance
            if (position.first < _PlaceRegionSet.at(i).xStart)
            {
                distance += std::fabs(position.first - _PlaceRegionSet.at(i).xStart);
            }
            else if (position.first > _PlaceRegionSet.at(i).xEnd)
            {
                distance += std::fabs(position.first - _PlaceRegionSet.at(i).xEnd);
            }
            // calculate the y distance
            if (position.second < _PlaceRegionSet.at(i).yStart)
            {
                distance += std::fabs(position.second - _PlaceRegionSet.at(i).yStart);
            }
            else if (position.second > _PlaceRegionSet.at(i).yEnd)
            {
                distance += std::fabs(position.second - _PlaceRegionSet.at(i).yEnd);
            }
            pair<double, size_t> initPair(distance, i);
            dis2Region.push_back(initPair);
        }
        // (2): sort the Region with the distance
        std::sort(dis2Region.begin(), dis2Region.end(), [](const pair<double, size_t> &a, const pair<double, size_t> &b)
                  {
                      return a.first < b.first; // sort in acsending order of distance
                  });

        for (const auto &disIDPair : dis2Region)
        {
            size_t RegionID = disIDPair.second;
            if (_ptrSolver->_ptr_Parser->_flipflopLib.at(fflibID).Hight <= _PlaceRegionSet.at(RegionID).yEnd - _PlaceRegionSet.at(RegionID).yStart)
            {
                if (_ptrSolver->_ptr_Parser->_flipflopLib.at(fflibID).Width <= _PlaceRegionSet.at(RegionID).xEnd - _PlaceRegionSet.at(RegionID).xStart)
                {
                    regionIndex = RegionID;
                    isInCertainRegion = true;
                    break;
                }
            }
        }
    }

    assert(regionIndex != -1);
    assert(isInCertainRegion != false);

    // step2: check how many points can be place
    //   define some constant for further judgement
    size_t numPlaceRow = int((_PlaceRegionSet.at(regionIndex).yEnd - _PlaceRegionSet.at(regionIndex).yStart) / _PlaceRegionSet.at(regionIndex).siteHight);
    double siteHeight = _PlaceRegionSet.at(regionIndex).siteHight;
    double siteWidth = _PlaceRegionSet.at(regionIndex).siteWidth;
    double siteNum = _PlaceRegionSet.at(regionIndex).totalNumOfSites;
    double LFx_pos = _PlaceRegionSet.at(regionIndex).xStart;
    double LFy_pos = _PlaceRegionSet.at(regionIndex).yStart;

    int libWidth = std::ceil(_ptrSolver->_ptr_Parser->_flipflopLib.at(fflibID).Width / siteWidth);
    int libHeight = std::ceil(_ptrSolver->_ptr_Parser->_flipflopLib.at(fflibID).Hight / siteHeight);
    int radius_x = std::ceil(radius / siteWidth);
    int radius_y = std::ceil(radius / siteHeight);
    int x_start = int(std::floor((position.first - LFx_pos) / siteWidth)) - radius_x;
    int x_end = x_start + 2 * radius_x;
    if (x_start < 0)
    {
        x_start = 0;
    }
    if (x_end + libWidth > int(siteNum))
    {
        x_end = int(siteNum) - libWidth;
    }

    int y_start = std::floor((position.second - LFy_pos) / siteHeight) - radius_y;
    int y_end = y_start + 2 * radius_y;
    if (y_start < 0)
    {
        y_start = 0;
    }
    if (y_end + libHeight > numPlaceRow)
    {
        y_end = int(numPlaceRow) - libHeight;
    }

    for (int y = y_start; y <= y_end; y++)
    {
        // construct height constrait for the row
        vector<int> heightConstraint;
        if (x_end - x_start + 1 < 1)
        {
            continue;
        }
        else
        {
            heightConstraint.clear();
            heightConstraint.shrink_to_fit();
            heightConstraint.reserve(x_end - x_start + 1);
            for (int g = 0; g < int(x_end - x_start + 1); g++)
            {
                heightConstraint.push_back(numPlaceRow - y);
            }
            // update the heightConstraint for  this placement row
            for (int k = 0; k <= int(x_end - x_start); k++)
            {
                // std::cout << "k = " << k << std::endl;
                for (int j = y; j < y_end + libHeight; j++)
                {
                    // std::cout << "j = " << j << std::endl;
                    if (_availPosTable.at(regionIndex).at(j).at(k))
                    {
                        continue;
                    }
                    else
                    {
                        heightConstraint.at(k) = j - y;
                        break;
                    }
                }
            }
        }

        // check whether the position is safe
        int breakIndex = x_start;
        for (int x = x_start; x <= x_end; x++)
        {
            if (x < breakIndex)
            {
                continue;
            }

            bool isSafe = true;
            for (int index = x; index < x + libWidth; index++)
            {
                if (heightConstraint.at(index - x_start) < libHeight)
                {
                    isSafe = false;
                    breakIndex = index;
                    break;
                }
            }
            if (isSafe)
            {
                pair<double, double> initPair;
                initPair.first = double(x) * siteWidth + LFx_pos;
                initPair.second = double(y) * siteHeight + LFy_pos;
                output.push_back(initPair);
            }
        }
    }
    return output;
}

void Solver::PositionChecker::initPositionChecker()
{
    if (!isInit)
    {
        // initialize the data structure for the placement row, we assume every placement row has one id
        _Pid2PlaceRow.clear();
        _Pid2PlaceRow.reserve(_ptrSolver->_PlaceRow.size());
        for (size_t i = 0; i < _ptrSolver->_PlaceRow.size(); i++)
        {
            struct PlacementRow_checker PR;
            PR.x = _ptrSolver->_PlaceRow.at(i).x;
            PR.y = _ptrSolver->_PlaceRow.at(i).y;
            PR.totalNumOfSites = _ptrSolver->_PlaceRow.at(i).totalNumOfSites;
            PR.siteHight = _ptrSolver->_PlaceRow.at(i).siteHight;
            PR.siteWidth = _ptrSolver->_PlaceRow.at(i).siteWidth;
            _Pid2PlaceRow.push_back(PR);
        }

        // initialize the set of Gate, we assume the whole gate would have only one id
        _id2gate.clear();
        _id2gate.reserve(_ptrSolver->_Gate_arr.size());
        for (size_t i = 0; i < _ptrSolver->_Gate_arr.size(); i++)
        {
            _id2gate.push_back(&_ptrSolver->_Gate_arr[i]);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////
        /*

        there sould be no need to modified the code below in this function

        */
        /////////////////////////////////////////////////////////////////////////////////////////////////

        categorizePlaceRegion();
        constructAvailableTable();
    }

    isInit = true;
}

void Solver::PositionChecker::categorizePlaceRegion()
{
    // connect every placement row with their base point
    vector<pair<pair<double, double>, size_t>> PRpos2Pid;
    PRpos2Pid.reserve(_Pid2PlaceRow.size());
    for (size_t i = 0; i < _Pid2PlaceRow.size(); i++)
    {
        pair<double, double> pos(_Pid2PlaceRow.at(i).x, _Pid2PlaceRow.at(i).y);
        pair<pair<double, double>, size_t> initPair(pos, i);
        PRpos2Pid.push_back(initPair);
    }

    // sorting the place row by their base point, making the place rows with same x coordinate would align together
    std::sort(PRpos2Pid.begin(), PRpos2Pid.end(), [](const pair<pair<double, double>, size_t> &a, const pair<pair<double, double>, size_t> &b)
              {
                  if (a.first.first != b.first.first)
                  {
                      return a.first.first < b.first.first; // sort in acsending order by x coordinate
                  }
                  return a.first.second < b.first.second; // if the x coordinate is the same, sort in acsending order by y coordinate
              });

    // by checking the base point sequentialy, we can distinguish the connected uniform placement row
    _PlaceRegionSet.clear();
    _PlaceRegionSet.reserve(128);
    if (true)
    {
        assert(PRpos2Pid.size() != 0);
        pair<double, double> basePosition(_Pid2PlaceRow.at(PRpos2Pid.at(0).second).x, _Pid2PlaceRow.at(PRpos2Pid.at(0).second).y);
        double siteWidth = _Pid2PlaceRow.at(PRpos2Pid.at(0).second).siteWidth;
        double totalNumOfSites = _Pid2PlaceRow.at(PRpos2Pid.at(0).second).totalNumOfSites;
        double siteHight = _Pid2PlaceRow.at(PRpos2Pid.at(0).second).siteHight;
        struct PlaceRegion_checker firstRegion;
        firstRegion.xStart = basePosition.first;
        firstRegion.xEnd = firstRegion.xStart + siteWidth * totalNumOfSites;
        firstRegion.yStart = basePosition.second;
        firstRegion.siteHight = siteHight;
        firstRegion.siteWidth = siteWidth;
        firstRegion.totalNumOfSites = totalNumOfSites;
        _PlaceRegionSet.push_back(firstRegion);
    }
    if (PRpos2Pid.size() < 2)
    {
        _PlaceRegionSet.at(0).yEnd = _PlaceRegionSet.at(0).yStart + _Pid2PlaceRow.at(PRpos2Pid.at(0).second).siteHight;
    }
    else
    {
        double ybefore = _PlaceRegionSet.at(0).yStart;
        for (size_t i = 1; i < PRpos2Pid.size(); i++)
        {
            pair<double, double> basePosition = PRpos2Pid.at(i).first;
            double siteWidth = _Pid2PlaceRow.at(PRpos2Pid.at(i).second).siteWidth;
            double totalNumOfSites = _Pid2PlaceRow.at(PRpos2Pid.at(i).second).totalNumOfSites;
            double siteHight = _Pid2PlaceRow.at(PRpos2Pid.at(i).second).siteHight;
            // if the next place row is uniform and connect to the former one
            if ((basePosition.first == _PlaceRegionSet.back().xStart) && siteWidth == _PlaceRegionSet.back().siteWidth && siteHight == _PlaceRegionSet.back().siteHight && totalNumOfSites == _PlaceRegionSet.back().totalNumOfSites && basePosition.second == ybefore + _PlaceRegionSet.back().siteHight)
            {
                if (i == PRpos2Pid.size() - 1)
                {
                    // if this is the last place row
                    _PlaceRegionSet.back().yEnd = basePosition.second + _PlaceRegionSet.back().siteHight;
                    continue;
                }
                else
                {
                    ybefore = basePosition.second;
                    continue;
                }
            }
            else
            {
                // the next place row is not uniform or connect to the former one
                // (1): close the former region
                _PlaceRegionSet.back().yEnd = ybefore + _PlaceRegionSet.back().siteHight;
                // (2): add an new region
                struct PlaceRegion_checker NewRegion;
                NewRegion.xStart = basePosition.first;
                NewRegion.xEnd = NewRegion.xStart + siteWidth * totalNumOfSites;
                NewRegion.yStart = basePosition.second;
                NewRegion.siteHight = siteHight;
                NewRegion.siteWidth = siteWidth;
                NewRegion.totalNumOfSites = totalNumOfSites;
                _PlaceRegionSet.push_back(NewRegion);
                // update ybefore
                ybefore = basePosition.second;

                // if it is the last item, close the region
                if (i == PRpos2Pid.size() - 1)
                {
                    _PlaceRegionSet.back().yEnd = basePosition.second + _PlaceRegionSet.back().siteHight;
                }
            }
        }
    }
}

void Solver::PositionChecker::constructAvailableTable()
{
    _availPosTable.clear();
    _availPosTable.shrink_to_fit();
    _availPosTable.reserve(_PlaceRegionSet.size());
    for (const auto &region : _PlaceRegionSet)
    {
        //  define some constant for construct available table
        size_t numPlaceRow = int((region.yEnd - region.yStart) / region.siteHight);
        double siteHeight = region.siteHight;
        double siteWidth = region.siteWidth;
        double siteNum = region.totalNumOfSites;
        double LFx_pos = region.xStart;
        double LFy_pos = region.yStart;

        vector<vector<bool>> vec_init;
        vec_init.reserve(numPlaceRow);
        _availPosTable.push_back(vec_init);

        //  initialize _availPosTable
        for (size_t i = 0; i < numPlaceRow; i++)
        {
            vector<bool> initVec(int(siteNum), true);
            _availPosTable.back().push_back(initVec);
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
            if (URx > siteNum)
            {
                URx = int(siteNum);
            }
            int URy = std::ceil((upperRight.second - LFy_pos) / siteHeight);
            if (URy > int(numPlaceRow))
            {
                URy = int(numPlaceRow);
            }
            // set grid related in _availPosTable to false
            for (int y = LLy; y < URy; y++)
            {
                for (int x = LLx; x < URx; x++)
                {
                    _availPosTable.back().at(y).at(x) = false;
                }
            }
        }
    }
}

pair<double, double> Solver::PositionChecker::getGateLF(Inst::Gate *ptr) const
{
    return _ptrSolver->getGateLF(ptr);
}

pair<double, double> Solver::PositionChecker::getGateUR(Inst::Gate *ptr) const
{
    return _ptrSolver->getGateUR(ptr);
}