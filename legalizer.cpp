#include <string>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
#include <algorithm>
#include <cmath>
#include <assert.h>

#include "inst.h"
#include "parser.h"
#include "legalizer.h"

bool Solver::legalizer::legalize()
{
    initLegalizer();
    bool isLegalSuccess = true;
    for (size_t i = 0; i < _PlaceRegionSet.size(); i++)
    {
        if (!legalRegion(_FFinPlaceRegion[i], &_PlaceRegionSet[i]))
        {
            isLegalSuccess = false;
        }
    }
    return isLegalSuccess;
}

void Solver::legalizer::initLegalizer() // initialize the basic data structure of legalizer
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /*

    below would be some data structure need to be adjusted to your own

    */
    /////////////////////////////////////////////////////////////////////////////////////////////////

    if (!isInit)
    {
        // initialize the data structure for the placement row, we assume every placement row has one id
        _Pid2PlaceRowPtr.clear();
        _Pid2PlaceRowPtr.reserve(_ptrSolver->_PlaceRow.size());
        for (size_t i = 0; i < _ptrSolver->_PlaceRow.size(); i++)
        {
            _Pid2PlaceRowPtr.push_back(&_ptrSolver->_PlaceRow[i]);
        }

        // initialize the data structure for FF, we assume every FF would have 1 instance and we will assign id for it
        _id2ffPtr.clear();
        _id2ffPtr.reserve(_ptrSolver->_FF_D_arr.size());
        for (size_t i = 0; i < _ptrSolver->_FF_D_arr.size(); i++)
        {
            _id2ffPtr.push_back(&_ptrSolver->_FF_D_arr[i]);
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
    }

    categorizeFF();
    isInit = true;
}

void Solver::legalizer::categorizePlaceRegion()
{
    // connect every placement row with their base point
    vector<pair<pair<double, double>, size_t>> PRpos2Pid;
    PRpos2Pid.reserve(_Pid2PlaceRowPtr.size());
    for (size_t i = 0; i < _Pid2PlaceRowPtr.size(); i++)
    {
        pair<double, double> pos = getLowerLeftPlaceRow(_Pid2PlaceRowPtr[i]);
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
        pair<double, double> basePosition = getLowerLeftPlaceRow(_Pid2PlaceRowPtr[PRpos2Pid[0].second]);
        double siteWidth = getSiteWidth(_Pid2PlaceRowPtr[PRpos2Pid[0].second]);
        double totalNumOfSites = getTotalNumOfSites(_Pid2PlaceRowPtr[PRpos2Pid[0].second]);
        double siteHight = getSiteHeight(_Pid2PlaceRowPtr[PRpos2Pid[0].second]);
        struct PlaceRegion firstRegion;
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
        _PlaceRegionSet[0].yEnd = _PlaceRegionSet[0].yStart + getSiteHeight(_Pid2PlaceRowPtr[PRpos2Pid[0].second]);
    }
    else
    {
        double ybefore = _PlaceRegionSet[0].yStart;
        for (size_t i = 1; i < PRpos2Pid.size(); i++)
        {
            pair<double, double> basePosition = PRpos2Pid[i].first;
            double siteWidth = getSiteWidth(_Pid2PlaceRowPtr[PRpos2Pid[i].second]);
            double totalNumOfSites = getTotalNumOfSites(_Pid2PlaceRowPtr[PRpos2Pid[i].second]);
            double siteHight = getSiteHeight(_Pid2PlaceRowPtr[PRpos2Pid[i].second]);
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
                struct PlaceRegion NewRegion;
                NewRegion.xStart = basePosition.first;
                NewRegion.xEnd = NewRegion.xStart + siteWidth * totalNumOfSites;
                NewRegion.yStart = basePosition.second;
                NewRegion.siteHight = siteHight;
                NewRegion.siteWidth = siteWidth;
                NewRegion.totalNumOfSites = totalNumOfSites;
                _PlaceRegionSet.push_back(NewRegion);

                // if it is the last item, close the region
                if (i == PRpos2Pid.size() - 1)
                {
                    _PlaceRegionSet.back().yEnd = basePosition.second + _PlaceRegionSet.back().siteHight;
                }
            }
        }
    }
}

void Solver::legalizer::categorizeFF()
{
    // initialize _FFinPlaceRegion
    _FFinPlaceRegion.clear();
    _FFinPlaceRegion.reserve(_PlaceRegionSet.size());
    for (size_t i = 0; i < _PlaceRegionSet.size(); i++)
    {
        list<size_t> initList;
        _FFinPlaceRegion.push_back(initList);
    }
    // sort the FF to the placeRegion
    vector<bool> isSorted(_id2ffPtr.size(), false);
    for (size_t ffID = 0; ffID < _id2ffPtr.size(); ffID++)
    {
        if (isSorted[ffID])
        {
            continue;
        }

        bool isInCertainRegion = false;
        pair<double, double> FFposition = getFFPosition(_id2ffPtr[ffID]);
        // check if the ff is in any of placeRegion
        for (size_t i = 0; i < _PlaceRegionSet.size(); i++)
        {
            if (FFposition.first >= _PlaceRegionSet[i].xStart && FFposition.first <= _PlaceRegionSet[i].xEnd)
            {
                if (FFposition.second >= _PlaceRegionSet[i].yStart && FFposition.second <= _PlaceRegionSet[i].yEnd)
                {
                    isInCertainRegion = true;
                    _FFinPlaceRegion[i].push_back(ffID);
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
                if (FFposition.first < _PlaceRegionSet[i].xStart)
                {
                    distance += std::fabs(FFposition.first - _PlaceRegionSet[i].xStart);
                }
                else if (FFposition.first > _PlaceRegionSet[i].xEnd)
                {
                    distance += std::fabs(FFposition.first - _PlaceRegionSet[i].xEnd);
                }
                // calculate the y distance
                if (FFposition.second < _PlaceRegionSet[i].yStart)
                {
                    distance += std::fabs(FFposition.second - _PlaceRegionSet[i].yStart);
                }
                else if (FFposition.second > _PlaceRegionSet[i].yEnd)
                {
                    distance += std::fabs(FFposition.second - _PlaceRegionSet[i].yEnd);
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
                if (getFFHeight(_id2ffPtr[ffID]) <= _PlaceRegionSet[RegionID].yEnd - _PlaceRegionSet[RegionID].yStart)
                {
                    if (getFFWidth(_id2ffPtr[ffID]) <= _PlaceRegionSet[RegionID].xEnd - _PlaceRegionSet[RegionID].xStart)
                    {
                        _FFinPlaceRegion[RegionID].push_back(ffID);
                        isInCertainRegion = true;
                        break;
                    }
                }
            }
            if (!isInCertainRegion)
            {
                std::cout << "Error: id " << ffID << " cannot be put in any placement region." << std::endl;
            }
        }

        // mark the related FF
        vector<size_t> groupMemberID = getGroupMem(_id2ffPtr[ffID]);
        for (const auto &id : groupMemberID)
        {
            isSorted[id] = true;
        }
    }
}

bool Solver::legalizer::legalRegion(const list<size_t> &_FFList, struct PlaceRegion *_ptr_placeRegion)
{
    int _MaxIteration = 30;
    std::cout << "start to get placement row information." << std::endl;
    // start the algorithm part
    // define some constant for the algorithm
    size_t numPlaceRow = int((_ptr_placeRegion->yEnd - _ptr_placeRegion->yStart) / _ptr_placeRegion->siteHight);
    double siteHeight = _ptr_placeRegion->siteHight;
    double siteWidth = _ptr_placeRegion->siteWidth;
    double siteNum = _ptr_placeRegion->totalNumOfSites;
    double LFx_pos = _ptr_placeRegion->xStart;
    double LFy_pos = _ptr_placeRegion->yStart;

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
    std::cout << "step 2: throw the FF to the nearests placement grid and sort it to placement row." << std::endl;
    // step 2: throw all FF to the nearests placement grid and sort it to placement row
    vector<bool> check_list4ffInit(_id2ffPtr.size(), false);
    for (const auto &id : _FFList)
    {
        check_list4ffInit[id] = true;
    }
    for (size_t i = 0; i < _id2ffPtr.size(); i++)
    {
        // only legal the FF in the _FFlist
        if (!check_list4ffInit[i])
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
                          return a.second <= b.second; // sort in acsending order
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

            //   as the assign for one row complete,
            //   we start to record the solution and start the legalize for next row
            int counter = 0;
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
                    // set all grouped instance's position to the result
                    vector<size_t> member = getGroupMem(_id2ffPtr[_legalist[i][j].first]);
                    for (const auto id : member)
                    {
                        setFFPosition(_id2ffPtr[id], position);
                    }
                }
            }
            std::cout << "legalization successfully. \n"
                      << std::endl;
            return true;
        }
    }
    std::cout << "legalization fail !!!!!!!!! \n"
              << std::endl;
    return false;
}

pair<double, double> Solver::legalizer::getLowerLeftPlaceRow(struct PlacementRow *ptr) const
{
    pair<double, double> position(ptr->x, ptr->y);
    return position;
}

pair<double, double> Solver::legalizer::getFFPosition(Inst::FF_D *ptr) const
{
    return _ptrSolver->getFFPosition(ptr);
}

vector<size_t> Solver::legalizer::getGroupMem(Inst::FF_D *ptr) const
{
    return _ptrSolver->getGroupMem(ptr);
}

double Solver::legalizer::getFFWidth(Inst::FF_D *ptr) const
{
    return _ptrSolver->getFFWidth(ptr);
}

double Solver::legalizer::getFFHeight(Inst::FF_D *ptr) const
{
    return _ptrSolver->getFFHeight(ptr);
}

pair<double, double> Solver::legalizer::getGateLF(Inst::Gate *ptr) const
{
    return _ptrSolver->getGateLF(ptr);
}

pair<double, double> Solver::legalizer::getGateUR(Inst::Gate *ptr) const
{
    return _ptrSolver->getGateUR(ptr);
}

void Solver::legalizer::setFFPosition(Inst::FF_D *ptr, pair<double, double> &pos)
{
    _ptrSolver->setFFPosition(ptr, pos);
    return;
}
