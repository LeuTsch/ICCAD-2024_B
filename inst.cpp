#include <string>
#include <vector>
#include <utility>
#include "inst.h"

using std::pair;
using std::string;
using std::vector;

void Inst::Inst::setPosition(double x, double y)
{
    _position.first = x;
    _position.second = y;
}

Inst::PIO::PIO(string name, const pair<double, double> pos) : _Name(name)
{
    setPosition(pos.first, pos.second);
}

Inst::Gate::Gate(string name, const pair<double, double> pos) : _Name(name)
{
    pinPosition.clear();
    pinName.clear();
    setPosition(pos.first, pos.second);
}

Inst::FF_D::FF_D(string name, const pair<double, double> pos) : _Name_origin(name)
{
    NetID_related.clear();
    faninCone.clear();
    grouped_member.clear();
    setPosition(pos.first, pos.second);
    _Name_now = "NEW" + name;
    grouped = false;
    _position_origin = pos;
}

Inst::FF_Q::FF_Q(string name, const pair<double, double> pos) : _Name_origin(name)
{
    fanoutCone.clear();
    setPosition(pos.first, pos.second);
    _Name_now = "NEW" + name;
    grouped = false;
    _position_origin = pos;
}
