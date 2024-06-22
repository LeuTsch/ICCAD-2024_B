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

Inst::Inst::~Inst()
{
    // delete &_position;
    NetID_related.clear();
    // delete &NetID_related;
}

Inst::PIO::PIO(string name, const pair<double, double> pos) : _Name(name)
{
    setPosition(pos.first, pos.second);
}
Inst::PIO::~PIO()
{
    _Name.clear();
}

Inst::Gate::Gate(string name, const pair<double, double> pos) : _Name(name)
{
    pinPosition.clear();
    pinName.clear();
    setPosition(pos.first, pos.second);
}
Inst::Gate::~Gate()
{
    _Name.clear();
    pinPosition.clear();
    pinName.clear();
}

Inst::FF_D::FF_D(string name, const pair<double, double> pos) : _Name_origin(name)
{
    faninCone.clear();
    grouped_member.clear();
    setPosition(pos.first, pos.second);
    _Name_now = "NEW" + name;
    grouped = false;
    _position_origin = pos;
    slack = 0;
}
Inst::FF_D::~FF_D()
{
    faninCone.clear();
    grouped_member.clear();
    _Name_now.clear();
}

Inst::FF_Q::FF_Q(string name, const pair<double, double> pos) : _Name_origin(name)
{
    fanoutCone.clear();
    setPosition(pos.first, pos.second);
    _Name_now = "NEW" + name;
    grouped = false;
    _position_origin = pos;
    slack = 0;
}

Inst::FF_Q::~FF_Q()
{
    fanoutCone.clear();
    _Name_now.clear();
}
