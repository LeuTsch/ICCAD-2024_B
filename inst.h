#ifndef INST_H
#define INST_H

#include <string>
#include <vector>
#include <utility>

using std::pair;
using std::string;
using std::vector;

namespace Inst
{
    enum InstType
    {
        INST_PIO,
        INST_FF_D,
        INST_FF_Q,
        INST_GATE
    };

    class Inst
    {
    public:
        Inst() { NetID_related.clear(); }
        ~Inst();
        virtual InstType getType() const = 0;
        pair<double, double> getPosition() const { return _position; };
        void setPosition(double, double);
        virtual string getName() const = 0;
        vector<size_t> getRelatedNet() { return NetID_related; };
        void addRelatedNet(size_t id) { NetID_related.push_back(id); };

    private:
        pair<double, double> _position;
        vector<size_t> NetID_related;
    };

    class PIO : public Inst
    {
    public:
        PIO() {}
        PIO(string, const pair<double, double>);
        ~PIO();
        InstType getType() const { return INST_PIO; };
        string getName() const { return _Name; };
        void setName(string s) { _Name = s; };

    private:
        string _Name;
    };

    class Gate : public Inst
    {
    public:
        Gate()
        {
            pinPosition.clear();
            pinName.clear();
        }
        Gate(string, const pair<double, double>);
        ~Gate();

        // data part
        size_t gate_type;
        vector<pair<double, double>> pinPosition;
        vector<string> pinName;
        size_t PIN_OFFSET; // ID - PINOFFSET = i-th pin

        // function part
        InstType getType() const { return INST_GATE; };
        string getName() const { return _Name; };
        void setName(string s) { _Name = s; };

    private:
        string _Name;
    };

    class FF_D : public Inst
    {
    public:
        FF_D()
        {
            grouped = false;
            faninCone.clear();
            slack = 0;
        }
        FF_D(string, const pair<double, double>);
        ~FF_D();

        // data part
        bool grouped;
        double slack;
        size_t FF_type;
        // size_t ID_to_Q;
        vector<size_t> faninCone;
        vector<size_t> grouped_member; // store the global ID for FF_D, not the position in FF_D_arr, and would include itself
                                       // this vector would have at least 1 element(itself)
        // function part
        InstType getType() const { return INST_FF_D; };
        string getName() const { return _Name_now; };
        string getOriName() const { return _Name_origin; };
        void setName(string s) { _Name_now = s; };
        double getOriSlack() const { return _Slack_origin; };
        void setOriSlack(double s) { _Slack_origin = s; };
        void setClk(size_t id) { _clkID = id; };
        size_t getClk() const { return _clkID; };
        pair<double, double> getOriPosition() const { return _position_origin; };

    private:
        string _Name_now;
        string _Name_origin;
        double _Slack_origin;
        size_t _clkID;
        pair<double, double> _position_origin;
    };

    class FF_Q : public Inst
    {
    public:
        FF_Q()
        {
            fanoutCone.clear();
            grouped = false;
            slack = 0;
        }
        FF_Q(string, const pair<double, double>);
        ~FF_Q();

        // data part
        bool grouped;
        double slack;
        size_t FF_type;
        // size_t ID_to_D;
        vector<size_t> fanoutCone;

        // function part
        InstType getType() const { return INST_FF_Q; };
        string getName() const { return _Name_now; };
        string getOriName() const { return _Name_origin; };
        void setName(string s) { _Name_now = s; };
        void setClk(size_t id) { _clkID = id; };
        size_t getClk() const { return _clkID; };
        pair<double, double> getOriPosition() const { return _position_origin; };

    private:
        string _Name_now;
        string _Name_origin;
        size_t _clkID;
        pair<double, double> _position_origin;
    };

}
#endif