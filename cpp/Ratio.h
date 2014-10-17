#ifndef __RATIO_H__
#define __RATIO_H__

#include <string>
#include <vector>

using namespace std;

class Ratio
{
public:
    Ratio();
    Ratio(vector<PhyloTreeEdge> eEdges, vector<PhyloTreeEdge> fEdges);
    Ratio(double eLength, double fLength);
    ~Ratio();
    Ratio combine(Ratio r1, Ratio r2);
    double geoAvg(double d1, double d2);
    double geoAvg(vector<PhyloTreeEdge> edges);
    vector<PhyloTreeEdge> getEEdges();
    void addEEdge(PhyloTreeEdge edge);
    void addAllEEdges(vector<PhyloTreeEdge> edges);
    double getELength();
    void setELength(double eLen);
    vector<PhyloTreeEdge> getFEdges();
    void addFEdge(PhyloTreeEdge edge);
    void addAllFEdges(vector<PhyloTreeEdge> edges);
    double getFLength();
    void setFLength(double fLen);
    double getRation();
    double getTime();
    Ratio reverse();
    bool containsOriginalEEdge(Bipartition edge);
    string toString();
    string toStringCombTypeAndValue();
    string toStringCombType();
    string toStringJustValue();
    string toStringVerbose(vector<string> leaf2NumMap);
    Ratio clone();
    
private:
    double eLength;
    double fLength;
    vector<PhyloTreeEdge> eEdges;
    vector<PhyloTreeEdge> fEdges;
};

#endif /* __RATIO_H__ */