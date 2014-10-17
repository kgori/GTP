#ifndef __GEODESIC_H__
#define __GEODESIC_H__

#include <string>
#include <vector>

using namespace std;

class Geodesic
{
public:
    Geodesic(RatioSequence rs);
    Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges);
    Geodesic(RatioSequence rs, Vector<PhyloTreeEdge> cEdges, double leafContributionSquared);
    ~Geodesic();
    PhyloTree getTreeAt(PhyloTree t1, PhyloTree t2, double position);
    RatioSequence getRS();
    void setRS(RatioSequence rs);
    double getDist();
    void addCommonEdge(PhyloTreeEdge e);
    Geodesic clone();
    string toString();
    vector<PhyloTreeEdge> getCommonEdges(PhyloTree t1, PhyloTree t2, double position);
    vector<PhyloTreeEdge> getCommonEdges();
    void setCommonEdges(vector<PhyloTreeEdge> commonEdges);
    int numCommonEdges();
    int numTopologies();
    Geodesic reverse();
    double getLeafContributionSquared();
    void setLeafContributionSquared(double leafContributionSquared);

private:
    RatioSequence rs;
    vector<PhyloTreeEdge> commonEdges;
    double leafContributionSquared = 0;

};
#endif /* __GEODESIC_H__ */