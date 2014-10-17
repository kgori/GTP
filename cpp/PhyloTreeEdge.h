#ifndef __PHYLOTREE_EDGE_H__
#define __PHYLOTREE_EDGE_H__

#include <vector>
#include <string>
#include <bitset>
#include "Bipartition"

using namespace std;

class PhyloTreeEdge : public Bipartition
{
public:
    PhyloTreeEdge();
    PhyloTreeEdge(bitset edge);
    PhyloTreeEdge(EdgeAttribute attrib);
    PhyloTreeEdge(EdgeAttribute attrib, int originalID);
    PhyloTreeEdge(EdgeAttribute attrib, Bipartition originalEdge, int originalID);
    PhyloTreeEdge(Bipartition edge, EdgeAttribute attrib, int originalID);
    PhyloTreeEdge(BitSet edge, EdgeAttribute attrib, BitSet originalEdge, int originalID);
    ~PhyloTreeEdge();
    string printEdgesVerbose(vector<PhyloTreeEdge> edges, vector<string> leaf2NumMap, bool originalEdges);
    double getLength();
    bool isZero();
    string toString();
    string toStringVerbose();
    PhyloTreeEdge clone();
    bool equals(PhyloTreeEdge e);
    bool sameBipartition(PhyloTreeEdge e);
    bool sameBipartition(Bipartition e);
    Bipartition asSplit();
    Bipartition getOriginalEdge();
    void setOriginalEdge(Bipartition originalEdge);
    int getOriginalID();
    void setOriginalID(int originalID);
    EdgeAttribute getAttribute();
    void setAttribute();

private:
    EdgeAttribute attribute;
    Bipartition originalEdge;
    int originalID;
    
};

#endif /* __PHYLOTREE_EDGE_H__ */
