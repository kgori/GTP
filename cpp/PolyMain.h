#ifndef __POLY_MAIN_H__
#define __POLY_MAIN_H__ value

#include <string>
#include <vector>
#include "PhyloTree.h"

class PolyMain
{
public:
    PolyMain();
    ~PolyMain();
    static vector<PhyloTree> aTreesNoCommonEdges;
    static vector<PhyloTree> bTreesNoCommonEdges;
    static bool normalize = false;
    static void splitOnCommonEdge(PhyloTree t1, PhyloTree t2);
    static double getRobinsonFouldsDistance(PhyloTree t1, PhyloTree t2, bool normalise);
    static double getWeightedRobinsonFouldsDistance(PhyloTree t1, PhyloTree t2, bool normalise);
    static double getEuclideanDistance(PhyloTree t1, PhyloTree t2, bool normalise);
    static double getGeodesicDistance(PhyloTree t1, PhyloTree t2, bool normalise);
    static Geodesic getGeodesic(PhyloTree t1, PhyloTree t2);
    static Geodesic getGeodesicNoCommonEdges(Phylotree t1, PhyloTree t2);
}; 
#endif /* __POLY_MAIN_H__ */