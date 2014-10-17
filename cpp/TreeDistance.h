#ifndef __TREE_DISTANCE_H__
#define __TREE_DISTANCE_H__

#include <vector>
#include <string>
#include <unorderedmap>
#include "Bipartition.h"

using namespace std;

class TreeDistance
{
public:
    TreeDistance();
    ~TreeDistance();
    static void resetTreeDistanceState();
    static vector<Bipartition> zeroCol(int col, vector<Bipartition> m);
    static bool checkForDuplicateEdges(vector<PhyloTreeEdge> edges);
    static Geodesic getGeodesic2(PhyloTree t1, PhyloTree t2, string algorithm, string geoFile);
    static vector<PhyloTree> splitOnCommonEdge(PhyloTree t1, PhyloTree t2);
    static Geodesic getGeodesicRecursive(PhyloTree t1, PhyloTree t2, string algorithm);
    static Geodesic getPruned1GeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2);
    static Geodesic getPruned2GeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2);
    static Geodesic getDivideAndConquerGeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2);
    static getDivideAndConquerRSNoCommonEdges(PhyloTree t1, PhyloTree t2);
    static Ratio calculatRatio(Bipartition minEl, vector<Bipartition> m, vector<PhyloTreeEdge> eEdges, vector<PhyloTreeEdge> fEdges);
    static void getMaxPathSpacesAsRatioSeqs(vector<Bipartition> m, RatioSequence ratioSeq, vector<PhyloTreeEdge> eEdges, vector<PhyloTreeEdge> fEdges);
    static void getPruned2MaxPathSpacesAsRatioSeqs(vector<Bipartition> m, RatioSequence ratioSeq, vector<PhyloTreeEdge> eEdges, vector<PhyloTreeEdge> fEdges);
    static vector<int> binaryToVector(bitset bin);
    static vector<Bipartition> getMinElements(vector<Bipartition> m);
    static vector<PhyloTreeEdge> deleteZeroEdges(vector<PhyloTreeEdge> v);
    static vector<PhyloTreeEdge> myVectorClonePhyloTreeEdge(vector<PhyloTreeEdge> v);
    static vector<RatioSequence> myVectorCloneRatioSequence(vector<RatioSequence> v);
    static vector<Bipartition> removeMinElFrom(vector<Bipartition> m, Bipartition minEl);
    static vector<string> myVectorCloneString(vector<string> v);
    static double truncate(double d, int p);
    static double round(double d, int p);
    static vector<Bipartition> myVectorCloneBipartition(vector<Bipartition> v);

public:
    static vector<RatioSequence> finalRatioSeqs;
    static double minTreeDist = -1;
    static RatioSequence minTreeDistRatioSeq;
    static long numMaxPaths = 0;
    static long numPrunes = 0;
    static long numnodes = 0;

    static int pathToSearch = 0;
    // stores pairs of trees with no common edges.  Should be reset at each new distance calculation
    static vector<PhyloTree> aTreesNoCommonEdges = new Vector<PhyloTree>();
    static vector<PhyloTree> bTreesNoCommonEdges = new Vector<PhyloTree>();
    static double firstTreeDist = -1;
    static boolean rooted = true;  //holds if the trees are rooted or not.
    static boolean normalize = false;  // holds if we should normalize the tree split lengths
    static int verbose = 0;
    static unorderedmap<string, double> nodeHashtable;
    static unorderedmap<string, Geodesic> subTreeHashtable;
    static unorderedmap<string, RatioSequence> subTreeRSHashtable;
    
private:
    static long numBaseCase = 0;
    
};

#endif /* __TREE_DISTANCE_H__ */