#ifndef __BIPARTITE_GRAPH__
#define __BIPARTITE_GRAPH__

#include <vector>
#include "Vertex.h"

using namespace std;

class BipartiteGraph
{
public:
	BipartiteGraph(vector<vector<bool>> IncidenceMatrix, vector<double> Aweight, vector<double> Bweight);
	~BipartiteGraph();
	vector<vector<int>> vertex_cover(vector<int> Aindex, vector<int> Bindex)

private:
	vector<vector<bool>> edge;
	int nA, nB, n, i, j;
	vector<Vertex> Avertex, Bvertex;
	bool debug = false;

};

#endif /* __BIPARTITE_GRAPH__ */