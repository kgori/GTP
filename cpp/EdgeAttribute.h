#ifndef __EDGE_ATTRIBUTE_H__
#define __EDGE_ATTRIBUTE_H__
#define TOLERANCE 0.000000000000001;

#include <string>
#include <vector>
using namespace std;

class EdgeAttribute
{
public:
	EdgeAttribute();
	EdgeAttribute(vector<double> vect);
	EdgeAttribute(string s);
	~EdgeAttribute();
	EdgeAttribute difference(EdgeAttribute a1, EdgeAttribute a2);
	EdgeAttribute add(EdgeAttribute a1, EdgeAttribute a2);
	EdgeAttribute weightedPairAverage(EdgeAttribute start, EdgeAttribute target, double position);
	EdgeAttribute zeroAttribute(int size);
	double getAttribute();
	void setEdgeAttribute(EdgeAttribute attrib);
	EdgeAttribute clone();
	string toString();
	bool equals(EdgeAttribute e);
	double norm();
	void scaleBy(double a);
	int size();
	void ensurePositive();
};

#endif /* __EDGE_ATTRIBUTE_H__ */