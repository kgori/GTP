#ifndef __RATIO_SEQUENCE_H__
#define __RATIO_SEQUENCE_H__

#include <vector>
#include <string>

using namespace std;

class RatioSequence : vector<Ratio>
{
public:
    RatioSequence();
    RatioSequence(string ptA, string ptB);
    ~RatioSequence();
    RatioSequence interleave(RatioSequence rs1, RatioSequence rs2);
    RatioSequence getRandomRS(int dim);
    int getCombineCode();
    void setCombineCode(int c);
    Ratio getRatio(int i);
    bool isAscending();
    double getDistance();
    double getMinNonDesRSDistance();
    RatioSequence clone();
    RatioSequence getCombinedRS(int combineCode);
    RatioSequence getNonDesRSWithMinDist();
    RatioSequence getAscRSWithMinDist();
    RatioSequence reverse();
    string toStringValueAndRatio();
    string toStringValue();
    string toStringVerbose(vector<string> leaf2NumMap);
    string toStringCombType();
    string toStringCombTypeAndValue();
    
private:
    int combineCode = 0;
};

#endif /* __RATIO_SEQUENCE_H__ */