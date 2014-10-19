#ifndef __BIPARTITION_H__
#define __BIPARTITION_H__

#include <bitset>
#include <string>
#include <vector>

using namespace std;

class Bipartition
{
public:
    Bipartition();
    Bipartition(bitset edge);
    Bipartition(string s);
    ~Bipartition();
    static string toStringVerbose(bitset edge, vector<string> leaf2NumMap);
    bitset getPartition();
    void setPartition(bitset edge);
    bool isEmpty() {
        return partition.none();
    }

    void addOne(int one) {
        partition.set(one);
    }

    void removeOne(int one) {
        partition.reset(one);
    }

    bool disjointFrom(Bipartition e);
    bool contains(Bipartition e);
    bool contains(int i);
    bool properlyContains(Bipartition e);
    bool crosses(Bipartition e);
    void complement(int numLeaves);
    string toString();
    bool equals(Bipartition e);
    Bipartition clone();
    bool isCompatibleWith(vector<Bipartition> splits);
    
protected:
    bitset partition;
};
#endif /* __BIPARTITION_H__ */
