#include "Bipartition.h"

void Bipartition::complement(int numLeaves) {
    for (int i = 0; i < numLeaves; i++) {
        partition.flip(i);
    }
}

bool Bipartition::contains(Bipartition e) {
}

bool Bipartition::contains(int i) {
    partition.get(i);
}