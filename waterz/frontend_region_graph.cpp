#include "frontend_region_graph.h"

void do_mapping_id(
     SegID* id1,
     SegID* id2,
     SegID  num_edge,
     SegID  num_id,
     SegID* mapping) {

    // assign to small seg ids
    SegID s1, s2;
    for (size_t i = 0; i < num_edge; ++i) {
        s1 = mapping[id1[i]];
        s2 = mapping[id2[i]];
        if (s1 > s2){
            mapping[id1[i]] = s2;
            mapping[id2[i]] = s2;
        } else {
            mapping[id1[i]] = s1;
            mapping[id2[i]] = s1;
        }
    }

    // skip 0th element: bg seg
    for (size_t i=1; i<num_id; ++i) {
        s1 = mapping[i];
        // root node: mapping to itself
        // is it mapped to the root node
        if (mapping[s1] == s1) continue;

        s1 = mapping[s1];
        while (mapping[s1] != s1){
            s1 = mapping[s1];
        }
        mapping[i] = s1;
    }
}
