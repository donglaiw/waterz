from libc.stdint cimport uint8_t, uint64_t, uint32_t
from libcpp cimport bool
import numpy as np
cimport numpy as np

def mapping_id(id1, id2):
    mapping = np.arange(max(id1.max(), id2.max()) + 1).astype(id1.dtype)
    __mapping_id(id1, id2, mapping)
    return mapping 

def __mapping_id(np.ndarray[np.uint32_t, ndim=1] id1,
                  np.ndarray[np.uint32_t, ndim=1] id2,
                  np.ndarray[np.uint32_t, ndim=1] mapping):
    '''Find the global mapping of IDs from the region graph without count constraints
    
    The region graph should be ordered by decreasing affinity and truncated
    at the affinity threshold.
    :param id1: a 1D array of the lefthand side of the two adjacent regions
    :param id2: a 1D array of the righthand side of the two adjacent regions
    :returns: a 1D array of the global IDs per local ID
    '''


    cdef uint32_t* id1_data;
    cdef uint32_t* id2_data;
    cdef uint32_t* mapping_data;
    id1_data = &id1[0];
    id2_data = &id2[0];
    mapping_data = &mapping[0];

    do_mapping_id(id1_data, id2_data, len(id1), len(mapping), mapping_data);

    return mapping


cdef extern from "frontend_region_graph.h":
    void do_mapping_id(
        uint32_t*      id1,
        uint32_t*      id2,
        uint32_t       num_edge,
        uint32_t       num_id,
        uint32_t*      mapping);
