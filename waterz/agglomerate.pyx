from libcpp.vector cimport vector
from libc.stdint cimport uint8_t, uint64_t, uint32_t
from libcpp cimport bool
import numpy as np
cimport numpy as np

def rgToArr(rg):
    # convert waterz rg into array
    num_rg = len(rg)
    rg_id = np.zeros([num_rg,2],np.uint32)
    rg_sc = np.zeros(num_rg, np.uint8)
    for i in range(num_rg):
        rg_id[i] = [rg[i]['u'],rg[i]['v']]
        rg_sc[i] = rg[i]['score']
    # smaller rg_sc, the closer 
    rg_sid = np.argsort(rg_sc)

    return rg_id[rg_sid], rg_sc[rg_sid]

def agglomerate(
    affs, 
    thresholds = [2], 
    gt = None, 
    fragments = None, 
    aff_threshold_low  = 1, 
    aff_threshold_high = 254, 
    rg_opt = 0,
    return_merge_history = False,
    return_region_graph=False):

    # the C++ part assumes contiguous memory, make sure we have it (and do 
    # nothing, if we do)
    if not affs.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous affinity arrray (avoid this by passing C_CONTIGUOUS arrays)")
        affs = np.ascontiguousarray(affs)
    if gt is not None and not gt.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous ground-truth arrray (avoid this by passing C_CONTIGUOUS arrays)")
        gt = np.ascontiguousarray(gt)
    if fragments is not None and not fragments.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous fragments arrray (avoid this by passing C_CONTIGUOUS arrays)")
        fragments = np.ascontiguousarray(fragments)

    print("Preparing segmentation volume...")

    if fragments is None:
        volume_shape = (affs.shape[1], affs.shape[2], affs.shape[3])
        segmentation = np.zeros(volume_shape, dtype=np.uint32)
        find_fragments = True
    else:
        segmentation = fragments.astype(np.uint32)
        find_fragments = False
    

    cdef WaterzState state;
    
    if rg_opt == 0: # do segmentation
        state = __initialize(affs, segmentation, gt, aff_threshold_low, aff_threshold_high, find_fragments)

        thresholds.sort()
        for threshold in thresholds:
            merge_history = mergeUntil(state, threshold)
            result = (segmentation,)
            if gt is not None:
                stats = {}
                stats['V_Rand_split'] = state.metrics.rand_split
                stats['V_Rand_merge'] = state.metrics.rand_merge
                stats['V_Info_split'] = state.metrics.voi_split
                stats['V_Info_merge'] = state.metrics.voi_merge
                result += (stats,)
            if return_merge_history:
                result += (merge_history,)
            if return_region_graph:
                result += (rgToArr(getRegionGraph(state)),)
            if len(result) == 1:
                yield result[0]
            else:
                yield result
        free(state)
    elif rg_opt in [1, 2, 3]: # segmentation -> rg
        yield __getRegionGraph(affs, segmentation, rg_opt)
    elif rg_opt == 4: # load rg
        # affs -> rg_score
        # segmentation -> rg_id
        # state = __initializeFromRg(segmentation[:,0], segmentation[:,1], affs)
        pass
    else:
        raise ValueError("Unknown region graph options")


def __getRegionGraph(np.ndarray[uint8_t, ndim=4] aff,
                np.ndarray[uint32_t, ndim=3] seg,
                np.uint32_t rg_opt):
    '''Return the initial region graph
    
    :param aff: the affinity predictions - an array of x, y, z, c where c == 0
                is the affinity prediction for x, c == 1 is the affinity
                prediction for y and c == 2 is the affinity prediction for z
    :param seg: the segmentation after finding basins
    :returns: a region graph as a 3-tuple of numpy 1-d arrays of affinity, 
              ID1 and ID2
    '''
    cdef uint8_t*    aff_data
    cdef uint32_t* seg_data
    aff_data = &aff[0,0,0,0]
    seg_data = &seg[0,0,0]

    width, height, depth = seg.shape[0], seg.shape[1], seg.shape[2]
    rg = rgFromSeg(width, height, depth, aff_data, seg_data, rg_opt)
    return rgToArr(rg);

def __initialize(
        np.ndarray[uint8_t, ndim=4] affs,
        np.ndarray[uint32_t, ndim=3]     segmentation,
        np.ndarray[uint32_t, ndim=3]     gt = None,
        aff_threshold_low  = 1,
        aff_threshold_high = 254,
        find_fragments = True):

    cdef uint8_t*    aff_data
    cdef uint32_t* segmentation_data
    cdef uint32_t* gt_data = NULL

    aff_data = &affs[0,0,0,0]
    segmentation_data = &segmentation[0,0,0]
    if gt is not None:
        gt_data = &gt[0,0,0]

    return initialize(
        affs.shape[1], affs.shape[2], affs.shape[3],
        aff_data,
        segmentation_data,
        gt_data,
        aff_threshold_low,
        aff_threshold_high,
        find_fragments)

"""
def __initializeFromRg(
        np.ndarray[uint32_t, ndim=1]      rg_id1,
        np.ndarray[uint32_t, ndim=1]      rg_id2,
        np.ndarray[uint8_t,  ndim=1]      rg_score):

    cdef uint8_t*    rg_score_data;
    cdef uint32_t*   rg_id1_data;
    cdef uint32_t*   rg_id2_data;

    rg_id1_data = &rg_id1[0]
    rg_id2_data = &rg_id2[0]
    rg_score_data = &rg_score[0]

    return initializeFromRg(
        max(rg_id1.max(), rg_id2.max()),
        np.uint32(len(rg_id1)),
        rg_id1_data,
        rg_id2_data,
        rg_score_data)
"""


cdef extern from "frontend_agglomerate.h":
    struct Metrics:
        double voi_split
        double voi_merge
        double rand_split
        double rand_merge

    struct Merge:
        uint32_t a
        uint32_t b
        uint32_t c
        uint8_t score

    struct ScoredEdge:
        uint32_t u
        uint32_t v
        uint8_t score

    struct WaterzState:
        int     context
        Metrics metrics

    WaterzState initialize(
            size_t          width,
            size_t          height,
            size_t          depth,
            const uint8_t*  affinity_data,
            uint32_t*       segmentation_data,
            const uint32_t* groundtruth_data,
            uint8_t         affThresholdLow,
            uint8_t         affThresholdHigh,
            bool            findFragments);

    WaterzState initializeFromRg(
            uint32_t        num_node,
            uint32_t        num_edge,
            uint32_t*       rg_id1,
            uint32_t*       rg_id2,
            uint8_t*        rg_score);

    vector[Merge] mergeUntil(
            WaterzState& state,
            uint8_t      threshold)

    vector[ScoredEdge] getRegionGraph(WaterzState& state)

    void free(WaterzState& state)
    
    # add region graph 
    vector[ScoredEdge] rgFromSeg(
		uint32_t     width,
		uint32_t     height,
		uint32_t     depth,
		const uint8_t* affinity_data,
		uint32_t*    segmentation_data,
		uint32_t     rg_opt);

