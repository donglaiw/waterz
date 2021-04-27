import numpy as np
from .region_graph import mapping_id

def somaBFS(ii,m0):
    # ii: Nx2
    # m0: soma ids (bigger than others)

    # remove redundant
    ii = ii[ii[:,0]!=ii[:,1]]

    # find m0-relevant ids
    out = mapping_id(ii[:,0], ii[:,1])
    ui, uc = np.unique()

    rl = np.zeros(ii.max()+1, np.uint8)
    rl[out[m0+1:]] = 1
    gg = rl[ii].max(axis=1) > 0
    

    ii2 = ii[gid]
    gid = np.ones(ii.shape[0]) > 0
    rll = np.arange(ii2.max()+1).astype(ii.dtype)

    while ii.max() > m0:
        print((ii>m0).sum())
        # can't have two soma ids in one row 
        bid = ii.min(axis=1) > m0
        gid[bid] = 0
        ii[bid] = 0

        # only use seg that connects to one soma
        # can't have two soma ids merge to the same id
        jj = ii.max(axis=1) > m0
        ii_j = ii[jj]
        ii_j_min = ii_j.min(axis=1)
        ii_j_max = ii_j.max(axis=1)
        gid_j = gid[jj]
        aa,bb = np.unique(ii_j_min, return_index=True)
        if len(bb) < len(ii_j_min):
            ii_j[:] = 0
            gid_j[:] = 0
            ii_j[bb] = ii[jj][bb]
            gid_j[bb] = 1

            ii[jj] = ii_j
            gid[jj] = gid_j

        # merge a step
        rll[ii_j_min[bb]] = ii_j_max[bb]
        ii = rll[ii]

        # remove redundant edges
        gid[ii[:,0]==ii[:,1]] = 0
        ii[ii[:,0]==ii[:,1]] = 0

    return ii2[gid]
