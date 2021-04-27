import numpy as np
from .region_graph import mapping_id

def findConnectEdge(ii, i0, i1):
    ii2 = ii.copy()
    sc = True
    g0 = []
    g1 = []
    while sc:
        g0 = np.unique(ii2[(ii2==i0).max(axis=1)])
        g1 = np.unique(ii2[(ii2==i1).max(axis=1)])
        sc = np.in1d(g0,g1).sum()==0
    return ii[sc==0]
        
def somaBFS(ii, m0, check_num=10):
    # ii: Nx2
    # m0: soma ids (bigger than others)

    # remove redundant 
    ii = ii[ii[:,0] != ii[:,1]]
    ii2 = ii.copy() # BFS
    gid = np.ones(ii.shape[0]) > 0

    rll = np.arange(ii2.max()+1).astype(ii.dtype)
    check_id = 0
    while ii2.max() > m0:
        if check_id == 0:
            # check if exist false merge in the end
            ii_g = ii[gid]
            out = mapping_id(ii_g[:,0], ii_g[:,1])
            ui, uc = np.unique(out[m0+1:], return_counts=True)
            # only keep the fm-affected ones
            rl = np.zeros(ii_g.max()+1, np.uint8)
            rl[ui[uc>1]] = 1
            gg = rl[out[ii_g]].max(axis=1)
            ii_g2 = ii2[gid]
            ii_g2[gg == 0] = 0
            ii2[gid] = ii_g2

            print('#soma left: ', (uc>1).sum())
        check_id = (check_id + 1) % check_num

        print((ii2>m0).sum())
        if (ii2>m0).sum() == 0:
            import pdb; pdb.set_trace()

        # can't have two soma ids in one row 
        bid = ii2.min(axis=1) > m0
        gid[bid] = 0
        ii2[bid] = 0

        # only use seg that connects to one soma
        # can't have two soma ids merge to the same id
        jj = ii2.max(axis=1) > m0
        ii_j = ii2[jj]
        gid_j = gid[jj]
        ii_j_min = ii_j.min(axis=1)
        ii_j_max = ii_j.max(axis=1)
        _, bb = np.unique(ii_j_min, return_index=True)
        if len(bb) < len(ii_j_min):
            ii_j[:] = 0
            gid_j[:] = 0
            ii_j[bb] = ii2[jj][bb]
            gid_j[bb] = 1
            ii2[jj] = ii_j
            gid[jj] = gid_j

        rll[ii_j_min[bb]] = ii_j_max[bb]

        # merge a step
        ii2 = rll[ii2]
        ii2[ii2[:,0]==ii2[:,1]] = 0

    return ii[gid]
