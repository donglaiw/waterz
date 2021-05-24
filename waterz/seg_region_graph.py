import numpy as np
from .region_graph import merge_id
from .seg_util import readh5

def findConnectEdge(ii, i0, i1):
    # TODO
    ii2 = ii.copy()
    sc = True
    g0 = []
    g1 = []
    while sc:
        g0 = np.unique(ii2[(ii2==i0).max(axis=1)])
        g1 = np.unique(ii2[(ii2==i1).max(axis=1)])
        sc = np.in1d(g0,g1).sum()==0
    return ii[sc==0]
        
def somaBFS(ii, m0=400000000, check_num=10):
    # ii: Nx2
    # m0: soma ids (bigger than others)

    # remove redundant 
    ii = ii[ii[:,0] != ii[:,1]]
    ii2 = ii.copy() # BFS
    ii3 = ii.copy() # global
    gid = np.ones(ii.shape[0]) > 0

    rll = np.arange(ii2.max()+1).astype(ii.dtype)
    check_id = 0
    while ii2.max() > m0:
        if check_id == 0:
            # check if exist false merge in the end
            ii_g = ii3[gid]
            out = merge_id(ii_g[:,0], ii_g[:,1])
            ui, uc = np.unique(out[m0+1:], return_counts=True)
            # only keep the fm-affected ones
            rl = np.zeros(ii_g.max()+1, np.uint8)
            rl[ui[uc>1]] = 1
            gg = rl[out[ii_g]].max(axis=1)
            ii_g2 = ii2[gid]
            ii_g2[gg == 0] = 0
            ii2[gid] = ii_g2
            
            # once unrelated, always unrelated
            # as we only remove edges
            ii_g[gg == 0] = 0
            ii3[gid] = ii_g
            print('#soma left: ', (uc>1).sum(), uc[uc>1].sum())
            if uc[uc>1].sum() == 0:
                return ii[gid]
        check_id = (check_id + 1) % check_num

        print((ii2>m0).sum())

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


def branchIoU(bbs, sid, rl, count_a, fn_count, fn_bb, fn_iou, chunk_sz = 228, thd_iou=0.5, m0=400000000, zmax=5699):
    # fn_iou[0]: z -> z-1
    # fn_iou[1]: z-1 -> z
    bb = bbs[bbs[:,0]==sid].squeeze()
    # step forward one slice to compute prev stat
    zran = range(max(0, bb[1]*4-4), min((bb[2]+1) * 4 + 4, zmax))
    gid = []
    # print('seg: %d (#z=%d, sz=%d)' % (sid, len(zran), bb[-1]))
    ii_pre = []
    for z in zran:
        """
        if z % 300 == 0:
            print(z)
        """
        cid = z // chunk_sz 
        zid = z - cid * chunk_sz
        count = np.loadtxt(fn_count % (z - zid)).astype(int)
        bb_z = readh5(fn_bb % z)[:,0]
        bb_z2 = bb_z.copy()
        bb_z2[bb_z2 < m0] += count_a[cid] + count[zid]
        # no connections among somas: >m0
        ii = list(bb_z[rl[bb_z2] == sid])

        if z>0 and len(ii) + len(ii_pre) > 0:
            # check upward link
            for kid in range(2):
                iou = readh5(fn_iou[kid] % z)
                iou = iou[np.in1d(iou[:, kid], ii) + np.in1d(iou[:, 1-kid], ii_pre)]
                # branch: small -> big
                # small seg: should be big enough, not diminishing to none
                # warning: exist size=0
                score = iou[:,-1]/(iou[:,-3:-1].min(axis=1).astype(float))
                ids = iou[score >= thd_iou, :2]
                
                ids[ids[:,kid]<m0, kid] += count_a[cid] + count[zid]
                if zid > 0:
                    ids[ids[:,1-kid]<m0,1-kid] += count_a[cid] + count[zid-1]
                else:
                    ids[ids[:,1-kid]<m0,1-kid] += count_a[cid-1] + count_pre[-2]

                ids = rl[ids]
                # not merge to other soma
                gid += list(ids[ids<m0])
                
        """
        if z == 3276:
        #if 179472776 in gid:
            import pdb; pdb.set_trace()
            db = readh5(fn_bb % (z-1))
            db[db[:,0]==46559]
        """
        count_pre, ii_pre = count, ii

    return np.unique(gid)

def branchIoUBFS(bbs, sid, rl, count_a, fn_count, fn_bb, fn_iou, chunk_sz = 228, thd_iou=0.5, m0=400000000, zmax=5699, thd_sz=200):
    did = []
    todo = [sid]
    rid = 0
    while len(todo) > 0:
        print('---------')
        print('round %d, #seg=%d' % (rid, len(todo)))
        did = did + todo
        todo2 = []
        for i, ii in enumerate(todo):
            gid = branchIoU(bbs, ii, rl, count_a, fn_count, fn_bb, fn_iou, chunk_sz = chunk_sz, thd_iou=thd_iou, m0=m0, zmax=zmax)
            gid = np.in1d(bbs[:,0], gid)
            aa, bb = bbs[gid,0], bbs[gid,-1]
            aa = aa[bb > thd_sz]
            todo2 += list(aa[np.in1d(aa, did, invert = True)])
            print('\t %d-th seg, # children=%d' % (i, len(todo2)))
        todo = list(np.unique(todo2))
        rid += 1
    return np.unique(did)
