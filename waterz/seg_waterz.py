import os, sys
import numpy as np
import h5py
import json

from .seg_watershed import watershed
from .seg_util import create_border_mask, writeh5, getScoreFunc
from . import agglomerate


def getRegionGraph(affs, fragments, rg_opt = 1, merge_function = None, rebuild = True):
    # rg_opt=1: all seg
    # rg_opt=2: skip first slice
    # rg_opt=3: only the connection
    for rg in agglomerate(
            affs,
            fragments = fragments,
            scoring_function = getScoreFunc(merge_function),
            rg_opt = rg_opt,
            force_rebuild=rebuild):
        return rg

def waterz(
        affs,
        thresholds,
        output_prefix = './',
        merge_function = None,
        gt = None,
        gt_border = 25/4.0,
        fragments = None,
        fragments_opt = 0,
        fragments_seed_nb = 5,
        discretize_queue = 0,
        fragments_mask = None,
        aff_threshold  = [1, 254],
        return_seg = True,
        bg_thres = 1,
        return_rg = False,
        rebuild = True):

    # affs shape: 3*z*y*x
    thresholds = list(thresholds)
    # print("waterz at thresholds " + str(thresholds))

    if fragments is None:
        print('initial watershed')
        if fragments_opt != 0: # mahotas
            fragments = watershed(affs, seed_method='maxima_distance', label_nb = np.ones([fragments_seed_nb,fragments_seed_nb]), bg_thres = bg_thres)
            if fragments_mask is not None:
                fragments[fragments_mask==False] = 0


    outs = []
    outs_rg = []
    if gt is not None and gt_border !=0:
        gt = create_border_mask(gt, gt_border, np.uint32(0))

    for i,out in enumerate(agglomerate(
            affs,
            thresholds,
            gt = gt,
            aff_threshold_low  = aff_threshold[0],
            aff_threshold_high = aff_threshold[1],
            fragments = fragments,
            scoring_function = getScoreFunc(merge_function),
            discretize_queue = discretize_queue,
            return_region_graph = return_rg,
            force_rebuild=rebuild)):

        rebuild = False
        threshold = thresholds[i]
        output_basename = output_prefix+merge_function+'_%.2f'%threshold
        seg = None
        rg = None
        if gt is not None:
            if return_seg:
                seg = out[0]
            if return_rg:
                rg = out[1]
        else:
            if return_seg:
                if return_rg:
                    seg, rg = out[0], out[1]
                else:
                    seg = out
        if return_seg:
            if seg is not None:
                outs.append(seg.copy())
            if rg is not None:
                outs_rg.append([rg[0].copy(), rg[1].copy()])
        else:
            print("Storing segmentation...")
            writeh5(output_basename + '.hdf', 'main', seg)

        if gt is not None:
            metrics = out[1]
            print("Storing record...")
            record = {
                'threshold': threshold,
                'merge_function': merge_function,
                'discretize_queue': discretize_queue,
                'voi_split': metrics['V_Info_split'],
                'voi_merge': metrics['V_Info_merge'],
                'rand_split': metrics['V_Rand_split'],
                'rand_merge': metrics['V_Rand_merge'],
            }
            with open(output_basename + '.json', 'w') as f:
                json.dump(record, f)
    if seg is not None:
        if rg is not None:
            return outs, outs_rg
        else:
            return outs
    else:
        return outs_rg
