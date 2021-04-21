import numpy as np
import mahotas
from scipy import ndimage

def get_seeds(boundary, method='grid', next_id = 1,
             seed_distance = 10, boundary_thres = 0.5, label_nb = None):
    if method == 'grid':
        height = boundary.shape[0]
        width  = boundary.shape[1]

        seed_positions = np.ogrid[0:height:seed_distance, 0:width:seed_distance]
        num_seeds_y = seed_positions[0].size
        num_seeds_x = seed_positions[1].size
        num_seeds = num_seeds_x*num_seeds_y
        seeds = np.zeros_like(boundary).astype(np.int32)
        seeds[seed_positions] = np.arange(next_id, next_id + num_seeds).reshape((num_seeds_y,num_seeds_x))

    if method in ['minima', 'maxima_distance']:
        if method == 'minima':
            peak = mahotas.regmin(boundary)
        elif method == 'maxima_distance':
            # distance = mahotas.distance(boundary < boundary_thres)
            distance = ndimage.morphology.distance_transform_cdt(boundary < boundary_thres)
            peak = mahotas.regmax(distance)
        if label_nb is None:
            seeds, num_seeds = mahotas.label(peak)
        else:
            seeds, num_seeds = mahotas.label(peak, label_nb)

        seeds[seeds > 0] += next_id

    return seeds, num_seeds

def watershed(affs, seed_method, boundary_thres = 0.5, label_nb = None, seg_bg = True):
    fragments = np.zeros_like(affs[0]).astype(np.uint64)
    depth  = fragments.shape[0]
    next_id = 0
    for z in range(depth):
        affs_xy = 1.0 - 0.5*(affs[1,z] + affs[2,z])
        seeds, num_seeds = get_seeds(affs_xy, next_id=next_id, method=seed_method, \
                                     boundary_thres=boundary_thres, label_nb = label_nb)
        fragments[z] = mahotas.cwatershed(affs_xy, seeds)
        next_id += num_seeds

    return fragments
