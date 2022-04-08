import numpy as np
import mahotas
from scipy import ndimage

def get_seeds(boundary, method='grid', \
             seed_distance = 10, boundary_thres = 0.5, label_nb = None):
    if method == 'grid':
        height = boundary.shape[0]
        width  = boundary.shape[1]

        seed_positions = np.ogrid[0:height:seed_distance, 0:width:seed_distance]
        num_seeds_y = seed_positions[0].size
        num_seeds_x = seed_positions[1].size
        num_seeds = num_seeds_x*num_seeds_y
        seeds = np.zeros_like(boundary).astype(np.int32)
        seeds[seed_positions] = np.arange(1, 1 + num_seeds).reshape((num_seeds_y,num_seeds_x))

    elif method in ['minima', 'maxima_distance']:
        if method == 'minima':
            peak = mahotas.regmin(boundary)
        elif method == 'maxima_distance':
            # edt: slow, maxima discontinous -> more strip-seg
            #distance = mahotas.distance(boundary < boundary_thres)
            distance = ndimage.morphology.distance_transform_cdt(boundary < boundary_thres)
            peak = mahotas.regmax(distance)
        if label_nb is None:
            seeds, num_seeds = mahotas.label(peak)
        else:
            seeds, num_seeds = mahotas.label(peak, label_nb)

    return seeds, num_seeds

def watershed(affs, seeds=None, seed_method='maxima_distance', boundary_thres = 0.5, label_nb = None, seg_bg = True):
    fragments = np.zeros_like(affs[0]).astype(np.uint32)
    depth  = fragments.shape[0]
    # 3D watershed is too slow -> do 2D watershed
    next_id = np.uint32(0)
    for z in range(depth):
        boundary = 1.0 - 0.5*(affs[1,z].astype(np.float32) + affs[2,z]) / 255.0
        if seeds is None:
            seeds_z, num_seeds = get_seeds(boundary, method = seed_method, boundary_thres = boundary_thres, label_nb = label_nb)
            fragments[z] = mahotas.cwatershed(boundary, seeds_z)
            fragments[z][fragments[z] > 0] += next_id
            next_id += num_seeds
        else:
            seeds_z = seeds[z]
            fragments[z] = mahotas.cwatershed(boundary, seeds_z)
        if seg_bg: # assign bg seg
            fragments[z][(affs[:,z]==0).max(axis=0)] = 0

        print('\tinit:', z, next_id)
    return fragments
