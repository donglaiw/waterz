import h5py
import numpy as np
import scipy

def mappingToList(mapping):
    # conver to sparse list for efficient i/o 
    ind = np.arange(len(mapping)).astype(mapping.dtype)
    ind = ind[ind!=mapping]
    return np.hstack([ind.reshape(-1,1), mapping[ind].reshape(-1,1)])

def getScoreFunc(scoreF):
    # aff50_his256
    config = {x[:3]: x[3:] for x in scoreF.split('_')}
    if 'aff' in config:
        if 'his' in config and config['his']!='0':
            if config['ran'] == '255':
                return 'One255Minus<HistogramQuantileAffinity<RegionGraphType, %s, ScoreValue, %s>>' % (config['aff'],config['his'])
            else:
                return 'OneMinus<HistogramQuantileAffinity<RegionGraphType, %s, ScoreValue, %s>>' % (config['aff'],config['his'])
        else:
            if config['ran'] == '255':
                return 'One255Minus<QuantileAffinity<RegionGraphType, '+config['aff']+', ScoreValue>>'
            else:
                return 'OneMinus<QuantileAffinity<RegionGraphType, '+config['aff']+', ScoreValue>>'
    elif 'max' in config:
            if config['ran'] == '255':
                return 'One255Minus<MeanMaxKAffinity<RegionGraphType, '+config['max']+', ScoreValue>>'
            else:
                return 'OneMinus<MeanMaxKAffinity<RegionGraphType, '+config['max']+', ScoreValue>>'

def writeh5(filename, datasetname, dtarray):                                                         
    fid=h5py.File(filename,'w')
    ds = fid.create_dataset(datasetname, dtarray.shape, compression="gzip", dtype=dtarray.dtype)
    ds[:] = dtarray
    fid.close()


def create_border_mask(input_data, max_dist, background_label,axis=0):
    """
    Overlay a border mask with background_label onto input data.
    A pixel is part of a border if one of its 4-neighbors has different label.
    
    Parameters
    ----------
    input_data : h5py.Dataset or numpy.ndarray - Input data containing neuron ids
    target : h5py.Datset or numpy.ndarray - Target which input data overlayed with border mask is written into.
    max_dist : int or float - Maximum distance from border for pixels to be included into the mask.
    background_label : int - Border mask will be overlayed using this label.
    axis : int - Axis of iteration (perpendicular to 2d images for which mask will be generated)
    """
    target = input_data.copy()
    sl = [slice(None) for d in xrange(len(target.shape))]

    for z in xrange(target.shape[axis]):
        sl[ axis ] = z
        border = create_border_mask_2d(input_data[tuple(sl)], max_dist)
        target_slice = input_data[tuple(sl)] if isinstance(input_data,h5py.Dataset) else np.copy(input_data[tuple(sl)]) 
        target_slice[border] = background_label
        target[tuple(sl)] = target_slice
    return target

def create_border_mask_2d(image, max_dist):
    """
    Create binary border mask for image.
    A pixel is part of a border if one of its 4-neighbors has different label.
    
    Parameters
    ----------
    image : numpy.ndarray - Image containing integer labels.
    max_dist : int or float - Maximum distance from border for pixels to be included into the mask.

    Returns
    -------
    mask : numpy.ndarray - Binary mask of border pixels. Same shape as image.
    """
    max_dist = max(max_dist, 0)
    
    padded = np.pad(image, 1, mode='edge')
    
    border_pixels = np.logical_and(
        np.logical_and( image == padded[:-2, 1:-1], image == padded[2:, 1:-1] ),
        np.logical_and( image == padded[1:-1, :-2], image == padded[1:-1, 2:] )
        )

    distances = scipy.ndimage.distance_transform_edt(
        border_pixels,
        return_distances=True,
        return_indices=False
        )

    return distances <= max_dist
    
