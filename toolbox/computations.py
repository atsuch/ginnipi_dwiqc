# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 17:30:27 2015
Functions for the REST pipeline
@author: herve
"""

def get_bval_per_frame(bvalsFilename, bval_list, delimiter=None, out_file=None):
    """
    Get intended bval for each frame and save with numpy savetxt.

    {frame: [frame index],
     bval: [intended bval]}

    """
    import os.path as op
    import numpy as np

    if not op.exists(bvalsFilename):
        raise RuntimeError('bvals file not exist:' + bvalsFilename)

    if out_file is None:
        out_file = op.join(op.split(bvalsFilename)[0], 'bval_per_frame.txt')

    bvals = np.loadtxt(bvalsFilename, delimiter=delimiter)
    bStep = np.array(bval_list)

    bval_per_frame = []
    for i in range(0, bvals.size):
        ind = np.argmin(abs(bvals[i] - bStep))
        bval_per_frame.append(bval_list[ind])

    # save as txt
    np.savetxt(out_file, np.asarray(bval_per_frame))

    return out_file


def compute_voxel_percentage(in_file, mask_file, out_file=None, thr=0):
    '''
    Compute the number of voxels above thr (default=0) and outputs the percentage relative to the total voxels of the mask

    Parameters
    ----------
    in_file : 3d image file
    mask_file : 3d binary mask file
    out_file : name of the csv file to be saved

    Returns
    -------
    out_file : a csv file containing the summary of voxel counts and percentage
    '''
    import nibabel as nb
    import numpy as np
    import pandas as pd
    import os

    if out_file is None:
        out_file = os.path.join(os.getcwd(), 'voxel_occupancy_percentage.csv')

    dat = np.asanyarray(nb.load(in_file).dataobj)
    mask_dat = np.asanyarray(nb.load(mask_file).dataobj)

    dat[(mask_dat == 0) | (dat <= thr)] = 0
    num_vox_dat = np.count_nonzero(dat)
    num_vox_mask = np.count_nonzero(mask_dat)
    percent_voxels = float(num_vox_dat)/float(num_vox_mask)

    df = pd.DataFrame({"in_file_num_vox": [num_vox_dat],
                       "mask_num_vox": [num_vox_mask],
                       "percent_voxels": [percent_voxels]})

    df[["in_file_num_vox", "mask_num_vox", "percent_voxels"]].to_csv(out_file, index=False)

    return out_file

