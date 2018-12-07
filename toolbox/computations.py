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


def identify_b0s(dcm_files, thr=6, return_series_num=False):
    '''
    From a set a of dicom files, this function will identify the series or
    mosaic frames with a b-value equal to zero and output two ordered lists
    of lists of dicom files for the b0 and actual dwi series.
    This function is SIEMENS based.

    Arguments:

    dcm_files:    A list of dicom files
    thr      :    The threshold for b-value=0 files (the DWI savvy
                  folks don't specify b=0, rather b0=5 because the scanners
                  can't reach a b-value of 0 anyway)

    Outputs:

    If return_series_num = False,

    b0_images: b0 DICOM files, in serial order
    dwi_series: dwi DICOM files, in serial order

    If return_series_num = True,

    b0_image_series number
    dwi_series_number
    '''

    import dicom

    bvalues = dict()
    series = dict()
    b0_series = dict()
    print('scanning DICOM instances')
    for dcm_file in dcm_files:
        dcm = dicom.read_file(dcm_file)
        try:
            bval = int(dcm[(0x0019,0x100c)].value)
        except:
            # if field is absent, SE-EPI and not quite a dwi b0 file
            bval = 0
            print('No b-value field!')
        series_number = dcm[(0x0020,0x0011)].value
        inst_number = dcm[(0x0020,0x0013)].value
        # build dictionaries for b=0 and b>0, per serie and instance number
        if bval <= thr:
            if series_number in b0_series.keys():
                b0_series[series_number][inst_number]=dcm_file
            else:
                b0_series[series_number]=dict()
                b0_series[series_number][inst_number]=dcm_file
        else:
            # That's a real diffusion weighted image
            if bval in bvalues.keys():
                bvalues[bval].append(dcm_file)
            else:
                bvalues[bval]=[dcm_file]
            if series_number in series.keys():
                series[series_number][inst_number]=dcm_file
            else:
                series[series_number]=dict()
                series[series_number][inst_number]=dcm_file

    print('reinjecting b0 instances in  dwi series')
    # reinject mosaic b0 in mosaic series
    for key in b0_series.keys():
        if key in series.keys():
            for instance in b0_series[key].keys():
                series[key][instance]=b0_series[key][instance]

    # now sort per series and instance number
    # list of lists for dwi mosaic series including b0s
    # list of lists for all b0s (individual 3D b0s + mosaic b0s)
    print('Sorting dwi series')
    dwi_series = []
    dwi_keys = list(map(int,series.keys()))
    dwi_keys.sort()
    for key in dwi_keys:
            dwi_keys2 = list(map(int,series[key].keys()))
            dwi_keys2.sort()
            dwi_series.append([series[key][key2] for key2 in dwi_keys2])

    print('Sorting b0 series')
    b0_images = []
    b0_keys = list(map(int,b0_series.keys()))
    b0_keys.sort()
    for key in b0_keys:
        b0_keys2 = list(map(int,b0_series[key].keys()))
        b0_keys2.sort()
        b0_images.append([b0_series[key][key2] for key2 in b0_keys2])
    print('FINISHED !')
    if return_series_num:
        return (b0_keys, dwi_keys)
    else:
        return (b0_images, dwi_series)


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

