#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 2018

Useful functions for extracting information from images.
Taken from plotting_tools by Ami Tsuchida (atsuch@gmail.com)


@author: tsuchida
"""

def get_img_dat(img):
    """
    convenience fxn to get img data as numpy
    """
    import numpy as np
    import nibabel as nib
    return np.asanyarray(nib.load(img).dataobj) if isinstance(img, str) else np.asanyarray(img.dataobj)


def get_1d_dat(img, mask=None, mask_label=None, l_thr=None, u_thr=None, nonzero=True, nonna=True):
    """
    Function to grab values from image data in img (NiftiImage or path to img)
    and retuns 1d array of the shape (num voxels, ).

    If a mask is provided, only values from voxels within the mask are extracted.
    Additionally, if a mask label is provided, only values from voxels with that
    mask label are extracted.
    
    If l_ or u_thr is provided, values below or above the values are discarded.

    If nonzero=True, only nonzero values are extracted, and if nonna=True,
    only non-na values are kept.

    """
    import numpy as np
    img_dat = get_img_dat(img)

    if mask is not None:
         mask_dat = get_img_dat(mask)
         if mask_label is None:
             dat = img_dat[mask_dat > 0]
         else:
             dat = img_dat[mask_dat == mask_label]
    else:
        dat = img_dat.flatten()
    
    if nonna:
        dat = dat[~np.isnan(dat)]

    if nonzero:
        if not nonna:
            dat = dat[np.array([e != 0 if ~np.isnan(e) else False for e in dat], dtype=bool)]
        else:
            dat = dat[dat != 0]
    
    if l_thr is not None:
        dat = dat[dat > l_thr]
        
    if u_thr is not None:
        dat = dat[dat < u_thr]
        
    return dat


def summarize_img_intensity_distribution(img, 
                                         stats=['mean', 'std', 'median', 
                                                '5p', '10p', '25p', '75p', '90p', '95p',
                                                'min', 'max', 'count'],
                                         custom_func_dict=None,
                                         mask=None, mask_label=None,
                                         l_thr=None, u_thr=None,
                                         nonzero=True, nonna=True,
                                         out=None):
    """
    Function to get a basic distribution stats summary from a img with or without mask.
    
    Available measures are;
    -mean and std
    -median
    -5th, 10th, 25th, 75th, 90th, 95th percentile
    -min and max
    -count of the voxels
    
    Returns a dict containing these summary measures. If out is not None, it will save
    the dict as a csv file.
    
    """
    import numpy as np
    import pandas as pd
    flat_dat = get_1d_dat(img, mask=mask, mask_label=mask_label,
                          l_thr=l_thr, u_thr=u_thr,
                          nonzero=nonzero, nonna=nonna)
    
    summary_dat = {}
    def_dict = {'mean': lambda x: np.mean(x),
                'std': lambda x: np.std(x),
                'median': lambda x: np.median(x),
                '5p': lambda x: np.percentile(x, 5),
                '10p': lambda x: np.percentile(x, 10),
                '25p': lambda x: np.percentile(x, 25),
                '75p': lambda x: np.percentile(x, 75),
                '90p': lambda x: np.percentile(x, 90),
                '95p': lambda x: np.percentile(x, 95),
                'min': lambda x: np.min(x),
                'max': lambda x: np.max(x),
                'count': lambda x: x.shape[0]}
    
    # Add any custom func
    stats_dict = def_dict.copy()
    if custom_func_dict is not None:
        stats_dict.update(custom_func_dict)
        stats += custom_func_dict.keys()
        
    for measure in stats:
        if flat_dat.size == 0:
            summary_dat[measure] = np.nan
        else:
            summary_dat[measure] = stats_dict[measure](flat_dat)
        
    if out is not None:
        summary_df = pd.DataFrame(summary_dat, columns=stats, index=[0])
        summary_df.to_csv(out, index=False)
   
    return summary_dat
        
    
def get_histmat_df(img_idx, im_list, bin_num=100, col_strform=None, density=False):
    """
    Function to get DataFrame containing a histogram values of images in im_list.

    The bins are determined by global distribution by taking 3 IQR +/- median.

    If density=True, histogram values are density of total voxels.

    Returns DataFrame of the length im_list/img_idx.
    """
    import numpy as np
    import pandas as pd
    
    # Grab nonzero (and also non-NA) data from im_path
    nonzero_data = []
    for im in im_list:
        im_dat = get_1d_dat(im)
        nonzero_data.append(im_dat)
    stacked = np.hstack(nonzero_data)

    # Get the range of plotting
    # Here we define the range by median +/- 3 IQR of the global distribution.
    # Global bins are defined by taking the linspace between the range.
    iqr = (np.percentile(stacked, 75) - np.percentile(stacked, 25))
    (l, u) = ((np.median(stacked) - 3*iqr), (np.median(stacked) + 3*iqr))
    bins = np.linspace(l, u, bin_num)
    bins = np.insert(bins, 0, 0)
    bins = np.append(bins, np.max(stacked))

    # Get histmat
    hists = []
    for dat in nonzero_data:
        hist, bins = np.histogram(dat, bins=bins)
        # calculate rather than using density option since it may cause problems
        # when using unequal bins (for min and max bins)
        if density:
            hist = np.divide(hist, np.sum(hist, dtype=np.float))
        hists.append(hist)
    histmat =  np.vstack([h.reshape((1, h.shape[0])) for h in hists])

    # Use str-formatted bins for columns
    if col_strform is None:
        col_strform = "%0.1f" if np.absolute(np.min(l)) > 0.1 else "%0.1e"
    columns = [col_strform % val for val in bins][:-1]
    columns[0] = ("<" + col_strform ) % l
    columns[-1] = (">" + col_strform ) % u
    cols = pd.Index(columns, name="Image Intensity")

    histmat_df = pd.DataFrame(data=histmat, index=img_idx, columns=cols)

    return histmat_df


def get_hist_df(img_idx, im_list, bins=None, num_bins=100, bin_strform=None, density=False):
    """
    Function to get DataFrame containing a histogram values of images in im_list.

    If the bins are not given, it will use the num_bins and set the range to the
    approximate  3 IQR +/- median of the global distribution, computed by taking 
    the mean of median and IQR values for each img.

    If density=True, histogram values are density of total voxels.

    Returns DataFrame of the length im_list/img_idx.
    """
    import numpy as np
    import pandas as pd
    
    # Get median and iqr from each img from im_path if bins are not provided 
    # to calculatethe bin
    if bins is None:
        stats_dict = {'median':[], 'iqr':[], 'min':[], 'max':[]}
        for im in im_list:
            im_summary = summarize_img_intensity_distribution( 
                            im, stats=['median', '25p', '75p', 'min', 'max'])
            stats_dict['median'].append(im_summary['median'])
            iqr = im_summary['75p'] - im_summary['25p']
            stats_dict['iqr'].append(iqr)
            stats_dict['min'].append(im_summary['min'])
            stats_dict['max'].append(im_summary['max'])
            
        # bins are calculated as the mean of median +/- mean of 3 IQR
        global_stats = {}
        for key, item in stats_dict.items():
            if key == 'max':
                global_stats[key] = np.max(np.asarray(item))
            elif key == 'min':
                global_stats[key] = np.min(np.asarray(item))
            else:
                global_stats[key] = np.mean(np.asarray(item))

        lower = global_stats['median'] - 3 * global_stats['iqr']
        upper = global_stats['median'] + 3 * global_stats['iqr']
        bins = np.linspace(lower, upper, num_bins)
        bins = np.insert(bins, 0, global_stats['min'])
        bins = np.append(bins, global_stats['max'])
 

    # Compute histogram values for each column in cols_to_plot and keep them in a df
    # For bins col, use str-formatted bins
    if bin_strform is None:
        bin_strform = "%0.1f" if np.absolute(np.min(global_stats['min'])) > 0.1 else "%0.1e"
    str_bins = ["-".join([bin_strform, bin_strform]) % (l, u) for l, u in zip(bins[:-1], bins[1:])]
    hist_df = pd.DataFrame({'bins': str_bins})
    for idx, img in zip(img_idx, im_list):
        dat = get_1d_dat(img)
        hist, bins = np.histogram(dat, bins=bins)
        # calculate rather than using density option since it may cause problems
        # when using unequal bins (for min and max bins)
        if density:
            hist = np.divide(hist, np.sum(hist, dtype=np.float))
        hist_df[idx] = hist
    return hist_df
