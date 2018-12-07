#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

def make_dwi_qc_traceplot(eddy_movement_rms, 
                          eddy_restricted_movement_rms,
                          eddy_outlier_map,
                          eddy_outlier_n_stdev_map,
                          eddy_outlier_n_sqr_stdev_map, 
                          afni_outliers_b0pre,
                          afni_outliers_nonb0pre_list,
                          afni_outliers_b0post,
                          afni_outliers_nonb0post_list,
                          bval_per_frame,
                          nonb0_bval_list,
                          out_dir=None):
    '''
    Creation of traceplot svg file for DWI QC across multiple shells, showing;
        - motion across time (from Eddy)
        - Eddy outlier across time (# slices per frame)
        - afni outlier across time (% voxels per frame)
      
    :param eddy_movement_rms: txt file containing the rms estimated by Eddy for each frame (rows) : 2nd col contains relative rms
    :type eddy_movement_rms: str
    :param eddy_restricted_movement_rms: txt file containing the restricted rms estimated by Eddy for each frame (rows) : 2nd col contains relative rms
    :type eddy_restricted_movement_rms: str
    :param eddy_outlier_map: txt file containing 1/0 identification for each slice (cols) per frame (rows)
    :type eddy_outlier_map: str
    :param eddy_outlier_n_stdev_map: txt file containing stdev for each slice (cols) per frame (rows)
    :type eddy_outlier_n_stdev_map: str
    :param eddy_outlier_n_sqr_stdev_map: txt file containing sqr stdev for each slice (cols) per frame (rows)
    :type eddy_outlier_n_sqr_stdev_map: str
    :param afni_outliers_b0pre: txt file containing % outlier voxels per frame (rows)
    :type afni_outliers_b0pre: str
    :param afni_outliers_nonb0pre_list: list of txt file containing % outlier voxels per frame (rows)
    :type afni_outliers_nonb0pre_list: list
    :param afni_outliers_b0post: txt file containing % outlier voxels per frame (rows)
    :type afni_outliers_b0post: str
    :param afni_outliers_nonb0post_list: list of txt file containing % outlier voxels per frame (rows)
    :type afni_outliers_nonb0post_list: list

    :param bval_per_frame: txt file containing intended bval for each frame
    :type bval_per_frame: str
    :param out_dir: output directory
    :type out_dir: str
    '''
    import os
    import os.path as op
    import numpy as np
    from ginnipi_dwiqc.toolbox.plotting_tools import QCtraceplot
    from ginnipi_dwiqc.toolbox.qualitycontrol_dwi import _get_dat_from_txt
    from seaborn import color_palette

    if out_dir is None:
        cwd = os.getcwd()
        out_dir = op.abspath(cwd)
    out_plot = op.join(out_dir, 'dwi_traceplot.svg')


    ## Plot 
    ############################################
    # initialization 
    figsize = (11, 8)

    title = 'DWI Eddy outlier and motion plot'
    myqc = QCtraceplot(out_dir=out_dir,
                       title=title,
                       figsize=figsize)

    # Set the color palette for plot backgrounds for bval > 0
    bpf = np.loadtxt(bval_per_frame)
    bvals = np.unique(bpf)
    nonb0_palette = color_palette('Pastel2', len(nonb0_bval_list))
    axvspan_d_list = []
    bval_change_pt = np.where(np.diff(bpf))[0]
    bval_change_pt = np.add(bval_change_pt, 1)
    bval_change_pt = np.insert(bval_change_pt, 0, 0)

    for i, pt in enumerate(bval_change_pt):
        x1 = pt
        x2 = bval_change_pt[i+1] if (i+1) < len(bval_change_pt) else (len(bpf) - 1)
        bval_at_pt = bpf[x1]
        if bval_at_pt > 0:
            axvspan_d = {'xmin': x1, 'xmax': x2, 'alpha': 0.5, 
                         'facecolor': nonb0_palette[nonb0_bval_list.index(bval_at_pt)]}
            axvspan_d_list.append(axvspan_d)
    # annotation for background colors
    base_annot_d = {'xycoords': 'figure fraction',
                    'ha':'center', 'va': 'center',
                    'fontsize': 9,
                    'bbox': {'boxstyle':'square', 'alpha':0.5, 'pad': 0.5, 'ec': 'grey'}}
    annot_d_list = []
    for i, val in enumerate(nonb0_bval_list):
        base_annot_d['s'] = 'b={}'.format(str(val))
        base_annot_d['xy'] = (0.7 + (0.1*i), 0.9)
        bbox_d = base_annot_d['bbox'].copy()
        bbox_d['fc'] = nonb0_palette[i]
        base_annot_d['bbox'] = bbox_d
        annot_d_list.append(base_annot_d.copy())

    # add motion rms traces
    rms = _get_dat_from_txt(eddy_movement_rms)
    res_rms = _get_dat_from_txt(eddy_restricted_movement_rms)
 
    if rms is not None and res_rms is not None:
        # input tseries should have shape (nsteps, ntrace)
        eddy_motion_tseries = np.transpose(np.vstack((rms[:, 1], res_rms[:, 1])))
        myqc.add_trace(eddy_motion_tseries,
                       'eddy_relative_rms',
                       True,
                       {'ylabel_name': 'Eddy relative rms',
                        'trace_name': ['rms', 'restricted_rms'],
                        'pcolor': color_palette('Spectral', 8)[0:2],
                        'hide_x': False,
                        'axvspan_d_list': axvspan_d_list,
                        'annot_d_list': annot_d_list})
            
    # add eddy outlier traces
    eddy_outlier_map_dat = _get_dat_from_txt(eddy_outlier_map, skiprows=1)
    eddy_outlier_n_stdev_map_dat = _get_dat_from_txt(eddy_outlier_n_stdev_map, skiprows=1)
    eddy_outlier_n_sqr_stdev_map_dat = _get_dat_from_txt(eddy_outlier_n_sqr_stdev_map, skiprows=1)

    if eddy_outlier_map_dat is not None:
        myqc.add_trace(eddy_outlier_map_dat.sum(axis=1),
                       'eddy_outlier_slices',
                       False,
                       {'ylabel_name': 'Number of eddy outlier slices',
                        'trace_name': ['eddy_outliers'],
                        'ylims': (0, None),
                        'pcolor': color_palette('Spectral_r', 8)[0:1],
                        'hide_x': False,
                        'axvspan_d_list': axvspan_d_list})

    if eddy_outlier_n_stdev_map_dat is not None:
        eddy_ol_stdev_tseries = np.transpose(
                                    np.vstack(
                                        (eddy_outlier_n_stdev_map_dat.mean(axis=1),
                                         eddy_outlier_n_stdev_map_dat.max(axis=1))
                                             )
                                            )
        myqc.add_trace(eddy_ol_stdev_tseries,
                       'eddy_outlier_stdev_stats',
                       True,
                       {'ylabel_name': 'Eddy obs-pred std',
                        'trace_name': ['mean', 'max'],
                        'pcolor': color_palette('PRGn', 11)[1:3],
                        'hide_x': False,
                        'axvspan_d_list': axvspan_d_list})

    if eddy_outlier_n_sqr_stdev_map_dat is not None:
        eddy_ol_sqr_stdev_tseries = np.transpose(
                                        np.vstack(
                                            (eddy_outlier_n_sqr_stdev_map_dat.mean(axis=1),
                                             eddy_outlier_n_sqr_stdev_map_dat.max(axis=1))
                                                 )
                                                )
        myqc.add_trace(eddy_ol_sqr_stdev_tseries,
                       'eddy_outlier_sqr_stdev_stats',
                       True,
                       {'ylabel_name': 'Eddy obs-pred sqr std',
                        'trace_name': ['mean', 'max'],
                        'pcolor': color_palette('PRGn_r', 11)[1:3],
                        'hide_x': False,
                        'axvspan_d_list': axvspan_d_list})

    # add afni outliers: we need to patch them together to create data per frame
    # first organize data into dict {bval:dat}
    afni_pre_dat_d, afni_post_dat_d = {}, {}
    for bval in bvals:
        if bval == 0.0:
            pre_dat = _get_dat_from_txt(afni_outliers_b0pre)
            post_dat = _get_dat_from_txt(afni_outliers_b0post)
        else:
            pre_dat = _get_dat_from_txt(afni_outliers_nonb0pre_list[nonb0_bval_list.index(bval)])
            post_dat = _get_dat_from_txt(afni_outliers_nonb0post_list[nonb0_bval_list.index(bval)])
        afni_pre_dat_d[bval] = pre_dat
        afni_post_dat_d[bval] = post_dat
        # keep track of any None data
        if pre_dat is None:
            afni_pre_dat_d['nonedata'] = afni_pre_dat_d.get('nonedata', 0) + 1
        if post_dat is None:
            afni_post_dat_d['nonedata'] = afni_post_dat_d.get('nonedata', 0) + 1

    # Now patch them up to create 1d dat if all the data are present
    afni_pre_dat, afni_post_dat = np.zeros(len(bpf)), np.zeros(len(bpf))
    for i, pt in enumerate(bval_change_pt):
        x2 = bval_change_pt[i+1] if (i+1) < len(bval_change_pt) else len(bpf)
        chunk_length = x2 - pt
        bval_at_pt = bpf[pt]
        if afni_pre_dat_d.get('nonedata', 0) == 0:
            afni_pre_dat[pt: pt+chunk_length], afni_pre_dat_d[bval_at_pt] = afni_pre_dat_d[bval_at_pt][:chunk_length], afni_pre_dat_d[bval_at_pt][chunk_length:]
        if afni_post_dat_d.get('nonedata', 0) == 0:
            afni_post_dat[pt: pt+chunk_length], afni_post_dat_d[bval_at_pt] = afni_post_dat_d[bval_at_pt][:chunk_length], afni_post_dat_d[bval_at_pt][chunk_length:]
    
    # add the trace
    if afni_pre_dat_d.get('nonedata', 0) == 0 or afni_post_dat_d.get('nonedata', 0) == 0:
        afni_ol_tseries = np.transpose(
                              np.vstack(
                                  (afni_pre_dat, afni_post_dat)
                                       )
                                      )
        myqc.add_trace(afni_ol_tseries,
                       'afni_percent_vox_outliers',
                       True,
                       {'ylabel_name': 'Afni percent vox outliers',
                        'trace_name': ['pre-processing', 'post-processing'],
                        'pcolor': color_palette('Accent_r', 8)[2:4],
                        'hide_x': False,
                        'axvspan_d_list': axvspan_d_list})

    # Finally: plot and save
    myqc.plot()
    myqc.fig.savefig(out_plot, bbox_inches='tight')
    myqc.fig.clf()

    return out_plot


def _get_dat_from_txt(in_file, skiprows=0):
    '''
    Get ndarray dat from in_file txt. Return None if the file does not exist.

    '''
    import os.path as op
    import numpy as np

    if op.exists(in_file):
        dat = np.loadtxt(in_file, skiprows=skiprows)
        return dat

    else:
        return None
