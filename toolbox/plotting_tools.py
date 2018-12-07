#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 2018

Useful functions for plotting using matplotlib and seaborn.
Some functions from plotting_tools by Ami Tsuchida (atsuch@gmail.com)

-plot_img_intensity_distribution: Image intensity histogram with mean, median, etc, for QC

-QCtraceplot class: a class of traceplot object that can be used to generate linked traceplots
                    -Initialized with;
                        -out_dir: where plotted traces are saved as text files,
                        -title (optional): title of the plot
                        -figsize (optional): size of the combined plot (default:11x12)
                    -Methods;
                        -add_trace(): to add one or more traces for plotting
                        -plot(): to plot all the traces added by add_trace()
    This code is largely based on MRIQC. When using, please cite the following:
        Esteban O, Birman D, Schaer M, Koyejo OO, Poldrack RA, Gorgolewski KJ; 
        MRIQC: Advancing the Automatic Prediction of Image Quality in MRI from Unseen Sites; 
        PLOS ONE 12(9):e0184661; doi:10.1371/journal.pone.0184661.

@author: tsuchida
"""
import os.path as op
import numpy as np
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import gridspec as mgs
import seaborn as sns
from seaborn import color_palette
sns.set_style("whitegrid")
from textwrap import wrap


def get_lims(arr, delta=0.01):
    max_val = np.nanmax(arr)
    min_val = np.nanmin(arr)
    margin = np.absolute((max_val - min_val))* delta
    upper = max_val + margin
    lower = min_val - margin
    return (lower, upper)


def str_format_val(val, num_digits=3, signed=False):
    """
    Function to return float or scientific format of the val.
    
    Returns float format if absolute(val) >= 0.1
    Returns scientifc format if absolute(val) < 0.1
    """
    
    if np.absolute(val) >= 0.1:
        return '{1:+0.{0:d}f}'.format(num_digits, val) if signed else '{1:0.{0:d}f}'.format(num_digits, val)
    
    else:
        return '{1:+0.{0:d}e}'.format(num_digits, val) if signed else '{1:0.{0:d}e}'.format(num_digits, val)
    
    
def plot_img_intensity_distribution(img, mask, mask_label=None,
                                    color='royalblue', title=None,
                                    stats=['mean', 'median', 'min', 'max', '5p', '95p'],
                                    stats_style={'mean': {'line_style': ':',
                                                          'color': 'm',
                                                          'y_relpos': 0.8},
                                                 'median': {'line_style': '--',
                                                            'color': 'royalblue',
                                                            'y_relpos': 0.9},
                                                 'min': {'line_style': ':',
                                                         'color': 'slategray',
                                                         'y_relpos': 0.1},
                                                 'max': {'line_style': ':',
                                                         'color': 'slategray',
                                                         'y_relpos': 0.1},
                                                 '5p': {'line_style': '--',
                                                        'color': 'royalblue',
                                                        'y_relpos': 0.2},
                                                 '95p': {'line_style': '--',
                                                         'color': 'royalblue',
                                                         'y_relpos': 0.2}},
                                    out_plot=None,
                                    out_stats=None):
    """
    Uses seaborn distplot to plot image intensity distribution within the mask.
    Along with the histogram and KDE plot, it will show image distribution stats
    that can be computed by summarize_img_intensity_distribution fxn. 
    
    By default it will display the following values;
    -mean
    -median
    -min, max
    -5th and 95th percentile
    
    If out_stats is not None, it will save all available summary stats as a csv.
    
    """
    import os
    import os.path as op
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from ginnipi_dwiqc.toolbox.image_utils import get_1d_dat, summarize_img_intensity_distribution
    sns.set_style("whitegrid")

    masked_im_dat = get_1d_dat(img, mask=mask, mask_label=mask_label)

    cwd = os.getcwd()
    out_dir = op.abspath(cwd)
    if out_plot is None:
        out_plot = op.join(out_dir, 'img_intensity_distribution.png')
    else:
        out_plot = op.join(out_dir, out_plot)
    if out_stats is not None:
        out_stats = op.join(out_dir, out_stats)


    ## Plot 
    ############################################
    # initialization
    title = title if title is not None else 'image intensity plot'

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    sns.distplot(masked_im_dat, hist=True, kde=True, rug=False, color=color, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Image intensity")
    ax.set_ylabel("Frequency")
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()

    # Get stat information of interest
    stat_summary = summarize_img_intensity_distribution(img,
                                                        mask=mask,
                                                        mask_label=mask_label,
                                                        out=out_stats)
    if not np.nan in stat_summary.values():
        # Check if the automatic xlim exceeds 25p - 3IQR or 75p + 6IQR. If it does,
        # reset the xlim so that the histogram is not overly skewed.
        iqr = stat_summary['75p'] - stat_summary['25p']
        l_lim = stat_summary['25p'] - 6*iqr
        u_lim = stat_summary['75p'] + 6*iqr
        new_l_xlim = l_lim if xlim[0] < l_lim else xlim[0]
        new_u_xlim = u_lim if xlim[1] > u_lim else xlim[1]
        ax.set_xlim(new_l_xlim, new_u_xlim)

        # Annotation parameters
        bbox_style = 'round'
        bbox_linewidth = 0
        bbox_edgecolor = 'none'

        arrow_relativepos = (0.2, 0.5)

        # Add vline and annotation for each stat
        for stat in stats:
            # vline
            if stat_summary[stat] > new_l_xlim and stat_summary[stat] < new_u_xlim:
                ax.axvline(x=stat_summary[stat],
                           linestyle=stats_style[stat]['line_style'],
                           color=stats_style[stat]['color'], lw=1)

            # Annotation
            if stat_summary['median'] < 0.1 or iqr < 0.1 or stat_summary['median'] > 10 or iqr > 10:
                label = '{}={:.3e}'.format(stat, stat_summary[stat])
            else:
                label = '{}={:.3f}'.format(stat, stat_summary[stat])
        
            posy = ylim[0] + stats_style[stat]['y_relpos']*(ylim[1] - ylim[0])
            if stat_summary[stat] < new_l_xlim:
                posx = new_l_xlim
                xytext = (12, 0)
                arrow_style = 'wedge,tail_width=0'
            elif stat_summary[stat] > new_u_xlim:
                posx = new_u_xlim
                xytext = (-12, 0)
                arrow_style = 'wedge,tail_width=0'
            else:
                posx = stat_summary[stat]
                xytext = (12, 0)
                arrow_style = 'wedge,tail_width=1.'

            ax.annotate(label,
                        xy=(posx, posy),
                        xytext=xytext,
                        textcoords='offset points',
                        va='center',
                        color='w',
                        size=10,
                        bbox=dict(boxstyle=bbox_style,
                                  fc=stats_style[stat]['color'],
                                  ec=bbox_edgecolor,
                                  lw=bbox_linewidth),
                        arrowprops=dict(arrowstyle=arrow_style,
                                        lw=bbox_linewidth,
                                        patchA=None,
                                        patchB=None,
                                        fc=stats_style[stat]['color'],
                                        ec=bbox_edgecolor,
                                        relpos=arrow_relativepos))

    fig.savefig(out_plot)

    return out_plot, out_stats


class QCtraceplot(object):
    """
    Generation of multiple linked traceplot
    """

    def __init__(self, out_dir, title=None, figsize=(11, 12)):
        '''
        :param bval_per_frame: txt file containing intended bval for each frame
        :type bval_per_frame: str
        :param out_dir: output directory
        :type out_dir: str
        :param title: title of the summary plot
        :type title: str
        :param figsize: size of the output png
        :type figsize: tuple of float
        '''
        ## Set attributes
        ############################################

        # Working directory 
        self.workingDir = out_dir

        ## Set plot
        ############################################
        # Figure size
        self.fig = plt.figure(figsize=figsize)

        # Font size
        mpl.rcParams['font.size'] = 7

        # Title
        if title is not None:
            self.fig.suptitle(title, fontsize=10)

        # Tab of traces that will be plotted
        # traceTab : one or more trace / plot

        self.traceTab = []
        self.plot_heightTab = []

    def add_trace(self, data, data_name, is_multi, kwargs):
        '''
        Add a trace to traceTab
        '''
        self.traceTab.append((data, kwargs))
        trace_height = 1 if is_multi else 0.3
        self.plot_heightTab.append(trace_height)
        out_dir = self.workingDir

        np.savetxt(op.join(out_dir, '{}.txt'.format(data_name)), data)

    def plot(self):
        ''' 
        Plot traces that are in traceTab.
        
        trace : plot one or more trace in one plot 
        '''

        # Length of the tab 
        ntraceTab = len(self.traceTab)

        # Grid that organize the figure
        nrows = ntraceTab 
        grid = mgs.GridSpec(nrows, 1, 
                            wspace=0.0, hspace=0.5,
                            height_ratios = self.plot_heightTab) 
        grid_id = 0

        # Plot trace ...
        for idx, (tseries, kwargs) in enumerate(self.traceTab):
            # plot
            self.traceplot(tseries, 
                           grid[grid_id], 
                           **kwargs) 

            # update grid idx
            grid_id += 1

        setattr(self, 'grid', grid)

    def traceplot(self, tseries, gridpos, gs_dist=None, ylabel_name=None,
                  trace_name=None, units=None, hide_x=True, pcolor=None,
                  cutoff=None, cutoff_name=True, ylims=None,
                  axvspan_d_list=None, annot_d_list=None):
        '''
        Plot one or more traces on a graph
        
        :param tseries: traces to plot
        :type tseries: list or ndarray 
        :param gridpos: position of this subplot the plot grid 
        :type gridpos: int
        :param gs_dist: 
        :type gs_dist: 
        :param name:  name of the graph
        :type name: str
        :param normalize: 
        :type normalize: boolean
        :param units: units of the tseries 
        :type units: str
        :param hide_x: hide (or not) X axis on the subplot 
        :type hide_x: boolean
        :param color: color of the trace
        :type color: color
        :param cutoff: cut-off value to plot
        :type cutoff: list of float 
        :param cutoff_name: display cut-off value
        :type cutoff_name: boolean
        :param ylims: limits of the Y axis
        :type ylims: tuple of floats
        '''
        ## Variable initialization 
        #############################################################################
        # Numpy array
        tseries = np.array(tseries)

         # reshape if tseries is one-dimensional
        if len(tseries.shape) == 1:
            tseries = tseries.reshape((-1, 1))

        (ntsteps, ntrace) = tseries.shape
 
        tseries_array = np.transpose(tseries)
        
        # Color palette
        if pcolor == None:
            pcolor = ['b'] * ntrace

        ## Plot
        #############################################################################
        # Define nested GridSpec
        gs = mgs.GridSpecFromSubplotSpec(1, 2, 
                                         subplot_spec=gridpos,
                                         width_ratios=[1, 100], 
                                         wspace=0.0)

        # Plot data 
        ax_ts = plt.subplot(gs[1])
        for trace in range(ntrace):
            ax_ts.plot(tseries_array[trace], linewidth = 1, color=pcolor[trace])

        # Set subplot
        ax_ts.grid(False)
        for side in ["top", "right"]:
            ax_ts.spines[side].set_color('none')
            ax_ts.spines[side].set_visible(False)

        # Set X axis
        ax_ts.set_xlim((0, ntsteps - 1))
        if not hide_x:
            xtick_step = ntsteps / 5
            xticks = list(range(0, ntsteps)[::int(xtick_step)])
            ax_ts.set_xticks(xticks)
            ax_ts.set_xlabel('frame #')
            ax_ts.spines["bottom"].set_position(('outward', 2))
            ax_ts.xaxis.set_ticks_position('bottom')
        else:
            ax_ts.set_xticklabels([])
            ax_ts.spines["bottom"].set_color('none')
            ax_ts.spines["bottom"].set_visible(False)

        # Set Y axis
        if ylabel_name is not None:
            var_label = ylabel_name
            if units is not None:
                var_label += (' [{}]').format(units)
            ax_ts.set_ylabel('\n'.join(wrap(var_label, 15)))
        ax_ts.spines["left"].set_position(('outward', 30))
        ax_ts.yaxis.set_ticks_position('left')
        def_ylims = [0.95 * tseries[~np.isnan(tseries)].min(),
                     1.1 * tseries[~np.isnan(tseries)].max()]
        if ylims is not None:
            if ylims[0] is not None:
                def_ylims[0] =  ylims[0]
            if ylims[1] is not None:
                def_ylims[1] =  ylims[1]
        ax_ts.set_ylim(def_ylims)
        yticks = sorted(def_ylims)
        ax_ts.set_yticks(yticks)
        ax_ts.set_yticklabels(['%.05f' % y for y in yticks])
        yrange = def_ylims[1] - def_ylims[0]

        # Add vspan if axvspan_d_list is provided
        if axvspan_d_list is not None:
            for axvspan_d in axvspan_d_list:
                ax_ts.axvspan(**axvspan_d)
        
        ## Annotate graph with cutoff and mean values
        #############################################################################
        # Annotation parameters
        fontsize = 7
        text_color = 'w'
        bbox_style = 'round'
        bbox_linewidth = 0
        bbox_edgecolor = 'none'
        arrow_style = 'wedge,tail_width=.9'
        arrow_linewidth = bbox_linewidth
        arrow_edgecolor = bbox_edgecolor

        # For each trace, add annoation and mean values
        ntrace_name = len(trace_name)
        if ntrace_name == ntrace:
            for trace in range(ntrace):
                # annotation
                xpos = trace*0.1
                label = trace_name[trace]
                position = (xpos,  1)
                text_offset = (xpos, 1)
                bbox_color = pcolor[trace]
                ax_ts.annotate(label, 
                               xy = position, 
                               xytext = text_offset, 
                               xycoords = 'axes fraction',
                               textcoords = 'offset points', 
                               va = 'center', 
                               color = text_color,
                               size = fontsize,
                               bbox=dict(boxstyle = bbox_style,
                                         fc = bbox_color, 
                                         ec = bbox_edgecolor,
                                         lw = bbox_linewidth))

                # mean value annotation and dotted line
                meanvalue = tseries_array[trace][~np.isnan(tseries_array)[trace]].mean() 
                mean_label = r'$\mu$=%.3f%s' % (meanvalue, units if units is not None else '')
                mean_position = (ntsteps - 1, meanvalue)
                mean_text_offset = (11, 0)
                ax_ts.annotate(mean_label, 
                               xy = mean_position, 
                               xytext = mean_text_offset,
                               textcoords = 'offset points', 
                               va = 'center', 
                               color = text_color, 
                               size = fontsize,
                               bbox=dict(boxstyle = bbox_style,
                                         fc = bbox_color, 
                                         ec = bbox_edgecolor, 
                                         lw = bbox_linewidth))
                # dotted line
                ax_ts.axhline(y=meanvalue,
                              linewidth=.75,
                              linestyle=':',
                              color=bbox_color)
                
                        
        else :
            print('Number of measure name : ' + str(ntrace_name) + ' is different from number of traces to plot : ' + str(ntrace))
        
        
        # Cutff values
        if cutoff is not None:
            for i, thr in enumerate(cutoff):
               
                # Plot as a dotted line
                x = (0, ntsteps - 1)
                y = [thr] * 2
                ax_ts.plot(x, y,
                           linewidth = 1.,
                           linestyle = '--',
                           color = 'k')
                
                # Annotate the graph
                if cutoff_name : 
                    # Compute offset position (avoid conflict with other annotation boxes)
                    y_off = [0.0, 0.0]
                    for pth in cutoff[:i]:
                        inc = abs(thr - pth)
                        if inc < yrange:
                            factor = (- (inc / yrange) + 1) ** 2
                            if (thr - pth) < 0.0:
                                y_off[0] -= factor * 20
                            else:
                                y_off[1] += factor * 20
                    offset = y_off[0] if abs(y_off[0]) > y_off[1] else y_off[1]

                    # Annotation parameters
                    label = '{:.5f}{}'.format(thr, units if units is not None else '')
                    position = (ntsteps - 1, thr)
                    text_offset = (11, offset)
                    bbox_color = 'dimgray'
                    arrow_color = bbox_color
                    arrow_relativepos = (.1, .5)

                    # Annotate
                    ax_ts.annotate(label, 
                                   xy = position, 
                                   xytext = text_offset,
                                   textcoords = 'offset points', 
                                   va = 'center',
                                   color = text_color, 
                                   size = fontsize,
                                   bbox=dict(boxstyle = bbox_style, 
                                             fc = bbox_color, 
                                             ec = bbox_edgecolor, 
                                             lw = bbox_linewidth),
                                   arrowprops=dict(arrowstyle = arrow_style, 
                                                   lw = arrow_linewidth, 
                                                   patchA = None, 
                                                   patchB = None,
                                                   fc = arrow_color, 
                                                   ec = arrow_edgecolor, 
                                                   relpos = arrow_relativepos))

        # Add any other annotation if provided:
        if annot_d_list is not None:
            for annot_d in annot_d_list:
                ax_ts.annotate(**annot_d)
            
            
        if not gs_dist is None:
            ax_dist = plt.subplot(gs_dist)
            sns.displot(tseries, vertical=True, ax=ax_dist)
            ax_dist.set_xlabel('Timesteps')
            ax_dist.set_ylim(ax_ts.get_ylim())
            ax_dist.set_yticklabels([])
            return [ax_ts, ax_dist], gs
        else:
            return ax_ts, gs
