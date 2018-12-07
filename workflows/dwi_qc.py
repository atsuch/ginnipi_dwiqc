# -*- coding: utf-8 -*
'''
QC pipe for DWI (and b0)

@author:    Ami Tsuchida 2018


'''
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.base import Undefined
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.mrtrix3.utils import DWIExtract
from nipype.interfaces.fsl.utils import Split
from nipype.interfaces.fsl.maths import ApplyMask
from nipype.interfaces.afni.preprocess import OutlierCount
from nipype.algorithms.confounds import TSNR
from ginnipi.interfaces.custom import SagQCplot, SNR_QCplot, MaskOverlayQCplot
from ginnipi.toolbox.flow import (getElementFromList, writeIdentity,
                                  makeFStringElementFromFnameList)
from ginnipi.toolbox.computations import compute_voxel_percentage
from ginnipi.toolbox.plotting_tools import plot_img_intensity_distribution
from ginnipi.toolbox.assessors import (createDWIshellNS_QCxml,
                                       createDWIshell_QCxml,
                                       createDTI_QCxml)
from ginnipi.toolbox.qualitycontrol_dwi import make_dwi_qc_traceplot


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


def genDWIqcPipeline(name,
                     numthreads=10,
                     bval_list=[0, 300.0, 1000.0, 2000.0],
                     meta=False):

    # Nonb0 bval list
    Nonb0_bval_list = [val for val in bval_list if val != 0]

    # Workflow
    dwiqc = Workflow(name)
    inputNode=Node(IdentityInterface(fields=[
                                             "dwi_merge",
                                             "preprocessed_dwi",
                                             "bvals",
                                             "bvecs",
                                             "eddy_cnr_maps",
                                             "eddy_outlier_map",
                                             "eddy_outlier_n_stdev_map",
                                             "eddy_outlier_n_sqr_stdev_map",
                                             "eddy_movement_rms",
                                             "eddy_restricted_movement_rms",
                                             "eddy_corrected_avg_b0",
                                             "dwi_b0_mask",
                                             "dwi_cropped_b0_mask",
                                             "single_dti_residual",
                                             "single_dti_implausible_mask",
                                             "multi_dti_residual",
                                             "multi_dti_implausible_mask"]),
                                             name='inputNode')

    ### DWI QC pipeline ###

    ###########################################################################
    ##
    ## In this pipeline, we will summarize;
    ##
    ## 1) Motion rms output from Eddy (restricted and non-restricted)
    ##
    ## 2) Slice-/voxel-level outlier detection by Eddy /afni 3dOutCount
    ##
    ## 3) Simple mosaic plots of saggital cuts (inspired by Eddy)
    ##
    ## 4) CNR/tSNR maps from Eddy/ nipype tSNR function
    ##
    ## 5) Scil DTI QC outputs summary
    ##      5-1) Mosaic plots of dti residual images and & implausible voxels
    ##      5-2) Summary statistics of these images
    ##
    ## Framewise plots of 1 & 2 are summarized with make_dwi_qc_traceplot.
    ## The summary statistics are aggregated to create dwi assessor xml.
    ##
    ###########################################################################


    # Get bval per frame so that we can use for extracing shell-specific dwi
    # and/or summarizing Eddy QC output by shell val.
    dwiGetbvals = Node(Function(input_names=['bvalsFilename', 'bval_list', 'delimiter', 'out_file'],
                                output_names=['out_file'],
                                function=get_bval_per_frame),
                       name="dwiGetbvals")
    dwiGetbvals.inputs.bval_list = bval_list
    dwiqc.connect(inputNode, "bvals", dwiGetbvals, "bvalsFilename")

    ### Extract pre- and post-processing b0 as well as non-b0 shells ##########

    # Create skull-stripped DWI (pre-Eddy/Denoise)
    dwiApplyMask = Node(ApplyMask(), name="dwiApplyMask")
    dwiApplyMask.inputs.out_file = "dwi_brain.nii.gz"
    dwiApplyMask.inputs.output_type = "NIFTI_GZ"
    dwiqc.connect(inputNode, "dwi_merge", dwiApplyMask, "in_file")
    dwiqc.connect(inputNode, "dwi_b0_mask", dwiApplyMask, "mask_file")

    # Node to extract b=0 frame pre-Eddy/denoising
    dwiGetAllb0Pre = Node(DWIExtract(), name='dwiGetAllb0Pre')
    dwiGetAllb0Pre.inputs.bzero = True
    dwiGetAllb0Pre.inputs.nthreads = numthreads
    dwiGetAllb0Pre.inputs.out_file = 'dwi_b0_pre.nii.gz'
    dwiqc.connect(dwiApplyMask, "out_file", dwiGetAllb0Pre, "in_file")
    dwiqc.connect(inputNode, "bvals", dwiGetAllb0Pre, "in_bval")
    dwiqc.connect(inputNode, "bvecs", dwiGetAllb0Pre, "in_bvec")

    # Node to extract b>0 frames for each bval separately
    dwiGetNonb0Pre = MapNode(DWIExtract(),
                             iterfield=["shell", "out_file"],
                             name='dwiGetNonb0Pre')
    dwiGetNonb0Pre.inputs.nthreads = numthreads
    dwiGetNonb0Pre.inputs.singleshell = True
    dwiGetNonb0Pre.inputs.shell = [[val] for val in Nonb0_bval_list]
    dwiGetNonb0Pre.inputs.out_file = ['dwi_b{}_pre.nii.gz'.format(str(int(val))) for val in Nonb0_bval_list]
    dwiqc.connect(dwiApplyMask, "out_file", dwiGetNonb0Pre, "in_file")
    dwiqc.connect(inputNode, "bvals", dwiGetNonb0Pre, "in_bval")
    dwiqc.connect(inputNode, "bvecs", dwiGetNonb0Pre, "in_bvec")

    # Node to extract b=0 frame post-Eddy/denoising
    dwiGetAllb0Post = Node(DWIExtract(), name='dwiGetAllb0Post')
    dwiGetAllb0Post.inputs.bzero = True
    dwiGetAllb0Post.inputs.nthreads = numthreads
    dwiGetAllb0Post.inputs.out_file = 'dwi_b0_post.nii.gz'
    dwiqc.connect(inputNode, "preprocessed_dwi", dwiGetAllb0Post, "in_file")
    dwiqc.connect(inputNode, "bvals", dwiGetAllb0Post, "in_bval")
    dwiqc.connect(inputNode, "bvecs", dwiGetAllb0Post, "in_bvec")


    # Node to extract b>0 frames for each bval separately for post Eddy/denoising
    dwiGetNonb0Post = MapNode(DWIExtract(),
                              iterfield=["shell", "out_file"],
                              name='dwiGetNonb0Post')
    dwiGetNonb0Post.inputs.nthreads = numthreads
    dwiGetNonb0Post.inputs.singleshell = True
    dwiGetNonb0Post.inputs.shell = [[val] for val in Nonb0_bval_list]
    dwiGetNonb0Post.inputs.out_file = ['dwi_b{}_post.nii.gz'.format(str(int(val))) for val in Nonb0_bval_list]
    dwiqc.connect(inputNode, "preprocessed_dwi", dwiGetNonb0Post, "in_file")
    dwiqc.connect(inputNode, "bvals", dwiGetNonb0Post, "in_bval")
    dwiqc.connect(inputNode, "bvecs", dwiGetNonb0Post, "in_bvec")

    ###########################################################################
    
    ### b0 mask QC to make sure it's a good mask ##############################
    dwiB0maskQC = Node(MaskOverlayQCplot(), name="dwiB0maskQC")
    dwiB0maskQC.inputs.transparency = 1
    dwiB0maskQC.inputs.bg_max = 99.99
    dwiB0maskQC.inputs.out_file = "dwiB0mask_overlay.png"
    dwiqc.connect(inputNode, "eddy_corrected_avg_b0", dwiB0maskQC, "bg_im_file")
    dwiqc.connect(inputNode, "dwi_b0_mask", dwiB0maskQC, "mask_file")

    ###########################################################################

    ### Nodes to get outlier counts using 3dToutcount #########################
    # b=0 Pre-Eddy/Denoise
    dwiCountOutb0Pre = Node(OutlierCount(), name="dwiCountOutb0Pre")
    dwiCountOutb0Pre.inputs.autoclip = Undefined
    dwiCountOutb0Pre.inputs.automask = Undefined
    dwiCountOutb0Pre.inputs.fraction = True
    dwiCountOutb0Pre.inputs.out_file = "dwi_b0_pre_afni_outliers.txt"
    dwiqc.connect(dwiGetAllb0Pre, "out_file", dwiCountOutb0Pre, "in_file")
    dwiqc.connect(inputNode, "dwi_b0_mask", dwiCountOutb0Pre, "mask")

    # b>0 Pre-Eddy/Denoise
    dwiCountOutNonb0Pre = MapNode(OutlierCount(),
                                  iterfield=["in_file", "out_file"],
                                  name="dwiCountOutNonb0Pre")
    dwiCountOutNonb0Pre.inputs.autoclip = Undefined
    dwiCountOutNonb0Pre.inputs.automask = Undefined
    dwiCountOutNonb0Pre.inputs.fraction = True
    dwiqc.connect(dwiGetNonb0Pre, "out_file", dwiCountOutNonb0Pre, "in_file")
    dwiqc.connect(inputNode, "dwi_b0_mask", dwiCountOutNonb0Pre, "mask")
    dwiqc.connect(dwiGetNonb0Pre,
                  ("out_file", makeFStringElementFromFnameList,
                   ".nii.gz", "_afni_outliers.txt"),
                  dwiCountOutNonb0Pre, "out_file")

    # b=0 Post-Eddy/Denoise
    dwiCountOutb0Post = Node(OutlierCount(), name="dwiCountOutb0Post")
    dwiCountOutb0Post.inputs.autoclip = Undefined
    dwiCountOutb0Post.inputs.automask = Undefined
    dwiCountOutb0Post.inputs.fraction = True
    dwiCountOutb0Post.inputs.out_file = "dwi_b0_post_afni_outliers.txt"
    dwiqc.connect(dwiGetAllb0Post, "out_file", dwiCountOutb0Post, "in_file")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask", dwiCountOutb0Post, "mask")

    # b>0 Post-Eddy/Denoise
    dwiCountOutNonb0Post = MapNode(OutlierCount(),
                                   iterfield=["in_file", "out_file"],
                                   name="dwiCountOutNonb0Post")
    dwiCountOutNonb0Post.inputs.autoclip = Undefined
    dwiCountOutNonb0Post.inputs.automask = Undefined
    dwiCountOutNonb0Post.inputs.fraction = True
    dwiqc.connect(dwiGetNonb0Post, "out_file", dwiCountOutNonb0Post, "in_file")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask",
                  dwiCountOutNonb0Post, "mask")
    dwiqc.connect(dwiGetNonb0Post,
                  ("out_file", makeFStringElementFromFnameList,
                   ".nii.gz", "_afni_outliers.txt"),
                  dwiCountOutNonb0Post, "out_file")

    ###########################################################################


    ### 1) and 2) Pass Eddy QC related output and afni outlier     ############
    ### output to make_dwi_qc_traceplot to get combined traceplot  ###########
    makeTraceplot = Node(Function(input_names=['eddy_movement_rms',
                                               'eddy_restricted_movement_rms',
                                               'eddy_outlier_map',
                                               'eddy_outlier_n_stdev_map',
                                               'eddy_outlier_n_sqr_stdev_map',
                                               'afni_outliers_b0pre',
                                               'afni_outliers_nonb0pre_list',
                                               'afni_outliers_b0post',
                                               'afni_outliers_nonb0post_list',
                                               'bval_per_frame',
                                               'nonb0_bval_list'],
                                  output_names=['out_plot'],
                                  function=make_dwi_qc_traceplot),
                         name="makeTraceplot")
    makeTraceplot.inputs.nonb0_bval_list = Nonb0_bval_list
    dwiqc.connect(inputNode, "eddy_movement_rms",
                  makeTraceplot, "eddy_movement_rms")
    dwiqc.connect(inputNode, "eddy_restricted_movement_rms",
                  makeTraceplot, "eddy_restricted_movement_rms")
    dwiqc.connect(inputNode, "eddy_outlier_map", makeTraceplot,
                  "eddy_outlier_map")
    dwiqc.connect(inputNode, "eddy_outlier_n_stdev_map",
                  makeTraceplot, "eddy_outlier_n_stdev_map")
    dwiqc.connect(inputNode, "eddy_outlier_n_sqr_stdev_map",
                  makeTraceplot, "eddy_outlier_n_sqr_stdev_map")
    dwiqc.connect(dwiCountOutb0Pre, "out_file",
                  makeTraceplot, "afni_outliers_b0pre")
    dwiqc.connect(dwiCountOutNonb0Pre, "out_file",
                  makeTraceplot, "afni_outliers_nonb0pre_list")
    dwiqc.connect(dwiCountOutb0Post, "out_file",
                  makeTraceplot, "afni_outliers_b0post")
    dwiqc.connect(dwiCountOutNonb0Post, "out_file",
                  makeTraceplot, "afni_outliers_nonb0post_list")
    dwiqc.connect(dwiGetbvals, "out_file",
                  makeTraceplot, "bval_per_frame")

    ### Also pass eddy motion and outlier info to DWI shell NS QC xml
    createShellNSAssessor = Node(Function(input_names=['eddy_movement_rms',
                                                       'eddy_restricted_movement_rms',
                                                       'eddy_outlier_map',
                                                       'eddy_outlier_n_stdev_map',
                                                       'eddy_outlier_n_sqr_stdev_map'],
                                          output_names=['xnat_assessor'],
                                          function=createDWIshellNS_QCxml),
                                 name="createShellNSAssessor")
    dwiqc.connect(inputNode, "eddy_movement_rms",
                  createShellNSAssessor, "eddy_movement_rms")
    dwiqc.connect(inputNode, "eddy_restricted_movement_rms",
                  createShellNSAssessor, "eddy_restricted_movement_rms")
    dwiqc.connect(inputNode, "eddy_outlier_map",
                  createShellNSAssessor, "eddy_outlier_map")
    dwiqc.connect(inputNode, "eddy_outlier_n_stdev_map",
                  createShellNSAssessor, "eddy_outlier_n_stdev_map")
    dwiqc.connect(inputNode, "eddy_outlier_n_sqr_stdev_map",
                  createShellNSAssessor, "eddy_outlier_n_sqr_stdev_map")

    ###########################################################################


    ### 3) Simple saggital plots of post-Eddy/Denoise DWI images per shell ####
    # Post Eddy/Denoise b0
    plotSagb0 = Node(SagQCplot(), name="plotSagb0")
    plotSagb0.inputs.out_file = "dwi_b0_post_sagplot.png"
    dwiqc.connect(dwiGetAllb0Post, "out_file", plotSagb0, "im_file")

    # Post Eddy/Denoise Nonb0
    plotSagNonb0 = MapNode(SagQCplot(),
                           iterfield=["im_file", "out_file"],
                           name="plotSagNonb0")
    plotSagNonb0.inputs.out_file = "postproc_dwi_b0_sagplot.png"
    dwiqc.connect(dwiGetNonb0Post, "out_file", plotSagNonb0, "im_file")
    dwiqc.connect(dwiGetNonb0Post,
                  ("out_file", makeFStringElementFromFnameList,
                   ".nii.gz", "_sagplot.png"),
                  plotSagNonb0, "out_file")

    ### 4) QC plots and summary from CNR/tSNR maps from Eddy and  #############
    ### Nipype TSNR function                                      #############
    # Node to split CNR map to get tSNR for b=0 and CNR for b>0 shells
    dwiSplitCNR = Node(Split(), name="dwiSplitCNR")
    dwiSplitCNR.inputs.dimension = 't'
    dwiSplitCNR.inputs.out_base_name = 'eddy_cnr'
    dwiqc.connect(inputNode, "eddy_cnr_maps", dwiSplitCNR, "in_file")

    # Node to get tSNR using nipype TSNR for b>0
    dwiNonb0TSNR = MapNode(TSNR(),
                           iterfield=["in_file", "tsnr_file"],
                           name="dwiNonb0TSNR")
    dwiqc.connect(dwiGetNonb0Post, "out_file", dwiNonb0TSNR, "in_file")
    dwiqc.connect(dwiGetNonb0Post,
                  ("out_file", makeFStringElementFromFnameList,
                   ".nii", "_tsnr.nii"),
                  dwiNonb0TSNR, "tsnr_file")

    # Plot eddy tSNR (for b=0) and CNR maps (for b>0)
    # Mosaic plots
    plotCNRmaps = MapNode(SNR_QCplot(),
                          iterfield=["snr_im_file", "out_file"],
                          name="plotCNRmaps")
    plotCNRmaps.inputs.out_file = ["eddy_cnr_b{}.png".format(str(int(bval))) for bval in bval_list]
    dwiqc.connect(inputNode, "dwi_b0_mask", plotCNRmaps, "bg_im_file")
    dwiqc.connect(inputNode, "dwi_b0_mask", plotCNRmaps, "mask_file")
    dwiqc.connect(dwiSplitCNR, "out_files", plotCNRmaps, "snr_im_file")

    # Image intensity histograms
    plotCNRmapHist = MapNode(Function(input_names=['img',
                                                   'mask',
                                                   'title',
                                                   'out_plot',
                                                   'out_stats'],
                                      output_names=['out_plot', 'out_stats'],
                                      function=plot_img_intensity_distribution),
                             iterfield=["img", "title", "out_plot", "out_stats"],
                             name="plotCNRmapHist")
    plotCNRmapHist.inputs.title = ["Image intensity distribution of Eddy CNR map: b={}".format(str(bval)) for bval in bval_list]
    dwiqc.connect(dwiSplitCNR, "out_files", plotCNRmapHist, "img")
    dwiqc.connect(inputNode, "dwi_b0_mask", plotCNRmapHist, "mask")
    dwiqc.connect(plotCNRmaps,
                  ("out_file", makeFStringElementFromFnameList,
                   ".png", "_hist.png"),
                  plotCNRmapHist, "out_plot")
    dwiqc.connect(plotCNRmaps,
                  ("out_file", makeFStringElementFromFnameList,
                   ".png", "_distribution_stats.csv"),
                  plotCNRmapHist, "out_stats")


    # Plot nipype tSNR maps for b>0
    # Mosaic plots
    plotTSNRmaps = MapNode(SNR_QCplot(),
                           iterfield=["snr_im_file", "out_file"],
                           name="plotTSNRmaps")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask", plotTSNRmaps, "bg_im_file")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask", plotTSNRmaps, "mask_file")
    dwiqc.connect(dwiNonb0TSNR, "tsnr_file", plotTSNRmaps, "snr_im_file")
    dwiqc.connect(dwiGetNonb0Post,
                  ("out_file", makeFStringElementFromFnameList,
                   ".nii.gz", "_tsnr.png"),
                  plotTSNRmaps, "out_file")

    # Image intensity histograms
    plotTSNRmapHist = MapNode(Function(input_names=['img',
                                                    'mask',
                                                    'title',
                                                    'out_plot',
                                                    'out_stats'],
                                        output_names=['out_plot', 'out_stats'],
                                        function=plot_img_intensity_distribution),
                              iterfield=["img", "title", "out_plot", "out_stats"],
                              name="plotTSNRmapHist")
    plotTSNRmapHist.inputs.title = ["Image intensity distribution of TSNR map: b={}".format(str(bval)) for bval in Nonb0_bval_list]
    dwiqc.connect(dwiNonb0TSNR, "tsnr_file", plotTSNRmapHist, "img")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask", plotTSNRmapHist, "mask")
    dwiqc.connect(plotTSNRmaps,
                  ("out_file", makeFStringElementFromFnameList,
                   ".png", "_hist.png"),
                  plotTSNRmapHist, "out_plot")
    dwiqc.connect(plotTSNRmaps,
                  ("out_file", makeFStringElementFromFnameList,
                   ".png", "_distribution_stats.csv"),
                  plotTSNRmapHist, "out_stats")

    ### Now pass shell-specific QC metrics to b0/nonb0 QC assessor
    # b0
    createDWIb0Assessor = Node(Function(input_names=['bval',
                                                     'afni_shell_specific_outlier_pre',
                                                     'afni_shell_specific_outlier_post',
                                                     'eddy_cnr_distribution_stats'],
                                        output_names=['xnat_assessor'],
                                        function=createDWIshell_QCxml),
                               name="createDWIb0Assessor")
    createDWIb0Assessor.inputs.bval = 0.0
    dwiqc.connect(dwiCountOutb0Pre, "out_file",
                  createDWIb0Assessor, "afni_shell_specific_outlier_pre")
    dwiqc.connect(dwiCountOutb0Post, "out_file",
                  createDWIb0Assessor, "afni_shell_specific_outlier_post")
    dwiqc.connect(plotCNRmapHist, ("out_stats", getElementFromList, 0),
                  createDWIb0Assessor, "eddy_cnr_distribution_stats")

    # b > 0
    createDWInonb0Assessor = MapNode(Function(input_names=['bval',
                                                           'afni_shell_specific_outlier_pre',
                                                           'afni_shell_specific_outlier_post',
                                                           'eddy_cnr_distribution_stats',
                                                           'post_tsnr_distribution_stats'],
                                              output_names=['xnat_assessor'],
                                              function=createDWIshell_QCxml),
                                     iterfield=["bval",
                                                "afni_shell_specific_outlier_pre",
                                                "afni_shell_specific_outlier_post",
                                                "eddy_cnr_distribution_stats",
                                                "post_tsnr_distribution_stats"],
                                     name="createDWInonb0Assessor")
    createDWInonb0Assessor.inputs.bval = Nonb0_bval_list
    dwiqc.connect(dwiCountOutNonb0Pre, "out_file",
                  createDWInonb0Assessor, "afni_shell_specific_outlier_pre")
    dwiqc.connect(dwiCountOutNonb0Post, "out_file",
                  createDWInonb0Assessor, "afni_shell_specific_outlier_post")
    dwiqc.connect(plotCNRmapHist, ("out_stats", getElementFromList, 1, -1),
                  createDWInonb0Assessor, "eddy_cnr_distribution_stats")
    dwiqc.connect(plotTSNRmapHist, "out_stats",
                  createDWInonb0Assessor, "post_tsnr_distribution_stats")

    ###########################################################################



    ###########################################################################

    ### 5) QC image and summary of scil DTI QC outputs      ###################

    # Plot of residual image with implausible voxels
    # for multi-shell DTI
    plotDTImultiQC = Node(MaskOverlayQCplot(), name="plotDTImultiQC")
    plotDTImultiQC.inputs.transparency = 1
    plotDTImultiQC.inputs.out_file = "multiDTIresidual_with_implausible_voxels.png"
    dwiqc.connect(inputNode, "multi_dti_residual",
                  plotDTImultiQC, "bg_im_file")
    dwiqc.connect(inputNode, "multi_dti_implausible_mask",
                  plotDTImultiQC, "mask_file")

    # for single-shell DTI
    plotDTIsingleQC = Node(MaskOverlayQCplot(), name="plotDTIsingleQC")
    plotDTIsingleQC.inputs.transparency = 1
    plotDTIsingleQC.inputs.out_file = "singleDTIresidual_with_implausible_voxels.png"
    dwiqc.connect(inputNode, "single_dti_residual",
                  plotDTIsingleQC, "bg_im_file")
    dwiqc.connect(inputNode, "single_dti_implausible_mask",
                  plotDTIsingleQC, "mask_file")

    # Image intensity analysis of residual image
    # for multi-shell DTI
    plotDTIresHist_multi = Node(Function(input_names=['img',
                                                      'mask',
                                                      'title',
                                                      'out_plot',
                                                      'out_stats'],
                                         output_names=['out_plot', 'out_stats'],
                                         function=plot_img_intensity_distribution),
                                name="plotDTIresHist_multi")
    plotDTIresHist_multi.inputs.title = "Image intensity distribution of multi-shell DTI residuals"
    plotDTIresHist_multi.inputs.out_plot = "multi_dti_residual_hist.png"
    plotDTIresHist_multi.inputs.out_stats = "multi_dti_residual_distribution_stats.csv"
    dwiqc.connect(inputNode, "multi_dti_residual",
                  plotDTIresHist_multi, "img")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask",
                  plotDTIresHist_multi, "mask")

    # for single-shell DTI
    plotDTIresHist_single = Node(Function(input_names=['img',
                                                       'mask',
                                                       'title',
                                                       'out_plot',
                                                       'out_stats'],
                                         output_names=['out_plot', 'out_stats'],
                                         function=plot_img_intensity_distribution),
                                name="plotDTIresHist_single")
    plotDTIresHist_single.inputs.title = "Image intensity distribution of single-shell DTI residuals"
    plotDTIresHist_single.inputs.out_plot = "single_dti_residual_hist.png"
    plotDTIresHist_single.inputs.out_stats = "single_dti_residual_distribution_stats.csv"
    dwiqc.connect(inputNode, "single_dti_residual",
                  plotDTIresHist_single, "img")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask",
                  plotDTIresHist_single, "mask")

    # Get percent voxels of physically implausible mask
    # for multi-shell DTI
    multiDTIimplausibleVoxelStats = Node(Function(input_names=['in_file',
                                                                'mask_file',
                                                                'out_file'],
                                                  output_names=['out_file'],
                                                  function=compute_voxel_percentage),
                                         name="multiDTIimplausibleVoxelStats")
    dwiqc.connect(inputNode, "multi_dti_implausible_mask",
                  multiDTIimplausibleVoxelStats, "in_file")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask",
                  multiDTIimplausibleVoxelStats, "mask_file")

    # for single-shell DTI
    singleDTIimplausibleVoxelStats = Node(Function(input_names=['in_file',
                                                                'mask_file',
                                                                'out_file'],
                                                   output_names=['out_file'],
                                                   function=compute_voxel_percentage),
                                          name="singleDTIimplausibleVoxelStats")
    dwiqc.connect(inputNode, "single_dti_implausible_mask",
                  singleDTIimplausibleVoxelStats, "in_file")
    dwiqc.connect(inputNode, "dwi_cropped_b0_mask",
                  singleDTIimplausibleVoxelStats, "mask_file")

    ### Now pass DTI QC metrics to DTI QC assessor
    # For multi-shell DTI
    createDTImultiAssessor = Node(Function(input_names=['dti_residual_distribution_stats',
                                                        'dti_implausible_mask_stats'],
                                           output_names=['xnat_assessor'],
                                           function=createDTI_QCxml),
                                  name="createDTImultiAssessor")
    dwiqc.connect(plotDTIresHist_multi, "out_stats",
                  createDTImultiAssessor, "dti_residual_distribution_stats")
    dwiqc.connect(multiDTIimplausibleVoxelStats, "out_file",
                  createDTImultiAssessor, "dti_implausible_mask_stats")
    
    # For single-shell DTI
    createDTIsingleAssessor = Node(Function(input_names=['dti_residual_distribution_stats',
                                                         'dti_implausible_mask_stats'],
                                            output_names=['xnat_assessor'],
                                            function=createDTI_QCxml),
                                   name="createDTIsingleAssessor")
    dwiqc.connect(plotDTIresHist_single, "out_stats",
                  createDTIsingleAssessor, "dti_residual_distribution_stats")
    dwiqc.connect(singleDTIimplausibleVoxelStats, "out_file",
                  createDTIsingleAssessor, "dti_implausible_mask_stats")


    ###########################################################################

    field_list = [
                  "dwiB0mask_QCplot",
                  "dwi_global_traceplot",
                  "dwi_global_assessor",
                  "b0_sagplot",
                  "Nonb0_sagplots",
                  "b0_eddy_tsnr_map",
                  "Nonb0_eddy_cnr_maps",
                  "b0_eddy_tsnr_distribution_plot",
                  "b0_eddy_tsnr_distribution_stat",
                  "Nonb0_eddy_cnr_distribution_plots",
                  "Nonb0_eddy_cnr_distribution_stats",
                  "Nonb0_post_tsnr_maps",
                  "Nonb0_post_tsnr_distribution_plots",
                  "Nonb0_post_tsnr_distribution_stats",
                  "b0_shell_specific_assessor",
                  "Nonb0_shell_specific_assessors",
                  "multiDTI_residual_plot",
                  "singleDTI_residual_plot",
                  "multiDTI_residual_distribution_plot",
                  "multiDTI_residual_distribution_stat",
                  "singleDTI_residual_distribution_plot",
                  "singleDTI_residual_distribution_stat",
                  "multiDTI_assessor",
                  "singleDTI_assessor",
                  ]

    outputNode = Node(IdentityInterface(fields=field_list), name="outputNode")
    dwiqc.connect(dwiB0maskQC, "out_file", outputNode, "dwiB0mask_QCplot")
    dwiqc.connect(makeTraceplot, "out_plot", outputNode, "dwi_global_traceplot")
    dwiqc.connect(createShellNSAssessor, "xnat_assessor", outputNode, "dwi_global_assessor")
    dwiqc.connect(plotSagb0, "out_file", outputNode, "b0_sagplot")
    dwiqc.connect(plotSagNonb0, "out_file", outputNode, "Nonb0_sagplots")
    dwiqc.connect(plotCNRmaps, ("out_file", getElementFromList, 0),
                  outputNode, "b0_eddy_tsnr_map")
    dwiqc.connect(plotCNRmaps, ("out_file", getElementFromList, 1, -1),
                  outputNode, "Nonb0_eddy_cnr_maps")
    dwiqc.connect(plotCNRmapHist, ("out_plot", getElementFromList, 0),
                  outputNode, "b0_eddy_tsnr_distribution_plot")
    dwiqc.connect(plotCNRmapHist, ("out_stats", getElementFromList, 0),
                  outputNode, "b0_eddy_tsnr_distribution_stat")
    dwiqc.connect(plotCNRmapHist, ("out_plot", getElementFromList, 1, -1),
                  outputNode, "Nonb0_eddy_cnr_distribution_plots")
    dwiqc.connect(plotCNRmapHist, ("out_stats", getElementFromList, 1, -1),
                  outputNode, "Nonb0_eddy_cnr_distribution_stats")
    dwiqc.connect(plotTSNRmaps, "out_file", outputNode, "Nonb0_post_tsnr_maps")
    dwiqc.connect(plotTSNRmapHist, "out_plot", outputNode, "Nonb0_post_tsnr_distribution_plots")
    dwiqc.connect(plotTSNRmapHist, "out_stats", outputNode, "Nonb0_post_tsnr_distribution_stats")
    dwiqc.connect(createDWIb0Assessor, "xnat_assessor", outputNode, "b0_shell_specific_assessor")
    dwiqc.connect(createDWInonb0Assessor, "xnat_assessor", outputNode, "Nonb0_shell_specific_assessors")
    dwiqc.connect(plotDTImultiQC, "out_file", outputNode, "multiDTI_residual_plot")
    dwiqc.connect(plotDTIsingleQC, "out_file", outputNode, "singleDTI_residual_plot")
    dwiqc.connect(plotDTIresHist_multi, "out_plot", outputNode, "multiDTI_residual_distribution_plot")
    dwiqc.connect(plotDTIresHist_multi, "out_stats", outputNode, "multiDTI_residual_distribution_stat")
    dwiqc.connect(plotDTIresHist_single, "out_plot", outputNode, "singleDTI_residual_distribution_plot")
    dwiqc.connect(plotDTIresHist_single, "out_stats", outputNode, "singleDTI_residual_distribution_stat")
    dwiqc.connect(createDTImultiAssessor, "xnat_assessor", outputNode, "multiDTI_assessor")
    dwiqc.connect(createDTIsingleAssessor, "xnat_assessor", outputNode, "singleDTI_assessor")

    # Node: outputDict
    if meta:
        outputDict = Node(Function(input_names=field_list,
                               output_names=['outDict'],
                               function=writeIdentity),
                      name='outputDict')
        for field in field_list:
            dwiqc.connect(outputNode,field, outputDict, field)

    dwiqc.write_graph(dotfilename='dwiqc', graph2use='flat', format='png', simple_form=True)
    return dwiqc
