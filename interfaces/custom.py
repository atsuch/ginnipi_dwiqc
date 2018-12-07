# -*- coding: utf-8 -*-
'''
ABC pipeline interfaces and customizations of the original nipype interfaces
@author: Pierre-Yves Hervé 2015-16
'''
import warnings
import os
from shutil import copyfile
import numpy as np
from scipy.io import savemat
import subprocess
from nipype.interfaces.base import (Directory, TraitedSpec, BaseInterface,
                                    DynamicTraitedSpec, traits, Undefined, 
                                    isdefined, File, InputMultiPath, BaseInterfaceInputSpec)
from nipype.interfaces.base import (CommandLineInputSpec, CommandLine,
                                    OutputMultiPath)
from nipype.interfaces.io import IOBase, add_traits
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.interfaces.fsl.utils import ImageMaths
from nipype.utils.filemanip import (split_filename, filename_to_list,
                                    list_to_filename)

from ginnipi_dwiqc.toolbox.computations import identify_b0s
warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)
from nipype.utils.functions import getsource, create_function_from_source
import ginnipi


class AbcEddyInputSpec(FSLCommandInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        argstr='--imain=%s',
        desc=('File containing all the images to estimate '
              'distortions for'))
    in_mask = File(
        exists=True,
        mandatory=True,
        argstr='--mask=%s',
        desc='Mask to indicate brain')
    in_index = File(
        exists=True,
        mandatory=True,
        argstr='--index=%s',
        desc=('File containing indices for all volumes in --imain '
              'into --acqp and --topup'))
    in_acqp = File(
        exists=True,
        mandatory=True,
        argstr='--acqp=%s',
        desc='File containing acquisition parameters')
    in_bvec = File(
        exists=True,
        mandatory=True,
        argstr='--bvecs=%s',
        desc=('File containing the b-vectors for all volumes in '
              '--imain'))
    in_bval = File(
        exists=True,
        mandatory=True,
        argstr='--bvals=%s',
        desc=('File containing the b-values for all volumes in '
              '--imain'))
    out_base = traits.Str(
        'eddy_corrected',
        argstr='--out=%s',
        usedefault=True,
        desc=('basename for output (warped) image'))
    session = File(
        exists=True,
        argstr='--session=%s',
        desc=('File containing session indices for all volumes in '
              '--imain'))
    in_topup_fieldcoef = File(
        exists=True,
        argstr="--topup=%s",
        requires=['in_topup_movpar'],
        desc=('topup file containing the field '
              'coefficients'))
    in_topup_movpar = File(
        exists=True,
        requires=['in_topup_fieldcoef'],
        desc='topup movpar.txt file')

    mb = traits.Int(
        argstr='--mb=%s',
        desc=('Multiband factor '))

    flm = traits.Enum(
        'linear',
        'quadratic',
        'cubic',
        argstr='--flm=%s',
        desc='First level EC model')

    slm = traits.Enum(
        'none',
        'linear',
        'quadratic',
        argstr='--slm=%s',
        desc='Second level EC model')

    fep = traits.Bool(
        False, argstr='--fep', desc='Fill empty planes in x- or y-directions')

    interp = traits.Enum(
        'spline',
        'trilinear',
        argstr='--interp=%s',
        desc='Interpolation model for estimation step')

    nvoxhp = traits.Int(
        1000,
        argstr='--nvoxhp=%s',
        desc=('# of voxels used to estimate the '
              'hyperparameters'))

    fudge_factor = traits.Float(
        10.0,
        argstr='--ff=%s',
        desc=('Fudge factor for hyperparameter '
              'error variance'))

    cnr_maps = traits.Bool(
        False,
        argstr='--cnr_maps',
        desc=('Write shell-wise cnr-maps '))

    dont_sep_offs_move = traits.Bool(
        False,
        argstr='--dont_sep_offs_move',
        desc=('Do NOT attempt to separate '
              'field offset from subject '
              'movement'))

    dont_peas = traits.Bool(
        False,
        argstr='--dont_peas',
        desc="Do NOT perform a post-eddy alignment of "
        "shells")

    fwhm = traits.Float(
        desc=('FWHM for conditioning filter when estimating '
              'the parameters'),
        argstr='--fwhm=%s')

    niter = traits.Int(5, argstr='--niter=%s', desc='Number of iterations')

    method = traits.Enum(
        'jac',
        'lsr',
        argstr='--resamp=%s',
        desc=('Final resampling method (jacobian/least '
              'squares)'))
    repol = traits.Bool(
        False, argstr='--repol', desc='Detect and replace outlier slices')

    ol_type = traits.Enum(
        'sw',
        'gw',
        'both',
        argstr='--ol_type=%s',
        desc=('Type of outliers, slicewise (sw), '
             'groupwise (gw) or both (both). (default sw)'))

    num_threads = traits.Int(
        1,
        usedefault=True,
        nohash=True,
        desc="Number of openmp threads to use")
    is_shelled = traits.Bool(
        False,
        argstr='--data_is_shelled',
        desc="Override internal check to ensure that "
        "date are acquired on a set of b-value "
        "shells")
    field = traits.Str(
        argstr='--field=%s',
        desc="NonTOPUP fieldmap scaled in Hz - filename has "
        "to be provided without an extension. TOPUP is "
        "strongly recommended")
    field_mat = File(
        exists=True,
        argstr='--field_mat=%s',
        desc="Matrix that specifies the relative locations of "
        "the field specified by --field and first volume "
        "in file --imain")
    use_cuda = traits.Bool(False, desc="Run eddy using cuda gpu")


class AbcEddyOutputSpec(TraitedSpec):
    out_corrected = File(
        exists=True, desc='4D image file containing all the corrected volumes')
    out_parameter = File(
        exists=True,
        desc=('text file with parameters definining the field and'
              'movement for each scan'))
    out_rotated_bvecs = File(
        exists=True, desc='File containing rotated b-values for all volumes')
    out_movement_rms = File(
        exists=True, desc='Summary of the "total movement" in each volume')
    out_restricted_movement_rms = File(
        exists=True,
        desc=('Summary of the "total movement" in each volume '
              'disregarding translation in the PE direction'))
    out_shell_alignment_parameters = File(
        exists=True,
        desc=('File containing rigid body movement parameters '
              'between the different shells as estimated by a '
              'post-hoc mutual information based registration'))
    out_outlier_report = File(
        exists=True,
        desc=('Text-file with a plain language report on what '
              'outlier slices eddy has found'))
    out_outlier_map = File(
        exists=True,
        desc=('Text-file with a numeric matrices in ASCII indicating '
              'whether the specific slice (col) in a specific frame (row) '
              'is an outlier (1) or not (0) '))
    out_outlier_n_stdev_map = File(
        exists=True,
        desc=('Text-file with a numeric matrices in ASCII indicating '
              'how many standard deviations off the mean difference '
              'between observation and prediction for each slice (col) '
              'in each frame (row) '))
    out_outlier_n_sqr_stdev_map = File(
        exists=True,
        desc=('Text-file with a numeric matrices in ASCII indicating '
              'how many standard deviations off the square root of '
              'the mean squared difference between observation and prediction '
              'for each slice (col) in each frame (row) '))
    out_cnr_maps = File(
        exists=True,
        desc=('This is a 4D image file with N+1 volumes where N is the number '
              'of non-zero b-value shells. The first volume contains '
              'the voxelwise SNR for the b=0 shell and the remaining volumes '
              'contain the voxelwise CNR (Contrast to Noise Ratio) for the '
              'non-zero b-shells in order of ascending b-value. For example '
              'if your data consists of 5 b=0, 48 b=1000 and 64 b=2000 volumes, '
              'my_eddy_output.eddy_cnr_maps will have three volumes where the '
              'first is the SNR for the b=0 volumes, followed by CNR maps for '
              'b=1000 and b=2000. The SNR for the b=0 shell is defined as '
              'mean(b0)/std(b0). The CNR for the DWI shells is defined as '
              'std(GP)/std(res) where std is the standard deviation of the '
              'Gaussian Process (GP) predictions and std(res) is the standard '
              'deviation of the residuals (the difference between the '
              'observations and the GP predictions). Useful for qualty control.'))

class AbcEddy(FSLCommand):
    """
    Interface for new FSL eddy (ver 5.0.11), a tool for estimating and correcting eddy
    currents induced distortions. `User guide
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Eddy/UsersGuide>`_ and
    `more info regarding acqp file
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file>`_.

    A temporary fix until nipype wrapping of Eddy is updated...
    Examples
    --------

    >>> from nipype.interfaces.fsl import Eddy
    >>> eddy = Eddy()
    >>> eddy.inputs.in_file = 'epi.nii'
    >>> eddy.inputs.in_mask  = 'epi_mask.nii'
    >>> eddy.inputs.in_index = 'epi_index.txt'
    >>> eddy.inputs.in_acqp  = 'epi_acqp.txt'
    >>> eddy.inputs.in_bvec  = 'bvecs.scheme'
    >>> eddy.inputs.in_bval  = 'bvals.scheme'
    >>> eddy.inputs.use_cuda = True
    >>> eddy.cmdline # doctest: +ELLIPSIS
    'eddy_cuda --acqp=epi_acqp.txt --bvals=bvals.scheme --bvecs=bvecs.scheme \
--imain=epi.nii --index=epi_index.txt --mask=epi_mask.nii \
--out=.../eddy_corrected'
    >>> eddy.inputs.use_cuda = False
    >>> eddy.cmdline # doctest: +ELLIPSIS
    'eddy_openmp --acqp=epi_acqp.txt --bvals=bvals.scheme \
--bvecs=bvecs.scheme --imain=epi.nii --index=epi_index.txt \
--mask=epi_mask.nii --out=.../eddy_corrected'
    >>> res = eddy.run() # doctest: +SKIP

    """
    _cmd = 'eddy_openmp'
    input_spec = AbcEddyInputSpec
    output_spec = AbcEddyOutputSpec

    _num_threads = 1

    def __init__(self, **inputs):
        super(AbcEddy, self).__init__(**inputs)
        self.inputs.on_trait_change(self._num_threads_update, 'num_threads')
        if not isdefined(self.inputs.num_threads):
            self.inputs.num_threads = self._num_threads
        else:
            self._num_threads_update()
        self.inputs.on_trait_change(self._use_cuda, 'use_cuda')
        if isdefined(self.inputs.use_cuda):
            self._use_cuda()

    def _num_threads_update(self):
        self._num_threads = self.inputs.num_threads
        if not isdefined(self.inputs.num_threads):
            if 'OMP_NUM_THREADS' in self.inputs.environ:
                del self.inputs.environ['OMP_NUM_THREADS']
        else:
            self.inputs.environ['OMP_NUM_THREADS'] = str(
                self.inputs.num_threads)

    def _use_cuda(self):
        self._cmd = 'eddy_cuda' if self.inputs.use_cuda else 'eddy_openmp'

    def _run_interface(self, runtime):
        # If 'eddy_openmp' is missing, use 'eddy'
        FSLDIR = os.getenv('FSLDIR', '')
        cmd = self._cmd
        if all((FSLDIR != '', cmd == 'eddy_openmp',
                not os.path.exists(os.path.join(FSLDIR, 'bin', cmd)))):
            self._cmd = 'eddy'
        runtime = super(AbcEddy, self)._run_interface(runtime)

        # Restore command to avoid side-effects
        self._cmd = cmd
        return runtime

    def _format_arg(self, name, spec, value):
        if name == 'in_topup_fieldcoef':
            return spec.argstr % value.split('_fieldcoef')[0]
        if name == 'out_base':
            return spec.argstr % os.path.abspath(value)
        return super(AbcEddy, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_corrected'] = os.path.abspath(
            '%s.nii.gz' % self.inputs.out_base)
        outputs['out_parameter'] = os.path.abspath(
            '%s.eddy_parameters' % self.inputs.out_base)

        # File generation might depend on the version of EDDY
        out_rotated_bvecs = os.path.abspath(
            '%s.eddy_rotated_bvecs' % self.inputs.out_base)
        out_movement_rms = os.path.abspath(
            '%s.eddy_movement_rms' % self.inputs.out_base)
        out_restricted_movement_rms = os.path.abspath(
            '%s.eddy_restricted_movement_rms' % self.inputs.out_base)
        out_shell_alignment_parameters = os.path.abspath(
            '%s.eddy_post_eddy_shell_alignment_parameters' %
            self.inputs.out_base)
        out_outlier_report = os.path.abspath(
            '%s.eddy_outlier_report' % self.inputs.out_base)
        out_outlier_map = os.path.abspath(
            '%s.eddy_outlier_map' % self.inputs.out_base)
        out_outlier_n_stdev_map = os.path.abspath(
            '%s.eddy_outlier_n_stdev_map' % self.inputs.out_base)
        out_outlier_n_sqr_stdev_map = os.path.abspath(
            '%s.eddy_outlier_n_sqr_stdev_map' % self.inputs.out_base)
        out_cnr_maps = os.path.abspath(
            '%s.eddy_cnr_maps.nii.gz' % self.inputs.out_base)

        if os.path.exists(out_rotated_bvecs):
            outputs['out_rotated_bvecs'] = out_rotated_bvecs
        if os.path.exists(out_movement_rms):
            outputs['out_movement_rms'] = out_movement_rms
        if os.path.exists(out_restricted_movement_rms):
            outputs['out_restricted_movement_rms'] = \
                out_restricted_movement_rms
        if os.path.exists(out_shell_alignment_parameters):
            outputs['out_shell_alignment_parameters'] = \
                out_shell_alignment_parameters
        if os.path.exists(out_outlier_report):
            outputs['out_outlier_report'] = out_outlier_report
        if os.path.exists(out_outlier_map):
            outputs['out_outlier_map'] = out_outlier_map
        if os.path.exists(out_outlier_n_stdev_map):
            outputs['out_outlier_n_stdev_map'] = out_outlier_n_stdev_map
        if os.path.exists(out_outlier_n_sqr_stdev_map):
            outputs['out_outlier_n_sqr_stdev_map'] = out_outlier_n_sqr_stdev_map
        if os.path.exists(out_cnr_maps):
            outputs['out_cnr_maps'] = out_cnr_maps

        return outputs


class MaskOverlayQCplotInputSpec(CommandLineInputSpec):
    bg_im_file = traits.Str(mandatory=True,
                            desc ='Background ref image file (typically T1 brain)',
                            argstr='%s',
                            position=1)
    mask_file = traits.Str(mandatory=True,
                           desc='Brain mask',
                           argstr='%s',
                           position=2)
    transparency = traits.Enum(0, 1,
                              argstr='%d',
                              desc='Set transparency (0: solid, 1:transparent)',
                              mandatory=True,
                              position=3)
    out_file = traits.Str('mask_overlay.png',
                           mandatory=False,
                           desc='Output png filename',
                           argstr='%s',
                           position=4)
    bg_max = traits.Float(argstr="%.3f",
                          mandatory=False,
                          desc='Optionally specifies the bg img intensity range as a percentile',
                          position=5)

class MaskOverlayQCplotOutputSpec(TraitedSpec):
    out_file = File(exists=True)

class MaskOverlayQCplot(CommandLine):
    """
    Creates multi-slice axial plot with Slicer showing mask overlaied on
    background image.
    """
    _cmd = 'mask_overlay_QC_images.sh'
    input_spec = MaskOverlayQCplotInputSpec
    output_spec = MaskOverlayQCplotOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

    
class SNR_QCplotInputSpec(CommandLineInputSpec):
    bg_im_file = traits.Str(mandatory=True,
                            desc ='Background image file',
                            argstr='%s',
                            position=1)
    mask_file = traits.Str(mandatory=True,
                           desc='Brain mask',
                           argstr='%s',
                           position=2)
    snr_im_file = traits.Str(mandatory=True,
                            desc='SNR image to be plotted',
                            argstr='%s',
                            position=3)
    out_file = traits.Str('snr_qc_plot.png',
                          mandatory=False,
                          usedefault=True,
                          desc='Output png filename',
                          argstr='%s',
                          position=4)


class SNR_QCplotOutputSpec(TraitedSpec):
    out_file = File(exists=True)

class SNR_QCplot(CommandLine):
    """
    Axial slices of SNR images in jet color scheme.
    """
    _cmd = 'snr_overlay_QC_images.sh'
    input_spec = SNR_QCplotInputSpec
    output_spec = SNR_QCplotOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs
    
class SagQCplotInputSpec(CommandLineInputSpec):
    im_file = traits.Str(mandatory=True,
                         desc ='image file to plot',
                         argstr='%s',
                         position=1)

    out_file = traits.Str('saggital_plot.png',
                          mandatory=False,
                          usedefault=True,
                          desc='Output png filename',
                          argstr='%s',
                          position=2)


class SagQCplotOutputSpec(TraitedSpec):
    out_file = File(exists=True)

class SagQCplot(CommandLine):
    """
    Plots mid-saggital slices per frame for viewing slice-to-volume
    movement artefacts in DWI data.
    """
    _cmd = 'plot_saggital.sh'
    input_spec = SagQCplotInputSpec
    output_spec = SagQCplotOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs
