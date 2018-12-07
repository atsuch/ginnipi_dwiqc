# -*- coding: utf-8 -*-
"""
Created on Sat June 1, 2018
Functions for creating assessors.
(Refactored from some functions originally included in computations.py)

@author: tsuchida
"""


def createDiffusionAssessorXml(atlas, coarseFS_files=None, detailedFS_files=None, JHU_files=None):
    '''
    Convert a csv file into xml format
    '''
    import lxml.etree as et
    import os
    from time import gmtime, strftime
    import csv

    coarse_fs_labels = ['left_wm','right_wm','corpus_callosum','left_cerebellum','right_cerebellum','left_ventral_diencephalon','right_ventral_diencephalon',
                        'left_lat_ventricle','right_lat_ventricle','brainstem']

    detailed_fs_labels = ['CC_Posterior','CC_Mid_Posterior','CC_Central','CC_Mid_Anterior','CC_Anterior','wm_lh_bankssts','wm_lh_caudalanteriorcingulate','wm_lh_caudalmiddlefrontal',
                          'wm_lh_cuneus','wm_lh_entorhinal','wm_lh_fusiform','wm_lh_inferiorparietal','wm_lh_inferiortemporal','wm_lh_isthmuscingulate','wm_lh_lateraloccipital','wm_lh_lateralorbitofrontal',
                          'wm_lh_lingual','wm_lh_medialorbitofrontal','wm_lh_middletemporal','wm_lh_parahippocampal','wm_lh_paracentral','wm_lh_parsopercularis','wm_lh_parsorbitalis','wm_lh_parstriangularis',
                          'wm_lh_pericalcarine','wm_lh_postcentral','wm_lh_posteriorcingulate','wm_lh_precentral','wm_lh_precuneus','wm_lh_rostralanteriorcingulate','wm_lh_rostralmiddlefrontal',
                          'wm_lh_superiorfrontal','wm_lh_superiorparietal','wm_lh_superiortemporal','wm_lh_supramarginal','wm_lh_frontalpole','wm_lh_temporalpole','wm_lh_transversetemporal','wm_lh_insula',
                          'wm_rh_bankssts','wm_rh_caudalanteriorcingulate','wm_rh_caudalmiddlefrontal','wm_rh_cuneus','wm_rh_entorhinal','wm_rh_fusiform','wm_rh_inferiorparietal','wm_rh_inferiortemporal',
                          'wm_rh_isthmuscingulate','wm_rh_lateraloccipital','wm_rh_lateralorbitofrontal','wm_rh_lingual','wm_rh_medialorbitofrontal','wm_rh_middletemporal','wm_rh_parahippocampal',
                          'wm_rh_paracentral','wm_rh_parsopercularis','wm_rh_parsorbitalis','wm_rh_parstriangularis','wm_rh_pericalcarine','wm_rh_postcentral','wm_rh_posteriorcingulate','wm_rh_precentral',
                          'wm_rh_precuneus','wm_rh_rostralanteriorcingulate','wm_rh_rostralmiddlefrontal','wm_rh_superiorfrontal','wm_rh_superiorparietal','wm_rh_superiortemporal','wm_rh_supramarginal',
                          'wm_rh_frontalpole','wm_rh_temporalpole','wm_rh_transversetemporal','wm_rh_insula','Left_UnsegmentedWhiteMatter','Right_UnsegmentedWhiteMatter']

    jhu_labels = ['middle_cerebellar_peduncle','pontine_crossing_tract','genu_cc','body_cc','splenium_cc','fornix','corticospinal_tract_R','corticospinal_tract_L',
                  'medial_lemniscus_R','medial_lemniscus_L','inf_cerebellar_peduncle_R','inf_cerebellar_peduncle_L','sup_cerebellar_peduncle_R','sup_cerebellar_peduncle_L',
                  'cerebral_peduncle_R','cerebral_peduncle_L','ant_limb_internal_capsule_R','ant_limb_internal_capsule_L','post_limb_internal_capsule_R','post_limb_internal_capsule_L',
                  'retrolenticular_part_internal_capsule_R','retrolenticular_part_internal_capsule_L','ant_corona_radiata_R','ant_corona_radiata_L','sup_corona_radiata_R',
                  'sup_corona_radiata_L','post_corona_radiata_R','post_corona_radiata_L','post_thalamic_radiation_R','post_thalamic_radiation_L','sagittal_stratum_R','sagittal_stratum_L',
                  'external_capsule_R','external_capsule_L','cingulum_cingulate_gyrus_R','cingulum_cingulate_gyrus_L','cingulum_hippocampus_R','cingulum_hippocampus_L',
                  'fornix_cres_or_stria_terminalis_R','fornix_cres_or_stria_terminalis_L','sup_longitudinal_fasciculus_R','sup_longitudinal_fasciculus_L','sup_fronto_occipital_fasciculus_R',
                  'sup_fronto_occipital_fasciculus_L','uncinate_fasciculus_R','uncinate_fasciculus_L','tapetum_R','tapetum_L']

    string = '<diffusion:{}Data xmlns:diffusion="http://nrg.wustl.edu/diffusion" xmlns:xnat="http://nrg.wustl.edu/xnat" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"></diffusion:{}Data>'.format(atlas,atlas)

    Assess = et.ElementTree(et.XML(string))

    namespaces={'diffusion':'http://nrg.wustl.edu/diffusion',
                'xnat':'http://nrg.wustl.edu/xnat',
                'xsi':'http://www.w3.org/2001/XMLSchema-instance'}

    if atlas == "fs" or atlas =="fs6":

        stats = ['FA','MD','AD','RD','vol']

        for fil in coarseFS_files:
            f = open(fil)
            f_read = csv.reader(f)
            data = []

            stat = stats[coarseFS_files.index(fil)]

            for row in f_read:
                data.append(row)

            f.close()

            for region in range(1,len(data[0])):
                field = et.Element(et.QName(namespaces["diffusion"],"_".join([stat,coarse_fs_labels[region-1]])),nsmap=namespaces)
                field.text = data[1][region]
                Assess.getroot().append(field)

        for fil in detailedFS_files:
            f = open(fil)
            f_read = csv.reader(f)
            data = []

            stat = stats[detailedFS_files.index(fil)]

            for row in f_read:
                data.append(row)

            f.close()

            for region in range(1,len(data[0])):
                field = et.Element(et.QName(namespaces["diffusion"],"_".join([stat,detailed_fs_labels[region-1]])),nsmap=namespaces)
                field.text = data[1][region]
                Assess.getroot().append(field)

    else:

        stats = ['FA','MD','AD','RD','wFA','wMD','wAD','wRD','vol']

        for fil in JHU_files:
            f = open(fil)
            f_read = csv.reader(f)
            data = []

            stat = stats[JHU_files.index(fil)]

            for row in f_read:
                data.append(row)

            f.close()

            for region in range(1,len(data[0])):
                field = et.Element(et.QName(namespaces["diffusion"],"_".join([stat,jhu_labels[region-1]])),nsmap=namespaces)
                field.text = data[1][region]
                Assess.getroot().append(field)

    date = et.Element(et.QName(namespaces["diffusion"],"date"),nsmap=namespaces)
    date.text = strftime("%Y-%m-%d", gmtime())
    Assess.getroot().append(date)

    # done generating XML, write it...
    fname=os.path.abspath('{}Data.xml'.format(atlas))
    with open(fname,'w') as f:
        Assess.write(f,pretty_print=True,encoding='utf-8', xml_declaration=True)
    f.close()
    return fname


def createDWIshellNS_QCxml(eddy_movement_rms,
                           eddy_restricted_movement_rms,
                           eddy_outlier_map,
                           eddy_outlier_n_stdev_map,
                           eddy_outlier_n_sqr_stdev_map):
    ''' 
    Write shell non-specific QC measures for DWI in a XML file (assessor), mainly
    from Eddy output.
    
    :param eddy_movement_rms: path to eddy motion rms file
    :param eddy_restricted_movement_rms: path to eddy restricted motion rms file
    :param eddy_outlier: path to eddy outlier file indicating outlier slices per frame
    :param eddy_outlier_n_stdev: path to eddy outlier stdev file
    :param eddy_outlier_n_sqr_stdev: path to eddy outlier sqr stdev file 
      
    '''

    import lxml.etree as et
    import os.path as op
    import numpy as np
    
    ## Load all the input files
    ############################################
    rms_dat = np.loadtxt(eddy_movement_rms)[:, 1] if op.exists(eddy_movement_rms) else None
    res_rms_dat = np.loadtxt(eddy_restricted_movement_rms)[:, 1] if op.exists(eddy_restricted_movement_rms) else None
    
    outlier_dat = np.loadtxt(eddy_outlier_map, skiprows=1) if op.exists(eddy_outlier_map) else None
    outlier_stdev_dat = np.loadtxt(eddy_outlier_n_stdev_map, skiprows=1) if op.exists(eddy_outlier_n_stdev_map) else None
    outlier_sqr_stdev_dat = np.loadtxt(eddy_outlier_n_sqr_stdev_map, skiprows=1) if op.exists(eddy_outlier_n_sqr_stdev_map) else None

    ## Compute QC values
    ############################################
    # mean RMS and mean restricted RMS
    mean_relative_rms = rms_dat.mean()
    mean_relative_res_rms = res_rms_dat.mean()

    # mean and max outlier slices
    mean_eddy_outlier_slices = outlier_dat.sum(axis=1).mean()
    max_eddy_outlier_slices = outlier_dat.sum(axis=1).max()

    # mean eddy outlier stdev and sqr_stdev
    mean_eddy_outlier_stdev = outlier_stdev_dat.mean()
    mean_eddy_outlier_sqr_stdev = outlier_sqr_stdev_dat.mean()

    ## Package into an assessor
    ############################################
    string = '<rest:restAssessorData xmlns:rest="http://nrg.wustl.edu/rest" xmlns:prov="http://www.nbirn.net/prov" xmlns:xnat="http://nrg.wustl.edu/xnat" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"></rest:restAssessorData>'
        
    Assess = et.ElementTree(et.XML(string))
        
    namespaces={'dwi':'http://nrg.wustl.edu/dwi',
                'xnat':'http://nrg.wustl.edu/xnat',
                'xsi':'http://www.w3.org/2001/XMLSchema-instance'}

    relative_rms = et.Element(et.QName(namespaces["dwi"], "QC_DWI_relative_rms"), nsmap=namespaces)
    relative_rms.text = str(mean_relative_rms)
    Assess.getroot().append(relative_rms)

    relative_res_rms = et.Element(et.QName(namespaces["dwi"], "QC_DWI_relative_restricted_rms"), nsmap=namespaces)
    relative_res_rms.text = str(mean_relative_res_rms)
    Assess.getroot().append(relative_res_rms)

    eddy_avg_ol_slices = et.Element(et.QName(namespaces["dwi"], "QC_eddy_mean_outlier_slices"), nsmap=namespaces)
    eddy_avg_ol_slices.text = str(mean_eddy_outlier_slices)
    Assess.getroot().append(eddy_avg_ol_slices)

    eddy_max_ol_slices = et.Element(et.QName(namespaces["dwi"], "QC_eddy_max_outlier_slices"), nsmap=namespaces)
    eddy_max_ol_slices.text = str(max_eddy_outlier_slices)
    Assess.getroot().append(eddy_max_ol_slices)

    eddy_avg_ol_std = et.Element(et.QName(namespaces["dwi"], "QC_eddy_mean_outlier_stdev"), nsmap=namespaces)
    eddy_avg_ol_std.text = str(mean_eddy_outlier_stdev)
    Assess.getroot().append(eddy_avg_ol_std)

    eddy_avg_ol_sqr_std = et.Element(et.QName(namespaces["dwi"], "QC_eddy_mean_outlier_sqr_stdev"), nsmap=namespaces)
    eddy_avg_ol_sqr_std.text = str(mean_eddy_outlier_sqr_stdev)
    Assess.getroot().append(eddy_avg_ol_sqr_std)
  
    fname = op.abspath('DWIshellNSAssessor.xml')
    with open(fname, 'wb') as f:
        Assess.write(f, pretty_print=True, encoding='utf-8', xml_declaration=True)

    return fname


def createDWIshell_QCxml(bval,
                         afni_shell_specific_outlier_pre,
                         afni_shell_specific_outlier_post,
                         eddy_cnr_distribution_stats,
                         post_tsnr_distribution_stats=None):
    ''' 
    Write shell-specific QC measures for DWI in a XML file (assessor)

    :param bval: bval for the shell
    :param afni_shell_specific_outlier_pre: path to the output file from afni 
                                            3dToutcount before preprocessing
    :param afni_shell_specific_outlier_post: path to the output file from afni 
                                             3dToutcount  after preprocessing
    :param eddy_cnr_distribution_stats: path to the eddy CNR map distribution stats
                                        (this is tSNR for b=0)
    :param post_tsnr_distribution_stats: path to the nipype TSNR distribution stats
                                         after preprocessing (only for b>0)
    '''
    import lxml.etree as et
    import os.path as op
    import numpy as np
    import pandas as pd
    
    ## Load all the input files
    ############################################
    afni_out_pre_dat = np.loadtxt(afni_shell_specific_outlier_pre) if op.exists(afni_shell_specific_outlier_pre) else None
    afni_out_post_dat = np.loadtxt(afni_shell_specific_outlier_post) if op.exists(afni_shell_specific_outlier_post) else None
    
    eddy_cnr_dist_df = pd.read_csv(eddy_cnr_distribution_stats) if op.exists(eddy_cnr_distribution_stats) else None
    if bval != 0:
        assert post_tsnr_distribution_stats is not None
        post_tsnr_dist_df = pd.read_csv(post_tsnr_distribution_stats) if op.exists(post_tsnr_distribution_stats) else None
    
    ## Compute QC values
    ############################################
     
    # mean and max afni outlier voxel proportions pre and post processing
    mean_afni_outlier_voxel_prop_pre = afni_out_pre_dat.mean() if afni_out_pre_dat is not None else 'NA'
    max_afni_outlier_voxel_prop_pre = afni_out_pre_dat.max() if afni_out_pre_dat is not None else 'NA'
    mean_afni_outlier_voxel_prop_post = afni_out_post_dat.mean() if afni_out_post_dat is not None else 'NA'
    max_afni_outlier_voxel_prop_post = afni_out_post_dat.max() if afni_out_post_dat is not None else 'NA'
    
    # inverse TSNR/CNRs
    inv_cnr = 1./(eddy_cnr_dist_df['median'].values[0]) if eddy_cnr_dist_df['median'].values[0] != 0.0 else 'INF'
    if bval != 0:
        if post_tsnr_dist_df is None:
            inv_tsnr = 'NA'
        else:
            inv_tsnr = 1./(post_tsnr_dist_df['median'].values[0]) if post_tsnr_dist_df['median'].values[0] != 0.0 else 'INF'

       
    ## Package into an assessor
    ############################################
    string = '<rest:restAssessorData xmlns:rest="http://nrg.wustl.edu/rest" xmlns:prov="http://www.nbirn.net/prov" xmlns:xnat="http://nrg.wustl.edu/xnat" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"></rest:restAssessorData>'
        
    Assess = et.ElementTree(et.XML(string))
        
    namespaces={'dwi':'http://nrg.wustl.edu/dwi',
                'xnat':'http://nrg.wustl.edu/xnat',
                'xsi':'http://www.w3.org/2001/XMLSchema-instance'}
        
    afni_avg_ol_pre = et.Element(et.QName(namespaces["dwi"], "QC_mean_afni_outlier_pre"), nsmap=namespaces)
    afni_avg_ol_pre.text = str(mean_afni_outlier_voxel_prop_pre)
    Assess.getroot().append(afni_avg_ol_pre)
        
    afni_max_ol_pre = et.Element(et.QName(namespaces["dwi"], "QC_max_afni_outlier_pre"), nsmap=namespaces)
    afni_max_ol_pre.text = str(max_afni_outlier_voxel_prop_pre)
    Assess.getroot().append(afni_max_ol_pre)
        
    afni_avg_ol_post = et.Element(et.QName(namespaces["dwi"], "QC_mean_afni_outlier_post"), nsmap=namespaces)
    afni_avg_ol_post.text = str(mean_afni_outlier_voxel_prop_post)
    Assess.getroot().append(afni_avg_ol_post)
        
    afni_max_ol_post = et.Element(et.QName(namespaces["dwi"], "QC_max_afni_outlier_post"), nsmap=namespaces)
    afni_max_ol_post.text = str(max_afni_outlier_voxel_prop_post)
    Assess.getroot().append(afni_max_ol_post)
    
    if bval == 0:
        inverse_tSNR = et.Element(et.QName(namespaces["dwi"],"QC_DWI_inverse_eddy_tSNR"), nsmap=namespaces)
        inverse_tSNR.text = str(inv_cnr)
        Assess.getroot().append(inverse_tSNR)
    
    else:
        inverse_CNR = et.Element(et.QName(namespaces["dwi"],"QC_DWI_inverse_eddy_CNR"), nsmap=namespaces)
        inverse_CNR.text = str(inv_cnr)
        Assess.getroot().append(inverse_CNR)
        
        inverse_tSNR = et.Element(et.QName(namespaces["dwi"],"QC_DWI_inverse_postproc_tSNR"), nsmap=namespaces)
        inverse_tSNR.text = str(inv_tsnr)
        Assess.getroot().append(inverse_tSNR)
        
    fname = op.abspath('dwiAssessor_b{}.xml'.format(str(int(bval))))
    with open(fname, 'wb') as f:
        Assess.write(f, pretty_print=True, encoding='utf-8', xml_declaration=True)
    
    return fname


def createDTI_QCxml(dti_residual_distribution_stats, dti_implausible_mask_stats, suffix=''):
    ''' 
    Write DTI QC measures in a XML file (assessor)

    :param dti_residual_distribution_stats: path to the scil dti residual image
                                            distribution stats
    :param dti_implausible_mask_stats: path to the csv summarizing % voxel occupation
                                       of the scil dti implausible voxels within the brainmask
    '''
    import lxml.etree as et
    import os.path as op
    import pandas as pd
    
    ## Load all the input files
    ############################################
    dti_res_dist_df = pd.read_csv(dti_residual_distribution_stats) if op.exists(dti_residual_distribution_stats) else None
    dti_perc_imp_df = pd.read_csv(dti_implausible_mask_stats) if op.exists(dti_implausible_mask_stats) else None
    
    ## Extract QC values
    ############################################
     
    # mean and max of dti residual image
    dti_res_mean = dti_res_dist_df['mean'].values[0] if dti_res_dist_df is not None else 'NA'
    dti_res_max = dti_res_dist_df['max'].values[0] if dti_res_dist_df is not None else 'NA'
    
    # % implausible voxels
    dti_perc_imp = dti_perc_imp_df['percent_voxels'].values[0] if dti_perc_imp_df is not None else 'NA'
    
       ## Package into an assessor
    ############################################
    string = '<rest:restAssessorData xmlns:rest="http://nrg.wustl.edu/rest" xmlns:prov="http://www.nbirn.net/prov" xmlns:xnat="http://nrg.wustl.edu/xnat" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"></rest:restAssessorData>'
        
    Assess = et.ElementTree(et.XML(string))
        
    namespaces={'dti':'http://nrg.wustl.edu/sti',
                'xnat':'http://nrg.wustl.edu/xnat',
                'xsi':'http://www.w3.org/2001/XMLSchema-instance'}
        
    avg_res = et.Element(et.QName(namespaces["dti"], "QC_mean_dti_residual"), nsmap=namespaces)
    avg_res.text = str(dti_res_mean)
    Assess.getroot().append(avg_res)
        
    max_res = et.Element(et.QName(namespaces["dti"], "QC_max_dti_residual"), nsmap=namespaces)
    max_res.text = str(dti_res_max)
    Assess.getroot().append(max_res)
        
    perc_imp = et.Element(et.QName(namespaces["dti"], "QC_percent_implausible_voxels"), nsmap=namespaces)
    perc_imp.text = str(dti_perc_imp)
    Assess.getroot().append(perc_imp)
    
    fname = op.abspath('dtiAssessor{}.xml'.format(suffix))
    with open(fname, 'wb') as f:
        Assess.write(f, pretty_print=True, encoding='utf-8', xml_declaration=True)
    
    return fname