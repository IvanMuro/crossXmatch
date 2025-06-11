from config import *
import numpy as np 

def get_query_des1_main(ra, dec, radius_deg):
    """
    This function returns the query to get darta from DES DR1 main. See data description: 
    https://datalab.noirlab.edu/query.php?name=des_dr1.main

    Input:
    ------
    - ra (float): decimal degrees
    - dec (float): decimal degrees
    - radius_deg (float): decimal degrees

    Output:
    -------  
    query
    """
    LS_COLS = """
            m.coadd_object_id, m.ALPHAWIN_J2000 , m.DELTAWIN_J2000, m.ebv_sfd98, m.mag_auto_g, m.mag_auto_r, m.mag_auto_i, m.mag_auto_z, m.mag_auto_y, m.spread_model_i, m.spreaderr_model_i, m.imaflags_iso_g,  m.imaflags_iso_r, m.imaflags_iso_i, m.imaflags_iso_z, m.imaflags_iso_y
            """
    query = f"""
            SELECT {LS_COLS}
            FROM des_dr1.main as m
            WHERE  Q3C_RADIAL_QUERY(ra, dec, {ra},{dec},{radius_deg})
            """
    return query

def get_query_des_gold(ra, dec, radius_deg):
    """
    This function returns the query to get darta from DES Y1GOLD DR1. See data description: 
    https://datalab.noirlab.edu/query.php?name=des_dr1.y3_gold
    """
    LS_COLS = """
            p.coadd_object_id, p.ALPHAWIN_J2000 , p.DELTAWIN_J2000, p.bpz_zmode_mof, p.bpz_zmode_sof, p.bpz_zsigma68_mof, p.bpz_zsigma68_sof, p.dnf_zmc_mof, p.dnf_zmc_sof
            """
            
    query = f"""
            SELECT {LS_COLS}
            FROM des_dr1.y3_gold as p
            WHERE  Q3C_RADIAL_QUERY(p.ra, p.dec, {ra},{dec},{radius_deg}) 
            """ 
    return query


def get_query_des_gold_main(ra, dec, radius_deg):
    """
    This function returns the query to get data from DES Y1GOLD DR1. See data description: 
    https://datalab.noirlab.edu/query.php?name=des_dr1.y3_gold
    """
    LS_COLS = """
            p.coadd_object_id, p.ALPHAWIN_J2000 , p.DELTAWIN_J2000, p.ebv_sfd98, p.mag_auto_g,
            p.mag_auto_r, p.mag_auto_i, p.mag_auto_z, p.mag_auto_y, p.spread_model_i,
            p.spreaderr_model_i, p.imaflags_iso_g,  p.imaflags_iso_r, p.imaflags_iso_i,
            p.imaflags_iso_z, p.imaflags_iso_y, p.coadd_object_id, p.ALPHAWIN_J2000,
            p.DELTAWIN_J2000, p.bpz_zmode_mof, p.bpz_zmode_sof, p.bpz_zsigma68_mof,
            p.bpz_zsigma68_sof, p.dnf_zmc_mof, p.dnf_zmc_sof
            """
            
    query = f"""
            SELECT {LS_COLS}
            FROM des_dr1.y3_gold as p
            WHERE  Q3C_RADIAL_QUERY(p.ra, p.dec, {ra},{dec},{radius_deg}) 
            """ 
    return query

def get_query_ls10(ra, dec, radius_deg):
    """
    This function returns the query to get darta from ls_dr9.tractor. See data description: 
    https://datalab.noirlab.edu/query.php?name=ls_dr9.tractor
    """
    LS_COLS = """
                ls.release, ls.ls_id, ls.brickid, ref_cat, ls.objid, ls.type,
                ls.ra, ls.dec, ls.ra_ivar, ls.dec_ivar, ls.glat, ls.glon,
                ls.mag_g, ls.mag_i, ls.mag_r, ls.mag_z, ls.mag_w1, ls.mag_w2, ls.mag_w3, ls.mag_w4,
                ls.dered_flux_g, ls.dered_flux_i, ls.dered_flux_r, ls.dered_flux_z,
                ls.dered_flux_w1, ls.dered_flux_w2, ls.dered_flux_w3, ls.dered_flux_w4,
                ls.snr_g, ls.snr_i, ls.snr_r, ls.snr_z, ls.snr_w1, ls.snr_w2, ls.snr_w3, ls.snr_w4,
                ls.gaia_phot_g_mean_mag, ls.gaia_phot_bp_mean_mag, ls.gaia_phot_rp_mean_mag,
                ls.gaia_phot_g_mean_flux_over_error, ls.gaia_phot_bp_mean_flux_over_error,
                ls.gaia_phot_rp_mean_flux_over_error,
                ls.pmra, pmra_ivar, ls.pmdec, ls.pmdec_ivar, ls.parallax,
                ls.parallax_ivar, ls.ref_epoch, ls.brick_primary,
                allmask_g, allmask_i, allmask_r, allmask_z,
                anymask_g, anymask_i, anymask_r, anymask_z,
                wisemask_w1, wisemask_w2, maskbits
            """
            
    query = f"""
            SELECT {LS_COLS}
            FROM ls_dr10.tractor as ls
            WHERE Q3C_RADIAL_QUERY(ls.ra, ls.dec, {ra},{dec},{radius_deg}) 
            """ 
    return query

def get_query_ls10Mara(ra, dec, radius_deg):
    """
    This function returns the query to get darta from ls_dr9.tractor. See data description: 
    https://datalab.noirlab.edu/query.php?name=ls_dr9.tractor
    """
    LS_COLS = """
                ls.release, ls.brickid, ls.brickname, ls.objid, ls.brick_primary, ls.maskbits, ls.fitbits, 
                ls.ra, ls.dec, ls.type, ls.ra_ivar, ls.dec_ivar, ls.bx, ls.by, ls.dchisq_1, ls.ebv, ls.mjd_min, ls.mjd_max, 
                ls.ref_cat, ls.ref_id, ls.pmra, ls.pmdec, ls.parallax, ls.pmra_ivar, ls.pmdec_ivar, 
                ls.parallax_ivar, ls.ref_epoch, ls.gaia_phot_g_mean_mag, ls.gaia_phot_g_mean_flux_over_error,
                ls.gaia_phot_g_n_obs, ls.gaia_phot_bp_mean_mag, ls.gaia_phot_bp_mean_flux_over_error, 
                ls.gaia_phot_bp_n_obs,gaia_phot_rp_mean_mag, ls.gaia_phot_rp_mean_flux_over_error, 
                ls.gaia_phot_rp_n_obs, ls.gaia_phot_variable_flag, ls.gaia_astrometric_excess_noise, 
                ls.gaia_astrometric_excess_noise_sig, ls.gaia_astrometric_n_obs_al, 
                ls.gaia_astrometric_n_good_obs_al, ls.gaia_astrometric_weight_al, ls.gaia_duplicated_source,
                ls.gaia_a_g_val, ls.gaia_e_bp_min_rp_val, ls.gaia_phot_bp_rp_excess_factor, 
                ls.gaia_astrometric_sigma5d_max, ls.gaia_astrometric_params_solved, ls.flux_g, ls.flux_r, 
                ls.flux_i, ls.flux_z, ls.flux_w1, ls.flux_w2, ls.flux_w3, ls.flux_w4, ls.flux_ivar_g,
                ls.flux_ivar_r, ls.flux_ivar_i, ls.flux_ivar_z, ls.flux_ivar_w1, ls.flux_ivar_w2,
                ls.flux_ivar_w3, ls.flux_ivar_w4, ls.fiberflux_g, ls.fiberflux_r, ls.fiberflux_i,
                ls.fiberflux_z, ls.fibertotflux_g, ls.fibertotflux_r, ls.fibertotflux_i, ls.fibertotflux_z, 
                ls.mw_transmission_g, ls.mw_transmission_r, ls.mw_transmission_i, ls.mw_transmission_z,
                ls.mw_transmission_w1, ls.mw_transmission_w2, ls.mw_transmission_w3, ls.mw_transmission_w4,
                ls.nobs_g, ls.nobs_r, ls.nobs_i, ls.nobs_z, ls.nobs_w1, ls.nobs_w2, ls.nobs_w3, ls.nobs_w4, 
                ls.rchisq_g, ls.rchisq_r, ls.rchisq_i, ls.rchisq_z, ls.rchisq_w1, ls.rchisq_w2, ls.rchisq_w3, 
                ls.rchisq_w4, ls.fracflux_g, ls.fracflux_r, ls.fracflux_i, ls.fracflux_z, ls.fracflux_w1, 
                ls.fracflux_w2, ls.fracflux_w3, ls.fracflux_w4, ls.fracmasked_g, ls.fracmasked_r, 
                ls.fracmasked_i, ls.fracmasked_z, ls.fracin_g, ls.fracin_r, ls.fracin_i, ls.fracin_z, 
                ls.ngood_g, ls.ngood_r, ls.ngood_i, ls.ngood_z, ls.anymask_g, ls.anymask_r, ls.anymask_i,
                ls.anymask_z, ls.allmask_g, ls.allmask_r, ls.allmask_i, ls.allmask_z, ls.wisemask_w1,
                ls.wisemask_w2, ls.psfsize_g, ls.psfsize_r, ls.psfsize_i, ls.psfsize_z, ls.psfdepth_g,
                ls.psfdepth_r, ls.psfdepth_i, ls.psfdepth_z, ls.galdepth_g, ls.galdepth_r, ls.galdepth_i,
                ls.galdepth_z, ls.nea_g, ls.nea_r, ls.nea_i, ls.nea_z, ls.blob_nea_g, ls.blob_nea_r,
                ls.blob_nea_i, ls.blob_nea_z, ls.psfdepth_w1, ls.psfdepth_w2, ls.psfdepth_w3, ls.psfdepth_w4,
                ls.wise_coadd_id, ls.wise_x, ls.wise_y, ls.sersic, ls.sersic_ivar, ls.shape_r, ls.shape_r_ivar, 
                ls.shape_e1, ls.shape_e1_ivar, ls.shape_e2, ls.shape_e2_ivar
            """ 
            
                #ls.apflux_g, ls.apflux_r, ls.apflux_i, ls.apflux_z, ls.apflux_resid_g, ls.apflux_resid_r, 
                #ls.apflux_resid_i, ls.apflux_resid_z, ls.apflux_blobresid_g, ls.apflux_blobresid_r, 
                #ls.apflux_blobresid_i, ls.apflux_blobresid_z, ls.apflux_ivar_g, ls.apflux_ivar_r, 
                #ls.apflux_ivar_i, ls.apflux_ivar_z, ls.apflux_masked_g, ls.apflux_masked_r, 
                #ls.apflux_masked_i, ls.apflux_masked_z, ls.apflux_w1, ls.apflux_w2, ls.apflux_w3, 
                #ls.apflux_w4, ls.apflux_resid_w1, ls.apflux_resid_w2, ls.apflux_resid_w3, ls.apflux_resid_w4,
                #ls.apflux_ivar_w1, ls.apflux_ivar_w2, ls.apflux_ivar_w3, ls.apflux_ivar_w4, 
                
                
                #ra_2,dec_2,Separation,
                #ls.lc_flux_w1, ls.lc_flux_w2, ls.lc_flux_ivar_w1, ls.lc_flux_ivar_w2, ls.lc_nobs_w1, 
                #ls.lc_nobs_w2, ls.lc_fracflux_w1, ls.lc_fracflux_w2, ls.lc_rchisq_w1, ls.lc_rchisq_w2, 
                #ls.lc_mjd_w1, ls.lc_mjd_w2, ls.lc_epoch_index_w1, ls.lc_epoch_index_w2,
            
    query = f"""
            SELECT {LS_COLS}
            FROM ls_dr10.tractor as ls
            WHERE Q3C_RADIAL_QUERY(ls.ra, ls.dec, {ra},{dec},{radius_deg}) 
            """ 
    return query

def get_query_ls10MaraAp(ra, dec, radius_deg):
    """
    This function returns the query to get darta from ls_dr9.tractor. See data description: 
    https://datalab.noirlab.edu/query.php?name=ls_dr9.tractor
    """
    LS_COLS = """
                ls.apflux_g_1, ls.apflux_r_1, ls.apflux_i_1, ls.apflux_z_1, 
                ls.apflux_resid_g_1, ls.apflux_resid_r_1, ls.apflux_resid_i_1, ls.apflux_resid_z_1,
                ls.apflux_blobresid_g_1, ls.apflux_blobresid_r_1, 
                ls.apflux_blobresid_i_1, ls.apflux_blobresid_z_1, ls.apflux_ivar_g_1, ls.apflux_ivar_r_1, 
                ls.apflux_ivar_i_1, ls.apflux_ivar_z_1, ls.apflux_masked_g_1, ls.apflux_masked_r_1, 
                ls.apflux_masked_i_1, ls.apflux_masked_z_1, ls.apflux_w1_1, ls.apflux_w2_1, ls.apflux_w3_1, 
                ls.apflux_w4_1, ls.apflux_resid_w1_1, ls.apflux_resid_w2_1, ls.apflux_resid_w3_1, 
                ls.apflux_resid_w4_1, ls.apflux_ivar_w1_1, ls.apflux_ivar_w2_1, ls.apflux_ivar_w3_1, ls.apflux_ivar_w4_1
            """    
                #ra_2,dec_2,Separation,
                #ls.lc_flux_w1, ls.lc_flux_w2, ls.lc_flux_ivar_w1, ls.lc_flux_ivar_w2, ls.lc_nobs_w1, 
                #ls.lc_nobs_w2, ls.lc_fracflux_w1, ls.lc_fracflux_w2, ls.lc_rchisq_w1, ls.lc_rchisq_w2, 
                #ls.lc_mjd_w1, ls.lc_mjd_w2, ls.lc_epoch_index_w1, ls.lc_epoch_index_w2,
            
    query = f"""
            SELECT {LS_COLS}
            FROM ls_dr10.apflux as ls
            WHERE Q3C_RADIAL_QUERY(ls.ra, ls.dec, {ra},{dec},{radius_deg}) 
            """ 
    return query

def get_query_vhs(ra, dec, radius_deg):
    """
    This function returns the query to get darta from DES Y1GOLD DR1. See data description: 
    https://datalab.noirlab.edu/query.php?name=des_dr1.y3_gold
    """  
    LS_COLS = """vhs.sourceid, vhs.framesetid, vhs.sourcename, vhs.ra2000, vhs.dec2000, vhs.l, vhs.b,
    vhs.mergedclass, vhs.psaturated, vhs.pstar, vhs.pgalaxy, vhs.pnoise, vhs.priorsec, vhs.primary_source,
    vhs.hppErrBits, vhs.jppErrBits, vhs.ksppErrBits, vhs.jppErrBits,
    vhs.hapermag4,  vhs.hapermag4err,  vhs.hapermagnoapercorr4,  vhs.hclass,  vhs.hpetromag,  vhs.hpetromagerr,
    vhs.japermag4,  vhs.japermag4err,  vhs.japermagnoapercorr4,  vhs.jclass,  vhs.jpetromag,  vhs.jpetromagerr,
    vhs.ksapermag4, vhs.ksapermag4err, vhs.ksapermagnoapercorr4, vhs.ksclass, vhs.kspetromag, vhs.kspetromagerr,
    vhs.yapermag4,  vhs.yapermag4err,  vhs.yapermagnoapercorr4,  vhs.yclass,  vhs.ypetromag,  vhs.ypetromagerr
    """
    query = f"""
            SELECT {LS_COLS}
            FROM vhs_dr5.vhs_cat_v3 as vhs
            WHERE  Q3C_RADIAL_QUERY(vhs.ra2000, vhs.dec2000, {ra},{dec},{radius_deg}) 
            """ 
    return query

#mergedclass Class flag, 1|0|-1|-2|-3|-9 = gal|noise|star|probStar|probGal|saturated

def _constrain_LS10(data):
    m = data['type'] == 'DUP'
    data = data[np.logical_not(m)]
        
    radec_err = np.sqrt(1.0/(data['ra_ivar']) + 1.0/(data['dec_ivar'])) * 3600.0
    merr = np.logical_or(np.isnan(radec_err), np.isinf(radec_err))
    radec_err[merr] = astrom_error['LS'] #associate error for sources w/o error in the catalogue
    data['e_radec'] = radec_err
    #c=[]
    mask   = radec_err > 0 
    unique = np.asarray(range(radec_err.size))
    lguid  =  (data['objid']).astype(np.int_) + \
                 ((data['brickid'].astype(np.int_))<<16) + \
                 ((data['release'].astype(np.int_))<<40)
    data['UOBJID_LS'] = unique
    data['LGUID']  = lguid
    return data
    #c.append(fits.Column(name='UOBJID', format='K', array=unique[mask]))
    #c.append(fits.Column(name='LGUID', format='K', array=lguid[mask]))


