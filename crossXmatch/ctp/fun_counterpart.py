"""

"""

import os

from astropy import units as u
from astropy.coordinates import SkyCoord, Angle, Longitude, Latitude
from astropy_healpix import HEALPix
from astropy.io import fits
from astropy.table import Table
from astropy.table import hstack, vstack, join

from astromatch import Catalogue
from astromatch import Match

#from dl import authClient as ac, queryClient as qc
#from dl.helpers.utils import convert
#from getpass import getpass

#import matplotlib.pyplot as plt

from mocpy import MOC

import numpy as np

import pyvo as vo


def xraymoc(xmap="DATA/pl26_soft_cell70.exp"):
    #from mocpy import MOC
    #from astropy.io import fits
    hdu = fits.open(xmap)
    
    data = hdu[0].data
    mocx = MOC.from_fits_image(hdu[0], max_norder=14, mask=data>0)
    hdu.close()
    return mocx


def runxmatch_des(
    xcatname="data/pl26_srclist_v1.fits", 
    xmocname="data/pl26_soft_cell70.moc", 
    ocatname="data/pl26_des1y3.fits", 
    omocname="data/pl26_des1y3.moc",
    exp_file="/storage2/ivan/chandra/pl26/pl26_soft_cell70.exp",
    corename='pl26',
    out_pth='/storage/ivan/output_matchDES/'
):
    """
    
    
    Input:
    ------
    xcatname (string):
    xmocname (string):
    ocatname (string):
    omocname (string):
    exp_file (string):
    corename (string):
    out_pth  (string):
    
    Output:
    -------
    
    """

    if(not os.path.isfile(xmocname)):
       xmoc = xraymoc(xmap=exp_file)
       xmoc.write(out_pth + f"data/{corename}_soft_cell70.moc", format='fits', overwrite=True)
       
    xmoc = MOC.from_fits(xmocname)       
    xcat = Table.read(xcatname)
    xcat['RAX'] = xcat['RA'].data*u.deg
    xcat['DECX'] = xcat['Dec'].data*u.deg
  
    omoc = MOC.from_fits(omocname)           
    ocat = Table.read(ocatname)
    ocat['Pos_error'] = 0.1
    un, m = np.unique(ar=ocat['coadd_object_id'].data, return_index=True)
    ocat = ocat[m]  

    xcat = Catalogue(
        xcat,
        name='C',
        id_col='name',
        coord_cols=['RAX', 'DECX'],
        poserr_cols=['Pos_error'],
        poserr_type='circle',
        area=omoc,
    )
        
    ocat = Catalogue(
        ocat,
        name='DES',
        id_col='coadd_object_id',
        coord_cols=['alphawin_j2000', 'deltawin_j2000'],
        poserr_cols=['Pos_error'],
        poserr_type='circle',
        area=omoc,
        mag_cols=['mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z', 'mag_auto_y']
    )

    xm = Match(xcat, ocat)
    match_results_lr = xm.run(method='lr', radius=3.5*u.arcsec, poserr_dist="normal")

    xm.set_best_matchs(cutoff=-9)
    all_matches = xm.get_matchs(match_type='all')
    dra, ddec = xm.offset('C', 'DES', match_type='all')
    dr = np.sqrt(dra * dra + ddec * ddec)
    all_matches = xm.get_matchs(match_type='all') #astropy table
    dra, ddec = xm.offset('C', 'DES', match_type='all')
    dr = np.sqrt(dra * dra + ddec * ddec)
    all_matches['CTP_DRA'] = dra
    all_matches['CTP_DDEC'] = ddec
    all_matches['CTP_DR'] = dr
    all_matches.rename_column('LR_BEST', 'CTP_LR') #CTP==counterpart
    all_matches.rename_column('REL_BEST', 'CTP_REL')#relaibility ==purity
    all_matches.rename_column('LR_BEST_MAG', 'CTP_PRIOR')
    all_matches.rename_column('match_flag', 'CTP_BEST')            
    all_matches.rename_column('prob_this_match', 'CTP_RELQUAL')

    all_matches.write(out_pth + 'tmp1.fits', format='fits', overwrite=True)
     
    #Up to here only information of the matches, then this table is joined with informaiton of the initial 
    #X-ray catalogue and optical.
    #For that first create something tha tonly has the ids to match
    xcat_all = xcat.select_by_id(all_matches['SRCID_C']) #xcat and ocat have the minimal information of the catalogue 
    ocat_all = ocat.select_by_id(all_matches['SRCID_DES'])
    
    
    from astropy.table import hstack
    all = hstack(
        [xcat_all.save(), ocat_all.save(),all_matches],  # Catalogues can be converted into Tables using the `save`  method
        table_names=[xcat_all.name, ocat_all.name,"match"],
        join_type='exact')
    all['UOBJID']=all['SRCID_DES_DES']
    #print(all)

    # enrich optical catalogues - add the additional  info into the best and all catalogues
    allowed_cols=['coadd_object_id', 'ALPHAWIN_J2000' ,'DELTAWIN_J2000', 'ebv_sfd98', 'mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z', 'm.mag_auto_y', 'm.spread_model_i', 'spreaderr_model_i', 'imaflags_iso_g',  'imaflags_iso_r', 'imaflags_iso_i', 'imaflags_iso_z', 'imaflags_iso_y',  'bpz_zmode_mof', 'bpz_zmode_sof', 'bpz_zsigma68_mof', 'bpz_zsigma68_sof', 'dnf_zmc_mof', 'dnf_zmc_sof']
                  
    FullOpt = Table.read(ocatname, format='fits')
    t1 = Table()
    t1['UOBJID'] = FullOpt['coadd_object_id'].astype(str)
    for col in FullOpt.columns:
        if((col in allowed_cols) and not(col in all.columns)):                
            t1[col]= FullOpt[col]            

    # this for all possible counterparts
    enocat_all = join(all, t1, keys='UOBJID', join_type='left', table_names=["1","2"])
    enocat_all.remove_columns('UOBJID')


    # enrich X-ray catalogue with additional columns. This is the same for both
    # best and all-possible counterparts
    FullX = Table.read(xcatname, format='fits')
    t = xcat.save()
    t['UOBJID']=t['SRCID_C']
    t1 = Table()
    t1['UOBJID'] = FullX['name'].astype(str)
    for col in FullX.columns:
        if(col  in ['name', 'Full_prob', 'Full_flux', 'Full_flux_Bay', 'Soft_prob', 'Soft_flux', 'Soft_flux_Bay', 'Hard_prob', 'Hard_flux', 'Hard_flux_Bay']):
            t1[col]= FullX[col]

    enxcat = join(t, t1, keys='UOBJID', join_type='left')
    enxcat.remove_columns('UOBJID')
   
    # match enxcat with enocat 
    enxcat['UX'] = enxcat['SRCID_C']        
    enocat_all['UX'] = enocat_all['SRCID_C_C']
    new = join(enxcat, enocat_all, keys='UX', join_type='left', table_names=["1","2"])
    new.remove_columns('UX')
    new.remove_columns(['SRCID_C', 'SRCID_C_C', 'SRCID_DES_DES', 'SRCID_DES_match', 'SRCID_C_match'])

    #Order sources by relaibility
    new['CTP_REL1'] = -new['CTP_REL']
    new['CTP_REL1'][new['CTP_REL'].mask] = 0
    g = new.group_by(['name', 'CTP_REL1'])
    g['NUM_CTP'] = np.zeros(len(g), dtype='int')
    g['CTP_RANK'] = np.zeros(len(g), dtype='int')
    g1 = new.group_by(['name'])
    a = g1.groups.indices
    n = a[1:]-a[0:-1]
    for i in range(len(a)-1):
        g['NUM_CTP'][a[i]:a[i+1]] = n[i]
        g['CTP_RANK'][a[i]:a[i+1]] = range(1,n[i]+1)
        #print(new['NUM_CTP', 'CTP_RANK'][a[i]:a[i+1]])
    g.remove_columns('CTP_REL1')
            
    g.write(out_pth + f'{corename}_des1y3_match.fits', format='fits', overwrite=True)
    
    return xm 