import os
import logging

from astromatch import Catalogue
from astromatch import Match

from astropy.table import Table
from astropy import units as u
from astropy.table import hstack, join

from mocpy import MOC
import numpy as np

from config import *
import crossXmatch.ctp.fun_counterpart as fun_cntrprt

logger = logging.getLogger(__name__)

def get_xcat(
    xmocname, xcatname, acorr=False, mlcoords=False, telescope='XMM'
    ):
    
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""
        
    logger.info(f'Reading X-ray MOC from {xmocname} ...')
    xmoc = MOC.from_fits(xmocname)       
    logger.info(f'Reading X-ray sources from {xcatname} ...')
    xcat_tab = Table.read(xcatname)
    
    logger.info(f'There are {len(xcat_tab)} X-ray sources')
    if telescope == 'XMM':
        id_col = 'UXID'
        if acorr == False:
            xcat_tab['RAX']  = xcat_tab['RA'].data  * u.deg
            xcat_tab['DECX'] = xcat_tab['DEC'].data * u.deg
            coord_cols = ['RAX', 'DECX']
        else:
            xcat_tab[f'RA_CORR_{acorr}X']  = xcat_tab[f'RA{mltag}_CORR_{acorr}'].data  * u.deg
            xcat_tab[f'DEC_CORR_{acorr}X'] = xcat_tab[f'DEC{mltag}_CORR_{acorr}'].data * u.deg
            coord_cols = [f'RA_CORR_{acorr}X', f'DEC_CORR_{acorr}X']

        if mlcoords :
            if type(xcat_tab.mask)!=type(None):
                xcat_tab = xcat_tab[~xcat_tab.mask['RADEC_ERR']]
                m = xcat_tab.mask['RADECML_ERR']
                xcat_tab = xcat_tab[~m]
            #print('hiiiiiiii')
            #radec_err = xcat_tab['RADECML_ERR']
            #merr = xcat_tab.mask['RADECML_ERR'] #np.logical_or(np.isnan(radec_err), np.isinf(radec_err))
            #print(xcat_tab[merr]['RADECML_ERR'] )
            #print(xcat_tab[merr]['RADECML_ERR'] )
            #print(xcat_tab_kk['RADECML_ERR'] )
            #kk = np.full(len(xcat_tab['RADECML_ERR']), 1.1)
            #xcat_tab_kk = np.ma.compress_rowcols(xcat_tab_kk)
            #m = xcat_tab_kk.mask['RADECML_ERR']
            #print(xcat_tab_kk[~m]['RADECML_ERR'])
            #median_err = np.median(xcat_tab[~merr]['RADECML_ERR'])
            #kk[merr]  = median_err
            #kk[~merr] = xcat_tab[~merr]['RADECML_ERR']
            #xcat_tab['RADECML_ERR_kk'] = kk
            #mfinite = np.isfinite(xcat_tab['RADECML_ERR'])
            #print( all(np.isfinite(xcat_tab['RADECML_ERR'])) )
            poserr_cols = ['RADECML_ERR']
        else:
            poserr_cols = ['RADEC_ERR']
            
    elif telescope == 'chandra':
        id_col = 'name'
        xcat_tab['RAX']  = xcat_tab['RA'].data  * u.deg
        xcat_tab['DECX'] = xcat_tab['Dec'].data * u.deg
        coord_cols = ['RAX', 'DECX']
        poserr_cols = ['Pos_error']
              
    xcat = Catalogue(
                     xcat_tab,
                     name='C',
                     id_col=id_col,
                     coord_cols=coord_cols,
                     poserr_cols=poserr_cols,
                     poserr_type='circle',
                     area=xmoc,
    )
    return xcat

def get_Optcat(ocatname, omocname, optCat, pos_error_optCat, match_code='official'):
    """_summary_

    Args:
        ocatname (_type_): _description_
        omocname (_type_): _description_
        optCat (_type_): _description_
        pos_error_optCat (_type_): _description_

    Returns:
        _type_: _description_
    """
    logger.info(f'Reading optical MOC from {omocname} ...')
    omoc = MOC.from_fits(omocname)       
    
    logger.info(f'Reading optical sources from {ocatname} ...')    
    ocat_tab = Table.read(ocatname)
    
    ocat_tab['Pos_error'] = pos_error_optCat
        
    name, id_col, coord_cols, poserr_cols, poserr_type, mag_cols = get_paramsOptCat(optCat)
    
    if match_code == 'age':
        mag_cols = ['MAG_W2',  'MAG_W1W2', 'PNT', 'MAG_G', 'MAG_RW2', 'GAIA_G']
    
    #two lines below needed bc some sources have two rows w same id and similar data
    #need to take only one of them if not is problematic later
    un, m = np.unique(ar=ocat_tab[id_col].data, return_index=True)
    ocat_tab = ocat_tab[m] 
    #print(mag_cols)
    #g(name)
    ocat = Catalogue(
                     ocat_tab,
                     name=name,
                     id_col=id_col,
                     coord_cols=coord_cols,
                     poserr_cols=poserr_cols,
                     poserr_type=poserr_type,
                     area=omoc,
                     mag_cols=mag_cols
    ) 
    return ocat
    
    
def get_paramsOptCat(optCat):
    """_summary_

    Args:
        optCat (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    name = cat_params[f'{optCat}']['name']
    id_col = cat_params[f'{optCat}']['id_col']
    coord_cols = cat_params[f'{optCat}']['coord_cols']
    poserr_cols = cat_params[f'{optCat}']['poserr_cols']
    poserr_type = cat_params[f'{optCat}']['poserr_type']
    mag_cols = cat_params[f'{optCat}']['mag_cols']
    
    return name, id_col, coord_cols, poserr_cols, poserr_type, mag_cols

def get_all_matches(xm, name_xcat='C', name_ocat=None, match_code='official'):
    """_summary_

    Args:
        xm (_type_): _description_
        name_xcat (str, optional): _description_. Defaults to 'C'.
        name_ocat (str, optional): _description_. Defaults to 'DES'.

    Returns:
        _type_: _description_
    """
    all_matches = xm.get_matchs(match_type='all')
    #print(len(all_matches))
    if match_code == 'official':
        dra, ddec = xm.offset(f'{name_xcat}', f'{name_ocat}', match_type='all')
    else:
        dra, ddec = xm.offset(f'{name_xcat}', f'{name_ocat}') #, match_type='all')
    dr = np.sqrt(dra * dra + ddec * ddec)

    all_matches['CTP_DRA']  = dra
    all_matches['CTP_DDEC'] = ddec
    all_matches['CTP_DR']   = dr

    all_matches.rename_column('LR_BEST', 'CTP_LR') 
    all_matches.rename_column('REL_BEST', 'CTP_REL') # Relaibility ~= purity
    all_matches.rename_column('LR_BEST_MAG', 'CTP_PRIOR')
    all_matches.rename_column('match_flag', 'CTP_BEST')            
    all_matches.rename_column('prob_this_match', 'CTP_RELQUAL')
    return all_matches

def get_enrichedOptcat(optCat, ocatname, in_cat):
    """
    Function that returns a catalogue from xmatch with properties
    from the initial optical catalogue. Allowed cols are defined in a 
    separate script config.py.

    Args:
        optCat (str): _description_
        ocatname (_type_): _description_
        in_cat (Table): Input catalogue

    Returns:
        _type_: _description_
    """
    allowed_cols = allowed_cols_optCat[f'{optCat}']
    print(allowed_cols)
    #if optCat == 'LS':
    #    print(allowed_cols)
    
    FullOpt = Table.read(ocatname, format='fits') #Initial optical catalogue
    temp_opt_tab = Table()
    
    id_col = cat_params[f'{optCat}']['id_col']
    temp_opt_tab['UOBJID'] = FullOpt[f'{id_col}'].astype(str)

    for col in FullOpt.columns:
        if np.logical_and( col in allowed_cols,
                           not(col in in_cat.columns) ):   
            temp_opt_tab[col] = FullOpt[col]            

    # this for all possible counterparts
    enocat_all = join(in_cat,
                      temp_opt_tab,
                      keys='UOBJID', 
                      join_type='left', 
                      table_names=["1", "2"])
    enocat_all.remove_columns('UOBJID')
    return enocat_all
    
def get_enrichedXcat(xcatname, xcat, telescope): 
    FullX = Table.read(xcatname, format='fits')
    tab_x = xcat.save() #Convert catalogue object to table
    tab_x['UOBJID'] = tab_x['SRCID_C']
    temp_x_tab = Table()
    if telescope =='XMM':
        temp_x_tab['UOBJID'] = FullX['UXID'].astype(str)
        allowed_cols = allowed_cols_XMM
    else:
        temp_x_tab['UOBJID'] = FullX['name'].astype(str)
        allowed_cols = allowed_cols_chandra

    for col in FullX.columns:
        if(col in allowed_cols):
            temp_x_tab[col] = FullX[col]

    enxcat = join(tab_x, temp_x_tab, keys='UOBJID', join_type='left')
    enxcat.remove_columns('UOBJID') 
    return enxcat

def get_join_enrich(enxcat, enocat_all, optCat):
    enxcat['UX'] = enxcat['SRCID_C']        
    enocat_all['UX'] = enocat_all['SRCID_C_C']

    joined_x_opt_cat = join(enxcat, 
                            enocat_all,
                            keys='UX',
                            join_type='left', 
                            table_names=["1", "2"])
    joined_x_opt_cat.remove_columns('UX')
    joined_x_opt_cat.remove_columns(['SRCID_C', 'SRCID_C_C', f'SRCID_{optCat}_{optCat}', 
                                    f'SRCID_{optCat}_match', 'SRCID_C_match'])
    return joined_x_opt_cat

def get_tab_orderRel(tab, telescope):
    tab['CTP_REL1'] = - tab['CTP_REL']
    tab['CTP_REL1'][tab['CTP_REL'].mask] = 0
    
    if telescope == 'XMM':
        id_xray = 'UXID'
    else:
        id_xray = 'name'

    g = tab.group_by([f'{id_xray}', 'CTP_REL1'])
    g['NUM_CTP'] = np.zeros(len(g), dtype='int')
    g['CTP_RANK'] = np.zeros(len(g), dtype='int')

    g1 = tab.group_by([f'{id_xray}'])
    a = g1.groups.indices
    n = a[1:] - a[0:-1]

    for i in range(len(a) - 1):
        g['NUM_CTP'][  a[i]:a[i+1] ] = n[i]
        g['CTP_RANK'][ a[i]:a[i+1] ] = range(1, n[i]+1)
        
    g.remove_columns('CTP_REL1')
    return g

def get_enriched_cat(clst, xcatname, ocatname, xcat, ocat, all_matches, optCat, out_pth, telescope):
    
    xcat_all = xcat.select_by_id( all_matches['SRCID_C'] ) 
    ocat_all = ocat.select_by_id( all_matches[f'SRCID_{optCat}'] )
    
    all_x_opt_ctlg = hstack( # Catalogues converted into Tables using the `save` method
                            [ xcat_all.save(), ocat_all.save(), all_matches ], 
                            table_names=[xcat_all.name, ocat_all.name, "match"],
                            join_type='exact'
    )

    all_x_opt_ctlg['UOBJID'] = all_x_opt_ctlg[f'SRCID_{optCat}_{optCat}']
    
    enocat_all = get_enrichedOptcat(optCat, ocatname, all_x_opt_ctlg)
    enxcat = get_enrichedXcat(xcatname, xcat, telescope)
    en_x_opt_ctlg = get_join_enrich(enxcat, enocat_all, optCat)
    
    g = get_tab_orderRel(en_x_opt_ctlg, telescope)
    
    return g
    
    
def do_match(#clst, 
             xcat,  
             ocat, 
             method='lr', 
             r_srch=3.5*u.arcsec,
             pos_error_optCat=0.1, 
             #optCat=None,
             #out_pth=None
             ): 
    """_summary_

    Args:
        clst (_type_): _description_
        xmocname (_type_): _description_
        xcatname (_type_): _description_
        omocname (_type_): _description_
        ocatname (_type_): _description_
        method (str, optional): _description_. Defaults to 'lr'.
        r_srch (_type_, optional): _description_. Defaults to 3.5*u.arcsec.
        pos_error_optCat (_type_, optional): _description_. Defaults to 0.1#arcsecDES=False.
        LS (bool, optional): _description_. Defaults to False.
        VHS (bool, optional): _description_. Defaults to False.
    """
    xm = Match(xcat, ocat)
    
    match_results_lr = xm.run(method=method, 
                              radius=r_srch, 
                              poserr_dist="normal")
    xm.set_best_matchs(cutoff=-9) #Gets all the objects because LR > 0

    #all_matches = get_all_matches(xm, name_xcat='C', name_ocat=f'{optCat}')
    #save_all_matches(out_pth, clst, optCat, all_matches)
        
    return xm

def do_match_age(xcat,
                 ocat,
                 radius=7*u.arcsec, 
                 priorfile="priors3D_LS10_4xmm_MW.fits"):

    from astromatch.priorsND import PriorND

    #mags = [['mag_r'], ['mag_z'], ['mag_w1']]
    #magmin = [ [10] ,[10], [5] ]
    #magmax = [ [30], [30], [25] ]
    #magbinsize = [ [0.2], [0.2], [0.2] ]

    priors = PriorND.from_table(priorfile)
    mags = priors.mags
        
    xm = Match( xcat, ocat )
    
    match_results = xm.run(method='lr',
                           radius=7*u.arcsec,
                           priors=priors,
                           mags=priors.mags,
                           magmin=priors.magmin, 
                           magmax=priors.magmax,
                           magbinsize=priors.magbinsize, 
                           poserr_dist="normal"
    )
    xm.set_best_matchs(cutoff=-9.)
    
    return xm 
    #lr_stats = xm.stats(mincutoff=-4, maxcutoff=1.5, cutoffstep=0.01)
    #lr_stats.write("stats.fits", format='fits', overwrite=True)
    #plot_stats(lr_stats, "legacy")
    
def save_all_matches(out_pth, clst, optCat, all_matches, match_code='oficial'):
    path_all_matches = out_pth + f'all_matches/'
    
    if not os.path.exists(path_all_matches):
        logger.info(f'Creating... {path_all_matches}...')
        os.makedirs(path_all_matches)
    logger.info('Here!')
    logger.info(f'Writing ALL matches in {path_all_matches}...')
    
    if match_code == 'official':
        all_matches.write(
            path_all_matches + f'{clst}_allMatch_{optCat}.fits', format='fits', overwrite=True
        )
    else:
        all_matches.write(
            path_all_matches + f'{clst}_allMatch_{optCat}_age.fits', format='fits', overwrite=True
        )
    
def save_enrich(out_pth, clst, optCat, ctlg, match_code='oficial'):
    path_out = out_pth + f'enrich/'
    if not os.path.exists(path_out):
        logger.info(f'Creating... {path_out}...')
        os.makedirs(path_out)
    logger.info(f'Writing enriched in {path_out}...')
    
    if match_code == 'official':
        ctlg.write(path_out+f'{clst}_enrich_{optCat}.fits', format='fits', overwrite=True)
    else:
        print(path_out+f'{clst}_enrich_{optCat}_age.fits')
        ctlg.write(path_out+f'{clst}_enrich_{optCat}_age.fits', format='fits', overwrite=True)
    
def write_all_matches(clst, 
                      xmocname, 
                      xcatname, 
                      path_optCat, 
                      method='lr', 
                      r_srch=3.5*u.arcsec,
                      pos_error_optCat=0.1, #arcsec
                      cat_lst=None,
                      pathOptCat=None,
                      out_pth=None,
                      save_matches=False,
                      xacorr=False,
                      mlcoords=False,
                      telescope='XMM',
                      match_code='official'
):

    logger.info(f'Looking for matches of {clst}...')
    
    xcat = get_xcat(xmocname, xcatname, acorr=xacorr, mlcoords=mlcoords, telescope=telescope)

    xm_dct = {}
    for optCat in cat_lst:
        logger.info(f'\tSearching in {optCat}...')

        #if telescope=='XMM'
        ocatname = pathOptCat + f'{clst}_{optCat}_rband_filter.fits'   
        omocname = pathOptCat + f'{clst}_{optCat}.moc'         
        #else:
        ocat = get_Optcat(ocatname, omocname, optCat, pos_error_optCat, match_code=match_code)  
        
        if match_code == 'official':
            xm = do_match(#clst, 
                           xcat,  
                           ocat, 
                           method='lr', 
                           r_srch=r_srch,
                           pos_error_optCat=0.1,
                           #optCat=optCat, 
                           #out_pth=out_pth,
            )
        else:
            xm = do_match_age(xcat,
                              ocat,
                              radius=r_srch, 
                              priorfile="./eso/ctp/priors3D_LS10_4xmm_MW.fits"
                 )
        
        xm_dct.update( {f'{optCat}':xm} )
        
        if save_matches == True:
            
            all_matches = get_all_matches(
                xm_dct[f'{optCat}'], name_xcat='C', name_ocat=f'{optCat}', match_code=match_code
            )
            save_all_matches(out_pth, clst, optCat, all_matches, match_code=match_code)

            enriched_cat = get_enriched_cat(
                clst, xcatname, ocatname, xcat, ocat, all_matches, optCat, out_pth, telescope
            )
            #print('kk', enriched_cat.column )
            save_enrich(out_pth, clst, optCat, enriched_cat, match_code=match_code)
    
    return xm_dct

        
        