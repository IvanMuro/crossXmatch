import os
import logging

from astropy.coordinates import SkyCoord, Angle, Longitude, Latitude
from astropy import units as u
from astropy_healpix import HEALPix
from astropy.table import Table
from astropy.table import hstack, vstack, join

from dl import authClient as ac, queryClient as qc
from dl.helpers.utils import convert
from getpass import getpass

from mocpy import MOC
import numpy as np
import pyvo as vo

import eso.ctp.get_queries as get_queries
from eso.config import *

logger = logging.getLogger(__name__)

def getCat(r, d, optCat, radius_deg=1): #DES=False, LS=False, VHS=False, extrcols,
    """
    Function that queries noirlab to get DES data (https://datalab.noirlab.edu/des/index.php). 
    This retrieves the sources around input ra0 and dec0. 
       
    Input:
    ------
    ra0 (float): RA coordinate
    dec0 (float): Dec coordinate
    extrcols: ???? 
    des1_main (bool): reads main DR1 catalogue of DES
    des3_photo (bool): reads gold Y3 DR1 catalogue of DES
    radius_deg (float): Radius of search around ra0, dec0
    
    Output:
    -------
    allobj (Table):
    """
    ac.login('ivan_rodriguez', "6zEHQEu4YAzgWA")
    #allobj = Table()
        
    if optCat == 'DES':
        query = get_queries.get_query_des_gold_main(r, d, radius_deg)
    if optCat == 'LS':
        query = get_queries.get_query_ls10(r, d, radius_deg)
    if optCat == 'LS_mara':
        query = get_queries.get_query_ls10Mara(r, d, radius_deg)
    if optCat == 'LS_maraAp':
        query = get_queries.get_query_ls10MaraAp(r, d, radius_deg)
    if optCat== 'VHS':
        query = get_queries.get_query_vhs(r, d, radius_deg)
    #print(query)
    
    #print (qc.get_timeout_request(), qc.isAlive())
    if(not qc.isAlive()):
        ac.login('ivan_rodriguez',"6zEHQEu4YAzgWA")
    
    #result = qc.query(sql=query, fmt='table', async=True, wait=True)
    #print(qc.isAlive())
    result = qc.query(sql=query, fmt='table')
    
    if optCat =='LS':
        result = get_queries._constrain_LS10(result)
        
    #print(result)
    #if(len(result)>0):
        #for k in extrcols.keys():
        #    result[k] = [extrcols[k][i]]*len(result)
        #if(i==0):
    #    allobj = result           
        #else:
        #    allobj = vstack([allobj, result])
    #print(result)
    ac.logout()
    
    return result #allobj

def makeMOC(ra, dec, radius_deg):
    moc_optical = MOC.from_cone(
        lon=ra*u.deg,
        lat=dec*u.deg,
        radius=radius_deg*u.deg,
        max_depth=14
    )
    return moc_optical 

def makeMOC_opt(table, optCat):
    """
    
    ** FUTURE UPDATE**:
    CHECK WHICH IS THE RESOLUTION OF max_norder, since 14 maybe it is too much and the MOC
    might have holes --> Try coarser grid w apropiate resolution

    Args:
        table (_type_): _description_
        optCat (_type_): _description_

    Returns:
        _type_: _description_
    """
    coords = cat_params[f'{optCat}']['coord_cols']
    ra, dec = coords[0], coords[1]
    moc = MOC.from_lonlat(
                          table[f"{ra}"].T * u.deg,
                          table[f"{dec}"].T * u.deg,
                          max_norder=12,#14,
        )
    return moc

def get_ids_cat(
    xfile="DATA/pl26_srclist_v1.fits", 
    corename='pl26', 
    radius_deg=1,
    optCat=False,
    writeCat=False,
    out_pth='/storage/ivan/output_matchESO/',
    logCat=None,
    xmoc=None,
    acorr=False,
    mlcoords=False
):
    
    """
    Function that calls get_DES(...) to get the correspondant DES catalogue (see inputs) and create a MOC.
    This function write the results into a file.
    
    Input:
    ------
    xfile (string): path to the X-ray source list
    corename (string) 
    des1_main (bool): 
    des3_photo (bool): 
    out_pth (string)
    
    Output:
    -------
    
    """
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""
        
    t = Table.read(xfile)
        
    if(not('seqnum' in t.keys())):
        t['seqnum'] = range(len(t))

    #if acorr == False:
    rax, decx = np.median(t['RA']), np.median(t['DEC'])
    #else:
    #    rax  = np.median(t[f'RA{mltag}_CORR_{acorr}']) 
    #    decx = np.median(t[f'DEC{mltag}_CORR_{acorr}'])
        
    logger.info(f"RA={rax:.2f}, DEC={decx:.2f}")
        
    logger.info(f"\t+ Searching in {optCat} ...")
    ids = getCat( rax , decx , optCat=optCat, radius_deg=radius_deg ) 
    logger.info(f'\t Number of objects within {radius_deg} deg = {len(ids)}')
    
    if len(ids) == 0 and logCat!=None:
        logger.info(f"\t ** NO SRCS WITHIN {radius_deg} DEG **, writting in {logCat}")
        f = open(out_pth + f"{logCat}", "a")
        f.write(f"{corename}\n")
        f.close() 
        
    else:
        omoc = makeMOC_opt(ids, optCat)
        file = out_pth + f"{corename}_{optCat}.moc"
        logger.info(f"\t Saving {optCat} MOC in file: {file}")
        
        omoc.save(file, format='fits', overwrite=True)

        moc_intrsec = xmoc.intersection(omoc)

        if not(moc_intrsec.sky_fraction>0) and logCat!=None: #len(ids) == 0:
            logger.info(f"\t ** NO INTERSECTION MOC *, writting in {logCat}")
            f = open(out_pth + f"{logCat}", "a")
            f.write(f"{corename}\n")
            f.close()        

        if writeCat and moc_intrsec.sky_fraction>0:#len(ids) > 0:
            file = out_pth + f'{corename}_{optCat}.fits'
            logger.info(f"\t Saving {optCat} catalogue in file: {file}")
            ids.write(file, overwrite=True) 
    logger.info("\t Finish!\n\n\n")
        #moc_optical = makeMOC(rax, decx, radius_deg)
    return ids

def call_get_ids(
                 clst, 
                 radius_deg=1, 
                 cat_lst=False, 
                 writeCat=False, 
                 out_pth=None, 
                 in_pth=None,
                 logDict=False,
                 loop_clst=None,
                 n_clst=None,
                 xacorr=False,
                 mlcoords=False
):
    logger.info(f'*************************************************')
    logger.info(f'*** {loop_clst}/{n_clst} --> CLUSTER: {clst} ***')
    logger.info(f'*************************************************')
    
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""
    
    if xacorr == False:
        xray_src_file = in_pth + f'{clst}/SRC/srclist.fits'
    else:
        xray_src_file = in_pth + f'{clst}/srclist_cor{mltag}.fits'
    #else:
    #    xray_src_file = in_pth + f'{clst}/SRC/srclist.fits'
    xmoc = MOC.from_fits(in_pth + f'{clst}/MOC.MOC')
    logger.info(f'Using X-ray catalogue: {xray_src_file}')
    
    if logDict == None:
        log=None
        [ get_ids_cat(xfile=xray_src_file, 
                  corename=clst, 
                  radius_deg=radius_deg,
                  optCat=optCat,
                  writeCat=writeCat,
                  out_pth=out_pth,
                  logCat=log,
                  xmoc=xmoc,
                  acorr=xacorr,
                  mlcoords=mlcoords
                  )
        for optCat in cat_lst ]
    else:
        [ get_ids_cat(xfile=xray_src_file, 
                  corename=clst, 
                  radius_deg=radius_deg,
                  optCat=optCat,
                  writeCat=writeCat,
                  out_pth=out_pth,
                  logCat=logDict[f'{optCat}'],
                  xmoc=xmoc,
                  acorr=xacorr,
                  mlcoords=mlcoords
                  )
        for optCat in cat_lst ]