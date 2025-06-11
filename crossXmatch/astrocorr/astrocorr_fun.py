import os
import logging

from astroquery.vizier import Vizier

from astropy.coordinates import SkyCoord, Angle, Longitude, Latitude
from astropy import units as u

from mocpy import MOC

import numpy as np

from crossXmatch.config import *

import crossXmatch.ctp.get_catalogues as get_catalogues
import crossXmatch.ctp.get_queries as get_queries

logger = logging.getLogger(__name__)

def _get_offset(RA_clst, dec_clst, RA_ctlg, dec_ctlg):
    """
    """
    from astropy import coordinates
    from astropy.coordinates import SkyCoord
    
    c_clust = SkyCoord(RA_clst*u.deg,
                       dec_clst*u.deg)
    c_clust = c_clust.transform_to(coordinates.FK5())

    c = SkyCoord(RA_ctlg*u.deg, 
                 dec_ctlg*u.deg)
    c = c.transform_to(coordinates.FK5())
    
    offset = c_clust.separation(c)

    return offset

def _get_qsoCat_region(RA_clst, dec_clst, ctlg, survey, r_search=1):
    """
    Function that returns objects within a given r_search
    """
    coords = QSO_cat_params[f'{survey}']['coord_cols']
    ra, dec = coords[0], coords[1]
    
    RA_ctlg, dec_ctlg = ctlg[ra], ctlg[dec]
    offset = _get_offset(RA_clst, dec_clst, RA_ctlg, dec_ctlg)
    mask_offset = offset.value < r_search
    
    return ctlg[mask_offset]

def _call_getCat(RA_clst, dec_clst, optCat='LS', radius_deg=1):
    ctlg = get_catalogues.getCat(RA_clst, dec_clst, optCat=optCat, radius_deg=radius_deg)
    ctlg = get_queries._constrain_LS10(ctlg)
    
    mask_starlike = ctlg['type'] == "PSF"
    
    return ctlg[mask_starlike]

def healpix_radius(nside=16, hpidx=1):
    hp = HEALPix(nside=nside, order="ring", frame=FK5())   
     
    center = hp.healpix_to_skycoord(hpidx)
    vertex = hp.boundaries_skycoord(hpidx, step=1)  
    
    return center.separation(vertex).max()

def contains(moc, ra, dec, unit=u.deg):
    # moc can be a file instead of a moc
    if not isinstance(moc, MOC):
        moc = read(moc)
    try:
        _ = ra.unit, dec.unit
    except AttributeError:
        ra = ra << unit
        dec = dec << unit

    return moc.contains(ra, dec)

def sources_in_moc(sources, moc_file, racol="RA", deccol="Dec", unit=u.deg):
    """
    Return the sources in the astropy table `sources` included in `moc_file`.
    """
    ra = sources[racol]
    dec = sources[deccol]

    if ra.unit is None:
        ra = ra * unit
    else:
        ra = ra.quantity

    if dec.unit is None:
        dec = dec * unit
    else:
        dec = dec.quantity

    mask = contains(moc_file, ra, dec)

    return sources[mask]

def _get_sdss_data(ra0, dec0, radius=1, survey='sdss', photometry=False):
    """
    Get all sources from `survey` in the observed area of `healpix`.
    If `photometry` is ``True``, in addition to their positions, the
    corresponding magnitudes for each source are also downloaded.
    """
    
    columns = ['objID', 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS']
    if photometry:
        # Model magnitudes
        columns += ['umag', 'e_umag', 'gmag', 'e_gmag', 'rmag', 'e_rmag',
                    'imag', 'e_imag', 'zmag', 'e_zmag', 'q_mode']
    
    # Get all sources within the healpix cell
    #center = SkyCoord(hpixel.ra, hpixel.dec, unit="deg")
    center = SkyCoord(ra0, dec0, unit='deg')
    #radius = healpix_radius(hpixel.nside, hpixel.hpix)
    
    data = _dwnl_sdss(center, radius, columns, catalog='V/147/sdss12')
    # Select sources in the observed area of the healpix cell
    #moc_file = hpixel.paths.heal.joinpath(str(hpixel.hpix), "MOC.MOC")

    #data = mocs.sources_in_moc(
    #    data, moc_file, racol=kwargs["ra"], deccol=kwargs["dec"], unit=u.deg
    #)
    
    return data

def _dwnl_sdss(coords, radius, columns, catalog='V/147/sdss12', **kwargs):
    # Selects only primary sources
    v = Vizier(columns=columns, #column_filters={"mode":"=1"},
               column_filters={"mode": "=1", "e_RA_ICRS": ">0", "e_DE_ICRS": ">0"},
               row_limit=np.inf, timeout=6000)

    vrsp = v.query_region(coords, radius=radius, catalog=catalog)
    src_table = vrsp[0]
    
    # Select sources with clean photometry
    if 'q_mode' in columns:
        msk_good = src_table['q_mode'] == '+'
        src_table = src_table[msk_good]

    return src_table

def _call_src_in_moc(catalogue, survey, r_search, xmoc):
    coords = QSO_cat_params[f"{survey}"]["coord_cols"]
    ra, dec = coords[0], coords[1]
    
    logger.info(f"Number of sources within {r_search} deg: {len(catalogue)}")  
    catalogue = sources_in_moc(
        catalogue, xmoc, racol=ra, deccol=dec, unit=u.deg
    )
    logger.info(f"Number of sources within MOC: {len(catalogue)}")
    
    return catalogue

def save_mw_data(xRA, xDec, xmoc, ctlg, path_core_mw, clst, survey, r_search=1):
    path_opt_save = os.path.join(path_core_mw, f'{survey}/')
    if not os.path.exists(path_opt_save):
        os.makedirs(path_opt_save)

    path_opt_save = os.path.join(path_opt_save, f'{clst}.fits')
    
    if survey == 'sdss':
        opt_clst_region = _get_sdss_data(xRA, xDec, radius=r_search*u.deg, survey='sdss', photometry=True)
    elif survey == 'gaia':
        opt_clst_region = _get_qsoCat_region(xRA, xDec, ctlg, survey, r_search=r_search)
    elif survey =='LS':
        opt_clst_region = _call_getCat(xRA, xDec, optCat='LS', radius_deg=r_search)
    
    opt_clst_region = _call_src_in_moc(opt_clst_region, survey, r_search, xmoc)
    
    logger.info(f'Saving mw data in {path_opt_save}')
    opt_clst_region.write(path_opt_save, format='fits', overwrite=True)

