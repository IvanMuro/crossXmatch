import os
import logging
import subprocess

from astropy.io import fits
from astropy.table import Table, hstack, vstack
from astropy import units as u

from mocpy import MOC
import numpy as np

import pxsas

from crossXmatch.config import *
import crossXmatch.utils as utils

logger = logging.getLogger(__name__)


def __astrocor(self, survey):
    """
    Astrometric correction of the catalogue using SAS's epscorr.
    """
    odata_file = _astrocor_get_optcat(survey)
    xdata_file = _astrocor_get_xcat()
    offsets = _astrocor_calc_offset(xdata_file, odata_file)


def _astrocor_get_optcat(clst, survey, odata_file, tmp_path, mag="G"):
    """This funciton gets the optical catalogue, change the name of some columns
    to match SAS requirements (I suppose) and save the enw table in a
    temporary directory

    Args:
        survey (str, optional): _description_. Defaults to 'sdss'.

    Returns:
        _type_: _description_
    """
    # from .crossmatch import MWData
    # odata_file = MWData.get_data(self, survey=survey, photometry=True)

    odata = Table.read(odata_file)
    if survey == "sdss":
        odata["CAT_RADEC_ERR"] = np.fmax(
            odata["e_RA_ICRS"], odata["e_DE_ICRS"]
        )  # This returns an array w the maximum value compared for each array position e.g. p.fmax([2, 3, 4], [1, 5, 2]) --> array([ 2.,  5.,  4.])
    elif survey == 'gaia':
        odata["CAT_RADEC_ERR"] = astrom_error[f"{survey}"]
    elif survey == 'LS':
        odata["CAT_RADEC_ERR"] = astrom_error[f"{survey}"]
        
    odata["CAT_RADEC_ERR"] = odata["CAT_RADEC_ERR"].astype(np.float32)

    coords = QSO_cat_params[f"{survey}"]["coord_cols"]
    ra, dec = coords[0], coords[1]

    odata.rename_column(f"{ra}", "CAT_RA")
    odata.rename_column(f"{dec}", "CAT_DEC")

    mag = QSO_cat_params[survey]["mag"]

    if survey != 'LS':
        odata["rmag"] = odata[mag].astype(np.float32)
    else:
        odata["rmag"] = np.array(odata[mag]).astype(np.float32)
    odata.keep_columns(["CAT_RA", "CAT_DEC", "CAT_RADEC_ERR", "rmag"])
    id_col = QSO_cat_params[f"{survey}"]["id_col"]

    odata.meta["EXTNAME"] = "SRCLIST"  # f'{id_col}'

    path_tmp_save = os.path.join(tmp_path, f"{clst}/")
    if not os.path.exists(path_tmp_save):
        os.makedirs(path_tmp_save)
    odata_file_TMP = os.path.join(path_tmp_save, f"{survey}_EPOSCORR.fits")
    odata.write(odata_file_TMP, format="fits", overwrite=True)
    logger.info(f"Writing optical cat for eposcorr in {odata_file_TMP}")
    return odata_file_TMP

def _astrocor_get_xcat(path_core, tmp_path, clst, survey, mlcoords=False):
    """This function gets out extended sources and save the new table in a temporary directory

    Returns:
        _type_: _description_
    """
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""

    path_clst = os.path.join(path_core, "hp/", f"{clst}")

    xdata_file = os.path.join(path_clst, "SRC/srclist.fits")
    # xmoc_file = os.path.join(path_clst, 'MOC.MOC')

    xdata = Table.read(xdata_file)
    #xdata = sources_inmoc(xdata, xmoc_file, racol='RA', deccol='DEC')

    mask = utils._get_mask_xray(xdata)

    xdata.keep_columns([f"RA{mltag}", f"DEC{mltag}", f"RADEC{mltag}_ERR"])
    xdata[f"RADEC{mltag}_ERR"] = xdata[f"RADEC{mltag}_ERR"].astype(np.float32)
    xdata.meta["EXTNAME"] = "SRCLIST"

    path_tmp_save = os.path.join(tmp_path, f"{clst}/")
    if not os.path.exists(path_tmp_save):
        os.makedirs(path_tmp_save)
    xdata_file_TMP = os.path.join(
        path_tmp_save, f"srclist_EPOSCORR{mltag}_{survey}.fits"
    )

    xdata[mask].write(xdata_file_TMP, format="fits", overwrite=True)

    logger.info(f"Writing X-ray cat for eposcorr in {xdata_file_TMP}")
    return xdata_file_TMP


def get_mag_param_epos_survey(survey, mlcoords=False):
    if survey == "sdss":
        # minbmagn = 14
        # maxbmagn = 22
        # usebmagn = "no"
        # bmag = "no"
        minrmagn = 14
        maxrmagn = 22
        usermagn = "yes"
        rmag = "rmag"
    elif survey == "gaia":
        minrmagn = 17
        maxrmagn = 20.2
        usermagn = "yes"
        rmag = "rmag"
    elif survey == "LS" and mlcoords==False:
        minrmagn = 10
        maxrmagn = 22
        usermagn = "yes"
        rmag = "rmag"
    elif survey == "LS" and mlcoords==True:
        minrmagn = 10
        maxrmagn = 24
        usermagn = "yes"
        rmag = "rmag"
    return minrmagn, maxrmagn, usermagn, rmag


def _astrocor_calc_offset(
    xdata_file, odata_file, clst, survey, nsrcsmin=10, mlcoords=False
):
    # subname = '{}:{}: '.format(type(self).__name__,
    #                           sys._getframe().f_code.co_name)
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""
        
    if fits_file_len(xdata_file) < nsrcsmin:
        # There are not enough X-ray sources in the observed area of the
        # healpix for an accurate astrometric correction, so none is applied
        message = (
            "The number of X-ray sources in HEALpix %d is below %d. "
            "No astrometric corrections will be calculated."
        )
        logging.warning(message, self.HEALpix, nsrcsmin)
        offsets = {"RAOFFSET": 0.0, "DEOFFSET": 0.0, "ROT_CORR": 0.0}

    else:
        minrmagn, maxrmagn, usermagn, rmag = get_mag_param_epos_survey(survey,mlcoords=mlcoords)
        pxsas.run(
            "eposcorr",
            xrayset=xdata_file,
            opticalset=odata_file,
            findrotation="yes",
            niter=3,
            maxoffset=6.0,
            maxrotation=1.0,
            maxdist=5.0,
            maxposnerr=1.0,
            usemaxdist="yes",
            usebmagn="no",
            minrmagn=minrmagn,
            maxrmagn=maxrmagn,
            usermagn=usermagn,
            rmag=rmag,
            makeimage="no",
            opticalra="CAT_RA",
            opticaldec="CAT_DEC",
            opticalradecerr="RADEC_ERR",
            opticaltableext="SRCLIST",
            xrayra=f"RA{mltag}",
            xraydec=f"DEC{mltag}",
            xrayradecerr=f"RADEC{mltag}_ERR",
            xraytableext="SRCLIST",
            calculateoffsets="yes",
            withmatchtable="no",
            rawxsyserr=1.0,
            minxsyserr=0.5,
            maxsig=5.0,
        )

        header = fits.getheader(xdata_file, 1)
        offsets = {
            f"RAOFFSET": header["RAOFFSET"],
            f"DEOFFSET": header["DEOFFSET"],
            f"ROT_CORR": header["ROT_CORR"],
        }

        # with fits.open(xdata_file, 'update') as f:
        #    f[1].header[f'RAOFFSET{mltag}_{survey}'] = header[f'RAOFFSET']
        #    f[1].header[f'DEOFFSET{mltag}_{survey}'] = header[f'DEOFFSET']
        #    f[1].header[f'ROT_CORR{mltag}_{survey}'] = header[f'ROT_CORR']
        #    f[1].header[f'NMATCHES{mltag}_{survey}'] = header[f'NMATCHES']
    return offsets


def _astrocor_apply_offset(
    offsets, odata_file, clst, path_core, tmp_path, survey, mlcoords=False
):
    if mlcoords:
        mltag = "ML"
    else:
        mltag = ""

    path_clst = os.path.join(path_core, "hp/", f"{clst}")
    xdata_file = os.path.join(path_clst, "SRC/srclist.fits")

    # xdata_file = os.path.join(self.hpdirs['src'], 'srclist.fits')
    path_tmp_save = os.path.join(tmp_path, f"{clst}/")
    xdata_file_TMP = os.path.join(path_tmp_save, f"srclist_EPOSCORR2{mltag}_{survey}.fits")

    xdata = Table.read(xdata_file)
    xdata.keep_columns([f"RA{mltag}", f"DEC{mltag}", f"RADEC{mltag}_ERR"])
    xdata[f"RADEC{mltag}_ERR"] = xdata[f"RADEC{mltag}_ERR"].astype(np.float32)

    xdata.meta["EXTNAME"] = "SRCLIST"
    xdata.write(xdata_file_TMP, format="fits", overwrite=True)

    pxsas.run(
        "eposcorr",
        xrayset=xdata_file_TMP,
        opticalset=odata_file,
        xrayra=f"RA{mltag}",
        xraydec=f"DEC{mltag}",
        xrayradecerr=f"RADEC{mltag}_ERR",
        xraytableext="SRCLIST",
        opticaltableext="SRCLIST",
        calculateoffsets="no",
        raoffset=offsets[f"RAOFFSET"],
        decoffset=offsets[f"DEOFFSET"],
        rotation=offsets[f"ROT_CORR"],
    )

    xdata_offsets = Table.read(xdata_file_TMP)

    xdata_offsets.rename_column("RA_CORR", f"RA_CORR_{survey}")
    xdata_offsets.rename_column("DEC_CORR", f"DEC_CORR_{survey}")

    # racor  = xdata_offsets['RA_CORR']
    # deccor = xdata_offsets['DEC_CORR']
    racor = xdata_offsets[f"RA_CORR_{survey}"]
    deccor = xdata_offsets[f"DEC_CORR_{survey}"]

    xdata_cor_file = os.path.join(path_clst, f"srclist_cor{mltag}.fits")

    if os.path.exists(xdata_cor_file):
        xdata_cor = Table.read(xdata_cor_file)
    else:
        xdata_cor = Table.read(xdata_file)

    # idx = xdata_cor.colnames.index('RADEC_ERR')
    # xdata_cor.add_columns([racor, deccor], indexes=[idx + 1, idx + 1])
    idx = xdata_cor.colnames.index(f"RADEC{mltag}_ERR")
    xdata_cor.add_columns(
        [racor, deccor],
        indexes=[idx + 1, idx + 1],
        names=[f"RA{mltag}_CORR_{survey}", f"DEC{mltag}_CORR_{survey}"],
    )
    # xdata_cor.replace_column([racor, deccor], indexes=[idx + 1, idx + 1])
    xdata_cor.meta[f"RAOFFSET{mltag}_{survey}"] = offsets[f"RAOFFSET"]
    xdata_cor.meta[f"DEOFFSET{mltag}_{survey}"] = offsets[f"DEOFFSET"]
    xdata_cor.meta[f"ROT_CORR{mltag}_{survey}"] = offsets[f"ROT_CORR"]

    xdata_cor.write(xdata_cor_file, format="fits", overwrite=True)


def fits_file_len(fname, hdu=1):
    try:
        n = fits.getval(fname, "NAXIS2", hdu)
    except IOError:
        logging.warning("Empty file!!! Returning 0")
        n = 0
    return n

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
############################################# OLD! ##################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

def _astrocor_calc_offset_old_old(xdata_file, odata_file, clst, survey, nsrcsmin=10):
    # subname = '{}:{}: '.format(type(self).__name__,
    #                           sys._getframe().f_code.co_name)

    if fits_file_len(xdata_file) < nsrcsmin:
        # There are not enough X-ray sources in the observed area of the
        # healpix for an accurate astrometric correction, so none is applied
        message = (
            "The number of X-ray sources in HEALpix %d is below %d. "
            "No astrometric corrections will be calculated."
        )
        logging.warning(message, self.HEALpix, nsrcsmin)
        offsets = {"RAOFFSET": 0.0, "DEOFFSET": 0.0, "ROT_CORR": 0.0}

    else:
        logfile = f"./log/sas_{clst}_calc_off.log"
        pysas(
            "eposcorr",
            logfile=logfile,
            xrayset=xdata_file,
            opticalset=odata_file,
            findrotation="yes",
            niter=3,
            maxoffset=6.0,
            maxrotation=1.0,
            maxdist=5.0,
            maxposnerr=1.0,
            usemaxdist="yes",
            usebmagn="no",
            minrmagn=14,
            maxrmagn=22,
            usermagn="yes",
            rmag="rmag",
            makeimage="no",
            opticalra="CAT_RA",
            opticaldec="CAT_DEC",
            opticalradecerr="RADEC_ERR",
            opticaltableext="SRCLIST",
            xrayra="RA",
            xraydec="DEC",
            xrayradecerr="RADEC_ERR",
            xraytableext="SRCLIST",
            calculateoffsets="yes",
            withmatchtable="no",
            rawxsyserr=1.0,
            minxsyserr=0.5,
            maxsig=5.0,
        )  # , message=subname

        header = fits.getheader(xdata_file, 1)
        offsets = {
            f"RAOFFSET": header["RAOFFSET"],
            f"DEOFFSET": header["DEOFFSET"],
            f"ROT_CORR": header["ROT_CORR"],
        }

    return offsets


def _astrocor_calc_offset_old(xdata_file, odata_file, clst, survey, nsrcsmin=10):
    # subname = '{}:{}: '.format(type(self).__name__,
    #                           sys._getframe().f_code.co_name)

    if fits_file_len(xdata_file) < nsrcsmin:
        # There are not enough X-ray sources in the observed area of the
        # healpix for an accurate astrometric correction, so none is applied
        message = (
            "The number of X-ray sources in HEALpix %d is below %d. "
            "No astrometric corrections will be calculated."
        )
        logging.warning(message, self.HEALpix, nsrcsmin)
        offsets = {"RAOFFSET": 0.0, "DEOFFSET": 0.0, "ROT_CORR": 0.0}

    else:
        logfile = f"./log/sas_{clst}_calc_off.log"
        # pysas(
        pxsas.run(
            "eposcorr",
            logfile=logfile,
            xrayset=xdata_file,
            opticalset=odata_file,
            findrotation="yes",
            niter=3,
            maxoffset=6.0,
            maxrotation=1.0,
            maxdist=5.0,
            maxposnerr=1.0,
            usemaxdist="yes",
            usebmagn="no",
            minrmagn=14,
            maxrmagn=22,
            usermagn="yes",
            rmag="rmag",
            makeimage="no",
            opticalra="CAT_RA",
            opticaldec="CAT_DEC",
            opticalradecerr="RADEC_ERR",
            opticaltableext="SRCLIST",
            xrayra="RA",
            xraydec="DEC",
            xrayradecerr="RADEC_ERR",
            xraytableext="SRCLIST",
            calculateoffsets="yes",
            withmatchtable="no",
            rawxsyserr=1.0,
            minxsyserr=0.5,
            maxsig=5.0,
        )  # , message=subname

        header = fits.getheader(xdata_file, 1)
        offsets = {
            f"RAOFFSET": header["RAOFFSET"],
            f"DEOFFSET": header["DEOFFSET"],
            f"ROT_CORR": header["ROT_CORR"],
        }

    return offsets


def pysas_old(task, message="", getoutput=False, logfile=None, **kwargs):
    """
    Wrapper for SAS task. 'task' must be the name of a SAS task, and kwargs
    contains all the parameters passed to the task, with the same name.
    """

    options = " ".join([f'{key}="{kwargs[key]}"' for key in kwargs.keys()])
    command = " ".join([task, options])

    # logging.debug(''.join([message, command]))

    if getoutput:
        log = subprocess.check_output(command, shell=True)
    else:
        if logfile is not None:
            command = "".join([command, f">> {logfile} 2>&1 "])
        # subprocess.check_output(command, shell=True)
        subprocess.check_call(command, shell=True)
        log = None

    return log


def pysas(task, message="", getoutput=False, logfile=None, **kwargs):
    """
    Wrapper for SAS task. 'task' must be the name of a SAS task, and kwargs
    contains all the parameters passed to the task, with the same name.
    """

    # options = ' '.join([f'{key}="{kwargs[key]}"'
    #                    for key in kwargs.keys()])

    options = " ".join([f'{k}="{v}"' for k, v in kwargs.items()])
    command = " ".join([task, options])

    # logging.debug(''.join([message, command]))

    if getoutput:
        log = subprocess.check_output(command, shell=True)
    else:
        if logfile is not None:
            command = "".join([command, f">> {logfile} 2>&1 "])
        # subprocess.check_output(command, shell=True)
        subprocess.check_call(command, shell=True)
        log = None

    return log


def _astrocor_apply_offset_old(offsets, odata_file, clst, path_core, tmp_path, survey):
    path_clst = os.path.join(path_core, "hp/", f"{clst}")
    xdata_file = os.path.join(path_clst, "SRC/srclist.fits")

    # xdata_file = os.path.join(self.hpdirs['src'], 'srclist.fits')
    path_tmp_save = os.path.join(tmp_path, f"{clst}/")
    xdata_file_TMP = os.path.join(path_tmp_save, f"srclist_EPOSCORR2_{survey}.fits")

    xdata = Table.read(xdata_file)
    xdata.keep_columns(["RA", "DEC", "RADEC_ERR"])
    xdata["RADEC_ERR"] = xdata["RADEC_ERR"].astype(np.float32)

    xdata.meta["EXTNAME"] = "SRCLIST"
    xdata.write(xdata_file_TMP, format="fits", overwrite=True)

    logfile = f"./log/sas_{clst}_apply_off.log"
    # pysas(
    pxsas.run(
        "eposcorr",
        logfile=logfile,
        xrayset=xdata_file_TMP,
        opticalset=odata_file,
        xrayra="RA",
        xraydec="DEC",
        xrayradecerr="RADEC_ERR",
        xraytableext="SRCLIST",
        opticaltableext="SRCLIST",
        calculateoffsets="no",
        raoffset=offsets[f"RAOFFSET"],
        decoffset=offsets[f"DEOFFSET"],
        rotation=offsets[f"ROT_CORR"],
    )  # message=subname,

    xdata_offsets = Table.read(xdata_file_TMP)

    xdata_offsets.rename_column("RA_CORR", f"RA_CORR_{survey}")
    xdata_offsets.rename_column("DEC_CORR", f"DEC_CORR_{survey}")

    # racor  = xdata_offsets['RA_CORR']
    # deccor = xdata_offsets['DEC_CORR']
    racor = xdata_offsets[f"RA_CORR_{survey}"]
    deccor = xdata_offsets[f"DEC_CORR_{survey}"]

    xdata_cor_file = os.path.join(path_clst, f"srclist_cor.fits")

    if os.path.exists(xdata_cor_file):
        xdata_cor = Table.read(xdata_cor_file)
    else:
        xdata_cor = Table.read(xdata_file)

    idx = xdata_cor.colnames.index("RADEC_ERR")
    xdata_cor.add_columns([racor, deccor], indexes=[idx + 1, idx + 1])
    # xdata_cor.replace_column([racor, deccor], indexes=[idx + 1, idx + 1])
    xdata_cor.meta[f"RAOFFSET_{survey}"] = offsets[f"RAOFFSET"]
    xdata_cor.meta[f"DEOFFSET_{survey}"] = offsets[f"DEOFFSET"]
    xdata_cor.meta[f"ROT_CORR_{survey}"] = offsets[f"ROT_CORR"]

    xdata_cor.write(xdata_cor_file, format="fits", overwrite=True)
