import astropy.units as u
from mocpy import MOC
import numpy as np

def _get_mask_xray(xdata, prob=1e-7):
    mask      = np.logical_and(xdata["DET_EXT_FULL"] <= 0, xdata["DET_ML_FULL"] > 0)
    mask_soft = np.logical_and(xdata["DET_EXT_SOFT"] <= 0, xdata["DET_ML_SOFT"] > 0)
    mask_prob = np.logical_and(xdata["PROB_FULL"] >= 0, xdata["PROB_FULL"] <= prob)

    mask = np.logical_and(mask, mask_soft)
    # mask = np.logical_and(mask, mask_prob)

    return mask

def get_best_hp(hp_lst, moc_clst, pth_moc):
    """
    Function which returns the hp which is more complete
    """
    #print(hp_lst)
    for loop_hp, hp in enumerate(hp_lst):
        #print(hp)
        moc_hp = MOC.from_fits(pth_moc + hp + "/MOC.MOC")
        #print(pth_moc + "/MOC.MOC")
        actual_moc_intersec = moc_clst.intersection(moc_hp).sky_fraction
        
        if loop_hp == 0:
            max_moc_intersec = actual_moc_intersec
            hp_max_moc_intersec = hp
            
        if actual_moc_intersec > max_moc_intersec:
            max_moc_intersec = actual_moc_intersec
            hp_max_moc_intersec = hp
            
    return hp_max_moc_intersec

def lamRef_2_lamObs(l_ref, z):
    l_obs = l_ref*(1+z)
    return l_obs

def get_clst_data(df):
    r_500    = float(df.r_500)*1e-3 #Mpc
    ra_clst  = float(df.ra)
    dec_clst = float(df.dec)
    redshift = float(df.z)
    return r_500, ra_clst, dec_clst, redshift

def calc_xKC(z, gamma):
    """
    rest frame/observed: F_rf/F_o
    """
    exp = gamma - 2 
    return (1+z)**exp

def xflux_2_lx(cosmo, flux, redshift, gamma=1.4):
    kc = calc_xKC(redshift, gamma)
    d_L = ((cosmo.luminosityDistance(redshift)*u.Mpc/cosmo.h).to(u.cm)).value
    d_L = d_L**2
    
    return 4*np.pi*d_L*flux*kc