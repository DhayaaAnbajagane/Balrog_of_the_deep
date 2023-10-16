import os
import functools
import collections
import numpy as np

import galsim
import ngmix
import fitsio

from galsiming import Our_params

WLDeblendData = collections.namedtuple(
    'WLDeblendData',
    [
        'cat', 'rand_rot', 'survey_name', 'bands', 'surveys',
        'builders', 'total_sky', 'noise', 'ngal_per_arcmin2',
        'psf_fwhm', 'pixel_scale',
    ],
)

DESDFData = collections.namedtuple(
    'DESDFData',
    [
        'cat', 'rand_rot',
    ],
)


def _cached_DESDF_catalog_read():
    fname = os.environ['CATDESDF_PATH']
    return fitsio.read(fname)


def init_desdf_catalog(*, rng):
    
    DESDF_cat = _cached_DESDF_catalog_read()
    
    #If rng not supplied then don't do random rotation
    if rng is None:
        angle = None
    else:
        angle = rng.uniform(low = 0, high = 1, size = len(DESDF_cat))*360
        
    return DESDFData(DESDF_cat, angle)


def get_desdf_galaxy(*, desdf_ind, rng, data, band = None):
    """Draw a galaxy from the DESDF model.

    Parameters
    ----------
    desdf_ind : int
        Index of galaxy in desdf catalog. Needed so galaxy in
        every band/exposure looks the same.
    rng : np.random.RandomState
        An RNG to use for galaxy orientation
    data : DESDF_cat data fitsio file
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.
    band : character
        A single character containing the band of the galaxy.
        If None then average over riz bands.

    Returns
    -------
    gal : galsim Object
        The galaxy as a galsim object.
    """
    
#     bulge_frac = data.cat['BDF_FRACDEV'][desdf_ind] #Fraction of bulge to total
    
    if band == None:
        flux = np.sum([data.cat['FLUX_%s' % i.upper()][desdf_ind] for i in ['r', 'i', 'z']])
        print("Why am I in here", band)
    else:
        flux = data.cat['FLUX_%s' % band.upper()][desdf_ind] #Assumed to be deredenned fluxes in DF catalog
        
#     disk  = galsim.Exponential(flux = flux,   half_light_radius = data.cat['BDF_HLR'][desdf_ind], gsparams = Our_params)
#     bulge = galsim.DeVaucouleurs(flux = flux, half_light_radius = data.cat['BDF_HLR'][desdf_ind], gsparams = Our_params)
    
#     prof  = bulge_frac*bulge + (1 - bulge_frac)*disk
#     prof  = prof.shear(g1 = data.cat['BDF_G1'][desdf_ind], g2 = data.cat['BDF_G2'][desdf_ind])
    
    f  = data.cat['BDF_FRACDEV'][desdf_ind]
    g1 = data.cat['BDF_G1'][desdf_ind]
    g2 = data.cat['BDF_G2'][desdf_ind]
    T  = data.cat['BDF_T'][desdf_ind]
    
    prof  = ngmix.GMixModel([0, 0, g1, g2, T, f, flux], "bdf").make_galsim_object(gsparams = Our_params)
    prof  = prof.rotate(data.rand_rot[desdf_ind]*galsim.degrees)

    return prof


def get_psf_config_wldeblend(*, data):
    """Get a config dict for a the PSF model for the weak lensing deblending
    objects.

    Parameters
    ----------
    data : WLDeblendData
        Namedtuple with data for making galaxies via the weak lesning
        deblending package.

    Returns
    -------
    gs_config : dict
        A dictionary with the PSF info.
    """
    gs_config = {}
    gs_config["type"] = "Kolmogorov"
    gs_config["fwhm"] = data.psf_fwhm
    return gs_config