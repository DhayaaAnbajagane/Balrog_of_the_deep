import os
import functools
import collections
import numpy as np

import galsim
import ngmix
import fitsio

from galsiming import Our_params


Dwarf_data = collections.namedtuple(
    'Dwarf_data',
    [
        'cat', 'rand_rot',
    ],
)

def init_dwarf_catalog(*, rng):
    
    fname = os.environ['CATDESDF_PATH']
    Dwarf_data = fitsio.read(fname)
    
    #If rng not supplied then don't do random rotation
    if rng is None:
        angle = None
    else:
        angle = rng.uniform(low = 0, high = 1, size = len(DESDF_cat))*360
        
    return (Dwarf_data, angle)


def get_dwarf_object(*, ind, rng, data, band = None):
    """Draw a galaxy from the DESDF model.

    Parameters
    ----------
    ind : int
        Index of object in dwarf catalog.
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
    
    flux = data.cat['FLUX_%s' % band.upper()][desdf_ind] #Assumed to be deredenned fluxes
    
    if data.cat['ISDWARF'][ind] == False:
            
        prof = galsim.DeltaFunction(flux=flux, gsparams = Our_params)
        
    elif data.cat['ISDWARF'][ind] == TRUE:
        
    
        f  = data.cat['BDF_FRACDEV'][desdf_ind]
        g1 = data.cat['BDF_G1'][desdf_ind]
        g2 = data.cat['BDF_G2'][desdf_ind]
        T  = data.cat['BDF_T'][desdf_ind]

        prof  = ngmix.GMixModel([0, 0, g1, g2, T, f, flux], "bdf").make_galsim_object(gsparams = Our_params)
        prof  = prof.rotate(data.rand_rot[desdf_ind]*galsim.degrees)

    return prof