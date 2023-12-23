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
    
    fname = os.environ['CAT_PATH']
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
    
    flux = data.cat['FLUX_%s' % band.upper()][ind] #Assumed to be deredenned fluxes
    
    if data.cat['ISDIFFUSE'][ind] == False:
        
        prof = galsim.DeltaFunction(flux = flux, gsparams = Our_params)
        
    
    elif data.cat['ISDIFFUSE'][ind] == True:
    
        hlr = data.cat['hlr'][ind]
        e1  = data.cat['e1'][ind]
        e2  = data.cat['e2'][ind]
        
        prof  = galsim.Exponential(half_light_radius = hlr, flux = flux).shear(g1 = e1, g2 = e2)
        prof  = prof.rotate(data.rand_rot[ind]*galsim.degrees)
    

    return prof