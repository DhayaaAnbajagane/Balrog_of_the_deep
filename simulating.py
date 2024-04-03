import logging
import shutil
import tempfile
import os

import numpy as np
import yaml
import joblib
import galsim
import fitsio
import healpy as hp #Needed to read extinction map
from esutil.ostools import StagedOutFile
from tqdm import tqdm

from files import (
    get_band_info_file,
    make_dirs_for_file,
    get_truth_catalog_path,
    expand_path)
from constants import MEDSCONF, R_SFD98
from truthing import make_coadd_grid_radec, make_coadd_random_radec, make_coadd_hexgrid_radec
from sky_bounding import get_rough_sky_bounds, radec_to_uv
from wcsing import get_esutil_wcs, get_galsim_wcs
from galsiming import render_sources_for_image, Our_params
from psf_wrapper import PSFWrapper
from realistic_dwarfing import init_dwarf_catalog, get_dwarf_object
from coadding import MakeSwarpCoadds

logger = logging.getLogger(__name__)

TMP_DIR = os.environ['TMPDIR']

class End2EndSimulation(object):
    """An end-to-end DES Y3 simulation.

    Parameters
    ----------
    seed : int
        The seed for the global RNG.
    output_meds_dir : str
        The output DEADATA/MEDS_DIR for the simulation data products.
    tilename : str
        The DES coadd tile to simulate.
    bands : str
        The bands to simulate.
    gal_kws : dict
        Keyword arguments to control the galaxy content of the simulation.
        Right now these should include:
            n_grid : int
                The galaxies will be put on a grid with `n_grid`
                on a side.
            g1 : float
                The true shear on the one-axis.
            g2 : float
                The true shear on the two-axis.
    psf_kws : dict
        Kyword arguments to control the PSF used for the simulation.
        Right now these should include:
            type : str
                One of 'gauss' and that's it.

    Methods
    -------
    run()
        Run the simulation, writing the data to disk.
    """
    def __init__(self, *,
                 seed, output_meds_dir, tilename, bands,
                 gal_kws, psf_kws, star_kws = None):
        
        self.output_meds_dir = output_meds_dir
        self.tilename = tilename
        self.bands = bands
        self.gal_kws  = gal_kws
        self.psf_kws  = psf_kws
        self.star_kws = star_kws
        self.seed = seed
        # any object within a 128 coadd pixel buffer of the edge of a CCD
        # will be rendered for that CCD
        self.bounds_buffer_uv = 128 * 0.263
       
        if self.psf_kws['type'] == 'psfex':
            self.draw_method = 'no_pixel'
        else:
            self.draw_method = 'phot'

        # make the RNGS. Extra initial seeds in case we need even more multiple random generators in future
        seeds = np.random.RandomState(seed=seed).randint(low=1, high=2**30, size=10)
        
        # one for galaxies in the truth catalog
        # one for noise in the images
        self.truth_cat_rng = np.random.RandomState(seed=seeds[0])
        self.noise_rng     = np.random.RandomState(seed=seeds[1])
        
        #one for drawing random objects from dwarf catalog
        self.dwarfsource_rng  = np.random.RandomState(seed=seeds[2])
        
        # load the image info for each band
        self.info = {}
        for band in bands:
            fname = get_band_info_file(
                meds_dir=self.output_meds_dir,
                medsconf=MEDSCONF,
                tilename=self.tilename,
                band=band)
            with open(fname, 'r') as fp:
                self.info[band] = yaml.load(fp, Loader=yaml.Loader)
        
    def run(self):
        """Run the simulation w/ galsim, writing the data to disk."""

        logger.info(' simulating coadd tile %s', self.tilename)

        # step 0 - Make coadd nwgint images
        self.info = MakeSwarpCoadds(tilename =  self.tilename, bands =  self.bands, output_meds_dir = self.output_meds_dir, config = np.NaN, n_files = None)._make_nwgint_files()
        
        # step 1 - Load simulated galaxy catalog if needed
        self.simulated_catalog = self._make_sim_catalog()
        
        # step 2 - make the truth catalog
        self.truth_catalog = self._make_truth_catalog()
        
        # step 3 - per band, write the images to a tile
        for band in self.bands:
            self._run_band(band=band)

    def _run_band(self, *, band):
        """Run a simulation of a truth cat for a given band."""

        logger.info(" rendering images in band %s", band)

        noise_seeds = self.noise_rng.randint(
            low=1, high=2**30, size=len(self.info[band]['src_info']))
        
        jobs = []
        for noise_seed, se_info in zip(noise_seeds, self.info[band]['src_info']):
            
            src_func = LazySourceCat(
                truth_cat=self.truth_catalog,
                wcs=get_galsim_wcs(
                    image_path=se_info['image_path'],
                    image_ext=se_info['image_ext']),
                psf=self._make_psf_wrapper(se_info=se_info),
                dwarfsource_rng = self.dwarfsource_rng,
                simulated_catalog = self.simulated_catalog,
                band = band)
            
            if self.gal_kws.get('galaxies', True):
                jobs.append(joblib.delayed(_render_se_image)(
                    se_info=se_info,
                    band=band,
                    truth_cat=self.truth_catalog,
                    bounds_buffer_uv=self.bounds_buffer_uv,
                    draw_method=self.draw_method,
                    noise_seed=noise_seed,
                    output_meds_dir=self.output_meds_dir,
                    src_func=src_func,
                    gal_kws = self.gal_kws))
            else:
                print("NO OBJECTS SIMULATED")
                jobs.append(joblib.delayed(_move_se_img_wgt_bkg)(se_info=se_info, output_meds_dir=self.output_meds_dir))

        #Have to do single threaded, cause dwarves are large objects (100s of arcsec)
        #so it is highly memory intensive to draw in galsim, on DECam images.
        with joblib.Parallel(n_jobs=40, backend='loky', verbose=50, max_nbytes=None) as p:
            p(jobs)

    def _make_psf_wrapper(self, *, se_info):
        
        wcs = get_galsim_wcs(image_path=se_info['image_path'], image_ext=se_info['image_ext'])

        if self.psf_kws['type'] == 'gauss':
            psf_model = galsim.Gaussian(fwhm=0.9)

        elif self.psf_kws['type'] == 'psfex':
            from galsim.des import DES_PSFEx
            psf_model = DES_PSFEx(expand_path(se_info['psfex_path']), wcs = wcs) #Need to pass wcs when reading file
            assert self.draw_method == 'phot'
        
        elif self.psf_kws['type'] == 'psfex_deconvolved':
            from psfex_deconvolved import PSFEx_Deconv
            psf_model = PSFEx_Deconv(expand_path(se_info['psfex_path']), wcs = wcs) #Need to pass wcs when reading file
            assert self.draw_method == 'phot' #Don't need no_pixel since psf already deconvolved
        
        else:
            raise ValueError(
                "psf type '%s' not recognized!" % self.psf_kws['type'])

        psf_wrap = PSFWrapper(psf_model, wcs)

        return psf_wrap

    def _make_truth_catalog(self):
        """Make the truth catalog."""
        # always done with first band
        band = self.bands[0]
        coadd_wcs = get_esutil_wcs(
            image_path=self.info[band]['image_path'],
            image_ext=self.info[band]['image_ext'])

        radius = 2*self.gal_kws['size_max']/0.263 #Radius of largest galaxy in pixel units. Factor of 2 to prevent overlap
        
        #These are the positions of the GALAXIES. We'll do the stars ourselves.
        ra_dwarf, dec_dwarf, x_dwarf, y_dwarf = make_coadd_hexgrid_radec(radius = radius,
            rng=self.truth_cat_rng, coadd_wcs=coadd_wcs,
            return_xy=True)
        
        #Don't inject a dwarf if it's center is in some masked part of the coadd
        if self.gal_kws.get('AvoidMaskedPix', True):
            
            #Get rid of galaxies in the masks.
            bit_mask = fitsio.read(self.info[band]['bmask_path'],  ext = self.info[band]['bmask_ext'])
            wgt      = fitsio.read(self.info[band]['weight_path'], ext = self.info[band]['weight_ext'])

            gal_mask = bit_mask[y_dwarf.astype(int), x_dwarf.astype(int)] == 0 #only select objects whose centers are unmasked
            gal_mask = gal_mask & (wgt[y_dwarf.astype(int), x_dwarf.astype(int)] != 0) #Do same thing but for wgt != 0 (nwgint sets wgt == 0 in some places)

            ra_dwarf, dec_dwarf = ra_dwarf[gal_mask], dec_dwarf[gal_mask]
            x_dwarf,  y_dwarf   = x_dwarf[gal_mask],  y_dwarf[gal_mask]
        
        print("TRUTH CATALOG HAS %d OBJECTS" % len(x_dwarf))
        
        
        #Find which dwarfs we will inject and subsample just the handful we need for this coadd
        mask_dwarf = self.simulated_catalog.cat['ISDIFFUSE'] == True
        mask_star  = self.simulated_catalog.cat['ISDIFFUSE'] == False
        inds_dwarf = self.dwarfsource_rng.choice(np.where(mask_dwarf)[0], len(ra_dwarf))
        
        #Positions and inds. The ind column also doubles as "DWARF or STAR" column since we can use it
        #to index into the simulated_catalog, which has this info.
        ra, dec, x, y, inds = [], [], [], [], []
        
        #Loop over each dwarf and build its ra and dec
        for d_i in tqdm(range(len(ra_dwarf)), desc = 'Building Truth catalog'):
            
            ind = inds_dwarf[d_i] #Find dwarf ind
            
            #Add dwarf properties to the inputs
            inds += [ind]
            ra   += [ra_dwarf[d_i]]
            dec  += [dec_dwarf[d_i]]
            x    += [x_dwarf[d_i]]
            y    += [y_dwarf[d_i]]
            
            dwarf_id = self.simulated_catalog.cat['PARENT_ID'][ind]
            inds_stars = np.where( (self.simulated_catalog.cat['PARENT_ID'] == dwarf_id) & mask_star)[0] #Find all stars associated with this
            Nstars   = len(inds_stars)
            #Get properties of the dwarf
            beta   = self.simulated_catalog.cat['beta'][ind]
            q      = self.simulated_catalog.cat['q'][ind]
            hlr    = self.simulated_catalog.cat['hlr'][ind]
            r0     = hlr / 1.6783469900166605/0.263 #Convert from hlr to scale radius of exponential, and in pixel units
            ang    = self.simulated_catalog.rand_rot[ind] * np.pi/180 #put it in radians
            
            
            if True:
                #Find star position using analytical inversion. We sample uniform dist. and convert to exponential dist.
                r = -r0 * np.log(1 - self.dwarfsource_rng.uniform(1e-16, 1 - 1e-4, Nstars)) #Use min != 0 else you will get r = infty error
                r = np.clip(r, 0, 5*hlr/0.263) #Prevent rare chances that star is placed crazy far away from the galaxy.

                #Get the x, y of the star. THIS IS W.R.T to galaxy. So x = 0 means center of dwarf. We will fix this later.
                t = self.dwarfsource_rng.uniform(0, 2*np.pi, Nstars) #Random angle. This is ang pos of star within galaxy
                x_tmp = r * np.cos(t)
                y_tmp = r * np.sin(t)

                #Now we change the position to account for galaxy ellipticity
                jac = galsim.Shear(beta = beta * galsim.degrees, q = q).getMatrix()

                x_star = x_tmp*jac[0, 0] + y_tmp*jac[0, 1]
                y_star = x_tmp*jac[1, 0] + y_tmp*jac[1, 1]

                #Finally, account for the new (random) galaxy rotation
                x_star_final = + x_star * np.cos(ang) + y_star * np.sin(ang)
                y_star_final = - x_star * np.sin(ang) + y_star * np.cos(ang)


                #Okay, now that all the rotation transforms are done, 
                #let's add back the host dwarf galaxy position (in image space)
                x_star_final += x_dwarf[d_i]
                y_star_final += y_dwarf[d_i]


                #Also need to get the corresponding RA and DEC now (skycoord). We already loaded the coadd_wcs before.
                #It's fine if a star goes "out of bounds" of coadd in x, y, ra, dec. We will remove out-of-bounds objs in later step
                ra_star_final, dec_star_final = coadd_wcs.image2sky(x_star_final, y_star_final)


                #Now just append everything to list of injection objects
                inds += list(inds_stars)
                ra   += list(ra_star_final)
                dec  += list(dec_star_final)
                x    += list(x_star_final)
                y    += list(y_star_final)
            
        
        inds = np.array(inds)
        ra   = np.array(ra)
        dec  = np.array(dec)
        x    = np.array(x)
        y    = np.array(y)
        
        dtype = [('number', 'i8'), ('ID', 'i8'), ('ind', 'i8'), 
                 ('ra',  'f8'), ('dec', 'f8'), ('x', 'f8'), ('y', 'f8'),
                 ('a_world', 'f8'), ('b_world', 'f8'), ('size', 'f8')]
        for b in self.bands:
            dtype += [('A%s'%b, 'f8')]
            
        truth_cat = np.zeros(len(ra), dtype = dtype)# + np.NaN
        
        
        #Now build truth catalog
        truth_cat['ind']    = inds
        truth_cat['number'] = np.arange(len(ra)).astype(np.int64) + 1 #This doesn't really matter
        truth_cat['ra']  = ra
        truth_cat['dec'] = dec
        truth_cat['x']   = x
        truth_cat['y']   = y
        
        truth_cat['ID']  = self.simulated_catalog.cat['ID'][truth_cat['ind']]
        
        if self.gal_kws['extinction'] == True:
            
            EBV   = hp.read_map(os.environ['EBV_PATH'])
            NSIDE = hp.npix2nside(EBV.size) 
            inds  = hp.ang2pix(NSIDE, ra, dec, lonlat = True)
            
            for b in self.bands: truth_cat['A%s' % b] = R_SFD98[b] * EBV[inds]
            
        
        #We shouldn't really need to use any of this since we won't be running
        #in "true-det" mode at any point.
        truth_cat['a_world'] = np.NaN
        truth_cat['b_world'] = np.NaN
        truth_cat['size']    = np.NaN
            
        truth_cat_path = get_truth_catalog_path(
            meds_dir=self.output_meds_dir,
            medsconf=MEDSCONF,
            tilename=self.tilename)

        make_dirs_for_file(truth_cat_path)
        fitsio.write(truth_cat_path, truth_cat, clobber=True)
        
        print("TRUTH CATALOG WRITTEN")

        return truth_cat
    

    def _make_sim_catalog(self):
        
        """Makes sim catalog"""
        
        self.simulated_catalog = init_dwarf_catalog(rng = self.dwarfsource_rng)
            
        return self.simulated_catalog


def _render_se_image(
        *, se_info, band, truth_cat, bounds_buffer_uv,
        draw_method, noise_seed, output_meds_dir, src_func, gal_kws):
    """Render an SE image.

    This function renders a full image and writes it to disk.

    Parameters
    ----------
    se_info : dict
        The entry from the `src_info` list for the coadd tile.
    band : str
        The band as a string.
    galaxy_truth_cat, star_truth_cat : np.ndarray
        A structured array (for galaxies and for stars) with the truth catalog. 
        Must at least have the columns 'ra' and 'dec' in degrees.
    bounds_buffer_uv : float
        The buffer in arcseconds for finding sources in the image. Any source
        whose center lies outside of this buffer area around the CCD will not
        be rendered for that CCD.
    draw_method : str
        The method used to draw the image. See the docs of `GSObject.drawImage`
        for details and options. Usually 'auto' is correct unless using a
        PSF with the pixel in which case 'no_pixel' is the right choice.
    noise_seed : int
        The RNG seed to use to generate the noise field for the image.
    output_meds_dir : str
        The output DEADATA/MEDS_DIR for the simulation data products.
    src_func : callable
        A function with signature `src_func(src_ind)` that
        returns the galsim object to be rendered and image position
        for a given index of the truth catalog.
    gal_kws : dict
        Dictionary containing the keywords passed to the
        the simulating code
    star_src_func : callable
        Similar to src_func, but for stars.
    """

    # step 1 - get the set of good objects for the CCD
    msk_inds = _cut_tuth_cat_to_se_image(
        truth_cat=truth_cat,
        se_info=se_info,
        bounds_buffer_uv=bounds_buffer_uv)

    # step 2 - render the objects
    im = _render_all_objects(
        msk_inds=msk_inds,
        truth_cat=truth_cat,
        se_info=se_info,
        band=band,
        src_func=src_func,
        draw_method=draw_method)
    
    # step 3 - add bkg and noise
    # also removes the zero point
    im, wgt, bkg, bmask = _add_noise_mask_background(
        image=im,
        se_info=se_info,
        noise_seed=noise_seed,
        gal_kws = gal_kws)

    # step 4 - write to disk
    _write_se_img_wgt_bkg(
        image=im,
        weight=wgt,
        background=bkg,
        bmask=bmask,
        se_info=se_info,
        output_meds_dir=output_meds_dir)


def _cut_tuth_cat_to_se_image(*, truth_cat, se_info, bounds_buffer_uv):
    """get the inds of the objects to render from the truth catalog"""
    wcs = get_esutil_wcs(
        image_path=se_info['image_path'],
        image_ext=se_info['image_ext'])
    sky_bnds, ra_ccd, dec_ccd = get_rough_sky_bounds(
        im_shape=se_info['image_shape'],
        wcs=wcs,
        position_offset=se_info['position_offset'],
        bounds_buffer_uv=bounds_buffer_uv,
        n_grid=4)
    u, v = radec_to_uv(truth_cat['ra'], truth_cat['dec'], ra_ccd, dec_ccd)
    sim_msk = sky_bnds.contains_points(u, v)
    msk_inds, = np.where(sim_msk)
    return msk_inds


def _render_all_objects(
        *, msk_inds, truth_cat, se_info, band, src_func, draw_method):
    gs_wcs = get_galsim_wcs(
        image_path=se_info['image_path'],
        image_ext=se_info['image_ext'])

    im = render_sources_for_image(
        image_shape=se_info['image_shape'],
        wcs=gs_wcs,
        draw_method=draw_method,
        src_inds=msk_inds,
        src_func=src_func,
        n_jobs=1)

    return im.array


def _add_noise_mask_background(*, image, se_info, noise_seed, gal_kws):
    """add noise, mask and background to an image, remove the zero point"""

    noise_rng = np.random.RandomState(seed=noise_seed)

    # first back to ADU units
    image /= se_info['scale']

    # now just read out these other images
    # in practice we just read out --> copy to other location
    # since balrog does not use wgt and bmask
    bkg   = fitsio.read(se_info['bkg_path'], ext=se_info['bkg_ext'])
#     wgt   = fitsio.read(se_info['weight_path'], ext=se_info['weight_ext'])
#     bmask = fitsio.read(se_info['bmask_path'], ext=se_info['bmask_ext'])

    wgt   = fitsio.read(se_info['nwgint_path'], ext=se_info['weight_ext'])
    bmask = fitsio.read(se_info['nwgint_path'], ext=se_info['bmask_ext'])
    
    
    #If we want Blank image, then we can't add original image
    if not gal_kws.get('BlankImage', False):
        # take the original image and add the simulated + original images together
        original_image = fitsio.read(se_info['nwgint_path'], ext=se_info['image_ext'])
        image += original_image

    else:
        
        # add the background
        image += bkg
        
        # now add noise
        img_std = 1.0 / np.sqrt(np.median(wgt[bmask == 0]))
        image += (noise_rng.normal(size=image.shape) * img_std)
    
    return image, wgt, bkg, bmask


def _write_se_img_wgt_bkg(
        *, image, weight, background, bmask, se_info, output_meds_dir):
    
    
#     # these should be the same
#     assert se_info['image_path'] == se_info['weight_path'], se_info
#     assert se_info['image_path'] == se_info['bmask_path'], se_info

#     # and not this
#     assert se_info['image_path'] != se_info['bkg_path']
    
    

    # get the final image file path and write
    image_file = se_info['nwgint_path'].replace(TMP_DIR, output_meds_dir)
    make_dirs_for_file(image_file)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        with StagedOutFile(image_file, tmpdir=tmpdir) as sf:
            # copy to the place we stage from
            shutil.copy(expand_path(se_info['nwgint_path']), sf.path)

            # open in read-write mode and replace the data
            with fitsio.FITS(sf.path, mode='rw') as fits:
                fits[se_info['image_ext']].write(image)
                fits[se_info['weight_ext']].write(weight)
                fits[se_info['bmask_ext']].write(bmask)

    # get the background file path and write
    bkg_file = se_info['bkg_path'].replace(TMP_DIR, output_meds_dir)
    make_dirs_for_file(bkg_file)
    with tempfile.TemporaryDirectory() as tmpdir:
        with StagedOutFile(bkg_file, tmpdir=tmpdir) as sf:
            # copy to the place we stage from
            shutil.copy(expand_path(se_info['bkg_path']), sf.path)

            # open in read-write mode and replace the data
            with fitsio.FITS(sf.path, mode='rw') as fits:
                fits[se_info['bkg_ext']].write(background)
                

def _move_se_img_wgt_bkg(*, se_info, output_meds_dir):
    '''
    Use this for blank image run where we do no source injection
    '''


    #Since nullweight is anyway made and transferred I dont think
    #we need any of this anymore
    
    '''
    # these should be the same
    assert se_info['image_path'] == se_info['weight_path'], se_info
    assert se_info['image_path'] == se_info['bmask_path'], se_info

    # and not this
    assert se_info['image_path'] != se_info['bkg_path']

    # get the final image file path and write
    image_file = se_info['image_path'].replace(TMP_DIR, output_meds_dir)
    make_dirs_for_file(image_file)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        with StagedOutFile(image_file, tmpdir=tmpdir) as sf:
            shutil.copy(expand_path(se_info['image_path']), sf.path)
    
    '''
    
    # get the background file path and write
    bkg_file = se_info['bkg_path'].replace(TMP_DIR, output_meds_dir)
    make_dirs_for_file(bkg_file)
    with tempfile.TemporaryDirectory() as tmpdir:
        with StagedOutFile(bkg_file, tmpdir=tmpdir) as sf:
            shutil.copy(expand_path(se_info['bkg_path']), sf.path)


class LazySourceCat(object):
    """A lazy source catalog that only builds objects to be rendered as they
    are needed.

    Parameters
    ----------
    truth_cat : structured np.array
        The truth catalog as a structured numpy array.
    wcs : galsim.GSFitsWCS
        A galsim WCS instance for the image to be rendered.
    psf : PSFWrapper
        A PSF wrapper object to use for the PSF.
    g1 : float
        The shear to apply on the 1-axis.
    g2 : float
        The shear to apply on the 2-axis.

    Methods
    -------
    __call__(ind)
        Returns the object to be rendered from the truth catalog at
        index `ind`.
    """
    def __init__(self, *, truth_cat, wcs, psf, band = None, dwarfsource_rng = None, simulated_catalog = None):
        self.truth_cat = truth_cat
        self.wcs = wcs
        self.psf = psf        
        
        self.dwarfsource_rng = dwarfsource_rng
        
        self.simulated_catalog = simulated_catalog
        self.band    = band
        
            

    def __call__(self, ind):
        pos = self.wcs.toImage(galsim.CelestialCoord(
            ra  = self.truth_cat['ra'][ind]  * galsim.degrees,
            dec = self.truth_cat['dec'][ind] * galsim.degrees))
        
        
        obj = get_dwarf_object(ind  = self.truth_cat['ind'][ind],
                               rng  = self.dwarfsource_rng, 
                               data = self.simulated_catalog,
                               band = self.band)

        #Now do extinction (the coefficients are just zero if we didnt set gal_kws['extinction'] = True)
        A_mag  = self.truth_cat[ind]['A%s' % self.band]
        A_flux = 10**(-A_mag/2.5)
        obj    = obj.withScaledFlux(A_flux)
        
        #Now psf
        psf = self.psf.getPSF(image_pos = pos)
        obj = galsim.Convolve([obj, psf], gsparams = Our_params)
        
        #For doing photon counting, need to do some workaround
        rng = galsim.BaseDeviate(self.dwarfsource_rng.randint(0, 2**16))
        
        return (obj, rng), pos