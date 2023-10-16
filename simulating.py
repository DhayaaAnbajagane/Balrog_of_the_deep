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
from realistic_galaxying import init_desdf_catalog, get_desdf_galaxy
from realistic_starsing import init_lsst_starsim_catalog
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
        
        #one for drawing random galaxies from descwl package
        self.galsource_rng  = np.random.RandomState(seed=seeds[2])
        self.starsource_rng = np.random.RandomState(seed=seeds[3])
        
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
        self.galaxy_simulated_catalog = self._make_sim_catalog()
        
        # step 2 - make the truth catalog
        self.galaxy_truth_catalog = self._make_truth_catalog()
        
        
        #step 2b - load simulated, truth star catalog
        if self.star_kws['stars'] == True:
            self.star_truth_catalog, self.star_simulated_catalog = self._make_star_catalogs()
        else:
            self.star_truth_catalog, self.star_simulated_catalog = None, None
        

        # step 3 - per band, write the images to a tile
        for band in self.bands:
            self._run_band(band=band)

    def _run_band(self, *, band):
        """Run a simulation of a truth cat for a given band."""

        logger.info(" rendering images in band %s", band)

        noise_seeds = self.noise_rng.randint(
            low=1, high=2**30, size=len(self.info[band]['src_info']))
        
        jobs = []
        for noise_seed, se_info in zip(
                noise_seeds, self.info[band]['src_info']):

            galaxy_src_func = LazySourceCat(
                truth_cat=self.galaxy_truth_catalog,
                wcs=get_galsim_wcs(
                    image_path=se_info['image_path'],
                    image_ext=se_info['image_ext']),
                psf=self._make_psf_wrapper(se_info=se_info),
                gal_mag = self.gal_kws['gal_mag'],
                gal_source = self.gal_kws['gal_source'],
                galsource_rng = self.galsource_rng,
                simulated_catalog = self.galaxy_simulated_catalog,
                band = band)
            
            if self.star_kws['stars'] == True:
                star_src_func = LazyStarSourceCat(
                    truth_cat=self.star_truth_catalog,
                    wcs=get_galsim_wcs(
                        image_path=se_info['image_path'],
                        image_ext=se_info['image_ext']),
                    psf=self._make_psf_wrapper(se_info=se_info),
                    star_mag = self.star_kws['star_mag'],
                    star_source = self.star_kws['star_source'],
                    starsource_rng = self.starsource_rng,
                    simulated_catalog = self.star_simulated_catalog,
                    band = band)
            else:
                star_src_func = None

            if self.gal_kws.get('galaxies', True):
                jobs.append(joblib.delayed(_render_se_image)(
                    se_info=se_info,
                    band=band,
                    galaxy_truth_cat=self.galaxy_truth_catalog,
                    star_truth_cat=self.star_truth_catalog,
                    bounds_buffer_uv=self.bounds_buffer_uv,
                    draw_method=self.draw_method,
                    noise_seed=noise_seed,
                    output_meds_dir=self.output_meds_dir,
                    galaxy_src_func=galaxy_src_func,
                    star_src_func = star_src_func,
                    gal_kws = self.gal_kws))
            else:
                print("NO GALAXY SIMULATED")
                jobs.append(joblib.delayed(_move_se_img_wgt_bkg)(se_info=se_info, output_meds_dir=self.output_meds_dir))

        with joblib.Parallel(
                n_jobs=-1, backend='loky', verbose=50, max_nbytes=None) as p:
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

        #Set what type of galaxy counts we use
        #Either constant counts per tile
        #or draw from poisson
        print(self.gal_kws)
        print(self.psf_kws)

        n_grid = self.gal_kws['n_grid']
        n_gal  = n_grid**2
        
        radius = 2*self.gal_kws['size_max']/0.263 #Radius of largest galaxy in pixel units. Factor of 2 to prevent overlap
        #Set what type of grid we use
        if self.gal_kws['truth_type'] in ['hexgrid', 'hexgrid-truedet']:
            ra, dec, x, y = make_coadd_hexgrid_radec(radius = radius,
                rng=self.truth_cat_rng, coadd_wcs=coadd_wcs,
                return_xy=True)
            
        elif self.gal_kws['truth_type'] in ['grid', 'grid-truedet']:
            ra, dec, x, y = make_coadd_grid_radec(
                rng=self.truth_cat_rng, coadd_wcs=coadd_wcs,
                return_xy=True, n_grid=n_grid)
            
        else:
            raise ValueError("Invalid option for `truth_type`. Use 'hexgrid', 'grid', 'random', 'grid-truedet', or 'random-truedet'.")
            
        if self.gal_kws.get('AvoidMaskedPix', True):
            
            #Get rid of galaxies in the masks.
            bit_mask = fitsio.read(self.info[band]['bmask_path'],  ext = self.info[band]['bmask_ext'])
            wgt      = fitsio.read(self.info[band]['weight_path'], ext = self.info[band]['weight_ext'])

            gal_mask = bit_mask[y.astype(int), x.astype(int)] == 0 #only select objects whose centers are unmasked
            gal_mask = gal_mask & (wgt[y.astype(int), x.astype(int)] != 0) #Do same thing but for wgt != 0 (nwgint sets wgt == 0 in some places)

            ra, dec = ra[gal_mask], dec[gal_mask]
            x,  y   = x[gal_mask],  y[gal_mask]
        
        print("TRUTH CATALOG HAS %d OBJECTS" % len(x))
        
        
        dtype = [('number', 'i8'), ('ID', 'i8'), ('ind', 'i8'), ('inj_class', 'i4'), 
                 ('ra',  'f8'), ('dec', 'f8'), ('x', 'f8'), ('y', 'f8'),
                 ('a_world', 'f8'), ('b_world', 'f8'), ('size', 'f8')]
        for b in self.bands:
            dtype += [('A%s'%b, 'f8')]
            
        truth_cat = np.zeros(len(ra), dtype = dtype)# + np.NaN
        
        
        #Use i-band as reference magnitude
        mag_ref = 30 - 2.5*np.log10(self.simulated_catalog.cat['FLUX_I'])
        
        # A value of 23.0 returns something similar to Balrog Y3
        # A value of 23.5 maybe is similar to the WL sample in Y6.
        # A value of 21.5 is maybe optimal for LSS samples in Y6.
        
        all_rand_inds = self.galsource_rng.randint(low=0, high=len(self.simulated_catalog.cat), size=len(ra))
        wl_rand_inds  = mock_balrog_sigmoid(mag_ref, 23.5, self.galsource_rng)
        
        #Loop until you get enough high-z galaxies
        #Just a safety loop so code doesn't fail because we
        #somehow selected too few galaxies.
        while True:
            wl_HQz_inds   = mock_balrog_sigmoid(mag_ref, 23.5, self.galsource_rng)
            wl_HQz_inds   = wl_HQz_inds[self.simulated_catalog.cat['Z_SOURCE'][wl_HQz_inds] != 0] #Select only indices with high-redshifts
            
            #Z_SOURCE HAS FOLLOWING INDICES
            # Z_SOURCE = 0 --- "NO REDSHIFT"
            # Z_SOURCE = 1 --- "SPEC-Z"
            # Z_SOURCE = 2 --- "COSMOS2020"
            # Z_SOURCE = 3 --- "PAUS+COSMOS"
            # Z_SOURCE = 4 ---  "C3R2"

            #If we have enough then break out
            if len(wl_HQz_inds) > len(ra)//4: 
                print("I HAVE ENOUGH HIGHQ GALAXIES. BREAKING OUT NOW")
                break
        
        
        inds = np.zeros_like(all_rand_inds) #Array to hold deep field index of galaxy
        cls  = np.zeros_like(all_rand_inds) #Array to hold injection class
        
        quarter = len(ra)//4 #Number of gals that make one quarter of required galaxies
        
        #First 1/2 is random balrog
        inds[:2*quarter] = all_rand_inds[:2*quarter];  cls[:2*quarter] = 0;
        
        #Next 1/4 is random WL specific
        inds[2*quarter:3*quarter] = wl_rand_inds[:quarter]; cls[2*quarter:3*quarter] = 1;
        
        #Final 1/4 is WL with high quality redshift data
        inds[3*quarter:4*quarter] = wl_HQz_inds[:quarter]; cls[3*quarter:4*quarter] = 2;
        
        
        #Now just randomly shuffle the inds. The ra/dec are in a uniform order.
        #If you don't shuffle inds then galaxies on one side of image will be 1/2 balrog
        #and remaining would be WL specific. We don't want that
        shuffle_inds = np.arange(len(inds))
        self.galsource_rng.shuffle(shuffle_inds)
        
        inds = inds[shuffle_inds]
        truth_cat['inj_class'] = cls[shuffle_inds]
        
        
        #Now build truth catalog
        truth_cat['ind']    = inds
        truth_cat['number'] = np.arange(len(ra)).astype(np.int64) + 1
        truth_cat['ra']  = ra
        truth_cat['dec'] = dec
        truth_cat['x'] = x
        truth_cat['y'] = y
        
        truth_cat['ID'] = self.simulated_catalog.cat['ID'][truth_cat['ind']]
        
        if self.gal_kws['extinction'] == True:
            
            EBV   = hp.read_map(os.environ['EBV_PATH'])
            NSIDE = hp.npix2nside(EBV.size) 
            inds  = hp.ang2pix(NSIDE, ra, dec, lonlat = True)
            
            for b in self.bands: truth_cat['A%s' % b] = R_SFD98[b] * EBV[inds]
            
        
        g1 = self.simulated_catalog.cat['BDF_G1'][truth_cat['ind']]
        g2 = self.simulated_catalog.cat['BDF_G2'][truth_cat['ind']]
        q  = np.sqrt(g1**2 + g2**2)

        truth_cat['a_world'] = 1
        truth_cat['b_world'] = q
        truth_cat['size']    = np.sqrt(self.simulated_catalog.cat['BDF_T'][truth_cat['ind']])
            
        truth_cat_path = get_truth_catalog_path(
            meds_dir=self.output_meds_dir,
            medsconf=MEDSCONF,
            tilename=self.tilename)

        make_dirs_for_file(truth_cat_path)
        fitsio.write(truth_cat_path, truth_cat, clobber=True)

        return truth_cat

    def _make_sim_catalog(self):
        
        """Makes sim catalog"""
        
        self.simulated_catalog = init_desdf_catalog(rng = self.galsource_rng)
            
        mag_i = 30 - 2.5*np.log10(self.simulated_catalog.cat['FLUX_I'])
        hlr   = np.sqrt(self.simulated_catalog.cat['BDF_T']) #This is only approximate
        Mask  = ((mag_i > self.gal_kws['mag_min']) &  (mag_i < self.gal_kws['mag_max']) &
                 (hlr > self.gal_kws['size_min'])  &  (hlr < self.gal_kws['size_max'])
                )

        self.simulated_catalog = self.simulated_catalog._replace(cat = self.simulated_catalog.cat[Mask])

        if self.gal_kws['circular']:
            #Temporarily remove all ellipticity
            self.simulated_catalog.cat['BDF_G1'] = 0
            self.simulated_catalog.cat['BDF_G2'] = 0

        return self.simulated_catalog
    
    
    def _make_star_catalogs(self):
        
        """Makes sim catalog and truth catalog at same time"""
        
        # always done with first band
        band = self.bands[0]
        coadd_wcs = get_esutil_wcs(
            image_path=self.info[band]['image_path'],
            image_ext=self.info[band]['image_ext'])
        coadd_info = self.info[band]
        
        if self.star_kws['star_source'] in ['lsst_sim']:

            star_catalog, binary_catalog = init_lsst_starsim_catalog(rng = self.starsource_rng)

            star_inds   = _cut_tuth_cat_to_se_image(truth_cat=star_catalog,   se_info=coadd_info, bounds_buffer_uv=self.bounds_buffer_uv)
            binary_inds = _cut_tuth_cat_to_se_image(truth_cat=binary_catalog, se_info=coadd_info, bounds_buffer_uv=self.bounds_buffer_uv)
            
            star_upsample   = 1000 #factor because we didn't download all stars
            binary_upsample = 10 * 100 #factor because we didn't download all binaries, and only 1/10th of binaries were simulated
            
            star_num    = self.star_kws['upscale'] * star_upsample   * len(star_inds)   * (1 - self.star_kws['f_bin'])
            binary_num  = self.star_kws['upscale'] * binary_upsample * len(binary_inds) * self.star_kws['f_bin']
            star_inds   = star_inds[self.starsource_rng.randint(len(star_inds),     size = int(star_num))]
            binary_inds = binary_inds[self.starsource_rng.randint(len(binary_inds), size = int(binary_num))]
            
            n_stars = len(star_inds)
            n_binar = len(binary_inds)
            n_tot   = n_stars + n_binar
            ra, dec, x, y = make_coadd_random_radec(rng=self.truth_cat_rng, coadd_wcs=coadd_wcs, return_xy=True, n_gal=n_tot)
            
            star_catalog   = star_catalog[star_inds]
            binary_catalog = binary_catalog[binary_inds]
            
            #Fill ra/dec for single star catalog
            star_catalog['ra']  = ra[:n_stars]
            star_catalog['dec'] = dec[:n_stars]
            
            #Offsets between two stars in binary system. 
            #Generated in flat-sky. Convert to curved sky for ra_offset. Dec offset is fine
            angles = self.starsource_rng.random(len(binary_inds))*np.pi
            sep    = (binary_catalog['a']*2.25461e-8) / (10**(1 + binary_catalog['mu0']/5)) * (180/np.pi) #conv. Rsun to pc, then rad to deg
            cos    = np.cos(angles) 
            sin    = np.sin(angles)
            
            ra_offset  = cos*sep / np.cos(dec[n_stars:]*np.pi/180)
            dec_offset = sin*sep
            
            
            binarystar1_catalog = np.zeros(len(binary_inds), dtype = star_catalog.dtype)
            binarystar1_catalog['mag_g'] = binary_catalog['mag_g_1']
            binarystar1_catalog['mag_r'] = binary_catalog['mag_r_1']
            binarystar1_catalog['mag_i'] = binary_catalog['mag_i_1']
            binarystar1_catalog['mag_z'] = binary_catalog['mag_z_1']
            binarystar1_catalog['ra']    = ra[n_stars:]
            binarystar1_catalog['dec']   = dec[n_stars:]
            
            binarystar2_catalog = np.zeros(len(binary_inds), dtype = star_catalog.dtype)
            binarystar2_catalog['mag_g'] = binary_catalog['mag_g_2']
            binarystar2_catalog['mag_r'] = binary_catalog['mag_r_2']
            binarystar2_catalog['mag_i'] = binary_catalog['mag_i_2']
            binarystar2_catalog['mag_z'] = binary_catalog['mag_z_2']
            binarystar2_catalog['ra']    = ra[n_stars:]  + ra_offset
            binarystar2_catalog['dec']   = dec[n_stars:] + dec_offset
            
            simulated_cat = np.concatenate([star_catalog, binarystar1_catalog, binarystar2_catalog])
            
            Mask = ((simulated_cat['mag_i'] > self.star_kws['mag_min']) & 
                    (simulated_cat['mag_i'] < self.star_kws['mag_max']))
            
            simulated_cat = simulated_cat[Mask]
        
        #NOW DO TRUTH CATALOG PART
        truth_cat = np.zeros(len(simulated_cat), dtype=[('number', 'i8'), ('ind', 'i8'), 
                                             ('ra',  'f8'), ('dec', 'f8'), 
                                             ('x', 'f8'), ('y', 'f8'),
                                             ('a_world', 'f8'), ('b_world', 'f8'), ('size', 'f8')])

        truth_cat['number'] = np.arange(len(simulated_cat)).astype(np.int64) + 1
        truth_cat['ra']     = simulated_cat['ra']
        truth_cat['dec']    = simulated_cat['dec']
        
        x, y = coadd_wcs.sky2image(simulated_cat['ra'], simulated_cat['dec'])
        truth_cat['x']   = x
        truth_cat['y']   = y
        truth_cat['ind'] = truth_cat['number']

        #We don't write the star catalog anywhere because we don't really use it as a data product
        #in the analysis. Otherwise would write the catalog in here
        
        return truth_cat, simulated_cat

def _render_se_image(
        *, se_info, band, galaxy_truth_cat, star_truth_cat, bounds_buffer_uv,
        draw_method, noise_seed, output_meds_dir, galaxy_src_func, gal_kws, star_src_func):
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
        truth_cat=galaxy_truth_cat,
        se_info=se_info,
        bounds_buffer_uv=bounds_buffer_uv)

    # step 2 - render the objects
    im = _render_all_objects(
        msk_inds=msk_inds,
        truth_cat=galaxy_truth_cat,
        se_info=se_info,
        band=band,
        src_func=galaxy_src_func,
        draw_method=draw_method)
    
    # step 2b - render the star objects
    if star_src_func is not None:
        
        msk_inds = _cut_tuth_cat_to_se_image(
        truth_cat=star_truth_cat,
        se_info=se_info,
        bounds_buffer_uv=bounds_buffer_uv)
        
        star_im = _render_all_objects(
                    msk_inds=msk_inds,
                    truth_cat=star_truth_cat,
                    se_info=se_info,
                    band=band,
                    src_func=star_src_func,
                    draw_method=draw_method)
        
        im += star_im

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

    #If we want Blank image, then we can't add original image
    if not gal_kws.get('BlankImage', False):
        # take the original image and add the simulated + original images together
        original_image = fitsio.read(se_info['nwgint_path'], ext=se_info['image_ext'])
        image += original_image

    # now just read out these other images
    # in practice we just read out --> copy to other location
    # since balrog does not use wgt and bmask
    bkg   = fitsio.read(se_info['bkg_path'], ext=se_info['bkg_ext'])
#     wgt   = fitsio.read(se_info['weight_path'], ext=se_info['weight_ext'])
#     bmask = fitsio.read(se_info['bmask_path'], ext=se_info['bmask_ext'])


    wgt   = fitsio.read(se_info['nwgint_path'], ext=se_info['weight_ext'])
    bmask = fitsio.read(se_info['nwgint_path'], ext=se_info['bmask_ext'])
    
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
    def __init__(self, *, truth_cat, wcs, psf, gal_mag, gal_source, band = None, galsource_rng = None, simulated_catalog = None):
        self.truth_cat = truth_cat
        self.wcs = wcs
        self.psf = psf        
        
        self.gal_source = gal_source
        self.galsource_rng = galsource_rng
        
        self.simulated_catalog = simulated_catalog
        
        self.gal_mag = gal_mag
        self.band    = band
        
            

    def __call__(self, ind):
        pos = self.wcs.toImage(galsim.CelestialCoord(
            ra  = self.truth_cat['ra'][ind]  * galsim.degrees,
            dec = self.truth_cat['dec'][ind] * galsim.degrees))
        
        
        obj = get_desdf_galaxy(desdf_ind = self.truth_cat['ind'][ind],
                                rng  = self.galsource_rng, 
                                data = self.simulated_catalog,
                                band = self.band)

        if self.gal_mag != 'custom':
            normalized_flux = 10**((30 - self.gal_mag)/2.5)
            obj = obj.withFlux(normalized_flux)
            
        #Now do extinction (the coefficients are just zero if we didnt set gal_kws['extinction'] = True)
        A_mag  = self.truth_cat[ind]['A%s' % self.band]
        A_flux = 10**(-A_mag/2.5)
        obj    = obj.withScaledFlux(A_flux)
        
        #Now psf
        psf = self.psf.getPSF(image_pos=pos)
        obj = galsim.Convolve([obj, psf], gsparams = Our_params)
        
        #For doing photon counting, need to do some silly work
        rng = galsim.BaseDeviate(self.galsource_rng.randint(0, 2**32))
        
        return (obj, rng), pos
    
    
class LazyStarSourceCat(object):
    """A lazy source catalog that only builds objects to be rendered as they
    are needed. But now just for stars.

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
    def __init__(self, *, truth_cat, wcs, psf, star_mag, star_source, starsource_rng = None, band = None, simulated_catalog = None):
        
        self.truth_cat = truth_cat
        self.wcs = wcs
        self.psf = psf
        
        self.star_source = star_source
        
        self.simulated_catalog = simulated_catalog
        
        self.star_mag  = star_mag
        self.band      = band
        
        self.starsource_rng = starsource_rng
            

    def __call__(self, ind):
        
        pos = self.wcs.toImage(galsim.CelestialCoord(ra  = self.truth_cat['ra'][ind]  * galsim.degrees,
                                                     dec = self.truth_cat['dec'][ind] * galsim.degrees))
        
        
        if self.star_mag == 'custom':
            mag = self.simulated_catalog['mag_%s'%self.band][ind]
        else:
            mag = self.star_mag

        normalized_flux = 10**((30 - mag)/2.5)

        #No extinction correction since the catalog (LSST sim) already has this included
        
        #Just PSF since stars ARE the point source
        obj = self.psf.getPSF(image_pos = pos).withFlux(normalized_flux)
        
        return obj, pos

    

def mock_balrog_sigmoid(mag_ref, sigmoid_x0, rng):
    """
    
    Function for selecting deep field galaxies at a rate that follows a sigmoid function that smoothly transitions from 1 for bright objects, to a value of 0 for faint objects. 
    Parameters
    ----------
    deep_data : pandas dataframe
        Pandas dataframe containing the deep field data.
    sigmoid_x0 : float
        Magnitude value at which the sigmoid function transitions from 1 to 0.
    N : int
        Number of galaxies to be drawn.
    ref_mag_col : string
        Column name of the reference magnitude in deep_data
    Returns
    -------
    deep_balrog_selection : pandas dataframe
        Pandas dataframe containing a list of N deep field objects to be injected by Balrog.
    """

    weights = 1.0 - 1.0 / (1.0 + np.exp(-4.0 * (mag_ref - sigmoid_x0)))
    weights /= np.sum(weights) #Need to normalize ourselves since numpy choice doesn't do this

    inds = rng.choice(len(mag_ref), len(mag_ref), p = weights, replace = True)
    
    return inds