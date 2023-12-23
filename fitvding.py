import numpy as np
import yaml
import os, glob, sys
from constants import MEDSCONF, PIFF_RUN
from files import get_meds_file_path, get_fitvd_path
import fitsio
import joblib
from constants import MEDSCONF, MAGZP_REF

TMP_DIR = os.environ['TMPDIR']


class MakeShredxCats(object):
    """
    Class to create shredx files from MEDS
    """
    
    def __init__(self, *, tilename, bands, output_meds_dir, config, shredx_config_path):
        
        self.output_meds_dir = output_meds_dir
        self.tilename = tilename
        self.bands = bands
        
        self.info = {}
        for band in bands:
            fname = get_band_info_file(
                meds_dir=self.output_meds_dir,
                medsconf=MEDSCONF,
                tilename=self.tilename,
                band=band)
            with open(fname, 'r') as fp:
                self.info[band] = yaml.load(fp, Loader=yaml.Loader)
                
                
        self.base_path   = get_fitvd_path(meds_dir = self.output_meds_dir, medsconf = MEDSCONF)
        self.fof_path    = self.base_path + '/shredx_fofs.fits'
        self.shredx_path = self.base_path + '/shredx.fits'
        self.fitvd_path  = self.base_path + '/fitvd.fits'
        self.config_path = shredx_config_path
        
    
    def run(self):
        
        self.get_file_list()
        self.make_fofs()
        self.run_shredx()
        self.collate_shredx()
        self.cleanup()
        
        print("FINISHED SHREDX")
    
    
    def get_file_list(self):
        
        segmap = {}
        coadd  = {}
        psfcat = {}
        meds   = {}
        srcext = {}
        
        for band in self.bands:
            
            segmap[band] = self.info[band]['seg_path'].replace(TMP_DIR, self.output_meds_dir)  
            coadd[band]  = self.info[band]['image_path'].replace(TMP_DIR, self.output_meds_dir)   
            psfcat[band] = self.info[band]['psf_path'].replace(TMP_DIR, self.output_meds_dir)  
            meds[band]   = get_meds_file_path(meds_dir = self.output_meds_dir, medsconf = MEDSCONF, tilename = self.timename, band = band)
            srcext[band] = self.info[band]['cat_path'].replace(TMP_DIR, self.output_meds_dir)  
        
        self.files = {'segmap' : segmap,
                      'coadd'  : coadd,
                      'psfcat' : psfcat,
                      'meds'   : meds,
                      'srcext' : srcext}
    
    def make_fofs(self):
       
        # Not in parallel:
        fof_command = "shredx-make-fofs --output " + self.fof_path + " --seg " + self.files['segmap']['r']

        os.system(fof_command)
        
    
    def run_shredx(self):
        
        # Chunk up the fofs groups:
        fofs   = fits.open(self.fof_path)[1].data
        n_fofs = len(np.unique(fofs['fof_id']))

        print("fof groups: ", n_fofs)

        n_chunks = os.cpu_count()
        n_objs_chunks = n_fofs // n_chunks

        print("n_chunks : ", n_chunks)
        print("n_objs_chunks : ", n_objs_chunks)

        starts_ends = []
        for i in range(n_chunks):
            start = i*n_objs_chunks
            end   = i*n_objs_chunks + n_objs_chunks - 1
            if i == n_chunks-1:
                end = n_fofs - 1
            starts_ends.append((start, end))
            
            
        print("MADE SHREDX CHUNKS")

        
        def run_shredx_per_chunk(start, end):
            
            os.makedirs(self.base + '/shredx_chunks/')
            shredx_chunk_path = self.base + '/shredx_chunks/' + self.tilename "_shredx-chunk_%09d-%09d.fits" % (start,end)

            args = {'SRCEXT_PATH' : self.files['srcext']['r'], #Just need one, so always use r-band
                    'START' : start,
                    'END' : end,
                    'COADD_IMAGES' : " ".join([self.files['coadd'][b] for b in bands]),
                    'COADD_PSFS' : " ".join([self.files['psf'][b] for b in bands]),
                    'CONFIG' : self.config_path,
                    'OUTFILE' : shredx_chunk_path,
                    'SEGMAP_PATH' : self.files['segmap']['i'] #Y6 uses i-band here so I do the same
                   }
                    
            SHREDX_COMMAND = "shred --cat %(SRCEXT_PATH)s \
                                    --start %(START)d \
                                    --end %(END)d \
                                    --seed 42 \
                                    --images %(COADD_IMAGES)s \
                                    --psf %(COADD_PSFS)s \
                                    --fofs %(FOFLIST)s \
                                    --config $(CONFIG)s \
                                    --outfile %(OUTFILE_PATH)s \
                                    --seg %(SEGMAP_PATH)s" % args
            
            os.system(SHREDX_COMMAND)
            
            
        with joblib.parallel_backend("loky"):
            jobs    = [joblib.delayed(run_shredx_per_chunk)(*starts_ends[i]) for i in range(len(start_ends))]
            outputs = joblib.Parallel(n_jobs = -1, verbose=10)(jobs)
            
        print("FINISHED SHREDX CHUNKS")
    
    
    
    def collate_shredx(self):
        
        file_list = sorted(glob.glob(self.base + '/shredx_chunks/' + self.tilename "_shredx-chunk_*"))
        
        args = {'TILENAME' : self.tilename,
                'OUTPUT' : self.shredx_path,
                'FLIST' : file_list}
        
        SHREDX_COLLATE_COMMAND = "shredx-collate --tilename %(TILENAME)s \
                                                 --output %(OUTPUT)s \
                                                 --flist %(FILELIST)s" % args
        
        os.system(SHREDX_COLLATE_COMMAND)
        
        print("FINISHED COLLATING SHREDX CHUNKS")
        
        
    
    def cleanup(self):
        
        os.system('rm -r %s' % (self.base + '/shredx_chunks/'))
        os.system('rm %s' self.fof_path)
        
        print("FINISHED CLEANING FOF FILES and SHREDX CHUNKS")
        
        

class MakeFitvdCats(object):
    """
    Class to create fitvd files from MEDS
    """
    
    def __init__(self, *, tilename, bands, output_meds_dir, config, fitvd_config_path):
        
        self.output_meds_dir = output_meds_dir
        self.tilename = tilename
        self.bands = bands
        
        self.info = {}
        for band in bands:
            fname = get_band_info_file(
                meds_dir=self.output_meds_dir,
                medsconf=MEDSCONF,
                tilename=self.tilename,
                band=band)
            with open(fname, 'r') as fp:
                self.info[band] = yaml.load(fp, Loader=yaml.Loader)
                
                
        self.base_path   = get_fitvd_path(meds_dir = self.output_meds_dir, medsconf = MEDSCONF)
        self.fof_path    = self.base_path + '/shredx_fofs.fits'
        self.shredx_path = self.base_path + '/shredx.fits'
        self.fitvd_path  = self.base_path + '/fitvd.fits'
        self.config_path = fitvd_config_path
    
    
    def run(self):
        
        self.get_file_list()
        self.run_fitvd()
        self.collate_fitvd()
        self.cleanup()
        
        print("FINISHED FITVD")
    
    
    def get_file_list(self):
        
        segmap = {}
        coadd  = {}
        psfcat = {}
        meds   = {}
        srcext = {}
        
        for band in self.bands:
            
            segmap[band] = self.info[band]['seg_path'].replace(TMP_DIR, self.output_meds_dir)  
            coadd[band]  = self.info[band]['image_path'].replace(TMP_DIR, self.output_meds_dir)   
            psfcat[band] = self.info[band]['psf_path'].replace(TMP_DIR, self.output_meds_dir)  
            meds[band]   = get_meds_file_path(meds_dir = self.output_meds_dir, medsconf = MEDSCONF, tilename = self.timename, band = band)
            srcext[band] = self.info[band]['cat_path'].replace(TMP_DIR, self.output_meds_dir)  
        
        self.files = {'segmap' : segmap,
                      'coadd'  : coadd,
                      'psfcat' : psfcat,
                      'meds'   : meds,
                      'srcext' : srcext}
    
    
    def run_fitvd(self):

        # Chunk up the fofs groups:
        shredx = fits.open(self.shredx_path)[1].data
        n_objs = len(shredx)

        print("n_objs groups: ", n_objs)

        n_chunks = os.cpu_count()
        n_objs_chunks = n_objs // n_chunks

        print("n_chunks : ", n_chunks)
        print("n_objs_chunks : ", n_objs_chunks)

        starts_ends = []
        for i in range(n_chunks):
            start = i*n_objs_chunks
            end   = i*n_objs_chunks + n_objs_chunks - 1
            if i == n_chunks-1:
                end = n_fofs - 1
            starts_ends.append((start, end))
            
        print("MADE FITVD CHUNKS")

        
        def run_fitvd_per_chunk(start, end):
            
            os.makedirs(self.base + '/fitvd_chunks/')
            fitvd_chunk_path = self.base + '/fitvd_chunks/' + self.tilename "_fitvd-chunk_%09d-%09d.fits" % (start,end)

            
            args = {'START' : start,
                    'END' : end,
                    'CONFIG' : self.config_path,
                    'OUTPUT' : fitvd_chunk_path,
                    'SHREDX_PATH' : self.shredx_path,
                    'MEDS_PATHS' : " ".join([self.files['meds'][b] for b in bands])
                   }
                    
            
            FITVD_COMMAND = "fitvd --start %(START)d \
                                   --end %(END)d \
                                   --seed 42 \
                                   --config $(CONFIG)s \
                                   --model-pars %(SHREDX_PATH)s \
                                   --output %(OUTPUT)s \
                                   %(MEDS_PATHS)s" % args
            
                
            os.system(FITVD_COMMAND)
            
            
        with joblib.parallel_backend("loky"):
            jobs    = [joblib.delayed(run_fitvd_per_chunk)(*starts_ends[i]) for i in range(len(start_ends))]
            outputs = joblib.Parallel(n_jobs = -1, verbose=10)(jobs)
            
        print("FINISHED FITVD CHUNKS")
    
    
    
    def collate_fitvd(self):
        
        file_list = sorted(glob.glob(self.base + '/fitvd_chunks/' + self.tilename "_fitvd-chunk_*"))
        
        args = {'TILENAME' : self.tilename,
                'MEDS_PATH' : self.files['meds']['r'],
                'FLIST' : file_list}
        
        "fitvd-collate --meds " + meds_paths[0] + " --output " + output + " " + " ".join(filelist)
        
        FITVD_COLLATE_COMMAND = "fitvd-collate --meds %(MEDS_PATH)s \
                                               --output %(OUTPUT)s \
                                               %(FILELIST)s" % args
        
        os.system(FITVD_COLLATE_COMMAND)
        
        print("FINISHED COLLATING FITVD CHUNKS")
        
        
    
    def cleanup(self):
        
        os.system('rm -r %s' % (self.base + '/fitvd_chunks/'))
        
        print("FINISHED CLEANING FOF FILES and FITVD CHUNKS")
        
        
        