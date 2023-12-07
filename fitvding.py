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
    
    def __init__(self, *, tilename, bands, output_meds_dir, config):
        
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
                    'CONFIG' : ,
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
    
    def __init__(self, *, tilename, bands, output_meds_dir, config):
        
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
                    'CONFIG' : ,
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
        
        



            
            
'''
# Want to chunk up individual tiles to run on 128 cores on Perlmutter

import multiprocessing as mp
from multiprocessing import Pool
import subprocess
import os
import argparse
from astropy.io import fits
from astropy.table import QTable
from astropy.table import vstack
import numpy as np

n_cpu = 120

parser = argparse.ArgumentParser(description='')
parser.add_argument("-t", '--tile', type=str, help='A tile')
args = parser.parse_args()

base_dir = os.environ['PSCRATCH'] + "/BalrogY6/Y6Running_NERSC/"
tilename = args.tile
fitvd_path = "fitvd"
    
    


tile_ext = tilename + "_" + req_num
shredx_fofslist = "fitvd/" + tile_ext + "_shredx-fofslist.fits"


#print(req_num)
#print(segmap_path)
#print(scat_path)
#print(coaddimage_paths)
#print(psfcat_paths)
#print(meds_paths)

    
    
def get_fofs_chunks(fofs_path):
    
    # Chunk up the fofs groups:
    with fits.open(fofs_path) as hdul:
        fofs = hdul[1].data
    n_fofs = len(np.unique(fofs['fof_id'])) 

    print("fof groups: ", n_fofs)

    n_chunks = n_cpu
    n_objs_chunks = n_fofs // n_chunks

    print("n_chunks : ", n_chunks)
    print("n_objs_chunks : ", n_objs_chunks)

    starts_ends = []
    for i in range(n_chunks):
        start = i*n_objs_chunks
        end = i*n_objs_chunks + n_objs_chunks - 1
        if i == n_chunks-1:
            end = n_fofs - 1
        starts_ends.append((start, end))
    
    return starts_ends
    

    
def get_fitvd_chunks(shredx_fits_path):
    
    with fits.open(shredx_fits_path) as hdul:
        shredx = hdul[1].data
    n_objs = len(shredx) 
    
    n_chunks = n_cpu
    n_objs_chunks = n_objs // n_chunks
    
    starts_ends = []
    for i in range(n_chunks):
        start = i*n_objs_chunks
        end = i*n_objs_chunks + n_objs_chunks - 1
        if i == n_chunks-1:
            end = n_objs - 1
        starts_ends.append((start, end))
    
    return starts_ends       
    
        
    
def run_make_fofs():
    
    # Not in parallel:
    fof_command = "shredx-make-fofs --output " + shredx_fofslist + " --seg " + segmap_path

    print(fof_command)
    print(shredx_fofslist)
    print(segmap_path)
    
    print(os.path.join(base_dir, tilename))
    print(fof_command)
    
    command = "cd " + os.path.join(base_dir, tilename) + "; " + fof_command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    print(stdout, stderr)

    return os.path.join(base_dir, tilename, shredx_fofslist)
    

# Old random seed: 206492073
def run_shredx(start, end):
    
    shredx_chunklist = "fitvd/" + tile_ext + "_shredx-chunk-" + str(start) + "-" + str(end) + ".fits"
    sof_output = "fitvd/" + tile_ext + "_sof-chunk-" + str(start) + "-" + str(end) + ".fits"
    
    shredx_command = ("shredx --cat " + scat_path + " --start " + str(start) + " --end " + str(end) + " --seed 206492073" +
                  " --images " + coaddimage_paths[0] + " " + coaddimage_paths[1] + " " + coaddimage_paths[2] + " " + coaddimage_paths[3] + 
                  " --psf " + psfcat_paths[0] + " " + psfcat_paths[1] + " " + psfcat_paths[2] + " " + psfcat_paths[3] + 
                  " --fofs " + shredx_fofslist + " --config %s/BalrogY6/inputs/Y6A1_v1_shredx-Y6A1v1.yaml " % os.environ['HOME'] + 
                  "--outfile " + shredx_chunklist + " --seg " + segmap_path)
    
    command = "cd " + os.path.join(base_dir, tilename) + "; " + shredx_command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    return stdout, stderr


# new seed: 234324
# Old random seed: 401349271
def run_fitvd(start, end):
    
    shredx_fits = tile_ext + "_shredx.fits"
    sof_output = "fitvd/" + tile_ext + "_sof-chunk-" + str(start) + "-" + str(end) + ".fits"
    
    fitvd_command = ("fitvd --start " + str(start) + " --end " + str(end) + " --seed 401349271" +
                 " --config %s/BalrogY6/inputs/Y6A1_v1_fitvd-Y6A1v5.yaml"  % os.environ['HOME'] + 
                 " --model-pars " + shredx_fits + 
                 " --output " + sof_output + 
                 " " + meds_paths[0] + " " + meds_paths[1] + " " + meds_paths[2] + " " + meds_paths[3])
    
    command = "cd " + os.path.join(base_dir, tilename) + "; " + fitvd_command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    return stdout, stderr



def join_shred_chunks():
    full_shred = QTable()
    for file in os.listdir(os.path.join(base_dir, tilename, fitvd_path)):
        if file[:34] == tile_ext + "_shredx-chunk":

            file_path = os.path.join(base_dir, tilename, fitvd_path, file)
            with fits.open(file_path) as hdul:
                data = hdul[1].data

            data = QTable(data)

            full_shred = vstack([full_shred, data])
            data = None

    file_name = tile_ext + "_shredx.fits"
    full_shred.write(os.path.join(base_dir, tilename, fitvd_path, file_name), overwrite=True)
    
    return os.path.join(base_dir, tilename, fitvd_path, file_name)
    

def join_sof_chunks():
    
    filelist = []
    file_name = tile_ext + "_sof.fits"
    output = os.path.join(base_dir, tilename, file_name)
    
    for file in os.listdir(os.path.join(base_dir, tilename, fitvd_path)):
        if tile_ext + "_sof-chunk" in file:
            file_path = os.path.join(base_dir, tilename, fitvd_path, file)
            
            filelist.append(file_path)
            
            print(file_path)
    
    fitvd_command = ("fitvd-collate --meds " + meds_paths[0] + " --output " + output + " " + " ".join(filelist))
    
    print(fitvd_command)
    command = "cd " + os.path.join(base_dir, tilename) + ";" + fitvd_command
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    
    print(stdout, stderr)
    return stdout, stderr

    
    


if __name__ == '__main__':
    
    # Make the fofs groups:
    shredx_fofslist = run_make_fofs()
    
    # Chunk the fofs groups for shredx:
    starts_ends = get_fofs_chunks(shredx_fofslist)
    
    # Run the groups through the shredder:
    with Pool(n_cpu+2) as p:
        print(p.starmap(run_shredx, starts_ends))
        p.close()
        p.join()
        
    # Combine the shredx files into one fits:
    shredx_fits_path = join_shred_chunks()
    
    # Chunk the objects for fitvd:
    starts_ends = get_fitvd_chunks(shredx_fits_path)
    
    # Run fitvd on the Chunks:
    with Pool(n_cpu+2) as p:
        print(p.starmap(run_fitvd, starts_ends))
        p.close()
        p.join()   
        
    # Join the sof chunks together:
    join_sof_chunks()
        
    print("Done with all")
'''