'''
Routines for moving files around, deleting remaining files
and saving relevant quantities.
'''

import shutil
import yaml
import os
from files import get_truth_catalog_path, get_balrog_file_path, get_band_info_file, get_mcal_file_path
from constants import MEDSCONF


def finalize_files(tilename, bands, output_desdata, config):
   
    try: 
        for b in bands:
            move_SrcExtractor_cat(tilename, b, output_desdata)
            move_OldSrcExtractor_cat(tilename, b, output_desdata)
            if config['files']['save_meds'] == True: 
                move_meds(tilename, b, output_desdata)
        
        move_metacal_cat(tilename, output_desdata)
        move_balrog_cat(tilename, output_desdata)
        move_Truth_cat(tilename, output_desdata)
    
    except:
        
        print("SOMETHING CRASHED. SIGH")

    if config['files']['clean_tmpdir'] == True:
        cleanup_tmpdir_files(tilename, output_desdata)
        

#Helper functions to run the above cleanup/re-organization code
def move_SrcExtractor_cat(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile,
            'band' : band}
    
    #This is always going to be path in our runs so just sort of hardcode this assumption
#     config_path = output_desdata + '/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
    with open(get_band_info_file(meds_dir=output_desdata, medsconf=MEDSCONF, tilename=tile, band = band), 'r') as fp:
        band_info = yaml.load(fp, Loader=yaml.Loader)

    cat_path = band_info['cat_path'].replace(os.environ['TMPDIR'], output_desdata)

    new_path = os.environ['BALROG_DIR'] + "/%(name)s/SrcExtractor_%(tile)s_%(band)s-cat.fits" % args

    #print(cat_path, new_path)
    shutil.move(cat_path, new_path)
        
    return True


def move_OldSrcExtractor_cat(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile,
            'band' : band}
    
    #This is always going to be path in our runs so just sort of hardcode this assumption
#     config_path = output_desdata + '/simple_des_y3_sims/y3v02/band_info_files/%(tile)s_%(band)s_info.yaml'%args
    with open(get_band_info_file(meds_dir=output_desdata, medsconf=MEDSCONF, tilename=tile, band = band), 'r') as fp:
        band_info = yaml.load(fp, Loader=yaml.Loader)

    cat_path = band_info['cat_path']

    new_path = os.environ['BALROG_DIR'] + "/%(name)s/OldSrcExtractor_%(tile)s_%(band)s-cat.fits" % args

    shutil.move(cat_path, new_path)
        
    return True


#Helper functions to run the above cleanup/re-organization code
def move_Truth_cat(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile}
    
    cat_path = get_truth_catalog_path(meds_dir = output_desdata, medsconf = MEDSCONF, tilename = tile)
    new_path = os.environ['BALROG_DIR'] + "/%(name)s/Input_%(tile)s-cat.fits" % args

    print(cat_path, new_path)
    shutil.move(cat_path, new_path)
        
    return True


def move_metacal_cat(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile}
    
    cat_path = get_mcal_file_path(meds_dir=output_desdata, medsconf =MEDSCONF, tilename = tile)
    new_path = os.environ['BALROG_DIR'] + "/%(name)s/metacal_%(tile)s.fits" % args

    shutil.move(cat_path, new_path)
         
    return True

def move_balrog_cat(tile, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile}
    
    cat_path = get_balrog_file_path(meds_dir=output_desdata, medsconf =MEDSCONF, tilename = tile)
    new_path = os.environ['BALROG_DIR'] + "/%(name)s/balrog_%(tile)s.fits" % args

    shutil.move(cat_path, new_path)
         
    return True

def move_meds(tile, band, output_desdata):
    
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile,
            'band' : band}
    
    meds_path = get_meds_file_path(meds_dir=output_desdata, medsconf =MEDSCONF, tilename = tile, band = band)
    new_path  = os.environ['BALROG_DIR'] + "/%(name)s/meds_%(tile)s_%(band)s-y3v02.fits.fz" % args

    #print(meds_path, new_path)
    shutil.move(meds_path, new_path)
        
    return True


def cleanup_tmpdir_files(tile, output_desdata):
    
    #Checks if both plus and minus measurements have been done, and deletes
    #input files accordingly
    args = {'dir'  : output_desdata,
            'name' : os.path.basename(os.path.dirname(output_desdata)),
            'tile' : tile}
    
    file  = os.environ['BALROG_DIR'] + "/%(name)s/metacal_%(tile)s.fits" % args

    if True: #os.path.isfile(file):
        file_paths = os.environ['PREP_DIR'] + "/%(name)s/*%(tile)s*" % args
        os.system("rm -rv %s" % file_paths)

        print(file_paths)
        
        file_paths = os.environ['TMPDIR'] + "/*%(tile)s*" % args
        os.system("rm -rv %s" % file_paths)

        print(file_paths)
        print(args)

        print(tile, output_desdata)
        
    return True
