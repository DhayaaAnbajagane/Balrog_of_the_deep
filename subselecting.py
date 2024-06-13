import numpy as np
import fitsio, yaml
from sklearn.neighbors import BallTree
import os

from files import get_mcal_file_path, get_truth_catalog_path, get_band_info_file, get_balrog_file_path
from constants import MEDSCONF

TMP_DIR = os.environ['TMPDIR']

def keep_injected_only(radius, output_desdata, tilename, bands):
    
    Truth_path = get_truth_catalog_path(meds_dir=output_desdata, medsconf=MEDSCONF, tilename=tilename)
    Binfo_path = {b : get_band_info_file(meds_dir=output_desdata, medsconf=MEDSCONF, tilename=tilename, band = b) for b in bands}
    info       = {band : yaml.load(open(Binfo_path[band], 'r'), Loader=yaml.Loader) for band in bands}
    BCat_path  = [info[band]['cat_path'].replace(TMP_DIR, output_desdata) for band in bands] #Path to new SrcExtractor
    OCat_path  = [B.replace('.fits', '_unmodded.fits') for B in BCat_path] #Path to modified SrcExtractor
    
    
    #Do the matching with just the first catalog. This works since the ra/dec is same for all bands
    fid   = fitsio.read(BCat_path[0])
    Truth = fitsio.read(Truth_path)
        
    tree = BallTree(
                    np.stack(
                             [fid['DELTAWIN_J2000'], fid['ALPHAWIN_J2000']], 
                             axis = 1) * np.pi/180, 
                    leaf_size=40, metric="haversine"
            )

    inds = tree.query_radius( 
                                np.stack([Truth['dec'], Truth['ra']], axis = 1) * np.pi/180,
                                radius * np.pi/180 / 60 #Extra 1/60 since radius is in arcmin
                        )
    
    inds = np.unique( np.concatenate(inds) )
    
    print(f"ONLY KEEPING {len(fid)} ---> {len(inds)} OBJECTS")

    
    for i in range(len(BCat_path)):
        
        fid = fitsio.read(BCat_path[i])
        new = fid[inds]
        
        fitsio.write(fid, OCat_path[i])
        fitsio.write(new, BCat_path[i])
        
        print("Finished modifying catalog", BCat_path[i])
        
        
    for band in bands:
        
        path   = info[band]['seg_path'].replace(TMP_DIR, output_desdata)
        segmap = fitsio.read(path)
        newmap = np.zeros_like(segmap)
        
        for i in range(inds.size):
            newmap[segmap == inds[i]] = i
            
        mask = (segmap > 0) & np.invert( np.isin(segmap, inds) )
        assert np.sum(newmap[mask]) == 0, "Some mixup in the segmentation map"
        
        newmap[mask] = int(999_999) #sentinel value so uberseg still knows about true objects
        
        fitsio.write(newmap, segmap)
        fitsio.write(segmap, path.replace('.fits', '_unmodded.fits'))
        

        
