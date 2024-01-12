import numpy as np
import h5py
import fitsio
import glob
import os
from tqdm import tqdm

if __name__ == "__main__":
    

    name     = os.path.basename(os.path.dirname(__file__))
    BROG_DIR = os.environ['BALROG_DIR']
    PATH     = BROG_DIR + '/' + name
    print('GETTING BALROG FILES FROM:', PATH)
    
    files = glob.glob(PATH + '/balrog*')
    
    FINAL_CAT = []
    tilenames = []
    for f in tqdm(files, desc = 'Concatenating fits files'):
        
        tile = os.path.basename(f).split('_')[1].split('.')[0]
        cat  = fitsio.read(f)
        tilenames.append([tile] * len(cat))
        FINAL_CAT.append(cat)
        
    FINAL_CAT = np.concatenate(FINAL_CAT, axis = 0)
    tilenames = np.concatenate(tilenames, axis = 0)
    
    with h5py.File(PATH + '/BalrogOfTheDECADE_Catalog.hdf5', 'w') as f:
    
        for i in tqdm(FINAL_CAT.dtype.names, desc = 'Making HDF5'):

            f.create_dataset(i, data = FINAL_CAT[i])
            
#         f.create_dataset('tilename', data = tilenames)
        
    