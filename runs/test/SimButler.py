import os, io, subprocess as sp
import jinja2
import yaml
import pandas as pd
import argparse
import datetime as dt
import numpy as np
import glob
import time

my_parser = argparse.ArgumentParser()

my_parser.add_argument('--Initialize', action='store_true', default = False)
my_parser.add_argument('--Maintain',   action='store_true', default = False)
my_parser.add_argument('--TileStart',  action='store', type = int, default = 0)
my_parser.add_argument('--TileCount',  action='store', type = int, default = 10)

my_parser.add_argument('--MaxConcurrentJobs', action='store', type = int, default = 10)
my_parser.add_argument('--MaxCutoffTime',     action='store', type = int, default = 3600*48) #Maxcutoff time in seconds

args = vars(my_parser.parse_args())

#Print args for debugging state
print('-------INPUT PARAMS----------')
for p in args.keys():
    print('%s : %s'%(p.upper(), args[p]))
print('-----------------------------')
print('-----------------------------')

TILENAME_SEED = 4242
NUM_OF_TILES  = 16800
tiles = pd.read_csv(os.environ['BALROG_RUN_DIR'] + '/data/Tilelist_DR3_1_All.csv')
tilenames = list(np.random.default_rng(TILENAME_SEED).choice(tiles['TILENAME'].values, size = NUM_OF_TILES, replace = False))[0:args['TileCount']]

print(tilenames)
if __name__ == '__main__':
    
    #Automatically get folder name (assuming certain folder structure)
    name = os.path.basename(os.path.dirname(__file__))

    #Create output directory for metacal
    BALROG_DIR = os.environ['BALROG_DIR']
    os.makedirs(BALROG_DIR +'/'+ name, exist_ok=True)

        
    ########################################################################3
    
    def is_finished(tilename):
        
        args  = {'name' : name, 'tile' : tilename}
        plus  = os.path.join(os.environ['BALROG_DIR'], "/%(name)s/metacal_%(tile)s.fits" % args)
        
        return os.path.isfile(plus)
    
    def is_job(tilename):
        
        return os.path.isfile('job_%s.sh' % tilename)
            
    def create_job(tilename, gal_seed):
        
        #Now create all job.sh files for running sims
        with open('job.sh.temp', 'r') as fp:
            tmp = jinja2.Template(fp.read())
        
        if not is_finished(tilename):
            with open('job_%s.sh' % tilename, 'w') as fp:
                fp.write(tmp.render(tilename=tilename, model_name = name, seed_galsim=gal_seed, seed_mcal=42))
            os.system('chmod u+x job_%s.sh' % tilename)
        
    def prep_job(tilename):
        
        os.system('python $BALROG_RUN_DIR/run_sims.py prep --tilename="%s" --bands="riz" --output-desdata="$PREP_DIR/%s/outputs_%s" --config-file="config.yaml"'%(tilename, name, tilename))
        
        
    def current_job_count():
        
        x = sp.check_output("squeue --format='%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R' --sort=+i -u dhayaa", shell = True)
        j = pd.read_csv(io.StringIO(x.decode("utf-8")), delim_whitespace=True)
        
        count = 0
        tiles = []
        
        #metacal_DES0849+0252_seed0_gminus
        for i in range(len(j)):
            
            condition1 = 'balrog' in j['NAME'].values[i]
            condition2 = j['STATE'].values[i] == 'RUNNING'
            
            if condition1:
                count += 1
                tiles.append(j['NAME'].values[i][-12:])
                
        for t in tiles:
            print("CURRENTLY RUNNING: %s" % t)
                
        return count
            
    def tilename_from_jobfilename(job):
        
        #job_DES0849+0252_minus.sh
            
        return job[4:4+12]
    
                
    ################################################################
    
    #In initial step we just create every job that we need
    if args['Initialize']:
        
        for i, tilename in enumerate(tilenames):
        
            if i < args['TileStart']: 
                continue
            else:
                create_job(tilename, i)
    
    #In next step we keep track of all jobs and add to queue when we can
    elif args['Maintain']:
        
        job_list = sorted(glob.glob('job_*.sh'))
        start = dt.datetime.now()

        while (dt.datetime.now() - start).seconds < 3600*48: #48 hours max
            
            if len(job_list) == 0:
                print("ALL JOBS HAVE BEEN STARTED. NOTHING TO MAINTAIN")
                break
                
            if current_job_count() >= args['MaxConcurrentJobs']:
                
                print("---------------------------")
                print(dt.datetime.now())
                print("---------------------------") 
                time.sleep(60*5)

            else:
                
                j = job_list[0]
                t = tilename_from_jobfilename(j)
                
                prep_job(t)
                
                os.system('sbatch %s' % j)
                os.system('rm %s' % j)
                
                print("SUBMITTED JOB %s" % j)
                
                #This just removes the top item from the list
                #which we do as we just ran that job
                job_list = job_list[1:]
                
        else:
            
            print("GONE OVERTIME. SHUTTING DOWN SCRIPT")            
                
    
