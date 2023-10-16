#!/bin/bash
#SBATCH --job-name balrog_concat
##SBATCH --partition=broadwl
#SBATCH --partition=kicp
#SBATCH --account=kicp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=30:00:00
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


if [ "$USER" == "dhayaa" ]
then
    cd /home/dhayaa/Desktop/DECADE/Balrog_of_the_decade/
    module load python
    conda activate shearDM
    source /home/dhayaa/Desktop/DECADE/Balrog_of_the_decade/bash_profile.sh
fi

cd $SLURM_SUBMIT_DIR

python -u Make_balrog_cat.py