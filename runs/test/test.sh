#!/bin/bash
#SBATCH --job-name balrog_test
#SBATCH --partition=chihway
#SBATCH --account=pi-chihway
#SBATCH --nodes=1
##SBATCH --nodelist=midway2-0696
#SBATCH --ntasks-per-node=40
#SBATCH --time=30:00:00
#SBATCH --output=/home/dhayaa/Desktop/DECADE/Balrog_of_the_deep/runs/test/%x.log
#SBATCH --mail-user=dhayaa@uchicago.edu
#SBATCH --mail-type=BEGIN,END


cd /home/dhayaa/Desktop/DECADE/Balrog_of_the_deep/
source /home/dhayaa/.bashrc
conda activate shearDM
source /home/dhayaa/Desktop/DECADE/Balrog_of_the_deep/bash_profile.sh

output="$PREP_DIR/test/outputs_test"
bands="griz"
tilename="DES1126+2543"

NOW=$( date "+%H:%M:%S" )
echo "Starting galsim stage at $NOW"


if [ -z "$SLURM_JOB_ID" ]
then
    echo "This script is running on a login node."
    python run_sims.py prep --tilename=$tilename --bands=$bands --output-desdata=$output --config-file="runs/test/config.yaml"
else
    python run_sims.py \
      galsim \
      --tilename="$tilename" \
      --bands="$bands" \
      --output-desdata="$output" \
      --config-file="runs/test/config.yaml" \
      --seed="42"

    NOW=$( date "+%H:%M:%S" )
    echo "Completed galsim and starting swarp at $NOW"


    python run_sims.py \
      swarp \
      --tilename="$tilename" \
      --bands="$bands" \
      --output-desdata="$output" \
      --config-file="runs/test/config.yaml"

    NOW=$( date "+%H:%M:%S" )
    echo "Completed swarp and starting source-extractor at $NOW"


    python run_sims.py \
      source-extractor \
      --tilename="$tilename" \
      --bands="$bands" \
      --output-desdata="$output" \
      --config-file="runs/test/config.yaml"

    NOW=$( date "+%H:%M:%S" )
    echo "Completed source-extractor and starting meds at $NOW"


    python run_sims.py \
      meds \
      --tilename="$tilename" \
      --bands="$bands" \
      --output-desdata="$output" \
      --config-file="runs/test/config.yaml" \
      --meds-config-file="runs/test/meds.yaml"
    NOW=$( date "+%H:%M:%S" )
    echo "Completed meds and starting mcal at $NOW"


    conda deactivate; conda activate des-y6-fitvd

    python run_sims.py \
      fitvd \
      --tilename="$tilename" \
      --bands="$bands" \
      --output-desdata="$output" \
      --seed="42" \
      --config-file="runs/test/config.yaml" \
      --shredx-config-file="/home/dhayaa/Desktop/DECADE/Balrog_of_the_deep/configs/Y6A1_v1_shredx-Y6A1v1.yaml" \
      --fitvd-config-file="/home/dhayaa/Desktop/DECADE/Balrog_of_the_deep/configs/Y6A1_v1_fitvd-Y6A1v5.yaml" \


    NOW=$( date "+%H:%M:%S" )
    echo "Completed meds and starting fitvd at $NOW"

#     conda deactivate; conda activate shearDM

#     python run_sims.py \
#       match \
#       --tilename="$tilename" \
#       --bands="$bands" \
#       --output-desdata="$output" \
#       --config-file="runs/test/config.yaml"

#     NOW=$( date "+%H:%M:%S" )
#     echo "Completed matching and starting finalize at $NOW"

#     python run_sims.py \
#       finalize \
#       --tilename="$tilename" \
#       --bands="$bands" \
#       --output-desdata="$output" \
#       --config-file="runs/testconfig.yaml"

#     NOW=$( date "+%H:%M:%S" )
#     echo "Finished moving meds, mcal, SrcExtractor output to $BALROG_DIR/test at $NOW"
fi
