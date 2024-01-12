
#source your conda environment first
#Then source this file as "source bash_profile.sh"

export PATH=$PATH:/home/dhayaa/.local/bin/

#Hacking to get cfitsio to work :/
PATH=$PATH:/home/dhayaa/Desktop/DECADE/cfitsio-4.0.0/

#PATH for files
export EBV_PATH=/project/chihway/dhayaa/DECADE/Imsim_Inputs/ebv_sfd98_fullres_nside_4096_ring_equatorial.fits
export CAT_PATH=/project/chihway/dhayaa/DECADE/Jonah_files/TMP.fits

export PREP_DIR=/scratch/midway2/dhayaa/PREP_OUTPUTS
export TMPDIR=/scratch/midway2/dhayaa/TMP_DIR
export MEDS_DIR=/scratch/midway2/dhayaa/MEDS_DIR
export BALROG_DIR=/project/chihway/dhayaa/DECADE/Balrog/
export DESDM_CONFIG=/home/dhayaa/Desktop/DECADE/Y6DESDM
export SWARP_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/swarp-2.40.1/
export SRCEXT_DIR=/home/dhayaa/Desktop/DECADE/Y6DESDM/sextractor-2.24.4

#DES-easyaccess
export DESREMOTE_RSYNC=desar2.cosmology.illinois.edu::ALLDESFiles/desarchive
export DESREMOTE_RSYNC_USER=dhayaa
export DES_RSYNC_PASSFILE=${HOME}/.desrsyncpass

export DECADEREMOTE_WGET=https://decade.ncsa.illinois.edu/deca_archive/

export EMAIL_ID=dhayaa@uchicago.edu
