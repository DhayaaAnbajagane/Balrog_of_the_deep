This directory contains the
versions of the terapix source code for

sextractor version 2.24.4

swarp version 2.40.1

pixcorrect version 0.5.5 (which contains the coadd_nwgint python call)

Along with the Y6 DESDM config files used for sextractor, swarp and coadd_nwgint

along with same command lines calls of each which have
parameters on the command line that override those
in the config file.

There are also product dependency lists for
these source codes,
for instance sextractor depends on openblas for
linear algebra libraries

pixcorrect depends on fitsio and a bunch of
other products -- one may wish to
reproduce the functionality of the coadd_nwgint
call in a small python code rather than
remake the whole pixcorrect stack  or one may
consider skip this step without too many consequences. 

The current stack for pixcorrect and dependencies
is python2 based, but there is a python3 version
also availale from DESDM EUPS.  

