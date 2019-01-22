#!/bin/bash

# Does a clean build of OPA including the new lib_ncdf I/O library
# Must be run from modipsl/modeles/NEMO
curdir=`pwd`
/bin/rm -rf ../../lib/oce/*
/bin/rm -rf ../../lib/*
#
#../UTIL/fait_config SS008_SPINUP
../UTIL/fait_config SS008
../UTIL/fait_AA_make 
../../util/clr_make
#
../../util/ins_make -p I4R8 -t lxiv8
#
./lib_ncdf_fix.sh
#
cd ../../../modipsl/config/SS008/
make
cd $curdir
