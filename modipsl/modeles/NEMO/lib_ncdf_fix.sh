#!/bin/bash
# Fixes OPA Makefile to build properly with lib_ncdf.F90
# Must be run from NEMO directory.
# Chris Nickerson, October 2007

TMP=$(pwd | grep NEMO | wc | gawk '{ print $1 }')
if [ $TMP -lt 1 ]
then
    echo "Not in NEMO directory."
    exit
fi

cd WORK
if [[ $? != "0" ]]
then echo "Couldn't find WORK directory."; exit
fi

mv Makefile old.Makefile
if [[ $? != "0" ]]
then echo "Couldn't find Makefile"; exit
fi
gawk '$0 !~ /netcdf.f90/' old.Makefile > temp1.Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi
gawk '$0 !~ /netcdf.o/' temp1.Makefile > temp2.Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi
gawk '$0 !~ /netcdf.mod/' temp2.Makefile > Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi

gawk '$0 !~ /calendar.f90/' temp2.Makefile > temp1.Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi
gawk '$0 !~ /calendar.o/' temp1.Makefile > temp2.Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi
gawk '$0 !~ /calendar.mod/' temp2.Makefile > Makefile
if [[ $? != "0" ]]
then echo "Error fixing Makefile"; exit
fi

/bin/rm -f temp1.Makefile temp2.Makefile&> /dev/null
if [[ $? != "0" ]]
then echo "Error cleaning up temp files"; exit
fi

echo 'Makefile fixed. Backup copy saved as old.Makefile'
