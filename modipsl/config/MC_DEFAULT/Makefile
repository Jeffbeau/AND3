#---------------------------------------------------------------------
#- Creation des elements relatifs a ORCA :
#-      ORCA  (lib, executable : opa)
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../../lib
BINDIR = ../../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
#-
#-
#- $Id: AA_make.gdef,v 2.24 2006/03/24 15:31:56 adm Exp $
#-
#
# Definitions fpr HP Proliant cluster running Red Hat Linux, HP MPI 
# and Intel compilers
# AD: Mar25-2010 -traceback -check
# I though we need -I/usr/include/ -I/opt/hpmpi/include to get mpif.h
#
#
#- Global definitions for BIO COMDA HPC
M_K = make
P_C =  cpp  
P_O = -P -C -traditional $(P_P) 
F_D =
F_P = -fdefault-real-8  
NCDF_INC =  -I/usr/include/
NCDF_LIB = -L/usr/lib64/openmpi/lib/ -lnetcdff -lnetcdf 
w_w = $(F_D) $(F_P) -I$(MODDIR) -I$(NCDF_INC) $(NCDF_LIB)
F_C = mpif90 -c $(w_w) -O2 
F_O = $(w_w) $(P_P) -J$(MODDIR) -fno-second-underscore -ffree-line-length-200
F_L = mpif90  -J$(MODDIR) -O2 
L_O = 
A_C = ar -r 
A_G = ar -x
C_C = mpicc -c
C_O =
C_L = mpicc
#-
#-
RM      = rm -f
STRIP   = strip
SIZE    = size
# $Header: /home/opalod/NEMOCVSROOT/CONFIG/ORCA2_LIM/AA_make,v 1.4 2005/03/23 10:07:37 opalod Exp $
all: ioipsl opa

ioipsl:
	(cd ../../modeles/IOIPSL/src; gmake;)


opa:
	(cd ../../modeles/NEMO/WORK;  gmake ;)


clean:
	(cd ../../modeles/NEMO/WORK;  gmake clean ;)
