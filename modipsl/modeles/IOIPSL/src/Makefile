#-
#- $Id: AA_make.ldef,v 2.6 2006/01/18 06:17:31 adm Exp $
#-
#---------------------------------------------------------------------
#- Local parameters for "libioipsl"
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../../lib
BINDIR = ../../../bin
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
D_P = I4R8
P_P =
MODEL_LIB = $(LIBDIR)/libioipsl.a
SXMODEL_LIB = $(MODEL_LIB)
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
#-
#- $Id: AA_make,v 2.16 2006/01/18 06:34:11 adm Exp $
#-
#---------------------------------------------------------------------
#- Create "libioipsl"
#---------------------------------------------------------------------

PRF_PRC = !\#-P-
DEF_PRC = $(PRF_PRC)$(D_P)
SRC_PRC   = defprec.f90
MODS1     = $(SRC_PRC)    \
            errioipsl.f90 \
            stringop.f90  \
            mathelp.f90   \
            getincom.f90  \
            calendar.f90  \
            fliocom.f90   \
            flincom.f90   \
            histcom.f90   \
            restcom.f90   \
            ioipsl.f90
OBJSMODS1 = $(MODS1:.f90=.o)

all: $(SRC_PRC) $(MODEL_LIB)($(OBJSMODS1))
	@echo IOIPSL is OK

$(SRC_PRC): def.prec
	(sed -e "s/^$(DEF_PRC) */  /g" def.prec | \
	 grep -v $(PRF_PRC) > $(SRC_PRC))

.PRECIOUS : $(MODEL_LIB)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) -I$(NCDF_INC) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	$(RM) $*.o

clean:
	$(RM) *.*~ Makefile~ core *.o *.mod i.*.L *.L i.*.f90
	$(RM) $(MODDIR)/*.mod $(MODDIR)/*.M *.M
	$(RM) $(MODEL_LIB) $(SRC_PRC)

#- Specific dependencies

$(MODEL_LIB)(errioipsl.o):  \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(stringop.o):  \
  $(MODEL_LIB)(defprec.o)

$(MODEL_LIB)(mathelp.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(getincom.o):  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(calendar.o):    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)

$(MODEL_LIB)(fliocom.o):     \
  $(MODEL_LIB)(defprec.o)    \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(flincom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(calendar.o)

$(MODEL_LIB)(histcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)  \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(restcom.o):     \
  $(MODEL_LIB)(errioipsl.o)  \
  $(MODEL_LIB)(stringop.o)   \
  $(MODEL_LIB)(mathelp.o)    \
  $(MODEL_LIB)(calendar.o)  \
  $(MODEL_LIB)(fliocom.o)

$(MODEL_LIB)(ioipsl.o):      \
  $(MODEL_LIB)(getincom.o)   \
  $(MODEL_LIB)(fliocom.o)    \
  $(MODEL_LIB)(flincom.o)    \
  $(MODEL_LIB)(histcom.o)    \
  $(MODEL_LIB)(restcom.o)

