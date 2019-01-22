#- $Id: AA_make.ldef,v 1.1.1.1 2006/03/10 17:58:35 opalod Exp $
#---------------------------------------------------------------------
#- Creation des elements relatifs a AGRIF (lib, config)
#---------------------------------------------------------------------
SHELL = /bin/sh
#---------------------------------------------------------------------
LIBDIR = ../../lib
MODDIR = $(LIBDIR)
#---------------------------------------------------------------------
#--- -DNC_REAL option is needed if compiled in simple precision
#--- In case you switch to double precision use -DNC_DOUBLE
#-
# Compilation:
FFLAGS_EXTEND=  -I..
#-
MODEL_LIB = $(LIBDIR)/libagrif.a
SXMODEL_LIB = $(MODEL_LIB)
#-
P_P=$(AGRIF_MPI)
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
#- $Id: AA_make,v 1.2 2006/04/26 10:37:06 opalod Exp $

TMP=../../../tmp

SRM=/bin/rm -f


AGRIFOBJS=AGRIF_OBJS
AGRIFFILES=AGRIF_FILES
AGRIFYOURFILES=AGRIF_YOURFILES

FILENAMES= modbc modmpp modcluster modinit modinitvars modinterp  \
	   modtypes modbcfunction modutil modcurgridfunctions \
	   modmask modsauv modupdate modinterpbasic\
	   modupdatebasic modlinktomodel modarrays
	   
OBJS1=$(addprefix $(AGRIFOBJS)/,$(FILENAMES))
OBJS=$(addsuffix .o,$(OBJS1))


AGRIF_FILES1=$(addprefix $(AGRIFFILES)/,$(FILENAMES))
AGRIF_FILES=$(addsuffix .F,$(AGRIF_FILES1))

all :   conv_exe key $(MODDIR)/libagrif.a $(SXMODEL_LIB)
	@echo AGRIF is OK

conv_exe : 
	@( cd LIB ; $(M_K) -f Makefile C_L="$(C_L)" C_O="$(C_O)" )

key : 
	-@if [ "`cat KEY_CPP`" != "$(P_P)" ] ; then \
	echo "$(P_P)" > KEY_CPP ; fi

$(MODDIR)/libagrif.a : $(OBJS)
	@$(A_C) $@ $(OBJS)

$(SXMODEL_LIB) : 

$(AGRIFOBJS)/%.o : $(AGRIFFILES)/%.F KEY_CPP
	@$(SRM) $(AGRIFYOURFILES)/$(*F).f
	@$(P_C) $(P_O) $(AGRIFFILES)/$(*F).F > $(AGRIFYOURFILES)/$(*F).f
	@$(F_C) $(F_O) $(FFLAGS_EXTEND) $(AGRIF_OPTS) -c $(AGRIFYOURFILES)/$(*F).f -o $(AGRIFOBJS)/$(*F).o


$(AGRIFOBJS)/modarrays.o: $(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modbc.o: $(AGRIFOBJS)/modinterp.o

$(AGRIFOBJS)/modbcfunction.o: $(AGRIFOBJS)/modupdate.o \
			 $(AGRIFOBJS)/modbc.o 

$(AGRIFOBJS)/modcluster.o: $(AGRIFOBJS)/modsauv.o \
			$(AGRIFOBJS)/modinitvars.o \
			$(AGRIFOBJS)/modcurgridfunctions.o

$(AGRIFOBJS)/modcurgridfunctions.o: $(AGRIFOBJS)/modinit.o

$(AGRIFOBJS)/modinit.o: $(AGRIFOBJS)/modlinktomodel.o \
			$(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modinitvars.o: $(AGRIFOBJS)/modlinktomodel.o \
			$(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modinterpbasic.o: $(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modinterp.o: $(AGRIFOBJS)/modcurgridfunctions.o \
			$(AGRIFOBJS)/modmask.o \
			$(AGRIFOBJS)/modarrays.o \
			$(AGRIFOBJS)/modmpp.o \
			$(AGRIFOBJS)/modinterpbasic.o

$(AGRIFOBJS)/modlinktomodel.o: $(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modmask.o: $(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modsauv.o: $(AGRIFOBJS)/modarrays.o \
			$(AGRIFOBJS)/modlinktomodel.o \
			$(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modupdatebasic.o: $(AGRIFOBJS)/modtypes.o

$(AGRIFOBJS)/modupdate.o: $(AGRIFOBJS)/modmask.o \
			$(AGRIFOBJS)/modcurgridfunctions.o \
			$(AGRIFOBJS)/modarrays.o \
			$(AGRIFOBJS)/modmpp.o \
			$(AGRIFOBJS)/modupdatebasic.o

$(AGRIFOBJS)/modutil.o: $(AGRIFOBJS)/modcurgridfunctions.o \
			$(AGRIFOBJS)/modcluster.o

$(AGRIFOBJS)/modmpp.o: $(AGRIFOBJS)/modtypes.o \
			$(AGRIFOBJS)/modarrays.o

clean : FORCE
	$(SRM) $(OBJS) $(MODOBJS)$(LIBDIR)/libagrif.a 
	
cleanall : FORCE
	$(SRM) $(OBJS) $(MODOBJS)$(LIBDIR)/libagrif.a 
	@(cd LIB ; $(M_K) -f Makefile clean )

FORCE :
