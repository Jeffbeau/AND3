MODULE sms
   !!======================================================================
   !!                        ***  sms  ***
   !! passive tracers :   set the passive tracers variables
   !!======================================================================
   !! History :
   !!   9.0  !  04-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/sms.F90,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
#if defined key_passivetrc

   USE par_oce
   USE par_trc
   USE par_sms
   IMPLICIT NONE

#if defined key_trc_lobster1
   !!----------------------------------------------------------------------
   !!   'key_trc_lobster1'                        LOBSTER1 biological model  
   !!----------------------------------------------------------------------
#  include "sms_lobster1.h90"

#elif defined key_trc_pisces
   !!----------------------------------------------------------------------
   !!   'key_trc_pisces'                            PISCES biological model                  
   !!----------------------------------------------------------------------
#  include "sms_pisces.h90"

#elif defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc  '                                          CFC model                  
   !!----------------------------------------------------------------------
#  include "sms_cfc.h90"

#endif

#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE sms
