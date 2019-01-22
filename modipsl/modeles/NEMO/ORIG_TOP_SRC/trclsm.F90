MODULE trclsm
   !!===============================================================
   !!
   !!                       *** MODULE trclsm ****
   !!
   !!  READS specific NAMELIST for sms terms
   !!
   !!=================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trclsm.F90,v 1.4 2005/11/14 15:44:41 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!-----------------------------------------------------------------
#if defined key_passivetrc
   !!-------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE sms


   IMPLICIT NONE                             
   PRIVATE

   !! * Accessibility
   PUBLIC trc_lsm


#if defined key_trc_lobster1
   !!----------------------------------------------------------------------
   !!   'key_trc_lobster1'                        LOBSTER1 biological model  
   !!----------------------------------------------------------------------
#  include "trclsm.lobster1.h90"

#elif defined key_trc_pisces
   !!----------------------------------------------------------------------
   !!   'key_trc_pisces'                            PISCES biological model                  
   !!----------------------------------------------------------------------
#  include "trclsm.pisces.h90"

#elif defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc  '                                          CFC model                  
   !!----------------------------------------------------------------------
#  include "trclsm.cfc.h90"

   !!----------------------------------------------------------------------
   !!   Default option                               
   !!----------------------------------------------------------------------
# endif

#else

CONTAINS

   SUBROUTINE trc_lsm
      !!================
      !! no passive tracers
   END  SUBROUTINE  trc_lsm

#endif  

END MODULE trclsm  
