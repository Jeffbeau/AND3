MODULE par_sms
   !!---------------------------------------------------------------------
   !!
   !!                         PARAMETER SMS
   !!                       *******************************
   !!
   !!  purpose :
   !!  ---------
   !!     INCLUDE PARAMETER FILE for SMS  models
   !!
   !!
   !!----------------------------------------------------------------------
   !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/par_sms.F90,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   USE par_trc_trp
   IMPLICIT NONE

#if defined key_trc_lobster1
   !!----------------------------------------------------------------------
   !!   'key_trc_lobster1'                        LOBSTER1 biological model  
   !!----------------------------------------------------------------------
#  include "par_sms_lobster1.h90"

#elif defined key_trc_pisces
   !!----------------------------------------------------------------------
   !!   'key_trc_pisces'                            PISCES biological model                  
   !!----------------------------------------------------------------------
#  include "par_sms_pisces.h90"

#elif defined key_cfc
   !!----------------------------------------------------------------------
   !!   'key_cfc  '                                          CFC model                  
   !!----------------------------------------------------------------------
#  include "par_sms_cfc.h90"

#else
   !!  purpose :
   !!  ---------
   !!     No SMS  models
#endif

END MODULE par_sms
