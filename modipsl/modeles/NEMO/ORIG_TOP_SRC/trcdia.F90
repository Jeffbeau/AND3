MODULE trcdia
   !!==========================================================================
   !!
   !!                       *** MODULE trcdia ***
   !! Output  for tracer concentration  
   !! O.Aumont and A.El Moussaoui 03/05 F90 
   !!==========================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcdia.F90,v 1.4 2005/11/14 15:44:40 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !! * Modules used
 
   USE trcdit

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_dia

   !! * Module variables

CONTAINS

   SUBROUTINE trc_dia(kt,kindic)  
      !!===========================================================================================
      !!
      !!                       ROUTINE trcdii_wr
      !!===========================================================================================

      INTEGER, INTENT( in ) :: kt, kindic

      ! outputs for tracer concentration
      ! -------------------------------- 

      CALL trcdit_wr(kt,kindic)

#if defined key_trc_diatrd

      ! outputs for dynamical trends
      ! ----------------------------

      CALL trcdid_wr(kt,kindic)

#endif
#if defined key_trc_diaadd

      ! outputs for additional arrays
      ! -----------------------------

      CALL trcdii_wr(kt,kindic)

#endif
#if defined key_trc_diabio

      ! outputs for biological trends
      ! -----------------------------

      CALL trcdib_wr(kt,kindic)

#endif

   END SUBROUTINE trc_dia

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE trc_dia
      
   END SUBROUTINE trc_dia   
#endif

END MODULE trcdia
