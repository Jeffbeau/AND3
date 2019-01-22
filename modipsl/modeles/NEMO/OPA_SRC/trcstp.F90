MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !!  Dummy module
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   Default key                                     No passive tracers
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/trcstp.F90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp

   !!======================================================================
END MODULE trcstp
