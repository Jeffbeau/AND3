MODULE initrc
   !!======================================================================
   !!                       ***  MODULE initrc  ***
   !!  Dummy module
   !! Initialization for the tracer model
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   Default key                                     No passive tracers
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/initrc.F90,v 1.1 2005/09/14 09:38:16 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ini_trc        ! Empty routine
      WRITE(*,*) 'ini_trc: You should not have seen this print! error?'
   END SUBROUTINE ini_trc

   !!======================================================================
END MODULE initrc
