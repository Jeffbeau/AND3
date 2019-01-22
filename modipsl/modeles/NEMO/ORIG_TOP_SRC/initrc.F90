MODULE initrc
   !!================================================
   !!
   !!                       *** MODULE initrc ***
   !! Initialisation the tracer model
   !!================================================
                                                                                                                            
#if defined key_passivetrc

   !!-------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/initrc.F90,v 1.5 2005/11/16 16:30:06 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!-------------------------------------------------------

   !!--------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE trcrst
   USE trcctl
   USE trclec
   USE trcdtr
   USE trcini
   USE prtctl_trc      ! Print control passive tracers (prt_ctl_trc_init routine)
   
   IMPLICIT NONE
   PRIVATE
   
   
   !! * Accessibility
   PUBLIC ini_trc
   
CONTAINS
   
   SUBROUTINE ini_trc
      !!---------------------------------------------------------------------
      !!
      !!                       ROUTINE ini_trc
      !!                     ******************
      !!
      !!  PURPOSE :
      !!  ---------
      !!     initialize the tracer model
      !!
      !!   METHOD :
      !!   -------
      !!
      !!
      !!   History:
      !!   -------
      !!      original  : 91-03 ()
      !!      additions : 92-01 (C. Levy)
      !!                  05-03 (O. Aumont and A. El Moussaoui) F90
      !!                  05-10 (C. Ethe ) print control initialization 
      !!----------------------------------------------------------------------

      !!---------------------------------------------------------------------
      !!  OPA.9, 03-2005
      !!---------------------------------------------------------------------

      !! 0.b PRINT the number of tracer
      !! ------------------------------

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
      IF(lwp) WRITE(numout,*) ' '

      ! 1. READ passive tracers namelists
      ! ---------------------------------

      CALL trc_lec

      ! 2. control consistency between parameters, cpp key and namelists
      ! ----------------------------------------------------------------

      CALL trc_ctl

      ! 3. computes some initializations
      ! --------------------------------

      CALL trc_ini

      ! 4. restart from a FILE (nutrst)
      ! ----------------------

      IF( lrsttr ) THEN


         CALL trc_rst

      ELSE

         ! start from anything ELSE

         CALL trc_dtr

      ENDIF

      ! 5. Print control
      !------------------

      IF( ln_ctl )    CALL prt_ctl_trc_init

   END SUBROUTINE ini_trc

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE ini_trc
      
   END SUBROUTINE ini_trc
#endif

END MODULE initrc 
