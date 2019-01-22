!!DB: 2009.09.09 -- deleted old DBG code 
!!and disable the routine. If you want it then dig up an older version
MODULE trcrad
   !!======================================================================
   !!                       ***  MODULE  trcrad  ***
   !! Ocean passive tracers:  correction of negative concentrations
   !!======================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_rad    : correction of negative concentrations
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and tracers variables
   USE trc                 ! ocean passive tracers variables
   USE lib_mpp
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_rad        ! routine called by trcstp.F90
   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcrad.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_rad( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad  ***
      !!
      !! ** Purpose : "crappy" routine to correct artificial negative
      !!      concentrations due to isopycnal scheme
      !!
      !! ** Method  : Set negative concentrations to zero
      !!              compute the corresponding mass added to the tracers
      !!              and remove it when possible 
      !!
      !! History :
      !!   8.2  !  01-01  (O. Aumont & E. Kestenare)  Original code
      !!   9.0  !  04-03  (C. Ethe)  free form F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::  ji, jj, jk, jn             ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'DB: trc_rad : ROUTINE DISABLED -- find an old one if you need it'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      
   END SUBROUTINE trc_rad

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rad (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rad
#endif
   
   !!======================================================================
END MODULE trcrad
