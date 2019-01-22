MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc              ! ocean passive tracers variables 
   USE trctrp           ! passive tracers transport
   USE trcsms           ! passive tracers sources and sinks
   USE prtctl_trc       ! Print control for debbuging
   USE trcdia
   USE trcdit
   USE trcrst

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_stp           ! called by step
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcstp.F90,v 1.11 2006/04/11 13:49:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_stp( kt, kindic )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe)  Original
      !!-------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      INTEGER, INTENT( in ) ::  kindic
      CHARACTER (len=25) :: charout

      ! this ROUTINE is called only every ndttrc time step
      IF( MOD( kt , ndttrc ) /= 0 ) RETURN

      ! tracers: sink and source 
      IF(ln_ctl) THEN
         WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
         CALL prt_ctl_trc_info(charout)
      ENDIF


      CALL trc_sms( kt )

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('SMS')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      ! transport of passive tracers
      CALL trc_trp( kt )

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('TRP')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      CALL trc_wri( kt )            ! outputs

      CALL trc_dia( kt, kindic )     ! diagnostics


   END SUBROUTINE trc_stp

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif

   !!======================================================================
END MODULE trcstp
