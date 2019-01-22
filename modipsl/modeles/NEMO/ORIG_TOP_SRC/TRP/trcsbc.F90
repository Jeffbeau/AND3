MODULE trcsbc
   !!==============================================================================
   !!                       ***  MODULE  trcsbc  ***
   !! Ocean passive tracers:  surface boundary condition
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean  passive tracers variables
   USE prtctl_trc          ! Print control for debbuging


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_sbc              ! routine called by step.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcsbc.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (concentration/dilution effect) and add it to the general 
      !!       trend of tracer equations.
      !!
      !! ** Method :
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer as :
      !!            tra = tra + emp * trn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation minus runoff) given in kg/m2/s is divided
      !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
      !!
      !! ** Action  : - Update the 1st level of tra with the trend associated
      !!                with the tracer surface boundary condition 
      !!
      !! History :
      !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
      !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jn           ! dummy loop indices
      REAL(wp) ::   ztra, zsrau, zse3t   ! temporary scalars
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_sbc : Passive tracers surface boundary condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! 0. initialization
      zsrau = 1. / rauw
#if ! defined key_s_coord
      zse3t = 1. / fse3t(1,1,1)
#endif

      DO jn = 1, jptra
         ! 1. Concentration dillution effect on tra
         DO jj = 2, jpj
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord
               zse3t = 1. / fse3t(ji,jj,1)
#endif
               ! concent./dilut. effect
               ztra = emps(ji,jj) * zsrau * trn(ji,jj,1,jn) * zse3t * tmask(ji,jj,1)
               
               ! add the trend to the general tracer trend
               tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra
            END DO
         END DO
         
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sbc')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_sbc

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_sbc (kt)              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sbc
#endif
   
   !!======================================================================
END MODULE trcsbc
