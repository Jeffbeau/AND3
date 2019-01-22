MODULE flxmod
   !!======================================================================
   !!                       ***  MODULE  flxmod  ***
   !! Ocean forcing:  thermohaline forcing of the ocean
   !!=====================================================================
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE daymod          ! calendar
   USE ocfzpt          ! ocean freezing point

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC flx       ! routine called by step.F90
   PUBLIC flx_init  ! routine called by opa.F90

   !! * local declarations
   REAL(wp), PUBLIC ::            & !!! surface fluxes namelist (namflx)
      q0    = 0.e0,               &  ! net heat flux
      qsr0  = 0.e0,               &  ! solar heat flux
      emp0  = 0.e0,               &  ! net freshwater flux
      dqdt0 = -40.,               &  ! coefficient for SST damping (W/m2/K)
      deds0 = 27.7                   ! coefficient for SSS damping (mm/day)
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxmod.F90,v 1.4 2006/04/19 14:43:16 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

#if defined key_flx_bulk_monthly
   !!----------------------------------------------------------------------
   !!   'key_flx_bulk_monthly'   and                           MONTHLY bulk
   !!   Default option                                         Net CDF file
   !!----------------------------------------------------------------------
#  include "flx_bulk_monthly.h90"

#elif defined key_flx_bulk_daily
   !!----------------------------------------------------------------------
   !!   'key_flx_bulk_daily'                                     DAILY bulk
   !!                                                          Net CDF file
   !!----------------------------------------------------------------------

!!DB 2008.06.27:
# if defined key_flx_bulk_cmc
   !!----------------------------------------------------------------------
   !!   'key_flx_bulk_cmc'                                     CMC bulk
   !!                                                          Net CDF file
   !!----------------------------------------------------------------------
#  include "flx_bulk_cmc.h90"        
!!DB 2009.08.06
#elif defined key_CORE_NY || defined key_CORE_ANNUAL
!! CORE Normal year forcing or 1958 ... 2006 forcing
!! DB modified flx_core.h90, using DB prepared SOPA12 files
#  include "db_flx_core.h90"
# else
#  include "flx_bulk_daily.h90"
# endif

#elif defined key_flx_forced_daily
   !!----------------------------------------------------------------------
   !!   'key_flx_forced_daily'                                 DAILY fluxes
   !!                                                          Net CDF file
   !!----------------------------------------------------------------------
#  include "flx_forced_daily.h90"

#elif defined key_flx_forced_monthly 
!byoung
#   include "flx_forced_monthly.h90"

#elif defined key_coupled
# if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_coupled'  and                          Coupled Ocan/Atmosphere
   !!   'key_ice_lim'                               with  LIM sea-ice model
   !!----------------------------------------------------------------------
#  include "flx_coupled_ice.h90"

# else
   !!----------------------------------------------------------------------
   !!   'key_flx_coupled'  and                      Coupled Ocan/Atmosphere
   !!   Default option                              without   sea-ice model
   !!----------------------------------------------------------------------
#  include "flx_coupled_noice.h90"

# endif
#else
   !!----------------------------------------------------------------------
   !!   Default option                                   Analytical forcing
   !!----------------------------------------------------------------------
   !!   flx          : define the thermohaline fluxes for the ocean
   !!----------------------------------------------------------------------
   !! * Module used
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flx ( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE flx  ***
      !!              
      !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
      !!      to the ocean at each time step.
      !!
      !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
      !!
      !! ** Action  : - q, qt, qsr, emp, emps, qrp, erp
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !!---------------------------------------------------------------------

      IF( kt == nit000 )  THEN
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)'flx     : Analytical/Constant surface fluxes'
         IF(lwp) WRITE(numout,*)'~~~~~~~ '
         IF(lwp) WRITE(numout,*)'          See the routine oce_sbc'     
         IF(lwp) WRITE(numout,*)' '
      ENDIF

   END SUBROUTINE flx

#endif


   SUBROUTINE flx_init
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE flx  ***
      !!              
      !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
      !!      to the ocean at each time step.
      !!
      !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
      !!
      !! ** Action  : - q, qt, qsr, emp, emps, qrp, erp
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      NAMELIST/namflx/ q0, qsr0, emp0, dqdt0, deds0
      !!---------------------------------------------------------------------

      ! Read Namelist namflx : surface thermohaline fluxes
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namflx )

      IF(lwp) THEN
         WRITE(numout,*)' '
         WRITE(numout,*)'flx_init : thermohaline forcing '
         WRITE(numout,*)'~~~~~~~~ '
         WRITE(numout,*)'           net heat flux                  q0   = ', q0  , ' W/m2'
         WRITE(numout,*)'           solar heat flux                qsr0 = ', qsr0, ' W/m2'
         WRITE(numout,*)'           net heat flux                  emp0 = ', emp0, ' W/m2'
         WRITE(numout,*)'           coefficient for SST damping   dqdt0 = ', dqdt0,' W/m2/K'
         WRITE(numout,*)'           coefficient for SSS damping   deds0 = ', deds0,' mm/day'
      ENDIF

   END SUBROUTINE flx_init

END MODULE flxmod
