!!DB: 2009.08.31 -- Eliminated GYRE config
MODULE taumod
   !!======================================================================
   !!                       ***  MODULE  taumod  ***
   !! Ocean forcing : stress at the the ocean surface
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   tau          : define the surface stress for the ocean
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar
   USE lbclnk          ! 

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tau                ! routine called by step.F90

   !! * Share modules variables

!!DB: NEW 2009.05.05 -- required for core forcing and (eventually) updated lim-model
! FD correction from JM Molines
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
      taux, tauy,  &      !: surface stress components in (i,j) referential
      tauxwo, tauywo,  &  !: surface stress components in (i,j) referential air-ocean
      tauxwi, tauywi,  &  !: surface stress components in (i,j) referential air-ice
      tauxg, tauyg        !: surface stress components in geographical
      !                   !  referential (used in output)

!!DB OLD
!   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
!      taux, tauy,      &  !: surface stress components in (i,j) referential
!      tauxg, tauyg        !: surface stress components in geographical
      !                   !  referential (used in output)

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/taumod.F90,v 1.7 2006/04/10 15:46:11 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

#if defined key_tau_monthly
   ! Monthly climatology in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
   !!----------------------------------------------------------------------
   !!   'key_tau_monthly'                        MONTHLY climatology stress
   !!   default case                                   NetCDF files
   !!----------------------------------------------------------------------
#   include "tau_forced_monthly.h90"

# elif defined key_tau_daily
   !!----------------------------------------------------------------------
   !!   'key_tau_daily'                                 DAILY stress
   !!                                                   NetCDF files
   !!----------------------------------------------------------------------
   ! Daily climatology/interannual in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
!#   include "tau_forced_daily.h90"
#   include "tau_forcing.h90"

# elif defined key_tau_cmc
   !!----------------------------------------------------------------------
   !!   'key_tau_cmc'                                   CMC WS stress
   !!                                                  default: 3 hourly for 6 days
   !!                                                   NetCDF files
   !!----------------------------------------------------------------------
   ! WS in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
!#   include "tau_forced_cmc.h90"
#   include "tau_forcing.h90"
#elif defined key_coupled
   ! Coupled case : stress at the coupling frequency
# if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_coupled'                              Coupled Ocean/Atmosphere
   !!   'key_ice_lim'                                   LIM sea-ice
   !!----------------------------------------------------------------------
   ! New way: 3D referential link to the earth (avoid north pole pb)
   ! (3 component stress defined at U- and V-points)
#  include "tau_coupled_ice.h90"
# else
   !!----------------------------------------------------------------------
   !!   'key_coupled'                              Coupled Ocean/Atmosphere
   !!   Default case                                  NO sea-ice
   !!----------------------------------------------------------------------
   ! old fashion: geographical referential
   ! (zonal and meridional stress defined at U- and V-points)
#  include "tau_coupled.h90"
# endif
#else
   !!----------------------------------------------------------------------
   !!   Default option                                     constant forcing
   !!----------------------------------------------------------------------
   !! * local modules variables
   INTEGER  ::       & !!! * Namelist numtau *
      ntau000 = 1       ! nb of time-step during which the surface stress
      !                 ! increase from 0 to its nominal value (taudta) (>0)
   REAL(wp) ::       & !!! * Namelist numtau *
      tau0x = 0.e0 , &  ! constant wind stress value in i-direction
      tau0y = 0.e0      ! constant wind stress value in j-direction
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tau( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE tau  ***
      !! 
      !! ** Purpose :   provide the ocean surface stress at each time step
      !!
      !! ** Method  :   Constant surface stress increasing from 0 to taudta 
      !!      value during the first ntau000 time-step (namelist)
      !!        CAUTION: never mask the surface stress field !
      !!
      !! ** Action  : - update taux , tauy the stress in (i,j) ref.
      !!              - update tauxg, tauyg the stress in geographic ref.
      !!
      !! History :
      !!   4.0  !  91-03  (G. Madec)  Original code
      !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  ) ::   kt    ! ocean time step
      REAL(wp) ::   ztau, ztau_sais, &  ! wind intensity and of the seasonal cycle
         ztime,                      &  ! time in hour
         ztimemax, ztimemin,         &  ! 21th June, and 21th decem. if date0 = 1st january
         ztaun                          ! intensity 
      INTEGER  ::   ji, jj              ! dummy loop indices

      INTEGER  ::           &
         zyear0,            &           ! initial year
         zmonth0,           &           ! initial month
         zday0,             &           ! initial day
         zday_year0                    ! initial day since january 1st
        

      !! * Local declarations
      REAL(wp) ::   zfacto              ! 

      NAMELIST/namtau/ ntau000, tau0x, tau0y
      !!---------------------------------------------------------------------

         IF( kt == nit000 ) THEN
   
            ! Read Namelist namtau : surface wind stress
            ! --------------------
            REWIND ( numnam )
            READ   ( numnam, namtau )
   
            IF(lwp) WRITE(numout,*)' '
            IF(lwp) WRITE(numout,*)' tau     : Constant surface wind stress read in namelist'
            IF(lwp) WRITE(numout,*)' ~~~~~~~ '
            IF(lwp) WRITE(numout,*)'           Namelist namtau: set the constant stress values'
            IF(lwp) WRITE(numout,*)'              spin up of the stress  ntau000 = ', ntau000, ' time-steps'
            IF(lwp) WRITE(numout,*)'              constant i-stress      tau0x   = ', tau0x  , ' N/m2'
            IF(lwp) WRITE(numout,*)'              constant j-stress      tau0y   = ', tau0y  , ' N/m2'
   
            ntau000 = MAX( ntau000, 1 )   ! must be >= 1
   
         ENDIF
   
         ! Increase the surface stress to its nominal value in ntau000 time-step
         
         IF( kt <= ntau000 ) THEN
            zfacto = 0.5 * (  1. - COS( rpi * FLOAT( kt ) / FLOAT( ntau000 ) )  )
            taux (:,:) = zfacto * tau0x
            tauy (:,:) = zfacto * tau0y
            tauxg(:,:) = zfacto * tau0x
            tauyg(:,:) = zfacto * tau0y
         ENDIF
      
   END SUBROUTINE tau
#endif
   !!======================================================================
END MODULE taumod
