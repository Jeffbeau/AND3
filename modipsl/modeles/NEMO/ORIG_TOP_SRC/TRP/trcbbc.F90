MODULE trcbbc
   !!==============================================================================
   !!                       ***  MODULE  trcbbc  ***
   !! Ocean passive tracers:  bottom boundary condition
   !!==============================================================================
#if   defined key_passivetrc && defined key_trcbbc
   !!----------------------------------------------------------------------
   !!   'key_trcbbc'                                  geothermal heat flux
   !!----------------------------------------------------------------------
   !!   trc_bbc      : update the tracer trend at ocean bottom 
   !!   trc_bbc_init : initialization of geothermal heat flux trend
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE prtctl_trc          ! Print control for debbuging
 
   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_bbc          ! routine called by trcstp.F90

   !! to be transfert in the namelist ???!   
   LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbc = .TRUE.   !: bbc flag

   !! * Module variables
   INTEGER ::                       & !!! ** bbc namelist (nambbc) **
      ngeo_trc_flux = 1                    ! Geothermal flux (0:no flux, 1:constant flux,
      !                                !                  2:read in file )
   REAL(wp) ::                      & !!! ** bbc namlist **
      ngeo_trc_flux_const = 86.4e-3        ! Constant value of geothermal heat flux

   INTEGER, DIMENSION(jpi,jpj) ::   &
      nbotlevt                         ! ocean bottom level index at T-pt
   REAL(wp), DIMENSION(jpi,jpj) ::  &
      qgh_trd                          ! geothermal heating trend
 
   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcbbc.F90,v 1.9 2005/12/12 14:18:10 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_bbc( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbc  ***
      !!
      !! ** Purpose :   Compute the bottom boundary contition on passive tracer 
      !!      associated with geothermal heating and add it to the general
      !!      trend of tracers equations.
      !!
      !! ** Method  :   The geothermal heat flux set to its constant value of 
      !!       86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
      !!       The temperature trend associated to this heat flux through the
      !!       ocean bottom can be computed once and is added to the temperature
      !!       trend juste above the bottom at each time step:
      !!            tra = tra + Qsf / (rau0 rcp e3T) for k= mbathy -1
      !!       Where Qsf is the geothermal heat flux.
      !!
      !! ** Action  : - update the temperature trends tra with the trend of
      !!                the ocean bottom boundary condition
      !!
      !! References :
      !!      Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
      !!
      !! History :
      !!   8.1  !  99-10  (G. Madec)  original code
      !!   8.5  !  02-08  (G. Madec)  free form + modules
      !!   9.0  !  04-03  (C. Ethe)  adpated for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index

      !! * Local declarations
#if defined key_vectopt_loop   &&   ! defined key_autotasking
      INTEGER ::   ji, jn                  ! dummy loop indices
#else
      INTEGER ::   ji, jj, jn              ! dummy loop indices
#endif
      REAL(wp) ::   ztra                ! temporary scalar
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      ! 0. Initialization
      IF( kt == nittrc000 )   CALL trc_bbc_init

      ! 1. Add the geothermal heat flux trend on temperature

      SELECT CASE ( ngeo_trc_flux )

      CASE ( 1:2 )                !  geothermal heat flux

         DO jn = 1, jptra
#if defined key_vectopt_loop   &&   ! defined key_autotasking
            DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
               tra(ji,1,nbotlevt(ji,1),jn) = tra(ji,1,nbotlevt(ji,1),jn) + qgh_trd(ji,1)
            END DO
#else
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  tra(ji,jj,nbotlevt(ji,jj),jn) = tra(ji,jj,nbotlevt(ji,jj),jn) + qgh_trd(ji,jj)
               END DO
            END DO
#endif
         END DO

         IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
            WRITE(charout, FMT="('bbc')")
            CALL prt_ctl_trc_info(charout)
            CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
         ENDIF
      END SELECT

   END SUBROUTINE trc_bbc


   SUBROUTINE trc_bbc_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbc_init  ***
      !!
      !! ** Purpose :   Compute once for all the trend associated with geo-
      !!      thermal heating that will be applied at each time step at the
      !!      bottom ocean level
      !!
      !! ** Method  :   Read the namtrabbc namelist and check the parameters.
      !!      called at the first time step (nittrc000)
      !!
      !! ** Input   : - Namlist namtrcbbc
      !!              - NetCDF file  : passivetrc_geothermal_heating.nc 
      !!                               ( if necessary )
      !!
      !! ** Action  : - compute the heat geothermal trend qgh_trd
      !!              - compute the bottom ocean level nbotlevt
      !!
      !! history :
      !!  8.5  ! 02-11 (A. Bozec) original code
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * local declarations
      CHARACTER (len=32) ::   clname
      INTEGER  ::   ji, jj              ! dummy loop indices
      INTEGER  ::   inum = 11           ! temporary logical unit
      INTEGER  ::   itime               ! temporary integers
      REAL(wp) ::   zdate0, zdt         ! temporary scalars
      REAL(wp), DIMENSION(1) :: zdept   ! temporary workspace
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
         zlamt, zphit, zdta   ! temporary workspace

      NAMELIST/namtrcbbc/ngeo_trc_flux, ngeo_trc_flux_const 
      !!----------------------------------------------------------------------

      ! Read Namelist nambbc : bottom momentum boundary condition
      REWIND ( numnamtra )
      READ   ( numnamtra, namtrcbbc )

      ! Control print
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'trc_bbc : Passive tracers Bottom Boundary Condition (bbc)'
      IF(lwp) WRITE(numout,*) '~~~~~~~   Geothermal heatflux'
      IF(lwp) WRITE(numout,*) '          Namelist namtrcbbc : set bbc parameters'
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '             Geothermal flux           ngeo_trc_flux       = ', ngeo_trc_flux
      IF(lwp) WRITE(numout,*) '             Constant geothermal flux  ngeo_trc_flux_const = ', ngeo_trc_flux_const
      IF(lwp) WRITE(numout,*)

      ! level of the ocean bottom at T-point

      DO jj = 1, jpj
         DO ji = 1, jpi
            nbotlevt(ji,jj) = MAX( mbathy(ji,jj)-1, 1 )
         END DO
      END DO

      ! initialization of geothermal heat flux

      SELECT CASE ( ngeo_trc_flux )

      CASE ( 0 )                ! no geothermal heat flux
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '             *** no geothermal heat flux'

      CASE ( 1 )                ! constant flux
         IF(lwp) WRITE(numout,*) '             *** constant heat flux  =   ', ngeo_trc_flux_const
         qgh_trd(:,:) = ngeo_trc_flux_const

      CASE ( 2 )                ! variable geothermal heat flux
         ! read the geothermal fluxes in mW/m2
         clname = 'passivetrc_geothermal_heating'
         itime = 1
         zlamt(:,:) = 0.
         zphit(:,:) = 0.
         IF(lwp) WRITE(numout,*) '             *** variable geothermal heat flux read in ', clname, ' file'
         CALL restini( clname, jpidta, jpjdta, zlamt, zphit, 1, zdept , clname,   &
                       itime, zdate0, zdt, inum , domain_id=nidom )
         CALL restget( inum, 'heatflow', jpidta, jpjdta, 1, 0, .FALSE., zdta )
         DO jj = 1, nlcj
            DO ji = 1, nlci
              qgh_trd(ji,jj) = zdta(mig(ji),mjg(jj))
            END DO
         END DO

         CALL restclo( inum )
         qgh_trd(:,:) = qgh_trd(:,:) * 1.e-3 ! conversion in W/m2

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '     bad flag value for ngeo_trc_flux = ', ngeo_trc_flux
         nstop = nstop + 1

      END SELECT

      ! geothermal heat flux trend

      SELECT CASE ( ngeo_trc_flux )

      CASE ( 1:2 )                !  geothermal heat flux

#if defined key_vectopt_loop   &&   ! defined key_autotasking
         DO ji = 1, jpij   ! vector opt. (forced unrolling)
            qgh_trd(ji,1) = ro0cpr * qgh_trd(ji,1) / fse3t(ji,1,nbotlevt(ji,1) )
         END DO
#else
         DO jj = 1, jpj
            DO ji = 1, jpi
               qgh_trd(ji,jj) = ro0cpr * qgh_trd(ji,jj) / fse3t(ji,jj,nbotlevt(ji,jj))
            END DO
         END DO
#endif

      END SELECT

   END SUBROUTINE trc_bbc_init

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbc = .FALSE.  !: bbc flag
CONTAINS
   SUBROUTINE trc_bbc( kt )           ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_bbc: You should not have seen this print! error?', kt
   END SUBROUTINE trc_bbc
#endif

   !!======================================================================
END MODULE trcbbc
