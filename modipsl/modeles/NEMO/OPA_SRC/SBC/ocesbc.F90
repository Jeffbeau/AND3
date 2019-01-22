!!DB: 2009.08.31 -- Eliminated GYRE config
MODULE ocesbc
   !!======================================================================
   !!                     ***  MODULE  ocesbc  ***
   !!                     Ocean surface boundary conditions
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   oce_sbc     : ???
   !!   oce_sbc_dmp : ???
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce            ! dynamics and tracers variables
   USE dom_oce        ! ocean space domain variables
   USE cpl_oce        ! coupled ocean-atmosphere variables
   USE ice_oce        ! sea-ice variable
   USE blk_oce        ! bulk variables
   USE flx_oce        ! sea-ice/ocean forcings variables
   USE phycst         ! Define parameters for the routines
   USE taumod         ! surface stress forcing
   USE flxmod         ! thermohaline fluxes
   USE flxrnf         ! runoffs forcing
   USE tradmp         ! damping salinity trend
   USE dtatem         ! ocean temperature data
   USE dtasal         ! ocean salinity data
   USE ocfzpt         ! surface ocean freezing point
   USE lbclnk         ! ocean lateral boundary condition
   USE lib_mpp        ! distribued memory computing library
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC oce_sbc    ! routine called by step

   !! * Shared module variables
   REAL(wp), PUBLIC ::   &  !:
      aplus, aminus,     &  !:
      empold = 0.e0         !: current year freshwater budget correction
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      qt  ,         &  !: total surface heat flux (w/m2)
      qsr ,         &  !: solar radiation (w/m2)
      emp ,         &  !: evaporation minus precipitation (kg/m2/s = mm/s)
      emps,         &  !: evaporation - precipitation (free surface)
      qrp ,         &  !: heat flux damping (w/m2)
      erp              !: evaporation damping (kg/m2/s = mm/s)
#if ! defined key_dynspg_rl
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      dmp              !: internal dampind term
#endif

#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/ocesbc.F90,v 1.16 2006/04/19 14:43:16 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
# if defined key_coupled
      !!----------------------------------------------------------------------
      !!   'key_coupled' :                            Coupled Ocean/Atmosphere
      !!----------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE oce_sbc  ***
      !!                    
      !! ** Purpose :   Ocean surface boundaries conditions with 
      !!        Louvain la Neuve Sea Ice Model in coupled mode
      !!
      !! History :
      !!   1.0  !  00-10  (O. Marti)  Original code
      !!   2.0  !  02-12  (G. Madec)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * Local declarations
      INTEGER ::   ji, jj                   ! dummy loop indices
      REAL(wp) ::   ztx, ztaux, zty, ztauy
      REAL(wp) ::   ztdta, ztgel, zqrp
      !!----------------------------------------------------------------------
 
      ! 1. initialization to zero at kt = nit000
      ! ---------------------------------------
      
      IF( kt == nit000 ) THEN     
         qsr   (:,:) = 0.e0
         freeze(:,:) = 0.e0
         qt    (:,:) = 0.e0
         qrp   (:,:) = 0.e0
         emp   (:,:) = 0.e0
         emps  (:,:) = 0.e0
         erp   (:,:) = 0.e0
#if ! defined key_dynspg_rl 
         dmp   (:,:) = 0.e0
#endif
      ENDIF

      IF( MOD( kt-1, nfice ) == 0 ) THEN 

         CALL oce_sbc_dmp   ! Computation of internal and evaporation damping terms       

         ! Surface heat flux (W/m2)
         ! -----------------------

         ! restoring heat flux
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztgel = fzptn(ji,jj)
#if defined key_dtasst
               ztdta = MAX( sst(ji,jj),    ztgel )
#else
               ztdta = MAX( t_dta(ji,jj,1), ztgel )
#endif
               zqrp = dqdt0 * ( tb(ji,jj,1) - ztdta )

               qrp(ji,jj) = (1.0-freeze(ji,jj) ) * zqrp
            END DO
         END DO

         ! non solar heat flux + solar flux + restoring
         qt (:,:) = fnsolar(:,:) + fsolar(:,:) + qrp(:,:)

         ! solar flux
         qsr(:,:) = fsolar(:,:)

#if ! defined key_dynspg_rl  
         ! total concentration/dilution effect (use on SSS)
         emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)

         ! total volume flux (use on sea-surface height)
         emp (:,:) = fmass(:,:)  -  dmp(:,:) + runoff(:,:) + erp(:,:)
#else
         ! Rigid-lid (emp=emps=E-P-R+Erp) 
         ! freshwater flux
         emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)
         emp (:,:) = emps(:,:)
#endif

         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vertor opt.
               ztx   = 0.5 * ( freeze(ji+1,jj) + freeze(ji+1,jj+1) )
               ztaux = 0.5 * ( ftaux (ji+1,jj) + ftaux (ji+1,jj+1) )
               taux(ji,jj) = (1.0-ztx) * taux(ji,jj) + ztx * ztaux

               zty   = 0.5 * ( freeze(ji,jj+1) + freeze(ji+1,jj+1) )
               ztauy = 0.5 * ( ftauy (ji,jj+1) + ftauy (ji+1,jj+1) )
               tauy(ji,jj) = (1.0-zty) * tauy(ji,jj) + zty * ztauy
            END DO
         END DO
         CALL lbc_lnk( taux, 'U', -1. )
         CALL lbc_lnk( tauy, 'V', -1. )    

         ! Re-initialization of fluxes
         sst_io(:,:) = 0.e0
         sss_io(:,:) = 0.e0
         u_io  (:,:) = 0.e0
         v_io  (:,:) = 0.e0
         gtaux (:,:) = 0.e0
         gtauy (:,:) = 0.e0

      ENDIF

   END SUBROUTINE oce_sbc

# elif defined key_flx_bulk_monthly || defined key_flx_bulk_daily
      !!----------------------------------------------------------------------
      !!   'key_ice_lim'                              with  LIM sea-ice model
      !!----------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc  ***
      !!                    
      !! ** Purpose : - Ocean surface boundary conditions with LIM sea-ice
      !!        model in forced mode using bulk formulea
      !!
      !! History :
      !!   1.0  !  99-11  (M. Imbard)  Original code
      !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
      !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * Local declarations
      INTEGER  ::   ji, jj                   ! dummy loop indices
      REAL(wp) ::   ztx, ztaux, zty, ztauy
      !!----------------------------------------------------------------------

      ! 1. initialization to zero at kt = nit000
      ! ---------------------------------------
      
      IF( kt == nit000 ) THEN     
         qsr    (:,:) = 0.e0
         qt     (:,:) = 0.e0
         qrp    (:,:) = 0.e0
         emp    (:,:) = 0.e0
         emps   (:,:) = 0.e0
         erp    (:,:) = 0.e0
#if ! defined key_dynspg_rl 
         dmp    (:,:) = 0.e0
#endif
      ENDIF

      IF( MOD( kt-1, nfice ) == 0 ) THEN

         CALL oce_sbc_dmp       ! Computation of internal and evaporation damping terms       

         ! Surface Ocean fluxes
         ! ====================

         ! Surface heat flux (W/m2)
         ! -----------------

         qt  (:,:) = fnsolar(:,:) + fsolar(:,:)     ! non solar heat flux + solar flux
         qsr (:,:) = fsolar(:,:)                     ! solar flux

#if ! defined key_dynspg_rl     
         ! total concentration/dilution effect (use on SSS)
         emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:) + empold

         ! total volume flux (use on sea-surface height)
         emp (:,:) = fmass(:,:) -   dmp(:,:) + runoff(:,:) + erp(:,:) + empold      
#else
         ! Rigid-lid (emp=emps=E-P-R+Erp)
         emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)     ! freshwater flux
         emp (:,:) = emps(:,:)

#endif

         ! Surface stress
         ! --------------

         ! update the stress beloww sea-ice area
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vertor opt.
               ztx         = MAX( freezn(ji,jj), freezn(ji,jj+1) )   ! ice/ocean indicator at U- and V-points
               zty         = MAX( freezn(ji,jj), freezn(ji+1,jj) )
               ztaux       = 0.5 *( ftaux(ji+1,jj) + ftaux(ji+1,jj+1) ) ! ice-ocean stress at U- and V-points
               ztauy       = 0.5 *( ftauy(ji,jj+1) + ftauy(ji+1,jj+1) )
               taux(ji,jj) = (1.-ztx) * taux(ji,jj) + ztx * ztaux    ! stress at the ocean surface
               tauy(ji,jj) = (1.-zty) * tauy(ji,jj) + zty * ztauy
            END DO
         END DO

         ! boundary condition on the stress (taux,tauy)
         CALL lbc_lnk( taux, 'U', -1. )
         CALL lbc_lnk( tauy, 'V', -1. )

         ! Re-initialization of fluxes
         sst_io(:,:) = 0.e0
         sss_io(:,:) = 0.e0
         u_io  (:,:) = 0.e0
         v_io  (:,:) = 0.e0

      ENDIF

   END SUBROUTINE oce_sbc

# else
      !!----------------------------------------------------------------------
      !!   Error option               LIM sea-ice model requires bulk formulea
      !!----------------------------------------------------------------------
      This line forced a compilation error
# endif

#else
   !!----------------------------------------------------------------------
   !!   Default option                                 NO LIM sea-ice model
   !!----------------------------------------------------------------------
# if defined key_coupled
      !!----------------------------------------------------------------------
      !!   'key_coupled' :                            Coupled Ocean/Atmosphere
      !!----------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE oce_sbc  ***
      !!                    
      !! ** Purpose :   Ocean surface boundaries conditions in 
      !!                coupled ocean/atmosphere case without sea-ice
      !!
      !! History :
      !!   1.0  !  00-10  (O. Marti)  Original code
      !!   2.0  !  02-12  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE cpl_oce       ! coupled ocean-atmosphere variables

      !! * Arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jf         ! dummy loop indices
      REAL(wp) ::   ztgel,          &  ! temporary scalars
         zice, zhemis, zqrp, zqri,  &  !    "         "
         zq, zqi, zerp, ze, zei, zro   !    "         "
      !!----------------------------------------------------------------------

      ! Compute fluxes
      ! --------------

      DO jj = 1, jpj
         DO ji = 1, jpi

            ztgel = fzptn(ji,jj)   ! local freezing temperature

            ! opa model ice freeze()

            zice = tmask(ji,jj,1)
            IF( tn(ji,jj,1) >=  ztgel ) zice = 0.
            freeze(ji,jj) = zice

            ! hemisphere indicator (=1 north, =-1 south)
            
            zhemis = float(isign(1, mjg(jj)-(jpjglo/2+1)))
            
            ! a) net downward radiative flux qsr()
            ! - AGCM qsrc if no ice
            ! - zero under ice (zice=1)

            qsr(ji,jj) = (1.-zice)*qsrc(ji,jj)*tmask(ji,jj,1)

            ! b) heat flux damping term qrp()
            ! - no damping if no  ice      (zice=0) 
            ! - gamma*min(0,t-tgel) if ice (zice=1)

            zqrp = 0.
            zqri = dqdt0*MIN( 0., tb(ji,jj,1)-ztgel )
            qrp(ji,jj) = ( ( 1. - zice ) * zqrp + zice * zqri ) * tmask(ji,jj,1)


            ! c) net downward heat flux q() = q0 + qrp()
            ! for q0
            ! - AGCM qc if no  ice (zice=0)
            ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if ice (zice=1)
            zq  = qc(ji,jj)
            zqi = -3. + zhemis
            qt(ji,jj) = ( (1.-zice) * zq + zice * zqi ) * tmask(ji,jj,1) + qrp(ji,jj)
            
            ! d) water flux damping term erp()
            ! - no damping
            zerp = 0.
            erp(ji,jj) = zerp
            
            ! e) net upward water flux e() = eo + runoff() + erp()
            ! for e0
            ! - AGCM if no ice (zice=0)
            ! - 1.mm/day if climatological and opa ice (zice=1)
            ze  = ec(ji,jj)
            zei = 1./rday
            zro = runoff(ji,jj)
            emp(ji,jj) = ( ( 1. - zice ) *  ze + zice * zei + zro ) * tmask(ji,jj,1) + erp(ji,jj)
            
            ! f) net upward water flux for the salinity surface 
            !    boundary condition
            emps(:,:) = emp(:,:)

         END DO
      END DO

   END SUBROUTINE oce_sbc

# elif defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_forced_daily
      !!-------------------------------------------------------------------------
      !!   'key_flx_bulk_monthly' or 'key_flx_bulk_daily' or        bulk formulea
      !!   'key_flx_forced_daily'                                or no bulk case 
      !!-------------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc  ***
      !!                    
      !! ** Purpose :   Ocean surface boundary conditions in forced mode
      !!      using either flux or bulk formulation.
      !!
      !! History :
      !!   1.0  !  99-11  (M. Imbard)  Original code
      !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
      !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE daymod                 ! calendar
#if ! defined key_dtasst
      USE dtasst, ONLY : rclice  ! sea surface temperature data
#endif
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
      USE blk_oce                ! bulk variables
#endif
#if defined key_flx_forced_daily
      USE flx_oce                ! sea-ice/ocean forcings variables
#endif

      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER ::   ji, jj        ! dummy loop arguments
      INTEGER ::   i15, ifreq             !      
      REAL(wp) ::  zxy
      REAL(wp) ::  zsice, zqri, zqrp, ztdta, zqrj
      REAL(wp) ::  zq, zqi, zhemis
      REAL(wp), DIMENSION(jpi,jpj) :: zeri, zerps, ziclim
      REAL(wp), DIMENSION(jpi,jpj) :: zqt, zqsr, zemp  
      !!----------------------------------------------------------------------
 
      ! 1. initialization to zero at kt = nit000
      ! ---------------------------------------
      
      IF( kt == nit000 ) THEN     
         qsr    (:,:) = 0.e0
         freeze (:,:) = 0.e0
         qt     (:,:) = 0.e0
         qrp    (:,:) = 0.e0
         emp    (:,:) = 0.e0
         emps   (:,:) = 0.e0
         erp    (:,:) = 0.e0
#if ! defined key_dynspg_rl 
         dmp    (:,:) = 0.e0
#endif
      ENDIF

#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
      ifreq      = nfbulk
      zqt (:,:)  = qsr_oce(:,:) + qnsr_oce(:,:)
      zqsr(:,:)  = qsr_oce(:,:)
      zemp(:,:)  = evap(:,:) - tprecip(:,:)
#endif 

#if defined key_flx_forced_daily
      ifreq      = 1
      zqt (:,:)  = p_qt (:,:)
      zqsr(:,:)  = p_qsr(:,:)
      zemp(:,:)  = p_emp(:,:)
#endif 

      IF( MOD( kt-1, ifreq) == 0 ) THEN
         ! Computation of internal and evaporation damping terms       
         CALL oce_sbc_dmp

         zsice = - 0.04 / 0.8    ! ratio of isohaline compressibility over isotherme compressibility 
                                 ! ( d rho / dt ) / ( d rho / ds )      ( s = 34, t = -1.8 )
         ! Flux computation
         DO jj = 1, jpj
            DO ji = 1, jpi      
               ! climatological ice 
               ziclim(ji,jj) = FLOAT( NINT( rclice(ji,jj,1) ) )

               ! avoid surfreezing point            
               tn(ji,jj,1) = MAX( tn(ji,jj,1), fzptn(ji,jj) )

               ! hemisphere indicator (=1 north, =-1 south)           
               zhemis = FLOAT( isign(1, mjg(jj) - (jpjdta/2+1) ) )

               ! restoring temperature (ztdta >= to local freezing temperature)            
#if defined key_dtasst
               ztdta = MAX( sst(ji,jj),    fzptn(ji,jj) )
#else
               ztdta = MAX( t_dta(ji,jj,1), fzptn(ji,jj) )
#endif

               ! a) net downward radiative flux qsr()           
               qsr(ji,jj) = (1.-ziclim(ji,jj)) * zqsr(ji,jj) * tmask(ji,jj,1)

               ! b) heat flux damping term qrp()
               ! - gamma*(t-tlevitus) if no  climatological ice (ziclim=0)
               ! - gamma*(t-(tgel-1.))  if climatological ice and no opa ice   (ziclim=1 zicopa=0)
               ! - gamma*min(0,t-tgel) if climatological and opa ice (ziclim=1 zicopa=1)

               zqri = dqdt0 * ( tb(ji,jj,1) - ( fzptn(ji,jj) - 1.) )
               zqrj = dqdt0 * MIN( 0., tb(ji,jj,1) - fzptn(ji,jj) )

               qrp(ji,jj) =  ( ziclim(ji,jj) * ( (1 - freeze(ji,jj)) * zqri    &
                 &                                  + freeze(ji,jj)  * zqrj ) ) * tmask(ji,jj,1)

#if ! defined key_flx_bulk_monthly || ! defined key_flx_bulk_daily
               zqrp = dqdt0 * ( tb(ji,jj,1) - ztdta )
               qrp(ji,jj) = qrp(ji,jj) + (1. - ziclim(ji,jj)) * zqrp
# endif

               ! c) net downward heat flux q() = q0 + qrp()
               ! for q0
               ! - ECMWF fluxes if no climatological ice      (ziclim=0)
               ! - qrp if climatological ice and no opa ice   (ziclim=1 zicopa=0)
               ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if climatological and opa ice 
               !                                              (ziclim=1 zicopa=1)
               zq  = zqt(ji,jj)
               zqi = -3. + zhemis
               qt (ji,jj) = ( (1.-ziclim(ji,jj)) * zq   &
                  +ziclim(ji,jj)  * freeze(ji,jj) * zqi )   &
                  * tmask(ji,jj,1)   &
                  + qrp(ji,jj)

            END DO
         END DO

#if ! defined key_dynspg_rl 
         ! Free-surface

         ! Water flux for zero buoyancy flux if no opa ice and ice clim
         zeri(:,:) = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
         zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )

         ! Contribution to sea level:
         ! net upward water flux emp() = e-p + runoff() + erp() + dmp + empold
         emp (:,:) = zemp(:,:)     &   ! e-p data
            &      + runoff(:,:)   &   ! runoff data
            &      + erp(:,:)      &   ! restoring term to SSS data
            &      + dmp(:,:)      &   ! freshwater flux associated with internal damping
            &      + empold            ! domain averaged annual mean correction

         ! Contribution to salinity:
         ! net upward water flux emps() = e-p + runoff() + erp() + zerps + empold
         emps(:,:) = zemp(:,:)     &
            &      + runoff(:,:)   &
            &      + erp(:,:)      &
            &      + zerps(:,:)    &
            &      + empold 
#else
         ! Rigid-lid (emp=emps=E-P-R+Erp)
         ! freshwater flux
         zeri(:,:)  = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
         zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
         emps (:,:) = zemp(:,:)     &
            &       + runoff(:,:)   &
            &       + erp(:,:)      &
            &       + zerps(:,:)
         emp (:,:) = emps(:,:)
#endif  

         ! Boundary condition on emp for free surface option
         ! -------------------------------------------------
         CALL lbc_lnk( emp, 'T', 1. )
      
      ENDIF

   END SUBROUTINE oce_sbc

!----------------------------------------------------------------!byoung
# elif defined key_flx_forced_monthly
      !!-------------------------------------------------------------------------
      !!   'key_flx_bulk_monthly' or 'key_flx_bulk_daily' or        bulk formulea
      !!   'key_flx_forced_daily'                                or no bulk case
      !!-------------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc  ***
      !!
      !! ** Purpose :   Ocean surface boundary conditions in forced mode
      !!      using either flux or bulk formulation.
      !!
      !! History :
      !!   1.0  !  99-11  (M. Imbard)  Original code
      !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
      !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE daymod                 ! calendar
#if ! defined key_dtasst
      USE dtasst, ONLY : rclice  ! sea surface temperature data
#endif
#if defined key_flx_forced_monthly
      USE flx_oce                ! sea-ice/ocean forcings variables
#endif

      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      INTEGER ::   ji, jj        ! dummy loop arguments
      INTEGER ::   i15, ifreq             !
      REAL(wp) ::  zxy
      REAL(wp) ::  zsice, zqri, zqrp, ztdta, zqrj
      REAL(wp) ::  zq, zqi, zhemis, ztrp
      REAL(wp), DIMENSION(jpi,jpj) :: zeri, zerps, ziclim
      REAL(wp), DIMENSION(jpi,jpj) :: zqt, zqsr, zemp
      !!----------------------------------------------------------------------

      ! 1. initialization to zero at kt = nit000
      ! ---------------------------------------

      IF( kt == nit000 ) THEN
         qsr    (:,:) = 0.e0
         freeze (:,:) = 0.e0
         qt     (:,:) = 0.e0
         qrp    (:,:) = 0.e0
         emp    (:,:) = 0.e0
         emps   (:,:) = 0.e0
         erp    (:,:) = 0.e0
#if ! defined key_dynspg_rl
         dmp    (:,:) = 0.e0
#endif
      ENDIF


#if defined key_flx_forced_monthly
      ifreq      = 1
      zqt (:,:)  = p_bqt (:,:)
      zqsr(:,:)  = p_bqsr(:,:)
      zemp(:,:)  = p_bemp(:,:)
#endif

      IF( MOD( kt-1, ifreq) == 0 ) THEN
         ! Computation of internal and evaporation damping terms
         CALL oce_sbc_dmp

         ztrp = -40.             ! restoring terme for temperature (w/m2/k)
         zsice = - 0.04 / 0.8    ! ratio of isohaline compressibility over isotherme compressibility
                                 ! ( d rho / dt ) / ( d rho / ds )      ( s = 34, t = -1.8 )
         ! Flux computation
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! climatological ice
               ziclim(ji,jj) = FLOAT( NINT( rclice(ji,jj,1) ) )

               ! avoid surfreezing point
               tn(ji,jj,1) = MAX( tn(ji,jj,1), fzptn(ji,jj) )

               ! hemisphere indicator (=1 north, =-1 south)
               zhemis = FLOAT( isign(1, mjg(jj) - (jpjdta/2+1) ) )

               ! restoring temperature (ztdta >= to local freezing temperature)
#if defined key_dtasst
               ztdta = MAX( sst(ji,jj),    fzptn(ji,jj) )
#else
               ztdta = MAX( t_dta(ji,jj,1), fzptn(ji,jj) )
#endif

               ! a) net downward radiative flux qsr()
               qsr(ji,jj) = (1.-ziclim(ji,jj)) * zqsr(ji,jj) * tmask(ji,jj,1)

               ! b) heat flux damping term qrp()
               ! - gamma*(t-tlevitus) if no  climatological ice (ziclim=0)
               ! - gamma*(t-(tgel-1.))  if climatological ice and no opa ice   (ziclim=1 zicopa=0)
               ! - gamma*min(0,t-tgel) if climatological and opa ice (ziclim=1 zicopa=1)

               zqri = ztrp * ( tb(ji,jj,1) - ( fzptn(ji,jj) - 1.) )
               zqrj = ztrp * MIN( 0., tb(ji,jj,1) - fzptn(ji,jj) )

               qrp(ji,jj) =  ( ziclim(ji,jj) * ( (1 - freeze(ji,jj)) * zqri    &
                 &                                  + freeze(ji,jj)  * zqrj ) ) * tmask(ji,jj,1)


               ! c) net downward heat flux q() = q0 + qrp()
               ! for q0
               ! - ECMWF fluxes if no climatological ice      (ziclim=0)
               ! - qrp if climatological ice and no opa ice   (ziclim=1 zicopa=0)
               ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if climatological and opa ice
               !                                              (ziclim=1 zicopa=1)
               zq  = zqt(ji,jj)
               zqi = -3. + zhemis
               qt (ji,jj) = ( (1.-ziclim(ji,jj)) * zq   &
                  +ziclim(ji,jj)  * freeze(ji,jj) * zqi )   &
                  * tmask(ji,jj,1)   &
                  + qrp(ji,jj)

            END DO
         END DO

#if ! defined key_dynspg_rl
         ! Free-surface

         ! Water flux for zero buoyancy flux if no opa ice and ice clim
         zeri(:,:) = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
         zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )

         ! Contribution to sea level:
         ! net upward water flux emp() = e-p + runoff() + erp() + dmp + empold
         emp (:,:) = zemp(:,:)     &   ! e-p data
            &      + runoff(:,:)   &   ! runoff data
            &      + erp(:,:)      &   ! restoring term to SSS data
            &      + dmp(:,:)      &   ! freshwater flux associated with internal damping
            &      + empold            ! domain averaged annual mean correction

         ! Contribution to salinity:
         ! net upward water flux emps() = e-p + runoff() + erp() + zerps + empold
         emps(:,:) = zemp(:,:)     &
            &      + runoff(:,:)   &
            &      + erp(:,:)      &
            &      + zerps(:,:)    &
            &      + empold
#else
         ! Rigid-lid (emp=emps=E-P-R+Erp)
         ! freshwater flux
         zeri(:,:)  = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
         zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
         emps (:,:) = zemp(:,:)     &
            &       + runoff(:,:)   &
            &       + erp(:,:)      &
            &       + zerps(:,:)
         emp (:,:) = emps(:,:)
#endif

         ! Boundary condition on emp for free surface option
         ! -------------------------------------------------
         CALL lbc_lnk( emp, 'T', 1. )

      ENDIF

   END SUBROUTINE oce_sbc





!--------------------------------BIO WZL END

# else
      !!----------------------------------------------------------------------
      !!   Default option :                                 Analytical forcing
      !!----------------------------------------------------------------------

   SUBROUTINE oce_sbc( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE oce_sbc  ***
      !!              
      !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
      !!                to the ocean at each time step.
      !!
      !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
      !!
      !! ** Action  : - qt, qsr, emp, emps, qrp, erp
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-05  (A. Koch-Larrouy) Add Gyre configuration 
      !!----------------------------------------------------------------------
      !! * Modules used
      USE flxrnf                       ! ocean runoffs
      USE daymod, ONLY : nyear         ! calendar
      USE dtasss                       ! sea surface salinity data

      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * local declarations
      REAL(wp) ::                   & !!! surface fluxes namelist (namflx)
         q0   = 0.e0,               &  ! net heat flux
         qsr0 = 0.e0,               &  ! solar heat flux
         emp0 = 0.e0                   ! net freshwater flux
      REAL(wp) ::   zsrp,           &
         zemp_S, zemp_N, zemp_sais, &
         zTstar, zcos_sais1, zconv, &
         zcos_sais2
      REAL(wp) ::           &
         zsumemp,           &          ! tampon used for the emp sum
         zsurf,             &          ! tampon used for the domain sum
         ztime,             &          ! time in hour
         ztimemax1, ztimemin1, &       ! 21th june,   and 21th december if date0 = 1st january
         ztimemax2, ztimemin2          ! 21th august, and 21th february if date0 = 1st january
      REAL(wp), DIMENSION(jpi,jpj) :: t_star
      INTEGER  ::   ji, jj             ! dummy loop indices

      INTEGER  ::           &
         zyear0,            &          ! initial year
         zmonth0,           &          ! initial month
         zday0,             &          ! initial day
         zday_year0                    ! initial day since january 1st

      NAMELIST/namflx/ q0, qsr0, emp0
      !!---------------------------------------------------------------------

      !same temperature, E-P as in HAZELEGER 2000

         ! Constant surface fluxes

         IF( kt == nit000 ) THEN
            IF(lwp) THEN
               WRITE(numout,*)' '
               WRITE(numout,*)' ocesbc  : Constant surface fluxes read in namelist'
               WRITE(numout,*)' ~~~~~~~ '
               WRITE(numout,*)'           Namelist namflx: set the constant flux values'
               WRITE(numout,*)'              net heat flux          q0   = ', q0  , ' W/m2'
               WRITE(numout,*)'              solar heat flux        qsr0 = ', qsr0, ' W/m2'
               WRITE(numout,*)'              net heat flux          emp0 = ', emp0, ' W/m2'
            ENDIF

            qt    (:,:) = q0
            qsr   (:,:) = qsr0
            emp   (:,:) = emp0
            emps  (:,:) = emp0
            qrp   (:,:) = 0.e0
            erp   (:,:) = 0.e0
   
            runoff(:,:) = 0.e0
         ENDIF


   END SUBROUTINE oce_sbc

# endif
#endif

#if defined key_dtasal
   !!----------------------------------------------------------------------
   !!   'key_dtasal'                                          salinity data
   !!----------------------------------------------------------------------
   SUBROUTINE oce_sbc_dmp
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc_dmp  ***
      !!                    
      !! ** Purpose : Computation of internal and evaporation damping terms 
      !!        for ocean surface boundary conditions 
      !!
      !! History :
      !!   9.0  !  04-01  (G. Madec, C. Ethe)  Original code
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj                   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj)  :: zsss, zfreeze
      REAL(wp) ::   zerp, zsrp
      CHARACTER (len=71) :: charout
#if ! defined key_dynspg_rl
      REAL(wp) ::   zwei
      REAL(wp) ::   zerpplus(jpi,jpj), zerpminus(jpi,jpj)
      REAL(wp) ::   zplus, zminus, zadefi
# if defined key_tradmp
      INTEGER jk
      REAL(wp), DIMENSION(jpi,jpj) ::   zstrdmp
# endif
#endif
      !!----------------------------------------------------------------------

#if defined key_ice_lim
      ! sea ice indicator (1 or 0)
      DO jj = 1, jpj
         DO ji = 1, jpi
            freezn(ji,jj) = MAX(0., SIGN(1., freeze(ji,jj)-rsmall) )
         END DO
      END DO
      zsss   (:,:) = sss_io(:,:)
      zfreeze(:,:) = freezn(:,:)
#else
      zsss   (:,:) = sb    (:,:,1)
      zfreeze(:,:) = freeze(:,:)
#endif

      ! Initialisation
      ! --------------
      ! Restoring coefficients on SST and SSS   
      zsrp = dqdt0 * ro0cpr * rauw   ! (Kg/m2/s) 

#if ! defined key_dynspg_rl 
      ! Free-surface
         
      ! Internal damping
# if defined key_tradmp
      ! Vertical mean of dampind trend (computed in tradmp module)
      zstrdmp(:,:) = 0.e0
      DO jk = 1, jpk
         zstrdmp(:,:) = zstrdmp(:,:) + strdmp(:,:,jk) * fse3t(:,:,jk)
      END DO
      ! volume flux associated to internal damping to climatology
      dmp(:,:) = zstrdmp(:,:) * rauw / ( zsss(:,:) + 1.e-20 )
# else
      dmp(:,:) = 0.e0            ! No internal damping
# endif
      
      !   evaporation damping term ( Surface restoring )
      zerpplus (:,:) = 0.e0
      zerpminus(:,:) = 0.e0
      zplus  =  15. / rday
      zminus = -15. / rday
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            zerp = ( 1. - 2.*upsrnfh(ji,jj) ) * zsrp   &
               & * ( zsss(ji,jj) - s_dta(ji,jj,1) )     &
               & / ( zsss(ji,jj) + 1.e-20        )
            
            zerp = MIN( zerp, zplus  )
            zerp = MAX( zerp, zminus )
            erp(ji,jj) = zerp
            zerpplus (ji,jj) = MAX( erp(ji,jj), 0.e0 )
            zerpminus(ji,jj) = MIN( erp(ji,jj), 0.e0 )
         END DO
      END DO

      aplus  = 0.e0
      aminus = 0.e0
      DO jj = 1, jpj
         DO ji = 1, jpi
            zwei   = e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            aplus  = aplus  + zerpplus (ji,jj) * zwei
            aminus = aminus - zerpminus(ji,jj) * zwei
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( aplus  )   ! sums over the global domain
      IF( lk_mpp )   CALL mpp_sum( aminus )

      IF(ln_ctl)   THEN
         WRITE(charout,FMT="('oce_sbc_dmp : a+ = ',D23.16, ' a- = ',D23.16)") aplus, aminus
         CALL prt_ctl_info(charout)
      ENDIF

      zadefi = MIN( aplus, aminus )
      IF( zadefi == 0.e0 ) THEN 
         erp(:,:) = 0.e0
      ELSE
         erp(:,:) = zadefi * ( zerpplus(:,:) / aplus + zerpminus(:,:) / aminus )
      ENDIF
      erp(:,:)=0 !byoung for 008
#else
      ! Rigid-lid (emp=emps=E-P-R+Erp)
      
      erp(:,:) = ( 1. - zfreeze(:,:) ) * zsrp    &   ! surface restoring term
         &     * ( zsss(:,:) - s_dta(:,:,1) )     &
         &     / ( zsss(:,:) + 1.e-20      )
#endif

   END SUBROUTINE oce_sbc_dmp

#else
   !!----------------------------------------------------------------------
   !!   Dummy routine                                      NO salinity data
   !!----------------------------------------------------------------------
   SUBROUTINE oce_sbc_dmp         ! Dummy routine
!      WRITE(*,*) 'oce_sbc_dmp: you should not have seen that print! error?'
   END SUBROUTINE oce_sbc_dmp
#endif

   !!======================================================================
END MODULE ocesbc
