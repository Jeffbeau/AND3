MODULE ldfeiv
   !!======================================================================
   !!                     ***  MODULE  ldfeiv  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
#if   defined key_traldf_eiv   &&   defined key_traldf_c2d
   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'      and                     eddy induced velocity
   !!   'key_traldf_c2d'                    2D tracer lateral  mixing coef.
   !!----------------------------------------------------------------------
   !!   ldf_eiv      : compute the eddy induced velocity coefficients
   !!                  Same results but not same routine if 'key_autotasking'
   !!                  is defined or not
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! iso-neutral slopes
   USE flxrnf          ! 
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE
   
   !! * Routine accessibility
   PUBLIC ldf_eiv               ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfeiv.F90,v 1.8 2006/04/28 12:24:20 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------

CONTAINS

# if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking' :                            autotasking (j-slab)
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_eiv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!      growth rate of baroclinic instability.
      !!
      !! ** Method :
      !!
      !! ** Action :   uslp(),   : i- and j-slopes of neutral surfaces
      !!               vslp()      at u- and v-points, resp.
      !!               wslpi(),  : i- and j-slopes of neutral surfaces
      !!               wslpj()     at w-points. 
      !!
      !! History :
      !!   8.1  !  99-03  (G. Madec, A. Jouzeau)  Original code
      !!   8.5  !  02-06  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   &
         zfw, ze3w, zn2, zf20,       &  ! temporary scalars
         zaht, zaht_min
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zn, zah, zhw, zross            ! workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) '~~~~~~~   key_autotasking'
      ENDIF
      
      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         
         ! 0. Local initialization
         ! -----------------------
         zn   (:,jj) = 0.e0
         zhw  (:,jj) = 5.e0
         zah  (:,jj) = 0.e0
         zross(:,jj) = 0.e0
         
         ! 1. Compute lateral diffusive coefficient 
         ! ----------------------------------------

!CDIR NOVERRCHK 
         DO jk = 1, jpk
!CDIR NOVERRCHK 
            DO ji = 2, jpim1
               ! Take the max of N^2 and zero then take the vertical sum 
               ! of the square root of the resulting N^2 ( required to compute 
               ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
               zn2 = MAX( rn2(ji,jj,jk), 0.e0 )
               ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
               zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
               ! Compute elements required for the inverse time scale of baroclinic
               ! eddies using the isopycnal slopes calculated in ldfslp.F : 
               ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
               zah(ji,jj) = zah(ji,jj) + zn2   &
                              * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)    &
                                + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) )   &
                              * ze3w
               zhw(ji,jj) = zhw(ji,jj) + ze3w
            END DO 
         END DO 
 
!CDIR NOVERRCHK 
         DO ji = 2, jpim1
            zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
            ! Rossby radius at w-point taken < 40km and  > 2km
            zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
            ! Compute aeiw by multiplying Ro^2 and T^-1
            aeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
            ! Take the minimum between aeiw and 1000m^2/s for depth levels
            ! lower than 20 (21 in w- point)
            IF( mbathy(ji,jj) <= 21. ) aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
         END DO

         ! Decrease the coefficient in the tropics (20N-20S) 
         zf20 = 2. * omega * sin( rad * 20. )
         DO ji = 2, jpim1
            aeiw(ji,jj) = MIN( 1., ABS( ff(ji,jj) / zf20 ) ) * aeiw(ji,jj)
         END DO
  
         ! ORCA R05: Take the minimum between aeiw  and 1000m2/s
         IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05
            DO ji = 2, jpim1
               aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
            END DO
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! lateral boundary condition on aeiw 
      CALL lbc_lnk( aeiw, 'W', 1. )

      ! Average the diffusive coefficient at u- v- points 
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji+1,jj  ))
            aeiv(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji  ,jj+1))
         END DO 
      END DO 
      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! lateral boundary condition on aeiu, aeiv 
      CALL lbc_lnk( aeiu, 'U', 1. )
      CALL lbc_lnk( aeiv, 'V', 1. )

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab2d_1=aeiu, clinfo1=' eiv  - u: ', ovlap=1)
         CALL prt_ctl(tab2d_1=aeiv, clinfo1=' eiv  - v: ', ovlap=1)
      ENDIF
      
      ! ORCA R05: add a space variation on aht (=aeiv except at the equator and river mouth)
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
         zf20     = 2. * omega * SIN( rad * 20. )
         zaht_min = 100.                              ! minimum value for aht
         DO jj = 1, jpj
            DO ji = 1, jpi
               zaht      = ( 1. -  MIN( 1., ABS( ff(ji,jj) / zf20 ) ) ) * ( aht0 - zaht_min )  &
                  &      + aht0 * upsrnfh(ji,jj)                          ! enhanced near river mouths
               ahtu(ji,jj) = MAX( MAX( zaht_min, aeiu(ji,jj) ) + zaht, aht0 )
               ahtv(ji,jj) = MAX( MAX( zaht_min, aeiv(ji,jj) ) + zaht, aht0 )
               ahtw(ji,jj) = MAX( MAX( zaht_min, aeiw(ji,jj) ) + zaht, aht0 )
            END DO
         END DO
         IF(ln_ctl) THEN
            CALL prt_ctl(tab2d_1=ahtu, clinfo1=' aht  - u: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtv, clinfo1=' aht  - v: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtw, clinfo1=' aht  - w: ', ovlap=1)
         ENDIF
      ENDIF

      IF( aeiv0 == 0.e0 ) THEN
         aeiu(:,:) = 0.e0
         aeiv(:,:) = 0.e0
         aeiw(:,:) = 0.e0
      ENDIF

   END SUBROUTINE ldf_eiv

# else
   !!----------------------------------------------------------------------
   !!   Default key                                             k-j-i loops
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_eiv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!      growth rate of baroclinic instability.
      !!
      !! ** Method :
      !!
      !! ** Action : - uslp(),  : i- and j-slopes of neutral surfaces
      !!             - vslp()      at u- and v-points, resp.
      !!             - wslpi(),  : i- and j-slopes of neutral surfaces
      !!             - wslpj()     at w-points. 
      !!
      !! History :
      !!   8.1  !  99-03  (G. Madec, A. Jouzeau)  Original code
      !!   8.5  !  02-06  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   &
         zfw, ze3w, zn2, zf20,       &  ! temporary scalars
         zaht, zaht_min
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zn, zah, zhw, zross            ! workspace
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF
      
      ! 0. Local initialization
      ! -----------------------
      zn   (:,:) = 0.e0
      zhw  (:,:) = 5.e0
      zah  (:,:) = 0.e0
      zross(:,:) = 0.e0


      ! 1. Compute lateral diffusive coefficient 
      ! ----------------------------------------

      DO jk = 1, jpk
#  if defined key_vectopt_loop  &&  ! defined key_autotasking
!CDIR NOVERRCHK 
         DO ji = 1, jpij   ! vector opt.
            ! Take the max of N^2 and zero then take the vertical sum
            ! of the square root of the resulting N^2 ( required to compute
            ! internal Rossby radius Ro = .5 * sum_jpk(N) / f
            zn2 = MAX( rn2(ji,1,jk), 0.e0 )
            zn(ji,1) = zn(ji,1) + SQRT( zn2 ) * fse3w(ji,1,jk)
            ! Compute elements required for the inverse time scale of baroclinic
            ! eddies using the isopycnal slopes calculated in ldfslp.F :
            ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
            ze3w = fse3w(ji,1,jk) * tmask(ji,1,jk)
               zah(ji,1) = zah(ji,1) + zn2   &
                              * ( wslpi(ji,1,jk) * wslpi(ji,1,jk)    &
                                + wslpj(ji,1,jk) * wslpj(ji,1,jk) )   &
                              * ze3w
            zhw(ji,1) = zhw(ji,1) + ze3w
         END DO
#  else
         DO jj = 2, jpjm1
!CDIR NOVERRCHK 
            DO ji = 2, jpim1
               ! Take the max of N^2 and zero then take the vertical sum 
               ! of the square root of the resulting N^2 ( required to compute 
               ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
               zn2 = MAX( rn2(ji,jj,jk), 0.e0 )
               zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
               ! Compute elements required for the inverse time scale of baroclinic
               ! eddies using the isopycnal slopes calculated in ldfslp.F : 
               ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
               ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
               zah(ji,jj) = zah(ji,jj) + zn2   &
                              * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)    &
                                + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) )  &
                              * ze3w
               zhw(ji,jj) = zhw(ji,jj) + ze3w
            END DO 
         END DO 
#  endif
      END DO 

      DO jj = 2, jpjm1
!CDIR NOVERRCHK 
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
            ! Rossby radius at w-point taken < 40km and  > 2km
            zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
            ! Compute aeiw by multiplying Ro^2 and T^-1
            aeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
            ! Take the minimum between aeiw and 1000m^2/s for depth levels
            ! lower than 20 (21 in w- point)
            IF( mbathy(ji,jj) <= 21. ) aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
         END DO
      END DO

      ! Decrease the coefficient in the tropics (20N-20S) 
         zf20 = 2. * omega * sin( rad * 20. )
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiw(ji,jj) = MIN( 1., ABS( ff(ji,jj) / zf20 ) ) * aeiw(ji,jj)
         END DO
      END DO

      ! ORCA R05: Take the minimum between aeiw  and 1000m2/s
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               aeiw(ji,jj) = MIN( aeiw(ji,jj), aeiv0 )
            END DO
         END DO
      ENDIF

      ! lateral boundary condition on aeiw 
      CALL lbc_lnk( aeiw, 'W', 1. )

      ! Average the diffusive coefficient at u- v- points 
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji+1,jj  ) )
            aeiv(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji  ,jj+1) )
         END DO 
      END DO 

      ! lateral boundary condition on aeiu, aeiv
      CALL lbc_lnk( aeiu, 'U', 1. )
      CALL lbc_lnk( aeiv, 'V', 1. )

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab2d_1=aeiu, clinfo1=' eiv  - u: ', ovlap=1)
         CALL prt_ctl(tab2d_1=aeiv, clinfo1=' eiv  - v: ', ovlap=1)
      ENDIF

      ! ORCA R05: add a space variation on aht (=aeiv except at the equator and river mouth)
      IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
         zf20     = 2. * omega * SIN( rad * 20. )
         zaht_min = 100.                              ! minimum value for aht
         DO jj = 1, jpj
            DO ji = 1, jpi
               zaht      = ( 1. -  MIN( 1., ABS( ff(ji,jj) / zf20 ) ) ) * ( aht0 - zaht_min )  &
                  &      + aht0 * upsrnfh(ji,jj)                          ! enhanced near river mouths
               ahtu(ji,jj) = MAX( MAX( zaht_min, aeiu(ji,jj) ) + zaht, aht0 )
               ahtv(ji,jj) = MAX( MAX( zaht_min, aeiv(ji,jj) ) + zaht, aht0 )
               ahtw(ji,jj) = MAX( MAX( zaht_min, aeiw(ji,jj) ) + zaht, aht0 )
            END DO
         END DO
         IF(ln_ctl) THEN
            CALL prt_ctl(tab2d_1=ahtu, clinfo1=' aht  - u: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtv, clinfo1=' aht  - v: ', ovlap=1)
            CALL prt_ctl(tab2d_1=ahtw, clinfo1=' aht  - w: ', ovlap=1)
         ENDIF
      ENDIF
      
      IF( aeiv0 == 0.e0 ) THEN
         aeiu(:,:) = 0.e0
         aeiv(:,:) = 0.e0
         aeiw(:,:) = 0.e0
      ENDIF

   END SUBROUTINE ldf_eiv

# endif

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE ldf_eiv( kt )       ! Empty routine
      if(lwp)WRITE(numout,*) 'ldf_eiv: You should not have seen this print! error?', kt
   END SUBROUTINE ldf_eiv
#endif

   !!======================================================================
END MODULE ldfeiv
