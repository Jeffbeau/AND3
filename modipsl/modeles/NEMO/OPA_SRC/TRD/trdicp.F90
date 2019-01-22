MODULE trdicp
   !!======================================================================
   !!                       ***  MODULE  trdicp  ***
   !! Ocean diagnostics:  ocean tracers and dynamic trends
   !!=====================================================================
#if  defined key_trdtra   ||   defined key_trddyn   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_trdtra'  or                  active tracers trends diagnostics
   !!   'key_trddyn'                            momentum trends diagnostics
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd              : verify the basin averaged properties for tra/dyn 
   !!   trd_dwr          : print dynmaic trends in ocean.output file
   !!   trd_twr          : print tracers trends in ocean.output file
   !!   trd_icp_init     : initialization step
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE trdmod_oce      ! ocean variables trends
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distibuted memory computing library
   USE eosbn2          ! equation of state
   USE phycst          ! physical constants

   IMPLICIT NONE
   PRIVATE

   !! * Interfaces
   INTERFACE trd
      MODULE PROCEDURE trd_2d, trd_3d
   END INTERFACE

   !! * Routine accessibility
   PUBLIC trd                   ! called by step.F90
   PUBLIC trd_dwr               ! called by step.F90
   PUBLIC trd_twr               ! called by step.F90
   PUBLIC trd_icp_init          ! called by opa.F90

   !! * Shared module variables
#if  defined key_trdtra   &&   defined key_trddyn    ||   defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .TRUE.    !: tracers  trend flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .TRUE.    !: momentum trend flag
#elif  defined key_trdtra
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .TRUE.    !: tracers  trend flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .FALSE.   !: momentum trend flag
#elif  defined key_trddyn
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .FALSE.   !: tracers  trend flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .TRUE.    !: momentum trend flag
#endif

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdicp.F90,v 1.3 2005/03/30 10:29:20 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trd_2d(ptrd2dx, ptrd2dy, ktrd , ctype)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_2d  ***
      !! 
      !! ** Purpose : verify the basin averaged properties of tracers and/or
      !!              momentum equations at every time step frequency ntrd.
      !!
      !! ** Method :
      !!
      !! History :
      !!        !  91-12 (G. Madec)
      !!        !  92-06 (M. Imbard) add time step frequency
      !!        !  96-01 (G. Madec)  terrain following coordinates
      !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08 (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         ptrd2dx,                      &   ! Temperature or U trend 
         ptrd2dy                           ! Salinity    or V trend

      INTEGER, INTENT( in ) ::   ktrd      ! tracer trend index

      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
         !                                 ! 'DYN' or 'TRA'

      !! * Local declarations
      INTEGER ::   ji, jj        ! loop indices
      REAL(wp) ::   &
         zbt, zbtu, zbtv,     &  ! temporary scalars
         zmsku, zmskv            !    "         "
      !!----------------------------------------------------------------------

      ! 1. Advective trends and forcing trend
      ! -------------------------------------

      ! 1.1 Mask the forcing trend and substract it from the vertical diffusion trend
      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,1)
               zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,1)
               ptrd2dx(ji,jj) = ptrd2dx(ji,jj) * zmsku
               ptrd2dy(ji,jj) = ptrd2dy(ji,jj) * zmskv
            END DO
         END DO
         ptrd2dx(jpi, : ) = 0.e0      ;      ptrd2dy(jpi, : ) = 0.e0
         ptrd2dx( : ,jpj) = 0.e0      ;      ptrd2dy( : ,jpj) = 0.e0

      CASE ('TRA')              ! Tracers
         ptrd2dx(:,:) = ptrd2dx(:,:) * tmask_i(:,:)
         ptrd2dy(:,:) = ptrd2dy(:,:) * tmask_i(:,:)

      END SELECT
      
      ! 2. Basin averaged tracer trend
      ! ------------------------------

      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum
         umo(ktrd) = 0.e0
         vmo(ktrd) = 0.e0

         SELECT CASE (ktrd)

         CASE (jpdtdswf)         ! surface forcing
            DO jj = 1, jpj
               DO ji = 1, jpi
                  umo(ktrd) = umo(ktrd) + ptrd2dx(ji,jj) * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,1)
                  vmo(ktrd) = vmo(ktrd) + ptrd2dy(ji,jj) * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,1)
               END DO
            END DO

         CASE (jpdtdbfr)         ! bottom friction fluxes
            DO jj = 1, jpj
               DO ji = 1, jpi
                  umo(ktrd) = umo(ktrd) + ptrd2dx(ji,jj)
                  vmo(ktrd) = vmo(ktrd) + ptrd2dy(ji,jj)
               END DO
            END DO

         END SELECT

      CASE ('TRA')              ! Tracers
         tmo(ktrd) = 0.e0
         smo(ktrd) = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1)
               tmo(ktrd) =  tmo(ktrd) + ptrd2dx(ji,jj) * zbt
               smo(ktrd) =  smo(ktrd) + ptrd2dy(ji,jj) * zbt
            END DO
         END DO

      END SELECT
      
      ! 3. Basin averaged tracer square trend
      ! -------------------------------------
      ! c a u t i o n: field now
      
      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum
         hke(ktrd) = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,1)
               zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,1)
               hke(ktrd) = hke(ktrd)   &
               &   + un(ji,jj,1) * ptrd2dx(ji,jj) * zbtu &
               &   + vn(ji,jj,1) * ptrd2dy(ji,jj) * zbtv
            END DO
         END DO

      CASE ('TRA')              ! Tracers
         t2(ktrd) = 0.e0
         s2(ktrd) = 0.e0
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1)
               t2(ktrd) = t2(ktrd) + ptrd2dx(ji,jj) * zbt * tn(ji,jj,1)
               s2(ktrd) = s2(ktrd) + ptrd2dy(ji,jj) * zbt * sn(ji,jj,1)
            END DO
         END DO
      
      END SELECT

   END SUBROUTINE trd_2d



   SUBROUTINE trd_3d(ptrd3dx, ptrd3dy, ktrd, ctype)
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_3d  ***
      !! 
      !! ** Purpose : verify the basin averaged properties of tracers and/or
      !!              momentum equations at every time step frequency ntrd.
      !!
      !! ** Method :
      !!
      !! History :
      !!        !  91-12 (G. Madec)
      !!        !  92-06 (M. Imbard) add time step frequency
      !!        !  96-01 (G. Madec)  terrain following coordinates
      !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08 (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
          ptrd3dx,                     &   ! Temperature or U trend 
          ptrd3dy                          ! Salinity    or V trend

      INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index

      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
         !                                 ! 'DYN' or 'TRA'

      !! * Local declarations
      INTEGER ::   ji, jj, jk
      REAL(wp) ::   &
         zbt, zbtu, zbtv,               &  ! temporary scalars
         zmsku, zmskv
      !!----------------------------------------------------------------------

      ! 1. Advective trends and forcing trend
      ! -------------------------------------

      ! Mask the trends
      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum        
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,jk)
                  zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)
                  ptrd3dx(ji,jj,jk) = ptrd3dx(ji,jj,jk) * zmsku
                  ptrd3dy(ji,jj,jk) = ptrd3dy(ji,jj,jk) * zmskv
               ENDDO
            ENDDO
         ENDDO

         ptrd3dx(jpi, : ,:) = 0.e0      ;      ptrd3dy(jpi, : ,:) = 0.e0
         ptrd3dx( : ,jpj,:) = 0.e0      ;      ptrd3dy( : ,jpj,:) = 0.e0

      CASE ('TRA')              ! Tracers
         DO jk = 1, jpk
            ptrd3dx(:,:,jk) = ptrd3dx(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
            ptrd3dy(:,:,jk) = ptrd3dy(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
         ENDDO

      END SELECT   

      ! 2. Basin averaged tracer/momentum trend
      ! ---------------------------------------
      
      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum
         umo(ktrd) = 0.e0
         vmo(ktrd) = 0.e0
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
                  zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
                  umo(ktrd) = umo(ktrd) + ptrd3dx(ji,jj,jk) * zbtu
                  vmo(ktrd) = vmo(ktrd) + ptrd3dy(ji,jj,jk) * zbtv
               END DO
            END DO
         END DO

      CASE ('TRA')              ! Tracers
         tmo(ktrd) = 0.e0
         smo(ktrd) = 0.e0
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) 
                  tmo(ktrd) =  tmo(ktrd) + ptrd3dx(ji,jj,jk) * zbt
                  smo(ktrd) =  smo(ktrd) + ptrd3dy(ji,jj,jk) * zbt
               END DO
            END DO
         END DO

      END SELECT

      ! 3. Basin averaged tracer/momentum square trend
      ! ----------------------------------------------
      ! c a u t i o n: field now
      
      SELECT CASE (ctype)

      CASE ('DYN')              ! Momentum
         hke(ktrd) = 0.e0
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
                  zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
                  hke(ktrd) = hke(ktrd)   &
                  &   + un(ji,jj,jk) * ptrd3dx(ji,jj,jk) * zbtu &
                  &   + vn(ji,jj,jk) * ptrd3dy(ji,jj,jk) * zbtv
               END DO
            END DO
         END DO

      CASE ('TRA')              ! Tracers
         t2(ktrd) = 0.e0
         s2(ktrd) = 0.e0
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
                  t2(ktrd) = t2(ktrd) + ptrd3dx(ji,jj,jk) * zbt * tn(ji,jj,jk)
                  s2(ktrd) = s2(ktrd) + ptrd3dy(ji,jj,jk) * zbt * sn(ji,jj,jk)
               END DO
            END DO
         END DO

      END SELECT

   END SUBROUTINE trd_3d



   SUBROUTINE trd_icp_init
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_icp_init  ***
      !! 
      !! ** Purpose :   
      !!
      !! ** Method  :
      !!
      !! History :
      !!   9.0  !  03-09 (G. Madec)  Original code
      !!        !  04-08 (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER :: ji, jj, jk

      REAL(wp) ::   zmskt
#if  defined key_trddyn
      REAL(wp) ::   zmsku,zmskv
#endif

      NAMELIST/namtrd/ ntrd, nctls
      !!----------------------------------------------------------------------

      ! namelist namtrd : trend diagnostic
      REWIND( numnam )
      READ  ( numnam, namtrd )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trd_icp_init : integral constraints properties trends'
         WRITE(numout,*) '~~~~~~~~~~~~~'
         WRITE(numout,*) ' '
         WRITE(numout,*) '          Namelist namtrd : '
         WRITE(numout,*) '             time step frequency trend       ntrd  = ', ntrd
      ENDIF

      ! initialisation of BBL tracers lateral diffusion to zero
      tldfbbl(:,:) = 0.e0   ;   sldfbbl(:,:) = 0.e0  
      ! initialisation of BBL tracers lateral advection to zero
      tladbbl(:,:) = 0.e0   ;   sladbbl(:,:) = 0.e0  
      ! initialisation of workspace
      tladi(:,:,:) = 0.e0  ;  tladj(:,:,:) = 0.e0
      sladi(:,:,:) = 0.e0  ;  sladj(:,:,:) = 0.e0

      ! Total volume at t-points:
      tvolt = 0.e0
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zmskt = tmask(ji,jj,jk) * tmask_i(ji,jj)
               tvolt = tvolt + zmskt * e1t(ji,jj) *e2t(ji,jj) * fse3t(ji,jj,jk)
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( tvolt )   ! sum over the global domain

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          total ocean volume at T-point   tvolt = ',tvolt
      ENDIF

#if  defined key_trddyn
      ! Initialization of potential to kinetic energy conversion
      rpktrd = 0.e0

      ! Total volume at u-, v- points:
      tvolu = 0.e0
      tvolv = 0.e0

      DO jk = 1, jpk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,jk)
               zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)
               tvolu = tvolu + zmsku * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
               tvolv = tvolv + zmskv * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( tvolu )   ! sums over the global domain
      IF( lk_mpp )   CALL mpp_sum( tvolv )

      IF(lwp) THEN
         WRITE(numout,*) '          total ocean volume at U-point   tvolu = ',tvolu
         WRITE(numout,*) '          total ocean volume at V-point   tvolv = ',tvolv
         WRITE(numout,*) ' '
      ENDIF
#endif

   END SUBROUTINE trd_icp_init



   SUBROUTINE trd_dwr( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dwr  ***
      !! 
      !! ** Purpose :  write dynamic trends in ocean.output 
      !!
      !! ** Method  :
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!        !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER ::   ji, jj, jk
      REAL(wp) ::   &
         ze1e2w,zcof,        &  !    "         "
         zbe1ru, zbe2rv,     &  !    "         "
         zbtr, ztz, zth 

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
               zkepe, zkx, zky, zkz              ! temporary arrays
      !!----------------------------------------------------------------------

      ! I. Momentum trends
      ! -------------------

      IF( MOD(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN

         ! I.1 Conversion potential energy - kinetic energy
         ! --------------------------------------------------
         ! c a u t i o n here, trends are computed at kt+1 (now , but after the swap)

         zkx(:,:,:) = 0.e0
         zky(:,:,:) = 0.e0
         zkz(:,:,:) = 0.e0
         zkepe(:,:,:) = 0.e0
   
         CALL eos( tn, sn, rhd, rhop )       ! now potential and in situ densities

         ! 4.1 Density flux at w-point
         DO jk = 2, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ze1e2w = 0.5 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk) * tmask_i(ji,jj)
                  zkz(ji,jj,jk) = ze1e2w / rau0 * ( rhop(ji,jj,jk) + rhop(ji,jj,jk-1) )
               END DO
            END DO
         END DO
         zkz  (:,:, 1 ) = 0.e0
         
         ! Density flux at u and v-points
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zcof   = 0.5 / rau0
                  zbe1ru = zcof * e2u(ji,jj) * fse3u(ji,jj,jk) * un(ji,jj,jk)
                  zbe2rv = zcof * e1v(ji,jj) * fse3v(ji,jj,jk) * vn(ji,jj,jk)
                  zkx(ji,jj,jk) = zbe1ru * ( rhop(ji,jj,jk) + rhop(ji+1,jj,jk) )
                  zky(ji,jj,jk) = zbe2rv * ( rhop(ji,jj,jk) + rhop(ji,jj+1,jk) )
               END DO
            END DO
         END DO
         
         ! Density flux divergence at t-point
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  ztz = - zbtr * (    zkz(ji,jj,jk) - zkz(ji,jj,jk+1) )
                  zth = - zbtr * (  ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk) )   &
                    &             + ( zky(ji,jj,jk) - zky(ji,jj-1,jk) )  )
                  zkepe(ji,jj,jk) = (zth + ztz) * tmask(ji,jj,jk) * tmask_i(ji,jj)
               END DO
            END DO
         END DO
         zkepe( : , : ,jpk) = 0.e0
         zkepe( : ,jpj, : ) = 0.e0
         zkepe(jpi, : , : ) = 0.e0

         ! I.2 Basin averaged kinetic energy trend
         ! ----------------------------------------
         peke = 0.e0
         DO jk = 1,jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  peke = peke + zkepe(ji,jj,jk) * grav * fsdept(ji,jj,jk)   &
                     &                     * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
               END DO
            END DO
         END DO

         ! I.3 Sums over the global domain
         ! ---------------------------------
         IF( lk_mpp ) THEN
               CALL mpp_sum( peke )
               CALL mpp_sum( umo , 11 )
               CALL mpp_sum( vmo , 11 )
               CALL mpp_sum( hke , 10 )
         END IF

         ! I.2 Print dynamic trends in the ocean.output file
         ! --------------------------------------------------

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9500) kt
            WRITE (numout,9501) umo( 1) / tvolu, vmo( 1) / tvolv
            WRITE (numout,9502) umo( 2) / tvolu, vmo( 2) / tvolv
            WRITE (numout,9503) umo( 3) / tvolu, vmo( 3) / tvolv
            WRITE (numout,9504) umo( 4) / tvolu, vmo( 4) / tvolv
            WRITE (numout,9505) umo( 5) / tvolu, vmo( 5) / tvolv
            WRITE (numout,9506) umo( 6) / tvolu, vmo( 6) / tvolv
            WRITE (numout,9507) umo( 7) / tvolu, vmo( 7) / tvolv
            WRITE (numout,9508) umo( 8) / tvolu, vmo( 8) / tvolv
            WRITE (numout,9509) umo(10) / tvolu, vmo(10) / tvolv
            WRITE (numout,9510) umo( 9) / tvolu, vmo( 9) / tvolv
            WRITE (numout,9511) umo(11) / tvolu, vmo(11) / tvolv
            WRITE (numout,9512)
            WRITE (numout,9513)                                                 &
            &     (  umo(1) + umo(2) + umo(3) + umo( 4) + umo( 5) + umo(6)    &
            &      + umo(7) + umo(8) + umo(9) + umo(10) + umo(11) ) / tvolu,   &
            &     (  vmo(1) + vmo(2) + vmo(3) + vmo( 4) + vmo( 5) + vmo(6)    &
            &      + vmo(7) + vmo(8) + vmo(9) + vmo(10) + vmo(11) ) / tvolv
         ENDIF

 9500    FORMAT(' momentum trend at it= ', i6, ' :', /' ==============================')
 9501    FORMAT(' pressure gradient          u= ', e20.13, '    v= ', e20.13)
 9502    FORMAT(' ke gradient                u= ', e20.13, '    v= ', e20.13)
 9503    FORMAT(' relative vorticity term    u= ', e20.13, '    v= ', e20.13)
 9504    FORMAT(' coriolis term              u= ', e20.13, '    v= ', e20.13)
 9505    FORMAT(' horizontal diffusion       u= ', e20.13, '    v= ', e20.13)
 9506    FORMAT(' vertical advection         u= ', e20.13, '    v= ', e20.13)
 9507    FORMAT(' vertical diffusion         u= ', e20.13, '    v= ', e20.13)
 9508    FORMAT(' surface pressure gradient  u= ', e20.13, '    v= ', e20.13)
 9509    FORMAT(' forcing term               u= ', e20.13, '    v= ', e20.13)
 9510    FORMAT(' dampimg term               u= ', e20.13, '    v= ', e20.13)
 9511    FORMAT(' bottom flux                u= ', e20.13, '    v= ', e20.13)
 9512    FORMAT(' -----------------------------------------------------------------------------')
 9513    FORMAT(' total trend                u= ', e20.13, '    v= ', e20.13)

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9520) kt
            WRITE (numout,9521) hke( 1) / tvolt
            WRITE (numout,9522) hke( 2) / tvolt
            WRITE (numout,9523) hke( 3) / tvolt
            WRITE (numout,9524) hke( 4) / tvolt
            WRITE (numout,9525) hke( 5) / tvolt
            WRITE (numout,9526) hke( 6) / tvolt
            WRITE (numout,9527) hke( 7) / tvolt
            WRITE (numout,9528) hke( 8) / tvolt
            WRITE (numout,9529) hke(10) / tvolt
            WRITE (numout,9530) hke( 9) / tvolt
            WRITE (numout,9531)
            WRITE (numout,9532)   &
            &     (  hke(1) + hke(2) + hke(3) + hke(4) + hke(5) + hke(6)   &
            &      + hke(7) + hke(8) + hke(9) + hke(10) ) / tvolt
         ENDIF

 9520    FORMAT(' kinetic energy trend at it= ', i6, ' :', /' ====================================')
 9521    FORMAT(' pressure gradient         u2= ', e20.13)
 9522    FORMAT(' ke gradient               u2= ', e20.13)
 9523    FORMAT(' relative vorticity term   u2= ', e20.13)
 9524    FORMAT(' coriolis term             u2= ', e20.13)
 9525    FORMAT(' horizontal diffusion      u2= ', e20.13)
 9526    FORMAT(' vertical advection        u2= ', e20.13)
 9527    FORMAT(' vertical diffusion        u2= ', e20.13)
 9528    FORMAT(' surface pressure gradient u2= ', e20.13)
 9529    FORMAT(' forcing term              u2= ', e20.13)
 9530    FORMAT(' dampimg term              u2= ', e20.13)
 9531    FORMAT(' --------------------------------------------------')
 9532    FORMAT(' total trend               u2= ', e20.13)

         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9540) kt
            WRITE (numout,9541) ( hke(2) + hke(3) + hke(6) ) / tvolt
            WRITE (numout,9542) ( hke(2) + hke(6) ) / tvolt
            WRITE (numout,9543) ( hke(4) ) / tvolt
            WRITE (numout,9544) ( hke(3) ) / tvolt
            WRITE (numout,9545) ( hke(8) ) / tvolt
            WRITE (numout,9546) ( hke(5) ) / tvolt
            WRITE (numout,9547) ( hke(7) ) / tvolt
            WRITE (numout,9548) ( hke(1) ) / tvolt, rpktrd / tvolt
         ENDIF

 9540    FORMAT(' energetic consistency at it= ', i6, ' :', /' =========================================')
 9541    FORMAT(' 0 = non linear term(true if key_vorenergy or key_combined): ', e20.13)
 9542    FORMAT(' 0 = ke gradient + vertical advection              : ', e20.13)
 9543    FORMAT(' 0 = coriolis term  (true if key_vorenergy or key_combined): ', e20.13)
 9544    FORMAT(' 0 = uh.( rot(u) x uh ) (true if enstrophy conser.)    : ', e20.13)
 9545    FORMAT(' 0 = surface pressure gradient                     : ', e20.13)
 9546    FORMAT(' 0 > horizontal diffusion                          : ', e20.13)
 9547    FORMAT(' 0 > vertical diffusion                            : ', e20.13)
 9548    FORMAT(' pressure gradient u2 = - 1/rau0 u.dz(rhop)        : ', e20.13, '  u.dz(rhop) =', e20.13)

         ! Save potential to kinetic energy conversion for next time step
         rpktrd = peke

      ENDIF

   END SUBROUTINE trd_dwr




   SUBROUTINE trd_twr( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_twr  ***
      !! 
      !! ** Purpose :  write active tracers trends in ocean.output 
      !!
      !! ** Method  :
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!        !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !!----------------------------------------------------------------------

      ! I. Tracers trends
      ! -----------------

      IF( MOD(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN

         ! I.1 Sums over the global domain
         ! -------------------------------
         IF( lk_mpp ) THEN
            CALL mpp_sum( tmo, 10 )   
            CALL mpp_sum( smo, 10 )
            CALL mpp_sum( t2 , 10 )
            CALL mpp_sum( s2 , 10 )
         ENDIF

         ! I.2 Print tracers trends in the ocean.output file
         ! --------------------------------------------------
         
         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9400) kt
            WRITE (numout,9401) tmo(1) / tvolt, smo(1) / tvolt
            WRITE (numout,9402) tmo(2) / tvolt, smo(2) / tvolt
            WRITE (numout,9403) tmo(3) / tvolt, smo(3) / tvolt
            WRITE (numout,9404) tmo(4) / tvolt, smo(4) / tvolt
            WRITE (numout,9405) tmo(5) / tvolt, smo(5) / tvolt
            WRITE (numout,9406) tmo(6) / tvolt, smo(6) / tvolt
            WRITE (numout,9407) tmo(7) / tvolt
            WRITE (numout,9408) tmo(8) / tvolt, smo(8) / tvolt
            WRITE (numout,9409)
            WRITE (numout,9410) (  tmo(1) + tmo(2) + tmo(3) + tmo(4)              &
            &                    + tmo(5) + tmo(6) + tmo(7) + tmo(8) ) / tvolt,   &
            &                   (  smo(1) + smo(2) + smo(3) + smo(4)              &
            &                    + smo(5) + smo(6)           + smo(8) ) / tvolt
         ENDIF

9400     FORMAT(' tracer trend at it= ',i6,' :     temperature',   &
              '              salinity',/' ============================')
9401     FORMAT(' horizontal advection        ',e20.13,'     ',e20.13)
9402     FORMAT(' vertical advection          ',e20.13,'     ',e20.13)
9403     FORMAT(' horizontal diffusion        ',e20.13,'     ',e20.13)
9404     FORMAT(' vertical diffusion          ',e20.13,'     ',e20.13)
9405     FORMAT(' STATIC instability mixing   ',e20.13,'     ',e20.13)
9406     FORMAT(' damping term                ',e20.13,'     ',e20.13)
9407     FORMAT(' penetrative qsr             ',e20.13,'     ',e20.13)
9408     FORMAT(' forcing term                ',e20.13,'     ',e20.13)
9409     FORMAT(' -------------------------------------------------------------------------')
9410     FORMAT(' total trend                 ',e20.13,'     ',e20.13)


         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9420) kt
            WRITE (numout,9421) t2(1) / tvolt, s2(1) / tvolt
            WRITE (numout,9422) t2(2) / tvolt, s2(2) / tvolt
            WRITE (numout,9423) t2(3) / tvolt, s2(3) / tvolt
            WRITE (numout,9424) t2(4) / tvolt, s2(4) / tvolt
            WRITE (numout,9425) t2(5) / tvolt, s2(5) / tvolt
            WRITE (numout,9426) t2(6) / tvolt, s2(6) / tvolt
            WRITE (numout,9427) t2(7) / tvolt
            WRITE (numout,9428) t2(8) / tvolt, s2(8) / tvolt
            WRITE (numout,9429)
            WRITE (numout,9430) (  t2(1) + t2(2) + t2(3) + t2(4)              &
            &                    + t2(5) + t2(6) + t2(7) + t2(8) ) / tvolt,   &
            &                   (  s2(1) + s2(2) + s2(3) + s2(4)              &
            &                    + s2(5) + s2(6)          + s2(8) ) / tvolt
         ENDIF

9420     FORMAT(' tracer**2 trend at it= ', i6, ' :      temperature',   &
            '               salinity', /, ' ===============================')
9421     FORMAT(' horizontal advection      * t   ', e20.13, '     ', e20.13)
9422     FORMAT(' vertical advection        * t   ', e20.13, '     ', e20.13)
9423     FORMAT(' horizontal diffusion      * t   ', e20.13, '     ', e20.13)
9424     FORMAT(' vertical diffusion        * t   ', e20.13, '     ', e20.13)
9425     FORMAT(' STATIC instability mixing * t   ', e20.13, '     ', e20.13)
9426     FORMAT(' damping term              * t   ', e20.13, '     ', e20.13)
9427     FORMAT(' penetrative qsr           * t   ', e20.13, '     ', e20.13)
9428     FORMAT(' forcing term              * t   ', e20.13, '     ', e20.13)
9429     FORMAT(' -----------------------------------------------------------------------------')
9430     FORMAT(' total trend                *t = ', e20.13, '  *s = ', e20.13)


         IF(lwp) THEN
            WRITE (numout,*)
            WRITE (numout,*)
            WRITE (numout,9440) kt
            WRITE (numout,9441) ( tmo(1)+tmo(2) )/tvolt, ( smo(1)+smo(2) )/tvolt
            WRITE (numout,9442)   tmo(3)/tvolt,  smo(3)/tvolt
            WRITE (numout,9443)   tmo(4)/tvolt,  smo(4)/tvolt
            WRITE (numout,9444)   tmo(5)/tvolt,  smo(5)/tvolt
            WRITE (numout,9445) ( t2(1)+t2(2) )/tvolt, ( s2(1)+s2(2) )/tvolt
            WRITE (numout,9446)   t2(3)/tvolt,   s2(3)/tvolt
            WRITE (numout,9447)   t2(4)/tvolt,   s2(4)/tvolt
            WRITE (numout,9448)   t2(5)/tvolt,   s2(5)/tvolt
         ENDIF

9440     FORMAT(' tracer consistency at it= ',i6,   &
            ' :         temperature','                salinity',/,   &
            ' ==================================')
9441     FORMAT(' 0 = horizontal+vertical advection      ',e20.13,'       ',e20.13)
9442     FORMAT(' 0 = horizontal diffusion               ',e20.13,'       ',e20.13)
9443     FORMAT(' 0 = vertical diffusion                 ',e20.13,'       ',e20.13)
9444     FORMAT(' 0 = static instability mixing          ',e20.13,'       ',e20.13)
9445     FORMAT(' 0 = horizontal+vertical advection * t  ',e20.13,'       ',e20.13)
9446     FORMAT(' 0 > horizontal diffusion          * t  ',e20.13,'       ',e20.13)
9447     FORMAT(' 0 > vertical diffusion            * t  ',e20.13,'       ',e20.13)
9448     FORMAT(' 0 > static instability mixing     * t  ',e20.13,'       ',e20.13)

      ENDIF

   END SUBROUTINE trd_twr

#   else
   !!----------------------------------------------------------------------
   !!   Default case :                                         Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .FALSE.   !: tracers  trend flag
   LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .FALSE.   !: momentum trend flag
CONTAINS
   SUBROUTINE trd_2d(ptrd2dx, ptrd2dy, ktrd , ctype)       ! Empty routine
      REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
          ptrd2dx,                     &   ! Temperature or U trend 
          ptrd2dy                          ! Salinity    or V trend
      INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
!      WRITE(*,*) 'trd_2d: You should not have seen this print! error ?', ptrd2dx(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd2dy(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
   END SUBROUTINE trd_2d
   SUBROUTINE trd_3d(ptrd3dx, ptrd3dy, ktrd , ctype)       ! Empty routine
      REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
          ptrd3dx,                     &   ! Temperature or U trend 
          ptrd3dy                          ! Salinity    or V trend
      INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
      CHARACTER(len=3), INTENT( in ) ::   &
         ctype                             ! momentum or tracers trends type
!      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
!      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
   END SUBROUTINE trd_3d
   SUBROUTINE trd_icp_init               ! Empty routine
   END SUBROUTINE trd_icp_init
   SUBROUTINE trd_dwr( kt )          ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trd_dwr: You should not have seen this print! error ?', kt
   END SUBROUTINE trd_dwr
   SUBROUTINE trd_twr( kt )          ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trd_twr: You should not have seen this print! error ?', kt
   END SUBROUTINE trd_twr
#endif

   !!======================================================================
END MODULE trdicp
