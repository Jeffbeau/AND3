MODULE dynvor
   !!======================================================================
   !!                       ***  MODULE  dynvor  ***
   !! Ocean dynamics: Update the momentum trend with the relative and
   !!                 planetary vorticity trends
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dyn_vor_enstrophy: enstrophy conserving scheme       (ln_dynvor_ens=T)
   !!   dyn_vor_energy   : energy conserving scheme          (ln_dynvor_ene=T)
   !!   dyn_vor_mixed    : mixed enstrophy/energy conserving (ln_dynvor_mix=T)
   !!   dyn_vor_ene_ens  : energy and enstrophy conserving   (ln_dynvor_een=T)
   !!   dyn_vor_ctl      : control of the different vorticity option
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_vor_enstrophy      ! routine called by step.F90
   PUBLIC dyn_vor_energy         ! routine called by step.F90
   PUBLIC dyn_vor_mixed          ! routine called by step.F90
   PUBLIC dyn_vor_ene_ens        ! routine called by step.F90
   PUBLIC dyn_vor_ctl            ! routine called by step.F90

   !! * Shared module variables
   LOGICAL, PUBLIC ::   ln_dynvor_ene = .FALSE.   !: energy conserving scheme
   LOGICAL, PUBLIC ::   ln_dynvor_ens = .TRUE.    !: enstrophy conserving scheme
   LOGICAL, PUBLIC ::   ln_dynvor_mix = .FALSE.   !: mixed scheme
   LOGICAL, PUBLIC ::   ln_dynvor_een = .FALSE.    !: energy and enstrophy conserving scheme

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynvor.F90,v 1.11 2006/03/21 08:26:27 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_vor_energy( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_vor_energy  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to 
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!      and the Sadourny (1975) flux form formulation : conserves the
      !!      horizontal kinetic energy.
      !!      The trend of the vorticity term is given by:
      !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f  mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f  mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ (rotn+f)  mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)  mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (utrd,vtrd) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References :
      !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !! History :
      !!   5.0  !  91-11  (G. Madec)  Original code
      !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
                         ztdva => sa        ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   &
         zfact2, zua, zva,               &  ! temporary scalars
         zx1, zx2, zy1, zy2                 !    "         "
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zwx, zwy, zwz                      ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zcu, zcv                           !    "         "
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_energy : vorticity term: energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ! Local constant initialization
      zfact2 = 0.5 * 0.5
      zcu (:,:,:) = 0.e0
      zcv (:,:,:) = 0.e0

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
         zcu(:,:,:) = 0.e0
         zcv(:,:,:) = 0.e0
      ENDIF
      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         IF( lk_sco ) THEN
            zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / fse3f(:,:,jk)
            zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
         ELSE
            zwz(:,:) = rotn(:,:,jk) + ff(:,:)
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
               zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
               zua = zfact2 / e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zva =-zfact2 / e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the relative & planetary vorticity trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ! Compute the planetary vorticity term trend
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
                  zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
                  zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
                  zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
# if defined key_s_coord
                 zcu(ji,jj,jk) = zfact2 / e1u(ji,jj) * ( ff(ji  ,jj-1) / fse3f(ji,jj-1,jk) * zy1  &
                   &                                   + ff(ji  ,jj  ) / fse3f(ji,jj  ,jk) * zy2 )
                 zcv(ji,jj,jk) =-zfact2 / e2v(ji,jj) * ( ff(ji-1,jj  ) / fse3f(ji-1,jj,jk) * zx1  &
                   &                                   + ff(ji  ,jj  ) / fse3f(ji  ,jj,jk) * zx2 )
# else
                 zcu(ji,jj,jk) = zfact2 / e1u(ji,jj) * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
                 zcv(ji,jj,jk) =-zfact2 / e2v(ji,jj) * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
# endif
               END DO  
            END DO  
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         ! Compute the relative vorticity term trend
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 

         CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
         CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
         CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_vor_energy


   SUBROUTINE dyn_vor_mixed( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_vor_mixed  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      Mixte formulation : conserves the potential enstrophy of a hori-
      !!      zontally non-divergent flow for (rotzu x uh), the relative vor-
      !!      ticity term and the horizontal kinetic energy for (f x uh), the
      !!      coriolis term. the now trend of the vorticity term is given by:
      !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1(rotn/e3f) mj-1[ mi(e1v*e3v vn) ]
      !!              +1/e1u  mj-1[ f/e3f          mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn/e3f) mi-1[ mj(e2u*e3u un) ]
      !!              +1/e2v  mi-1[ f/e3f          mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1(rotn) mj-1[ mi(e1v vn) ]
      !!              +1/e1u  mj-1[ f          mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn) mi-1[ mj(e2u un) ]
      !!              +1/e2v  mi-1[ f          mj(e2u un) ]
      !!      Add this now trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!             - Save the trends in (utrd,vtrd) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References :
      !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !! History :
      !!   5.0  !  91-11  (G. Madec) Original code, enstrophy-energy-combined schemes
      !!   6.0  !  96-01  (G. Madec) s-coord, suppress work arrays
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,     & ! use ta as 3D workspace   
                         ztdva => sa        ! use sa as 3D workspace   
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean timestep index

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   &
         zfact1, zfact2, zua, zva,       &  ! temporary scalars
         zcua, zcva, zx1, zx2, zy1, zy2
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zwx, zwy, zwz, zww                 ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zcu, zcv                           !    "         "
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_mixed : vorticity term: mixed energy/enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ! Local constant initialization
      zfact1 = 0.5 * 0.25
      zfact2 = 0.5 * 0.5

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
         zcu(:,:,:) = 0.e0
         zcv(:,:,:) = 0.e0
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Relative and planetary potential vorticity and horizontal fluxes
         ! ----------------------------------------------------------------
         IF( lk_sco ) THEN        
            zwz(:,:) = ff  (:,:)    / fse3f(:,:,jk)
            zww(:,:) = rotn(:,:,jk) / fse3f(:,:,jk)
            zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
         ELSE
            zwz(:,:) = ff(:,:)
            zww(:,:) = rotn(:,:,jk)
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               ! enstrophy conserving formulation for relative vorticity term
               zua = zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1 + zy2 )
               zva =-zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1 + zx2 )
               ! energy conserving formulation for planetary vorticity term
               zcua = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zcva =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )

               ua(ji,jj,jk) = ua(ji,jj,jk) + zcua + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zcva + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the relative & planetary vorticity trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ! Compute the planetary vorticity term trend
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
                  zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
                  zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
                  zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
                  ! energy conserving formulation for planetary vorticity term
                  zcu(ji,jj,jk) = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
                  zcv(ji,jj,jk) =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               END DO  
            END DO  
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         ! Compute the relative vorticity term trend
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 

         CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
         CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
         CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_vor_mixed


   SUBROUTINE dyn_vor_enstrophy( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dyn_vor_enstrophy  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!       * s-coordinate (lk_sco=T), the e3. are inside the derivative:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!             - Save the trends in (utrd,vtrd) in 2 parts (relative 
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References :
      !!      Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !! History :
      !!   5.0  !  91-11  (G. Madec)  Original code
      !!   6.0  !  96-01  (G. Madec)  s-coord, suppress work arrays
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY:   zwx  => ta,   & ! use ta as 3D workspace
                       zwy  => sa      ! use sa as 3D workspace
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean timestep

      !! * Local declarations
      INTEGER ::   ji, jj, jk          ! dummy loop indices
      REAL(wp) ::   &
         zfact1, zua, zva, zuav, zvau  ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zcu, zcv, zwz,              & ! temporary workspace
         ztdua, ztdva                  ! temporary workspace
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_enstrophy : vorticity term: enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~'
      ENDIF

      ! Local constant initialization
      zfact1 = 0.5 * 0.25

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
         zcu(:,:,:) = 0.e0
         zcv(:,:,:) = 0.e0
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         IF( lk_sco ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwz(ji,jj,jk) = ( rotn(ji,jj,jk) + ff(ji,jj) ) / fse3f(ji,jj,jk)
                  zwx(ji,jj,jk) = e2u(ji,jj) * fse3u(ji,jj,jk) * un(ji,jj,jk)
                  zwy(ji,jj,jk) = e1v(ji,jj) * fse3v(ji,jj,jk) * vn(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop 
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwz(ji,jj,jk) = rotn(ji,jj,jk) + ff(ji,jj)
                  zwx(ji,jj,jk) = e2u(ji,jj) * un(ji,jj,jk)
                  zwy(ji,jj,jk) = e1v(ji,jj) * vn(ji,jj,jk)
               END DO
            END DO
         ENDIF


         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
                                            + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
               zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
                                            + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )

               zua  = zuav * ( zwz(ji  ,jj-1,jk) + zwz(ji,jj,jk) )
               zva  = zvau * ( zwz(ji-1,jj  ,jk) + zwz(ji,jj,jk) )

               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============


      ! save the relative & planetary vorticity trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ! Compute the planetary vorticity term trend
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1,jk) + zwy(ji+1,jj-1,jk)   &
                     &                         + zwy(ji  ,jj  ,jk) + zwy(ji+1,jj  ,jk) )
                  zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ,jk) + zwx(ji-1,jj+1,jk)   &
                     &                         + zwx(ji  ,jj  ,jk) + zwx(ji  ,jj+1,jk) )
# if defined key_s_coord
                  zcu(ji,jj,jk)  = zuav * ( ff(ji  ,jj-1) / fse3f(ji  ,jj-1,jk)   &
                    &                     + ff(ji  ,jj  ) / fse3f(ji  ,jj  ,jk) )
                  zcv(ji,jj,jk)  = zvau * ( ff(ji-1,jj  ) / fse3f(ji-1,jj  ,jk)   &
                    &                     + ff(ji  ,jj  ) / fse3f(ji  ,jj  ,jk) )
# else
                  zcu(ji,jj,jk) = zuav * ( ff(ji  ,jj-1) + ff(ji,jj) )
                  zcv(ji,jj,jk) = zvau * ( ff(ji-1,jj  ) + ff(ji,jj) )
# endif
               END DO  
            END DO  
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         ! Compute the relative vorticity term trend
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 

         CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
         CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
         CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_vor_enstrophy


   SUBROUTINE dyn_vor_ene_ens( kt )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE dyn_vor_ene_ens  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to 
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!      and the Arakawa and Lamb (19XX) flux form formulation : conserves 
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero.
      !!      The trend of the vorticity term is given by:
      !!       * s-coordinate (lk_sco=T), the e3. are inside the derivatives:
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (utrd,vtrd) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References :
      !!      Arakawa and Lamb 1980, A potential enstrophy and energy conserving
      !!                             scheme for the Shallow water equations, 
      !!                             Monthly Weather Review, vol. 109, p 18-36
      !!
      !! History :
      !!   8.5  !  04-02  (G. Madec)  Original code
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
                         ztdva => sa       ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk              ! dummy loop indices
      REAL(wp) ::   &
         zfac12, zua, zva                  ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zwx, zwy, zwz,                 &  ! temporary workspace
         ztnw, ztne, ztsw, ztse,        &  !    "           "
         zcor                              ! potential planetary vorticity (f/e3)
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zcu, zcv                          ! temporary workspace  
      REAL(wp), DIMENSION(jpi,jpj,jpk), SAVE ::   &
         ze3f
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_ene_ens : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'

         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ze3f(ji,jj,jk) = ( fse3t(ji,jj+1,jk)*tmask(ji,jj+1,jk) + fse3t(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                     &             + fse3t(ji,jj  ,jk)*tmask(ji,jj  ,jk) + fse3t(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) ) * 0.25
!!!               ze3f(ji,jj,jk) = MAX( ze3f(ji,jj,jk) , 1.e-20)
                  IF( ze3f(ji,jj,jk) /= 0.e0 )   ze3f(ji,jj,jk) = 1.e0 / ze3f(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( ze3f, 'F', 1. )
      ENDIF

      ! Local constant initialization
      zfac12 = 1.e0 / 12.e0

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
         zcu(:,:,:) = 0.e0
         zcv(:,:,:) = 0.e0
      ENDIF
      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
!!!bug   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) / fse3f(:,:,jk)
         zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f(:,:,jk)
         zwx(:,:) = e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
         zwy(:,:) = e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
         zcor(:,:) = ff(:,:) * ze3f(:,:,jk)

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         jj=2
         ztne(1,:) = 0 ; ztnw(1,:) = 0 ; ztse(1,:) = 0 ; ztsw(1,:) = 0
         DO ji = 2, jpi   
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END DO
         DO jj = 3, jpj
            DO ji = fs_2, jpi   ! vector opt.
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
            END DO
         END DO
         
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zua = + zfac12 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                  &                           + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zva = - zfac12 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                  &                           + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the relative & planetary vorticity trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ! Compute the planetary vorticity term trend
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zcu(ji,jj,jk) = + zfac12 / e1u(ji,jj) * (  zcor(ji,jj  ) * zwy(ji  ,jj  ) + zcor(ji+1,jj) * zwy(ji+1,jj  )   &
                    &                                      + zcor(ji,jj  ) * zwy(ji  ,jj-1) + zcor(ji+1,jj) * zwy(ji+1,jj-1) )
                  zcv(ji,jj,jk) = - zfac12 / e2v(ji,jj) * (  zcor(ji,jj+1) * zwx(ji-1,jj+1) + zcor(ji,jj+1) * zwx(ji  ,jj+1)   &
                    &                                      + zcor(ji,jj  ) * zwx(ji-1,jj  ) + zcor(ji,jj  ) * zwx(ji  ,jj  ) )
               END DO  
            END DO  
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         ! Compute the relative vorticity term trend
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:) - zcu(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:) - zcv(:,:,:) 

         CALL trd_mod(zcu  , zcv  , jpdtdpvo, 'DYN', kt)
         CALL trd_mod(zcu  , zcv  , jpdtddat, 'DYN', kt)
         CALL trd_mod(ztdua, ztdva, jpdtdrvo, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' vor  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_vor_ene_ens


   SUBROUTINE dyn_vor_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_vor_ctl  ***
      !!
      !! ** Purpose :   Control the consistency between cpp options for
      !!      tracer advection schemes
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ioptio          ! temporary integer

      NAMELIST/nam_dynvor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_mix, ln_dynvor_een
      !!----------------------------------------------------------------------

      ! Read Namelist nam_dynvor : Vorticity scheme options
      ! ------------------------
      REWIND ( numnam )
      READ   ( numnam, nam_dynvor )

      ! Control of vorticity scheme options
      ! -----------------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_vor_ctl : vorticity term : read namelist and control the consistency'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '              Namelist nam_dynvor : oice of the vorticity term scheme'
         WRITE(numout,*) '                 enstrophy conserving scheme                ln_dynvor_ens = ', ln_dynvor_ens
         WRITE(numout,*) '                 energy    conserving scheme                ln_dynvor_ene = ', ln_dynvor_ene
         WRITE(numout,*) '                 mixed enstrophy/energy conserving scheme   ln_dynvor_mix = ', ln_dynvor_mix
         WRITE(numout,*) '                 enstrophy and energy conserving scheme     ln_dynvor_een = ', ln_dynvor_een
      ENDIF

      ioptio = 0

      IF( ln_dynvor_ens ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              vorticity term : enstrophy conserving scheme'
         ioptio = ioptio + 1
      ENDIF
      IF( ln_dynvor_ene ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              vorticity term : energy conserving scheme'
         ioptio = ioptio + 1
      ENDIF
      IF( ln_dynvor_mix ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              vorticity term : mixed enstrophy/energy conserving scheme'
         ioptio = ioptio + 1
      ENDIF
      IF( ln_dynvor_een ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              vorticity term : energy and enstrophy conserving scheme'
         ioptio = ioptio + 1
      ENDIF
      IF ( ioptio /= 1 .AND. .NOT. lk_esopa ) THEN
          if(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) ' use ONE and ONLY one vorticity scheme'
          nstop = nstop + 1
      ENDIF

   END SUBROUTINE dyn_vor_ctl

!!==============================================================================
END MODULE dynvor
