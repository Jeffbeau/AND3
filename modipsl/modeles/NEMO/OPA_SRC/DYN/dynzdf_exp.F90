MODULE dynzdf_exp
   !!==============================================================================
   !!                     ***  MODULE  dynzdf_exp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dyn_zdf_exp  : update the momentum trend with the vertical diffu-
   !!                  sion using an explicit time-stepping scheme.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE taumod          ! surface ocean stress
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_zdf_exp    ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_exp.F90,v 1.5 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_zdf_exp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!
      !! ** Method  :   Explicit forward time stepping with a time splitting
      !!      technique. The vertical diffusion of momentum is given by:
      !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
      !!      Surface boundary conditions: wind stress input
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) with the vertical diffusive trend
      !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
      !!
      !! History :
      !!        !  90-10  (B. Blanke)  Original code
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!---------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
                         ztdva => sa       ! use sa as 3D workspace   
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index

      !! * Local declarations
      INTEGER ::   &
         ji, jj, jk, jl,                 & ! dummy loop indices
         ikbu, ikbum1 , ikbv, ikbvm1       ! temporary integers
      REAL(wp) ::   &
         zrau0r, zlavmr, z2dt, zua, zva    ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpk) ::    &
         zwx, zwy, zwz, zww                ! temporary workspace arrays
      REAL(wp), DIMENSION(jpi,jpj) ::    &
         ztsx, ztsy, ztbx, ztby            ! temporary workspace arrays
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_exp : vertical momentum diffusion explicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! Local constant initialization
      ! -----------------------------
      zrau0r = 1. / rau0                                   ! inverse of the reference density
      zlavmr = 1. / float( n_zdfexp )                      ! inverse of the number of sub time step
      z2dt = 2. * rdt                                      ! Leap-frog environnement
      ztsx(:,:) = 0.e0
      ztsy(:,:) = 0.e0 
      ztbx(:,:) = 0.e0
      ztby(:,:) = 0.e0

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt    ! Euler time stepping when starting from rest

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

         ! Surface boundary condition
         DO ji = 2, jpim1
            zwy(ji,1) = taux(ji,jj) * zrau0r
            zww(ji,1) = tauy(ji,jj) * zrau0r
         END DO  

         ! Initialization of x, z and contingently trends array
         DO jk = 1, jpk
            DO ji = 2, jpim1
               zwx(ji,jk) = ub(ji,jj,jk)
               zwz(ji,jk) = vb(ji,jj,jk)
            END DO  
         END DO  

         ! Time splitting loop
         DO jl = 1, n_zdfexp

            ! First vertical derivative
            DO jk = 2, jpk
               DO ji = 2, jpim1
                  zwy(ji,jk) = avmu(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) / fse3uw(ji,jj,jk) 
                  zww(ji,jk) = avmv(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) / fse3vw(ji,jj,jk)
               END DO  
            END DO  

            ! Second vertical derivative and trend estimation at kt+l*rdt/n_zdfexp
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  zua = zlavmr*( zwy(ji,jk) - zwy(ji,jk+1) ) / fse3u(ji,jj,jk)
                  zva = zlavmr*( zww(ji,jk) - zww(ji,jk+1) ) / fse3v(ji,jj,jk)
                  ua(ji,jj,jk) = ua(ji,jj,jk) + zua
                  va(ji,jj,jk) = va(ji,jj,jk) + zva

                  zwx(ji,jk) = zwx(ji,jk) + z2dt*zua*umask(ji,jj,jk)
                  zwz(ji,jk) = zwz(ji,jk) + z2dt*zva*vmask(ji,jj,jk)
               END DO  
            END DO  

         END DO  

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the vertical diffusion trends for diagnostic
      ! momentum trends
      IF( l_trddyn )  THEN 
         ! save the total vertical momentum diffusive trend
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
 
         ! subtract and save surface and momentum fluxes
         !                                                ! ===============
         DO jj = 2, jpjm1                                 !  Horizontal slab
            !                                             ! ===============
            DO ji = 2, jpim1
               ! save the surface momentum fluxes 
               ztsx(ji,jj) = zwy(ji,1) / fse3u(ji,jj,1)
               ztsy(ji,jj) = zww(ji,1) / fse3v(ji,jj,1)
               ! save bottom friction momentum fluxes 
               ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
               ikbum1 = MAX( ikbu-1, 1 )
               ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
               ikbvm1 = MAX( ikbv-1, 1 )
               ztbx(ji,jj) = avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
                               / ( fse3u(ji,jj,ikbum1) * fse3uw(ji,jj,ikbu) )
               ztby(ji,jj) = avmv(ji,jj,ikbv) * zwz(ji,ikbvm1)   &
                               / ( fse3v(ji,jj,ikbvm1) * fse3vw(ji,jj,ikbv) )
               ! subtract surface forcing and bottom friction trend from vertical
               ! diffusive momentum trend
               ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
               ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
               ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
               ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
         ztdua(:,:,:) = 0.e0
         ztdva(:,:,:) = 0.e0
         ztdua(:,:,1) = ztsx(:,:)
         ztdva(:,:,1) = ztsy(:,:)
         CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
         ztdua(:,:,:) = 0.e0
         ztdva(:,:,:) = 0.e0
         ztdua(:,:,1) = ztbx(:,:)
         ztdva(:,:,1) = ztby(:,:)
         CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_zdf_exp

   !!==============================================================================
END MODULE dynzdf_exp
