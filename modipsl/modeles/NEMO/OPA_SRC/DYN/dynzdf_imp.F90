MODULE dynzdf_imp
   !!==============================================================================
   !!                    ***  MODULE  dynzdf_imp  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dyn_zdf_imp  : update the momentum trend with the vertical diffu-
   !!                  sion using a implicit time-stepping.
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
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
   PUBLIC dyn_zdf_imp    ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS


   SUBROUTINE dyn_zdf_imp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_imp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
      !!      and the surface forcing, and add it to the general trend of 
      !!      the momentum equations.
      !!
      !! ** Method  :   The vertical momentum mixing trend is given by :
      !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
      !!      backward time stepping
      !!      Surface boundary conditions: wind stress input
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
      !!               mixing trend.
      !!             - Save the trends in (ztdua,ztdva) ('l_trddyn')
      !!
      !! History :
      !!        !  90-10  (B. Blanke)  Original code
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce, ONLY :  zwd   => ta,   &  ! use ta as workspace
                       zws   => sa       ! use sa as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER ::   &
         ji, jj, jk,                  &  ! dummy loop indices
         ikbu, ikbum1, ikbv, ikbvm1      ! temporary integers
      REAL(wp) ::   &
         zrau0r, z2dt,                &  ! temporary scalars
         z2dtf, zcoef, zzws, zrhs        !    "         "
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         ztsx, ztsy, ztbx, ztby          ! temporary workspace arrays
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwi, ztdua, ztdva               ! temporary workspace arrays
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      zrau0r = 1. / rau0      ! inverse of the reference density
      z2dt   = 2. * rdt       ! Leap-frog environnement
      ztsx(:,:) = 0.e0
      ztsy(:,:) = 0.e0 
      ztbx(:,:) = 0.e0
      ztby(:,:) = 0.e0
      ! Euler time stepping when starting from rest
      IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt

      ! Save previous ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      ! 1. Vertical diffusion on u
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: only zws must be masked as avmu can take
      ! non zero value at the ocean bottom depending on the bottom friction
      ! used (see zdfmix.F)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1 
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef = - z2dt / fse3u(ji,jj,jk)
               zwi(ji,jj,jk) = zcoef * avmu(ji,jj,jk  ) / fse3uw(ji,jj,jk  )
               zzws          = zcoef * avmu(ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
               zws(ji,jj,jk) = zzws * umask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
            END DO
         END DO
      END DO

      ! Surface boudary conditions
      DO jj = 2, jpjm1 
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0.
            zwd(ji,jj,1) = 1. - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in ua
      !-----------------------------------------------------------------------

      ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
      DO jj = 2, jpjm1   
         DO ji = fs_2, fs_jpim1   ! vector opt.
!!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
!!!         ua(ji,jj,1) = ub(ji,jj,1)  &
!!!                      + z2dt * ( ua(ji,jj,1) + taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 ) )
            z2dtf = z2dt / ( fse3u(ji,jj,1)*rau0 )
            ua(ji,jj,1) = ub(ji,jj,1)  &
                         + z2dt *  ua(ji,jj,1) + z2dtf * taux(ji,jj)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zrhs = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)   ! zrhs=right hand side
               ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk
      DO jj = 2, jpjm1   
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) =( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( l_trddyn )  THEN 
         ! diagnose surface and bottom momentum fluxes
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! save the surface forcing momentum fluxes
               ztsx(ji,jj) = taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 )
               ! save bottom friction momentum fluxes
               ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
               ikbum1 = MAX( ikbu-1, 1 )
               ztbx(ji,jj) = - avmu(ji,jj,ikbu) * ua(ji,jj,ikbum1)   &
                  / ( fse3u(ji,jj,ikbum1)*fse3uw(ji,jj,ikbu) )
               ! subtract surface forcing and bottom friction trend from vertical
               ! diffusive momentum trend
               ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
               ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
            END DO
         END DO
      ENDIF

      ! Normalization to obtain the general momentum trend ua
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ua(ji,jj,jk) - ub(ji,jj,jk) ) / z2dt
            END DO
         END DO
      END DO


      ! 2. Vertical diffusion on v
      ! ---------------------------
      ! Matrix and second member construction
      ! bottom boundary condition: only zws must be masked as avmv can take
      ! non zero value at the ocean bottom depending on the bottom friction
      ! used (see zdfmix.F)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef = -z2dt / fse3v(ji,jj,jk)
               zwi(ji,jj,jk) = zcoef * avmv(ji,jj,jk  ) / fse3vw(ji,jj,jk  )
               zzws       = zcoef * avmv(ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
               zws(ji,jj,jk) =  zzws * vmask(ji,jj,jk+1)
               zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
            END DO
         END DO
      END DO

      ! Surface boudary conditions
      DO jj = 2, jpjm1   
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0.e0
            zwd(ji,jj,1) = 1. - zws(ji,jj,1)
         END DO
      END DO

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular
      !   matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------

      ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
!!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
!!!         va(ji,jj,1) = vb(ji,jj,1)  &
!!!                      + z2dt * ( va(ji,jj,1) + tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 ) )
            z2dtf = z2dt / ( fse3v(ji,jj,1)*rau0 )
            va(ji,jj,1) = vb(ji,jj,1)  &
                         + z2dt * va(ji,jj,1) + z2dtf * tauy(ji,jj)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zrhs = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)   ! zrhs=right hand side
               va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! thrid recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk
      DO jj = 2, jpjm1   
         DO ji = fs_2, fs_jpim1   ! vector opt.
            va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) =( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( l_trddyn )  THEN 
         ! diagnose surface and bottom momentum fluxes
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! save the surface forcing momentum fluxes
               ztsy(ji,jj) = tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 )
               ! save bottom friction momentum fluxes
               ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
               ikbvm1 = MAX( ikbv-1, 1 )
               ztby(ji,jj) = - avmv(ji,jj,ikbv) * va(ji,jj,ikbvm1)   &
                  / ( fse3v(ji,jj,ikbvm1)*fse3vw(ji,jj,ikbv) )
               ! subtract surface forcing and bottom friction trend from vertical
               ! diffusive momentum trend
               ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
               ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
            END DO
         END DO
      ENDIF

      ! Normalization to obtain the general momentum trend va
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1   
            DO ji = fs_2, fs_jpim1   ! vector opt.
               va(ji,jj,jk) = ( va(ji,jj,jk) - vb(ji,jj,jk) ) / z2dt
            END DO
         END DO
      END DO

      ! save the vertical diffusion trends for diagnostic
      ! momentum trends
      IF( l_trddyn )  THEN 
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

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

   END SUBROUTINE dyn_zdf_imp

   !!==============================================================================
END MODULE dynzdf_imp
