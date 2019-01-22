MODULE dynzad
   !!======================================================================
   !!                       ***  MODULE  dynzad  ***
   !! Ocean dynamics : vertical advection trend
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   dyn_zad      : vertical advection momentum trend
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE flxrnf          ! ocean runoffs
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE
   
   !! * Accessibility
   PUBLIC dyn_zad                ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzad.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking'                              j-k-i loops (j-slab)
   !!----------------------------------------------------------------------

   SUBROUTINE dyn_zad( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad  ***
      !!
      !! ** Purpose :   Compute the now vertical momentum advection trend and 
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method  :   Use j-slab (j-k-i loops) for auto-tasking
      !!      The now vertical advection of momentum is given by:
      !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
      !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
      !!      Add this trend to the general trend (ua,va):
      !!         (ua,va) = (ua,va) + w dz(u,v)
      !!
      !! ** Action  : - Update (ua,va) with the vert. momentum advection trends
      !!              - Save the trends in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   6.0  !  91-01  (G. Madec) Original code
      !!   7.0  !  91-11  (G. Madec)
      !!   7.5  !  96-01  (G. Madec) statement function for e3
      !!   8.5  !  02-07  (G. Madec) Free form, F90
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY:   zwuw => ta,   & ! use ta as 3D workspace
                       zwvw => sa      ! use sa as 3D workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step inedx
      
      !! * Local declarations
      INTEGER  ::   ji, jj, jk         ! dummy loop indices
      REAL(wp) ::   zvn, zua, zva      ! temporary scalars
      REAL(wp), DIMENSION(jpi) ::   &
         zww                           ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztdua, ztdva                  ! temporary workspace
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zad : arakawa advection scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~   Auto-tasking case, j-slab, no vector opt.'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

         ! Vertical momentum advection at level w and u- and v- vertical
         ! ----------------------------------------------------------------
         DO jk = 2, jpkm1
            ! vertical fluxes 
            DO ji = 2, jpi
               zww(ji) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
            END DO
            ! vertical momentum advection at w-point
            DO ji = 2, jpim1
               zvn = 0.25 * e1t(ji,jj+1) * e2t(ji,jj+1) * wn(ji,jj+1,jk)
               zwuw(ji,jj,jk) = ( zww(ji+1) + zww(ji) ) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
               zwvw(ji,jj,jk) = ( zvn       + zww(ji) ) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
            END DO  
         END DO   

         ! Surface and bottom values set to zero
         DO ji = 2, jpim1
            zwuw(ji,jj, 1 ) = 0.e0
            zwvw(ji,jj, 1 ) = 0.e0
            zwuw(ji,jj,jpk) = 0.e0
            zwvw(ji,jj,jpk) = 0.e0
         END DO  

         ! Vertical momentum advection at u- and v-points
         ! ----------------------------------------------
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ! vertical momentum advective trends
               zua = - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) )
               zva = - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) )
               ! add the trends to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO  
         END DO  
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the vertical advection trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(ztdua, ztdva, jpdtdzad, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_zad

#else
   !!----------------------------------------------------------------------
   !!   Default option                             k-j-i loop (vector opt.)
   !!----------------------------------------------------------------------

   SUBROUTINE dyn_zad ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad  ***
      !! 
      !! ** Purpose :   Compute the now vertical momentum advection trend and 
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The now vertical advection of momentum is given by:
      !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
      !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
      !!      Add this trend to the general trend (ua,va):
      !!         (ua,va) = (ua,va) + w dz(u,v)
      !!
      !! ** Action  : - Update (ua,va) with the vert. momentum adv. trends
      !!              - Save the trends in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   8.5  !  02-07  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY:   zwuw => ta,   & ! use ta as 3D workspace
                       zwvw => sa      ! use sa as 3D workspace
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt    ! ocean time-step inedx
      
      !! * Local declarations
      INTEGER  ::   ji, jj, jk         ! dummy loop indices
      REAL(wp) ::   zua, zva           ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zww                           ! temporary  workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztdua, ztdva                  ! temporary workspace
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'dyn_zad : arakawa advection scheme'
         IF(lwp)WRITE(numout,*) '~~~~~~~   vector optimization k-j-i loop'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF
      
      ! Vertical momentum advection at level w and u- and v- vertical
      ! -------------------------------------------------------------
      DO jk = 2, jpkm1
         ! vertical fluxes 
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               zww(ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
            END DO
         END DO
         ! vertical momentum advection at w-point
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwuw(ji,jj,jk) = ( zww(ji+1,jj  ) + zww(ji,jj) ) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
               zwvw(ji,jj,jk) = ( zww(ji  ,jj+1) + zww(ji,jj) ) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
            END DO  
         END DO   
      END DO

      ! Surface and bottom values set to zero
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwuw(ji,jj, 1 ) = 0.e0
            zwvw(ji,jj, 1 ) = 0.e0
            zwuw(ji,jj,jpk) = 0.e0
            zwvw(ji,jj,jpk) = 0.e0
         END DO  
      END DO


      ! Vertical momentum advection at u- and v-points
      ! ----------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! vertical momentum advective trends
               zua = - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) )
               zva = - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) )
               ! add the trends to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO  
         END DO  
      END DO

      ! save the vertical advection trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(ztdua, ztdva, jpdtdzad, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_zad
#endif

!!======================================================================
END MODULE dynzad
