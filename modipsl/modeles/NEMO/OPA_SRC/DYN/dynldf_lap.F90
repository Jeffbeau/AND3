MODULE dynldf_lap
   !!======================================================================
   !!                       ***  MODULE  dynldf_lap  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap  : update the momentum trend with the lateral diffusion
   !!                  using an iso-level harmonic operator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE ldfslp          ! iso-neutral slopes 
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_ldf_lap  ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_lap.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_ldf_lap( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is an
      !!      harmonic operator (laplacian type) which separates the divergent
      !!      and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!         difu = 1/e1u di[ahmt hdivb] - 1/(e2u*e3u) dj-1[e3f ahmf rotb]
      !!         difv = 1/e2v dj[ahmt hdivb] + 1/(e1v*e3v) di-1[e3f ahmf rotb]
      !!      If 'key_s_coord' key is not activated, the vertical scale factor
      !!      is simplified in the rotational part of the diffusion.
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' activated: the two components of the horizontal
      !!                                 diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level harmonic 
      !!               mixing trend.
      !!             - Save in (ztdua,ztdva) arrays the trends ('key_trddyn')
      !!
      !! History :
      !!        !  90-09 (G. Madec) Original code
      !!        !  91-11 (G. Madec)
      !!        !  96-01 (G. Madec) statement function for e3 and ahm
      !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08 (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
                         ztdva => sa      ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   &
         zua, zva, ze2u, ze1v             ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacien) operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               ze2u = rotb (ji,jj,jk)*fsahmf(ji,jj,jk)*fse3f(ji,jj,jk)
               ze1v = hdivb(ji,jj,jk)*fsahmt(ji,jj,jk)
               ! horizontal diffusive trends
               zua = - ( ze2u - rotb (ji,jj-1,jk)*fsahmf(ji,jj-1,jk)*fse3f(ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     + ( hdivb(ji+1,jj,jk)*fsahmt(ji+1,jj,jk) - ze1v                   ) / e1u(ji,jj)

               zva = + ( ze2u - rotb (ji-1,jj,jk)*fsahmf(ji-1,jj,jk)*fse3f(ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     + ( hdivb(ji,jj+1,jk)*fsahmt(ji,jj+1,jk) - ze1v                   ) / e2v(ji,jj)
#else
               ! horizontal diffusive trends
               ze2u = rotb (ji,jj,jk)*fsahmf(ji,jj,jk)
               ze1v = hdivb(ji,jj,jk)*fsahmt(ji,jj,jk)
               zua = - (                ze2u                  - rotb (ji,jj-1,jk)*fsahmf(ji,jj-1,jk) ) / e2u(ji,jj)   &
                     + ( hdivb(ji+1,jj,jk)*fsahmt(ji+1,jj,jk) -                ze1v                  ) / e1u(ji,jj)

               zva = + (                ze2u                  - rotb (ji-1,jj,jk)*fsahmf(ji-1,jj,jk) ) / e1v(ji,jj)   &
                     + ( hdivb(ji,jj+1,jk)*fsahmt(ji,jj+1,jk) -                ze1v                  ) / e2v(ji,jj)
#endif

               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the lateral diffusion trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(ztdua, ztdva, jpdtdldf, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_ldf_lap

   !!======================================================================
END MODULE dynldf_lap
