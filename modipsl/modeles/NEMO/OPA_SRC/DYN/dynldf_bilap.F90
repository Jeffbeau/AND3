MODULE dynldf_bilap
   !!======================================================================
   !!                     ***  MODULE  dynldf_bilap  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dyn_ldf_bilap : update the momentum trend with the lateral diffusion
   !!                   using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_ldf_bilap  ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_bilap.F90,v 1.8 2005/09/02 15:45:23 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_ldf_bilap( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before trend of the lateral momentum
      !!      diffusion and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is a 
      !!      bi-harmonic operator (bilaplacian type) which separates the
      !!      divergent and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!      laplacian:
      !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
      !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
      !!      third derivative:
      !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
      !!          zlu = ahmu * zlu
      !!          zlv = ahmv * zlv
      !!       * curl and divergence of the laplacian
      !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
      !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
      !!      bilaplacian:
      !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
      !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
      !!      If lk_sco=F and lk_zps=F, the vertical scale factors in the
      !!      rotational part of the diffusion are simplified
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' defined: the two components of the horizontal
      !!                               diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level biharmonic
      !!               mixing trend.
      !!             - Save in (ztdua,ztdva) the trends ('key_trddyn')
      !!
      !! History :
      !!        !  90-09  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-07  (G. Madec)  lbc calls
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,      & ! use ta as 3D workspace   
                         ztdva => sa         ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   zua, zva, zbt, ze2u, ze2v ! temporary scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zuf, zut, zlu, zlv, zcu, zcv         ! temporary workspace
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_bilap : iso-level bilaplacian operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      zuf(:,:) = 0.e0
      zut(:,:) = 0.e0
      zlu(:,:) = 0.e0
      zlv(:,:) = 0.e0

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Laplacian
         ! ---------

         IF( lk_sco .OR. lk_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
            zuf(:,:) = rotb(:,:,jk) * fse3f(:,:,jk)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlu(ji,jj) = - ( zuf(ji,jj) - zuf(ji,jj-1) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) )   &
                     &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj,jk) ) / e1u(ji,jj)
   
                  zlv(ji,jj) = + ( zuf(ji,jj) - zuf(ji-1,jj) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) )   &
                     &         + ( hdivb(ji,jj+1,jk) - hdivb(ji,jj,jk) ) / e2v(ji,jj)
               END DO
            END DO
         ELSE                            ! z-coordinate
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlu(ji,jj) = - ( rotb (ji  ,jj,jk) - rotb (ji,jj-1,jk) ) / e2u(ji,jj)   &
                     &         + ( hdivb(ji+1,jj,jk) - hdivb(ji,jj  ,jk) ) / e1u(ji,jj)
   
                  zlv(ji,jj) = + ( rotb (ji,jj  ,jk) - rotb (ji-1,jj,jk) ) / e1v(ji,jj)   &
                     &         + ( hdivb(ji,jj+1,jk) - hdivb(ji  ,jj,jk) ) / e2v(ji,jj)
               END DO  
            END DO  
         ENDIF

         ! Boundary conditions on the laplacian  (zlu,zlv)
         CALL lbc_lnk( zlu, 'U', -1. )
         CALL lbc_lnk( zlv, 'V', -1. )
         
         
         ! Third derivative
         ! ----------------
         
         ! Multiply by the eddy viscosity coef. (at u- and v-points)
         zlu(:,:) = zlu(:,:) * fsahmu(:,:,jk)
         zlv(:,:) = zlv(:,:) * fsahmv(:,:,jk)
         
         ! Contravariant "laplacian"
         zcu(:,:) = e1u(:,:) * zlu(:,:)
         zcv(:,:) = e2v(:,:) * zlv(:,:)
         
         ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zuf(ji,jj) = fmask(ji,jj,jk) * (  zcv(ji+1,jj  ) - zcv(ji,jj)      &
                  &                            - zcu(ji  ,jj+1) + zcu(ji,jj)  )   &
#if defined key_s_coord || defined key_partial_steps
                  &       * fse3f(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
#else
                  &                         / ( e1f(ji,jj)*e2f(ji,jj) )
#endif
            END DO  
         END DO  

         ! Laplacian Horizontal fluxes
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               zlu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) * zlu(ji,jj)
               zlv(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) * zlv(ji,jj)
#else
               zlu(ji,jj) = e2u(ji,jj) * zlu(ji,jj)
               zlv(ji,jj) = e1v(ji,jj) * zlv(ji,jj)
#endif
            END DO
         END DO

         ! Laplacian divergence
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
#else
               zbt = e1t(ji,jj) * e2t(ji,jj)
#endif
               zut(ji,jj) = (  zlu(ji,jj) - zlu(ji-1,jj  )   &
                  &          + zlv(ji,jj) - zlv(ji  ,jj-1)  ) / zbt
            END DO
         END DO


      ! boundary conditions on the laplacian curl and div (zuf,zut)
      CALL lbc_lnk( zuf, 'F', 1. )
      CALL lbc_lnk( zut, 'T', 1. )

         
         ! Bilaplacian
         ! -----------

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               ze2u = e2u(ji,jj) * fse3u(ji,jj,jk)
               ze2v = e1v(ji,jj) * fse3v(ji,jj,jk)
#else
               ze2u = e2u(ji,jj)
               ze2v = e1v(ji,jj)
#endif
               ! horizontal biharmonic diffusive trends
               zua = - ( zuf(ji  ,jj) - zuf(ji,jj-1) ) / ze2u   &
                  &  + ( zut(ji+1,jj) - zut(ji,jj  ) ) / e1u(ji,jj)

               zva = + ( zuf(ji,jj  ) - zuf(ji-1,jj) ) / ze2v   &
                  &  + ( zut(ji,jj+1) - zut(ji  ,jj) ) / e2v(ji,jj)
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

   END SUBROUTINE dyn_ldf_bilap

   !!======================================================================
END MODULE dynldf_bilap
