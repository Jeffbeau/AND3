MODULE traldf_bilap
   !!==============================================================================
   !!                   ***  MODULE  traldf_bilap  ***
   !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_ldf_bilap : update the tracer trend with the horizontal diffusion
   !!                   using a iso-level biharmonic operator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE ldfslp          ! iso-neutral slopes 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE diaptr          ! poleward transport diagnostics
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_ldf_bilap   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "ldfeiv_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_bilap.F90,v 1.6 2005/09/02 15:45:33 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE tra_ldf_bilap( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces 
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends of temperature (idem for salinity) is given by:
      !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
      !!      factors e3. are inside the derivatives:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(tb) ]
      !!                              + dj-1[ e1v/e2v dj(tb) ] }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
      !!                              + dj-1[ e1v/e2v dj(zlt) ]  }
      !!
      !!      Add this trend to the general trend (ta,sa):
      !!         (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action : - Update (ta,sa) arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!             - Save the trends in (ztdta,ztdsa) ('key_trdtra')
      !!
      !! History :
      !!        !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  96-01  (M. Imbard)  mpp exchange
      !!        !  97-07  (G. Madec)  optimization, and ahtt
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce           , ztu => ua,  &  ! use ua as workspace
         &                ztv => va      ! use va as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk             ! dummy loop indices
#if defined key_partial_steps
      INTEGER ::   iku, ikv               ! temporary integers
#endif
      REAL(wp) ::   zta, zsa              ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zeeu, zeev, zbtr,              & ! workspace
         zlt, zls
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
         zsu, zsv,                          & ! workspace arrays
         ztdta, ztdsa
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_bilap : iso-level biharmonic operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! 0. Initialization of metric arrays (for z- or s-coordinates)
         ! ----------------------------------

         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               ! s-coordinates, vertical scale factor are used
               zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
               zeeu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
               zeev(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
#else
               ! z-coordinates, no vertical scale factors
               zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
               zeeu(ji,jj) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
               zeev(ji,jj) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
#endif
            END DO
         END DO


         ! 1. Laplacian
         ! ------------

         ! First derivative (gradient)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ztu(ji,jj,jk) = zeeu(ji,jj) * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
               zsu(ji,jj,jk) = zeeu(ji,jj) * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
               ztv(ji,jj,jk) = zeev(ji,jj) * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
               zsv(ji,jj,jk) = zeev(ji,jj) * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
            END DO
         END DO
#if defined key_partial_steps
         DO jj = 1, jpj-1
            DO ji = 1, jpi-1
               ! last level
               iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
               ikv = MIN ( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
               IF( iku == jk ) THEN
                  ztu(ji,jj,jk) = zeeu(ji,jj) * gtu(ji,jj)
                  zsu(ji,jj,jk) = zeeu(ji,jj) * gsu(ji,jj)
               ENDIF
               IF( ikv == jk ) THEN
                  ztv(ji,jj,jk) = zeev(ji,jj) * gtv(ji,jj)
                  zsv(ji,jj,jk) = zeev(ji,jj) * gsv(ji,jj)
               ENDIF
            END DO
         END DO
#endif

         ! Second derivative (divergence)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zlt(ji,jj) = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               zls(ji,jj) = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
            END DO
         END DO

         ! Multiply by the eddy diffusivity coefficient
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zlt(ji,jj) = fsahtt(ji,jj,jk) * zlt(ji,jj)
               zls(ji,jj) = fsahtt(ji,jj,jk) * zls(ji,jj)
            END DO
         END DO

         ! Lateral boundary conditions on the laplacian (zlt,zls)   (unchanged sgn)
         CALL lbc_lnk( zlt, 'T', 1. )   ;    CALL lbc_lnk( zls, 'T', 1. )

         ! 2. Bilaplacian
         ! --------------

         ! third derivative (gradient)
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
               zsu(ji,jj,jk) = zeeu(ji,jj) * ( zls(ji+1,jj  ) - zls(ji,jj) )
               ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
               zsv(ji,jj,jk) = zeev(ji,jj) * ( zls(ji  ,jj+1) - zls(ji,jj) )
            END DO
         END DO

         ! fourth derivative (divergence) and add to the general tracer trend
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! horizontal diffusive trends
               zta = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               zsa = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
               ! add it to the general tracer trends
               ta(ji,jj,jk) = ta(ji,jj,jk) + zta
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
            END DO
         END DO
         !                                             ! ===============
      END DO                                           ! Horizontal slab
      !                                                ! ===============

      ! save the trends for diagnostic
      ! save the horizontal diffusive trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)

         CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      ! "zonal" mean lateral diffusive heat and salt transport
      IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
# if defined key_s_coord || defined key_partial_steps
         pht_ldf(:) = ptr_vj( ztv(:,:,:) )
         pst_ldf(:) = ptr_vj( zsv(:,:,:) )
# else
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztv(ji,jj,jk) = ztv(ji,jj,jk) * fse3v(ji,jj,jk)
                  zsv(ji,jj,jk) = zsv(ji,jj,jk) * fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
         pht_ldf(:) = ptr_vj( ztv(:,:,:) )
         pst_ldf(:) = ptr_vj( zsv(:,:,:) )
# endif
      ENDIF

   END SUBROUTINE tra_ldf_bilap

   !!==============================================================================
END MODULE traldf_bilap
