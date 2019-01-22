MODULE traldf_lap
   !!==============================================================================
   !!                       ***  MODULE  traldf_lap  ***
   !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_ldf_lap  : update the tracer trend with the horizontal diffusion
   !!                 using a iso-level harmonic (laplacien) operator.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE diaptr          ! poleward transport diagnostics
   USE prtctl          ! Print control


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_ldf_lap  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_lap.F90,v 1.7 2005/09/02 15:45:33 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE tra_ldf_lap( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_lap  ***
      !!                   
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of 
      !!      temperature (idem for salinity) is given by:
      !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
      !!      factors e3. are inside the derivatives:
      !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
      !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(tb) ]
      !!                               + dj-1[ aht e1v/e2v dj(tb) ] }
      !!      Add this trend to the general tracer trend (ta,sa):
      !!          (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action  : - Update (ta,sa) arrays with the before iso-level 
      !!                harmonic mixing trend.
      !!              - Save the trends in (ztdta,ztdsa) ('key_trdtra')
      !!
      !! History :
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      USE oce              , ztu => ua,  &  ! use ua as workspace
         &                   zsu => va      ! use va as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      
      !! * Local save
      REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
         ze1ur, ze2vr, zbtr2              ! scale factor coefficients
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   &
         zabe1, zabe2, zbtr, zta, zsa     ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztv, zsv,                      &  ! temporary workspace arrays
         ztdta, ztdsa                      !    "         "
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_lap : iso-level laplacian diffusion'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         ze1ur(:,:) = e2u(:,:) / e1u(:,:)
         ze2vr(:,:) = e1v(:,:) / e2v(:,:)
         zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
      ENDIF
      
      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      !                                                  ! =============
      DO jk = 1, jpkm1                                   ! Vertical slab
         !                                               ! =============
         ! 1. First derivative (gradient)
         ! -------------------
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
#if defined key_s_coord
               zabe1 = fsahtu(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj) * fse3u(ji,jj,jk)
               zabe2 = fsahtv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj) * fse3v(ji,jj,jk)
#else
               zabe1 = fsahtu(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj)
               zabe2 = fsahtv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj)
#endif
               ztu(ji,jj,jk) = zabe1 * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
               zsu(ji,jj,jk) = zabe1 * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
               ztv(ji,jj,jk) = zabe2 * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
               zsv(ji,jj,jk) = zabe2 * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
            END DO  
         END DO  
         
         
         ! 2. Second derivative (divergence)
         ! --------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord
               zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
#else
               zbtr = zbtr2(ji,jj)
#endif
               ! horizontal diffusive trends
               zta = zbtr * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk)   &
                  &          + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               zsa = zbtr * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk)   &
                  &          + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
               ! add it to the general tracer trends
               ta(ji,jj,jk) = ta(ji,jj,jk) + zta
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
            END DO  
         END DO  
         !                                               ! =============
      END DO                                             !  End of slab  
      !                                                  ! =============

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

   END SUBROUTINE tra_ldf_lap

   !!==============================================================================
END MODULE traldf_lap
