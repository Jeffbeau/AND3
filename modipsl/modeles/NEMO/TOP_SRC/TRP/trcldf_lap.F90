MODULE trcldf_lap
   !!==============================================================================
   !!                       ***  MODULE  trcldf_lap  ***
   !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_ldf_lap  : update the tracer trend with the horizontal diffusion
   !!                 using a iso-level harmonic (laplacien) operator.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_ldf_lap  ! routine called by step.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_lap.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE trc_ldf_lap( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ldf_lap  ***
      !!                   
      !! ** Purpose :   Compute the before horizontal tracer diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of 
      !!      the passive tracer is given by:
      !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
      !!      factors e3. are inside the derivatives:
      !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(trb) ]
      !!                                   + dj-1[ aht e1v*e3v/e2v dj(trb) ] }
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(trb) ]
      !!                               + dj-1[ aht e1v/e2v dj(trb) ] }
      !!      Add this trend to the general tracer trend tra :
      !!          tra = tra + difft
      !!
      !! ** Action  : - Update tra arrays with the before iso-level 
      !!                harmonic mixing trend.
      !!              - Save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  95-02  (M. Levy)    passive tracers
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)   passive tracer
      !!----------------------------------------------------------------------
      USE oce_trc          , ztu => ua,  &  ! use ua as workspace
         &                   ztv => va      ! use va as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      
      !! * Local save
      REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
         ze1ur, ze2vr, zbtr2              ! scale factor coefficients
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn         ! dummy loop indices
      REAL(wp) ::   &
         zabe1, zabe2, zbtr              ! temporary scalars

      REAL(wp) ::   &
         ztra, ztrax, ztray              ! workspace
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------
      
      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_ldf_lap : iso-level laplacian diffusion'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         ze1ur(:,:) = e2u(:,:) / e1u(:,:)
         ze2vr(:,:) = e1v(:,:) / e2v(:,:)
         zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
      ENDIF
 
      DO jn = 1, jptra
         
         !                                                  ! =============
         DO jk = 1, jpkm1                                   ! Vertical slab
            !                                               ! =============
            ! 1. First derivative (gradient)
            ! -------------------
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
#if defined key_s_coord
                  zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj) * fse3u(ji,jj,jk)
                  zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj) * fse3v(ji,jj,jk)
#else
                  zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj)
                  zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj)
#endif
                  ztu(ji,jj,jk) = zabe1 * ( trb(ji+1,jj  ,jk,jn) - trb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zabe2 * ( trb(ji  ,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
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
                  ztrax = zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )
                  ztray = zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )

                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztrax + ztray

#if defined key_trc_diatrd
                  ! save the horizontal diffusive trends
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),4) = ztrax
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),5) = ztray
#endif
               END DO
            END DO
            !                                               ! =============
         END DO                                             !  End of slab  
         !                                                  ! =============

      END DO

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ldf - lap')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_ldf_lap

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf_lap( kt )  
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_ldf_lap: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf_lap
#endif

   !!==============================================================================
END MODULE trcldf_lap
