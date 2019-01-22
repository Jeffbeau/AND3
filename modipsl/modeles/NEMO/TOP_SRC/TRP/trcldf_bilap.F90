MODULE trcldf_bilap
   !!==============================================================================
   !!                   ***  MODULE  trcldf_bilap  ***
   !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_ldf_bilap : update the tracer trend with the horizontal diffusion
   !!                   using a iso-level biharmonic operator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_ldf_bilap   ! routine called by step.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_bilap.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE trc_ldf_bilap( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer tra diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces 
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends of passive tracer is given by:
      !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
      !!      factors e3. are inside the derivatives:
      !!      Laplacian of trb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(trb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(trb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!      Laplacian of trb:
      !!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(trb) ]
      !!                              + dj-1[ e1v/e2v dj(trb) ] }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
      !!                              + dj-1[ e1v/e2v dj(zlt) ]  }
      !!
      !!      Add this trend to the general trend tra :
      !!         tra = tra + difft 
      !!
      !! ** Action : - Update tra arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!             - Save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!        !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  96-01  (M. Imbard)  mpp exchange
      !!        !  97-07  (G. Madec)  optimization, and ahtt
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-10  (MA Foujols E. Kestenare) use passive tracer coefficient
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
#if defined key_partial_steps
      INTEGER ::   iku, ikv                   ! temporary integers
#endif
      REAL(wp) ::   ztra     ! temporary scalars

      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zeeu, zeev, zbtr, zlt                 ! workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
         ztu, ztv                              ! workspace
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_ldf_bilap : iso-level biharmonic operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      ! 

      DO jn = 1, jptra
                                                          ! ===============
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
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( trb(ji+1,jj  ,jk,jn) - trb(ji,jj,jk,jn) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( trb(ji  ,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
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
                  ENDIF
                  IF( ikv == jk ) THEN
                     ztv(ji,jj,jk) = zeev(ji,jj) * gtv(ji,jj)
                  ENDIF
               END DO
            END DO
#endif

            ! Second derivative (divergence)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlt(ji,jj) = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
               END DO
            END DO

            ! Multiply by the eddy diffusivity coefficient
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zlt(ji,jj) = fsahtrt(ji,jj,jk) * zlt(ji,jj)
               END DO
            END DO

            ! Lateral boundary conditions on the laplacian zlt   (unchanged sgn)
            CALL lbc_lnk( zlt, 'T', 1. )  

            ! 2. Bilaplacian
            ! --------------

            ! third derivative (gradient)
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
                  ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
               END DO
            END DO

            ! fourth derivative (divergence) and add to the general tracer trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! horizontal diffusive trends
                  ztra = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
                  ! save the horizontal diffusive trends
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),4) = (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) ) * zbtr(ji,jj)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),5) = (  ztv(ji,jj,jk) - ztv(ji-1,jj,jk) ) * zbtr(ji,jj)
#endif
               END DO
            END DO
            !                                             ! ===============
         END DO                                           ! Horizontal slab
         !                                                ! ===============
#if defined key_trc_diatrd
         ! Lateral boundary conditions on the laplacian zlt   (unchanged sgn)
         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),5), 'T', 1. )  
#endif
      END DO

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ldf - bilap')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_ldf_bilap

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf_bilap( kt )  
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_ldf_bilap: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf_bilap
#endif
   !!==============================================================================
END MODULE trcldf_bilap
