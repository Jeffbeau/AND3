MODULE trcldf_bilapg
   !!==============================================================================
   !!                       ***  MODULE  trcldf_bilapg  ***
   !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
#if key_passivetrc && defined key_ldfslp
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   trc_ldf_bilapg : update the tracer trend with the horizontal diffusion
   !!                    using an horizontal biharmonic operator in s-coordinate 
   !!   ldfght         :  ???
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and tracers variables
   USE trc                 ! ocean passive tracers variables
   USE lbclnk              ! ocean lateral boundary condition (or mpp link)
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_ldf_bilapg    ! routine called by step.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_bilapg.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE trc_ldf_bilapg( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_ldf_bilapg  ***
      !!                    
      !! ** Purpose :   Compute the before horizontal passive tracer diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The lateral diffusive trends is provided by a 4th order
      !!      operator rotated along geopotential surfaces. It is computed
      !!      using before fields (forward in time) and geopotential slopes
      !!      computed in routine inildf.
      !!         -1- compute the geopotential harmonic operator applied to
      !!      trb and multiply it by the eddy diffusivity coefficient
      !!      (done by a call to ldfght routine, result in wk1 array).
      !!      Applied the domain lateral boundary conditions by call to lbc_lnk
      !!         -2- compute the geopotential harmonic operator applied to
      !!      wk1 by a second call to ldfght routine (result in wk3 array).
      !!         -3- Add this trend to the general trend (ta,sa):
      !!            tra = tra + wk3
      !!
      !! ** Action : - Update tra arrays with the before geopotential 
      !!               biharmonic mixing trend.
      !!             - Save the trends  in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   8.0  !  97-07  (G. Madec)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn              ! dummy loop indices
      REAL(wp) ::  ztra                       ! workspace   
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   &
         wk1, wk2               ! work array used for rotated biharmonic
                                ! operator on tracers and/or momentum
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF

      ! 1. Laplacian of passive tracers trb * aht
      ! ----------------------------- 
      ! rotated harmonic operator applied to trb 
      ! and multiply by aht (output in wk1 )

      CALL ldfght ( trb, wk1, 1 )

      DO jn = 1, jptra
      ! Lateral boundary conditions on wk1   (unchanged sign)
         CALL lbc_lnk( wk1(:,:,:,jn) , 'T', 1. )
      END DO

      ! 2. Bilaplacian of trb
      ! -------------------------
      ! rotated harmonic operator applied to wk1 (output in wk2 )

      CALL ldfght ( wk1, wk2, 2 )


      DO jn = 1, jptra
         ! 3. Update the tracer trends                    (j-slab :   2, jpj-1)
         ! ---------------------------
         !                                                ! ===============
         DO jj = 2, jpjm1                                 !  Vertical slab
            !                                             ! ===============
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + wk2(ji,jj,jk,jn)
#if defined key_trc_diatrd
                  ! save the horizontal diffusive trends
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = wk2(ji,jj,jk,jn)
#endif
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      END DO

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ldf - bilapg')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_ldf_bilapg


   SUBROUTINE ldfght ( pt, plt, kaht )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldfght  ***
      !!          
      !! ** Purpose :   Apply a geopotential harmonic operator to pt and 
      !!      multiply it by the eddy diffusivity coefficient (if kaht=1).
      !!      Routine only used in s-coordinates (l_sco=T) with bilaplacian
      !!      operator (ln_traldf_bilap=T) acting along geopotential surfaces
      !!      (ln_traldf_hor).
      !!
      !! ** Method  :   The harmonic operator rotated along geopotential 
      !!      surfaces is applied to pt using the slopes of geopotential
      !!      surfaces computed in inildf routine. The result is provided in
      !!      plt arrays. It is computed in 2 steps:
      !!
      !!      First step: horizontal part of the operator. It is computed on
      !!      ==========  pt as follows (idem on ps)
      !!      horizontal fluxes :
      !!         zftu = e2u*e3u/e1u di[ pt ] - e2u*uslp dk[ mi(mk(pt)) ]
      !!         zftv = e1v*e3v/e2v dj[ pt ] - e1v*vslp dk[ mj(mk(pt)) ]
      !!      take the horizontal divergence of the fluxes (no divided by
      !!      the volume element :
      !!         plt  = di-1[ zftu ] +  dj-1[ zftv ]
      !!
      !!      Second step: vertical part of the operator. It is computed on
      !!      ===========  pt as follows (idem on ps)
      !!      vertical fluxes :
      !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pt ]
      !!              -     e2t     *       wslpi        di[ mi(mk(pt)) ]
      !!              -     e1t     *       wslpj        dj[ mj(mk(pt)) ]
      !!      take the vertical divergence of the fluxes add it to the hori-
      !!      zontal component, divide the result by the volume element and
      !!      if kaht=1, multiply by the eddy diffusivity coefficient:
      !!         plt  = aht / (e1t*e2t*e3t) { plt + dk[ zftw ] }
      !!      else:
      !!         plt  =  1  / (e1t*e2t*e3t) { plt + dk[ zftw ] }
      !!
      !! * Action :
      !!      'key_trdtra' defined: the trend is saved for diagnostics.
      !!
      !! History :
      !!   8.0  !  97-07  (G. Madec)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( in  ) ::   &
         pt               ! tracer fields before for 1st call
      !                   ! and laplacian of these fields for 2nd call.
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( out ) ::   &
         plt              ! partial harmonic operator applied to
      !                   ! pt components except
      !                   ! second order vertical derivative term)
      INTEGER, INTENT( in ) ::   &
         kaht             ! =1 multiply the laplacian by the eddy diffusivity coeff.
      !                   ! =2 no multiplication

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn             ! dummy loop indices
      REAL(wp) ::   &
         zabe1, zabe2, zmku, zmkv,     &  ! temporary scalars
         zbtr, ztah, ztav, &
         zcof0, zcof1, zcof2,          &
         zcof3, zcof4

      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zftu, zftv ,                   &  ! workspace
         zdkt, zdk1t

      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zftw,                          &  ! workspace
         zdit, zdjt, zdj1t

      !!----------------------------------------------------------------------


      DO jn = 1, jptra
         !                               ! ********** !   ! ===============
         DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
            !                            ! ********** !   ! ===============

            ! I.1 Vertical gradient of pt and ps at level jk and jk+1
            ! -------------------------------------------------------
            !     surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

            zdk1t(:,:) = ( pt(:,:,jk,jn) - pt(:,:,jk+1,jn) ) * tmask(:,:,jk+1)

            IF( jk == 1 ) THEN
               zdkt(:,:) = zdk1t(:,:)
            ELSE
               zdkt(:,:) = ( pt(:,:,jk-1,jn) - pt(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF


            ! I.2 Horizontal fluxes
            ! ---------------------

            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  zabe1 = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
                  zabe2 = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)

                  zmku=1./MAX( tmask(ji+1,jj,jk  )+tmask(ji,jj,jk+1)   &
                              +tmask(ji+1,jj,jk+1)+tmask(ji,jj,jk  ),1. )
                  zmkv=1./MAX( tmask(ji,jj+1,jk  )+tmask(ji,jj,jk+1)   &
                              +tmask(ji,jj+1,jk+1)+tmask(ji,jj,jk  ),1. )

                  zcof1= -e2u(ji,jj) * uslp(ji,jj,jk) * zmku
                  zcof2= -e1v(ji,jj) * vslp(ji,jj,jk) * zmkv

                  zftu(ji,jj)= umask(ji,jj,jk) *   &
                     (  zabe1 *( pt(ji+1,jj,jk,jn) - pt(ji,jj,jk,jn) )   &
                      + zcof1 *( zdkt (ji+1,jj) + zdk1t(ji,jj)           &
                                +zdk1t(ji+1,jj) + zdkt (ji,jj) )  )

                  zftv(ji,jj)= vmask(ji,jj,jk) *   &
                     (  zabe2 *( pt(ji,jj+1,jk,jn) - pt(ji,jj,jk,jn) )   &
                      + zcof2 *( zdkt (ji,jj+1) + zdk1t(ji,jj)           &
                                +zdk1t(ji,jj+1) + zdkt (ji,jj) )  )
               END DO
            END DO


            ! I.3 Second derivative (divergence) (not divided by the volume)
            ! ---------------------

            DO jj = 2 , jpjm1
               DO ji = 2 , jpim1
                  ztah = zftu(ji,jj) - zftu(ji-1,jj) + zftv(ji,jj) - zftv(ji,jj-1)
                  plt(ji,jj,jk,jn) = ztah
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


         !                             ! ************ !   ! ===============
         DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
            !                          ! ************ !   ! ===============

            ! II.1 horizontal tracer gradient
            ! -------------------------------

            DO jk = 1, jpk
               DO ji = 1, jpim1
                  zdit (ji,jk) = ( pt(ji+1,jj  ,jk,jn) - pt(ji,jj  ,jk,jn) ) * umask(ji,jj  ,jk)
                  zdjt (ji,jk) = ( pt(ji  ,jj+1,jk,jn) - pt(ji,jj  ,jk,jn) ) * vmask(ji,jj  ,jk)
                  zdj1t(ji,jk) = ( pt(ji  ,jj  ,jk,jn) - pt(ji,jj-1,jk,jn) ) * vmask(ji,jj-1,jk)
               END DO
            END DO


            ! II.2 Vertical fluxes
            ! --------------------

            ! Surface and bottom vertical fluxes set to zero
            zftw(:, 1 ) = 0.e0
            zftw(:,jpk) = 0.e0

            ! interior (2=<jk=<jpk-1)
            DO jk = 2, jpkm1
               DO ji = 2, jpim1
                  zcof0 = e1t(ji,jj) * e2t(ji,jj) / fse3w(ji,jj,jk)   &
                     * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)        &
                        + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )

                  zmku =1./MAX(  umask(ji  ,jj,jk-1)+umask(ji-1,jj,jk)   &
                                +umask(ji-1,jj,jk-1)+umask(ji  ,jj,jk), 1. )

                  zmkv =1./MAX(  vmask(ji,jj  ,jk-1)+vmask(ji,jj-1,jk)   &
                                +vmask(ji,jj-1,jk-1)+vmask(ji,jj  ,jk), 1. )

                  zcof3 = - e2t(ji,jj) * wslpi (ji,jj,jk) * zmku
                  zcof4 = - e1t(ji,jj) * wslpj (ji,jj,jk) * zmkv

                  zftw(ji,jk) = tmask(ji,jj,jk) *   &
                     (  zcof0 * ( pt  (ji,jj,jk-1,jn) - pt  (ji,jj,jk,jn) )   &
                      + zcof3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)           &
                                 +zdit (ji-1,jk-1) + zdit (ji  ,jk) )         &
                      + zcof4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)           &
                                 +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )  )
               END DO
            END DO


            ! II.3 Divergence of vertical fluxes added to the horizontal divergence
            ! ---------------------------------------------------------------------

            IF( kaht == 1 ) THEN
               ! multiply the laplacian by the eddy diffusivity coefficient
               DO jk = 1, jpkm1
                  DO ji = 2, jpim1
                     ! eddy coef. divided by the volume element
                     zbtr = fsahtrt(ji,jj,jk) / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                     ! vertical divergence
                     ztav = zftw(ji,jk) - zftw(ji,jk+1)
                     ! harmonic operator applied to (pt,ps) and multiply by aht
                     plt(ji,jj,jk,jn) = ( plt(ji,jj,jk,jn) + ztav ) * zbtr
                  END DO
               END DO
            ELSEIF( kaht == 2 ) THEN
               ! second call, no multiplication
               DO jk = 1, jpkm1
                  DO ji = 2, jpim1
                     ! inverse of the volume element
                     zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                     ! vertical divergence
                     ztav = zftw(ji,jk) - zftw(ji,jk+1)
                     ! harmonic operator applied to (pt,ps) 
                     plt(ji,jj,jk,jn) = ( plt(ji,jj,jk,jn) + ztav ) * zbtr
                  END DO
               END DO
            ELSE
               IF(lwp) WRITE(numout,*) ' ldfght: kaht= 1 or 2, here =', kaht
               IF(lwp) WRITE(numout,*) '         We stop'
               STOP 'ldfght'
            ENDIF
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      END DO

   END SUBROUTINE ldfght

#else  
   !!----------------------------------------------------------------------
   !!   Dummy module :             NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf_bilapg( kt )               ! Dummy routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_ldf_bilapg: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf_bilapg
#endif

   !!==============================================================================
END MODULE trcldf_bilapg
