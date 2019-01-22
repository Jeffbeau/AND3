MODULE trcldf_iso
   !!==============================================================================
   !!                    ***  MODULE  trcldf_iso  ***
   !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
#if key_passivetrc && defined key_ldfslp 
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   trc_ldf_iso : update the tracer trend with the horizontal component
   !!                 of iso neutral laplacian operator or horizontal 
   !!                 laplacian operator in s-coordinate
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc      ! ocean dynamics and tracers variables
   USE trc          ! ocean passive tracers variables
   USE prtctl_trc   ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_ldf_iso  ! routine called by step.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_iso.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_ldf_iso( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ldf_iso  ***
      !! 
      !! ** Purpose :   Compute the before horizontal tracer  diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The horizontal component of the lateral diffusive trends 
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      horizontal fluxes associated with the rotated lateral mixing:
      !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
      !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
      !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
      !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
      !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
      !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( tb ) 
      !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( tb ) 
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
      !!      Add this trend to the general trend tra :
      !!         tra = tra + difft
      !!
      !! ** Action  : - Update tra arrays with the before isopycnal or
      !!                geopotential s-coord harmonic mixing trend.
      !!              - Save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!        !  94-08  (G. Madec, M. Imbard)
      !!        !  97-05  (G. Madec)  split into traldf and trazdf
      !!        !  98-03  (L. Bopp, MA Foujols) passive tracer generalisation
      !!        !  00-10  (MA Foujols E Kestenare) USE passive tracer coefficient
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!   9.0  !  04-03  (C. Ethe)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce_trc       , zftu => ua,  &  ! use ua as workspace
         &                zfsu => va      ! use va as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn             ! dummy loop indices
      REAL(wp) ::   &
         zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
         zmsku, zmskv, zbtr,           &
#if defined key_trcldf_eiv
         zcg1, zcg2, zuwk, zvwk,       &
         zuwk1, zvwk1,                 &
#endif
         ztra

      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zdkt, zdk1t            ! workspace

#if defined key_trcldf_eiv
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zftug, zftvg
#endif

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zftv                       ! workspace
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_ldf_iso : iso neutral lateral diffusion or'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   horizontal laplacian diffusion in s-coordinate'
#if defined key_trcldf_eiv && defined key_diaeiv
         u_trc_eiv(:,:,:) = 0.e0
         v_trc_eiv(:,:,:) = 0.e0
#endif
      ENDIF


      DO jn = 1, jptra

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            ! 1. Vertical tracer gradient at level jk and jk+1
            ! ------------------------------------------------
            ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

            zdk1t(:,:) = ( trb(:,:,jk,jn) - trb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)

            IF( jk == 1 ) THEN
               zdkt(:,:) = zdk1t(:,:)
            ELSE
               zdkt(:,:) = ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF


            ! 2. Horizontal fluxes
            ! --------------------

            DO jj = 1 , jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zabe1 = ( fsahtru(ji,jj,jk) + ahtrb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
                  zabe2 = ( fsahtrv(ji,jj,jk) + ahtrb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)

                  zmsku = 1. / MAX(   tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
                     + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )

                  zmskv = 1. / MAX(   tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
                     + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )

                  zcof1 = -fsahtru(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
                  zcof2 = -fsahtrv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv

                  zftu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn)  )   &
                     &                              + zcof1 * (   zdkt (ji+1,jj) + zdk1t(ji,jj)      &
                     &                                          + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )

                  zftv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn)  )   &
                     &                              + zcof2 * (   zdkt (ji,jj+1) + zdk1t(ji,jj)      &
                     &                                          + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )

               END DO
            END DO

#   if defined key_trcldf_eiv
            !                              ! ---------------------------------------!
            !                              ! Eddy induced vertical advective fluxes !
            !                              ! ---------------------------------------!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk  ) ) * fsaeitru(ji,jj,jk  ) * umask(ji,jj,jk  )
                  zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) ) * fsaeitru(ji,jj,jk+1) * umask(ji,jj,jk+1)
                  zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk  ) ) * fsaeitrv(ji,jj,jk  ) * vmask(ji,jj,jk  )
                  zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) ) * fsaeitrv(ji,jj,jk+1) * vmask(ji,jj,jk+1)

                  zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
                  zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )

                  zftug(ji,jj) = zcg1 * ( trb(ji+1,jj,jk,jn) + trb(ji,jj,jk,jn) )
                  zftvg(ji,jj) = zcg2 * ( trb(ji,jj+1,jk,jn) + trb(ji,jj,jk,jn) )

                  zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
                  zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)

#   if defined key_diaeiv
                  u_trc_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
                  v_trc_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
#   endif
               END DO
            END DO
#   endif

            ! II.4 Second derivative (divergence) and add to the general trend
            ! ----------------------------------------------------------------

            DO jj = 2 , jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  ztra = zbtr * (  zftu(ji,jj,jk) - zftu(ji-1,jj  ,jk)   &
                     &          + zftv(ji,jj,jk) - zftv(ji  ,jj-1,jk)  )
                  tra (ji,jj,jk,jn) = tra (ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd
                  IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),4) = ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk  ) ) * zbtr
                  IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),5) = ( zftv(ji,jj,jk) - zftv(ji,jj-1,jk  ) ) * zbtr
#endif
               END DO
            END DO
            !                                          ! ===============
         END DO                                        !   End of slab  
         !                                             ! ===============

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('ldf - iso')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_ldf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :             No rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_ldf_iso( kt )               ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_ldf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE trc_ldf_iso
#endif

   !!==============================================================================
END MODULE trcldf_iso
