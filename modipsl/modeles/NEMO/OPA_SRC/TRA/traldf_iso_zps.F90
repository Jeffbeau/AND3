MODULE traldf_iso_zps
   !!==============================================================================
   !!                   ***  MODULE  traldf_iso_zps  ***
   !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
   !!==============================================================================
#if ( defined key_ldfslp   &&   defined key_partial_steps )   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'               slope of the lateral diffusive direction
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso_zps : update the tracer trend with the horizontal 
   !!                     component of a iso-neutral laplacian operator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE ldfslp          ! iso-neutral slopes
   USE diaptr          ! poleward transport diagnostics
   USE prtctl          ! Print control


   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC tra_ldf_iso_zps  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "ldfeiv_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_iso_zps.F90,v 1.9 2006/03/20 16:52:22 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_ldf_iso_zps( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso_zps  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
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
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      'key_trdtra' defined: the trend is saved for diagnostics.
      !!
      !!      macro-tasked on horizontal slab (jk-loop).
      !!
      !! ** Action :
      !!         Update (ta,sa) arrays with the before along level biharmonic
      !!      mixing trend.
      !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
      !!
      !! History :
      !!        !  94-08  (G. Madec, M. Imbard)
      !!        !  97-05  (G. Madec)  split into traldf and trazdf
      !!   8.5  !  02-08  (G. Madec)  Free form, F90
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce           , zftu => ua,  &  ! use ua as workspace
         &                zfsu => va      ! use va as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk             ! dummy loop indices
      INTEGER ::   iku, ikv               ! temporary integer
      REAL(wp) ::   &
         zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
         zmsku, zmskv, zbtr, zta, zsa     !    "           "
      REAL(wp), DIMENSION(jpi,jpj) ::   & ! temporary workspace
         zdkt , zdk1t, zdks , zdk1s       !    "           "
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
         zftv, zgtbu, zgtbv,                &  ! temporary workspace
         zfsv, zgsbu, zgsbv,                &  !    "           "
         ztdta, ztdsa
         
#if defined key_traldf_eiv 
      REAL(wp) ::   &
         zcg1, zcg2, zuwk, zvwk,            &  ! temporary scalars
         zuwk1, zvwk1                          !    "           "
      REAL(wp), DIMENSION(jpi,jpj) ::       &  ! temporary workspace
         zftug, zftvg, zfsug, zfsvg            !     "        "     
#endif
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso_zps : iso neutral laplacian diffusion in '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinates with partial steps'
#if defined key_diaeiv
         u_eiv(:,:,:) = 0.e0
         v_eiv(:,:,:) = 0.e0
#endif
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      ! Horizontal temperature and salinity gradient 
      DO jk = 1, jpk
         DO jj = 1, jpj-1
            DO ji = 1, fs_jpim1   ! vector opt.
               zgtbu(ji,jj,jk) = tmask(ji,jj,jk) * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
               zgsbu(ji,jj,jk) = tmask(ji,jj,jk) * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
               zgtbv(ji,jj,jk) = tmask(ji,jj,jk) * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
               zgsbv(ji,jj,jk) = tmask(ji,jj,jk) * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
            END DO
         END DO
      END DO
      ! partial steps correction at the last level 
      DO jj = 1, jpj-1
         DO ji = 1, jpi-1
            ! last level
            iku = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
            ikv = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
            zgtbu(ji,jj,iku) = gtu(ji,jj) 
            zgsbu(ji,jj,iku) = gsu(ji,jj)               
            zgtbv(ji,jj,ikv) = gtv(ji,jj) 
            zgsbv(ji,jj,ikv) = gsv(ji,jj)               
         END DO
      END DO
      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! 1. Vertical tracer gradient at level jk and jk+1
         ! ------------------------------------------------
         ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

         zdk1t(:,:) = ( tb(:,:,jk) - tb(:,:,jk+1) ) * tmask(:,:,jk+1)
         zdk1s(:,:) = ( sb(:,:,jk) - sb(:,:,jk+1) ) * tmask(:,:,jk+1)

         IF( jk == 1 ) THEN
            zdkt(:,:) = zdk1t(:,:)
            zdks(:,:) = zdk1s(:,:)
         ELSE
            zdkt(:,:) = ( tb(:,:,jk-1) - tb(:,:,jk) ) * tmask(:,:,jk)
            zdks(:,:) = ( sb(:,:,jk-1) - sb(:,:,jk) ) * tmask(:,:,jk)
         ENDIF


         ! 2. Horizontal fluxes
         ! --------------------

         DO jj = 1 , jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zabe1 = ( fsahtu(ji,jj,jk) + ahtb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
               zabe2 = ( fsahtv(ji,jj,jk) + ahtb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)

               zmsku = 1. / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
                                + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )

               zmskv = 1. / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
                                + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )

               zcof1 = -fsahtu(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
               zcof2 = -fsahtv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv

               zftu(ji,jj,jk) = umask(ji,jj,jk) * (  zabe1 * zgtbu(ji,jj,jk)   &
                  &                                + zcof1 * (  zdkt (ji+1,jj) + zdk1t(ji,jj)      &
                  &                                           + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
               zftv(ji,jj,jk) = vmask(ji,jj,jk) * (  zabe2 * zgtbv(ji,jj,jk)   &
                  &                                + zcof2 * (  zdkt (ji,jj+1) + zdk1t(ji,jj)      &
                  &                                           + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
               zfsu(ji,jj,jk) = umask(ji,jj,jk) * (  zabe1 * zgsbu(ji,jj,jk)   &
                  &                                + zcof1 * (  zdks (ji+1,jj) + zdk1s(ji,jj)      &
                  &                                           + zdk1s(ji+1,jj) + zdks (ji,jj)  )  )
               zfsv(ji,jj,jk) = vmask(ji,jj,jk) * (  zabe2 * zgsbv(ji,jj,jk)   &
                  &                                + zcof2 * (  zdks (ji,jj+1) + zdk1s(ji,jj)      &
                  &                                           + zdk1s(ji,jj+1) + zdks (ji,jj)  )  )
            END DO
         END DO

#if defined key_traldf_eiv
                                        ! ---------------------------------------!
                                        ! Eddy induced vertical advective fluxes !
                                        ! ---------------------------------------!
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj  ,jk  ) ) * fsaeiu(ji,jj,jk  ) * umask(ji,jj,jk  )
                  zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj  ,jk+1) ) * fsaeiu(ji,jj,jk+1) * umask(ji,jj,jk+1)
                  zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji  ,jj+1,jk  ) ) * fsaeiv(ji,jj,jk  ) * vmask(ji,jj,jk  )
                  zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji  ,jj+1,jk+1) ) * fsaeiv(ji,jj,jk+1) * vmask(ji,jj,jk+1)

                  zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
                  zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )

                  zftug(ji,jj) = zcg1 * ( tb(ji+1,jj,jk) + tb(ji,jj,jk) )
                  zftvg(ji,jj) = zcg2 * ( tb(ji,jj+1,jk) + tb(ji,jj,jk) )
                  zfsug(ji,jj) = zcg1 * ( sb(ji+1,jj,jk) + sb(ji,jj,jk) )
                  zfsvg(ji,jj) = zcg2 * ( sb(ji,jj+1,jk) + sb(ji,jj,jk) )

                  zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
                  zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)
                  zfsu(ji,jj,jk) = zfsu(ji,jj,jk) + zfsug(ji,jj)
                  zfsv(ji,jj,jk) = zfsv(ji,jj,jk) + zfsvg(ji,jj)
#   if defined key_diaeiv
                  u_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
                  v_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
#   endif
               END DO
            END DO
#endif

         ! II.4 Second derivative (divergence) and add to the general trend
         ! ----------------------------------------------------------------

         DO jj = 2 , jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
               zta = zbtr * ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)  )
               zsa = zbtr * ( zfsu(ji,jj,jk) - zfsu(ji-1,jj,jk) + zfsv(ji,jj,jk) - zfsv(ji,jj-1,jk)  )
               ta (ji,jj,jk) = ta (ji,jj,jk) + zta
               sa (ji,jj,jk) = sa (ji,jj,jk) + zsa
            END DO
         END DO
         !                                          ! ===============
      END DO                                        !   End of slab  
      !                                             ! ===============

      ! save the trends for diagnostic
      ! save the horizontal diffusive trends
      IF( l_trdtra )   THEN
#   if defined key_traldf_eiv
         DO jk = 1 , jpkm1
            DO jj = 2 , jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  tladi(ji,jj,jk) = ( zftug(ji,jj) - zftug(ji-1,jj  ) ) * zbtr
                  tladj(ji,jj,jk) = ( zftvg(ji,jj) - zftvg(ji  ,jj-1) ) * zbtr
                  sladi(ji,jj,jk) = ( zfsug(ji,jj) - zfsug(ji-1,jj  ) ) * zbtr
                  sladj(ji,jj,jk) = ( zfsvg(ji,jj) - zfsvg(ji  ,jj-1) ) * zbtr
               END DO
            END DO
         END DO
#   else
         tladi(:,:,:) = 0.e0
         tladj(:,:,:) = 0.e0
         sladi(:,:,:) = 0.e0
         sladj(:,:,:) = 0.e0
#   endif

         ! Substract the eddy induced velocity for T/S
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tladi(:,:,:) - tladj(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sladi(:,:,:) - sladj(:,:,:) 

         CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF


      !!bug  no separation of diff iso and eiv
      IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
         ! "zonal" mean lateral diffusive heat and salt transports
         pht_ldf(:) = ptr_vj( zftv(:,:,:) )
         pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
         ! "zonal" mean lateral eddy induced velocity heat and salt transports
#if defined key_diaeiv
         pht_eiv(:) = ptr_vj( zftv(:,:,:) )
         pst_eiv(:) = ptr_vj( zfsv(:,:,:) )
#endif
      ENDIF

   END SUBROUTINE tra_ldf_iso_zps

#else
   !!----------------------------------------------------------------------
   !!   default option :   Dummy code   NO rotation of the diffusive tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_ldf_iso_zps( kt )               ! Empty routine
!      WRITE(*,*) 'tra_ldf_iso_zps: You should not have seen this print! error?', kt
   END SUBROUTINE tra_ldf_iso_zps
#endif

   !!==============================================================================
END MODULE traldf_iso_zps
