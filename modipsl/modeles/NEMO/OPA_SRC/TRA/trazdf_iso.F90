MODULE trazdf_iso
   !!==============================================================================
   !!                    ***  MODULE  trazdf_iso  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
#if defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   tra_zdf_iso  : update the tracer trend with the vertical part of 
   !!                  the isopycnal or geopotential s-coord. operator and
   !!                  the vertical diffusion
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE ldfslp          ! Make iso-neutral slopes available
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdf_oce         ! ocean vertical physics 
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE zdfkpp          ! KPP parameterisation
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC tra_zdf_iso    ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "ldfeiv_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso.F90,v 1.13 2005/09/22 15:43:10 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf_iso( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_iso  ***
      !!
      !! ** Purpose :
      !!     Compute the trend due to the vertical tracer diffusion inclu-
      !!     ding the vertical component of lateral mixing (only for second
      !!     order operator, for fourth order it is already computed and
      !!     add to the general trend in traldf.F) and add it to the general
      !!     trend of the tracer equations.
      !!
      !! ** Method :
      !!         The vertical component of the lateral diffusive trends is
      !!      provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      First part: vertical trends associated with the lateral mixing
      !!      ==========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      save avt coef. resulting from vertical physics alone in zavt:
      !!         zavt = avt
      !!      update and save in zavt the vertical eddy viscosity coefficient:
      !!         avt = avt + wslpi^2+wslj^2
      !!      add vertical Eddy Induced advective fluxes ('lk_traldf_eiv=T):
      !!         zftw = zftw + { di[aht e2u mi(wslpi)]
      !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers (t & s) is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme, t=ta.
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!         ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in zdfmxl.F90).
      !!         avt = zavt
      !!         (avs = zavs if lk_zdfddm=T )
      !!
      !!      'key_trdtra' defined: trend saved for futher diagnostics.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action :
      !!         Update (ta,sa) arrays with the before vertical diffusion trend
      !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
      !!
      !! History :
      !!   7.0  !  91-11  (G. Madec)  Original code
      !!        !  92-06  (M. Imbard)  correction on tracer trend loops
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!   9.0  !  05-06  (C. Ethe )  non-local flux in KPP vertical mixing scheme
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce                ,   &
# if defined key_zdfddm
         zavs => va,  &  ! use va as workspace
# endif
         zavt => ua      ! use ua as workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index

      !! * Local save
      REAL(wp), DIMENSION(jpk), SAVE ::  &
         z2dt

      !! * Local declarations
      INTEGER ::   ji, jj, jk                 ! dummy loop indices
      INTEGER ::   ikst, ikenm2, ikstp1       ! temporary integers
#if defined key_partial_steps
      INTEGER ::   iku, ikv, ikv1             ! temporary integers
#endif
      REAL(wp) ::   &
         zcoef0, zcoef3,        &  ! ???
         zcoef4, zavi,          &  ! ???
         zbtr, zmku, zmkv,      &  !
         ztav, zsav
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwd, zws, zwi,         &  ! ???
         zwx, zwy, zwz, zwt        ! ???
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         ztfw, zdit, zdjt, zdj1t,   &
         zsfw, zdis, zdjs, zdj1s
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztavg, zsavg,                      & ! workspace arrays
         ztdta, ztdsa                         ! workspace arrays
#if defined key_traldf_eiv   ||   defined key_esopa
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztfwg, zsfwg 
      REAL(wp) ::                       &
         zcoeg3,                        &
         zuwk, zvwk,                    &
         zuwki, zvwki
#endif
      !!---------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!---------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_zdf_iso : vertical mixing (including isopycnal component)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
#if defined key_diaeiv
         w_eiv(:,:,:) = 0.e0
#endif
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      ztavg(:,:,:) = 0.e0
      zsavg(:,:,:) = 0.e0

      ! time step = 2 rdttra ex
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
      ELSEIF( kt <= nit000 + 1) THEN
         z2dt(:) = 2. * rdttra(:)          ! leapfrog
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

         ! I. vertical trends associated with the lateral mixing
         ! =====================================================
         !  (excluding the vertical flux proportional to dk[t]


         ! I.1 horizontal tracer gradient
         ! ------------------------------

         DO jk = 1, jpkm1
            DO ji = 1, jpim1
               ! i-gradient of T and S at jj
               zdit (ji,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
               zdis (ji,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
               ! j-gradient of T and S at jj
               zdjt (ji,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
               zdjs (ji,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
               ! j-gradient of T and S at jj+1
               zdj1t(ji,jk) = ( tb(ji,jj,jk)-tb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
               zdj1s(ji,jk) = ( sb(ji,jj,jk)-sb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
            END DO
         END DO
#  if defined key_partial_steps
         ! partial steps correction at the bottom ocean level 
         DO ji = 1, jpim1
            ! last ocean level
            iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
            ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
            ikv1 = MIN( mbathy(ji,jj), mbathy(ji  ,jj-1) ) - 1
            ! i-gradient of T and S at jj
            zdit (ji,iku) = gtu(ji,jj)
            zdis (ji,iku) = gsu(ji,jj)
            ! j-gradient of T and S at jj
            zdjt (ji,ikv) = gtv(ji,jj) 
            zdjs (ji,ikv) = gsv(ji,jj) 
            ! j-gradient of T and S at jj+1
            zdj1t(ji,ikv1)= gtv(ji,jj-1)
            zdj1s(ji,ikv1)= gsv(ji,jj-1)
         END DO
#endif


         ! I.2 Vertical fluxes
         ! -------------------

         ! Surface and bottom vertical fluxes set to zero
         ztfw(:, 1 ) = 0.e0
         zsfw(:, 1 ) = 0.e0
         ztfw(:,jpk) = 0.e0
         zsfw(:,jpk) = 0.e0
#if defined key_traldf_eiv
         ztfwg(:,:, 1 ) = 0.e0
         zsfwg(:,:, 1 ) = 0.e0
         ztfwg(:,:,jpk) = 0.e0
         zsfwg(:,:,jpk) = 0.e0
#endif

         ! interior (2=<jk=<jpk-1)
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)

               zmku = 1./MAX( umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)   &
                  &          +umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1. )

               zmkv = 1./MAX( vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)   &
                  &          +vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1. )

               zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
               zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)

               ztfw(ji,jk) = zcoef3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
                  &                    +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
                  &        + zcoef4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
                  &                    +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )

               zsfw(ji,jk) = zcoef3 * ( zdis (ji  ,jk-1) + zdis (ji-1,jk)     &
                  &                    +zdis (ji-1,jk-1) + zdis (ji  ,jk) )   &
                  &        + zcoef4 * ( zdjs (ji  ,jk-1) + zdj1s(ji  ,jk)     &
                  &                    +zdj1s(ji  ,jk-1) + zdjs (ji  ,jk) )
            END DO
         END DO


         ! I.3  update and save of avt (and avs if double diffusive mixing)
         ! ---------------------------

         DO jk = 2, jpkm1
            DO ji = 2, jpim1

               zavi = fsahtw(ji,jj,jk)*( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)   &
                  &                     +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )

               ! save avt in zavt to recover avt for mixed layer depth diag.
               zavt(ji,jj,jk) = avt(ji,jj,jk)
               ! add isopycnal vertical coeff. to avt
               avt(ji,jj,jk) = avt(ji,jj,jk) + zavi
               ! same procedure on avs if necessary
#if defined key_zdfddm
               ! save avs in zavs to recover avs in output files
               zavs(ji,jj,jk) = fsavs(ji,jj,jk)
               ! add isopycnal vertical coeff. to avs
               fsavs(ji,jj,jk) = fsavs(ji,jj,jk) + zavi
#endif
            END DO
         END DO

#if defined key_traldf_eiv
         !                              ! ---------------------------------------!
         !                              ! Eddy induced vertical advective fluxes !
         !                              ! ---------------------------------------!
#if defined key_traldf_c2d || defined key_traldf_c3d
            DO jk = 2, jpkm1
               DO ji = 2, jpim1
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * fsaeiu(ji-1,jj,jk) * e2u(ji-1,jj)*umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * fsaeiu(ji  ,jj,jk) * e2u(ji  ,jj)*umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * fsaeiv(ji,jj-1,jk) * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * fsaeiv(ji,jj  ,jk) * e1v(ji  ,jj)*vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
   
                  ztfwg(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
                  zsfwg(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )
   
                  ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jj,jk)
                  zsfw(ji,jk) = zsfw(ji,jk) + zsfwg(ji,jj,jk)
# if defined key_diaeiv
                  w_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
# endif
               END DO
            END DO

#else
            DO jk = 2, jpkm1
               DO ji = 2, jpim1
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * e2u(ji-1,jj)*umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * e2u(ji  ,jj)*umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
   
                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
                     &            * ( zuwk - zuwki + zvwk - zvwki )

                  ztfwg(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
                  zsfwg(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )

                  ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jj,jk)
                  zsfw(ji,jk) = zsfw(ji,jk) + zsfwg(ji,jj,jk)
# if defined key_diaeiv
                  w_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
# endif
               END DO
            END DO
#endif

#endif

         ! I.5 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------

         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
               ztav = (  ztfw(ji,jk) - ztfw(ji,jk+1)  ) * zbtr
               zsav = (  zsfw(ji,jk) - zsfw(ji,jk+1)  ) * zbtr
               ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the trends for diagnostic
      !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
      IF( l_trdtra )   THEN
#   if defined key_traldf_eiv
         ! Compute the vertical Gent velocity trend
         !                                                ! ===============
         DO jj = 2, jpjm1                                 !  Vertical slab
            !                                             ! ===============
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  ztavg(ji,jj,jk) = ( ztfwg(ji,jj,jk) - ztfwg(ji,jj,jk+1) ) * zbtr
                  zsavg(ji,jj,jk) = ( zsfwg(ji,jj,jk) - zsfwg(ji,jj,jk+1) ) * zbtr
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         CALL trd_mod(ztavg, zsavg, jpttddoe, 'TRA', kt)
#   endif
         ! Recompute the divergence of vertical fluxes ztav & zsav trends 
         ! computed at step 1.5 above in making the difference between the new 
         ! trend ta()/sa() and the previous one ztdta()/ztdsa() and substract 
         ! the vertical Gent velocity trend ztavg()/zsavg() (zero if not used)
         ztavg(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztavg(:,:,:) 
         zsavg(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsavg(:,:,:) 

         ! Save the new ta and sa trends
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

         ! II. Vertical trend associated with the vertical physics
         ! =======================================================
         !     (including the vertical flux proportional to dk[t] associated
         !      with the lateral mixing, through the avt update)
         !     dk[ avt dk[ (t,s) ] ] diffusive trends


         ! II.0 Matrix construction
         ! ------------------------

         ! Diagonal, inferior, superior
         ! (including the bottom boundary condition via avt masked)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwi(ji,jk) = - z2dt(jk) * avt(ji,jj,jk  )   &
                                       / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jk) = - z2dt(jk) * avt(ji,jj,jk+1)   &
                                       / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
            END DO
         END DO

         ! Surface boudary conditions
         DO ji = 2, jpim1
            zwi(ji,1) = 0.e0
            zwd(ji,1) = 1. - zws(ji,1)
         END DO


         ! II.1. Vertical diffusion on t
         ! ---------------------------

         ! Second member construction
#if defined key_zdfkpp
         ! add non-local temperature flux ( in convective case only)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ! zrhs=right hand side 
               zwy(ji,jk) = tb(ji,jj,jk) + z2dt(jk) * ta(ji,jj,jk) &
                  &  - z2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
                  &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
            END DO
         END DO
#else
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwy(ji,jk) = tb(ji,jj,jk) + z2dt(jk) * ta(ji,jj,jk)
            END DO
         END DO
#endif

         ! Matrix inversion from the first level
         ikst = 1
#   include "zdf.matrixsolver.h90"

         ! Save the masked temperature after in ta
         ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
         !                  it will not be done in tranxt)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
            END DO
         END DO


         ! II.2 Vertical diffusion on s
         ! ---------------------------

#if defined key_zdfddm
         ! Rebuild the Matrix as avt /= avs

         ! Diagonal, inferior, superior
         ! (including the bottom boundary condition via avs masked)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwi(ji,jk) = - z2dt(jk) * fsavs(ji,jj,jk  )   &
                                       /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jk) = - z2dt(jk) * fsavs(ji,jj,jk+1)   &
                                       /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
            END DO
         END DO

         ! Surface boudary conditions
         DO ji = 2, jpim1
            zwi(ji,1) = 0.e0
            zwd(ji,1) = 1. - zws(ji,1)
         END DO
#endif
         ! Second member construction

#if defined key_zdfkpp
         ! add non-local temperature flux ( in convective case only)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ! zrhs=right hand side 
               zwy(ji,jk) = sb(ji,jj,jk) + z2dt(jk) * sa(ji,jj,jk) &
                  &  - z2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
                  &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
            END DO
         END DO
#else
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwy(ji,jk) = sb(ji,jj,jk) + z2dt(jk) * sa(ji,jj,jk)
            END DO
         END DO
#endif

         ! Matrix inversion from the first level
         ikst = 1
#   include "zdf.matrixsolver.h90"

         ! Save the masked salinity after in sa
         ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
         !                  it will not be done in tranxt)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
            END DO
         END DO


         ! III. recover the avt (avs) resulting from vertical physics only
         ! ===============================================================

         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               avt(ji,jj,jk) = zavt(ji,jj,jk)
#if defined key_zdfddm
               fsavs(ji,jj,jk) = zavs(ji,jj,jk)
#endif
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the trends for diagnostic
      ! compute the vertical diffusive trends in substracting the previous 
      ! trends ztdta()/ztdsa() to the new one computed via dT/dt or dS/dt 
      ! i.e. with the new temperature/salinity ta/sa computed above
      IF( l_trdtra )   THEN
         IF( ln_traldf_iso)   THEN
            DO jk = 1, jpkm1
               ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk) + ztavg(:,:,jk) 
               ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk) + zsavg(:,:,jk) 
            END DO
         ELSE
            DO jk = 1, jpkm1
               ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk)                             
               ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk)                             
            END DO
         ENDIF

         CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)
      ENDIF

   END SUBROUTINE tra_zdf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module               NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_zdf_iso( kt )              ! empty routine
!      WRITE(*,*) 'tra_zdf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE tra_zdf_iso
#endif

   !!==============================================================================
END MODULE trazdf_iso
