MODULE traadv_cen2
   !!==============================================================================
   !!                       ***  MODULE  traadv_cen2  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!==============================================================================
 
   !!----------------------------------------------------------------------
   !!   tra_adv_cen2 : update the tracer trend with the horizontal
   !!                  and vertical advection trends using a 2nd order 
   !!                  centered finite difference scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE flxrnf          !
   USE trabbl          ! advective term in the BBL
   USE ocfzpt          !
   USE lib_mpp
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)
   USE in_out_manager  ! I/O manager
   USE diaptr          ! poleward transport diagnostics
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC tra_adv_cen2    ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_cen2.F90,v 1.11 2005/12/28 09:25:10 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking' :      2nd order centered scheme (k- and j-slabs)
   !!----------------------------------------------------------------------
#  include "traadv_cen2_atsk.h90"

#else
   !!----------------------------------------------------------------------
   !!   Default option :             2nd order centered scheme (k-j-i loop)
   !!----------------------------------------------------------------------

   SUBROUTINE tra_adv_cen2( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_cen2  ***
      !!                 
      !! ** Purpose :   Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The advection is evaluated by a second order centered
      !!      scheme using now fields (leap-frog scheme). In specific areas
      !!      (vicinity of major river mouths, some straits, or where tn is
      !!      is approaching the freezing point) it is mixed with an upstream
      !!      scheme for stability reasons.
      !!        Part 0 : compute the upstream / centered flag
      !!                 (3D array, zind, defined at T-point (0<zind<1))
      !!        Part I : horizontal advection
      !!      * centered flux:
      !!         * s-coordinate (lk_sco=T) or
      !!         * z-coordinate with partial steps (lk_zps=T),
      !!        the vertical scale factors e3. are inside the derivatives:
      !!               zcenu = e2u*e3u  un  mi(tn)
      !!               zcenv = e1v*e3v  vn  mj(tn)
      !!         * z-coordinate (default key), e3t=e3u=e3v:
      !!               zcenu = e2u  un  mi(tn)
      !!               zcenv = e1v  vn  mj(tn)
      !!      * upstream flux:
      !!         * s-coordinate (lk_sco=T) or
      !!         * z-coordinate with partial steps (lk_zps=T)
      !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
      !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
      !!         * z-coordinate (default key)
      !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
      !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
      !!      * mixed upstream / centered horizontal advection scheme
      !!               zcofi = max(zind(i+1), zind(i))
      !!               zcofj = max(zind(j+1), zind(j))
      !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
      !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
      !!      * horizontal advective trend (divergence of the fluxes)
      !!         * s-coordinate (lk_sco=T) or
      !!         * z-coordinate with partial steps (lk_zps=T)
      !!               zta = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
      !!         * z-coordinate (default key), e3t=e3u=e3v:
      !!               zta = 1/(e1t*e2t) { di-1[zwx] + dj-1[zwy] }
      !!      * Add this trend now to the general trend of tracer (ta,sa):
      !!              (ta,sa) = (ta,sa) + ( zta , zsa )
      !!      * trend diagnostic ('key_trdtra'): the trend is saved
      !!      for diagnostics. The trends saved is expressed as
      !!      Uh.gradh(T), (save trend = zta + tn divn).
      !!         In addition, the advective trend in the two horizontal direc-
      !!      tion is also re-computed as Uh gradh(T). Indeed hadt+tn divn is
      !!      equal to (in s-coordinates, and similarly in z-coord.):
      !!         zta+tn*divn=1/(e1t*e2t*e3t) { mi-1( e2u*e3u  un  di[tn] )
      !!                                      +mj-1( e1v*e3v  vn  mj[tn] )  }
      !!         C A U T I O N : the trend saved is the centered trend only.
      !!      It doesn't take into account the upstream part of the scheme.
      !!
      !!         Part II : vertical advection
      !!      For temperature (idem for salinity) the advective trend is com-
      !!      puted as follows :
      !!            zta = 1/e3t dk+1[ zwz ]
      !!      where the vertical advective flux, zwz, is given by :
      !!            zwz = zcofk * zupst + (1-zcofk) * zcent
      !!      with 
      !!        zupsv = upstream flux = wn * (tb(k) or tb(k-1) ) [wn>0 or <0]
      !!        zcenu = centered flux = wn * mk(tn)
      !!         The surface boundary condition is : 
      !!      rigid-lid (default option) : zero advective flux
      !!      free-surf ("key_fresurf_cstvol") : wn(:,:,1) * tn(:,:,1)
      !!         Add this trend now to the general trend of tracer (ta,sa):
      !!            (ta,sa) = (ta,sa) + ( zta , zsa )
      !!         Trend diagnostic ('key_trdtra'): the trend is saved for
      !!      diagnostics. The trends saved is expressed as :
      !!             save trend =  w.gradz(T) = zta - tn divn.
      !!
      !! ** Action : - update (ta,sa) with the now advective tracer trends
      !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
      !!
      !! History :
      !!   8.2  !  01-08  (G. Madec, E. Durand)  trahad+trazad = traadv 
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce                , zwx => ua,  &  ! use ua as workspace
         &                     zwy => va      ! use va as workspace
#if defined key_trabbl_adv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
         &         zun, zvn, zwn
#else
      USE oce                , zun => un,  &  ! When no bbl, zun == un
         &                     zvn => vn,  &  ! When no bbl, zvn == vn
         &                     zwn => wn      ! When no bbl, zwn == wn
#endif

 
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
 
      !! * Local save
      REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
         zbtr2
 
      !! * Local declarations
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::                           &
         zbtr, zta, zsa, zfui, zfvj,        &  ! temporary scalars
         zhw, ze3tr, zcofi, zcofj,          &  !    "         "
         zupsut, zupsvt, zupsus, zupsvs,    &  !    "         "
         zfp_ui, zfp_vj, zfm_ui, zfm_vj,    &  !    "         "
         zcofk, zupst, zupss, zcent,        &  !    "         "
         zcens, zfp_w, zfm_w,               &  !    "         "
         zcenut, zcenvt, zcenus, zcenvs        !    "         "
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwz, zww, zind,                    &  ! temporary workspace arrays
         ztdta, ztdsa                          !    "         "
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_cen2 : 2nd order centered advection scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
         IF(lwp) WRITE(numout,*)
   
         zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
         l_adv = 'ce2'
      ENDIF

#if defined key_trabbl_adv
      ! Advective bottom boundary layer 
      ! -------------------------------
      zun(:,:,:) = un(:,:,:) - u_bbl(:,:,:)
      zvn(:,:,:) = vn(:,:,:) - v_bbl(:,:,:)
      zwn(:,:,:) = wn(:,:,:) + w_bbl(:,:,:)
#endif

      ! Upstream / centered scheme indicator
      ! ------------------------------------
 
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zind(ji,jj,jk) =  MAX ( upsrnfh(ji,jj) * upsrnfz(jk),     &  ! changing advection scheme near runoff
                  &                    upsadv(ji,jj)                     &  ! in the vicinity of some straits
#if defined key_ice_lim
                  &                  , tmask(ji,jj,jk)                   &  ! half upstream tracer fluxes
                  &                  * MAX( 0., SIGN( 1., fzptn(ji,jj)   &  ! if tn < ("freezing"+0.1 )
                  &                                +0.1-tn(ji,jj,jk) ) ) &
#endif
                  &                  )
            END DO
         END DO
      END DO


      ! I. Horizontal advective fluxes
      ! ------------------------------

      ! 1. Second order centered tracer flux at u and v-points
      ! -------------------------------------------------------

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ! upstream indicator
               zcofi = MAX( zind(ji+1,jj,jk), zind(ji,jj,jk) )
               zcofj = MAX( zind(ji,jj+1,jk), zind(ji,jj,jk) )
               ! volume fluxes * 1/2
#if defined key_s_coord || defined key_partial_steps
               zfui = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
               zfvj = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
#else
               zfui = 0.5 * e2u(ji,jj) * zun(ji,jj,jk)
               zfvj = 0.5 * e1v(ji,jj) * zvn(ji,jj,jk)
#endif
               ! upstream scheme
               zfp_ui = zfui + ABS( zfui )
               zfp_vj = zfvj + ABS( zfvj )
               zfm_ui = zfui - ABS( zfui )
               zfm_vj = zfvj - ABS( zfvj )
               zupsut = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
               zupsvt = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
               zupsus = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
               zupsvs = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
               ! centered scheme
               zcenut = zfui * ( tn(ji,jj,jk) + tn(ji+1,jj  ,jk) )
               zcenvt = zfvj * ( tn(ji,jj,jk) + tn(ji  ,jj+1,jk) )
               zcenus = zfui * ( sn(ji,jj,jk) + sn(ji+1,jj  ,jk) )
               zcenvs = zfvj * ( sn(ji,jj,jk) + sn(ji  ,jj+1,jk) )
               ! mixed centered / upstream scheme
               zwx(ji,jj,jk) = zcofi * zupsut + (1.-zcofi) * zcenut
               zwy(ji,jj,jk) = zcofj * zupsvt + (1.-zcofj) * zcenvt
               zww(ji,jj,jk) = zcofi * zupsus + (1.-zcofi) * zcenus
               zwz(ji,jj,jk) = zcofj * zupsvs + (1.-zcofj) * zcenvs
            END DO
         END DO

         ! 2. Tracer flux divergence at t-point added to the general trend
         ! ---------------------------------------------------------------

         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
#else
               zbtr = zbtr2(ji,jj) 
#endif
               ! horizontal advective trends
               zta = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
                  &            + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
               zsa = - zbtr * (  zww(ji,jj,jk) - zww(ji-1,jj  ,jk)   &
                  &            + zwz(ji,jj,jk) - zwz(ji  ,jj-1,jk)  )
               ! add it to the general tracer trends
               ta(ji,jj,jk) = ta(ji,jj,jk) + zta
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! 3. Save the horizontal advective trends for diagnostic
      ! ------------------------------------------------------
      IF( l_trdtra )   THEN
         ! Recompute the hoizontal advection zta & zsa trends computed 
         ! at the step 2. above in making the difference between the new 
         ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
         ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)

         CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)

         ! Save the new ta and sa trends
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 

      ENDIF

      IF(ln_ctl)   THEN
          CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 had  - Ta: ', mask1=tmask, &
             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      ! "zonal" mean advective heat and salt transport 
      IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
# if defined key_s_coord || defined key_partial_steps
         pht_adv(:) = ptr_vj( zwy(:,:,:) )
         pst_adv(:) = ptr_vj( zwz(:,:,:) )
# else
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwy(ji,jj,jk) = zwy(ji,jj,jk) * fse3v(ji,jj,jk)
                  zwz(ji,jj,jk) = zwz(ji,jj,jk) * fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
         pht_adv(:) = ptr_vj( zwy(:,:,:) )
         pst_adv(:) = ptr_vj( zwz(:,:,:) )
# endif
      ENDIF


      ! II. Vertical advection
      ! ----------------------

      ! Bottom value : flux set to zero
      zwx(:,:,jpk) = 0.e0     ;    zwy(:,:,jpk) = 0.e0

      ! Surface value
      IF( lk_dynspg_rl ) THEN
         ! rigid lid : flux set to zero
         zwx(:,:, 1 ) = 0.e0    ;    zwy(:,:, 1 ) = 0.e0
      ELSE
         ! free surface
         zwx(:,:, 1 ) = zwn(:,:,1) * tn(:,:,1)
         zwy(:,:, 1 ) = zwn(:,:,1) * sn(:,:,1)
      ENDIF

      ! 1. Vertical advective fluxes
      ! ----------------------------

      ! Second order centered tracer flux at w-point

      DO jk = 2, jpk
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! upstream indicator
               zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) )
               ! velocity * 1/2
               zhw = 0.5 * zwn(ji,jj,jk)
               ! upstream scheme
               zfp_w = zhw + ABS( zhw )
               zfm_w = zhw - ABS( zhw )
               zupst = zfp_w * tb(ji,jj,jk) + zfm_w * tb(ji,jj,jk-1)
               zupss = zfp_w * sb(ji,jj,jk) + zfm_w * sb(ji,jj,jk-1)
               ! centered scheme
               zcent = zhw * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) )
               zcens = zhw * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) )
               ! mixed centered / upstream scheme
               zwx(ji,jj,jk) = zcofk * zupst + (1.-zcofk) * zcent
               zwy(ji,jj,jk) = zcofk * zupss + (1.-zcofk) * zcens
            END DO
         END DO
      END DO


      ! 2. Tracer flux divergence at t-point added to the general trend
      ! -------------------------
  
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ze3tr = 1. / fse3t(ji,jj,jk)
               ! vertical advective trends
               zta = - ze3tr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
               zsa = - ze3tr * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) )
               ! add it to the general tracer trends
               ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
               sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
            END DO
         END DO
      END DO

      ! 3. Save the vertical advective trends for diagnostic
      ! ----------------------------------------------------
      IF( l_trdtra )   THEN
         ! Recompute the vertical advection zta & zsa trends computed 
         ! at the step 2. above in making the difference between the new 
         ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
         ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)

         CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
      ENDIF

      IF(ln_ctl)   THEN
          CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 zad  - Ta: ', mask1=tmask, &
             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

   END SUBROUTINE tra_adv_cen2

#endif

   !!======================================================================
END MODULE traadv_cen2
