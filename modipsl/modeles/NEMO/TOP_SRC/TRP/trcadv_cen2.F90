MODULE trcadv_cen2
   !!==============================================================================
   !!                       ***  MODULE  trcadv_cen2  ***
   !! Ocean passive tracers:  horizontal & vertical advective tracer trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_adv_cen2 : update the tracer trend with the horizontal
   !!                  and vertical advection trends using a 2nd order 
   !!                  centered finite difference scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE trcbbl              ! advective passive tracers in the BBL
   USE prtctl_trc
#ifdef key_RIVER_INPUT
   USE rivers
#endif

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_adv_cen2    ! routine called by trcstp.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_cen2.F90,v 1.11 2006/04/10 15:38:53 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option :             2nd order centered scheme (k-j-i loop)
   !!----------------------------------------------------------------------

   SUBROUTINE trc_adv_cen2( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_cen2  ***
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
      !!      * horizontal advective trend (divergence of the fluxes)
      !!         * s-coordinate (lk_sco=T) or
      !!         * z-coordinate with partial steps (lk_zps=T)
      !!               ztra = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
      !!         * z-coordinate (default key), e3t=e3u=e3v:
      !!               ztra = 1/(e1t*e2t) { di-1[zwx] + dj-1[zwy] }
      !!      * Add this trend now to the general trend of tracer tra:
      !!              tra = tra + ztra
      !!      * trend diagnostic ('key_trc_diatrd'): the trend is saved
      !!      for diagnostics. The trends saved is expressed as
      !!      Uh.gradh(T)
      !!
      !!         Part II : vertical advection
      !!      For any tracer  the advective trend is computed as follows :
      !!            ztra = 1/e3t dk+1[ zwz ]
      !!      where the vertical advective flux, zwz, is given by :
      !!            zwz = zcofk * zupst + (1-zcofk) * zcent
      !!      with 
      !!        zupsv = upstream flux = wn * (trb(k) or trb(k-1) ) [wn>0 or <0]
      !!        zcenu = centered flux = wn * mk(trn)
      !!         The surface boundary condition is : 
      !!      rigid-lid (default option) : zero advective flux
      !!      free-surf ("key_fresurf_cstvol") : wn(:,:,1) * trn(:,:,1)
      !!         Add this trend now to the general trend of tracer tra :
      !!            tra = tra + ztra
      !!         Trend diagnostic ('key_trc_diatrd'): the trend is saved for
      !!      diagnostics. The trends saved is expressed as :
      !!             save trend =  w.gradz(T) = ztra - trn divn.
      !!
      !! ** Action : - update tra with the now advective tracer trends
      !!             - save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   8.2  !  01-08  (M-A Filiberti, and M.Levy)  trahad+trazad = traadv 
      !!   8.5  !  02-06  (G. Madec, C. Ethe)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce_trc            , zwx => ua,  &  ! use ua as workspace
         &                     zwy => va      ! use va as workspace
#if defined key_trcbbl_adv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
         &         zun, zvn, zwn
#else
      USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
         &                     zvn => vn,  &  ! When no bbl, zvn == vn
         &                     zwn => wn      ! When no bbl, zwn == wn
#endif
!!DBG
      USE obc_oce

 
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
 
      !! * Local save
      REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
         zbtr2
 
      !! * Local declarations
      INTEGER  ::   ji, jj, jk, jn             ! dummy loop indices
      REAL(wp) ::                           &
         zbtr, ztra, zfui, zfvj,            &  ! temporary scalars
         zhw, ze3tr, zcofi, zcofj,          &  !    "         "
         zupsut, zupsvt,                    &  !    "         "
         zfp_ui, zfp_vj, zfm_ui, zfm_vj,    &  !    "         "
         zcofk, zupst, zcent,               &  !    "         "
         zfp_w, zfm_w,                      &  !    "         "
         zcenut, zcenvt                        ! 

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zind                              ! temporary workspace arrays
#if defined key_trc_diatrd
      REAL(wp) ::                           &
         ztai, ztaj,                        &  ! temporary scalars
         zfui1, zfvj1                          !    "         "
#endif
      CHARACTER (len=22) :: charout
!!DB
      INTEGER :: ii, ij, ik, ilast, ifirst
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_adv_cen2 : 2nd order centered advection scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
         IF(lwp) WRITE(numout,*)
   
         zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
      ENDIF

!!DBG
      zwx(:,:,:) = 0.0; zwy(:,:,:) = 0.0


#if defined key_trcbbl_adv

      ! Advective bottom boundary layer 
      ! -------------------------------
      zun(:,:,:) = un(:,:,:) - u_trc_bbl(:,:,:)
      zvn(:,:,:) = vn(:,:,:) - v_trc_bbl(:,:,:)
      zwn(:,:,:) = wn(:,:,:) + w_trc_bbl(:,:,:)
#endif

#ifdef key_RIVER_INPUT
!NL#1 (same call as the routine tra_adv_tvd )
      call river_set_mask
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

!!DBG: change indices on east & west boundaries
      ilast = fs_jpim1
!      do ji = nie0, nie1
!         if(fs_jpim1 > nie0) ilast = nie0
!      enddo
      ifirst = 1
!      do ji = niw0, niw1
!         ifirst = 2
!      enddo


      DO jn = 1, jptra
         ! I. Horizontal advective fluxes
         ! ------------------------------
         
         ! Second order centered tracer flux at u and v-points
         
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpjm1
!               DO ji = 1, ilast   ! vector opt.
!!DBG
               DO ji = ifirst, ilast   ! vector opt.
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
                  zupsut = zfp_ui * trb(ji,jj,jk,jn) + zfm_ui * trb(ji+1,jj  ,jk,jn)
                  zupsvt = zfp_vj * trb(ji,jj,jk,jn) + zfm_vj * trb(ji  ,jj+1,jk,jn)
                  ! centered scheme
                  zcenut = zfui * ( trn(ji,jj,jk,jn) + trn(ji+1,jj  ,jk,jn) )
                  zcenvt = zfvj * ( trn(ji,jj,jk,jn) + trn(ji  ,jj+1,jk,jn) )
                  ! mixed centered / upstream scheme
                  zwx(ji,jj,jk) = zcofi * zupsut + (1.-zcofi) * zcenut
                  zwy(ji,jj,jk) = zcofj * zupsvt + (1.-zcofj) * zcenvt                
               END DO
            END DO


            ! 2. Tracer flux divergence at t-point added to the general trend
            ! -------------------------

            DO jj = 2, jpjm1
               DO ji = fs_2, ilast   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
                  zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
#else
                  zbtr = zbtr2(ji,jj) 
#endif
                  ! horizontal advective trends
                  ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
                     &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )

                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra

#if defined key_trc_diatrd 
                  ! recompute the trends in i- and j-direction as Uh gradh(T)
# if defined key_s_coord || defined key_partial_steps
                  zfui = 0.5 * e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)
                  zfui1= 0.5 * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
                  zfvj = 0.5 * e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)
                  zfvj1= 0.5 * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
# else
                  zfui = 0.5 * e2u(ji  ,jj) * zun(ji,  jj,jk)
                  zfui1= 0.5 * e2u(ji-1,jj) * zun(ji-1,jj,jk)
                  zfvj = 0.5 * e1v(ji,jj  ) * zvn(ji,jj  ,jk)
                  zfvj1= 0.5 * e1v(ji,jj-1) * zvn(ji,jj-1,jk)
# endif
                  ztai = - zbtr * ( zfui  * ( trn(ji+1,jj  ,jk,jn) - trn(ji,  jj,jk,jn) )   &
                     &                + zfui1 * ( trn(ji,  jj,  jk,jn) - trn(ji-1,jj,jk,jn) ) )
                  ztaj = - zbtr * ( zfvj  * ( trn(ji  ,jj+1,jk,jn) - trn(ji,jj  ,jk,jn) )    &
                     &                + zfvj1 * ( trn(ji  ,jj  ,jk,jn) - trn(ji,jj-1,jk,jn) ) )
                  ! save i- and j- advective trends computed as Uh gradh(T)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = ztai
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = ztaj
#endif
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============
      ENDDO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('centered2 - had')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF
      
      ! II. Vertical advection
      ! ----------------------
      DO jn = 1, jptra

         ! Bottom value : flux set to zero
         zwx(:,:,jpk) = 0.e0 

         ! Surface value
         IF ( lk_dynspg_rl ) THEN       ! rigid lid : flux set to zero
            zwx(:,:, 1 ) = 0.e0  
         ELSE                           ! free surface-constant volume
            zwx(:,:, 1 ) = zwn(:,:,1) * trn(:,:,1,jn)
         ENDIF

         ! 1. Vertical advective fluxes
         ! ----------------------------

         ! Second order centered tracer flux at w-point

         DO jk = 2, jpk
            DO jj = 2, jpjm1
               DO ji = fs_2, ilast   ! vector opt.
                  ! upstream indicator
                  zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) )
                  ! velocity * 1/2
                  zhw = 0.5 * zwn(ji,jj,jk)
                  ! upstream scheme
                  zfp_w = zhw + ABS( zhw )
                  zfm_w = zhw - ABS( zhw )
                  zupst = zfp_w * trb(ji,jj,jk,jn) + zfm_w * trb(ji,jj,jk-1,jn)
                  ! centered scheme
                  zcent = zhw * ( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) )
                  ! centered scheme
                  zwx(ji,jj,jk) = zcofk * zupst + (1.-zcofk) * zcent
               END DO
            END DO
         END DO


         ! 2. Tracer flux divergence at t-point added to the general trend
         ! -------------------------

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, ilast   ! vector opt.
                  ze3tr = 1. / fse3t(ji,jj,jk)
                  ! vertical advective trends
                  ztra = - ze3tr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
                  ! add it to the general tracer trends
                  tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra
#if defined key_trc_diatrd 
                  ! save the vertical advective trends computed as w gradz(T)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = ztra - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
#endif
               END DO
            END DO
         END DO

      END DO

#ifdef key_RIVER_INPUT
!NL#1 (same call as the routine tra_adv_tvd )
      call river_unset_mask
#endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('centered - zad')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_adv_cen2
#else

   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv_cen2( kt )  
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_adv_cen2: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv_cen2
#endif
   !!======================================================================
END MODULE trcadv_cen2
