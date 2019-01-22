MODULE ldfslp
   !!======================================================================
   !!                       ***  MODULE  ldfslp  ***
   !! Ocean physics: slopes of neutral surfaces
   !!======================================================================
#if   defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                      Rotation of lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   ldf_slp      : compute the slopes of neutral surface
   !!   ldf_slp_mxl  : compute the slopes of iso-neutral surface
   !!   ldf_slp_init : initialization of the slopes computation
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce
   USE phycst          ! physical constants
   USE zdfmxl          ! mixed layer depth
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC ldf_slp        ! routine called by step.F90

   !! * Share module variables
   LOGICAL , PUBLIC, PARAMETER ::   lk_ldfslp = .TRUE.     !: slopes flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      uslp, wslpi,         &  !: i_slope at U- and W-points
      vslp, wslpj             !: j-slope at V- and W-points
   
   !! * Module variables
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
      omlmask                 ! mask of the surface mixed layer at T-pt
   REAL(wp), DIMENSION(jpi,jpj) ::   &
      uslpml, wslpiml,     &  ! i_slope at U- and W-points just below
   !                          ! the surface mixed layer
      vslpml, wslpjml         ! j_slope at V- and W-points just below
   !                          ! the surface mixed layer

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/LDF/ldfslp.F90,v 1.2 2005/11/16 16:13:32 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ldf_slp( kt, prd, pn2 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE ldf_slp  ***
      !! 
      !! ** Purpose :   Compute the slopes of neutral surface (slope of iso-
      !!      pycnal surfaces referenced locally) ('key_traldfiso').
      !!
      !! ** Method  :   The slope in the i-direction is computed at U- and 
      !!      W-points (uslp, wslpi) and the slope in the j-direction is 
      !!      computed at V- and W-points (vslp, wslpj).
      !!      They are bounded by 1/100 over the whole ocean, and within the
      !!      surface layer they are bounded by the distance to the surface
      !!      ( slope<= depth/l  where l is the length scale of horizontal
      !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
      !!      of 10cm/s)
      !!        A horizontal shapiro filter is applied to the slopes
      !!        'key_s_coord' defined: add to the previously computed slopes
      !!      the slope of the model level surface.
      !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
      !!      [slopes already set to zero at level 1, and to zero or the ocean
      !!      bottom slope ('key_s_coord' defined) at level jpk in inildf]
      !!
      !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes 
      !!               of now neutral surfaces at u-, w- and v- w-points, resp.
      !!
      !! History :
      !!   7.0  !  94-12  (G. Madec, M. Imbard)  Original code
      !!   8.0  !  97-06  (G. Madec)  optimization, lbc
      !!   8.1  !  99-10  (A. Jouzeau)  NEW profile
      !!   8.5  !  99-10  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce            , zgru  => ua,  &  ! use ua as workspace
                           zgrv  => va,  &  ! use va as workspace
                           zwy   => ta,  &  ! use ta as workspace
                           zwz   => sa      ! use sa as workspace
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         prd,                            &  ! in situ density
         pn2                                ! Brunt-Vaisala frequency (locally ref.)

      !! * Local declarations
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1      ! temporary integer
#if defined key_partial_steps
      INTEGER  ::   iku, ikv  ! temporary integers
#endif
      REAL(wp) ::   &
         zeps, zmg, zm05g, zcoef1, zcoef2,   &  ! temporary scalars
         zau, zbu, zav, zbv,                 &
         zai, zbi, zaj, zbj,                &
         zcofu, zcofv, zcofw,                &
         z1u, z1v, z1wu, z1wv,               &
         zalpha
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zww
      !!----------------------------------------------------------------------

      
      ! 0. Initialization (first time-step only)
      !    --------------
      
      IF( kt == nit000 ) CALL ldf_slp_init
      

      ! 0. Local constant initialization
      ! --------------------------------
      
      zeps  =  1.e-20
      zmg   = -1.0 / grav
      zm05g = -0.5 / grav

      zww(:,:,:) = 0.e0
      zwz(:,:,:) = 0.e0

      ! horizontal density gradient computation
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zgru(ji,jj,jk) = umask(ji,jj,jk) * ( prd(ji+1,jj  ,jk) - prd(ji,jj,jk) ) 
               zgrv(ji,jj,jk) = vmask(ji,jj,jk) * ( prd(ji  ,jj+1,jk) - prd(ji,jj,jk) ) 
            END DO
         END DO
      END DO

#if defined key_partial_steps
      ! partial steps correction at the bottom ocean level (zps_hde routine)
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
# else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
# endif
            ! last ocean level
            iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj) ) - 1
            ikv = MIN ( mbathy(ji,jj), mbathy(ji,jj+1) ) - 1
            zgru(ji,jj,iku) = gru(ji,jj) 
            zgrv(ji,jj,ikv) = grv(ji,jj)               
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO
#endif

      ! Slopes of isopycnal surfaces just below the mixed layer
      ! -------------------------------------------------------
      
      CALL ldf_slp_mxl( prd, pn2 )
      
      !-------------------synchro---------------------------------------------

      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! I.  slopes at u and v point
         ! ===========================
         
         
         ! I.1. Slopes of isopycnal surfaces
         ! ---------------------------------
         ! uslp = d/di( prd ) / d/dz( prd )
         ! vslp = d/dj( prd ) / d/dz( prd )
         
         ! Local vertical density gradient evaluated from N^2
         ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point

         DO jj = 1, jpj
            DO ji = 1, jpi
               zwy(ji,jj,jk) = zmg * ( prd(ji,jj,jk) + 1. )                &
                  &          * ( pn2(ji,jj,jk) + pn2(ji,jj,jk+1) )       &
                  &          / MAX( tmask(ji,jj,jk) + tmask (ji,jj,jk+1), 1. )
            END DO
         END DO
         
         ! Slope at u and v points
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! horizontal and vertical density gradient at u- and v-points
               zau = 1. / e1u(ji,jj) * zgru(ji,jj,jk)
               zav = 1. / e2v(ji,jj) * zgrv(ji,jj,jk)
               zbu = 0.5 * ( zwy(ji,jj,jk) + zwy(ji+1,jj  ,jk) )
               zbv = 0.5 * ( zwy(ji,jj,jk) + zwy(ji  ,jj+1,jk) )
               ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
               !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
               zbu = MIN( zbu, -100.*ABS( zau ), -7.e+3/fse3u(ji,jj,jk)*ABS( zau ) )
               zbv = MIN( zbv, -100.*ABS( zav ), -7.e+3/fse3v(ji,jj,jk)*ABS( zav ) )
               ! uslp and vslp output in zwz and zww, resp.
               zalpha = MAX( omlmask(ji,jj,jk), omlmask(ji+1,jj,jk) )
#if defined key_s_coord 
               zwz (ji,jj,jk) = ( zau / ( zbu - zeps ) * ( 1. - zalpha)   &
                  &        + zalpha * uslpml(ji,jj)   &
                  &        * ( fsdepu(ji,jj,jk) - .5*fse3u(ji,jj,1) )   &
                  &        / MAX( hmlpt(ji,jj), hmlpt(ji+1,jj), 5. ) )   &
                  &        * umask(ji,jj,jk)
               zalpha = MAX( omlmask(ji,jj,jk), omlmask(ji,jj+1,jk) )
               zww (ji,jj,jk) = ( zav / ( zbv - zeps ) * ( 1. - zalpha)   &
                  &        + zalpha * vslpml(ji,jj)   &
                  &        * ( fsdepv(ji,jj,jk) - .5*fse3v(ji,jj,1) )   &
                  &         / MAX( hmlpt(ji,jj), hmlpt(ji,jj+1), 5. ) )   &
                  &        * vmask(ji,jj,jk)
#else
               ! z-coord and partial steps slope computed in the same way
               zwz (ji,jj,jk) = ( zau / ( zbu - zeps ) * ( 1. - zalpha)    &
                  &        + zalpha * uslpml(ji,jj)    &
                  &        * ( fsdept(ji,jj,jk) - .5*fse3u(ji,jj,1))    &
                  &        / MAX (hmlpt(ji,jj),hmlpt(ji+1,jj),5.) )    &
                  &        * umask (ji,jj,jk)
               zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj+1,jk))
               zww (ji,jj,jk) = ( zav / ( zbv - zeps ) * ( 1. - zalpha)    &
                  &        + zalpha * vslpml(ji,jj)    &
                  &        * ( fsdept(ji,jj,jk) - .5*fse3v(ji,jj,1))    &
                  &        / MAX(hmlpt(ji,jj),hmlpt(ji,jj+1),5.) )    &
                  &        * vmask (ji,jj,jk)
#endif
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   end of slab
      !                                                ! ===============


         ! lateral boundary conditions on zww and zwz
         CALL lbc_lnk( zwz, 'U', -1. )
         CALL lbc_lnk( zww, 'V', -1. )

      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         
         ! Shapiro filter applied in the horizontal direction
         zcofu = 1. / 16.
         zcofv = 1. / 16.
         DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1 only
            DO ji = 2, jpim1  
               !uslop
               uslp(ji,jj,jk) = zcofu * (       zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
                  &                       +     zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
                  &                       + 2.*(zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
                  &                       +     zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
                  &                       + 4.* zwz(ji  ,jj  ,jk)                       )
               ! vslop
               vslp(ji,jj,jk) = zcofv * (       zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
                  &                       +     zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
                  &                       + 2.*(zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
                  &                       +     zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
                  &                       + 4.* zww(ji,jj    ,jk)                       )
            END DO
         END DO

         DO jj = 3, jpj-2
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! uslop
               uslp(ji,jj,jk) = zcofu * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
                  &                       +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
                  &                       + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
                  &                       +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zwz(ji  ,jj  ,jk)                       )
               ! vslop
               vslp(ji,jj,jk) = zcofv * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
                  &                       +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
                  &                       + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
                  &                       +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
                  &                       + 4.*  zww(ji,jj    ,jk)                       )
            END DO
         END DO

         ! decrease along coastal boundaries
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z1u  = ( umask(ji,jj+1,jk) + umask(ji,jj-1,jk) )*.5
               z1v  = ( vmask(ji+1,jj,jk) + vmask(ji-1,jj,jk) )*.5
               z1wu = ( umask(ji,jj,jk)   + umask(ji,jj,jk+1) )*.5
               z1wv = ( vmask(ji,jj,jk)   + vmask(ji,jj,jk+1) )*.5
               uslp(ji,jj,jk) = uslp(ji,jj,jk) * z1u * z1wu
               vslp(ji,jj,jk) = vslp(ji,jj,jk) * z1v * z1wv
            END DO
         END DO


         IF( lk_sco ) THEN
            ! Add the slope of level surfaces
            ! -----------------------------------
            ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
            ! in inildf, ldfslp never called
            ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
            ! is added to the slope of isopycnal surfaces.
            ! c a u t i o n : minus sign as fsdep has positive value 
         
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  uslp(ji,jj,jk) = uslp(ji,jj,jk) - 1. / e1u(ji,jj)   &
                     &           * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )
                  vslp(ji,jj,jk) = vslp(ji,jj,jk) - 1. / e2v(ji,jj)   &
                     &           * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
               END DO
            END DO
         ENDIF


         ! II. Computation of slopes at w point
         ! ====================================
         
         
         ! II.1 Slopes of isopycnal surfaces
         ! ---------------------------------
         ! wslpi = mij( d/di( prd ) / d/dz( prd )
         ! wslpj = mij( d/dj( prd ) / d/dz( prd )
         
         
         ! Local vertical density gradient evaluated from N^2
         !     zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
         DO jj = 1, jpj
            DO ji = 1, jpi
               zwy (ji,jj,jk) = zm05g * pn2 (ji,jj,jk) *   &
                  &                     ( prd (ji,jj,jk) + prd (ji,jj,jk-1) + 2. )
            END DO
         END DO
         
         ! Slope at w point
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! horizontal density i-gradient at w-points
               zcoef1 = MAX( zeps, umask(ji-1,jj,jk  )+umask(ji,jj,jk  )    &
                  &               +umask(ji-1,jj,jk-1)+umask(ji,jj,jk-1) )
               zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
               zai = zcoef1 * (  zgru(ji  ,jj,jk  ) + zgru(ji  ,jj,jk-1)   &
                  &            + zgru(ji-1,jj,jk-1) + zgru(ji-1,jj,jk  ) ) * tmask (ji,jj,jk)
               ! horizontal density j-gradient at w-points
               zcoef2 = MAX( zeps, vmask(ji,jj-1,jk  )+vmask(ji,jj,jk-1)   &
                  &               +vmask(ji,jj-1,jk-1)+vmask(ji,jj,jk  ) )
               zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
               zaj = zcoef2 * (  zgrv(ji,jj  ,jk  ) + zgrv(ji,jj  ,jk-1)   &
                  &            + zgrv(ji,jj-1,jk-1) + zgrv(ji,jj-1,jk  ) ) * tmask (ji,jj,jk)
               ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
               !                   static instability:
               !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
               zbi = MIN( zwy (ji,jj,jk),- 100.*ABS(zai), -7.e+3/fse3w(ji,jj,jk)*ABS(zai) )
               zbj = MIN( zwy (ji,jj,jk), -100.*ABS(zaj), -7.e+3/fse3w(ji,jj,jk)*ABS(zaj) )
               ! wslpi and wslpj output in zwz and zww, resp.
               zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj,jk-1))
               zwz(ji,jj,jk) = ( zai / ( zbi - zeps) * ( 1. - zalpha )   &
                  &            + zalpha * wslpiml(ji,jj)   &
                  &            * fsdepw(ji,jj,jk) / MAX( hmlp(ji,jj),10. ) )   &
                  &            * tmask (ji,jj,jk)
               zww(ji,jj,jk) = ( zaj / ( zbj - zeps) * ( 1. - zalpha )   &
                  &            + zalpha * wslpjml(ji,jj)   &
                  &            * fsdepw(ji,jj,jk) / MAX( hmlp(ji,jj),10. ) )   &
                  &            * tmask (ji,jj,jk)
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   end of slab
      !                                                ! ===============


      ! lateral boundary conditions on zwz and zww
      CALL lbc_lnk( zwz, 'T', -1. )
      CALL lbc_lnk( zww, 'T', -1. )

      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Shapiro filter applied in the horizontal direction

         DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1
            DO ji = 2, jpim1
               zcofw = tmask(ji,jj,jk)/16.
               wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
                  &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
                  &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
                  &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
                  &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw

               wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
                  &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
                  &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
                  &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
                  &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
            END DO
         END DO
         
         DO jj = 3, jpj-2
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcofw = tmask(ji,jj,jk)/16.
               wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
                  &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
                  &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
                  &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
                  &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw

               wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
                  &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
                  &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
                  &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
                  &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
            END DO
         END DO
         
         ! decrease the slope along the boundaries
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               z1u = ( umask(ji,jj,jk) + umask(ji-1,jj,jk) ) *.5
               z1v = ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) *.5
               wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * z1u * z1v
               wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * z1u * z1v
            END DO
         END DO
         
         IF( lk_sco ) THEN
         
            ! Slope of level surfaces
            ! -----------------------
            ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
            ! in inildf, ldfslp never called
            ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
            ! is added to the slope of isopycnal surfaces.
         
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  wslpi(ji,jj,jk) = wslpi(ji,jj,jk) - 1. / e1t(ji,jj)   &
                     &                                   * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
                  wslpj(ji,jj,jk) = wslpj(ji,jj,jk) - 1. / e2t(ji,jj)   &
                     &                                   * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
               END DO
            END DO
         ENDIF
         
         ! III. Specific grid points
         ! -------------------------

         IF( cp_cfg == "orca" .AND. jp_cfg == 4 ) THEN
            !                                        ! =======================
            ! Horizontal diffusion in                !  ORCA_R4 configuration 
            ! specific area                          ! =======================
            !
            !                                             ! Gibraltar Strait
            ij0 =  50   ;   ij1 =  53
            ii0 =  69   ;   ii1 =  71   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ij0 =  51   ;   ij1 =  53
            ii0 =  68   ;   ii1 =  71   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ii0 =  69   ;   ii1 =  71   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ii0 =  69   ;   ii1 =  71   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0

            !                                             ! Mediterrannean Sea
            ij0 =  49   ;   ij1 =  56
            ii0 =  71   ;   ii1 =  90   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ij0 =  50   ;   ij1 =  56
            ii0 =  70   ;   ii1 =  90   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ii0 =  71   ;   ii1 =  90   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
            ii0 =  71   ;   ii1 =  90   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
         ENDIF
         !                                             ! ===============
      END DO                                           !   end of slab
      !                                                ! ===============       

      
      ! III Lateral boundary conditions on all slopes (uslp , vslp, 
      ! -------------------------------                wslpi, wslpj )
      CALL lbc_lnk( uslp , 'U', -1. )
      CALL lbc_lnk( vslp , 'V', -1. )
      CALL lbc_lnk( wslpi, 'W', -1. )
      CALL lbc_lnk( wslpj, 'W', -1. )

   END SUBROUTINE ldf_slp


   SUBROUTINE ldf_slp_mxl( prd, pn2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_mxl  ***
      !! ** Purpose :
      !!     Compute the slopes of iso-neutral surface (slope of isopycnal
      !!   surfaces referenced locally) just above the mixed layer.
      !!
      !! ** Method :
      !!      The slope in the i-direction is computed at u- and w-points
      !!   (uslp, wslpi) and the slope in the j-direction is computed at
      !!   v- and w-points (vslp, wslpj).
      !!   They are bounded by 1/100 over the whole ocean, and within the
      !!   surface layer they are bounded by the distance to the surface
      !!   ( slope<= depth/l  where l is the length scale of horizontal
      !!   diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
      !!   of 10cm/s)
      !!
      !! ** Action :
      !!      Compute uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
      !!   of now neutral surfaces at u-, w- and v- w-points, resp.
      !!
      !! History :
      !!   8.1  !  99-10  (A. Jouzeau)  Original code
      !!   8.5  !  99-10  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE oce           , zgru  => ua,  &  ! ua, va used as workspace and set to hor. 
                          zgrv  => va      ! density gradient in ldf_slp

      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         prd,                           &  ! in situ density
         pn2                              ! Brunt-Vaisala frequency (locally ref.)

      !! * Local declarations
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      INTEGER  ::   ik, ikm1               ! temporary integers
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zwy                               ! temporary workspace
      REAL(wp) ::   &
         zeps, zmg, zm05g,              &  ! temporary scalars
         zcoef1, zcoef2,                &  !    "         "
         zau, zbu, zav, zbv,            &  !    "         "
         zai, zbi, zaj, zbj                !    "         "
      !!----------------------------------------------------------------------


      ! 0. Local constant initialization
      ! --------------------------------

      zeps  =  1.e-20
      zmg   = -1.0 / grav
      zm05g = -0.5 / grav


      uslpml (1,:) = 0.e0      ;      uslpml (jpi,:) = 0.e0
      vslpml (1,:) = 0.e0      ;      vslpml (jpi,:) = 0.e0
      wslpiml(1,:) = 0.e0      ;      wslpiml(jpi,:) = 0.e0
      wslpjml(1,:) = 0.e0      ;      wslpjml(jpi,:) = 0.e0

      ! surface mixed layer mask

      ! mask for mixed layer
      DO jk = 1, jpk
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = 1, jpij   ! vector opt. (forced unrolling)
# else
         DO jj = 1, jpj
            DO ji = 1, jpi
# endif
               ! mixed layer interior (mask = 1) and exterior (mask = 0)
               ik = nmln(ji,jj) - 1
               IF( jk <= ik ) THEN
                  omlmask(ji,jj,jk) = 1.e0
               ELSE
                  omlmask(ji,jj,jk) = 0.e0
               ENDIF
# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO
      END DO


      ! Slopes of isopycnal surfaces just before bottom of mixed layer
      ! --------------------------------------------------------------
      ! uslpml = d/di( prd ) / d/dz( prd )
      ! vslpml = d/dj( prd ) / d/dz( prd )

      ! Local vertical density gradient evaluated from N^2
      ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point

      !-----------------------------------------------------------------------
      zwy(:,jpj) = 0.e0
      zwy(jpi,:) = 0.e0
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
# else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
# endif
            ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
            ! if ik = jpk take jpkm1 values
            ik = MIN( ik,jpkm1 )
            zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
               &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
               &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO
      ! lateral boundary conditions on zwy
      CALL lbc_lnk( zwy, 'U', -1. )

      ! Slope at u points
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
# endif
            ! horizontal and vertical density gradient at u-points
            ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
            ik = MIN( ik,jpkm1 )
            zau = 1./ e1u(ji,jj) * zgru(ji,jj,ik)
            zbu = 0.5*( zwy(ji,jj) + zwy(ji+1,jj) )
            ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
            !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbu = MIN( zbu, -100.*ABS(zau), -7.e+3/fse3u(ji,jj,ik)*ABS(zau) )
            ! uslpml
            uslpml (ji,jj) = zau / ( zbu - zeps ) * umask (ji,jj,ik)
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! lateral boundary conditions on uslpml
      CALL lbc_lnk( uslpml, 'U', -1. )

      ! Local vertical density gradient evaluated from N^2
      !     zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
      zwy ( :, jpj) = 0.e0
      zwy ( jpi, :) = 0.e0
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
# else
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
# endif
            ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
            ik = MIN( ik,jpkm1 )
            zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
               &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
               &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! lateral boundary conditions on zwy
      CALL lbc_lnk( zwy, 'V', -1. )

      ! Slope at v points
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
# endif
            ! horizontal and vertical density gradient at v-points
            ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
            ik = MIN( ik,jpkm1 )
            zav = 1./ e2v(ji,jj) * zgrv(ji,jj,ik)
            zbv = 0.5*( zwy(ji,jj) + zwy(ji,jj+1) )
            ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
            !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbv = MIN( zbv, -100.*ABS(zav), -7.e+3/fse3v(ji,jj,ik)*ABS( zav ) )
            ! vslpml
            vslpml (ji,jj) = zav / ( zbv - zeps ) * vmask (ji,jj,ik)
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! lateral boundary conditions on vslpml
      CALL lbc_lnk( vslpml, 'V', -1. )

      ! wslpiml = mij( d/di( prd ) / d/dz( prd )
      ! wslpjml = mij( d/dj( prd ) / d/dz( prd )


      ! Local vertical density gradient evaluated from N^2
      ! zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = 1, jpij   ! vector opt. (forced unrolling)
# else
      DO jj = 1, jpj
         DO ji = 1, jpi
# endif
            ik = nmln(ji,jj)+1
            ik = MIN( ik,jpk )
            ikm1 = MAX ( 1, ik-1)
            zwy (ji,jj) = zm05g * pn2 (ji,jj,ik) *     &
               &             ( prd (ji,jj,ik) + prd (ji,jj,ikm1) + 2. )
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! Slope at w point
# if defined key_vectopt_loop   &&   ! defined key_autotasking
      jj = 1
      DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
# endif
            ik = nmln(ji,jj)+1
            ik = MIN( ik,jpk )
            ikm1 = MAX ( 1, ik-1)
            ! horizontal density i-gradient at w-points
            zcoef1 = MAX( zeps, umask(ji-1,jj,ik  )+umask(ji,jj,ik  )   &
               &               +umask(ji-1,jj,ikm1)+umask(ji,jj,ikm1) )
            zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
            zai = zcoef1 * (  zgru(ji  ,jj,ik  ) + zgru(ji  ,jj,ikm1)   &
               &            + zgru(ji-1,jj,ikm1) + zgru(ji-1,jj,ik  ) ) * tmask (ji,jj,ik)
            ! horizontal density j-gradient at w-points
            zcoef2 = MAX( zeps, vmask(ji,jj-1,ik  )+vmask(ji,jj,ikm1)    &
               &               +vmask(ji,jj-1,ikm1)+vmask(ji,jj,ik  ) )
            zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
            zaj = zcoef2 * (  zgrv(ji,jj  ,ik  ) + zgrv(ji,jj  ,ikm1)   &
               &            + zgrv(ji,jj-1,ikm1) + zgrv(ji,jj-1,ik  ) ) * tmask (ji,jj,ik)
            ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
            !                   static instability:
            !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
            zbi = MIN ( zwy (ji,jj),- 100.*ABS(zai), -7.e+3/fse3w(ji,jj,ik)*ABS(zai) )
            zbj = MIN ( zwy (ji,jj), -100.*ABS(zaj), -7.e+3/fse3w(ji,jj,ik)*ABS(zaj) )
            ! wslpiml and wslpjml
            wslpiml (ji,jj) = zai / ( zbi - zeps) * tmask (ji,jj,ik)
            wslpjml (ji,jj) = zaj / ( zbj - zeps) * tmask (ji,jj,ik)
# if ! defined key_vectopt_loop   ||   defined key_autotasking
         END DO
# endif
      END DO

      ! lateral boundary conditions on wslpiml and wslpjml
      CALL lbc_lnk( wslpiml, 'W', -1. )
      CALL lbc_lnk( wslpjml, 'W', -1. )

   END SUBROUTINE ldf_slp_mxl


   SUBROUTINE ldf_slp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_slp_init  ***
      !!
      !! ** Purpose :   Initialization for the isopycnal slopes computation
      !!
      !! ** Method  :   read the nammbf namelist and check the parameter 
      !!      values called by tra_dmp at the first timestep (nit000)
      !!
      !! History :
      !!  8.5  ! 02-06 (G. Madec) original code
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      
      
      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_slp : direction of lateral mixing'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Direction of lateral diffusion (tracers and/or momentum)
      ! ------------------------------
      ! set the slope to zero (even in s-coordinates)

      uslp (:,:,:) = 0.e0
      vslp (:,:,:) = 0.e0
      wslpi(:,:,:) = 0.e0
      wslpj(:,:,:) = 0.e0

      uslpml (:,:) = 0.e0
      vslpml (:,:) = 0.e0
      wslpiml(:,:) = 0.e0
      wslpjml(:,:) = 0.e0

      IF( ln_traldf_hor ) THEN

         ! geopotential diffusion in s-coordinates on tracers and/or momentum
         ! The slopes of s-surfaces are computed once (no call to ldfslp in step)
         ! The slopes for momentum diffusion are i- or j- averaged of those on tracers

         ! set the slope of diffusion to the slope of s-surfaces
         !      ( c a u t i o n : minus sign as fsdep has positive value )
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  uslp (ji,jj,jk) = -1. / e1u(ji,jj) * umask(ji,jj,jk)   &
                     &                               * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )
                  vslp (ji,jj,jk) = -1. / e2v(ji,jj) * vmask(ji,jj,jk)   &
                     &                               * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
                  wslpi(ji,jj,jk) = -1. / e1t(ji,jj) * tmask(ji,jj,jk)   &
                     &                               * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
                  wslpj(ji,jj,jk) = -1. / e2t(ji,jj) * tmask(ji,jj,jk)   &
                     &                               * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
               END DO
            END DO
         END DO

         ! Lateral boundary conditions on the slopes
         CALL lbc_lnk( uslp , 'U', -1. )
         CALL lbc_lnk( vslp , 'V', -1. )
         CALL lbc_lnk( wslpi, 'W', -1. )
         CALL lbc_lnk( wslpj, 'W', -1. )
      ENDIF

   END SUBROUTINE ldf_slp_init

#else
   !!------------------------------------------------------------------------
   !!   Dummy module :                 NO Rotation of lateral mixing tensor
   !!------------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_ldfslp = .FALSE.    !: slopes flag
CONTAINS
   SUBROUTINE ldf_slp( kt, prd, pn2 )        ! Dummy routine
      INTEGER, INTENT(in) :: kt 
      REAL,DIMENSION(:,:,:), INTENT(in) :: prd, pn2
      WRITE(*,*) 'ldf_slp: You should not have seen this print! error?', kt, prd(1,1,1), pn2(1,1,1)
   END SUBROUTINE ldf_slp
#endif

   !!======================================================================
END MODULE ldfslp
