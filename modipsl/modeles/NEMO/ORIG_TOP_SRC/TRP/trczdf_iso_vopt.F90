MODULE trczdf_iso_vopt
   !!==============================================================================
   !!                 ***  MODULE  trczdf_iso_vopt  ***
   !! Ocean passive tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
#if defined key_passivetrc && ( defined key_ldfslp   ||   defined key_esopa )
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   trc_zdf_iso_vopt : Update the tracer trend with the vertical part of 
   !!                  the isopycnal or geopotential s-coord. operator and
   !!                  the vertical diffusion. vector optimization, use
   !!                  k-j-i loops.
   !!   trc_zdf_iso  :
   !!   trc_zdf_zdf  :
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE trctrp_lec      ! passive tracers transport
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_zdf_iso_vopt   !  routine called by step.F90

   !! * Module variables
   REAL(wp), DIMENSION(jpk) ::  &
      rdttrc                          ! vertical profile of 2 x time-step

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_iso_vopt.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE trc_zdf_iso_vopt( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_iso_vopt  ***
      !!
      !! ** Purpose :
      !! ** Method  :
      !! ** Action  :
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe)   adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      CHARACTER (len=22) :: charout
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'trc_zdf_iso_vopt : vertical mixing computation'
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~~   is  iso-neutral diffusion : implicit vertical time stepping'
#if defined key_trcldf_eiv && defined key_diaeiv 
         w_trc_eiv(:,:,:) = 0.e0
#endif
      ENDIF


      ! I. vertical extra-diagonal part of the rotated tensor
      ! -----------------------------------------------------

      CALL trc_zdf_iso

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - 1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

      ! II. vertical diffusion (including the vertical diagonal part of the rotated tensor)
      ! ----------------------

      CALL trc_zdf_zdf( kt )

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - 2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_zdf_iso_vopt


   SUBROUTINE trc_zdf_zdf( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE trc_zdf_zdf  ***
      !!                    
      !! ** Purpose :   Compute the trend due to the vertical tracer diffusion
      !!     including the vertical component of lateral mixing (only for 2nd
      !!     order operator, for fourth order it is already computed and add
      !!     to the general trend in traldf.F) and add it to the general trend
      !!     of the tracer equations.
      !!
      !! ** Method  :   The vertical component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geo-
      !!      potential surfaces to which an eddy induced advection can be 
      !!      added. It is computed using before fields (forward in time) and 
      !!      isopycnal or geopotential slopes computed in routine ldfslp.
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers  is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=tra).
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend tra :
      !!         tra = tra + dz( avt dz(t) )
      !!         (tra = tra + dz( avs dz(t) ) if lk_trc_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in diamld).
      !!         avt = zavt
      !!         (avs = zavs if lk_trc_zdfddm=T )
      !!
      !!      'key_trdtra' defined: trend saved for futher diagnostics.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action  : - Update tra with before vertical diffusion trend
      !!              - Save the trend in trtrd  ('key_trc_diatrd')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11 (G. Madec)
      !!        !  92-06 (M. Imbard) correction on tracer trend loops
      !!        !  96-01 (G. Madec) statement function for e3
      !!        !  97-05 (G. Madec) vertical component of isopycnal
      !!        !  97-07 (G. Madec) geopotential diffusion in s-coord
      !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
      !!        !                     avt multiple correction
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce_trc, ONLY :   zwd   => ua,  &  ! ua, va used as
                            zws   => va      ! workspace
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn                ! dummy loop indices
      REAL(wp) ::   &
         zavi, zrhs                          ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwi, zwt, zavsi                     ! temporary workspace arrays
      REAL(wp) ::    ztra              !temporary scalars
#  if defined key_trc_diatrd
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztrd
#  endif
      !!---------------------------------------------------------------------


      ! I. Local constant initialization
      ! --------------------------------
      ! ... time step = 2 rdttra ex
      IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
         ! time step = 2 rdttra with Arakawa or TVD advection scheme
         IF( neuler == 0 .AND. kt == nittrc000 ) THEN
            rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
         ELSEIF( kt <= nittrc000 + 1 ) THEN
            rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
         ENDIF
      ELSE
         rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
      ENDIF

      DO jn = 1, jptra
         
         zwd( 1 ,:,:)=0.e0     ;     zwd(jpi,:,:)=0.e0
         zws( 1 ,:,:)=0.e0     ;     zws(jpi,:,:)=0.e0
         zwi( 1 ,:,:)=0.e0     ;     zwi(jpi,:,:)=0.e0

         zwt( 1 ,:,:)=0.e0     ;     zwt(jpi,:,:)=0.e0
         zwt(  :,:,1)=0.e0     ;     zwt(:,:,jpk)= 0.e0
         zavsi( 1 ,:,:)=0.e0   ;     zavsi(jpi,:,:)=0.e0 
         zavsi(  :,:,1)=0.e0   ;     zavsi(:,:,jpk)=0.e0

#  if defined key_trc_diatrd
         ! save the tra trend
         ztrd(:,:,:) = tra(:,:,:,jn)
#  endif

         ! II. Vertical trend associated with the vertical physics
         ! =======================================================
         !     (including the vertical flux proportional to dk[t] associated
         !      with the lateral mixing, through the avt update)
         !     dk[ avt dk[ (t,s) ] ] diffusive trends


         ! II.0 Matrix construction
         ! ------------------------        
         ! update and save of avt (and avs if double diffusive mixing)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zavi = fsahtw(ji,jj,jk) * (                 &   ! vertical mixing coef. due to lateral mixing
                     &                           wslpi(ji,jj,jk) * wslpi(ji,jj,jk)      &
                     &                         + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
                  zavsi(ji,jj,jk) = fstravs(ji,jj,jk) + zavi        ! dd mixing: zavsi = total vertical mixing coef. on tracer

               END DO
            END DO
         END DO


         ! II.2 Vertical diffusion on tracer
         ! ---------------------------========

         ! Rebuild the Matrix as avt /= avs

         ! Diagonal, inferior, superior  (including the bottom boundary condition via avs masked)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwi(ji,jj,jk) = - rdttrc(jk) * zavsi(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
                  zws(ji,jj,jk) = - rdttrc(jk) * zavsi(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
                  zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Surface boudary conditions
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwi(ji,jj,1) = 0.e0
               zwd(ji,jj,1) = 1. - zws(ji,jj,1)
            END DO
         END DO

         !! Matrix inversion from the first level
         !!----------------------------------------------------------------------
         !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
         !
         !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
         !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
         !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
         !        (        ...               )( ...  ) ( ...  )
         !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
         !
         !   m is decomposed in the product of an upper and lower triangular
         !   matrix
         !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
         !   The second member is in 2d array zwy
         !   The solution is in 2d array zwx
         !   The 3d arry zwt is a work space array
         !   zwy is used and then used as a work space array : its value is modified!

         ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwt(ji,jj,1) = zwd(ji,jj,1)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
               END DO
            END DO
         END DO

         ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               tra(ji,jj,1,jn) = trb(ji,jj,1,jn) + rdttrc(1) * tra(ji,jj,1,jn)
            END DO
         END DO
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  zrhs = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)   ! zrhs=right hand side
                  tra(ji,jj,jk,jn) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * tra(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO

         ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
         ! Save the masked passive tracer after in tra
         ! (c a u t i o n: passive tracer not its trend, Leap-frog scheme done it will not be done in tranxt)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               tra(ji,jj,jpkm1,jn) = tra(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
            END DO
         END DO
         DO jk = jpk-2, 1, -1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  tra(ji,jj,jk,jn) = ( tra(ji,jj,jk,jn) - zws(ji,jj,jk) * tra(ji,jj,jk+1,jn) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

#if defined key_trc_diatrd
         ! Compute and save the vertical diffusive passive tracer trends
#  if defined key_trc_ldfiso 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztra = ( tra(ji,jj,jk,jn) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - ztrd(ji,jj,jk) + trtrd(ji,jj,jk,ikeep(jn),6)
               END DO
            END DO
         END DO
#  else
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztra = ( tra(ji,jj,jk,jn) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - ztrd(ji,jj,jk)
               END DO
            END DO
         END DO
#  endif
#endif

      END DO

   END SUBROUTINE trc_zdf_zdf


   SUBROUTINE trc_zdf_iso
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_iso  ***
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
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(trb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(trb)) ]  }
      !!      save avt coef. resulting from vertical physics alone in zavt:
      !!         zavt = avt
      !!      update and save in zavt the vertical eddy viscosity coefficient:
      !!         avt = avt + wslpi^2+wslj^2
      !!      add vertical Eddy Induced advective fluxes (lk_traldf_eiv=T):
      !!         zftw = zftw + { di[aht e2u mi(wslpi)]
      !!                    +dj[aht e1v mj(wslpj)] } mk(trb)
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
      !!      Add this trend to the general trend tra :
      !!         tra = tra + difft
      !!
      !! ** Action :
      !!         Update tra arrays with the before vertical diffusion trend
      !!         Save in trtrd arrays the trends if 'key_trc_diatrd' defined
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard) correction on tracer trend loops
      !!        !  96-01  (G. Madec) statement function for e3
      !!        !  97-05  (G. Madec) vertical component of isopycnal
      !!        !  97-07  (G. Madec) geopotential diffusion in s-coord
      !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
      !!        !                     avt multiple correction
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce_trc, ONLY :   zwx => ua,  &  ! use ua, va as
                            zwy => va      ! workspace arrays

      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn       ! dummy loop indices
#if defined key_partial_steps
      INTEGER ::   iku, ikv
#endif
      REAL(wp) ::   &
         ztavg,                  &  ! temporary scalars
         zcoef0, zcoef3,         &  !    "         "
         zcoef4,                 &  !    "         "
         zbtr, zmku, zmkv,       &  !    "         "
#if defined key_trcldf_eiv
         zcoeg3,                 &  !    "         "
         zuwki, zvwki,           &  !    "         "
         zuwk, zvwk,             &  !    "         "
#endif
         ztav
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwz, zwt, ztfw             ! temporary workspace arrays
      !!---------------------------------------------------------------------

      DO jn = 1, jptra

         ! 0. Local constant initialization
         ! --------------------------------
         ztavg = 0.e0

         zwx( 1 ,:,:)=0.e0     ;     zwx(jpi,:,:)=0.e0
         zwy( 1 ,:,:)=0.e0     ;     zwy(jpi,:,:)=0.e0
         zwz( 1 ,:,:)=0.e0     ;     zwz(jpi,:,:)=0.e0
         zwt( 1 ,:,:)=0.e0     ;     zwt(jpi,:,:)=0.e0
         ztfw( 1 ,:,:)=0.e0    ;     ztfw(jpi,:,:)=0.e0

         ! I. Vertical trends associated with lateral mixing
         ! -------------------------------------------------
         !    (excluding the vertical flux proportional to dk[t] )


         ! I.1 horizontal tracer gradient
         ! ------------------------------

         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  ! i-gradient of passive tracer at ji
                  zwx (ji,jj,jk) = ( trb(ji+1,jj,jk,jn)-trb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  ! j-gradient of passive tracer at jj
                  zwy (ji,jj,jk) = ( trb(ji,jj+1,jk,jn)-trb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
#  if defined key_partial_steps
         ! partial steps correction at the bottom ocean level 
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ! last ocean level
               iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
               ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
               ! i-gradient of passive tracer
               zwx (ji,jj,iku) = gtru(ji,jj,jn)
               ! j-gradient of passive tracer
               zwy (ji,jj,ikv) = gtrv(ji,jj,jn)  
            END DO
         END DO
#endif


         ! I.2 Vertical fluxes
         ! -------------------

         ! Surface and bottom vertical fluxes set to zero
         ztfw(:,:, 1 ) = 0.e0
         ztfw(:,:,jpk) = 0.e0

         ! interior (2=<jk=<jpk-1)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)

                  zmku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
                     &           + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )

                  zmkv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
                     &           + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )

                  zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
                  zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)

                  ztfw(ji,jj,jk) = zcoef3 * (   zwx(ji  ,jj  ,jk-1) + zwx(ji-1,jj  ,jk)      &
                     &                        + zwx(ji-1,jj  ,jk-1) + zwx(ji  ,jj  ,jk)  )   &
                     &           + zcoef4 * (   zwy(ji  ,jj  ,jk-1) + zwy(ji  ,jj-1,jk)      &
                     &                        + zwy(ji  ,jj-1,jk-1) + zwy(ji  ,jj  ,jk)  )
               END DO
            END DO
         END DO

#if defined key_trcldf_eiv
         !                              ! ---------------------------------------!
         !                              ! Eddy induced vertical advective fluxes !
         !                              ! ---------------------------------------!
         zwx(:,:, 1 ) = 0.e0
         zwx(:,:,jpk) = 0.e0

         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
#   if defined key_traldf_c2d || defined key_traldf_c3d
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * fsaeitru(ji-1,jj,jk) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * fsaeitru(ji  ,jj,jk) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * fsaeitrv(ji,jj-1,jk) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * fsaeitrv(ji,jj  ,jk) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
#   else
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * e1v(ji  ,jj) * vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
                     &            * ( zuwk - zuwki + zvwk - zvwki )
#   endif
                  zwx(ji,jj,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) )

                  ztfw(ji,jj,jk) = ztfw(ji,jj,jk) + zwx(ji,jj,jk)
#   if defined key_diaeiv
                  w_trc_eiv(ji,jj,jk) = -2. * zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
#   endif
               END DO
            END DO
         END DO
#endif

         ! I.5 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  ztav = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztav
#if defined key_trc_diatrd
#   if defined key_trcldf_eiv
                  ztavg = ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) ) * zbtr
                  !  WARNING trtrd(ji,jj,jk,7) used for vertical gent velocity trend  not for damping !!!
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),7) = ztavg
#   endif
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztav - ztavg
#endif
               END DO
            END DO
         END DO

      END DO

   END SUBROUTINE trc_zdf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :             NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf_iso_vopt( kt )              ! empty routine
      WRITE(*,*) 'trc_zdf_iso_vopt: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf_iso_vopt
#endif

   !!==============================================================================
END MODULE trczdf_iso_vopt
