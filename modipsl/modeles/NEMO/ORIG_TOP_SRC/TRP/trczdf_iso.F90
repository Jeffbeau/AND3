MODULE trczdf_iso
   !!==============================================================================
   !!                    ***  MODULE  trczdf_iso  ***
   !! Ocean passive tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
#if defined key_passivetrc && ( defined key_ldfslp   ||   defined key_esopa )
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   trc_zdf_iso  : update the tracer trend with the vertical part of 
   !!                  the isopycnal or geopotential s-coord. operator and
   !!                  the vertical diffusion
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc          ! ocean dynamics and tracers variables
   USE trc              ! ocean passive tracers variables 
   USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
   USE trctrp_lec       ! passive tracers transport
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_zdf_iso    ! routine called by step.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      rdttrc                     ! vertical profile of 2 x tracer time-step

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_iso.F90,v 1.10 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_zdf_iso( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_iso  ***
      !!
      !! ** Purpose :
      !!     Compute the trend due to the vertical tracer diffusion inclu-
      !!     ding the vertical component of lateral mixing (only for second
      !!     order operator, for fourth order it is already computed and
      !!     add to the general trend in trcldf.F) and add it to the general
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
      !!      add vertical Eddy Induced advective fluxes ('lk_trcldf_eiv=T):
      !!         zftw = zftw + { di[aht e2u mi(wslpi)]
      !!                    +dj[aht e1v mj(wslpj)] } mk(trb)
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
      !!      Add this trend to the general trend tra :
      !!         tra = tra + difft
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers tra  is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme, t=ta.
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend tra  :
      !!         tra = tra + dz( avt dz(t) )
      !!         (tra = tra + dz( avs dz(t) ) if lk_trc_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in diamld).
      !!         avt = zavt
      !!         (avs = zavs if lk_trc_zdfddm=T )
      !!
      !!      'key_trc_diatrd' defined: trend saved for futher diagnostics.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action :
      !!         Update tra arrays with the before vertical diffusion trend
      !!         Save in trtrd arrays the trends if 'key_trc_diatrd' defined
      !!
      !! History :
      !!   7.0  !  91-11  (G. Madec)  Original code
      !!        !  92-06  (M. Imbard)  correction on tracer trend loops
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
      !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
      !!        !                     avt multiple correction
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce_trc               ,   &
         zavs => va

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn                 ! dummy loop indices
      INTEGER ::   ikst, ikenm2, ikstp1       ! temporary integers
#if defined key_partial_steps
      INTEGER ::   iku, ikv, ikv1             ! temporary integers
#endif
      REAL(wp) ::   ztra
      REAL(wp) ::   &
         ztavg,                 &  ! ???
         zcoef0, zcoef3,        &  ! ???
         zcoef4, zavi,          &  ! ???
         zbtr, zmku, zmkv,      &  !
         ztav
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwd, zws, zwi,         &  ! ???
         zwx, zwy, zwz, zwt        ! ???
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         ztfw, zdit, zdjt, zdj1t
#if defined key_trcldf_eiv   ||   defined key_esopa
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         ztfwg

      REAL(wp) ::         &
         zcoeg3,          &
         zuwk, zvwk,      &
         zuwki, zvwki
#endif
      CHARACTER (len=22) :: charout
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_zdf_iso : vertical mixing (including isopycnal component)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
#if defined key_trcldf_eiv && defined key_diaeiv
         w_trc_eiv(:,:,:) = 0.e0
#endif
      ENDIF

      ! 0.0  Local constant initialization
      ! --------------------------------
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



      ! 0.1 Save avs in zavs to recover avs in output files
      !---------------------------------------------------
      zavs(:,:,:) = fstravs(:,:,:)



      DO jn = 1, jptra

         ztavg = 0.e0

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
                  ! i-gradient of passive tracer at jj
                  zdit (ji,jk) = ( trb(ji+1,jj,jk,jn)-trb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
                  ! j-gradient of passive tracer at jj
                  zdjt (ji,jk) = ( trb(ji,jj+1,jk,jn)-trb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
                  ! j-gradient of passive tracer at jj+1
                  zdj1t(ji,jk) = ( trb(ji,jj,jk,jn)-trb(ji,jj-1,jk,jn) ) * vmask(ji,jj-1,jk)
               END DO
            END DO
#  if defined key_partial_steps
            ! partial steps correction at the bottom ocean level 
            DO ji = 1, jpim1
               ! last ocean level
               iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
               ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
               ikv1 = MIN( mbathy(ji,jj), mbathy(ji  ,jj-1) ) - 1
               ! i-gradient of of passive tracer at jj
               zdit (ji,iku) = gtru(ji,jj,jn)
               ! j-gradient of of passive tracer at jj
               zdjt (ji,ikv) = gtrv(ji,jj,jn) 
               ! j-gradient of of passive tracer at jj+1
               zdj1t(ji,ikv1)= gtrv(ji,jj-1,jn)
            END DO
#endif


            ! I.2 Vertical fluxes
            ! -------------------

            ! Surface and bottom vertical fluxes set to zero
            ztfw(:, 1 ) = 0.e0
            ztfw(:,jpk) = 0.e0

#if defined key_trcldf_eiv
            ztfwg(:, 1 ) = 0.e0
            ztfwg(:,jpk) = 0.e0
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

               END DO
            END DO


            ! I.3  update and save of avt (and avs if double diffusive mixing)
            ! ---------------------------

            DO jk = 2, jpkm1
               DO ji = 2, jpim1

                  zavi = fsahtw(ji,jj,jk)*( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)   &
                     &                     +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )

                  ! add isopycnal vertical coeff. to avs
                  fstravs(ji,jj,jk) = fstravs(ji,jj,jk) + zavi

               END DO
            END DO

#if defined key_trcldf_eiv
            !                              ! ---------------------------------------!
            !                              ! Eddy induced vertical advective fluxes !
            !                              ! ---------------------------------------!
#if defined key_traldf_c2d || defined key_traldf_c3d
            DO jk = 2, jpkm1
               DO ji = 2, jpim1
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * fsaeitru(ji-1,jj,jk) * e2u(ji-1,jj)*umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * fsaeitru(ji  ,jj,jk) * e2u(ji  ,jj)*umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * fsaeitrv(ji,jj-1,jk) * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * fsaeitrv(ji,jj  ,jk) * e1v(ji  ,jj)*vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )

                  ztfwg(ji,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) ) 
                  ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jk)

# if defined key_diaeiv
                  w_trc_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
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

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeitrw(ji,jj,jk)   &
                     &            * ( zuwk - zuwki + zvwk - zvwki )

                  ztfwg(ji,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) )
                  ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jk)

# if defined key_diaeiv
                  w_trc_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
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
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztav

#if defined key_trc_diatrd
#   if defined key_trcldf_eiv
                  ztavg = ( ztfwg(ji,jk) - ztfwg(ji,jk+1) ) * zbtr
                  !  WARNING trtrd(ji,jj,jk,6) used for vertical gent velocity trend
                  !                           not for damping !!!
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztavg
#   endif
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztav - ztavg
#endif
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - 1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

      DO jn = 1, jptra
         !                                                ! ===============
         DO jj = 2, jpjm1                                 !  Vertical slab
            !                                             ! ===============

            ! II. Vertical trend associated with the vertical physics
            ! =======================================================
            !     (including the vertical flux proportional to dk[t] associated
            !      with the lateral mixing, through the avt update)
            !     dk[ avt dk[ (t,s) ] ] diffusive trends


            ! Diagonal, inferior, superior
            ! (including the bottom boundary condition via avs masked)
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  zwi(ji,jk) = - rdttrc(jk) * fstravs(ji,jj,jk  )   &
                     /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
                  zws(ji,jk) = - rdttrc(jk) * fstravs(ji,jj,jk+1)   &
                     /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
                  zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
               END DO
            END DO

            ! Surface boudary conditions
            DO ji = 2, jpim1
               zwi(ji,1) = 0.e0
               zwd(ji,1) = 1. - zws(ji,1)
            END DO

            ! Second member construction
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  zwy(ji,jk) = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)
               END DO
            END DO

            ! Matrix inversion from the first level
            ikst = 1
#   include "zdf.matrixsolver.h90"
#if defined key_trc_diatrd
            ! Compute and save the vertical diffusive of tracers trends
#  if defined key_trc_ldfiso
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ztra = ( zwx(ji,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn) + trtrd(ji,jj,jk,ikeep(jn),6)
               END DO
            END DO
#  else
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ztra = ( zwx(ji,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn)
               END DO
            END DO
#  endif
#endif
            ! Save the masked passive tracer after in tra
            ! (c a u t i o n:  tracer not its trend, Leap-frog scheme done
            !                  it will not be done in tranxt)
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  tra(ji,jj,jk,jn) = zwx(ji,jk)  * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      END DO



      ! III. recover the avt (avs) resulting from vertical physics only
      !---------------------------------------------------------------
      fstravs(:,:,:) = zavs(:,:,:)

      
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - 2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_zdf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module               NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf_iso( kt )              ! empty routine
      WRITE(*,*) 'trc_zdf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf_iso
#endif

   !!==============================================================================
END MODULE trczdf_iso
