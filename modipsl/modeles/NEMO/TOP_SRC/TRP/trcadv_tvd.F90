MODULE trcadv_tvd
   !!==============================================================================
   !!                       ***  MODULE  trcadv_tvd  ***
   !! Ocean passive tracers:  horizontal & vertical advective trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_adv_tvd  : update the passive tracer trend with the horizontal
   !!                  and vertical advection trends using a TVD scheme
   !!   nonosc       : compute monotonic tracer fluxes by a nonoscillatory
   !!                  algorithm 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE lbclnk              ! ocean lateral boundary conditions (or mpp link)
   USE trcbbl              ! advective passive tracers in the BBL
   USE prtctl_trc      ! Print control for debbuging
#ifdef key_RIVER_INPUT
   USE rivers
#endif

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_adv_tvd    ! routine called by trcstp.F90

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_tvd.F90,v 1.12 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_adv_tvd( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_tvd  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction)
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - update tra with the now advective tracer trends
      !!             - save the trends in trtrd ('key_trc_diatrd)
      !!
      !! History :
      !!        !  95-12  (L. Mortier)  Original code
      !!        !  00-01  (H. Loukos)  adapted to ORCA 
      !!        !  00-10  (MA Foujols E.Kestenare)  include file not routine
      !!        !  00-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
      !!        !  01-07  (E. Durand G. Madec)  adaptation to ORCA config
      !!   9.0  !  02-06  (C. Ethe, G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
#if defined key_trcbbl_adv
      USE oce_trc            , zun => ua,  &  ! use ua as workspace
         &                     zvn => va      ! use va as workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
#else
      USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
                               zvn => vn,  &  !             zvn == vn
                               zwn => wn      !             zwn == wn
#endif
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER  ::   ji, jj, jk,jn           ! dummy loop indices
      REAL(wp) ::   ztra                    ! temporary scalar

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zti, ztu, ztv, ztw                ! temporary workspace

      REAL(wp) ::   &
         z2dtt, zbtr, zeu, zev, zew, z2, &  ! temporary scalar
         zfp_ui, zfp_vj, zfp_wk,         &  !    "         "
         zfm_ui, zfm_vj, zfm_wk             !    "         "

#if defined key_trc_diatrd
       REAL(wp) :: &
          zgm, zgz
#endif

      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      zti(:,:,:) = 0.e0

      IF( kt == nittrc000  .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_adv_tvd : TVD advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      IF( neuler == 0 .AND. kt == nittrc000 ) THEN
         z2=1.
      ELSE
         z2=2.
      ENDIF

#if defined key_trcbbl_adv
      ! Advective Bottom boundary layer: add the velocity
      ! -------------------------------------------------
      zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
      zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
      zwn(:,:,:) = wn (:,:,:) + w_trc_bbl(:,:,:)
#endif

      DO jn = 1, jptra

         ! 1. Bottom value : flux set to zero
         ! ---------------
         ztu(:,:,jpk) = 0.e0
         ztv(:,:,jpk) = 0.e0
         ztw(:,:,jpk) = 0.e0
         zti(:,:,jpk) = 0.e0


         ! 2. upstream advection with initial mass fluxes & intermediate update
         ! --------------------------------------------------------------------
         ! upstream tracer flux in the i and j direction
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
                  zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
                  ! upstream scheme
                  zfp_ui = zeu + ABS( zeu )
                  zfm_ui = zeu - ABS( zeu )
                  zfp_vj = zev + ABS( zev )
                  zfm_vj = zev - ABS( zev )
                  ztu(ji,jj,jk) = zfp_ui * trb(ji,jj,jk,jn) + zfm_ui * trb(ji+1,jj  ,jk,jn)
                  ztv(ji,jj,jk) = zfp_vj * trb(ji,jj,jk,jn) + zfm_vj * trb(ji  ,jj+1,jk,jn)
               END DO
            END DO
         END DO

         ! upstream tracer flux in the k direction
         ! Surface value
         IF( lk_dynspg_rl ) THEN   ! rigid lid : flux set to zero
            ztw(:,:,1) = 0.e0
         ELSE                      ! free surface 
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zew = e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,1)
                  ztw(ji,jj,1) = zew * trb(ji,jj,1,jn)
               END DO
            END DO
         ENDIF

         ! Interior value
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
                  zfp_wk = zew + ABS( zew )
                  zfm_wk = zew - ABS( zew )
                  ztw(ji,jj,jk) = zfp_wk * trb(ji,jj,jk,jn) + zfm_wk * trb(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO

         ! total advective trend
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  zti(ji,jj,jk) = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )   &
                     &              + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )   &
                     &              + ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr

#if defined key_trc_diatrd
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
                        &                          zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )                     
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
                        &                          zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )                     
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
                        &                          zbtr * ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) )
#endif
               END DO
            END DO
         END DO

#ifdef key_RIVER_INPUT
!NL#2 (same call as the routine tra_adv_tvd )
      call river_set_mask
#endif

         ! update and guess with monotonic sheme
         DO jk = 1, jpkm1
            z2dtt = z2 * rdttra(jk) * FLOAT(ndttrc)
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + zti(ji,jj,jk)
                  zti (ji,jj,jk) = ( trb(ji,jj,jk,jn) + z2dtt * zti(ji,jj,jk) ) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Lateral boundary conditions on zti, zsi   (unchanged sign)
         CALL lbc_lnk( zti, 'T', 1. )

         ! 3. antidiffusive flux : high order minus low order
         ! --------------------------------------------------
         ! antidiffusive flux on i and j
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
                  zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
                  zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
                  ztu(ji,jj,jk) = zeu * ( trn(ji,jj,jk,jn) + trn(ji+1,jj,jk,jn) ) - ztu(ji,jj,jk)
                  ztv(ji,jj,jk) = zev * ( trn(ji,jj,jk,jn) + trn(ji,jj+1,jk,jn) ) - ztv(ji,jj,jk)
               END DO
            END DO
         END DO

         ! antidiffusive flux on k
         ! Surface value
         ztw(:,:,1) = 0.

         ! Interior value
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
                  ztw(ji,jj,jk) = zew * ( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) ) - ztw(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Lateral bondary conditions
         CALL lbc_lnk( ztu, 'U', -1. )
         CALL lbc_lnk( ztv, 'V', -1. )
         CALL lbc_lnk( ztw, 'W',  1. )

         ! 4. monotonicity algorithm
         ! -------------------------
         CALL nonosc( trb(:,:,:,jn), ztu, ztv, ztw, zti, z2 )
#ifdef key_RIVER_INPUT
!NL#2 (same call as the routine tra_adv_tvd )
      call river_unset_mask
#endif


         ! 5. final trend with corrected fluxes
         ! ------------------------------------
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
#if defined key_trc_diatrd
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
                        &                          zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )                     
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
                        &                          zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )                     
                  IF ( luttrd(jn) ) &
                     trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
                        &                          zbtr * ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) )
#endif
                  tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn)   &
                     &         - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )   &
                     &           + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )   &
                     &           + ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
               END DO
            END DO
         END DO
         ! 6.0 convert the transport trend into advection trend
         ! ----------------------------------------------------
         
#if defined key_trc_diatrd
         DO jk = 1,jpk
            DO jj = 2,jpjm1
               DO  ji = 2,jpim1
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  zgm = zbtr * trn(ji,jj,jk,jn) *     &
                     &         (  zun(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)    &
                     &          - zun(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) )
                  
                  zgz = zbtr * trn(ji,jj,jk,jn) *     &
                     &         (  zvn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)    &
                     &          - zvn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) )
                  
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) + zgm
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) + zgz
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3)    &
                     &            - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
               END DO
            END DO
         END DO
         
         ! Lateral boundary conditions on trtrd:
         
         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),1), 'T', 1. )
         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),2), 'T', 1. )
         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),3), 'T', 1. )
#endif

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('tvd - adv')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_adv_tvd


   SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, prdt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE nonosc  ***
      !!     
      !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
      !!       scheme and the before field by a nonoscillatory algorithm 
      !!
      !! **  Method  :   ... ???
      !!       warning : pbef and paft must be masked, but the boundaries
      !!       conditions on the fluxes are not necessary zalezak (1979)
      !!       drange (1995) multi-dimensional forward-in-time and upstream-
      !!       in-space based differencing for fluid
      !!
      !! History :
      !!        !  97-04  (L. Mortier) Original code
      !!        !  00-02  (H. Loukos)  rewritting for opa8
      !!        !  00-10  (M.A Foujols, E. Kestenare)  lateral b.c.
      !!        !  01-03  (E. Kestenare)  add key_passivetrc
      !!        !  01-07  (E. Durand G. Madec)  adapted for T & S
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT( in ) ::   &
         prdt                               ! ???
      REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT( inout ) ::   &
         pbef,                            & ! before field
         paft,                            & ! after field
         paa,                             & ! monotonic flux in the i direction
         pbb,                             & ! monotonic flux in the j direction
         pcc                                ! monotonic flux in the k direction

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      INTEGER ::   ikm1
      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   zbetup, zbetdo
      REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt
      !!----------------------------------------------------------------------

      zbig = 1.e+40
      zrtrn = 1.e-15
      zbetup(:,:,:) = 0.e0   ;   zbetdo(:,:,:) = 0.e0

      ! Search local extrema
      ! --------------------
      ! large negative value (-zbig) inside land
      ! large negative value (-zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
      ! search maximum in neighbourhood
      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
                  &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
                  &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO
      ! large positive value (+zbig) inside land
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
      ! search minimum in neighbourhood
      DO jk = 1, jpkm1
         ikm1 = MAX(jk-1,1)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
                  &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
                  &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
                  &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
                  &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
                  &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
                  &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
            END DO
         END DO
      END DO

      ! restore masked values to zero
      pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
      paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)
 

      ! 2. Positive and negative part of fluxes and beta terms
      ! ------------------------------------------------------

      DO jk = 1, jpkm1
         z2dtt = prdt * rdttra(jk) * FLOAT(ndttrc)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! positive & negative part of the flux
               zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
               zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
                  & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
                  & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
               ! up & down beta terms
               zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / z2dtt
               zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
               zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
            END DO
         END DO
      END DO

      ! lateral boundary condition on zbetup & zbetdo   (unchanged sign)
      CALL lbc_lnk( zbetup, 'T', 1. )
      CALL lbc_lnk( zbetdo, 'T', 1. )


      ! 3. monotonic flux in the i direction, i.e. paa
      ! ----------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zc = paa(ji,jj,jk)
               IF( zc >= 0. ) THEN
                  za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
                  paa(ji,jj,jk) = za * zc
               ELSE
                  zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
                  paa(ji,jj,jk) = zb * zc
               ENDIF
            END DO
         END DO
      END DO

      ! lateral boundary condition on paa   (changed sign)
      CALL lbc_lnk( paa, 'U', -1. )


      ! 4. monotonic flux in the j direction, i.e. pbb
      ! ----------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zc = pbb(ji,jj,jk)
               IF( zc >= 0. ) THEN
                  za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
                  pbb(ji,jj,jk) = za * zc
               ELSE
                  zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
                  pbb(ji,jj,jk) = zb * zc
               ENDIF
            END DO
         END DO
      END DO

      ! lateral boundary condition on pbb   (changed sign)
      CALL lbc_lnk( pbb, 'V', -1. )


      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zc = pcc(ji,jj,jk)
               IF( zc >= 0. ) THEN
                  za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
                  pcc(ji,jj,jk) = za * zc
               ELSE
                  zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
                  pcc(ji,jj,jk) = zb * zc
               ENDIF
            END DO
         END DO
      END DO

      ! lateral boundary condition on pcc   (unchanged sign)
      CALL lbc_lnk( pcc, 'W', 1. )

   END SUBROUTINE nonosc

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_adv_tvd( kt )  
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_adv_tvd: You should not have seen this print! error?', kt
   END SUBROUTINE trc_adv_tvd
#endif

   !!======================================================================
END MODULE trcadv_tvd
