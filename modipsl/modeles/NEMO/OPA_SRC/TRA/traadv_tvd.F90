MODULE traadv_tvd
   !!==============================================================================
   !!                       ***  MODULE  traadv_tvd  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_adv_tvd  : update the tracer trend with the horizontal
   !!                  and vertical advection trends using a TVD scheme
   !!   nonosc       : compute monotonic tracer fluxes by a nonoscillatory
   !!                  algorithm 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE trabbl          ! Advective term of BBL
   USE lib_mpp
   USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
   USE diaptr          ! poleward transport diagnostics
   USE prtctl          ! Print control
#ifdef key_RIVER_INPUT
!!DB
   USE rivers
#endif

   IMPLICIT NONE
   PRIVATE


   !! * Accessibility
   PUBLIC tra_adv_tvd    ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_tvd.F90,v 1.12 2006/03/20 16:34:20 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_adv_tvd( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_tvd  ***
      !! 
      !! **  Purpose :   Compute the now trend due to total advection of 
      !!       tracers and add it to the general trend of tracer equations
      !!
      !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
      !!       corrected flux (monotonic correction)
      !!       note: - this advection scheme needs a leap-frog time scheme
      !!
      !! ** Action : - update (ta,sa) with the now advective tracer trends
      !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
      !!
      !! History :
      !!        !  95-12  (L. Mortier)  Original code
      !!        !  00-01  (H. Loukos)  adapted to ORCA 
      !!        !  00-10  (MA Foujols E.Kestenare)  include file not routine
      !!        !  00-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
      !!        !  01-07  (E. Durand G. Madec)  adaptation to ORCA config
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-01  (A. de Miranda, G. Madec, J.M. Molines ): advective bbl
      !!    "   !  08-04  (S. Cravatte) add the i-, j- & k- trends computation
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE trdmod_oce         , ztay => tladj,  &  ! use tladj latter
         &                     zsay => sladj,  &  ! use sladj latter
         &                     ztaz => tladi,  &  ! use ua as workspace
         &                     zsaz => sladi      ! use ua as workspace
#if defined key_trabbl_adv
      USE oce                , zun => ua,  &  ! use ua as workspace
         &                     zvn => va      ! use va as workspace

      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
#else
      USE oce                , zun => un,  &  ! When no bbl, zun == un
                               zvn => vn,  &  !             zvn == vn
                               zwn => wn      !             zwn == wn
#endif



      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      REAL(wp) ::                        &  ! temporary scalar
         ztai, ztaj, ztak,               &  !    "         "   
         zsai, zsaj, zsak                   !    "         "   
      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
         zti, ztu, ztv, ztw,             &  ! temporary workspace
         zsi, zsu, zsv, zsw,             &  !    "           "
         ztdta, ztdsa                       !    "           "
      REAL(wp) ::   &
         z2dtt, zbtr, zeu, zev, zew, z2, &  ! temporary scalar
         zfp_ui, zfp_vj, zfp_wk,         &  !    "         "
         zfm_ui, zfm_vj, zfm_wk             !    "         "
      !!----------------------------------------------------------------------
!!ylu DB -- removed the delta_adv code as it did not work for us
!      REAL(wp), DIMENSION (jpi,jpj) ::   &
!         delta_adv  ! temporary workspace

!!----------------------------------------------------------------------

      zti(:,:,:) = 0.e0   ;   zsi(:,:,:) = 0.e0

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_tvd : TVD advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2=1.
      ELSE
         z2=2.
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
         l_adv = 'tvd'
      ENDIF

#if defined key_trabbl_adv
      ! Advective Bottom boundary layer: add the velocity
      ! -------------------------------------------------
      zun(:,:,:) = un (:,:,:) - u_bbl(:,:,:)
      zvn(:,:,:) = vn (:,:,:) - v_bbl(:,:,:)
      zwn(:,:,:) = wn (:,:,:) + w_bbl(:,:,:)
#endif

      ! 1. Bottom value : flux set to zero
      ! ---------------
      ztu(:,:,jpk) = 0.e0   ;   zsu(:,:,jpk) = 0.e0
      ztv(:,:,jpk) = 0.e0   ;   zsv(:,:,jpk) = 0.e0
      ztw(:,:,jpk) = 0.e0   ;   zsw(:,:,jpk) = 0.e0
      zti(:,:,jpk) = 0.e0   ;   zsi(:,:,jpk) = 0.e0


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
               ztu(ji,jj,jk) = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
               ztv(ji,jj,jk) = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
               zsu(ji,jj,jk) = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
               zsv(ji,jj,jk) = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
            END DO
         END DO
      END DO

      ! upstream tracer flux in the k direction
      ! Surface value
      IF( lk_dynspg_rl ) THEN				! rigid lid : flux set to zero
         ztw(:,:,1) = 0.e0
         zsw(:,:,1) = 0.e0
      ELSE						! free surface
         DO jj = 1, jpj
            DO ji = 1, jpi
               zew = e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,1)
               ztw(ji,jj,1) = zew * tb(ji,jj,1)
               zsw(ji,jj,1) = zew * sb(ji,jj,1)
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
               ztw(ji,jj,jk) = zfp_wk * tb(ji,jj,jk) + zfm_wk * tb(ji,jj,jk-1)
               zsw(ji,jj,jk) = zfp_wk * sb(ji,jj,jk) + zfm_wk * sb(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! total advective trend
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )

               ! i- j- horizontal & k- vertical advective trends
               ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  ) ) * zbtr
               ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
               ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
               zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  ) ) * zbtr
               zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
               zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr

               ! total intermediate advective trends
               zti(ji,jj,jk) = ztai + ztaj + ztak
               zsi(ji,jj,jk) = zsai + zsaj + zsak
            END DO
         END DO
      END DO


     ! Save the intermediate vertical & j- horizontal advection trends
     IF( l_trdtra )   THEN
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr
                  zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr
                  ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
                  zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr
               END DO
            END DO
         END DO
      ENDIF
!
#ifdef key_RIVER_INPUT
!JC:
      call river_set_mask
#endif

      ! update and guess with monotonic sheme
      DO jk = 1, jpkm1
         z2dtt = z2 * rdttra(jk)
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ta(ji,jj,jk) =  ta(ji,jj,jk) + zti(ji,jj,jk)
               sa(ji,jj,jk) =  sa(ji,jj,jk) + zsi(ji,jj,jk)
               zti (ji,jj,jk) = ( tb(ji,jj,jk) + z2dtt * zti(ji,jj,jk) ) * tmask(ji,jj,jk)
               zsi (ji,jj,jk) = ( sb(ji,jj,jk) + z2dtt * zsi(ji,jj,jk) ) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Lateral boundary conditions on zti, zsi   (unchanged sign)
      CALL lbc_lnk( zti, 'T', 1. )
      CALL lbc_lnk( zsi, 'T', 1. )


      ! 3. antidiffusive flux : high order minus low order
      ! --------------------------------------------------
      ! antidiffusive flux on i and j
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
               zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
               ztu(ji,jj,jk) = zeu * ( tn(ji,jj,jk) + tn(ji+1,jj,jk) ) - ztu(ji,jj,jk)
               zsu(ji,jj,jk) = zeu * ( sn(ji,jj,jk) + sn(ji+1,jj,jk) ) - zsu(ji,jj,jk)
               ztv(ji,jj,jk) = zev * ( tn(ji,jj,jk) + tn(ji,jj+1,jk) ) - ztv(ji,jj,jk)
               zsv(ji,jj,jk) = zev * ( sn(ji,jj,jk) + sn(ji,jj+1,jk) ) - zsv(ji,jj,jk)
            END DO
         END DO
      END DO
      
      ! antidiffusive flux on k
      ! Surface value
      ztw(:,:,1) = 0.
      zsw(:,:,1) = 0.

      ! Interior value
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
               ztw(ji,jj,jk) = zew * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) ) - ztw(ji,jj,jk)
               zsw(ji,jj,jk) = zew * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) ) - zsw(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Lateral bondary conditions
      CALL lbc_lnk( ztu, 'U', -1. )   ;   CALL lbc_lnk( zsu, 'U', -1. )
      CALL lbc_lnk( ztv, 'V', -1. )   ;   CALL lbc_lnk( zsv, 'V', -1. )
      CALL lbc_lnk( ztw, 'W',  1. )   ;   CALL lbc_lnk( zsw, 'W',  1. )

      ! 4. monotonicity algorithm
      ! -------------------------
      CALL nonosc( tb, ztu, ztv, ztw, zti, z2 )
      CALL nonosc( sb, zsu, zsv, zsw, zsi, z2 )
!
#ifdef key_RIVER_INPUT
!JC:
      call river_unset_mask
#endif


      ! 5. final trend with corrected fluxes
      ! ------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.  
               zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
               ! i- j- horizontal & k- vertical advective trends
               ztai = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )) * zbtr
               ztaj = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )) * zbtr
               ztak = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1)) * zbtr
               zsai = - ( zsu(ji,jj,jk) - zsu(ji-1,jj  ,jk  )) * zbtr
               zsaj = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  )) * zbtr
               zsak = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1)) * zbtr

               ! add them to the general tracer trends
               ta(ji,jj,jk) = ta(ji,jj,jk) + ztai + ztaj + ztak
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsai + zsaj + zsak
            END DO
         END DO
      END DO

      ! save the advective trends for diagnostic
      ! tracers trends
      IF( l_trdtra )   THEN
         ! Compute the final vertical & j- horizontal advection trends
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
                  ztay(ji,jj,jk) = - ( ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  ) ) * zbtr   &
                     &             + ztay(ji,jj,jk) 
                  zsay(ji,jj,jk) = - ( zsv(ji,jj,jk) - zsv(ji  ,jj-1,jk  ) ) * zbtr   &
                     &             + zsay(ji,jj,jk) 
                  ztaz(ji,jj,jk) = - ( ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr   &
                     &             + ztaz(ji,jj,jk) 
                  zsaz(ji,jj,jk) = - ( zsw(ji,jj,jk) - zsw(ji  ,jj  ,jk+1) ) * zbtr   &
                     &             + zsaz(ji,jj,jk) 
               END DO
            END DO
         END DO

         ! horizontal advection: 
         ! make the difference between the new trends ta()/sa() and the 
         ! previous one ztdta()/ztdsa() to have the total advection trends 
         ! to which we substract the vertical trends ztaz()/zsaz()
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztaz(:,:,:)
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsaz(:,:,:)

         ! Add the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
         ztdta(:,:,:) = ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:)
         ztdsa(:,:,:) = ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)

         CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)

         ! vertical advection: 
         ! Substract the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
         ztaz(:,:,:) = ztaz(:,:,:) - tn(:,:,:) * hdivn(:,:,:)
         zsaz(:,:,:) = zsaz(:,:,:) - sn(:,:,:) * hdivn(:,:,:)

         CALL trd_mod(ztaz, zsaz, jpttdzad, 'TRA', kt)

      ENDIF

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=ta, clinfo1=' tvd adv  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      ! "zonal" mean advective heat and salt transport
      IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
         pht_adv(:) = ptr_vj( ztv(:,:,:) )
         pst_adv(:) = ptr_vj( zsv(:,:,:) )
      ENDIF

   END SUBROUTINE tra_adv_tvd


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
         z2dtt = prdt * rdttra(jk)
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


      ! 3. monotonic flux in the i & j direction (paa & pbb)
      ! ----------------------------------------
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
               zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
               zc = 0.5 * ( 1.e0 + SIGN( 1.e0, paa(ji,jj,jk) ) )
               paa(ji,jj,jk) = paa(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )

               za = MIN( 1.e0, zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )
               zb = MIN( 1.e0, zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
               zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pbb(ji,jj,jk) ) )
               pbb(ji,jj,jk) = pbb(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
            END DO
         END DO
      END DO


      ! monotonic flux in the k direction, i.e. pcc
      ! -------------------------------------------
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

               za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
               zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
               zc = 0.5 * ( 1.e0 + SIGN( 1.e0, pcc(ji,jj,jk) ) )
               pcc(ji,jj,jk) = pcc(ji,jj,jk) * ( zc * za + ( 1.e0 - zc) * zb )
            END DO
         END DO
      END DO

      ! lateral boundary condition on paa, pbb, pcc
      CALL lbc_lnk( paa, 'U', -1. )      ! changed sign
      CALL lbc_lnk( pbb, 'V', -1. )      ! changed sign
      CALL lbc_lnk( pcc, 'W',  1. )      ! NO changed sign

   END SUBROUTINE nonosc

   !!======================================================================
END MODULE traadv_tvd
