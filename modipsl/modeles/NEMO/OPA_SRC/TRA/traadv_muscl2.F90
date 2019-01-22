MODULE traadv_muscl2
   !!==============================================================================
   !!                       ***  MODULE  traadv_muscl2  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_adv_muscl2 : update the tracer trend with the horizontal
   !!                    and vertical advection trends using MUSCL2 scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE trabbl          ! tracers: bottom boundary layer
   USE lib_mpp
   USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
   USE diaptr          ! poleward transport diagnostics
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC tra_adv_muscl2        ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_muscl2.F90,v 1.10 2005/12/28 09:25:10 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_adv_muscl2( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_adv_muscl2  ***
      !!
      !! ** Purpose :   Compute the now trend due to total advection of T and 
      !!      S using a MUSCL scheme (Monotone Upstream-centered Scheme for
      !!      Conservation Laws) and add it to the general tracer trend.
      !!
      !! ** Method  : MUSCL scheme plus centered scheme at ocean boundaries
      !!
      !! ** Action  : - update (ta,sa) with the now advective tracer trends
      !!              - save trends in (ztdta,ztdsa) ('key_trdtra')
      !!
      !! References :                
      !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
      !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
      !!
      !! History :
      !!        !  06-00  (A.Estublier)  for passive tracers
      !!        !  01-08  (E.Durand G.Madec)  adapted for T & S
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * modules used
#if defined key_trabbl_adv
      USE oce                , zun => ua,  &  ! use ua as workspace
         &                     zvn => va      ! use va as workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
#else
      USE oce                , zun => un,  &  ! When no bbl, zun == un
                               zvn => vn,  &  !              zvn == vn
                               zwn => wn      !              zwn == wn
#endif

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   &
         zu, zv, zw, zeu, zev,           &  
         zew, zbtr, zstep,               &
         z0u, z0v, z0w,                  &
         zzt1, zzt2, zalpha,             &
         zzs1, zzs2, z2,                 &
         zta, zsa
      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
         zt1, zt2, ztp1, ztp2,   &
         zs1, zs2, zsp1, zsp2,   &
         ztdta, ztdsa
      !!----------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_muscl2 : MUSCL2 advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
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
         l_adv = 'mu2'
      ENDIF

#if defined key_trabbl_adv
      ! Advective bottom boundary layer
      ! -------------------------------
      zun(:,:,:) = un (:,:,:) - u_bbl(:,:,:)
      zvn(:,:,:) = vn (:,:,:) - v_bbl(:,:,:)
      zwn(:,:,:) = wn (:,:,:) + w_bbl( :,:,:)
#endif


      ! I. Horizontal advective fluxes
      ! ------------------------------

      ! first guess of the slopes
      ! interior values
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1      
            DO ji = 1, fs_jpim1   ! vector opt.
               zt1(ji,jj,jk) = umask(ji,jj,jk) * ( tb(ji+1,jj,jk) - tb(ji,jj,jk) )
               zs1(ji,jj,jk) = umask(ji,jj,jk) * ( sb(ji+1,jj,jk) - sb(ji,jj,jk) )
               zt2(ji,jj,jk) = vmask(ji,jj,jk) * ( tb(ji,jj+1,jk) - tb(ji,jj,jk) )
               zs2(ji,jj,jk) = vmask(ji,jj,jk) * ( sb(ji,jj+1,jk) - sb(ji,jj,jk) )
            END DO
         END DO
      END DO
      ! bottom values
      zt1(:,:,jpk) = 0.e0    ;    zt2(:,:,jpk) = 0.e0
      zs1(:,:,jpk) = 0.e0    ;    zs2(:,:,jpk) = 0.e0

      ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
      CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. )
      CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )

      ! Slopes
      ! interior values
      DO jk = 1, jpkm1
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
                  &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
               zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji-1,jj  ,jk) )   &
                  &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji-1,jj  ,jk) ) )
               ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
                  &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )
               zsp2(ji,jj,jk) =                    ( zs2(ji,jj,jk) + zs2(ji  ,jj-1,jk) )   &
                  &           * ( 0.25 + SIGN( 0.25, zs2(ji,jj,jk) * zs2(ji  ,jj-1,jk) ) )
            END DO
         END DO
      END DO
      ! bottom values
      ztp1(:,:,jpk) = 0.e0    ;    ztp2(:,:,jpk) = 0.e0
      zsp1(:,:,jpk) = 0.e0    ;    zsp2(:,:,jpk) = 0.e0

      ! Slopes limitation
      DO jk = 1, jpkm1
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
                  &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
                  &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
                  &                  2.*ABS( zt1 (ji  ,jj,jk) ) )
               zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
                  &           * MIN(    ABS( zsp1(ji  ,jj,jk) ),   &
                  &                  2.*ABS( zs1 (ji-1,jj,jk) ),   &
                  &                  2.*ABS( zs1 (ji  ,jj,jk) ) )
               ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
                  &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
                  &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
                  &                  2.*ABS( zt2 (ji,jj  ,jk) ) )
               zsp2(ji,jj,jk) = SIGN( 1., zsp2(ji,jj,jk) )   &
                  &           * MIN(    ABS( zsp2(ji,jj  ,jk) ),   &
                  &                  2.*ABS( zs2 (ji,jj-1,jk) ),   &
                  &                  2.*ABS( zs2 (ji,jj  ,jk) ) )
            END DO
         END DO
      END DO        

      ! Advection terms
      ! interior values
      DO jk = 1, jpkm1
         zstep  = z2 * rdttra(jk)
         DO jj = 2, jpjm1      
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! volume fluxes
#if defined key_s_coord || defined key_partial_steps
               zeu = e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
               zev = e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
#else
               zeu = e2u(ji,jj) * zun(ji,jj,jk)
               zev = e1v(ji,jj) * zvn(ji,jj,jk)
#endif
               ! MUSCL fluxes
               z0u = SIGN( 0.5, zun(ji,jj,jk) )            
               zalpha = 0.5 - z0u
               zu  = z0u - 0.5 * zun(ji,jj,jk) * zstep / e1u(ji,jj)
               zzt1 = tb(ji+1,jj,jk) + zu*ztp1(ji+1,jj,jk)
               zzt2 = tb(ji  ,jj,jk) + zu*ztp1(ji  ,jj,jk)
               zzs1 = sb(ji+1,jj,jk) + zu*zsp1(ji+1,jj,jk)
               zzs2 = sb(ji  ,jj,jk) + zu*zsp1(ji  ,jj,jk)
               zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
               zs1(ji,jj,jk) = zeu * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )

               z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
               zalpha = 0.5 - z0v
               zv  = z0v - 0.5 * zvn(ji,jj,jk) * zstep / e2v(ji,jj)
               zzt1 = tb(ji,jj+1,jk) + zv*ztp2(ji,jj+1,jk)
               zzt2 = tb(ji,jj  ,jk) + zv*ztp2(ji,jj  ,jk)
               zzs1 = sb(ji,jj+1,jk) + zv*zsp2(ji,jj+1,jk)
               zzs2 = sb(ji,jj  ,jk) + zv*zsp2(ji,jj  ,jk)
               zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
               zs2(ji,jj,jk) = zev * ( zalpha * zzs1 + (1.-zalpha) * zzs2 )
            END DO
         END DO
      END DO

      !!!!  centered scheme at lateral b.C. if off-shore velocity
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               zev = e1v(ji,jj) * fse3v(ji,jj,jk)
               IF( umask(ji,jj,jk) == 0. ) THEN
                  IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
                     zt1(ji+1,jj,jk) = e2u(ji+1,jj)* fse3u(ji+1,jj,jk)   &
                        &            * zun(ji+1,jj,jk) * ( tb(ji+1,jj,jk) + tb(ji+2,jj,jk) ) * 0.5
                     zs1(ji+1,jj,jk) = e2u(ji+1,jj)* fse3u(ji+1,jj,jk)   &
                        &            * zun(ji+1,jj,jk) * ( sb(ji+1,jj,jk) + sb(ji+2,jj,jk) ) * 0.5
                  ENDIF
                  IF( zun(ji-1,jj,jk) < 0. ) THEN
                     zt1(ji-1,jj,jk) = e2u(ji-1,jj)* fse3u(ji-1,jj,jk)   &
                        &            * zun(ji-1,jj,jk) * ( tb(ji-1,jj,jk) + tb(ji  ,jj,jk) ) * 0.5
                     zs1(ji-1,jj,jk) = e2u(ji-1,jj)* fse3u(ji-1,jj,jk)   &
                        &            * zun(ji-1,jj,jk) * ( sb(ji-1,jj,jk) + sb(ji  ,jj,jk) ) * 0.5
                  ENDIF
               ENDIF
               IF( vmask(ji,jj,jk) == 0. ) THEN
                  IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
                     zt2(ji,jj+1,jk) = e1v(ji,jj+1) * fse3v(ji,jj+1,jk)   &
                        &            * zvn(ji,jj+1,jk) * ( tb(ji,jj+1,jk) + tb(ji,jj+2,jk) ) * 0.5
                     zs2(ji,jj+1,jk) = e1v(ji,jj+1) * fse3v(ji,jj+1,jk)   &
                        &            * zvn(ji,jj+1,jk) * ( sb(ji,jj+1,jk) + sb(ji,jj+2,jk) ) * 0.5
                  ENDIF
                  IF( zvn(ji,jj-1,jk) < 0. ) THEN
                     zt2(ji,jj-1,jk) = e1v(ji,jj-1)* fse3v(ji,jj-1,jk)   &
                        &            * zvn(ji,jj-1,jk) * ( tb(ji,jj-1,jk) + tb(ji  ,jj,jk) ) * 0.5
                     zs2(ji,jj-1,jk) = e1v(ji,jj-1)* fse3v(ji,jj-1,jk)   &
                        &            * zvn(ji,jj-1,jk) * ( sb(ji,jj-1,jk) + sb(ji  ,jj,jk) ) * 0.5
                  ENDIF
               ENDIF

#else
               IF( umask(ji,jj,jk) == 0. ) THEN
                  IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
                     zt1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( tb(ji+1,jj,jk) + tb(ji+2,jj,jk) ) * 0.5
                     zs1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( sb(ji+1,jj,jk) + sb(ji+2,jj,jk) ) * 0.5
                  ENDIF
                  IF( zun(ji-1,jj,jk) < 0. ) THEN
                     zt1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( tb(ji-1,jj,jk) + tb(ji  ,jj,jk) ) * 0.5
                     zs1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( sb(ji-1,jj,jk) + sb(ji  ,jj,jk) ) * 0.5
                  ENDIF
               ENDIF
               IF( vmask(ji,jj,jk) == 0. ) THEN
                  IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
                     zt2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( tb(ji,jj+1,jk) + tb(ji,jj+2,jk) ) * 0.5
                     zs2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( sb(ji,jj+1,jk) + sb(ji,jj+2,jk) ) * 0.5
                  ENDIF
                  IF( zvn(ji,jj-1,jk) < 0. ) THEN
                     zt2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( tb(ji,jj-1,jk) + tb(ji  ,jj,jk) ) * 0.5
                     zs2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( sb(ji,jj-1,jk) + sb(ji  ,jj,jk) ) * 0.5
                  ENDIF
               ENDIF
#endif
            END DO
         END DO
      END DO

      ! lateral boundary conditions on zt1, zt2 ; zs1, zs2   (changed sign)
      CALL lbc_lnk( zt1, 'U', -1. )   ;   CALL lbc_lnk( zs1, 'U', -1. ) 
      CALL lbc_lnk( zt2, 'V', -1. )   ;   CALL lbc_lnk( zs2, 'V', -1. )

      ! Save MUSCL fluxes to compute i- & j- horizontal 
      ! advection trends in the MLD
      IF( l_trdtra )   THEN
         ! save i- terms
         tladi(:,:,:) = zt1(:,:,:) 
         sladi(:,:,:) = zs1(:,:,:) 
         ! save j- terms
         tladj(:,:,:) = zt2(:,:,:) 
         sladj(:,:,:) = zs2(:,:,:) 
      ENDIF

      ! Compute & add the horizontal advective trend

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1      
            DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord || defined key_partial_steps
               zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
#else
               zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
#endif
               ! horizontal advective trends
               zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
                  &           + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
               zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji-1,jj  ,jk  )   &
                  &           + zs2(ji,jj,jk) - zs2(ji  ,jj-1,jk  ) ) 
               ! add it to the general tracer trends
               ta(ji,jj,jk) = ta(ji,jj,jk) + zta
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
            END DO
        END DO
      END DO        

      ! Save the horizontal advective trends for diagnostic

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

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 had  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      ! "zonal" mean advective heat and salt transport
      IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
# if defined key_s_coord || defined key_partial_steps
         pht_adv(:) = ptr_vj( zt2(:,:,:) )
         pst_adv(:) = ptr_vj( zs2(:,:,:) )
# else
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                 zt2(ji,jj,jk) = zt2(ji,jj,jk) * fse3v(ji,jj,jk)
                 zs2(ji,jj,jk) = zs2(ji,jj,jk) * fse3v(ji,jj,jk)
               END DO
            END DO
         END DO
         pht_adv(:) = ptr_vj( zt2(:,:,:) )
         pst_adv(:) = ptr_vj( zs2(:,:,:) )
# endif
      ENDIF

      ! II. Vertical advective fluxes
      ! -----------------------------
      
      ! First guess of the slope
      ! interior values
      DO jk = 2, jpkm1
         zt1(:,:,jk) = tmask(:,:,jk) * ( tb(:,:,jk-1) - tb(:,:,jk) )
         zs1(:,:,jk) = tmask(:,:,jk) * ( sb(:,:,jk-1) - sb(:,:,jk) )
      END DO
      ! surface & bottom boundary conditions
      zt1 (:,:, 1 ) = 0.e0    ;    zt1 (:,:,jpk) = 0.e0
      zs1 (:,:, 1 ) = 0.e0    ;    zs1 (:,:,jpk) = 0.e0

      ! Slopes
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
                  &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
               zsp1(ji,jj,jk) =                    ( zs1(ji,jj,jk) + zs1(ji,jj,jk+1) )   &
                  &           * ( 0.25 + SIGN( 0.25, zs1(ji,jj,jk) * zs1(ji,jj,jk+1) ) )
            END DO
         END DO
      END DO

      ! Slopes limitation
      ! interior values
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
                  &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
                  &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
                  &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
               zsp1(ji,jj,jk) = SIGN( 1., zsp1(ji,jj,jk) )   &
                  &           * MIN(    ABS( zsp1(ji,jj,jk  ) ),   &
                  &                  2.*ABS( zs1 (ji,jj,jk+1) ),   &
                  &                  2.*ABS( zs1 (ji,jj,jk  ) ) )
            END DO
         END DO
      END DO
      ! surface values
      ztp1(:,:,1) = 0.e0
      zsp1(:,:,1) = 0.e0

      ! vertical advective flux
      ! interior values
      DO jk = 1, jpkm1
         zstep  = z2 * rdttra(jk)
         DO jj = 2, jpjm1      
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zew = zwn(ji,jj,jk+1)
               z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
               zalpha = 0.5 + z0w
               zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*zstep / fse3w(ji,jj,jk+1)
               zzt1 = tb(ji,jj,jk+1) + zw*ztp1(ji,jj,jk+1)
               zzt2 = tb(ji,jj,jk  ) + zw*ztp1(ji,jj,jk  )
               zzs1 = sb(ji,jj,jk+1) + zw*zsp1(ji,jj,jk+1)
               zzs2 = sb(ji,jj,jk  ) + zw*zsp1(ji,jj,jk  )
               zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
               zs1(ji,jj,jk+1) = zew * ( zalpha * zzs1 + (1.-zalpha)*zzs2 )
            END DO
         END DO
      END DO
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               IF( tmask(ji,jj,jk+1) == 0. ) THEN
                  IF( zwn(ji,jj,jk) > 0. ) THEN
                     zt1(ji,jj,jk) = zwn(ji,jj,jk) * ( tb(ji,jj,jk-1) + tb(ji,jj,jk) ) * 0.5
                     zs1(ji,jj,jk) = zwn(ji,jj,jk) * ( sb(ji,jj,jk-1) + sb(ji,jj,jk) ) * 0.5
                  ENDIF
               ENDIF
            END DO
         END DO
      END DO

      ! surface values
      IF( lk_dynspg_rl ) THEN                           ! rigid lid : flux set to zero
         zt1(:,:, 1 ) = 0.e0
         zs1(:,:, 1 ) = 0.e0
      ELSE                                              ! free surface
         zt1(:,:, 1 ) = zwn(:,:,1) * tb(:,:,1)
         zs1(:,:, 1 ) = zwn(:,:,1) * sb(:,:,1)
      ENDIF

      ! bottom values
      zt1(:,:,jpk) = 0.e0
      zs1(:,:,jpk) = 0.e0


      ! Compute & add the vertical advective trend

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1      
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr = 1. / fse3t(ji,jj,jk)
               ! horizontal advective trends
               zta = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
               zsa = - zbtr * ( zs1(ji,jj,jk) - zs1(ji,jj,jk+1) )
               ! add it to the general tracer trends
               ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
               sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
            END DO
         END DO
      END DO

      ! Save the vertical advective trends for diagnostic

      IF( l_trdtra )   THEN
         ! Recompute the vertical advection zta & zsa trends computed 
         ! at the step 2. above in making the difference between the new 
         ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
         ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)

         CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN
         CALL prt_ctl(tab3d_1=ta, clinfo1=' muscl2 zad  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')

      ENDIF

   END SUBROUTINE tra_adv_muscl2

   !!======================================================================
END MODULE traadv_muscl2
