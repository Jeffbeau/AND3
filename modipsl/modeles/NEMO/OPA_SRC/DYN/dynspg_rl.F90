MODULE dynspg_rl
   !!======================================================================
   !!                    ***  MODULE  dynspg_rl  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
#if   defined key_dynspg_rl   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_rl'                                           rigid lid
   !!----------------------------------------------------------------------
   !!   dyn_spg_rl   : update the momentum trend with the surface pressure
   !!                  for the rigid-lid case.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE sol_oce         ! ocean elliptic solver
   USE solpcg          ! preconditionned conjugate gradient solver
   USE solsor          ! Successive Over-relaxation solver
   USE solfet          ! FETI solver
   USE solsor_e        ! Successive Over-relaxation solver with MPP optimization
   USE solisl          ! ???
   USE obc_oce         ! Lateral open boundary condition
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_rl   ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
#  include "obc_vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_rl.F90,v 1.8 2005/12/21 10:46:38 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_rl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_rl  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!      gradient for the rigid-lid case,  add it to the general trend of 
      !!      momentum equation.
      !!
      !! ** Method  :   Rigid-lid appromimation: the surface pressure gradient 
      !!      is given by:
      !!         spgu = 1/rau0 d/dx(ps) = Mu + 1/(hu e2u) dj-1(bsfd)
      !!         spgv = 1/rau0 d/dy(ps) = Mv - 1/(hv e1v) di-1(bsfd)
      !!      where (Mu,Mv) is the vertically averaged momentum trend (i.e.
      !!      the vertical ponderated sum of the general momentum trend),
      !!      and bsfd is the barotropic streamfunction trend.
      !!      The trend is computed as follows:
      !!         -1- compute the vertically averaged momentum trend (Mu,Mv)
      !!         -2- compute the barotropic streamfunction trend by solving an
      !!      ellipic equation using a diagonal preconditioned conjugate
      !!      gradient or a successive-over-relaxation method (depending
      !!      on nsolv, a namelist parameter).
      !!         -3- add to bsfd the island trends if lk_isl=T
      !!         -4- compute the after streamfunction is for further diagnos-
      !!      tics using a leap-frog scheme.
      !!         -5- add the momentum trend associated with the surface pres-
      !!      sure gradient to the general trend.
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References :
      !!      Madec et al. 1988, ocean modelling, issue 78, 1-6.
      !!
      !! History :
      !!        !  96-05  (G. Madec, M. Imbard, M. Guyon)  rewitting in 1
      !!                  routine, without pointers, and with the same matrix
      !!                  for sor and pcg, mpp exchange, and symmetric conditions
      !!        !  96-07  (A. Weaver)  Euler forward step
      !!        !  96-11  (A. Weaver)  correction to preconditioning:
      !!        !  98-02  (M. Guyon)  FETI method
      !!        !  98-05  (G. Roullet)  free surface
      !!        !  97-09  (J.-M. Molines)  Open boundaries
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( out ) ::   kindic   ! solver flag, take a negative value
      !                                    ! when the solver doesnot converge
      !! * Local declarations
      INTEGER ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::  zbsfa, zgcx, z2dt
# if defined key_obc
      INTEGER ::   ip, ii, ij
      INTEGER ::   iii, ijj, jip, jnic
      INTEGER ::   it, itm, itm2, ib, ibm, ibm2
      REAL(wp) ::   z2dtr
# endif
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_rl : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

         ! set to zero rigid-lid specific arrays
         spgu(:,:) = 0.e0      ! surface pressure gradient (i-direction) 
         spgv(:,:) = 0.e0      ! surface pressure gradient (j-direction)
      ENDIF

      ! 0. Initializations:
      ! -------------------
# if defined key_obc
      ! space index on boundary arrays
      ib = 1
      ibm = 2
      ibm2 = 3
      ! time index on boundary arrays
      it = 1
      itm = 2
      itm2 = 3
# endif

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

         ! 1. Vertically averaged momentum trend
         ! -------------------------------------
         ! initialization to zero
         spgu(:,jj) = 0.
         spgv(:,jj) = 0.
         ! vertical sum
         DO jk = 1, jpk
            DO ji = 2, jpim1
               spgu(ji,jj) = spgu(ji,jj) + ua(ji,jj,jk) * fse3u(ji,jj,jk) * umask(ji,jj,jk)
               spgv(ji,jj) = spgv(ji,jj) + va(ji,jj,jk) * fse3v(ji,jj,jk) * vmask(ji,jj,jk)
            END DO 
         END DO 
         ! divide by the depth
         spgu(:,jj) = spgu(:,jj) * hur(:,jj)
         spgv(:,jj) = spgv(:,jj) * hvr(:,jj)

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Boundary conditions on (spgu,spgv)
      CALL lbc_lnk( spgu, 'U', -1. )
      CALL lbc_lnk( spgv, 'V', -1. )

      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! 2. Barotropic streamfunction trend (bsfd)
      ! ----------------------------------

      ! Curl of the vertically averaged velocity
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            gcb(ji,jj) = -gcdprc(ji,jj)   &
                       *(  ( e2v(ji+1,jj  )*spgv(ji+1,jj  ) - e2v(ji,jj)*spgv(ji,jj) )   &
                          -( e1u(ji  ,jj+1)*spgu(ji  ,jj+1) - e1u(ji,jj)*spgu(ji,jj) )   ) 
         END DO
      END DO

# if defined key_obc
      ! Open boundary contribution
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
           gcb(ji,jj) = gcb(ji,jj) - gcdprc(ji,jj) * gcbob(ji,jj)
         END DO
      END DO
# else
      ! No open boundary contribution
# endif

      ! First guess using previous solution of the elliptic system and
      ! not bsfd since the system is solved with 0 as coastal boundary
      ! condition. Also include a swap array (gcx,gxcb)
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zgcx        = gcx(ji,jj)
            gcx (ji,jj) = 2.*zgcx - gcxb(ji,jj)
            gcxb(ji,jj) =    zgcx
         END DO
      END DO
      ! applied the lateral boundary conditions
      IF( nsolv == 4)   CALL lbc_lnk_e( gcb, c_solver_pt, 1. )   

      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! Relative precision (computation on one processor)
      rnorme = 0.e0
      rnorme = SUM( gcb(1:nlci,1:nlcj) * gcdmat(1:nlci,1:nlcj) * gcb(1:nlci,1:nlcj) * bmask(:,:) )
      IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain

      epsr = eps*eps*rnorme
      ncut = 0
      ! if rnorme is 0, the solution is 0, the solver isn't called
      IF( rnorme == 0.e0 ) THEN
         bsfd (:,:) = 0.e0
         res   = 0.e0
         niter = 0
         ncut  = 999
      ENDIF

      kindic = 0
      ! solve the bsf system  ===> solution in gcx array
      IF( ncut == 0 ) THEN
         SELECT CASE ( nsolv )
         CASE ( 1 )                    ! diagonal preconditioned conjuguate gradient
            CALL sol_pcg( kindic )
         CASE( 2 )                     ! successive-over-relaxation
            CALL sol_sor( kindic )
         CASE( 3 )                     ! FETI solver
            CALL sol_fet( kindic )
         CASE( 4 )                     ! successive-over-relaxation with extra outer halo
            CALL sol_sor_e( kindic )
         CASE DEFAULT                  ! e r r o r in nsolv namelist parameter
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' dyn_spg_rl : e r r o r, nsolv = 1, 2 ,3 or 4'
            IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~                not = ', nsolv
            nstop = nstop + 1
         END SELECT
      ENDIF


      ! bsf trend update
      ! ----------------

      bsfd(1:nlci,1:nlcj) = gcx(1:nlci,1:nlcj)

      
      ! update bsf trend with islands trend
      ! -----------------------------------

      IF( lk_isl )   CALL isl_dyn_spg         ! update bsfd


# if defined key_obc
      ! Compute bsf trend for OBC points (not masked)

      IF( lp_obc_east ) THEN
      ! compute bsf trend at the boundary from bsfeob, computed in obc_spg
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtr = 1. / rdt
      ELSE
         z2dtr = 1. / (2. * rdt )
      ENDIF
      ! (jped,jpefm1),nieob
      DO ji = fs_nie0, fs_nie1   ! vector opt.
         DO jj = nje0m1, nje1m1
            bsfd(ji,jj) = ( bsfeob(jj) - bsfb(ji,jj) ) * z2dtr
         END DO
      END DO
      ENDIF

      IF( lp_obc_west ) THEN
      ! compute bsf trend at the boundary from bsfwob, computed in obc_spg
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtr = 1. / rdt
      ELSE
         z2dtr = 1. / ( 2. * rdt )
      ENDIF
      ! (jpwd,jpwfm1),niwob
      DO ji = fs_niw0, fs_niw1   ! vector opt.
         DO jj = njw0m1, njw1m1
            bsfd(ji,jj) = ( bsfwob(jj) - bsfb(ji,jj) ) * z2dtr
         END DO
      END DO
      ENDIF

      IF( lp_obc_north ) THEN
      ! compute bsf trend at the boundary from bsfnob, computed in obc_spg
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtr = 1. / rdt
      ELSE
         z2dtr = 1. / ( 2. * rdt )
      ENDIF
      ! njnob,(jpnd,jpnfm1)
      DO jj = fs_njn0, fs_njn1   ! vector opt.
         DO ji = nin0m1, nin1m1
            bsfd(ji,jj) = ( bsfnob(ji) - bsfb(ji,jj) ) * z2dtr
         END DO
      END DO
      ENDIF

      IF( lp_obc_south ) THEN
      ! compute bsf trend at the boundary from bsfsob, computed in obc_spg
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtr = 1. / rdt
      ELSE
         z2dtr = 1. / ( 2. * rdt )
      ENDIF
      ! njsob,(jpsd,jpsfm1)
      DO jj = fs_njs0, fs_njs1   ! vector opt.
         DO ji = nis0m1, nis1m1
            bsfd(ji,jj) = ( bsfsob(ji) - bsfb(ji,jj) ) * z2dtr
         END DO
      END DO
      ENDIF

      ! compute bsf trend for isolated coastline points

      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dtr = 1. / rdt
      ELSE
         z2dtr = 1. /( 2. * rdt )
      ENDIF

      IF( nbobc > 1 ) THEN
         DO jnic = 1,nbobc - 1
            ip = mnic(0,jnic)
            DO jip = 1,ip
               ii = miic(jip,0,jnic)
               ij = mjic(jip,0,jnic)
               IF( ii >= 1 + nimpp - 1 .AND. ii <= jpi + nimpp -1 .AND. &
                   ij >= 1 + njmpp - 1 .AND. ij <= jpj + njmpp -1 ) THEN
                  iii = ii - nimpp + 1
                  ijj = ij - njmpp + 1
                  bsfd(iii,ijj) = ( bsfic(jnic) - bsfb(iii,ijj) ) * z2dtr
               ENDIF
            END DO
         END DO
      ENDIF
# endif

      ! 4. Barotropic stream function and array swap
      ! --------------------------------------------

      ! Leap-frog time scheme, time filter and array swap
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         ! Euler time stepping (first time step, starting from rest)
         z2dt = rdt
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbsfa       = bsfb(ji,jj) + z2dt * bsfd(ji,jj)
               bsfb(ji,jj) = bsfn(ji,jj)
               bsfn(ji,jj) = zbsfa 
            END DO
         END DO
      ELSE
         ! Leap-frog time stepping - Asselin filter
         z2dt = 2.*rdt
         DO jj = 1, jpj
            DO ji = 1, jpi
               zbsfa       = bsfb(ji,jj) + z2dt * bsfd(ji,jj)
               bsfb(ji,jj) = atfp * ( bsfb(ji,jj) + zbsfa ) + atfp1 * bsfn(ji,jj)
               bsfn(ji,jj) = zbsfa
            END DO
         END DO
      ENDIF

# if defined key_obc
      ! Swap of boundary arrays
      IF( lp_obc_east ) THEN
      ! (jped,jpef),nieob
      IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
         DO jj = nje0m1, nje1
            ! fields itm2 <== itm
            bebnd(jj,ib  ,itm2) = bebnd(jj,ib  ,itm)
            bebnd(jj,ibm ,itm2) = bebnd(jj,ibm ,itm)
            bebnd(jj,ibm2,itm2) = bebnd(jj,ibm2,itm)
            bebnd(jj,ib  ,itm ) = bebnd(jj,ib  ,it )
         END DO
      ELSE
         ! fields itm <== it  plus time filter at the boundary
         DO ji = fs_nie0, fs_nie1   ! vector opt.
            DO jj = nje0m1, nje1
               bebnd(jj,ib  ,itm2) = bebnd(jj,ib  ,itm)
               bebnd(jj,ibm ,itm2) = bebnd(jj,ibm ,itm)
               bebnd(jj,ibm2,itm2) = bebnd(jj,ibm2,itm)
               bebnd(jj,ib  ,itm ) = atfp * ( bebnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bebnd(jj,ib,it)
               bebnd(jj,ibm ,itm ) = bebnd(jj,ibm ,it )
               bebnd(jj,ibm2,itm ) = bebnd(jj,ibm2,it )
            END DO
         END DO
      ENDIF
      ! fields it <== now (kt+1)
      DO ji = fs_nie0, fs_nie1   ! vector opt.
         DO jj = nje0m1, nje1
            bebnd(jj,ib  ,it  ) = bsfn (ji  ,jj)
            bebnd(jj,ibm ,it  ) = bsfn (ji-1,jj)
            bebnd(jj,ibm2,it  ) = bsfn (ji-2,jj)
         END DO
      END DO
      IF( lk_mpp )   CALL mppobc( bebnd, jpjed, jpjef, jpieob, 3*3, 2, jpj )
      ENDIF

      IF( lp_obc_west ) THEN
      ! (jpwd,jpwf),niwob
      IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
         DO jj = njw0m1, njw1
            ! fields itm2 <== itm
            bwbnd(jj,ib  ,itm2) = bwbnd(jj,ib  ,itm)
            bwbnd(jj,ibm ,itm2) = bwbnd(jj,ibm ,itm)
            bwbnd(jj,ibm2,itm2) = bwbnd(jj,ibm2,itm)
            bwbnd(jj,ib  ,itm ) = bwbnd(jj,ib  ,it )
         END DO
      ELSE
         DO ji = fs_niw0, fs_niw1   ! Vector opt.
            DO jj = njw0m1, njw1
               bwbnd(jj,ib  ,itm2) = bwbnd(jj,ib  ,itm)
               bwbnd(jj,ibm ,itm2) = bwbnd(jj,ibm ,itm)
               bwbnd(jj,ibm2,itm2) = bwbnd(jj,ibm2,itm)
               ! fields itm <== it  plus time filter at the boundary
               bwbnd(jj,ib  ,itm ) = atfp * ( bwbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bwbnd(jj,ib,it)
               bwbnd(jj,ibm ,itm ) = bwbnd(jj,ibm ,it )
               bwbnd(jj,ibm2,itm ) = bwbnd(jj,ibm2,it )
            END DO
         END DO
      ENDIF
      ! fields it <== now (kt+1)
      DO ji = fs_niw0, fs_niw1   ! Vector opt.
         DO jj = njw0m1, njw1
            bwbnd(jj,ib  ,it  ) = bsfn (ji  ,jj)
            bwbnd(jj,ibm ,it  ) = bsfn (ji+1,jj)
            bwbnd(jj,ibm2,it  ) = bsfn (ji+2,jj)
         END DO
      END DO
      IF( lk_mpp )   CALL mppobc( bwbnd, jpjwd, jpjwf, jpiwob, 3*3, 2, jpj )
      ENDIF

      IF( lp_obc_north ) THEN
   ! njnob,(jpnd,jpnf)
      IF( kt < nit000 + 3 .AND. .NOT.ln_rstart ) THEN
         DO ji = nin0m1, nin1
            ! fields itm2 <== itm
            bnbnd(ji,ib  ,itm2) = bnbnd(ji,ib  ,itm)
            bnbnd(ji,ibm ,itm2) = bnbnd(ji,ibm ,itm)
            bnbnd(ji,ibm2,itm2) = bnbnd(ji,ibm2,itm)
            bnbnd(ji,ib  ,itm ) = bnbnd(ji,ib  ,it )
         END DO
      ELSE
         DO jj = fs_njn0, fs_njn1   ! Vector opt.
            DO ji = nin0m1, nin1
               bnbnd(ji,ib  ,itm2) = bnbnd(ji,ib  ,itm)
               bnbnd(ji,ibm ,itm2) = bnbnd(ji,ibm ,itm)
               bnbnd(ji,ibm2,itm2) = bnbnd(ji,ibm2,itm)
               ! fields itm <== it  plus time filter at the boundary
               bnbnd(jj,ib  ,itm ) = atfp * ( bnbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bnbnd(jj,ib,it)
               bnbnd(ji,ibm ,itm ) = bnbnd(ji,ibm ,it )
               bnbnd(ji,ibm2,itm ) = bnbnd(ji,ibm2,it )
            END DO
          END DO
      ENDIF
      ! fields it <== now (kt+1)
      DO jj = fs_njn0, fs_njn1   ! Vector opt.
         DO ji = nin0m1, nin1
            bnbnd(ji,ib  ,it  ) = bsfn (ji,jj  )
            bnbnd(ji,ibm ,it  ) = bsfn (ji,jj-1)
            bnbnd(ji,ibm2,it  ) = bsfn (ji,jj-2)
         END DO
      END DO
      IF( lk_mpp )   CALL mppobc( bnbnd, jpind, jpinf, jpjnob, 3*3, 1, jpi )
      ENDIF

      IF( lp_obc_south ) THEN
         ! njsob,(jpsd,jpsf)
         IF( kt < nit000+3 .AND. .NOT.ln_rstart ) THEN
            DO ji = nis0m1, nis1
               ! fields itm2 <== itm
               bsbnd(ji,ib  ,itm2) = bsbnd(ji,ib  ,itm)
               bsbnd(ji,ibm ,itm2) = bsbnd(ji,ibm ,itm)
               bsbnd(ji,ibm2,itm2) = bsbnd(ji,ibm2,itm)
               bsbnd(ji,ib  ,itm ) = bsbnd(ji,ib  ,it )
            END DO
         ELSE
            DO jj = fs_njs0, fs_njs1   ! vector opt.
               DO ji = nis0m1, nis1
                  bsbnd(ji,ib  ,itm2) = bsbnd(ji,ib  ,itm)
                  bsbnd(ji,ibm ,itm2) = bsbnd(ji,ibm ,itm)
                  bsbnd(ji,ibm2,itm2) = bsbnd(ji,ibm2,itm)
                  ! fields itm <== it  plus time filter at the boundary
                  bsbnd(jj,ib  ,itm ) = atfp * ( bsbnd(jj,ib,itm) + bsfn(ji,jj) ) + atfp1 * bsbnd(jj,ib,it)
                  bsbnd(ji,ibm ,itm ) = bsbnd(ji,ibm ,it )
                  bsbnd(ji,ibm2,itm ) = bsbnd(ji,ibm2,it )
               END DO
            END DO
         ENDIF
         DO jj = fs_njs0, fs_njs1   ! vector opt.
            DO ji = nis0m1, nis1 
               ! fields it <== now (kt+1)
               bsbnd(ji,ib  ,it  ) = bsfn (ji,jj  )
               bsbnd(ji,ibm ,it  ) = bsfn (ji,jj+1)
               bsbnd(ji,ibm2,it  ) = bsfn (ji,jj+2)
            END DO
         END DO
         IF( lk_mpp )   CALL mppobc( bsbnd, jpisd, jpisf, jpjsob, 3*3, 1, jpi )
      ENDIF
# endif
      !
      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      
      !  add the surface pressure trend to the general trend
      ! -----------------------------------------------------
      
      DO jj = 2, jpjm1

         ! update the surface pressure gradient with the barotropic trend
         DO ji = 2, jpim1
            spgu(ji,jj) = spgu(ji,jj) + hur(ji,jj) / e2u(ji,jj) * ( bsfd(ji,jj) - bsfd(ji  ,jj-1) )
            spgv(ji,jj) = spgv(ji,jj) - hvr(ji,jj) / e1v(ji,jj) * ( bsfd(ji,jj) - bsfd(ji-1,jj  ) )
         END DO
         ! add the surface pressure gradient trend to the general trend
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ua(ji,jj,jk) - spgu(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - spgv(ji,jj)
            END DO
         END DO

      END DO

   END SUBROUTINE dyn_spg_rl

#else
   !!----------------------------------------------------------------------
   !!   'key_dynspg_rl'                                        NO rigid lid
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_spg_rl( kt, kindic )       ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_spg_rl: You should not have seen this print! error?', kt, kindic
   END SUBROUTINE dyn_spg_rl
#endif
   
   !!======================================================================
END MODULE dynspg_rl
