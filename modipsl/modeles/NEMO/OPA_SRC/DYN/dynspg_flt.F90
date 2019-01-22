MODULE dynspg_flt
   !!======================================================================
   !!                   ***  MODULE  dynspg_flt  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
#if ( defined key_dynspg_flt && ! defined key_autotasking )   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_flt'                              filtered free surface
   !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
   !!----------------------------------------------------------------------
   !!   dyn_spg_flt  : update the momentum trend with the surface pressure
   !!                  gradient in the filtered free surface case with
   !!                  vector optimization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE zdf_oce         ! ocean vertical physics
   USE phycst          ! physical constants
   USE ocesbc          ! ocean surface boundary condition
   USE flxrnf          ! ocean runoffs
   USE sol_oce         ! ocean elliptic solver
   USE solpcg          ! preconditionned conjugate gradient solver
   USE solsor          ! Successive Over-relaxation solver
   USE solfet          ! FETI solver
   USE solsor_e        ! Successive Over-relaxation solver with MPP optimization
   USE obc_oce         ! Lateral open boundary condition
   USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
   USE obcvol          ! ocean open boundary condition (obc_vol routines)
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE cla_dynspg      ! cross land advection
   USE prtctl          ! Print control
   USE in_out_manager  ! I/O manager
   USE solmat          ! matrix construction for elliptic solvers
   USE agrif_opa_interp
#ifdef key_RIVER_INPUT
!JC:
   USE rivers
#endif


   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_flt  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_flt.F90,v 1.4 2006/03/20 17:25:49 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_flt( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_flt  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure 
      !!      gradient in case of filtered free surface formulation  and add
      !!      it to the general trend of momentum equation.
      !!
      !! ** Method  :   Filtered free surface formulation. The surface
      !!      pressure gradient is given by:
      !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( sshn + rnu btda )
      !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( sshn + rnu btda )
      !!      where sshn is the free surface elevation and btda is the after
      !!      of the free surface elevation
      !!       -1- compute the after sea surface elevation from the kinematic
      !!      surface boundary condition:
      !!              zssha = sshb + 2 rdt ( wn - emp )
      !!           Time filter applied on now sea surface elevation to avoid 
      !!      the divergence of two consecutive time-steps and swap of free
      !!      surface arrays to start the next time step:
      !!              sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
      !!              sshn = zssha
      !!       -2- evaluate the surface presure trend (including the addi-
      !!      tional force) in three steps:
      !!        a- compute the right hand side of the elliptic equation:
      !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
      !!         where (spgu,spgv) are given by:
      !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
      !!                 - grav 2 rdt hu /e1u di[sshn + emp]
      !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
      !!                 - grav 2 rdt hv /e2v dj[sshn + emp]
      !!         and define the first guess from previous computation :
      !!            zbtd = btda
      !!            btda = 2 zbtd - btdb
      !!            btdb = zbtd
      !!        b- compute the relative accuracy to be reached by the
      !!         iterative solver
      !!        c- apply the solver by a call to sol... routine
      !!       -3- compute and add the free surface pressure gradient inclu-
      !!      ding the additional force used to stabilize the equation.
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References :
      !!      Roullet and Madec 1999, JGR.
      !!
      !! History :
      !!        !  98-05 (G. Roullet)  Original code
      !!        !  98-10 (G. Madec, M. Imbard)  release 8.2
      !!   8.5  !  02-08 (G. Madec)  F90: Free form and module
      !!        !  02-11 (C. Talandier, A-M Treguier) Open boundaries
      !!   9.0  !  04-08 (C. Talandier) New trends organization
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in )  ::   kt         ! ocean time-step index
      INTEGER, INTENT( out ) ::   kindic     ! solver convergence flag
                                             ! if the solver doesn t converge
                                             ! the flag is < 0
      !! * Local declarations
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::                         & 
         z2dt, z2dtg, zraur, znugdt,      &  ! temporary scalars
         znurau, zssha, zgcb, zbtd,       &  !   "          "
         ztdgu, ztdgv                        !   "          "
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_flt : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   (free surface constant volume case)'
       
         ! set to zero free surface specific arrays
         spgu(:,:) = 0.e0                     ! surface pressure gradient (i-direction)
         spgv(:,:) = 0.e0                     ! surface pressure gradient (j-direction)
      ENDIF

      ! Local constant initialization
      z2dt = 2. * rdt                                       ! time step: leap-frog
      IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt       ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000+1 ) CALL sol_mat(kt)
      z2dtg  = grav * z2dt
      zraur  = 1. / rauw
      znugdt =  rnu * grav * z2dt
      znurau =  znugdt * zraur

      ! Surface pressure gradient (now)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
            spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
         END DO 
      END DO 

      ! Add the surface pressure trend to the general trend
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
            END DO
         END DO
      END DO

      ! Evaluate the masked next velocity (effect of the additional force not included)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
               va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO

#if defined key_obc
      ! Update velocities on each open boundary with the radiation algorithm
      CALL obc_dyn( kt )
      ! Correction of the barotropic componant velocity to control the volume of the system
      CALL obc_vol( kt )
#endif
#ifdef key_RIVER_INPUT
!JC:
   call riv_dyna( kt ) !JC: Same as FD  for ua/va velocities river points
#endif

#if defined key_agrif
      ! Update velocities on each coarse/fine interfaces

      CALL Agrif_dyn( kt )
#endif

!!DB: deleted ORCA
!#if defined key_orca_r2
!      IF( n_cla == 1 )   CALL dyn_spg_cla( kt )      ! Cross Land Advection (update (ua,va))
!#endif

      ! compute the next vertically averaged velocity (effect of the additional force not included)
      ! ---------------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            spgu(ji,jj) = 0.e0
            spgv(ji,jj) = 0.e0
         END DO
      END DO

      ! vertical sum
!CDIR NOLOOPCHG
      IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
         DO jk = 1, jpkm1
            DO ji = 1, jpij
               spgu(ji,1) = spgu(ji,1) + fse3u(ji,1,jk) * ua(ji,1,jk)
               spgv(ji,1) = spgv(ji,1) + fse3v(ji,1,jk) * va(ji,1,jk)
            END DO
         END DO
      ELSE                        ! No  vector opt.
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  spgu(ji,jj) = spgu(ji,jj) + fse3u(ji,jj,jk) * ua(ji,jj,jk)
                  spgv(ji,jj) = spgv(ji,jj) + fse3v(ji,jj,jk) * va(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! transport: multiplied by the horizontal scale factor
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            spgu(ji,jj) = spgu(ji,jj) * e2u(ji,jj)
            spgv(ji,jj) = spgv(ji,jj) * e1v(ji,jj)
         END DO
      END DO

      ! Boundary conditions on (spgu,spgv)
      CALL lbc_lnk( spgu, 'U', -1. )
      CALL lbc_lnk( spgv, 'V', -1. )

      ! Right hand side of the elliptic equation and first guess
      ! -----------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! Divergence of the after vertically averaged velocity
            zgcb =  spgu(ji,jj) - spgu(ji-1,jj)   &
                  + spgv(ji,jj) - spgv(ji,jj-1)
            gcb(ji,jj) = gcdprc(ji,jj) * zgcb
            ! First guess of the after barotropic transport divergence
            zbtd = gcx(ji,jj)
            gcx (ji,jj) = 2. * zbtd   - gcxb(ji,jj)
            gcxb(ji,jj) =      zbtd
         END DO
      END DO
      ! applied the lateral boundary conditions
      IF( nsolv == 4 )   CALL lbc_lnk_e( gcb, c_solver_pt, 1. )   

#if defined key_agrif
      IF( .NOT. AGRIF_ROOT() ) THEN
         ! add contribution of gradient of after barotropic transport divergence 
         IF( (nbondi == -1) .OR. (nbondi == 2) ) gcb(3,:) = gcb(3,:) &
            &            -znugdt * z2dt*laplacu(2,:)*gcdprc(3,:)*hu(2,:)*e2u(2,:)
         IF( (nbondi ==  1) .OR. (nbondi == 2) )  gcb(nlci-2,:) = gcb(nlci-2,:) &
            &           +znugdt * z2dt*laplacu(nlci-2,:)*gcdprc(nlci-2,:)*hu(nlci-2,:)*e2u(nlci-2,:)
         IF( (nbondj == -1) .OR. (nbondj == 2) ) gcb(:,3) = gcb(:,3) &
            &           -znugdt * z2dt*laplacv(:,2)*gcdprc(:,3)*hv(:,2)*e1v(:,2)
         IF( (nbondj ==  1) .OR. (nbondj == 2) )  gcb(:,nlcj-2) = gcb(:,nlcj-2) &
            &           +znugdt * z2dt*laplacv(:,nlcj-2)*gcdprc(:,nlcj-2)*hv(:,nlcj-2)*e1v(:,nlcj-2)
      ENDIF
#endif


      ! Relative precision (computation on one processor)
      ! ------------------
      rnorme =0.
      rnorme = SUM( gcb(1:jpi,1:jpj) * gcdmat(1:jpi,1:jpj) * gcb(1:jpi,1:jpj) * bmask(:,:) )
      IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain

      epsr = eps * eps * rnorme
      ncut = 0
      ! if rnorme is 0, the solution is 0, the solver isn't called
      IF( rnorme == 0.e0 ) THEN
         gcx(:,:) = 0.e0
         res   = 0.e0
         niter = 0
         ncut  = 999
      ENDIF

      ! Evaluate the next transport divergence
      ! --------------------------------------
      !    Iterarive solver for the elliptic equation (except IF sol.=0)
      !    (output in gcx with boundary conditions applied)
      kindic = 0
      IF( ncut == 0 ) THEN
         IF( nsolv == 1 ) THEN         ! diagonal preconditioned conjuguate gradient
            CALL sol_pcg( kindic )
         ELSEIF( nsolv == 2 ) THEN     ! successive-over-relaxation
            CALL sol_sor( kindic )
         ELSEIF( nsolv == 3 ) THEN     ! FETI solver
            CALL sol_fet( kindic )
         ELSEIF( nsolv == 4 ) THEN     ! successive-over-relaxation with extra outer halo
            CALL sol_sor_e( kindic )
         ELSE                          ! e r r o r in nsolv namelist parameter
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' dyn_spg_flt : e r r o r, nsolv = 1, 2 ,3 or 4'
            IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~                not = ', nsolv
            nstop = nstop + 1
         ENDIF
      ENDIF

      ! Transport divergence gradient multiplied by z2dt
      ! --------------------------------------------====
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            ! trend of Transport divergence gradient
            ztdgu = znugdt * (gcx(ji+1,jj  ) - gcx(ji,jj) ) / e1u(ji,jj)
            ztdgv = znugdt * (gcx(ji  ,jj+1) - gcx(ji,jj) ) / e2v(ji,jj)
            ! multiplied by z2dt
#if defined key_obc
            ! caution : grad D = 0 along open boundaries
            spgu(ji,jj) = z2dt * ztdgu * obcumask(ji,jj)
            spgv(ji,jj) = z2dt * ztdgv * obcvmask(ji,jj)
#else
            spgu(ji,jj) = z2dt * ztdgu
            spgv(ji,jj) = z2dt * ztdgv
#endif
         END DO
      END DO

#if defined key_agrif      
      IF( .NOT. Agrif_Root() ) THEN
         ! caution : grad D (fine) = grad D (coarse) at coarse/fine interface
         IF( (nbondi == -1) .OR. (nbondi == 2) ) spgu(2,:) = znugdt * z2dt * laplacu(2,:) * umask(2,:,1)
         IF( (nbondi ==  1) .OR. (nbondi == 2) ) spgu(nlci-2,:) = znugdt * z2dt * laplacu(nlci-2,:) * umask(nlci-2,:,1)
         IF( (nbondj == -1) .OR. (nbondj == 2) ) spgv(:,2) = znugdt * z2dt * laplacv(:,2) * vmask(:,2,1)
         IF( (nbondj ==  1) .OR. (nbondj == 2) ) spgv(:,nlcj-2) = znugdt * z2dt * laplacv(:,nlcj-2) * vmask(:,nlcj-2,1)
      ENDIF
#endif      
      ! 7.  Add the trends multiplied by z2dt to the after velocity
      ! -----------------------------------------------------------
      !     ( c a u t i o n : (ua,va) here are the after velocity not the
      !                       trend, the leap-frog time stepping will not
      !                       be done in dynnxt.F routine)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = (ua(ji,jj,jk) + spgu(ji,jj)) * umask(ji,jj,jk)
               va(ji,jj,jk) = (va(ji,jj,jk) + spgv(ji,jj)) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO


      ! Sea surface elevation time stepping
      ! -----------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN      ! Euler (forward) time stepping, no time filter
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! after free surface elevation
               zssha = sshb(ji,jj) + rdt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
               ! swap of arrays
               sshb(ji,jj) = sshn(ji,jj)
               sshn(ji,jj) = zssha
            END DO
         END DO
      ELSE                                           ! Leap-frog time stepping and time filter
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! after free surface elevation
               zssha = sshb(ji,jj) + z2dt * ( wn(ji,jj,1) - emp(ji,jj) * zraur ) * tmask(ji,jj,1)
               ! time filter and array swap
               sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
               sshn(ji,jj) = zssha
            END DO
         END DO
      ENDIF

      !                       ! print sum trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=sshn, clinfo1=' spg  - ssh: ', mask1=tmask )

   END SUBROUTINE dyn_spg_flt

#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart free surface cst volume
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_spg_flt( kt, kindic )       ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_spg_flt: You should not have seen this print! error?', kt, kindic
   END SUBROUTINE dyn_spg_flt
#endif
   
   !!======================================================================
END MODULE dynspg_flt
