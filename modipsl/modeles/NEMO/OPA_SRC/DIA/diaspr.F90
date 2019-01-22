MODULE diaspr
   !!======================================================================
   !!                       ***  MODULE  diaspr  ***
   !! Ocean diagnostics:  surface pressure (rigid-lid case) 
   !!=====================================================================
#if   defined key_diaspr   &&   defined key_dynspg_rl
   !!----------------------------------------------------------------------
   !!   'key_diaspr'        and                surface pressure diagnostics
   !!   'key_dynspg_rl'                                      rigid-lid case
   !!----------------------------------------------------------------------
   !!   dia_spr      : update momentum and tracer Kz from a tke scheme
   !!   sprmat       : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE sol_oce         ! ocean elliptic solver
   USE solpcg          ! preconditionned conjugate gradient solver
   USE solsor          ! Successive Over-relaxation solver
   USE solfet          ! FETI solver
   USE lib_mpp         ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dia_spr   ! routine called by step.F90

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_diaspr = .TRUE.    !: surface pressure diag. flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   gps         !: surface pressure

   !! * Module variables
   INTEGER ::                 &
      nmoyps,                 &  ! time step for average
      nindic,                 &  ! indicator of convergence of the solver
      !                          ! namspr  surface pressure diagnostic
      nmaxp ,                 &  ! maximum of iterations for the solver
      niterp                     ! number of iteration done by the solver

   REAL(wp) ::     &
      ! namspr  surface pressure diagnostic
      epsp                       ! absolute precision of the solver

      !! * Namelist
      NAMELIST/namspr/ nmaxp, epsp, niterp

   REAL(wp) ::     &
      e1e2t                      ! ???

   REAL(wp), PUBLIC DIMENSION(jpi,jpj) ::   &
      spgum, spgvm,           &  ! average value of the surface pressure gradients
      gpsuu, gpsvv,           &  ! surface pressure gradients computed from comp. PS
      gcdpsc,                 &  ! inverse diagonal preconditioning matrix
      gcsmat,                 &  ! diagonal preconditioning matrix
      spmsk                      ! surface pressure Mask

   REAL(wp), DIMENSION(jpi,jpj,4) ::   &
      gcps                       ! extra-diagonal elements of SPG matrix
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diaspr.F90,v 1.4 2005/12/21 10:46:34 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_spr( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_spr  ***
      !!
      !! ** Purpose :   compute the surface pressure from its gradient
      !!
      !! ** Method  :   rigid-lid appromimation: the surface pressure 
      !!      gradient is given by:
      !!           spgu = 1/rau0 d/dx(ps) = Mu + 1/(hu e2u) dj-1(bsfd)
      !!           spgv = 1/rau0 d/dy(ps) = Mv - 1/(hv e1v) di-1(bsfd)
      !!
      !!      where (Mu,Mv) is the vertically averaged momentum trend, i.e.
      !!      the vertical ponderated sum of the general momentum trend.
      !!      where bsfd is the trend of the barotropic stream function.
      !!
      !!       taking the divergence of the surface pressure gradient provides
      !!      an elliptic equation for ps which is solved using either a
      !!      diagonal preconditioned conjugate gradient method (solpcg.f) or
      !!      an successive-over-relaxation method (solsor.f) or FETI method
      !!      (solfet.F).
      !!
      !!      n.b. this resolution is valid with topography, cyclic east-west
      !!      boundary conditions and islands.
      !!
      !! History :
      !!        !  98-01  (G. Madec & M. Ioualalen)  Original code
      !!        !  98-02  (M. Guyon)  FETI method
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER  ::  ji, jj
      INTEGER  ::   imax, ijt, iju
      REAL(wp) ::   zpsmea, zeps, zmoyr
      REAL(wp) ::   ztab(jpi,jpj,8)
      REAL(wp) ::   zemin1, zemax1, zemin2, zemax2, zgwgt
      REAL(wp) ::   z1, z2, zcompt,z3,z4
      REAL(wp) ::   zdif1, zdif2, zvar1, zvar2
      !!----------------------------------------------------------------------


      ! 0. initialisation (the first time step)
      ! ---------------------------------------
      
      IF( kt == nit000 ) THEN

         ! Namelist namspr : surface pressure

         nmaxp  = 2000
         epsp   = 1.e-6
         niterp = 16

         ! Read Namelist namspr : surface pressure diagnostics
         REWIND ( numnam )
         READ(numnam,namspr)

         IF(lwp) THEN
            WRITE(numout,*) 'dia_spr : surface pressure diagnostic (rigid-lid case)'
            WRITE(numout,*) '~~~~~~~'
            WRITE(numout,*)
            WRITE(numout,*) '          Namelist namspr : set solver parameters'
            WRITE(numout,*)
            WRITE(numout,*) '             maximum iterations for solver  nmaxp  = ', nmaxp
            WRITE(numout,*) '             absolute precision of solver   epsp   = ', epsp
            WRITE(numout,*) '             number of solver iterations    niterp = ', niterp
            WRITE(numout,*) '             frequeny of averaged output    nwrite = ', nwrite
            WRITE(numout,*)
         ENDIF

         ! control
# if ! defined key_dynspg_rl
      IF(lwp) WRITE(numout,cform_err)
      IF(lwp) WRITE(numout,*) '          surface pressure already explicitly computed !!'
      nstop = nstop + 1
# endif

         ! compute the ocean surface
         e1e2t = 0.e0
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               e1e2t = e1e2t + e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
            END DO
         END DO
         IF( lk_mpp )   CALL  mpp_sum( e1e2t )   ! sum over the global domain
         
         ! build the matrix for the surface pressure
         CALL sprmat
         
         ! set to zero the mean surface pressure gradient
         nmoyps = 0
         spgum(:,:) = 0.e0
         spgvm(:,:) = 0.e0

      ENDIF

      ! 1. cumulate the surface pressure gradient (at each time step)
      ! -----------------------------------------

      nmoyps = nmoyps + 1
      spgum(:,:) = spgum(:,:) + spgu(:,:)
      spgvm(:,:) = spgvm(:,:) + spgv(:,:)
      

      ! 2. ps computation each nwrite time step
      ! ---------------------------------------
      
      ! RETURN IF not the right time to compute ps
      IF ( MOD(kt-nit000+1,nwrite) /= 0 ) RETURN
      
      
      ! mean surface pressure gradient
      !   averaging and mask
      zmoyr = 1./float(nmoyps)
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            spgum(ji,jj) = spgum(ji,jj) * zmoyr * umask(ji,jj,1)
            spgvm(ji,jj) = spgvm(ji,jj) * zmoyr * vmask(ji,jj,1)
         END DO
      END DO

      CALL  lbc_lnk(spgum, 'U', -1. )
      CALL  lbc_lnk(spgvm, 'V', -1. )

      
      ! SAVE in local arrays and variables of solver informations
      zeps   = eps
      imax  = nmax 
      ztab(:,:,1) = gcp   (:,:,1)
      ztab(:,:,2) = gcp   (:,:,2)
      ztab(:,:,3) = gcp   (:,:,3)
      ztab(:,:,4) = gcp   (:,:,4)
      ztab(:,:,5) = gcdprc(:,:  )
      ztab(:,:,6) = gcdmat(:,:  )
      ztab(:,:,7) = gcx   (:,:  )
      ztab(:,:,8) = bmask (:,:  )

      ! replace bsf solver informations by ps solver one
      eps    = epsp
      nmax   = nmaxp
      gcp   (:,:,1) = gcps  (:,:,1)
      gcp   (:,:,2) = gcps  (:,:,2)
      gcp   (:,:,3) = gcps  (:,:,3)
      gcp   (:,:,4) = gcps  (:,:,4)
      gcdprc(:,:  ) = gcdpsc(:,:  )
      gcdmat(:,:  ) = gcsmat(:,:  )
      bmask (:,:  ) = spmsk (:,:  )
      !    first guess: ps
      gcx   (:,:  ) = gps   (:,:  )

      !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,

      ! right hand side: 2d div. of the surface pressure gradient
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            gcb(ji,jj) = -gcdpsc(ji,jj)*   &
               (  e2u(ji,jj)*spgum(ji,jj) - e2u(ji-1,jj)*spgum(ji-1,jj)   &
               + e1v(ji,jj)*spgvm(ji,jj) - e1v(ji,jj-1)*spgvm(ji,jj-1) )
         END DO
      END DO
      
      !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
      
      ! relative PRECISION
      rnorme = 0.
      DO jj = 1, jpj
         DO ji = 1, jpi
            rnorme = rnorme + gcb(ji,jj) * gcsmat(ji,jj) * gcb(ji,jj)
         END DO
      END DO
      IF( lk_mpp )   CALL  mpp_sum( rnorme )   ! sum over the global domain

      epsr=eps*eps*rnorme
      ncut=0
      !   IF the second member is 0 the solution is 0, solpcg isn't called
      IF ( rnorme == 0.e0 ) THEN
         gps(:,:) = 0.e0
         res   = 0.e0
         niter = 0
         ncut  = 999
      ENDIF
      
      !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
      
      nindic = 0

      ! iterarive solver of the spg system (except IF sol.=0)
      !     (OUTPUT in gcx with boundary conditions applied)
      IF ( ncut == 0 ) THEN
         IF ( nsolv == 1 ) THEN
            CALL sol_pcg( nindic )         !   diagonal preconditioned conjuguate gradient
         ELSE IF ( nsolv == 2 ) THEN
            CALL sol_sor( nindic )     !   successive-over-relaxation
         ELSE IF(nsolv == 3) THEN
            CALL sol_fet( nindic )         !   FETI solver
         ELSE
            !   e r r o r  in nsolv namelist PARAMETER
            IF(lwp) THEN
               WRITE(numout,*) ' dia_spr: e r r o r, nsolv = 1 or 2'
               WRITE(numout,*) ' ******               not = ',nsolv
            ENDIF
            STOP 'dia_spr'
         ENDIF
      ENDIF
      
      !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,

      
      ! sp solver statistics  (i.e. problem for the solver) 
      IF ( epsr < 0.) THEN
         IF(lwp) THEN 
            WRITE(numout,*)'rrrrrrrrrrrrrrr'
            IF(lwp)WRITE(numout,*)'dia_spr-1:',epsr
            IF(lwp)WRITE(numout,*)'rrrrrrrrrrrrrrr'
         ENDIF
      ENDIF
      IF(lwp)WRITE(numout,9300) kt, niter, res, SQRT(epsr)/eps
      IF (nindic < 0) THEN 
         IF(lwp) THEN 
            WRITE(numout,9100)
            WRITE(numout,*) ' dia_spr : the surface pressure solver DO not converge'
            WRITE(numout,*) ' ====== ' 
            WRITE(numout,*) 
         ENDIF
      ENDIF
9100  FORMAT( /,' ===>>>> : w a r n i n g',/,'          ===============',/ )
9300  FORMAT(' it :', i8, ' niter :', i4, ' res :',e20.10,' b :', e20.10)

      ! recover bsf solver informations and SAVE ps for next computation
      eps    = zeps
      nmax   = imax 
      gps   (:,:  ) = gcx (:,:)
      gcp   (:,:,1) = ztab(:,:,1)
      gcp   (:,:,2) = ztab(:,:,2)
      gcp   (:,:,3) = ztab(:,:,3)
      gcp   (:,:,4) = ztab(:,:,4)
      gcdprc(:,:  ) = ztab(:,:,5)
      gcdmat(:,:  ) = ztab(:,:,6)
      gcx   (:,:  ) = ztab(:,:,7)
      bmask (:,:  ) = ztab(:,:,8)
      
      ! compute and substract the mean value
      
      zpsmea = 0.e0
      DO jj=2,jpjm1
         DO ji=2,jpim1
            zpsmea = zpsmea + gps(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask(ji,jj,1)
         END DO
      END DO
      IF( lk_mpp )   CALL  mpp_sum( zpsmea )   ! sum over the global domain

      zpsmea = zpsmea / e1e2t
      gps(:,:) = ( gps(:,:) - zpsmea ) * tmask(:,:,1)
 
      IF(lwp)WRITE(numout,*) ' mean value of ps = ',zpsmea,' is substracted'
      ! ----------------------------------------
      ! i. compute the surface pressure gradient
      !    from the computed surface pressure
      ! ----------------------------------------

      DO jj=2,jpjm1
         DO ji=2,jpim1
            gpsuu(ji,jj)=(gps(ji+1,jj)-gps(ji,jj))/e1u(ji,jj) * umask(ji,jj,1)
            gpsvv(ji,jj)=(gps(ji,jj+1)-gps(ji,jj))/e2v(ji,jj) * vmask(ji,jj,1)
         END DO
      END DO
      
      ! compute the max and min error
      
      zemax1 = 0.e0
      zemin1 = 0.e0
      zemax2 = 0.e0
      zemin2 = 0.e0
      DO jj = 2,jpj-1
         DO ji = 2,jpi-1
            z1 = ABS( spgum(ji,jj)-gpsuu(ji,jj) )*umask(ji,jj,1)
            z2 = ABS( spgvm(ji,jj)-gpsvv(ji,jj) )*vmask(ji,jj,1)
            z3 = MAX ( ABS( spgum(ji,jj) ), ABS( spgvm(ji,jj) ) )
            z4 = MAX ( ABS( gpsuu(ji,jj) ), ABS( gpsvv(ji,jj) ) )
            zemax1 = MAX(z1,zemax1)
            zemax2 = MAX(z2,zemax2)
            zemin1 = MAX(z3,zemin1)
            zemin2 = MAX(z4,zemin2)
         END DO
      END DO
      IF( lk_mpp )   CALL  mpp_sum( zemax1 )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zemax2 )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zemin1 )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zemin2 )   ! sum over the global domain

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'pserro : time step = ', kt
         WRITE(numout,*) '******** ------------------'
         WRITE(numout,*)
         WRITE(numout,*) '         gpsx error  max=',zemax1
         WRITE(numout,*) '         gpsy error  max=',zemax2
         WRITE(numout,*) '         gps max =',zemin1
         WRITE(numout,*) '         gpsc max =',zemin2
         WRITE(numout,*)
      ENDIF

      ! compute the norme and variance of this error

      zcompt = 0.e0
      zdif1  = 0.e0
      zdif2  = 0.e0
      zvar1  = 0.e0
      zvar2  = 0.e0
      DO jj = 2, jpj-1
         DO ji = 2, jpi-1
            z1 = ( spgum(ji,jj)-gpsuu(ji,jj) ) * umask(ji,jj,1)
            z2 = ( spgvm(ji,jj)-gpsvv(ji,jj) ) * vmask(ji,jj,1)
            zcompt=zcompt+tmask(ji,jj,1)
            zdif1=zdif1+z1
            zdif2=zdif2+z2
            zvar1=zvar1+z1*z1
            zvar2=zvar2+z2*z2
         END DO
      END DO
      IF( lk_mpp )   CALL  mpp_sum( zcompt )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zdif1  )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zdif2  )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zvar1  )   ! sum over the global domain
      IF( lk_mpp )   CALL  mpp_sum( zvar2  )   ! sum over the global domain

      IF(lwp) WRITE(numout,*) '        zcompt = ',zcompt
      zdif1=zdif1/zcompt
      zdif2=zdif2/zcompt
      IF( zvar1 < 0.) THEN 
         IF(lwp) THEN
            WRITE(numout,*)'rrrrrrrrrrrrrrr'
            WRITE(numout,*)'dia_spr-2:',zvar1
            WRITE(numout,*)'rrrrrrrrrrrrrrr'
         ENDIF
      ENDIF
      zvar1 = SQRT(zvar1)/zcompt
      IF( zvar2 < 0. ) THEN 
         IF(lwp)THEN
            WRITE(numout,*)'rrrrrrrrrrrrrrr'
            WRITE(numout,*)'dia_spr-3:',zvar2
            WRITE(numout,*)'rrrrrrrrrrrrrrr'
         ENDIF
      ENDIF
      zvar2 = SQRT(zvar2)/zcompt
      
      IF(lwp) THEN 
         WRITE(numout,*)
         WRITE(numout,*) '         gpsx mean error = ',zdif1
         WRITE(numout,*) '         gpsy mean error = ',zdif2
         WRITE(numout,*)
         WRITE(numout,*) '         gpsx var. error = ',zvar1
         WRITE(numout,*) '         gpsy var. error = ',zvar2
         WRITE(numout,*)
         WRITE(numout,*)
      ENDIF
      
      ! reset to zero nmoyps and the mean surface pressure gradient
      nmoyps = 0
      spgum(:,:) = 0.e0
      spgvm(:,:) = 0.e0
      
   END SUBROUTINE dia_spr


   SUBROUTINE sprmat
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sprmat  ***
      !!               
      !! ** Purpose :   construction of the matrix of the surface pressure
      !!      system and the diagonal preconditioning matrix.
      !!
      !! ** Method :
      !!
      !! History :
      !!        !  98-01  (G. Madec & M. Ioualalen)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jl                ! dummy loop indices
      REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn
      !!----------------------------------------------------------------------
      
      
      ! 0. ocean/land mask at ps-point (computed from surface tmask): spmsk
      ! --------------------------------------------------------------------
      
      ! computation
      spmsk(:,:) = tmask(:,:,1)
      
      ! boundary conditions
      ! south symmetry: psmsk must be set to 0. on 1
      IF( nperio == 2 ) THEN
         spmsk(:, 1 ) = 0.e0
      ENDIF
      
      ! east-west cyclic: spmsk must be set to 0. on 1 and jpi
      IF( nperio == 1 .OR. nperio == 4 .OR.nperio == 6) THEN 
         spmsk( 1 ,:) = 0.e0
         spmsk(jpi,:) = 0.e0
      ENDIF
      
      ! north fold: spmsk must be set to 0. on ligne jpj and on half
      !                   ligne jpj-1
      ! T-point pivot
      IF( nperio == 3 .OR. nperio == 4 ) THEN
         spmsk(:,jpj) = 0.e0
         DO ji = jpi/2+1, jpi
            spmsk(ji,jpjm1) = 0.e0
         END DO
      ENDIF
      ! F-point pivot
      IF( nperio == 5 .OR. nperio == 6 ) THEN
         spmsk(:,jpj) = 0.e0
      ENDIF
      
      ! mpp boundary cond.: spmsk is initialized at zero on the overlap
      ! region for both the preconjugate gradient and the sor algorithms

      IF( nbondi /= -1 .AND. nbondi /= 2 ) THEN
         DO jl = 1, jpreci
            spmsk(jl,:) = 0.e0
         END DO
      ENDIF
      IF( nbondi /= 1 .AND. nbondi /= 2 ) THEN
         DO ji = nlci, jpi
            spmsk(ji,:) = 0.e0
         END DO
      ENDIF
      IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN
         DO jl=1,jprecj
            spmsk(:,jl) = 0.e0
         END DO
      ENDIF
      IF( nbondj /= 1 .AND. nbondj /= 2 ) THEN
         DO jj = nlcj, jpj
            spmsk(:,jj) = 0.e0
         END DO
      ENDIF
      
      ! 1. construction of the matrix
      ! -----------------------------
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            
            IF( spmsk(ji,jj) == 0. ) THEN
               ! land points
               gcps  (ji,jj,1) = 0.e0
               gcps  (ji,jj,2) = 0.e0
               gcps  (ji,jj,3) = 0.e0
               gcps  (ji,jj,4) = 0.e0
               gcdpsc(ji,jj  ) = 0.e0
               gcsmat(ji,jj  ) = 0.e0
            ELSE
               ! south coefficient
               zcoefs = -e1v(ji,jj-1) / e2v(ji,jj-1) * vmask(ji,jj-1,1)
               gcps(ji,jj,1) = zcoefs
               ! west coefficient
               zcoefw = -e2u(ji-1,jj) / e1u(ji-1,jj) * umask(ji-1,jj,1)
               gcps(ji,jj,2) = zcoefw
               ! east coefficient
               zcoefe = -e2u(ji  ,jj) / e1u(ji  ,jj) * umask(ji  ,jj,1)
               gcps(ji,jj,3) = zcoefe
               ! north coefficient
               zcoefn = -e1v(ji,jj  ) / e2v(ji,jj  ) * vmask(ji,jj  ,1)
               gcps(ji,jj,4) = zcoefn
               
               ! diagonal coefficient
               gcsmat(ji,jj) = -zcoefs-zcoefw-zcoefe-zcoefn
            ENDIF
         END DO
      END DO
      

      ! 2. boundary conditions 
      ! ----------------------
      
      ! cyclic east-west boundary conditions
      ! ji=2 is the column east of ji=jpim1 and reciprocally,
      ! ji=jpim1 is the column west of ji=2
      ! all the coef are already set to zero as spmask is initialized to
      ! zero for ji=1 and ji=jpj.
      
      ! symetrical conditions
      ! the diagonal coefficient of the southern grid points must be modify to
      ! account for the existence of the south symmetric bassin.
      IF( nperio == 2 ) THEN
         DO ji = 1, jpi
            IF( spmsk(ji,2) /= 0 ) THEN
               zcoefs = e1v(ji,1) / e2v(ji,1)
               gcsmat(ji,2) = gcsmat(ji,2) - zcoefs
            ENDIF
         END DO
      ENDIF
      
      ! North fold boundary condition
      ! all the coef are already set to zero as bmask is initialized to
      ! zero on duplicated lignes and portion of lignes
      
      
      ! 3. preconditioned matrix
      ! ------------------------
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( spmsk(ji,jj) /= 0. ) gcdpsc(ji,jj) = 1.e0 / gcsmat(ji,jj)
         END DO
      END DO
      
      gcps(:,:,1) = gcps(:,:,1) * gcdpsc(:,:)
      gcps(:,:,2) = gcps(:,:,2) * gcdpsc(:,:)
      gcps(:,:,3) = gcps(:,:,3) * gcdpsc(:,:)
      gcps(:,:,4) = gcps(:,:,4) * gcdpsc(:,:)
      
      
      ! 3. initialization the arrays used in sp solver
      ! ----------------------------------------------
      
      gps  (:,:) = 0.e0
      gpsuu(:,:) = 0.e0
      gpsvv(:,:) = 0.e0
      
   END SUBROUTINE sprmat

#else
   !!----------------------------------------------------------------------
   !!   Default option :                    NO surface pressure diagnostics
   !!----------------------------------------------------------------------
   USE in_out_manager  
   LOGICAL, PUBLIC, PARAMETER ::   lk_diaspr = .FALSE.   !: surface pressure diag. flag
CONTAINS
   SUBROUTINE dia_spr( kt )      ! Empty routine
      if(lwp) WRITE(numout,*) 'dia_spr: You should not have seen this print! error?', kt
   END SUBROUTINE dia_spr
#endif

   !!======================================================================
END MODULE diaspr
