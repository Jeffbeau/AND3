MODULE solmat
   !!======================================================================
   !!                       ***  MODULE  solmat  ***
   !! solver       : construction of the matrix 
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   sol_mat       : Construction of the matrix of used by the elliptic solvers
   !!   fetsch        :
   !!   fetmat        :
   !!   fetstr        :
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE sol_oce         ! ocean solver
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lbclnk          ! lateral boudary conditions
   USE lib_mpp         ! distributed memory computing
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sol_mat     ! routine called by inisol.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solmat.F90,v 1.11 2006/03/20 17:27:14 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sol_mat( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_mat  ***
      !!
      !! ** Purpose :   Construction of the matrix of used by the elliptic 
      !!      solvers (either sor, pcg or feti methods).
      !!
      !! ** Method  :   The matrix depends on the type of free surface:
      !!       * lk_dynspg_rl=T: rigid lid formulation
      !!      The matrix is built for the barotropic stream function system.
      !!      a diagonal preconditioning matrix is also defined.
      !!       * lk_dynspg_flt=T: free surface formulation
      !!      The matrix is built for the divergence of the transport system
      !!      a diagonal preconditioning matrix is also defined.
      !!        Note that for feti solver (nsolv=3) a specific initialization 
      !!      is required (call to fetstr.F) for memory allocation and inter-
      !!      face definition.
      !! 
      !! ** Action  : - gcp    : extra-diagonal elements of the matrix
      !!              - gcdmat : preconditioning matrix (diagonal elements)
      !!              - gcdprc : inverse of the preconditioning matrix
      !!
      !! History :
      !!   1.0  !  88-04  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  93-06  (M. Guyon)  suppress pointers
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!        !  98-02  (M. Guyon)  FETI method
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!        !  02-11  (C. Talandier, A-M. Treguier) Free surface & Open boundaries
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in) :: kt

      !! * Local declarations
      INTEGER ::   ji, jj                    ! dummy loop indices
      INTEGER ::   ii, ij, iiend, ijend      ! temporary integers
      REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn  ! temporary scalars
      REAL(wp) ::   z2dt, zcoef
      !!----------------------------------------------------------------------

      ! FETI method ( nsolv = 3)
      ! memory allocation and interface definition for the solver

      IF( nsolv == 3 )   CALL fetstr

      
      ! 1. Construction of the matrix
      ! -----------------------------
      
      ! initialize to zero
      zcoef = 0.e0
      gcp(:,:,1) = 0.e0
      gcp(:,:,2) = 0.e0
      gcp(:,:,3) = 0.e0
      gcp(:,:,4) = 0.e0
      
      gcdprc(:,:) = 0.e0
      gcdmat(:,:) = 0.e0
      
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         z2dt = rdt
      ELSE
         z2dt = 2. * rdt
      ENDIF

#if defined key_dynspg_flt && ! defined key_obc
!!cr      IF( lk_dynspg_flt .AND. .NOT.lk_obc ) THEN   !bug missing lk_dynspg_flt_atsk

      ! defined the coefficients for free surface elliptic system

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
            zcoefs = -zcoef * hv(ji  ,jj-1) * e1v(ji  ,jj-1) / e2v(ji  ,jj-1)    ! south coefficient
            zcoefw = -zcoef * hu(ji-1,jj  ) * e2u(ji-1,jj  ) / e1u(ji-1,jj  )    ! west coefficient
            zcoefe = -zcoef * hu(ji  ,jj  ) * e2u(ji  ,jj  ) / e1u(ji  ,jj  )    ! east coefficient
            zcoefn = -zcoef * hv(ji  ,jj  ) * e1v(ji  ,jj  ) / e2v(ji  ,jj  )    ! north coefficient
            gcp(ji,jj,1) = zcoefs
            gcp(ji,jj,2) = zcoefw
            gcp(ji,jj,3) = zcoefe
            gcp(ji,jj,4) = zcoefn
            gcdmat(ji,jj) = e1t(ji,jj) * e2t(ji,jj) * bmask(ji,jj)    &          ! diagonal coefficient
               &          - zcoefs -zcoefw -zcoefe -zcoefn
         END DO
      END DO
      
#  elif defined key_dynspg_flt && defined key_obc
!!cr      ELSEIF( lk_dynspg_flt .AND. lk_obc ) THEN     !bug missing lk_dynspg_flt_atsk 

      !   defined gcdmat in the case of open boundaries

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
            !  south coefficient
            IF( lp_obc_south .AND. ( jj == njs0p1 ) ) THEN
               zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)*(1.-vsmsk(ji,1))
            ELSE
               zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)
            END IF
            gcp(ji,jj,1) = zcoefs

            !  west coefficient
            IF( lp_obc_west  .AND. ( ji == niw0p1 ) ) THEN
               zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)*(1.-uwmsk(jj,1))
            ELSE
               zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)
            END IF
            gcp(ji,jj,2) = zcoefw

            !   east coefficient
            IF( lp_obc_east  .AND. ( ji == nie0 ) ) THEN
               zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)*(1.-uemsk(jj,1))
            ELSE
               zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)
            END IF
            gcp(ji,jj,3) = zcoefe

            !   north coefficient
            IF( lp_obc_north .AND. ( jj == njn0 ) ) THEN
               zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)*(1.-vnmsk(ji,1))
            ELSE
               zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)
            END IF
            gcp(ji,jj,4) = zcoefn

            ! diagonal coefficient
            gcdmat(ji,jj) = e1t(ji,jj)*e2t(ji,jj)*bmask(ji,jj) &
                            - zcoefs -zcoefw -zcoefe -zcoefn
         END DO
      END DO

#  else
!!cr      ELSE

      !   defined the coefficients for bsf elliptic system
      
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zcoefs = -hur(ji  ,jj  ) * e1u(ji  ,jj  ) / e2u(ji  ,jj  ) * bmask(ji,jj)   ! south coefficient
            zcoefw = -hvr(ji  ,jj  ) * e2v(ji  ,jj  ) / e1v(ji  ,jj  ) * bmask(ji,jj)   ! west coefficient
            zcoefe = -hvr(ji+1,jj  ) * e2v(ji+1,jj  ) / e1v(ji+1,jj  ) * bmask(ji,jj)   ! east coefficient
            zcoefn = -hur(ji  ,jj+1) * e1u(ji  ,jj+1) / e2u(ji  ,jj+1) * bmask(ji,jj)   ! north coefficient
            gcp(ji,jj,1) = zcoefs
            gcp(ji,jj,2) = zcoefw
            gcp(ji,jj,3) = zcoefe
            gcp(ji,jj,4) = zcoefn
            gcdmat(ji,jj) = -zcoefs -zcoefw -zcoefe -zcoefn                             ! diagonal coefficient
         END DO
      END DO
      
!!cr  ENDIF
#endif
#if defined key_agrif
       IF (.NOT.AGRIF_ROOT()) THEN
       
       IF ( (nbondi == -1)  .OR. (nbondi == 2) ) bmask(2,:)=0.
       IF ( (nbondi ==  1)  .OR. (nbondi == 2) ) bmask(nlci-1,:)=0.
       IF ( (nbondj == -1)  .OR. (nbondj == 2) ) bmask(:,2)=0.
       IF ( (nbondj ==  1)  .OR. (nbondj == 2) ) bmask(:,nlcj-1)=0.

      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zcoef = z2dt * z2dt * grav * rnu * bmask(ji,jj)
            !  south coefficient
            IF( ((nbondj == -1)  .OR. (nbondj == 2)) .AND. ( jj == 3 ) ) THEN
               zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)*(1.-vmask(ji,jj-1,1))
            ELSE
               zcoefs = -zcoef * hv(ji,jj-1) * e1v(ji,jj-1)/e2v(ji,jj-1)
            END IF
            gcp(ji,jj,1) = zcoefs

            !  west coefficient
	    IF( ( (nbondi == -1)  .OR. (nbondi == 2) ) .AND. ( ji == 3 )  ) THEN
               zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)*(1.-umask(ji-1,jj,1))
            ELSE
               zcoefw = -zcoef * hu(ji-1,jj) * e2u(ji-1,jj)/e1u(ji-1,jj)
            END IF
            gcp(ji,jj,2) = zcoefw

            !   east coefficient
            IF( ((nbondi == 1)  .OR. (nbondi == 2)) .AND. ( ji == nlci-2 ) ) THEN
               zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)*(1.-umask(ji,jj,1))
            ELSE
               zcoefe = -zcoef * hu(ji,jj) * e2u(ji,jj)/e1u(ji,jj)
            END IF
            gcp(ji,jj,3) = zcoefe

            !   north coefficient
            IF( ((nbondj == 1)  .OR. (nbondj == 2)) .AND. ( jj == nlcj-2 ) ) THEN
               zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)*(1.-vmask(ji,jj,1))
            ELSE
               zcoefn = -zcoef * hv(ji,jj) * e1v(ji,jj)/e2v(ji,jj)
            END IF
            gcp(ji,jj,4) = zcoefn

            ! diagonal coefficient
            gcdmat(ji,jj) = e1t(ji,jj)*e2t(ji,jj)*bmask(ji,jj) &
                            - zcoefs -zcoefw -zcoefe -zcoefn
         END DO
      END DO
      
       ENDIF
#endif

      ! 2. Boundary conditions 
      ! ----------------------
      
      ! Cyclic east-west boundary conditions
      !     ji=2 is the column east of ji=jpim1 and reciprocally,
      !     ji=jpim1 is the column west of ji=2
      !     all the coef are already set to zero as bmask is initialized to
      !     zero for ji=1 and ji=jpj in dommsk.
      
      ! Symetrical conditions
      ! free surface: no specific action
      ! bsf system: n-s gradient of bsf = 0 along j=2 (perhaps a bug !!!!!!)
      ! the diagonal coefficient of the southern grid points must be modify to
      ! account for the existence of the south symmetric bassin.
      
!!cr      IF( .NOT.lk_dynspg_flt ) THEN   !bug missing lk_dynspg_flt_atsk
#if ! defined key_dynspg_flt
      IF( nperio == 2 ) THEN
         DO ji = 1, jpi
            IF( bmask(ji,2) /= 0.e0 ) THEN
               zcoefs = - hur(ji,2)*e1u(ji,2)/e2u(ji,2)
               gcdmat(ji,2) = gcdmat(ji,2) - zcoefs
            ENDIF
         END DO
      ENDIF
!!cr      ENDIF
#endif
      
      ! North fold boundary condition
      !     all the coef are already set to zero as bmask is initialized to
      !     zero on duplicated lignes and portion of lignes
      
      ! 3. Preconditioned matrix
      ! ------------------------
      
      IF( nsolv /= 3 ) THEN
         
         ! SOR and PCG solvers
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( bmask(ji,jj) /= 0.e0 )   gcdprc(ji,jj) = 1.e0 / gcdmat(ji,jj)
            END DO
         END DO
         
         gcp(:,:,1) = gcp(:,:,1) * gcdprc(:,:)
         gcp(:,:,2) = gcp(:,:,2) * gcdprc(:,:)
         gcp(:,:,3) = gcp(:,:,3) * gcdprc(:,:)
         gcp(:,:,4) = gcp(:,:,4) * gcdprc(:,:)
         IF( ( nsolv == 2 ) .OR. ( nsolv == 4 ) )  gccd(:,:) = sor * gcp(:,:,2)

         IF( nsolv == 4 ) THEN
            CALL lbc_lnk_e( gcp   (:,:,1), c_solver_pt, 1. )   ! lateral boundary conditions
            CALL lbc_lnk_e( gcp   (:,:,2), c_solver_pt, 1. )   ! lateral boundary conditions
            CALL lbc_lnk_e( gcp   (:,:,3), c_solver_pt, 1. )   ! lateral boundary conditions
            CALL lbc_lnk_e( gcp   (:,:,4), c_solver_pt, 1. )   ! lateral boundary conditions
            CALL lbc_lnk_e( gcdprc(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions
            CALL lbc_lnk_e( gcdmat(:,:)  , c_solver_pt, 1. )   ! lateral boundary conditions         
            IF( npolj /= 0 ) CALL sol_exd( gcp , c_solver_pt ) ! switch northernelements
         END IF

      ELSE
         
         ! FETI method
         ! if feti solver : gcdprc is a mask for the non-overlapping
         !   data structuring
         
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( bmask(ji,jj) /= 0.e0 ) THEN
                  gcdprc(ji,jj) = 1.e0
               ELSE
                  gcdprc(ji,jj) = 0.e0
               ENDIF
            END DO
         END DO
         
         ! so "common" line & "common" column have to be !=0 except on global
         !   domain boundaries
         ! pbs with nbondi if nperio != 2 ?
         !   ii = nldi-1
         ! pb with nldi value if jperio==1 : nbondi modifyed at the end
         !   of inimpp.F => pb
         ! pb with periodicity conditions : iiend, ijend
         
         ijend = nlej
         iiend = nlei
         IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )   iiend = nlci - jpreci
         ii = jpreci
         
         ! case number 1
         
         IF( nbondi /= -1 .AND. nbondi /= 2 ) THEN
            DO jj = 1, ijend
               IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
            END DO
         ENDIF
         
         ! case number 2
         
         IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
            DO jj = 1, ijend
               IF( fmask(ii,jj,1) == 1. ) gcdprc(ii,jj) = 1.
            END DO
         ENDIF
         
         !      ij = nldj-1
         ! pb with nldi value if jperio==1 : nbondi modifyed at the end
         !   of inimpp.F => pb, here homogeneisation...
         
         ij = jprecj
         IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN
            DO ji = 1, iiend
               IF( fmask(ji,ij,1) == 1. ) gcdprc(ji,ij) = 1.
            END DO
         ENDIF
      ENDIF
      
      
      ! 4. Initialization the arrays used in pcg
      ! ----------------------------------------
      gcx  (:,:) = 0.e0
      gcxb (:,:) = 0.e0
      gcb  (:,:) = 0.e0
      gcr  (:,:) = 0.e0
      gcdes(:,:) = 0.e0
      gccd (:,:) = 0.e0
      
      ! FETI method
      IF( nsolv == 3 ) THEN
         CALL fetmat       ! Matrix treatment : Neumann condition, inverse computation
         CALL fetsch       ! data framework for the Schur Dual solver
      ENDIF
      
   END SUBROUTINE sol_mat


   SUBROUTINE sol_exd( pt3d, cd_type )
      !!----------------------------------------------------------------------
      !!                  ***  routine sol_exd  ***
      !!                  
      !! ** Purpose :   Reorder gcb coefficient on the extra outer  halo 
      !!                at north fold in case of T or F pivot
      !!
      !! ** Method  :   Perform a circular permutation of the coefficients on 
      !!                the total area strictly above the pivot point,
      !!                and on the semi-row of the pivot point   
      !!                
      !! History :
      !!   9.0  !  05-09  (R. Benshila)  original routine
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=1) , INTENT( in ) ::   &
         cd_type       ! define the nature of pt2d array grid-points
         !             !  = T , U , V , F , W 
         !             !  = S : T-point, north fold treatment
         !             !  = G : F-point, north fold treatment
         !             !  = I : sea-ice velocity at F-point with index shift
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4), INTENT( inout ) ::   &
         pt3d          ! 2D array on which the boundary condition is applied

      !! * Local variables
      INTEGER  ::   ji, jk      ! dummy loop indices
      INTEGER  ::   iloc                ! temporary integers
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &
         ztab          ! 2D array on which the boundary condition is applied
      !!----------------------------------------------------------------------

      ztab = pt3d

      ! north fold treatment
      ! -----------------------
  
      SELECT CASE ( npolj )
         
         CASE ( 3 , 4 )   !  T pivot
         iloc = jpiglo/2 +1 
            
            SELECT CASE ( cd_type )
  
            CASE ( 'T', 'S', 'U', 'W' )
               DO jk =1, 4
                  DO ji = 1-jpr2di, nlci+jpr2di
                     pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
                  ENDDO
               ENDDO

              DO jk =1, 4
                  DO ji = nlci+jpr2di, 1-jpr2di,  -1
                     IF( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
                       & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
                     pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
                  ENDDO
               ENDDO

            CASE ( 'F' ,'G' , 'I', 'V' )
               DO jk =1, 4
                  DO ji = 1-jpr2di, nlci+jpr2di
                     pt3d(ji,nlcj-1:nlcj+jpr2dj,jk) = ztab(ji,nlcj-1:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
                  ENDDO
               ENDDO

            END SELECT   ! cd_type
  
         CASE ( 5 , 6 )                 ! F pivot
          iloc=jpiglo/2

            SELECT CASE (cd_type )

            CASE ( 'T'  ,'S', 'U', 'W')
               DO jk =1, 4
                  DO ji = 1-jpr2di, nlci+jpr2di
                     pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
                  ENDDO
               ENDDO

            CASE ( 'F' ,'G' , 'I', 'V' )
               DO jk =1, 4
                  DO ji = 1-jpr2di, nlci+jpr2di
                     pt3d(ji,nlcj:nlcj+jpr2dj,jk) = ztab(ji,nlcj:nlcj+jpr2dj,jk+3-2*MOD(jk+3,4))           
                  ENDDO
               ENDDO
               DO jk =1, 4
                  DO ji = nlci+jpr2di, 1-jpr2di,  -1
                    IF ( ( ji .LT. mi0(iloc) .AND. mi0(iloc) /= 1 ) &
                       & .OR. ( mi0(iloc) == jpi+1 ) ) EXIT
                    pt3d(ji,nlcj-1,jk) = ztab(ji,nlcj-1,jk+3-2*MOD(jk+3,4))
                  ENDDO
               ENDDO

            END SELECT   ! cd_type

         END SELECT   ! npolj
  
   END SUBROUTINE sol_exd

#if defined key_feti

   SUBROUTINE fetstr
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE fetstr  ***
      !!               
      !! ** Purpose :   Construction of the matrix of the barotropic stream 
      !!      function system.
      !!      Finite Elements Tearing & Interconnecting (FETI) approach
      !!      Memory allocation and interface definition for the solver
      !!
      !! ** Method :
      !!
      !! References :
      !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
      !!      A domain decomposition solver to compute the barotropic 
      !!      component of an OGCM in the parallel processing field.
      !!      Ocean Modelling, issue 105, december 94.
      !!
      !! History :
      !!        !  98-02 (M. Guyon)  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE lib_feti            ! feti librairy
      !! * Local declarations
      INTEGER ::   iiend, ijend, iperio   ! temporary integers
      !!---------------------------------------------------------------------
      
      
      ! Preconditioning technics of the Dual Schur Operator
      ! <= definition of the Coarse Grid solver
      ! <= dimension of the nullspace of the local operators
      ! <= Neumann boundaries conditions
      
      ! 0. Initializations
      ! ------------------
      
      ndkerep = 1

      ! initialization of the superstructures management

      malxm = 1
      malim = 1

      ! memory space for the pcpg associated with the FETI dual formulation
      ! ndkerep is associated to the list of rigid modes, 
      ! ndkerep == 1 because the Dual Operator
      ! is a first order operator due to SPG elliptic Operator is a
      ! second order operator

      nim = 50
      nim = nim + ndkerep
      nim = nim + 2*jpi + 2*jpj
      nim = nim + jpi*jpj

      nxm = 33
      nxm = nxm + 4*jpnij
      nxm = nxm + 19*(jpi+jpj)
      nxm = nxm + 13*jpi*jpj
      nxm = nxm + jpi*jpi*jpj

      ! krylov space memory

      iperio = 0
      IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6) iperio = 1
      nxm = nxm + 3*(jpnij-jpni)*jpi
      nxm = nxm + 3*(jpnij-jpnj+iperio)*jpj
      nxm = nxm + 2*(jpi+jpj)*(jpnij-jpni)*jpi
      nxm = nxm + 2*(jpi+jpj)*(jpnij-jpnj+iperio)*jpj

      ! Resolution with the Schur dual Method ( frontal and local solver by
      ! blocks
      ! Case with a local symetrical matrix
      ! The local matrix is stored in a multi-column form
      ! The total number of nodes for this subdomain is named "noeuds"

      noeuds = jpi*jpj
      nifmat = jpi-1
      njfmat = jpj-1
      nelem = nifmat*njfmat
      npe = 4
      nmorse = 5*noeuds
      
      ! 1. mesh building
      ! ----------------
      
      ! definition of specific information for a subdomain
      !  narea           : subdomain number = processor number +1 
      !  ninterf         : neighbour subdomain number
      !  nni             : interface point number
      !  ndvois array    : neighbour subdomain list 
      !  maplistin array : node pointer at each interface
      !  maplistin array : concatened list of interface nodes
      
      !  messag coding is necessary by interface type for avoid collision
      !  if nperio == 1 
      
      !  lint  array     : indoor interface list / type
      !  lext  array     : outdoor interface list / type
      
      !  domain with jpniXjpnj subdomains
      
      CALL feti_inisub(nifmat,njfmat,nbondi,nbondj,nperio,   &
          nbsw,nbnw,nbse,nbne,ninterf,ninterfc,nni,nnic)

      CALL feti_creadr(malim,malimax,nim,3*ninterf ,mandvois ,'ndvois' )
      CALL feti_creadr(malim,malimax,nim,3*ninterfc,mandvoisc,'ndvoisc')
      CALL feti_creadr(malim,malimax,nim,ninterfc+1,maplistin,'plistin')
      CALL feti_creadr(malim,malimax,nim,nnic      ,malistin ,'listin' )

      ! pb with periodicity conditions : iiend, ijend

      ijend = nlej
      iiend = nlei
      IF (jperio == 1) iiend = nlci - jpreci

      CALL feti_subound(nifmat,njfmat,nldi,iiend,nldj,ijend,   &
          narea,nbondi,nbondj,nperio,   &
          ninterf,ninterfc,   &
          nowe,noea,noso,nono,   &
          nbsw,nbnw,nbse,nbne,   &
          npsw,npnw,npse,npne,   &
          mfet(mandvois),mfet(mandvoisc),   &
          mfet(maplistin),nnic,mfet(malistin) )

   END SUBROUTINE fetstr


   SUBROUTINE fetmat
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE fetmat  ***
      !!            
      !! ** Purpose :   Construction of the matrix of the barotropic stream 
      !!      function system.
      !!      Finite Elements Tearing & Interconnecting (FETI) approach
      !!      Matrix treatment : Neumann condition, inverse computation
      !!
      !! ** Method :
      !!
      !! References :
      !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
      !!      A domain decomposition solver to compute the barotropic 
      !!      component of an OGCM in the parallel processing field.
      !!      Ocean Modelling, issue 105, december 94.
      !!
      !! History :
      !!        !  98-02 (M. Guyon)  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE lib_feti            ! feti librairy
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jl
      INTEGER ::   iimask(jpi,jpj)
      INTEGER ::   iiend, ijend
      REAL(wp) ::   zres, zres2, zdemi
      !!---------------------------------------------------------------------

      ! Matrix computation
      ! ------------------

      CALL feti_creadr(malxm,malxmax,nxm,nmorse,maan,'matrice a')

      nnitot = nni

      CALL mpp_sum( nnitot, 1, numit0ete )
      CALL feti_creadr(malxm,malxmax,nxm,npe*npe,maae,'ae')

      ! initialisation of the local barotropic matrix
      ! local boundary conditions on the halo

      CALL lbc_lnk( gcp(:,:,1), 'F', 1)
      CALL lbc_lnk( gcp(:,:,2), 'F', 1) 
      CALL lbc_lnk( gcp(:,:,3), 'F', 1) 
      CALL lbc_lnk( gcp(:,:,4), 'F', 1) 
      CALL lbc_lnk( gcdmat    , 'T', 1)

      ! Neumann conditions
      ! initialisation of the integer Neumann Mask

      CALL feti_iclr(jpi*jpj,iimask)
      DO jj = 1, jpj
         DO ji = 1, jpi
            iimask(ji,jj) = INT( gcdprc(ji,jj) )
         END DO
      END DO

      ! regularization of the local matrix

      DO jj = 1, jpj
         DO ji = 1, jpi
            gcdmat(ji,jj) = gcdmat(ji,jj) * gcdprc(ji,jj) + 1. - gcdprc(ji,jj)
         END DO
      END DO

      DO jk = 1, 4
         DO jj = 1, jpj
            DO ji = 1, jpi
               gcp(ji,jj,jk) = gcp(ji,jj,jk) * gcdprc(ji,jj)
            END DO
         END DO
      END DO
      
      ! implementation of the west, east, north & south Neumann conditions

      zdemi  = 0.5

      ! pb with periodicity conditions : iiend, ijend

      ijend = nlej
      iiend = nlei
      IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) iiend = nlci - jpreci

      IF( nbondi == 2 .AND. (nperio /= 1 .OR. nperio /= 4 .OR. nperio == 6) ) THEN

         ! with the periodicity : no east/west interface if nbondi = 2
         ! and nperio != 1

      ELSE 
         ! west
         IF( nbondi /= -1 ) THEN 
            DO jj = 1, jpj
               IF( iimask(1,jj) /= 0 ) THEN
                  gcp(1,jj,2) = 0.e0
                  gcp(1,jj,1) = zdemi * gcp(1,jj,1)
                  gcp(1,jj,4) = zdemi * gcp(1,jj,4)
               ENDIF
            END DO
            DO jj = 1, jpj
               IF( iimask(1,jj) /= 0 ) THEN
                  gcdmat(1,jj) = - ( gcp(1,jj,1) + gcp(1,jj,2) + gcp(1,jj,3) + gcp(1,jj,4) )
               ENDIF
            END DO
         ENDIF
         ! east
         IF( nbondi /= 1 ) THEN 
            DO jj = 1, jpj
               IF( iimask(iiend,jj) /= 0 ) THEN
                  gcp(iiend,jj,3) = 0.e0
                  gcp(iiend,jj,1) = zdemi * gcp(iiend,jj,1)
                  gcp(iiend,jj,4) = zdemi * gcp(iiend,jj,4)
               ENDIF
            END DO
            DO jj = 1, jpj
               IF( iimask(iiend,jj) /= 0 ) THEN
                  gcdmat(iiend,jj) = - ( gcp(iiend,jj,1) + gcp(iiend,jj,2)   &
                                       + gcp(iiend,jj,3) + gcp(iiend,jj,4) )
               ENDIF
            END DO
         ENDIF
      ENDIF

      ! south
      IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN 
         DO ji = 1, jpi
            IF( iimask(ji,1) /= 0 ) THEN
               gcp(ji,1,1) = 0.e0
               gcp(ji,1,2) = zdemi * gcp(ji,1,2)
               gcp(ji,1,3) = zdemi * gcp(ji,1,3)
            ENDIF
         END DO
         DO ji = 1, jpi
            IF( iimask(ji,1) /= 0 ) THEN
               gcdmat(ji,1) = - ( gcp(ji,1,1) + gcp(ji,1,2) + gcp(ji,1,3) + gcp(ji,1,4) )
            ENDIF
         END DO
      ENDIF
      
      ! north
      IF( nbondj /= 1 .AND. nbondj /= 2 ) THEN 
         DO ji = 1, jpi
            IF( iimask(ji,ijend) /= 0 ) THEN
               gcp(ji,ijend,4) = 0.e0
               gcp(ji,ijend,2) = zdemi * gcp(ji,ijend,2) 
               gcp(ji,ijend,3) = zdemi * gcp(ji,ijend,3)
            ENDIF
         END DO
         DO ji = 1, jpi
            IF( iimask(ji,ijend) /= 0 ) THEN
               gcdmat(ji,ijend) = - ( gcp(ji,ijend,1) + gcp(ji,ijend,2)   &
                                    + gcp(ji,ijend,3) + gcp(ji,ijend,4) )
            ENDIF
         END DO
      ENDIF

      ! matrix terms are  saved in FETI solver arrays
      CALL feti_vmov(noeuds,gcp(1,1,1),wfeti(maan))
      CALL feti_vmov(noeuds,gcp(1,1,2),wfeti(maan+noeuds))
      CALL feti_vmov(noeuds,gcdmat,wfeti(maan+2*noeuds))
      CALL feti_vmov(noeuds,gcp(1,1,3),wfeti(maan+3*noeuds))
      CALL feti_vmov(noeuds,gcp(1,1,4),wfeti(maan+4*noeuds))

      ! construction of Dirichlet liberty degrees array
      CALL feti_subdir(nifmat,njfmat,noeuds,ndir,iimask)
      CALL feti_creadr(malim,malimax,nim,ndir,malisdir,'lisdir')
      CALL feti_listdir(jpi,jpj,iimask,ndir,mfet(malisdir))

      ! stop onto  matrix term for Dirichlet conditions
      CALL feti_blomat(nifmat+1,njfmat+1,wfeti(maan),ndir,mfet(malisdir))

      ! reservation of factorized diagonal blocs and temporary array for
      ! factorization
      npblo = (njfmat+1) * (nifmat+1) * (nifmat+1)
      ndimax = nifmat+1

      CALL feti_creadr(malxm,malxmax,nxm,npblo,mablo,'blo')
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,madia,'dia')
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,mav,'v')
      CALL feti_creadr(malxm,malxmax,nxm,ndimax*ndimax,mautil,'util')

      ! stoping the rigid modes

      ! the number of rigid modes =< Max [dim(Ker(Ep))]
      !                                p=1,Np

      CALL feti_creadr(malim,malimax,nim,ndkerep,malisblo,'lisblo')

      ! Matrix factorization

      CALL feti_front(noeuds,nifmat+1,njfmat+1,wfeti(maan),npblo,   &
          wfeti(mablo),wfeti(madia),   &
          wfeti(mautil),wfeti(mav),ndlblo,mfet(malisblo),ndkerep)
      CALL feti_prext(noeuds,wfeti(madia))

      ! virtual dealloc => we have to see for a light f90 version
      ! the super structure is removed to clean the coarse grid
      ! solver structure

      malxm = madia
      CALL feti_vclr(noeuds,wfeti(madia))
      CALL feti_vclr(noeuds,wfeti(mav))
      CALL feti_vclr(ndimax*ndimax,wfeti(mautil))

      ! ndlblo is the dimension of the local nullspace .=<. the size of the 
      ! memory of the superstructure associated to the nullspace : ndkerep
      ! ndkerep is introduced to avoid messages "out of bounds" when memory
      ! is checked

      ! copy matrix for Dirichlet condition

      CALL feti_creadr(malxm,malxmax,nxm,noeuds,miax,'x')
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,may,'y')
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,maz,'z')

      ! stoping the rigid modes

      ! ndlblo is the dimension of the local nullspace .=<. the size of the 
      ! memory of the superstructure associated to the nullspace : ndkerep
      ! ndkerep is introduced to avoid messages "out of bounds" when memory
      ! is checked

      CALL feti_creadr(malxm,malxmax,nxm,ndkerep*noeuds,mansp,'nsp')
      CALL feti_blomat1(nifmat+1,njfmat+1,wfeti(maan),ndlblo,   &
          mfet(malisblo),wfeti(mansp))      

      ! computation of operator kernel

      CALL feti_nullsp(noeuds,nifmat+1,njfmat+1,npblo,wfeti(mablo),   &
          wfeti(maan),ndlblo,mfet(malisblo),wfeti(mansp),   &
          wfeti(maz))

      ! test of the factorisation onto each sub domain

      CALL feti_init(noeuds,wfeti(may))
      CALL feti_blodir(noeuds,wfeti(may),ndir,mfet(malisdir))
      CALL feti_blodir(noeuds,wfeti(may),ndlblo,mfet(malisblo))
      CALL feti_vclr(noeuds,wfeti(miax))
      CALL feti_resloc(noeuds,nifmat+1,njfmat+1,wfeti(maan),npblo,   &
          wfeti(mablo),wfeti(may),wfeti(miax),wfeti(maz)) 
      CALL feti_proax(noeuds,nifmat+1,njfmat+1,wfeti(maan),wfeti(miax),   &
          wfeti(maz))
      CALL feti_blodir(noeuds,wfeti(maz),ndlblo,mfet(malisblo))
      CALL feti_vsub(noeuds,wfeti(may),wfeti(maz),wfeti(maz))

      zres2 = 0.e0
      DO jl = 1, noeuds
         zres2 = zres2 + wfeti(may+jl-1) * wfeti(may+jl-1)
      END DO
      CALL mpp_sum(zres2,1,zres)

      res2 = 0.e0
      DO jl = 1, noeuds
         res2 = res2 + wfeti(maz+jl-1) * wfeti(maz+jl-1)
      END DO 
      res2 = res2 / zres2
      CALL mpp_sum(res2,1,zres)

      res2 = SQRT(res2)
      IF(lwp) WRITE(numout,*) 'global residu : sqrt((Ax-b,Ax-b)/(b.b)) =', res2

      IF( res2 > (eps/100.) ) THEN 
         IF(lwp) WRITE (numout,*) 'eps is :',eps
         IF(lwp) WRITE (numout,*) 'factorized matrix precision :',res2
         STOP 
      ENDIF

   END SUBROUTINE fetmat


   SUBROUTINE fetsch
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE fetsch  ***
      !!        
      !! ** Purpose :
      !!     Construction of the matrix of the barotropic stream function
      !!     system.
      !!     Finite Elements Tearing & Interconnecting (FETI) approach
      !!     Data framework for the Schur Dual solve
      !!
      !! ** Method :
      !!
      !! References :
      !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
      !!      A domain decomposition solver to compute the barotropic 
      !!      component of an OGCM in the parallel processing field.
      !!      Ocean Modelling, issue 105, december 94.
      !!
      !! History :
      !!        !  98-02 (M. Guyon)  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE lib_feti            ! feti librairy
      !! * Local declarations
      !!---------------------------------------------------------------------

      !  computing weights for the conform construction

      CALL feti_creadr(malxm,malxmax,nxm,noeuds,mapoids ,'poids' )
      CALL feti_creadr(malxm,malxmax,nxm,nnic  ,mabufin ,'bufin' )
      CALL feti_creadr(malxm,malxmax,nxm,nnic  ,mabufout,'bufout')

!!    CALL feti_poids(ninterfc,mfet(mandvoisc),mfet(maplistin),nnic,   &
!!        mfet(malistin),narea,noeuds,wfeti(mapoids),wfeti(mabufin),   &
!!        wfeti(mabufout) )
      CALL feti_poids(ninterfc,                                nnic,   &
          mfet(malistin),      noeuds,wfeti(mapoids)                )


      ! Schur dual arrays
 
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,mabitw,'bitw') 
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,mautilu,'utilu') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,malambda,'lambda') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,mag,'g') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,mapg,'pg') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,mamg,'mg') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,maw,'w') 
      CALL feti_creadr(malxm,malxmax,nxm,nni,madw,'dw')

      !  coarse grid solver dimension and arrays

      nitmaxete = ndlblo
      CALL  mpp_sum(nitmaxete,1,numit0ete)

      nitmaxete = nitmaxete + 1
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maxnul,'xnul')
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maynul,'ynul')
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maeteg,'eteg')
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep,maeteag,'eteag')
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep*nitmaxete,maeted,'eted')
      CALL feti_creadr(malxm,malxmax,nxm,ndkerep*nitmaxete,maetead,'etead')
      CALL feti_creadr(malxm,malxmax,nxm,nitmaxete,maeteadd,'eteadd')
      CALL feti_creadr(malxm,malxmax,nxm,nitmaxete,maetegamm,'etegamm')
      CALL feti_creadr(malxm,malxmax,nxm,nni,maetew,'etew') 
      CALL feti_creadr(malxm,malxmax,nxm,noeuds,maetev,'etev') 

      ! construction of semi interface arrays

      CALL feti_creadr(malim,malimax,nim,ninterf+1,maplistih,'plistih')
!!    CALL feti_halfint(ninterf,mfet(mandvois),mfet(maplistin),nni,   &
!!        mfet(maplistih),nnih,narea)
      CALL feti_halfint(ninterf,mfet(mandvois),mfet(maplistin),       & 
          mfet(maplistih),nnih      )

      CALL feti_creadr(malxm,malxmax,nxm,nnih,magh,'gh') 

      ! Schur Dual Method

      nmaxd = nnitot / 2

      !  computation of the remain array for descent directions

      nmaxd = min0(nmaxd,(nxm-nitmaxete-malxm)/(2*nnih+3))
      CALL mpp_min(nmaxd,1,numit0ete)

      nitmax = nnitot/2
      epsilo = eps
      ntest = 0

      ! Krylov space construction

      CALL feti_creadr(malxm,malxmax,nxm,nnih*nmaxd,mawj,'wj') 
      CALL feti_creadr(malxm,malxmax,nxm,nnih*nmaxd,madwj,'dwj') 
      CALL feti_creadr(malxm,malxmax,nxm,nmaxd,madwwj,'dwwj') 
      CALL feti_creadr(malxm,malxmax,nxm,nmaxd,magamm,'gamm') 
      CALL feti_creadr(malxm,malxmax,nxm,max0(nmaxd,nitmaxete),mawork,'work')
      mjj0 = 0 
      numit0ete = 0 

   END SUBROUTINE fetsch

#else
   SUBROUTINE fetstr                 ! Empty routine
   END SUBROUTINE fetstr
   SUBROUTINE fetmat                 ! Empty routine
   END SUBROUTINE fetmat
   SUBROUTINE fetsch                 ! Empty routine
   END SUBROUTINE fetsch
#endif

   !!======================================================================
END MODULE solmat
