MODULE solsor
   !!======================================================================
   !!                     ***  MODULE  solsor  ***
   !! Ocean solver :  Successive Over-Relaxation solver
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sol_sor     : Red-Black Successive Over-Relaxation solver
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE sol_oce         ! solver variables
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sol_sor              ! ???

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor.F90,v 1.8 2005/12/21 10:46:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
      
   SUBROUTINE sol_sor( kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_sor  ***
      !!                 
      !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
      !!      function system (lk_dynspg_rl=T) or the transport divergence 
      !!      system (lk_dynspg_flt=T) using a red-black successive-over-
      !!      relaxation method.
      !!       In the former case, the barotropic stream function trend has a
      !!     zero boundary condition along all coastlines (i.e. continent
      !!     as well as islands) while in the latter the boundary condition
      !!     specification is not required.
      !!
      !! ** Method  :   Successive-over-relaxation method using the red-black 
      !!      technique. The former technique used was not compatible with 
      !!      the north-fold boundary condition used in orca configurations.
      !!
      !! References :
      !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
      !!
      !! History :
      !!        !  90-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!   7.1  !  93-04  (G. Madec)  time filter
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step
      !                                      ! set to zero before the call of solsor
      !! * Local declarations
      INTEGER  ::   ji, jj, jn               ! dummy loop indices
      INTEGER  ::   ishift
      REAL(wp) ::   ztmp, zres, zres2

      INTEGER  ::   ijmppodd, ijmppeven
      !!----------------------------------------------------------------------
      
      ijmppeven = MOD(nimpp+njmpp  ,2)
      ijmppodd  = MOD(nimpp+njmpp+1,2)
      !                                                       ! ==============
      DO jn = 1, nmax                                         ! Iterative loop 
         !                                                    ! ==============

         CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boundary conditions
         
         ! Residus
         ! -------

         ! Guess black update
         DO jj = 2, jpjm1
            ishift = MOD( jj-ijmppodd, 2 )
            DO ji = 2+ishift, jpim1, 2
               ztmp =                  gcb(ji  ,jj  )   &
                  &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
                  &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
                  &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
                  &   - gcp(ji,jj,4) * gcx(ji  ,jj+1)
               ! Estimate of the residual
               zres = ztmp - gcx(ji,jj)
               gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
               ! Guess update
               gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
            END DO
         END DO

         CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! applied the lateral boubary conditions

         ! Guess red update
         DO jj = 2, jpjm1
            ishift = MOD( jj-ijmppeven, 2 )
            DO ji = 2+ishift, jpim1, 2
               ztmp =                  gcb(ji  ,jj  )   &
                  &   - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
                  &   - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
                  &   - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
                  &   - gcp(ji,jj,4) * gcx(ji  ,jj+1) 
               ! Estimate of the residual
               zres = ztmp - gcx(ji,jj)
               gcr(ji,jj) = zres * gcdmat(ji,jj) * zres
               ! Guess update
               gcx(ji,jj) = sor * ztmp + (1-sor) * gcx(ji,jj)
            END DO
         END DO

         ! test of convergence
         IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then

            SELECT CASE ( nsol_arp )
            CASE ( 0 )                 ! absolute precision (maximum value of the residual)
               zres2 = MAXVAL( gcr(2:jpim1,2:jpjm1) )
               IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
               ! test of convergence
               IF( zres2 < resmax .OR. jn == nmax ) THEN
                  res = SQRT( zres2 )
                  niter = jn
                  ncut = 999
               ENDIF
            CASE ( 1 )                 ! relative precision
               rnorme = SUM( gcr(2:jpim1,2:jpjm1) )
               IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
               ! test of convergence
               IF( rnorme < epsr .OR. jn == nmax ) THEN
                  res = SQRT( rnorme )
                  niter = jn
                  ncut = 999
               ENDIF
            END SELECT
         
         !****
9300     FORMAT('          niter :',i4,' res :',e20.10,' b :',e20.10)
         !****
         
         ENDIF
         ! indicator of non-convergence or explosion
         IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999
         
         !                                                 ! =====================
      END DO                                               ! END of iterative loop
      !                                                    ! =====================
      
999   CONTINUE
      
      
      !  Output in gcx
      !  -------------

      CALL lbc_lnk( gcx, c_solver_pt, 1. )    ! boundary conditions

      
   END SUBROUTINE sol_sor

   !!=====================================================================
END MODULE solsor
