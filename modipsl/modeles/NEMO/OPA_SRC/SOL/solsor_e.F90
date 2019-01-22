MODULE solsor_e
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
   PUBLIC sol_sor_e              ! ???

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solsor_e.F90,v 1.2 2005/12/21 10:46:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS
      
   SUBROUTINE sol_sor_e( kindic )
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
      !!       This routine provides a MPI optimization to the existing solsor
      !!     by reducing the number of call to lbc.
      !! 
      !! ** Method  :   Successive-over-relaxation method using the red-black 
      !!      technique. The former technique used was not compatible with 
      !!      the north-fold boundary condition used in orca configurations.
      !!      Compared to the classical sol_sor, this routine provides a 
      !!      mpp optimization by reducing the number of calls to lnc_lnk
      !!      The solution is computed on a larger area and the boudary
      !!      conditions only when the inside domain is reached.
      !! 
      !! References :
      !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
      !!      Beare and Stevens 1997 Ann. Geophysicae 15, 1369-1377
      !!
      !! History :
      !!        !  90-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!   7.1  !  93-04  (G. Madec)  time filter
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
      !!   9.0  !  05-09  (R. Benshila, G. Madec)  MPI optimization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step
      !                                      ! set to zero before the call of solsor
      !! * Local declarations
      INTEGER  ::   ji, jj, jn               ! dummy loop indices
      INTEGER  ::   ishift, icount
      REAL(wp) ::   ztmp, zres, zres2

      INTEGER  ::   ijmppodd, ijmppeven
      INTEGER  ::   ijpr2d
      !!----------------------------------------------------------------------
      
      ijmppeven = MOD(nimpp+njmpp+jpr2di+jpr2dj,2)
      ijmppodd  = MOD(nimpp+njmpp+jpr2di+jpr2dj+1,2)
      ijpr2d = MAX(jpr2di,jpr2dj)
      icount = 0
      !                                                       ! ==============
      DO jn = 1, nmax                                         ! Iterative loop 
         !                                                    ! ==============

         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )   
        
         ! Residus
         ! -------

         ! Guess black update
         DO jj = 2-jpr2dj, nlcj-1+jpr2dj
            ishift = MOD( jj-ijmppodd-jpr2dj, 2 )
            DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
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
         icount = icount + 1 
 
         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx, c_solver_pt, 1. )  

         ! Guess red update
         DO jj = 2-jpr2dj, nlcj-1+jpr2dj
            ishift = MOD( jj-ijmppeven-jpr2dj, 2 )
            DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
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
         icount = icount + 1

         ! test of convergence
         IF ( jn > nmin .AND. MOD( jn-nmin, nmod ) == 0 ) then

            SELECT CASE ( nsol_arp )
            CASE ( 0 )                 ! absolute precision (maximum value of the residual)
               zres2 = MAXVAL( gcr(2:nlci-1,2:nlcj-1) )
               IF( lk_mpp )   CALL mpp_max( zres2 )   ! max over the global domain
               ! test of convergence
               IF( zres2 < resmax .OR. jn == nmax ) THEN
                  res = SQRT( zres2 )
                  niter = jn
                  ncut = 999
               ENDIF
            CASE ( 1 )                 ! relative precision
               rnorme = SUM( gcr(2:nlci-1,2:nlcj-1) )
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

      CALL lbc_lnk_e( gcx, c_solver_pt, 1. )    ! boundary conditions

      
   END SUBROUTINE sol_sor_e

   !!=====================================================================
END MODULE solsor_e
