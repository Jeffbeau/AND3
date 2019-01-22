MODULE solpcg
   !!======================================================================
   !!                     ***  MODULE  solfet
   !! Ocean solver :  preconditionned conjugate gradient solver
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sol_pcg    : preconditionned conjugate gradient solver
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean solver variables
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sol_pcg              ! ???

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solpcg.F90,v 1.4 2005/12/21 10:46:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sol_pcg( kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_pcg  ***
      !!                    
      !! ** Purpose :   Solve the ellipic equation for the barotropic stream 
      !!      function system (lk_dynspg_rl=T) or the transport divergence 
      !!      system (lk_dynspg_flt=T) using a diagonal preconditionned
      !!      conjugate gradient method.
      !!      In the former case, the barotropic stream function trend has a
      !!      zero boundary condition along all coastlines (i.e. continent
      !!      as well as islands) while in the latter the boundary condition
      !!      specification is not required.
      !!
      !! ** Method  :   Diagonal preconditionned conjugate gradient method.
      !!      the algorithm is multitasked. (case of 5 points matrix)
      !!      define              pa  = q^-1 * a
      !!                        pgcb  = q^-1 * gcb
      !!                 < . ; . >_q  = ( . )^t q ( . )
      !!      where q is the preconditioning matrix = diagonal matrix of the
      !!                                              diagonal elements of a
      !!      Initialization:
      !!         x(o) = gcx
      !!         r(o) = d(o) = pgcb - pa.x(o)
      !!         rr(o)= < r(o) , r(o) >_q
      !!      Iteration n   :
      !!         z(n)   = pa.d(n)
      !!         alp(n) = rr(n) / < z(n) , d(n) >_q
      !!         x(n+1) = x(n) + alp(n) d(n)
      !!         r(n+1) = r(n) - alp(n) z(n)
      !!         rr(n+1)= < r(n+1) , r(n+1) >_q
      !!         bet(n) = rr(n+1) / rr(n)
      !!         r(n+1) = r(n+1) + bet(n+1) d(n)
      !!      Convergence test :
      !!         rr(n+1) / < gcb , gcb >_q   =< epsr
      !!
      !! ** Action : - niter  : solver number of iteration done
      !!             - res    : solver residu reached
      !!             - gcx()  : solution of the elliptic system
      !!
      !! References :
      !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
      !!
      !! History :
      !!        !  90-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-04  (M. Guyon)  loops and suppress pointers
      !!        !  95-09  (M. Imbard, J. Escobar)  mpp exchange 
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step
      !                                      ! set to zero before the call of solpcg

      !! * Local declarations
      INTEGER ::   ji, jj, jn                ! dummy loop indices
      REAL(wp) ::   zgcad                    ! temporary scalars
      !!----------------------------------------------------------------------

      !                                                !================
      DO jn = 1, nmax                                  ! Iterative loop
         !                                             !================

         IF( jn == 1 ) THEN           ! Initialization of the algorithm

            CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! lateral boundary condition

            ! gcr   = gcb-a.gcx
            ! gcdes = gsr
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zgcad = bmask(ji,jj) * ( gcb(ji,jj  ) -                gcx(ji  ,jj  )   &
                     &                                  - gcp(ji,jj,1) * gcx(ji  ,jj-1)   &
                     &                                  - gcp(ji,jj,2) * gcx(ji-1,jj  )   &
                     &                                  - gcp(ji,jj,3) * gcx(ji+1,jj  )   &
                     &                                  - gcp(ji,jj,4) * gcx(ji  ,jj+1)   )
                  gcr  (ji,jj) = zgcad
                  gcdes(ji,jj) = zgcad
               END DO
            END DO
            
            rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
            IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
            rr = rnorme

         ENDIF
        
         !                             ! Algorithm
        
         CALL lbc_lnk( gcdes, c_solver_pt, 1. )   ! lateral boundary condition
        
         ! ... gccd = matrix . gcdes
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gccd(ji,jj) = bmask(ji,jj)*( gcdes(ji,jj)   &
                  &        +gcp(ji,jj,1)*gcdes(ji,jj-1)+gcp(ji,jj,2)*gcdes(ji-1,jj)   &
                  &        +gcp(ji,jj,4)*gcdes(ji,jj+1)+gcp(ji,jj,3)*gcdes(ji+1,jj)   )
            END DO
         END DO
 
         ! alph = (gcr,gcr)/(gcdes,gccd)
         radd = SUM(  gcdes(:,:) * gcdmat(:,:) * gccd(:,:)  )
         IF( lk_mpp )   CALL mpp_sum( radd )   ! sum over the global domain
         alph = rr / radd
         
         ! gcx = gcx + alph * gcdes
         ! gcr = gcr - alph * gccd
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gcx(ji,jj) = bmask(ji,jj) * ( gcx(ji,jj) + alph * gcdes(ji,jj) )
               gcr(ji,jj) = bmask(ji,jj) * ( gcr(ji,jj) - alph * gccd (ji,jj) )
            END DO
         END DO
        
         ! rnorme = (gcr,gcr)
         rnorme = SUM(  gcr(:,:) * gcdmat(:,:) * gcr(:,:)  )
         IF( lk_mpp )   CALL  mpp_sum( rnorme )   ! sum over the global domain
        
         ! test of convergence
         IF( rnorme < epsr .OR. jn == nmax ) THEN
            res = SQRT( rnorme )
            niter = jn
            ncut = 999
         ENDIF
        
         ! beta = (rk+1,rk+1)/(rk,rk)
         beta = rnorme / rr
         rr   = rnorme

         ! indicator of non-convergence or explosion
         IF( jn == nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999

         ! gcdes = gcr + beta * gcdes
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               gcdes(ji,jj) = bmask(ji,jj)*( gcr(ji,jj) + beta * gcdes(ji,jj) )
            END DO
         END DO
        
         !                                             !================
      END DO                                           !    End Loop
      !                                                !================
     
999   CONTINUE
     
     
      ! Output in gcx with lateral b.c. applied
      ! ---------------------------------------
     
      CALL lbc_lnk( gcx, c_solver_pt, 1. )
     
   END SUBROUTINE sol_pcg

   !!=====================================================================
END MODULE solpcg
