MODULE solfet
   !!======================================================================
   !!                     ***  MODULE  solfet
   !! Ocean solver :  Finite Elements Tearing & Interconnecting solver
   !!=====================================================================
#if defined key_feti
   !!----------------------------------------------------------------------
   !!   'key_feti' :                                            FETI solver
   !!----------------------------------------------------------------------
   !!   sol_fet     : FETI solver
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean solver 
   USE lib_mpp         ! distribued memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC sol_fet              ! ???
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solfet.F90,v 1.4 2005/12/21 10:46:41 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sol_fet( kindic )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sol_fet  ***
      !!
      !! ** Purpose :
      !!     Solve the ellipic equation for the barotropic stream function
      !!     system (default option) or the transport divergence system
      !!     (lk_dynspg_flt=T) using a Finite Elements Tearing and 
      !!      Interconnecting (FETI) approach.
      !!     In the former case, the barotropic stream function trend has a
      !!     zero boundary condition along all coastlines (i.e. continent
      !!     as well as islands) while in the latter the boundary condition
      !!     specification is not required.
      !!
      !! ** Method :
      !!      Resolution of the elliptic equation by a Dual formulation of
      !!      the Schur Complement Method or Finite Elements Tearing & 
      !!      Interconnecting (FETI) approach
      !!
      !! ** Action :
      !!
      !! ** References :
      !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
      !!      A domain decomposition solver to compute the barotropic 
      !!      component of an OGCM in the parallel processing field.
      !!      Ocean Modelling, issue 105, december 94.
      !!
      !! History :
      !!        !  97-02  (M. Guyon)  original code 
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step

      !! * Local declarations
      INTEGER ::   ji, jj                    ! dummy loop indices
      INTEGER ::   iit0, inn, iju 
      REAL(wp) ::   zgcad, zgwgt


      ! Norme(gcr)   = (gcb, gcb)
      ! gcdes = gsr
      
      ! bmask field is filtering the differents contribution on the 
      ! non-overlapping interface
      
      gcb(:,:) = bmask(:,:) * gcb(:,:)
      
      ! Mpp: sum over all the global domain
      
      ! copy the right hand member
      
      CALL feti_vmov(noeuds,gcb(1,1),wfeti(may))
      
      ! conservation of descent direction if ntest = 0
      
      IF(ntest /= 0) mjj0=0
      iit0 = mjj0
      !                      --->
      !    resolution of the Grad(PS) equation by a Dual formulation of the
      !    Schur Complement Method or Finite Elements Tearing & Interconnecting
      !    (FETI) approach
      !    interface problem (Lagrange  multiplier) : PCPG algorithm
      !    local problem (trend of the 2D potential field) : LU factorization
      !    preconditioner : lumped
      !    optimisation : Krylov initialisation + Krylov correction

      CALL feti_dualschur(noeuds,nifmat+1,njfmat+1,wfeti(maan),   &
          npblo,wfeti(mablo),ninterf,   &
          ninterfc,nni,nnic,mfet(mandvois),mfet(mandvoisc),   &
          mfet(maplistin),mfet(malistin),   &
          wfeti(mapoids),wfeti(miax),   &
          wfeti(maz),wfeti(may),   &
          wfeti(mabitw),wfeti(mautilu),   &
          wfeti(malambda),wfeti(mag),   &
          wfeti(mapg),wfeti(mamg),nitmax,nmaxd,mjj0,   &
          wfeti(mawj),   &
          wfeti(madwj),wfeti(madwwj),wfeti(magamm),   &
          wfeti(mawork),   &
          wfeti(mabufin),wfeti(mabufout),narea,epsilo,   &
          ndlblo,mfet(malisblo),ndkerep,   &
          wfeti(maxnul),wfeti(maynul),numit0ete,nitmaxete,   &
          wfeti(maeteg),wfeti(maeteag),wfeti(maeted),   &
          wfeti(maetead),wfeti(maeteadd),wfeti(maetegamm),   &
          wfeti(mansp),   &
          wfeti(maetev),wfeti(maetew),nnih,mfet(maplistih),   &
          wfeti(magh),   &
          wfeti(maw),wfeti(madw),   &
          res,kindic,inn)

      ! number of iteration of the pcg to solve the interface pb

      inn =  mjj0 - iit0

      ! test of convergence
      IF( res < epsilo .OR. inn == nmax ) THEN
          niter = inn
          ncut  = 999
      ENDIF

      ! indicator of non-convergence or explosion
      IF( inn == nmax .OR. rr > 1.e+20 ) kindic = -2
      IF( ncut == 999 ) GOTO 999

999   CONTINUE

      !  2. Output in gcx
      !  -----------------
      
      CALL feti_vmov( noeuds, wfeti(miax), gcx )

      CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! lateral boundary condition

   END SUBROUTINE sol_fet

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE sol_fet( kindic )        ! Empty routine
      INTEGER, INTENT( inout ) ::   kindic  ! solver problem 
      kindic = -100
   END SUBROUTINE sol_fet
#endif

   !!=====================================================================
END MODULE solfet
