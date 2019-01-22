MODULE sol_oce
   !!======================================================================
   !!                    ***  MODULE  sol_oce  ***
   !! Ocean solver :  solver variables defined in memory 
   !!=====================================================================
   !!
   !! ** Purpose :   Define in memory solver variables
   !!
   !! History :
   !!   9.0  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/sol_oce.F90,v 1.8 2006/03/10 10:55:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   !!-----------------------------------
   !! elliptic solver: SOR, PCG or FETI
   !! ----------------------------------
   INTEGER , PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
      nsolv    =    1 ,     &  !: = 1/2/3/4 type of elliptic solver
      nsol_arp =    0 ,     &  !: = 0/1 absolute/relative precision convergence test
      nmin     =  300 ,     &  !: minimum of iterations for the SOR solver
      nmax     =  800 ,     &  !: maximum of iterations for the SOR solver
      nmod     =   10 ,     &  !: frequency of test for the SOR solver
      nmisl    = 4000          !: maximum pcg iterations for island
     
   REAL(wp), PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
      eps    =  1.e-6_wp ,  &  !: absolute precision of the solver
      resmax = 1.e-14_wp ,  &  !: absolute precision for the SOR solver
      sor    =   1.92_wp ,  &  !: optimal coefficient for the SOR solver
      epsisl = 1.e-10_wp ,  &  !: absolute precision on stream function solver
      rnu    =    1.0_wp       !: strength of the additional force used in free surface

   CHARACTER(len=1), PUBLIC ::   &  !:
      c_solver_pt = 'T'        !: nature of grid-points T (S) for free surface case
      !                        !                        F (G) for rigid-lid case

   INTEGER , PUBLIC ::   &  !:
      ncut,         &  !: indicator of solver convergence
      niter            !: number of iteration done by the solver

   REAL(wp), PUBLIC ::   &  !:
      epsr,         &  !: relative precision for SOR & PCG solvers
      epsilo,       &  !: precision for the FETI solver
      rnorme, res,  &  !: intermediate modulus, solver residu
      alph,         &  !: coefficient  =(gcr,gcr)/(gcx,gccd)
      beta,         &  !: coefficient  =(rn+1,rn+1)/(rn,rn)
      radd,         &  !: coefficient  =(gccd,gcdes)
      rr               !: coefficient  =(rn,rn)

   REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &  !:
      gcp              !: barotropic matrix extra-diagonal elements

   REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj) ::   &  !:
      gcx, gcxb,    &  !: now, before solution of the elliptic equation
      gcdprc,       &  !: inverse diagonal preconditioning matrix
      gcdmat,       &  !: diagonal preconditioning matrix
      gcb,          &  !: second member of the barotropic linear system
      gcr,          &  !: residu =b-a.x
      gcdes,        &  !: vector descente
      gccd             !: vector such that ca.gccd=a.d (ca-1=gcdprc)

#if defined key_agrif
      REAL(wp), DIMENSION(jpi,jpj) :: laplacu, laplacv
#endif

#if defined key_feti
   !!----------------------------------------------------------------------
   !!   'key_feti' :                                            FETI solver
   !!----------------------------------------------------------------------
   !!      noeuds           : total number of nodes for a subdomnain
   !!      ninterf          : neighbour subdomain number
   !!      nni              : interface point number
   !!      ndvois()         : neighbour subdomain list
   !!      maplistin()      : node pointer at each interface
   !!      malistin()       : concatened list of interface nodes

   INTEGER, PUBLIC :: nim,nxm,   &
       malxm,malim,malxmax,malimax,   &
       nifmat,njfmat,nelem,npe,matopo,   &
       noeuds,nmorse,maan,   &
       ninterf,ninterfc,nni,nnic,nnih,nnitot,ndir,   &
       mandvois,mandvoisc,   &
       maplistin,maplistih,malistin,   &
       malisdir,npblo,mablo,ndlblo,ndkerep,   &
       malisblo,mansp,ndimax,   &
       miax,may,maz,mapoids,   &
       nmaxd,mjj0,nitmax,ntest,   &
       mabitw,mautilu,malambda,mag,mamg,mapg,mawj,madwj,madwwj,   &
       magh,maw,madw,   &
       mautil,mav,madia,mabufin,mabufout,mawork,maae,magamm,   &
       maxnul,maynul,numit0ete,nitmaxete,maeteg,maeteag,   &
       maeted,maetead,maeteadd,maetegamm,maetev,maetew,   &
       madwork

   INTEGER, PUBLIC :: mfet(jpi*jpj+2*jpi+2*jpj+51)

   REAL(wp), PUBLIC ::   &  !:
      wfeti(jpj*jpi*jpi+13*jpi*jpj+19*(jpi+jpj)   &
       +4*jpnij+33   &
       +2*(jpi+jpj)*(jpnij-jpni)*jpi   &
       +2*(jpi+jpj)*(jpnij-jpnj+jperio)*jpj   &
       +3*(jpnij-jpni)*jpi   &
       +3*(jpnij-jpnj+jperio)*jpj) 

   REAL(wp), PUBLIC ::   res2, rcompt

#endif

   !!----------------------------------------------------------------------
END MODULE sol_oce
