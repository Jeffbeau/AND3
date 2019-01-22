


PGF90 (Version     10.5)          09/23/2010  09:49:24      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: sol_oce.F90

(    1) MODULE sol_oce
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  sol_oce  ***
(    4)    !! Ocean solver :  solver variables defined in memory 
(    5)    !!=====================================================================
(    6)    !!
(    7)    !! ** Purpose :   Define in memory solver variables
(    8)    !!
(    9)    !! History :
(   10)    !!   9.0  !  02-11  (G. Madec)  F90: Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/sol_oce.F90,v 1.8 2006/03/10 10:55:43 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !!-----------------------------------
(   23)    !! elliptic solver: SOR, PCG or FETI
(   24)    !! ----------------------------------
(   25)    INTEGER , PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
(   26)       nsolv    =    1 ,     &  !: = 1/2/3/4 type of elliptic solver
(   27)       nsol_arp =    0 ,     &  !: = 0/1 absolute/relative precision convergence test
(   28)       nmin     =  300 ,     &  !: minimum of iterations for the SOR solver
(   29)       nmax     =  800 ,     &  !: maximum of iterations for the SOR solver
(   30)       nmod     =   10 ,     &  !: frequency of test for the SOR solver
(   31)       nmisl    = 4000          !: maximum pcg iterations for island
(   32)      
(   33)    REAL(wp), PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
(   34)       eps    =  1.e-6_wp ,  &  !: absolute precision of the solver
(   35)       resmax = 1.e-14_wp ,  &  !: absolute precision for the SOR solver
(   36)       sor    =   1.92_wp ,  &  !: optimal coefficient for the SOR solver
(   37)       epsisl = 1.e-10_wp ,  &  !: absolute precision on stream function solver
(   38)       rnu    =    1.0_wp       !: strength of the additional force used in free surface
(   39) 
(   40)    CHARACTER(len=1), PUBLIC ::   &  !:
(   41)       c_solver_pt = 'T'        !: nature of grid-points T (S) for free surface case
(   42)       !                        !                        F (G) for rigid-lid case
(   43) 
(   44)    INTEGER , PUBLIC ::   &  !:
(   45)       ncut,         &  !: indicator of solver convergence
(   46)       niter            !: number of iteration done by the solver






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 2

(   47) 
(   48)    REAL(wp), PUBLIC ::   &  !:
(   49)       epsr,         &  !: relative precision for SOR & PCG solvers
(   50)       epsilo,       &  !: precision for the FETI solver
(   51)       rnorme, res,  &  !: intermediate modulus, solver residu
(   52)       alph,         &  !: coefficient  =(gcr,gcr)/(gcx,gccd)
(   53)       beta,         &  !: coefficient  =(rn+1,rn+1)/(rn,rn)
(   54)       radd,         &  !: coefficient  =(gccd,gcdes)
(   55)       rr               !: coefficient  =(rn,rn)
(   56) 
(   57)    REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &  !:
(   58)       gcp              !: barotropic matrix extra-diagonal elements
(   59) 
(   60)    REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj) ::   &  !:
(   61)       gcx, gcxb,    &  !: now, before solution of the elliptic equation
(   62)       gcdprc,       &  !: inverse diagonal preconditioning matrix
(   63)       gcdmat,       &  !: diagonal preconditioning matrix
(   64)       gcb,          &  !: second member of the barotropic linear system
(   65)       gcr,          &  !: residu =b-a.x
(   66)       gcdes,        &  !: vector descente
(   67)       gccd             !: vector such that ca.gccd=a.d (ca-1=gcdprc)
(   68) 
(   69) #if defined key_agrif
(   70)       REAL(wp), DIMENSION(jpi,jpj) :: laplacu, laplacv
(   71) #endif
(   72) 
(   73) #if defined key_feti
(   74)    !!----------------------------------------------------------------------
(   75)    !!   'key_feti' :                                            FETI solver
(   76)    !!----------------------------------------------------------------------
(   77)    !!      noeuds           : total number of nodes for a subdomnain
(   78)    !!      ninterf          : neighbour subdomain number
(   79)    !!      nni              : interface point number
(   80)    !!      ndvois()         : neighbour subdomain list
(   81)    !!      maplistin()      : node pointer at each interface
(   82)    !!      malistin()       : concatened list of interface nodes
(   83) 
(   84)    INTEGER, PUBLIC :: nim,nxm,   &
(   85)        malxm,malim,malxmax,malimax,   &
(   86)        nifmat,njfmat,nelem,npe,matopo,   &
(   87)        noeuds,nmorse,maan,   &
(   88)        ninterf,ninterfc,nni,nnic,nnih,nnitot,ndir,   &
(   89)        mandvois,mandvoisc,   &
(   90)        maplistin,maplistih,malistin,   &
(   91)        malisdir,npblo,mablo,ndlblo,ndkerep,   &
(   92)        malisblo,mansp,ndimax,   &
(   93)        miax,may,maz,mapoids,   &
(   94)        nmaxd,mjj0,nitmax,ntest,   &
(   95)        mabitw,mautilu,malambda,mag,mamg,mapg,mawj,madwj,madwwj,   &
(   96)        magh,maw,madw,   &
(   97)        mautil,mav,madia,mabufin,mabufout,mawork,maae,magamm,   &
(   98)        maxnul,maynul,numit0ete,nitmaxete,maeteg,maeteag,   &
(   99)        maeted,maetead,maeteadd,maetegamm,maetev,maetew,   &
(  100)        madwork
(  101) 
(  102)    INTEGER, PUBLIC :: mfet(jpi*jpj+2*jpi+2*jpj+51)
(  103) 
(  104)    REAL(wp), PUBLIC ::   &  !:






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 3

(  105)       wfeti(jpj*jpi*jpi+13*jpi*jpj+19*(jpi+jpj)   &
(  106)        +4*jpnij+33   &
(  107)        +2*(jpi+jpj)*(jpnij-jpni)*jpi   &
(  108)        +2*(jpi+jpj)*(jpnij-jpnj+jperio)*jpj   &
(  109)        +3*(jpnij-jpni)*jpi   &
(  110)        +3*(jpnij-jpnj+jperio)*jpj) 
(  111) 
(  112)    REAL(wp), PUBLIC ::   res2, rcompt
(  113) 
(  114) #endif
(  115) 
(  116)    !!----------------------------------------------------------------------
(  117) END MODULE sol_oce



















































PGF90 (Version     10.5)          09/23/2010  09:49:24      page 4

(    1) # 1 "sol_oce.F90"
(    1) MODULE sol_oce
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  sol_oce  ***
(    4)    !! Ocean solver :  solver variables defined in memory 
(    5)    !!=====================================================================
(    6)    !!
(    7)    !! ** Purpose :   Define in memory solver variables
(    8)    !!
(    9)    !! History :
(   10)    !!   9.0  !  02-11  (G. Madec)  F90: Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/sol_oce.F90,v 1.8 2006/03/10 10:55:43 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !!-----------------------------------
(   23)    !! elliptic solver: SOR, PCG or FETI
(   24)    !! ----------------------------------
(   25)    INTEGER , PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
(   26)       nsolv    =    1 ,     &  !: = 1/2/3/4 type of elliptic solver
(   27)       nsol_arp =    0 ,     &  !: = 0/1 absolute/relative precision convergence test
(   28)       nmin     =  300 ,     &  !: minimum of iterations for the SOR solver
(   29)       nmax     =  800 ,     &  !: maximum of iterations for the SOR solver
(   30)       nmod     =   10 ,     &  !: frequency of test for the SOR solver
(   31)       nmisl    = 4000          !: maximum pcg iterations for island
(   32)      
(   33)    REAL(wp), PUBLIC ::      & !!: namsol   elliptic solver / island / free surface
(   34)       eps    =  1.e-6_wp ,  &  !: absolute precision of the solver
(   35)       resmax = 1.e-14_wp ,  &  !: absolute precision for the SOR solver
(   36)       sor    =   1.92_wp ,  &  !: optimal coefficient for the SOR solver
(   37)       epsisl = 1.e-10_wp ,  &  !: absolute precision on stream function solver
(   38)       rnu    =    1.0_wp       !: strength of the additional force used in free surface
(   39) 
(   40)    CHARACTER(len=1), PUBLIC ::   &  !:
(   41)       c_solver_pt = 'T'        !: nature of grid-points T (S) for free surface case
(   42)       !                        !                        F (G) for rigid-lid case
(   43) 
(   44)    INTEGER , PUBLIC ::   &  !:
(   45)       ncut,         &  !: indicator of solver convergence
(   46)       niter            !: number of iteration done by the solver
(   47) 
(   48)    REAL(wp), PUBLIC ::   &  !:
(   49)       epsr,         &  !: relative precision for SOR & PCG solvers
(   50)       epsilo,       &  !: precision for the FETI solver
(   51)       rnorme, res,  &  !: intermediate modulus, solver residu
(   52)       alph,         &  !: coefficient  =(gcr,gcr)/(gcx,gccd)
(   53)       beta,         &  !: coefficient  =(rn+1,rn+1)/(rn,rn)
(   54)       radd,         &  !: coefficient  =(gccd,gcdes)
(   55)       rr               !: coefficient  =(rn,rn)
(   56) 
(   57)    REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj,4) ::   &  !:






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 5

(   58)       gcp              !: barotropic matrix extra-diagonal elements
(   59) 
(   60)    REAL(wp), PUBLIC, DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj) ::   &  !:
(   61)       gcx, gcxb,    &  !: now, before solution of the elliptic equation
(   62)       gcdprc,       &  !: inverse diagonal preconditioning matrix
(   63)       gcdmat,       &  !: diagonal preconditioning matrix
(   64)       gcb,          &  !: second member of the barotropic linear system
(   65)       gcr,          &  !: residu =b-a.x
(   66)       gcdes,        &  !: vector descente
(   67)       gccd             !: vector such that ca.gccd=a.d (ca-1=gcdprc)
(   68) 
(   69) 
(   70) 
(   71) # 116
(  116)    !!----------------------------------------------------------------------
(  117) END MODULE sol_oce
