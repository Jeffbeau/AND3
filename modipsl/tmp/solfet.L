


PGF90 (Version     10.5)          09/23/2010  09:49:25      page 1

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

Filename: solfet.F90

(    1) MODULE solfet
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solfet
(    4)    !! Ocean solver :  Finite Elements Tearing & Interconnecting solver
(    5)    !!=====================================================================
(    6) #if defined key_feti
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_feti' :                                            FETI solver
(    9)    !!----------------------------------------------------------------------
(   10)    !!   sol_fet     : FETI solver
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers variables
(   14)    USE dom_oce         ! ocean space and time domain variables 
(   15)    USE sol_oce         ! ocean solver 
(   16)    USE lib_mpp         ! distribued memory computing
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Routine accessibility
(   23)    PUBLIC sol_fet              ! ???
(   24)    !!----------------------------------------------------------------------
(   25)    !!----------------------------------------------------------------------
(   26)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solfet.F90,v 1.4 2005/12/21 10:46:41 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) CONTAINS
(   31) 
(   32)    SUBROUTINE sol_fet( kindic )
(   33)       !!---------------------------------------------------------------------
(   34)       !!                  ***  ROUTINE sol_fet  ***
(   35)       !!
(   36)       !! ** Purpose :
(   37)       !!     Solve the ellipic equation for the barotropic stream function
(   38)       !!     system (default option) or the transport divergence system
(   39)       !!     (lk_dynspg_flt=T) using a Finite Elements Tearing and 
(   40)       !!      Interconnecting (FETI) approach.
(   41)       !!     In the former case, the barotropic stream function trend has a
(   42)       !!     zero boundary condition along all coastlines (i.e. continent
(   43)       !!     as well as islands) while in the latter the boundary condition
(   44)       !!     specification is not required.
(   45)       !!
(   46)       !! ** Method :






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 2

(   47)       !!      Resolution of the elliptic equation by a Dual formulation of
(   48)       !!      the Schur Complement Method or Finite Elements Tearing & 
(   49)       !!      Interconnecting (FETI) approach
(   50)       !!
(   51)       !! ** Action :
(   52)       !!
(   53)       !! ** References :
(   54)       !!      Guyon, M, Roux, F-X, Chartier, M and Fraunie, P, 1994 :
(   55)       !!      A domain decomposition solver to compute the barotropic 
(   56)       !!      component of an OGCM in the parallel processing field.
(   57)       !!      Ocean Modelling, issue 105, december 94.
(   58)       !!
(   59)       !! History :
(   60)       !!        !  97-02  (M. Guyon)  original code 
(   61)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   62)       !!----------------------------------------------------------------------
(   63)       !! * Arguments
(   64)       INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
(   65)       !                                      ! gence is not reached: the model is
(   66)       !                                      ! stopped in step
(   67) 
(   68)       !! * Local declarations
(   69)       INTEGER ::   ji, jj                    ! dummy loop indices
(   70)       INTEGER ::   iit0, inn, iju 
(   71)       REAL(wp) ::   zgcad, zgwgt
(   72) 
(   73) 
(   74)       ! Norme(gcr)   = (gcb, gcb)
(   75)       ! gcdes = gsr
(   76)       
(   77)       ! bmask field is filtering the differents contribution on the 
(   78)       ! non-overlapping interface
(   79)       
(   80)       gcb(:,:) = bmask(:,:) * gcb(:,:)
(   81)       
(   82)       ! Mpp: sum over all the global domain
(   83)       
(   84)       ! copy the right hand member
(   85)       
(   86)       CALL feti_vmov(noeuds,gcb(1,1),wfeti(may))
(   87)       
(   88)       ! conservation of descent direction if ntest = 0
(   89)       
(   90)       IF(ntest /= 0) mjj0=0
(   91)       iit0 = mjj0
(   92)       !                      --->
(   93)       !    resolution of the Grad(PS) equation by a Dual formulation of the
(   94)       !    Schur Complement Method or Finite Elements Tearing & Interconnecting
(   95)       !    (FETI) approach
(   96)       !    interface problem (Lagrange  multiplier) : PCPG algorithm
(   97)       !    local problem (trend of the 2D potential field) : LU factorization
(   98)       !    preconditioner : lumped
(   99)       !    optimisation : Krylov initialisation + Krylov correction
(  100) 
(  101)       CALL feti_dualschur(noeuds,nifmat+1,njfmat+1,wfeti(maan),   &
(  102)           npblo,wfeti(mablo),ninterf,   &
(  103)           ninterfc,nni,nnic,mfet(mandvois),mfet(mandvoisc),   &
(  104)           mfet(maplistin),mfet(malistin),   &






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 3

(  105)           wfeti(mapoids),wfeti(miax),   &
(  106)           wfeti(maz),wfeti(may),   &
(  107)           wfeti(mabitw),wfeti(mautilu),   &
(  108)           wfeti(malambda),wfeti(mag),   &
(  109)           wfeti(mapg),wfeti(mamg),nitmax,nmaxd,mjj0,   &
(  110)           wfeti(mawj),   &
(  111)           wfeti(madwj),wfeti(madwwj),wfeti(magamm),   &
(  112)           wfeti(mawork),   &
(  113)           wfeti(mabufin),wfeti(mabufout),narea,epsilo,   &
(  114)           ndlblo,mfet(malisblo),ndkerep,   &
(  115)           wfeti(maxnul),wfeti(maynul),numit0ete,nitmaxete,   &
(  116)           wfeti(maeteg),wfeti(maeteag),wfeti(maeted),   &
(  117)           wfeti(maetead),wfeti(maeteadd),wfeti(maetegamm),   &
(  118)           wfeti(mansp),   &
(  119)           wfeti(maetev),wfeti(maetew),nnih,mfet(maplistih),   &
(  120)           wfeti(magh),   &
(  121)           wfeti(maw),wfeti(madw),   &
(  122)           res,kindic,inn)
(  123) 
(  124)       ! number of iteration of the pcg to solve the interface pb
(  125) 
(  126)       inn =  mjj0 - iit0
(  127) 
(  128)       ! test of convergence
(  129)       IF( res < epsilo .OR. inn == nmax ) THEN
(  130)           niter = inn
(  131)           ncut  = 999
(  132)       ENDIF
(  133) 
(  134)       ! indicator of non-convergence or explosion
(  135)       IF( inn == nmax .OR. rr > 1.e+20 ) kindic = -2
(  136)       IF( ncut == 999 ) GOTO 999
(  137) 
(  138) 999   CONTINUE
(  139) 
(  140)       !  2. Output in gcx
(  141)       !  -----------------
(  142)       
(  143)       CALL feti_vmov( noeuds, wfeti(miax), gcx )
(  144) 
(  145)       CALL lbc_lnk( gcx, c_solver_pt, 1. )   ! lateral boundary condition
(  146) 
(  147)    END SUBROUTINE sol_fet
(  148) 
(  149) #else
(  150)    !!----------------------------------------------------------------------
(  151)    !!   Default option :                                       Empty module
(  152)    !!----------------------------------------------------------------------
(  153) CONTAINS
(  154)    SUBROUTINE sol_fet( kindic )        ! Empty routine
(  155)       INTEGER, INTENT( inout ) ::   kindic  ! solver problem 
(  156)       kindic = -100
(  157)    END SUBROUTINE sol_fet
(  158) #endif
(  159) 
(  160)    !!=====================================================================
(  161) END MODULE solfet







PGF90 (Version     10.5)          09/23/2010  09:49:25      page 4

(    1) # 1 "solfet.F90"
(    1) MODULE solfet
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solfet
(    4)    !! Ocean solver :  Finite Elements Tearing & Interconnecting solver
(    5)    !!=====================================================================
(    6) # 150
(  150)    !!----------------------------------------------------------------------
(  151)    !!   Default option :                                       Empty module
(  152)    !!----------------------------------------------------------------------
(  153) CONTAINS
(  154)    SUBROUTINE sol_fet( kindic )        ! Empty routine
(  155)       INTEGER, INTENT( inout ) ::   kindic  ! solver problem 
(  156)       kindic = -100
(  157)    END SUBROUTINE sol_fet
(  158) 
(  159) # 160
(  160)    !!=====================================================================
(  161) END MODULE solfet
