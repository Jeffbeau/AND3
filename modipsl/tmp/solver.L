


PGF90 (Version     10.5)          09/23/2010  09:50:21      page 1

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

Filename: solver.F90

(    1) MODULE solver
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solver  ***
(    4)    !! Ocean solver :  initialization of ocean solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   solver_init: solver initialization
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE sol_oce         ! solver variables
(   15)    USE solmat          ! ???
(   16)    USE solisl          ! ???
(   17)    USE obc_oce         ! Lateral open boundary condition
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE lib_mpp
(   21)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   22) 
(   23)    IMPLICIT NONE
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solver.F90,v 1.13 2006/03/20 17:27:15 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE solver_init( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE solver_init  ***
(   36)       !!                   
(   37)       !! ** Purpose :   Initialization for the solver of the elliptic equation:
(   38)       !!       * default option: barotropic stream function system
(   39)       !!         and islands initialization (if lk_isl=T)
(   40)       !!       * lk_dynspg_flt = T : transport divergence system. No specific
(   41)       !!         treatment of islands.
(   42)       !!      
(   43)       !! ** Method :
(   44)       !!       - Compute the local depth of the water column at u- and v-point
(   45)       !!      (lk_dynspg_flt = T) or its inverse (lk_dynspg_rl = T).
(   46)       !!      The local depth of the water column is computed by summing 






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 2

(   47)       !!      the vertical scale factors. For its inverse, the thickness of
(   48)       !!      the first model level is imposed as lower bound. The inverse of
(   49)       !!      this depth is THEN taken and masked, so that the inverse of the
(   50)       !!      local depth is zero when the local depth is zero.
(   51)       !!       - Construct the matrix of the elliptic system by a call to
(   52)       !!      solmat.F routine.
(   53)       !!       - island (if lk_isl=T)
(   54)       !!            isl_dom: find islands from the bathymetry file
(   55)       !!            isl_bsf: compute the island barotropic stream function
(   56)       !!            isl_mat: compute the inverse island matrix
(   57)       !!            set mbathy to the number of non-zero w-levels of a water
(   58)       !!            column (the minimum value of mbathy is 2):
(   59)       !!                  mbathy = min( mbathy, 1 ) + 1
(   60)       !!
(   61)       !! ** Action : - hur, hvr : masked inverse of the local depth at
(   62)       !!                                u- and v-point. (lk_dynspg_rl = T)
(   63)       !!             - hu, hv   : masked local depth at u- and v- points
(   64)       !!                                (lk_dynspg_flt = T)
(   65)       !!             - c_solver_pt : nature of the gridpoint at which the
(   66)       !!                                solver is applied
(   67)       !! References :
(   68)       !!      Jensen, 1986: adv. phys. oceanogr. num. mod.,ed. o brien,87-110.
(   69)       !!      Madec & Marti, 1990: internal rep. LODYC, 90/03., 29pp.
(   70)       !!
(   71)       !! History :
(   72)       !!        !  90-10  (G. Madec)  Original code           
(   73)       !!        !  93-02  (O. Marti)                         
(   74)       !!        !  97-02  (G. Madec)  local depth inverse computation
(   75)       !!        !  98-10  (G. Roullet, G. Madec)  free surface 
(   76)       !!   9.0  !  03-07  (G. Madec)  free form, F90
(   77)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   78)       !!----------------------------------------------------------------------
(   79)       !! * Arguments
(   80)       INTEGER, INTENT(in) :: kt
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER :: ji, jj   ! dummy loop indices
(   84)       CHARACTER(len=80) :: clname
(   85) 
(   86)       NAMELIST/namsol/ nsolv, nsol_arp, nmin, nmax, nmod, eps, resmax, sor, epsisl, nmisl, rnu
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89)       IF(lwp) WRITE(numout,*)
(   90)       IF(lwp) WRITE(numout,*) 'solver_init : solver to compute the surface pressure gradient'
(   91)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(   92) 
(   93)       ! open elliptic solver statistics file
(   94)       clname = 'solver.stat'
(   95)       CALL ctlopn( numsol, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
(   96)                    1, numout, lwp, 1 )
(   97) 
(   98) 
(   99)       ! 0. Define the solver parameters
(  100)       !    ----------------------------
(  101)       ! Namelist namsol : elliptic solver / islands / free surface
(  102)       REWIND( numnam )
(  103)       READ  ( numnam, namsol )
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 3

(  105) #if defined key_feti
(  106)       ! FETI algorithm, we force nsolv at 3
(  107)       nsolv = 3
(  108) #endif
(  109) 
(  110) 
(  111)       ! 0. Parameter control and print
(  112)       !    ---------------------------
(  113) 
(  114)       ! Control print
(  115)       IF(lwp) WRITE(numout,*) '          Namelist namsol : set solver parameters'
(  116) 
(  117)       IF(lwp) THEN
(  118)          WRITE(numout,*) '             type of elliptic solver            nsolv    = ', nsolv
(  119)          WRITE(numout,*) '             absolute/relative (0/1) precision  nsol_arp = ', nsol_arp
(  120)          WRITE(numout,*) '             minimum iterations for solver      nmin     = ', nmin
(  121)          WRITE(numout,*) '             maximum iterations for solver      nmax     = ', nmax
(  122)          WRITE(numout,*) '             frequency for test                 nmod     = ', nmod
(  123)          WRITE(numout,*) '             absolute precision of solver       eps      = ', eps
(  124)          WRITE(numout,*) '             absolute precision for SOR solver  resmax   = ', resmax
(  125)          WRITE(numout,*) '             optimal coefficient of sor         sor      = ', sor
(  126)          IF(lk_isl) WRITE(numout,*) '             absolute precision stream fct    epsisl   = ', epsisl
(  127)          IF(lk_isl) WRITE(numout,*) '             maximum pcg iterations island    nmisl    = ', nmisl
(  128)          WRITE(numout,*) '             free surface parameter         rnu    = ', rnu
(  129)          WRITE(numout,*)
(  130)       ENDIF
(  131) 
(  132)       IF( lk_dynspg_flt ) THEN
(  133)          IF(lwp) WRITE(numout,*)
(  134)          IF(lwp) WRITE(numout,*) '          free surface formulation'
(  135)          IF( lk_isl ) THEN
(  136)             IF(lwp) WRITE(numout,cform_err)
(  137)             IF(lwp) WRITE(numout,*) ' key_islands inconsistent with key_dynspg_flt'
(  138)             nstop = nstop + 1
(  139)          ENDIF
(  140)       ELSEIF( lk_dynspg_rl ) THEN
(  141)          IF(lwp) WRITE(numout,*)
(  142)          IF(lwp) WRITE(numout,*) '          Rigid lid formulation'
(  143)       ELSE
(  144)          IF(lwp) WRITE(numout,cform_err)
(  145)          IF(lwp) WRITE(numout,*) '          Choose only one surface pressure gradient calculation: filtered or rigid-lid'
(  146)          IF(lwp) WRITE(numout,*) '          Should not call this routine if dynspg_exp or dynspg_ts has been chosen'
(  147)          nstop = nstop + 1
(  148)       ENDIF
(  149)       IF( lk_dynspg_flt .AND. lk_dynspg_rl ) THEN
(  150)          IF(lwp) WRITE(numout,cform_err)
(  151)          IF(lwp) WRITE(numout,*) '          Chose between free surface or rigid-lid, not both'
(  152)          nstop = nstop + 1
(  153)       ENDIF
(  154) 
(  155)       SELECT CASE ( nsolv )
(  156) 
(  157)       CASE ( 1 )                ! preconditioned conjugate gradient solver
(  158)          IF(lwp) WRITE(numout,*) '          a preconditioned conjugate gradient solver is used'
(  159)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  160)             IF(lwp) WRITE(numout,cform_err)
(  161)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  162)             nstop = nstop + 1






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 4

(  163)          ENDIF
(  164) 
(  165)       CASE ( 2 )                ! successive-over-relaxation solver
(  166)          IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver is used'
(  167)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  168)             IF(lwp) WRITE(numout,cform_err)
(  169)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  170)             nstop = nstop + 1
(  171)          ENDIF
(  172) 
(  173)       CASE ( 3 )                ! FETI solver
(  174)          IF(lwp) WRITE(numout,*) '          the FETI solver is used'
(  175)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  176)             IF(lwp) WRITE(numout,cform_err)
(  177)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  178)             nstop = nstop + 1
(  179)          ENDIF
(  180)          IF( .NOT.lk_mpp ) THEN
(  181)             IF(lwp) WRITE(numout,cform_err)
(  182)             IF(lwp) WRITE(numout,*) ' The FETI algorithm is used only with the key_mpp_... option'
(  183)             nstop = nstop + 1
(  184)          ELSE
(  185)             IF( jpnij == 1 ) THEN
(  186)                IF(lwp) WRITE(numout,cform_err)
(  187)                IF(lwp) WRITE(numout,*) ' The FETI algorithm needs more than one processor'
(  188)                nstop = nstop + 1
(  189)             ENDIF
(  190)          ENDIF
(  191)          
(  192)       CASE ( 4 )                ! successive-over-relaxation solver with extra outer halo
(  193)          IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver with extra outer halo is used'
(  194)          IF(lwp) WRITE(numout,*) '          with jpr2di =', jpr2di, ' and  jpr2dj =', jpr2dj
(  195)          IF( .NOT. lk_mpp .AND. jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  196)             IF(lwp) WRITE(numout,cform_err)
(  197)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj are not equal to zero'
(  198)             IF(lwp) WRITE(numout,*) ' In this case this algorithm should be used only with the key_mpp_... option'
(  199)             nstop = nstop + 1
(  200)          ELSE
(  201)             IF( ( ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) .OR. ( jpni /= 1 ) ) &
(  202)               &  .AND. ( jpr2di /= jpr2dj ) ) THEN  
(  203)                IF(lwp) WRITE(numout,cform_err)
(  204)                IF(lwp) WRITE(numout,*) '          jpr2di should be equal to jpr2dj'
(  205)                nstop = nstop + 1
(  206)             ENDIF
(  207)          ENDIF
(  208) 
(  209)       CASE DEFAULT
(  210)          IF(lwp) WRITE(numout,cform_err)
(  211)          IF(lwp) WRITE(numout,*) '          bad flag value for nsolv = ', nsolv
(  212)          nstop = nstop + 1
(  213)          
(  214)       END SELECT
(  215) 
(  216)       ! Grid-point at which the solver is applied
(  217)       ! -----------------------------------------
(  218) 
(  219)       IF( lk_dynspg_rl ) THEN       ! rigid-lid
(  220)          IF( lk_mpp ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 5

(  221)             c_solver_pt = 'G'   ! G= F with special staff ??? which one?
(  222)          ELSE
(  223)             c_solver_pt = 'F'
(  224)          ENDIF
(  225)       ELSE                          ! free surface T-point
(  226)          IF( lk_mpp ) THEN
(  227)             c_solver_pt = 'S'   ! S=T with special staff ??? which one?
(  228)          ELSE
(  229)             c_solver_pt = 'T'
(  230)          ENDIF
(  231)       ENDIF
(  232) 
(  233) 
(  234)       ! Construction of the elliptic system matrix
(  235)       ! ------------------------------------------
(  236) 
(  237)       CALL sol_mat( kt )
(  238) 
(  239) 
(  240)       IF( lk_isl ) THEN
(  241)       
(  242)          ! Islands in the domain
(  243)          ! ---------------------
(  244) 
(  245)          IF ( jpisl == 0 ) THEN
(  246)              IF(lwp)WRITE(numout,cform_err)
(  247)              IF(lwp)WRITE(numout,*) ' bad islands parameter jpisl =', jpisl
(  248)              nstop = nstop + 1
(  249)          ENDIF
(  250) 
(  251)          ! open Island streamfunction statistic file
(  252)          CALL ctlopn( numisp, 'islands.stat', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
(  253)             &         1     , numout        , lwp      , 1                         )
(  254)    
(  255)          CALL isl_dom       ! Island identification
(  256) 
(  257)          CALL isl_bsf       ! Island barotropic stream function
(  258) 
(  259)          CALL isl_mat       ! Comput and invert the island matrix
(  260) 
(  261)          ! mbathy set to the number of w-level (minimum value 2)
(  262)          DO jj = 1, jpj
(  263)             DO ji = 1, jpi
(  264)                mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
(  265)             END DO
(  266)          END DO
(  267) 
(  268)       ENDIF
(  269) 
(  270)    END SUBROUTINE solver_init
(  271) 
(  272)    !!======================================================================
(  273) END MODULE solver











PGF90 (Version     10.5)          09/23/2010  09:50:21      page 6

(    1) # 1 "solver.F90"
(    1) MODULE solver
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  solver  ***
(    4)    !! Ocean solver :  initialization of ocean solver
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   solver_init: solver initialization
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE sol_oce         ! solver variables
(   15)    USE solmat          ! ???
(   16)    USE solisl          ! ???
(   17)    USE obc_oce         ! Lateral open boundary condition
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE lib_mpp
(   21)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   22) 
(   23)    IMPLICIT NONE
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solver.F90,v 1.13 2006/03/20 17:27:15 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE solver_init( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE solver_init  ***
(   36)       !!                   
(   37)       !! ** Purpose :   Initialization for the solver of the elliptic equation:
(   38)       !!       * default option: barotropic stream function system
(   39)       !!         and islands initialization (if lk_isl=T)
(   40)       !!       * lk_dynspg_flt = T : transport divergence system. No specific
(   41)       !!         treatment of islands.
(   42)       !!      
(   43)       !! ** Method :
(   44)       !!       - Compute the local depth of the water column at u- and v-point
(   45)       !!      (lk_dynspg_flt = T) or its inverse (lk_dynspg_rl = T).
(   46)       !!      The local depth of the water column is computed by summing 
(   47)       !!      the vertical scale factors. For its inverse, the thickness of
(   48)       !!      the first model level is imposed as lower bound. The inverse of
(   49)       !!      this depth is THEN taken and masked, so that the inverse of the
(   50)       !!      local depth is zero when the local depth is zero.
(   51)       !!       - Construct the matrix of the elliptic system by a call to
(   52)       !!      solmat.F routine.
(   53)       !!       - island (if lk_isl=T)
(   54)       !!            isl_dom: find islands from the bathymetry file
(   55)       !!            isl_bsf: compute the island barotropic stream function
(   56)       !!            isl_mat: compute the inverse island matrix
(   57)       !!            set mbathy to the number of non-zero w-levels of a water






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 7

(   58)       !!            column (the minimum value of mbathy is 2):
(   59)       !!                  mbathy = min( mbathy, 1 ) + 1
(   60)       !!
(   61)       !! ** Action : - hur, hvr : masked inverse of the local depth at
(   62)       !!                                u- and v-point. (lk_dynspg_rl = T)
(   63)       !!             - hu, hv   : masked local depth at u- and v- points
(   64)       !!                                (lk_dynspg_flt = T)
(   65)       !!             - c_solver_pt : nature of the gridpoint at which the
(   66)       !!                                solver is applied
(   67)       !! References :
(   68)       !!      Jensen, 1986: adv. phys. oceanogr. num. mod.,ed. o brien,87-110.
(   69)       !!      Madec & Marti, 1990: internal rep. LODYC, 90/03., 29pp.
(   70)       !!
(   71)       !! History :
(   72)       !!        !  90-10  (G. Madec)  Original code           
(   73)       !!        !  93-02  (O. Marti)                         
(   74)       !!        !  97-02  (G. Madec)  local depth inverse computation
(   75)       !!        !  98-10  (G. Roullet, G. Madec)  free surface 
(   76)       !!   9.0  !  03-07  (G. Madec)  free form, F90
(   77)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   78)       !!----------------------------------------------------------------------
(   79)       !! * Arguments
(   80)       INTEGER, INTENT(in) :: kt
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER :: ji, jj   ! dummy loop indices
(   84)       CHARACTER(len=80) :: clname
(   85) 
(   86)       NAMELIST/namsol/ nsolv, nsol_arp, nmin, nmax, nmod, eps, resmax, sor, epsisl, nmisl, rnu
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89)       IF(lwp) WRITE(numout,*)
(   90)       IF(lwp) WRITE(numout,*) 'solver_init : solver to compute the surface pressure gradient'
(   91)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(   92) 
(   93)       ! open elliptic solver statistics file
(   94)       clname = 'solver.stat'
(   95)       CALL ctlopn( numsol, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
(   96)                    1, numout, lwp, 1 )
(   97) 
(   98) 
(   99)       ! 0. Define the solver parameters
(  100)       !    ----------------------------
(  101)       ! Namelist namsol : elliptic solver / islands / free surface
(  102)       REWIND( numnam )
(  103)       READ  ( numnam, namsol )
(  104) 
(  105) 
(  106) 
(  107) # 111
(  111)       ! 0. Parameter control and print
(  112)       !    ---------------------------
(  113) 
(  114)       ! Control print
(  115)       IF(lwp) WRITE(numout,*) '          Namelist namsol : set solver parameters'
(  116) 
(  117)       IF(lwp) THEN
(  118)          WRITE(numout,*) '             type of elliptic solver            nsolv    = ', nsolv






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 8

(  119)          WRITE(numout,*) '             absolute/relative (0/1) precision  nsol_arp = ', nsol_arp
(  120)          WRITE(numout,*) '             minimum iterations for solver      nmin     = ', nmin
(  121)          WRITE(numout,*) '             maximum iterations for solver      nmax     = ', nmax
(  122)          WRITE(numout,*) '             frequency for test                 nmod     = ', nmod
(  123)          WRITE(numout,*) '             absolute precision of solver       eps      = ', eps
(  124)          WRITE(numout,*) '             absolute precision for SOR solver  resmax   = ', resmax
(  125)          WRITE(numout,*) '             optimal coefficient of sor         sor      = ', sor
(  126)          IF(lk_isl) WRITE(numout,*) '             absolute precision stream fct    epsisl   = ', epsisl
(  127)          IF(lk_isl) WRITE(numout,*) '             maximum pcg iterations island    nmisl    = ', nmisl
(  128)          WRITE(numout,*) '             free surface parameter         rnu    = ', rnu
(  129)          WRITE(numout,*)
(  130)       ENDIF
(  131) 
(  132)       IF( lk_dynspg_flt ) THEN
(  133)          IF(lwp) WRITE(numout,*)
(  134)          IF(lwp) WRITE(numout,*) '          free surface formulation'
(  135)          IF( lk_isl ) THEN
(  136)             IF(lwp) WRITE(numout,cform_err)
(  137)             IF(lwp) WRITE(numout,*) ' key_islands inconsistent with key_dynspg_flt'
(  138)             nstop = nstop + 1
(  139)          ENDIF
(  140)       ELSEIF( lk_dynspg_rl ) THEN
(  141)          IF(lwp) WRITE(numout,*)
(  142)          IF(lwp) WRITE(numout,*) '          Rigid lid formulation'
(  143)       ELSE
(  144)          IF(lwp) WRITE(numout,cform_err)
(  145)          IF(lwp) WRITE(numout,*) '          Choose only one surface pressure gradient calculation: filtered or rigid-lid'
(  146)          IF(lwp) WRITE(numout,*) '          Should not call this routine if dynspg_exp or dynspg_ts has been chosen'
(  147)          nstop = nstop + 1
(  148)       ENDIF
(  149)       IF( lk_dynspg_flt .AND. lk_dynspg_rl ) THEN
(  150)          IF(lwp) WRITE(numout,cform_err)
(  151)          IF(lwp) WRITE(numout,*) '          Chose between free surface or rigid-lid, not both'
(  152)          nstop = nstop + 1
(  153)       ENDIF
(  154) 
(  155)       SELECT CASE ( nsolv )
(  156) 
(  157)       CASE ( 1 )                ! preconditioned conjugate gradient solver
(  158)          IF(lwp) WRITE(numout,*) '          a preconditioned conjugate gradient solver is used'
(  159)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  160)             IF(lwp) WRITE(numout,cform_err)
(  161)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  162)             nstop = nstop + 1
(  163)          ENDIF
(  164) 
(  165)       CASE ( 2 )                ! successive-over-relaxation solver
(  166)          IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver is used'
(  167)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  168)             IF(lwp) WRITE(numout,cform_err)
(  169)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  170)             nstop = nstop + 1
(  171)          ENDIF
(  172) 
(  173)       CASE ( 3 )                ! FETI solver
(  174)          IF(lwp) WRITE(numout,*) '          the FETI solver is used'
(  175)          IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  176)             IF(lwp) WRITE(numout,cform_err)






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 9

(  177)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
(  178)             nstop = nstop + 1
(  179)          ENDIF
(  180)          IF( .NOT.lk_mpp ) THEN
(  181)             IF(lwp) WRITE(numout,cform_err)
(  182)             IF(lwp) WRITE(numout,*) ' The FETI algorithm is used only with the key_mpp_... option'
(  183)             nstop = nstop + 1
(  184)          ELSE
(  185)             IF( jpnij == 1 ) THEN
(  186)                IF(lwp) WRITE(numout,cform_err)
(  187)                IF(lwp) WRITE(numout,*) ' The FETI algorithm needs more than one processor'
(  188)                nstop = nstop + 1
(  189)             ENDIF
(  190)          ENDIF
(  191)          
(  192)       CASE ( 4 )                ! successive-over-relaxation solver with extra outer halo
(  193)          IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver with extra outer halo is used'
(  194)          IF(lwp) WRITE(numout,*) '          with jpr2di =', jpr2di, ' and  jpr2dj =', jpr2dj
(  195)          IF( .NOT. lk_mpp .AND. jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
(  196)             IF(lwp) WRITE(numout,cform_err)
(  197)             IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj are not equal to zero'
(  198)             IF(lwp) WRITE(numout,*) ' In this case this algorithm should be used only with the key_mpp_... option'
(  199)             nstop = nstop + 1
(  200)          ELSE
(  201)             IF( ( ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) .OR. ( jpni /= 1 ) ) &
(  202)               &  .AND. ( jpr2di /= jpr2dj ) ) THEN  
(  203)                IF(lwp) WRITE(numout,cform_err)
(  204)                IF(lwp) WRITE(numout,*) '          jpr2di should be equal to jpr2dj'
(  205)                nstop = nstop + 1
(  206)             ENDIF
(  207)          ENDIF
(  208) 
(  209)       CASE DEFAULT
(  210)          IF(lwp) WRITE(numout,cform_err)
(  211)          IF(lwp) WRITE(numout,*) '          bad flag value for nsolv = ', nsolv
(  212)          nstop = nstop + 1
(  213)          
(  214)       END SELECT
(  215) 
(  216)       ! Grid-point at which the solver is applied
(  217)       ! -----------------------------------------
(  218) 
(  219)       IF( lk_dynspg_rl ) THEN       ! rigid-lid
(  220)          IF( lk_mpp ) THEN
(  221)             c_solver_pt = 'G'   ! G= F with special staff ??? which one?
(  222)          ELSE
(  223)             c_solver_pt = 'F'
(  224)          ENDIF
(  225)       ELSE                          ! free surface T-point
(  226)          IF( lk_mpp ) THEN
(  227)             c_solver_pt = 'S'   ! S=T with special staff ??? which one?
(  228)          ELSE
(  229)             c_solver_pt = 'T'
(  230)          ENDIF
(  231)       ENDIF
(  232) 
(  233) 
(  234)       ! Construction of the elliptic system matrix






PGF90 (Version     10.5)          09/23/2010  09:50:21      page 10

(  235)       ! ------------------------------------------
(  236) 
(  237)       CALL sol_mat( kt )
(  238) 
(  239) 
(  240)       IF( lk_isl ) THEN
(  241)       
(  242)          ! Islands in the domain
(  243)          ! ---------------------
(  244) 
(  245)          IF ( jpisl == 0 ) THEN
(  246)              IF(lwp)WRITE(numout,cform_err)
(  247)              IF(lwp)WRITE(numout,*) ' bad islands parameter jpisl =', jpisl
(  248)              nstop = nstop + 1
(  249)          ENDIF
(  250) 
(  251)          ! open Island streamfunction statistic file
(  252)          CALL ctlopn( numisp, 'islands.stat', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
(  253)             &         1     , numout        , lwp      , 1                         )
(  254)    
(  255)          CALL isl_dom       ! Island identification
(  256) 
(  257)          CALL isl_bsf       ! Island barotropic stream function
(  258) 
(  259)          CALL isl_mat       ! Comput and invert the island matrix
(  260) 
(  261)          ! mbathy set to the number of w-level (minimum value 2)
(  262)          DO jj = 1, jpj
(  263)             DO ji = 1, jpi
(  264)                mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
(  265)             END DO
(  266)          END DO
(  267) 
(  268)       ENDIF
(  269) 
(  270)    END SUBROUTINE solver_init
(  271) 
(  272)    !!======================================================================
(  273) END MODULE solver
