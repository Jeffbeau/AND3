


PGF90 (Version     10.5)          09/23/2010  09:49:09      page 1

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

Filename: trdicp.F90

(    1) MODULE trdicp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdicp  ***
(    4)    !! Ocean diagnostics:  ocean tracers and dynamic trends
(    5)    !!=====================================================================
(    6) #if  defined key_trdtra   ||   defined key_trddyn   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_trdtra'  or                  active tracers trends diagnostics
(    9)    !!   'key_trddyn'                            momentum trends diagnostics
(   10)    !!----------------------------------------------------------------------
(   11) 
(   12)    !!----------------------------------------------------------------------
(   13)    !!   trd              : verify the basin averaged properties for tra/dyn 
(   14)    !!   trd_dwr          : print dynmaic trends in ocean.output file
(   15)    !!   trd_twr          : print tracers trends in ocean.output file
(   16)    !!   trd_icp_init     : initialization step
(   17)    !!----------------------------------------------------------------------
(   18)    !! * Modules used
(   19)    USE oce             ! ocean dynamics and tracers variables
(   20)    USE dom_oce         ! ocean space and time domain variables
(   21)    USE trdmod_oce      ! ocean variables trends
(   22)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   23)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   24)    USE zdf_oce         ! ocean vertical physics
(   25)    USE in_out_manager  ! I/O manager
(   26)    USE lib_mpp         ! distibuted memory computing library
(   27)    USE eosbn2          ! equation of state
(   28)    USE phycst          ! physical constants
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Interfaces
(   34)    INTERFACE trd
(   35)       MODULE PROCEDURE trd_2d, trd_3d
(   36)    END INTERFACE
(   37) 
(   38)    !! * Routine accessibility
(   39)    PUBLIC trd                   ! called by step.F90
(   40)    PUBLIC trd_dwr               ! called by step.F90
(   41)    PUBLIC trd_twr               ! called by step.F90
(   42)    PUBLIC trd_icp_init          ! called by opa.F90
(   43) 
(   44)    !! * Shared module variables
(   45) #if  defined key_trdtra   &&   defined key_trddyn    ||   defined key_esopa
(   46)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .TRUE.    !: tracers  trend flag






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 2

(   47)    LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .TRUE.    !: momentum trend flag
(   48) #elif  defined key_trdtra
(   49)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .TRUE.    !: tracers  trend flag
(   50)    LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .FALSE.   !: momentum trend flag
(   51) #elif  defined key_trddyn
(   52)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .FALSE.   !: tracers  trend flag
(   53)    LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .TRUE.    !: momentum trend flag
(   54) #endif
(   55) 
(   56)    !! * Substitutions
(   57) #  include "domzgr_substitute.h90"
(   58) #  include "vectopt_loop_substitute.h90"
(   59)    !!----------------------------------------------------------------------
(   60)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   61)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdicp.F90,v 1.3 2005/03/30 10:29:20 opalod Exp $ 
(   62)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   63)    !!----------------------------------------------------------------------
(   64) 
(   65) CONTAINS
(   66) 
(   67)    SUBROUTINE trd_2d(ptrd2dx, ptrd2dy, ktrd , ctype)
(   68)       !!---------------------------------------------------------------------
(   69)       !!                  ***  ROUTINE trd_2d  ***
(   70)       !! 
(   71)       !! ** Purpose : verify the basin averaged properties of tracers and/or
(   72)       !!              momentum equations at every time step frequency ntrd.
(   73)       !!
(   74)       !! ** Method :
(   75)       !!
(   76)       !! History :
(   77)       !!        !  91-12 (G. Madec)
(   78)       !!        !  92-06 (M. Imbard) add time step frequency
(   79)       !!        !  96-01 (G. Madec)  terrain following coordinates
(   80)       !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Arguments
(   84)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
(   85)          ptrd2dx,                      &   ! Temperature or U trend 
(   86)          ptrd2dy                           ! Salinity    or V trend
(   87) 
(   88)       INTEGER, INTENT( in ) ::   ktrd      ! tracer trend index
(   89) 
(   90)       CHARACTER(len=3), INTENT( in ) ::   &
(   91)          ctype                             ! momentum or tracers trends type
(   92)          !                                 ! 'DYN' or 'TRA'
(   93) 
(   94)       !! * Local declarations
(   95)       INTEGER ::   ji, jj        ! loop indices
(   96)       REAL(wp) ::   &
(   97)          zbt, zbtu, zbtv,     &  ! temporary scalars
(   98)          zmsku, zmskv            !    "         "
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       ! 1. Advective trends and forcing trend
(  102)       ! -------------------------------------
(  103) 
(  104)       ! 1.1 Mask the forcing trend and substract it from the vertical diffusion trend






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 3

(  105)       SELECT CASE (ctype)
(  106) 
(  107)       CASE ('DYN')              ! Momentum
(  108)          DO jj = 1, jpjm1
(  109)             DO ji = 1, jpim1
(  110)                zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,1)
(  111)                zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,1)
(  112)                ptrd2dx(ji,jj) = ptrd2dx(ji,jj) * zmsku
(  113)                ptrd2dy(ji,jj) = ptrd2dy(ji,jj) * zmskv
(  114)             END DO
(  115)          END DO
(  116)          ptrd2dx(jpi, : ) = 0.e0      ;      ptrd2dy(jpi, : ) = 0.e0
(  117)          ptrd2dx( : ,jpj) = 0.e0      ;      ptrd2dy( : ,jpj) = 0.e0
(  118) 
(  119)       CASE ('TRA')              ! Tracers
(  120)          ptrd2dx(:,:) = ptrd2dx(:,:) * tmask_i(:,:)
(  121)          ptrd2dy(:,:) = ptrd2dy(:,:) * tmask_i(:,:)
(  122) 
(  123)       END SELECT
(  124)       
(  125)       ! 2. Basin averaged tracer trend
(  126)       ! ------------------------------
(  127) 
(  128)       SELECT CASE (ctype)
(  129) 
(  130)       CASE ('DYN')              ! Momentum
(  131)          umo(ktrd) = 0.e0
(  132)          vmo(ktrd) = 0.e0
(  133) 
(  134)          SELECT CASE (ktrd)
(  135) 
(  136)          CASE (jpdtdswf)         ! surface forcing
(  137)             DO jj = 1, jpj
(  138)                DO ji = 1, jpi
(  139)                   umo(ktrd) = umo(ktrd) + ptrd2dx(ji,jj) * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,1)
(  140)                   vmo(ktrd) = vmo(ktrd) + ptrd2dy(ji,jj) * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,1)
(  141)                END DO
(  142)             END DO
(  143) 
(  144)          CASE (jpdtdbfr)         ! bottom friction fluxes
(  145)             DO jj = 1, jpj
(  146)                DO ji = 1, jpi
(  147)                   umo(ktrd) = umo(ktrd) + ptrd2dx(ji,jj)
(  148)                   vmo(ktrd) = vmo(ktrd) + ptrd2dy(ji,jj)
(  149)                END DO
(  150)             END DO
(  151) 
(  152)          END SELECT
(  153) 
(  154)       CASE ('TRA')              ! Tracers
(  155)          tmo(ktrd) = 0.e0
(  156)          smo(ktrd) = 0.e0
(  157)          DO jj = 1, jpj
(  158)             DO ji = 1, jpi
(  159)                zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1)
(  160)                tmo(ktrd) =  tmo(ktrd) + ptrd2dx(ji,jj) * zbt
(  161)                smo(ktrd) =  smo(ktrd) + ptrd2dy(ji,jj) * zbt
(  162)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 4

(  163)          END DO
(  164) 
(  165)       END SELECT
(  166)       
(  167)       ! 3. Basin averaged tracer square trend
(  168)       ! -------------------------------------
(  169)       ! c a u t i o n: field now
(  170)       
(  171)       SELECT CASE (ctype)
(  172) 
(  173)       CASE ('DYN')              ! Momentum
(  174)          hke(ktrd) = 0.e0
(  175)          DO jj = 1, jpj
(  176)             DO ji = 1, jpi
(  177)                zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,1)
(  178)                zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,1)
(  179)                hke(ktrd) = hke(ktrd)   &
(  180)                &   + un(ji,jj,1) * ptrd2dx(ji,jj) * zbtu &
(  181)                &   + vn(ji,jj,1) * ptrd2dy(ji,jj) * zbtv
(  182)             END DO
(  183)          END DO
(  184) 
(  185)       CASE ('TRA')              ! Tracers
(  186)          t2(ktrd) = 0.e0
(  187)          s2(ktrd) = 0.e0
(  188)          DO jj = 1, jpj
(  189)             DO ji = 1, jpi
(  190)                zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,1)
(  191)                t2(ktrd) = t2(ktrd) + ptrd2dx(ji,jj) * zbt * tn(ji,jj,1)
(  192)                s2(ktrd) = s2(ktrd) + ptrd2dy(ji,jj) * zbt * sn(ji,jj,1)
(  193)             END DO
(  194)          END DO
(  195)       
(  196)       END SELECT
(  197) 
(  198)    END SUBROUTINE trd_2d
(  199) 
(  200) 
(  201) 
(  202)    SUBROUTINE trd_3d(ptrd3dx, ptrd3dy, ktrd, ctype)
(  203)       !!---------------------------------------------------------------------
(  204)       !!                  ***  ROUTINE trd_3d  ***
(  205)       !! 
(  206)       !! ** Purpose : verify the basin averaged properties of tracers and/or
(  207)       !!              momentum equations at every time step frequency ntrd.
(  208)       !!
(  209)       !! ** Method :
(  210)       !!
(  211)       !! History :
(  212)       !!        !  91-12 (G. Madec)
(  213)       !!        !  92-06 (M. Imbard) add time step frequency
(  214)       !!        !  96-01 (G. Madec)  terrain following coordinates
(  215)       !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
(  216)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(  217)       !!----------------------------------------------------------------------
(  218)       !! * Arguments
(  219)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(  220)           ptrd3dx,                     &   ! Temperature or U trend 






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 5

(  221)           ptrd3dy                          ! Salinity    or V trend
(  222) 
(  223)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  224) 
(  225)       CHARACTER(len=3), INTENT( in ) ::   &
(  226)          ctype                             ! momentum or tracers trends type
(  227)          !                                 ! 'DYN' or 'TRA'
(  228) 
(  229)       !! * Local declarations
(  230)       INTEGER ::   ji, jj, jk
(  231)       REAL(wp) ::   &
(  232)          zbt, zbtu, zbtv,               &  ! temporary scalars
(  233)          zmsku, zmskv
(  234)       !!----------------------------------------------------------------------
(  235) 
(  236)       ! 1. Advective trends and forcing trend
(  237)       ! -------------------------------------
(  238) 
(  239)       ! Mask the trends
(  240)       SELECT CASE (ctype)
(  241) 
(  242)       CASE ('DYN')              ! Momentum        
(  243)          DO jk = 1, jpk
(  244)             DO jj = 1, jpjm1
(  245)                DO ji = 1, jpim1
(  246)                   zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,jk)
(  247)                   zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)
(  248)                   ptrd3dx(ji,jj,jk) = ptrd3dx(ji,jj,jk) * zmsku
(  249)                   ptrd3dy(ji,jj,jk) = ptrd3dy(ji,jj,jk) * zmskv
(  250)                ENDDO
(  251)             ENDDO
(  252)          ENDDO
(  253) 
(  254)          ptrd3dx(jpi, : ,:) = 0.e0      ;      ptrd3dy(jpi, : ,:) = 0.e0
(  255)          ptrd3dx( : ,jpj,:) = 0.e0      ;      ptrd3dy( : ,jpj,:) = 0.e0
(  256) 
(  257)       CASE ('TRA')              ! Tracers
(  258)          DO jk = 1, jpk
(  259)             ptrd3dx(:,:,jk) = ptrd3dx(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
(  260)             ptrd3dy(:,:,jk) = ptrd3dy(:,:,jk) * tmask(:,:,jk) * tmask_i(:,:)
(  261)          ENDDO
(  262) 
(  263)       END SELECT   
(  264) 
(  265)       ! 2. Basin averaged tracer/momentum trend
(  266)       ! ---------------------------------------
(  267)       
(  268)       SELECT CASE (ctype)
(  269) 
(  270)       CASE ('DYN')              ! Momentum
(  271)          umo(ktrd) = 0.e0
(  272)          vmo(ktrd) = 0.e0
(  273)          DO jk = 1, jpk
(  274)             DO jj = 1, jpj
(  275)                DO ji = 1, jpi
(  276)                   zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
(  277)                   zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
(  278)                   umo(ktrd) = umo(ktrd) + ptrd3dx(ji,jj,jk) * zbtu






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 6

(  279)                   vmo(ktrd) = vmo(ktrd) + ptrd3dy(ji,jj,jk) * zbtv
(  280)                END DO
(  281)             END DO
(  282)          END DO
(  283) 
(  284)       CASE ('TRA')              ! Tracers
(  285)          tmo(ktrd) = 0.e0
(  286)          smo(ktrd) = 0.e0
(  287)          DO jk = 1, jpkm1
(  288)             DO jj = 1, jpj
(  289)                DO ji = 1, jpi
(  290)                   zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) 
(  291)                   tmo(ktrd) =  tmo(ktrd) + ptrd3dx(ji,jj,jk) * zbt
(  292)                   smo(ktrd) =  smo(ktrd) + ptrd3dy(ji,jj,jk) * zbt
(  293)                END DO
(  294)             END DO
(  295)          END DO
(  296) 
(  297)       END SELECT
(  298) 
(  299)       ! 3. Basin averaged tracer/momentum square trend
(  300)       ! ----------------------------------------------
(  301)       ! c a u t i o n: field now
(  302)       
(  303)       SELECT CASE (ctype)
(  304) 
(  305)       CASE ('DYN')              ! Momentum
(  306)          hke(ktrd) = 0.e0
(  307)          DO jk = 1, jpk
(  308)             DO jj = 1, jpj
(  309)                DO ji = 1, jpi
(  310)                   zbtu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
(  311)                   zbtv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
(  312)                   hke(ktrd) = hke(ktrd)   &
(  313)                   &   + un(ji,jj,jk) * ptrd3dx(ji,jj,jk) * zbtu &
(  314)                   &   + vn(ji,jj,jk) * ptrd3dy(ji,jj,jk) * zbtv
(  315)                END DO
(  316)             END DO
(  317)          END DO
(  318) 
(  319)       CASE ('TRA')              ! Tracers
(  320)          t2(ktrd) = 0.e0
(  321)          s2(ktrd) = 0.e0
(  322)          DO jk = 1, jpk
(  323)             DO jj = 1, jpj
(  324)                DO ji = 1, jpi
(  325)                   zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
(  326)                   t2(ktrd) = t2(ktrd) + ptrd3dx(ji,jj,jk) * zbt * tn(ji,jj,jk)
(  327)                   s2(ktrd) = s2(ktrd) + ptrd3dy(ji,jj,jk) * zbt * sn(ji,jj,jk)
(  328)                END DO
(  329)             END DO
(  330)          END DO
(  331) 
(  332)       END SELECT
(  333) 
(  334)    END SUBROUTINE trd_3d
(  335) 
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 7

(  337) 
(  338)    SUBROUTINE trd_icp_init
(  339)       !!---------------------------------------------------------------------
(  340)       !!                  ***  ROUTINE trd_icp_init  ***
(  341)       !! 
(  342)       !! ** Purpose :   
(  343)       !!
(  344)       !! ** Method  :
(  345)       !!
(  346)       !! History :
(  347)       !!   9.0  !  03-09 (G. Madec)  Original code
(  348)       !!        !  04-08 (C. Talandier) New trends organization
(  349)       !!----------------------------------------------------------------------
(  350)       !! * Local declarations
(  351)       INTEGER :: ji, jj, jk
(  352) 
(  353)       REAL(wp) ::   zmskt
(  354) #if  defined key_trddyn
(  355)       REAL(wp) ::   zmsku,zmskv
(  356) #endif
(  357) 
(  358)       NAMELIST/namtrd/ ntrd, nctls
(  359)       !!----------------------------------------------------------------------
(  360) 
(  361)       ! namelist namtrd : trend diagnostic
(  362)       REWIND( numnam )
(  363)       READ  ( numnam, namtrd )
(  364) 
(  365)       IF(lwp) THEN
(  366)          WRITE(numout,*)
(  367)          WRITE(numout,*) 'trd_icp_init : integral constraints properties trends'
(  368)          WRITE(numout,*) '~~~~~~~~~~~~~'
(  369)          WRITE(numout,*) ' '
(  370)          WRITE(numout,*) '          Namelist namtrd : '
(  371)          WRITE(numout,*) '             time step frequency trend       ntrd  = ', ntrd
(  372)       ENDIF
(  373) 
(  374)       ! initialisation of BBL tracers lateral diffusion to zero
(  375)       tldfbbl(:,:) = 0.e0   ;   sldfbbl(:,:) = 0.e0  
(  376)       ! initialisation of BBL tracers lateral advection to zero
(  377)       tladbbl(:,:) = 0.e0   ;   sladbbl(:,:) = 0.e0  
(  378)       ! initialisation of workspace
(  379)       tladi(:,:,:) = 0.e0  ;  tladj(:,:,:) = 0.e0
(  380)       sladi(:,:,:) = 0.e0  ;  sladj(:,:,:) = 0.e0
(  381) 
(  382)       ! Total volume at t-points:
(  383)       tvolt = 0.e0
(  384)       DO jk = 1, jpkm1
(  385)          DO jj = 2, jpjm1
(  386)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  387)                zmskt = tmask(ji,jj,jk) * tmask_i(ji,jj)
(  388)                tvolt = tvolt + zmskt * e1t(ji,jj) *e2t(ji,jj) * fse3t(ji,jj,jk)
(  389)             END DO
(  390)          END DO
(  391)       END DO
(  392)       IF( lk_mpp )   CALL mpp_sum( tvolt )   ! sum over the global domain
(  393) 
(  394)       IF(lwp) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 8

(  395)          WRITE(numout,*)
(  396)          WRITE(numout,*) '          total ocean volume at T-point   tvolt = ',tvolt
(  397)       ENDIF
(  398) 
(  399) #if  defined key_trddyn
(  400)       ! Initialization of potential to kinetic energy conversion
(  401)       rpktrd = 0.e0
(  402) 
(  403)       ! Total volume at u-, v- points:
(  404)       tvolu = 0.e0
(  405)       tvolv = 0.e0
(  406) 
(  407)       DO jk = 1, jpk
(  408)          DO jj = 2, jpjm1
(  409)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  410)                zmsku = tmask_i(ji+1,jj  ) * tmask_i(ji,jj) * umask(ji,jj,jk)
(  411)                zmskv = tmask_i(ji  ,jj+1) * tmask_i(ji,jj) * vmask(ji,jj,jk)
(  412)                tvolu = tvolu + zmsku * e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
(  413)                tvolv = tvolv + zmskv * e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
(  414)             END DO
(  415)          END DO
(  416)       END DO
(  417)       IF( lk_mpp )   CALL mpp_sum( tvolu )   ! sums over the global domain
(  418)       IF( lk_mpp )   CALL mpp_sum( tvolv )
(  419) 
(  420)       IF(lwp) THEN
(  421)          WRITE(numout,*) '          total ocean volume at U-point   tvolu = ',tvolu
(  422)          WRITE(numout,*) '          total ocean volume at V-point   tvolv = ',tvolv
(  423)          WRITE(numout,*) ' '
(  424)       ENDIF
(  425) #endif
(  426) 
(  427)    END SUBROUTINE trd_icp_init
(  428) 
(  429) 
(  430) 
(  431)    SUBROUTINE trd_dwr( kt )
(  432)       !!---------------------------------------------------------------------
(  433)       !!                  ***  ROUTINE trd_dwr  ***
(  434)       !! 
(  435)       !! ** Purpose :  write dynamic trends in ocean.output 
(  436)       !!
(  437)       !! ** Method  :
(  438)       !!
(  439)       !! History :
(  440)       !!   9.0  !  03-09  (G. Madec)  Original code
(  441)       !!        !  04-08  (C. Talandier)  New trends organization
(  442)       !!----------------------------------------------------------------------
(  443)       !! * Arguments
(  444)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  445)       INTEGER ::   ji, jj, jk
(  446)       REAL(wp) ::   &
(  447)          ze1e2w,zcof,        &  !    "         "
(  448)          zbe1ru, zbe2rv,     &  !    "         "
(  449)          zbtr, ztz, zth 
(  450) 
(  451)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  452)                zkepe, zkx, zky, zkz              ! temporary arrays






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 9

(  453)       !!----------------------------------------------------------------------
(  454) 
(  455)       ! I. Momentum trends
(  456)       ! -------------------
(  457) 
(  458)       IF( MOD(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN
(  459) 
(  460)          ! I.1 Conversion potential energy - kinetic energy
(  461)          ! --------------------------------------------------
(  462)          ! c a u t i o n here, trends are computed at kt+1 (now , but after the swap)
(  463) 
(  464)          zkx(:,:,:) = 0.e0
(  465)          zky(:,:,:) = 0.e0
(  466)          zkz(:,:,:) = 0.e0
(  467)          zkepe(:,:,:) = 0.e0
(  468)    
(  469)          CALL eos( tn, sn, rhd, rhop )       ! now potential and in situ densities
(  470) 
(  471)          ! 4.1 Density flux at w-point
(  472)          DO jk = 2, jpk
(  473)             DO jj = 1, jpj
(  474)                DO ji = 1, jpi
(  475)                   ze1e2w = 0.5 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk) * tmask_i(ji,jj)
(  476)                   zkz(ji,jj,jk) = ze1e2w / rau0 * ( rhop(ji,jj,jk) + rhop(ji,jj,jk-1) )
(  477)                END DO
(  478)             END DO
(  479)          END DO
(  480)          zkz  (:,:, 1 ) = 0.e0
(  481)          
(  482)          ! Density flux at u and v-points
(  483)          DO jk = 1, jpk
(  484)             DO jj = 1, jpjm1
(  485)                DO ji = 1, jpim1
(  486)                   zcof   = 0.5 / rau0
(  487)                   zbe1ru = zcof * e2u(ji,jj) * fse3u(ji,jj,jk) * un(ji,jj,jk)
(  488)                   zbe2rv = zcof * e1v(ji,jj) * fse3v(ji,jj,jk) * vn(ji,jj,jk)
(  489)                   zkx(ji,jj,jk) = zbe1ru * ( rhop(ji,jj,jk) + rhop(ji+1,jj,jk) )
(  490)                   zky(ji,jj,jk) = zbe2rv * ( rhop(ji,jj,jk) + rhop(ji,jj+1,jk) )
(  491)                END DO
(  492)             END DO
(  493)          END DO
(  494)          
(  495)          ! Density flux divergence at t-point
(  496)          DO jk = 1, jpkm1
(  497)             DO jj = 2, jpjm1
(  498)                DO ji = 2, jpim1
(  499)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  500)                   ztz = - zbtr * (    zkz(ji,jj,jk) - zkz(ji,jj,jk+1) )
(  501)                   zth = - zbtr * (  ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk) )   &
(  502)                     &             + ( zky(ji,jj,jk) - zky(ji,jj-1,jk) )  )
(  503)                   zkepe(ji,jj,jk) = (zth + ztz) * tmask(ji,jj,jk) * tmask_i(ji,jj)
(  504)                END DO
(  505)             END DO
(  506)          END DO
(  507)          zkepe( : , : ,jpk) = 0.e0
(  508)          zkepe( : ,jpj, : ) = 0.e0
(  509)          zkepe(jpi, : , : ) = 0.e0
(  510) 






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 10

(  511)          ! I.2 Basin averaged kinetic energy trend
(  512)          ! ----------------------------------------
(  513)          peke = 0.e0
(  514)          DO jk = 1,jpk
(  515)             DO jj = 1, jpj
(  516)                DO ji = 1, jpi
(  517)                   peke = peke + zkepe(ji,jj,jk) * grav * fsdept(ji,jj,jk)   &
(  518)                      &                     * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk)
(  519)                END DO
(  520)             END DO
(  521)          END DO
(  522) 
(  523)          ! I.3 Sums over the global domain
(  524)          ! ---------------------------------
(  525)          IF( lk_mpp ) THEN
(  526)                CALL mpp_sum( peke )
(  527)                CALL mpp_sum( umo , 11 )
(  528)                CALL mpp_sum( vmo , 11 )
(  529)                CALL mpp_sum( hke , 10 )
(  530)          END IF
(  531) 
(  532)          ! I.2 Print dynamic trends in the ocean.output file
(  533)          ! --------------------------------------------------
(  534) 
(  535)          IF(lwp) THEN
(  536)             WRITE (numout,*)
(  537)             WRITE (numout,*)
(  538)             WRITE (numout,9500) kt
(  539)             WRITE (numout,9501) umo( 1) / tvolu, vmo( 1) / tvolv
(  540)             WRITE (numout,9502) umo( 2) / tvolu, vmo( 2) / tvolv
(  541)             WRITE (numout,9503) umo( 3) / tvolu, vmo( 3) / tvolv
(  542)             WRITE (numout,9504) umo( 4) / tvolu, vmo( 4) / tvolv
(  543)             WRITE (numout,9505) umo( 5) / tvolu, vmo( 5) / tvolv
(  544)             WRITE (numout,9506) umo( 6) / tvolu, vmo( 6) / tvolv
(  545)             WRITE (numout,9507) umo( 7) / tvolu, vmo( 7) / tvolv
(  546)             WRITE (numout,9508) umo( 8) / tvolu, vmo( 8) / tvolv
(  547)             WRITE (numout,9509) umo(10) / tvolu, vmo(10) / tvolv
(  548)             WRITE (numout,9510) umo( 9) / tvolu, vmo( 9) / tvolv
(  549)             WRITE (numout,9511) umo(11) / tvolu, vmo(11) / tvolv
(  550)             WRITE (numout,9512)
(  551)             WRITE (numout,9513)                                                 &
(  552)             &     (  umo(1) + umo(2) + umo(3) + umo( 4) + umo( 5) + umo(6)    &
(  553)             &      + umo(7) + umo(8) + umo(9) + umo(10) + umo(11) ) / tvolu,   &
(  554)             &     (  vmo(1) + vmo(2) + vmo(3) + vmo( 4) + vmo( 5) + vmo(6)    &
(  555)             &      + vmo(7) + vmo(8) + vmo(9) + vmo(10) + vmo(11) ) / tvolv
(  556)          ENDIF
(  557) 
(  558)  9500    FORMAT(' momentum trend at it= ', i6, ' :', /' ==============================')
(  559)  9501    FORMAT(' pressure gradient          u= ', e20.13, '    v= ', e20.13)
(  560)  9502    FORMAT(' ke gradient                u= ', e20.13, '    v= ', e20.13)
(  561)  9503    FORMAT(' relative vorticity term    u= ', e20.13, '    v= ', e20.13)
(  562)  9504    FORMAT(' coriolis term              u= ', e20.13, '    v= ', e20.13)
(  563)  9505    FORMAT(' horizontal diffusion       u= ', e20.13, '    v= ', e20.13)
(  564)  9506    FORMAT(' vertical advection         u= ', e20.13, '    v= ', e20.13)
(  565)  9507    FORMAT(' vertical diffusion         u= ', e20.13, '    v= ', e20.13)
(  566)  9508    FORMAT(' surface pressure gradient  u= ', e20.13, '    v= ', e20.13)
(  567)  9509    FORMAT(' forcing term               u= ', e20.13, '    v= ', e20.13)
(  568)  9510    FORMAT(' dampimg term               u= ', e20.13, '    v= ', e20.13)






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 11

(  569)  9511    FORMAT(' bottom flux                u= ', e20.13, '    v= ', e20.13)
(  570)  9512    FORMAT(' -----------------------------------------------------------------------------')
(  571)  9513    FORMAT(' total trend                u= ', e20.13, '    v= ', e20.13)
(  572) 
(  573)          IF(lwp) THEN
(  574)             WRITE (numout,*)
(  575)             WRITE (numout,*)
(  576)             WRITE (numout,9520) kt
(  577)             WRITE (numout,9521) hke( 1) / tvolt
(  578)             WRITE (numout,9522) hke( 2) / tvolt
(  579)             WRITE (numout,9523) hke( 3) / tvolt
(  580)             WRITE (numout,9524) hke( 4) / tvolt
(  581)             WRITE (numout,9525) hke( 5) / tvolt
(  582)             WRITE (numout,9526) hke( 6) / tvolt
(  583)             WRITE (numout,9527) hke( 7) / tvolt
(  584)             WRITE (numout,9528) hke( 8) / tvolt
(  585)             WRITE (numout,9529) hke(10) / tvolt
(  586)             WRITE (numout,9530) hke( 9) / tvolt
(  587)             WRITE (numout,9531)
(  588)             WRITE (numout,9532)   &
(  589)             &     (  hke(1) + hke(2) + hke(3) + hke(4) + hke(5) + hke(6)   &
(  590)             &      + hke(7) + hke(8) + hke(9) + hke(10) ) / tvolt
(  591)          ENDIF
(  592) 
(  593)  9520    FORMAT(' kinetic energy trend at it= ', i6, ' :', /' ====================================')
(  594)  9521    FORMAT(' pressure gradient         u2= ', e20.13)
(  595)  9522    FORMAT(' ke gradient               u2= ', e20.13)
(  596)  9523    FORMAT(' relative vorticity term   u2= ', e20.13)
(  597)  9524    FORMAT(' coriolis term             u2= ', e20.13)
(  598)  9525    FORMAT(' horizontal diffusion      u2= ', e20.13)
(  599)  9526    FORMAT(' vertical advection        u2= ', e20.13)
(  600)  9527    FORMAT(' vertical diffusion        u2= ', e20.13)
(  601)  9528    FORMAT(' surface pressure gradient u2= ', e20.13)
(  602)  9529    FORMAT(' forcing term              u2= ', e20.13)
(  603)  9530    FORMAT(' dampimg term              u2= ', e20.13)
(  604)  9531    FORMAT(' --------------------------------------------------')
(  605)  9532    FORMAT(' total trend               u2= ', e20.13)
(  606) 
(  607)          IF(lwp) THEN
(  608)             WRITE (numout,*)
(  609)             WRITE (numout,*)
(  610)             WRITE (numout,9540) kt
(  611)             WRITE (numout,9541) ( hke(2) + hke(3) + hke(6) ) / tvolt
(  612)             WRITE (numout,9542) ( hke(2) + hke(6) ) / tvolt
(  613)             WRITE (numout,9543) ( hke(4) ) / tvolt
(  614)             WRITE (numout,9544) ( hke(3) ) / tvolt
(  615)             WRITE (numout,9545) ( hke(8) ) / tvolt
(  616)             WRITE (numout,9546) ( hke(5) ) / tvolt
(  617)             WRITE (numout,9547) ( hke(7) ) / tvolt
(  618)             WRITE (numout,9548) ( hke(1) ) / tvolt, rpktrd / tvolt
(  619)          ENDIF
(  620) 
(  621)  9540    FORMAT(' energetic consistency at it= ', i6, ' :', /' =========================================')
(  622)  9541    FORMAT(' 0 = non linear term(true if key_vorenergy or key_combined): ', e20.13)
(  623)  9542    FORMAT(' 0 = ke gradient + vertical advection              : ', e20.13)
(  624)  9543    FORMAT(' 0 = coriolis term  (true if key_vorenergy or key_combined): ', e20.13)
(  625)  9544    FORMAT(' 0 = uh.( rot(u) x uh ) (true if enstrophy conser.)    : ', e20.13)
(  626)  9545    FORMAT(' 0 = surface pressure gradient                     : ', e20.13)






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 12

(  627)  9546    FORMAT(' 0 > horizontal diffusion                          : ', e20.13)
(  628)  9547    FORMAT(' 0 > vertical diffusion                            : ', e20.13)
(  629)  9548    FORMAT(' pressure gradient u2 = - 1/rau0 u.dz(rhop)        : ', e20.13, '  u.dz(rhop) =', e20.13)
(  630) 
(  631)          ! Save potential to kinetic energy conversion for next time step
(  632)          rpktrd = peke
(  633) 
(  634)       ENDIF
(  635) 
(  636)    END SUBROUTINE trd_dwr
(  637) 
(  638) 
(  639) 
(  640) 
(  641)    SUBROUTINE trd_twr( kt )
(  642)       !!---------------------------------------------------------------------
(  643)       !!                  ***  ROUTINE trd_twr  ***
(  644)       !! 
(  645)       !! ** Purpose :  write active tracers trends in ocean.output 
(  646)       !!
(  647)       !! ** Method  :
(  648)       !!
(  649)       !! History :
(  650)       !!   9.0  !  03-09  (G. Madec)  Original code
(  651)       !!        !  04-08  (C. Talandier)  New trends organization
(  652)       !!----------------------------------------------------------------------
(  653)       !! * Arguments
(  654)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  655) 
(  656)       !!----------------------------------------------------------------------
(  657) 
(  658)       ! I. Tracers trends
(  659)       ! -----------------
(  660) 
(  661)       IF( MOD(kt,ntrd) == 0 .OR. kt == nit000 .OR. kt == nitend ) THEN
(  662) 
(  663)          ! I.1 Sums over the global domain
(  664)          ! -------------------------------
(  665)          IF( lk_mpp ) THEN
(  666)             CALL mpp_sum( tmo, 10 )   
(  667)             CALL mpp_sum( smo, 10 )
(  668)             CALL mpp_sum( t2 , 10 )
(  669)             CALL mpp_sum( s2 , 10 )
(  670)          ENDIF
(  671) 
(  672)          ! I.2 Print tracers trends in the ocean.output file
(  673)          ! --------------------------------------------------
(  674)          
(  675)          IF(lwp) THEN
(  676)             WRITE (numout,*)
(  677)             WRITE (numout,*)
(  678)             WRITE (numout,9400) kt
(  679)             WRITE (numout,9401) tmo(1) / tvolt, smo(1) / tvolt
(  680)             WRITE (numout,9402) tmo(2) / tvolt, smo(2) / tvolt
(  681)             WRITE (numout,9403) tmo(3) / tvolt, smo(3) / tvolt
(  682)             WRITE (numout,9404) tmo(4) / tvolt, smo(4) / tvolt
(  683)             WRITE (numout,9405) tmo(5) / tvolt, smo(5) / tvolt
(  684)             WRITE (numout,9406) tmo(6) / tvolt, smo(6) / tvolt






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 13

(  685)             WRITE (numout,9407) tmo(7) / tvolt
(  686)             WRITE (numout,9408) tmo(8) / tvolt, smo(8) / tvolt
(  687)             WRITE (numout,9409)
(  688)             WRITE (numout,9410) (  tmo(1) + tmo(2) + tmo(3) + tmo(4)              &
(  689)             &                    + tmo(5) + tmo(6) + tmo(7) + tmo(8) ) / tvolt,   &
(  690)             &                   (  smo(1) + smo(2) + smo(3) + smo(4)              &
(  691)             &                    + smo(5) + smo(6)           + smo(8) ) / tvolt
(  692)          ENDIF
(  693) 
(  694) 9400     FORMAT(' tracer trend at it= ',i6,' :     temperature',   &
(  695)               '              salinity',/' ============================')
(  696) 9401     FORMAT(' horizontal advection        ',e20.13,'     ',e20.13)
(  697) 9402     FORMAT(' vertical advection          ',e20.13,'     ',e20.13)
(  698) 9403     FORMAT(' horizontal diffusion        ',e20.13,'     ',e20.13)
(  699) 9404     FORMAT(' vertical diffusion          ',e20.13,'     ',e20.13)
(  700) 9405     FORMAT(' STATIC instability mixing   ',e20.13,'     ',e20.13)
(  701) 9406     FORMAT(' damping term                ',e20.13,'     ',e20.13)
(  702) 9407     FORMAT(' penetrative qsr             ',e20.13,'     ',e20.13)
(  703) 9408     FORMAT(' forcing term                ',e20.13,'     ',e20.13)
(  704) 9409     FORMAT(' -------------------------------------------------------------------------')
(  705) 9410     FORMAT(' total trend                 ',e20.13,'     ',e20.13)
(  706) 
(  707) 
(  708)          IF(lwp) THEN
(  709)             WRITE (numout,*)
(  710)             WRITE (numout,*)
(  711)             WRITE (numout,9420) kt
(  712)             WRITE (numout,9421) t2(1) / tvolt, s2(1) / tvolt
(  713)             WRITE (numout,9422) t2(2) / tvolt, s2(2) / tvolt
(  714)             WRITE (numout,9423) t2(3) / tvolt, s2(3) / tvolt
(  715)             WRITE (numout,9424) t2(4) / tvolt, s2(4) / tvolt
(  716)             WRITE (numout,9425) t2(5) / tvolt, s2(5) / tvolt
(  717)             WRITE (numout,9426) t2(6) / tvolt, s2(6) / tvolt
(  718)             WRITE (numout,9427) t2(7) / tvolt
(  719)             WRITE (numout,9428) t2(8) / tvolt, s2(8) / tvolt
(  720)             WRITE (numout,9429)
(  721)             WRITE (numout,9430) (  t2(1) + t2(2) + t2(3) + t2(4)              &
(  722)             &                    + t2(5) + t2(6) + t2(7) + t2(8) ) / tvolt,   &
(  723)             &                   (  s2(1) + s2(2) + s2(3) + s2(4)              &
(  724)             &                    + s2(5) + s2(6)          + s2(8) ) / tvolt
(  725)          ENDIF
(  726) 
(  727) 9420     FORMAT(' tracer**2 trend at it= ', i6, ' :      temperature',   &
(  728)             '               salinity', /, ' ===============================')
(  729) 9421     FORMAT(' horizontal advection      * t   ', e20.13, '     ', e20.13)
(  730) 9422     FORMAT(' vertical advection        * t   ', e20.13, '     ', e20.13)
(  731) 9423     FORMAT(' horizontal diffusion      * t   ', e20.13, '     ', e20.13)
(  732) 9424     FORMAT(' vertical diffusion        * t   ', e20.13, '     ', e20.13)
(  733) 9425     FORMAT(' STATIC instability mixing * t   ', e20.13, '     ', e20.13)
(  734) 9426     FORMAT(' damping term              * t   ', e20.13, '     ', e20.13)
(  735) 9427     FORMAT(' penetrative qsr           * t   ', e20.13, '     ', e20.13)
(  736) 9428     FORMAT(' forcing term              * t   ', e20.13, '     ', e20.13)
(  737) 9429     FORMAT(' -----------------------------------------------------------------------------')
(  738) 9430     FORMAT(' total trend                *t = ', e20.13, '  *s = ', e20.13)
(  739) 
(  740) 
(  741)          IF(lwp) THEN
(  742)             WRITE (numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 14

(  743)             WRITE (numout,*)
(  744)             WRITE (numout,9440) kt
(  745)             WRITE (numout,9441) ( tmo(1)+tmo(2) )/tvolt, ( smo(1)+smo(2) )/tvolt
(  746)             WRITE (numout,9442)   tmo(3)/tvolt,  smo(3)/tvolt
(  747)             WRITE (numout,9443)   tmo(4)/tvolt,  smo(4)/tvolt
(  748)             WRITE (numout,9444)   tmo(5)/tvolt,  smo(5)/tvolt
(  749)             WRITE (numout,9445) ( t2(1)+t2(2) )/tvolt, ( s2(1)+s2(2) )/tvolt
(  750)             WRITE (numout,9446)   t2(3)/tvolt,   s2(3)/tvolt
(  751)             WRITE (numout,9447)   t2(4)/tvolt,   s2(4)/tvolt
(  752)             WRITE (numout,9448)   t2(5)/tvolt,   s2(5)/tvolt
(  753)          ENDIF
(  754) 
(  755) 9440     FORMAT(' tracer consistency at it= ',i6,   &
(  756)             ' :         temperature','                salinity',/,   &
(  757)             ' ==================================')
(  758) 9441     FORMAT(' 0 = horizontal+vertical advection      ',e20.13,'       ',e20.13)
(  759) 9442     FORMAT(' 0 = horizontal diffusion               ',e20.13,'       ',e20.13)
(  760) 9443     FORMAT(' 0 = vertical diffusion                 ',e20.13,'       ',e20.13)
(  761) 9444     FORMAT(' 0 = static instability mixing          ',e20.13,'       ',e20.13)
(  762) 9445     FORMAT(' 0 = horizontal+vertical advection * t  ',e20.13,'       ',e20.13)
(  763) 9446     FORMAT(' 0 > horizontal diffusion          * t  ',e20.13,'       ',e20.13)
(  764) 9447     FORMAT(' 0 > vertical diffusion            * t  ',e20.13,'       ',e20.13)
(  765) 9448     FORMAT(' 0 > static instability mixing     * t  ',e20.13,'       ',e20.13)
(  766) 
(  767)       ENDIF
(  768) 
(  769)    END SUBROUTINE trd_twr
(  770) 
(  771) #   else
(  772)    !!----------------------------------------------------------------------
(  773)    !!   Default case :                                         Empty module
(  774)    !!----------------------------------------------------------------------
(  775)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .FALSE.   !: tracers  trend flag
(  776)    LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .FALSE.   !: momentum trend flag
(  777) CONTAINS
(  778)    SUBROUTINE trd_2d(ptrd2dx, ptrd2dy, ktrd , ctype)       ! Empty routine
(  779)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  780)           ptrd2dx,                     &   ! Temperature or U trend 
(  781)           ptrd2dy                          ! Salinity    or V trend
(  782)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  783)       CHARACTER(len=3), INTENT( in ) ::   &
(  784)          ctype                             ! momentum or tracers trends type
(  785) !      WRITE(*,*) 'trd_2d: You should not have seen this print! error ?', ptrd2dx(1,1,1)
(  786) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd2dy(1,1,1)
(  787) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  788) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
(  789)    END SUBROUTINE trd_2d
(  790)    SUBROUTINE trd_3d(ptrd3dx, ptrd3dy, ktrd , ctype)       ! Empty routine
(  791)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  792)           ptrd3dx,                     &   ! Temperature or U trend 
(  793)           ptrd3dy                          ! Salinity    or V trend
(  794)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  795)       CHARACTER(len=3), INTENT( in ) ::   &
(  796)          ctype                             ! momentum or tracers trends type
(  797) !      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
(  798) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
(  799) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  800) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 15

(  801)    END SUBROUTINE trd_3d
(  802)    SUBROUTINE trd_icp_init               ! Empty routine
(  803)    END SUBROUTINE trd_icp_init
(  804)    SUBROUTINE trd_dwr( kt )          ! Empty routine
(  805)       INTEGER, INTENT(in) :: kt
(  806) !      WRITE(*,*) 'trd_dwr: You should not have seen this print! error ?', kt
(  807)    END SUBROUTINE trd_dwr
(  808)    SUBROUTINE trd_twr( kt )          ! Empty routine
(  809)       INTEGER, INTENT(in) :: kt
(  810) !      WRITE(*,*) 'trd_twr: You should not have seen this print! error ?', kt
(  811)    END SUBROUTINE trd_twr
(  812) #endif
(  813) 
(  814)    !!======================================================================
(  815) END MODULE trdicp

















































PGF90 (Version     10.5)          09/23/2010  09:49:09      page 16

(    1) # 1 "trdicp.F90"
(    1) MODULE trdicp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdicp  ***
(    4)    !! Ocean diagnostics:  ocean tracers and dynamic trends
(    5)    !!=====================================================================
(    6) # 772
(  772)    !!----------------------------------------------------------------------
(  773)    !!   Default case :                                         Empty module
(  774)    !!----------------------------------------------------------------------
(  775)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdtra = .FALSE.   !: tracers  trend flag
(  776)    LOGICAL, PUBLIC, PARAMETER ::   lk_trddyn = .FALSE.   !: momentum trend flag
(  777) CONTAINS
(  778)    SUBROUTINE trd_2d(ptrd2dx, ptrd2dy, ktrd , ctype)       ! Empty routine
(  779)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  780)           ptrd2dx,                     &   ! Temperature or U trend 
(  781)           ptrd2dy                          ! Salinity    or V trend
(  782)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  783)       CHARACTER(len=3), INTENT( in ) ::   &
(  784)          ctype                             ! momentum or tracers trends type
(  785) !      WRITE(*,*) 'trd_2d: You should not have seen this print! error ?', ptrd2dx(1,1,1)
(  786) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd2dy(1,1,1)
(  787) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  788) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
(  789)    END SUBROUTINE trd_2d
(  790)    SUBROUTINE trd_3d(ptrd3dx, ptrd3dy, ktrd , ctype)       ! Empty routine
(  791)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  792)           ptrd3dx,                     &   ! Temperature or U trend 
(  793)           ptrd3dy                          ! Salinity    or V trend
(  794)       INTEGER, INTENT( in ) ::   ktrd      ! momentum or tracer trend index
(  795)       CHARACTER(len=3), INTENT( in ) ::   &
(  796)          ctype                             ! momentum or tracers trends type
(  797) !      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1)
(  798) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ptrd3dy(1,1,1)
(  799) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd
(  800) !      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ctype
(  801)    END SUBROUTINE trd_3d
(  802)    SUBROUTINE trd_icp_init               ! Empty routine
(  803)    END SUBROUTINE trd_icp_init
(  804)    SUBROUTINE trd_dwr( kt )          ! Empty routine
(  805)       INTEGER, INTENT(in) :: kt
(  806) !      WRITE(*,*) 'trd_dwr: You should not have seen this print! error ?', kt
(  807)    END SUBROUTINE trd_dwr
(  808)    SUBROUTINE trd_twr( kt )          ! Empty routine
(  809)       INTEGER, INTENT(in) :: kt
(  810) !      WRITE(*,*) 'trd_twr: You should not have seen this print! error ?', kt
(  811)    END SUBROUTINE trd_twr
(  812) 
(  813) # 814
(  814)    !!======================================================================
(  815) END MODULE trdicp
