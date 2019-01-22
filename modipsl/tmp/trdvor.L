


PGF90 (Version     10.5)          09/23/2010  09:49:10      page 1

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

Filename: trdvor.F90

(    1) MODULE trdvor
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdvor  ***
(    4)    !! Ocean diagnostics:  momentum trends
(    5)    !!=====================================================================
(    6)    
(    7) #if defined key_trdvor   ||   defined key_esopa
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_trdvor'   : momentum trend diagnostics
(   10)    !!----------------------------------------------------------------------
(   11)    !!   trd_vor      : momentum trends averaged over the depth
(   12)    !!   trd_vor_zint : vorticity vertical integration
(   13)    !!   trd_vor_init : initialization step
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers variables
(   17)    USE dom_oce         ! ocean space and time domain variables
(   18)    USE trdmod_oce      ! ocean variables trends
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE phycst          ! Define parameters for the routines
(   22)    USE ldfdyn_oce      ! ocean active tracers: lateral physics
(   23)    USE daymod          ! calandar
(   24)    USE dianam          ! build the name of file (routine)
(   25)    USE ldfslp          ! iso-neutral slopes
(   26)    USE zdfmxl          ! mixed layer depth
(   27)    USE ioipsl          ! NetCDF library
(   28)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   29) 
(   30) 
(   31)    IMPLICIT NONE
(   32)    PRIVATE
(   33) 
(   34)    !! * Interfaces
(   35)    INTERFACE trd_vor_zint
(   36)       MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
(   37)    END INTERFACE
(   38) 
(   39)    !! * Accessibility
(   40)    PUBLIC trd_vor        ! routine called by step.F90
(   41)    PUBLIC trd_vor_zint   ! routine called by dynamics routines
(   42)    PUBLIC trd_vor_init   ! routine called by opa.F90
(   43) 
(   44)    !! * Shared module variables
(   45)    LOGICAL, PUBLIC ::   lk_trdvor = .TRUE.   ! momentum trend flag
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 2

(   47)    !! * Module variables
(   48)    INTEGER ::                &
(   49)       nh_t, nmoydpvor  ,     &
(   50)       nidvor, nhoridvor,     &
(   51)       ndexvor1(jpi*jpj),     &
(   52)       ndimvor1, icount,      &
(   53)       idebug                    ! (0/1) set it to 1 in case of problem to have more print
(   54) 
(   55)    REAL(wp), DIMENSION(jpi,jpj) ::  &
(   56)      vor_avr    ,     &  ! average
(   57)      vor_avrb   ,     &  ! before vorticity (kt-1)
(   58)      vor_avrbb  ,     &  ! vorticity at begining of the nwrite-1 timestep averaging period
(   59)      vor_avrbn  ,     &  ! after vorticity at time step after the
(   60)      rotot      ,     &  ! begining of the NWRITE-1 timesteps
(   61)      vor_avrtot ,     &
(   62)      vor_avrres
(   63) 
(   64)    REAL(wp), DIMENSION(jpi,jpj,jplvor)::   &  !: curl of trends
(   65)       vortrd   
(   66) 
(   67)    CHARACTER(len=12) ::   cvort
(   68) 
(   69)    !! * Substitutions
(   70) #  include "domzgr_substitute.h90"
(   71) #  include "ldfdyn_substitute.h90"
(   72) #  include "vectopt_loop_substitute.h90"
(   73) 
(   74)    !!----------------------------------------------------------------------
(   75)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   76)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdvor.F90,v 1.5 2005/12/12 14:18:08 opalod Exp $ 
(   77)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   78)    !!----------------------------------------------------------------------
(   79)   
(   80) CONTAINS
(   81) 
(   82)    SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
(   83)       !!----------------------------------------------------------------------------
(   84)       !!                  ***  ROUTINE trd_vor_zint  ***
(   85)       !!
(   86)       !! ** Purpose :   computation of vertically integrated vorticity budgets
(   87)       !!      from ocean surface down to control surface (NetCDF output)
(   88)       !!
(   89)       !! ** Method/usage :
(   90)       !!      integration done over nwrite-1 time steps
(   91)       !!
(   92)       !!
(   93)       !! ** Action :
(   94)       !!            /comvor/   :
(   95)       !!                         vor_avr          average
(   96)       !!                         vor_avrb         vorticity at kt-1
(   97)       !!                         vor_avrbb        vorticity at begining of the NWRITE-1
(   98)       !!                                          time steps averaging period
(   99)       !!                         vor_avrbn         vorticity at time step after the
(  100)       !!                                          begining of the NWRITE-1 time
(  101)       !!                                          steps averaging period
(  102)       !!
(  103)       !!                 trends :
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 3

(  105)       !!                  vortrd (,,1) = Pressure Gradient Trend
(  106)       !!                  vortrd (,,2) = KE Gradient Trend
(  107)       !!                  vortrd (,,3) = Relative Vorticity Trend
(  108)       !!                  vortrd (,,4) = Coriolis Term Trend
(  109)       !!                  vortrd (,,5) = Horizontal Diffusion Trend
(  110)       !!                  vortrd (,,6) = Vertical Advection Trend
(  111)       !!                  vortrd (,,7) = Vertical Diffusion Trend
(  112)       !!                  vortrd (,,8) = Surface Pressure Grad. Trend
(  113)       !!                  vortrd (,,9) = Beta V
(  114)       !!                  vortrd (,,10) = forcing term
(  115)       !!		  vortrd (,,11) = bottom friction term
(  116)       !!                  rotot(,) : total cumulative trends over nwrite-1 time steps
(  117)       !!                  vor_avrtot(,) : first membre of vrticity equation
(  118)       !!                  vor_avrres(,) : residual = dh/dt entrainment
(  119)       !!
(  120)       !!      trends output in netCDF format using ioipsl
(  121)       !!
(  122)       !! History :
(  123)       !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
(  124)       !!        !  04-08  (C. Talandier) New trends organization
(  125)       !!----------------------------------------------------------------------
(  126)       !! * Arguments
(  127)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  128) 
(  129)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
(  130)          putrdvor,                         &  ! u vorticity trend 
(  131)          pvtrdvor                             ! v vorticity trend
(  132) 
(  133)       !! * Local declarations
(  134)       INTEGER ::   ji, jj
(  135)       INTEGER ::   ikbu, ikbum1, ikbv, ikbvm1
(  136)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  137)          zudpvor,                       &  ! total cmulative trends
(  138)          zvdpvor                           !   "      "        "
(  139)       !!----------------------------------------------------------------------
(  140) 
(  141)       ! Initialization
(  142)       zudpvor(:,:) = 0.e0
(  143)       zvdpvor(:,:) = 0.e0
(  144) 
(  145)       CALL lbc_lnk( putrdvor,  'U' , -1. )
(  146)       CALL lbc_lnk( pvtrdvor,  'V' , -1. )
(  147) 
(  148)       !  =====================================
(  149)       !  I vertical integration of 2D trends
(  150)       !  =====================================
(  151) 
(  152)       SELECT CASE (ktrd) 
(  153) 
(  154)       CASE (jpvorbfr)        ! bottom friction
(  155) 
(  156)          DO jj = 2, jpjm1
(  157)             DO ji = fs_2, fs_jpim1 
(  158)                ikbu   = min( mbathy(ji+1,jj), mbathy(ji,jj) )
(  159)                ikbum1 = max( ikbu-1, 1 )
(  160)                ikbv   = min( mbathy(ji,jj+1), mbathy(ji,jj) )
(  161)                ikbvm1 = max( ikbv-1, 1 )
(  162)             






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 4

(  163)                zudpvor(ji,jj) = putrdvor(ji,jj) * fse3u(ji,jj,ikbum1) * e1u(ji,jj) * umask(ji,jj,ikbum1)
(  164)                zvdpvor(ji,jj) = pvtrdvor(ji,jj) * fse3v(ji,jj,ikbvm1) * e2v(ji,jj) * vmask(ji,jj,ikbvm1)
(  165)             END DO
(  166)          END DO
(  167) 
(  168)       CASE (jpvorswf)        ! wind stress
(  169) 
(  170)          zudpvor(:,:) = putrdvor(:,:) * fse3u(:,:,1) * e1u(:,:) * umask(:,:,1)
(  171)          zvdpvor(:,:) = pvtrdvor(:,:) * fse3v(:,:,1) * e2v(:,:) * vmask(:,:,1)
(  172) 
(  173)       END SELECT
(  174) 
(  175)       ! Average except for Beta.V
(  176)       zudpvor(:,:) = zudpvor(:,:) * hur(:,:)
(  177)       zvdpvor(:,:) = zvdpvor(:,:) * hvr(:,:)
(  178)    
(  179)       ! Curl
(  180)       DO ji=1,jpim1
(  181)          DO jj=1,jpjm1
(  182)             vortrd(ji,jj,ktrd) = ( zvdpvor(ji+1,jj) - zvdpvor(ji,jj)        &
(  183)                  &                - ( zudpvor(ji,jj+1) - zudpvor(ji,jj) ) ) &
(  184)                  &               / ( e1f(ji,jj) * e2f(ji,jj) )
(  185)          END DO
(  186)       END DO
(  187) 
(  188)       ! Surface mask
(  189)       vortrd(:,:,ktrd) = vortrd(:,:,ktrd) * fmask(:,:,1)
(  190) 
(  191)       IF( idebug /= 0 ) THEN
(  192)          IF(lwp) WRITE(numout,*) ' debuging trd_vor_zint: I done'
(  193)          CALL FLUSH(numout)
(  194)       ENDIF
(  195) 
(  196)    END SUBROUTINE trd_vor_zint_2d
(  197) 
(  198) 
(  199) 
(  200)    SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
(  201)       !!----------------------------------------------------------------------------
(  202)       !!                  ***  ROUTINE trd_vor_zint  ***
(  203)       !!
(  204)       !! ** Purpose :   computation of vertically integrated vorticity budgets
(  205)       !!      from ocean surface down to control surface (NetCDF output)
(  206)       !!
(  207)       !! ** Method/usage :
(  208)       !!      integration done over nwrite-1 time steps
(  209)       !!
(  210)       !!
(  211)       !! ** Action :
(  212)       !!            /comvor/   :
(  213)       !!                         vor_avr          average
(  214)       !!                         vor_avrb         vorticity at kt-1
(  215)       !!                         vor_avrbb        vorticity at begining of the NWRITE-1
(  216)       !!                                          time steps averaging period
(  217)       !!                         vor_avrbn         vorticity at time step after the
(  218)       !!                                          begining of the NWRITE-1 time
(  219)       !!                                          steps averaging period
(  220)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 5

(  221)       !!                 trends :
(  222)       !!
(  223)       !!                  vortrd (,,1) = Pressure Gradient Trend
(  224)       !!                  vortrd (,,2) = KE Gradient Trend
(  225)       !!                  vortrd (,,3) = Relative Vorticity Trend
(  226)       !!                  vortrd (,,4) = Coriolis Term Trend
(  227)       !!                  vortrd (,,5) = Horizontal Diffusion Trend
(  228)       !!                  vortrd (,,6) = Vertical Advection Trend
(  229)       !!                  vortrd (,,7) = Vertical Diffusion Trend
(  230)       !!                  vortrd (,,8) = Surface Pressure Grad. Trend
(  231)       !!                  vortrd (,,9) = Beta V
(  232)       !!                  vortrd (,,10) = forcing term
(  233)       !!		  vortrd (,,11) = bottom friction term
(  234)       !!                  rotot(,) : total cumulative trends over nwrite-1 time steps
(  235)       !!                  vor_avrtot(,) : first membre of vrticity equation
(  236)       !!                  vor_avrres(,) : residual = dh/dt entrainment
(  237)       !!
(  238)       !!      trends output in netCDF format using ioipsl
(  239)       !!
(  240)       !! History :
(  241)       !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
(  242)       !!        !  04-08  (C. Talandier) New trends organization
(  243)       !!----------------------------------------------------------------------
(  244)       !! * Arguments
(  245)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  246) 
(  247)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(  248)          putrdvor,                         &  ! u vorticity trend 
(  249)          pvtrdvor                             ! v vorticity trend
(  250) 
(  251)       !! * Local declarations
(  252)       INTEGER ::   ji, jj, jk
(  253) 
(  254)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  255)          zubet,                         &  ! u Beta.V case
(  256)          zvbet,                         &  ! v Beta.V case
(  257)          zudpvor,                       &  ! total cmulative trends
(  258)          zvdpvor                           !   "      "        "
(  259)       !!----------------------------------------------------------------------
(  260)      
(  261)       ! Initialization
(  262)       zubet(:,:) = 0.e0
(  263)       zvbet(:,:) = 0.e0
(  264)       zudpvor(:,:) = 0.e0
(  265)       zvdpvor(:,:) = 0.e0
(  266) 
(  267)       !  =====================================
(  268)       !  I vertical integration of 3D trends
(  269)       !  =====================================
(  270) 
(  271)       CALL lbc_lnk( putrdvor, 'U' , -1. )
(  272)       CALL lbc_lnk( pvtrdvor, 'V' , -1. )
(  273) 
(  274)       ! putrdvor and pvtrdvor terms
(  275)       DO jk = 1,jpk
(  276)         zudpvor(:,:) = zudpvor(:,:) + putrdvor(:,:,jk) * fse3u(:,:,jk) * e1u(:,:) * umask(:,:,jk)
(  277)         zvdpvor(:,:) = zvdpvor(:,:) + pvtrdvor(:,:,jk) * fse3v(:,:,jk) * e2v(:,:) * vmask(:,:,jk)
(  278)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 6

(  279) 
(  280)       ! Save Beta.V term to avoid average before Curl
(  281)       ! Beta.V : intergration, no average
(  282)       IF( ktrd == jpvorbev ) THEN 
(  283)          zubet(:,:) = zudpvor(:,:)
(  284)          zvbet(:,:) = zvdpvor(:,:)
(  285)       ENDIF
(  286) 
(  287)       ! Average except for Beta.V
(  288)       zudpvor(:,:) = zudpvor(:,:) * hur(:,:)
(  289)       zvdpvor(:,:) = zvdpvor(:,:) * hvr(:,:)
(  290)    
(  291)       ! Curl
(  292)       DO ji=1,jpim1
(  293)          DO jj=1,jpjm1
(  294)             vortrd(ji,jj,ktrd) = (  zvdpvor(ji+1,jj) - zvdpvor(ji,jj) -   &
(  295)                  &                ( zudpvor(ji,jj+1) - zudpvor(ji,jj) ) ) &
(  296)                  &               / ( e1f(ji,jj) * e2f(ji,jj) )
(  297)          END DO
(  298)       END DO
(  299) 
(  300)       ! Surface mask
(  301)       vortrd(:,:,ktrd) = vortrd(:,:,ktrd) * fmask(:,:,1)
(  302) 
(  303)       ! Special treatement for the Beta.V term
(  304)       ! Compute the Curl of the Beta.V term which is not averaged
(  305)       IF( ktrd == jpvorbev ) THEN
(  306)          DO ji=1,jpim1
(  307)             DO jj=1,jpjm1
(  308)                vortrd(ji,jj,jpvorbev) = (  zvbet(ji+1,jj) - zvbet(ji,jj) -   &
(  309)                     &                    ( zubet(ji,jj+1) - zubet(ji,jj) ) ) &
(  310)                     &                   / ( e1f(ji,jj) * e2f(ji,jj) )
(  311)             END DO
(  312)          END DO
(  313) 
(  314)          ! Average on the Curl
(  315)          vortrd(:,:,jpvorbev) = vortrd(:,:,jpvorbev) * hur(:,:)
(  316) 
(  317)          ! Surface mask
(  318)          vortrd(:,:,jpvorbev) = vortrd(:,:,jpvorbev) * fmask(:,:,1)
(  319)       ENDIF
(  320)    
(  321)       IF( idebug /= 0 ) THEN
(  322)          IF(lwp) WRITE(numout,*) ' debuging trd_vor_zint: I done'
(  323)          CALL FLUSH(numout)
(  324)       ENDIF
(  325) 
(  326)    END SUBROUTINE trd_vor_zint_3d
(  327) 
(  328) 
(  329) 
(  330)    SUBROUTINE trd_vor( kt )
(  331)       !!----------------------------------------------------------------------
(  332)       !!                  ***  ROUTINE trd_vor  ***
(  333)       !! 
(  334)       !! ** Purpose :  computation of cumulated trends over analysis period
(  335)       !!               and make outputs (NetCDF or DIMG format)
(  336)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 7

(  337)       !! ** Method/usage :
(  338)       !!
(  339)       !! History :
(  340)       !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
(  341)       !!        !  04-08  (C. Talandier) New trends organization
(  342)       !!----------------------------------------------------------------------
(  343)       !! * Arguments
(  344)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(  345) 
(  346)       !! * Local declarations
(  347)       INTEGER :: ji, jj, jk, jl, it
(  348) 
(  349)       REAL(wp) :: zmean
(  350) 
(  351)       REAL(wp) ,DIMENSION(jpi,jpj) ::   &
(  352)          zun, zvn
(  353)       !!----------------------------------------------------------------------
(  354) 
(  355)       !  =================
(  356)       !  I. Initialization
(  357)       !  =================
(  358)      
(  359)      
(  360)       ! I.1 set before values of vertically average u and v
(  361)       ! ---------------------------------------------------
(  362) 
(  363)       IF( kt > nit000 ) THEN
(  364)          vor_avrb(:,:) = vor_avr(:,:)
(  365)       ENDIF
(  366) 
(  367)        IF( idebug /= 0 ) THEN
(  368)           if(lwp)WRITE(numout,*) ' debuging trd_vor: I.1 done '
(  369)           CALL FLUSH(numout)
(  370)       ENDIF
(  371) 
(  372)       ! I.2 vertically integrated vorticity
(  373)       !  ----------------------------------
(  374) 
(  375)       vor_avr(:,:) = 0.
(  376)       zun(:,:)=0
(  377)       zvn(:,:)=0
(  378)       vor_avrtot(:,:)=0
(  379)       vor_avrres(:,:)=0
(  380)       
(  381)       ! Vertically averaged velocity
(  382)       DO jk = 1, jpk - 1
(  383)          zun(:,:)=zun(:,:) + e1u(:,:)*un(:,:,jk)*fse3u(:,:,jk)
(  384)          zvn(:,:)=zvn(:,:) + e2v(:,:)*vn(:,:,jk)*fse3v(:,:,jk)
(  385)       END DO
(  386)  
(  387)       zun(:,:)=zun(:,:)*hur(:,:)
(  388)       zvn(:,:)=zvn(:,:)*hvr(:,:)
(  389) 
(  390)       ! Curl
(  391)       DO ji=1,jpim1
(  392)          DO jj=1,jpjm1
(  393)             vor_avr(ji,jj) = ((zvn(ji+1,jj)-zvn(ji,jj))-   &
(  394)                               (zun(ji,jj+1)-zun(ji,jj)))   &






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 8

(  395)                              /( e1f(ji,jj) * e2f(ji,jj) )
(  396)             vor_avr(ji,jj) = vor_avr(ji,jj)*fmask(ji,jj,1)
(  397)          END DO
(  398)       END DO
(  399)       
(  400)       IF(idebug /= 0) THEN
(  401)          if(lwp) WRITE(numout,*) ' debuging trd_vor: I.2 done'
(  402)          CALL FLUSH(numout)
(  403)       ENDIF
(  404) 
(  405)       !  =================================
(  406)       !   II. Cumulated trends
(  407)       !  =================================
(  408) 
(  409)       ! II.1 set `before' mixed layer values for kt = nit000+1
(  410)       ! ------------------------------------------------------
(  411)       IF( kt == nit000+1 ) THEN
(  412)          vor_avrbb(:,:) = vor_avrb(:,:)
(  413)          vor_avrbn(:,:) = vor_avr (:,:)
(  414)       ENDIF
(  415) 
(  416)       IF( idebug /= 0 ) THEN
(  417)          if(lwp)  WRITE(numout,*) ' debuging trd_vor: I1.1 done'
(  418)          CALL FLUSH(numout)
(  419)       ENDIF
(  420) 
(  421)       ! II.2 cumulated trends over analysis period (kt=2 to nwrite)
(  422)       ! ----------------------
(  423)       ! trends cumulated over nwrite-2 time steps
(  424) 
(  425)       IF( kt >= nit000+2 ) THEN
(  426)          nmoydpvor = nmoydpvor + 1
(  427)          DO jl = 1, jplvor
(  428)             IF( jl /= 9 ) THEN
(  429)                rotot(:,:) = rotot(:,:) + vortrd(:,:,jl)
(  430)             ENDIF
(  431)          END DO
(  432)       ENDIF
(  433) 
(  434)       IF( idebug /= 0 ) THEN
(  435)          if(lwp)  WRITE(numout,*) ' debuging trd_vor: II.2 done'
(  436)          CALL FLUSH(numout)
(  437)       ENDIF
(  438) 
(  439)       !  =============================================
(  440)       !   III. Output in netCDF + residual computation
(  441)       !  =============================================
(  442) 
(  443)       IF( MOD( kt - nit000+1, ntrd ) == 0 ) THEN
(  444) 
(  445)          ! III.1 compute total trend
(  446)          ! ------------------------
(  447)          zmean = float(nmoydpvor)
(  448) 
(  449)          vor_avrtot(:,:) = ( vor_avr(:,:) - vor_avrbn(:,:) + vor_avrb(:,:) - &
(  450)                              vor_avrbb(:,:) ) /  (zmean * 2. * rdt)
(  451) 
(  452)          IF( idebug /= 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 9

(  453)              if(lwp)  WRITE(numout,*) ' zmean = ',zmean
(  454)              if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.1 done'
(  455)              CALL FLUSH(numout)
(  456)          ENDIF
(  457) 
(  458)          ! III.2 compute residual
(  459)          ! ---------------------
(  460)          vor_avrres(:,:) = vor_avrtot(:,:) - rotot(:,:) / zmean
(  461) 
(  462)          ! Boundary conditions
(  463)          CALL lbc_lnk( vor_avrtot, 'F', 1. )
(  464)          CALL lbc_lnk( vor_avrres, 'F', 1. )
(  465) 
(  466)          IF( idebug /= 0 ) THEN
(  467)             if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.2 done'
(  468)             CALL FLUSH(numout)
(  469)          ENDIF
(  470) 
(  471)          ! III.3 time evolution array swap
(  472)          ! ------------------------------
(  473)          vor_avrbb(:,:) = vor_avrb(:,:)
(  474)          vor_avrbn(:,:) = vor_avr(:,:)
(  475) 
(  476)          IF( idebug /= 0 ) THEN
(  477)             if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.3 done'
(  478)             CALL FLUSH(numout)
(  479)          ENDIF
(  480) 
(  481)          nmoydpvor=0
(  482) 
(  483)       ENDIF
(  484) 
(  485)       ! III.4 write trends to output
(  486)       ! ---------------------------
(  487) 
(  488)       IF( kt >=  nit000+1 ) THEN
(  489) 
(  490)          ! define time axis
(  491)          it= kt-nit000+1
(  492)          IF( lwp .AND. MOD( kt, ntrd ) == 0 ) THEN
(  493)             WRITE(numout,*) '     trdvor_ncwrite : write NetCDF fields'
(  494)          ENDIF
(  495)  
(  496)          CALL histwrite( nidvor,"sovortPh",it,vortrd(:,:,1),ndimvor1,ndexvor1)  ! grad Ph
(  497)          CALL histwrite( nidvor,"sovortEk",it,vortrd(:,:,2),ndimvor1,ndexvor1)  ! Energy
(  498)          CALL histwrite( nidvor,"sovozeta",it,vortrd(:,:,3),ndimvor1,ndexvor1)  ! rel vorticity
(  499)          CALL histwrite( nidvor,"sovortif",it,vortrd(:,:,4),ndimvor1,ndexvor1)  ! coriolis
(  500)          CALL histwrite( nidvor,"sovodifl",it,vortrd(:,:,5),ndimvor1,ndexvor1)  ! lat diff
(  501)          CALL histwrite( nidvor,"sovoadvv",it,vortrd(:,:,6),ndimvor1,ndexvor1)  ! vert adv
(  502)          CALL histwrite( nidvor,"sovodifv",it,vortrd(:,:,7),ndimvor1,ndexvor1)  ! vert diff
(  503)          CALL histwrite( nidvor,"sovortPs",it,vortrd(:,:,8),ndimvor1,ndexvor1)  ! grad Ps
(  504)          CALL histwrite( nidvor,"sovortbv",it,vortrd(:,:,9),ndimvor1,ndexvor1)  ! beta.V
(  505)          CALL histwrite( nidvor,"sovowind",it,vortrd(:,:,10),ndimvor1,ndexvor1) ! wind stress
(  506)          CALL histwrite( nidvor,"sovobfri",it,vortrd(:,:,11),ndimvor1,ndexvor1) ! bottom friction
(  507)          CALL histwrite( nidvor,"1st_mbre",it,vor_avrtot    ,ndimvor1,ndexvor1) ! First membre
(  508)          CALL histwrite( nidvor,"sovorgap",it,vor_avrres    ,ndimvor1,ndexvor1) ! gap between 1st and 2 nd mbre
(  509) 
(  510)          IF( idebug /= 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 10

(  511)             if(lwp)  WRITE(numout,*) ' debuging trd_vor: III.4 done'
(  512)             CALL FLUSH(numout)
(  513)          ENDIF
(  514) 
(  515)       ENDIF
(  516) 
(  517)       IF( MOD( kt - nit000+1, ntrd ) == 0 ) rotot(:,:)=0
(  518) 
(  519)       IF( kt == nitend )   CALL histclo( nidvor )
(  520) 
(  521)    END SUBROUTINE trd_vor
(  522) 
(  523) 
(  524) 
(  525)    SUBROUTINE trd_vor_init
(  526)       !!----------------------------------------------------------------------
(  527)       !!                  ***  ROUTINE trd_vor_init  ***
(  528)       !! 
(  529)       !! ** Purpose :   computation of vertically integrated T and S budgets
(  530)       !!      from ocean surface down to control surface (NetCDF output)
(  531)       !!
(  532)       !! ** Method/usage :
(  533)       !!
(  534)       !! History :
(  535)       !!   9.0  !  04-06  (L. Brunier, A-M. Treguier) Original code 
(  536)       !!        !  04-08  (C. Talandier) New trends organization
(  537)       !!----------------------------------------------------------------------
(  538)       !! * Local declarations
(  539)       REAL(wp) :: zjulian, zsto, zout
(  540) 
(  541)       CHARACTER (len=40) ::   clhstnam
(  542)       CHARACTER (len=40) ::   clop
(  543) 
(  544)       NAMELIST/namtrd/ ntrd,nctls
(  545)       !!----------------------------------------------------------------------
(  546) 
(  547)       !  ===================
(  548)       !   I. initialization
(  549)       !  ===================
(  550) 
(  551)       cvort='averaged-vor'
(  552) 
(  553)       ! Open specifier
(  554)       idebug = 0      ! set it to 1 in case of problem to have more Print
(  555) 
(  556)       ! namelist namtrd : trend diagnostic
(  557)       REWIND( numnam )
(  558)       READ  ( numnam, namtrd )
(  559) 
(  560)       IF(lwp) THEN
(  561)          WRITE(numout,*) ' '
(  562)          WRITE(numout,*) 'trd_vor_init: vorticity trends'
(  563)          WRITE(numout,*) '~~~~~~~~~~~~~'
(  564)          WRITE(numout,*) ' '
(  565)          WRITE(numout,*) '          Namelist namtrd : '
(  566)          WRITE(numout,*) '             time step frequency trend       ntrd  = ',ntrd
(  567)          WRITE(numout,*) ' '
(  568)          WRITE(numout,*) '##########################################################################'






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 11

(  569)          WRITE(numout,*) ' CAUTION: The interpretation of the vorticity trends is'
(  570)          WRITE(numout,*) ' not obvious, please contact Anne-Marie TREGUIER at: treguier@ifremer.fr '
(  571)          WRITE(numout,*) '##########################################################################'
(  572)          WRITE(numout,*) ' '
(  573)       ENDIF
(  574) 
(  575)       ! cumulated trends array init
(  576)       nmoydpvor = 0
(  577)       rotot(:,:)=0
(  578)       vor_avrtot(:,:)=0
(  579)       vor_avrres(:,:)=0
(  580) 
(  581)       IF( idebug /= 0 ) THEN
(  582)          if(lwp)  WRITE(numout,*) ' debuging trd_vor_init: I. done'
(  583)          CALL FLUSH(numout)
(  584)       ENDIF
(  585) 
(  586)       !  =================================
(  587)       !   II. netCDF output initialization
(  588)       !  =================================
(  589) 
(  590)       !-----------------------------------------
(  591)       ! II.1 Define frequency of output and means
(  592)       ! -----------------------------------------
(  593) #if defined key_diainstant
(  594)       zsto = nwrite*rdt
(  595)       clop ="inst(x)"
(  596) #else
(  597)       zsto = rdt
(  598)       clop ="ave(x)"
(  599) #endif
(  600)       zout = ntrd*rdt
(  601) 
(  602)       IF(lwp) WRITE (numout,*) ' trdvor_ncinit: netCDF initialization'
(  603) 
(  604)       ! II.2 Compute julian date from starting date of the run
(  605)       ! ------------------------
(  606)       CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  607)       IF (lwp) WRITE(numout,*)' '  
(  608)       IF (lwp) WRITE(numout,*)' Date 0 used :',nit000         &
(  609)            ,' YEAR ', nyear,' MONTH ', nmonth,' DAY ', nday   &
(  610)            ,'Julian day : ', zjulian
(  611) 
(  612)       ! II.3 Define the T grid trend file (nidvor)
(  613)       ! ---------------------------------
(  614)       CALL dia_nam( clhstnam, ntrd, 'vort' )                  ! filename
(  615)       IF(lwp) WRITE(numout,*) ' Name of NETCDF file ', clhstnam
(  616)       CALL histbeg( clhstnam, jpi, glamf, jpj, gphif,1, jpi,   &  ! Horizontal grid : glamt and gphit
(  617)          &          1, jpj, 0, zjulian, rdt, nh_t, nidvor, domain_id=nidom )
(  618)       CALL wheneq( jpi*jpj, fmask, 1, 1., ndexvor1, ndimvor1 )    ! surface
(  619) 
(  620)       ! Declare output fields as netCDF variables
(  621)       CALL histdef( nidvor, "sovortPh", cvort//"grad Ph" , "s-2",        & ! grad Ph
(  622)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  623)       CALL histdef( nidvor, "sovortEk", cvort//"Energy", "s-2",          & ! Energy
(  624)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  625)       CALL histdef( nidvor, "sovozeta", cvort//"rel vorticity", "s-2",   & ! rel vorticity
(  626)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 12

(  627)       CALL histdef( nidvor, "sovortif", cvort//"coriolis", "s-2",        & ! coriolis
(  628)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  629)       CALL histdef( nidvor, "sovodifl", cvort//"lat diff ", "s-2",       & ! lat diff
(  630)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  631)       CALL histdef( nidvor, "sovoadvv", cvort//"vert adv", "s-2",        & ! vert adv
(  632)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  633)       CALL histdef( nidvor, "sovodifv", cvort//"vert diff" , "s-2",      & ! vert diff
(  634)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  635)       CALL histdef( nidvor, "sovortPs", cvort//"grad Ps", "s-2",         & ! grad Ps
(  636)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  637)       CALL histdef( nidvor, "sovortbv", cvort//"Beta V", "s-2",          & ! beta.V
(  638)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  639)       CALL histdef( nidvor, "sovowind", cvort//"wind stress", "s-2",     & ! wind stress
(  640)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  641)       CALL histdef( nidvor, "sovobfri", cvort//"bottom friction", "s-2", & ! bottom friction
(  642)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  643)       CALL histdef( nidvor, "1st_mbre", cvort//"1st mbre", "s-2",        & ! First membre
(  644)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  645)       CALL histdef( nidvor, "sovorgap", cvort//"gap", "s-2",             & ! gap between 1st and 2 nd mbre
(  646)          &          jpi,jpj,nh_t,1,1,1,-99,32,clop,zsto,zout)
(  647)       CALL histend( nidvor )
(  648) 
(  649)       IF( idebug /= 0 ) THEN
(  650)          if(lwp)  WRITE(numout,*) ' debuging trd_vor_init: II. done'
(  651)          CALL FLUSH(numout)
(  652)       ENDIF
(  653) 
(  654)    END SUBROUTINE trd_vor_init
(  655) 
(  656) #else
(  657)    !!----------------------------------------------------------------------
(  658)    !!   Default option :                                       Empty module
(  659)    !!----------------------------------------------------------------------
(  660)    LOGICAL, PUBLIC ::   lk_trdvor = .FALSE.   ! momentum trend flag
(  661) 
(  662)    !! * Interfaces
(  663)    INTERFACE trd_vor_zint
(  664)       MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
(  665)    END INTERFACE
(  666) 
(  667) CONTAINS
(  668)    SUBROUTINE trd_vor( kt )        ! Empty routine
(  669) !      WRITE(*,*) 'trd_vor: You should not have seen this print! error?', kt
(  670)    END SUBROUTINE trd_vor
(  671)    SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
(  672)       REAL, DIMENSION(:,:), INTENT( inout ) ::   &
(  673)          putrdvor, pvtrdvor                  ! U and V momentum trends
(  674)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  675) !      WRITE(*,*) 'trd_vor_zint_2d: You should not have seen this print! error?', putrdvor(1,1)
(  676) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1)
(  677) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
(  678)    END SUBROUTINE trd_vor_zint_2d
(  679)    SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
(  680)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  681)          putrdvor, pvtrdvor                  ! U and V momentum trends
(  682)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  683) !      WRITE(*,*) 'trd_vor_zint_3d: You should not have seen this print! error?', putrdvor(1,1,1)
(  684) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1,1)






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 13

(  685) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
(  686)    END SUBROUTINE trd_vor_zint_3d
(  687)    SUBROUTINE trd_vor_init              ! Empty routine
(  688) !      WRITE(*,*) 'trd_vor_init: You should not have seen this print! error?'
(  689)    END SUBROUTINE trd_vor_init
(  690) #endif
(  691)    !!======================================================================
(  692) END MODULE trdvor
























































PGF90 (Version     10.5)          09/23/2010  09:49:10      page 14

(    1) # 1 "trdvor.F90"
(    1) MODULE trdvor
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdvor  ***
(    4)    !! Ocean diagnostics:  momentum trends
(    5)    !!=====================================================================
(    6)    
(    7) # 657
(  657)    !!----------------------------------------------------------------------
(  658)    !!   Default option :                                       Empty module
(  659)    !!----------------------------------------------------------------------
(  660)    LOGICAL, PUBLIC ::   lk_trdvor = .FALSE.   ! momentum trend flag
(  661) 
(  662)    !! * Interfaces
(  663)    INTERFACE trd_vor_zint
(  664)       MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
(  665)    END INTERFACE
(  666) 
(  667) CONTAINS
(  668)    SUBROUTINE trd_vor( kt )        ! Empty routine
(  669) !      WRITE(*,*) 'trd_vor: You should not have seen this print! error?', kt
(  670)    END SUBROUTINE trd_vor
(  671)    SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
(  672)       REAL, DIMENSION(:,:), INTENT( inout ) ::   &
(  673)          putrdvor, pvtrdvor                  ! U and V momentum trends
(  674)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  675) !      WRITE(*,*) 'trd_vor_zint_2d: You should not have seen this print! error?', putrdvor(1,1)
(  676) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1)
(  677) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
(  678)    END SUBROUTINE trd_vor_zint_2d
(  679)    SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
(  680)       REAL, DIMENSION(:,:,:), INTENT( inout ) ::   &
(  681)          putrdvor, pvtrdvor                  ! U and V momentum trends
(  682)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  683) !      WRITE(*,*) 'trd_vor_zint_3d: You should not have seen this print! error?', putrdvor(1,1,1)
(  684) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', pvtrdvor(1,1,1)
(  685) !      WRITE(*,*) '  "      "     : You should not have seen this print! error?', ktrd
(  686)    END SUBROUTINE trd_vor_zint_3d
(  687)    SUBROUTINE trd_vor_init              ! Empty routine
(  688) !      WRITE(*,*) 'trd_vor_init: You should not have seen this print! error?'
(  689)    END SUBROUTINE trd_vor_init
(  690) # 691
(  691)    !!======================================================================
(  692) END MODULE trdvor
