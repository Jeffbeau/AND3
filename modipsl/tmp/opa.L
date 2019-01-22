


PGF90 (Version     10.5)          09/23/2010  09:50:22      page 1

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

Filename: opa.F90

(    1) MODULE opa
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE opa   ***
(    4)    !! Ocean system   : OPA ocean dynamics (including on-line tracers and sea-ice)
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   opa_model      : solve ocean dynamics, tracer and/or sea-ice
(    9)    !!   opa_flg        : initialisation of algorithm flag 
(   10)    !!   opa_closefile  : close remaining files
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE cpl_oce         ! ocean-atmosphere-sea ice coupled exchanges
(   14)    USE dom_oce         ! ocean space domain variables
(   15)    USE oce             ! dynamics and tracers variables
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE daymod          ! calendar
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lib_mpp         ! distributed memory computing
(   20) 
(   21)    USE domcfg          ! domain configuration               (dom_cfg routine)
(   22)    USE mppini          ! shared/distributed memory setting (mpp_init routine)
(   23)    USE domain          ! domain initialization             (dom_init routine)
(   24)    USE obc_par         ! open boundary cond. parameters
(   25)    USE obcini          ! open boundary cond. initialization (obc_ini routine)
(   26)    USE solver          ! solver initialization          (solver_init routine)
(   27)    USE istate          ! initial state setting          (istate_init routine)
(   28)    USE eosbn2          ! equation of state            (eos bn2 routine)
(   29)    USE zpshde          ! partial step: hor. derivative (zps_hde routine)
(   30) 
(   31)    ! ocean physics
(   32)    USE traqsr          ! solar radiation penetration   (tra_qsr_init routine)
(   33)    USE ldfdyn          ! lateral viscosity setting      (ldfdyn_init routine)
(   34)    USE ldftra          ! lateral diffusivity setting    (ldftra_init routine)
(   35)    USE zdfini
(   36) 
(   37)    USE phycst          ! physical constant                  (par_cst routine)
(   38)    USE iceini          ! initialization of sea-ice         (ice_init routine)
(   39)    USE cpl             ! coupled ocean/atmos.              (cpl_init routine)
(   40)    USE ocfzpt          ! ocean freezing point              (oc_fz_pt routine)
(   41)    USE trdicp          ! momentum/tracers trends       (trd_icp_init routine)
(   42)    USE trdvor          ! vorticity trends              (trd_vor_init routine)
(   43)    USE trdmld          ! tracer mixed layer trends     (trd_mld_init routine)
(   44)    USE flxfwb          ! freshwater budget correction  (flx_fwb_init routine)
(   45)    USE flxmod          ! thermohaline forcing of the ocean (flx_init routine)
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 2

(   47)    USE diaptr          ! poleward transports           (dia_ptr_init routine)
(   48) 
(   49)    USE step            ! OPA time-stepping                  (stp     routine)
(   50)    USE dynspg_oce      ! Control choice of surface pressure gradient schemes
(   51)    USE prtctl          ! Print control                 (prt_ctl_init routine)
(   52)    USE ini1d           ! re-initialization of u-v mask for the 1D configuration
(   53)    USE dyncor1d        ! Coriolis factor at T-point
(   54)    USE step1d          ! Time stepping loop for the 1D configuration
(   55) 
(   56)    USE initrc          ! Initialization of the passive tracers
(   57) 
(   58)    IMPLICIT NONE
(   59)    PRIVATE
(   60) 
(   61)    !! * Routine accessibility
(   62)    PUBLIC opa_model      ! called by model.F90
(   63)    PUBLIC opa_init
(   64)    !!----------------------------------------------------------------------
(   65)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   66)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/opa.F90,v 1.23 2006/04/19 14:43:14 opalod Exp $ 
(   67)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   68)    !!----------------------------------------------------------------------
(   69) 
(   70) CONTAINS
(   71) 
(   72)    SUBROUTINE opa_model
(   73)       !!----------------------------------------------------------------------
(   74)       !!                     ***  ROUTINE opa  ***
(   75)       !!
(   76)       !! ** Purpose :   opa solves the primitive equations on an orthogonal 
(   77)       !!      curvilinear mesh on the sphere.
(   78)       !!
(   79)       !! ** Method  : - model general initialization
(   80)       !!              - launch the time-stepping (stp routine)
(   81)       !!
(   82)       !! References :
(   83)       !!      Madec, Delecluse,Imbard, and Levy, 1997: reference manual.
(   84)       !!              internal report, IPSL.
(   85)       !!
(   86)       !! History :
(   87)       !!   4.0  !  90-10  (C. Levy, G. Madec)  Original code
(   88)       !!   7.0  !  91-11  (M. Imbard, C. Levy, G. Madec)
(   89)       !!   7.1  !  93-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
(   90)       !!                   M. Guyon, A. Lazar, P. Delecluse, C. Perigaud,
(   91)       !!                   G. Caniaux, B. Colot, C. Maes ) release 7.1 
(   92)       !!        !  92-06  (L.Terray) coupling implementation
(   93)       !!        !  93-11  (M.A. Filiberti) IGLOO sea-ice 
(   94)       !!   8.0  !  96-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
(   95)       !!                   M. Guyon, A. Lazar, P. Delecluse, L.Terray,
(   96)       !!                   M.A. Filiberti, J. Vialar, A.M. Treguier,
(   97)       !!                   M. Levy)  release 8.0
(   98)       !!   8.1  !  97-06  (M. Imbard, G. Madec)
(   99)       !!   8.2  !  99-11  (M. Imbard, H. Goosse)  LIM sea-ice model 
(  100)       !!        !  99-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP 
(  101)       !!        !  00-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
(  102)       !!   9.0  !  02-08  (G. Madec)  F90: Free form and modules
(  103)       !!    "   !  04-08  (C. Talandier) New trends organization
(  104)       !!    "   !  05-06  (C. Ethe) Add the 1D configuration possibility






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 3

(  105)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  106)       !!----------------------------------------------------------------------
(  107)       !! * Local declarations
(  108)       INTEGER ::   istp       ! time step index
(  109)       CHARACTER (len=64) ::        &
(  110)          cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
(  111)       !!----------------------------------------------------------------------
(  112) 
(  113) #if defined key_agrif
(  114) 
(  115)       Call Agrif_Init_Grids()
(  116) #endif
(  117)       
(  118)       Call opa_init  ! Initializations
(  119) 
(  120)       IF( lk_cfg_1d  )  THEN 
(  121)          istp = nit000
(  122)          DO WHILE ( istp <= nitend .AND. nstop == 0 )
(  123) #if defined key_agrif
(  124)             CALL Agrif_Step(stp_1d)
(  125) #else
(  126)             CALL stp_1d( istp )
(  127) #endif
(  128)             istp = istp + 1
(  129)          END DO
(  130)       ELSE
(  131)          istp = nit000
(  132)          DO WHILE ( istp <= nitend .AND. nstop == 0 )
(  133) #if defined key_agrif
(  134)             CALL Agrif_Step(stp)
(  135) #else
(  136)             CALL stp( istp )
(  137) #endif
(  138)             istp = istp + 1
(  139)          END DO
(  140)       ENDIF
(  141)       !                                     ! ========= !
(  142)       !                                     !  Job end  !
(  143)       !                                     ! ========= !
(  144) 
(  145)       IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  146) 
(  147)       IF( nstop /= 0 ) THEN                 ! error print
(  148)       IF(lwp) WRITE(numout,cform_err)
(  149)       IF(lwp) WRITE(numout,*) nstop, ' error have been found' 
(  150)       ENDIF
(  151) 
(  152)       CALL opa_closefile
(  153)       IF( lk_mpp )   CALL mppstop                          ! Close all files (mpp)
(  154) 
(  155)    END SUBROUTINE opa_model
(  156) 
(  157) 
(  158)    SUBROUTINE opa_flg
(  159)       !!----------------------------------------------------------------------
(  160)       !!                     ***  ROUTINE opa  ***
(  161)       !!
(  162)       !! ** Purpose :   Initialize logical flags that control the choice of






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 4

(  163)       !!      some algorithm or control print
(  164)       !!
(  165)       !! ** Method  :    Read in namilist namflg logical flags
(  166)       !!
(  167)       !! History :
(  168)       !!   9.0  !  03-11  (G. Madec)  Original code
(  169)       !!----------------------------------------------------------------------
(  170)       !! * Local declarations
(  171) 
(  172)       NAMELIST/namflg/ ln_dynhpg_imp
(  173)       !!----------------------------------------------------------------------
(  174) 
(  175)       ! Read Namelist namflg : algorithm FLaG 
(  176)       ! --------------------
(  177)       REWIND ( numnam )
(  178)       READ   ( numnam, namflg )
(  179) 
(  180)       ! Parameter control and print
(  181)       ! ---------------------------
(  182)       ! Control print
(  183)       IF(lwp) THEN
(  184)          WRITE(numout,*)
(  185)          WRITE(numout,*) 'opa_flg : algorithm flag initialization'
(  186)          WRITE(numout,*) '~~~~~~~'
(  187)          WRITE(numout,*) '          Namelist namflg : set algorithm flags'
(  188)          WRITE(numout,*)
(  189)          WRITE(numout,*) '             centered (F) or semi-implicit (T)   ln_dynhpg_imp = ', ln_dynhpg_imp
(  190)          WRITE(numout,*) '             hydrostatic pressure gradient'
(  191)       ENDIF
(  192) 
(  193)    END SUBROUTINE opa_flg
(  194) 
(  195)    SUBROUTINE opa_closefile
(  196)       !!----------------------------------------------------------------------
(  197)       !!                     ***  ROUTINE opa_closefile  ***
(  198)       !!
(  199)       !! ** Purpose :   Close the files
(  200)       !!           
(  201)       !! ** Method  : 
(  202)       !!
(  203)       !! History :
(  204)       !!   9.0  !  05-01  (O. Le Galloudec)  Original code
(  205)       !!----------------------------------------------------------------------
(  206)       !! * Modules used
(  207)       USE dtatem        ! temperature data
(  208)       USE dtasal        ! salinity data
(  209)       USE dtasst        ! sea surface temperature data
(  210)       !!----------------------------------------------------------------------
(  211) 
(  212)       IF ( lk_mpp ) CALL mppsync
(  213) 
(  214)       ! 1. Unit close
(  215)       ! -------------
(  216) 
(  217)       CLOSE( numnam )       ! namelist
(  218)       CLOSE( numout )       ! standard model output file
(  219)       CLOSE( numstp )       ! time-step file
(  220)       CLOSE( numwrs )       ! ocean restart file






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 5

(  221) 
(  222)       IF( lk_dtatem )   CLOSE( numtdt )
(  223)       IF( lk_dtasal )   CLOSE( numsdt )
(  224)       IF( lk_dtasst )   CLOSE( numsst )
(  225) 
(  226)       IF(lwp) CLOSE( numsol )
(  227) 
(  228)       IF( lk_cpl ) THEN
(  229)          CLOSE( numlhf )
(  230)          CLOSE( numlts )
(  231)       ENDIF
(  232) 
(  233)       CLOSE( numwri )
(  234) 
(  235)    END SUBROUTINE opa_closefile
(  236) 
(  237)    !!======================================================================
(  238)    SUBROUTINE opa_init
(  239)       !!----------------------------------------------------------------------
(  240)       !!                     ***  ROUTINE opa_init  ***
(  241)       !!
(  242)       !! ** Purpose :   initialization of the opa model
(  243)       !!
(  244)       !! ** Method  : 
(  245)       !!
(  246)       !! References :
(  247)       !!----------------------------------------------------------------------
(  248)       !! * Local declarations
(  249) 
(  250) #if defined key_coupled
(  251)       INTEGER ::   itro, istp0        ! ???
(  252) #endif
(  253)       CHARACTER (len=64) ::        &
(  254)          cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
(  255)       CHARACTER (len=20) :: namelistname
(  256)       CHARACTER (len=28) :: file_out
(  257)       !!----------------------------------------------------------------------
(  258) 
(  259)       ! Initializations
(  260)       ! ===============
(  261) 
(  262)       file_out = 'ocean.output'
(  263)       
(  264)       ! open listing and namelist units
(  265)       IF ( numout /= 0 .AND. numout /= 6 ) THEN 
(  266)          CALL ctlopn(numout,file_out,'UNKNOWN', 'FORMATTED',   &
(  267)                       'SEQUENTIAL',1,numout,.FALSE.,1)
(  268) !         OPEN( UNIT=numout, FILE=TRIM(file_out), FORM='FORMATTED' )
(  269)       ENDIF
(  270) 
(  271)       namelistname = 'namelist'
(  272)       CALL ctlopn(numnam,namelistname,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
(  273)                      1,numout,.FALSE.,1)
(  274) !!!!      OPEN( UNIT=numnam, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )
(  275) 
(  276) !!DB 2009.09.09 -- restrict write to lwp
(  277) !      WRITE(numout,*)
(  278) !      WRITE(numout,*) '                 L O D Y C - I P S L'






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 6

(  279) !      WRITE(numout,*) '                     O P A model'
(  280) !      WRITE(numout,*) '            Ocean General Circulation Model'
(  281) !      WRITE(numout,*) '               version OPA 9.0  (2005) '
(  282) !      WRITE(numout,*)
(  283) !      WRITE(numout,*)
(  284) 
(  285)       ! Nodes selection
(  286)       narea = mynode()
(  287)       narea = narea + 1    ! mynode return the rank of proc (0 --> jpnij -1 )
(  288)       lwp   = narea == 1
(  289) 
(  290)       !                                     ! ============================== !
(  291)       !                                     !  Model general initialization  !
(  292)       !                                     ! ============================== !
(  293) 
(  294)       IF(lwp) then
(  295)          WRITE(numout,*)
(  296)          WRITE(numout,*) '                 L O D Y C - I P S L'
(  297)          WRITE(numout,*) '                     O P A model'
(  298)          WRITE(numout,*) '            Ocean General Circulation Model'
(  299)          WRITE(numout,*) '               version OPA 9.0  (2005) '
(  300)          WRITE(numout,*)
(  301)          WRITE(numout,*)
(  302)          WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  303)       endif
(  304)       
(  305) !!DB 2009.10.01 open run_info.out file as numout2 (in_out_manager)
(  306)       if(lwp) then 
(  307)          open(numout2, file = 'sopa_run_info.out',STATUS='replace')
(  308)          write(numout2,*) ' SOPA MC runtime information file '
(  309)          write(numout2,*) ' -------------------------------- '
(  310)          write(numout2,*) '                                  '
(  311)          write(numout2,*) '                 L O D Y C - I P S L'
(  312)          write(numout2,*) '                     O P A model'
(  313)          write(numout2,*) '            Ocean General Circulation Model'
(  314)          write(numout2,*) '               version OPA 9.0  (2005) '
(  315)          write(numout2,*) '                                  '
(  316)          write(numout2,*) '                                  '
(  317)       endif
(  318) 
(  319)                                             ! Domain decomposition
(  320)       IF( jpni*jpnj == jpnij ) THEN
(  321)          CALL mpp_init                          ! standard cutting out
(  322)       ELSE
(  323)          CALL mpp_init2                         ! eliminate land processors
(  324)       ENDIF
(  325)       
(  326)       CALL phy_cst                          ! Physical constants
(  327) 
(  328)       CALL dom_cfg                          ! Domain configuration
(  329)       
(  330)       CALL dom_init                         ! Domain
(  331)       IF( ln_ctl )      CALL prt_ctl_init   ! Print control
(  332) 
(  333)       IF( lk_cfg_1d )   CALL fcorio_1d      ! redefine Coriolis at T-point
(  334) 
(  335)       IF( lk_obc    )   CALL obc_init       ! Open boundaries 
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 7

(  337)       CALL day( nit000 )                    ! Calendar
(  338) 
(  339)       CALL istate_init                      ! ocean initial state (Dynamics and tracers)
(  340) 
(  341)       IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
(  342)          CALL solver_init( nit000 )         ! Elliptic solver
(  343)       ENDIF
(  344) 
(  345) !!add
(  346)                        CALL eos( tb, sb, rhd, rhop )        ! before potential and in situ densities
(  347)                        
(  348)                        CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
(  349) 
(  350)       IF( lk_zps .AND. .NOT. lk_cfg_1d )   &
(  351)          &             CALL zps_hde( nit000, tb, sb, rhd,  &  ! Partial steps: before Horizontal DErivative
(  352)                                             gtu, gsu, gru, &  ! of t, s, rd at the bottom ocean level
(  353)                                             gtv, gsv, grv )
(  354) 
(  355) !!add
(  356) 
(  357)       CALL oc_fz_pt                         ! Surface freezing point
(  358) 
(  359) #if defined key_ice_lim
(  360)       CALL ice_init                         ! Sea ice model
(  361) #endif
(  362) 
(  363)       !                                     ! Ocean scheme
(  364) 
(  365)       CALL opa_flg                              ! Choice of algorithms
(  366) 
(  367)       !                                     ! Ocean physics
(  368) 
(  369)       CALL tra_qsr_init                         ! Solar radiation penetration
(  370) 
(  371)       CALL ldf_dyn_init                         ! Lateral ocean momentum physics
(  372) 
(  373)       CALL ldf_tra_init                         ! Lateral ocean tracer physics
(  374) 
(  375)       CALL zdf_init                             ! Vertical ocean physics
(  376) 
(  377)       !                                     ! Ocean trends
(  378)       ! Control parameters 
(  379)       IF( lk_trdtra .OR. lk_trdmld )   l_trdtra = .TRUE.
(  380)       IF( lk_trddyn .OR. lk_trdvor )   l_trddyn = .TRUE.
(  381) 
(  382)       IF( lk_trddyn .OR. lk_trdtra )   &
(  383)          &            CALL trd_icp_init         ! active tracers and/or momentum
(  384) 
(  385)       IF( lk_trdmld ) CALL trd_mld_init         ! mixed layer
(  386) 
(  387)       IF( lk_trdvor ) CALL trd_vor_init         ! vorticity
(  388) 
(  389) #if defined key_passivetrc
(  390)       CALL ini_trc                           ! Passive tracers
(  391) #endif
(  392) 
(  393) #if defined key_coupled
(  394)       itro  = nitend - nit000 + 1           ! Coupled






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 8

(  395)       istp0 = NINT( rdt )
(  396)       CALL cpl_init( itro, nexco, istp0 )   ! Signal processing and process id exchange
(  397) #endif
(  398) 
(  399)       CALL flx_init                         ! Thermohaline forcing initialization
(  400) 
(  401)       CALL flx_fwb_init                     ! FreshWater Budget correction
(  402) 
(  403)       CALL dia_ptr_init                     ! Poleward TRansports initialization
(  404) 
(  405)       !                                     ! =============== !
(  406)       !                                     !  time stepping  !
(  407)       !                                     ! =============== !
(  408) 
(  409) 
(  410) !!DB: 2007.12.06
(  411) !!Code for additional control parameters
(  412) !!Look for a file called run_params and get run params
(  413) !!This is a HARDWIRED list. At this time, the default values
(  414) !!should be found at the bottom of oce.F90
(  415) !!Exit gracefully if file not found
(  416)       if(lwp)write(numout2,*)'DBG Default perpetual_forcing,M2_ave,ioutput_ave: ', &
(  417)            perpetual_forcing,M2_ave,ioutput_ave
(  418)       open(66,file='run_params',STATUS='OLD',ERR=666)
(  419)       read(66,*,ERR=666)perpetual_forcing
(  420)       read(66,*,ERR=666)M2_ave
(  421)       read(66,*,ERR=666)ioutput_ave
(  422)       close(66)
(  423) !!DBG
(  424)       if(lwp) then
(  425)          write(numout2,*)'DBG: run_params --> perpetual_forcing,M2_ave,ioutput_ave: ', &
(  426)               perpetual_forcing,M2_ave, ioutput_ave
(  427)       endif
(  428) 666   continue
(  429) 
(  430) !!DB 2009.10.01
(  431)       if(ioutput_ave /= 0) then
(  432)          ioutput_ave = nwrite
(  433)          if(lwp) write(numout2,*)'DBG: ioutput_ave /= 0 ; setting to nwrite  = ', nwrite
(  434)       endif
(  435) 
(  436) 
(  437)       IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  438) 
(  439)       IF( lk_cfg_1d  )  THEN 
(  440)          CALL init_1d
(  441)       ENDIF
(  442)    END SUBROUTINE opa_init
(  443)    !!======================================================================
(  444) END MODULE opa














PGF90 (Version     10.5)          09/23/2010  09:50:22      page 9

(    1) # 1 "opa.F90"
(    1) MODULE opa
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE opa   ***
(    4)    !! Ocean system   : OPA ocean dynamics (including on-line tracers and sea-ice)
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   opa_model      : solve ocean dynamics, tracer and/or sea-ice
(    9)    !!   opa_flg        : initialisation of algorithm flag 
(   10)    !!   opa_closefile  : close remaining files
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE cpl_oce         ! ocean-atmosphere-sea ice coupled exchanges
(   14)    USE dom_oce         ! ocean space domain variables
(   15)    USE oce             ! dynamics and tracers variables
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE daymod          ! calendar
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lib_mpp         ! distributed memory computing
(   20) 
(   21)    USE domcfg          ! domain configuration               (dom_cfg routine)
(   22)    USE mppini          ! shared/distributed memory setting (mpp_init routine)
(   23)    USE domain          ! domain initialization             (dom_init routine)
(   24)    USE obc_par         ! open boundary cond. parameters
(   25)    USE obcini          ! open boundary cond. initialization (obc_ini routine)
(   26)    USE solver          ! solver initialization          (solver_init routine)
(   27)    USE istate          ! initial state setting          (istate_init routine)
(   28)    USE eosbn2          ! equation of state            (eos bn2 routine)
(   29)    USE zpshde          ! partial step: hor. derivative (zps_hde routine)
(   30) 
(   31)    ! ocean physics
(   32)    USE traqsr          ! solar radiation penetration   (tra_qsr_init routine)
(   33)    USE ldfdyn          ! lateral viscosity setting      (ldfdyn_init routine)
(   34)    USE ldftra          ! lateral diffusivity setting    (ldftra_init routine)
(   35)    USE zdfini
(   36) 
(   37)    USE phycst          ! physical constant                  (par_cst routine)
(   38)    USE iceini          ! initialization of sea-ice         (ice_init routine)
(   39)    USE cpl             ! coupled ocean/atmos.              (cpl_init routine)
(   40)    USE ocfzpt          ! ocean freezing point              (oc_fz_pt routine)
(   41)    USE trdicp          ! momentum/tracers trends       (trd_icp_init routine)
(   42)    USE trdvor          ! vorticity trends              (trd_vor_init routine)
(   43)    USE trdmld          ! tracer mixed layer trends     (trd_mld_init routine)
(   44)    USE flxfwb          ! freshwater budget correction  (flx_fwb_init routine)
(   45)    USE flxmod          ! thermohaline forcing of the ocean (flx_init routine)
(   46) 
(   47)    USE diaptr          ! poleward transports           (dia_ptr_init routine)
(   48) 
(   49)    USE step            ! OPA time-stepping                  (stp     routine)
(   50)    USE dynspg_oce      ! Control choice of surface pressure gradient schemes
(   51)    USE prtctl          ! Print control                 (prt_ctl_init routine)
(   52)    USE ini1d           ! re-initialization of u-v mask for the 1D configuration
(   53)    USE dyncor1d        ! Coriolis factor at T-point
(   54)    USE step1d          ! Time stepping loop for the 1D configuration
(   55) 
(   56)    USE initrc          ! Initialization of the passive tracers
(   57) 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 10

(   58)    IMPLICIT NONE
(   59)    PRIVATE
(   60) 
(   61)    !! * Routine accessibility
(   62)    PUBLIC opa_model      ! called by model.F90
(   63)    PUBLIC opa_init
(   64)    !!----------------------------------------------------------------------
(   65)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   66)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/opa.F90,v 1.23 2006/04/19 14:43:14 opalod Exp $ 
(   67)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   68)    !!----------------------------------------------------------------------
(   69) 
(   70) CONTAINS
(   71) 
(   72)    SUBROUTINE opa_model
(   73)       !!----------------------------------------------------------------------
(   74)       !!                     ***  ROUTINE opa  ***
(   75)       !!
(   76)       !! ** Purpose :   opa solves the primitive equations on an orthogonal 
(   77)       !!      curvilinear mesh on the sphere.
(   78)       !!
(   79)       !! ** Method  : - model general initialization
(   80)       !!              - launch the time-stepping (stp routine)
(   81)       !!
(   82)       !! References :
(   83)       !!      Madec, Delecluse,Imbard, and Levy, 1997: reference manual.
(   84)       !!              internal report, IPSL.
(   85)       !!
(   86)       !! History :
(   87)       !!   4.0  !  90-10  (C. Levy, G. Madec)  Original code
(   88)       !!   7.0  !  91-11  (M. Imbard, C. Levy, G. Madec)
(   89)       !!   7.1  !  93-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
(   90)       !!                   M. Guyon, A. Lazar, P. Delecluse, C. Perigaud,
(   91)       !!                   G. Caniaux, B. Colot, C. Maes ) release 7.1 
(   92)       !!        !  92-06  (L.Terray) coupling implementation
(   93)       !!        !  93-11  (M.A. Filiberti) IGLOO sea-ice 
(   94)       !!   8.0  !  96-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
(   95)       !!                   M. Guyon, A. Lazar, P. Delecluse, L.Terray,
(   96)       !!                   M.A. Filiberti, J. Vialar, A.M. Treguier,
(   97)       !!                   M. Levy)  release 8.0
(   98)       !!   8.1  !  97-06  (M. Imbard, G. Madec)
(   99)       !!   8.2  !  99-11  (M. Imbard, H. Goosse)  LIM sea-ice model 
(  100)       !!        !  99-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP 
(  101)       !!        !  00-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
(  102)       !!   9.0  !  02-08  (G. Madec)  F90: Free form and modules
(  103)       !!    "   !  04-08  (C. Talandier) New trends organization
(  104)       !!    "   !  05-06  (C. Ethe) Add the 1D configuration possibility
(  105)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  106)       !!----------------------------------------------------------------------
(  107)       !! * Local declarations
(  108)       INTEGER ::   istp       ! time step index
(  109)       CHARACTER (len=64) ::        &
(  110)          cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
(  111)       !!----------------------------------------------------------------------
(  112) 
(  113) # 117
(  117)       
(  118)       Call opa_init  ! Initializations






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 11

(  119) 
(  120)       IF( lk_cfg_1d  )  THEN 
(  121)          istp = nit000
(  122)          DO WHILE ( istp <= nitend .AND. nstop == 0 )
(  123) # 126
(  126)             CALL stp_1d( istp )
(  127) # 128
(  128)             istp = istp + 1
(  129)          END DO
(  130)       ELSE
(  131)          istp = nit000
(  132)          DO WHILE ( istp <= nitend .AND. nstop == 0 )
(  133) # 136
(  136)             CALL stp( istp )
(  137) # 138
(  138)             istp = istp + 1
(  139)          END DO
(  140)       ENDIF
(  141)       !                                     ! ========= !
(  142)       !                                     !  Job end  !
(  143)       !                                     ! ========= !
(  144) 
(  145)       IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  146) 
(  147)       IF( nstop /= 0 ) THEN                 ! error print
(  148)       IF(lwp) WRITE(numout,cform_err)
(  149)       IF(lwp) WRITE(numout,*) nstop, ' error have been found' 
(  150)       ENDIF
(  151) 
(  152)       CALL opa_closefile
(  153)       IF( lk_mpp )   CALL mppstop                          ! Close all files (mpp)
(  154) 
(  155)    END SUBROUTINE opa_model
(  156) 
(  157) 
(  158)    SUBROUTINE opa_flg
(  159)       !!----------------------------------------------------------------------
(  160)       !!                     ***  ROUTINE opa  ***
(  161)       !!
(  162)       !! ** Purpose :   Initialize logical flags that control the choice of
(  163)       !!      some algorithm or control print
(  164)       !!
(  165)       !! ** Method  :    Read in namilist namflg logical flags
(  166)       !!
(  167)       !! History :
(  168)       !!   9.0  !  03-11  (G. Madec)  Original code
(  169)       !!----------------------------------------------------------------------
(  170)       !! * Local declarations
(  171) 
(  172)       NAMELIST/namflg/ ln_dynhpg_imp
(  173)       !!----------------------------------------------------------------------
(  174) 
(  175)       ! Read Namelist namflg : algorithm FLaG 
(  176)       ! --------------------
(  177)       REWIND ( numnam )
(  178)       READ   ( numnam, namflg )
(  179) 
(  180)       ! Parameter control and print






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 12

(  181)       ! ---------------------------
(  182)       ! Control print
(  183)       IF(lwp) THEN
(  184)          WRITE(numout,*)
(  185)          WRITE(numout,*) 'opa_flg : algorithm flag initialization'
(  186)          WRITE(numout,*) '~~~~~~~'
(  187)          WRITE(numout,*) '          Namelist namflg : set algorithm flags'
(  188)          WRITE(numout,*)
(  189)          WRITE(numout,*) '             centered (F) or semi-implicit (T)   ln_dynhpg_imp = ', ln_dynhpg_imp
(  190)          WRITE(numout,*) '             hydrostatic pressure gradient'
(  191)       ENDIF
(  192) 
(  193)    END SUBROUTINE opa_flg
(  194) 
(  195)    SUBROUTINE opa_closefile
(  196)       !!----------------------------------------------------------------------
(  197)       !!                     ***  ROUTINE opa_closefile  ***
(  198)       !!
(  199)       !! ** Purpose :   Close the files
(  200)       !!           
(  201)       !! ** Method  : 
(  202)       !!
(  203)       !! History :
(  204)       !!   9.0  !  05-01  (O. Le Galloudec)  Original code
(  205)       !!----------------------------------------------------------------------
(  206)       !! * Modules used
(  207)       USE dtatem        ! temperature data
(  208)       USE dtasal        ! salinity data
(  209)       USE dtasst        ! sea surface temperature data
(  210)       !!----------------------------------------------------------------------
(  211) 
(  212)       IF ( lk_mpp ) CALL mppsync
(  213) 
(  214)       ! 1. Unit close
(  215)       ! -------------
(  216) 
(  217)       CLOSE( numnam )       ! namelist
(  218)       CLOSE( numout )       ! standard model output file
(  219)       CLOSE( numstp )       ! time-step file
(  220)       CLOSE( numwrs )       ! ocean restart file
(  221) 
(  222)       IF( lk_dtatem )   CLOSE( numtdt )
(  223)       IF( lk_dtasal )   CLOSE( numsdt )
(  224)       IF( lk_dtasst )   CLOSE( numsst )
(  225) 
(  226)       IF(lwp) CLOSE( numsol )
(  227) 
(  228)       IF( lk_cpl ) THEN
(  229)          CLOSE( numlhf )
(  230)          CLOSE( numlts )
(  231)       ENDIF
(  232) 
(  233)       CLOSE( numwri )
(  234) 
(  235)    END SUBROUTINE opa_closefile
(  236) 
(  237)    !!======================================================================
(  238)    SUBROUTINE opa_init






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 13

(  239)       !!----------------------------------------------------------------------
(  240)       !!                     ***  ROUTINE opa_init  ***
(  241)       !!
(  242)       !! ** Purpose :   initialization of the opa model
(  243)       !!
(  244)       !! ** Method  : 
(  245)       !!
(  246)       !! References :
(  247)       !!----------------------------------------------------------------------
(  248)       !! * Local declarations
(  249) 
(  250) # 253
(  253)       CHARACTER (len=64) ::        &
(  254)          cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
(  255)       CHARACTER (len=20) :: namelistname
(  256)       CHARACTER (len=28) :: file_out
(  257)       !!----------------------------------------------------------------------
(  258) 
(  259)       ! Initializations
(  260)       ! ===============
(  261) 
(  262)       file_out = 'ocean.output'
(  263)       
(  264)       ! open listing and namelist units
(  265)       IF ( numout /= 0 .AND. numout /= 6 ) THEN 
(  266)          CALL ctlopn(numout,file_out,'UNKNOWN', 'FORMATTED',   &
(  267)                       'SEQUENTIAL',1,numout,.FALSE.,1)
(  268) !         OPEN( UNIT=numout, FILE=TRIM(file_out), FORM='FORMATTED' )
(  269)       ENDIF
(  270) 
(  271)       namelistname = 'namelist'
(  272)       CALL ctlopn(numnam,namelistname,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
(  273)                      1,numout,.FALSE.,1)
(  274) !!!!      OPEN( UNIT=numnam, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )
(  275) 
(  276) !!DB 2009.09.09 -- restrict write to lwp
(  277) !      WRITE(numout,*)
(  278) !      WRITE(numout,*) '                 L O D Y C - I P S L'
(  279) !      WRITE(numout,*) '                     O P A model'
(  280) !      WRITE(numout,*) '            Ocean General Circulation Model'
(  281) !      WRITE(numout,*) '               version OPA 9.0  (2005) '
(  282) !      WRITE(numout,*)
(  283) !      WRITE(numout,*)
(  284) 
(  285)       ! Nodes selection
(  286)       narea = mynode()
(  287)       narea = narea + 1    ! mynode return the rank of proc (0 --> jpnij -1 )
(  288)       lwp   = narea == 1
(  289) 
(  290)       !                                     ! ============================== !
(  291)       !                                     !  Model general initialization  !
(  292)       !                                     ! ============================== !
(  293) 
(  294)       IF(lwp) then
(  295)          WRITE(numout,*)
(  296)          WRITE(numout,*) '                 L O D Y C - I P S L'
(  297)          WRITE(numout,*) '                     O P A model'
(  298)          WRITE(numout,*) '            Ocean General Circulation Model'






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 14

(  299)          WRITE(numout,*) '               version OPA 9.0  (2005) '
(  300)          WRITE(numout,*)
(  301)          WRITE(numout,*)
(  302)          WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  303)       endif
(  304)       
(  305) !!DB 2009.10.01 open run_info.out file as numout2 (in_out_manager)
(  306)       if(lwp) then 
(  307)          open(numout2, file = 'sopa_run_info.out',STATUS='replace')
(  308)          write(numout2,*) ' SOPA MC runtime information file '
(  309)          write(numout2,*) ' -------------------------------- '
(  310)          write(numout2,*) '                                  '
(  311)          write(numout2,*) '                 L O D Y C - I P S L'
(  312)          write(numout2,*) '                     O P A model'
(  313)          write(numout2,*) '            Ocean General Circulation Model'
(  314)          write(numout2,*) '               version OPA 9.0  (2005) '
(  315)          write(numout2,*) '                                  '
(  316)          write(numout2,*) '                                  '
(  317)       endif
(  318) 
(  319)                                             ! Domain decomposition
(  320)       IF( jpni*jpnj == jpnij ) THEN
(  321)          CALL mpp_init                          ! standard cutting out
(  322)       ELSE
(  323)          CALL mpp_init2                         ! eliminate land processors
(  324)       ENDIF
(  325)       
(  326)       CALL phy_cst                          ! Physical constants
(  327) 
(  328)       CALL dom_cfg                          ! Domain configuration
(  329)       
(  330)       CALL dom_init                         ! Domain
(  331)       IF( ln_ctl )      CALL prt_ctl_init   ! Print control
(  332) 
(  333)       IF( lk_cfg_1d )   CALL fcorio_1d      ! redefine Coriolis at T-point
(  334) 
(  335)       IF( lk_obc    )   CALL obc_init       ! Open boundaries 
(  336) 
(  337)       CALL day( nit000 )                    ! Calendar
(  338) 
(  339)       CALL istate_init                      ! ocean initial state (Dynamics and tracers)
(  340) 
(  341)       IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
(  342)          CALL solver_init( nit000 )         ! Elliptic solver
(  343)       ENDIF
(  344) 
(  345) !!add
(  346)                        CALL eos( tb, sb, rhd, rhop )        ! before potential and in situ densities
(  347)                        
(  348)                        CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
(  349) 
(  350)       IF( lk_zps .AND. .NOT. lk_cfg_1d )   &
(  351)          &             CALL zps_hde( nit000, tb, sb, rhd,  &  ! Partial steps: before Horizontal DErivative
(  352)                                             gtu, gsu, gru, &  ! of t, s, rd at the bottom ocean level
(  353)                                             gtv, gsv, grv )
(  354) 
(  355) !!add
(  356) 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 15

(  357)       CALL oc_fz_pt                         ! Surface freezing point
(  358) 
(  359) # 360
(  360)       CALL ice_init                         ! Sea ice model
(  361) 
(  362) # 363
(  363)       !                                     ! Ocean scheme
(  364) 
(  365)       CALL opa_flg                              ! Choice of algorithms
(  366) 
(  367)       !                                     ! Ocean physics
(  368) 
(  369)       CALL tra_qsr_init                         ! Solar radiation penetration
(  370) 
(  371)       CALL ldf_dyn_init                         ! Lateral ocean momentum physics
(  372) 
(  373)       CALL ldf_tra_init                         ! Lateral ocean tracer physics
(  374) 
(  375)       CALL zdf_init                             ! Vertical ocean physics
(  376) 
(  377)       !                                     ! Ocean trends
(  378)       ! Control parameters 
(  379)       IF( lk_trdtra .OR. lk_trdmld )   l_trdtra = .TRUE.
(  380)       IF( lk_trddyn .OR. lk_trdvor )   l_trddyn = .TRUE.
(  381) 
(  382)       IF( lk_trddyn .OR. lk_trdtra )   &
(  383)          &            CALL trd_icp_init         ! active tracers and/or momentum
(  384) 
(  385)       IF( lk_trdmld ) CALL trd_mld_init         ! mixed layer
(  386) 
(  387)       IF( lk_trdvor ) CALL trd_vor_init         ! vorticity
(  388) 
(  389) 
(  390) 
(  391) # 399
(  399)       CALL flx_init                         ! Thermohaline forcing initialization
(  400) 
(  401)       CALL flx_fwb_init                     ! FreshWater Budget correction
(  402) 
(  403)       CALL dia_ptr_init                     ! Poleward TRansports initialization
(  404) 
(  405)       !                                     ! =============== !
(  406)       !                                     !  time stepping  !
(  407)       !                                     ! =============== !
(  408) 
(  409) 
(  410) !!DB: 2007.12.06
(  411) !!Code for additional control parameters
(  412) !!Look for a file called run_params and get run params
(  413) !!This is a HARDWIRED list. At this time, the default values
(  414) !!should be found at the bottom of oce.F90
(  415) !!Exit gracefully if file not found
(  416)       if(lwp)write(numout2,*)'DBG Default perpetual_forcing,M2_ave,ioutput_ave: ', &
(  417)            perpetual_forcing,M2_ave,ioutput_ave
(  418)       open(66,file='run_params',STATUS='OLD',ERR=666)
(  419)       read(66,*,ERR=666)perpetual_forcing
(  420)       read(66,*,ERR=666)M2_ave
(  421)       read(66,*,ERR=666)ioutput_ave






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 16

(  422)       close(66)
(  423) !!DBG
(  424)       if(lwp) then
(  425)          write(numout2,*)'DBG: run_params --> perpetual_forcing,M2_ave,ioutput_ave: ', &
(  426)               perpetual_forcing,M2_ave, ioutput_ave
(  427)       endif
(  428) 666   continue
(  429) 
(  430) !!DB 2009.10.01
(  431)       if(ioutput_ave /= 0) then
(  432)          ioutput_ave = nwrite
(  433)          if(lwp) write(numout2,*)'DBG: ioutput_ave /= 0 ; setting to nwrite  = ', nwrite
(  434)       endif
(  435) 
(  436) 
(  437)       IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA
(  438) 
(  439)       IF( lk_cfg_1d  )  THEN 
(  440)          CALL init_1d
(  441)       ENDIF
(  442)    END SUBROUTINE opa_init
(  443)    !!======================================================================
(  444) END MODULE opa
