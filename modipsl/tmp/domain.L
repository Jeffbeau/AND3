


PGF90 (Version     10.5)          09/23/2010  09:49:34      page 1

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

Filename: domain.F90

(    1) MODULE domain
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domain   ***
(    4)    !! Ocean initialization : domain initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_init       : initialize the space and time domain
(    9)    !!   dom_nam        : read and contral domain namelists
(   10)    !!   dom_ctl        : control print for the ocean domain
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! 
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE ice_oce         ! ice variables
(   18)    USE blk_oce         ! bulk variables
(   19)    USE flxrnf          ! runoffs
(   20)    USE daymod          ! calendar
(   21)    USE lib_mpp         ! distributed memory computing library
(   22) 
(   23)    USE domhgr          ! domain: set the horizontal mesh
(   24)    USE domzgr          ! domain: set the vertical mesh
(   25)    USE domstp          ! domain: set the time-step
(   26)    USE dommsk          ! domain: set the mask system
(   27)    USE domwri          ! domain: write the meshmask file
(   28)    USE closea          ! closed sea or lake              (dom_clo routine)
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Routine accessibility
(   34)    PUBLIC dom_init       ! called by opa.F90
(   35) 
(   36)    !! * Substitutions
(   37) #  include "domzgr_substitute.h90"
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domain.F90,v 1.17 2006/03/10 10:55:38 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE dom_init






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !!                  ***  ROUTINE dom_init  ***
(   49)       !!                    
(   50)       !! ** Purpose :   Domain initialization. Call the routines that are 
(   51)       !!      required to create the arrays which define the space and time
(   52)       !!      domain of the ocean model.
(   53)       !!
(   54)       !! ** Method  :
(   55)       !!      - dom_msk: compute the masks from the bathymetry file
(   56)       !!      - dom_hgr: compute or read the horizontal grid-point position and
(   57)       !!                scale factors, and the coriolis factor
(   58)       !!      - dom_zgr: define the vertical coordinate system and the bathymetry
(   59)       !!      - dom_stp: defined the model time step
(   60)       !!      - dom_wri: create the meshmask file if nmsh=1
(   61)       !!
(   62)       !! History :
(   63)       !!        !  90-10  (C. Levy - G. Madec)  Original code
(   64)       !!        !  91-11  (G. Madec)
(   65)       !!        !  92-01  (M. Imbard) insert time step initialization
(   66)       !!        !  96-06  (G. Madec) generalized vertical coordinate 
(   67)       !!        !  97-02  (G. Madec) creation of domwri.F
(   68)       !!        !  01-05  (E.Durand - G. Madec) insert closed sea
(   69)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   70)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   71)       !!----------------------------------------------------------------------
(   72)       !! * Local declarations
(   73)       INTEGER ::   jk                ! dummy loop argument
(   74)       INTEGER ::   iconf = 0         ! temporary integers
(   75)       !!----------------------------------------------------------------------
(   76) 
(   77)       IF(lwp) THEN
(   78)          WRITE(numout,*)
(   79)          WRITE(numout,*) 'dom_init : domain initialization'
(   80)          WRITE(numout,*) '~~~~~~~~'
(   81)       ENDIF
(   82) 
(   83)       CALL dom_nam                        ! read namelist ( namrun, namdom, namcla )
(   84) 
(   85)       CALL dom_clo                        ! Closed seas and lake
(   86) 
(   87)       CALL dom_hgr                        ! Horizontal mesh
(   88) 
(   89)       CALL dom_zgr                        ! Vertical mesh and bathymetry
(   90) 
(   91)       CALL dom_msk                        ! Masks
(   92) 
(   93) 
(   94)       ! Local depth or Inverse of the local depth of the water column at u- and v-points
(   95)       ! ------------------------------
(   96)       ! Ocean depth at U- and V-points
(   97)       hu(:,:) = 0.
(   98)       hv(:,:) = 0.
(   99) 
(  100)       DO jk = 1, jpk
(  101)          hu(:,:) = hu(:,:) + fse3u(:,:,jk) * umask(:,:,jk)
(  102)          hv(:,:) = hv(:,:) + fse3v(:,:,jk) * vmask(:,:,jk)
(  103)       END DO
(  104)       ! Inverse of the local depth






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 3

(  105)       hur(:,:) = fse3u(:,:,1)             ! Lower bound : thickness of the first model level
(  106)       hvr(:,:) = fse3v(:,:,1)
(  107)       
(  108)       DO jk = 2, jpk                      ! Sum of the vertical scale factors
(  109)          hur(:,:) = hur(:,:) + fse3u(:,:,jk) * umask(:,:,jk)
(  110)          hvr(:,:) = hvr(:,:) + fse3v(:,:,jk) * vmask(:,:,jk)
(  111)       END DO
(  112) 
(  113)       ! Compute and mask the inverse of the local depth
(  114)       hur(:,:) = 1. / hur(:,:) * umask(:,:,1)
(  115)       hvr(:,:) = 1. / hvr(:,:) * vmask(:,:,1)
(  116) 
(  117) 
(  118)       CALL dom_stp                        ! Time step
(  119) 
(  120)       IF( nmsh /= 0 )   CALL dom_wri      ! Create a domain file
(  121) 
(  122)       IF( .NOT.ln_rstart )   CALL dom_ctl    ! Domain control
(  123) 
(  124)    END SUBROUTINE dom_init
(  125) 
(  126) 
(  127)    SUBROUTINE dom_nam
(  128)       !!----------------------------------------------------------------------
(  129)       !!                     ***  ROUTINE dom_nam  ***
(  130)       !!                    
(  131)       !! ** Purpose :   read domaine namelists and print the variables.
(  132)       !!
(  133)       !! ** input   : - namrun namelist
(  134)       !!              - namdom namelist
(  135)       !!              - namcla namelist
(  136)       !!
(  137)       !! History :
(  138)       !!   9.0  !  03-08  (G. Madec)  Original code
(  139)       !!----------------------------------------------------------------------
(  140)       !! * Modules used
(  141)       USE ioipsl
(  142)       NAMELIST/namrun/ no    , cexper   , ln_rstart , nrstdt , nit000,         &
(  143)          &             nitend, ndate0   , nleapy   , ninist , nstock,          &
(  144)          &             nprint, nwrite   , nrunoff  , ln_ctl , nictls, nictle,   &
(  145)          &             njctls, njctle   , nbench   , isplt  , jsplt
(  146) 
(  147)       NAMELIST/namdom/ ntopo , e3zps_min, e3zps_rat, ngrid  , nmsh  ,   &
(  148)          &             nacc  , atfp     , rdt      , rdtmin , rdtmax,   &
(  149)          &             rdth  , rdtbt    , nfice    , nfbulk , nclosea
(  150)       NAMELIST/namcla/ n_cla
(  151)       !!----------------------------------------------------------------------
(  152) 
(  153)       IF(lwp) THEN
(  154)          WRITE(numout,*)
(  155)          WRITE(numout,*) 'dom_nam  : domain initialization through namelist read'
(  156)          WRITE(numout,*) '~~~~~~~ '
(  157)       ENDIF
(  158) 
(  159)       ! Namelist namrun : parameters of the run
(  160)       REWIND( numnam )
(  161)       READ  ( numnam, namrun )
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 4

(  163)       IF(lwp) THEN
(  164)          WRITE(numout,*) '        Namelist namrun'
(  165)          WRITE(numout,*) '           job number                      no        = ', no
(  166)          WRITE(numout,*) '           experiment name for output      cexper    = ', cexper
(  167)          WRITE(numout,*) '           restart logical                 ln_rstart = ', ln_rstart
(  168)          WRITE(numout,*) '           control of time step            nrstdt    = ', nrstdt
(  169)          WRITE(numout,*) '           number of the first time step   nit000    = ', nit000
(  170)          WRITE(numout,*) '           number of the last time step    nitend    = ', nitend
(  171)          WRITE(numout,*) '           initial calendar date aammjj    ndate0    = ', ndate0
(  172)          WRITE(numout,*) '           leap year calendar (0/1)        nleapy    = ', nleapy
(  173)          WRITE(numout,*) '           initial state output            ninist    = ', ninist
(  174)          WRITE(numout,*) '           level of print                  nprint    = ', nprint
(  175)          WRITE(numout,*) '           frequency of restart file       nstock    = ', nstock
(  176)          WRITE(numout,*) '           frequency of output file        nwrite    = ', nwrite
(  177)          WRITE(numout,*) '           runoff option                   nrunoff   = ', nrunoff
(  178)          WRITE(numout,*) '           run control (for debugging)     ln_ctl    = ', ln_ctl
(  179)          WRITE(numout,*) '           Start i indice for SUM control  nictls    = ', nictls
(  180)          WRITE(numout,*) '           End i indice for SUM control    nictle    = ', nictle
(  181)          WRITE(numout,*) '           Start j indice for SUM control  njctls    = ', njctls
(  182)          WRITE(numout,*) '           End j indice for SUM control    njctle    = ', njctle
(  183)          WRITE(numout,*) '           number of proc. following i     isplt     = ', isplt
(  184)          WRITE(numout,*) '           number of proc. following j     jsplt     = ', jsplt
(  185)          WRITE(numout,*) '           benchmark parameter (0/1)       nbench    = ', nbench
(  186)       ENDIF
(  187) 
(  188)       ndastp = ndate0                ! Assign initial date to current date
(  189) 
(  190) ! ... Control the sub-domain area indices for the print control
(  191)       IF(ln_ctl)   THEN
(  192)          IF( lk_mpp ) THEN
(  193)             ! the domain is forced to the real splitted domain in MPI
(  194)             isplt = jpni ; jsplt = jpnj ; ijsplt = jpni*jpnj
(  195)          ELSE
(  196)             IF( isplt == 1 .AND. jsplt == 1  ) THEN
(  197)                IF(lwp) WRITE(numout,cform_war)
(  198)                IF(lwp) WRITE(numout,*)'          - isplt & jsplt are equal to 1'
(  199)                IF(lwp) WRITE(numout,*)'          - the print control will be done over the whole domain'
(  200)                IF(lwp) WRITE(numout,*)
(  201)             ENDIF
(  202) 
(  203)             ! compute the total number of processors ijsplt
(  204)             ijsplt = isplt*jsplt
(  205)          ENDIF
(  206) 
(  207)          IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
(  208)          IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
(  209) 
(  210)          ! Control the indices used for the SUM control
(  211)          IF( nictls+nictle+njctls+njctle == 0 )   THEN
(  212)             ! the print control is done over the default area
(  213)             lsp_area = .FALSE.
(  214)          ELSE
(  215)             ! the print control is done over a specific  area
(  216)             lsp_area = .TRUE.
(  217)             IF( nictls < 1 .OR. nictls > jpiglo )   THEN
(  218)                IF(lwp) WRITE(numout,cform_war)
(  219)                IF(lwp) WRITE(numout,*)'          - nictls must be 1<=nictls>=jpiglo, it is forced to 1'
(  220)                IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 5

(  221)                nwarn = nwarn + 1
(  222)                nictls = 1
(  223)             ENDIF
(  224) 
(  225)             IF( nictle < 1 .OR. nictle > jpiglo )   THEN
(  226)                IF(lwp) WRITE(numout,cform_war)
(  227)                IF(lwp) WRITE(numout,*)'          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo'
(  228)                IF(lwp) WRITE(numout,*)
(  229)                nwarn = nwarn + 1
(  230)                nictle = jpjglo
(  231)             ENDIF
(  232) 
(  233)             IF( njctls < 1 .OR. njctls > jpjglo )   THEN
(  234)                IF(lwp) WRITE(numout,cform_war)
(  235)                IF(lwp) WRITE(numout,*)'          - njctls must be 1<=njctls>=jpjglo, it is forced to 1'
(  236)                IF(lwp) WRITE(numout,*)
(  237)                nwarn = nwarn + 1
(  238)                njctls = 1
(  239)             ENDIF
(  240) 
(  241)             IF( njctle < 1 .OR. njctle > jpjglo )   THEN
(  242)                IF(lwp) WRITE(numout,cform_war)
(  243)                IF(lwp) WRITE(numout,*)'          - njctle must be 1<=njctle>= jpjglo, it is forced to jpjglo'
(  244)                IF(lwp) WRITE(numout,*)
(  245)                nwarn = nwarn + 1
(  246)                njctle = jpjglo
(  247)             ENDIF
(  248) 
(  249)          ENDIF          ! IF( nictls+nictle+njctls+njctle == 0 )
(  250)        ENDIF            ! IF(ln_ctl)
(  251) 
(  252) ! ... Control of output frequency
(  253)       IF ( nstock == 0 ) THEN
(  254)           IF(lwp)WRITE(numout,cform_war)
(  255)           IF(lwp)WRITE(numout,*) '           nstock = ', nstock, ' it is forced to ', nitend
(  256)           nstock = nitend
(  257)           nwarn = nwarn + 1
(  258)       ENDIF
(  259)       IF ( nwrite == 0 ) THEN
(  260)           IF(lwp)WRITE(numout,cform_war)
(  261)           IF(lwp)WRITE(numout,*) '           nwrite = ', nwrite, ' it is forced to ', nitend
(  262)           nwrite = nitend
(  263)           nwarn = nwarn + 1
(  264)       ENDIF
(  265) 
(  266) #if defined key_agrif
(  267)       if ( Agrif_Root() ) then
(  268) #endif
(  269)       SELECT CASE ( nleapy )   ! Choose calendar for IOIPSL
(  270)       CASE (  1 ) 
(  271)          CALL ioconf_calendar('gregorian')
(  272)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "gregorian", i.e. leap year'
(  273)       CASE (  0 )
(  274)          CALL ioconf_calendar('noleap')
(  275)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "noleap", i.e. no leap year'
(  276)       CASE ( 30 )
(  277)          CALL ioconf_calendar('360d')
(  278)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "360d", i.e. 360 days in a year'






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 6

(  279)       END SELECT
(  280) #if defined key_agrif
(  281)       endif
(  282) #endif
(  283) 
(  284)       SELECT CASE ( nleapy )   ! year=raajj*days day=rjjhh*hours hour=rhhmm*minutes etc ...
(  285)       CASE ( 1 )
(  286)          raajj = 365.25
(  287)          raass = raajj * rjjss
(  288)          rmoss = raass/raamo
(  289)       CASE ( 0 )
(  290)          raajj = 365.
(  291)          raass = raajj * rjjss
(  292)          rmoss = raass/raamo
(  293)       CASE DEFAULT
(  294)          raajj = FLOAT( nleapy ) * raamo
(  295)          raass =        raajj    * rjjss
(  296)          rmoss = FLOAT( nleapy ) * rjjss
(  297)       END SELECT
(  298)       IF(lwp) THEN
(  299)          WRITE(numout,*)
(  300)          WRITE(numout,*) '           nb of days per year      raajj = ', raajj,' days'
(  301)          WRITE(numout,*) '           nb of seconds per year   raass = ', raass, ' s'
(  302)          WRITE(numout,*) '           nb of seconds per month  rmoss = ', rmoss, ' s'
(  303)       ENDIF
(  304) 
(  305)       ! Namelist namdom : space/time domain (bathymetry, mesh, timestep)
(  306)       REWIND( numnam )
(  307)       READ  ( numnam, namdom )
(  308) 
(  309)       IF(lwp) THEN
(  310)          WRITE(numout,*)
(  311)          WRITE(numout,*) '        Namelist namdom'
(  312)          WRITE(numout,*) '           flag read/compute bathymetry   ntopo     = ', ntopo
(  313)          WRITE(numout,*) '           minimum thickness of partial   e3zps_min = ', e3zps_min, ' (m)'
(  314)          WRITE(numout,*) '              step level                  e3zps_rat = ', e3zps_rat
(  315)          WRITE(numout,*) '           flag read/compute coordinates  ngrid     = ', ngrid
(  316)          WRITE(numout,*) '           flag write mesh/mask file(s)   nmsh      = ', nmsh
(  317)          WRITE(numout,*) '                = 0   no file created                 '
(  318)          WRITE(numout,*) '                = 1   mesh_mask                       '
(  319)          WRITE(numout,*) '                = 2   mesh and mask                   '
(  320)          WRITE(numout,*) '                = 3   mesh_hgr, msh_zgr and mask      '
(  321)          WRITE(numout,*) '           acceleration of converge       nacc      = ', nacc
(  322)          WRITE(numout,*) '           asselin time filter parameter  atfp      = ', atfp
(  323)          WRITE(numout,*) '           time step                      rdt       = ', rdt
(  324)          WRITE(numout,*) '           minimum time step on tracers   rdtmin    = ', rdtmin
(  325)          WRITE(numout,*) '           maximum time step on tracers   rdtmax    = ', rdtmax
(  326)          WRITE(numout,*) '           depth variation tracer step    rdth      = ', rdth
(  327)          WRITE(numout,*) '           barotropic time step           rdtbt     = ', rdtbt
(  328)       ENDIF
(  329) 
(  330)       IF( lk_ice_lim ) THEN
(  331)          IF(lwp) WRITE(numout,*) '           ice model coupling frequency      nfice  = ', nfice
(  332)          nfbulk = nfice
(  333)          IF( MOD( rday, nfice*rdt ) /= 0 ) THEN 
(  334)             IF(lwp) WRITE(numout,*) ' '
(  335)             IF(lwp) WRITE(numout,*) 'W A R N I N G :  nfice is NOT a multiple of the number of time steps in a day'
(  336)             IF(lwp) WRITE(numout,*) ' '






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 7

(  337)          ENDIF
(  338)          IF(lwp) WRITE(numout,*) '           bulk computation frequency       nfbulk  = ', nfbulk, ' = nfice if ice model us
(  339)          IF(lwp) WRITE(numout,*) '           flag closed sea or not           nclosea = ', nclosea
(  340)       ENDIF
(  341) 
(  342)       ! Default values
(  343)       n_cla = 0
(  344) 
(  345)       ! Namelist cross land advection
(  346)       REWIND( numnam )
(  347)       READ  ( numnam, namcla )
(  348)       IF(lwp) THEN
(  349)          WRITE(numout,*)
(  350)          WRITE(numout,*) '        Namelist namcla'
(  351)          WRITE(numout,*) '           cross land advection        n_cla        = ',n_cla
(  352)       ENDIF
(  353) 
(  354)    END SUBROUTINE dom_nam
(  355) 
(  356) 
(  357)    SUBROUTINE dom_ctl
(  358)       !!----------------------------------------------------------------------
(  359)       !!                     ***  ROUTINE dom_ctl  ***
(  360)       !!
(  361)       !! ** Purpose :   Domain control.
(  362)       !!
(  363)       !! ** Method  :   compute and print extrema of masked scale factors
(  364)       !!
(  365)       !! History :
(  366)       !!   8.5  !  02-08  (G. Madec)    Original code
(  367)       !!----------------------------------------------------------------------
(  368)       !! * Local declarations
(  369)       INTEGER ::   iimi1, ijmi1, iimi2, ijmi2, iima1, ijma1, iima2, ijma2
(  370)       INTEGER, DIMENSION(2) ::   iloc      ! 
(  371)       REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
(  372)       !!----------------------------------------------------------------------
(  373) 
(  374)       ! Extrema of the scale factors
(  375) 
(  376)       IF(lwp)WRITE(numout,*)
(  377)       IF(lwp)WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
(  378)       IF(lwp)WRITE(numout,*) '~~~~~~~'
(  379) 
(  380)       IF (lk_mpp) THEN
(  381)          CALL mpp_minloc( e1t(:,:), tmask(:,:,1), ze1min, iimi1,ijmi1 )
(  382)          CALL mpp_minloc( e2t(:,:), tmask(:,:,1), ze2min, iimi2,ijmi2 )
(  383)          CALL mpp_maxloc( e1t(:,:), tmask(:,:,1), ze1max, iima1,ijma1 )
(  384)          CALL mpp_maxloc( e2t(:,:), tmask(:,:,1), ze2max, iima2,ijma2 )
(  385)       ELSE
(  386)          ze1min = MINVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  387)          ze2min = MINVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  388)          ze1max = MAXVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  389)          ze2max = MAXVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  390) 
(  391)          iloc  = MINLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  392)          iimi1 = iloc(1) + nimpp - 1
(  393)          ijmi1 = iloc(2) + njmpp - 1
(  394)          iloc  = MINLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 8

(  395)          iimi2 = iloc(1) + nimpp - 1
(  396)          ijmi2 = iloc(2) + njmpp - 1
(  397)          iloc  = MAXLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  398)          iima1 = iloc(1) + nimpp - 1
(  399)          ijma1 = iloc(2) + njmpp - 1
(  400)          iloc  = MAXLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  401)          iima2 = iloc(1) + nimpp - 1
(  402)          ijma2 = iloc(2) + njmpp - 1
(  403)       ENDIF
(  404) 
(  405)       IF(lwp) THEN
(  406)          WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, iima1, ijma1
(  407)          WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, iimi1, ijmi1
(  408)          WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, iima2, ijma2
(  409)          WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, iimi2, ijmi2
(  410)       ENDIF
(  411) 
(  412)    END SUBROUTINE dom_ctl
(  413) 
(  414)    !!======================================================================
(  415) END MODULE domain











































PGF90 (Version     10.5)          09/23/2010  09:49:34      page 9

(    1) # 1 "domain.F90"
(    1) MODULE domain
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domain   ***
(    4)    !! Ocean initialization : domain initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_init       : initialize the space and time domain
(    9)    !!   dom_nam        : read and contral domain namelists
(   10)    !!   dom_ctl        : control print for the ocean domain
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! 
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE ice_oce         ! ice variables
(   18)    USE blk_oce         ! bulk variables
(   19)    USE flxrnf          ! runoffs
(   20)    USE daymod          ! calendar
(   21)    USE lib_mpp         ! distributed memory computing library
(   22) 
(   23)    USE domhgr          ! domain: set the horizontal mesh
(   24)    USE domzgr          ! domain: set the vertical mesh
(   25)    USE domstp          ! domain: set the time-step
(   26)    USE dommsk          ! domain: set the mask system
(   27)    USE domwri          ! domain: write the meshmask file
(   28)    USE closea          ! closed sea or lake              (dom_clo routine)
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Routine accessibility
(   34)    PUBLIC dom_init       ! called by opa.F90
(   35) 
(   36)    !! * Substitutions
(   37) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 38 "domain.F90"
(   38) # 38
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domain.F90,v 1.17 2006/03/10 10:55:38 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 10

(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE dom_init
(   47)       !!----------------------------------------------------------------------
(   48)       !!                  ***  ROUTINE dom_init  ***
(   49)       !!                    
(   50)       !! ** Purpose :   Domain initialization. Call the routines that are 
(   51)       !!      required to create the arrays which define the space and time
(   52)       !!      domain of the ocean model.
(   53)       !!
(   54)       !! ** Method  :
(   55)       !!      - dom_msk: compute the masks from the bathymetry file
(   56)       !!      - dom_hgr: compute or read the horizontal grid-point position and
(   57)       !!                scale factors, and the coriolis factor
(   58)       !!      - dom_zgr: define the vertical coordinate system and the bathymetry
(   59)       !!      - dom_stp: defined the model time step
(   60)       !!      - dom_wri: create the meshmask file if nmsh=1
(   61)       !!
(   62)       !! History :
(   63)       !!        !  90-10  (C. Levy - G. Madec)  Original code
(   64)       !!        !  91-11  (G. Madec)
(   65)       !!        !  92-01  (M. Imbard) insert time step initialization
(   66)       !!        !  96-06  (G. Madec) generalized vertical coordinate 
(   67)       !!        !  97-02  (G. Madec) creation of domwri.F
(   68)       !!        !  01-05  (E.Durand - G. Madec) insert closed sea
(   69)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   70)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   71)       !!----------------------------------------------------------------------
(   72)       !! * Local declarations
(   73)       INTEGER ::   jk                ! dummy loop argument
(   74)       INTEGER ::   iconf = 0         ! temporary integers
(   75)       !!----------------------------------------------------------------------
(   76) 
(   77)       IF(lwp) THEN
(   78)          WRITE(numout,*)
(   79)          WRITE(numout,*) 'dom_init : domain initialization'
(   80)          WRITE(numout,*) '~~~~~~~~'
(   81)       ENDIF
(   82) 
(   83)       CALL dom_nam                        ! read namelist ( namrun, namdom, namcla )
(   84) 
(   85)       CALL dom_clo                        ! Closed seas and lake
(   86) 
(   87)       CALL dom_hgr                        ! Horizontal mesh
(   88) 
(   89)       CALL dom_zgr                        ! Vertical mesh and bathymetry
(   90) 
(   91)       CALL dom_msk                        ! Masks
(   92) 
(   93) 
(   94)       ! Local depth or Inverse of the local depth of the water column at u- and v-points
(   95)       ! ------------------------------
(   96)       ! Ocean depth at U- and V-points
(   97)       hu(:,:) = 0.
(   98)       hv(:,:) = 0.
(   99) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 11

(  100)       DO jk = 1, jpk
(  101)          hu(:,:) = hu(:,:) + e3t(jk) * umask(:,:,jk)
(  102)          hv(:,:) = hv(:,:) + e3t(jk) * vmask(:,:,jk)
(  103)       END DO
(  104)       ! Inverse of the local depth
(  105)       hur(:,:) = e3t(1)             ! Lower bound : thickness of the first model level
(  106)       hvr(:,:) = e3t(1)
(  107)       
(  108)       DO jk = 2, jpk                      ! Sum of the vertical scale factors
(  109)          hur(:,:) = hur(:,:) + e3t(jk) * umask(:,:,jk)
(  110)          hvr(:,:) = hvr(:,:) + e3t(jk) * vmask(:,:,jk)
(  111)       END DO
(  112) 
(  113)       ! Compute and mask the inverse of the local depth
(  114)       hur(:,:) = 1. / hur(:,:) * umask(:,:,1)
(  115)       hvr(:,:) = 1. / hvr(:,:) * vmask(:,:,1)
(  116) 
(  117) 
(  118)       CALL dom_stp                        ! Time step
(  119) 
(  120)       IF( nmsh /= 0 )   CALL dom_wri      ! Create a domain file
(  121) 
(  122)       IF( .NOT.ln_rstart )   CALL dom_ctl    ! Domain control
(  123) 
(  124)    END SUBROUTINE dom_init
(  125) 
(  126) 
(  127)    SUBROUTINE dom_nam
(  128)       !!----------------------------------------------------------------------
(  129)       !!                     ***  ROUTINE dom_nam  ***
(  130)       !!                    
(  131)       !! ** Purpose :   read domaine namelists and print the variables.
(  132)       !!
(  133)       !! ** input   : - namrun namelist
(  134)       !!              - namdom namelist
(  135)       !!              - namcla namelist
(  136)       !!
(  137)       !! History :
(  138)       !!   9.0  !  03-08  (G. Madec)  Original code
(  139)       !!----------------------------------------------------------------------
(  140)       !! * Modules used
(  141)       USE ioipsl
(  142)       NAMELIST/namrun/ no    , cexper   , ln_rstart , nrstdt , nit000,         &
(  143)          &             nitend, ndate0   , nleapy   , ninist , nstock,          &
(  144)          &             nprint, nwrite   , nrunoff  , ln_ctl , nictls, nictle,   &
(  145)          &             njctls, njctle   , nbench   , isplt  , jsplt
(  146) 
(  147)       NAMELIST/namdom/ ntopo , e3zps_min, e3zps_rat, ngrid  , nmsh  ,   &
(  148)          &             nacc  , atfp     , rdt      , rdtmin , rdtmax,   &
(  149)          &             rdth  , rdtbt    , nfice    , nfbulk , nclosea
(  150)       NAMELIST/namcla/ n_cla
(  151)       !!----------------------------------------------------------------------
(  152) 
(  153)       IF(lwp) THEN
(  154)          WRITE(numout,*)
(  155)          WRITE(numout,*) 'dom_nam  : domain initialization through namelist read'
(  156)          WRITE(numout,*) '~~~~~~~ '
(  157)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 12

(  158) 
(  159)       ! Namelist namrun : parameters of the run
(  160)       REWIND( numnam )
(  161)       READ  ( numnam, namrun )
(  162) 
(  163)       IF(lwp) THEN
(  164)          WRITE(numout,*) '        Namelist namrun'
(  165)          WRITE(numout,*) '           job number                      no        = ', no
(  166)          WRITE(numout,*) '           experiment name for output      cexper    = ', cexper
(  167)          WRITE(numout,*) '           restart logical                 ln_rstart = ', ln_rstart
(  168)          WRITE(numout,*) '           control of time step            nrstdt    = ', nrstdt
(  169)          WRITE(numout,*) '           number of the first time step   nit000    = ', nit000
(  170)          WRITE(numout,*) '           number of the last time step    nitend    = ', nitend
(  171)          WRITE(numout,*) '           initial calendar date aammjj    ndate0    = ', ndate0
(  172)          WRITE(numout,*) '           leap year calendar (0/1)        nleapy    = ', nleapy
(  173)          WRITE(numout,*) '           initial state output            ninist    = ', ninist
(  174)          WRITE(numout,*) '           level of print                  nprint    = ', nprint
(  175)          WRITE(numout,*) '           frequency of restart file       nstock    = ', nstock
(  176)          WRITE(numout,*) '           frequency of output file        nwrite    = ', nwrite
(  177)          WRITE(numout,*) '           runoff option                   nrunoff   = ', nrunoff
(  178)          WRITE(numout,*) '           run control (for debugging)     ln_ctl    = ', ln_ctl
(  179)          WRITE(numout,*) '           Start i indice for SUM control  nictls    = ', nictls
(  180)          WRITE(numout,*) '           End i indice for SUM control    nictle    = ', nictle
(  181)          WRITE(numout,*) '           Start j indice for SUM control  njctls    = ', njctls
(  182)          WRITE(numout,*) '           End j indice for SUM control    njctle    = ', njctle
(  183)          WRITE(numout,*) '           number of proc. following i     isplt     = ', isplt
(  184)          WRITE(numout,*) '           number of proc. following j     jsplt     = ', jsplt
(  185)          WRITE(numout,*) '           benchmark parameter (0/1)       nbench    = ', nbench
(  186)       ENDIF
(  187) 
(  188)       ndastp = ndate0                ! Assign initial date to current date
(  189) 
(  190) ! ... Control the sub-domain area indices for the print control
(  191)       IF(ln_ctl)   THEN
(  192)          IF( lk_mpp ) THEN
(  193)             ! the domain is forced to the real splitted domain in MPI
(  194)             isplt = jpni ; jsplt = jpnj ; ijsplt = jpni*jpnj
(  195)          ELSE
(  196)             IF( isplt == 1 .AND. jsplt == 1  ) THEN
(  197)                IF(lwp) WRITE(numout,cform_war)
(  198)                IF(lwp) WRITE(numout,*)'          - isplt & jsplt are equal to 1'
(  199)                IF(lwp) WRITE(numout,*)'          - the print control will be done over the whole domain'
(  200)                IF(lwp) WRITE(numout,*)
(  201)             ENDIF
(  202) 
(  203)             ! compute the total number of processors ijsplt
(  204)             ijsplt = isplt*jsplt
(  205)          ENDIF
(  206) 
(  207)          IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
(  208)          IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt
(  209) 
(  210)          ! Control the indices used for the SUM control
(  211)          IF( nictls+nictle+njctls+njctle == 0 )   THEN
(  212)             ! the print control is done over the default area
(  213)             lsp_area = .FALSE.
(  214)          ELSE
(  215)             ! the print control is done over a specific  area






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 13

(  216)             lsp_area = .TRUE.
(  217)             IF( nictls < 1 .OR. nictls > jpiglo )   THEN
(  218)                IF(lwp) WRITE(numout,cform_war)
(  219)                IF(lwp) WRITE(numout,*)'          - nictls must be 1<=nictls>=jpiglo, it is forced to 1'
(  220)                IF(lwp) WRITE(numout,*)
(  221)                nwarn = nwarn + 1
(  222)                nictls = 1
(  223)             ENDIF
(  224) 
(  225)             IF( nictle < 1 .OR. nictle > jpiglo )   THEN
(  226)                IF(lwp) WRITE(numout,cform_war)
(  227)                IF(lwp) WRITE(numout,*)'          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo'
(  228)                IF(lwp) WRITE(numout,*)
(  229)                nwarn = nwarn + 1
(  230)                nictle = jpjglo
(  231)             ENDIF
(  232) 
(  233)             IF( njctls < 1 .OR. njctls > jpjglo )   THEN
(  234)                IF(lwp) WRITE(numout,cform_war)
(  235)                IF(lwp) WRITE(numout,*)'          - njctls must be 1<=njctls>=jpjglo, it is forced to 1'
(  236)                IF(lwp) WRITE(numout,*)
(  237)                nwarn = nwarn + 1
(  238)                njctls = 1
(  239)             ENDIF
(  240) 
(  241)             IF( njctle < 1 .OR. njctle > jpjglo )   THEN
(  242)                IF(lwp) WRITE(numout,cform_war)
(  243)                IF(lwp) WRITE(numout,*)'          - njctle must be 1<=njctle>= jpjglo, it is forced to jpjglo'
(  244)                IF(lwp) WRITE(numout,*)
(  245)                nwarn = nwarn + 1
(  246)                njctle = jpjglo
(  247)             ENDIF
(  248) 
(  249)          ENDIF          ! IF( nictls+nictle+njctls+njctle == 0 )
(  250)        ENDIF            ! IF(ln_ctl)
(  251) 
(  252) ! ... Control of output frequency
(  253)       IF ( nstock == 0 ) THEN
(  254)           IF(lwp)WRITE(numout,cform_war)
(  255)           IF(lwp)WRITE(numout,*) '           nstock = ', nstock, ' it is forced to ', nitend
(  256)           nstock = nitend
(  257)           nwarn = nwarn + 1
(  258)       ENDIF
(  259)       IF ( nwrite == 0 ) THEN
(  260)           IF(lwp)WRITE(numout,cform_war)
(  261)           IF(lwp)WRITE(numout,*) '           nwrite = ', nwrite, ' it is forced to ', nitend
(  262)           nwrite = nitend
(  263)           nwarn = nwarn + 1
(  264)       ENDIF
(  265) 
(  266) # 269
(  269)       SELECT CASE ( nleapy )   ! Choose calendar for IOIPSL
(  270)       CASE (  1 ) 
(  271)          CALL ioconf_calendar('gregorian')
(  272)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "gregorian", i.e. leap year'
(  273)       CASE (  0 )
(  274)          CALL ioconf_calendar('noleap')
(  275)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "noleap", i.e. no leap year'






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 14

(  276)       CASE ( 30 )
(  277)          CALL ioconf_calendar('360d')
(  278)          IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "360d", i.e. 360 days in a year'
(  279)       END SELECT
(  280) 
(  281) # 284
(  284)       SELECT CASE ( nleapy )   ! year=raajj*days day=rjjhh*hours hour=rhhmm*minutes etc ...
(  285)       CASE ( 1 )
(  286)          raajj = 365.25
(  287)          raass = raajj * rjjss
(  288)          rmoss = raass/raamo
(  289)       CASE ( 0 )
(  290)          raajj = 365.
(  291)          raass = raajj * rjjss
(  292)          rmoss = raass/raamo
(  293)       CASE DEFAULT
(  294)          raajj = FLOAT( nleapy ) * raamo
(  295)          raass =        raajj    * rjjss
(  296)          rmoss = FLOAT( nleapy ) * rjjss
(  297)       END SELECT
(  298)       IF(lwp) THEN
(  299)          WRITE(numout,*)
(  300)          WRITE(numout,*) '           nb of days per year      raajj = ', raajj,' days'
(  301)          WRITE(numout,*) '           nb of seconds per year   raass = ', raass, ' s'
(  302)          WRITE(numout,*) '           nb of seconds per month  rmoss = ', rmoss, ' s'
(  303)       ENDIF
(  304) 
(  305)       ! Namelist namdom : space/time domain (bathymetry, mesh, timestep)
(  306)       REWIND( numnam )
(  307)       READ  ( numnam, namdom )
(  308) 
(  309)       IF(lwp) THEN
(  310)          WRITE(numout,*)
(  311)          WRITE(numout,*) '        Namelist namdom'
(  312)          WRITE(numout,*) '           flag read/compute bathymetry   ntopo     = ', ntopo
(  313)          WRITE(numout,*) '           minimum thickness of partial   e3zps_min = ', e3zps_min, ' (m)'
(  314)          WRITE(numout,*) '              step level                  e3zps_rat = ', e3zps_rat
(  315)          WRITE(numout,*) '           flag read/compute coordinates  ngrid     = ', ngrid
(  316)          WRITE(numout,*) '           flag write mesh/mask file(s)   nmsh      = ', nmsh
(  317)          WRITE(numout,*) '                = 0   no file created                 '
(  318)          WRITE(numout,*) '                = 1   mesh_mask                       '
(  319)          WRITE(numout,*) '                = 2   mesh and mask                   '
(  320)          WRITE(numout,*) '                = 3   mesh_hgr, msh_zgr and mask      '
(  321)          WRITE(numout,*) '           acceleration of converge       nacc      = ', nacc
(  322)          WRITE(numout,*) '           asselin time filter parameter  atfp      = ', atfp
(  323)          WRITE(numout,*) '           time step                      rdt       = ', rdt
(  324)          WRITE(numout,*) '           minimum time step on tracers   rdtmin    = ', rdtmin
(  325)          WRITE(numout,*) '           maximum time step on tracers   rdtmax    = ', rdtmax
(  326)          WRITE(numout,*) '           depth variation tracer step    rdth      = ', rdth
(  327)          WRITE(numout,*) '           barotropic time step           rdtbt     = ', rdtbt
(  328)       ENDIF
(  329) 
(  330)       IF( lk_ice_lim ) THEN
(  331)          IF(lwp) WRITE(numout,*) '           ice model coupling frequency      nfice  = ', nfice
(  332)          nfbulk = nfice
(  333)          IF( MOD( rday, nfice*rdt ) /= 0 ) THEN 
(  334)             IF(lwp) WRITE(numout,*) ' '
(  335)             IF(lwp) WRITE(numout,*) 'W A R N I N G :  nfice is NOT a multiple of the number of time steps in a day'






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 15

(  336)             IF(lwp) WRITE(numout,*) ' '
(  337)          ENDIF
(  338)          IF(lwp) WRITE(numout,*) '           bulk computation frequency       nfbulk  = ', nfbulk, ' = nfice if ice model used'
(  339)          IF(lwp) WRITE(numout,*) '           flag closed sea or not           nclosea = ', nclosea
(  340)       ENDIF
(  341) 
(  342)       ! Default values
(  343)       n_cla = 0
(  344) 
(  345)       ! Namelist cross land advection
(  346)       REWIND( numnam )
(  347)       READ  ( numnam, namcla )
(  348)       IF(lwp) THEN
(  349)          WRITE(numout,*)
(  350)          WRITE(numout,*) '        Namelist namcla'
(  351)          WRITE(numout,*) '           cross land advection        n_cla        = ',n_cla
(  352)       ENDIF
(  353) 
(  354)    END SUBROUTINE dom_nam
(  355) 
(  356) 
(  357)    SUBROUTINE dom_ctl
(  358)       !!----------------------------------------------------------------------
(  359)       !!                     ***  ROUTINE dom_ctl  ***
(  360)       !!
(  361)       !! ** Purpose :   Domain control.
(  362)       !!
(  363)       !! ** Method  :   compute and print extrema of masked scale factors
(  364)       !!
(  365)       !! History :
(  366)       !!   8.5  !  02-08  (G. Madec)    Original code
(  367)       !!----------------------------------------------------------------------
(  368)       !! * Local declarations
(  369)       INTEGER ::   iimi1, ijmi1, iimi2, ijmi2, iima1, ijma1, iima2, ijma2
(  370)       INTEGER, DIMENSION(2) ::   iloc      ! 
(  371)       REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
(  372)       !!----------------------------------------------------------------------
(  373) 
(  374)       ! Extrema of the scale factors
(  375) 
(  376)       IF(lwp)WRITE(numout,*)
(  377)       IF(lwp)WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
(  378)       IF(lwp)WRITE(numout,*) '~~~~~~~'
(  379) 
(  380)       IF (lk_mpp) THEN
(  381)          CALL mpp_minloc( e1t(:,:), tmask(:,:,1), ze1min, iimi1,ijmi1 )
(  382)          CALL mpp_minloc( e2t(:,:), tmask(:,:,1), ze2min, iimi2,ijmi2 )
(  383)          CALL mpp_maxloc( e1t(:,:), tmask(:,:,1), ze1max, iima1,ijma1 )
(  384)          CALL mpp_maxloc( e2t(:,:), tmask(:,:,1), ze2max, iima2,ijma2 )
(  385)       ELSE
(  386)          ze1min = MINVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  387)          ze2min = MINVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  388)          ze1max = MAXVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  389)          ze2max = MAXVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    
(  390) 
(  391)          iloc  = MINLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  392)          iimi1 = iloc(1) + nimpp - 1
(  393)          ijmi1 = iloc(2) + njmpp - 1






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 16

(  394)          iloc  = MINLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  395)          iimi2 = iloc(1) + nimpp - 1
(  396)          ijmi2 = iloc(2) + njmpp - 1
(  397)          iloc  = MAXLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  398)          iima1 = iloc(1) + nimpp - 1
(  399)          ijma1 = iloc(2) + njmpp - 1
(  400)          iloc  = MAXLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )
(  401)          iima2 = iloc(1) + nimpp - 1
(  402)          ijma2 = iloc(2) + njmpp - 1
(  403)       ENDIF
(  404) 
(  405)       IF(lwp) THEN
(  406)          WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, iima1, ijma1
(  407)          WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, iimi1, ijmi1
(  408)          WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, iima2, ijma2
(  409)          WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, iimi2, ijmi2
(  410)       ENDIF
(  411) 
(  412)    END SUBROUTINE dom_ctl
(  413) 
(  414)    !!======================================================================
(  415) END MODULE domain
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (domain.F90: 294)
