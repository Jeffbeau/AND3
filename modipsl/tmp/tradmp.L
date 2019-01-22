


PGF90 (Version     10.5)          09/23/2010  09:49:13      page 1

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

Filename: tradmp.F90

(    1) MODULE tradmp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  tradmp  ***
(    4)    !! Ocean physics: internal restoring trend on active tracers (T and S)
(    5)    !!======================================================================
(    6) #if   defined key_tradmp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   key_tradmp                                         internal damping
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_dmp      : update the tracer trend with the internal damping
(   11)    !!   tra_dmp_init : initialization, namlist read, parameters control
(   12)    !!   dtacof_zoom  : restoring coefficient for zoom domain
(   13)    !!   dtacof       : restoring coefficient for global domain
(   14)    !!   cofdis       : compute the distance to the coastline
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers variables
(   18)    USE dom_oce         ! ocean space and time domain variables
(   19)    USE trdmod          ! ocean active tracers trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE zdf_oce         ! ocean vertical physics
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE phycst          ! Define parameters for the routines
(   24)    USE dtatem          ! temperature data
(   25)    USE dtasal          ! salinity data
(   26)    USE zdfmxl          ! mixed layer depth
(   27)    USE lib_mpp         ! distribued memory computing
(   28)    USE prtctl          ! Print control
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Routine accessibility
(   34)    PUBLIC tra_dmp   ! routine called by step.F90
(   35) 
(   36)    !! * Shared module variables
(   37)    LOGICAL , PUBLIC &
(   38) #if ! defined key_agrif
(   39)    , PARAMETER  &
(   40) #endif
(   41)    ::   lk_tradmp = .TRUE.    !: internal damping flag
(   42) 
(   43)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &
(   44)       strdmp,              &  ! damping salinity trend (psu/s)
(   45)       resto                   ! restoring coeff. on T and S (s-1)
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 2

(   47)    !! * Module variables
(   48)    INTEGER  ::             & !!! * newtonian damping namelist (mandmp) *
(   49)       ndmp   =   -1 ,      &  ! = 0/-1/'latitude' for damping over T and S
(   50)       ndmpf  =    2 ,      &  ! = 1 create a damping.coeff NetCDF file 
(   51)       nmldmp =    0           ! = 0/1/2 flag for damping in the mixed layer
(   52)    REAL(wp) ::             & !!!  * newtonian damping namelist *
(   53)       sdmp   =   50.,      &  ! surface time scale for internal damping (days)
(   54)       bdmp   =  360.,      &  ! bottom time scale for internal damping (days)
(   55)       hdmp   =  800.          ! depth of transition between sdmp and bdmp (meters)
(   56) 
(   57)    !! * Substitutions
(   58) #  include "domzgr_substitute.h90"
(   59) #  include "vectopt_loop_substitute.h90"
(   60)    !!----------------------------------------------------------------------
(   61)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   62)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tradmp.F90,v 1.15 2006/04/19 14:43:17 opalod Exp $ 
(   63)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   64)    !!----------------------------------------------------------------------
(   65) 
(   66) CONTAINS
(   67) 
(   68)    SUBROUTINE tra_dmp( kt )
(   69)       !!----------------------------------------------------------------------
(   70)       !!                   ***  ROUTINE tra_dmp  ***
(   71)       !!                  
(   72)       !! ** Purpose :   Compute the tracer trend due to a newtonian damping
(   73)       !!      of the tracer field towards given data field and add it to the
(   74)       !!      general tracer trends.
(   75)       !!
(   76)       !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
(   77)       !!      and add to the general tracer trends:
(   78)       !!                     ta = ta + resto * (t_dta - tb)
(   79)       !!                     sa = sa + resto * (s_dta - sb)
(   80)       !!         The trend is computed either throughout the water column
(   81)       !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
(   82)       !!      below the well mixed layer (nlmdmp=2)
(   83)       !!
(   84)       !! ** Action  : - update the tracer trends (ta,sa) with the newtonian 
(   85)       !!                damping trends.
(   86)       !!              - save the trends in (ttrd,strd) ('key_trdtra')
(   87)       !!
(   88)       !! History :
(   89)       !!   7.0  !         (G. Madec)  Original code
(   90)       !!        !  96-01  (G. Madec) 
(   91)       !!        !  97-05  (G. Madec)  macro-tasked on jk-slab
(   92)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(   93)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   94)       !!----------------------------------------------------------------------
(   95)       !! * Modules used     
(   96)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   97)                          ztdsa => va      ! use va as 3D workspace   
(   98) 
(   99)       !! * Arguments
(  100)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(  101) 
(  102)       !! * Local declarations
(  103)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(  104)       REAL(wp) ::   ztest, zta, zsa       ! temporary scalars






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 3

(  105)       !!----------------------------------------------------------------------
(  106) 
(  107)       ! 0. Initialization (first time-step only)
(  108)       !    --------------
(  109)       IF( kt == nit000 ) CALL tra_dmp_init
(  110) 
(  111)       ! Save ta and sa trends
(  112)       IF( l_trdtra )   THEN
(  113)          ztdta(:,:,:) = ta(:,:,:) 
(  114)          ztdsa(:,:,:) = sa(:,:,:) 
(  115)       ENDIF
(  116) 
(  117)       ! 1. Newtonian damping trends on tracer fields
(  118)       ! --------------------------------------------
(  119)       !    compute the newtonian damping trends depending on nmldmp
(  120) 
(  121)       SELECT CASE ( nmldmp )
(  122) 
(  123)       CASE( 0 )                ! newtonian damping throughout the water column
(  124) !!DBG: 2009.06.24 -- HARDWIRE strong sfce S restoring
(  125) !         DO jk = 1, 1
(  126) !            DO jj = 2, jpjm1
(  127) !               DO ji = fs_2, fs_jpim1   ! vector opt.
(  128) !                  zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  129) !!                  zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  130) !                  zsa = (1./(3.*rday)) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  131) !                  ! add the trends to the general tracer trends
(  132) !                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  133) !                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  134) !                  ! save the salinity trend (used in flx to close the salt budget)
(  135) !                  strdmp(ji,jj,jk) = zsa
(  136) !               END DO
(  137) !            END DO
(  138) !         END DO
(  139) !         DO jk = 2, jpkm1
(  140) 
(  141)          DO jk = 1, jpkm1
(  142)             DO jj = 2, jpjm1
(  143)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  144)                   zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  145)                   zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  146)                   ! add the trends to the general tracer trends
(  147)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  148)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  149)                   ! save the salinity trend (used in flx to close the salt budget)
(  150)                   strdmp(ji,jj,jk) = zsa
(  151)                END DO
(  152)             END DO
(  153)          END DO
(  154) 
(  155)       CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
(  156)          DO jk = 1, jpkm1
(  157)             DO jj = 2, jpjm1
(  158)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  159)                   ztest = avt(ji,jj,jk) - 5.e-4
(  160)                   IF( ztest < 0. ) THEN
(  161)                      zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  162)                      zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 4

(  163)                   ELSE
(  164)                      zta = 0.e0
(  165)                      zsa = 0.e0
(  166)                   ENDIF
(  167)                   ! add the trends to the general tracer trends
(  168)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  169)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  170)                   ! save the salinity trend (used in flx to close the salt budget)
(  171)                   strdmp(ji,jj,jk) = zsa
(  172)                END DO
(  173)             END DO
(  174)          END DO
(  175) 
(  176)       CASE ( 2 )                ! no damping in the mixed layer 
(  177)          DO jk = 1, jpkm1
(  178)             DO jj = 2, jpjm1
(  179)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  180)                   IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
(  181)                      zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  182)                      zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  183)                   ELSE
(  184)                      zta = 0.e0
(  185)                      zsa = 0.e0
(  186)                   ENDIF
(  187)                   ! add the trends to the general tracer trends
(  188)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  189)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  190)                   ! save the salinity trend (used in flx to close the salt budget)
(  191)                   strdmp(ji,jj,jk) = zsa
(  192)                END DO
(  193)             END DO
(  194)          END DO
(  195) 
(  196)       END SELECT
(  197) 
(  198)       ! save the trends for diagnostic
(  199)       ! damping salinity trends
(  200)       IF( l_trdtra )   THEN
(  201)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  202)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  203)          CALL trd_mod(ztdta, ztdsa, jpttddoe, 'TRA', kt)
(  204)       ENDIF
(  205) 
(  206)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  207)          CALL prt_ctl(tab3d_1=ta, clinfo1=' dmp  - Ta: ', mask1=tmask, &
(  208)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  209)       ENDIF
(  210) 
(  211) 
(  212)    END SUBROUTINE tra_dmp
(  213) 
(  214) 
(  215)    SUBROUTINE tra_dmp_init
(  216)       !!----------------------------------------------------------------------
(  217)       !!                  ***  ROUTINE tra_dmp_init  ***
(  218)       !! 
(  219)       !! ** Purpose :   Initialization for the newtonian damping 
(  220)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 5

(  221)       !! ** Method  :   read the nammbf namelist and check the parameters
(  222)       !!      called by tra_dmp at the first timestep (nit000)
(  223)       !!
(  224)       !! History :
(  225)       !!   8.5  !  02-08  (G. Madec)  Original code
(  226)       !!----------------------------------------------------------------------
(  227)       !! * Local declarations
(  228)       NAMELIST/namtdp/ ndmp, ndmpf, nmldmp, sdmp, bdmp, hdmp
(  229)       !!----------------------------------------------------------------------
(  230) 
(  231)       ! Read Namelist namtdp : temperature and salinity damping term
(  232)       ! --------------------
(  233)       REWIND ( numnam )
(  234)       READ   ( numnam, namtdp )
(  235)       IF( lzoom )   nmldmp = 0           ! restoring to climatology at closed north or south boundaries
(  236) 
(  237)       ! Parameter control and print
(  238)       ! ---------------------------
(  239)       IF(lwp) THEN
(  240)          WRITE(numout,*)
(  241)          WRITE(numout,*) 'tra_dmp : T and S newtonian damping'
(  242)          WRITE(numout,*) '~~~~~~~'
(  243)          WRITE(numout,*) '          Namelist namtdp : set damping parameter'
(  244)          WRITE(numout,*)
(  245)          WRITE(numout,*) '             T and S damping option         ndmp   = ', ndmp
(  246)          WRITE(numout,*) '             create a damping.coeff file    ndmpf  = ', ndmpf
(  247)          WRITE(numout,*) '             mixed layer damping option     nmldmp = ', nmldmp, '(zoom: forced to 0)'
(  248)          WRITE(numout,*) '             surface time scale (days)      sdmp   = ', sdmp
(  249)          WRITE(numout,*) '             bottom time scale (days)       bdmp   = ', bdmp
(  250)          WRITE(numout,*) '             depth of transition (meters)   hdmp   = ', hdmp
(  251)          WRITE(numout,*)
(  252)       ENDIF
(  253) 
(  254)       SELECT CASE ( ndmp )
(  255) 
(  256) !byoung------------------------------------------------------------------------
(  257)       CASE ( 0 )               ! SS: damping in the Scotian Shelf only
(  258)          IF(lwp) WRITE(numout,*) '          tracer damping in the Scotian Shelf sea only'
(  259) !----------------------------------------------------------------------------	 
(  260)       CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
(  261)          IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'
(  262) 
(  263)       CASE ( 1:90 )             ! Damping poleward of 'ndmp' degrees
(  264)          IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmp, ' degrees'
(  265) 
(  266)       CASE DEFAULT
(  267)          IF(lwp) WRITE(numout,cform_err)
(  268)          IF(lwp) WRITE(numout,*) '          bad flag value for ndmp = ', ndmp
(  269)          nstop = nstop + 1
(  270) 
(  271)       END SELECT
(  272) 
(  273) 
(  274)       SELECT CASE ( nmldmp )
(  275) 
(  276)       CASE ( 0 )                ! newtonian damping throughout the water column
(  277)          IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 6

(  279)       CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
(  280)          IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'
(  281) 
(  282)       CASE ( 2 )                ! no damping in the mixed layer 
(  283)          IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'
(  284) 
(  285)       CASE DEFAULT
(  286)          IF(lwp) WRITE(numout,cform_err)
(  287)          IF(lwp) WRITE(numout,*) '          bad flag value for nmldmp = ', nmldmp
(  288)          nstop = nstop + 1
(  289) 
(  290)       END SELECT
(  291) 
(  292)       IF( .NOT.lk_dtasal .OR. .NOT.lk_dtatem ) THEN
(  293)          IF(lwp) WRITE(numout,cform_err)
(  294)          IF(lwp) WRITE(numout,*) '          no temperature and/or salinity data '
(  295)          IF(lwp) WRITE(numout,*) '          define key_dtatem and key_dtasal'
(  296)          nstop = nstop + 1
(  297)       ENDIF
(  298) 
(  299) 
(  300)       strdmp(:,:,:) = 0.e0       ! internal damping salinity trend (used in ocesbc)
(  301) 
(  302)       ! Damping coefficients initialization
(  303)       ! -----------------------------------
(  304) 
(  305)       IF( lzoom ) THEN
(  306)          CALL dtacof_zoom
(  307)       ELSE
(  308)          CALL dtacof
(  309)       ENDIF
(  310) 
(  311)    END SUBROUTINE tra_dmp_init
(  312) 
(  313) 
(  314)    SUBROUTINE dtacof_zoom
(  315)       !!----------------------------------------------------------------------
(  316)       !!                  ***  ROUTINE dtacof_zoom  ***
(  317)       !!
(  318)       !! ** Purpose :   Compute the damping coefficient for zoom domain
(  319)       !!
(  320)       !! ** Method  : - set along closed boundary due to zoom a damping over
(  321)       !!      6 points with a max time scale of 5 days.
(  322)       !!              - ORCA arctic/antarctic zoom: set the damping along
(  323)       !!      south/north boundary over a latitude strip.
(  324)       !!
(  325)       !! ** Action  : - resto, the damping coeff. for T and S
(  326)       !!
(  327)       !! History :
(  328)       !!   9.0  !  03-09  (G. Madec)  Original code
(  329)       !!----------------------------------------------------------------------
(  330)       !! * Local declarations
(  331)       INTEGER ::   ji, jj, jk, jn      ! dummy loop indices
(  332)       REAL(wp) ::   &
(  333)          zlat, zlat0, zlat1, zlat2     ! temporary scalar
(  334)       REAL(wp), DIMENSION(6)  ::   &
(  335)          zfact                         ! temporary workspace
(  336)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 7

(  337) 
(  338)       zfact(1) =  1.
(  339)       zfact(2) =  1. 
(  340)       zfact(3) = 11./12.
(  341)       zfact(4) =  8./12.
(  342)       zfact(5) =  4./12.
(  343)       zfact(6) =  1./12.
(  344)       zfact(:) = zfact(:) / ( 5. * rday )    ! 5 days max restoring time scale
(  345) 
(  346)       resto(:,:,:) = 0.e0
(  347) 
(  348)       ! damping along the forced closed boundary over 6 grid-points
(  349)       DO jn = 1, 6
(  350)          IF( lzoom_w )   resto( mi0(jn+jpizoom):mi1(jn+jpizoom), : , : ) = zfact(jn) ! west  closed
(  351)          IF( lzoom_s )   resto( : , mj0(jn+jpjzoom):mj1(jn+jpjzoom), : ) = zfact(jn) ! south closed 
(  352)          IF( lzoom_e )   resto( mi0(jpiglo+jpizoom-1-jn):mi1(jpiglo+jpizoom-1-jn) , : , : ) &
(  353)                        &              = zfact(jn)                                 ! east  closed 
(  354)          IF( lzoom_n )   resto( : , mj0(jpjglo+jpjzoom-1-jn):mj1(jpjglo+jpjzoom-1-jn) , : ) &
(  355)                        &              = zfact(jn)                                 ! north closed
(  356)       END DO
(  357) 
(  358) 
(  359)       IF( lzoom_arct .AND. lzoom_anta ) THEN
(  360) 
(  361)          ! ====================================================
(  362)          !  ORCA configuration : arctic zoom or antarctic zoom
(  363)          ! ====================================================
(  364) 
(  365)          IF(lwp) WRITE(numout,*)
(  366)          IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA    Arctic zoom'
(  367)          IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA Antarctic zoom'
(  368)          IF(lwp) WRITE(numout,*)
(  369) 
(  370)          ! ... Initialization : 
(  371)          !     zlat0 : latitude strip where resto decreases
(  372)          !     zlat1 : resto = 1 before zlat1
(  373)          !     zlat2 : resto decreases from 1 to 0 between zlat1 and zlat2
(  374)          resto(:,:,:) = 0.e0
(  375)          zlat0 = 10.
(  376)          zlat1 = 30.
(  377)          zlat2 = zlat1 + zlat0
(  378) 
(  379)          ! ... Compute arrays resto ; value for internal damping : 5 days
(  380)          DO jk = 2, jpkm1
(  381)             DO jj = 1, jpj
(  382)                DO ji = 1, jpi
(  383)                   zlat = ABS( gphit(ji,jj) )
(  384)                   IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
(  385)                      resto(ji,jj,jk) = 0.5 * ( 1./(5.*rday) ) *   &
(  386)                         ( 1. - cos(rpi*(zlat2-zlat)/zlat0) ) 
(  387)                   ELSE IF ( zlat < zlat1 ) THEN
(  388)                      resto(ji,jj,jk) = 1./(5.*rday)
(  389)                   ENDIF
(  390)                END DO
(  391)             END DO
(  392)          END DO
(  393) 
(  394)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 8

(  395) 
(  396)       ! ... Mask resto array
(  397)       resto(:,:,:) = resto(:,:,:) * tmask(:,:,:)
(  398) 
(  399)    END SUBROUTINE dtacof_zoom
(  400) 
(  401)    SUBROUTINE dtacof
(  402)       !!----------------------------------------------------------------------
(  403)       !!                  ***  ROUTINE dtacof  ***
(  404)       !!
(  405)       !! ** Purpose :   Compute the damping coefficient
(  406)       !!
(  407)       !! ** Method  :   Arrays defining the damping are computed for each grid
(  408)       !!      point for temperature and salinity (resto)
(  409)       !!      Damping depends on distance to coast, depth and latitude
(  410)       !!
(  411)       !! ** Action  : - resto, the damping coeff. for T and S
(  412)       !!
(  413)       !! History :
(  414)       !!   5.0  !  91-03  (O. Marti, G. Madec)  Original code
(  415)       !!        !  92-06  (M. Imbard)  doctor norme
(  416)       !!        !  96-01  (G. Madec) statement function for e3
(  417)       !!        !  98-07  (M. Imbard, G. Madec) ORCA version
(  418)       !!        !  00-08  (G. Madec, D. Ludicone) 
(  419)       !!----------------------------------------------------------------------
(  420)       !! * Modules used
(  421)       USE ioipsl
(  422) 
(  423)       !! * Local declarations
(  424)       INTEGER ::   ji, jj, jk, je      ! dummy loop indices
(  425)       INTEGER, PARAMETER ::   jpmois=1
(  426)       INTEGER ::   ipi, ipj, ipk       ! temporary integers
(  427)       INTEGER ::   ii0, ii1, ij0, ij1  !    "          "
(  428)       INTEGER ::   &
(  429)          idmp,     &  ! logical unit for file restoring damping term
(  430)          icot         ! logical unit for file distance to the coast
(  431)       INTEGER :: itime, istep(jpmois), ie
(  432)       LOGICAL :: llbon
(  433)       CHARACTER (len=32) ::  clname, clname2, clname3
(  434)       REAL(wp) ::   &
(  435)          zdate0, zinfl, zlon,         & ! temporary scalars
(  436)          zlat, zlat0, zlat1, zlat2,   & !    "         "
(  437)          zsdmp, zbdmp                   !    "         "
(  438)       REAL(wp), DIMENSION(jpk) ::   &
(  439)          zdept, zhfac
(  440)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  441)          zmrs, zlamt, zphit
(  442)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  443)          zdct
(  444) !!DB
(  445)       INTEGER :: k_hdmp
(  446) 
(  447)       !!----------------------------------------------------------------------
(  448) 
(  449)       ! ====================================
(  450)       !  ORCA configuration : global domain
(  451)       ! ====================================
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 9

(  453)       IF(lwp) WRITE(numout,*)
(  454)       IF(lwp) WRITE(numout,*) '              dtacof : Global domain of ORCA'
(  455)       IF(lwp) WRITE(numout,*) '              ------------------------------'
(  456) 
(  457)       ! ... Initialization : 
(  458)       !   zdct()      : distant to the coastline
(  459)       !   resto()     : array of restoring coeff. on T and S
(  460) 
(  461)       zdct (:,:,:) = 0.e0
(  462)       resto(:,:,:) = 0.e0
(  463) 
(  464)       IF ( ndmp > 0 ) THEN
(  465) 
(  466)          !    ------------------------------------
(  467)          !     Damping poleward of 'ndmp' degrees
(  468)          !    ------------------------------------
(  469) 
(  470)          IF(lwp) WRITE(numout,*)
(  471)          IF(lwp) WRITE(numout,*) '              Damping poleward of ', ndmp,' deg.'
(  472)          IF(lwp) WRITE(numout,*)
(  473) 
(  474)          ! ... Distance to coast (zdct)
(  475) 
(  476)          !   ... Test the existance of distance-to-coast file
(  477)          itime = jpmois
(  478)          ipi = jpiglo
(  479)          ipj = jpjglo
(  480)          ipk = jpk
(  481)          clname = 'dist.coast'
(  482)          DO je = 1,32
(  483)             IF( clname(je:je) == ' ' ) go to 140
(  484)          END DO
(  485) 140      CONTINUE
(  486)          ie = je
(  487)          clname2 = clname(1:ie-1)//".nc"
(  488)          inquire( FILE = clname2, EXIST = llbon )
(  489) 
(  490)          IF ( llbon ) THEN
(  491) 
(  492)             !   ... Read file distance to coast if possible
(  493)             CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj, .false.,   &
(  494)                ipi, ipj, ipk, zlamt, zphit, zdept, jpmois,   &
(  495)                istep, zdate0, rdt, icot )
(  496)             CALL flinget( icot, 'Tcoast', jpidta, jpjdta, jpk,    &
(  497)                jpmois, 1, 1, mig(1), nlci, mjg(1), nlcj, zdct(1:nlci,1:nlcj,1:jpk) )
(  498)             CALL flinclo( icot )
(  499)             IF(lwp)WRITE(numout,*) '    ** : File dist.coast.nc read'
(  500) 
(  501)          ELSE
(  502) 
(  503)             !   ... Compute and save the distance-to-coast array (output in zdct)
(  504)             CALL cofdis ( zdct )
(  505) 
(  506)          ENDIF
(  507) 
(  508)          ! ... Compute arrays resto 
(  509)          !      zinfl : distance of influence for damping term
(  510)          !      zlat0 : latitude strip where resto decreases






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 10

(  511)          !      zlat1 : resto = 0 between -zlat1 and zlat1
(  512)          !      zlat2 : resto increases from 0 to 1 between |zlat1| and |zlat2|
(  513)          !          and resto = 1 between |zlat2| and 90 deg.
(  514)          zinfl = 1000.e3
(  515)          zlat0 = 10
(  516)          zlat1 = ndmp
(  517)          zlat2 = zlat1 + zlat0
(  518) 
(  519)          DO jj = 1, jpj
(  520)             DO ji = 1, jpi
(  521)                zlat = ABS( gphit(ji,jj) )
(  522)                IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
(  523)                   resto(ji,jj,1) = 0.5 * ( 1. - cos(rpi*(zlat-zlat1)/zlat0 ) )
(  524)                ELSEIF ( zlat > zlat2 ) THEN
(  525)                   resto(ji,jj,1) = 1.
(  526)                ENDIF
(  527)             END DO
(  528)          END DO
(  529) 
(  530)          !   ... North Indian ocean (20N/30N x 45E/100E) : resto=0
(  531)          IF ( ndmp == 20 ) THEN
(  532)             DO jj = 1, jpj
(  533)                DO ji = 1, jpi
(  534)                   zlat = gphit(ji,jj)
(  535)                   zlon = MOD( glamt(ji,jj), 360. )
(  536)                   IF ( zlat1 < zlat .AND. zlat < zlat2 .AND.   &
(  537)                      45.  < zlon .AND. zlon < 100. ) THEN
(  538)                      resto(ji,jj,1) = 0.
(  539)                   ENDIF
(  540)                END DO
(  541)             END DO
(  542)          ENDIF
(  543) 
(  544)          zsdmp = 1./(sdmp * rday)
(  545)          zbdmp = 1./(bdmp * rday)
(  546)          DO jk = 2, jpkm1
(  547)             DO jj = 1, jpj
(  548)                DO ji = 1, jpi
(  549)                   zdct(ji,jj,jk) = MIN( zinfl, zdct(ji,jj,jk) )
(  550) 
(  551)                   !   ... Decrease the value in the vicinity of the coast
(  552)                   resto(ji,jj,jk) = resto(ji,jj,1)*0.5   &
(  553)                      &            * ( 1. - COS( rpi*zdct(ji,jj,jk)/zinfl) )
(  554) 
(  555)                   !   ... Vertical variation from zsdmp (sea surface) to zbdmp (bottom)
(  556)                   resto(ji,jj,jk) = resto(ji,jj,jk)   &
(  557)                      &            * ( zbdmp + (zsdmp-zbdmp)*EXP(-fsdept(ji,jj,jk)/hdmp) )
(  558)                END DO
(  559)             END DO
(  560)          END DO
(  561) 
(  562)       ENDIF
(  563) 
(  564) 
(  565) !!DB: deleted ORCA
(  566) !      IF( cp_cfg == "orca" .AND. ( ndmp > 0 .OR. ndmp == -1 ) ) THEN
(  567) !!DB -----------------------------------------------------------------------
(  568)       if( cp_cfg == "SS" .AND. ndmp == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 11

(  569)       
(  570)          resto(:,:,:) = 0.e0
(  571) 
(  572) !!DB 2008.06.30 
(  573) !Determine depth at which hdmp applies
(  574) !ifdef defaults to old system if key_flux_bulk_* not defined
(  575) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  576)          k_hdmp = jpk
(  577)          do jk = 1, jpk
(  578)             if(hdmp > gdept(jk)) k_hdmp = jk
(  579)          enddo
(  580) !!DBG
(  581)          if(lwp) write(numout2,*)'DBG: hdmp, k_hdmp: ', hdmp, k_hdmp 
(  582) #else
(  583)          k_hdmp = 0
(  584) #endif
(  585)          
(  586)          DO jk = 1, k_hdmp
(  587)             DO jj = 1, jpj
(  588)                DO ji = 1, jpi
(  589)                   resto(ji,jj,jk) = 1./(sdmp * rday)   
(  590)                END DO
(  591)             END DO
(  592)          END DO
(  593)          DO jk = k_hdmp+1, jpk
(  594)             DO jj = 1, jpj
(  595)                DO ji = 1, jpi
(  596)                   resto(ji,jj,jk) = 1./(bdmp * rday)   
(  597)                END DO
(  598)             END DO
(  599)          END DO
(  600)          
(  601)          resto(:,:, : ) = resto(:,:,:) * tmask(:,:,:)
(  602) !!DB
(  603)          resto(:,:, 1 ) = 1./(sdmp *rday) 
(  604)          
(  605)          !-------------------------------------------------------------------------------         	 
(  606)       ELSE
(  607)          !    ------------
(  608)          !     No damping
(  609)          !    ------------
(  610)          IF(lwp) WRITE(numout,cform_err)
(  611)          IF(lwp) WRITE(numout,*) 'Choose a correct value of ndmp or DO NOT defined key_tradmp'
(  612)          nstop = nstop + 1
(  613)       ENDIF
(  614) 
(  615)       !    ----------------------------
(  616)       !     Create Print damping array
(  617)       !    ----------------------------
(  618) 
(  619)       ! ndmpf   : = 1 create a damping.coeff NetCDF file
(  620) 
(  621)       IF( ndmpf == 1 ) THEN
(  622)          IF(lwp) WRITE(numout,*) '              create damping.coeff.nc file'
(  623)          itime   = 0
(  624)          clname3 = 'damping.coeff'
(  625)          CALL ymds2ju( 0     , 1     , 1      , 0.e0 , zdate0 )
(  626)          CALL restini( 'NONE', jpi   , jpj    , glamt, gphit,    &






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 12

(  627)                        jpk   , gdept , clname3, itime, zdate0,   &
(  628)                        rdt   , idmp, domain_id=nidom )
(  629)          CALL restput( idmp, 'Resto', jpi, jpj, jpk,   &
(  630)                        0   , resto  )
(  631)          CALL restclo( idmp )
(  632)       ENDIF
(  633) 
(  634)    END SUBROUTINE dtacof
(  635) 
(  636) 
(  637)    SUBROUTINE cofdis ( pdct )
(  638)       !!----------------------------------------------------------------------
(  639)       !!                 ***  ROUTINE cofdis  ***
(  640)       !!
(  641)       !! ** Purpose :   Compute the distance between ocean T-points and the
(  642)       !!      ocean model coastlines. Save the distance in a NetCDF file.
(  643)       !!
(  644)       !! ** Method  :   For each model level, the distance-to-coast is 
(  645)       !!      computed as follows : 
(  646)       !!       - The coastline is defined as the serie of U-,V-,F-points
(  647)       !!      that are at the ocean-land bound.
(  648)       !!       - For each ocean T-point, the distance-to-coast is then 
(  649)       !!      computed as the smallest distance (on the sphere) between the 
(  650)       !!      T-point and all the coastline points.
(  651)       !!       - For land T-points, the distance-to-coast is set to zero.
(  652)       !!      C A U T I O N : Computation not yet implemented in mpp case.
(  653)       !!
(  654)       !! ** Action  : - pdct, distance to the coastline (argument)
(  655)       !!              - NetCDF file 'dist.coast.nc' 
(  656)       !!        
(  657)       !! History :
(  658)       !!   7.0  !  01-02  (M. Imbard)  Original code
(  659)       !!   8.1  !  01-02  (G. Madec, E. Durand)
(  660)       !!   8.5  !  02-08  (G. Madec, E. Durand)  Free form, F90
(  661)       !!----------------------------------------------------------------------
(  662)       !! * Modules used
(  663)       USE ioipsl
(  664) 
(  665)       !! * Arguments
(  666)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  667)          pdct                     ! distance to the coastline
(  668) 
(  669)       !! * local declarations
(  670)       INTEGER :: ji, jj, jk, jl      ! dummy loop indices
(  671)       INTEGER :: iju, ijt            ! temporary integers
(  672)       INTEGER :: icoast, itime
(  673)       INTEGER ::   &
(  674)          icot         ! logical unit for file distance to the coast
(  675)       LOGICAL, DIMENSION(jpi,jpj) ::   &
(  676)          llcotu, llcotv, llcotf   ! ???
(  677)       CHARACTER (len=32) ::   clname
(  678)       REAL(wp) ::   zdate0
(  679)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  680)          zxt, zyt, zzt,                 &  ! cartesian coordinates for T-points
(  681)          zmask                             
(  682)       REAL(wp), DIMENSION(3*jpi*jpj) ::   &
(  683)          zxc, zyc, zzc, zdis      ! temporary workspace
(  684)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 13

(  685) 
(  686)       ! 0. Initialization
(  687)       ! -----------------
(  688)       IF(lwp) WRITE(numout,*)
(  689)       IF(lwp) WRITE(numout,*) 'cofdis : compute the distance to coastline'
(  690)       IF(lwp) WRITE(numout,*) '~~~~~~'
(  691)       IF(lwp) WRITE(numout,*)
(  692)       IF( lk_mpp ) THEN
(  693)          IF(lwp) WRITE(numout,cform_err)
(  694)          IF(lwp) WRITE(numout,*) '         Computation not yet implemented with key_mpp_...'
(  695)          IF(lwp) WRITE(numout,*) '         Rerun the code on another computer or '
(  696)          IF(lwp) WRITE(numout,*) '         create the "dist.coast.nc" file using IDL'
(  697)          nstop = nstop + 1
(  698)       ENDIF
(  699) 
(  700)       pdct(:,:,:) = 0.e0
(  701)       zxt(:,:) = cos( rad * gphit(:,:) ) * cos( rad * glamt(:,:) )
(  702)       zyt(:,:) = cos( rad * gphit(:,:) ) * sin( rad * glamt(:,:) )
(  703)       zzt(:,:) = sin( rad * gphit(:,:) )
(  704) 
(  705) 
(  706)       ! 1. Loop on vertical levels
(  707)       ! --------------------------
(  708)       !                                                ! ===============
(  709)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  710)          !                                             ! ===============
(  711)          ! Define the coastline points (U, V and F)
(  712)          DO jj = 2, jpjm1
(  713)             DO ji = 2, jpim1
(  714)                zmask(ji,jj) =  ( tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) &
(  715)                    &           + tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk) )
(  716)                llcotu(ji,jj) = ( tmask(ji,jj,  jk) + tmask(ji+1,jj  ,jk) == 1. ) 
(  717)                llcotv(ji,jj) = ( tmask(ji,jj  ,jk) + tmask(ji  ,jj+1,jk) == 1. ) 
(  718)                llcotf(ji,jj) = ( zmask(ji,jj) > 0. ) .AND. ( zmask(ji,jj) < 4. )
(  719)             END DO
(  720)          END DO
(  721) 
(  722)          ! Lateral boundaries conditions
(  723)          llcotu(:, 1 ) = umask(:,  2  ,jk) == 1
(  724)          llcotu(:,jpj) = umask(:,jpjm1,jk) == 1
(  725)          llcotv(:, 1 ) = vmask(:,  2  ,jk) == 1
(  726)          llcotv(:,jpj) = vmask(:,jpjm1,jk) == 1
(  727)          llcotf(:, 1 ) = fmask(:,  2  ,jk) == 1
(  728)          llcotf(:,jpj) = fmask(:,jpjm1,jk) == 1
(  729) 
(  730)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  731)             llcotu( 1 ,:) = llcotu(jpim1,:)
(  732)             llcotu(jpi,:) = llcotu(  2  ,:)
(  733)             llcotv( 1 ,:) = llcotv(jpim1,:)
(  734)             llcotv(jpi,:) = llcotv(  2  ,:)
(  735)             llcotf( 1 ,:) = llcotf(jpim1,:)
(  736)             llcotf(jpi,:) = llcotf(  2  ,:)
(  737)          ELSE
(  738)             llcotu( 1 ,:) = umask(  2  ,:,jk) == 1
(  739)             llcotu(jpi,:) = umask(jpim1,:,jk) == 1
(  740)             llcotv( 1 ,:) = vmask(  2  ,:,jk) == 1
(  741)             llcotv(jpi,:) = vmask(jpim1,:,jk) == 1
(  742)             llcotf( 1 ,:) = fmask(  2  ,:,jk) == 1






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 14

(  743)             llcotf(jpi,:) = fmask(jpim1,:,jk) == 1
(  744)          ENDIF
(  745)          IF( nperio == 3 .OR. nperio == 4 ) THEN
(  746)             DO ji = 1, jpim1
(  747)                iju = jpi - ji + 1
(  748)                llcotu(ji,jpj  ) = llcotu(iju,jpj-2)
(  749)                llcotf(ji,jpj-1) = llcotf(iju,jpj-2)
(  750)                llcotf(ji,jpj  ) = llcotf(iju,jpj-3)
(  751)             END DO
(  752)             DO ji = jpi/2, jpi-1
(  753)                iju = jpi - ji + 1
(  754)                llcotu(ji,jpjm1) = llcotu(iju,jpjm1)
(  755)             END DO
(  756)             DO ji = 2, jpi
(  757)                ijt = jpi - ji + 2
(  758)                llcotv(ji,jpj-1) = llcotv(ijt,jpj-2)
(  759)                llcotv(ji,jpj  ) = llcotv(ijt,jpj-3)
(  760)             END DO
(  761)          ENDIF
(  762)          IF( nperio == 5 .OR. nperio == 6 ) THEN
(  763)             DO ji = 1, jpim1
(  764)                iju = jpi - ji
(  765)                llcotu(ji,jpj  ) = llcotu(iju,jpj-1)
(  766)                llcotf(ji,jpj  ) = llcotf(iju,jpj-2)
(  767)             END DO
(  768)             DO ji = jpi/2, jpi-1
(  769)                iju = jpi - ji
(  770)                llcotf(ji,jpjm1) = llcotf(iju,jpjm1)
(  771)             END DO
(  772)             DO ji = 1, jpi
(  773)                ijt = jpi - ji + 1
(  774)                llcotv(ji,jpj  ) = llcotv(ijt,jpj-1)
(  775)             END DO
(  776)             DO ji = jpi/2+1, jpi
(  777)                ijt = jpi - ji + 1
(  778)                llcotv(ji,jpjm1) = llcotv(ijt,jpjm1)
(  779)             END DO
(  780)          ENDIF
(  781) 
(  782)          ! Compute cartesian coordinates of coastline points
(  783)          ! and the number of coastline points
(  784) 
(  785)          icoast = 0
(  786)          DO jj = 1, jpj
(  787)             DO ji = 1, jpi
(  788)                IF( llcotf(ji,jj) ) THEN
(  789)                   icoast = icoast + 1
(  790)                   zxc(icoast) = COS( rad*gphif(ji,jj) ) * COS( rad*glamf(ji,jj) )
(  791)                   zyc(icoast) = COS( rad*gphif(ji,jj) ) * SIN( rad*glamf(ji,jj) )
(  792)                   zzc(icoast) = SIN( rad*gphif(ji,jj) )
(  793)                ENDIF
(  794)                IF( llcotu(ji,jj) ) THEN
(  795)                   icoast = icoast+1
(  796)                   zxc(icoast) = COS( rad*gphiu(ji,jj) ) * COS( rad*glamu(ji,jj) )
(  797)                   zyc(icoast) = COS( rad*gphiu(ji,jj) ) * SIN( rad*glamu(ji,jj) )
(  798)                   zzc(icoast) = SIN( rad*gphiu(ji,jj) )
(  799)                ENDIF
(  800)                IF( llcotv(ji,jj) ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 15

(  801)                   icoast = icoast+1
(  802)                   zxc(icoast) = COS( rad*gphiv(ji,jj) ) * COS( rad*glamv(ji,jj) )
(  803)                   zyc(icoast) = COS( rad*gphiv(ji,jj) ) * SIN( rad*glamv(ji,jj) )
(  804)                   zzc(icoast) = SIN( rad*gphiv(ji,jj) )
(  805)                ENDIF
(  806)             END DO
(  807)          END DO
(  808) 
(  809)          ! Distance for the T-points
(  810) 
(  811)          DO jj = 1, jpj
(  812)             DO ji = 1, jpi
(  813)                IF( tmask(ji,jj,jk) == 0. ) THEN
(  814)                   pdct(ji,jj,jk) = 0.
(  815)                ELSE
(  816)                   DO jl = 1, icoast
(  817)                      zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
(  818)                               + ( zyt(ji,jj) - zyc(jl) )**2   &
(  819)                               + ( zzt(ji,jj) - zzc(jl) )**2
(  820)                   END DO
(  821)                   pdct(ji,jj,jk) = ra * SQRT( MINVAL( zdis(1:icoast) ) )
(  822)                ENDIF
(  823)             END DO
(  824)          END DO
(  825)          !                                                ! ===============
(  826)       END DO                                              !   End of slab
(  827)       !                                                   ! ===============
(  828) 
(  829) 
(  830)       ! 2. Create the  distance to the coast file in NetCDF format
(  831)       ! ----------------------------------------------------------    
(  832)       clname = 'dist.coast'
(  833)       itime = 0
(  834)       CALL ymds2ju( 0     , 1     , 1     , 0.e0 , zdate0 )
(  835)       CALL restini( 'NONE', jpi   , jpj   , glamt, gphit ,   &
(  836)                     jpk   , gdept , clname, itime, zdate0,   &
(  837)                     rdt   , icot                         )
(  838)       CALL restput( icot, 'Tcoast', jpi, jpj, jpk, 0, pdct )
(  839)       CALL restclo( icot )
(  840) 
(  841)    END SUBROUTINE cofdis
(  842) 
(  843) #else
(  844)    !!----------------------------------------------------------------------
(  845)    !!   Default key                                     NO internal damping
(  846)    !!----------------------------------------------------------------------
(  847)    LOGICAL , PUBLIC, PARAMETER ::   lk_tradmp = .FALSE.    !: internal damping flag
(  848) CONTAINS
(  849)    SUBROUTINE tra_dmp( kt )        ! Empty routine
(  850) !      WRITE(*,*) 'tra_dmp: You should not have seen this print! error?', kt
(  851)    END SUBROUTINE tra_dmp
(  852) #endif
(  853) 
(  854)    !!======================================================================
(  855) END MODULE tradmp









PGF90 (Version     10.5)          09/23/2010  09:49:13      page 16

(    1) # 1 "tradmp.F90"
(    1) MODULE tradmp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  tradmp  ***
(    4)    !! Ocean physics: internal restoring trend on active tracers (T and S)
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   1                                         internal damping
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_dmp      : update the tracer trend with the internal damping
(   11)    !!   tra_dmp_init : initialization, namlist read, parameters control
(   12)    !!   dtacof_zoom  : restoring coefficient for zoom domain
(   13)    !!   dtacof       : restoring coefficient for global domain
(   14)    !!   cofdis       : compute the distance to the coastline
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers variables
(   18)    USE dom_oce         ! ocean space and time domain variables
(   19)    USE trdmod          ! ocean active tracers trends 
(   20)    USE trdmod_oce      ! ocean variables trends
(   21)    USE zdf_oce         ! ocean vertical physics
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE phycst          ! Define parameters for the routines
(   24)    USE dtatem          ! temperature data
(   25)    USE dtasal          ! salinity data
(   26)    USE zdfmxl          ! mixed layer depth
(   27)    USE lib_mpp         ! distribued memory computing
(   28)    USE prtctl          ! Print control
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Routine accessibility
(   34)    PUBLIC tra_dmp   ! routine called by step.F90
(   35) 
(   36)    !! * Shared module variables
(   37)    LOGICAL , PUBLIC &
(   39)    , PARAMETER  &
(   41)    ::   lk_tradmp = .TRUE.    !: internal damping flag
(   42) 
(   43)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &
(   44)       strdmp,              &  ! damping salinity trend (psu/s)
(   45)       resto                   ! restoring coeff. on T and S (s-1)
(   46) 
(   47)    !! * Module variables
(   48)    INTEGER  ::             & !!! * newtonian damping namelist (mandmp) *
(   49)       ndmp   =   -1 ,      &  ! = 0/-1/'latitude' for damping over T and S
(   50)       ndmpf  =    2 ,      &  ! = 1 create a damping.coeff NetCDF file 
(   51)       nmldmp =    0           ! = 0/1/2 flag for damping in the mixed layer
(   52)    REAL(wp) ::             & !!!  * newtonian damping namelist *
(   53)       sdmp   =   50.,      &  ! surface time scale for internal damping (days)
(   54)       bdmp   =  360.,      &  ! bottom time scale for internal damping (days)
(   55)       hdmp   =  800.          ! depth of transition between sdmp and bdmp (meters)
(   56) 
(   57)    !! * Substitutions
(   58) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 17

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
(   48) # 59 "tradmp.F90"
(   59) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 60 "tradmp.F90"
(   60) # 60
(   60)    !!----------------------------------------------------------------------
(   61)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   62)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tradmp.F90,v 1.15 2006/04/19 14:43:17 opalod Exp $ 
(   63)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   64)    !!----------------------------------------------------------------------
(   65) 
(   66) CONTAINS
(   67) 
(   68)    SUBROUTINE tra_dmp( kt )
(   69)       !!----------------------------------------------------------------------
(   70)       !!                   ***  ROUTINE tra_dmp  ***
(   71)       !!                  
(   72)       !! ** Purpose :   Compute the tracer trend due to a newtonian damping
(   73)       !!      of the tracer field towards given data field and add it to the
(   74)       !!      general tracer trends.
(   75)       !!
(   76)       !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
(   77)       !!      and add to the general tracer trends:
(   78)       !!                     ta = ta + resto * (t_dta - tb)
(   79)       !!                     sa = sa + resto * (s_dta - sb)
(   80)       !!         The trend is computed either throughout the water column
(   81)       !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
(   82)       !!      below the well mixed layer (nlmdmp=2)
(   83)       !!
(   84)       !! ** Action  : - update the tracer trends (ta,sa) with the newtonian 
(   85)       !!                damping trends.
(   86)       !!              - save the trends in (ttrd,strd) ('key_trdtra')
(   87)       !!
(   88)       !! History :
(   89)       !!   7.0  !         (G. Madec)  Original code






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 18

(   90)       !!        !  96-01  (G. Madec) 
(   91)       !!        !  97-05  (G. Madec)  macro-tasked on jk-slab
(   92)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(   93)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   94)       !!----------------------------------------------------------------------
(   95)       !! * Modules used     
(   96)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   97)                          ztdsa => va      ! use va as 3D workspace   
(   98) 
(   99)       !! * Arguments
(  100)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(  101) 
(  102)       !! * Local declarations
(  103)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(  104)       REAL(wp) ::   ztest, zta, zsa       ! temporary scalars
(  105)       !!----------------------------------------------------------------------
(  106) 
(  107)       ! 0. Initialization (first time-step only)
(  108)       !    --------------
(  109)       IF( kt == nit000 ) CALL tra_dmp_init
(  110) 
(  111)       ! Save ta and sa trends
(  112)       IF( l_trdtra )   THEN
(  113)          ztdta(:,:,:) = ta(:,:,:) 
(  114)          ztdsa(:,:,:) = sa(:,:,:) 
(  115)       ENDIF
(  116) 
(  117)       ! 1. Newtonian damping trends on tracer fields
(  118)       ! --------------------------------------------
(  119)       !    compute the newtonian damping trends depending on nmldmp
(  120) 
(  121)       SELECT CASE ( nmldmp )
(  122) 
(  123)       CASE( 0 )                ! newtonian damping throughout the water column
(  124) !!DBG: 2009.06.24 -- HARDWIRE strong sfce S restoring
(  125) !         DO jk = 1, 1
(  126) !            DO jj = 2, jpjm1
(  127) !               DO ji = 2, jpim1   ! vector opt.
(  128) !                  zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  129) !!                  zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  130) !                  zsa = (1./(3.*rday)) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  131) !                  ! add the trends to the general tracer trends
(  132) !                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  133) !                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  134) !                  ! save the salinity trend (used in flx to close the salt budget)
(  135) !                  strdmp(ji,jj,jk) = zsa
(  136) !               END DO
(  137) !            END DO
(  138) !         END DO
(  139) !         DO jk = 2, jpkm1
(  140) 
(  141)          DO jk = 1, jpkm1
(  142)             DO jj = 2, jpjm1
(  143)                DO ji = 2, jpim1   ! vector opt.
(  144)                   zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  145)                   zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  146)                   ! add the trends to the general tracer trends
(  147)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 19

(  148)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  149)                   ! save the salinity trend (used in flx to close the salt budget)
(  150)                   strdmp(ji,jj,jk) = zsa
(  151)                END DO
(  152)             END DO
(  153)          END DO
(  154) 
(  155)       CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
(  156)          DO jk = 1, jpkm1
(  157)             DO jj = 2, jpjm1
(  158)                DO ji = 2, jpim1   ! vector opt.
(  159)                   ztest = avt(ji,jj,jk) - 5.e-4
(  160)                   IF( ztest < 0. ) THEN
(  161)                      zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  162)                      zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  163)                   ELSE
(  164)                      zta = 0.e0
(  165)                      zsa = 0.e0
(  166)                   ENDIF
(  167)                   ! add the trends to the general tracer trends
(  168)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  169)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  170)                   ! save the salinity trend (used in flx to close the salt budget)
(  171)                   strdmp(ji,jj,jk) = zsa
(  172)                END DO
(  173)             END DO
(  174)          END DO
(  175) 
(  176)       CASE ( 2 )                ! no damping in the mixed layer 
(  177)          DO jk = 1, jpkm1
(  178)             DO jj = 2, jpjm1
(  179)                DO ji = 2, jpim1   ! vector opt.
(  180)                   IF( gdept(jk) >= hmlp (ji,jj) ) THEN
(  181)                      zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
(  182)                      zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
(  183)                   ELSE
(  184)                      zta = 0.e0
(  185)                      zsa = 0.e0
(  186)                   ENDIF
(  187)                   ! add the trends to the general tracer trends
(  188)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  189)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  190)                   ! save the salinity trend (used in flx to close the salt budget)
(  191)                   strdmp(ji,jj,jk) = zsa
(  192)                END DO
(  193)             END DO
(  194)          END DO
(  195) 
(  196)       END SELECT
(  197) 
(  198)       ! save the trends for diagnostic
(  199)       ! damping salinity trends
(  200)       IF( l_trdtra )   THEN
(  201)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  202)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  203)          CALL trd_mod(ztdta, ztdsa, jpttddoe, 'TRA', kt)
(  204)       ENDIF
(  205) 






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 20

(  206)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  207)          CALL prt_ctl(tab3d_1=ta, clinfo1=' dmp  - Ta: ', mask1=tmask, &
(  208)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  209)       ENDIF
(  210) 
(  211) 
(  212)    END SUBROUTINE tra_dmp
(  213) 
(  214) 
(  215)    SUBROUTINE tra_dmp_init
(  216)       !!----------------------------------------------------------------------
(  217)       !!                  ***  ROUTINE tra_dmp_init  ***
(  218)       !! 
(  219)       !! ** Purpose :   Initialization for the newtonian damping 
(  220)       !!
(  221)       !! ** Method  :   read the nammbf namelist and check the parameters
(  222)       !!      called by tra_dmp at the first timestep (nit000)
(  223)       !!
(  224)       !! History :
(  225)       !!   8.5  !  02-08  (G. Madec)  Original code
(  226)       !!----------------------------------------------------------------------
(  227)       !! * Local declarations
(  228)       NAMELIST/namtdp/ ndmp, ndmpf, nmldmp, sdmp, bdmp, hdmp
(  229)       !!----------------------------------------------------------------------
(  230) 
(  231)       ! Read Namelist namtdp : temperature and salinity damping term
(  232)       ! --------------------
(  233)       REWIND ( numnam )
(  234)       READ   ( numnam, namtdp )
(  235)       IF( lzoom )   nmldmp = 0           ! restoring to climatology at closed north or south boundaries
(  236) 
(  237)       ! Parameter control and print
(  238)       ! ---------------------------
(  239)       IF(lwp) THEN
(  240)          WRITE(numout,*)
(  241)          WRITE(numout,*) 'tra_dmp : T and S newtonian damping'
(  242)          WRITE(numout,*) '~~~~~~~'
(  243)          WRITE(numout,*) '          Namelist namtdp : set damping parameter'
(  244)          WRITE(numout,*)
(  245)          WRITE(numout,*) '             T and S damping option         ndmp   = ', ndmp
(  246)          WRITE(numout,*) '             create a damping.coeff file    ndmpf  = ', ndmpf
(  247)          WRITE(numout,*) '             mixed layer damping option     nmldmp = ', nmldmp, '(zoom: forced to 0)'
(  248)          WRITE(numout,*) '             surface time scale (days)      sdmp   = ', sdmp
(  249)          WRITE(numout,*) '             bottom time scale (days)       bdmp   = ', bdmp
(  250)          WRITE(numout,*) '             depth of transition (meters)   hdmp   = ', hdmp
(  251)          WRITE(numout,*)
(  252)       ENDIF
(  253) 
(  254)       SELECT CASE ( ndmp )
(  255) 
(  256) !byoung------------------------------------------------------------------------
(  257)       CASE ( 0 )               ! SS: damping in the Scotian Shelf only
(  258)          IF(lwp) WRITE(numout,*) '          tracer damping in the Scotian Shelf sea only'
(  259) !----------------------------------------------------------------------------	 
(  260)       CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
(  261)          IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'
(  262) 
(  263)       CASE ( 1:90 )             ! Damping poleward of 'ndmp' degrees






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 21

(  264)          IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmp, ' degrees'
(  265) 
(  266)       CASE DEFAULT
(  267)          IF(lwp) WRITE(numout,cform_err)
(  268)          IF(lwp) WRITE(numout,*) '          bad flag value for ndmp = ', ndmp
(  269)          nstop = nstop + 1
(  270) 
(  271)       END SELECT
(  272) 
(  273) 
(  274)       SELECT CASE ( nmldmp )
(  275) 
(  276)       CASE ( 0 )                ! newtonian damping throughout the water column
(  277)          IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'
(  278) 
(  279)       CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
(  280)          IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'
(  281) 
(  282)       CASE ( 2 )                ! no damping in the mixed layer 
(  283)          IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'
(  284) 
(  285)       CASE DEFAULT
(  286)          IF(lwp) WRITE(numout,cform_err)
(  287)          IF(lwp) WRITE(numout,*) '          bad flag value for nmldmp = ', nmldmp
(  288)          nstop = nstop + 1
(  289) 
(  290)       END SELECT
(  291) 
(  292)       IF( .NOT.lk_dtasal .OR. .NOT.lk_dtatem ) THEN
(  293)          IF(lwp) WRITE(numout,cform_err)
(  294)          IF(lwp) WRITE(numout,*) '          no temperature and/or salinity data '
(  295)          IF(lwp) WRITE(numout,*) '          define key_dtatem and key_dtasal'
(  296)          nstop = nstop + 1
(  297)       ENDIF
(  298) 
(  299) 
(  300)       strdmp(:,:,:) = 0.e0       ! internal damping salinity trend (used in ocesbc)
(  301) 
(  302)       ! Damping coefficients initialization
(  303)       ! -----------------------------------
(  304) 
(  305)       IF( lzoom ) THEN
(  306)          CALL dtacof_zoom
(  307)       ELSE
(  308)          CALL dtacof
(  309)       ENDIF
(  310) 
(  311)    END SUBROUTINE tra_dmp_init
(  312) 
(  313) 
(  314)    SUBROUTINE dtacof_zoom
(  315)       !!----------------------------------------------------------------------
(  316)       !!                  ***  ROUTINE dtacof_zoom  ***
(  317)       !!
(  318)       !! ** Purpose :   Compute the damping coefficient for zoom domain
(  319)       !!
(  320)       !! ** Method  : - set along closed boundary due to zoom a damping over
(  321)       !!      6 points with a max time scale of 5 days.






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 22

(  322)       !!              - ORCA arctic/antarctic zoom: set the damping along
(  323)       !!      south/north boundary over a latitude strip.
(  324)       !!
(  325)       !! ** Action  : - resto, the damping coeff. for T and S
(  326)       !!
(  327)       !! History :
(  328)       !!   9.0  !  03-09  (G. Madec)  Original code
(  329)       !!----------------------------------------------------------------------
(  330)       !! * Local declarations
(  331)       INTEGER ::   ji, jj, jk, jn      ! dummy loop indices
(  332)       REAL(wp) ::   &
(  333)          zlat, zlat0, zlat1, zlat2     ! temporary scalar
(  334)       REAL(wp), DIMENSION(6)  ::   &
(  335)          zfact                         ! temporary workspace
(  336)       !!----------------------------------------------------------------------
(  337) 
(  338)       zfact(1) =  1.
(  339)       zfact(2) =  1. 
(  340)       zfact(3) = 11./12.
(  341)       zfact(4) =  8./12.
(  342)       zfact(5) =  4./12.
(  343)       zfact(6) =  1./12.
(  344)       zfact(:) = zfact(:) / ( 5. * rday )    ! 5 days max restoring time scale
(  345) 
(  346)       resto(:,:,:) = 0.e0
(  347) 
(  348)       ! damping along the forced closed boundary over 6 grid-points
(  349)       DO jn = 1, 6
(  350)          IF( lzoom_w )   resto( mi0(jn+jpizoom):mi1(jn+jpizoom), : , : ) = zfact(jn) ! west  closed
(  351)          IF( lzoom_s )   resto( : , mj0(jn+jpjzoom):mj1(jn+jpjzoom), : ) = zfact(jn) ! south closed 
(  352)          IF( lzoom_e )   resto( mi0(jpiglo+jpizoom-1-jn):mi1(jpiglo+jpizoom-1-jn) , : , : ) &
(  353)                        &              = zfact(jn)                                 ! east  closed 
(  354)          IF( lzoom_n )   resto( : , mj0(jpjglo+jpjzoom-1-jn):mj1(jpjglo+jpjzoom-1-jn) , : ) &
(  355)                        &              = zfact(jn)                                 ! north closed
(  356)       END DO
(  357) 
(  358) 
(  359)       IF( lzoom_arct .AND. lzoom_anta ) THEN
(  360) 
(  361)          ! ====================================================
(  362)          !  ORCA configuration : arctic zoom or antarctic zoom
(  363)          ! ====================================================
(  364) 
(  365)          IF(lwp) WRITE(numout,*)
(  366)          IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA    Arctic zoom'
(  367)          IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA Antarctic zoom'
(  368)          IF(lwp) WRITE(numout,*)
(  369) 
(  370)          ! ... Initialization : 
(  371)          !     zlat0 : latitude strip where resto decreases
(  372)          !     zlat1 : resto = 1 before zlat1
(  373)          !     zlat2 : resto decreases from 1 to 0 between zlat1 and zlat2
(  374)          resto(:,:,:) = 0.e0
(  375)          zlat0 = 10.
(  376)          zlat1 = 30.
(  377)          zlat2 = zlat1 + zlat0
(  378) 
(  379)          ! ... Compute arrays resto ; value for internal damping : 5 days






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 23

(  380)          DO jk = 2, jpkm1
(  381)             DO jj = 1, jpj
(  382)                DO ji = 1, jpi
(  383)                   zlat = ABS( gphit(ji,jj) )
(  384)                   IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
(  385)                      resto(ji,jj,jk) = 0.5 * ( 1./(5.*rday) ) *   &
(  386)                         ( 1. - cos(rpi*(zlat2-zlat)/zlat0) ) 
(  387)                   ELSE IF ( zlat < zlat1 ) THEN
(  388)                      resto(ji,jj,jk) = 1./(5.*rday)
(  389)                   ENDIF
(  390)                END DO
(  391)             END DO
(  392)          END DO
(  393) 
(  394)       ENDIF
(  395) 
(  396)       ! ... Mask resto array
(  397)       resto(:,:,:) = resto(:,:,:) * tmask(:,:,:)
(  398) 
(  399)    END SUBROUTINE dtacof_zoom
(  400) 
(  401)    SUBROUTINE dtacof
(  402)       !!----------------------------------------------------------------------
(  403)       !!                  ***  ROUTINE dtacof  ***
(  404)       !!
(  405)       !! ** Purpose :   Compute the damping coefficient
(  406)       !!
(  407)       !! ** Method  :   Arrays defining the damping are computed for each grid
(  408)       !!      point for temperature and salinity (resto)
(  409)       !!      Damping depends on distance to coast, depth and latitude
(  410)       !!
(  411)       !! ** Action  : - resto, the damping coeff. for T and S
(  412)       !!
(  413)       !! History :
(  414)       !!   5.0  !  91-03  (O. Marti, G. Madec)  Original code
(  415)       !!        !  92-06  (M. Imbard)  doctor norme
(  416)       !!        !  96-01  (G. Madec) statement function for e3
(  417)       !!        !  98-07  (M. Imbard, G. Madec) ORCA version
(  418)       !!        !  00-08  (G. Madec, D. Ludicone) 
(  419)       !!----------------------------------------------------------------------
(  420)       !! * Modules used
(  421)       USE ioipsl
(  422) 
(  423)       !! * Local declarations
(  424)       INTEGER ::   ji, jj, jk, je      ! dummy loop indices
(  425)       INTEGER, PARAMETER ::   jpmois=1
(  426)       INTEGER ::   ipi, ipj, ipk       ! temporary integers
(  427)       INTEGER ::   ii0, ii1, ij0, ij1  !    "          "
(  428)       INTEGER ::   &
(  429)          idmp,     &  ! logical unit for file restoring damping term
(  430)          icot         ! logical unit for file distance to the coast
(  431)       INTEGER :: itime, istep(jpmois), ie
(  432)       LOGICAL :: llbon
(  433)       CHARACTER (len=32) ::  clname, clname2, clname3
(  434)       REAL(wp) ::   &
(  435)          zdate0, zinfl, zlon,         & ! temporary scalars
(  436)          zlat, zlat0, zlat1, zlat2,   & !    "         "
(  437)          zsdmp, zbdmp                   !    "         "






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 24

(  438)       REAL(wp), DIMENSION(jpk) ::   &
(  439)          zdept, zhfac
(  440)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  441)          zmrs, zlamt, zphit
(  442)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  443)          zdct
(  444) !!DB
(  445)       INTEGER :: k_hdmp
(  446) 
(  447)       !!----------------------------------------------------------------------
(  448) 
(  449)       ! ====================================
(  450)       !  ORCA configuration : global domain
(  451)       ! ====================================
(  452) 
(  453)       IF(lwp) WRITE(numout,*)
(  454)       IF(lwp) WRITE(numout,*) '              dtacof : Global domain of ORCA'
(  455)       IF(lwp) WRITE(numout,*) '              ------------------------------'
(  456) 
(  457)       ! ... Initialization : 
(  458)       !   zdct()      : distant to the coastline
(  459)       !   resto()     : array of restoring coeff. on T and S
(  460) 
(  461)       zdct (:,:,:) = 0.e0
(  462)       resto(:,:,:) = 0.e0
(  463) 
(  464)       IF ( ndmp > 0 ) THEN
(  465) 
(  466)          !    ------------------------------------
(  467)          !     Damping poleward of 'ndmp' degrees
(  468)          !    ------------------------------------
(  469) 
(  470)          IF(lwp) WRITE(numout,*)
(  471)          IF(lwp) WRITE(numout,*) '              Damping poleward of ', ndmp,' deg.'
(  472)          IF(lwp) WRITE(numout,*)
(  473) 
(  474)          ! ... Distance to coast (zdct)
(  475) 
(  476)          !   ... Test the existance of distance-to-coast file
(  477)          itime = jpmois
(  478)          ipi = jpiglo
(  479)          ipj = jpjglo
(  480)          ipk = jpk
(  481)          clname = 'dist.coast'
(  482)          DO je = 1,32
(  483)             IF( clname(je:je) == ' ' ) go to 140
(  484)          END DO
(  485) 140      CONTINUE
(  486)          ie = je
(  487)          clname2 = clname(1:ie-1)//".nc"
(  488)          inquire( FILE = clname2, EXIST = llbon )
(  489) 
(  490)          IF ( llbon ) THEN
(  491) 
(  492)             !   ... Read file distance to coast if possible
(  493)             CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj, .false.,   &
(  494)                ipi, ipj, ipk, zlamt, zphit, zdept, jpmois,   &
(  495)                istep, zdate0, rdt, icot )






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 25

(  496)             CALL flinget( icot, 'Tcoast', jpidta, jpjdta, jpk,    &
(  497)                jpmois, 1, 1, mig(1), nlci, mjg(1), nlcj, zdct(1:nlci,1:nlcj,1:jpk) )
(  498)             CALL flinclo( icot )
(  499)             IF(lwp)WRITE(numout,*) '    ** : File dist.coast.nc read'
(  500) 
(  501)          ELSE
(  502) 
(  503)             !   ... Compute and save the distance-to-coast array (output in zdct)
(  504)             CALL cofdis ( zdct )
(  505) 
(  506)          ENDIF
(  507) 
(  508)          ! ... Compute arrays resto 
(  509)          !      zinfl : distance of influence for damping term
(  510)          !      zlat0 : latitude strip where resto decreases
(  511)          !      zlat1 : resto = 0 between -zlat1 and zlat1
(  512)          !      zlat2 : resto increases from 0 to 1 between |zlat1| and |zlat2|
(  513)          !          and resto = 1 between |zlat2| and 90 deg.
(  514)          zinfl = 1000.e3
(  515)          zlat0 = 10
(  516)          zlat1 = ndmp
(  517)          zlat2 = zlat1 + zlat0
(  518) 
(  519)          DO jj = 1, jpj
(  520)             DO ji = 1, jpi
(  521)                zlat = ABS( gphit(ji,jj) )
(  522)                IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
(  523)                   resto(ji,jj,1) = 0.5 * ( 1. - cos(rpi*(zlat-zlat1)/zlat0 ) )
(  524)                ELSEIF ( zlat > zlat2 ) THEN
(  525)                   resto(ji,jj,1) = 1.
(  526)                ENDIF
(  527)             END DO
(  528)          END DO
(  529) 
(  530)          !   ... North Indian ocean (20N/30N x 45E/100E) : resto=0
(  531)          IF ( ndmp == 20 ) THEN
(  532)             DO jj = 1, jpj
(  533)                DO ji = 1, jpi
(  534)                   zlat = gphit(ji,jj)
(  535)                   zlon = MOD( glamt(ji,jj), 360. )
(  536)                   IF ( zlat1 < zlat .AND. zlat < zlat2 .AND.   &
(  537)                      45.  < zlon .AND. zlon < 100. ) THEN
(  538)                      resto(ji,jj,1) = 0.
(  539)                   ENDIF
(  540)                END DO
(  541)             END DO
(  542)          ENDIF
(  543) 
(  544)          zsdmp = 1./(sdmp * rday)
(  545)          zbdmp = 1./(bdmp * rday)
(  546)          DO jk = 2, jpkm1
(  547)             DO jj = 1, jpj
(  548)                DO ji = 1, jpi
(  549)                   zdct(ji,jj,jk) = MIN( zinfl, zdct(ji,jj,jk) )
(  550) 
(  551)                   !   ... Decrease the value in the vicinity of the coast
(  552)                   resto(ji,jj,jk) = resto(ji,jj,1)*0.5   &
(  553)                      &            * ( 1. - COS( rpi*zdct(ji,jj,jk)/zinfl) )






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 26

(  554) 
(  555)                   !   ... Vertical variation from zsdmp (sea surface) to zbdmp (bottom)
(  556)                   resto(ji,jj,jk) = resto(ji,jj,jk)   &
(  557)                      &            * ( zbdmp + (zsdmp-zbdmp)*EXP(-gdept(jk)/hdmp) )
(  558)                END DO
(  559)             END DO
(  560)          END DO
(  561) 
(  562)       ENDIF
(  563) 
(  564) 
(  565) !!DB: deleted ORCA
(  566) !      IF( cp_cfg == "orca" .AND. ( ndmp > 0 .OR. ndmp == -1 ) ) THEN
(  567) !!DB -----------------------------------------------------------------------
(  568)       if( cp_cfg == "SS" .AND. ndmp == 0 ) THEN
(  569)       
(  570)          resto(:,:,:) = 0.e0
(  571) 
(  572) !!DB 2008.06.30 
(  573) !Determine depth at which hdmp applies
(  574) !ifdef defaults to old system if key_flux_bulk_* not defined
(  575) # 576
(  576)          k_hdmp = jpk
(  577)          do jk = 1, jpk
(  578)             if(hdmp > gdept(jk)) k_hdmp = jk
(  579)          enddo
(  580) !!DBG
(  581)          if(lwp) write(numout2,*)'DBG: hdmp, k_hdmp: ', hdmp, k_hdmp 
(  582) # 585
(  585)          
(  586)          DO jk = 1, k_hdmp
(  587)             DO jj = 1, jpj
(  588)                DO ji = 1, jpi
(  589)                   resto(ji,jj,jk) = 1./(sdmp * rday)   
(  590)                END DO
(  591)             END DO
(  592)          END DO
(  593)          DO jk = k_hdmp+1, jpk
(  594)             DO jj = 1, jpj
(  595)                DO ji = 1, jpi
(  596)                   resto(ji,jj,jk) = 1./(bdmp * rday)   
(  597)                END DO
(  598)             END DO
(  599)          END DO
(  600)          
(  601)          resto(:,:, : ) = resto(:,:,:) * tmask(:,:,:)
(  602) !!DB
(  603)          resto(:,:, 1 ) = 1./(sdmp *rday) 
(  604)          
(  605)          !-------------------------------------------------------------------------------         	 
(  606)       ELSE
(  607)          !    ------------
(  608)          !     No damping
(  609)          !    ------------
(  610)          IF(lwp) WRITE(numout,cform_err)
(  611)          IF(lwp) WRITE(numout,*) 'Choose a correct value of ndmp or DO NOT defined key_tradmp'
(  612)          nstop = nstop + 1
(  613)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 27

(  614) 
(  615)       !    ----------------------------
(  616)       !     Create Print damping array
(  617)       !    ----------------------------
(  618) 
(  619)       ! ndmpf   : = 1 create a damping.coeff NetCDF file
(  620) 
(  621)       IF( ndmpf == 1 ) THEN
(  622)          IF(lwp) WRITE(numout,*) '              create damping.coeff.nc file'
(  623)          itime   = 0
(  624)          clname3 = 'damping.coeff'
(  625)          CALL ymds2ju( 0     , 1     , 1      , 0.e0 , zdate0 )
(  626)          CALL restini( 'NONE', jpi   , jpj    , glamt, gphit,    &
(  627)                        jpk   , gdept , clname3, itime, zdate0,   &
(  628)                        rdt   , idmp, domain_id=nidom )
(  629)          CALL restput( idmp, 'Resto', jpi, jpj, jpk,   &
(  630)                        0   , resto  )
(  631)          CALL restclo( idmp )
(  632)       ENDIF
(  633) 
(  634)    END SUBROUTINE dtacof
(  635) 
(  636) 
(  637)    SUBROUTINE cofdis ( pdct )
(  638)       !!----------------------------------------------------------------------
(  639)       !!                 ***  ROUTINE cofdis  ***
(  640)       !!
(  641)       !! ** Purpose :   Compute the distance between ocean T-points and the
(  642)       !!      ocean model coastlines. Save the distance in a NetCDF file.
(  643)       !!
(  644)       !! ** Method  :   For each model level, the distance-to-coast is 
(  645)       !!      computed as follows : 
(  646)       !!       - The coastline is defined as the serie of U-,V-,F-points
(  647)       !!      that are at the ocean-land bound.
(  648)       !!       - For each ocean T-point, the distance-to-coast is then 
(  649)       !!      computed as the smallest distance (on the sphere) between the 
(  650)       !!      T-point and all the coastline points.
(  651)       !!       - For land T-points, the distance-to-coast is set to zero.
(  652)       !!      C A U T I O N : Computation not yet implemented in mpp case.
(  653)       !!
(  654)       !! ** Action  : - pdct, distance to the coastline (argument)
(  655)       !!              - NetCDF file 'dist.coast.nc' 
(  656)       !!        
(  657)       !! History :
(  658)       !!   7.0  !  01-02  (M. Imbard)  Original code
(  659)       !!   8.1  !  01-02  (G. Madec, E. Durand)
(  660)       !!   8.5  !  02-08  (G. Madec, E. Durand)  Free form, F90
(  661)       !!----------------------------------------------------------------------
(  662)       !! * Modules used
(  663)       USE ioipsl
(  664) 
(  665)       !! * Arguments
(  666)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  667)          pdct                     ! distance to the coastline
(  668) 
(  669)       !! * local declarations
(  670)       INTEGER :: ji, jj, jk, jl      ! dummy loop indices
(  671)       INTEGER :: iju, ijt            ! temporary integers






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 28

(  672)       INTEGER :: icoast, itime
(  673)       INTEGER ::   &
(  674)          icot         ! logical unit for file distance to the coast
(  675)       LOGICAL, DIMENSION(jpi,jpj) ::   &
(  676)          llcotu, llcotv, llcotf   ! ???
(  677)       CHARACTER (len=32) ::   clname
(  678)       REAL(wp) ::   zdate0
(  679)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  680)          zxt, zyt, zzt,                 &  ! cartesian coordinates for T-points
(  681)          zmask                             
(  682)       REAL(wp), DIMENSION(3*jpi*jpj) ::   &
(  683)          zxc, zyc, zzc, zdis      ! temporary workspace
(  684)       !!----------------------------------------------------------------------
(  685) 
(  686)       ! 0. Initialization
(  687)       ! -----------------
(  688)       IF(lwp) WRITE(numout,*)
(  689)       IF(lwp) WRITE(numout,*) 'cofdis : compute the distance to coastline'
(  690)       IF(lwp) WRITE(numout,*) '~~~~~~'
(  691)       IF(lwp) WRITE(numout,*)
(  692)       IF( lk_mpp ) THEN
(  693)          IF(lwp) WRITE(numout,cform_err)
(  694)          IF(lwp) WRITE(numout,*) '         Computation not yet implemented with key_mpp_...'
(  695)          IF(lwp) WRITE(numout,*) '         Rerun the code on another computer or '
(  696)          IF(lwp) WRITE(numout,*) '         create the "dist.coast.nc" file using IDL'
(  697)          nstop = nstop + 1
(  698)       ENDIF
(  699) 
(  700)       pdct(:,:,:) = 0.e0
(  701)       zxt(:,:) = cos( rad * gphit(:,:) ) * cos( rad * glamt(:,:) )
(  702)       zyt(:,:) = cos( rad * gphit(:,:) ) * sin( rad * glamt(:,:) )
(  703)       zzt(:,:) = sin( rad * gphit(:,:) )
(  704) 
(  705) 
(  706)       ! 1. Loop on vertical levels
(  707)       ! --------------------------
(  708)       !                                                ! ===============
(  709)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  710)          !                                             ! ===============
(  711)          ! Define the coastline points (U, V and F)
(  712)          DO jj = 2, jpjm1
(  713)             DO ji = 2, jpim1
(  714)                zmask(ji,jj) =  ( tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) &
(  715)                    &           + tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk) )
(  716)                llcotu(ji,jj) = ( tmask(ji,jj,  jk) + tmask(ji+1,jj  ,jk) == 1. ) 
(  717)                llcotv(ji,jj) = ( tmask(ji,jj  ,jk) + tmask(ji  ,jj+1,jk) == 1. ) 
(  718)                llcotf(ji,jj) = ( zmask(ji,jj) > 0. ) .AND. ( zmask(ji,jj) < 4. )
(  719)             END DO
(  720)          END DO
(  721) 
(  722)          ! Lateral boundaries conditions
(  723)          llcotu(:, 1 ) = umask(:,  2  ,jk) == 1
(  724)          llcotu(:,jpj) = umask(:,jpjm1,jk) == 1
(  725)          llcotv(:, 1 ) = vmask(:,  2  ,jk) == 1
(  726)          llcotv(:,jpj) = vmask(:,jpjm1,jk) == 1
(  727)          llcotf(:, 1 ) = fmask(:,  2  ,jk) == 1
(  728)          llcotf(:,jpj) = fmask(:,jpjm1,jk) == 1
(  729) 






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 29

(  730)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  731)             llcotu( 1 ,:) = llcotu(jpim1,:)
(  732)             llcotu(jpi,:) = llcotu(  2  ,:)
(  733)             llcotv( 1 ,:) = llcotv(jpim1,:)
(  734)             llcotv(jpi,:) = llcotv(  2  ,:)
(  735)             llcotf( 1 ,:) = llcotf(jpim1,:)
(  736)             llcotf(jpi,:) = llcotf(  2  ,:)
(  737)          ELSE
(  738)             llcotu( 1 ,:) = umask(  2  ,:,jk) == 1
(  739)             llcotu(jpi,:) = umask(jpim1,:,jk) == 1
(  740)             llcotv( 1 ,:) = vmask(  2  ,:,jk) == 1
(  741)             llcotv(jpi,:) = vmask(jpim1,:,jk) == 1
(  742)             llcotf( 1 ,:) = fmask(  2  ,:,jk) == 1
(  743)             llcotf(jpi,:) = fmask(jpim1,:,jk) == 1
(  744)          ENDIF
(  745)          IF( nperio == 3 .OR. nperio == 4 ) THEN
(  746)             DO ji = 1, jpim1
(  747)                iju = jpi - ji + 1
(  748)                llcotu(ji,jpj  ) = llcotu(iju,jpj-2)
(  749)                llcotf(ji,jpj-1) = llcotf(iju,jpj-2)
(  750)                llcotf(ji,jpj  ) = llcotf(iju,jpj-3)
(  751)             END DO
(  752)             DO ji = jpi/2, jpi-1
(  753)                iju = jpi - ji + 1
(  754)                llcotu(ji,jpjm1) = llcotu(iju,jpjm1)
(  755)             END DO
(  756)             DO ji = 2, jpi
(  757)                ijt = jpi - ji + 2
(  758)                llcotv(ji,jpj-1) = llcotv(ijt,jpj-2)
(  759)                llcotv(ji,jpj  ) = llcotv(ijt,jpj-3)
(  760)             END DO
(  761)          ENDIF
(  762)          IF( nperio == 5 .OR. nperio == 6 ) THEN
(  763)             DO ji = 1, jpim1
(  764)                iju = jpi - ji
(  765)                llcotu(ji,jpj  ) = llcotu(iju,jpj-1)
(  766)                llcotf(ji,jpj  ) = llcotf(iju,jpj-2)
(  767)             END DO
(  768)             DO ji = jpi/2, jpi-1
(  769)                iju = jpi - ji
(  770)                llcotf(ji,jpjm1) = llcotf(iju,jpjm1)
(  771)             END DO
(  772)             DO ji = 1, jpi
(  773)                ijt = jpi - ji + 1
(  774)                llcotv(ji,jpj  ) = llcotv(ijt,jpj-1)
(  775)             END DO
(  776)             DO ji = jpi/2+1, jpi
(  777)                ijt = jpi - ji + 1
(  778)                llcotv(ji,jpjm1) = llcotv(ijt,jpjm1)
(  779)             END DO
(  780)          ENDIF
(  781) 
(  782)          ! Compute cartesian coordinates of coastline points
(  783)          ! and the number of coastline points
(  784) 
(  785)          icoast = 0
(  786)          DO jj = 1, jpj
(  787)             DO ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 30

(  788)                IF( llcotf(ji,jj) ) THEN
(  789)                   icoast = icoast + 1
(  790)                   zxc(icoast) = COS( rad*gphif(ji,jj) ) * COS( rad*glamf(ji,jj) )
(  791)                   zyc(icoast) = COS( rad*gphif(ji,jj) ) * SIN( rad*glamf(ji,jj) )
(  792)                   zzc(icoast) = SIN( rad*gphif(ji,jj) )
(  793)                ENDIF
(  794)                IF( llcotu(ji,jj) ) THEN
(  795)                   icoast = icoast+1
(  796)                   zxc(icoast) = COS( rad*gphiu(ji,jj) ) * COS( rad*glamu(ji,jj) )
(  797)                   zyc(icoast) = COS( rad*gphiu(ji,jj) ) * SIN( rad*glamu(ji,jj) )
(  798)                   zzc(icoast) = SIN( rad*gphiu(ji,jj) )
(  799)                ENDIF
(  800)                IF( llcotv(ji,jj) ) THEN
(  801)                   icoast = icoast+1
(  802)                   zxc(icoast) = COS( rad*gphiv(ji,jj) ) * COS( rad*glamv(ji,jj) )
(  803)                   zyc(icoast) = COS( rad*gphiv(ji,jj) ) * SIN( rad*glamv(ji,jj) )
(  804)                   zzc(icoast) = SIN( rad*gphiv(ji,jj) )
(  805)                ENDIF
(  806)             END DO
(  807)          END DO
(  808) 
(  809)          ! Distance for the T-points
(  810) 
(  811)          DO jj = 1, jpj
(  812)             DO ji = 1, jpi
(  813)                IF( tmask(ji,jj,jk) == 0. ) THEN
(  814)                   pdct(ji,jj,jk) = 0.
(  815)                ELSE
(  816)                   DO jl = 1, icoast
(  817)                      zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
(  818)                               + ( zyt(ji,jj) - zyc(jl) )**2   &
(  819)                               + ( zzt(ji,jj) - zzc(jl) )**2
(  820)                   END DO
(  821)                   pdct(ji,jj,jk) = ra * SQRT( MINVAL( zdis(1:icoast) ) )
(  822)                ENDIF
(  823)             END DO
(  824)          END DO
(  825)          !                                                ! ===============
(  826)       END DO                                              !   End of slab
(  827)       !                                                   ! ===============
(  828) 
(  829) 
(  830)       ! 2. Create the  distance to the coast file in NetCDF format
(  831)       ! ----------------------------------------------------------    
(  832)       clname = 'dist.coast'
(  833)       itime = 0
(  834)       CALL ymds2ju( 0     , 1     , 1     , 0.e0 , zdate0 )
(  835)       CALL restini( 'NONE', jpi   , jpj   , glamt, gphit ,   &
(  836)                     jpk   , gdept , clname, itime, zdate0,   &
(  837)                     rdt   , icot                         )
(  838)       CALL restput( icot, 'Tcoast', jpi, jpj, jpk, 0, pdct )
(  839)       CALL restclo( icot )
(  840) 
(  841)    END SUBROUTINE cofdis
(  842) 
(  843) 
(  844) # 854
(  854)    !!======================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:13      page 31

(  855) END MODULE tradmp
