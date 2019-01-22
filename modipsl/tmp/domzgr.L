


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

Filename: domzgr.F90

(    1) MODULE domzgr
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domzgr   ***
(    4)    !! Ocean initialization : domain initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_zgr     : defined the ocean vertical coordinate system
(    9)    !!       zgr_bat      : bathymetry fields (levels and meters)
(   10)    !!       zgr_bat_zoom : modify the bathymetry field if zoom domain
(   11)    !!       zgr_bat_ctl  : check the bathymetry files
(   12)    !!       zgr_z        : reference z-coordinate 
(   13)    !!       zgr_zps      : z-coordinate with partial steps
(   14)    !!       zgr_s        : s-coordinate
(   15)    !!---------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distributed memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE closea
(   23)    USE solisl
(   24)    USE ini1d           ! initialization of the 1D configuration
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dom_zgr        ! called by dom_init.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "vectopt_loop_substitute.h90"
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr.F90,v 1.14 2006/04/28 12:24:19 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS       
(   42) 
(   43)    SUBROUTINE dom_zgr
(   44)       !!----------------------------------------------------------------------
(   45)       !!                ***  ROUTINE dom_zgr  ***
(   46)       !!                   






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 2

(   47)       !! ** Purpose :  set the depth of model levels and the resulting 
(   48)       !!      vertical scale factors.
(   49)       !!
(   50)       !! ** Method  reference vertical coordinate
(   51)       !!        Z-coordinates : The depth of model levels is defined
(   52)       !!      from an analytical function the derivative of which gives
(   53)       !!      the vertical scale factors.
(   54)       !!      both depth and scale factors only depend on k (1d arrays).
(   55)       !!              w-level: gdepw  = fsdep(k)
(   56)       !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
(   57)       !!              t-level: gdept  = fsdep(k+0.5)
(   58)       !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
(   59)       !!
(   60)       !! ** Action : - gdept, gdepw : depth of T- and W-point (m)
(   61)       !!             -  e3t, e3w    : scale factors at T- and W-levels (m)
(   62)       !!
(   63)       !! Reference :
(   64)       !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
(   65)       !!
(   66)       !! History :
(   67)       !!   9.0  !  03-08  (G. Madec)  original code
(   68)       !!----------------------------------------------------------------------
(   69)       INTEGER ::   ioptio = 0      ! temporary integer
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72)       ! Check Vertical coordinate options
(   73)       ! ---------------------------------
(   74)       ioptio = 0
(   75)       IF( lk_sco ) THEN
(   76)          IF(lwp) WRITE(numout,*)
(   77)          IF(lwp) WRITE(numout,*) 'dom_zgr : s-coordinate'
(   78)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   79)          ioptio = ioptio + 1
(   80)       ENDIF
(   81)       IF( lk_zps ) THEN
(   82)          IF(lwp) WRITE(numout,*)
(   83)          IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate with partial steps'
(   84)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   85)          ioptio = ioptio + 1
(   86)       ENDIF
(   87)       IF( ioptio == 0 ) THEN
(   88)          IF(lwp) WRITE(numout,*)
(   89)          IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate'
(   90)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   91)       ENDIF
(   92) 
(   93)       IF ( ioptio > 1 ) THEN
(   94)           IF(lwp) WRITE(numout,cform_err)
(   95)           IF(lwp) WRITE(numout,*) ' several vertical coordinate options used'
(   96)           nstop = nstop + 1
(   97)       ENDIF
(   98) 
(   99)       ! Build the vertical coordinate system
(  100)       ! ------------------------------------
(  101) 
(  102)       CALL zgr_z                       ! Reference z-coordinate system
(  103) 
(  104)       CALL zgr_bat                     ! Bathymetry fields (levels and meters)






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 3

(  105) 
(  106)       CALL zgr_zps                     ! Partial step z-coordinate
(  107) 
(  108)       CALL zgr_s                       ! s-coordinate
(  109) 
(  110)    END SUBROUTINE dom_zgr
(  111) 
(  112) 
(  113)    SUBROUTINE zgr_z
(  114)       !!----------------------------------------------------------------------
(  115)       !!                   ***  ROUTINE zgr_z  ***
(  116)       !!                   
(  117)       !! ** Purpose :   set the depth of model levels and the resulting 
(  118)       !!      vertical scale factors.
(  119)       !!
(  120)       !! ** Method  :   z-coordinate system (use in all type of coordinate)
(  121)       !!        The depth of model levels is defined from an analytical
(  122)       !!      function the derivative of which gives the scale factors.
(  123)       !!        both depth and scale factors only depend on k (1d arrays).
(  124)       !!              w-level: gdepw  = fsdep(k)
(  125)       !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
(  126)       !!              t-level: gdept  = fsdep(k+0.5)
(  127)       !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
(  128)       !!
(  129)       !! ** Action  : - gdept, gdepw : depth of T- and W-point (m)
(  130)       !!              -  e3t, e3w    : scale factors at T- and W-levels (m)
(  131)       !!
(  132)       !! Reference :
(  133)       !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
(  134)       !!
(  135)       !! History :
(  136)       !!   9.0  !  03-08  (G. Madec)  F90: Free form and module
(  137)       !!----------------------------------------------------------------------
(  138)       !! * Local declarations
(  139)       INTEGER  ::   jk                     ! dummy loop indices
(  140)       REAL(wp) ::   zt, zw                 ! temporary scalars
(  141)       REAL(wp) ::   &
(  142)       zsur , za0, za1, zkth, zacr,      &  ! Values set from parameters in
(  143)       zdzmin, zhmax                        ! par_CONFIG_Rxx.h90
(  144)       !!----------------------------------------------------------------------
(  145) 
(  146)       ! Set variables from parameters
(  147)       ! ------------------------------
(  148)        zkth = ppkth       ;   zacr = ppacr
(  149)        zdzmin = ppdzmin   ;   zhmax = pphmax
(  150) 
(  151)       ! If ppa1 and ppa0 and ppsur are et to pp_to_be_computed
(  152)       !  za0, za1, zsur are computed from ppdzmin , pphmax, ppkth, ppacr
(  153)       !
(  154)        IF(  ppa1  == pp_to_be_computed  .AND.  &
(  155)          &  ppa0  == pp_to_be_computed  .AND.  &
(  156)          &  ppsur == pp_to_be_computed           ) THEN
(  157)          za1 = ( ppdzmin - pphmax / FLOAT(jpk-1) )          &
(  158)              / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) &
(  159)              &                         *  (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
(  160)              &                             - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )
(  161) 
(  162)          za0  = ppdzmin - za1 * TANH( (1-ppkth) / ppacr )






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 4

(  163)          zsur = - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )
(  164) 
(  165)        ELSE
(  166)          za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
(  167)        ENDIF
(  168) 
(  169) 
(  170)       ! Parameter print
(  171)       ! ---------------
(  172)       IF(lwp) THEN
(  173)          WRITE(numout,*)
(  174)          WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates'
(  175)          WRITE(numout,*) '    ~~~~~~~'
(  176)          IF (  ppkth == 0. ) THEN              
(  177)               WRITE(numout,*) '            Uniform grid with ',jpk-1,' layers'
(  178)               WRITE(numout,*) '            Total depth    :', zhmax
(  179)               WRITE(numout,*) '            Layer thickness:', zhmax/(jpk-1)
(  180)          ELSE
(  181)             IF ( ppa1 == 0. .AND. ppa0 == 0. .AND. ppsur == 0. ) THEN
(  182)                WRITE(numout,*) '         zsur, za0, za1 computed from '
(  183)                WRITE(numout,*) '                 zdzmin = ', zdzmin
(  184)                WRITE(numout,*) '                 zhmax  = ', zhmax
(  185)             ENDIF
(  186)             WRITE(numout,*) '           Value of coefficients for vertical mesh:'
(  187)             WRITE(numout,*) '                 zsur = ', zsur
(  188)             WRITE(numout,*) '                 za0  = ', za0
(  189)             WRITE(numout,*) '                 za1  = ', za1
(  190)             WRITE(numout,*) '                 zkth = ', zkth
(  191)             WRITE(numout,*) '                 zacr = ', zacr
(  192)          ENDIF
(  193)       ENDIF
(  194) 
(  195) 
(  196)       ! Reference z-coordinate (depth - scale factor at T- and W-points)
(  197)       ! ======================
(  198)       IF (  ppkth == 0. ) THEN            !  uniform vertical grid       
(  199) 
(  200)          za1 = zhmax/FLOAT(jpk-1) 
(  201)          DO jk = 1, jpk
(  202)             zw = FLOAT( jk )
(  203)             zt = FLOAT( jk ) + 0.5
(  204)             gdepw(jk) = ( zw - 1 ) * za1
(  205)             gdept(jk) = ( zt - 1 ) * za1
(  206)             e3w  (jk) =  za1
(  207)             e3t  (jk) =  za1
(  208)          END DO
(  209) 
(  210)       ELSE
(  211) 
(  212)          DO jk = 1, jpk
(  213)             zw = FLOAT( jk )
(  214)             zt = FLOAT( jk ) + 0.5
(  215)             gdepw(jk) = ( zsur + za0 * zw + za1 * zacr * LOG( COSH( (zw-zkth)/zacr ) )  )
(  216)             gdept(jk) = ( zsur + za0 * zt + za1 * zacr * LOG( COSH( (zt-zkth)/zacr ) )  )
(  217)             e3w  (jk) =          za0      + za1        * TANH(      (zw-zkth)/zacr   )
(  218)             e3t  (jk) =          za0      + za1        * TANH(      (zt-zkth)/zacr   )
(  219)          END DO
(  220)          gdepw(1) = 0.e0   ! force first w-level to be exactly at zero






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 5

(  221) 
(  222) !!DBG: 2008.11.17: the params in par_SS008.h90 yield for the first 4 levels
(  223) !            level   gdept    gdepw     e3t      e3w  
(  224) !             1     3.05     0.00     6.19     6.00
(  225) !             2     9.45     6.19     6.64     6.40
(  226) !             3    16.36    12.84     7.20     6.90
(  227) !             4    23.90    20.04     7.89     7.53
(  228) !!Since the min #-levels = 2, this means that the min btm depth should be 12.84m
(  229) !!which we want to keep the same to start.
(  230) !!REM: e3t=thickness centered on T-pt; e3w=thickness centered on w-pt; 
(  231) !!e3w(1)=??? Pretend that e3w(1) ~ 2*gdept(1)
(  232) !!e.g. of change to a 3m thick top layer, while keeping gdepw(3) unchanged
(  233) !         gdepw(2) = 3.0 
(  234) !         gdept(1) = gdepw(2)/2.0
(  235) !         e3w(1) = 2.0*gdept(1)
(  236) !         e3t(1) = gdepw(2)-gdepw(1)
(  237) !         e3t(2) = gdepw(3)-gdepw(2) 
(  238) !         gdept(2) = gdepw(2) + e3t(2)/2.0 
(  239) !         e3w(2) = gdept(2)-gdept(1)
(  240) !         e3w(3) = gdept(3)-gdept(2)
(  241) 
(  242) 
(  243)       ENDIF
(  244) 
(  245)       ! Control and  print
(  246)       ! ==================
(  247) 
(  248)       IF(lwp) THEN
(  249)          WRITE(numout,*)
(  250)          WRITE(numout,*) '              Reference z-coordinate depth and scale factors:'
(  251)          WRITE(numout, "(9x,' level   gdept    gdepw     e3t      e3w  ')" )
(  252)          WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, gdept(jk), gdepw(jk), e3t(jk), e3w(jk), jk = 1, jpk )
(  253)       ENDIF
(  254) 
(  255)       DO jk = 1, jpk
(  256)          IF( e3w(jk) <= 0. .OR. e3t(jk) <= 0. ) THEN
(  257)             IF(lwp) WRITE(numout,cform_err)
(  258)             IF(lwp) WRITE(numout,*) ' e3w or e3t =< 0 '
(  259)             nstop = nstop + 1
(  260)          ENDIF
(  261)          IF( gdepw(jk) < 0. .OR. gdept(jk) < 0.) THEN
(  262)             IF(lwp) WRITE(numout,cform_err)
(  263)             IF(lwp) WRITE(numout,*) ' gdepw or gdept < 0 '
(  264)             nstop = nstop + 1
(  265)          ENDIF
(  266)       END DO
(  267) 
(  268)    END SUBROUTINE zgr_z
(  269) 
(  270) 
(  271)    SUBROUTINE zgr_bat
(  272)       !!----------------------------------------------------------------------
(  273)       !!                    ***  ROUTINE zgr_bat  ***
(  274)       !! 
(  275)       !! ** Purpose :   set bathymetry both in levels and meters
(  276)       !!
(  277)       !! ** Method  :   read or define mbathy and bathy arrays
(  278)       !!       * level bathymetry:






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 6

(  279)       !!      The ocean basin geometry is given by a two-dimensional array,
(  280)       !!      mbathy, which is defined as follow :
(  281)       !!            mbathy(ji,jj) = 1, ..., jpk-1, the number of ocean level
(  282)       !!                              at t-point (ji,jj).
(  283)       !!                            = 0  over the continental t-point.
(  284)       !!                            = -n over the nth island t-point.
(  285)       !!      The array mbathy is checked to verified its consistency with
(  286)       !!      model option. in particular:
(  287)       !!            mbathy must have at least 1 land grid-points (mbathy<=0)
(  288)       !!                  along closed boundary.
(  289)       !!            mbathy must be cyclic IF jperio=1.
(  290)       !!            mbathy must be lower or equal to jpk-1.
(  291)       !!            isolated ocean grid points are suppressed from mbathy
(  292)       !!                  since they are only connected to remaining
(  293)       !!                  ocean through vertical diffusion.
(  294)       !!      ntopo=-1 :   rectangular channel or bassin with a bump 
(  295)       !!      ntopo= 0 :   flat rectangular channel or basin 
(  296)       !!      ntopo= 1 :   mbathy is read in 'bathy_level.nc' NetCDF file
(  297)       !!                   bathy  is read in 'bathy_meter.nc' NetCDF file
(  298)       !!      C A U T I O N : mbathy will be modified during the initializa-
(  299)       !!      tion phase to become the number of non-zero w-levels of a water
(  300)       !!      column, with a minimum value of 1.
(  301)       !!
(  302)       !! ** Action  : - mbathy: level bathymetry (in level index)
(  303)       !!              - bathy : meter bathymetry (in meters)
(  304)       !!
(  305)       !! History :
(  306)       !!   9.0  !  03-08  (G. Madec)  Original code
(  307)       !!----------------------------------------------------------------------
(  308)       !! * Modules used
(  309)       USE ioipsl
(  310) 
(  311)       !! * Local declarations
(  312)       CHARACTER (len=18) ::   clname    ! temporary characters
(  313)       LOGICAL ::   llbon                ! check the existence of bathy files
(  314)       INTEGER ::   ji, jj, jl, jk       ! dummy loop indices
(  315)       INTEGER ::   inum = 11            ! temporary logical unit
(  316)       INTEGER  ::   &
(  317)          ipi, ipj, ipk,              &  ! temporary integers
(  318)          itime,                      &  !    "          "
(  319)          ii_bump, ij_bump               ! bump center position
(  320)       INTEGER, DIMENSION (1) ::   istep
(  321)       INTEGER , DIMENSION(jpidta,jpjdta) ::   &
(  322)          idta                           ! global domain integer data
(  323)       REAL(wp) ::   &
(  324)          r_bump, h_bump, h_oce,      &  ! bump characteristics 
(  325)          zi, zj, zdate0, zdt            ! temporary scalars
(  326)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  327)          zlamt, zphit,               &  ! temporary workspace (NetCDF read)
(  328)          zdta                           ! global domain scalar data
(  329)       REAL(wp), DIMENSION(jpk) ::   &
(  330)          zdept                          ! temporary workspace (NetCDF read)
(  331)       !!----------------------------------------------------------------------
(  332) 
(  333)       IF(lwp) WRITE(numout,*)
(  334)       IF(lwp) WRITE(numout,*) '    zgr_bat : defines level and meter bathymetry'
(  335)       IF(lwp) WRITE(numout,*) '    ~~~~~~~'
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 7

(  337)       ! ========================================
(  338)       ! global domain level and meter bathymetry (idta,zdta)
(  339)       ! ========================================
(  340)       !                                               ! =============== ! 
(  341)       IF( ntopo == 0 .OR. ntopo == -1 ) THEN          ! defined by hand !
(  342)          !                                            ! =============== !
(  343) 
(  344)          IF( ntopo == 0 ) THEN                        ! flat basin
(  345) 
(  346)             IF(lwp) WRITE(numout,*)
(  347)             IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin'
(  348) 
(  349)             idta(:,:) = jpkm1                            ! flat basin 
(  350)             zdta(:,:) = gdepw(jpk)
(  351) 
(  352)          ELSE                                         ! bump
(  353)             IF(lwp) WRITE(numout,*)
(  354)             IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin with a bump'
(  355) 
(  356)             ii_bump = jpidta / 3 + 3       ! i-index of the bump center
(  357)             ij_bump = jpjdta / 2           ! j-index of the bump center
(  358)             r_bump  =    6              ! bump radius (index)       
(  359)             h_bump  =  240.e0           ! bump height (meters)
(  360)             h_oce   = gdepw(jpk)        ! background ocean depth (meters)
(  361)             IF(lwp) WRITE(numout,*) '            bump characteristics: '
(  362)             IF(lwp) WRITE(numout,*) '               bump center (i,j)   = ', ii_bump, ii_bump
(  363)             IF(lwp) WRITE(numout,*) '               bump height         = ', h_bump , ' meters'
(  364)             IF(lwp) WRITE(numout,*) '               bump radius         = ', r_bump , ' index'
(  365)             IF(lwp) WRITE(numout,*) '            background ocean depth = ', h_oce  , ' meters'
(  366)             ! zdta :
(  367)             DO jj = 1, jpjdta
(  368)                DO ji = 1, jpidta
(  369)                   zi = FLOAT( ji - ii_bump ) / r_bump      
(  370)                   zj = FLOAT( jj - ij_bump ) / r_bump       
(  371)                   zdta(ji,jj) = h_oce - h_bump * EXP( -( zi*zi + zj*zj ) )
(  372)                END DO
(  373)             END DO
(  374)             ! idta :
(  375)             idta(:,:) = jpkm1
(  376)             DO jk = 1, jpkm1
(  377)                DO jj = 1, jpjdta
(  378)                   DO ji = 1, jpidta
(  379)                      IF( gdept(jk) < zdta(ji,jj) .AND. zdta(ji,jj) <= gdept(jk+1) )   idta(ji,jj) = jk
(  380)                   END DO
(  381)                END DO
(  382)             END DO
(  383)          ENDIF
(  384) 
(  385)          ! set boundary conditions (caution, idta on the global domain: use of jperio, not nperio)
(  386)          IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
(  387)             idta( :    , 1    ) = -1                ;      zdta( :    , 1    ) = -1.e0
(  388)             idta( :    ,jpjdta) =  0                ;      zdta( :    ,jpjdta) =  0.e0
(  389)          ELSEIF( jperio == 2 ) THEN
(  390)             idta( :    , 1    ) = idta( : ,  3  )   ;      zdta( :    , 1    ) = zdta( : ,  3  )
(  391)             idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
(  392)             idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
(  393)             idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
(  394)          ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 8

(  395)             idta( :    , 1    ) = 0                 ;      zdta( :    , 1    ) =  0.e0
(  396)             idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
(  397)             idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
(  398)             idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
(  399)          ENDIF
(  400) 
(  401)          !  EEL R5 configuration with east and west open boundaries.
(  402)          !  Two rows of zeroes are needed at the south and north for OBCs
(  403)          !  This is for compatibility with the rigid lid option. 
(  404)           
(  405) !         IF( cp_cfg == "eel" .AND. jp_cfg == 5 ) THEN
(  406) !            idta( : , 2      ) = 0                 ;      zdta( : , 2      ) =  0.e0
(  407) !            idta( : ,jpjdta-1) = 0                 ;      zdta( : ,jpjdta-1) =  0.e0
(  408) !         ENDIF
(  409) 
(  410)          !                                            ! =============== !
(  411)       ELSEIF( ntopo == 1 ) THEN                       !   read in file  !
(  412)          !                                            ! =============== !
(  413)          IF( lk_zco ) THEN
(  414)             clname = 'bathy_level.nc'                       ! Level bathymetry
(  415) #if defined key_agrif
(  416)             IF( .NOT. Agrif_Root() ) THEN
(  417)                clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  418)             ENDIF
(  419) #endif         
(  420)             INQUIRE( FILE=clname, EXIST=llbon )
(  421)             IF( llbon ) THEN
(  422)                IF(lwp) WRITE(numout,*)
(  423)                IF(lwp) WRITE(numout,*) '         read level bathymetry in ', clname
(  424)                IF(lwp) WRITE(numout,*)
(  425)                itime = 1
(  426)                ipi = jpidta
(  427)                ipj = jpjdta
(  428)                ipk = 1
(  429)                zdt = rdt
(  430)                CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &
(  431)                               ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
(  432)                CALL flinget( inum, 'Bathy_level', jpidta, jpjdta, 1,   &
(  433)                              itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) )
(  434)                idta(:,:) = zdta(:,:)
(  435)                CALL flinclo( inum )
(  436) 
(  437)             ELSE
(  438)                IF(lwp) WRITE(numout,cform_err)
(  439)                IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname
(  440)                nstop = nstop + 1
(  441)             ENDIF   
(  442)   
(  443)          ELSEIF( lk_zps ) THEN
(  444)             clname = 'bathy_meter.nc'                       ! meter bathymetry
(  445) #if defined key_agrif
(  446)             IF( .NOT. Agrif_Root() ) THEN
(  447)                clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  448)             ENDIF
(  449) #endif   	 
(  450)             INQUIRE( FILE=clname, EXIST=llbon )
(  451)             IF( llbon ) THEN
(  452)                IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 9

(  453)                IF(lwp) WRITE(numout,*) '         read meter bathymetry in ', clname
(  454)                IF(lwp) WRITE(numout,*)
(  455)                itime = 1
(  456)                ipi = jpidta
(  457)                ipj = jpjdta
(  458)                ipk = 1
(  459)                zdt = rdt
(  460)                CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &    
(  461)                               ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
(  462)                CALL flinget( inum, 'Bathymetry', jpidta, jpjdta, 1,   &
(  463)                              itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) ) 
(  464)                CALL flinclo( inum )
(  465)                idta(:,:) = jpkm1      ! initialisation
(  466)             ELSE
(  467)                IF(lwp) WRITE(numout,cform_err)       
(  468)                IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname
(  469)                nstop = nstop + 1
(  470)             ENDIF
(  471)          ENDIF
(  472)          !                                            ! =============== !
(  473)       ELSE                                            !      error      !
(  474)          !                                            ! =============== !
(  475)          IF(lwp) WRITE(numout,cform_err)
(  476)          IF(lwp) WRITE(numout,*) '          parameter , ntopo = ', ntopo
(  477)          nstop = nstop + 1
(  478)       ENDIF
(  479) 
(  480) 
(  481)       ! =======================================
(  482)       ! local domain level and meter bathymetry (mbathy,bathy)
(  483)       ! =======================================
(  484) 
(  485)       mbathy(:,:) = 0                                 ! set to zero extra halo points
(  486)       bathy (:,:) = 0.e0                              ! (require for mpp case)
(  487) 
(  488)       DO jj = 1, nlcj                                 ! interior values
(  489)          DO ji = 1, nlci
(  490)             mbathy(ji,jj) = idta( mig(ji), mjg(jj) )
(  491)             bathy (ji,jj) = zdta( mig(ji), mjg(jj) )
(  492)          END DO
(  493)       END DO
(  494) 
(  495)       ! =======================
(  496)       ! NO closed seas or lakes
(  497)       ! =======================
(  498) 
(  499)       IF( nclosea == 0 ) THEN
(  500)          DO jl = 1, jpncs
(  501)             DO jj = ncsj1(jl), ncsj2(jl)
(  502)                DO ji = ncsi1(jl), ncsi2(jl)
(  503)                   mbathy(ji,jj) = 0                   ! suppress closed seas
(  504)                   bathy (ji,jj) = 0.e0                ! and lakes
(  505)                END DO
(  506)             END DO
(  507)          END DO
(  508)       ENDIF
(  509) 
(  510)       ! ===========






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 10

(  511)       ! Zoom domain 
(  512)       ! ===========
(  513) 
(  514)       IF( lzoom )   CALL zgr_bat_zoom
(  515) 
(  516)       ! ================
(  517)       ! Bathymetry check
(  518)       ! ================
(  519) 
(  520)       CALL zgr_bat_ctl
(  521) 
(  522)    END SUBROUTINE zgr_bat
(  523) 
(  524) 
(  525)    SUBROUTINE zgr_bat_zoom
(  526)       !!----------------------------------------------------------------------
(  527)       !!                    ***  ROUTINE zgr_bat_zoom  ***
(  528)       !!
(  529)       !! ** Purpose : - Close zoom domain boundary if necessary
(  530)       !!              - Suppress Med Sea from ORCA R2 and R05 arctic zoom
(  531)       !!
(  532)       !! ** Method  : 
(  533)       !!
(  534)       !! ** Action  : - update mbathy: level bathymetry (in level index)
(  535)       !!
(  536)       !! History :
(  537)       !!   9.0  !  03-08  (G. Madec)  Original code
(  538)       !!----------------------------------------------------------------------
(  539)       !! * Local variables
(  540)       INTEGER ::   ii0, ii1, ij0, ij1   ! temporary integers
(  541)       !!----------------------------------------------------------------------
(  542) 
(  543)       IF(lwp) WRITE(numout,*)
(  544)       IF(lwp) WRITE(numout,*) '    zgr_bat_zoom : modify the level bathymetry for zoom domain'
(  545)       IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~'
(  546) 
(  547)       ! Zoom domain
(  548)       ! ===========
(  549) 
(  550)       ! Forced closed boundary if required
(  551)       IF( lzoom_w )   mbathy( mi0(jpizoom):mi1(jpizoom) , :  ) = 0
(  552)       IF( lzoom_s )   mbathy(  : , mj0(jpjzoom):mj1(jpjzoom) ) = 0
(  553)       IF( lzoom_e )   mbathy( mi0(jpiglo+jpizoom-1):mi1(jpiglo+jpizoom-1) , :  ) = 0
(  554)       IF( lzoom_n )   mbathy(  : , mj0(jpjglo+jpjzoom-1):mj1(jpjglo+jpjzoom-1) ) = 0
(  555) 
(  556) !!DB: delete ORCA
(  557)       ! Configuration specific domain modifications
(  558)       ! (here, ORCA arctic configuration: suppress Med Sea)
(  559) !      IF( cp_cfg == "orca" .AND. lzoom_arct ) THEN
(  560) 
(  561)    END SUBROUTINE zgr_bat_zoom
(  562) 
(  563) 
(  564)    SUBROUTINE zgr_bat_ctl
(  565)       !!----------------------------------------------------------------------
(  566)       !!                    ***  ROUTINE zgr_bat_ctl  ***
(  567)       !!
(  568)       !! ** Purpose :   check the bathymetry in levels






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 11

(  569)       !!
(  570)       !! ** Method  :   The array mbathy is checked to verified its consistency
(  571)       !!      with the model options. in particular:
(  572)       !!            mbathy must have at least 1 land grid-points (mbathy<=0)
(  573)       !!                  along closed boundary.
(  574)       !!            mbathy must be cyclic IF jperio=1.
(  575)       !!            mbathy must be lower or equal to jpk-1.
(  576)       !!            isolated ocean grid points are suppressed from mbathy
(  577)       !!                  since they are only connected to remaining
(  578)       !!                  ocean through vertical diffusion.
(  579)       !!      C A U T I O N : mbathy will be modified during the initializa-
(  580)       !!      tion phase to become the number of non-zero w-levels of a water
(  581)       !!      column, with a minimum value of 1.
(  582)       !!
(  583)       !! ** Action  : - update mbathy: level bathymetry (in level index)
(  584)       !!              - update bathy : meter bathymetry (in meters)
(  585)       !!
(  586)       !! History :
(  587)       !!   9.0  !  03-08  (G. Madec)  Original code
(  588)       !!----------------------------------------------------------------------
(  589)       !! * Local declarations
(  590)       INTEGER ::   ji, jj, jl           ! dummy loop indices
(  591)       INTEGER ::   &
(  592)          icompt, ibtest, ikmax          ! temporary integers
(  593)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  594)          zbathy                         ! temporary workspace
(  595)       !!----------------------------------------------------------------------
(  596) 
(  597)       IF(lwp) WRITE(numout,*)
(  598)       IF(lwp) WRITE(numout,*) '    zgr_bat_ctl : check the bathymetry'
(  599)       IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
(  600) 
(  601)       ! ================
(  602)       ! Bathymetry check
(  603)       ! ================
(  604) 
(  605)       IF( .NOT. lk_cfg_1d )   THEN
(  606) 
(  607)          ! Suppress isolated ocean grid points
(  608) 
(  609)          IF(lwp) WRITE(numout,*)
(  610)          IF(lwp) WRITE(numout,*)'                   suppress isolated ocean grid points'
(  611)          IF(lwp) WRITE(numout,*)'                   -----------------------------------'
(  612) 
(  613)          icompt = 0
(  614) 
(  615)          DO jl = 1, 2
(  616) 
(  617)             IF( nperio == 1 .OR. nperio  ==  4 .OR. nperio  ==  6 ) THEN
(  618)                mbathy( 1 ,:) = mbathy(jpim1,:)
(  619)                mbathy(jpi,:) = mbathy(  2  ,:)
(  620)             ENDIF
(  621)             DO jj = 2, jpjm1
(  622)                DO ji = 2, jpim1
(  623)                   ibtest = MAX( mbathy(ji-1,jj), mbathy(ji+1,jj),   &
(  624)                      mbathy(ji,jj-1),mbathy(ji,jj+1) )
(  625)                   IF( ibtest < mbathy(ji,jj) ) THEN
(  626)                      IF(lwp) WRITE(numout,*) ' the number of ocean level at ',   &






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 12

(  627)                         'grid-point (i,j) =  ',ji,jj,' is changed from ',   &
(  628)                         mbathy(ji,jj),' to ', ibtest
(  629)                      mbathy(ji,jj) = ibtest
(  630)                      icompt = icompt + 1
(  631)                   ENDIF
(  632)                END DO
(  633)             END DO
(  634) 
(  635)          END DO
(  636)          IF( icompt == 0 ) THEN
(  637)             IF(lwp) WRITE(numout,*)'     no isolated ocean grid points'
(  638)          ELSE
(  639)             IF(lwp) WRITE(numout,*)'    ',icompt,' ocean grid points suppressed'
(  640)          ENDIF
(  641)          IF( lk_mpp ) THEN
(  642)             zbathy(:,:) = FLOAT( mbathy(:,:) )
(  643)             CALL lbc_lnk( zbathy, 'T', 1. )
(  644)             mbathy(:,:) = INT( zbathy(:,:) )
(  645)          ENDIF
(  646) 
(  647)          ! 3.2 East-west cyclic boundary conditions
(  648) 
(  649)          IF( nperio == 0 ) THEN
(  650)             IF(lwp) WRITE(numout,*) ' mbathy set to 0 along east and west',   &
(  651)                ' boundary: nperio = ', nperio
(  652)             IF( lk_mpp ) THEN
(  653)                IF( nbondi == -1 .OR. nbondi == 2 ) THEN
(  654)                   IF( jperio /= 1 )   mbathy(1,:) = 0
(  655)                ENDIF
(  656)                IF( nbondi == 1 .OR. nbondi == 2 ) THEN
(  657)                   IF( jperio /= 1 )   mbathy(nlci,:) = 0
(  658)                ENDIF
(  659)             ELSE
(  660)                mbathy( 1 ,:) = 0
(  661)                mbathy(jpi,:) = 0
(  662)             ENDIF
(  663)          ELSEIF( nperio == 1 .OR. nperio == 4 .OR. nperio ==  6 ) THEN
(  664)             IF(lwp) WRITE(numout,*)' east-west cyclic boundary conditions',   &
(  665)                ' on mbathy: nperio = ', nperio
(  666)             mbathy( 1 ,:) = mbathy(jpim1,:)
(  667)             mbathy(jpi,:) = mbathy(  2  ,:)
(  668)          ELSEIF( nperio == 2 ) THEN
(  669)             IF(lwp) WRITE(numout,*) '   equatorial boundary conditions',   &
(  670)                ' on mbathy: nperio = ', nperio
(  671)          ELSE
(  672)             IF(lwp) WRITE(numout,*) '    e r r o r'
(  673)             IF(lwp) WRITE(numout,*) '    parameter , nperio = ', nperio
(  674)             !         STOP 'dom_mba'
(  675)          ENDIF
(  676) 
(  677)          ! Set to zero mbathy over islands if necessary  (lk_isl=F)
(  678)          IF( .NOT. lk_isl ) THEN    ! No island
(  679)             IF(lwp) WRITE(numout,*)
(  680)             IF(lwp) WRITE(numout,*) '         mbathy set to 0 over islands'
(  681)             IF(lwp) WRITE(numout,*) '         ----------------------------'
(  682) 
(  683)             mbathy(:,:) = MAX( 0, mbathy(:,:) )
(  684) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 13

(  685)             !  Boundary condition on mbathy
(  686)             IF( .NOT.lk_mpp ) THEN 
(  687)                
(  688)                !!bug ???  y reflechir!
(  689)                !   ... mono- or macro-tasking: T-point, >0, 2D array, no slab
(  690)                
(  691)                zbathy(:,:) = FLOAT( mbathy(:,:) )
(  692)                CALL lbc_lnk( zbathy, 'T', 1. )
(  693)                mbathy(:,:) = INT( zbathy(:,:) )
(  694)             ENDIF
(  695) 
(  696)          ENDIF
(  697) 
(  698)       ENDIF
(  699) 
(  700)       ! Number of ocean level inferior or equal to jpkm1
(  701) 
(  702)       ikmax = 0
(  703)       DO jj = 1, jpj
(  704)          DO ji = 1, jpi
(  705)             ikmax = MAX( ikmax, mbathy(ji,jj) )
(  706)          END DO
(  707)       END DO
(  708)       !!! test a faire:   ikmax = MAX( mbathy(:,:) )   ???
(  709) 
(  710)       IF( ikmax > jpkm1 ) THEN
(  711)          IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' >  jpk-1'
(  712)          IF(lwp) WRITE(numout,*) ' change jpk to ',ikmax+1,' to use the exact ead bathymetry'
(  713)       ELSE IF( ikmax < jpkm1 ) THEN
(  714)          IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' < jpk-1' 
(  715)          IF(lwp) WRITE(numout,*) ' you can decrease jpk to ', ikmax+1
(  716)       ENDIF
(  717) 
(  718)       IF( lwp .AND. nprint == 1 ) THEN
(  719)          WRITE(numout,*)
(  720)          WRITE(numout,*) ' bathymetric field '
(  721)          WRITE(numout,*) ' ------------------'
(  722)          WRITE(numout,*) ' number of non-zero T-levels '
(  723)          CALL prihin( mbathy, jpi, jpj, 1, jpi,   &
(  724)                       1     , 1  , jpj, 1, 3  ,   &
(  725)                       numout )
(  726)          WRITE(numout,*)
(  727)       ENDIF
(  728) 
(  729)    END SUBROUTINE zgr_bat_ctl
(  730) 
(  731) 
(  732) #  include "domzgr_zps.h90"
(  733) 
(  734) 
(  735) #  include "domzgr_s.h90"
(  736) 
(  737) 
(  738)    !!======================================================================
(  739) END MODULE domzgr









PGF90 (Version     10.5)          09/23/2010  09:49:34      page 14

(    1) # 1 "domzgr.F90"
(    1) MODULE domzgr
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE domzgr   ***
(    4)    !! Ocean initialization : domain initialization
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_zgr     : defined the ocean vertical coordinate system
(    9)    !!       zgr_bat      : bathymetry fields (levels and meters)
(   10)    !!       zgr_bat_zoom : modify the bathymetry field if zoom domain
(   11)    !!       zgr_bat_ctl  : check the bathymetry files
(   12)    !!       zgr_z        : reference z-coordinate 
(   13)    !!       zgr_zps      : z-coordinate with partial steps
(   14)    !!       zgr_s        : s-coordinate
(   15)    !!---------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distributed memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE closea
(   23)    USE solisl
(   24)    USE ini1d           ! initialization of the 1D configuration
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dom_zgr        ! called by dom_init.F90
(   31) 
(   32)    !! * Substitutions
(   33) # 1 "./domzgr_substitute.h90"
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
(   48) # 34 "domzgr.F90"
(   34) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 15

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 35 "domzgr.F90"
(   35) # 35
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr.F90,v 1.14 2006/04/28 12:24:19 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS       
(   42) 
(   43)    SUBROUTINE dom_zgr
(   44)       !!----------------------------------------------------------------------
(   45)       !!                ***  ROUTINE dom_zgr  ***
(   46)       !!                   
(   47)       !! ** Purpose :  set the depth of model levels and the resulting 
(   48)       !!      vertical scale factors.
(   49)       !!
(   50)       !! ** Method  reference vertical coordinate
(   51)       !!        Z-coordinates : The depth of model levels is defined
(   52)       !!      from an analytical function the derivative of which gives
(   53)       !!      the vertical scale factors.
(   54)       !!      both depth and scale factors only depend on k (1d arrays).
(   55)       !!              w-level: gdepw  = fsdep(k)
(   56)       !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
(   57)       !!              t-level: gdept  = fsdep(k+0.5)
(   58)       !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
(   59)       !!
(   60)       !! ** Action : - gdept, gdepw : depth of T- and W-point (m)
(   61)       !!             -  e3t, e3w    : scale factors at T- and W-levels (m)
(   62)       !!
(   63)       !! Reference :
(   64)       !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
(   65)       !!
(   66)       !! History :
(   67)       !!   9.0  !  03-08  (G. Madec)  original code
(   68)       !!----------------------------------------------------------------------
(   69)       INTEGER ::   ioptio = 0      ! temporary integer
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72)       ! Check Vertical coordinate options
(   73)       ! ---------------------------------
(   74)       ioptio = 0
(   75)       IF( lk_sco ) THEN
(   76)          IF(lwp) WRITE(numout,*)
(   77)          IF(lwp) WRITE(numout,*) 'dom_zgr : s-coordinate'
(   78)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   79)          ioptio = ioptio + 1
(   80)       ENDIF
(   81)       IF( lk_zps ) THEN
(   82)          IF(lwp) WRITE(numout,*)
(   83)          IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate with partial steps'
(   84)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   85)          ioptio = ioptio + 1
(   86)       ENDIF
(   87)       IF( ioptio == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 16

(   88)          IF(lwp) WRITE(numout,*)
(   89)          IF(lwp) WRITE(numout,*) 'dom_zgr : z-coordinate'
(   90)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   91)       ENDIF
(   92) 
(   93)       IF ( ioptio > 1 ) THEN
(   94)           IF(lwp) WRITE(numout,cform_err)
(   95)           IF(lwp) WRITE(numout,*) ' several vertical coordinate options used'
(   96)           nstop = nstop + 1
(   97)       ENDIF
(   98) 
(   99)       ! Build the vertical coordinate system
(  100)       ! ------------------------------------
(  101) 
(  102)       CALL zgr_z                       ! Reference z-coordinate system
(  103) 
(  104)       CALL zgr_bat                     ! Bathymetry fields (levels and meters)
(  105) 
(  106)       CALL zgr_zps                     ! Partial step z-coordinate
(  107) 
(  108)       CALL zgr_s                       ! s-coordinate
(  109) 
(  110)    END SUBROUTINE dom_zgr
(  111) 
(  112) 
(  113)    SUBROUTINE zgr_z
(  114)       !!----------------------------------------------------------------------
(  115)       !!                   ***  ROUTINE zgr_z  ***
(  116)       !!                   
(  117)       !! ** Purpose :   set the depth of model levels and the resulting 
(  118)       !!      vertical scale factors.
(  119)       !!
(  120)       !! ** Method  :   z-coordinate system (use in all type of coordinate)
(  121)       !!        The depth of model levels is defined from an analytical
(  122)       !!      function the derivative of which gives the scale factors.
(  123)       !!        both depth and scale factors only depend on k (1d arrays).
(  124)       !!              w-level: gdepw  = fsdep(k)
(  125)       !!                       e3w(k) = dk(fsdep)(k)     = fse3(k)
(  126)       !!              t-level: gdept  = fsdep(k+0.5)
(  127)       !!                       e3t(k) = dk(fsdep)(k+0.5) = fse3(k+0.5)
(  128)       !!
(  129)       !! ** Action  : - gdept, gdepw : depth of T- and W-point (m)
(  130)       !!              -  e3t, e3w    : scale factors at T- and W-levels (m)
(  131)       !!
(  132)       !! Reference :
(  133)       !!      Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
(  134)       !!
(  135)       !! History :
(  136)       !!   9.0  !  03-08  (G. Madec)  F90: Free form and module
(  137)       !!----------------------------------------------------------------------
(  138)       !! * Local declarations
(  139)       INTEGER  ::   jk                     ! dummy loop indices
(  140)       REAL(wp) ::   zt, zw                 ! temporary scalars
(  141)       REAL(wp) ::   &
(  142)       zsur , za0, za1, zkth, zacr,      &  ! Values set from parameters in
(  143)       zdzmin, zhmax                        ! par_CONFIG_Rxx.h90
(  144)       !!----------------------------------------------------------------------
(  145) 






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 17

(  146)       ! Set variables from parameters
(  147)       ! ------------------------------
(  148)        zkth = ppkth       ;   zacr = ppacr
(  149)        zdzmin = ppdzmin   ;   zhmax = pphmax
(  150) 
(  151)       ! If ppa1 and ppa0 and ppsur are et to pp_to_be_computed
(  152)       !  za0, za1, zsur are computed from ppdzmin , pphmax, ppkth, ppacr
(  153)       !
(  154)        IF(  ppa1  == pp_to_be_computed  .AND.  &
(  155)          &  ppa0  == pp_to_be_computed  .AND.  &
(  156)          &  ppsur == pp_to_be_computed           ) THEN
(  157)          za1 = ( ppdzmin - pphmax / FLOAT(jpk-1) )          &
(  158)              / ( TANH((1-ppkth)/ppacr) - ppacr/FLOAT(jpk-1) &
(  159)              &                         *  (  LOG( COSH( (jpk - ppkth) / ppacr) )      &
(  160)              &                             - LOG( COSH( ( 1  - ppkth) / ppacr) )  )  )
(  161) 
(  162)          za0  = ppdzmin - za1 * TANH( (1-ppkth) / ppacr )
(  163)          zsur = - za0 - za1 * ppacr * LOG( COSH( (1-ppkth) / ppacr )  )
(  164) 
(  165)        ELSE
(  166)          za1 = ppa1 ;       za0 = ppa0 ;          zsur = ppsur
(  167)        ENDIF
(  168) 
(  169) 
(  170)       ! Parameter print
(  171)       ! ---------------
(  172)       IF(lwp) THEN
(  173)          WRITE(numout,*)
(  174)          WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates'
(  175)          WRITE(numout,*) '    ~~~~~~~'
(  176)          IF (  ppkth == 0. ) THEN              
(  177)               WRITE(numout,*) '            Uniform grid with ',jpk-1,' layers'
(  178)               WRITE(numout,*) '            Total depth    :', zhmax
(  179)               WRITE(numout,*) '            Layer thickness:', zhmax/(jpk-1)
(  180)          ELSE
(  181)             IF ( ppa1 == 0. .AND. ppa0 == 0. .AND. ppsur == 0. ) THEN
(  182)                WRITE(numout,*) '         zsur, za0, za1 computed from '
(  183)                WRITE(numout,*) '                 zdzmin = ', zdzmin
(  184)                WRITE(numout,*) '                 zhmax  = ', zhmax
(  185)             ENDIF
(  186)             WRITE(numout,*) '           Value of coefficients for vertical mesh:'
(  187)             WRITE(numout,*) '                 zsur = ', zsur
(  188)             WRITE(numout,*) '                 za0  = ', za0
(  189)             WRITE(numout,*) '                 za1  = ', za1
(  190)             WRITE(numout,*) '                 zkth = ', zkth
(  191)             WRITE(numout,*) '                 zacr = ', zacr
(  192)          ENDIF
(  193)       ENDIF
(  194) 
(  195) 
(  196)       ! Reference z-coordinate (depth - scale factor at T- and W-points)
(  197)       ! ======================
(  198)       IF (  ppkth == 0. ) THEN            !  uniform vertical grid       
(  199) 
(  200)          za1 = zhmax/FLOAT(jpk-1) 
(  201)          DO jk = 1, jpk
(  202)             zw = FLOAT( jk )
(  203)             zt = FLOAT( jk ) + 0.5






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 18

(  204)             gdepw(jk) = ( zw - 1 ) * za1
(  205)             gdept(jk) = ( zt - 1 ) * za1
(  206)             e3w  (jk) =  za1
(  207)             e3t  (jk) =  za1
(  208)          END DO
(  209) 
(  210)       ELSE
(  211) 
(  212)          DO jk = 1, jpk
(  213)             zw = FLOAT( jk )
(  214)             zt = FLOAT( jk ) + 0.5
(  215)             gdepw(jk) = ( zsur + za0 * zw + za1 * zacr * LOG( COSH( (zw-zkth)/zacr ) )  )
(  216)             gdept(jk) = ( zsur + za0 * zt + za1 * zacr * LOG( COSH( (zt-zkth)/zacr ) )  )
(  217)             e3w  (jk) =          za0      + za1        * TANH(      (zw-zkth)/zacr   )
(  218)             e3t  (jk) =          za0      + za1        * TANH(      (zt-zkth)/zacr   )
(  219)          END DO
(  220)          gdepw(1) = 0.e0   ! force first w-level to be exactly at zero
(  221) 
(  222) !!DBG: 2008.11.17: the params in par_SS008.h90 yield for the first 4 levels
(  223) !            level   gdept    gdepw     e3t      e3w  
(  224) !             1     3.05     0.00     6.19     6.00
(  225) !             2     9.45     6.19     6.64     6.40
(  226) !             3    16.36    12.84     7.20     6.90
(  227) !             4    23.90    20.04     7.89     7.53
(  228) !!Since the min #-levels = 2, this means that the min btm depth should be 12.84m
(  229) !!which we want to keep the same to start.
(  230) !!REM: e3t=thickness centered on T-pt; e3w=thickness centered on w-pt; 
(  231) !!e3w(1)=??? Pretend that e3w(1) ~ 2*gdept(1)
(  232) !!e.g. of change to a 3m thick top layer, while keeping gdepw(3) unchanged
(  233) !         gdepw(2) = 3.0 
(  234) !         gdept(1) = gdepw(2)/2.0
(  235) !         e3w(1) = 2.0*gdept(1)
(  236) !         e3t(1) = gdepw(2)-gdepw(1)
(  237) !         e3t(2) = gdepw(3)-gdepw(2) 
(  238) !         gdept(2) = gdepw(2) + e3t(2)/2.0 
(  239) !         e3w(2) = gdept(2)-gdept(1)
(  240) !         e3w(3) = gdept(3)-gdept(2)
(  241) 
(  242) 
(  243)       ENDIF
(  244) 
(  245)       ! Control and  print
(  246)       ! ==================
(  247) 
(  248)       IF(lwp) THEN
(  249)          WRITE(numout,*)
(  250)          WRITE(numout,*) '              Reference z-coordinate depth and scale factors:'
(  251)          WRITE(numout, "(9x,' level   gdept    gdepw     e3t      e3w  ')" )
(  252)          WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, gdept(jk), gdepw(jk), e3t(jk), e3w(jk), jk = 1, jpk )
(  253)       ENDIF
(  254) 
(  255)       DO jk = 1, jpk
(  256)          IF( e3w(jk) <= 0. .OR. e3t(jk) <= 0. ) THEN
(  257)             IF(lwp) WRITE(numout,cform_err)
(  258)             IF(lwp) WRITE(numout,*) ' e3w or e3t =< 0 '
(  259)             nstop = nstop + 1
(  260)          ENDIF
(  261)          IF( gdepw(jk) < 0. .OR. gdept(jk) < 0.) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 19

(  262)             IF(lwp) WRITE(numout,cform_err)
(  263)             IF(lwp) WRITE(numout,*) ' gdepw or gdept < 0 '
(  264)             nstop = nstop + 1
(  265)          ENDIF
(  266)       END DO
(  267) 
(  268)    END SUBROUTINE zgr_z
(  269) 
(  270) 
(  271)    SUBROUTINE zgr_bat
(  272)       !!----------------------------------------------------------------------
(  273)       !!                    ***  ROUTINE zgr_bat  ***
(  274)       !! 
(  275)       !! ** Purpose :   set bathymetry both in levels and meters
(  276)       !!
(  277)       !! ** Method  :   read or define mbathy and bathy arrays
(  278)       !!       * level bathymetry:
(  279)       !!      The ocean basin geometry is given by a two-dimensional array,
(  280)       !!      mbathy, which is defined as follow :
(  281)       !!            mbathy(ji,jj) = 1, ..., jpk-1, the number of ocean level
(  282)       !!                              at t-point (ji,jj).
(  283)       !!                            = 0  over the continental t-point.
(  284)       !!                            = -n over the nth island t-point.
(  285)       !!      The array mbathy is checked to verified its consistency with
(  286)       !!      model option. in particular:
(  287)       !!            mbathy must have at least 1 land grid-points (mbathy<=0)
(  288)       !!                  along closed boundary.
(  289)       !!            mbathy must be cyclic IF jperio=1.
(  290)       !!            mbathy must be lower or equal to jpk-1.
(  291)       !!            isolated ocean grid points are suppressed from mbathy
(  292)       !!                  since they are only connected to remaining
(  293)       !!                  ocean through vertical diffusion.
(  294)       !!      ntopo=-1 :   rectangular channel or bassin with a bump 
(  295)       !!      ntopo= 0 :   flat rectangular channel or basin 
(  296)       !!      ntopo= 1 :   mbathy is read in 'bathy_level.nc' NetCDF file
(  297)       !!                   bathy  is read in 'bathy_meter.nc' NetCDF file
(  298)       !!      C A U T I O N : mbathy will be modified during the initializa-
(  299)       !!      tion phase to become the number of non-zero w-levels of a water
(  300)       !!      column, with a minimum value of 1.
(  301)       !!
(  302)       !! ** Action  : - mbathy: level bathymetry (in level index)
(  303)       !!              - bathy : meter bathymetry (in meters)
(  304)       !!
(  305)       !! History :
(  306)       !!   9.0  !  03-08  (G. Madec)  Original code
(  307)       !!----------------------------------------------------------------------
(  308)       !! * Modules used
(  309)       USE ioipsl
(  310) 
(  311)       !! * Local declarations
(  312)       CHARACTER (len=18) ::   clname    ! temporary characters
(  313)       LOGICAL ::   llbon                ! check the existence of bathy files
(  314)       INTEGER ::   ji, jj, jl, jk       ! dummy loop indices
(  315)       INTEGER ::   inum = 11            ! temporary logical unit
(  316)       INTEGER  ::   &
(  317)          ipi, ipj, ipk,              &  ! temporary integers
(  318)          itime,                      &  !    "          "
(  319)          ii_bump, ij_bump               ! bump center position






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 20

(  320)       INTEGER, DIMENSION (1) ::   istep
(  321)       INTEGER , DIMENSION(jpidta,jpjdta) ::   &
(  322)          idta                           ! global domain integer data
(  323)       REAL(wp) ::   &
(  324)          r_bump, h_bump, h_oce,      &  ! bump characteristics 
(  325)          zi, zj, zdate0, zdt            ! temporary scalars
(  326)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  327)          zlamt, zphit,               &  ! temporary workspace (NetCDF read)
(  328)          zdta                           ! global domain scalar data
(  329)       REAL(wp), DIMENSION(jpk) ::   &
(  330)          zdept                          ! temporary workspace (NetCDF read)
(  331)       !!----------------------------------------------------------------------
(  332) 
(  333)       IF(lwp) WRITE(numout,*)
(  334)       IF(lwp) WRITE(numout,*) '    zgr_bat : defines level and meter bathymetry'
(  335)       IF(lwp) WRITE(numout,*) '    ~~~~~~~'
(  336) 
(  337)       ! ========================================
(  338)       ! global domain level and meter bathymetry (idta,zdta)
(  339)       ! ========================================
(  340)       !                                               ! =============== ! 
(  341)       IF( ntopo == 0 .OR. ntopo == -1 ) THEN          ! defined by hand !
(  342)          !                                            ! =============== !
(  343) 
(  344)          IF( ntopo == 0 ) THEN                        ! flat basin
(  345) 
(  346)             IF(lwp) WRITE(numout,*)
(  347)             IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin'
(  348) 
(  349)             idta(:,:) = jpkm1                            ! flat basin 
(  350)             zdta(:,:) = gdepw(jpk)
(  351) 
(  352)          ELSE                                         ! bump
(  353)             IF(lwp) WRITE(numout,*)
(  354)             IF(lwp) WRITE(numout,*) '         bathymetry field: flat basin with a bump'
(  355) 
(  356)             ii_bump = jpidta / 3 + 3       ! i-index of the bump center
(  357)             ij_bump = jpjdta / 2           ! j-index of the bump center
(  358)             r_bump  =    6              ! bump radius (index)       
(  359)             h_bump  =  240.e0           ! bump height (meters)
(  360)             h_oce   = gdepw(jpk)        ! background ocean depth (meters)
(  361)             IF(lwp) WRITE(numout,*) '            bump characteristics: '
(  362)             IF(lwp) WRITE(numout,*) '               bump center (i,j)   = ', ii_bump, ii_bump
(  363)             IF(lwp) WRITE(numout,*) '               bump height         = ', h_bump , ' meters'
(  364)             IF(lwp) WRITE(numout,*) '               bump radius         = ', r_bump , ' index'
(  365)             IF(lwp) WRITE(numout,*) '            background ocean depth = ', h_oce  , ' meters'
(  366)             ! zdta :
(  367)             DO jj = 1, jpjdta
(  368)                DO ji = 1, jpidta
(  369)                   zi = FLOAT( ji - ii_bump ) / r_bump      
(  370)                   zj = FLOAT( jj - ij_bump ) / r_bump       
(  371)                   zdta(ji,jj) = h_oce - h_bump * EXP( -( zi*zi + zj*zj ) )
(  372)                END DO
(  373)             END DO
(  374)             ! idta :
(  375)             idta(:,:) = jpkm1
(  376)             DO jk = 1, jpkm1
(  377)                DO jj = 1, jpjdta






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 21

(  378)                   DO ji = 1, jpidta
(  379)                      IF( gdept(jk) < zdta(ji,jj) .AND. zdta(ji,jj) <= gdept(jk+1) )   idta(ji,jj) = jk
(  380)                   END DO
(  381)                END DO
(  382)             END DO
(  383)          ENDIF
(  384) 
(  385)          ! set boundary conditions (caution, idta on the global domain: use of jperio, not nperio)
(  386)          IF( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) THEN
(  387)             idta( :    , 1    ) = -1                ;      zdta( :    , 1    ) = -1.e0
(  388)             idta( :    ,jpjdta) =  0                ;      zdta( :    ,jpjdta) =  0.e0
(  389)          ELSEIF( jperio == 2 ) THEN
(  390)             idta( :    , 1    ) = idta( : ,  3  )   ;      zdta( :    , 1    ) = zdta( : ,  3  )
(  391)             idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
(  392)             idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
(  393)             idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
(  394)          ELSE
(  395)             idta( :    , 1    ) = 0                 ;      zdta( :    , 1    ) =  0.e0
(  396)             idta( :    ,jpjdta) = 0                 ;      zdta( :    ,jpjdta) =  0.e0
(  397)             idta( 1    , :    ) = 0                 ;      zdta( 1    , :    ) =  0.e0
(  398)             idta(jpidta, :    ) = 0                 ;      zdta(jpidta, :    ) =  0.e0
(  399)          ENDIF
(  400) 
(  401)          !  EEL R5 configuration with east and west open boundaries.
(  402)          !  Two rows of zeroes are needed at the south and north for OBCs
(  403)          !  This is for compatibility with the rigid lid option. 
(  404)           
(  405) !         IF( cp_cfg == "eel" .AND. jp_cfg == 5 ) THEN
(  406) !            idta( : , 2      ) = 0                 ;      zdta( : , 2      ) =  0.e0
(  407) !            idta( : ,jpjdta-1) = 0                 ;      zdta( : ,jpjdta-1) =  0.e0
(  408) !         ENDIF
(  409) 
(  410)          !                                            ! =============== !
(  411)       ELSEIF( ntopo == 1 ) THEN                       !   read in file  !
(  412)          !                                            ! =============== !
(  413)          IF( lk_zco ) THEN
(  414)             clname = 'bathy_level.nc'                       ! Level bathymetry
(  415) # 420
(  420)             INQUIRE( FILE=clname, EXIST=llbon )
(  421)             IF( llbon ) THEN
(  422)                IF(lwp) WRITE(numout,*)
(  423)                IF(lwp) WRITE(numout,*) '         read level bathymetry in ', clname
(  424)                IF(lwp) WRITE(numout,*)
(  425)                itime = 1
(  426)                ipi = jpidta
(  427)                ipj = jpjdta
(  428)                ipk = 1
(  429)                zdt = rdt
(  430)                CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &
(  431)                               ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
(  432)                CALL flinget( inum, 'Bathy_level', jpidta, jpjdta, 1,   &
(  433)                              itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) )
(  434)                idta(:,:) = zdta(:,:)
(  435)                CALL flinclo( inum )
(  436) 
(  437)             ELSE
(  438)                IF(lwp) WRITE(numout,cform_err)
(  439)                IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 22

(  440)                nstop = nstop + 1
(  441)             ENDIF   
(  442)   
(  443)          ELSEIF( lk_zps ) THEN
(  444)             clname = 'bathy_meter.nc'                       ! meter bathymetry
(  445) # 450
(  450)             INQUIRE( FILE=clname, EXIST=llbon )
(  451)             IF( llbon ) THEN
(  452)                IF(lwp) WRITE(numout,*)
(  453)                IF(lwp) WRITE(numout,*) '         read meter bathymetry in ', clname
(  454)                IF(lwp) WRITE(numout,*)
(  455)                itime = 1
(  456)                ipi = jpidta
(  457)                ipj = jpjdta
(  458)                ipk = 1
(  459)                zdt = rdt
(  460)                CALL flinopen( clname, 1, jpidta, 1, jpjdta, .FALSE.,   &    
(  461)                               ipi, ipj, ipk, zlamt, zphit, zdept, itime, istep, zdate0, zdt, inum )
(  462)                CALL flinget( inum, 'Bathymetry', jpidta, jpjdta, 1,   &
(  463)                              itime, 1, 1, 1, jpidta, 1, jpjdta, zdta(:,:) ) 
(  464)                CALL flinclo( inum )
(  465)                idta(:,:) = jpkm1      ! initialisation
(  466)             ELSE
(  467)                IF(lwp) WRITE(numout,cform_err)       
(  468)                IF(lwp) WRITE(numout,*)'    zgr_bat : unable to read the file', clname
(  469)                nstop = nstop + 1
(  470)             ENDIF
(  471)          ENDIF
(  472)          !                                            ! =============== !
(  473)       ELSE                                            !      error      !
(  474)          !                                            ! =============== !
(  475)          IF(lwp) WRITE(numout,cform_err)
(  476)          IF(lwp) WRITE(numout,*) '          parameter , ntopo = ', ntopo
(  477)          nstop = nstop + 1
(  478)       ENDIF
(  479) 
(  480) 
(  481)       ! =======================================
(  482)       ! local domain level and meter bathymetry (mbathy,bathy)
(  483)       ! =======================================
(  484) 
(  485)       mbathy(:,:) = 0                                 ! set to zero extra halo points
(  486)       bathy (:,:) = 0.e0                              ! (require for mpp case)
(  487) 
(  488)       DO jj = 1, nlcj                                 ! interior values
(  489)          DO ji = 1, nlci
(  490)             mbathy(ji,jj) = idta( mig(ji), mjg(jj) )
(  491)             bathy (ji,jj) = zdta( mig(ji), mjg(jj) )
(  492)          END DO
(  493)       END DO
(  494) 
(  495)       ! =======================
(  496)       ! NO closed seas or lakes
(  497)       ! =======================
(  498) 
(  499)       IF( nclosea == 0 ) THEN
(  500)          DO jl = 1, jpncs
(  501)             DO jj = ncsj1(jl), ncsj2(jl)






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 23

(  502)                DO ji = ncsi1(jl), ncsi2(jl)
(  503)                   mbathy(ji,jj) = 0                   ! suppress closed seas
(  504)                   bathy (ji,jj) = 0.e0                ! and lakes
(  505)                END DO
(  506)             END DO
(  507)          END DO
(  508)       ENDIF
(  509) 
(  510)       ! ===========
(  511)       ! Zoom domain 
(  512)       ! ===========
(  513) 
(  514)       IF( lzoom )   CALL zgr_bat_zoom
(  515) 
(  516)       ! ================
(  517)       ! Bathymetry check
(  518)       ! ================
(  519) 
(  520)       CALL zgr_bat_ctl
(  521) 
(  522)    END SUBROUTINE zgr_bat
(  523) 
(  524) 
(  525)    SUBROUTINE zgr_bat_zoom
(  526)       !!----------------------------------------------------------------------
(  527)       !!                    ***  ROUTINE zgr_bat_zoom  ***
(  528)       !!
(  529)       !! ** Purpose : - Close zoom domain boundary if necessary
(  530)       !!              - Suppress Med Sea from ORCA R2 and R05 arctic zoom
(  531)       !!
(  532)       !! ** Method  : 
(  533)       !!
(  534)       !! ** Action  : - update mbathy: level bathymetry (in level index)
(  535)       !!
(  536)       !! History :
(  537)       !!   9.0  !  03-08  (G. Madec)  Original code
(  538)       !!----------------------------------------------------------------------
(  539)       !! * Local variables
(  540)       INTEGER ::   ii0, ii1, ij0, ij1   ! temporary integers
(  541)       !!----------------------------------------------------------------------
(  542) 
(  543)       IF(lwp) WRITE(numout,*)
(  544)       IF(lwp) WRITE(numout,*) '    zgr_bat_zoom : modify the level bathymetry for zoom domain'
(  545)       IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~'
(  546) 
(  547)       ! Zoom domain
(  548)       ! ===========
(  549) 
(  550)       ! Forced closed boundary if required
(  551)       IF( lzoom_w )   mbathy( mi0(jpizoom):mi1(jpizoom) , :  ) = 0
(  552)       IF( lzoom_s )   mbathy(  : , mj0(jpjzoom):mj1(jpjzoom) ) = 0
(  553)       IF( lzoom_e )   mbathy( mi0(jpiglo+jpizoom-1):mi1(jpiglo+jpizoom-1) , :  ) = 0
(  554)       IF( lzoom_n )   mbathy(  : , mj0(jpjglo+jpjzoom-1):mj1(jpjglo+jpjzoom-1) ) = 0
(  555) 
(  556) !!DB: delete ORCA
(  557)       ! Configuration specific domain modifications
(  558)       ! (here, ORCA arctic configuration: suppress Med Sea)
(  559) !      IF( cp_cfg == "orca" .AND. lzoom_arct ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 24

(  560) 
(  561)    END SUBROUTINE zgr_bat_zoom
(  562) 
(  563) 
(  564)    SUBROUTINE zgr_bat_ctl
(  565)       !!----------------------------------------------------------------------
(  566)       !!                    ***  ROUTINE zgr_bat_ctl  ***
(  567)       !!
(  568)       !! ** Purpose :   check the bathymetry in levels
(  569)       !!
(  570)       !! ** Method  :   The array mbathy is checked to verified its consistency
(  571)       !!      with the model options. in particular:
(  572)       !!            mbathy must have at least 1 land grid-points (mbathy<=0)
(  573)       !!                  along closed boundary.
(  574)       !!            mbathy must be cyclic IF jperio=1.
(  575)       !!            mbathy must be lower or equal to jpk-1.
(  576)       !!            isolated ocean grid points are suppressed from mbathy
(  577)       !!                  since they are only connected to remaining
(  578)       !!                  ocean through vertical diffusion.
(  579)       !!      C A U T I O N : mbathy will be modified during the initializa-
(  580)       !!      tion phase to become the number of non-zero w-levels of a water
(  581)       !!      column, with a minimum value of 1.
(  582)       !!
(  583)       !! ** Action  : - update mbathy: level bathymetry (in level index)
(  584)       !!              - update bathy : meter bathymetry (in meters)
(  585)       !!
(  586)       !! History :
(  587)       !!   9.0  !  03-08  (G. Madec)  Original code
(  588)       !!----------------------------------------------------------------------
(  589)       !! * Local declarations
(  590)       INTEGER ::   ji, jj, jl           ! dummy loop indices
(  591)       INTEGER ::   &
(  592)          icompt, ibtest, ikmax          ! temporary integers
(  593)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  594)          zbathy                         ! temporary workspace
(  595)       !!----------------------------------------------------------------------
(  596) 
(  597)       IF(lwp) WRITE(numout,*)
(  598)       IF(lwp) WRITE(numout,*) '    zgr_bat_ctl : check the bathymetry'
(  599)       IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
(  600) 
(  601)       ! ================
(  602)       ! Bathymetry check
(  603)       ! ================
(  604) 
(  605)       IF( .NOT. lk_cfg_1d )   THEN
(  606) 
(  607)          ! Suppress isolated ocean grid points
(  608) 
(  609)          IF(lwp) WRITE(numout,*)
(  610)          IF(lwp) WRITE(numout,*)'                   suppress isolated ocean grid points'
(  611)          IF(lwp) WRITE(numout,*)'                   -----------------------------------'
(  612) 
(  613)          icompt = 0
(  614) 
(  615)          DO jl = 1, 2
(  616) 
(  617)             IF( nperio == 1 .OR. nperio  ==  4 .OR. nperio  ==  6 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 25

(  618)                mbathy( 1 ,:) = mbathy(jpim1,:)
(  619)                mbathy(jpi,:) = mbathy(  2  ,:)
(  620)             ENDIF
(  621)             DO jj = 2, jpjm1
(  622)                DO ji = 2, jpim1
(  623)                   ibtest = MAX( mbathy(ji-1,jj), mbathy(ji+1,jj),   &
(  624)                      mbathy(ji,jj-1),mbathy(ji,jj+1) )
(  625)                   IF( ibtest < mbathy(ji,jj) ) THEN
(  626)                      IF(lwp) WRITE(numout,*) ' the number of ocean level at ',   &
(  627)                         'grid-point (i,j) =  ',ji,jj,' is changed from ',   &
(  628)                         mbathy(ji,jj),' to ', ibtest
(  629)                      mbathy(ji,jj) = ibtest
(  630)                      icompt = icompt + 1
(  631)                   ENDIF
(  632)                END DO
(  633)             END DO
(  634) 
(  635)          END DO
(  636)          IF( icompt == 0 ) THEN
(  637)             IF(lwp) WRITE(numout,*)'     no isolated ocean grid points'
(  638)          ELSE
(  639)             IF(lwp) WRITE(numout,*)'    ',icompt,' ocean grid points suppressed'
(  640)          ENDIF
(  641)          IF( lk_mpp ) THEN
(  642)             zbathy(:,:) = FLOAT( mbathy(:,:) )
(  643)             CALL lbc_lnk( zbathy, 'T', 1. )
(  644)             mbathy(:,:) = INT( zbathy(:,:) )
(  645)          ENDIF
(  646) 
(  647)          ! 3.2 East-west cyclic boundary conditions
(  648) 
(  649)          IF( nperio == 0 ) THEN
(  650)             IF(lwp) WRITE(numout,*) ' mbathy set to 0 along east and west',   &
(  651)                ' boundary: nperio = ', nperio
(  652)             IF( lk_mpp ) THEN
(  653)                IF( nbondi == -1 .OR. nbondi == 2 ) THEN
(  654)                   IF( jperio /= 1 )   mbathy(1,:) = 0
(  655)                ENDIF
(  656)                IF( nbondi == 1 .OR. nbondi == 2 ) THEN
(  657)                   IF( jperio /= 1 )   mbathy(nlci,:) = 0
(  658)                ENDIF
(  659)             ELSE
(  660)                mbathy( 1 ,:) = 0
(  661)                mbathy(jpi,:) = 0
(  662)             ENDIF
(  663)          ELSEIF( nperio == 1 .OR. nperio == 4 .OR. nperio ==  6 ) THEN
(  664)             IF(lwp) WRITE(numout,*)' east-west cyclic boundary conditions',   &
(  665)                ' on mbathy: nperio = ', nperio
(  666)             mbathy( 1 ,:) = mbathy(jpim1,:)
(  667)             mbathy(jpi,:) = mbathy(  2  ,:)
(  668)          ELSEIF( nperio == 2 ) THEN
(  669)             IF(lwp) WRITE(numout,*) '   equatorial boundary conditions',   &
(  670)                ' on mbathy: nperio = ', nperio
(  671)          ELSE
(  672)             IF(lwp) WRITE(numout,*) '    e r r o r'
(  673)             IF(lwp) WRITE(numout,*) '    parameter , nperio = ', nperio
(  674)             !         STOP 'dom_mba'
(  675)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 26

(  676) 
(  677)          ! Set to zero mbathy over islands if necessary  (lk_isl=F)
(  678)          IF( .NOT. lk_isl ) THEN    ! No island
(  679)             IF(lwp) WRITE(numout,*)
(  680)             IF(lwp) WRITE(numout,*) '         mbathy set to 0 over islands'
(  681)             IF(lwp) WRITE(numout,*) '         ----------------------------'
(  682) 
(  683)             mbathy(:,:) = MAX( 0, mbathy(:,:) )
(  684) 
(  685)             !  Boundary condition on mbathy
(  686)             IF( .NOT.lk_mpp ) THEN 
(  687)                
(  688)                !!bug ???  y reflechir!
(  689)                !   ... mono- or macro-tasking: T-point, >0, 2D array, no slab
(  690)                
(  691)                zbathy(:,:) = FLOAT( mbathy(:,:) )
(  692)                CALL lbc_lnk( zbathy, 'T', 1. )
(  693)                mbathy(:,:) = INT( zbathy(:,:) )
(  694)             ENDIF
(  695) 
(  696)          ENDIF
(  697) 
(  698)       ENDIF
(  699) 
(  700)       ! Number of ocean level inferior or equal to jpkm1
(  701) 
(  702)       ikmax = 0
(  703)       DO jj = 1, jpj
(  704)          DO ji = 1, jpi
(  705)             ikmax = MAX( ikmax, mbathy(ji,jj) )
(  706)          END DO
(  707)       END DO
(  708)       !!! test a faire:   ikmax = MAX( mbathy(:,:) )   ???
(  709) 
(  710)       IF( ikmax > jpkm1 ) THEN
(  711)          IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' >  jpk-1'
(  712)          IF(lwp) WRITE(numout,*) ' change jpk to ',ikmax+1,' to use the exact ead bathymetry'
(  713)       ELSE IF( ikmax < jpkm1 ) THEN
(  714)          IF(lwp) WRITE(numout,*) ' maximum number of ocean level = ', ikmax,' < jpk-1' 
(  715)          IF(lwp) WRITE(numout,*) ' you can decrease jpk to ', ikmax+1
(  716)       ENDIF
(  717) 
(  718)       IF( lwp .AND. nprint == 1 ) THEN
(  719)          WRITE(numout,*)
(  720)          WRITE(numout,*) ' bathymetric field '
(  721)          WRITE(numout,*) ' ------------------'
(  722)          WRITE(numout,*) ' number of non-zero T-levels '
(  723)          CALL prihin( mbathy, jpi, jpj, 1, jpi,   &
(  724)                       1     , 1  , jpj, 1, 3  ,   &
(  725)                       numout )
(  726)          WRITE(numout,*)
(  727)       ENDIF
(  728) 
(  729)    END SUBROUTINE zgr_bat_ctl
(  730) 
(  731) 
(  732) # 1 "./domzgr_zps.h90"
(    1)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:34      page 27

(    2)    !!                      ***  domzgr_zps.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !!----------------------------------------------------------------------
(    5)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    6)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_zps.h90,v 1.5 2006/03/09 17:21:56 opalod Exp $ 
(    7)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(    8)    !!----------------------------------------------------------------------
(    9) 
(   10) # 436
(  436)    !!----------------------------------------------------------------------
(  437)    !!   Default option :                                      Empty routine
(  438)    !!----------------------------------------------------------------------
(  439)    SUBROUTINE zgr_zps
(  440)    END SUBROUTINE zgr_zps
(  441) # 733 "domzgr.F90"
(  733) 
(  734) 
(  735) # 1 "./domzgr_s.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                     ***  domzgr.s.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !!----------------------------------------------------------------------
(    5)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    6)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_s.h90,v 1.3 2005/08/30 13:26:10 opalod Exp $ 
(    7)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(    8)    !!----------------------------------------------------------------------
(    9) 
(   10) # 259
(  259)    !!----------------------------------------------------------------------
(  260)    !!   Default option :                                      Empty routine
(  261)    !!----------------------------------------------------------------------
(  262)    SUBROUTINE zgr_s
(  263)    END SUBROUTINE zgr_s
(  264) # 736 "domzgr.F90"
(  736) 
(  737) 
(  738) # 738
(  738)    !!======================================================================
(  739) END MODULE domzgr
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (domzgr.F90: 157)
