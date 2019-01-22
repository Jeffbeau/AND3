


PGF90 (Version     10.5)          09/23/2010  09:49:51      page 1

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

Filename: limrst.F90

(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE limrst
(    3)    !!======================================================================
(    4)    !!                     ***  MODULE  limrst  ***
(    5)    !! Ice restart :  write the ice restart file
(    6)    !!======================================================================
(    7) #if defined key_ice_lim
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   10)    !!----------------------------------------------------------------------
(   11)    !!   lim_rst_write   : write of the restart file 
(   12)    !!   lim_rst_read    : read  the restart file 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE in_out_manager
(   16)    USE ice
(   17)    USE ioipsl
(   18)    USE dom_oce
(   19)    USE ice_oce         ! ice variables
(   20)    USE daymod
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Accessibility
(   26)    PUBLIC lim_rst_write  ! routine called by lim_step.F90
(   27)    PUBLIC lim_rst_read   ! routine called by lim_init.F90
(   28) 
(   29)    !!----------------------------------------------------------------------
(   30)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrst.F90,v 1.8 2006/03/10 10:35:43 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    !!----------------------------------------------------------------------
(   38)    !!   Default option                                          NetCDF file
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41)    SUBROUTINE lim_rst_write( niter )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                    ***  lim_rst_write  ***
(   44)       !!
(   45)       !! ** purpose  :   output of sea-ice variable in a netcdf file
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       ! Arguments
(   49)       INTEGER  ::    niter        ! number of iteration
(   50) 
(   51)       !- dummy variables :
(   52)       LOGICAL :: &
(   53)          llbon
(   54)       INTEGER :: &
(   55)          ji, jj
(   56)       INTEGER :: &
(   57)          inumwrs, it0, itime
(   58)       REAL(wp), DIMENSION(1) :: &
(   59)          zdept
(   60)       REAL(wp), DIMENSION(2) :: &
(   61)          zinfo
(   62)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(   63)          zmoment
(   64)       REAL(wp) :: &
(   65)          zsec, zdate0, zdt
(   66) 
(   67)       CHARACTER(len=45)  ::  ccfile
(   68) 
(   69) !sujie add-----------------------------------
(   70)       CHARACTER (len=8) ::   cln
(   71)      ! Name of the new restart file
(   72)        WRITE(cln,'(i4.4,i2.2,i2.2)') nyear, nmonth,nday
(   73)        ccfile = 'restart_ice_'//cln//'_out.nc'
(   74) !      ccfile = 'restart_ice_out.nc'
(   75) !--------------------------------------------
(   76) 
(   77) #if defined key_agrif
(   78)       if ( .NOT. Agrif_Root() ) then
(   79)          ccfile= TRIM(Agrif_CFixed())//'_'//TRIM(ccfile)
(   80)       endif
(   81) #endif
(   82) 
(   83)       inumwrs  = 61
(   84)       INQUIRE ( FILE = ccfile, EXIST = llbon )
(   85)       IF( llbon ) THEN
(   86)          OPEN ( UNIT = inumwrs , FILE = ccfile, STATUS = 'old' )
(   87)          CLOSE( inumwrs , STATUS = 'delete' )
(   88)       ENDIF
(   89) 
(   90) 
(   91)       it0      = niter
(   92)       zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
(   93)       zinfo(2) = FLOAT( it0   )   ! iteration number
(   94) 
(   95)       zsec     = 0.e0
(   96)       itime    = 0
(   97)       zdept(1) = 0.e0
(   98)       zdt      = rdt_ice * nstock
(   99) 
(  100)       ! Write in inumwrs
(  101) 
(  102)       DO jj = 1, jpj              ! 3D array: 10 time faster than 35 restput
(  103)          DO ji = 1, jpi
(  104)             zmoment(ji,jj,1)  = sxice(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 3

(  105)             zmoment(ji,jj,2)  = syice(ji,jj)
(  106)             zmoment(ji,jj,3)  = sxxice(ji,jj)
(  107)             zmoment(ji,jj,4)  = syyice(ji,jj)
(  108)             zmoment(ji,jj,5)  = sxyice(ji,jj)
(  109)             zmoment(ji,jj,6)  = sxsn(ji,jj)
(  110)             zmoment(ji,jj,7)  = sysn(ji,jj)
(  111)             zmoment(ji,jj,8)  = sxxsn(ji,jj)
(  112)             zmoment(ji,jj,9)  = syysn(ji,jj)
(  113)             zmoment(ji,jj,10) = sxysn(ji,jj)
(  114)             zmoment(ji,jj,11) = sxa(ji,jj)
(  115)             zmoment(ji,jj,12) = sya(ji,jj)
(  116)             zmoment(ji,jj,13) = sxxa(ji,jj)
(  117)             zmoment(ji,jj,14) = syya(ji,jj)
(  118)             zmoment(ji,jj,15) = sxya(ji,jj)
(  119)             zmoment(ji,jj,16) = sxc0(ji,jj)
(  120)             zmoment(ji,jj,17) = syc0(ji,jj)
(  121)             zmoment(ji,jj,18) = sxxc0(ji,jj)
(  122)             zmoment(ji,jj,19) = syyc0(ji,jj)
(  123)             zmoment(ji,jj,20) = sxyc0(ji,jj)
(  124)             zmoment(ji,jj,21) = sxc1(ji,jj)
(  125)             zmoment(ji,jj,22) = syc1(ji,jj)
(  126)             zmoment(ji,jj,23) = sxxc1(ji,jj)
(  127)             zmoment(ji,jj,24) = syyc1(ji,jj)
(  128)             zmoment(ji,jj,25) = sxyc1(ji,jj)
(  129)             zmoment(ji,jj,26) = sxc2(ji,jj)
(  130)             zmoment(ji,jj,27) = syc2(ji,jj)
(  131)             zmoment(ji,jj,28) = sxxc2(ji,jj)
(  132)             zmoment(ji,jj,29) = syyc2(ji,jj)
(  133)             zmoment(ji,jj,30) = sxyc2(ji,jj)
(  134)             zmoment(ji,jj,31) = sxst(ji,jj)
(  135)             zmoment(ji,jj,32) = syst(ji,jj)
(  136)             zmoment(ji,jj,33) = sxxst(ji,jj)
(  137)             zmoment(ji,jj,34) = syyst(ji,jj)
(  138)             zmoment(ji,jj,35) = sxyst(ji,jj)
(  139)          END DO
(  140)       END DO
(  141) 
(  142)       CALL ymds2ju( nyear, nmonth, nday, zsec, zdate0 )
(  143)       CALL restini( 'NONE', jpi, jpj, glamt, gphit, 1 , zdept, ccfile, itime, zdate0, zdt, &
(  144)          &         inumwrs, domain_id=nidom )
(  145)       
(  146)       CALL restput( inumwrs, 'info'   ,   1,   1, 2 , 0, zinfo   )  ! restart informations
(  147)        
(  148)       CALL restput( inumwrs, 'hicif'  , jpi, jpj, 1 , 0, hicif   )  ! prognostic variables 
(  149)       CALL restput( inumwrs, 'hsnif'  , jpi, jpj, 1 , 0, hsnif   )
(  150)       CALL restput( inumwrs, 'frld'   , jpi, jpj, 1 , 0, frld    )
(  151)       CALL restput( inumwrs, 'sist'   , jpi, jpj, 1 , 0, sist    )
(  152) # if defined key_coupled
(  153)       CALL restput( inumwrs, 'albege' , jpi, jpj, 1 , 0, albege  )
(  154) # endif
(  155)       CALL restput( inumwrs, 'tbif'   , jpi, jpj, 3 , 0, tbif    )
(  156)       CALL restput( inumwrs, 'u_ice'  , jpi, jpj, 1 , 0, u_ice   )
(  157)       CALL restput( inumwrs, 'v_ice'  , jpi, jpj, 1 , 0, v_ice   )
(  158)       CALL restput( inumwrs, 'gtaux'  , jpi, jpj, 1 , 0, gtaux  )
(  159)       CALL restput( inumwrs, 'gtauy'  , jpi, jpj, 1 , 0, gtauy  )
(  160)       CALL restput( inumwrs, 'qstoif' , jpi, jpj, 1 , 0, qstoif  )
(  161)       CALL restput( inumwrs, 'fsbbq'  , jpi, jpj, 1 , 0, fsbbq   )
(  162)       CALL restput( inumwrs, 'moment' , jpi, jpj, 35, 0, zmoment )






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 4

(  163) 
(  164)       
(  165)       CALL restclo( inumwrs )
(  166) 
(  167)    END SUBROUTINE lim_rst_write
(  168) 
(  169) 
(  170)    SUBROUTINE lim_rst_read( niter )
(  171)       !-----------------------------------------------------------------------
(  172)       !  restart from a state defined in a binary file
(  173)       !-----------------------------------------------------------------------
(  174)       ! Arguments
(  175)       INTEGER  ::   niter        ! number of iteration
(  176) 
(  177)       !- dummy variables :
(  178)       CHARACTER(len=45)  ::  ccfile
(  179)       INTEGER :: &
(  180)         ji, jj
(  181)       INTEGER :: &
(  182)          inumrst, it0, it1, itime, ibvar, ifice
(  183)       LOGICAL :: &
(  184)          llog
(  185)       REAL(wp),DIMENSION(jpi,jpj) :: &
(  186)          zlamt, zphit
(  187)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(  188)          zmoment
(  189)       REAL(wp),DIMENSION(1) :: &
(  190)          zdept
(  191)       REAL(wp),DIMENSION(2) :: &
(  192)          zinfo
(  193)       REAL(wp) :: &
(  194)          zdate0, zdt
(  195)       CHARACTER ( len = 10 ) ::  &
(  196)          clvnames(60)       
(  197) 
(  198) 
(  199)        ccfile = 'restart_ice_in.nc'
(  200) #if defined key_agrif
(  201)       if ( .NOT. Agrif_Root() ) then
(  202)          ccfile= TRIM(Agrif_CFixed())//'_'//TRIM(ccfile)
(  203)       endif
(  204) #endif
(  205) 
(  206)       !Initialisations
(  207)       inumrst    = 71
(  208)       it0        = nit000
(  209)       itime      = 0
(  210)       llog       = .FALSE.
(  211)       zlamt(:,:) = 0.
(  212)       zphit(:,:) = 0.
(  213)       zdept(1)   = 0.
(  214) 
(  215)       CALL restini(ccfile , jpi, jpj, zlamt, zphit, 1 , zdept, 'NONE', itime, zdate0, zdt, inumrst, &
(  216)          &         domain_id=nidom )      
(  217)       CALL ioget_vname( inumrst, ibvar, clvnames )
(  218) 
(  219)       CALL restget    ( inumrst,'info', 1, 1 , 2, 0, llog, zinfo )
(  220)  






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 5

(  221)       ifice   = INT( zinfo(1) )
(  222)       it1     = INT( zinfo(2) )
(  223) 
(  224)       IF(lwp) THEN
(  225)          WRITE(numout,*)
(  226)          WRITE(numout,*) 'lim_rst_read : READ restart file name ', ccfile,' at time step : ', it1
(  227)          WRITE(numout,*) '~~~~~~~~~~~~   number of variables   : ', ibvar
(  228)          WRITE(numout,*) '               NetCDF variables      : ', clvnames(1:ibvar)
(  229)       ENDIF
(  230) 
(  231)       
(  232)       !Control of date
(  233)       
(  234)       IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
(  235)          IF(lwp) THEN
(  236)             WRITE(numout,cform_err)
(  237)             WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
(  238)             WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
(  239)             WRITE(numout,*) '   control of time parameter  nrstdt'
(  240)             nstop = nstop + 1
(  241)          ENDIF
(  242)       ENDIF
(  243) 
(  244)       CALL restget( inumrst, 'hicif'  , jpi, jpj, 1 , 0, llog, hicif   )
(  245)       CALL restget( inumrst, 'hsnif'  , jpi, jpj, 1 , 0, llog, hsnif   )
(  246)       CALL restget( inumrst, 'frld'   , jpi, jpj, 1 , 0, llog, frld    )
(  247)       CALL restget( inumrst, 'sist'   , jpi, jpj, 1 , 0, llog, sist    )
(  248) # if defined key_coupled 
(  249)       CALL restget( inumrst, 'albege' , jpi, jpj, 1 , 0, llog, albege  )
(  250) # endif
(  251)       CALL restget( inumrst, 'tbif'   , jpi, jpj, 3 , 0, llog, tbif    )
(  252)       CALL restget( inumrst, 'u_ice'  , jpi, jpj, 1 , 0, llog, u_ice   )
(  253)       CALL restget( inumrst, 'v_ice'  , jpi, jpj, 1 , 0, llog, v_ice   )
(  254)       CALL restget( inumrst, 'gtaux'  , jpi, jpj, 1 , 0, llog, gtaux  )
(  255)       CALL restget( inumrst, 'gtauy'  , jpi, jpj, 1 , 0, llog, gtauy  )
(  256)       CALL restget( inumrst, 'qstoif' , jpi, jpj, 1 , 0, llog, qstoif  )
(  257)       CALL restget( inumrst, 'fsbbq'  , jpi, jpj, 1 , 0, llog, fsbbq   )
(  258)       CALL restget( inumrst, 'moment' , jpi, jpj, 35, 0, llog, zmoment )
(  259) 
(  260)       CALL restclo( inumrst )
(  261) 
(  262)       niter = it1
(  263)       DO jj = 1, jpj
(  264)          DO ji = 1, jpi
(  265)             sxice(ji,jj)  = zmoment(ji,jj,1)
(  266)             syice(ji,jj)  = zmoment(ji,jj,2)
(  267)             sxxice(ji,jj) = zmoment(ji,jj,3)
(  268)             syyice(ji,jj) = zmoment(ji,jj,4)
(  269)             sxyice(ji,jj) = zmoment(ji,jj,5)
(  270)             sxsn(ji,jj)   = zmoment(ji,jj,6)
(  271)             sysn(ji,jj)   = zmoment(ji,jj,7)
(  272)             sxxsn(ji,jj)  = zmoment(ji,jj,8)
(  273)             syysn(ji,jj)  = zmoment(ji,jj,9)
(  274)             sxysn(ji,jj)  = zmoment(ji,jj,10)
(  275)             sxa(ji,jj)    = zmoment(ji,jj,11)
(  276)             sya(ji,jj)    = zmoment(ji,jj,12)
(  277)             sxxa(ji,jj)   = zmoment(ji,jj,13)
(  278)             syya(ji,jj)   = zmoment(ji,jj,14)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 6

(  279)             sxya(ji,jj)   = zmoment(ji,jj,15)
(  280)             sxc0(ji,jj)   = zmoment(ji,jj,16)
(  281)             syc0(ji,jj)   = zmoment(ji,jj,17)
(  282)             sxxc0(ji,jj)  = zmoment(ji,jj,18)
(  283)             syyc0(ji,jj)  = zmoment(ji,jj,19)
(  284)             sxyc0(ji,jj)  = zmoment(ji,jj,20)
(  285)             sxc1(ji,jj)   = zmoment(ji,jj,21)
(  286)             syc1(ji,jj)   = zmoment(ji,jj,22)
(  287)             sxxc1(ji,jj)  = zmoment(ji,jj,23)
(  288)             syyc1(ji,jj)  = zmoment(ji,jj,24)
(  289)             sxyc1(ji,jj)  = zmoment(ji,jj,25)
(  290)             sxc2(ji,jj)   = zmoment(ji,jj,26)
(  291)             syc2(ji,jj)   = zmoment(ji,jj,27)
(  292)             sxxc2(ji,jj)  = zmoment(ji,jj,28)
(  293)             syyc2(ji,jj)  = zmoment(ji,jj,29)
(  294)             sxyc2(ji,jj)  = zmoment(ji,jj,30)
(  295)             sxst(ji,jj)   = zmoment(ji,jj,31)
(  296)             syst(ji,jj)   = zmoment(ji,jj,32)
(  297)             sxxst(ji,jj)  = zmoment(ji,jj,33)
(  298)             syyst(ji,jj)  = zmoment(ji,jj,34)
(  299)             sxyst(ji,jj)  = zmoment(ji,jj,35)
(  300)          END DO
(  301)       END DO
(  302) 
(  303)       
(  304)    END SUBROUTINE lim_rst_read
(  305) 
(  306) 
(  307) #else
(  308)    !!----------------------------------------------------------------------
(  309)    !!   Default option :       Empty module            NO LIM sea-ice model
(  310)    !!----------------------------------------------------------------------
(  311) CONTAINS
(  312)    SUBROUTINE lim_rst_read             ! Empty routine
(  313)    END SUBROUTINE lim_rst_read
(  314)    SUBROUTINE lim_rst_write            ! Empty routine
(  315)    END SUBROUTINE lim_rst_write
(  316) #endif
(  317) 
(  318)    !!======================================================================
(  319) END MODULE limrst























PGF90 (Version     10.5)          09/23/2010  09:49:51      page 7

(    1) # 1 "limrst.F90"
(    1) !!DB -- 2009.09.04 -- key_diadimg eliminated
(    2) MODULE limrst
(    3)    !!======================================================================
(    4)    !!                     ***  MODULE  limrst  ***
(    5)    !! Ice restart :  write the ice restart file
(    6)    !!======================================================================
(    7) # 8
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   10)    !!----------------------------------------------------------------------
(   11)    !!   lim_rst_write   : write of the restart file 
(   12)    !!   lim_rst_read    : read  the restart file 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE in_out_manager
(   16)    USE ice
(   17)    USE ioipsl
(   18)    USE dom_oce
(   19)    USE ice_oce         ! ice variables
(   20)    USE daymod
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Accessibility
(   26)    PUBLIC lim_rst_write  ! routine called by lim_step.F90
(   27)    PUBLIC lim_rst_read   ! routine called by lim_init.F90
(   28) 
(   29)    !!----------------------------------------------------------------------
(   30)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrst.F90,v 1.8 2006/03/10 10:35:43 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    !!----------------------------------------------------------------------
(   38)    !!   Default option                                          NetCDF file
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41)    SUBROUTINE lim_rst_write( niter )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                    ***  lim_rst_write  ***
(   44)       !!
(   45)       !! ** purpose  :   output of sea-ice variable in a netcdf file
(   46)       !!
(   47)       !!----------------------------------------------------------------------
(   48)       ! Arguments
(   49)       INTEGER  ::    niter        ! number of iteration
(   50) 
(   51)       !- dummy variables :
(   52)       LOGICAL :: &
(   53)          llbon
(   54)       INTEGER :: &
(   55)          ji, jj
(   56)       INTEGER :: &
(   57)          inumwrs, it0, itime






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 8

(   58)       REAL(wp), DIMENSION(1) :: &
(   59)          zdept
(   60)       REAL(wp), DIMENSION(2) :: &
(   61)          zinfo
(   62)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(   63)          zmoment
(   64)       REAL(wp) :: &
(   65)          zsec, zdate0, zdt
(   66) 
(   67)       CHARACTER(len=45)  ::  ccfile
(   68) 
(   69) !sujie add-----------------------------------
(   70)       CHARACTER (len=8) ::   cln
(   71)      ! Name of the new restart file
(   72)        WRITE(cln,'(i4.4,i2.2,i2.2)') nyear, nmonth,nday
(   73)        ccfile = 'restart_ice_'//cln//'_out.nc'
(   74) !      ccfile = 'restart_ice_out.nc'
(   75) !--------------------------------------------
(   76) 
(   77) 
(   78) # 83
(   83)       inumwrs  = 61
(   84)       INQUIRE ( FILE = ccfile, EXIST = llbon )
(   85)       IF( llbon ) THEN
(   86)          OPEN ( UNIT = inumwrs , FILE = ccfile, STATUS = 'old' )
(   87)          CLOSE( inumwrs , STATUS = 'delete' )
(   88)       ENDIF
(   89) 
(   90) 
(   91)       it0      = niter
(   92)       zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
(   93)       zinfo(2) = FLOAT( it0   )   ! iteration number
(   94) 
(   95)       zsec     = 0.e0
(   96)       itime    = 0
(   97)       zdept(1) = 0.e0
(   98)       zdt      = rdt_ice * nstock
(   99) 
(  100)       ! Write in inumwrs
(  101) 
(  102)       DO jj = 1, jpj              ! 3D array: 10 time faster than 35 restput
(  103)          DO ji = 1, jpi
(  104)             zmoment(ji,jj,1)  = sxice(ji,jj)
(  105)             zmoment(ji,jj,2)  = syice(ji,jj)
(  106)             zmoment(ji,jj,3)  = sxxice(ji,jj)
(  107)             zmoment(ji,jj,4)  = syyice(ji,jj)
(  108)             zmoment(ji,jj,5)  = sxyice(ji,jj)
(  109)             zmoment(ji,jj,6)  = sxsn(ji,jj)
(  110)             zmoment(ji,jj,7)  = sysn(ji,jj)
(  111)             zmoment(ji,jj,8)  = sxxsn(ji,jj)
(  112)             zmoment(ji,jj,9)  = syysn(ji,jj)
(  113)             zmoment(ji,jj,10) = sxysn(ji,jj)
(  114)             zmoment(ji,jj,11) = sxa(ji,jj)
(  115)             zmoment(ji,jj,12) = sya(ji,jj)
(  116)             zmoment(ji,jj,13) = sxxa(ji,jj)
(  117)             zmoment(ji,jj,14) = syya(ji,jj)
(  118)             zmoment(ji,jj,15) = sxya(ji,jj)
(  119)             zmoment(ji,jj,16) = sxc0(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 9

(  120)             zmoment(ji,jj,17) = syc0(ji,jj)
(  121)             zmoment(ji,jj,18) = sxxc0(ji,jj)
(  122)             zmoment(ji,jj,19) = syyc0(ji,jj)
(  123)             zmoment(ji,jj,20) = sxyc0(ji,jj)
(  124)             zmoment(ji,jj,21) = sxc1(ji,jj)
(  125)             zmoment(ji,jj,22) = syc1(ji,jj)
(  126)             zmoment(ji,jj,23) = sxxc1(ji,jj)
(  127)             zmoment(ji,jj,24) = syyc1(ji,jj)
(  128)             zmoment(ji,jj,25) = sxyc1(ji,jj)
(  129)             zmoment(ji,jj,26) = sxc2(ji,jj)
(  130)             zmoment(ji,jj,27) = syc2(ji,jj)
(  131)             zmoment(ji,jj,28) = sxxc2(ji,jj)
(  132)             zmoment(ji,jj,29) = syyc2(ji,jj)
(  133)             zmoment(ji,jj,30) = sxyc2(ji,jj)
(  134)             zmoment(ji,jj,31) = sxst(ji,jj)
(  135)             zmoment(ji,jj,32) = syst(ji,jj)
(  136)             zmoment(ji,jj,33) = sxxst(ji,jj)
(  137)             zmoment(ji,jj,34) = syyst(ji,jj)
(  138)             zmoment(ji,jj,35) = sxyst(ji,jj)
(  139)          END DO
(  140)       END DO
(  141) 
(  142)       CALL ymds2ju( nyear, nmonth, nday, zsec, zdate0 )
(  143)       CALL restini( 'NONE', jpi, jpj, glamt, gphit, 1 , zdept, ccfile, itime, zdate0, zdt, &
(  144)          &         inumwrs, domain_id=nidom )
(  145)       
(  146)       CALL restput( inumwrs, 'info'   ,   1,   1, 2 , 0, zinfo   )  ! restart informations
(  147)        
(  148)       CALL restput( inumwrs, 'hicif'  , jpi, jpj, 1 , 0, hicif   )  ! prognostic variables 
(  149)       CALL restput( inumwrs, 'hsnif'  , jpi, jpj, 1 , 0, hsnif   )
(  150)       CALL restput( inumwrs, 'frld'   , jpi, jpj, 1 , 0, frld    )
(  151)       CALL restput( inumwrs, 'sist'   , jpi, jpj, 1 , 0, sist    )
(  152) # 155
(  155)       CALL restput( inumwrs, 'tbif'   , jpi, jpj, 3 , 0, tbif    )
(  156)       CALL restput( inumwrs, 'u_ice'  , jpi, jpj, 1 , 0, u_ice   )
(  157)       CALL restput( inumwrs, 'v_ice'  , jpi, jpj, 1 , 0, v_ice   )
(  158)       CALL restput( inumwrs, 'gtaux'  , jpi, jpj, 1 , 0, gtaux  )
(  159)       CALL restput( inumwrs, 'gtauy'  , jpi, jpj, 1 , 0, gtauy  )
(  160)       CALL restput( inumwrs, 'qstoif' , jpi, jpj, 1 , 0, qstoif  )
(  161)       CALL restput( inumwrs, 'fsbbq'  , jpi, jpj, 1 , 0, fsbbq   )
(  162)       CALL restput( inumwrs, 'moment' , jpi, jpj, 35, 0, zmoment )
(  163) 
(  164)       
(  165)       CALL restclo( inumwrs )
(  166) 
(  167)    END SUBROUTINE lim_rst_write
(  168) 
(  169) 
(  170)    SUBROUTINE lim_rst_read( niter )
(  171)       !-----------------------------------------------------------------------
(  172)       !  restart from a state defined in a binary file
(  173)       !-----------------------------------------------------------------------
(  174)       ! Arguments
(  175)       INTEGER  ::   niter        ! number of iteration
(  176) 
(  177)       !- dummy variables :
(  178)       CHARACTER(len=45)  ::  ccfile
(  179)       INTEGER :: &






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 10

(  180)         ji, jj
(  181)       INTEGER :: &
(  182)          inumrst, it0, it1, itime, ibvar, ifice
(  183)       LOGICAL :: &
(  184)          llog
(  185)       REAL(wp),DIMENSION(jpi,jpj) :: &
(  186)          zlamt, zphit
(  187)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(  188)          zmoment
(  189)       REAL(wp),DIMENSION(1) :: &
(  190)          zdept
(  191)       REAL(wp),DIMENSION(2) :: &
(  192)          zinfo
(  193)       REAL(wp) :: &
(  194)          zdate0, zdt
(  195)       CHARACTER ( len = 10 ) ::  &
(  196)          clvnames(60)       
(  197) 
(  198) 
(  199)        ccfile = 'restart_ice_in.nc'
(  200) 
(  201) # 206
(  206)       !Initialisations
(  207)       inumrst    = 71
(  208)       it0        = nit000
(  209)       itime      = 0
(  210)       llog       = .FALSE.
(  211)       zlamt(:,:) = 0.
(  212)       zphit(:,:) = 0.
(  213)       zdept(1)   = 0.
(  214) 
(  215)       CALL restini(ccfile , jpi, jpj, zlamt, zphit, 1 , zdept, 'NONE', itime, zdate0, zdt, inumrst, &
(  216)          &         domain_id=nidom )      
(  217)       CALL ioget_vname( inumrst, ibvar, clvnames )
(  218) 
(  219)       CALL restget    ( inumrst,'info', 1, 1 , 2, 0, llog, zinfo )
(  220)  
(  221)       ifice   = INT( zinfo(1) )
(  222)       it1     = INT( zinfo(2) )
(  223) 
(  224)       IF(lwp) THEN
(  225)          WRITE(numout,*)
(  226)          WRITE(numout,*) 'lim_rst_read : READ restart file name ', ccfile,' at time step : ', it1
(  227)          WRITE(numout,*) '~~~~~~~~~~~~   number of variables   : ', ibvar
(  228)          WRITE(numout,*) '               NetCDF variables      : ', clvnames(1:ibvar)
(  229)       ENDIF
(  230) 
(  231)       
(  232)       !Control of date
(  233)       
(  234)       IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
(  235)          IF(lwp) THEN
(  236)             WRITE(numout,cform_err)
(  237)             WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
(  238)             WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
(  239)             WRITE(numout,*) '   control of time parameter  nrstdt'
(  240)             nstop = nstop + 1
(  241)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 11

(  242)       ENDIF
(  243) 
(  244)       CALL restget( inumrst, 'hicif'  , jpi, jpj, 1 , 0, llog, hicif   )
(  245)       CALL restget( inumrst, 'hsnif'  , jpi, jpj, 1 , 0, llog, hsnif   )
(  246)       CALL restget( inumrst, 'frld'   , jpi, jpj, 1 , 0, llog, frld    )
(  247)       CALL restget( inumrst, 'sist'   , jpi, jpj, 1 , 0, llog, sist    )
(  248) # 251
(  251)       CALL restget( inumrst, 'tbif'   , jpi, jpj, 3 , 0, llog, tbif    )
(  252)       CALL restget( inumrst, 'u_ice'  , jpi, jpj, 1 , 0, llog, u_ice   )
(  253)       CALL restget( inumrst, 'v_ice'  , jpi, jpj, 1 , 0, llog, v_ice   )
(  254)       CALL restget( inumrst, 'gtaux'  , jpi, jpj, 1 , 0, llog, gtaux  )
(  255)       CALL restget( inumrst, 'gtauy'  , jpi, jpj, 1 , 0, llog, gtauy  )
(  256)       CALL restget( inumrst, 'qstoif' , jpi, jpj, 1 , 0, llog, qstoif  )
(  257)       CALL restget( inumrst, 'fsbbq'  , jpi, jpj, 1 , 0, llog, fsbbq   )
(  258)       CALL restget( inumrst, 'moment' , jpi, jpj, 35, 0, llog, zmoment )
(  259) 
(  260)       CALL restclo( inumrst )
(  261) 
(  262)       niter = it1
(  263)       DO jj = 1, jpj
(  264)          DO ji = 1, jpi
(  265)             sxice(ji,jj)  = zmoment(ji,jj,1)
(  266)             syice(ji,jj)  = zmoment(ji,jj,2)
(  267)             sxxice(ji,jj) = zmoment(ji,jj,3)
(  268)             syyice(ji,jj) = zmoment(ji,jj,4)
(  269)             sxyice(ji,jj) = zmoment(ji,jj,5)
(  270)             sxsn(ji,jj)   = zmoment(ji,jj,6)
(  271)             sysn(ji,jj)   = zmoment(ji,jj,7)
(  272)             sxxsn(ji,jj)  = zmoment(ji,jj,8)
(  273)             syysn(ji,jj)  = zmoment(ji,jj,9)
(  274)             sxysn(ji,jj)  = zmoment(ji,jj,10)
(  275)             sxa(ji,jj)    = zmoment(ji,jj,11)
(  276)             sya(ji,jj)    = zmoment(ji,jj,12)
(  277)             sxxa(ji,jj)   = zmoment(ji,jj,13)
(  278)             syya(ji,jj)   = zmoment(ji,jj,14)
(  279)             sxya(ji,jj)   = zmoment(ji,jj,15)
(  280)             sxc0(ji,jj)   = zmoment(ji,jj,16)
(  281)             syc0(ji,jj)   = zmoment(ji,jj,17)
(  282)             sxxc0(ji,jj)  = zmoment(ji,jj,18)
(  283)             syyc0(ji,jj)  = zmoment(ji,jj,19)
(  284)             sxyc0(ji,jj)  = zmoment(ji,jj,20)
(  285)             sxc1(ji,jj)   = zmoment(ji,jj,21)
(  286)             syc1(ji,jj)   = zmoment(ji,jj,22)
(  287)             sxxc1(ji,jj)  = zmoment(ji,jj,23)
(  288)             syyc1(ji,jj)  = zmoment(ji,jj,24)
(  289)             sxyc1(ji,jj)  = zmoment(ji,jj,25)
(  290)             sxc2(ji,jj)   = zmoment(ji,jj,26)
(  291)             syc2(ji,jj)   = zmoment(ji,jj,27)
(  292)             sxxc2(ji,jj)  = zmoment(ji,jj,28)
(  293)             syyc2(ji,jj)  = zmoment(ji,jj,29)
(  294)             sxyc2(ji,jj)  = zmoment(ji,jj,30)
(  295)             sxst(ji,jj)   = zmoment(ji,jj,31)
(  296)             syst(ji,jj)   = zmoment(ji,jj,32)
(  297)             sxxst(ji,jj)  = zmoment(ji,jj,33)
(  298)             syyst(ji,jj)  = zmoment(ji,jj,34)
(  299)             sxyst(ji,jj)  = zmoment(ji,jj,35)
(  300)          END DO
(  301)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 12

(  302) 
(  303)       
(  304)    END SUBROUTINE lim_rst_read
(  305) 
(  306) 
(  307) 
(  308) # 318
(  318)    !!======================================================================
(  319) END MODULE limrst
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (limrst.F90: 92)
