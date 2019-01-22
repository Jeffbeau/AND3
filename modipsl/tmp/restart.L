


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

Filename: restart.F90

(    1) !!DB - 2007.12.04
(    2) !!This version has eliminated the IOIPSL code
(    3) 
(    4) MODULE restart
(    5)    !!======================================================================
(    6)    !!                     ***  MODULE  restart  ***
(    7)    !! Ocean restart :  write the ocean restart file
(    8)    !!=====================================================================
(    9) 
(   10)    !!----------------------------------------------------------------------
(   11)    !!   rst_write  : write of the restart file
(   12)    !!   rst_read   : read the restart file
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE oce             ! ocean dynamics and tracers 
(   17)    USE phycst          ! physical constants
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE daymod          ! calendar
(   20)    USE sol_oce         ! ocean elliptic solver
(   21)    USE zdf_oce         ! ???
(   22)    USE zdftke          ! turbulent kinetic energy scheme
(   23)    USE ice_oce         ! ice variables
(   24)    USE blk_oce         ! bulk variables
(   25)    USE flx_oce         ! sea-ice/ocean forcings variables
(   26)    USE dynspg_oce      ! free surface time splitting scheme variables
(   27)    USE cpl_oce,         ONLY : lk_cpl              !
(   28) 
(   29)    USE lib_ncdf          ! netCDF I/O library
(   30) #if defined key_ice_lim
(   31)    USE ice               ! for ice related I/O
(   32) #endif
(   33)    IMPLICIT NONE
(   34)    PRIVATE
(   35) 
(   36)    !! * Routine accessibility
(   37)    PUBLIC rst_write  ! routine called by step.F90
(   38)    PUBLIC rst_read   ! routine called by inidtr.F90
(   39) #if defined key_ice_lim
(   40) !!DB for ice
(   41)    PUBLIC rst_ice_write  ! routine called by LIM_SRC/icestp.F90
(   42)    PUBLIC rst_ice_read   ! routine called by LIM_SRC/iceini
(   43) #endif 
(   44) 
(   45) 
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 2

(   47)    !! * Module variables
(   48)    CHARACTER (len=48) ::   &
(   49)       crestart = 'initial.nc'   ! restart file name
(   50)    INTEGER :: write_count
(   51) #if defined key_ice_lim
(   52)    INTEGER ::  ice_write_count = 0
(   53) #endif
(   54)    !!----------------------------------------------------------------------
(   55)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   56)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/restart.F90,v 1.15 2006/03/10 10:55:35 opalod Exp $ 
(   57)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   58)    !!----------------------------------------------------------------------
(   59) 
(   60) 
(   61) CONTAINS
(   62) 
(   63) 
(   64)    !!----------------------------------------------------------------------
(   65)    !!   Default option                                          NetCDF file
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68)    SUBROUTINE rst_write( kt )
(   69)       !!---------------------------------------------------------------------
(   70)       !!                   ***  ROUTINE rstwrite  ***
(   71)       !!                     
(   72)       !! ** Purpose :   Write restart fields in NetCDF format
(   73)       !!
(   74)       !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
(   75)       !!      file, save fields which are necessary for restart
(   76)       !!
(   77)       !! History :
(   78)       !!        !  99-11  (M. Imbard)  Original code
(   79)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   80)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   81)       !!----------------------------------------------------------------------
(   82)       !! * Modules used
(   83)       USE ioipsl
(   84) 
(   85)       !! * Arguments 
(   86)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   87) 
(   88)       !! * Local declarations
(   89)       LOGICAL ::   llbon
(   90)       CHARACTER (len=50) ::   clname, cln
(   91)       INTEGER ::   ic, jc, itime
(   92)       REAL(wp) ::   zdate0
(   93)       REAL(wp), DIMENSION( 1) ::   zfice, zfblk   ! used only in case of ice & bulk
(   94)       REAL(wp), DIMENSION(10) ::   zinfo(10)
(   95)       REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(   96)       
(   97)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
(   98)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
(   99)       INTEGER :: lncdf_stat                    ! lib_ncdf call return status
(  100)       CHARACTER(LEN=80) :: restfilename
(  101)       CHARACTER(LEN=4) :: chnum
(  102) 
(  103) #if defined key_agrif
(  104)        Integer :: knum






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 3

(  105) #endif
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108)       IF( kt == nit000 ) THEN
(  109)          IF(lwp) WRITE(numout,*)
(  110)          IF(lwp) WRITE(numout,*) 'rst_wri : write restart.output NetCDF file'
(  111)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  112)          zfice(1) = 1.e0   ;   zfblk(1) = 1.e0
(  113)       ENDIF
(  114)       IF (kt == nit000) THEN
(  115)          write_count = 0
(  116)       END IF
(  117) 
(  118) 
(  119)       IF( MOD( kt, nstock ) == 0 .OR. kt == nitend ) THEN
(  120)          
(  121)          ! 0. Initializations
(  122)          ! ------------------
(  123) 
(  124)          IF(lwp) WRITE(numout,*) ' '
(  125)          IF(lwp) WRITE(numout,*) 'rst_write : write the restart file in NetCDF format ',   &
(  126)                                               'at it= ',kt,' date= ',ndastp
(  127)          IF(lwp) WRITE(numout,*) '~~~~~~~~~'
(  128) 
(  129)          ! Job informations
(  130)          zinfo(1) = FLOAT( no        )   ! job number
(  131)          zinfo(2) = FLOAT( kt        )   ! time-step
(  132)          zinfo(3) = FLOAT( 2 - nsolv )   ! pcg solver
(  133)          zinfo(4) = FLOAT( nsolv - 1 )   ! sor solver
(  134)          IF( lk_zdftke ) THEN
(  135)             zinfo(5) = 1.e0              ! TKE 
(  136)          ELSE
(  137)             zinfo(5) = 0.e0              ! no TKE 
(  138)          ENDIF
(  139)          zinfo(6) = FLOAT( ndastp )      ! date
(  140)          zinfo(7) = adatrj               ! number of elapsed days since the begining of the run
(  141) 
(  142) 
(  143) !!DB - code deleted ...
(  144) !         IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL?
(  145)       ! LIB_NCDF CALLS
(  146) 
(  147) !DB 
(  148) !      IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf?
(  149) !         WRITE(UNIT=restfilename,FMT="('rstfile_',I4.4,'.nc')") write_count
(  150) !!DB
(  151)          write(chnum,'(I4.4)')write_count
(  152)          restfilename = trim(cexper)//'_rstfile_'//chnum//'.nc'
(  153) 
(  154)          CALL ncdf_create_restart(restfilename, lncdf_stat)
(  155)          CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   
(  156)          CALL ncdf_write(restfilename, 'ub', ub, nwrite, lncdf_stat) 
(  157)          CALL ncdf_write(restfilename, 'vb', vb, nwrite, lncdf_stat)  
(  158)          CALL ncdf_write(restfilename, 'tb', tb, nwrite, lncdf_stat) 
(  159)          CALL ncdf_write(restfilename, 'sb', sb, nwrite, lncdf_stat)
(  160)          CALL ncdf_write(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
(  161)          CALL ncdf_write(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
(  162)          CALL ncdf_write(restfilename, 'un', un, nwrite, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 4

(  163)          CALL ncdf_write(restfilename, 'vn', vn, nwrite, lncdf_stat)
(  164)          CALL ncdf_write(restfilename, 'tn', tn, nwrite, lncdf_stat) 
(  165)          CALL ncdf_write(restfilename, 'sn', sn, nwrite, lncdf_stat)
(  166)          CALL ncdf_write(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
(  167)          CALL ncdf_write(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
(  168)          ztab(:,:) = gcx(1:jpi,1:jpj)
(  169)          CALL ncdf_write(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
(  170)          ztab(:,:) = gcxb(1:jpi,1:jpj)
(  171)          CALL ncdf_write(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)       
(  172) #if defined key_dynspg_rl
(  173)          CALL ncdf_write(restfilename, 'bsfb', bsfb, nwrite, lncdf_stat)! Rigid-lid formulation (bsf)
(  174)          CALL ncdf_write(restfilename, 'bsfn', bsfn, nwrite, lncdf_stat)
(  175)          CALL ncdf_write(restfilename, 'bsfd', bsfd, nwrite, lncdf_stat)
(  176) #else
(  177)          CALL ncdf_write(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
(  178)          CALL ncdf_write(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
(  179) # if defined key_dynspg_ts      
(  180)          CALL ncdf_write(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
(  181)          CALL ncdf_write(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
(  182)          CALL ncdf_write(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
(  183)          CALL ncdf_write(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
(  184) # endif
(  185) #endif
(  186) #if defined key_zdftke   ||   defined key_esopa
(  187)          IF( lk_zdftke ) THEN
(  188)             CALL ncdf_write(restfilename, 'en', en, nwrite, lncdf_stat)   ! TKE arrays
(  189)          ENDIF
(  190) #endif
(  191) #if defined key_ice_lim
(  192)           zfice(1) = FLOAT( nfice )                                      ! Louvain La Neuve Sea Ice Model
(  193) !ADbug          CALL ncdf_write(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
(  194)          CALL ncdf_write(restfilename, 'nfice', zfice, lncdf_stat)
(  195)          CALL ncdf_write(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
(  196)          CALL ncdf_write(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
(  197)          CALL ncdf_write(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
(  198)          CALL ncdf_write(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
(  199) # if defined key_coupled
(  200)          CALL ncdf_write(restfilename, 'alb_ice', alb_ice, nwrite, lncdf_stat)
(  201) # endif
(  202) #endif
(  203) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  204)          zfblk(1) = FLOAT( nfbulk )                                 ! Bulk
(  205) ! ADbug         CALL ncdf_write(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
(  206)          CALL ncdf_write(restfilename, 'nfbulk', zfblk, lncdf_stat)
(  207)          CALL ncdf_write(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
(  208) #endif
(  209) 
(  210)          CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
(  211)          CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
(  212)          CALL ncdf_write(restfilename, 'nav_lev', gdept, lncdf_stat)
(  213)          CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  214)          CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  215)          write_count = write_count + 1
(  216) 
(  217) !      END IF ! End lib_ncdf block
(  218) 
(  219)       ! END LIB_NCDF CALLS
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 5

(  221)    ENDIF
(  222) 
(  223)    END SUBROUTINE rst_write
(  224) 
(  225) 
(  226) 
(  227)    SUBROUTINE rst_read
(  228)      !!---------------------------------------------------------------------- 
(  229)      !!                   ***  ROUTINE rst_read  ***
(  230)      !! 
(  231)      !! ** Purpose :   Read files for restart
(  232)      !! 
(  233)      !! ** Method  :   Read the previous fields on the NetCDF file 
(  234)      !!      the first record indicates previous characterics
(  235)      !!      after control with the present run, we read :
(  236)      !!      - prognostic variables on the second record
(  237)      !!      - elliptic solver arrays 
(  238)      !!      - barotropic stream function arrays ("key_dynspg_rl" defined)
(  239)      !!        or free surface arrays 
(  240)      !!      - tke arrays (lk_zdftke=T)
(  241)      !!      for this last three records,  the previous characteristics 
(  242)      !!      could be different with those used in the present run. 
(  243)      !!
(  244)      !!   According to namelist parameter nrstdt,
(  245)      !!       nrstdt = 0  no control on the date (nit000 is  arbitrary).
(  246)      !!       nrstdt = 1  we verify that nit000 is equal to the last
(  247)      !!                   time step of previous run + 1.
(  248)      !!       In both those options, the  exact duration of the experiment
(  249)      !!       since the beginning (cumulated duration of all previous restart runs)
(  250)      !!       is not stored in the restart and is assumed to be (nit000-1)*rdt.
(  251)      !!       This is valid is the time step has remained constant.
(  252)      !!
(  253)      !!       nrstdt = 2  the duration of the experiment in days (adatrj)
(  254)      !!                    has been stored in the restart file.
(  255)      !!
(  256)      !! History :
(  257)      !!        !  99-05  (M. Imbard)  Original code
(  258)      !!   8.5  !  02-09  (G. Madec)  F90: Free form
(  259)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  260)      !!----------------------------------------------------------------------
(  261)      !! * Modules used
(  262)      USE ioipsl
(  263)      
(  264)      !! * Local declarations
(  265)      LOGICAL ::   llog
(  266)      CHARACTER (len=8 ) ::   clvnames(50)
(  267)      CHARACTER (len=32) ::   clname
(  268)      INTEGER  ::   &
(  269)           itime, ibvar,     &  !
(  270)           inum                 ! temporary logical unit
(  271)      REAL(wp) ::   zdate0, zdt, zinfo(10)
(  272)      REAL(wp) ::   zdept(jpk), zlamt(jpi,jpj), zphit(jpi,jpj)
(  273)      REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(  274) #   if defined key_ice_lim
(  275)      INTEGER  ::   ios1, ji, jj, jn
(  276)      REAL(wp) ::   zfice(1)
(  277) #   endif
(  278) # if defined key_flx_bulk_monthly || defined key_flx_bulk_daily






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 6

(  279)      INTEGER  ::   ios2, jk
(  280)      REAL(wp) ::   zfblk(1)
(  281) #   endif
(  282)      
(  283)      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.   ! Use IOIPSL subroutines for restart input
(  284)      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.  ! Use lib_ncdf subroutines for restart input
(  285)      INTEGER :: lncdf_stat                     ! lib_ncdf routine return status
(  286)      REAL :: tstp
(  287)      CHARACTER(LEN=40) :: restfilename="restart.nc" !Default filename 
(  288)      
(  289)      
(  290)      !!----------------------------------------------------------------------
(  291)      !!  OPA 8.5, LODYC-IPSL (2002)
(  292)      !!----------------------------------------------------------------------
(  293)      clname = 'restart'
(  294) #if defined key_agrif       
(  295)      inum = Agrif_Get_Unit()
(  296)      If(.NOT. Agrif_root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  297) #endif 
(  298)      
(  299)      IF(lwp) WRITE(numout,*)
(  300)      IF(lwp) WRITE(numout,*) 'rst_read : read the NetCDF restart file'
(  301)      IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  302)      
(  303)      IF(lwp) WRITE(numout,*) ' Info on the present job : '
(  304)      IF(lwp) WRITE(numout,*) '   job number          : ', no
(  305)      IF(lwp) WRITE(numout,*) '   time-step           : ', nit000
(  306)      IF(lwp) WRITE(numout,*) '   solver type         : ', nsolv
(  307)      IF( lk_zdftke ) THEN
(  308)         IF(lwp) WRITE(numout,*) '   tke option          : 1 '
(  309)      ELSE
(  310)         IF(lwp) WRITE(numout,*) '   tke option          : 0 '
(  311)      ENDIF
(  312)      IF(lwp) WRITE(numout,*) '   date ndastp         : ', ndastp
(  313)      IF(lwp) WRITE(numout,*)
(  314)      
(  315)      ! Time domain : restart
(  316)      ! -------------------------
(  317)      
(  318)      IF(lwp) WRITE(numout,*)
(  319)      IF(lwp) WRITE(numout,*)
(  320)      IF(lwp) WRITE(numout,*) ' *** restart option'
(  321)      SELECT CASE ( nrstdt )
(  322)      CASE ( 0 ) 
(  323)         IF(lwp) WRITE(numout,*) ' nrstdt = 0 no control of nit000'
(  324)      CASE ( 1 ) 
(  325)         IF(lwp) WRITE(numout,*) ' nrstdt = 1 we control the date of nit000'
(  326)      CASE ( 2 )
(  327)         IF(lwp) WRITE(numout,*) ' nrstdt = 2 the date adatrj is read in restart file'
(  328)      CASE DEFAULT
(  329)         IF(lwp) WRITE(numout,*) '  ===>>>> nrstdt not equal 0, 1 or 2 : no control of the date'
(  330)         IF(lwp) WRITE(numout,*) ' =======                   ========='
(  331)      END SELECT
(  332)      
(  333)      itime = 0
(  334)      llog  = .FALSE.
(  335)      zlamt(:,:) = 0.e0
(  336)      zphit(:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 7

(  337)      zdept(:)   = 0.e0
(  338) 
(  339) !!DB -- code deleted ...
(  340) !     IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Do we want to use IOIPSL?
(  341) 
(  342)      ! LIB_NCDF CALLS
(  343) !!DB
(  344) !     IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Do we want to use lib_ncdf?
(  345) !        restfilename="restart.nc"
(  346)      CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
(  347) !!DB -- check if found file
(  348)      if(lwp .AND. lncdf_stat /= 0)then
(  349)         write(numout,*)'                '
(  350)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  351)         write(numout,*)'STOP: Problem reading restart file ', restfilename
(  352)         write(numout,*)'narea = ', narea
(  353)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  354)         write(numout,*)'                '
(  355)         stop
(  356)      endif
(  357) 
(  358)      
(  359)      ! NON LIB_NCDF RELATED CODE
(  360)      IF(lwp) WRITE(numout,*)
(  361)      IF(lwp) WRITE(numout,*) ' Info on the restart file read : '
(  362)      IF(lwp) WRITE(numout,*) '   job number          : ', NINT( zinfo(1) )
(  363)      IF(lwp) WRITE(numout,*) '   time-step           : ', NINT( zinfo(2) )
(  364)      IF(lwp) WRITE(numout,*) '   solver type         : ', NINT( zinfo(4) ) + 1
(  365)      IF(lwp) WRITE(numout,*) '   tke option          : ', NINT( zinfo(5) )
(  366)      IF(lwp) WRITE(numout,*) '   date ndastp         : ', NINT( zinfo(6) )
(  367)      IF(lwp) WRITE(numout,*)
(  368)      
(  369)      ! Control of date
(  370)      IF( nit000 - NINT( zinfo(2) )  /= 1 .AND. nrstdt /= 0 ) THEN
(  371)         IF(lwp) WRITE(numout,cform_err)
(  372)         IF(lwp) WRITE(numout,*) ' ===>>>> : problem with nit000 for the restart'
(  373)         IF(lwp) WRITE(numout,*) ' verify the restart file or rerun with nrstdt = 0 (namelist)'
(  374)         nstop = nstop + 1
(  375)      ENDIF
(  376)      
(  377)      ! re-initialisation of  adatrj0
(  378)      adatrj0 =  ( FLOAT( nit000-1 ) * rdttra(1) ) / rday
(  379)      
(  380)      IF ( nrstdt == 2 ) THEN
(  381)         !                             by default ndatsp has been set to ndate0 in dom_nam
(  382)         !                             ndate0 has been read in the namelist (standard OPA 8)
(  383)         !                             here when nrstdt=2 we keep the  final date of previous run
(  384)         ndastp = NINT( zinfo(6) )
(  385)         adatrj0 =  zinfo(7)
(  386)      ENDIF
(  387)      ! END NON LIB_NCDF RELATED CODE
(  388)      
(  389)      ! End of job info stuff
(  390)      CALL ncdf_read(restfilename, 'ub', ub, nwrite, lncdf_stat)   
(  391)      if(lncdf_stat /= 0) stop
(  392)      CALL ncdf_read(restfilename, 'vb', vb, nwrite, lncdf_stat)  
(  393)      if(lncdf_stat /= 0) stop
(  394)      CALL ncdf_read(restfilename, 'tb', tb, nwrite, lncdf_stat) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 8

(  395)      if(lncdf_stat /= 0) stop
(  396)      CALL ncdf_read(restfilename, 'sb', sb, nwrite, lncdf_stat)
(  397)      if(lncdf_stat /= 0) stop
(  398)      CALL ncdf_read(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
(  399)      if(lncdf_stat /= 0) stop
(  400)      CALL ncdf_read(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
(  401)      if(lncdf_stat /= 0) stop
(  402)      CALL ncdf_read(restfilename, 'un', un, nwrite, lncdf_stat)
(  403)      if(lncdf_stat /= 0) stop
(  404)      CALL ncdf_read(restfilename, 'vn', vn, nwrite, lncdf_stat)
(  405)      if(lncdf_stat /= 0) stop
(  406)      CALL ncdf_read(restfilename, 'tn', tn, nwrite, lncdf_stat) 
(  407)      if(lncdf_stat /= 0) stop
(  408)      CALL ncdf_read(restfilename, 'sn', sn, nwrite, lncdf_stat)
(  409)      if(lncdf_stat /= 0) stop
(  410)      CALL ncdf_read(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
(  411)      if(lncdf_stat /= 0) stop
(  412)      CALL ncdf_read(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
(  413)      if(lncdf_stat /= 0) stop
(  414)      CALL ncdf_read(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)    
(  415)      if(lncdf_stat /= 0) stop
(  416)      gcxb(1:jpi,1:jpj) = ztab(:,:) 
(  417)      CALL ncdf_read(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
(  418)      if(lncdf_stat /= 0) stop
(  419)      gcx(1:jpi,1:jpj) = ztab(:,:)   
(  420) #if defined key_dynspg_rl
(  421)      CALL ncdf_read(restfilename, 'bsfb', bsfb, nwrite, lncdf_stat)! Rigid-lid formulation (bsf)
(  422)      if(lncdf_stat /= 0) stop
(  423)      CALL ncdf_read(restfilename, 'bsfn', bsfn, nwrite, lncdf_stat)
(  424)      if(lncdf_stat /= 0) stop
(  425)      CALL ncdf_read(restfilename, 'bsfd', bsfd, nwrite, lncdf_stat)
(  426)      if(lncdf_stat /= 0) stop
(  427) #else
(  428)      CALL ncdf_read(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
(  429)      if(lncdf_stat /= 0) stop
(  430)      CALL ncdf_read(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
(  431)      if(lncdf_stat /= 0) stop
(  432) # if defined key_dynspg_ts      
(  433)      CALL ncdf_read(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
(  434)      if(lncdf_stat /= 0) stop
(  435)      CALL ncdf_read(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
(  436)      if(lncdf_stat /= 0) stop
(  437)      CALL ncdf_read(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
(  438)      if(lncdf_stat /= 0) stop
(  439)      CALL ncdf_read(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
(  440)      if(lncdf_stat /= 0) stop
(  441) # endif
(  442) #endif
(  443) #if defined key_zdftke   ||   defined key_esopa
(  444)      IF( lk_zdftke ) THEN
(  445)         IF( NINT( zinfo(5) ) == 1 ) THEN                                ! Read tke arrays
(  446)            CALL ncdf_read(restfilename, 'en', en, nwrite, lncdf_stat)
(  447)            ln_rstke = .FALSE.
(  448)         ELSE
(  449)            IF(lwp) WRITE(numout,*) ' ===>>>> : the previous restart file didnot used  tke scheme'
(  450)            IF(lwp) WRITE(numout,*) ' =======                ======='
(  451)            nrstdt = 2
(  452)            ln_rstke = .TRUE.






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 9

(  453)         ENDIF
(  454)      ENDIF
(  455) #endif
(  456) # if defined key_ice_lim
(  457)       ! Louvain La Neuve Sea Ice Model
(  458)       ios1 = 0
(  459)       DO jn = 1, 30
(  460)          IF( clvnames(jn) == 'nfice' )  ios1 = 1
(  461)       END DO
(  462)       IF( ios1 == 1 ) THEN
(  463) !ADbug          CALL ncdf_read(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
(  464)          CALL ncdf_read(restfilename, 'nfice', zfice, lncdf_stat)
(  465)          CALL ncdf_read(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
(  466)          CALL ncdf_read(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
(  467)          CALL ncdf_read(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
(  468)          CALL ncdf_read(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
(  469) # if defined key_coupled
(  470)          CALL ncdf_read(restfilename, 'alb_ice', alb_ice, nwrite, lncdf_stat)
(  471) # endif
(  472)       ENDIF
(  473)       IF( zfice(1) /= FLOAT(nfice) .OR. ios1 == 0 ) THEN
(  474)          IF(lwp) WRITE(numout,*)
(  475)          IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
(  476)          IF(lwp) WRITE(numout,*)
(  477)          sst_io(:,:) = ( nfice-1 )*( tn(:,:,1) + rt0 )          !!bug a explanation is needed here!
(  478)          sss_io(:,:) = ( nfice-1 )*  sn(:,:,1)
(  479)          DO jj = 2, jpj
(  480)             DO ji = 2, jpi
(  481)                u_io(ji,jj) = ( nfice-1 ) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  482)                v_io(ji,jj) = ( nfice-1 ) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  483)             END DO
(  484)          END DO
(  485) #    if defined key_coupled
(  486)          alb_ice(:,:) = 0.8 * tmask(:,:,1)
(  487) #    endif
(  488)       ENDIF
(  489) # endif
(  490) # if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  491)       ! Louvain La Neuve Sea Ice Model
(  492)       ios2 = 0
(  493)       DO jk = 1, 30
(  494)          IF( clvnames(jk) == 'nfbulk' )  ios2 = 1
(  495)       END DO
(  496)       IF( ios2 == 1 ) THEN
(  497) !ADbug          CALL ncdf_read(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
(  498)          CALL ncdf_read(restfilename, 'nfbulk', zfblk, lncdf_stat)
(  499)          CALL ncdf_read(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
(  500)       ENDIF
(  501)       IF( zfblk(1) /= FLOAT(nfbulk) .OR. ios2 == 0 ) THEN
(  502)          IF(lwp) WRITE(numout,*)
(  503)          IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
(  504)          IF(lwp) WRITE(numout,*)
(  505)          gsst(:,:) = 0.
(  506)          gsst(:,:) = gsst(:,:) + ( nfbulk-1 )*( tn(:,:,1) + rt0 )
(  507)       ENDIF
(  508) # endif
(  509) 
(  510)      CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 10

(  511)      if(lncdf_stat /= 0) stop
(  512)      CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
(  513)      if(lncdf_stat /= 0) stop
(  514)      CALL ncdf_read(restfilename, 'nav_lev', gdept, lncdf_stat)
(  515)      if(lncdf_stat /= 0) stop
(  516)      tstp = REAL(itime)
(  517)      CALL ncdf_read(restfilename, 'time', tstp, 1, lncdf_stat)
(  518)      if(lncdf_stat /= 0) stop
(  519)      itime = NINT(tstp)
(  520)      CALL ncdf_read(restfilename, 'time_steps', zdate0, 1, lncdf_stat)
(  521)      if(lncdf_stat /= 0) stop
(  522)   ! END LIB_NCDF CALLS
(  523)   
(  524)   ! In case of restart with neuler = 0 then put all before fields = to now fields
(  525)   IF ( neuler == 0 ) THEN
(  526)      tb(:,:,:)=tn(:,:,:)
(  527)      sb(:,:,:)=sn(:,:,:)
(  528)      ub(:,:,:)=un(:,:,:)
(  529)      vb(:,:,:)=vn(:,:,:)
(  530)      rotb(:,:,:)=rotn(:,:,:)
(  531)      hdivb(:,:,:)=hdivn(:,:,:)
(  532) #if defined key_dynspg_rl
(  533)      ! rigid lid
(  534)      bsfb(:,:)=bsfn(:,:)
(  535) #else
(  536)      ! free surface formulation (eta)
(  537)      sshb(:,:)=sshn(:,:)
(  538) #endif
(  539)   ENDIF
(  540)   
(  541) END SUBROUTINE rst_read
(  542) 
(  543) 
(  544) 
(  545) 
(  546) 
(  547) #ifdef key_ice_lim
(  548) !!DB 2008.05.22 
(  549) !!NB: called when needed by icestp
(  550) !!         IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
(  551) !            CALL lim_rst_write( numit )                              ! Ice restart file !
(  552) !         ENDIF
(  553) 
(  554)    SUBROUTINE rst_ice_write( niter )
(  555)       !!---------------------------------------------------------------------
(  556)       !!                   ***  ROUTINE rst ice write  ***
(  557)       !!                     
(  558)       !! ** Purpose :   Write lim ice restart fields in NetCDF format
(  559)       !!
(  560)       !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
(  561)       !!      file, save fields which are necessary for restart
(  562)       !!
(  563)       !! History : DB  2008.05.22
(  564)       !!
(  565)       !!----------------------------------------------------------------------
(  566)       !! * Modules used
(  567)       USE ioipsl
(  568) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 11

(  569)       !! * Arguments 
(  570)       INTEGER, INTENT( in ) ::   niter         ! ice iteration 
(  571) 
(  572)       !! * Local declarations
(  573)       LOGICAL ::   llbon
(  574)       CHARACTER (len=50) ::   clname, cln
(  575)       INTEGER ::   ic, jc, itime, it0
(  576)       REAL(wp), DIMENSION(2) ::   zinfo(2)
(  577)       REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(  578)       REAL(wp) ::  zsec, zdate0, zdt
(  579)       
(  580)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
(  581)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
(  582)       INTEGER :: lncdf_stat                    ! lib_ncdf call return status
(  583)       CHARACTER(LEN=80) :: restfilename
(  584)       CHARACTER(LEN=4) :: chnum
(  585) 
(  586)       REAL(wp),DIMENSION(jpi,jpj,35) :: zmoment
(  587)       INTEGER :: ji, jj
(  588) 
(  589) #if defined key_agrif
(  590)        Integer :: knum
(  591) #endif
(  592)       !!----------------------------------------------------------------------
(  593) 
(  594) 
(  595) ! Job informations
(  596)       it0      = niter
(  597)       zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
(  598)       zinfo(2) = FLOAT( it0   )   ! iteration number
(  599) 
(  600)       zsec     = 0.e0
(  601)       itime    = 0
(  602) !      zdept(1) = 0.e0   !! this is default z-level -- handled in ncdf_create_ice_restart()
(  603)       zdt      = rdt_ice * nstock  !!does not seem to enter restart file
(  604) 
(  605)       IF(lwp) WRITE(numout,*) ' '
(  606)       IF(lwp) WRITE(numout,*) 'ice_rst_write : write the ice restart file in NetCDF format ',   &
(  607)            'at it= ',niter,' date= ',ndastp
(  608)       IF(lwp) WRITE(numout,*) '~~~~~~~~~'
(  609) 
(  610) 
(  611) !!DB
(  612)       write(chnum,'(I4.4)')ice_write_count
(  613)       restfilename = trim(cexper)//'_ice_rstfile_'//chnum//'.nc'
(  614) 
(  615)       CALL ncdf_create_ice_restart(restfilename, lncdf_stat)
(  616)       CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   
(  617) 
(  618)       CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
(  619)       CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
(  620) !      CALL ncdf_write(restfilename, 'nav_lev', real(0), nwrite, lncdf_stat)   !!nav_lev treated as a scalar
(  621)       CALL ncdf_write(restfilename, 'nav_lev', REAL(0), 1, lncdf_stat)   !!nav_lev treated as a scalar
(  622)       CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  623)       CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  624) 
(  625)       CALL ncdf_write(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
(  626)       CALL ncdf_write(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 12

(  627)       CALL ncdf_write(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
(  628)       CALL ncdf_write(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
(  629) # if defined key_coupled
(  630)       CALL ncdf_write(restfilename, 'albege', albege   , nwrite, lncdf_stat)
(  631) # endif, nwrite, lncdf_stat)
(  632)       CALL ncdf_write(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
(  633)       CALL ncdf_write(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
(  634)       CALL ncdf_write(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
(  635)       CALL ncdf_write(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)
(  636)       CALL ncdf_write(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
(  637)       CALL ncdf_write(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)
(  638)       CALL ncdf_write(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)
(  639) 
(  640) 
(  641) !!DB load zmoment
(  642)       DO jj = 1, jpj  
(  643)          DO ji = 1, jpi
(  644)             zmoment(ji,jj,1)  = sxice(ji,jj)
(  645)             zmoment(ji,jj,2)  = syice(ji,jj)
(  646)             zmoment(ji,jj,3)  = sxxice(ji,jj)
(  647)             zmoment(ji,jj,4)  = syyice(ji,jj)
(  648)             zmoment(ji,jj,5)  = sxyice(ji,jj)
(  649)             zmoment(ji,jj,6)  = sxsn(ji,jj)
(  650)             zmoment(ji,jj,7)  = sysn(ji,jj)
(  651)             zmoment(ji,jj,8)  = sxxsn(ji,jj)
(  652)             zmoment(ji,jj,9)  = syysn(ji,jj)
(  653)             zmoment(ji,jj,10) = sxysn(ji,jj)
(  654)             zmoment(ji,jj,11) = sxa(ji,jj)
(  655)             zmoment(ji,jj,12) = sya(ji,jj)
(  656)             zmoment(ji,jj,13) = sxxa(ji,jj)
(  657)             zmoment(ji,jj,14) = syya(ji,jj)
(  658)             zmoment(ji,jj,15) = sxya(ji,jj)
(  659)             zmoment(ji,jj,16) = sxc0(ji,jj)
(  660)             zmoment(ji,jj,17) = syc0(ji,jj)
(  661)             zmoment(ji,jj,18) = sxxc0(ji,jj)
(  662)             zmoment(ji,jj,19) = syyc0(ji,jj)
(  663)             zmoment(ji,jj,20) = sxyc0(ji,jj)
(  664)             zmoment(ji,jj,21) = sxc1(ji,jj)
(  665)             zmoment(ji,jj,22) = syc1(ji,jj)
(  666)             zmoment(ji,jj,23) = sxxc1(ji,jj)
(  667)             zmoment(ji,jj,24) = syyc1(ji,jj)
(  668)             zmoment(ji,jj,25) = sxyc1(ji,jj)
(  669)             zmoment(ji,jj,26) = sxc2(ji,jj)
(  670)             zmoment(ji,jj,27) = syc2(ji,jj)
(  671)             zmoment(ji,jj,28) = sxxc2(ji,jj)
(  672)             zmoment(ji,jj,29) = syyc2(ji,jj)
(  673)             zmoment(ji,jj,30) = sxyc2(ji,jj)
(  674)             zmoment(ji,jj,31) = sxst(ji,jj)
(  675)             zmoment(ji,jj,32) = syst(ji,jj)
(  676)             zmoment(ji,jj,33) = sxxst(ji,jj)
(  677)             zmoment(ji,jj,34) = syyst(ji,jj)
(  678)             zmoment(ji,jj,35) = sxyst(ji,jj)
(  679)          END DO
(  680)       END DO
(  681) !!DB probably should call mppsync here ??? NB: Works if I do not
(  682) !      CALL MPI_BARRIER(MPI_COMM_WORLD, lncdf_stat)
(  683)       CALL ncdf_write(restfilename, 'moment', zmoment, nwrite, lncdf_stat)
(  684) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 13

(  685)       ice_write_count = ice_write_count + 1
(  686) 
(  687) 
(  688)    END SUBROUTINE rst_ice_write
(  689) 
(  690) 
(  691) !!DB: 2008.05.27 read ice restart 
(  692) !!Tested to run OK, but there is no proof that it did the right thing
(  693) !!NB; In general, there is no proof that ice restarts are perfect
(  694) !!Similar to old lim_rst_read()
(  695) !!NB: default restfilename = restart_ice.nc
(  696) !!NB: Agrif stuff not done yet
(  697)    SUBROUTINE rst_ice_read( niter )
(  698)       !-----------------------------------------------------------------------
(  699)       ! Arguments
(  700)       INTEGER  ::   niter        ! number of iteration
(  701) 
(  702)       !- dummy variables :
(  703)       CHARACTER(len=45)  ::  restfilename = 'restart_ice.nc'
(  704)       INTEGER :: &
(  705)         ji, jj, lncdf_stat
(  706)       INTEGER :: &
(  707)          inumrst, it0, it1, itime, ibvar, ifice
(  708)       LOGICAL :: &
(  709)          llog
(  710)       REAL(wp),DIMENSION(jpi,jpj) :: &
(  711)          zlamt, zphit
(  712)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(  713)          zmoment
(  714)       REAL(wp),DIMENSION(1) :: &
(  715)          zdept
(  716)       REAL(wp),DIMENSION(2) :: &
(  717)          zinfo
(  718)       REAL(wp) :: &
(  719)          zdate0, zdt
(  720)       CHARACTER ( len = 10 ) ::  &
(  721)          clvnames(60)       
(  722) 
(  723) 
(  724) 
(  725) #if defined key_agrif
(  726)       if ( .NOT. Agrif_Root() ) then
(  727)          restfilename= TRIM(Agrif_CFixed())//'_'//TRIM(restfilename)
(  728)       endif
(  729) #endif
(  730) 
(  731)      CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
(  732) !!DB -- check if found file
(  733)      if(lwp .AND. lncdf_stat /= 0)then
(  734)         write(numout,*)'                '
(  735)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  736)         write(numout,*)'STOP: Problem reading ice restart file ', restfilename
(  737)         write(numout,*)'narea = ', narea
(  738)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  739)         write(numout,*)'                '
(  740)         stop
(  741)      endif
(  742) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 14

(  743) 
(  744) 
(  745)       !Initialisations
(  746)       inumrst    = 71
(  747)       it0        = nit000
(  748)       itime      = 0
(  749)       llog       = .FALSE.
(  750)       zlamt(:,:) = 0.
(  751)       zphit(:,:) = 0.
(  752)       zdept(1)   = 0.
(  753) 
(  754)       ifice   = INT( zinfo(1) )
(  755)       it1     = INT( zinfo(2) )
(  756) 
(  757) !!DB Note that old code gets varnames but otherwise does not use this info 
(  758)       IF(lwp) THEN
(  759)          WRITE(numout,*)
(  760)          WRITE(numout,*) 'lim_rst_read : READ restart file name ', restfilename,' at time step : ', it1
(  761)       ENDIF
(  762)       
(  763)       !Control of date
(  764) !!DB: 2009.09.30
(  765)       IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
(  766)          IF(lwp) THEN
(  767)             WRITE(numout,cform_err)
(  768)             WRITE(numout,*) 'DBG: lim_rst_read nrstdt=1 and nit000 /= zinfo(2)'
(  769)             WRITE(numout,*) 'DBG: This case covered in restart() ====> Continue'
(  770) !            WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
(  771) !            WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
(  772) !            WRITE(numout,*) '   control of time parameter  nrstdt'
(  773) !            nstop = nstop + 1
(  774)          ENDIF
(  775)       ENDIF
(  776) 
(  777) !!DBG: writes from rst_ice_write()
(  778)      CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)
(  779)      if(lncdf_stat /= 0) stop
(  780)      CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
(  781)      if(lncdf_stat /= 0) stop
(  782) 
(  783) !!DBG ignore these vars to start
(  784) !      CALL ncdf_read(restfilename, 'nav_lev', zdept, 1, lncdf_stat)   !!nav_lev treated as a scalar
(  785) !      CALL ncdf_read(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  786) !      CALL ncdf_read(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  787) 
(  788)       CALL ncdf_read(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
(  789)       CALL ncdf_read(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)
(  790)       CALL ncdf_read(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
(  791)       CALL ncdf_read(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
(  792) # if defined key_coupled
(  793)       CALL ncdf_read(restfilename, 'albege', albege   , nwrite, lncdf_stat)
(  794) # endif, nwrite, lncdf_stat)
(  795)       CALL ncdf_read(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
(  796)       CALL ncdf_read(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
(  797)       CALL ncdf_read(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
(  798)       CALL ncdf_read(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)
(  799)       CALL ncdf_read(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
(  800)       CALL ncdf_read(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 15

(  801)       CALL ncdf_read(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)
(  802)       CALL ncdf_read(restfilename, 'moment', zmoment, nwrite, lncdf_stat)
(  803) 
(  804) !!DB: 2009.09.30
(  805) !      niter = it1
(  806)       niter = nit000
(  807) 
(  808)       DO jj = 1, jpj
(  809)          DO ji = 1, jpi
(  810)             sxice(ji,jj)  = zmoment(ji,jj,1)
(  811)             syice(ji,jj)  = zmoment(ji,jj,2)
(  812)             sxxice(ji,jj) = zmoment(ji,jj,3)
(  813)             syyice(ji,jj) = zmoment(ji,jj,4)
(  814)             sxyice(ji,jj) = zmoment(ji,jj,5)
(  815)             sxsn(ji,jj)   = zmoment(ji,jj,6)
(  816)             sysn(ji,jj)   = zmoment(ji,jj,7)
(  817)             sxxsn(ji,jj)  = zmoment(ji,jj,8)
(  818)             syysn(ji,jj)  = zmoment(ji,jj,9)
(  819)             sxysn(ji,jj)  = zmoment(ji,jj,10)
(  820)             sxa(ji,jj)    = zmoment(ji,jj,11)
(  821)             sya(ji,jj)    = zmoment(ji,jj,12)
(  822)             sxxa(ji,jj)   = zmoment(ji,jj,13)
(  823)             syya(ji,jj)   = zmoment(ji,jj,14)
(  824)             sxya(ji,jj)   = zmoment(ji,jj,15)
(  825)             sxc0(ji,jj)   = zmoment(ji,jj,16)
(  826)             syc0(ji,jj)   = zmoment(ji,jj,17)
(  827)             sxxc0(ji,jj)  = zmoment(ji,jj,18)
(  828)             syyc0(ji,jj)  = zmoment(ji,jj,19)
(  829)             sxyc0(ji,jj)  = zmoment(ji,jj,20)
(  830)             sxc1(ji,jj)   = zmoment(ji,jj,21)
(  831)             syc1(ji,jj)   = zmoment(ji,jj,22)
(  832)             sxxc1(ji,jj)  = zmoment(ji,jj,23)
(  833)             syyc1(ji,jj)  = zmoment(ji,jj,24)
(  834)             sxyc1(ji,jj)  = zmoment(ji,jj,25)
(  835)             sxc2(ji,jj)   = zmoment(ji,jj,26)
(  836)             syc2(ji,jj)   = zmoment(ji,jj,27)
(  837)             sxxc2(ji,jj)  = zmoment(ji,jj,28)
(  838)             syyc2(ji,jj)  = zmoment(ji,jj,29)
(  839)             sxyc2(ji,jj)  = zmoment(ji,jj,30)
(  840)             sxst(ji,jj)   = zmoment(ji,jj,31)
(  841)             syst(ji,jj)   = zmoment(ji,jj,32)
(  842)             sxxst(ji,jj)  = zmoment(ji,jj,33)
(  843)             syyst(ji,jj)  = zmoment(ji,jj,34)
(  844)             sxyst(ji,jj)  = zmoment(ji,jj,35)
(  845)          END DO
(  846)       END DO
(  847) 
(  848)       
(  849)    END SUBROUTINE rst_ice_read
(  850) 
(  851) 
(  852) 
(  853) #endif   !!key_ice_lim
(  854) 
(  855) 
(  856) !!=====================================================================
(  857) END MODULE restart
(  858) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 16

(    1) # 1 "restart.F90"
(    1) !!DB - 2007.12.04
(    2) !!This version has eliminated the IOIPSL code
(    3) 
(    4) MODULE restart
(    5)    !!======================================================================
(    6)    !!                     ***  MODULE  restart  ***
(    7)    !! Ocean restart :  write the ocean restart file
(    8)    !!=====================================================================
(    9) 
(   10)    !!----------------------------------------------------------------------
(   11)    !!   rst_write  : write of the restart file
(   12)    !!   rst_read   : read the restart file
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE oce             ! ocean dynamics and tracers 
(   17)    USE phycst          ! physical constants
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE daymod          ! calendar
(   20)    USE sol_oce         ! ocean elliptic solver
(   21)    USE zdf_oce         ! ???
(   22)    USE zdftke          ! turbulent kinetic energy scheme
(   23)    USE ice_oce         ! ice variables
(   24)    USE blk_oce         ! bulk variables
(   25)    USE flx_oce         ! sea-ice/ocean forcings variables
(   26)    USE dynspg_oce      ! free surface time splitting scheme variables
(   27)    USE cpl_oce,         ONLY : lk_cpl              !
(   28) 
(   29)    USE lib_ncdf          ! netCDF I/O library
(   30) # 31
(   31)    USE ice               ! for ice related I/O
(   32) # 33
(   33)    IMPLICIT NONE
(   34)    PRIVATE
(   35) 
(   36)    !! * Routine accessibility
(   37)    PUBLIC rst_write  ! routine called by step.F90
(   38)    PUBLIC rst_read   ! routine called by inidtr.F90
(   39) # 40
(   40) !!DB for ice
(   41)    PUBLIC rst_ice_write  ! routine called by LIM_SRC/icestp.F90
(   42)    PUBLIC rst_ice_read   ! routine called by LIM_SRC/iceini
(   43) 
(   44) 
(   45) 
(   46) # 47
(   47)    !! * Module variables
(   48)    CHARACTER (len=48) ::   &
(   49)       crestart = 'initial.nc'   ! restart file name
(   50)    INTEGER :: write_count
(   51) # 52
(   52)    INTEGER ::  ice_write_count = 0
(   53) # 54
(   54)    !!----------------------------------------------------------------------
(   55)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   56)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/restart.F90,v 1.15 2006/03/10 10:55:35 opalod Exp $ 
(   57)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 17

(   58)    !!----------------------------------------------------------------------
(   59) 
(   60) 
(   61) CONTAINS
(   62) 
(   63) 
(   64)    !!----------------------------------------------------------------------
(   65)    !!   Default option                                          NetCDF file
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68)    SUBROUTINE rst_write( kt )
(   69)       !!---------------------------------------------------------------------
(   70)       !!                   ***  ROUTINE rstwrite  ***
(   71)       !!                     
(   72)       !! ** Purpose :   Write restart fields in NetCDF format
(   73)       !!
(   74)       !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
(   75)       !!      file, save fields which are necessary for restart
(   76)       !!
(   77)       !! History :
(   78)       !!        !  99-11  (M. Imbard)  Original code
(   79)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   80)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   81)       !!----------------------------------------------------------------------
(   82)       !! * Modules used
(   83)       USE ioipsl
(   84) 
(   85)       !! * Arguments 
(   86)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   87) 
(   88)       !! * Local declarations
(   89)       LOGICAL ::   llbon
(   90)       CHARACTER (len=50) ::   clname, cln
(   91)       INTEGER ::   ic, jc, itime
(   92)       REAL(wp) ::   zdate0
(   93)       REAL(wp), DIMENSION( 1) ::   zfice, zfblk   ! used only in case of ice & bulk
(   94)       REAL(wp), DIMENSION(10) ::   zinfo(10)
(   95)       REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(   96)       
(   97)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
(   98)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
(   99)       INTEGER :: lncdf_stat                    ! lib_ncdf call return status
(  100)       CHARACTER(LEN=80) :: restfilename
(  101)       CHARACTER(LEN=4) :: chnum
(  102) 
(  103) # 106
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108)       IF( kt == nit000 ) THEN
(  109)          IF(lwp) WRITE(numout,*)
(  110)          IF(lwp) WRITE(numout,*) 'rst_wri : write restart.output NetCDF file'
(  111)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  112)          zfice(1) = 1.e0   ;   zfblk(1) = 1.e0
(  113)       ENDIF
(  114)       IF (kt == nit000) THEN
(  115)          write_count = 0
(  116)       END IF
(  117) 






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 18

(  118) 
(  119)       IF( MOD( kt, nstock ) == 0 .OR. kt == nitend ) THEN
(  120)          
(  121)          ! 0. Initializations
(  122)          ! ------------------
(  123) 
(  124)          IF(lwp) WRITE(numout,*) ' '
(  125)          IF(lwp) WRITE(numout,*) 'rst_write : write the restart file in NetCDF format ',   &
(  126)                                               'at it= ',kt,' date= ',ndastp
(  127)          IF(lwp) WRITE(numout,*) '~~~~~~~~~'
(  128) 
(  129)          ! Job informations
(  130)          zinfo(1) = FLOAT( no        )   ! job number
(  131)          zinfo(2) = FLOAT( kt        )   ! time-step
(  132)          zinfo(3) = FLOAT( 2 - nsolv )   ! pcg solver
(  133)          zinfo(4) = FLOAT( nsolv - 1 )   ! sor solver
(  134)          IF( lk_zdftke ) THEN
(  135)             zinfo(5) = 1.e0              ! TKE 
(  136)          ELSE
(  137)             zinfo(5) = 0.e0              ! no TKE 
(  138)          ENDIF
(  139)          zinfo(6) = FLOAT( ndastp )      ! date
(  140)          zinfo(7) = adatrj               ! number of elapsed days since the begining of the run
(  141) 
(  142) 
(  143) !!DB - code deleted ...
(  144) !         IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL?
(  145)       ! LIB_NCDF CALLS
(  146) 
(  147) !DB 
(  148) !      IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf?
(  149) !         WRITE(UNIT=restfilename,FMT="('rstfile_',I4.4,'.nc')") write_count
(  150) !!DB
(  151)          write(chnum,'(I4.4)')write_count
(  152)          restfilename = trim(cexper)//'_rstfile_'//chnum//'.nc'
(  153) 
(  154)          CALL ncdf_create_restart(restfilename, lncdf_stat)
(  155)          CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   
(  156)          CALL ncdf_write(restfilename, 'ub', ub, nwrite, lncdf_stat) 
(  157)          CALL ncdf_write(restfilename, 'vb', vb, nwrite, lncdf_stat)  
(  158)          CALL ncdf_write(restfilename, 'tb', tb, nwrite, lncdf_stat) 
(  159)          CALL ncdf_write(restfilename, 'sb', sb, nwrite, lncdf_stat)
(  160)          CALL ncdf_write(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
(  161)          CALL ncdf_write(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
(  162)          CALL ncdf_write(restfilename, 'un', un, nwrite, lncdf_stat)
(  163)          CALL ncdf_write(restfilename, 'vn', vn, nwrite, lncdf_stat)
(  164)          CALL ncdf_write(restfilename, 'tn', tn, nwrite, lncdf_stat) 
(  165)          CALL ncdf_write(restfilename, 'sn', sn, nwrite, lncdf_stat)
(  166)          CALL ncdf_write(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
(  167)          CALL ncdf_write(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
(  168)          ztab(:,:) = gcx(1:jpi,1:jpj)
(  169)          CALL ncdf_write(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
(  170)          ztab(:,:) = gcxb(1:jpi,1:jpj)
(  171)          CALL ncdf_write(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)       
(  172) # 177
(  177)          CALL ncdf_write(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
(  178)          CALL ncdf_write(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
(  179) # 180






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 19

(  180)          CALL ncdf_write(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
(  181)          CALL ncdf_write(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
(  182)          CALL ncdf_write(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
(  183)          CALL ncdf_write(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
(  184) # 187
(  187)          IF( lk_zdftke ) THEN
(  188)             CALL ncdf_write(restfilename, 'en', en, nwrite, lncdf_stat)   ! TKE arrays
(  189)          ENDIF
(  190) # 192
(  192)           zfice(1) = FLOAT( nfice )                                      ! Louvain La Neuve Sea Ice Model
(  193) !ADbug          CALL ncdf_write(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
(  194)          CALL ncdf_write(restfilename, 'nfice', zfice, lncdf_stat)
(  195)          CALL ncdf_write(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
(  196)          CALL ncdf_write(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
(  197)          CALL ncdf_write(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
(  198)          CALL ncdf_write(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
(  199) # 204
(  204)          zfblk(1) = FLOAT( nfbulk )                                 ! Bulk
(  205) ! ADbug         CALL ncdf_write(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
(  206)          CALL ncdf_write(restfilename, 'nfbulk', zfblk, lncdf_stat)
(  207)          CALL ncdf_write(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
(  208) 
(  209) # 210
(  210)          CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
(  211)          CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
(  212)          CALL ncdf_write(restfilename, 'nav_lev', gdept, lncdf_stat)
(  213)          CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  214)          CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  215)          write_count = write_count + 1
(  216) 
(  217) !      END IF ! End lib_ncdf block
(  218) 
(  219)       ! END LIB_NCDF CALLS
(  220) 
(  221)    ENDIF
(  222) 
(  223)    END SUBROUTINE rst_write
(  224) 
(  225) 
(  226) 
(  227)    SUBROUTINE rst_read
(  228)      !!---------------------------------------------------------------------- 
(  229)      !!                   ***  ROUTINE rst_read  ***
(  230)      !! 
(  231)      !! ** Purpose :   Read files for restart
(  232)      !! 
(  233)      !! ** Method  :   Read the previous fields on the NetCDF file 
(  234)      !!      the first record indicates previous characterics
(  235)      !!      after control with the present run, we read :
(  236)      !!      - prognostic variables on the second record
(  237)      !!      - elliptic solver arrays 
(  238)      !!      - barotropic stream function arrays ("key_dynspg_rl" defined)
(  239)      !!        or free surface arrays 
(  240)      !!      - tke arrays (lk_zdftke=T)
(  241)      !!      for this last three records,  the previous characteristics 
(  242)      !!      could be different with those used in the present run. 
(  243)      !!
(  244)      !!   According to namelist parameter nrstdt,






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 20

(  245)      !!       nrstdt = 0  no control on the date (nit000 is  arbitrary).
(  246)      !!       nrstdt = 1  we verify that nit000 is equal to the last
(  247)      !!                   time step of previous run + 1.
(  248)      !!       In both those options, the  exact duration of the experiment
(  249)      !!       since the beginning (cumulated duration of all previous restart runs)
(  250)      !!       is not stored in the restart and is assumed to be (nit000-1)*rdt.
(  251)      !!       This is valid is the time step has remained constant.
(  252)      !!
(  253)      !!       nrstdt = 2  the duration of the experiment in days (adatrj)
(  254)      !!                    has been stored in the restart file.
(  255)      !!
(  256)      !! History :
(  257)      !!        !  99-05  (M. Imbard)  Original code
(  258)      !!   8.5  !  02-09  (G. Madec)  F90: Free form
(  259)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  260)      !!----------------------------------------------------------------------
(  261)      !! * Modules used
(  262)      USE ioipsl
(  263)      
(  264)      !! * Local declarations
(  265)      LOGICAL ::   llog
(  266)      CHARACTER (len=8 ) ::   clvnames(50)
(  267)      CHARACTER (len=32) ::   clname
(  268)      INTEGER  ::   &
(  269)           itime, ibvar,     &  !
(  270)           inum                 ! temporary logical unit
(  271)      REAL(wp) ::   zdate0, zdt, zinfo(10)
(  272)      REAL(wp) ::   zdept(jpk), zlamt(jpi,jpj), zphit(jpi,jpj)
(  273)      REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(  274) # 275
(  275)      INTEGER  ::   ios1, ji, jj, jn
(  276)      REAL(wp) ::   zfice(1)
(  277) # 279
(  279)      INTEGER  ::   ios2, jk
(  280)      REAL(wp) ::   zfblk(1)
(  281) # 282
(  282)      
(  283)      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.   ! Use IOIPSL subroutines for restart input
(  284)      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.  ! Use lib_ncdf subroutines for restart input
(  285)      INTEGER :: lncdf_stat                     ! lib_ncdf routine return status
(  286)      REAL :: tstp
(  287)      CHARACTER(LEN=40) :: restfilename="restart.nc" !Default filename 
(  288)      
(  289)      
(  290)      !!----------------------------------------------------------------------
(  291)      !!  OPA 8.5, LODYC-IPSL (2002)
(  292)      !!----------------------------------------------------------------------
(  293)      clname = 'restart'
(  294) # 298
(  298)      
(  299)      IF(lwp) WRITE(numout,*)
(  300)      IF(lwp) WRITE(numout,*) 'rst_read : read the NetCDF restart file'
(  301)      IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  302)      
(  303)      IF(lwp) WRITE(numout,*) ' Info on the present job : '
(  304)      IF(lwp) WRITE(numout,*) '   job number          : ', no
(  305)      IF(lwp) WRITE(numout,*) '   time-step           : ', nit000
(  306)      IF(lwp) WRITE(numout,*) '   solver type         : ', nsolv






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 21

(  307)      IF( lk_zdftke ) THEN
(  308)         IF(lwp) WRITE(numout,*) '   tke option          : 1 '
(  309)      ELSE
(  310)         IF(lwp) WRITE(numout,*) '   tke option          : 0 '
(  311)      ENDIF
(  312)      IF(lwp) WRITE(numout,*) '   date ndastp         : ', ndastp
(  313)      IF(lwp) WRITE(numout,*)
(  314)      
(  315)      ! Time domain : restart
(  316)      ! -------------------------
(  317)      
(  318)      IF(lwp) WRITE(numout,*)
(  319)      IF(lwp) WRITE(numout,*)
(  320)      IF(lwp) WRITE(numout,*) ' *** restart option'
(  321)      SELECT CASE ( nrstdt )
(  322)      CASE ( 0 ) 
(  323)         IF(lwp) WRITE(numout,*) ' nrstdt = 0 no control of nit000'
(  324)      CASE ( 1 ) 
(  325)         IF(lwp) WRITE(numout,*) ' nrstdt = 1 we control the date of nit000'
(  326)      CASE ( 2 )
(  327)         IF(lwp) WRITE(numout,*) ' nrstdt = 2 the date adatrj is read in restart file'
(  328)      CASE DEFAULT
(  329)         IF(lwp) WRITE(numout,*) '  ===>>>> nrstdt not equal 0, 1 or 2 : no control of the date'
(  330)         IF(lwp) WRITE(numout,*) ' =======                   ========='
(  331)      END SELECT
(  332)      
(  333)      itime = 0
(  334)      llog  = .FALSE.
(  335)      zlamt(:,:) = 0.e0
(  336)      zphit(:,:) = 0.e0
(  337)      zdept(:)   = 0.e0
(  338) 
(  339) !!DB -- code deleted ...
(  340) !     IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Do we want to use IOIPSL?
(  341) 
(  342)      ! LIB_NCDF CALLS
(  343) !!DB
(  344) !     IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Do we want to use lib_ncdf?
(  345) !        restfilename="restart.nc"
(  346)      CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
(  347) !!DB -- check if found file
(  348)      if(lwp .AND. lncdf_stat /= 0)then
(  349)         write(numout,*)'                '
(  350)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  351)         write(numout,*)'STOP: Problem reading restart file ', restfilename
(  352)         write(numout,*)'narea = ', narea
(  353)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  354)         write(numout,*)'                '
(  355)         stop
(  356)      endif
(  357) 
(  358)      
(  359)      ! NON LIB_NCDF RELATED CODE
(  360)      IF(lwp) WRITE(numout,*)
(  361)      IF(lwp) WRITE(numout,*) ' Info on the restart file read : '
(  362)      IF(lwp) WRITE(numout,*) '   job number          : ', NINT( zinfo(1) )
(  363)      IF(lwp) WRITE(numout,*) '   time-step           : ', NINT( zinfo(2) )
(  364)      IF(lwp) WRITE(numout,*) '   solver type         : ', NINT( zinfo(4) ) + 1






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 22

(  365)      IF(lwp) WRITE(numout,*) '   tke option          : ', NINT( zinfo(5) )
(  366)      IF(lwp) WRITE(numout,*) '   date ndastp         : ', NINT( zinfo(6) )
(  367)      IF(lwp) WRITE(numout,*)
(  368)      
(  369)      ! Control of date
(  370)      IF( nit000 - NINT( zinfo(2) )  /= 1 .AND. nrstdt /= 0 ) THEN
(  371)         IF(lwp) WRITE(numout,cform_err)
(  372)         IF(lwp) WRITE(numout,*) ' ===>>>> : problem with nit000 for the restart'
(  373)         IF(lwp) WRITE(numout,*) ' verify the restart file or rerun with nrstdt = 0 (namelist)'
(  374)         nstop = nstop + 1
(  375)      ENDIF
(  376)      
(  377)      ! re-initialisation of  adatrj0
(  378)      adatrj0 =  ( FLOAT( nit000-1 ) * rdttra(1) ) / rday
(  379)      
(  380)      IF ( nrstdt == 2 ) THEN
(  381)         !                             by default ndatsp has been set to ndate0 in dom_nam
(  382)         !                             ndate0 has been read in the namelist (standard OPA 8)
(  383)         !                             here when nrstdt=2 we keep the  final date of previous run
(  384)         ndastp = NINT( zinfo(6) )
(  385)         adatrj0 =  zinfo(7)
(  386)      ENDIF
(  387)      ! END NON LIB_NCDF RELATED CODE
(  388)      
(  389)      ! End of job info stuff
(  390)      CALL ncdf_read(restfilename, 'ub', ub, nwrite, lncdf_stat)   
(  391)      if(lncdf_stat /= 0) stop
(  392)      CALL ncdf_read(restfilename, 'vb', vb, nwrite, lncdf_stat)  
(  393)      if(lncdf_stat /= 0) stop
(  394)      CALL ncdf_read(restfilename, 'tb', tb, nwrite, lncdf_stat) 
(  395)      if(lncdf_stat /= 0) stop
(  396)      CALL ncdf_read(restfilename, 'sb', sb, nwrite, lncdf_stat)
(  397)      if(lncdf_stat /= 0) stop
(  398)      CALL ncdf_read(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
(  399)      if(lncdf_stat /= 0) stop
(  400)      CALL ncdf_read(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
(  401)      if(lncdf_stat /= 0) stop
(  402)      CALL ncdf_read(restfilename, 'un', un, nwrite, lncdf_stat)
(  403)      if(lncdf_stat /= 0) stop
(  404)      CALL ncdf_read(restfilename, 'vn', vn, nwrite, lncdf_stat)
(  405)      if(lncdf_stat /= 0) stop
(  406)      CALL ncdf_read(restfilename, 'tn', tn, nwrite, lncdf_stat) 
(  407)      if(lncdf_stat /= 0) stop
(  408)      CALL ncdf_read(restfilename, 'sn', sn, nwrite, lncdf_stat)
(  409)      if(lncdf_stat /= 0) stop
(  410)      CALL ncdf_read(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
(  411)      if(lncdf_stat /= 0) stop
(  412)      CALL ncdf_read(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
(  413)      if(lncdf_stat /= 0) stop
(  414)      CALL ncdf_read(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)    
(  415)      if(lncdf_stat /= 0) stop
(  416)      gcxb(1:jpi,1:jpj) = ztab(:,:) 
(  417)      CALL ncdf_read(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
(  418)      if(lncdf_stat /= 0) stop
(  419)      gcx(1:jpi,1:jpj) = ztab(:,:)   
(  420) # 428
(  428)      CALL ncdf_read(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
(  429)      if(lncdf_stat /= 0) stop






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 23

(  430)      CALL ncdf_read(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
(  431)      if(lncdf_stat /= 0) stop
(  432) # 433
(  433)      CALL ncdf_read(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
(  434)      if(lncdf_stat /= 0) stop
(  435)      CALL ncdf_read(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
(  436)      if(lncdf_stat /= 0) stop
(  437)      CALL ncdf_read(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
(  438)      if(lncdf_stat /= 0) stop
(  439)      CALL ncdf_read(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
(  440)      if(lncdf_stat /= 0) stop
(  441) # 444
(  444)      IF( lk_zdftke ) THEN
(  445)         IF( NINT( zinfo(5) ) == 1 ) THEN                                ! Read tke arrays
(  446)            CALL ncdf_read(restfilename, 'en', en, nwrite, lncdf_stat)
(  447)            ln_rstke = .FALSE.
(  448)         ELSE
(  449)            IF(lwp) WRITE(numout,*) ' ===>>>> : the previous restart file didnot used  tke scheme'
(  450)            IF(lwp) WRITE(numout,*) ' =======                ======='
(  451)            nrstdt = 2
(  452)            ln_rstke = .TRUE.
(  453)         ENDIF
(  454)      ENDIF
(  455) # 457
(  457)       ! Louvain La Neuve Sea Ice Model
(  458)       ios1 = 0
(  459)       DO jn = 1, 30
(  460)          IF( clvnames(jn) == 'nfice' )  ios1 = 1
(  461)       END DO
(  462)       IF( ios1 == 1 ) THEN
(  463) !ADbug          CALL ncdf_read(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
(  464)          CALL ncdf_read(restfilename, 'nfice', zfice, lncdf_stat)
(  465)          CALL ncdf_read(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
(  466)          CALL ncdf_read(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
(  467)          CALL ncdf_read(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
(  468)          CALL ncdf_read(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
(  469) # 472
(  472)       ENDIF
(  473)       IF( zfice(1) /= FLOAT(nfice) .OR. ios1 == 0 ) THEN
(  474)          IF(lwp) WRITE(numout,*)
(  475)          IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
(  476)          IF(lwp) WRITE(numout,*)
(  477)          sst_io(:,:) = ( nfice-1 )*( tn(:,:,1) + rt0 )          !!bug a explanation is needed here!
(  478)          sss_io(:,:) = ( nfice-1 )*  sn(:,:,1)
(  479)          DO jj = 2, jpj
(  480)             DO ji = 2, jpi
(  481)                u_io(ji,jj) = ( nfice-1 ) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  482)                v_io(ji,jj) = ( nfice-1 ) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  483)             END DO
(  484)          END DO
(  485) # 488
(  488)       ENDIF
(  489) # 491
(  491)       ! Louvain La Neuve Sea Ice Model
(  492)       ios2 = 0
(  493)       DO jk = 1, 30
(  494)          IF( clvnames(jk) == 'nfbulk' )  ios2 = 1
(  495)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 24

(  496)       IF( ios2 == 1 ) THEN
(  497) !ADbug          CALL ncdf_read(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
(  498)          CALL ncdf_read(restfilename, 'nfbulk', zfblk, lncdf_stat)
(  499)          CALL ncdf_read(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
(  500)       ENDIF
(  501)       IF( zfblk(1) /= FLOAT(nfbulk) .OR. ios2 == 0 ) THEN
(  502)          IF(lwp) WRITE(numout,*)
(  503)          IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
(  504)          IF(lwp) WRITE(numout,*)
(  505)          gsst(:,:) = 0.
(  506)          gsst(:,:) = gsst(:,:) + ( nfbulk-1 )*( tn(:,:,1) + rt0 )
(  507)       ENDIF
(  508) 
(  509) # 510
(  510)      CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)
(  511)      if(lncdf_stat /= 0) stop
(  512)      CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
(  513)      if(lncdf_stat /= 0) stop
(  514)      CALL ncdf_read(restfilename, 'nav_lev', gdept, lncdf_stat)
(  515)      if(lncdf_stat /= 0) stop
(  516)      tstp = REAL(itime)
(  517)      CALL ncdf_read(restfilename, 'time', tstp, 1, lncdf_stat)
(  518)      if(lncdf_stat /= 0) stop
(  519)      itime = NINT(tstp)
(  520)      CALL ncdf_read(restfilename, 'time_steps', zdate0, 1, lncdf_stat)
(  521)      if(lncdf_stat /= 0) stop
(  522)   ! END LIB_NCDF CALLS
(  523)   
(  524)   ! In case of restart with neuler = 0 then put all before fields = to now fields
(  525)   IF ( neuler == 0 ) THEN
(  526)      tb(:,:,:)=tn(:,:,:)
(  527)      sb(:,:,:)=sn(:,:,:)
(  528)      ub(:,:,:)=un(:,:,:)
(  529)      vb(:,:,:)=vn(:,:,:)
(  530)      rotb(:,:,:)=rotn(:,:,:)
(  531)      hdivb(:,:,:)=hdivn(:,:,:)
(  532) # 536
(  536)      ! free surface formulation (eta)
(  537)      sshb(:,:)=sshn(:,:)
(  538) # 539
(  539)   ENDIF
(  540)   
(  541) END SUBROUTINE rst_read
(  542) 
(  543) 
(  544) 
(  545) 
(  546) 
(  547) # 548
(  548) !!DB 2008.05.22 
(  549) !!NB: called when needed by icestp
(  550) !!         IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
(  551) !            CALL lim_rst_write( numit )                              ! Ice restart file !
(  552) !         ENDIF
(  553) 
(  554)    SUBROUTINE rst_ice_write( niter )
(  555)       !!---------------------------------------------------------------------
(  556)       !!                   ***  ROUTINE rst ice write  ***






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 25

(  557)       !!                     
(  558)       !! ** Purpose :   Write lim ice restart fields in NetCDF format
(  559)       !!
(  560)       !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
(  561)       !!      file, save fields which are necessary for restart
(  562)       !!
(  563)       !! History : DB  2008.05.22
(  564)       !!
(  565)       !!----------------------------------------------------------------------
(  566)       !! * Modules used
(  567)       USE ioipsl
(  568) 
(  569)       !! * Arguments 
(  570)       INTEGER, INTENT( in ) ::   niter         ! ice iteration 
(  571) 
(  572)       !! * Local declarations
(  573)       LOGICAL ::   llbon
(  574)       CHARACTER (len=50) ::   clname, cln
(  575)       INTEGER ::   ic, jc, itime, it0
(  576)       REAL(wp), DIMENSION(2) ::   zinfo(2)
(  577)       REAL(wp), DIMENSION(jpi,jpj) :: ztab 
(  578)       REAL(wp) ::  zsec, zdate0, zdt
(  579)       
(  580)       LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
(  581)       LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
(  582)       INTEGER :: lncdf_stat                    ! lib_ncdf call return status
(  583)       CHARACTER(LEN=80) :: restfilename
(  584)       CHARACTER(LEN=4) :: chnum
(  585) 
(  586)       REAL(wp),DIMENSION(jpi,jpj,35) :: zmoment
(  587)       INTEGER :: ji, jj
(  588) 
(  589) # 592
(  592)       !!----------------------------------------------------------------------
(  593) 
(  594) 
(  595) ! Job informations
(  596)       it0      = niter
(  597)       zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
(  598)       zinfo(2) = FLOAT( it0   )   ! iteration number
(  599) 
(  600)       zsec     = 0.e0
(  601)       itime    = 0
(  602) !      zdept(1) = 0.e0   !! this is default z-level -- handled in ncdf_create_ice_restart()
(  603)       zdt      = rdt_ice * nstock  !!does not seem to enter restart file
(  604) 
(  605)       IF(lwp) WRITE(numout,*) ' '
(  606)       IF(lwp) WRITE(numout,*) 'ice_rst_write : write the ice restart file in NetCDF format ',   &
(  607)            'at it= ',niter,' date= ',ndastp
(  608)       IF(lwp) WRITE(numout,*) '~~~~~~~~~'
(  609) 
(  610) 
(  611) !!DB
(  612)       write(chnum,'(I4.4)')ice_write_count
(  613)       restfilename = trim(cexper)//'_ice_rstfile_'//chnum//'.nc'
(  614) 
(  615)       CALL ncdf_create_ice_restart(restfilename, lncdf_stat)
(  616)       CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 26

(  617) 
(  618)       CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
(  619)       CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
(  620) !      CALL ncdf_write(restfilename, 'nav_lev', real(0), nwrite, lncdf_stat)   !!nav_lev treated as a scalar
(  621)       CALL ncdf_write(restfilename, 'nav_lev', REAL(0), 1, lncdf_stat)   !!nav_lev treated as a scalar
(  622)       CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  623)       CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  624) 
(  625)       CALL ncdf_write(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
(  626)       CALL ncdf_write(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)
(  627)       CALL ncdf_write(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
(  628)       CALL ncdf_write(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
(  629) # 632
(  632)       CALL ncdf_write(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
(  633)       CALL ncdf_write(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
(  634)       CALL ncdf_write(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
(  635)       CALL ncdf_write(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)
(  636)       CALL ncdf_write(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
(  637)       CALL ncdf_write(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)
(  638)       CALL ncdf_write(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)
(  639) 
(  640) 
(  641) !!DB load zmoment
(  642)       DO jj = 1, jpj  
(  643)          DO ji = 1, jpi
(  644)             zmoment(ji,jj,1)  = sxice(ji,jj)
(  645)             zmoment(ji,jj,2)  = syice(ji,jj)
(  646)             zmoment(ji,jj,3)  = sxxice(ji,jj)
(  647)             zmoment(ji,jj,4)  = syyice(ji,jj)
(  648)             zmoment(ji,jj,5)  = sxyice(ji,jj)
(  649)             zmoment(ji,jj,6)  = sxsn(ji,jj)
(  650)             zmoment(ji,jj,7)  = sysn(ji,jj)
(  651)             zmoment(ji,jj,8)  = sxxsn(ji,jj)
(  652)             zmoment(ji,jj,9)  = syysn(ji,jj)
(  653)             zmoment(ji,jj,10) = sxysn(ji,jj)
(  654)             zmoment(ji,jj,11) = sxa(ji,jj)
(  655)             zmoment(ji,jj,12) = sya(ji,jj)
(  656)             zmoment(ji,jj,13) = sxxa(ji,jj)
(  657)             zmoment(ji,jj,14) = syya(ji,jj)
(  658)             zmoment(ji,jj,15) = sxya(ji,jj)
(  659)             zmoment(ji,jj,16) = sxc0(ji,jj)
(  660)             zmoment(ji,jj,17) = syc0(ji,jj)
(  661)             zmoment(ji,jj,18) = sxxc0(ji,jj)
(  662)             zmoment(ji,jj,19) = syyc0(ji,jj)
(  663)             zmoment(ji,jj,20) = sxyc0(ji,jj)
(  664)             zmoment(ji,jj,21) = sxc1(ji,jj)
(  665)             zmoment(ji,jj,22) = syc1(ji,jj)
(  666)             zmoment(ji,jj,23) = sxxc1(ji,jj)
(  667)             zmoment(ji,jj,24) = syyc1(ji,jj)
(  668)             zmoment(ji,jj,25) = sxyc1(ji,jj)
(  669)             zmoment(ji,jj,26) = sxc2(ji,jj)
(  670)             zmoment(ji,jj,27) = syc2(ji,jj)
(  671)             zmoment(ji,jj,28) = sxxc2(ji,jj)
(  672)             zmoment(ji,jj,29) = syyc2(ji,jj)
(  673)             zmoment(ji,jj,30) = sxyc2(ji,jj)
(  674)             zmoment(ji,jj,31) = sxst(ji,jj)
(  675)             zmoment(ji,jj,32) = syst(ji,jj)
(  676)             zmoment(ji,jj,33) = sxxst(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 27

(  677)             zmoment(ji,jj,34) = syyst(ji,jj)
(  678)             zmoment(ji,jj,35) = sxyst(ji,jj)
(  679)          END DO
(  680)       END DO
(  681) !!DB probably should call mppsync here ??? NB: Works if I do not
(  682) !      CALL MPI_BARRIER(MPI_COMM_WORLD, lncdf_stat)
(  683)       CALL ncdf_write(restfilename, 'moment', zmoment, nwrite, lncdf_stat)
(  684) 
(  685)       ice_write_count = ice_write_count + 1
(  686) 
(  687) 
(  688)    END SUBROUTINE rst_ice_write
(  689) 
(  690) 
(  691) !!DB: 2008.05.27 read ice restart 
(  692) !!Tested to run OK, but there is no proof that it did the right thing
(  693) !!NB; In general, there is no proof that ice restarts are perfect
(  694) !!Similar to old lim_rst_read()
(  695) !!NB: default restfilename = restart_ice.nc
(  696) !!NB: Agrif stuff not done yet
(  697)    SUBROUTINE rst_ice_read( niter )
(  698)       !-----------------------------------------------------------------------
(  699)       ! Arguments
(  700)       INTEGER  ::   niter        ! number of iteration
(  701) 
(  702)       !- dummy variables :
(  703)       CHARACTER(len=45)  ::  restfilename = 'restart_ice.nc'
(  704)       INTEGER :: &
(  705)         ji, jj, lncdf_stat
(  706)       INTEGER :: &
(  707)          inumrst, it0, it1, itime, ibvar, ifice
(  708)       LOGICAL :: &
(  709)          llog
(  710)       REAL(wp),DIMENSION(jpi,jpj) :: &
(  711)          zlamt, zphit
(  712)       REAL(wp),DIMENSION(jpi,jpj,35) :: &
(  713)          zmoment
(  714)       REAL(wp),DIMENSION(1) :: &
(  715)          zdept
(  716)       REAL(wp),DIMENSION(2) :: &
(  717)          zinfo
(  718)       REAL(wp) :: &
(  719)          zdate0, zdt
(  720)       CHARACTER ( len = 10 ) ::  &
(  721)          clvnames(60)       
(  722) 
(  723) 
(  724) 
(  725) 
(  726) # 731
(  731)      CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
(  732) !!DB -- check if found file
(  733)      if(lwp .AND. lncdf_stat /= 0)then
(  734)         write(numout,*)'                '
(  735)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
(  736)         write(numout,*)'STOP: Problem reading ice restart file ', restfilename
(  737)         write(numout,*)'narea = ', narea
(  738)         write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 28

(  739)         write(numout,*)'                '
(  740)         stop
(  741)      endif
(  742) 
(  743) 
(  744) 
(  745)       !Initialisations
(  746)       inumrst    = 71
(  747)       it0        = nit000
(  748)       itime      = 0
(  749)       llog       = .FALSE.
(  750)       zlamt(:,:) = 0.
(  751)       zphit(:,:) = 0.
(  752)       zdept(1)   = 0.
(  753) 
(  754)       ifice   = INT( zinfo(1) )
(  755)       it1     = INT( zinfo(2) )
(  756) 
(  757) !!DB Note that old code gets varnames but otherwise does not use this info 
(  758)       IF(lwp) THEN
(  759)          WRITE(numout,*)
(  760)          WRITE(numout,*) 'lim_rst_read : READ restart file name ', restfilename,' at time step : ', it1
(  761)       ENDIF
(  762)       
(  763)       !Control of date
(  764) !!DB: 2009.09.30
(  765)       IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
(  766)          IF(lwp) THEN
(  767)             WRITE(numout,cform_err)
(  768)             WRITE(numout,*) 'DBG: lim_rst_read nrstdt=1 and nit000 /= zinfo(2)'
(  769)             WRITE(numout,*) 'DBG: This case covered in restart() ====> Continue'
(  770) !            WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
(  771) !            WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
(  772) !            WRITE(numout,*) '   control of time parameter  nrstdt'
(  773) !            nstop = nstop + 1
(  774)          ENDIF
(  775)       ENDIF
(  776) 
(  777) !!DBG: writes from rst_ice_write()
(  778)      CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)
(  779)      if(lncdf_stat /= 0) stop
(  780)      CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
(  781)      if(lncdf_stat /= 0) stop
(  782) 
(  783) !!DBG ignore these vars to start
(  784) !      CALL ncdf_read(restfilename, 'nav_lev', zdept, 1, lncdf_stat)   !!nav_lev treated as a scalar
(  785) !      CALL ncdf_read(restfilename, 'time', REAL(0), 1, lncdf_stat)
(  786) !      CALL ncdf_read(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
(  787) 
(  788)       CALL ncdf_read(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
(  789)       CALL ncdf_read(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)
(  790)       CALL ncdf_read(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
(  791)       CALL ncdf_read(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
(  792) # 795
(  795)       CALL ncdf_read(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
(  796)       CALL ncdf_read(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
(  797)       CALL ncdf_read(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
(  798)       CALL ncdf_read(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 29

(  799)       CALL ncdf_read(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
(  800)       CALL ncdf_read(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)
(  801)       CALL ncdf_read(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)
(  802)       CALL ncdf_read(restfilename, 'moment', zmoment, nwrite, lncdf_stat)
(  803) 
(  804) !!DB: 2009.09.30
(  805) !      niter = it1
(  806)       niter = nit000
(  807) 
(  808)       DO jj = 1, jpj
(  809)          DO ji = 1, jpi
(  810)             sxice(ji,jj)  = zmoment(ji,jj,1)
(  811)             syice(ji,jj)  = zmoment(ji,jj,2)
(  812)             sxxice(ji,jj) = zmoment(ji,jj,3)
(  813)             syyice(ji,jj) = zmoment(ji,jj,4)
(  814)             sxyice(ji,jj) = zmoment(ji,jj,5)
(  815)             sxsn(ji,jj)   = zmoment(ji,jj,6)
(  816)             sysn(ji,jj)   = zmoment(ji,jj,7)
(  817)             sxxsn(ji,jj)  = zmoment(ji,jj,8)
(  818)             syysn(ji,jj)  = zmoment(ji,jj,9)
(  819)             sxysn(ji,jj)  = zmoment(ji,jj,10)
(  820)             sxa(ji,jj)    = zmoment(ji,jj,11)
(  821)             sya(ji,jj)    = zmoment(ji,jj,12)
(  822)             sxxa(ji,jj)   = zmoment(ji,jj,13)
(  823)             syya(ji,jj)   = zmoment(ji,jj,14)
(  824)             sxya(ji,jj)   = zmoment(ji,jj,15)
(  825)             sxc0(ji,jj)   = zmoment(ji,jj,16)
(  826)             syc0(ji,jj)   = zmoment(ji,jj,17)
(  827)             sxxc0(ji,jj)  = zmoment(ji,jj,18)
(  828)             syyc0(ji,jj)  = zmoment(ji,jj,19)
(  829)             sxyc0(ji,jj)  = zmoment(ji,jj,20)
(  830)             sxc1(ji,jj)   = zmoment(ji,jj,21)
(  831)             syc1(ji,jj)   = zmoment(ji,jj,22)
(  832)             sxxc1(ji,jj)  = zmoment(ji,jj,23)
(  833)             syyc1(ji,jj)  = zmoment(ji,jj,24)
(  834)             sxyc1(ji,jj)  = zmoment(ji,jj,25)
(  835)             sxc2(ji,jj)   = zmoment(ji,jj,26)
(  836)             syc2(ji,jj)   = zmoment(ji,jj,27)
(  837)             sxxc2(ji,jj)  = zmoment(ji,jj,28)
(  838)             syyc2(ji,jj)  = zmoment(ji,jj,29)
(  839)             sxyc2(ji,jj)  = zmoment(ji,jj,30)
(  840)             sxst(ji,jj)   = zmoment(ji,jj,31)
(  841)             syst(ji,jj)   = zmoment(ji,jj,32)
(  842)             sxxst(ji,jj)  = zmoment(ji,jj,33)
(  843)             syyst(ji,jj)  = zmoment(ji,jj,34)
(  844)             sxyst(ji,jj)  = zmoment(ji,jj,35)
(  845)          END DO
(  846)       END DO
(  847) 
(  848)       
(  849)    END SUBROUTINE rst_ice_read
(  850) 
(  851) 
(  852) 
(  853) 
(  854) 
(  855) # 856
(  856) !!=====================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:51      page 30

(  857) END MODULE restart
(  858) 
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (restart.F90: 130)
