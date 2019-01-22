


PGF90 (Version     10.5)          09/23/2010  09:49:07      page 1

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

Filename: dtatem.F90

(    1) !!DB: 2008.10.01 -- Routine modified 
(    2) !!-1- deleted reference to orca
(    3) !!-2- replaced IOIPSL with lib_ncdf 
(    4) 
(    5) 
(    6) MODULE dtatem
(    7)    !!======================================================================
(    8)    !!                     ***  MODULE  dtatem  ***
(    9)    !! Ocean data  :  read ocean temperature data from monthly atlas data
(   10)    !!=====================================================================
(   11) #if defined key_dtatem   ||   defined key_esopa
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_dtatem'                              3D temperature data field
(   14)    !!----------------------------------------------------------------------
(   15)    !!   dta_tem      : read ocean temperature data
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce             ! ocean dynamics and tracers
(   19)    USE dom_oce         ! ocean space and time domain
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE daymod          ! calendar
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dta_tem   ! called by step.F90 and inidta.F90
(   28) 
(   29)    !! * Shared module variables
(   30)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .TRUE.   !: temperature data flag
(   31)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   32)       t_dta             !: temperature data at given time-step
(   33) 
(   34)    !! * Module variables
(   35)    CHARACTER (len=45) ::   &
(   36)       cl_tdata
(   37)    INTEGER ::   &
(   38)       nlecte =  0,   &  ! switch for the first read
(   39)       ntem1      ,   &  ! first record used
(   40)       ntem2             ! second record used
(   41)    REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
(   42)       temdta            ! temperature data at two consecutive times
(   43) 
(   44)    !! * Substitutions
(   45) #  include "domzgr_substitute.h90"
(   46)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 2

(   47)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   48)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtatem.F90,v 1.12 2006/04/19 14:43:15 opalod Exp $ 
(   49)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)    !!----------------------------------------------------------------------
(   55)    !!   Default case                                            NetCDF file
(   56)    !!----------------------------------------------------------------------
(   57) 
(   58)    SUBROUTINE dta_tem( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                   ***  ROUTINE dta_tem  ***
(   61)       !!                    
(   62)       !! ** Purpose :   Reads monthly temperature data 
(   63)       !! 
(   64)       !! ** Method  :   Read on unit numtdt the interpolated temperature 
(   65)       !!      onto the model grid.
(   66)       !!      Data begin at january. 
(   67)       !!      The value is centered at the middle of month.
(   68)       !!      In the opa model, kt=1 agree with january 1.
(   69)       !!      At each time step, a linear interpolation is applied between 
(   70)       !!      two monthly values.
(   71)       !!      Read on unit numtdt
(   72)       !!
(   73)       !! ** Action  :   define t_dta array at time-step kt
(   74)       !!
(   75)       !! History :
(   76)       !!        !  91-03  ()  Original code
(   77)       !!        !  92-07  (M. Imbard)
(   78)       !!        !  99-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
(   79)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   80)       !!----------------------------------------------------------------------
(   81)       !! * Modules used
(   82) !!DB
(   83)       USE lib_ncdf
(   84) 
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER, PARAMETER ::   &
(   90)          jpmois = 12                    ! number of month
(   91)       INTEGER ::   ji, jj, jl           ! dummy loop indicies
(   92)       INTEGER ::   &
(   93)          imois, iman, itime, ik ,    &  ! temporary integers
(   94)          i15, ipi, ipj, ipk             !    "          "
(   95) #  if defined key_tradmp
(   96)       INTEGER ::   &
(   97)          il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
(   98) # endif
(   99) 
(  100)       INTEGER, DIMENSION(jpmois) ::   istep
(  101)       REAL(wp) ::   zxy, zl, zdate0
(  102)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
(  103)       REAL(wp), DIMENSION(jpk) ::   zlev
(  104) !!DB






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 3

(  105)       INTEGER :: len, status 
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108)       ! 0. Initialization
(  109)       ! -----------------
(  110) 
(  111)       iman  = jpmois
(  112)       i15   = nday / 16
(  113)       imois = nmonth + i15 - 1
(  114)       IF( imois == 0 )   imois = iman
(  115) 
(  116)       itime = jpmois
(  117)       ipi = jpiglo
(  118)       ipj = jpjglo
(  119)       ipk = jpk
(  120) 
(  121)       ! 1. First call kt=nit000
(  122)       ! -----------------------
(  123) 
(  124)       IF( kt == nit000 .AND. nlecte == 0 ) THEN
(  125)          ntem1 = 0
(  126)          cl_tdata = 'data_1m_potential_temperature_nomask.nc'
(  127)          IF(lwp) WRITE(numout,*)
(  128)          IF(lwp) WRITE(numout,*) ' dtatem : Monthly Temperature fields'
(  129)          IF(lwp) WRITE(numout,*) ' ~~~~~~'
(  130)          IF(lwp) WRITE(numout,*) '             NetCDF File'
(  131)          IF(lwp) WRITE(numout,*)cl_tdata
(  132)          IF(lwp) WRITE(numout,*)
(  133)          
(  134) #if defined key_agrif
(  135)          if ( .NOT. Agrif_Root() ) then
(  136)             cl_tdata = TRIM(Agrif_CFixed())//'_'//TRIM(cl_tdata)
(  137)          endif
(  138) #endif         
(  139) !!DB: 
(  140)          call ncdf_get_dim_size(cl_tdata, 'time_counter', itime, status)
(  141)          call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
(  142)          call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
(  143)          call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
(  144)          if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
(  145)             if(lwp) then
(  146)                write(numout,*)
(  147)                write(numout,*) 'problem with dimensions'
(  148)                write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  149)                write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  150)                write(numout,*) ' ipk ',ipk,' jpk ',jpk
(  151)             endif
(  152)             stop 'dtatem'
(  153)          endif
(  154) !!DB: 
(  155)          if( itime /= jpmois ) then
(  156)             if(lwp) then
(  157)                write(numout,*)
(  158)                write(numout,*) 'problem with time coordinates'
(  159)                write(numout,*) ' itime ',itime,' jpmois ',jpmois
(  160)             endif
(  161)             stop 'dta_sal'
(  162)          endif






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 4

(  163) 
(  164)       ENDIF
(  165) 
(  166) 
(  167)       ! 2. Read monthly file
(  168)       ! -------------------
(  169) 
(  170)       IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
(  171)          nlecte = 1
(  172) 
(  173)          ! Calendar computation
(  174)          
(  175)          ntem1 = imois        ! first file record used 
(  176)          ntem2 = ntem1 + 1    ! last  file record used
(  177)          ntem1 = MOD( ntem1, iman )
(  178)          IF( ntem1 == 0 )   ntem1 = iman
(  179)          ntem2 = MOD( ntem2, iman )
(  180)          IF( ntem2 == 0 )   ntem2 = iman
(  181)          if(lwp) write(numout,*) 'dtatem reading records ',ntem1, ntem2
(  182)          
(  183)          ! Read monthly temperature data
(  184) !!DB
(  185)          call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,1),-ntem1,status)
(  186)          call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,2),-ntem2,status)
(  187) 
(  188) 
(  189)          IF(lwp) WRITE(numout,*)
(  190)          IF(lwp) WRITE(numout,*) ' Monthly temperature input OK'
(  191)          IF(lwp) WRITE(numout,*)
(  192)          
(  193) #if defined key_tradmp
(  194) !!DB -- orca-related code deleted
(  195) #endif
(  196) 
(  197)          !                                  ! Mask
(  198)          DO jl = 1, 2
(  199)             temdta(:,:,:,jl) = temdta(:,:,:,jl) * tmask(:,:,:)
(  200)             temdta(:,:,jpk,jl) = 0.
(  201)             IF( lk_zps ) THEN                ! z-coord. with partial steps
(  202)                DO jj = 1, jpj                  ! interpolation of temperature at the last level
(  203)                   DO ji = 1, jpi
(  204)                      ik = mbathy(ji,jj) - 1
(  205)                      IF( ik > 2 ) THEN
(  206)                         zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
(  207)                         temdta(ji,jj,ik,jl) = (1.-zl) * temdta(ji,jj,ik,jl) + zl * temdta(ji,jj,ik-1,jl) 
(  208)                      ENDIF
(  209)                   END DO
(  210)                END DO
(  211)             ENDIF
(  212)          END DO
(  213) 
(  214) !         IF(lwp) THEN
(  215) !            WRITE(numout,*) ' temperature month ', ntem1, ntem2
(  216) !            WRITE(numout,*)
(  217) !            WRITE(numout,*) '  month = ', ntem1, '  level = 1'
(  218) !            CALL prihre( temdta(:,:,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  219) !            WRITE(numout,*) '  month = ', ntem1, '  level = ', jpk/2
(  220) !            CALL prihre( temdta(:,:,jpk/2,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 5

(  221) !            WRITE(numout,*) '  month = ',ntem1,'  level = ', jpkm1
(  222) !            CALL prihre( temdta(:,:,jpkm1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  223) !         ENDIF
(  224)       ENDIF
(  225) 
(  226)  
(  227)       ! 2. At every time step compute temperature data
(  228)       ! ----------------------------------------------
(  229)       zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  230) !!DB: Note that for monthly data it is not necessary to do this every kt.
(  231) !!    So in the future use a code fragment like: 
(  232) !      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
(  233) !         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
(  234) !      endif
(  235) 
(  236)       t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
(  237) 
(  238) 
(  239)    END SUBROUTINE dta_tem
(  240) 
(  241) #else
(  242)    !!----------------------------------------------------------------------
(  243)    !!   Default case                           NO 3D temperature data field
(  244)    !!----------------------------------------------------------------------
(  245)    USE in_out_manager
(  246)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .FALSE.   !: temperature data flag
(  247) CONTAINS
(  248)    SUBROUTINE dta_tem( kt )        ! Empty routine
(  249)       if(lwp) WRITE(numout,*) 'dta_tem: You should not have seen this print! error?', kt
(  250)    END SUBROUTINE dta_tem
(  251) #endif
(  252)    !!======================================================================
(  253) END MODULE dtatem































PGF90 (Version     10.5)          09/23/2010  09:49:07      page 6

(    1) # 1 "dtatem.F90"
(    1) !!DB: 2008.10.01 -- Routine modified 
(    2) !!-1- deleted reference to orca
(    3) !!-2- replaced IOIPSL with lib_ncdf 
(    4) 
(    5) 
(    6) MODULE dtatem
(    7)    !!======================================================================
(    8)    !!                     ***  MODULE  dtatem  ***
(    9)    !! Ocean data  :  read ocean temperature data from monthly atlas data
(   10)    !!=====================================================================
(   11) # 12
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_dtatem'                              3D temperature data field
(   14)    !!----------------------------------------------------------------------
(   15)    !!   dta_tem      : read ocean temperature data
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce             ! ocean dynamics and tracers
(   19)    USE dom_oce         ! ocean space and time domain
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE daymod          ! calendar
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC dta_tem   ! called by step.F90 and inidta.F90
(   28) 
(   29)    !! * Shared module variables
(   30)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .TRUE.   !: temperature data flag
(   31)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   32)       t_dta             !: temperature data at given time-step
(   33) 
(   34)    !! * Module variables
(   35)    CHARACTER (len=45) ::   &
(   36)       cl_tdata
(   37)    INTEGER ::   &
(   38)       nlecte =  0,   &  ! switch for the first read
(   39)       ntem1      ,   &  ! first record used
(   40)       ntem2             ! second record used
(   41)    REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
(   42)       temdta            ! temperature data at two consecutive times
(   43) 
(   44)    !! * Substitutions
(   45) # 1 "./domzgr_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 7

(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 46 "dtatem.F90"
(   46) # 46
(   46)    !!----------------------------------------------------------------------
(   47)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   48)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtatem.F90,v 1.12 2006/04/19 14:43:15 opalod Exp $ 
(   49)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)    !!----------------------------------------------------------------------
(   55)    !!   Default case                                            NetCDF file
(   56)    !!----------------------------------------------------------------------
(   57) 
(   58)    SUBROUTINE dta_tem( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                   ***  ROUTINE dta_tem  ***
(   61)       !!                    
(   62)       !! ** Purpose :   Reads monthly temperature data 
(   63)       !! 
(   64)       !! ** Method  :   Read on unit numtdt the interpolated temperature 
(   65)       !!      onto the model grid.
(   66)       !!      Data begin at january. 
(   67)       !!      The value is centered at the middle of month.
(   68)       !!      In the opa model, kt=1 agree with january 1.
(   69)       !!      At each time step, a linear interpolation is applied between 
(   70)       !!      two monthly values.
(   71)       !!      Read on unit numtdt
(   72)       !!
(   73)       !! ** Action  :   define t_dta array at time-step kt
(   74)       !!
(   75)       !! History :
(   76)       !!        !  91-03  ()  Original code
(   77)       !!        !  92-07  (M. Imbard)
(   78)       !!        !  99-10  (M.A. Foujols, M. Imbard)  NetCDF FORMAT 
(   79)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   80)       !!----------------------------------------------------------------------
(   81)       !! * Modules used
(   82) !!DB
(   83)       USE lib_ncdf
(   84) 
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER, PARAMETER ::   &
(   90)          jpmois = 12                    ! number of month
(   91)       INTEGER ::   ji, jj, jl           ! dummy loop indicies
(   92)       INTEGER ::   &
(   93)          imois, iman, itime, ik ,    &  ! temporary integers
(   94)          i15, ipi, ipj, ipk             !    "          "
(   95) # 96
(   96)       INTEGER ::   &
(   97)          il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
(   98) 
(   99) # 100






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 8

(  100)       INTEGER, DIMENSION(jpmois) ::   istep
(  101)       REAL(wp) ::   zxy, zl, zdate0
(  102)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
(  103)       REAL(wp), DIMENSION(jpk) ::   zlev
(  104) !!DB
(  105)       INTEGER :: len, status 
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108)       ! 0. Initialization
(  109)       ! -----------------
(  110) 
(  111)       iman  = jpmois
(  112)       i15   = nday / 16
(  113)       imois = nmonth + i15 - 1
(  114)       IF( imois == 0 )   imois = iman
(  115) 
(  116)       itime = jpmois
(  117)       ipi = jpiglo
(  118)       ipj = jpjglo
(  119)       ipk = jpk
(  120) 
(  121)       ! 1. First call kt=nit000
(  122)       ! -----------------------
(  123) 
(  124)       IF( kt == nit000 .AND. nlecte == 0 ) THEN
(  125)          ntem1 = 0
(  126)          cl_tdata = 'data_1m_potential_temperature_nomask.nc'
(  127)          IF(lwp) WRITE(numout,*)
(  128)          IF(lwp) WRITE(numout,*) ' dtatem : Monthly Temperature fields'
(  129)          IF(lwp) WRITE(numout,*) ' ~~~~~~'
(  130)          IF(lwp) WRITE(numout,*) '             NetCDF File'
(  131)          IF(lwp) WRITE(numout,*)cl_tdata
(  132)          IF(lwp) WRITE(numout,*)
(  133)          
(  134) # 139
(  139) !!DB: 
(  140)          call ncdf_get_dim_size(cl_tdata, 'time_counter', itime, status)
(  141)          call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
(  142)          call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
(  143)          call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
(  144)          if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
(  145)             if(lwp) then
(  146)                write(numout,*)
(  147)                write(numout,*) 'problem with dimensions'
(  148)                write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  149)                write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  150)                write(numout,*) ' ipk ',ipk,' jpk ',jpk
(  151)             endif
(  152)             stop 'dtatem'
(  153)          endif
(  154) !!DB: 
(  155)          if( itime /= jpmois ) then
(  156)             if(lwp) then
(  157)                write(numout,*)
(  158)                write(numout,*) 'problem with time coordinates'
(  159)                write(numout,*) ' itime ',itime,' jpmois ',jpmois
(  160)             endif
(  161)             stop 'dta_sal'






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 9

(  162)          endif
(  163) 
(  164)       ENDIF
(  165) 
(  166) 
(  167)       ! 2. Read monthly file
(  168)       ! -------------------
(  169) 
(  170)       IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
(  171)          nlecte = 1
(  172) 
(  173)          ! Calendar computation
(  174)          
(  175)          ntem1 = imois        ! first file record used 
(  176)          ntem2 = ntem1 + 1    ! last  file record used
(  177)          ntem1 = MOD( ntem1, iman )
(  178)          IF( ntem1 == 0 )   ntem1 = iman
(  179)          ntem2 = MOD( ntem2, iman )
(  180)          IF( ntem2 == 0 )   ntem2 = iman
(  181)          if(lwp) write(numout,*) 'dtatem reading records ',ntem1, ntem2
(  182)          
(  183)          ! Read monthly temperature data
(  184) !!DB
(  185)          call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,1),-ntem1,status)
(  186)          call ncdf_read(cl_tdata,'votemper',temdta(:,:,:,2),-ntem2,status)
(  187) 
(  188) 
(  189)          IF(lwp) WRITE(numout,*)
(  190)          IF(lwp) WRITE(numout,*) ' Monthly temperature input OK'
(  191)          IF(lwp) WRITE(numout,*)
(  192)          
(  193) # 194
(  194) !!DB -- orca-related code deleted
(  195) 
(  196) # 197
(  197)          !                                  ! Mask
(  198)          DO jl = 1, 2
(  199)             temdta(:,:,:,jl) = temdta(:,:,:,jl) * tmask(:,:,:)
(  200)             temdta(:,:,jpk,jl) = 0.
(  201)             IF( lk_zps ) THEN                ! z-coord. with partial steps
(  202)                DO jj = 1, jpj                  ! interpolation of temperature at the last level
(  203)                   DO ji = 1, jpi
(  204)                      ik = mbathy(ji,jj) - 1
(  205)                      IF( ik > 2 ) THEN
(  206)                         zl = ( gdept(ik) - gdept(ik) ) / ( gdept(ik) - gdept(ik-1) )
(  207)                         temdta(ji,jj,ik,jl) = (1.-zl) * temdta(ji,jj,ik,jl) + zl * temdta(ji,jj,ik-1,jl) 
(  208)                      ENDIF
(  209)                   END DO
(  210)                END DO
(  211)             ENDIF
(  212)          END DO
(  213) 
(  214) !         IF(lwp) THEN
(  215) !            WRITE(numout,*) ' temperature month ', ntem1, ntem2
(  216) !            WRITE(numout,*)
(  217) !            WRITE(numout,*) '  month = ', ntem1, '  level = 1'
(  218) !            CALL prihre( temdta(:,:,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  219) !            WRITE(numout,*) '  month = ', ntem1, '  level = ', jpk/2






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 10

(  220) !            CALL prihre( temdta(:,:,jpk/2,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  221) !            WRITE(numout,*) '  month = ',ntem1,'  level = ', jpkm1
(  222) !            CALL prihre( temdta(:,:,jpkm1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  223) !         ENDIF
(  224)       ENDIF
(  225) 
(  226)  
(  227)       ! 2. At every time step compute temperature data
(  228)       ! ----------------------------------------------
(  229)       zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  230) !!DB: Note that for monthly data it is not necessary to do this every kt.
(  231) !!    So in the future use a code fragment like: 
(  232) !      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
(  233) !         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
(  234) !      endif
(  235) 
(  236)       t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
(  237) 
(  238) 
(  239)    END SUBROUTINE dta_tem
(  240) 
(  241) # 252
(  252)    !!======================================================================
(  253) END MODULE dtatem
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dtatem.F90: 229)
