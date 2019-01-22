


PGF90 (Version     10.5)          09/23/2010  09:49:06      page 1

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

Filename: dtasal.F90

(    1) !!DB: 2008.10.01 -- Routine modified 
(    2) !!-1- deleted reference to orca
(    3) !!-2- replaced IOIPSL with lib_ncdf 
(    4) 
(    5) MODULE dtasal
(    6)    !!======================================================================
(    7)    !!                     ***  MODULE  dtasal  ***
(    8)    !! Ocean data  :  read ocean salinity data from monthly atlas data
(    9)    !!=====================================================================
(   10) #if defined key_dtasal   ||   defined key_esopa
(   11)    !!----------------------------------------------------------------------
(   12)    !!   'key_dtasal'                                          salinity data
(   13)    !!----------------------------------------------------------------------
(   14)    !!   dta_sal      : read ocean salinity data
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE daymod          ! calendar
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dta_sal   ! called by step.F90 and inidta.F90
(   27)    
(   28)    !! * Shared module variables
(   29)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .TRUE.    !: salinity data flag
(   30)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   31)       s_dta       !: salinity data at given time-step
(   32) 
(   33)    !! * Module variables
(   34)    CHARACTER (len=32) ::   clname
(   35)    INTEGER ::   &
(   36)       nlecsa = 0,   &  ! switch for the first read
(   37)       nsal1     ,   &  ! first record used
(   38)       nsal2            ! second record used
(   39)    REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
(   40)       saldta    ! salinity data at two consecutive times
(   41) 
(   42)    !! * Substitutions
(   43) #  include "domzgr_substitute.h90"
(   44)    !!----------------------------------------------------------------------
(   45)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtasal.F90,v 1.12 2006/04/19 14:43:14 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 2

(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    !!----------------------------------------------------------------------
(   53)    !!   Default option:                                         NetCDF file
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56)    SUBROUTINE dta_sal( kt )
(   57)       !!----------------------------------------------------------------------
(   58)       !!                   ***  ROUTINE dta_sal  ***
(   59)       !!        
(   60)       !! ** Purpose :   Reads monthly salinity data
(   61)       !!             
(   62)       !! ** Method  : - Read on unit numsdt the monthly salinity data interpo-
(   63)       !!     lated onto the model grid.
(   64)       !!              - At each time step, a linear interpolation is applied
(   65)       !!     between two monthly values.
(   66)       !!
(   67)       !! History :
(   68)       !!        !  91-03  ()  Original code
(   69)       !!        !  92-07  (M. Imbard)
(   70)       !!   9.0  !  02-06  (G. Madec)  F90: Free form and module 
(   71)       !!----------------------------------------------------------------------
(   72)       !! * Modules used
(   73) !!DB
(   74)       USE lib_ncdf
(   75) 
(   76)       !! * Arguments
(   77)       INTEGER, INTENT(in) ::   kt             ! ocean time step
(   78) 
(   79)       !! * Local declarations
(   80) 
(   81) 
(   82)       INTEGER, PARAMETER ::   jpmois = 12, jpf = 1
(   83)       INTEGER ::   ji, jj, jl           ! dummy loop indicies
(   84)       INTEGER ::   &
(   85)          imois, iman, ik, i15,       &  ! temporary integers
(   86)          ipi, ipj, ipk, itime           !    "          "
(   87) #if defined key_tradmp
(   88)       INTEGER ::   &
(   89)          jk, il0, il1,               &  ! temporary integers
(   90)          ii0, ii1, ij0, ij1             !    "          "
(   91) #endif
(   92) 
(   93)       INTEGER, DIMENSION(jpmois) ::   istep
(   94)       REAL(wp) ::   &
(   95)          zxy, zl, zdate0
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
(   97)       REAL(wp), DIMENSION(jpk) ::   zlev
(   98) !!DB
(   99)       INTEGER :: len, status 
(  100) 
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103)       ! 0. Initialization
(  104)       ! -----------------






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 3

(  105) 
(  106)       iman  = jpmois
(  107)       i15   = nday / 16
(  108) 
(  109)       imois = nmonth + i15 - 1
(  110)       IF( imois == 0 ) imois = iman
(  111) 
(  112)       itime = jpmois
(  113)       ipi=jpiglo
(  114)       ipj=jpjglo
(  115)       ipk = jpk
(  116) 
(  117)       ! 1. First call kt=nit000
(  118)       ! -----------------------
(  119) 
(  120)       IF( kt == nit000 .AND. nlecsa == 0 ) THEN
(  121)          nsal1 = 0
(  122)          IF(lwp) THEN
(  123)             WRITE(numout,*)
(  124)             WRITE(numout,*) ' dta_sal : monthly salinity data in NetCDF file'
(  125)             WRITE(numout,*) ' ~~~~~~~'
(  126)             WRITE(numout,*)
(  127)          ENDIF
(  128) 
(  129)          ! open file
(  130)          clname = 'data_1m_salinity_nomask.nc'
(  131) #if defined key_agrif
(  132)          if ( .NOT. Agrif_Root() ) then
(  133)             clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  134)          endif
(  135) #endif          
(  136) !!DB: 
(  137)          call ncdf_get_dim_size(clname, 'time_counter', itime, status)
(  138)          call ncdf_get_dim_size(clname, 'x', ipi, status)
(  139)          call ncdf_get_dim_size(clname, 'y', ipj, status)
(  140)          call ncdf_get_dim_size(clname, 'z', ipk, status)
(  141)          if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
(  142)             if(lwp) then
(  143)                write(numout,*)
(  144)                write(numout,*) 'problem with dimensions'
(  145)                write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  146)                write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  147)                write(numout,*) ' ipk ',ipk,' jpk ',jpk
(  148)             endif
(  149)             stop 'dtasal'
(  150)          endif
(  151) !!DB: 
(  152)          if( itime /= jpmois ) then
(  153)             if(lwp) then
(  154)                write(numout,*)
(  155)                write(numout,*) 'problem with time coordinates'
(  156)                write(numout,*) ' itime ',itime,' jpmois ',jpmois
(  157)             endif
(  158)             stop 'dta_sal'
(  159)          endif
(  160) 
(  161)       ENDIF
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 4

(  163) 
(  164)       ! 2. Read monthly file
(  165)       ! -------------------
(  166) 
(  167)       IF( ( kt == nit000 .AND. nlecsa == 0) .OR. imois /= nsal1 ) THEN
(  168)          nlecsa = 1
(  169)          
(  170)          ! 2.1 Calendar computation
(  171)          
(  172)          nsal1 = imois        ! first file record used 
(  173)          nsal2 = nsal1 + 1    ! last  file record used
(  174)          nsal1 = MOD( nsal1, iman )
(  175)          IF( nsal1 == 0 ) nsal1 = iman
(  176)          nsal2 = MOD( nsal2, iman )
(  177)          IF( nsal2 == 0 ) nsal2 = iman
(  178)          if(lwp) write(numout,*) 'dtasal reading records ',nsal1, nsal2
(  179)          
(  180)          ! 2.3 Read monthly salinity data
(  181) 
(  182) !!DB
(  183)          call ncdf_read(clname,'vosaline',saldta(:,:,:,1),-nsal1,status)
(  184)          call ncdf_read(clname,'vosaline',saldta(:,:,:,2),-nsal2,status)
(  185) 
(  186)          IF(lwp) WRITE(numout,*)
(  187)          IF(lwp) WRITE(numout,*) ' Monthly salinity input OK'
(  188)          IF(lwp) WRITE(numout,*)
(  189)          
(  190) #if defined key_tradmp
(  191) !!DB -- orca code deleted
(  192) #endif
(  193)          
(  194)          !                                     ! Mask
(  195)          DO jl = 1, 2
(  196)             saldta(:,:,:,jl) = saldta(:,:,:,jl)*tmask(:,:,:)
(  197)             saldta(:,:,jpk,jl) = 0.
(  198)             IF( lk_zps ) THEN          ! z-coord. partial steps
(  199)                DO jj = 1, jpj          ! interpolation of salinity at the last ocean level (i.e. the partial step)
(  200)                   DO ji = 1, jpi
(  201)                      ik = mbathy(ji,jj) - 1
(  202)                      IF( ik > 2 ) THEN
(  203)                         zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
(  204)                         saldta(ji,jj,ik,jl) = (1.-zl) * saldta(ji,jj,ik,jl) +zl * saldta(ji,jj,ik-1,jl) 
(  205)                      ENDIF
(  206)                   END DO
(  207)                END DO
(  208)             ENDIF
(  209)          END DO
(  210)          
(  211) 
(  212) !         IF(lwp) THEN
(  213) !            WRITE(numout,*)' salinity month ',nsal1,nsal2
(  214) !            WRITE(numout,*)
(  215) !            WRITE(numout,*) ' month = ',nsal1,'  level = 1'
(  216) !            CALL prihre(saldta(:,:,1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
(  217) !            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpk/2
(  218) !            CALL prihre(saldta(:,:,jpk/2,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
(  219) !            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpkm1
(  220) !            CALL prihre(saldta(:,:,jpkm1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 5

(  221) !         ENDIF
(  222)       ENDIF
(  223)       
(  224)  
(  225)       ! 3. At every time step compute salinity data
(  226)       ! -------------------------------------------
(  227)       zxy = FLOAT(nday + 15 - 30*i15)/30.
(  228) !!DB: Note that for monthly data it is not necessary to do this every kt.
(  229) !!    So in the future use a code fragment like: 
(  230) !      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
(  231) !         s_dta(:,:,:) = (1.-zxy) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)
(  232) !      endif
(  233) 
(  234)       s_dta(:,:,:) = ( 1.- zxy ) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)
(  235) 
(  236)    END SUBROUTINE dta_sal
(  237) 
(  238) #else
(  239)    !!----------------------------------------------------------------------
(  240)    !!   Default option:                                    NO salinity data
(  241)    !!----------------------------------------------------------------------
(  242)    USE in_out_manager
(  243)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .FALSE.   !: salinity data flag
(  244) CONTAINS
(  245)    SUBROUTINE dta_sal( kt )        ! Empty routine
(  246)       if(lwp) WRITE(numout,*) 'dta_sal: You should not have seen this print! error?', kt
(  247)    END SUBROUTINE dta_sal
(  248) #endif
(  249)    !!======================================================================
(  250) END MODULE dtasal


































PGF90 (Version     10.5)          09/23/2010  09:49:06      page 6

(    1) # 1 "dtasal.F90"
(    1) !!DB: 2008.10.01 -- Routine modified 
(    2) !!-1- deleted reference to orca
(    3) !!-2- replaced IOIPSL with lib_ncdf 
(    4) 
(    5) MODULE dtasal
(    6)    !!======================================================================
(    7)    !!                     ***  MODULE  dtasal  ***
(    8)    !! Ocean data  :  read ocean salinity data from monthly atlas data
(    9)    !!=====================================================================
(   10) # 11
(   11)    !!----------------------------------------------------------------------
(   12)    !!   'key_dtasal'                                          salinity data
(   13)    !!----------------------------------------------------------------------
(   14)    !!   dta_sal      : read ocean salinity data
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE daymod          ! calendar
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dta_sal   ! called by step.F90 and inidta.F90
(   27)    
(   28)    !! * Shared module variables
(   29)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasal = .TRUE.    !: salinity data flag
(   30)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   31)       s_dta       !: salinity data at given time-step
(   32) 
(   33)    !! * Module variables
(   34)    CHARACTER (len=32) ::   clname
(   35)    INTEGER ::   &
(   36)       nlecsa = 0,   &  ! switch for the first read
(   37)       nsal1     ,   &  ! first record used
(   38)       nsal2            ! second record used
(   39)    REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
(   40)       saldta    ! salinity data at two consecutive times
(   41) 
(   42)    !! * Substitutions
(   43) # 1 "./domzgr_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 7

(   48) # 44 "dtasal.F90"
(   44) # 44
(   44)    !!----------------------------------------------------------------------
(   45)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DTA/dtasal.F90,v 1.12 2006/04/19 14:43:14 opalod Exp $ 
(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   48)    !!----------------------------------------------------------------------
(   49) 
(   50) CONTAINS
(   51) 
(   52)    !!----------------------------------------------------------------------
(   53)    !!   Default option:                                         NetCDF file
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56)    SUBROUTINE dta_sal( kt )
(   57)       !!----------------------------------------------------------------------
(   58)       !!                   ***  ROUTINE dta_sal  ***
(   59)       !!        
(   60)       !! ** Purpose :   Reads monthly salinity data
(   61)       !!             
(   62)       !! ** Method  : - Read on unit numsdt the monthly salinity data interpo-
(   63)       !!     lated onto the model grid.
(   64)       !!              - At each time step, a linear interpolation is applied
(   65)       !!     between two monthly values.
(   66)       !!
(   67)       !! History :
(   68)       !!        !  91-03  ()  Original code
(   69)       !!        !  92-07  (M. Imbard)
(   70)       !!   9.0  !  02-06  (G. Madec)  F90: Free form and module 
(   71)       !!----------------------------------------------------------------------
(   72)       !! * Modules used
(   73) !!DB
(   74)       USE lib_ncdf
(   75) 
(   76)       !! * Arguments
(   77)       INTEGER, INTENT(in) ::   kt             ! ocean time step
(   78) 
(   79)       !! * Local declarations
(   80) 
(   81) 
(   82)       INTEGER, PARAMETER ::   jpmois = 12, jpf = 1
(   83)       INTEGER ::   ji, jj, jl           ! dummy loop indicies
(   84)       INTEGER ::   &
(   85)          imois, iman, ik, i15,       &  ! temporary integers
(   86)          ipi, ipj, ipk, itime           !    "          "
(   87) # 88
(   88)       INTEGER ::   &
(   89)          jk, il0, il1,               &  ! temporary integers
(   90)          ii0, ii1, ij0, ij1             !    "          "
(   91) 
(   92) # 93
(   93)       INTEGER, DIMENSION(jpmois) ::   istep
(   94)       REAL(wp) ::   &
(   95)          zxy, zl, zdate0
(   96)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon, zlat
(   97)       REAL(wp), DIMENSION(jpk) ::   zlev
(   98) !!DB
(   99)       INTEGER :: len, status 






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 8

(  100) 
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103)       ! 0. Initialization
(  104)       ! -----------------
(  105) 
(  106)       iman  = jpmois
(  107)       i15   = nday / 16
(  108) 
(  109)       imois = nmonth + i15 - 1
(  110)       IF( imois == 0 ) imois = iman
(  111) 
(  112)       itime = jpmois
(  113)       ipi=jpiglo
(  114)       ipj=jpjglo
(  115)       ipk = jpk
(  116) 
(  117)       ! 1. First call kt=nit000
(  118)       ! -----------------------
(  119) 
(  120)       IF( kt == nit000 .AND. nlecsa == 0 ) THEN
(  121)          nsal1 = 0
(  122)          IF(lwp) THEN
(  123)             WRITE(numout,*)
(  124)             WRITE(numout,*) ' dta_sal : monthly salinity data in NetCDF file'
(  125)             WRITE(numout,*) ' ~~~~~~~'
(  126)             WRITE(numout,*)
(  127)          ENDIF
(  128) 
(  129)          ! open file
(  130)          clname = 'data_1m_salinity_nomask.nc'
(  131) # 136
(  136) !!DB: 
(  137)          call ncdf_get_dim_size(clname, 'time_counter', itime, status)
(  138)          call ncdf_get_dim_size(clname, 'x', ipi, status)
(  139)          call ncdf_get_dim_size(clname, 'y', ipj, status)
(  140)          call ncdf_get_dim_size(clname, 'z', ipk, status)
(  141)          if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
(  142)             if(lwp) then
(  143)                write(numout,*)
(  144)                write(numout,*) 'problem with dimensions'
(  145)                write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  146)                write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  147)                write(numout,*) ' ipk ',ipk,' jpk ',jpk
(  148)             endif
(  149)             stop 'dtasal'
(  150)          endif
(  151) !!DB: 
(  152)          if( itime /= jpmois ) then
(  153)             if(lwp) then
(  154)                write(numout,*)
(  155)                write(numout,*) 'problem with time coordinates'
(  156)                write(numout,*) ' itime ',itime,' jpmois ',jpmois
(  157)             endif
(  158)             stop 'dta_sal'
(  159)          endif
(  160) 
(  161)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 9

(  162) 
(  163) 
(  164)       ! 2. Read monthly file
(  165)       ! -------------------
(  166) 
(  167)       IF( ( kt == nit000 .AND. nlecsa == 0) .OR. imois /= nsal1 ) THEN
(  168)          nlecsa = 1
(  169)          
(  170)          ! 2.1 Calendar computation
(  171)          
(  172)          nsal1 = imois        ! first file record used 
(  173)          nsal2 = nsal1 + 1    ! last  file record used
(  174)          nsal1 = MOD( nsal1, iman )
(  175)          IF( nsal1 == 0 ) nsal1 = iman
(  176)          nsal2 = MOD( nsal2, iman )
(  177)          IF( nsal2 == 0 ) nsal2 = iman
(  178)          if(lwp) write(numout,*) 'dtasal reading records ',nsal1, nsal2
(  179)          
(  180)          ! 2.3 Read monthly salinity data
(  181) 
(  182) !!DB
(  183)          call ncdf_read(clname,'vosaline',saldta(:,:,:,1),-nsal1,status)
(  184)          call ncdf_read(clname,'vosaline',saldta(:,:,:,2),-nsal2,status)
(  185) 
(  186)          IF(lwp) WRITE(numout,*)
(  187)          IF(lwp) WRITE(numout,*) ' Monthly salinity input OK'
(  188)          IF(lwp) WRITE(numout,*)
(  189)          
(  190) # 191
(  191) !!DB -- orca code deleted
(  192) # 193
(  193)          
(  194)          !                                     ! Mask
(  195)          DO jl = 1, 2
(  196)             saldta(:,:,:,jl) = saldta(:,:,:,jl)*tmask(:,:,:)
(  197)             saldta(:,:,jpk,jl) = 0.
(  198)             IF( lk_zps ) THEN          ! z-coord. partial steps
(  199)                DO jj = 1, jpj          ! interpolation of salinity at the last ocean level (i.e. the partial step)
(  200)                   DO ji = 1, jpi
(  201)                      ik = mbathy(ji,jj) - 1
(  202)                      IF( ik > 2 ) THEN
(  203)                         zl = ( gdept(ik) - gdept(ik) ) / ( gdept(ik) - gdept(ik-1) )
(  204)                         saldta(ji,jj,ik,jl) = (1.-zl) * saldta(ji,jj,ik,jl) +zl * saldta(ji,jj,ik-1,jl) 
(  205)                      ENDIF
(  206)                   END DO
(  207)                END DO
(  208)             ENDIF
(  209)          END DO
(  210)          
(  211) 
(  212) !         IF(lwp) THEN
(  213) !            WRITE(numout,*)' salinity month ',nsal1,nsal2
(  214) !            WRITE(numout,*)
(  215) !            WRITE(numout,*) ' month = ',nsal1,'  level = 1'
(  216) !            CALL prihre(saldta(:,:,1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
(  217) !            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpk/2
(  218) !            CALL prihre(saldta(:,:,jpk/2,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
(  219) !            WRITE(numout,*) ' month = ',nsal1,'  level = ',jpkm1






PGF90 (Version     10.5)          09/23/2010  09:49:06      page 10

(  220) !            CALL prihre(saldta(:,:,jpkm1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
(  221) !         ENDIF
(  222)       ENDIF
(  223)       
(  224)  
(  225)       ! 3. At every time step compute salinity data
(  226)       ! -------------------------------------------
(  227)       zxy = FLOAT(nday + 15 - 30*i15)/30.
(  228) !!DB: Note that for monthly data it is not necessary to do this every kt.
(  229) !!    So in the future use a code fragment like: 
(  230) !      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
(  231) !         s_dta(:,:,:) = (1.-zxy) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)
(  232) !      endif
(  233) 
(  234)       s_dta(:,:,:) = ( 1.- zxy ) * saldta(:,:,:,1) + zxy * saldta(:,:,:,2)
(  235) 
(  236)    END SUBROUTINE dta_sal
(  237) 
(  238) # 249
(  249)    !!======================================================================
(  250) END MODULE dtasal
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dtasal.F90: 227)
