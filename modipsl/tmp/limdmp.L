


PGF90 (Version     10.5)          09/23/2010  09:50:01      page 1

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

Filename: limdmp.F90

(    1) MODULE limdmp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limdmp   ***
(    4)    !!  Ice model : restoring Ice thickness and Fraction leads
(    5)    !!======================================================================
(    6) !byoung #if defined key_ice_lim
(    7) #if defined key_ice_lim && defined key_tradmp
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   10)    !!----------------------------------------------------------------------
(   11)    !!   lim_dmp      : ice model damping
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE ice
(   16)    USE ice_oce
(   17)    USE tradmp
(   18)    USE dom_oce
(   19)    USE oce
(   20)    USE daymod          ! calendar
(   21)    
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC lim_dmp     ! called by ice_step
(   27)    
(   28)    !! * Shared module variables
(   29)    CHARACTER (len=38) ::   &
(   30)       cl_icedata = 'ice_damping_ATL4.nc'
(   31)    INTEGER ::   &
(   32)         nice1      ,   &  ! first record used
(   33)         nice2             ! second record used
(   34)    
(   35)     REAL(wp), DIMENSION(jpi,jpj,2) ::   &
(   36)          hicif_data ,   & ! ice thickness data at two consecutive times
(   37)          frld_data        ! fraction lead data at two consecutive times
(   38) 
(   39)     REAL(wp), DIMENSION(jpi,jpj) ::   &
(   40)          hicif_dta ,   &  ! ice thickness at a given time
(   41)          frld_dta         ! fraction lead at a given time
(   42) 
(   43)    !! * Substitution
(   44) #  include "vectopt_loop_substitute.h90"
(   45)    !!----------------------------------------------------------------------
(   46)    !!   LIM 2.0 , UCL-LOCEAN-IPSL  (2005)






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 2

(   47)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdmp.F90,v 1.1 2006/03/21 08:42:21 opalod Exp $
(   48)    !! This software is governed by the CeCILL licence see !modipsl/doc/NEMO_CeCILL.txt
(   49)    !!----------------------------------------------------------------------
(   50) 
(   51) CONTAINS
(   52) 
(   53)    SUBROUTINE lim_dmp(kt)
(   54)       !!-------------------------------------------------------------------
(   55)       !!                   ***  ROUTINE lim_dmp ***
(   56)       !!
(   57)       !! ** purpose : ice model damping : restoring ice thickness and 
(   58)       !!              fraction leads
(   59)       !!
(   60)       !! ** method  : the key_tradmp must be used to compute resto(:,:) coef.
(   61)       !!     
(   62)       !! ** action :
(   63)       !!
(   64)       !! History :
(   65)       !!
(   66)       !!   2.0  !  04-04 (S. Theetten) Original
(   67)       !!---------------------------------------------------------------------
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(   70) 
(   71)       !! * Local Variables
(   72)       INTEGER  ::   ji, jj         ! dummy loop indices
(   73)       !!---------------------------------------------------------------------
(   74)     
(   75)       CALL dta_lim(kt)
(   76) 
(   77)       DO jj = 2, jpjm1
(   78)          DO ji = fs_2, fs_jpim1   ! vector opt.
(   79) 
(   80)             hicif(ji,jj) = hicif(ji,jj) - rdt_ice * resto(ji,jj,1) * ( hicif(ji,jj) -  hicif_dta(ji,jj))
(   81)             frld(ji,jj)  = frld(ji,jj)  - rdt_ice * resto(ji,jj,1) * ( frld(ji,jj)  - frld_dta(ji,jj))  
(   82) 
(   83)          ENDDO
(   84)       ENDDO
(   85) 
(   86)    END SUBROUTINE lim_dmp
(   87) 
(   88) 
(   89) 
(   90)    SUBROUTINE dta_lim( kt ) 
(   91)       !!----------------------------------------------------------------------
(   92)       !!                   ***  ROUTINE dta_lim  ***
(   93)       !!
(   94)       !! ** Purpose :   Reads monthly ice thickness and fraction lead  data
(   95)       !!
(   96)       !! ** Method  :   Read on unit numicedt the interpolated ice variable
(   97)       !!      onto the model grid.
(   98)       !!      Data begin at january.
(   99)       !!      The value is centered at the middle of month.
(  100)       !!      In the opa model, kt=1 agree with january 1.
(  101)       !!      At each time step, a linear interpolation is applied between
(  102)       !!      two monthly values.
(  103)       !!      
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 3

(  105)       !! ** Action  :   define hicif_dta and frld_dta arrays at time-step kt
(  106)       !!
(  107)       !! History :
(  108)       !!   2.0   !   04-04 (S. Theetten) Original
(  109)       !!----------------------------------------------------------------------
(  110)       !! * Modules used
(  111)       USE ioipsl
(  112) 
(  113)       !! * Arguments
(  114)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(  115) 
(  116)       !! * Local declarations
(  117)       INTEGER, PARAMETER ::   jpmois = 12       ! number of month
(  118)      
(  119)       INTEGER ::   &
(  120)          imois, iman, itime ,    &  ! temporary integers
(  121)          i15, ipi, ipj, ipk         !    "          "
(  122) 
(  123)       INTEGER, DIMENSION(jpmois) ::   istep
(  124)       REAL(wp) ::   zxy, zdate0, zdt
(  125)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
(  126)       REAL(wp), DIMENSION(jpk) ::   zlev
(  127)       !!----------------------------------------------------------------------
(  128) 
(  129)       ! 0. Initialization
(  130)       ! -----------------
(  131)       iman  = jpmois
(  132)       i15   = nday / 16
(  133)       imois = nmonth + i15 - 1
(  134)       IF( imois == 0 )   imois = iman
(  135) 
(  136)       itime = jpmois
(  137)       ipi=jpiglo
(  138)       ipj=jpjglo
(  139)       ipk=1
(  140)       zdt=rdt
(  141) 
(  142)       ! 1. First call kt=nit000
(  143)       ! -----------------------
(  144) 
(  145)       IF( kt == nit000 ) THEN
(  146)          nice1 = 0
(  147)          IF(lwp) WRITE(numout,*)
(  148)          IF(lwp) WRITE(numout,*) 'dtalim : Ice thickness and lead fraction  monthly fields'
(  149)          IF(lwp) WRITE(numout,*) '~~~~~~'
(  150)          IF(lwp) WRITE(numout,*) '             NetCDF FORMAT'
(  151)          IF(lwp) WRITE(numout,*)
(  152)          
(  153)          ! open file
(  154)          
(  155)          CALL flinopen( TRIM(cl_icedata), mig(1), nlci , mjg(1),  nlcj, .FALSE.,  &
(  156)             &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, numice_dmp )
(  157) 
(  158)           ! title, dimensions and tests
(  159)          IF( itime /= jpmois ) THEN
(  160)             IF(lwp) THEN
(  161)                WRITE(numout,*)
(  162)                WRITE(numout,*) 'problem with time coordinates'






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 4

(  163)                WRITE(numout,*) ' itime ',itime,' jpmois ',jpmois
(  164)             ENDIF
(  165)             STOP 'dta_lim'
(  166)          ENDIF
(  167)          IF( ipi /= jpidta .OR. ipj /= jpjdta ) THEN
(  168)             IF(lwp) THEN
(  169)                WRITE(numout,*)
(  170)                WRITE(numout,*) 'problem with dimensions'
(  171)                WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  172)                WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  173)             ENDIF
(  174)             STOP 'dta_lim'
(  175)          ENDIF
(  176)          IF(lwp) WRITE(numout,*) itime,istep,zdate0,zdt,numice_dmp
(  177) 
(  178)       ENDIF
(  179) 
(  180) 
(  181)       ! 2. Read monthly file
(  182)       ! -------------------
(  183) 
(  184)       IF( ( kt == nit000 ) .OR. imois /= nice1 ) THEN
(  185) 
(  186)          ! Calendar computation
(  187)          
(  188)          nice1 = imois        ! first file record used 
(  189)          nice2 = nice1 + 1    ! last  file record used
(  190)          nice1 = MOD( nice1, iman )
(  191)          IF( nice1 == 0 )   nice1 = iman
(  192)          nice2 = MOD( nice2, iman )
(  193)          IF( nice2 == 0 )   nice2 = iman
(  194)          IF(lwp) WRITE(numout,*) 'first record file used nice1 ', nice1
(  195)          IF(lwp) WRITE(numout,*) 'last  record file used nice2 ', nice2
(  196)          
(  197)          ! Read monthly ice thickness Levitus 
(  198)          
(  199)          CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
(  200)             &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,1) )
(  201)          CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
(  202)             &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,2) )
(  203)          
(  204)          IF(lwp) WRITE(numout,*)
(  205)          IF(lwp) WRITE(numout,*) ' read ice thickness ok'
(  206)          IF(lwp) WRITE(numout,*)
(  207) 
(  208)          ! Read monthly ice thickness Levitus 
(  209)          
(  210)          CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
(  211)             &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,1) )
(  212)          CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
(  213)             &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,2) )
(  214)          
(  215)          ! The fraction lead read in the file is in fact the 
(  216)          ! ice concentration which is 1 - the fraction lead
(  217)          frld_data = 1 - frld_data          
(  218)          
(  219)          IF(lwp) WRITE(numout,*)
(  220)          IF(lwp) WRITE(numout,*) ' read fraction lead ok'






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 5

(  221)          IF(lwp) WRITE(numout,*)
(  222) 
(  223) 
(  224)          IF(lwp) THEN
(  225)             WRITE(numout,*) ' Ice thickness month ', nice1,' and ', nice2
(  226)             WRITE(numout,*)
(  227)             WRITE(numout,*) ' Ice thickness month = ', nice1
(  228)             CALL prihre( hicif_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  229)             WRITE(numout,*)
(  230)             WRITE(numout,*) ' Fraction lead months ', nice1,' and ', nice2
(  231)             WRITE(numout,*)
(  232)             WRITE(numout,*) ' Fraction lead month = ', nice1
(  233)             CALL prihre( frld_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  234)          ENDIF
(  235)          
(  236)          ! 2. At every time step compute ice thickness and fraction lead data
(  237)          ! ------------------------------------------------------------------
(  238)          
(  239)          zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  240)          hicif_dta(:,:) = (1.-zxy) * hicif_data(:,:,1) + zxy * hicif_data(:,:,2)
(  241)          frld_dta(:,:) = (1.-zxy) * frld_data(:,:,1) + zxy * frld_data(:,:,2)
(  242) 
(  243)       ENDIF
(  244) 
(  245) 
(  246)    END SUBROUTINE dta_lim
(  247) 
(  248) #else
(  249)    !!----------------------------------------------------------------------
(  250)    !!   Default option         Empty Module                  No ice damping
(  251)    !!----------------------------------------------------------------------
(  252) CONTAINS
(  253)    SUBROUTINE lim_dmp(kt)        ! Dummy routine
(  254)       WRITE(*,*) 'lim_dmp: You should not see this print! error? ', kt
(  255)    END SUBROUTINE lim_dmp
(  256) #endif
(  257) 
(  258)    !!======================================================================
(  259) 
(  260) END MODULE limdmp
























PGF90 (Version     10.5)          09/23/2010  09:50:01      page 6

(    1) # 1 "limdmp.F90"
(    1) MODULE limdmp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limdmp   ***
(    4)    !!  Ice model : restoring Ice thickness and Fraction leads
(    5)    !!======================================================================
(    6) !byoung #if defined 1
(    7) # 8
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   10)    !!----------------------------------------------------------------------
(   11)    !!   lim_dmp      : ice model damping
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE ice
(   16)    USE ice_oce
(   17)    USE tradmp
(   18)    USE dom_oce
(   19)    USE oce
(   20)    USE daymod          ! calendar
(   21)    
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC lim_dmp     ! called by ice_step
(   27)    
(   28)    !! * Shared module variables
(   29)    CHARACTER (len=38) ::   &
(   30)       cl_icedata = 'ice_damping_ATL4.nc'
(   31)    INTEGER ::   &
(   32)         nice1      ,   &  ! first record used
(   33)         nice2             ! second record used
(   34)    
(   35)     REAL(wp), DIMENSION(jpi,jpj,2) ::   &
(   36)          hicif_data ,   & ! ice thickness data at two consecutive times
(   37)          frld_data        ! fraction lead data at two consecutive times
(   38) 
(   39)     REAL(wp), DIMENSION(jpi,jpj) ::   &
(   40)          hicif_dta ,   &  ! ice thickness at a given time
(   41)          frld_dta         ! fraction lead at a given time
(   42) 
(   43)    !! * Substitution
(   44) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 45 "limdmp.F90"
(   45) # 45






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 7

(   45)    !!----------------------------------------------------------------------
(   46)    !!   LIM 2.0 , UCL-LOCEAN-IPSL  (2005)
(   47)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdmp.F90,v 1.1 2006/03/21 08:42:21 opalod Exp $
(   48)    !! This software is governed by the CeCILL licence see !modipsl/doc/NEMO_CeCILL.txt
(   49)    !!----------------------------------------------------------------------
(   50) 
(   51) CONTAINS
(   52) 
(   53)    SUBROUTINE lim_dmp(kt)
(   54)       !!-------------------------------------------------------------------
(   55)       !!                   ***  ROUTINE lim_dmp ***
(   56)       !!
(   57)       !! ** purpose : ice model damping : restoring ice thickness and 
(   58)       !!              fraction leads
(   59)       !!
(   60)       !! ** method  : the 1 must be used to compute resto(:,:) coef.
(   61)       !!     
(   62)       !! ** action :
(   63)       !!
(   64)       !! History :
(   65)       !!
(   66)       !!   2.0  !  04-04 (S. Theetten) Original
(   67)       !!---------------------------------------------------------------------
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(   70) 
(   71)       !! * Local Variables
(   72)       INTEGER  ::   ji, jj         ! dummy loop indices
(   73)       !!---------------------------------------------------------------------
(   74)     
(   75)       CALL dta_lim(kt)
(   76) 
(   77)       DO jj = 2, jpjm1
(   78)          DO ji = 2, jpim1   ! vector opt.
(   79) 
(   80)             hicif(ji,jj) = hicif(ji,jj) - rdt_ice * resto(ji,jj,1) * ( hicif(ji,jj) -  hicif_dta(ji,jj))
(   81)             frld(ji,jj)  = frld(ji,jj)  - rdt_ice * resto(ji,jj,1) * ( frld(ji,jj)  - frld_dta(ji,jj))  
(   82) 
(   83)          ENDDO
(   84)       ENDDO
(   85) 
(   86)    END SUBROUTINE lim_dmp
(   87) 
(   88) 
(   89) 
(   90)    SUBROUTINE dta_lim( kt ) 
(   91)       !!----------------------------------------------------------------------
(   92)       !!                   ***  ROUTINE dta_lim  ***
(   93)       !!
(   94)       !! ** Purpose :   Reads monthly ice thickness and fraction lead  data
(   95)       !!
(   96)       !! ** Method  :   Read on unit numicedt the interpolated ice variable
(   97)       !!      onto the model grid.
(   98)       !!      Data begin at january.
(   99)       !!      The value is centered at the middle of month.
(  100)       !!      In the opa model, kt=1 agree with january 1.
(  101)       !!      At each time step, a linear interpolation is applied between
(  102)       !!      two monthly values.






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 8

(  103)       !!      
(  104)       !!
(  105)       !! ** Action  :   define hicif_dta and frld_dta arrays at time-step kt
(  106)       !!
(  107)       !! History :
(  108)       !!   2.0   !   04-04 (S. Theetten) Original
(  109)       !!----------------------------------------------------------------------
(  110)       !! * Modules used
(  111)       USE ioipsl
(  112) 
(  113)       !! * Arguments
(  114)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step
(  115) 
(  116)       !! * Local declarations
(  117)       INTEGER, PARAMETER ::   jpmois = 12       ! number of month
(  118)      
(  119)       INTEGER ::   &
(  120)          imois, iman, itime ,    &  ! temporary integers
(  121)          i15, ipi, ipj, ipk         !    "          "
(  122) 
(  123)       INTEGER, DIMENSION(jpmois) ::   istep
(  124)       REAL(wp) ::   zxy, zdate0, zdt
(  125)       REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
(  126)       REAL(wp), DIMENSION(jpk) ::   zlev
(  127)       !!----------------------------------------------------------------------
(  128) 
(  129)       ! 0. Initialization
(  130)       ! -----------------
(  131)       iman  = jpmois
(  132)       i15   = nday / 16
(  133)       imois = nmonth + i15 - 1
(  134)       IF( imois == 0 )   imois = iman
(  135) 
(  136)       itime = jpmois
(  137)       ipi=jpiglo
(  138)       ipj=jpjglo
(  139)       ipk=1
(  140)       zdt=rdt
(  141) 
(  142)       ! 1. First call kt=nit000
(  143)       ! -----------------------
(  144) 
(  145)       IF( kt == nit000 ) THEN
(  146)          nice1 = 0
(  147)          IF(lwp) WRITE(numout,*)
(  148)          IF(lwp) WRITE(numout,*) 'dtalim : Ice thickness and lead fraction  monthly fields'
(  149)          IF(lwp) WRITE(numout,*) '~~~~~~'
(  150)          IF(lwp) WRITE(numout,*) '             NetCDF FORMAT'
(  151)          IF(lwp) WRITE(numout,*)
(  152)          
(  153)          ! open file
(  154)          
(  155)          CALL flinopen( TRIM(cl_icedata), mig(1), nlci , mjg(1),  nlcj, .FALSE.,  &
(  156)             &           ipi, ipj, ipk, zlon, zlat, zlev, itime, istep, zdate0, zdt, numice_dmp )
(  157) 
(  158)           ! title, dimensions and tests
(  159)          IF( itime /= jpmois ) THEN
(  160)             IF(lwp) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 9

(  161)                WRITE(numout,*)
(  162)                WRITE(numout,*) 'problem with time coordinates'
(  163)                WRITE(numout,*) ' itime ',itime,' jpmois ',jpmois
(  164)             ENDIF
(  165)             STOP 'dta_lim'
(  166)          ENDIF
(  167)          IF( ipi /= jpidta .OR. ipj /= jpjdta ) THEN
(  168)             IF(lwp) THEN
(  169)                WRITE(numout,*)
(  170)                WRITE(numout,*) 'problem with dimensions'
(  171)                WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
(  172)                WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  173)             ENDIF
(  174)             STOP 'dta_lim'
(  175)          ENDIF
(  176)          IF(lwp) WRITE(numout,*) itime,istep,zdate0,zdt,numice_dmp
(  177) 
(  178)       ENDIF
(  179) 
(  180) 
(  181)       ! 2. Read monthly file
(  182)       ! -------------------
(  183) 
(  184)       IF( ( kt == nit000 ) .OR. imois /= nice1 ) THEN
(  185) 
(  186)          ! Calendar computation
(  187)          
(  188)          nice1 = imois        ! first file record used 
(  189)          nice2 = nice1 + 1    ! last  file record used
(  190)          nice1 = MOD( nice1, iman )
(  191)          IF( nice1 == 0 )   nice1 = iman
(  192)          nice2 = MOD( nice2, iman )
(  193)          IF( nice2 == 0 )   nice2 = iman
(  194)          IF(lwp) WRITE(numout,*) 'first record file used nice1 ', nice1
(  195)          IF(lwp) WRITE(numout,*) 'last  record file used nice2 ', nice2
(  196)          
(  197)          ! Read monthly ice thickness Levitus 
(  198)          
(  199)          CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
(  200)             &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,1) )
(  201)          CALL flinget( numice_dmp, 'icethic', jpidta, jpjdta, jpk,  &
(  202)             &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, hicif_data(1:nlci,1:nlcj,2) )
(  203)          
(  204)          IF(lwp) WRITE(numout,*)
(  205)          IF(lwp) WRITE(numout,*) ' read ice thickness ok'
(  206)          IF(lwp) WRITE(numout,*)
(  207) 
(  208)          ! Read monthly ice thickness Levitus 
(  209)          
(  210)          CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
(  211)             &          jpmois, nice1, nice1, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,1) )
(  212)          CALL flinget( numice_dmp, 'ileadfra', jpidta, jpjdta, jpk,  &
(  213)             &          jpmois, nice2, nice2, mig(1), nlci, mjg(1), nlcj, frld_data(1:nlci,1:nlcj,2) )
(  214)          
(  215)          ! The fraction lead read in the file is in fact the 
(  216)          ! ice concentration which is 1 - the fraction lead
(  217)          frld_data = 1 - frld_data          
(  218)          






PGF90 (Version     10.5)          09/23/2010  09:50:01      page 10

(  219)          IF(lwp) WRITE(numout,*)
(  220)          IF(lwp) WRITE(numout,*) ' read fraction lead ok'
(  221)          IF(lwp) WRITE(numout,*)
(  222) 
(  223) 
(  224)          IF(lwp) THEN
(  225)             WRITE(numout,*) ' Ice thickness month ', nice1,' and ', nice2
(  226)             WRITE(numout,*)
(  227)             WRITE(numout,*) ' Ice thickness month = ', nice1
(  228)             CALL prihre( hicif_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  229)             WRITE(numout,*)
(  230)             WRITE(numout,*) ' Fraction lead months ', nice1,' and ', nice2
(  231)             WRITE(numout,*)
(  232)             WRITE(numout,*) ' Fraction lead month = ', nice1
(  233)             CALL prihre( frld_data(1,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
(  234)          ENDIF
(  235)          
(  236)          ! 2. At every time step compute ice thickness and fraction lead data
(  237)          ! ------------------------------------------------------------------
(  238)          
(  239)          zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  240)          hicif_dta(:,:) = (1.-zxy) * hicif_data(:,:,1) + zxy * hicif_data(:,:,2)
(  241)          frld_dta(:,:) = (1.-zxy) * frld_data(:,:,1) + zxy * frld_data(:,:,2)
(  242) 
(  243)       ENDIF
(  244) 
(  245) 
(  246)    END SUBROUTINE dta_lim
(  247) 
(  248) 
(  249) # 258
(  258)    !!======================================================================
(  259) 
(  260) END MODULE limdmp
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (limdmp.F90: 239)
