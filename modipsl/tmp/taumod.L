


PGF90 (Version     10.5)          09/23/2010  09:49:12      page 1

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

Filename: taumod.F90

(    1) !!DB: 2009.08.31 -- Eliminated GYRE config
(    2) MODULE taumod
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  taumod  ***
(    5)    !! Ocean forcing : stress at the the ocean surface
(    6)    !!=====================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   tau          : define the surface stress for the ocean
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE phycst          ! physical constants
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE daymod          ! calendar
(   16)    USE lbclnk          ! 
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC tau                ! routine called by step.F90
(   23) 
(   24)    !! * Share modules variables
(   25) 
(   26) !!DB: NEW 2009.05.05 -- required for core forcing and (eventually) updated lim-model
(   27) ! FD correction from JM Molines
(   28)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
(   29)       taux, tauy,  &      !: surface stress components in (i,j) referential
(   30)       tauxwo, tauywo,  &  !: surface stress components in (i,j) referential air-ocean
(   31)       tauxwi, tauywi,  &  !: surface stress components in (i,j) referential air-ice
(   32)       tauxg, tauyg        !: surface stress components in geographical
(   33)       !                   !  referential (used in output)
(   34) 
(   35) !!DB OLD
(   36) !   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
(   37) !      taux, tauy,      &  !: surface stress components in (i,j) referential
(   38) !      tauxg, tauyg        !: surface stress components in geographical
(   39)       !                   !  referential (used in output)
(   40) 
(   41)    !!----------------------------------------------------------------------
(   42)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   43)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/taumod.F90,v 1.7 2006/04/10 15:46:11 opalod Exp $ 
(   44)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   45)    !!----------------------------------------------------------------------
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 2

(   47) #if defined key_tau_monthly
(   48)    ! Monthly climatology in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
(   49)    !!----------------------------------------------------------------------
(   50)    !!   'key_tau_monthly'                        MONTHLY climatology stress
(   51)    !!   default case                                   NetCDF files
(   52)    !!----------------------------------------------------------------------
(   53) #   include "tau_forced_monthly.h90"
(   54) 
(   55) # elif defined key_tau_daily
(   56)    !!----------------------------------------------------------------------
(   57)    !!   'key_tau_daily'                                 DAILY stress
(   58)    !!                                                   NetCDF files
(   59)    !!----------------------------------------------------------------------
(   60)    ! Daily climatology/interannual in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
(   61) !#   include "tau_forced_daily.h90"
(   62) #   include "tau_forcing.h90"
(   63) 
(   64) # elif defined key_tau_cmc
(   65)    !!----------------------------------------------------------------------
(   66)    !!   'key_tau_cmc'                                   CMC WS stress
(   67)    !!                                                  default: 3 hourly for 6 days
(   68)    !!                                                   NetCDF files
(   69)    !!----------------------------------------------------------------------
(   70)    ! WS in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
(   71) !#   include "tau_forced_cmc.h90"
(   72) #   include "tau_forcing.h90"
(   73) #elif defined key_coupled
(   74)    ! Coupled case : stress at the coupling frequency
(   75) # if defined key_ice_lim
(   76)    !!----------------------------------------------------------------------
(   77)    !!   'key_coupled'                              Coupled Ocean/Atmosphere
(   78)    !!   'key_ice_lim'                                   LIM sea-ice
(   79)    !!----------------------------------------------------------------------
(   80)    ! New way: 3D referential link to the earth (avoid north pole pb)
(   81)    ! (3 component stress defined at U- and V-points)
(   82) #  include "tau_coupled_ice.h90"
(   83) # else
(   84)    !!----------------------------------------------------------------------
(   85)    !!   'key_coupled'                              Coupled Ocean/Atmosphere
(   86)    !!   Default case                                  NO sea-ice
(   87)    !!----------------------------------------------------------------------
(   88)    ! old fashion: geographical referential
(   89)    ! (zonal and meridional stress defined at U- and V-points)
(   90) #  include "tau_coupled.h90"
(   91) # endif
(   92) #else
(   93)    !!----------------------------------------------------------------------
(   94)    !!   Default option                                     constant forcing
(   95)    !!----------------------------------------------------------------------
(   96)    !! * local modules variables
(   97)    INTEGER  ::       & !!! * Namelist numtau *
(   98)       ntau000 = 1       ! nb of time-step during which the surface stress
(   99)       !                 ! increase from 0 to its nominal value (taudta) (>0)
(  100)    REAL(wp) ::       & !!! * Namelist numtau *
(  101)       tau0x = 0.e0 , &  ! constant wind stress value in i-direction
(  102)       tau0y = 0.e0      ! constant wind stress value in j-direction
(  103)    !!----------------------------------------------------------------------
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 3

(  105) CONTAINS
(  106) 
(  107)    SUBROUTINE tau( kt )
(  108)       !!---------------------------------------------------------------------
(  109)       !!                    ***  ROUTINE tau  ***
(  110)       !! 
(  111)       !! ** Purpose :   provide the ocean surface stress at each time step
(  112)       !!
(  113)       !! ** Method  :   Constant surface stress increasing from 0 to taudta 
(  114)       !!      value during the first ntau000 time-step (namelist)
(  115)       !!        CAUTION: never mask the surface stress field !
(  116)       !!
(  117)       !! ** Action  : - update taux , tauy the stress in (i,j) ref.
(  118)       !!              - update tauxg, tauyg the stress in geographic ref.
(  119)       !!
(  120)       !! History :
(  121)       !!   4.0  !  91-03  (G. Madec)  Original code
(  122)       !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(  123)       !!----------------------------------------------------------------------
(  124)       !! * Arguments
(  125)       INTEGER, INTENT( in  ) ::   kt    ! ocean time step
(  126)       REAL(wp) ::   ztau, ztau_sais, &  ! wind intensity and of the seasonal cycle
(  127)          ztime,                      &  ! time in hour
(  128)          ztimemax, ztimemin,         &  ! 21th June, and 21th decem. if date0 = 1st january
(  129)          ztaun                          ! intensity 
(  130)       INTEGER  ::   ji, jj              ! dummy loop indices
(  131) 
(  132)       INTEGER  ::           &
(  133)          zyear0,            &           ! initial year
(  134)          zmonth0,           &           ! initial month
(  135)          zday0,             &           ! initial day
(  136)          zday_year0                    ! initial day since january 1st
(  137)         
(  138) 
(  139)       !! * Local declarations
(  140)       REAL(wp) ::   zfacto              ! 
(  141) 
(  142)       NAMELIST/namtau/ ntau000, tau0x, tau0y
(  143)       !!---------------------------------------------------------------------
(  144) 
(  145)          IF( kt == nit000 ) THEN
(  146)    
(  147)             ! Read Namelist namtau : surface wind stress
(  148)             ! --------------------
(  149)             REWIND ( numnam )
(  150)             READ   ( numnam, namtau )
(  151)    
(  152)             IF(lwp) WRITE(numout,*)' '
(  153)             IF(lwp) WRITE(numout,*)' tau     : Constant surface wind stress read in namelist'
(  154)             IF(lwp) WRITE(numout,*)' ~~~~~~~ '
(  155)             IF(lwp) WRITE(numout,*)'           Namelist namtau: set the constant stress values'
(  156)             IF(lwp) WRITE(numout,*)'              spin up of the stress  ntau000 = ', ntau000, ' time-steps'
(  157)             IF(lwp) WRITE(numout,*)'              constant i-stress      tau0x   = ', tau0x  , ' N/m2'
(  158)             IF(lwp) WRITE(numout,*)'              constant j-stress      tau0y   = ', tau0y  , ' N/m2'
(  159)    
(  160)             ntau000 = MAX( ntau000, 1 )   ! must be >= 1
(  161)    
(  162)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 4

(  163)    
(  164)          ! Increase the surface stress to its nominal value in ntau000 time-step
(  165)          
(  166)          IF( kt <= ntau000 ) THEN
(  167)             zfacto = 0.5 * (  1. - COS( rpi * FLOAT( kt ) / FLOAT( ntau000 ) )  )
(  168)             taux (:,:) = zfacto * tau0x
(  169)             tauy (:,:) = zfacto * tau0y
(  170)             tauxg(:,:) = zfacto * tau0x
(  171)             tauyg(:,:) = zfacto * tau0y
(  172)          ENDIF
(  173) 
(  174)       
(  175)    END SUBROUTINE tau
(  176) #endif
(  177)    !!======================================================================
(  178) END MODULE taumod
















































PGF90 (Version     10.5)          09/23/2010  09:49:12      page 5

(    1) # 1 "taumod.F90"
(    1) !!DB: 2009.08.31 -- Eliminated GYRE config
(    2) MODULE taumod
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  taumod  ***
(    5)    !! Ocean forcing : stress at the the ocean surface
(    6)    !!=====================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   tau          : define the surface stress for the ocean
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE phycst          ! physical constants
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE daymod          ! calendar
(   16)    USE lbclnk          ! 
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC tau                ! routine called by step.F90
(   23) 
(   24)    !! * Share modules variables
(   25) 
(   26) !!DB: NEW 2009.05.05 -- required for core forcing and (eventually) updated lim-model
(   27) ! FD correction from JM Molines
(   28)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
(   29)       taux, tauy,  &      !: surface stress components in (i,j) referential
(   30)       tauxwo, tauywo,  &  !: surface stress components in (i,j) referential air-ocean
(   31)       tauxwi, tauywi,  &  !: surface stress components in (i,j) referential air-ice
(   32)       tauxg, tauyg        !: surface stress components in geographical
(   33)       !                   !  referential (used in output)
(   34) 
(   35) !!DB OLD
(   36) !   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
(   37) !      taux, tauy,      &  !: surface stress components in (i,j) referential
(   38) !      tauxg, tauyg        !: surface stress components in geographical
(   39)       !                   !  referential (used in output)
(   40) 
(   41)    !!----------------------------------------------------------------------
(   42)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   43)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/taumod.F90,v 1.7 2006/04/10 15:46:11 opalod Exp $ 
(   44)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   45)    !!----------------------------------------------------------------------
(   46) 
(   47) # 56
(   56)    !!----------------------------------------------------------------------
(   57)    !!   'key_tau_daily'                                 DAILY stress
(   58)    !!                                                   NetCDF files
(   59)    !!----------------------------------------------------------------------
(   60)    ! Daily climatology/interannual in (i,j) referential  (i-comp. at U-pt and j-comp. at V-pt)
(   61) !#   include "tau_forced_daily.h90"
(   62) # 1 "./tau_forcing.h90"
(    1) !!DB: 2009.05.07
(    2) !!General routine to input wind stress forcing and perform
(    3) !!time-interpolation






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 6

(    4) !!Routine accounts for the case where the time interval between input
(    5) !!records is NOT an integer # of timesteps
(    6) !!NB: I/O routines use lib_ncdf
(    7) !!NB: idy variable retained, but routine does not take account of
(    8) !!leap year (it is not clear that it needs to)
(    9) 
(   10) !!NB: Routine clamps maximum input record to max-record-in-file if
(   11) !!is_cyclical is FALSE
(   12) !!NB: 2009.05.07  Routine  modified if normal-year daily forcing 
(   13) !!is used to account for the cyclical nature of the forcing.
(   14) !!NB: if is_cyclical = .true. must have is_daily=.true. as well
(   15) !!NB: Routine does not check if WS files have consistent dimensions wrt model domain
(   16) 
(   17)    !!----------------------------------------------------------------------
(   18)    !!                   ***  tau_forcing.h90  ***
(   19)    !!----------------------------------------------------------------------
(   20) 
(   21)    !!----------------------------------------------------------------------
(   22)    !!   tau     :   update the surface stress - daily fields in NetCDF file
(   23)    !!----------------------------------------------------------------------
(   24)    !! * local modules variables
(   25)    INTEGER ::   &
(   26)       ntau1, ntau2 ,  &  ! index of the first and second record used
(   27)       ndt_in,ndt_interp,  &  ! related to updating fields
(   28)       numrecs_tau, iread_tau
(   29) 
(   30)    CHARACTER (len=34) ::   &      !!! 
(   31)       cl_taux ,  & ! generic name of the i-component monthly NetCDF file
(   32)       cl_tauy      ! generic name of the j-component monthly NetCDF file
(   33) 
(   34) !!DB
(   35)    REAL(wp), DIMENSION(jpi,jpj,2) ::   &
(   36)       taux_dta,    &  ! i- and j-components of the surface stress (Pascal)
(   37)       tauy_dta        ! at 2 consecutive frames in (i,j) referential
(   38) 
(   39)    !!----------------------------------------------------------------------
(   40)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE tau( kt )
(   46)       !!---------------------------------------------------------------------
(   47)       !!                    ***  ROUTINE tau  ***
(   48)       !!     
(   49)       !! DB: 2009.02
(   50)       !! ** Purpose :   provide to the ocean the stress at each time step
(   51)       !!
(   52)       !! ** Method  :   Read the surface stress components from NetCDF 
(   53)       !!      file in  (i,j) coordinates. 
(   54)       !!        The i-component is given at U-point (INTERP package)
(   55)       !!        The j-component is given at V-point (INTERP package)
(   56)       !!
(   57)       !!    CAUTION: never mask the surface stress field !
(   58)       !!
(   59)       !! ** Action  :  update every ndt_interp time-steps the two components of the 
(   60)       !!      surface stress in (i,j) coords
(   61)       !!      NB:  geographical coords not used






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 7

(   62)       !!
(   63)       !!----------------------------------------------------------------------
(   64)       !! * Modules used
(   65) 
(   66) !!DB
(   67)       USE oce, only :  perpetual_forcing, ramp
(   68)       USE lib_ncdf
(   69) 
(   70)       !! * Arguments
(   71)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER :: idy
(   75) !!DB
(   76)       REAL(wp) :: zxy, rdt_in, d_rec, d_interp
(   77)       INTEGER :: i, j, status, rec1, rec2
(   78)       INTEGER, SAVE :: frame0
(   79)       LOGICAL :: is_daily, is_cmc, is_monthly, is_cyclical
(   80) 
(   81) 
(   82)       !!---------------------------------------------------------------------
(   83) 
(   84) 
(   85) 
(   86)       ! -------------- !
(   87)       ! Initialization !
(   88)       ! -------------- !
(   89) 
(   90) !!DB: 
(   91)       is_cyclical = .true. !!NB: hardwired
(   92) # 93
(   93)       is_daily = .true.  
(   94)       is_cmc = .false. 
(   95) 
(   96) # 102
(  102) !!time between input file records in hours; can be determined from file itself
(  103)       if(is_daily) then
(  104)          d_rec = 24.0
(  105)          d_interp = 6.0       ! time interpolation in hours = default
(  106)       endif
(  107)       if(is_cmc) then
(  108)          d_rec = 3.0
(  109)          d_interp = 1.0       ! time interpolation in hours = default
(  110)       endif
(  111)       rdt_in = d_rec*3600./rdt       !# of timesteps between input records (float)
(  112)       ndt_in = nint(d_rec*3600./rdt) !# of timesteps between input records (int; not used)
(  113)       ndt_interp = nint(d_interp*3600./rdt) !# of timesteps between time interpolations (int)
(  114) 
(  115) 
(  116) !!DB: OLD code 
(  117)       idy   = 365
(  118)       IF ( nleapy == 1 ) idy = 366 
(  119) 
(  120) !!DBG
(  121)       if(kt==nit000 .AND. lwp) write(numout2,*)'DBG: generic tau_forcing kt, ramp = ', kt, ramp
(  122) 
(  123) 
(  124)       ! -------------------- !






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 8

(  125)       ! First call kt=nit000 !
(  126)       ! -------------------- !
(  127) 
(  128)       if( kt == nit000 ) then
(  129) !!DB: set frame0 == offset into the files. 
(  130) !!NB: frame0 when is_daily=true is a guess based on assuming a cyclical year of forcing
(  131)          if(is_daily) then
(  132)             if(is_cyclical) then
(  133)                frame0 = nday_year     
(  134)             else
(  135)                frame0 = 1               ! for non-cyclical daily forcing
(  136)             endif
(  137)             cl_taux = 'taux_1d.nc'
(  138)             cl_tauy = 'tauy_1d.nc'
(  139)          endif
(  140)          if(is_cmc) then
(  141)             !!DB In this case can be determined from namelist (ndate0?)
(  142)             frame0 = 1  !default
(  143)             cl_taux = 'taux_3h.nc'
(  144)             cl_tauy = 'tauy_3h.nc'
(  145)          endif
(  146)          
(  147)          if(lwp) then
(  148)             write(numout2,*) ' '
(  149)             write(numout2,*) ' tau    : generic wind stress in NetCDF files'
(  150)             write(numout2,*) ' ~~~~~~~'
(  151)             write(numout2,*) ' rdt_in, d_rec, ndt_interp = ', rdt_in, d_rec, ndt_interp 
(  152)             write(numout2,*) ' HARDWIRED: is_cyclical = ', is_cyclical
(  153)             write(numout2,*) ' frame0 = ', frame0
(  154)             write(numout2,*) ' taux & tauy files: ',cl_taux, cl_tauy
(  155)          endif
(  156)          ! title, dimensions and tests
(  157) # 162
(  162)          
(  163) 
(  164) 
(  165) # 170
(  170) !!DB -- crucial assignments
(  171)          call ncdf_get_dim_size(cl_taux, 'time_counter',numrecs_tau,status)
(  172)          if(status == 0) then
(  173)             if(lwp) write(numout2,*)'DBG: # records in ',cl_taux, ' = ', numrecs_tau
(  174)          else
(  175)             if(lwp) write(numout2,*)'DBG: prob with ncdf_get_dim_size ====> STOP'
(  176)             call flush(numout2)
(  177)             stop
(  178)          endif
(  179) 
(  180)          if(frame0 .gt. numrecs_tau) frame0 = numrecs_tau  !to be safe, but should not occur
(  181)          ntau1 =  frame0 - 1 !!to get the first input correct (see below)
(  182)          if(ntau1 .gt. numrecs_tau) ntau1 = numrecs_tau
(  183)          ntau2 = ntau1 + 1
(  184)          if(ntau2 .gt. numrecs_tau) ntau2 = numrecs_tau
(  185)          iread_tau = 0 ! next input timestep rel nit000; set to 0 to start
(  186) 
(  187)       endif
(  188) 
(  189)       ! ----------------- !
(  190)       ! Read file   !






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 9

(  191)       ! ----------------- !
(  192) !!Account for cyclical nature of forcing
(  193) !!Note that ntau1 must be continually incremented if is_cyclical=.true. 
(  194) !!and rec1 & rec2 point to the appropriate record to input
(  195)       if(kt-nit000 == iread_tau) then
(  196) 
(  197)          ntau1 = ntau1 + 1
(  198)          if(is_cyclical) then    !roll over to next cycle
(  199)             rec1 = mod(ntau1,numrecs_tau)
(  200)             if(rec1 == 0) rec1 = numrecs_tau
(  201)          else                    !clamp if necessary
(  202)             if(ntau1 .gt. numrecs_tau) ntau1 = numrecs_tau
(  203)             rec1 = ntau1
(  204)          endif
(  205) 
(  206)          ntau2 = ntau1 + 1
(  207)          if(is_cyclical) then    !roll over to next cycle
(  208)             rec2 = mod(ntau2,numrecs_tau)
(  209)             if(rec2 == 0) rec2 = numrecs_tau
(  210)          else                    !clamp if necessary
(  211)             if(ntau2 .gt. numrecs_tau) ntau2 = numrecs_tau
(  212)             rec2 = ntau2
(  213)          endif
(  214)          iread_tau = nint((ntau1-(frame0-1))*rdt_in) !! next input timestep rel nit000
(  215) 
(  216) !!DB: perpetual_forcing from oce; 
(  217)          if(perpetual_forcing > 0) then
(  218)             ntau1 = frame0
(  219)             ntau2 = ntau1 
(  220)             rec1 = ntau1
(  221)             rec2 = rec1
(  222)             iread_tau = -1   !ensures (kt-nit000) /= iread_tau, thus only 1 input of WS files
(  223)             if(lwp .AND. kt==nit000) then 
(  224)                write(numout2,*)'DBG: In tau_forcing, perpetual_forcing = ', perpetual_forcing
(  225)                write(numout2,*)'DBG: In tau_forcing, ntau2 = ntau1 = ', frame0
(  226)             endif
(  227)          endif
(  228)             
(  229) ! Read wind stress data using lib_ncdf; NB: use -ve index to point directly to file frame-#
(  230)          call ncdf_read(cl_taux, 'sozotaux', taux_dta(:,:,1), -rec1, status)
(  231)          if(status /= 0) then
(  232)             if(lwp) then
(  233)                write(numout2,*)'DBG: Prob reading taux file, status = ', status
(  234)                call flush(numout2)
(  235)             endif
(  236)          endif
(  237) 
(  238)          call ncdf_read(cl_taux, 'sozotaux', taux_dta(:,:,2), -rec2, status)
(  239)          if(status /= 0) then
(  240)             if(lwp) then
(  241)                write(numout2,*)'DBG: Prob reading taux file, status = ', status
(  242)             endif
(  243)          endif
(  244) 
(  245)          call ncdf_read(cl_tauy, 'sometauy', tauy_dta(:,:,1), -rec1, status)
(  246)          if(status /= 0) then
(  247)             if(lwp) then
(  248)                write(numout2,*)'DBG: Prob reading tauy file'






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 10

(  249)             endif
(  250)          endif
(  251)          call ncdf_read(cl_tauy, 'sometauy', tauy_dta(:,:,2), -rec2, status)
(  252)          if(status /= 0) then
(  253)             if(lwp) then
(  254)                write(numout2,*)'DBG: Prob reading tauy file'
(  255)             endif
(  256)          endif
(  257) 
(  258) !!DBG:
(  259)          if(lwp) write(numout2,*)'DBG: Input WS fields ', rec1, rec2, ' kt = ', kt
(  260)          if(lwp)write(numout2,*)'DBG: Next WS input at iread_tau = ', iread_tau
(  261) 
(  262)          CALL FLUSH(numout2)
(  263) 
(  264) 
(  265)       ENDIF
(  266) 
(  267) !!DB: apply ramp here for perpetual_forcing case
(  268)       if(perpetual_forcing > 0) then
(  269)          taux(:,:) = ramp*taux_dta(:,:,1)
(  270)          tauy(:,:) = ramp*tauy_dta(:,:,1)
(  271) !!DBG
(  272)          if(lwp .AND. mod(kt-nit000,int(0.25*86400.0/rdt))==0) then
(  273)             write(3200+narea,'(a40,i7,2x,(3f10.6,1x))')'WS: perpetual forcing on -- kt, ramp ', & 
(  274)                  kt,ramp,taux_dta(20,20,1),taux(20,20)
(  275)          endif
(  276) 
(  277)       endif
(  278) 
(  279) 
(  280) 
(  281) !!DB Use rdt_in & ndt_interp in a formula that accounts for the nature of 
(  282) !!the ndt_* variables which may not be integers (depending on rdt)
(  283) !!and also gives a consistent zxy upon restart
(  284) !!NB: perpetual_forcing case is handled above, and although the below should work
(  285) !!if perpetual_forcing /= 0 it is safer to avoid this
(  286)       if(mod(kt-nit000,ndt_interp) == 0  .AND.  perpetual_forcing == 0) then
(  287)          zxy = (float(kt-nit000) - (ntau1-1-(frame0-1))*rdt_in)/rdt_in
(  288)          zxy = min(zxy,1.0)  !to be safe
(  289)          taux(:,:) = (1.-zxy) * taux_dta(:,:,1) + zxy * taux_dta(:,:,2)
(  290)          tauy(:,:) = (1.-zxy) * tauy_dta(:,:,1) + zxy * tauy_dta(:,:,2)
(  291)          taux(:,:) = ramp*taux(:,:)
(  292)          tauy(:,:) = ramp*tauy(:,:)
(  293) !!DBG
(  294) !         if(lwp) write(3200+narea,'(a35,i7,2x,(5f10.6,1x))')'WS: time interp -- kt, ramp, zxy ', & 
(  295) !              kt,ramp,zxy,taux_dta(20,20,2),taux_dta(20,20,1),taux(20,20)
(  296)       endif
(  297) 
(  298) 
(  299)    END SUBROUTINE tau
(  300) # 63 "taumod.F90"
(   63) 
(   64) # 177
(  177)    !!======================================================================
(  178) END MODULE taumod
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (taumod.F90: 287)
