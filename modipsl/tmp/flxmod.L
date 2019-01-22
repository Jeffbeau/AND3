


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

Filename: flxmod.F90

(    1) MODULE flxmod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  flxmod  ***
(    4)    !! Ocean forcing:  thermohaline forcing of the ocean
(    5)    !!=====================================================================
(    6)    !! * Modules used
(    7)    USE oce             ! ocean dynamics and tracers
(    8)    USE dom_oce         ! ocean space and time domain
(    9)    USE phycst          ! physical constants
(   10)    USE in_out_manager  ! I/O manager
(   11)    USE lib_mpp         ! distribued memory computing library
(   12)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   13)    USE daymod          ! calendar
(   14)    USE ocfzpt          ! ocean freezing point
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC flx       ! routine called by step.F90
(   21)    PUBLIC flx_init  ! routine called by opa.F90
(   22) 
(   23)    !! * local declarations
(   24)    REAL(wp), PUBLIC ::            & !!! surface fluxes namelist (namflx)
(   25)       q0    = 0.e0,               &  ! net heat flux
(   26)       qsr0  = 0.e0,               &  ! solar heat flux
(   27)       emp0  = 0.e0,               &  ! net freshwater flux
(   28)       dqdt0 = -40.,               &  ! coefficient for SST damping (W/m2/K)
(   29)       deds0 = 27.7                   ! coefficient for SSS damping (mm/day)
(   30)    
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxmod.F90,v 1.4 2006/04/19 14:43:16 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) #if defined key_flx_bulk_monthly
(   40)    !!----------------------------------------------------------------------
(   41)    !!   'key_flx_bulk_monthly'   and                           MONTHLY bulk
(   42)    !!   Default option                                         Net CDF file
(   43)    !!----------------------------------------------------------------------
(   44) #  include "flx_bulk_monthly.h90"
(   45) 
(   46) #elif defined key_flx_bulk_daily






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 2

(   47)    !!----------------------------------------------------------------------
(   48)    !!   'key_flx_bulk_daily'                                     DAILY bulk
(   49)    !!                                                          Net CDF file
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) !!DB 2008.06.27:
(   53) # if defined key_flx_bulk_cmc
(   54)    !!----------------------------------------------------------------------
(   55)    !!   'key_flx_bulk_cmc'                                     CMC bulk
(   56)    !!                                                          Net CDF file
(   57)    !!----------------------------------------------------------------------
(   58) #  include "flx_bulk_cmc.h90"        
(   59) !!DB 2009.08.06
(   60) #elif defined key_CORE_NY || defined key_CORE_ANNUAL
(   61) !! CORE Normal year forcing or 1958 ... 2006 forcing
(   62) !! DB modified flx_core.h90, using DB prepared SOPA12 files
(   63) #  include "db_flx_core.h90"
(   64) # else
(   65) #  include "flx_bulk_daily.h90"
(   66) # endif
(   67) 
(   68) #elif defined key_flx_forced_daily
(   69)    !!----------------------------------------------------------------------
(   70)    !!   'key_flx_forced_daily'                                 DAILY fluxes
(   71)    !!                                                          Net CDF file
(   72)    !!----------------------------------------------------------------------
(   73) #  include "flx_forced_daily.h90"
(   74) 
(   75) #elif defined key_flx_forced_monthly 
(   76) !byoung
(   77) #   include "flx_forced_monthly.h90"
(   78) 
(   79) #elif defined key_coupled
(   80) # if defined key_ice_lim
(   81)    !!----------------------------------------------------------------------
(   82)    !!   'key_coupled'  and                          Coupled Ocan/Atmosphere
(   83)    !!   'key_ice_lim'                               with  LIM sea-ice model
(   84)    !!----------------------------------------------------------------------
(   85) #  include "flx_coupled_ice.h90"
(   86) 
(   87) # else
(   88)    !!----------------------------------------------------------------------
(   89)    !!   'key_flx_coupled'  and                      Coupled Ocan/Atmosphere
(   90)    !!   Default option                              without   sea-ice model
(   91)    !!----------------------------------------------------------------------
(   92) #  include "flx_coupled_noice.h90"
(   93) 
(   94) # endif
(   95) #else
(   96)    !!----------------------------------------------------------------------
(   97)    !!   Default option                                   Analytical forcing
(   98)    !!----------------------------------------------------------------------
(   99)    !!   flx          : define the thermohaline fluxes for the ocean
(  100)    !!----------------------------------------------------------------------
(  101)    !! * Module used
(  102)    !!----------------------------------------------------------------------
(  103) 
(  104) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 3

(  105) 
(  106)    SUBROUTINE flx ( kt )
(  107)       !!---------------------------------------------------------------------
(  108)       !!                    ***  ROUTINE flx  ***
(  109)       !!              
(  110)       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(  111)       !!      to the ocean at each time step.
(  112)       !!
(  113)       !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
(  114)       !!
(  115)       !! ** Action  : - q, qt, qsr, emp, emps, qrp, erp
(  116)       !!
(  117)       !! History :
(  118)       !!        !  91-03  ()  Original code
(  119)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  120)       !!----------------------------------------------------------------------
(  121)       !! * arguments
(  122)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  123)       !!---------------------------------------------------------------------
(  124) 
(  125)       IF( kt == nit000 )  THEN
(  126)          IF(lwp) WRITE(numout,*)' '
(  127)          IF(lwp) WRITE(numout,*)'flx     : Analytical/Constant surface fluxes'
(  128)          IF(lwp) WRITE(numout,*)'~~~~~~~ '
(  129)          IF(lwp) WRITE(numout,*)'          See the routine oce_sbc'     
(  130)          IF(lwp) WRITE(numout,*)' '
(  131)       ENDIF
(  132) 
(  133)    END SUBROUTINE flx
(  134) 
(  135) #endif
(  136) 
(  137) 
(  138)    SUBROUTINE flx_init
(  139)       !!---------------------------------------------------------------------
(  140)       !!                    ***  ROUTINE flx  ***
(  141)       !!              
(  142)       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(  143)       !!      to the ocean at each time step.
(  144)       !!
(  145)       !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
(  146)       !!
(  147)       !! ** Action  : - q, qt, qsr, emp, emps, qrp, erp
(  148)       !!
(  149)       !! History :
(  150)       !!        !  91-03  ()  Original code
(  151)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  152)       !!----------------------------------------------------------------------
(  153)       NAMELIST/namflx/ q0, qsr0, emp0, dqdt0, deds0
(  154)       !!---------------------------------------------------------------------
(  155) 
(  156)       ! Read Namelist namflx : surface thermohaline fluxes
(  157)       ! --------------------
(  158)       REWIND ( numnam )
(  159)       READ   ( numnam, namflx )
(  160) 
(  161)       IF(lwp) THEN
(  162)          WRITE(numout,*)' '






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 4

(  163)          WRITE(numout,*)'flx_init : thermohaline forcing '
(  164)          WRITE(numout,*)'~~~~~~~~ '
(  165)          WRITE(numout,*)'           net heat flux                  q0   = ', q0  , ' W/m2'
(  166)          WRITE(numout,*)'           solar heat flux                qsr0 = ', qsr0, ' W/m2'
(  167)          WRITE(numout,*)'           net heat flux                  emp0 = ', emp0, ' W/m2'
(  168)          WRITE(numout,*)'           coefficient for SST damping   dqdt0 = ', dqdt0,' W/m2/K'
(  169)          WRITE(numout,*)'           coefficient for SSS damping   deds0 = ', deds0,' mm/day'
(  170)       ENDIF
(  171) 
(  172)    END SUBROUTINE flx_init
(  173) 
(  174) END MODULE flxmod




















































PGF90 (Version     10.5)          09/23/2010  09:49:12      page 5

(    1) # 1 "flxmod.F90"
(    1) MODULE flxmod
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  flxmod  ***
(    4)    !! Ocean forcing:  thermohaline forcing of the ocean
(    5)    !!=====================================================================
(    6)    !! * Modules used
(    7)    USE oce             ! ocean dynamics and tracers
(    8)    USE dom_oce         ! ocean space and time domain
(    9)    USE phycst          ! physical constants
(   10)    USE in_out_manager  ! I/O manager
(   11)    USE lib_mpp         ! distribued memory computing library
(   12)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   13)    USE daymod          ! calendar
(   14)    USE ocfzpt          ! ocean freezing point
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC flx       ! routine called by step.F90
(   21)    PUBLIC flx_init  ! routine called by opa.F90
(   22) 
(   23)    !! * local declarations
(   24)    REAL(wp), PUBLIC ::            & !!! surface fluxes namelist (namflx)
(   25)       q0    = 0.e0,               &  ! net heat flux
(   26)       qsr0  = 0.e0,               &  ! solar heat flux
(   27)       emp0  = 0.e0,               &  ! net freshwater flux
(   28)       dqdt0 = -40.,               &  ! coefficient for SST damping (W/m2/K)
(   29)       deds0 = 27.7                   ! coefficient for SSS damping (mm/day)
(   30)    
(   31)    !! * Substitutions
(   32) # 1 "./domzgr_substitute.h90"
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
(   48) # 33 "flxmod.F90"
(   33) # 33
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxmod.F90,v 1.4 2006/04/19 14:43:16 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) # 47
(   47)    !!----------------------------------------------------------------------
(   48)    !!   'key_flx_bulk_daily'                                     DAILY bulk






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 6

(   49)    !!                                                          Net CDF file
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) !!DB 2008.06.27:
(   53) # 61
(   61) !! CORE Normal year forcing or 1958 ... 2006 forcing
(   62) !! DB modified flx_core.h90, using DB prepared SOPA12 files
(   63) # 1 "./db_flx_core.h90"
(    1) !!DB 2009.05.05
(    2) 
(    3)    !!----------------------------------------------------------------------
(    4)    !!                    ***  flx_core.h90  ***
(    5)    !!----------------------------------------------------------------------
(    6)    !!   flx     : update surface thermohaline fluxes from the NCAR data set
(    7)    !!             read in a NetCDF file
(    8)    !!----------------------------------------------------------------------
(    9)    !! * Modules used     C A U T I O N  already defined in flxmod.F90
(   10) 
(   11)    !! * Shared module variables
(   12)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &
(   13)       gsss,            &  !: SSS mean on nfbulk ocean time step
(   14)       gu  , gv ,       &  !: (un,vn)(jk=1) mean over nfbulk ocean time step
(   15)       !                   ! these variables are used for output in diawri
(   16)       qlw_oce  ,       &  !: long wave flx over ocean
(   17)       qla_oce  ,       &  !: latent heat flx over ocean
(   18)       qsb_oce  ,       &  !: sensible heat flx over ocean
(   19)       qlw_ice  ,       &  !: long wave flx over ice
(   20)       qsb_ice             !: sensible heat flx over ice
(   21) 
(   22)    !! * Module variables
(   23)    INTEGER, PARAMETER ::   jpfile  = 8  ! maximum number of files to read 
(   24)    CHARACTER (LEN=32), DIMENSION (jpfile) ::  clvarname  
(   25)    CHARACTER (LEN=50), DIMENSION (jpfile) ::  clname 
(   26) 
(   27)    INTEGER   ::      isnap
(   28)    INTEGER, DIMENSION(jpfile) ::   &
(   29)       numflxall,           &  ! logical units for surface fluxes data
(   30)       nrecflx  , nrecflx2     ! first and second record to be read in flux file
(   31)    REAL(wp), DIMENSION(jpfile) ::   &
(   32)       freqh                 ! Frequency for each forcing file (hours)
(   33)       !                     ! a negative value of -12 corresponds to monthly
(   34)    REAL(wp), DIMENSION(jpi,jpj,jpfile,2) ::   &
(   35)       flxdta                !:  set of NCAR 6hourly/daily/monthly fluxes
(   36) 
(   37) !!DB -- deleted this option
(   38) 
(   39) !!DB -- keep for now
(   40)    LOGICAL  :: &
(   41)         &   ln_2m   = .FALSE.     !: logical flag for height of air temp. and hum.
(   42) 
(   43)    REAL(wp) :: alpha_precip=1.    !: multiplication factor for precipitation
(   44) 
(   45) 
(   46) !!DB -- related to I/O
(   47)    INTEGER :: bflx1, bflx2, num_recsd, ireadd, frame_daily
(   48)    INTEGER :: nyear0
(   49) 
(   50)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 7

(   51)    !! History :
(   52)    !!   9.0  !  04-08  (U. Schweckendiek)  Original code
(   53)    !!        !  05-04  (L. Brodeau, A.M. Treguier)  severals modifications:
(   54)    !!        !         - new bulk routine for efficiency
(   55)    !!        !         - WINDS ARE NOW ASSUMED TO BE AT T POINTS in input files
(   56)    !!        !         - file names and file characteristics in namelist 
(   57)    !!        !         - Implement reading of 6-hourly fields   
(   58)    !!        !  06-02  (S. Masson, G. Madec)  IOM read + NEMO "style"
(   59)    !!        !  12-06  (L. Brodeau) handle both 2m and 10m input fields
(   60)    !!----------------------------------------------------------------------
(   61)    !!   OPA 9.0 , LOCEAN-IPSL (2006)
(   62)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flx_core.h90,v 1.4 2007/05/30 08:32:50 opalod Exp $
(   63)    !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
(   64)    !!----------------------------------------------------------------------
(   65)    !
(   66) CONTAINS
(   67)    !
(   68)    SUBROUTINE flx( kt )
(   69)       !!---------------------------------------------------------------------
(   70)       !!                     ***  ROUTINE flx  ***
(   71)       !!                   
(   72)       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(   73)       !!                and momentum fluxes (tau) 
(   74)       !!      to the ocean at each time step.
(   75)       !!
(   76)       !! ** Method  :   Read NCAR data in a NetCDF file
(   77)       !!               (file names and frequency of inputs specified in namelist)
(   78)       !!          precipitation:  total (rain+snow)
(   79)       !!          precipitation:  snow only 
(   80)       !!          u10,v10 -> scalar wind at 10m in m/s -  ON 'T' GRID POINTS!!!
(   81)       !!          solar radiation (short wave)  in W/m2
(   82)       !!          thermal radiation (long wave) in W/m2
(   83)       !!          specific humidity             in %
(   84)       !!          temperature at 10m            in degrees K
(   85)       !!
(   86)       !!  ** Action  :   
(   87)       !!          call flx_blk_albedo to compute ocean and ice albedo
(   88)       !!          Calculates forcing fluxes to input into ice model 
(   89)       !          or to be used directly for the ocean in ocesbc. 
(   90)       !!
(   91)       !!  ** Outputs 
(   92)       !!          COMMENTS TO BE ADDED -units to be verified!
(   93)       !!          taux    : zonal      wind stress on "u" points  (N/m2)
(   94)       !!          tauy    : meridional wind stress on "v" points  (N/m2)
(   95)       !!          qsr_oce : Solar flux over the ocean (W/m2)
(   96)       !!          qnsr_oce: longwave flux over the ocean (W/m2)
(   97)       !!          qsr_ice : solar flux over the ice (W/m2)
(   98)       !!          qnsr_ice: longwave flux over the ice (W/m2)
(   99)       !!          qla_ice : latent flux over sea-ice (W/m2)
(  100)       !!          dqns_ice: total heat fluxes sensitivity over ice (dQ/dT)
(  101)       !!          dqla_ice: latent heat flux sensitivity over ice (dQla/dT)
(  102)       !!          fr1_i0
(  103)       !!          fr2_i0
(  104)       !!          tprecip
(  105)       !!          sprecip
(  106)       !!          evap
(  107)       !!          
(  108)       !! caution : now, in the opa global model, the net upward water flux is






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 8

(  109)       !! -------   with mm/day unit,i.e. Kg/m2/s.
(  110)       !!---------------------------------------------------------------------
(  111)       !! * modules used
(  112) !!DB
(  113) !      USE iom
(  114) !      USE restart         
(  115) !! Will need (some of) these later
(  116) !!DB
(  117)       USE oce, only :  perpetual_forcing, ramp
(  118)       USE lib_ncdf
(  119) 
(  120) !
(  121)       USE par_oce
(  122)       USE flx_oce
(  123)       USE blk_oce         ! bulk variable
(  124)       USE taumod
(  125)       USE bulk
(  126)       USE phycst
(  127)       USE lbclnk
(  128)       USE dtatem          ! FOR Ce = F(SST(levitus)):
(  129)       USE daymod 
(  130) 
(  131) !!DB: -- in case no ice model
(  132) # 133
(  133) ! FD add for wind relative to ice velocity
(  134)       USE ice, only    : u_ice , v_ice ! : ice velocity
(  135) 
(  136) # 137
(  137)       !! * arguments
(  138)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  139) 
(  140) !!DB -- will change
(  141)       !! * Local declarations      
(  142) 
(  143) # 155
(  155)       CHARACTER(len=45)  ::  &
(  156)          clname_1 = 'core_t10_1d.nc',        &
(  157)          clname_2 = 'core_precip_1d.nc',   &
(  158)          clname_3 = 'core_q10_1d.nc',        &
(  159)          clname_4 = 'core_u10_v10_1d.nc',        &
(  160)          clname_5 = 'core_SWLW_1d.nc'
(  161) 
(  162)       LOGICAL ::  &
(  163)          is_cyclical = .true.
(  164) 
(  165) # 166
(  166)       REAL(wp) :: zxy, update_frac, rdt_in  !!NB: not used
(  167)       INTEGER :: i, j,  status, ipi, ipj, idy, ndt, rec1, rec2
(  168) 
(  169) 
(  170)       !! * Local declarations
(  171)       INTEGER ::   &
(  172)          jj   , ji                     !  Loop indices 
(  173)       REAL(wp) ::   &
(  174)          zcof , zzu , zzv              ! scalar 
(  175)       REAL(wp), DIMENSION(jpi,jpj, jpfile) ::   &
(  176)          flxnow             ! input flux values at current time step 
(  177)       REAL(wp), DIMENSION(jpi,jpj) ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 9

(  178)          dqlw_ice ,      &  ! long wave heat flx sensitivity over ice
(  179)          dqsb_ice ,      &  ! sensible heat flx sensitivity over ice
(  180)          alb_oce_os,     &  ! albedo of the ocean under overcast sky
(  181)          alb_ice_os,     &  ! albedo of the ice under overcast sky
(  182)          alb_ice_cs,     &  ! albedo of ice under clear sky
(  183)          alb_oce_cs,     &  ! albedo of the ocean under clear sky
(  184)          zsst,           &  ! nfbulk : mean SST
(  185)          zsss,           &  ! nfbulk : mean tn_ice(:,:,1)
(  186)          zut,            &  ! nfbulk : mean U at T-point
(  187)          zvt,            &  ! nfbulk : mean V at T-point
(  188)          dUnormt,        &  ! scalar wind (norm) on T points
(  189) ! FD modifications from JM Molines
(  190)          tauxtwo, tauytwo,   &  ! wind stress on ocean computed at T-point
(  191)          tauxtwi, tauytwi,   &  ! wind stress on ice computed at T-point
(  192) ! FD different Unormt for ice from Steven Alderson
(  193)          Unormt,         &  ! scalar wind (norm) on T points
(  194)          zus,            &  ! nfbulk : mean u_ice at T-point
(  195)          zvs,            &  ! nfbulk : mean v_ice at T-point
(  196) ! FD end 
(  197)          tauxt, tauyt,   &  ! wind stress computed at T-point
(  198)          qsatw,          &  ! specific humidity at zSST
(  199)          qsat,           &  ! specific humidity at zsss 
(  200)          Ch,             &  ! coefficient for sensible heat flux
(  201)          Ce,             &  ! coefficient for latent heat flux
(  202)          Cd,             &  ! coefficient for wind stress
(  203)          zt_zu,          &  !: air temperature at wind speed height
(  204)          zq_zu              !: air spec. hum.  at wind speed height
(  205)       REAL(wp), PARAMETER ::   &
(  206)          rhoa  = 1.22,    &  ! air density
(  207)          cp    = 1000.5,  &  ! specific heat of air
(  208)          Lv    = 2.5e6,   &  ! latent heat of vaporization
(  209)          Ls    = 2.839e6, &  ! latent heat of sublimation
(  210)          Stef  = 5.67e-8, &  ! Stefan Boltzmann constant
(  211)          Cice = 1.63e-3      ! transfert coefficient over ice
(  212)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  213)          catm1               ! 
(  214) 
(  215)       !!---------------------------------------------------------------------
(  216) 
(  217) 
(  218) !!DB: adapt flx_bulk_daily routine here
(  219) !!NB: I do not do time-interpolation yet, just input a new frame every day
(  220) !!DB: 
(  221) 
(  222) 
(  223)       rdt_in = rday/rdt    !# of timesteps between input records (float)
(  224)       ndt = int(rday/rdt)  !int version (REM: daily input routine)
(  225)       idy = 365
(  226)       IF(nleapy == 1) idy = 366
(  227) 
(  228) 
(  229) 
(  230)                   
(  231) 
(  232) 
(  233) 
(  234)       !! =============================================== !!
(  235)       !!   PART A:   READING FLUX FILES WHEN NECESSARY   !!






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 10

(  236)       !! =============================================== !!
(  237) 
(  238)       !                                         ! ------------------------ !
(  239)       IF( kt == nit000 ) THEN                   !   first call kt=nit000   !
(  240) 
(  241)          !--- Initializes default values of file names, frequency of forcings 
(  242)          !    and variable to be read in the files, before reading namelist 
(  243) !         clname(1) = 'precip.nc'   ; freqh(1) = -12 ; clvarname(1) = 'precip'     ! monthly
(  244) !         clname(2) = 'u10.nc'      ; freqh(2) =   6 ; clvarname(2) = 'u10'        ! 6 hourly          
(  245) !         clname(3) = 'v10.nc'      ; freqh(3) =   6 ; clvarname(3) = 'v10'        ! 6 hourly    
(  246) !         clname(4) = 'q10.nc'      ; freqh(4) =   6 ; clvarname(4) = 'q10'        ! 6 hourly   
(  247) !         clname(5) = 'radsw.nc'    ; freqh(5) =  24 ; clvarname(5) = 'radsw'      ! daily 
(  248) !         clname(6) = 'radlw.nc'    ; freqh(6) =  24 ; clvarname(6) = 'radlw'      ! daily 
(  249) !         clname(7) = 't10.nc'      ; freqh(7) =   6 ; clvarname(7) = 't10'        ! 6 hourly 
(  250) !         clname(8) = 'snow.nc'     ; freqh(8) = -12 ; clvarname(8) = 'snow'       ! monthly
(  251) 
(  252)          IF(lwp) THEN
(  253)             WRITE(numout2,*)' '
(  254)             WRITE(numout2,*)' DB flx_core : global CORE fields in NetCDF format'
(  255)             WRITE(numout2,*)' ~~~~~~~~~~ '
(  256)             WRITE(numout2,*) 'DB: #-dt between inputs  rdt_in        = ', rdt_in
(  257)             WRITE(numout2,*) 'DB: is_cyclical            = ', is_cyclical
(  258)             WRITE(numout2,*) '           ln_2m        = ', ln_2m
(  259)             WRITE(numout2,*) '           alpha_precip = ', alpha_precip
(  260) !            WRITE(numout2,*) '           list of files and frequency (hour), or monthly (-12) '
(  261)          ENDIF
(  262) 
(  263)          !--- Initialization to zero
(  264)          qsr_oce (:,:) = 0.e0
(  265)          qsb_oce (:,:) = 0.e0
(  266)          qla_oce (:,:) = 0.e0
(  267)          qlw_oce (:,:) = 0.e0
(  268)          qnsr_oce(:,:) = 0.e0
(  269)          qsr_ice (:,:) = 0.e0
(  270)          qnsr_ice(:,:) = 0.e0
(  271)          qla_ice (:,:) = 0.e0
(  272)          qlw_ice (:,:) = 0.e0
(  273)          qsb_ice (:,:) = 0.e0
(  274) 
(  275)          dqns_ice(:,:) = 0.e0
(  276)          dqla_ice(:,:) = 0.e0
(  277)          tprecip (:,:) = 0.e0
(  278)          sprecip (:,:) = 0.e0
(  279)          evap    (:,:) = 0.e0
(  280) 
(  281)          flxnow  (:,:,:) = 1.e-15         !RB
(  282)          flxdta  (:,:,:,:) = 1.e-15       !RB
(  283) 
(  284) 
(  285) 
(  286) 
(  287)          if(is_cyclical) then
(  288) !DB  frame_daily is the offset into the files. 
(  289)             frame_daily = nday_year     
(  290)          else
(  291) !AD:            frame_daily = 1               ! for non-cyclical daily forcing
(  292) !AD:  non-cyclical should start on nday_year in ndate0
(  293)             frame_daily = nday_year     






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 11

(  294) !AD:  insert nyear
(  295) 	    write(clname_1(18:21),'(I4.4)') nyear
(  296) 	    write(clname_2(21:24),'(I4.4)') nyear
(  297) 	    write(clname_3(18:21),'(I4.4)') nyear
(  298) 	    write(clname_4(22:25),'(I4.4)') nyear
(  299) 	    write(clname_5(19:22),'(I4.4)') nyear  
(  300)             nyear0=nyear	                          
(  301)          endif
(  302) 
(  303) !!DB -- crucial assignments
(  304)          call ncdf_get_dim_size(clname_1, 'time_counter',num_recsd,status)
(  305)          if(status == 0) then
(  306)             if(lwp) write(numout2,*)'db_flx_core: # records in ',clname_1, ' = ', num_recsd
(  307)          else
(  308)             if(lwp) write(numout2,*)'DBG: prob with ncdf_get_dim_size ====> STOP'
(  309)             stop
(  310)          endif
(  311)          if(num_recsd /= idy) then
(  312)             if(lwp) then
(  313)                write(numout2,*)'DBG -- db_flx_core: # records /= ', idy, ' Could be an error'
(  314)             endif
(  315)          endif
(  316)          call ncdf_get_dim_size(clname_1, 'x', ipi, status)
(  317)          call ncdf_get_dim_size(clname_1, 'y', ipj, status)
(  318)          if( ipi /= jpidta .OR. ipj /= jpjdta) then
(  319)             if(lwp) then
(  320)                write(numout2,*)
(  321)                write(numout2,*) 'problem with file ', clname_1
(  322)                write(numout2,*) 'problem with dimensions'
(  323)                write(numout2,*) ' ipi ',ipi,' jpidta ',jpidta
(  324)                write(numout2,*) ' ipj ',ipj,' jpjdta ',jpjdta
(  325)             endif
(  326)             stop 'db_flx_core'
(  327)          endif
(  328) 
(  329) 
(  330)          bflx1 =  0 + (frame_daily-1)  !!NB: to get 1st record correct, see below
(  331)          if(bflx1 .gt. num_recsd) bflx1 = num_recsd
(  332)          bflx2 = bflx1 + 1
(  333)          if(bflx2 .gt. num_recsd) bflx2 = num_recsd
(  334)          ireadd = 0 ! next input timestep rel nit000 = 0 to start
(  335) 
(  336) 
(  337)       ENDIF   !! kt==nit000
(  338) 
(  339) !!DB version here:
(  340) 
(  341)     ! Read daily data
(  342)     ! ----------------------------
(  343) 
(  344)     if(kt-nit000 == ireadd) then   
(  345) 
(  346)        bflx1 = bflx1 + 1
(  347)        if(is_cyclical) then    !roll over to next cycle
(  348)           rec1 = mod(bflx1,num_recsd)
(  349)           if(rec1 == 0) rec1 = num_recsd
(  350) 
(  351)        else          






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 12

(  352) 
(  353)           if(nyear .gt. nyear0) then
(  354) !AD:  insert nyear if calendar rolled over
(  355) !   ---> should have a check here to see if file exists
(  356) !        and if not it should use the last value as before
(  357)  	     write(clname_1(18:21),'(I4.4)') nyear
(  358) 	     write(clname_2(21:24),'(I4.4)') nyear
(  359) 	     write(clname_3(18:21),'(I4.4)') nyear
(  360) 	     write(clname_4(22:25),'(I4.4)') nyear
(  361) 	     write(clname_5(19:22),'(I4.4)') nyear                        
(  362) 	     nyear0=nyear
(  363) 	   endif   
(  364) 
(  365) !AD:   Leave just in case	             !clamp if necessary
(  366)           rec1 = nday_year !AD: leave bflx1 alone for ireadd to work
(  367)           if(rec1 .gt. num_recsd) rec1 = num_recsd
(  368) 
(  369) 
(  370)        endif
(  371)        
(  372)        bflx2 = bflx1 + 1
(  373)        if(is_cyclical) then    !roll over to next cycle
(  374)           rec2 = mod(bflx2,num_recsd)
(  375)           if(rec2 == 0) rec2 = num_recsd
(  376)        else                    !clamp if necessary
(  377) 
(  378)            if(nyear .gt. nyear0) then
(  379) !AD:  insert nyear if calendar rolled over
(  380) !   ---> should have a check here to see if file exists
(  381) !        and if not it should use the last value as before
(  382) 
(  383) !         will need to keep 2 file names!!!           
(  384) ! 	     write(clname_1(18:21),'(I4.4)') nyear
(  385) !	     write(clname_2(21:24),'(I4.4)') nyear
(  386) !	     write(clname_3(18:21),'(I4.4)') nyear
(  387) !	     write(clname_4(22:25),'(I4.4)') nyear
(  388) !	     write(clname_5(19:22),'(I4.4)') nyear                        
(  389) !	     nyear0=nyear
(  390) 
(  391) 	   endif   
(  392) !AD:   Leave just in case	             !clamp if necessary
(  393)           rec2 = rec1 + 1 
(  394)           if(rec2 .gt. num_recsd) rec2 = num_recsd 
(  395)        endif
(  396) 
(  397) !!DB: better formula if #-timesteps-per-day is NOT an integer (see tau_forcing.h90)
(  398)        ireadd = nint((bflx1-(frame_daily-1))*rdt_in) !! next input timestep rel nit000       
(  399) !       ireadd = ireadd + ndt  !! next input timestep rel nit000
(  400)        
(  401) !!DB: perpetual_forcing from oce; 
(  402)        if(perpetual_forcing > 0) then
(  403)           bflx1 = frame_daily
(  404)           if(bflx1 .gt. num_recsd) bflx1 = num_recsd
(  405)           bflx2 = bflx1 
(  406)           rec1 = bflx1
(  407)           rec2 = rec1
(  408)           ireadd = -1   !ensures (kt-nit000) /= ireadd, thus only 1 input of files
(  409) !!DBG






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 13

(  410)           if(lwp .AND. kt==nit000) then 
(  411)              write(numout2,*)'DBG: In db_flx_core, perpetual_forcing = ', perpetual_forcing
(  412)              write(numout2,*)'DBG: In db_flx_core, bflx2 = bflx1 = 1'
(  413)           endif
(  414)        endif
(  415) 
(  416) !!DB
(  417) 
(  418)        if(lwp) write(numout2,*)'DBG: In db_flx_core -- reading field ',rec1, 'at kt = ', kt 
(  419)        if(lwp) write(numout2,*)'DBG: In db_flx_core -- next input at ireadd = ', ireadd
(  420)        if(lwp) write(numout2,*)'AD:  In db_flx_core -- nyear,nday_year = ',nyear,nday_year
(  421)        if(lwp) write(numout2,*)'--->',bflx1,frame_daily,rdt_in
(  422)        ! Read data
(  423) !!Adjust using frame_daily 
(  424) 
(  425)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_1
(  426)        call ncdf_read(clname_1, 'T_10_MOD', tatm, -(rec1), status)
(  427)        if(status /= 0) then
(  428)           if(lwp) then
(  429)              write(numout2,*)'DBG: Prob reading tair file, status = ', status
(  430)              !               nstop = nstop + 1
(  431)           endif
(  432)        endif
(  433) !!DB -- done elsewhere or not necessary???
(  434)        ! conversion of temperature Kelvin --> Celsius  [rt0=273.15]
(  435) !       tatm(:,:) = ( tatm(:,:) - rt0 ) 
(  436)        flxdta(:,:,7,1) = tatm(:,:)
(  437) 
(  438)        ! read rain and snow
(  439)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_2
(  440)        call ncdf_read(clname_2, 'RAIN', vatm, -(rec1), status)
(  441)        if(status /= 0) then
(  442)           if(lwp) then
(  443)              write(numout2,*)'DBG: Prob reading wspd file, status = ', status
(  444)              !               nstop = nstop + 1
(  445)           endif
(  446)        endif
(  447)        flxdta(:,:,1,1) = vatm(:,:)       
(  448) 
(  449)        ! read rain and snow
(  450)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_2
(  451)        call ncdf_read(clname_2, 'SNOW', vatm, -(rec1), status)
(  452)        if(status /= 0) then
(  453)           if(lwp) then
(  454)              write(numout2,*)'DBG: Prob reading wspd file, status = ', status
(  455)              !               nstop = nstop + 1
(  456)           endif
(  457)        endif
(  458)        flxdta(:,:,8,1) = vatm(:,:)       
(  459) 
(  460)        ! read humidity
(  461)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_3
(  462)        call ncdf_read(clname_3, 'Q_10_MOD', vatm, -(rec1), status)
(  463)        if(status /= 0) then
(  464)           if(lwp) then
(  465)              write(numout2,*)'DBG: Prob reading hum_cloud file, status = ', status
(  466)              !               nstop = nstop + 1
(  467)           endif






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 14

(  468)        endif
(  469)        flxdta(:,:,4,1) = vatm(:,:)       
(  470) 
(  471)        ! read U and V atm
(  472)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_4
(  473)        call ncdf_read(clname_4, 'U10_MOD', vatm, -(rec1), status)
(  474)        if(status /= 0) then
(  475)           if(lwp) then
(  476)              write(numout2,*)'DBG: Prob reading UV file, status = ', status
(  477)              !               nstop = nstop + 1
(  478)           endif
(  479)        endif
(  480)        flxdta(:,:,2,1) = vatm(:,:)       
(  481) 
(  482)        ! read U and V atm
(  483)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_4
(  484)        call ncdf_read(clname_4, 'V10_MOD', vatm, -(rec1), status)
(  485)        if(status /= 0) then
(  486)           if(lwp) then
(  487)              write(numout2,*)'DBG: Prob reading UV file, status = ', status
(  488)              !               nstop = nstop + 1
(  489)           endif
(  490)        endif
(  491)        flxdta(:,:,3,1) = vatm(:,:)       
(  492) 
(  493)        ! read SW & LW
(  494)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_5
(  495)        call ncdf_read(clname_5, 'SWDN_MOD', vatm, -(rec1), status)
(  496)        if(status /= 0) then
(  497)           if(lwp) then
(  498)              write(numout2,*)'DBG: Prob reading SWLW file, status = ', status
(  499)              !               nstop = nstop + 1
(  500)           endif
(  501)        endif
(  502)        flxdta(:,:,5,1) = vatm(:,:)       
(  503) 
(  504)        ! read SW & LW
(  505)        if(lwp) write(numout2,*)'db_flx_core reading ---> ',clname_5
(  506)        call ncdf_read(clname_5, 'LWDN_MOD', vatm, -(rec1), status)
(  507)        if(status /= 0) then
(  508)           if(lwp) then
(  509)              write(numout2,*)'DBG: Prob reading SWLW file, status = ', status
(  510)              !               nstop = nstop + 1
(  511)           endif
(  512)        endif
(  513)        flxdta(:,:,6,1) = vatm(:,:)       
(  514) 
(  515)        call flush(numout2)
(  516)          
(  517)     ENDIF
(  518) 
(  519) 
(  520) !!DB: basic structure to do time interpolation every ndt_1h
(  521) !!See tau_forced_cmc for more details
(  522) !!NB: would have to have input the bflx2 fields as well (not done above)
(  523) !!DB Use ndt_3h & ndt_1h
(  524) !!Want a formula that also accounts for nature of ndt_* variables which may not
(  525) !!be integers (depending on rdt), and still gives a consistent zxy upon restart






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 15

(  526) !!NB: formula below is incorrect if perpetual_forcing is on so avoid this possibility 
(  527) !      if(mod(kt-nit000,ndt_1h) == 0  .AND.  perpetual_forcing == 0) then
(  528) !         zxy = (float(kt-nit000) - (bflx1-1)*rdt_3h)/rdt_3h
(  529) !         zxy = min(zxy,1.0)  !to be safe
(  530) !         taux(:,:) = (1.-zxy) * taux_dta(:,:,1) + zxy * taux_dta(:,:,2)
(  531) !         tauy(:,:) = (1.-zxy) * tauy_dta(:,:,1) + zxy * tauy_dta(:,:,2)
(  532) !      endif
(  533) 
(  534) !!DB
(  535)     !   Interpolates in time   !
(  536) !! adapt existing routine
(  537) !!OLD
(  538) !      DO ji = 1, jpfile       
(  539) !         IF( freqh(ji) == -12 ) THEN
(  540) !            zxy = REAL(nday) / REAL( nobis(imois) ) + 0.5 - i15         !!! <== Caution nobis hard coded !!!
(  541) !            flxnow(:,:,ji) = ( ( 1. - zxy) * flxdta(:,:,ji,1) + zxy * flxdta(:,:,ji,2) )
(  542) !         ELSE
(  543) !            flxnow(:,:,ji) = flxdta(:,:,ji,1)
(  544) !         ENDIF
(  545) !      END DO
(  546) !      ! JMM : add vatm needed in tracer routines  ==> wind module ???
(  547) !      vatm(:,:) = SQRT( flxnow(:,:,2) * flxnow(:,:,2) + flxnow(:,:,3) * flxnow(:,:,3) )
(  548) 
(  549) !!DB: assign to flxnow here as flxnow is local and flxdta is global and thus
(  550) !!still exist if perpetual_forcing
(  551)     do ji = 1, jpfile
(  552)        flxnow(:,:,ji) = flxdta(:,:,ji,1)
(  553)     enddo
(  554) 
(  555) 
(  556) 
(  557) 
(  558)       !! =============================================== !!
(  559)       !!   PART B:   CORE BULK CALCULATION               !!
(  560)       !! =============================================== !!
(  561) 
(  562)       !  for other forcing cases this is done in modules bulk.F90 and flxblk
(  563) 
(  564)       !                                         ! ------------------------ !
(  565)       !                                         !   Bulk initialisation    !
(  566)       !                                         ! ------------------------ !
(  567) 
(  568) !!DB from flxblk -- in case no ice model
(  569) 
(  570) 
(  571) # 574
(  574)       ! code part from bulk.F90 :
(  575) 
(  576)       IF( kt == nit000) THEN
(  577)          !--- computation of rdtbs2     ===>  !gm  is it really usefull ????
(  578)          IF( nacc == 1 ) THEN
(  579)             rdtbs2 = nfbulk * rdtmin * 0.5
(  580)          ELSE
(  581)            rdtbs2 = nfbulk * rdt * 0.5
(  582)          ENDIF
(  583)          IF( .NOT.ln_rstart ) THEN
(  584)             zcof =  REAL( nfbulk - 1, wp )
(  585)             gsst(:,:) =  zcof * ( tn(:,:,1)  + rt0 )






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 16

(  586)             gsss(:,:) =  zcof * tn_ice(:,:)
(  587)             gu  (:,:) =  zcof * un(:,:,1)
(  588)             gv  (:,:) =  zcof * vn(:,:,1)
(  589)          ENDIF
(  590)       ENDIF
(  591) 
(  592)       ! ----------------------------------------------------------------------------- !
(  593)       !      0   Mean SST, ice temperature and ocean currents over nfbulk pdt         !
(  594)       ! ----------------------------------------------------------------------------- !
(  595)       !
(  596)       gsst(:,:) = gsst(:,:) + (tn(:,:,1) + rt0 )  
(  597)       gsss(:,:) = gsss(:,:) + tn_ice(:,:)
(  598)       gu  (:,:) = gu  (:,:) + un    (:,:,1)
(  599)       gv  (:,:) = gv  (:,:) + vn    (:,:,1)
(  600) 
(  601)       IF( MOD( kt - 1 , nfbulk ) == 0 ) THEN
(  602)          !
(  603) 
(  604) !!gm  better coded:
(  605)          ! Caution set to rt0 over land, not 0 Kelvin (otherwise floating point exception in bulk computation) 
(  606)          zcof =  1. / REAL( nfbulk, wp )
(  607)          zsst(:,:) = gsst(:,:) * zcof * tmask(:,:,1) + rt0 * (1.- tmask(:,:,1))    ! mean sst in K
(  608)          zsss(:,:) = gsss(:,:) * zcof * tmask(:,:,1) + rt0 * (1.- tmask(:,:,1))   ! mean tn_ice in K
(  609) !!gm
(  610) 
(  611)          zut(:,:) = 0.e0      !!gm not necessary but at least first and last column
(  612)          zvt(:,:) = 0.e0
(  613)          !            lb
(  614)          ! Interpolation of surface current at T-point, zut and zvt :
(  615) !!gm better coded
(  616)          zcof =  0.5 / REAL( nfbulk, wp )
(  617)          DO jj = 2, jpjm1
(  618)             DO ji = 2, jpim1
(  619)                zut(ji,jj) = ( gu(ji-1,jj  ) + gu(ji,jj) ) * zcof
(  620)                zvt(ji,jj) = ( gv(ji  ,jj-1) + gv(ji,jj) ) * zcof
(  621)             END DO
(  622)          END DO
(  623) !!gm
(  624) ! FD for ice velocity: B-grid (ice_vel-point) to T-grid
(  625) !!DB -- in case no ice model
(  626) # 627
(  627)          zus(:,:) = 0.
(  628)          zus(:,:) = 0.
(  629)          DO jj = 2, jpjm1
(  630)             DO ji = 2, jpim1
(  631)                zus(ji,jj) = ( u_ice(ji+1,jj  ) + u_ice(ji,jj  ) &
(  632) 	                    + u_ice(ji+1,jj+1) + u_ice(ji,jj+1) ) * 0.25
(  633)                zvs(ji,jj) = ( v_ice(ji+1,jj  ) + v_ice(ji,jj  ) &
(  634) 	                    + v_ice(ji+1,jj+1) + v_ice(ji,jj+1) ) * 0.25
(  635)             END DO
(  636)          END DO
(  637) 
(  638) 
(  639) # 643
(  643)          CALL lbc_lnk( zut(:,:), 'T', -1. )         
(  644)          CALL lbc_lnk( zvt(:,:), 'T', -1. )         
(  645) 
(  646)          ! ----------------------------------------------------------------------------- !






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 17

(  647)          !      I   Radiative FLUXES                                                     !
(  648)          ! ----------------------------------------------------------------------------- !
(  649) 
(  650)          !--- Ocean & Ice Albedo
(  651)          alb_oce_os(:,:) = 0.   ;   alb_oce_cs(:,:) = 0.   !gm is it necessary ???
(  652)          alb_ice_os(:,:) = 0.   ;   alb_ice_cs(:,:) = 0.
(  653)          CALL flx_blk_albedo( alb_ice_os, alb_oce_os, alb_ice_cs, alb_oce_cs )
(  654) 
(  655)          !--- Radiative fluxes over ocean
(  656)          qsr_oce(:,:) = ( 1. - 0.066 ) * flxnow(:,:,5)                                 ! Solar Radiation
(  657)          qlw_oce(:,:) = flxnow(:,:,6) - Stef*zsst(:,:)*zsst(:,:)*zsst(:,:)*zsst(:,:)   ! Long Waves (Infra-red)
(  658) 
(  659)          !--- Radiative fluxes over ice
(  660)          qsr_ice(:,:) = ( 1. - alb_ice_cs(:,:) ) * flxnow(:,:,5)                                  ! Solar Radiation
(  661)          qlw_ice(:,:) = 0.95 * ( flxnow(:,:,6) - Stef*zsss(:,:)*zsss(:,:)*zsss(:,:)*zsss(:,:) )   ! Long Waves (Infra-red)
(  662)          !-------------------------------------------------------------------------------!
(  663) !	 ji=175;jj=21;write(*,*) 'rad',qsr_oce(ji,jj),qlw_oce(ji,jj),qsr_ice(ji,jj),qlw_ice(ji,jj)
(  664) 
(  665)          ! ----------------------------------------------------------------------------- !
(  666)          !     II    Turbulent FLUXES                                                    !
(  667)          ! ----------------------------------------------------------------------------- !
(  668)          !
(  669)          ! scalar wind ( = | U10m - SSU | )
(  670)          ! It is important to take into account the sea surface courant
(  671)          !                lb
(  672)          ! Now, wind components are provided on T-points within the netcdf input file.
(  673)          ! Thus, the wind module is computded at T-points taking into account the sea
(  674) 
(  675)          !--- Module of relative wind
(  676) !        dUnormt(:,:) = sqrt( (flxnow(:,:,2) - zut(:,:))*(flxnow(:,:,2) - zut(:,:))  &
(  677) !             &             + (flxnow(:,:,3) - zvt(:,:))*(flxnow(:,:,3) - zvt(:,:)) )
(  678) !!gm  more efficient coding:
(  679)          DO jj = 1, jpj
(  680)             DO ji = 1, jpi
(  681)                zzu = flxnow(ji,jj,2) - zut(ji,jj)
(  682)                zzv = flxnow(ji,jj,3) - zvt(ji,jj)
(  683)                dUnormt(ji,jj) = SQRT( zzu*zzu + zzv*zzv )
(  684)             END DO 
(  685)          END DO 
(  686) ! FD Unormt for ice from Steven Alderson + difference from ice vel
(  687)          DO jj = 1, jpj
(  688)             DO ji = 1, jpi
(  689)                zzu = flxnow(ji,jj,2) - zus(ji,jj)
(  690)                zzv = flxnow(ji,jj,3) - zvs(ji,jj)
(  691)                Unormt(ji,jj) = SQRT( zzu*zzu + zzv*zzv )
(  692)             END DO 
(  693)          END DO 
(  694) !!gm
(  695)          !--- specific humidity at temp SST
(  696)          qsatw(:,:) = 0.98*640380*exp(-5107.4/zsst(:,:))/rhoa
(  697)          !
(  698)          !--- specific humidity at temp tn_ice
(  699)          qsat(:,:)  = 11637800*exp(-5897.8/zsss(:,:))/rhoa
(  700)          !
(  701)          ! CORE iterartive algo for computation of Cd, Ch, Ce at T-point :
(  702)          ! ===============================================================
(  703)          IF( ln_2m ) THEN
(  704)             IF( kt == nit000 .AND. lwp )   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 18

(  705)                WRITE(numout,*) 
(  706)                WRITE(numout,*) '           Calling TURB_CORE_2Z for bulk transfert coefficients'
(  707)             ENDIF
(  708)             !! If air temp. and spec. hum. are given at different height (2m) than wind (10m) :
(  709)             CALL TURB_CORE_2Z(2.,10.,zsst(:,:), flxnow(:,:,7), qsatw(:,:), flxnow(:,:,4), &
(  710)                &              dUnormt(:,:), Cd(:,:), Ch(:,:), Ce(:,:), zt_zu(:,:), zq_zu(:,:))
(  711)          ELSE
(  712)             IF( kt == nit000 .AND. lwp )    THEN
(  713)                WRITE(numout,*) 
(  714)                WRITE(numout,*) '           Calling TURB_CORE_1Z for bulk transfert coefficients'
(  715)             ENDIF
(  716)             !! If air temp. and spec. hum. are given at same height than wind (10m) :
(  717)             CALL TURB_CORE_1Z(10., zsst(:,:), flxnow(:,:,7), qsatw(:,:), flxnow(:,:,4),  &
(  718)                &              dUnormt(:,:), Cd(:,:), Ch(:,:), Ce(:,:) )
(  719)          ENDIF
(  720) 
(  721) ! FD correction from JM Molines + Zeliang Wang suggestion
(  722)          !   II.1) Momentum over ocean and ice
(  723)          !   ---------------------------------
(  724)          ! Tau_x and Tau_y at T-point
(  725)          tauxtwo(:,:) = rhoa*dUnormt(:,:)*Cd(:,:)*(flxnow(:,:,2) - zut(:,:))
(  726)          tauytwo(:,:) = rhoa*dUnormt(:,:)*Cd(:,:)*(flxnow(:,:,3) - zvt(:,:))
(  727) 
(  728) ! FD Unormt for ice from Steven Alderson
(  729)          ! Tau_x and Tau_y at T-point
(  730) 	 !
(  731) ! FD        tauxtwi(:,:) = rhoa*Unormt(:,:)*Cice*flxnow(:,:,2)
(  732) ! FD        tauytwi(:,:) = rhoa*Unormt(:,:)*Cice*flxnow(:,:,3)
(  733) ! FD do the diff from ice vel
(  734)          tauxtwi(:,:) = rhoa*Unormt(:,:)*Cice*(flxnow(:,:,2)-zus(:,:))
(  735)          tauytwi(:,:) = rhoa*Unormt(:,:)*Cice*(flxnow(:,:,3)-zvs(:,:))
(  736) ! FD        tauxtwi(:,:) = rhoa*dUnormt(:,:)*Cice*flxnow(:,:,2)
(  737) ! FD        tauytwi(:,:) = rhoa*dUnormt(:,:)*Cice*flxnow(:,:,3) 
(  738)          !
(  739)          CALL lbc_lnk( tauxtwo(:,:), 'T', -1. )       !!gm seems not decessary at this point....
(  740)          CALL lbc_lnk( tauxtwi(:,:), 'T', -1. )       !!gm seems not decessary at this point....
(  741)          CALL lbc_lnk( tauytwo(:,:), 'T', -1. )
(  742)          CALL lbc_lnk( tauytwi(:,:), 'T', -1. )
(  743) 
(  744)          !
(  745)          !Tau_x at U-point
(  746)          DO ji = 1, jpim1
(  747)             tauxwo(ji,:) = 0.5*(tauxtwo(ji,:) + tauxtwo(ji+1,:))
(  748)             tauxwi(ji,:) = 0.5*(tauxtwi(ji,:) + tauxtwi(ji+1,:))
(  749)          END DO
(  750)          !
(  751)          ! Tau_y at V-point
(  752)          DO jj = 1, jpjm1
(  753)             tauywo(:,jj) = 0.5*(tauytwo(:,jj) + tauytwo(:,jj+1))
(  754)             tauywi(:,jj) = 0.5*(tauytwi(:,jj) + tauytwi(:,jj+1))
(  755)          END DO
(  756) 
(  757) 
(  758) ! FD what follows was the original lines in flx_core.h90
(  759) ! FD         !   II.1) Momentum over ocean and ice
(  760) ! FD         !   ---------------------------------
(  761) ! FD         ! Tau_x at T-point
(  762) ! FD         tauxt(:,:) = rhoa*dUnormt(:,:)*( (1. - freeze(:,:))*Cd(:,:)*(flxnow(:,:,2) - zut(:,:))    &






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 19

(  763) ! FD              &                          + freeze(:,:)*Cice*flxnow(:,:,2) ) !lb correct pour glace
(  764) ! FD         ! Tau_y at T-point
(  765) ! FD         tauyt(:,:) = rhoa*dUnormt(:,:)*( (1. - freeze(:,:))*Cd(:,:)*(flxnow(:,:,3) - zvt(:,:))    &
(  766) ! FD              &                          + freeze(:,:)*Cice*flxnow(:,:,3) )
(  767) ! FD         !
(  768) ! FD         CALL lbc_lnk( tauxt(:,:), 'T', -1. )       !!gm seems not decessary at this point....
(  769) ! FD         CALL lbc_lnk( tauyt(:,:), 'T', -1. )
(  770) ! FD
(  771) ! FD         !
(  772) ! FD         !Tau_x at U-point
(  773) ! FD         DO ji = 1, jpim1
(  774) ! FD            taux(ji,:) = 0.5*(tauxt(ji,:) + tauxt(ji+1,:))
(  775) ! FD         END DO
(  776) ! FD         !
(  777) ! FD         ! Tau_y at V-point
(  778) ! FD         DO jj = 1, jpjm1
(  779) ! FD            tauy(:,jj) = 0.5*(tauyt(:,jj) + tauyt(:,jj+1))
(  780) ! FD         END DO
(  781) 
(  782) ! FD add just in case
(  783)          taux = tauxwo
(  784)          tauy = tauywo
(  785) 
(  786)          !
(  787)          !   II.2) Turbulent fluxes over ocean
(  788)          !   ---------------------------------
(  789)          !
(  790)          IF( ln_2m ) THEN
(  791)             !!
(  792)             !! Values of temp. and hum. adjusted to 10m must be used instead of 2m values
(  793)             !! Sensible Heat :    ! right sign for ocean
(  794)             qsb_oce(:,:) = rhoa*cp*Ch(:,:)*(zt_zu(:,:) - zsst(:,:))*dUnormt(:,:)
(  795)             !!
(  796)             !! Latent Heat :      ! wrong sign for ocean
(  797)             evap(:,:)    = rhoa*Ce(:,:)*(qsatw(:,:) - zq_zu(:,:))*dUnormt(:,:)
(  798)             !!
(  799)          ELSE
(  800)             !!
(  801)             !! Sensible Heat :    ! right sign for ocean
(  802) !!DBG: 2009.07.27 -- jack up sensible heat flux to see if more ice is formed
(  803) !!ANS=yes; Ch*=2 ----> too much ice (spatial distribution not bad)
(  804) !!      ====> try Ch*=1.5
(  805) !            qsb_oce(:,:) = rhoa*cp*Ch(:,:)*(flxnow(:,:,7) - zsst(:,:))*dUnormt(:,:)
(  806)             qsb_oce(:,:) = (1.5) * rhoa*cp*Ch(:,:)*(flxnow(:,:,7) - zsst(:,:))*dUnormt(:,:)
(  807)             !!
(  808)             !! Latent Heat :      ! wrong sign for ocean
(  809)             evap(:,:)    = rhoa*Ce(:,:)*(qsatw(:,:) - flxnow(:,:,4))*dUnormt(:,:)
(  810)             !!
(  811)          END IF
(  812)          !!
(  813)          qla_oce(:,:) =  -Lv * evap(:,:) ! right sign for ocean
(  814)          
(  815) !	 ji=175;jj=21;write(*,*) 'evap',evap(ji,jj)
(  816) !	 ji=175;jj=21;write(*,*) 'qlatent',qla_oce(ji,jj)
(  817)          ! 
(  818)          !   II.3) Turbulent fluxes over ice
(  819)          !   -------------------------------
(  820)          ! 






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 20

(  821)          ! Sensible Heat :
(  822) ! FD Unormt for ice from Steven Alderson
(  823) ! FD         qsb_ice(:,:) = rhoa*cp*Cice*( flxnow(:,:,7) - zsss(:,:) )*dUnormt(:,:) !lb     use
(  824)          qsb_ice(:,:) = rhoa*cp*Cice*( flxnow(:,:,7) - zsss(:,:) )*Unormt(:,:) !lb     use
(  825)          !                                                                    !lb   dUnormt???
(  826)          ! Latent Heat :                                                      !lb or rather Unormt?
(  827) ! FD Unormt for ice from Steven Alderson
(  828) ! FD         qla_ice(:,:) = Ls*rhoa*Cice*( flxnow(:,:,4) - qsat(:,:) )*dUnormt(:,:)
(  829)          qla_ice(:,:) = Ls*rhoa*Cice*( flxnow(:,:,4) - qsat(:,:) )*Unormt(:,:)
(  830)          !
(  831)          !-------------------------------------------------------------------------------------
(  832) 
(  833) 
(  834)          ! ----------------------------------------------------------------------------- !
(  835)          !     III    Total FLUXES                                                       !
(  836)          ! ----------------------------------------------------------------------------- !
(  837)          !
(  838)          !   III.1) Downward Non Solar flux over ocean
(  839)          !   -----------------------------------------
(  840)          qnsr_oce(:,:) = qlw_oce(:,:) + qsb_oce(:,:) + qla_oce(:,:)
(  841)          !
(  842)          !   III.1) Downward Non Solar flux over ice
(  843)          !   ---------------------------------------
(  844)          qnsr_ice(:,:) = qlw_ice(:,:) + qsb_ice(:,:) + qla_ice(:,:) 
(  845)          !
(  846)          !-------------------------------------------------------------------------------!
(  847) !	 ji=175;jj=21;write(*,*) 'total rad',qnsr_oce(ji,jj),qnsr_ice(ji,jj)
(  848) 
(  849)          ! 6. TOTAL NON SOLAR SENSITIVITY 
(  850) 
(  851)          dqlw_ice(:,:)= 4.0*0.95*Stef*zsss(:,:)*zsss(:,:)*zsss(:,:)
(  852) 
(  853)          ! d qla_ice/ d zsss
(  854) ! FD correction from Steven Alderson
(  855) ! FD         dqla_ice(:,:) = -Ls*Cice*0.98*11637800/(rhoa*rhoa) &
(  856) ! FD              * (-5897.8)/(zsss(:,:)*zsss(:,:)) &
(  857) ! FD              * exp(-5897.8/zsss(:,:)) * dUnormt(:,:)
(  858) ! FD check on derivation:
(  859) ! FD we start from qla_ice = Ls*rhoa*Cice*( flxnow(:,:,4) - qsat )*Unormt
(  860) ! FD 
(  861) ! FD therefore d qla_ice / d qsat = - Ls*rhoa*Cice*Unormt
(  862) ! FD
(  863) ! FD but qsat  = 11637800*exp(-5897.8/zsss)/rhoa
(  864) ! FD
(  865) ! FD therefore d qsat / d zsss = 
(  866) ! FD          11637800/rhoa*(-5897.8) * (-1)/(zsss)**2 * exp(-5897.8/zsss)
(  867) ! FD and therefore
(  868) ! FD d qla_ice / d zsss = - Ls*Cice*Unormt * 
(  869) ! FD              11637800*(+5897.8) * 1/(zsss)**2 * exp(-5897.8/zsss)
(  870) ! FD
(  871) ! FD put sign back in 5897.8 because of the final minus sign below for total dqns_ice
(  872)          dqla_ice(:,:) = - Ls*Cice*11637800 &
(  873)               * (-5897.8)/(zsss(:,:)*zsss(:,:)) &
(  874)               * exp(-5897.8/zsss(:,:)) * Unormt(:,:) 
(  875) 
(  876)          ! d qsb_ice/ d zsss
(  877) ! FD Unormt for ice from Steven Alderson
(  878) ! FD         dqsb_ice(:,:) = rhoa * cp * Cice * dUnormt(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 21

(  879)          dqsb_ice(:,:) = rhoa * cp * Cice * dUnormt(:,:)
(  880) 
(  881)          dqns_ice(:,:) = - ( dqlw_ice(:,:) + dqsb_ice(:,:) + dqla_ice(:,:) )
(  882) 
(  883)          !--------------------------------------------------------------------
(  884)          ! FRACTION of net shortwave radiation which is not absorbed in the
(  885)          ! thin surface layer and penetrates inside the ice cover
(  886)          ! ( Maykut and Untersteiner, 1971 ; Elbert and Curry, 1993 )
(  887) 
(  888) ! FD correction from Steven Alderson
(  889)          catm1(:,:) = 1.0  - 0.3  ! flxnow(:,:,8)
(  890) ! FD         fr1_i0(:,:) = & 
(  891) ! FD              (0.18 * catm1(:,:) + 0.35 * flxnow(:,:,8))
(  892) ! FD         fr2_i0(:,:) = &
(  893) ! FD              (0.82 * catm1(:,:) + 0.65 * flxnow(:,:,8))
(  894)          fr1_i0(:,:) = (0.18 * (1.0-catm1(:,:)) + 0.35 * catm1(:,:) ) 
(  895)          fr2_i0(:,:) = (0.82 * (1.0-catm1(:,:)) + 0.65 * catm1(:,:) ) 
(  896) 
(  897) 
(  898)          !  Total  PRECIPITATION (kg/m2/s)
(  899)          tprecip(:,:) = alpha_precip*flxnow(:,:,1) 
(  900)          ! rename precipitation for freshwater budget calculations:
(  901) !         watm(:,:) =  flxnow(:,:,1)*1000
(  902)          watm(:,:) =  flxnow(:,:,1)*rday
(  903)          !
(  904) 
(  905)          !  SNOW PRECIPITATION  (kg/m2/s)
(  906)          sprecip(:,:) =  alpha_precip*flxnow(:,:,8) 
(  907) 
(  908)          !---------------------------------------------------------------------
(  909) 
(  910)          CALL lbc_lnk( taux (:,:)    , 'U', -1. )
(  911)          CALL lbc_lnk( tauy (:,:)    , 'V', -1. )
(  912)          CALL lbc_lnk( qsr_oce (:,:) , 'T', 1. )
(  913)          CALL lbc_lnk( qnsr_oce(:,:) , 'T', 1. )
(  914)          CALL lbc_lnk( qsr_ice (:,:) , 'T', 1. )
(  915)          CALL lbc_lnk( qnsr_ice(:,:) , 'T', 1. )
(  916)          CALL lbc_lnk( qla_ice (:,:) , 'T', 1. )
(  917)          CALL lbc_lnk( dqns_ice(:,:) , 'T', 1. )
(  918)          CALL lbc_lnk( dqla_ice(:,:) , 'T', 1. )
(  919)          CALL lbc_lnk( fr1_i0  (:,:) , 'T', 1. )
(  920)          CALL lbc_lnk( fr2_i0  (:,:) , 'T', 1. )
(  921)          CALL lbc_lnk( tprecip (:,:) , 'T', 1. )
(  922)          CALL lbc_lnk( sprecip (:,:) , 'T', 1. )
(  923)          CALL lbc_lnk( evap    (:,:) , 'T', 1. )
(  924)          !!
(  925)          !!
(  926)          !!   NEVER mask the windstress!!
(  927)          qsr_oce (:,:) = qsr_oce (:,:)*tmask(:,:,1)
(  928)          qnsr_oce(:,:) = qnsr_oce(:,:)*tmask(:,:,1)
(  929)          qsr_ice (:,:) = qsr_ice (:,:)*tmask(:,:,1)
(  930)          qnsr_ice(:,:) = qnsr_ice(:,:)*tmask(:,:,1)
(  931)          qla_ice (:,:) = qla_ice (:,:)*tmask(:,:,1)
(  932)          dqns_ice(:,:) = dqns_ice(:,:)*tmask(:,:,1)
(  933)          dqla_ice(:,:) = dqla_ice(:,:)*tmask(:,:,1)
(  934)          fr1_i0  (:,:) = fr1_i0  (:,:)*tmask(:,:,1)
(  935)          fr2_i0  (:,:) = fr2_i0  (:,:)*tmask(:,:,1)
(  936)          tprecip (:,:) = tprecip (:,:)*tmask(:,:,1)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 22

(  937)          sprecip (:,:) = sprecip (:,:)*tmask(:,:,1)
(  938)          evap    (:,:) = evap    (:,:)*tmask(:,:,1)
(  939)          gsst(:,:) = 0.e0
(  940)          gsss(:,:) = 0.e0
(  941)          gu  (:,:) = 0.e0
(  942)          gv  (:,:) = 0.e0
(  943) 
(  944)       ENDIF
(  945) 
(  946)       ! ------------------- !
(  947)       ! Last call kt=nitend ! not required in SOPA code
(  948)       ! ------------------- !
(  949) 
(  950)       ! Closing of the numflx file (required in mpp)
(  951) 
(  952) 
(  953)    END SUBROUTINE flx
(  954) 
(  955)    SUBROUTINE TURB_CORE_1Z(zu, sst, T_a, q_sat, q_a,   &
(  956)       &                        dU, Cd, Ch, Ce   )
(  957)       !!----------------------------------------------------------------------
(  958)       !!                      ***  ROUTINE  turb_core  ***
(  959)       !!
(  960)       !! ** Purpose :   Computes turbulent transfert coefficients of surface
(  961)       !!                fluxes according to Large & Yeager (2004)
(  962)       !!
(  963)       !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
(  964)       !!      Momentum, Latent and sensible heat exchange coefficients
(  965)       !!      Caution: this procedure should only be used in cases when air
(  966)       !!      temperature (T_air), air specific humidity (q_air) and wind (dU)
(  967)       !!      are provided at the same height 'zzu'!
(  968)       !!
(  969)       !! References :
(  970)       !!      Large & Yeager, 2004 : ???
(  971)       !! History :
(  972)       !!        !  XX-XX  (???     )  Original code
(  973)       !!   9.0  !  05-08  (L. Brodeau) Rewriting and optimization
(  974)       !!----------------------------------------------------------------------
(  975)       !! * Arguments
(  976) 
(  977)       REAL(wp), INTENT(in) :: zu                 ! altitude of wind measurement       [m]
(  978)       REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::  &
(  979)          sst,       &       ! sea surface temperature         [Kelvin]
(  980)          T_a,       &       ! potential air temperature       [Kelvin]
(  981)          q_sat,     &       ! sea surface specific humidity   [kg/kg]
(  982)          q_a,       &       ! specific air humidity           [kg/kg]
(  983)          dU                 ! wind module |U(zu)-U(0)|        [m/s]
(  984)       REAL(wp), intent(out), DIMENSION(jpi,jpj)  :: &
(  985)          Cd,    &                ! transfert coefficient for momentum       (tau)
(  986)          Ch,    &                ! transfert coefficient for temperature (Q_sens)
(  987)          Ce                      ! transfert coefficient for evaporation  (Q_lat)
(  988) 
(  989)       !! * Local declarations
(  990)       REAL(wp), DIMENSION(jpi,jpj)  ::   &
(  991)          dU10,        &       ! dU                                   [m/s]
(  992)          dT,          &       ! air/sea temperature differeence      [K]
(  993)          dq,          &       ! air/sea humidity difference          [K]
(  994)          Cd_n10,      &       ! 10m neutral drag coefficient






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 23

(  995)          Ce_n10,      &       ! 10m neutral latent coefficient
(  996)          Ch_n10,      &       ! 10m neutral sensible coefficient
(  997)          sqrt_Cd_n10, &       ! root square of Cd_n10
(  998)          sqrt_Cd,     &       ! root square of Cd
(  999)          T_vpot,      &       ! virtual potential temperature        [K]
( 1000)          T_star,      &       ! turbulent scale of tem. fluct.
( 1001)          q_star,      &       ! turbulent humidity of temp. fluct.
( 1002)          U_star,      &       ! turb. scale of velocity fluct.
( 1003)          L,           &       ! Monin-Obukov length                  [m]
( 1004)          zeta,        &       ! stability parameter at height zu
( 1005)          U_n10,       &       ! neutral wind velocity at 10m         [m]   
( 1006)          xlogt, xct, zpsi_h, zpsi_m
( 1007)       !!
( 1008)       INTEGER :: j_itt
( 1009)       INTEGER, PARAMETER :: nb_itt = 3
( 1010)       INTEGER, DIMENSION(jpi,jpj)  ::   &
( 1011)          stab         ! 1st guess stability test integer
( 1012) 
( 1013)       REAL(wp), PARAMETER ::                        &
( 1014)          grav   = 9.8,          &  ! gravity                       
( 1015)          kappa  = 0.4              ! von Karman s constant
( 1016) 
( 1017)       !! * Start
( 1018)       !! Air/sea differences
( 1019)       dU10 = max(0.5, dU)     ! we don't want to fall under 0.5 m/s
( 1020)       dT = T_a - sst          ! assuming that T_a is allready the potential temp. at zzu
( 1021)       dq = q_a - q_sat
( 1022)       !!    
( 1023)       !! Virtual potential temperature
( 1024)       T_vpot = T_a*(1. + 0.608*q_a)
( 1025)       !!
( 1026)       !! Neutral Drag Coefficient
( 1027)       stab    = 0.5 + sign(0.5,dT)    ! stable : stab = 1 ; unstable : stab = 0 
( 1028)       Cd_n10  = 1E-3 * ( 2.7/dU10 + 0.142 + dU10/13.09 )    !   L & Y eq. (6a)
( 1029)       sqrt_Cd_n10 = sqrt(Cd_n10)
( 1030)       Ce_n10  = 1E-3 * ( 34.6 * sqrt_Cd_n10 )               !   L & Y eq. (6b)
( 1031)       Ch_n10  = 1E-3*sqrt_Cd_n10*(18*stab + 32.7*(1-stab)) !   L & Y eq. (6c), (6d)
( 1032)       !!
( 1033)       !! Initializing transfert coefficients with their first guess neutral equivalents :
( 1034)       Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)
( 1035) 
( 1036)       !! * Now starting iteration loop
( 1037)       DO j_itt=1, nb_itt
( 1038)          !! Turbulent scales :
( 1039)          U_star  = sqrt_Cd*dU10                !   L & Y eq. (7a)
( 1040)          T_star  = Ch/sqrt_Cd*dT               !   L & Y eq. (7b)
( 1041)          q_star  = Ce/sqrt_Cd*dq               !   L & Y eq. (7c)
( 1042) 
( 1043)          !! Estimate the Monin-Obukov length :
( 1044)          L  = (U_star**2)/( kappa*grav*(T_star/T_vpot + q_star/(q_a + 1./0.608)) )
( 1045) 
( 1046)          !! Stability parameters :
( 1047)          zeta = zu/L ;  zeta   = sign( min(abs(zeta),10.0), zeta )
( 1048)          zpsi_h = psi_h(zeta)
( 1049)          zpsi_m = psi_m(zeta)
( 1050) 
( 1051)          !! Shifting the wind speed to 10m and neutral stability :
( 1052)          U_n10 = dU10*1./(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)) !  L & Y eq. (9a)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 24

( 1053) 
( 1054)          !! Updating the neutral 10m transfer coefficients :
( 1055)          Cd_n10  = 1E-3 * (2.7/U_n10 + 0.142 + U_n10/13.09)              !  L & Y eq. (6a)
( 1056)          sqrt_Cd_n10 = sqrt(Cd_n10)
( 1057)          Ce_n10  = 1E-3 * (34.6 * sqrt_Cd_n10)                           !  L & Y eq. (6b)
( 1058)          stab    = 0.5 + sign(0.5,zeta)
( 1059)          Ch_n10  = 1E-3*sqrt_Cd_n10*(18.*stab + 32.7*(1-stab))           !  L & Y eq. (6c), (6d)
( 1060) 
( 1061)          !! Shifting the neutral  10m transfer coefficients to ( zu , zeta ) :
( 1062)          !!
( 1063)          xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10) - zpsi_m)
( 1064)          Cd  = Cd_n10/(xct*xct) ;  sqrt_Cd = sqrt(Cd)
( 1065)          !!
( 1066)          xlogt = log(zu/10.) - zpsi_h
( 1067)          !!
( 1068)          xct = 1. + Ch_n10*xlogt/kappa/sqrt_Cd_n10
( 1069)          Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
( 1070)          !!
( 1071)          xct = 1. + Ce_n10*xlogt/kappa/sqrt_Cd_n10
( 1072)          Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
( 1073)          !!
( 1074)       END DO
( 1075)       !!
( 1076)     END SUBROUTINE TURB_CORE_1Z
( 1077) 
( 1078) 
( 1079)     SUBROUTINE TURB_CORE_2Z(zt, zu, sst, T_zt, q_sat, q_zt, dU, Cd, Ch, Ce, T_zu, q_zu)
( 1080)       !!----------------------------------------------------------------------
( 1081)       !!                      ***  ROUTINE  turb_core  ***
( 1082)       !!
( 1083)       !! ** Purpose :   Computes turbulent transfert coefficients of surface 
( 1084)       !!                fluxes according to Large & Yeager (2004).
( 1085)       !!
( 1086)       !! ** Method  :   I N E R T I A L   D I S S I P A T I O N   M E T H O D
( 1087)       !!      Momentum, Latent and sensible heat exchange coefficients
( 1088)       !!      Caution: this procedure should only be used in cases when air
( 1089)       !!      temperature (T_air) and air specific humidity (q_air) are at 2m
( 1090)       !!      whereas wind (dU) is at 10m.
( 1091)       !!
( 1092)       !! References :
( 1093)       !!      Large & Yeager, 2004 : ???
( 1094)       !! History :
( 1095)       !!   9.0  !  06-12  (L. Brodeau) Original code for 2Z
( 1096)       !!----------------------------------------------------------------------
( 1097)       !! * Arguments
( 1098)       REAL(wp), INTENT(in)   :: &
( 1099)          zt,      &     ! height for T_zt and q_zt                   [m]
( 1100)          zu             ! height for dU                              [m]
( 1101)       REAL(wp), INTENT(in), DIMENSION(jpi,jpj) ::  &
( 1102)          sst,      &     ! sea surface temperature              [Kelvin]
( 1103)          T_zt,     &     ! potential air temperature            [Kelvin]
( 1104)          q_sat,    &     ! sea surface specific humidity         [kg/kg]
( 1105)          q_zt,     &     ! specific air humidity                 [kg/kg]
( 1106)          dU              ! relative wind module |U(zu)-U(0)|       [m/s]
( 1107)       REAL(wp), INTENT(out), DIMENSION(jpi,jpj)  ::  &
( 1108)          Cd,       &     ! transfer coefficient for momentum         (tau)
( 1109)          Ch,       &     ! transfer coefficient for sensible heat (Q_sens)
( 1110)          Ce,       &     ! transfert coefficient for evaporation   (Q_lat)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 25

( 1111)          T_zu,     &     ! air temp. shifted at zu                     [K]
( 1112)          q_zu            ! spec. hum.  shifted at zu               [kg/kg]
( 1113) 
( 1114)       !! * Local declarations
( 1115)       REAL(wp), DIMENSION(jpi,jpj) ::  &
( 1116)          dU10,        &     ! dU                                [m/s]
( 1117)          dT,          &     ! air/sea temperature differeence   [K]
( 1118)          dq,          &     ! air/sea humidity difference       [K]
( 1119)          Cd_n10,      &     ! 10m neutral drag coefficient
( 1120)          Ce_n10,      &     ! 10m neutral latent coefficient
( 1121)          Ch_n10,      &     ! 10m neutral sensible coefficient
( 1122)          sqrt_Cd_n10, &     ! root square of Cd_n10
( 1123)          sqrt_Cd,     &     ! root square of Cd
( 1124)          T_vpot_u,    &     ! virtual potential temperature        [K]
( 1125)          T_star,      &     ! turbulent scale of tem. fluct.
( 1126)          q_star,      &     ! turbulent humidity of temp. fluct.
( 1127)          U_star,      &     ! turb. scale of velocity fluct.
( 1128)          L,           &     ! Monin-Obukov length                  [m]
( 1129)          zeta_u,      &     ! stability parameter at height zu
( 1130)          zeta_t,      &     ! stability parameter at height zt
( 1131)          U_n10,       &     ! neutral wind velocity at 10m        [m]
( 1132)          xlogt, xct, zpsi_hu, zpsi_ht, zpsi_m
( 1133) 
( 1134)       INTEGER :: j_itt
( 1135)       INTEGER, PARAMETER :: nb_itt = 3   ! number of itterations
( 1136)       INTEGER, DIMENSION(jpi,jpj) :: &
( 1137)            &     stab                ! 1st stability test integer
( 1138)       REAL(wp), PARAMETER ::                        &
( 1139)          grav   = 9.8,      &  ! gravity                       
( 1140)          kappa  = 0.4          ! von Karman's constant
( 1141) 
( 1142)       !!  * Start
( 1143) 
( 1144)       !! Initial air/sea differences
( 1145)       dU10 = max(0.5, dU)      !  we don't want to fall under 0.5 m/s
( 1146)       dT = T_zt - sst 
( 1147)       dq = q_zt - q_sat
( 1148) 
( 1149)       !! Neutral Drag Coefficient :
( 1150)       stab = 0.5 + sign(0.5,dT)                 ! stab = 1  if dT > 0  -> STABLE
( 1151)       Cd_n10  = 1E-3*( 2.7/dU10 + 0.142 + dU10/13.09 ) 
( 1152)       sqrt_Cd_n10 = sqrt(Cd_n10)
( 1153)       Ce_n10  = 1E-3*( 34.6 * sqrt_Cd_n10 )
( 1154)       Ch_n10  = 1E-3*sqrt_Cd_n10*(18*stab + 32.7*(1 - stab))
( 1155) 
( 1156)       !! Initializing transf. coeff. with their first guess neutral equivalents :
( 1157)       Cd = Cd_n10 ;  Ce = Ce_n10 ;  Ch = Ch_n10 ;  sqrt_Cd = sqrt(Cd)
( 1158) 
( 1159)       !! Initializing z_u values with z_t values :
( 1160)       T_zu = T_zt ;  q_zu = q_zt
( 1161) 
( 1162)       !!  * Now starting iteration loop
( 1163)       DO j_itt=1, nb_itt
( 1164)          dT = T_zu - sst ;  dq = q_zu - q_sat ! Updating air/sea differences
( 1165)          T_vpot_u = T_zu*(1. + 0.608*q_zu)    ! Updating virtual potential temperature at zu
( 1166)          U_star = sqrt_Cd*dU10                ! Updating turbulent scales :   (L & Y eq. (7))
( 1167)          T_star  = Ch/sqrt_Cd*dT              !
( 1168)          q_star  = Ce/sqrt_Cd*dq              !






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 26

( 1169)          !!
( 1170)          L = (U_star*U_star) &                ! Estimate the Monin-Obukov length at height zu
( 1171)               & / (kappa*grav/T_vpot_u*(T_star*(1.+0.608*q_zu) + 0.608*T_zu*q_star))
( 1172)          !! Stability parameters :
( 1173)          zeta_u  = zu/L  ;  zeta_u = sign( min(abs(zeta_u),10.0), zeta_u )
( 1174)          zeta_t  = zt/L  ;  zeta_t = sign( min(abs(zeta_t),10.0), zeta_t )
( 1175)          zpsi_hu = psi_h(zeta_u)
( 1176)          zpsi_ht = psi_h(zeta_t)
( 1177)          zpsi_m  = psi_m(zeta_u)
( 1178)          !!
( 1179)          !! Shifting the wind speed to 10m and neutral stability : (L & Y eq.(9a))
( 1180) !        U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u)))
( 1181)          U_n10 = dU10/(1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m))
( 1182)          !!
( 1183)          !! Shifting temperature and humidity at zu :          (L & Y eq. (9b-9c))
( 1184) !        T_zu = T_zt - T_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
( 1185)          T_zu = T_zt - T_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
( 1186) !        q_zu = q_zt - q_star/kappa*(log(zt/zu) + psi_h(zeta_u) - psi_h(zeta_t))
( 1187)          q_zu = q_zt - q_star/kappa*(log(zt/zu) + zpsi_hu - zpsi_ht)
( 1188)          !!
( 1189)          !! q_zu cannot have a negative value : forcing 0
( 1190)          stab = 0.5 + sign(0.5,q_zu) ;  q_zu = stab*q_zu
( 1191)          !!
( 1192)          !! Updating the neutral 10m transfer coefficients :
( 1193)          Cd_n10  = 1E-3 * (2.7/U_n10 + 0.142 + U_n10/13.09)    ! L & Y eq. (6a)
( 1194)          sqrt_Cd_n10 = sqrt(Cd_n10)
( 1195)          Ce_n10  = 1E-3 * (34.6 * sqrt_Cd_n10)                 ! L & Y eq. (6b)
( 1196)          stab    = 0.5 + sign(0.5,zeta_u)
( 1197)          Ch_n10  = 1E-3*sqrt_Cd_n10*(18.*stab + 32.7*(1-stab)) ! L & Y eq. (6c-6d)
( 1198)          !!
( 1199)          !!
( 1200)          !! Shifting the neutral 10m transfer coefficients to (zu,zeta_u) :
( 1201) !        xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - psi_m(zeta_u))
( 1202)          xct = 1. + sqrt_Cd_n10/kappa*(log(zu/10.) - zpsi_m)
( 1203)          Cd = Cd_n10/(xct*xct) ; sqrt_Cd = sqrt(Cd)
( 1204)          !!
( 1205) !        xlogt = log(zu/10.) - psi_h(zeta_u)
( 1206)          xlogt = log(zu/10.) - zpsi_hu
( 1207)          !!
( 1208)          xct = 1. + Ch_n10*xlogt/kappa/sqrt_Cd_n10
( 1209)          Ch  = Ch_n10*sqrt_Cd/sqrt_Cd_n10/xct
( 1210)          !!
( 1211)          xct = 1. + Ce_n10*xlogt/kappa/sqrt_Cd_n10
( 1212)          Ce  = Ce_n10*sqrt_Cd/sqrt_Cd_n10/xct
( 1213)          !!
( 1214)          !!
( 1215)       END DO
( 1216)       !!
( 1217)     END SUBROUTINE TURB_CORE_2Z
( 1218) 
( 1219)     FUNCTION psi_m(zta)   !! Psis, L & Y eq. (8c), (8d), (8e)
( 1220)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: zta
( 1221) 
( 1222)       REAL(wp), PARAMETER :: pi = 3.14159
( 1223)       REAL(wp), DIMENSION(jpi,jpj)             :: psi_m
( 1224)       REAL(wp), DIMENSION(jpi,jpj)             :: X2, X, stabit
( 1225)       X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.0) ;  X  = sqrt(X2)
( 1226)       stabit    = 0.5 + sign(0.5,zta)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 27

( 1227)       psi_m = -5.*zta*stabit  &                                                  ! Stable
( 1228)            & + (1. - stabit)*(2*log((1. + X)/2) + log((1. + X2)/2) - 2*atan(X) + pi/2)  ! Unstable 
( 1229)     END FUNCTION psi_m
( 1230) 
( 1231)     FUNCTION psi_h(zta)    !! Psis, L & Y eq. (8c), (8d), (8e)
( 1232)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: zta
( 1233) 
( 1234)       REAL(wp), DIMENSION(jpi,jpj)             :: psi_h
( 1235)       REAL(wp), DIMENSION(jpi,jpj)             :: X2, X, stabit
( 1236)       X2 = sqrt(abs(1. - 16.*zta))  ;  X2 = max(X2 , 1.) ;  X  = sqrt(X2)
( 1237)       stabit    = 0.5 + sign(0.5,zta)
( 1238)       psi_h = -5.*zta*stabit  &                                       ! Stable
( 1239)            & + (1. - stabit)*(2.*log( (1. + X2)/2. ))                 ! Unstable
( 1240)     END FUNCTION psi_h
( 1241) 
( 1242) 
( 1243)   SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
( 1244)     !!----------------------------------------------------------------------
( 1245)     !!               ***  ROUTINE flx_blk_albedo  ***
( 1246)     !!
( 1247)     !! ** Purpose :   Computation of the albedo of the snow/ice system
( 1248)     !!      as well as the ocean one
( 1249)     !!
( 1250)     !! ** Method  : - Computation of the albedo of snow or ice (choose the
( 1251)     !!      right one by a large number of tests
( 1252)     !!              - Computation of the albedo of the ocean
( 1253)     !!
( 1254)     !! References :
( 1255)     !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
( 1256)     !!
( 1257)     !! History :
( 1258)     !!  8.0   !  01-04  (LIM 1.0)
( 1259)     !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
( 1260)     !!----------------------------------------------------------------------
( 1261)     !! * Modules used
( 1262)     USE ice                   ! ???
( 1263) 
( 1264)     !! * Arguments
( 1265)     REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
( 1266)          palb         , & !  albedo of ice under overcast sky
( 1267)          palcn        , & !  albedo of ocean under overcast sky
( 1268)          palbp        , & !  albedo of ice under clear sky
( 1269)          palcnp                  !  albedo of ocean under clear sky
( 1270) 
( 1271)     !! * Local variables
( 1272)     INTEGER ::    &
( 1273)          ji, jj                   ! dummy loop indices
( 1274)     REAL(wp) ::   &
( 1275)          c1     = 0.05  , & ! constants values
( 1276)          c2     = 0.1   , &
( 1277)          albice = 0.50  , & !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
( 1278)          cgren  = 0.06  , & !  correction of the snow or ice albedo to take into account
( 1279)          !  effects of cloudiness (Grenfell & Perovich, 1984)
( 1280)          alphd  = 0.80  , & !  coefficients for linear interpolation used to compute
( 1281)          alphdi = 0.72  , & !  albedo between two extremes values (Pyane, 1972)
( 1282)          alphc  = 0.65  , &
( 1283)          zmue   = 0.4   , &     !  cosine of local solar altitude
( 1284)          zzero   = 0.0     ,  &






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 28

( 1285)          zone    = 1.0
( 1286) 
( 1287)     REAL(wp) ::   &
( 1288)          zmue14         , & !  zmue**1.4
( 1289)          zalbpsnm       , & !  albedo of ice under clear sky when snow is melting
( 1290)          zalbpsnf       , & !  albedo of ice under clear sky when snow is freezing
( 1291)          zalbpsn        , & !  albedo of snow/ice system when ice is coverd by snow
( 1292)          zalbpic        , & !  albedo of snow/ice system when ice is free of snow
( 1293)          zithsn         , & !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
( 1294)          zitmlsn        , & !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
( 1295)          zihsc1         , & !  = 1 hsn <= c1 ; = 0 hsn > c1
( 1296)          zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2
( 1297)     REAL(wp), DIMENSION(jpi,jpj) ::  &
( 1298)          zalbfz         , & !  ( = alphdi for freezing ice ; = albice for melting ice )
( 1299)          zficeth                  !  function of ice thickness
( 1300)     LOGICAL , DIMENSION(jpi,jpj) ::  &
( 1301)          llmask
( 1302)     !!    to be included for without seaice
( 1303)     !!      REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
( 1304)     !!      sist   ,   &  !: Sea-Ice Surface Temperature (Kelvin )
( 1305)     !!      hsnif  ,   &  !: Snow thickness
( 1306)     !!      hicif       !: Ice thickness
( 1307)     
( 1308)     !!---------------------------------------------------------------------
( 1309) 
( 1310)     !-------------------------
( 1311)     !  Computation of  zficeth
( 1312)     !--------------------------
( 1313) 
( 1314)     llmask = (hsnif == 0.0) .AND. ( sist >= rt0_ice )
( 1315)     WHERE ( llmask )   !  ice free of snow and melts
( 1316)        zalbfz = albice
( 1317)     ELSEWHERE
( 1318)        zalbfz = alphdi
( 1319)     END WHERE
( 1320) 
( 1321)     DO jj = 1, jpj
( 1322)        DO ji = 1, jpi
( 1323)           IF( hicif(ji,jj) > 1.5 ) THEN
( 1324)              zficeth(ji,jj) = zalbfz(ji,jj)
( 1325)           ELSEIF( hicif(ji,jj) > 1.0  .AND. hicif(ji,jj) <= 1.5 ) THEN
( 1326)              zficeth(ji,jj) = 0.472 + 2.0 * ( zalbfz(ji,jj) - 0.472 ) * ( hicif(ji,jj) - 1.0 )
( 1327)           ELSEIF( hicif(ji,jj) > 0.05 .AND. hicif(ji,jj) <= 1.0 ) THEN
( 1328)              zficeth(ji,jj) = 0.2467 + 0.7049 * hicif(ji,jj)                                &
( 1329)                   &                    - 0.8608 * hicif(ji,jj) * hicif(ji,jj)                 &
( 1330)                   &                    + 0.3812 * hicif(ji,jj) * hicif(ji,jj) * hicif (ji,jj)
( 1331)           ELSE
( 1332)              zficeth(ji,jj) = 0.1 + 3.6 * hicif(ji,jj)
( 1333)           ENDIF
( 1334)        END DO
( 1335)     END DO
( 1336) 
( 1337)     !-----------------------------------------------
( 1338)     !    Computation of the snow/ice albedo system
( 1339)     !-------------------------- ---------------------
( 1340) 
( 1341)     !    Albedo of snow-ice for clear sky.
( 1342)     !-----------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 29

( 1343)      DO jj = 1, jpj
( 1344)         DO ji = 1, jpi
( 1345)            !  Case of ice covered by snow.
( 1346)            !  melting snow
( 1347)            zihsc1       = 1.0 - MAX ( zzero , SIGN ( zone , - ( hsnif(ji,jj) - c1 ) ) )
( 1348)            zalbpsnm     = ( 1.0 - zihsc1 ) & 
( 1349)                 * ( zficeth(ji,jj) + hsnif(ji,jj) * ( alphd - zficeth(ji,jj) ) / c1 ) &
( 1350)                 &                 + zihsc1   * alphd
( 1351)           !  freezing snow
( 1352)            zihsc2       = MAX ( zzero , SIGN ( zone , hsnif(ji,jj) - c2 ) )
( 1353)            zalbpsnf     = ( 1.0 - zihsc2 ) * & 
( 1354)                 ( albice + hsnif(ji,jj) * ( alphc - albice ) / c2 ) &
( 1355)                 &                 + zihsc2   * alphc
( 1356) 
( 1357)            zitmlsn      =  MAX ( zzero , SIGN ( zone , sist(ji,jj) - rt0_snow ) )
( 1358)            zalbpsn      =  zitmlsn * zalbpsnf + ( 1.0 - zitmlsn ) * zalbpsnm
( 1359) 
( 1360)            !  Case of ice free of snow.
( 1361)            zalbpic      = zficeth(ji,jj)
( 1362) 
( 1363)            ! albedo of the system
( 1364)            zithsn       = 1.0 - MAX ( zzero , SIGN ( zone , - hsnif(ji,jj) ) )
( 1365)            palbp(ji,jj) =  zithsn * zalbpsn + ( 1.0 - zithsn ) *  zalbpic
( 1366)         END DO
( 1367)      END DO
( 1368) 
( 1369)      !    Albedo of snow-ice for overcast sky.
( 1370)      !----------------------------------------------
( 1371)      palb(:,:)   = palbp(:,:) + cgren
( 1372) 
( 1373)      !--------------------------------------------
( 1374)      !    Computation of the albedo of the ocean
( 1375)      !-------------------------- -----------------
( 1376) 
( 1377) 
( 1378)      !  Parameterization of Briegled and Ramanathan, 1982
( 1379)      zmue14      = zmue**1.4
( 1380)      palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )
( 1381) 
( 1382)      !  Parameterization of Kondratyev, 1969 and Payne, 1972
( 1383)      palcn(:,:)  = 0.06
( 1384) 
( 1385)    END SUBROUTINE flx_blk_albedo
( 1386) 
( 1387) 
( 1388) # 64 "flxmod.F90"
(   64) 
(   65) 
(   66) 
(   67) # 138
(  138)    SUBROUTINE flx_init
(  139)       !!---------------------------------------------------------------------
(  140)       !!                    ***  ROUTINE flx  ***
(  141)       !!              
(  142)       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(  143)       !!      to the ocean at each time step.
(  144)       !!
(  145)       !! ** Method  :   Constant surface fluxes (read in namelist (namflx))






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 30

(  146)       !!
(  147)       !! ** Action  : - q, qt, qsr, emp, emps, qrp, erp
(  148)       !!
(  149)       !! History :
(  150)       !!        !  91-03  ()  Original code
(  151)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  152)       !!----------------------------------------------------------------------
(  153)       NAMELIST/namflx/ q0, qsr0, emp0, dqdt0, deds0
(  154)       !!---------------------------------------------------------------------
(  155) 
(  156)       ! Read Namelist namflx : surface thermohaline fluxes
(  157)       ! --------------------
(  158)       REWIND ( numnam )
(  159)       READ   ( numnam, namflx )
(  160) 
(  161)       IF(lwp) THEN
(  162)          WRITE(numout,*)' '
(  163)          WRITE(numout,*)'flx_init : thermohaline forcing '
(  164)          WRITE(numout,*)'~~~~~~~~ '
(  165)          WRITE(numout,*)'           net heat flux                  q0   = ', q0  , ' W/m2'
(  166)          WRITE(numout,*)'           solar heat flux                qsr0 = ', qsr0, ' W/m2'
(  167)          WRITE(numout,*)'           net heat flux                  emp0 = ', emp0, ' W/m2'
(  168)          WRITE(numout,*)'           coefficient for SST damping   dqdt0 = ', dqdt0,' W/m2/K'
(  169)          WRITE(numout,*)'           coefficient for SSS damping   deds0 = ', deds0,' mm/day'
(  170)       ENDIF
(  171) 
(  172)    END SUBROUTINE flx_init
(  173) 
(  174) END MODULE flxmod
