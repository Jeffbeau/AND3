


PGF90 (Version     10.5)          09/23/2010  09:49:23      page 1

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

Filename: flxblk_solar.F90

(    1) !!DB 2008.12.04
(    2) !!Modified routine to remove possible conflicts with lim-calls and
(    3) !!for easier use in BGCM_01
(    4) !!-1- changed qsr_oce? to ad_qsr_oce?
(    5) !!-2- the variable ad_qsr_oce() now contains a true daily ave
(    6) !!    (for BGCM_01 use)
(    7) !!-3- #ifdef key_lim_ice ==> ad_qsr_oce(:,:)=frld(:,:)*ad_qsr_oce(:,:)
(    8) !!    where frld is the fraction of open water. Note that this is not
(    9) !!    meant to be a perfect accounting of ice coverage, but rather
(   10) !!    something likely better than nothing. It will set qsr=0 for
(   11) !!    completely ice-covered cells. 
(   12) !!DB 2009.08.17 
(   13) !!Re -2- above: I adjusted the already-daily-averaged insolation value to be
(   14) !!daily-averaged which reduces the insolation. This turns out to be not
(   15) !!so bad as it looks like an average cloudiness correction ===> Keep as is
(   16) !!
(   17) 
(   18) MODULE flxblk_solar
(   19)    !!======================================================================
(   20)    !!                       ***  MODULE  flxblk_solar  ***
(   21)    !!
(   22)    !!=================== SUMMARY =========================================
(   23)    !!Purpose:         Stand alone Short wave flux calcuation
(   24)    !!
(   25)    !!Inputs:         None / CMC humidity cloud cover and air temp
(   26)    !!
(   27)    !!Outputs:             qsr_oce,     &        !: solar flux over ocean (daily average)
(   28)    !!                     qsr_oce1              !: solar flux over ocean (instantaneous) &
(   29)    !!                                                available  only with cmc
(   30)    !! 
(   31)    !! Usage:   -compilie with key_bulk_solar. No other keys required
(   32)    !!          -SUBROUTINE flx_blk_solar( kt) must be called at each timestep
(   33)    !!          -Access qsr_oce and qsr_oce1 with "USE flxblk_solar"
(   34)    !!================= DESCRIPTION =======================================
(   35)    !!
(   36)    !!    SUBROUTINE flx_blk_solar( kt) computes shortwave radiation based
(   37)    !!     on the Zillman (1972). It has been stripped out of flxblk.F90
(   38)    !!     which uses this formula for the open ocean and something else for
(   39)    !!     ice (and as a result is tied to runing opa with lim and bulk
(   40)    !!     forcing).
(   41)    !!   
(   42)    !!    The calculation from flxblk.F90 computes a daily average (kept
(   43)    !!    here as qsr_oce) and is independent of local time.  This module
(   44)    !!    also offers an instantanous flux based on the "real_time"
(   45)    !!    varialble read from the cmc forcing.  read_cmc_flx opens the cmc
(   46)    !!    netcdf forcing files (based on flx_bulk_cmc.h90) and reads the






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47)    !!    nearest 3h forcast).  For this part, the nit000 is assumed to
(   48)    !!    coincide with the first fields in the cmc files. The correct
(   49)    !!    yearday is provided through ndate0 in namelist.  The reading of
(   50)    !!    the "real_time" in the cmc files is accomplished with
(   51)    !!    ncdf_readdate, a new subroutine in lib_ncdf, designed to read a
(   52)    !!    string of length 19 characters from an array of unlimited
(   53)    !!    length. The only time information from the real_time variable that
(   54)    !!    is actually used, is jhour which is the local gmt time. jhour is
(   55)    !!    used to calculate the appropriate time offset based on the
(   56)    !!    longitude of each grid point.
(   57)    !!   
(   58)    !!    If do_read_cmc logical parameter is set to false. The cmc files
(   59)    !!    are not read. The Zillman formula is computed from constant values
(   60)    !!    for humidity, cloud cover and air temperate. These are hardwired
(   61)    !!    in the code (search for do_read_cmc).  In this scenario qsr_oce1
(   62)    !!    is meaningless and only qsr_oce should be used.
(   63)    !!   
(   64)    !!    Last Modify: 28Sep2008 by A. Drozdowski
(   65)    !!   =====================================================================
(   66) 
(   67) 
(   68) 
(   69) #if defined key_bulk_solar || defined key_BGCM_01
(   70)    !!----------------------------------------------------------------------
(   71)    !! * Modules used
(   72)    USE oce             ! ocean dynamics and tracers
(   73)    USE dom_oce         ! ocean space and time domain
(   74)    USE cpl_oce         ! ???
(   75)    USE phycst          ! physical constants
(   76)    USE daymod
(   77) 
(   78)    USE ocfzpt          ! ???
(   79)    USE in_out_manager
(   80)    USE lbclnk
(   81)    USE albedo
(   82)    USE prtctl          ! Print control
(   83) #ifdef key_ice_lim
(   84)    USE ice, ONLY : frld
(   85) #endif
(   86) 
(   87) 
(   88) 
(   89)    IMPLICIT NONE
(   90)    PRIVATE
(   91) 
(   92)    !! * Accessibility
(   93)    PUBLIC flx_blk_solar 
(   94) 
(   95)    !! * Module variables
(   96)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
(   97)       ad_qsr_oce,     &        !: solar flux over ocean (daily average)
(   98)       ad_qsr_oce1,    &          !: solar flux over ocean (instantanous)    
(   99)       frac_day_length            !DB
(  100)    
(  101)    INTEGER, PARAMETER  ::   &
(  102)       jpintsr = 24          ! number of time step between sunrise and sunset
(  103)       !                     ! uses for heat flux computation
(  104)    LOGICAL ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105)       lbulk_init = .TRUE.   ! flag, bulk initialization done or not)
(  106) 
(  107)    !! * constants for bulk computation (flx_blk)
(  108)    REAL(wp)  ::            &  ! constant values
(  109)       zeps    = 1.e-20  ,  &
(  110)       zeps0   = 1.e-13  ,  &
(  111)       zeps1   = 1.e-06  ,  &
(  112)       zzero   = 0.e0    ,  &
(  113)       zone    = 1.0
(  114) 
(  115) !! * constants for solar declinaison computation (flx_blk_declin)
(  116)    REAL(wp) ::                &
(  117)       a0  =  0.39507671   ,   &  ! coefficients
(  118)       a1  = 22.85684301   ,   &
(  119)       a2  = -0.38637317   ,   &
(  120)       a3  =  0.15096535   ,   &
(  121)       a4  = -0.00961411   ,   &
(  122)       b1  = -4.29692073   ,   &
(  123)       b2  =  0.05702074   ,   &
(  124)       b3  = -0.09028607   ,   &
(  125)       b4  =  0.00592797
(  126) 
(  127) !*** Control reading from cmc files
(  128)    LOGICAL do_read_cmc
(  129)    parameter (do_read_cmc=.false.)
(  130)       
(  131)             
(  132) !!! import from blk_oce
(  133)    REAL(wp) ::        &
(  134)       yearday  ,      &  !: number of days per year
(  135)       rdtbs2             !: bulk time step divide by 2
(  136) 
(  137)    INTEGER ::         & !!: namdom : space/time domain (namlist)
(  138)       nfbulk =  5        !: bulk computation frequency 
(  139)       
(  140)             
(  141)    REAL(wp), PRIVATE, DIMENSION(jpi,jpj) ::   &  !:
(  142) !      watm     ,      &  !: precipitation
(  143)       tatm     ,      &  !: atmospheric temperature
(  144)       hatm     ,      &  !: relative humidity
(  145) !      vatm     ,      &  !: wind speed
(  146)       catm               !: percent of cloud cover
(  147) 
(  148) 
(  149) !import from flx_bulk_cmc: subroutine flx was change to read_cmc_flux
(  150) !!DB
(  151)    INTEGER :: nflx1, nflx2, num_recsf, ireadf 
(  152)       
(  153) !A.D: date string to be read from cmc file
(  154)    CHARACTER(LEN=19) :: cmc_date      
(  155)    REAL(wp) :: jhour  !julian hour of the day               
(  156)    
(  157)    !!----------------------------------------------------------------------
(  158)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(  159)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxblk.F90,v 1.8 2005/09/02 15:45:31 opalod Exp $ 
(  160)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(  161)    !!----------------------------------------------------------------------
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163) CONTAINS
(  164) 
(  165)    SUBROUTINE flx_blk_solar( kt)
(  166)       !!---------------------------------------------------------------------------
(  167)       !!                     ***  ROUTINE flx_blk  ***
(  168)       !!                 
(  169)       !!  ** Purpose :   Computation of the soloar heat flux over ocean: Zillman
(  170)       !!
(  171)       !! History :
(  172)       !!   AD: 11Sep08 Modified flx_blk
(  173)       !!----------------------------------------------------------------------
(  174)       !! * Arguments
(  175)       
(  176)       !! * Local variables
(  177)       INTEGER  ::             &
(  178)          kt,ji, jj, jt        ,  &  ! dummy loop indices
(  179)          indaet            ,  &  !  = -1, 0, 1 for odd, normal and leap years resp.
(  180)          iday              ,  &  ! integer part of day
(  181)          indxb             ,  &  ! index for budyko coefficient        
(  182)          indxc             , &     ! index for cloud depth coefficient
(  183) !A.D:         
(  184)          lhour         
(  185) 
(  186)       REAL(wp)  ::            & 
(  187)          zalat , zclat     ,  &  ! latitude in degrees 
(  188)          zmt1, zmt2, zmt3  ,  &  ! tempory air temperatures variables
(  189)          ztatm3, ztatm4    ,  &  ! power 3 and 4 of air temperature
(  190)          z4tatm3           ,  &  ! 4 * zftatm3
(  191)          zcmue             ,  &  ! cosine of local solar altitude
(  192)          zcmue2            ,  &  ! root of zcmue1 
(  193)          zpcmue            ,  &  ! zcmue1**1.4
(  194)          zdecl             ,  &  ! solar declination
(  195)          zsdecl , zcdecl   ,  &  ! sine and cosine of solar declination 
(  196)          zalbo             ,  &  ! albedo of sea-water
(  197)          zalbi             ,  &  ! albedo of ice
(  198)          ztamr             ,  &  ! air temperature minus triple point of water (rtt)
(  199)          ztaevbk           ,  &  ! part of net longwave radiation
(  200)          zevi , zevo       ,  &  ! vapour pressure of ice and ocean 
(  201)          zpis2             ,  &  ! pi / 2
(  202)          z2pi                    ! 2 * pi 
(  203) 
(  204)       REAL(wp)  ::            & 
(  205)          zxday             ,  &  ! day of year
(  206)          zdist             ,  &  ! distance between the sun and the earth during the year
(  207)          zdaycor           ,  &  ! corr. factor to take into account the variation of 
(  208)          !                       ! zday when calc. the solar rad.    
(  209)          zesi, zeso        ,  &  ! vapour pressure of ice and ocean at saturation
(  210)          zesi2             ,  &  ! root of zesi 
(  211)          zqsato            ,  &  ! humidity close to the ocean surface (at saturation)   
(  212)          zqsati            ,  &  ! humidity close to the ice surface (at saturation) 
(  213)          zqsati2           ,  &  ! root of  zqsati 
(  214)          zdesidt           ,  &  ! derivative of zesi, function of ice temperature
(  215)          zdteta            ,  &  ! diff. betw. sst and air temperature
(  216)          zdeltaq           ,  &  ! diff. betw. spec. hum. and hum. close to the surface
(  217)          ztvmoy, zobouks   ,  &  ! tempory scalars
(  218)          zpsims, zpsihs, zpsils, zobouku, zxins, zpsimu ,  & 
(  219)          zpsihu, zpsilu, zstab,zpsim, zpsih, zpsil      ,  & 
(  220)          zvatmg, zcmn, zchn, zcln, zcmcmn, zdenum       ,  & 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(  221)          zdtetar, ztvmoyr, zlxins, zcmn2, zchcm, zclcm , zcoef
(  222) 
(  223)       REAL(wp)  ::            & 
(  224)          zrhova            ,  &  ! air density per wind speed
(  225)          zcsho , zcleo     ,  &  ! transfer coefficient over ocean
(  226)          zcshi , zclei     ,  &  ! transfer coefficient over ice-free
(  227)          zrhovacleo        ,  &  ! air density per wind speed per transfer coef.
(  228)          zrhovacsho, zrhovaclei, zrhovacshi, & 
(  229)          ztice3            ,  &  ! power 3 of ice temperature
(  230)          zticemb, zticemb2 ,  &  ! tempory air temperatures variables
(  231)          zdl, zdr                ! fractionnal part of latitude
(  232)          
(  233)       REAL(wp), DIMENSION(jpi,jpj) :: & 
(  234)          zpatm            ,  &   ! atmospheric pressure
(  235)          zqatm            ,  &   ! specific humidity
(  236)          zes              ,  &   ! vapour pressure at saturation
(  237)          zev, zevsqr      ,  &   ! vapour pressure and his square-root
(  238)          zrhoa            ,  &   ! air density
(  239)          ztatm            ,  &   ! air temperature in Kelvins
(  240)          zcatm1                  ! fraction of cloud
(  241) 
(  242)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  243)          zalboos          ,  &   ! albedo of ocean under overcast sky
(  244)          zalbics          ,  &   ! albedo of ice under clear sky
(  245)          zalbios          ,  &   ! albedo of ice under overcast sky
(  246)          zalbomu          ,  &   ! albedo of ocean when zcmue is 0.4
(  247)          zqsro            ,  &   ! solar radiation over ocean
(  248)          zcldcor          ,  &   ! cloud correction
(  249)          zlsrise, zlsset  ,  &   ! sunrise and sunset
(  250)          zlmunoon         ,  &   ! local noon solar altitude
(  251)          zdlha            ,  &   ! length of the ninstr segments of the solar day
(  252)          zps              ,  &   ! sine of latitude per sine of solar decli.
(  253)          zpc                     ! cosine of latitude per cosine of solar decli. 
(  254) 
(  255)       REAL(wp), DIMENSION(jpi,jpj,jpintsr) ::    &
(  256)          zlha             ,  &   ! local hour angle
(  257)          zalbocs          ,  &   ! tempory var. of ocean albedo under clear sky
(  258)          zsqsro                  ! tempory var. of solar rad. over ocean 
(  259) 
(  260)       !!---------------------------------------------------------------------
(  261) 
(  262)       !---------------------
(  263)       !   Initilization    !
(  264)       !---------------------
(  265) 
(  266)       !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
(  267)       !  and the correction factor for taking into account  the effect of clouds 
(  268)       !------------------------------------------------------
(  269)       IF( lbulk_init ) THEN
(  270)          DO jj = 1, jpj  
(  271)             DO ji = 1 , jpi
(  272)                zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
(  273)                zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
(  274)                indxb          = 1 + INT( zalat ) 
(  275)                indxc          = 1 + INT( zclat )  
(  276)                zdl            = zclat - INT( zclat ) 
(  277)                zdr            = 1.0 - zdl
(  278)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 6

(  279)          END DO
(  280)          IF( nleapy == 1 ) THEN
(  281)             yearday = 366.e0
(  282)          ELSE IF( nleapy == 0 ) THEN
(  283)             yearday = 365.e0
(  284)          ELSEIF( nleapy == 30) THEN
(  285)             yearday = 360.e0
(  286)          ENDIF
(  287) 
(  288) !!! import from blk_oce                  
(  289)       ! computation of rdtbs2
(  290)         IF( nacc == 1 ) THEN
(  291)            rdtbs2 = nfbulk * rdtmin * 0.5
(  292)         ELSE
(  293)            rdtbs2 = nfbulk * rdt * 0.5
(  294)         ENDIF
(  295)                  
(  296)          
(  297)          lbulk_init = .FALSE.
(  298)       ENDIF
(  299) 
(  300) 
(  301)       
(  302) 
(  303) !#if defined key_flx_bulk_cmc
(  304)       IF (do_read_cmc)  then
(  305)            call read_cmc_flx( kt )
(  306)       ELSE
(  307)    !*** Setting Default values       
(  308)          tatm(:,:)=10.
(  309)          catm(:,:)=0.10      
(  310)          hatm(:,:)=0.4     
(  311)       ENDIF
(  312) !#endif      
(  313)         
(  314)             
(  315)       IF( MOD( kt - nit000 , nfbulk ) == 0 ) THEN
(  316)          zpis2       = rpi / 2.
(  317)          z2pi        = 2. * rpi
(  318)    
(  319)         
(  320)          
(  321)          DO jj = 1, jpj
(  322)             DO ji = 1, jpi
(  323)   
(  324)                ztatm (ji,jj) = 273.15 + tatm  (ji,jj)  !  air temperature in Kelvins 
(  325)                zcatm1(ji,jj) = 1.0    - catm  (ji,jj)  !  fractional cloud cover
(  326)                zpatm(ji,jj)  = 101000.               !  pressure 
(  327)          
(  328)                !  Computation of air density, obtained from the equation of state for dry air. 
(  329)                zrhoa(ji,jj) = zpatm(ji,jj) / ( 287.04 * ztatm(ji,jj) )
(  330)          
(  331)                !  zes : Saturation water vapour
(  332)                ztamr = ztatm(ji,jj) - rtt
(  333)                zmt1  = SIGN( 17.269, ztamr )
(  334)                zmt2  = SIGN( 21.875, ztamr )
(  335)                zmt3  = SIGN( 28.200, -ztamr )
(  336)                zes(ji,jj) = 611.0 * EXP (  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 7

(  337)                   &                      / ( ztatm(ji,jj) - 35.86  + MAX( zzero, zmt3 ) ) )
(  338)    
(  339)                !  zev : vapour pressure  (hatm is relative humidity)  
(  340)                zev(ji,jj)   = hatm(ji,jj) * zes(ji,jj) 
(  341)                !  square-root of vapour pressure
(  342)                zevsqr(ji,jj) = SQRT( zev(ji,jj) * 0.01 )
(  343)                !  zqapb  : specific humidity 
(  344)                zqatm(ji,jj) = 0.622 * zev(ji,jj) / ( zpatm(ji,jj) - 0.378 * zev(ji,jj) )
(  345)    
(  346)    
(  347)                !----------------------------------------------------
(  348)                !   Computation of snow precipitation (Ledley, 1985) |
(  349)                !----------------------------------------------------
(  350)    
(  351)                zmt1  =   253.0 - ztatm(ji,jj)
(  352)                zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0 
(  353)                zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0
(  354)   
(  355)             END DO
(  356)          END DO
(  357)   
(  358)         !----------------------------------------------------------
(  359)         !   Computation of albedo (need to calculates heat fluxes)|
(  360)         !-----------------------------------------------------------
(  361)         
(  362)          CALL flx_blk_albedo( zalbios, zalboos, zalbics, zalbomu )
(  363)   
(  364)         !-------------------------------------
(  365)         !   Computation of solar irradiance. |
(  366)         !----------------------------------------
(  367)          indaet   = 1  
(  368)         !  compution of the day of the year at which the fluxes have to be calculate 
(  369)         !--The date corresponds to the middle of the time step.
(  370)          zxday=nday_year + rdtbs2/rday
(  371)   
(  372)          iday   = INT( zxday )
(  373)    
(  374)          IF(ln_ctl) CALL prt_ctl_info('declin : iday ', ivar1=iday, clinfo2=' nfbulk= ', ivar2=nfbulk)
(  375)    
(  376)          !   computation of the solar declination, his sine and his cosine
(  377)          CALL flx_blk_declin( indaet, iday, zdecl )
(  378)          
(  379)          zdecl    = zdecl * rad
(  380)          zsdecl   = SIN( zdecl )
(  381)          zcdecl   = COS( zdecl )
(  382)          
(  383)          !  correction factor added for computation of shortwave flux to take into account the variation of
(  384)          !  the distance between the sun and the earth during the year (Oberhuber 1988)
(  385)          zdist    = zxday * z2pi / yearday
(  386)          zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )
(  387)    
(  388)          DO jj = 1, jpj
(  389)             DO ji = 1, jpi
(  390)                !  product of sine of latitude and sine of solar declination
(  391)                zps     (ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
(  392)                !  product of cosine of latitude and cosine of solar declination
(  393)                zpc     (ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
(  394)                !  computation of the both local time of sunrise and sunset






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 8

(  395)                zlsrise (ji,jj) = ACOS( - SIGN( zone, zps(ji,jj) ) * MIN( zone, SIGN( zone, zps(ji,jj) )  &
(  396)                   &                     * ( zps(ji,jj) / zpc(ji,jj) ) ) ) 
(  397)                zlsset  (ji,jj) = - zlsrise(ji,jj)
(  398)                !  dividing the solar day into jpintsr segments of length zdlha
(  399)                zdlha   (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jpintsr )
(  400)                !  computation of the local noon solar altitude
(  401)                zlmunoon(ji,jj) = ASIN ( ( zps(ji,jj) + zpc(ji,jj) ) ) / rad
(  402)                
(  403)                !  cloud correction taken from Reed (1977) (imposed lower than 1)
(  404)                zcldcor (ji,jj) = MIN( zone, ( 1.e0 - 0.62 * catm(ji,jj) + 0.0019 * zlmunoon(ji,jj) ) )
(  405)             END DO
(  406)          END DO
(  407)    
(  408)             !  Computation of solar heat flux at each time of the day between sunrise and sunset. 
(  409)             !  We do this to a better optimisation of the code 
(  410)             !------------------------------------------------------       
(  411)    
(  412)          DO jt = 1, jpintsr   
(  413)             zcoef = FLOAT( jt ) - 0.5
(  414)             DO jj = 1, jpj
(  415)                DO ji = 1, jpi
(  416)                   !  local hour angle
(  417)                   zlha (ji,jj,jt) = COS ( zlsrise(ji,jj) - zcoef * zdlha(ji,jj) )
(  418)    
(  419)                   ! cosine of local solar altitude
(  420)                   zcmue              = MAX ( zzero ,   zps(ji,jj) + zpc(ji,jj) * zlha (ji,jj,jt)  )
(  421)                   zcmue2             = 1368.0 * zcmue * zcmue
(  422)                   zpcmue             = zcmue**1.4
(  423)                   ! computation of sea-water albedo (Payne, 1972)
(  424)                   zalbocs(ji,jj,jt)  = 0.05 / ( 1.1 * zpcmue + 0.15 )
(  425)                   zalbo              = zcatm1(ji,jj) * zalbocs(ji,jj,jt) + catm(ji,jj) * zalboos(ji,jj)
(  426)                   
(  427)                   ! solar heat flux absorbed at ocean surfaces (Zillman, 1972)
(  428)                   zevo               = zev(ji,jj) * 1.0e-05
(  429)                   zsqsro(ji,jj,jt)   =  ( 1.0 - zalbo ) * zdlha(ji,jj) * zcmue2                &
(  430)                                        / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
(  431) !A.D Warning!!!! not appling albedo for testing purposes                  
(  432) !                  zsqsro(ji,jj,jt)   =  zdlha(ji,jj) * zcmue2                &
(  433) !                                      / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
(  434)                
(  435)                END DO
(  436)             END DO
(  437)          END DO
(  438)    
(  439)   
(  440)          !  Computation of daily (between sunrise and sunset) solar heat flux absorbed 
(  441)          !  at the ocean and snow/ice surfaces.
(  442)          !--------------------------------------------------------------------
(  443)    
(  444)          zqsro   (:,:) = 0.e0
(  445)          DO jt = 1, jpintsr 
(  446)             DO jj = 1, jpj
(  447)                DO ji = 1, jpi  
(  448)                   zqsro  (ji,jj)  = zqsro   (ji,jj) + zsqsro  (ji,jj,jt)
(  449)                END DO
(  450)             END DO
(  451)          END DO
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 9

(  453)  
(  454)          
(  455)                  
(  456)          DO jj = 1, jpj
(  457)             DO ji = 1, jpi 
(  458)    
(  459)                !------------------------------------------- 
(  460)                !  Computation of shortwave radiation.
(  461)                !-------------------------------------------
(  462)                ! For ocean
(  463)                lhour=mod(int(jhour+glamt(ji,jj)/360.*24.),24)  !compute local hour  
(  464) 
(  465)                ad_qsr_oce(ji,jj)  = srgamma*zcldcor(ji,jj)*zqsro(ji,jj) / z2pi
(  466)                ad_qsr_oce1(ji,jj) = REAL( jpintsr )*srgamma*zcldcor(ji,jj)*zsqsro(ji,jj,lhour) / z2pi
(  467) !!DB -- NB: 
(  468) !               frac_day_length(ji,jj) = 24.0*(zlsrise(ji,jj)-zlsset(ji,jj))/z2pi
(  469)                frac_day_length(ji,jj) = (zlsrise(ji,jj)-zlsset(ji,jj))/z2pi
(  470)                               
(  471) !Diagnostic 1
(  472) !                if (jj .eq. 10 .AND. ji .eq. 10) then               
(  473) !                write(9090+narea,'(I7,F8.4,F8.1,2(F8.4))')kt,jhour,tatm(ji,jj),gphit(ji,jj),glamt(ji,jj)
(  474) !                endif
(  475)                                 
(  476)                if ((jj -1 + njmpp) .eq. 80 .AND. (ji -1 + nimpp) .eq. 135) then
(  477) !                  write(9081,'(I7,30(F6.1))') kt,srgamma,zcldcor(ji,jj),z2pi,zone,zzero,zdaycor
(  478) !                  write(9080,'(I7,30(F6.1))') kt,ad_qsr_oce(ji,jj),zqsro(ji,jj),zsqsro(ji,jj,:)
(  479) !DB
(  480)                   write(9080,'(I7,I3,2(F9.4),30(F6.1,1x))') kt,lhour,gphit(ji,jj),glamt(ji,jj) &
(  481)                    ,ad_qsr_oce(ji,jj),ad_qsr_oce1(ji,jj), zlsrise(ji,jj),zlsset(ji,jj),frac_day_length(ji,jj), &
(  482)                    frac_day_length(ji,jj)*ad_qsr_oce(ji,jj)
(  483)                endif
(  484) !!DB: convert qsr to average over 24h instead of average over day_length
(  485) !!2009.08.17 --NB: This correction is in fact an error (see note at top of file) BUT
(  486) !!                 it serves as a sort-of cloudiness correction (reduces insolation)
(  487) !!                 ====> keep as is
(  488)                ad_qsr_oce(ji,jj)  =  ad_qsr_oce(ji,jj)*frac_day_length(ji,jj)
(  489)                
(  490)    
(  491)             END DO
(  492)          END DO
(  493) 
(  494) 
(  495) #ifdef key_ice_lim
(  496) !!DB: use frac open water to adjust insolation
(  497) !!Note the adjustment in case frld=0==> ice covered cell
(  498) !!This guards against zero insolation with a 1% error otherwise
(  499)          ad_qsr_oce(:,:) = ad_qsr_oce(:,:) * (frld(:,:)+0.01)
(  500) #endif
(  501)    
(  502) !!DBG: this causes some zero values which are not desired. 
(  503) !!===> do not call
(  504) !         CALL lbc_lnk( ad_qsr_oce (:,:) , 'T', 1. )
(  505) 
(  506) !!DB -- NB: DO NOT MASK
(  507) !         ad_qsr_oce (:,:) = ad_qsr_oce (:,:)*tmask(:,:,1)
(  508) !         write(3500+narea,*)-kt, minval(ad_qsr_oce(:,:)),minloc(ad_qsr_oce(:,:))
(  509) 
(  510)       ENDIF  






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 10

(  511)   
(  512)    END SUBROUTINE flx_blk_solar
(  513)   
(  514)    
(  515)    
(  516)   
(  517)    SUBROUTINE flx_blk_declin( ky, kday, pdecl )
(  518)       !!---------------------------------------------------------------------------
(  519)       !!               ***  ROUTINE flx_blk_declin  ***
(  520)       !!          
(  521)       !! ** Purpose :   Computation of the solar declination for the day
(  522)       !!         kday ( in decimal degrees ).
(  523)       !!       
(  524)       !! ** Method  :
(  525)       !!
(  526)       !! History :
(  527)       !!         original    : 01-04 (LIM)
(  528)       !!         addition    : 02-08 (C. Ethe, G. Madec)
(  529)       !!---------------------------------------------------------------------
(  530)       !! * Arguments
(  531)       INTEGER, INTENT( in ) ::   &
(  532)          ky  ,        &  ! = -1, 0, 1 for odd, normal and leap years resp.
(  533)          kday            ! day of the year ( kday = 1 on january 1)
(  534)       REAL(wp), INTENT(out) ::  &
(  535)          pdecl            ! solar declination
(  536) 
(  537)       !! * Local variables
(  538)       REAL(wp) ::                & 
(  539)          zday                ,   &  ! corresponding day of type year (cf. ky)
(  540)          zp1, zp2, zp3, zp4         ! temporary scalars
(  541)       !!---------------------------------------------------------------------
(  542)       
(  543)       zday = FLOAT( kday ) 
(  544)       
(  545)       IF( ky == 1 )  THEN 
(  546)          zday = zday - 0.5
(  547)       ELSEIF( ky == 3 )  THEN
(  548)          zday = zday - 1.
(  549)       ELSE 
(  550)          zday = REAL( kday )
(  551)       ENDIF
(  552)       
(  553)       zp1 = rpi * ( 2.0 * zday - 367.0 ) / yearday
(  554)       zp2 = 2. * zp1
(  555)       zp3 = 3. * zp1
(  556)       zp4 = 4. * zp1
(  557)       
(  558)       pdecl  = a0                                                                      &
(  559)          &   + a1 * COS( zp1 ) + a2 * COS( zp2 ) + a3 * COS( zp3 ) + a4 * COS( zp4 )   &
(  560)          &   + b1 * SIN( zp1 ) + b2 * SIN( zp2 ) + b3 * SIN( zp3 ) + b4 * SIN( zp4 )
(  561)         
(  562)    END SUBROUTINE flx_blk_declin
(  563)   
(  564)      
(  565)    SUBROUTINE read_cmc_flx( kt )
(  566)        !!---------------------------------------------------------------------
(  567)        !!                    ***  ROUTINE flx  ***
(  568)        !!        






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 11

(  569)        !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(  570)        !!      to the ocean at each time step.
(  571)        !!
(  572)        !!       Read several cmc fluxes file
(  573)        !!             temperature at 2m   tatm   (K)
(  574)        !!             relative humidite   hatm   (%)
(  575)  !      !!             wind speed          vatm   (m/s)
(  576)  !      !!             precip              watm   (kg/m2/day)
(  577)        !!             clouds              catm   (%)
(  578)        !!
(  579)        !! caution : now, in the opa global model, the net upward water flux is
(  580)        !! -------   with mm/day unit.
(  581)        !!
(  582)        !! History :
(  583)        !!        !  91-03  (O. Marti and Ph Dandin)  Original code
(  584)        !!        !  92-07  (M. Imbard)
(  585)        !!        !  96-11  (E. Guilyardi)  Daily AGCM input files
(  586)        !!        !  99-11  (M. Imbard)  NetCDF FORMAT with ioipsl
(  587)        !!        !  00-05  (K. Rodgers) Daily Netcdf
(  588)        !!   8.5  !  02-09  (C. Ethe and G. Madec)  F90: Free form and MODULE 
(  589)        !!----------------------------------------------------------------------
(  590)        !! * modules used
(  591)  !!DB
(  592)        USE oce, only :  perpetual_forcing, ramp
(  593)        USE lib_ncdf
(  594)  
(  595)        !! * arguments
(  596)        INTEGER, INTENT( in  ) ::   kt ! ocean time step
(  597)  
(  598)        !! * Local declarations      
(  599)        CHARACTER(len=45)  ::  &
(  600)           clname_n = 'tair_3h.nc',        &
(  601)           clname_c = 'hum_cloud_3h.nc',   &
(  602)           clname_x = 'rain_3h.nc',        &
(  603)           clname_w = 'wspd_3h.nc'
(  604)  !!DB
(  605)        REAL(wp) :: zxy, update_frac, rdt_3h
(  606)        INTEGER :: i, j, frame0, status, ndt_3h, ndt_1h
(  607)  
(  608)  
(  609)        !!---------------------------------------------------------------------
(  610)  
(  611)  
(  612)        ! Initialization
(  613)        ! --------------
(  614)        rdt_3h = 3.0*3600./rdt
(  615)        ndt_3h = nint(3.0*3600./rdt)
(  616)        ndt_1h = nint(1.0*3600./rdt) 
(  617)  !!DB  frame0 is a hardwired offset into the files. Use for the (likely restart) 
(  618)  !!    case where  ntau1 = (kt-nit000)/ndt_3h + 1 points to the wrong frame
(  619)        frame0 = 0
(  620)  
(  621)  
(  622)        ! 1. first call kt = nit000
(  623)        ! -------------------------
(  624)        
(  625)        IF( kt == nit000 ) THEN
(  626)  






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 12

(  627)  !!DB -- crucial assignments
(  628)  !!Use new routine to get #-records
(  629)           call ncdf_get_dim_size(clname_x, 'time_counter',num_recsf,status)
(  630)           if(status == 0) then
(  631)              if(lwp) write(numout,*)'read_cmc_flx: # records in ',clname_x, ' = ', num_recsf
(  632)           else
(  633)              if(lwp) write(numout,*)'read_cmc_flx: prob with ncdf_get_dim_size ====> STOP'
(  634)              stop
(  635)           endif
(  636)  
(  637)           nflx1 =  0
(  638)           if(nflx1 .gt. num_recsf) nflx1 = num_recsf
(  639)           nflx2 = nflx1 + 1
(  640)           if(nflx2 .gt. num_recsf) nflx2 = num_recsf
(  641)           ireadf = 0 ! next input timestep rel nit000 = 0 to start
(  642)  
(  643)  
(  644)      ENDIF
(  645)  
(  646)  
(  647)      ! 2. Read (3h) cmc data
(  648)      ! ----------------------------
(  649)  
(  650)      if(kt-nit000 == ireadf) then   !!orig DB statement
(  651)  
(  652)  
(  653)         nflx1 = nflx1 + 1
(  654)  !!DBG if using daily OMIP must turn off these checks
(  655)         if(nflx1 .gt. num_recsf) nflx1 = num_recsf
(  656)         nflx2 = nflx1 + 1
(  657)         if(nflx2 .gt. num_recsf) nflx2 = num_recsf
(  658)         ireadf = nint((nflx1)*rdt_3h) !! next input timestep rel nit000
(  659)         
(  660)  !!DB: perpetual_forcing from oce; 
(  661)         if(perpetual_forcing > 0) then
(  662)            nflx1 = 1
(  663)            nflx2 = nflx1 
(  664)  !!DBG
(  665)            if(lwp .AND. kt==nit000) then 
(  666)               write(numout,*)'read_cmc_flx: perpetual_forcing = ', perpetual_forcing
(  667)               write(numout,*)'read_cmc_flx: nflx2 = nflx1 = 1'
(  668)            endif
(  669)         endif
(  670)  
(  671)  !!DB
(  672)         if(lwp) write(numout,*)'read_cmc_flx:  -- reading field ',nflx1, 'at kt = ', kt 
(  673)         if(lwp) write(numout,*)'read_cmc_flx:  -- next input at ireadf = ', ireadf
(  674)         
(  675)         ! Read data
(  676)  !!Adjust using frame0 
(  677)         ! read temp at 2m (in K)
(  678)         call ncdf_read(clname_n, 'air', tatm, -(nflx1+frame0), status)
(  679)         if(status /= 0) then
(  680)            if(lwp) then
(  681)               write(numout,*)'read_cmc_flx: Prob reading tair file, status = ', status
(  682)               !               nstop = nstop + 1
(  683)            endif
(  684)         endif






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 13

(  685)         ! conversion of temperature Kelvin --> Celsius  [rt0=273.15]
(  686)         tatm(:,:) = ( tatm(:,:) - rt0 ) 
(  687)  
(  688)  !A.D: read date from tair file
(  689)         call ncdf_readdate(clname_n, 'real_time', cmc_date, (nflx1+frame0), status)        
(  690)         read(cmc_date(12:13),*) jhour                
(  691)            if(lwp) then
(  692)               write(numout,*)'read_cmc_flx: cmc_date = ', cmc_date,(nflx1+frame0)              
(  693)               write(numout,*)'read_cmc_flx: cmc_date (hour) = ', cmc_date(12:13)           
(  694)               write(numout,*)'read_cmc_flx: julian hour = ', jhour
(  695)            endif
(  696) 
(  697)         ! read humidity
(  698)         call ncdf_read(clname_c, 'socliohu', hatm, -(nflx1+frame0), status)
(  699)         if(status /= 0) then
(  700)            if(lwp) then
(  701)               write(numout,*)'read_cmc_flx: Prob reading hum_cloud file, status = ', status
(  702)               !               nstop = nstop + 1
(  703)            endif
(  704)         endif
(  705)         
(  706)         ! read clouds
(  707)         call ncdf_read(clname_c, 'socliocl', catm, -(nflx1+frame0), status)
(  708)         if(status /= 0) then
(  709)            if(lwp) then
(  710)               write(numout,*)'read_cmc_flx: Prob reading hum_cloud file, status = ', status
(  711)               !               nstop = nstop + 1
(  712)            endif
(  713)         endif
(  714)         
(  715)  
(  716)         call flush(numout)
(  717)           
(  718)      ENDIF
(  719)  
(  720)  
(  721)  !!DB: basic structure to do time interpolation every ndt_1h
(  722)  !!See tau_forced_cmc for more details
(  723)  !!NB: would have to have input the nflx2 fields as well (not done above)
(  724)  !!DB Use ndt_3h & ndt_1h
(  725)  !!Want a formula that also accounts for nature of ndt_* variables which may not
(  726)  !!be integers (depending on rdt), and still gives a consistent zxy upon restart
(  727)  !!NB: formula below is incorrect if perpetual_forcing is on so avoid this possibility 
(  728)  !      if(mod(kt-nit000,ndt_1h) == 0  .AND.  perpetual_forcing == 0) then
(  729)  !         zxy = (float(kt-nit000) - (nflx1-1)*rdt_3h)/rdt_3h
(  730)  !         zxy = min(zxy,1.0)  !to be safe
(  731)  !         taux(:,:) = (1.-zxy) * taux_dta(:,:,1) + zxy * taux_dta(:,:,2)
(  732)  !         tauy(:,:) = (1.-zxy) * tauy_dta(:,:,1) + zxy * tauy_dta(:,:,2)
(  733)  !      endif
(  734)  
(  735)  !       CALL blk(kt)
(  736)        
(  737)        CALL FLUSH(numout)
(  738)        
(  739)    END SUBROUTINE read_cmc_flx
(  740)       
(  741)    
(  742) #else






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 14

(  743)    !!----------------------------------------------------------------------
(  744)    !!   Default option :           Empty module                     NO bulk
(  745)    !!----------------------------------------------------------------------
(  746) CONTAINS
(  747)    SUBROUTINE flx_blk_solar            ! Empty routine
(  748)    END SUBROUTINE flx_blk_solar
(  749) #endif
(  750) 
(  751)    !!======================================================================
(  752) END MODULE flxblk_solar






















































PGF90 (Version     10.5)          09/23/2010  09:49:23      page 15

(    1) # 1 "flxblk_solar.F90"
(    1) !!DB 2008.12.04
(    2) !!Modified routine to remove possible conflicts with lim-calls and
(    3) !!for easier use in BGCM_01
(    4) !!-1- changed qsr_oce? to ad_qsr_oce?
(    5) !!-2- the variable ad_qsr_oce() now contains a true daily ave
(    6) !!    (for BGCM_01 use)
(    7) !!-3- #ifdef key_lim_ice ==> ad_qsr_oce(:,:)=frld(:,:)*ad_qsr_oce(:,:)
(    8) !!    where frld is the fraction of open water. Note that this is not
(    9) !!    meant to be a perfect accounting of ice coverage, but rather
(   10) !!    something likely better than nothing. It will set qsr=0 for
(   11) !!    completely ice-covered cells. 
(   12) !!DB 2009.08.17 
(   13) !!Re -2- above: I adjusted the already-daily-averaged insolation value to be
(   14) !!daily-averaged which reduces the insolation. This turns out to be not
(   15) !!so bad as it looks like an average cloudiness correction ===> Keep as is
(   16) !!
(   17) 
(   18) MODULE flxblk_solar
(   19)    !!======================================================================
(   20)    !!                       ***  MODULE  flxblk_solar  ***
(   21)    !!
(   22)    !!=================== SUMMARY =========================================
(   23)    !!Purpose:         Stand alone Short wave flux calcuation
(   24)    !!
(   25)    !!Inputs:         None / CMC humidity cloud cover and air temp
(   26)    !!
(   27)    !!Outputs:             qsr_oce,     &        !: solar flux over ocean (daily average)
(   28)    !!                     qsr_oce1              !: solar flux over ocean (instantaneous) &
(   29)    !!                                                available  only with cmc
(   30)    !! 
(   31)    !! Usage:   -compilie with key_bulk_solar. No other keys required
(   32)    !!          -SUBROUTINE flx_blk_solar( kt) must be called at each timestep
(   33)    !!          -Access qsr_oce and qsr_oce1 with "USE flxblk_solar"
(   34)    !!================= DESCRIPTION =======================================
(   35)    !!
(   36)    !!    SUBROUTINE flx_blk_solar( kt) computes shortwave radiation based
(   37)    !!     on the Zillman (1972). It has been stripped out of flxblk.F90
(   38)    !!     which uses this formula for the open ocean and something else for
(   39)    !!     ice (and as a result is tied to runing opa with lim and bulk
(   40)    !!     forcing).
(   41)    !!   
(   42)    !!    The calculation from flxblk.F90 computes a daily average (kept
(   43)    !!    here as qsr_oce) and is independent of local time.  This module
(   44)    !!    also offers an instantanous flux based on the "real_time"
(   45)    !!    varialble read from the cmc forcing.  read_cmc_flx opens the cmc
(   46)    !!    netcdf forcing files (based on flx_bulk_cmc.h90) and reads the
(   47)    !!    nearest 3h forcast).  For this part, the nit000 is assumed to
(   48)    !!    coincide with the first fields in the cmc files. The correct
(   49)    !!    yearday is provided through ndate0 in namelist.  The reading of
(   50)    !!    the "real_time" in the cmc files is accomplished with
(   51)    !!    ncdf_readdate, a new subroutine in lib_ncdf, designed to read a
(   52)    !!    string of length 19 characters from an array of unlimited
(   53)    !!    length. The only time information from the real_time variable that
(   54)    !!    is actually used, is jhour which is the local gmt time. jhour is
(   55)    !!    used to calculate the appropriate time offset based on the
(   56)    !!    longitude of each grid point.
(   57)    !!   






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 16

(   58)    !!    If do_read_cmc logical parameter is set to false. The cmc files
(   59)    !!    are not read. The Zillman formula is computed from constant values
(   60)    !!    for humidity, cloud cover and air temperate. These are hardwired
(   61)    !!    in the code (search for do_read_cmc).  In this scenario qsr_oce1
(   62)    !!    is meaningless and only qsr_oce should be used.
(   63)    !!   
(   64)    !!    Last Modify: 28Sep2008 by A. Drozdowski
(   65)    !!   =====================================================================
(   66) 
(   67) 
(   68) 
(   69) # 743
(  743)    !!----------------------------------------------------------------------
(  744)    !!   Default option :           Empty module                     NO bulk
(  745)    !!----------------------------------------------------------------------
(  746) CONTAINS
(  747)    SUBROUTINE flx_blk_solar            ! Empty routine
(  748)    END SUBROUTINE flx_blk_solar
(  749) 
(  750) # 751
(  751)    !!======================================================================
(  752) END MODULE flxblk_solar
