


PGF90 (Version     10.5)          09/23/2010  09:49:11      page 1

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

Filename: flxblk.F90

(    1) !!DB: original flxblk as of 2009.05.04-12:49
(    2) 
(    3) MODULE flxblk
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  flxblk  ***
(    6)    !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
(    7)    !!=====================================================================
(    8) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(    9)    !!----------------------------------------------------------------------
(   10)    !!   'key_flx_bulk_monthly'   or                            MONTHLY bulk
(   11)    !!   'key_flx_bulk_daily'                                     DAILY bulk
(   12)    !!----------------------------------------------------------------------
(   13)    !!   flx_blk        : thermohaline fluxes from bulk
(   14)    !!   flx_blk_declin : solar declinaison
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE cpl_oce         ! ???
(   20)    USE phycst          ! physical constants
(   21)    USE daymod
(   22)    USE blk_oce         ! bulk variables
(   23)    USE flx_oce         ! forcings variables
(   24)    USE ocfzpt          ! ???
(   25)    USE in_out_manager
(   26)    USE lbclnk
(   27)    USE albedo
(   28)    USE prtctl          ! Print control
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Accessibility
(   34)    PUBLIC flx_blk        ! routine called by flx.F90 
(   35) 
(   36)    !! * Module variables
(   37)    INTEGER, PARAMETER  ::   &
(   38)       jpintsr = 24          ! number of time step between sunrise and sunset
(   39)       !                     ! uses for heat flux computation
(   40)    LOGICAL ::   &
(   41)       lbulk_init = .TRUE.   ! flag, bulk initialization done or not)
(   42) 
(   43)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   44)       stauc            ,  &   ! cloud optical depth 
(   45)       sbudyko   
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 2

(   47)    !! * constants for bulk computation (flx_blk)
(   48)    REAL(wp), DIMENSION(19)  ::  &
(   49)       budyko                  ! BUDYKO's coefficient
(   50)    ! BUDYKO's coefficient (cloudiness effect on LW radiation):
(   51)    DATA budyko / 1.00, 0.98, 0.95, 0.92, 0.89, 0.86, 0.83, 0.80, 0.78, 0.75,  &
(   52)       &          0.72, 0.69, 0.67, 0.64, 0.61, 0.58, 0.56, 0.53, 0.50 /
(   53)    REAL(wp), DIMENSION(20)  :: &
(   54)       tauco                  ! cloud optical depth coefficient
(   55)    ! Cloud optical depth coefficient
(   56)    DATA tauco / 6.6, 6.6, 7.0, 7.2, 7.1, 6.8, 6.5, 6.6, 7.1, 7.6,   &
(   57)       &         6.6, 6.1, 5.6, 5.5, 5.8, 5.8, 5.6, 5.6, 5.6, 5.6 /
(   58)    REAL(wp)  ::            &  ! constant values
(   59)       zeps    = 1.e-20  ,  &
(   60)       zeps0   = 1.e-13  ,  &
(   61)       zeps1   = 1.e-06  ,  &
(   62)       zzero   = 0.e0    ,  &
(   63)       zone    = 1.0
(   64) 
(   65)    !! * constants for solar declinaison computation (flx_blk_declin)
(   66)    REAL(wp) ::                &
(   67)       a0  =  0.39507671   ,   &  ! coefficients
(   68)       a1  = 22.85684301   ,   &
(   69)       a2  = -0.38637317   ,   &
(   70)       a3  =  0.15096535   ,   &
(   71)       a4  = -0.00961411   ,   &
(   72)       b1  = -4.29692073   ,   &
(   73)       b2  =  0.05702074   ,   &
(   74)       b3  = -0.09028607   ,   &
(   75)       b4  =  0.00592797
(   76)    !!----------------------------------------------------------------------
(   77)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   78)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxblk.F90,v 1.8 2005/09/02 15:45:31 opalod Exp $ 
(   79)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   80)    !!----------------------------------------------------------------------
(   81) 
(   82) CONTAINS
(   83) 
(   84)    SUBROUTINE flx_blk( psst )
(   85)       !!---------------------------------------------------------------------------
(   86)       !!                     ***  ROUTINE flx_blk  ***
(   87)       !!                 
(   88)       !!  ** Purpose :   Computation of the heat fluxes at ocean and snow/ice
(   89)       !!       surface the solar heat at ocean and snow/ice surfaces and the 
(   90)       !!       sensitivity of total heat fluxes to the SST variations
(   91)       !!         
(   92)       !!  ** Method  :   The flux of heat at the ice and ocean surfaces are derived
(   93)       !!       from semi-empirical ( or bulk ) formulae which relate the flux to 
(   94)       !!       the properties of the surface and of the lower atmosphere. Here, we
(   95)       !!       follow the work of Oberhuber, 1988   
(   96)       !!
(   97)       !!  ** Action  :   call flx_blk_albedo to compute ocean and ice albedo 
(   98)       !!          computation of snow precipitation
(   99)       !!          computation of solar flux at the ocean and ice surfaces
(  100)       !!          computation of the long-wave radiation for the ocean and sea/ice
(  101)       !!          computation of turbulent heat fluxes over water and ice
(  102)       !!          computation of evaporation over water
(  103)       !!          computation of total heat fluxes sensitivity over ice (dQ/dT)
(  104)       !!          computation of latent heat flux sensitivity over ice (dQla/dT)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 3

(  105)       !!
(  106)       !! History :
(  107)       !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
(  108)       !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
(  109)       !!----------------------------------------------------------------------
(  110)       !! * Arguments
(  111)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  112)          &                          psst      ! Sea Surface Temperature 
(  113) 
(  114)       !! * Local variables
(  115)       INTEGER  ::             &
(  116)          ji, jj, jt        ,  &  ! dummy loop indices
(  117)          indaet            ,  &  !  = -1, 0, 1 for odd, normal and leap years resp.
(  118)          iday              ,  &  ! integer part of day
(  119)          indxb             ,  &  ! index for budyko coefficient
(  120)          indxc                   ! index for cloud depth coefficient
(  121) 
(  122)       REAL(wp)  ::            & 
(  123)          zalat , zclat     ,  &  ! latitude in degrees 
(  124)          zmt1, zmt2, zmt3  ,  &  ! tempory air temperatures variables
(  125)          ztatm3, ztatm4    ,  &  ! power 3 and 4 of air temperature
(  126)          z4tatm3           ,  &  ! 4 * ztatm3
(  127)          zcmue             ,  &  ! cosine of local solar altitude
(  128)          zcmue2            ,  &  ! root of zcmue1 
(  129)          zscmue            ,  &  ! square-root of zcmue1 
(  130)          zpcmue            ,  &  ! zcmue1**1.4
(  131)          zdecl             ,  &  ! solar declination
(  132)          zsdecl , zcdecl   ,  &  ! sine and cosine of solar declination 
(  133)          zalbo             ,  &  ! albedo of sea-water
(  134)          zalbi             ,  &  ! albedo of ice
(  135)          ztamr             ,  &  ! air temperature minus triple point of water (rtt)
(  136)          ztaevbk           ,  &  ! part of net longwave radiation
(  137)          zevi , zevo       ,  &  ! vapour pressure of ice and ocean 
(  138)          zind1,zind2,zind3 ,  &  ! switch for testing the values of air temperature
(  139)          zinda             ,  &  ! switch for testing the values of sea ice cover
(  140)          zpis2             ,  &  ! pi / 2
(  141)          z2pi                    ! 2 * pi 
(  142) 
(  143)       REAL(wp)  ::            & 
(  144)          zxday             ,  &  ! day of year
(  145)          zdist             ,  &  ! distance between the sun and the earth during the year
(  146)          zdaycor           ,  &  ! corr. factor to take into account the variation of 
(  147)          !                       ! zday when calc. the solar rad.    
(  148)          zesi, zeso        ,  &  ! vapour pressure of ice and ocean at saturation
(  149)          zesi2             ,  &  ! root of zesi 
(  150)          zqsato            ,  &  ! humidity close to the ocean surface (at saturation)   
(  151)          zqsati            ,  &  ! humidity close to the ice surface (at saturation) 
(  152)          zqsati2           ,  &  ! root of  zqsati 
(  153)          zdesidt           ,  &  ! derivative of zesi, function of ice temperature
(  154)          zdteta            ,  &  ! diff. betw. sst and air temperature
(  155)          zdeltaq           ,  &  ! diff. betw. spec. hum. and hum. close to the surface
(  156)          ztvmoy, zobouks   ,  &  ! tempory scalars
(  157)          zpsims, zpsihs, zpsils, zobouku, zxins, zpsimu ,  & 
(  158)          zpsihu, zpsilu, zstab,zpsim, zpsih, zpsil      ,  & 
(  159)          zvatmg, zcmn, zchn, zcln, zcmcmn, zdenum       ,  & 
(  160)          zdtetar, ztvmoyr, zlxins, zcmn2, zchcm, zclcm , zcoef
(  161) 
(  162)       REAL(wp)  ::            & 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 4

(  163)          zrhova            ,  &  ! air density per wind speed
(  164)          zcsho , zcleo     ,  &  ! transfer coefficient over ocean
(  165)          zcshi , zclei     ,  &  ! transfer coefficient over ice-free
(  166)          zrhovacleo        ,  &  ! air density per wind speed per transfer coef.
(  167)          zrhovacsho, zrhovaclei, zrhovacshi, & 
(  168)          ztice3            ,  &  ! power 3 of ice temperature
(  169)          zticemb, zticemb2 ,  &  ! tempory air temperatures variables
(  170)          zdqlw_ice         ,  &  ! sensitivity of long-wave flux over ice
(  171)          zdqsb_ice         ,  &  ! sensitivity of sensible heat flux over ice
(  172)          zdqla_ice         ,  &  ! sensitivity of latent heat flux over ice
(  173)          zdl, zdr                ! fractionnal part of latitude
(  174)       REAL(wp), DIMENSION(jpi,jpj) :: & 
(  175)          zpatm            ,  &   ! atmospheric pressure
(  176)          zqatm            ,  &   ! specific humidity
(  177)          zes              ,  &   ! vapour pressure at saturation
(  178)          zev, zevsqr      ,  &   ! vapour pressure and his square-root
(  179)          zrhoa            ,  &   ! air density
(  180)          ztatm            ,  &   ! air temperature in Kelvins
(  181)          zfrld            ,  &   ! fraction of sea ice cover 
(  182)          zcatm1           ,  &   ! fraction of cloud
(  183)          zcldeff                 ! correction factor to account cloud effect
(  184)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  185)          zalbocsd         ,  &   ! albedo of ocean
(  186)          zalboos          ,  &   ! albedo of ocean under overcast sky
(  187)          zalbics          ,  &   ! albedo of ice under clear sky
(  188)          zalbios          ,  &   ! albedo of ice under overcast sky
(  189)          zalbomu          ,  &   ! albedo of ocean when zcmue is 0.4
(  190)          zqsro            ,  &   ! solar radiation over ocean
(  191)          zqsrics          ,  &   ! solar radiation over ice under clear sky
(  192)          zqsrios          ,  &   ! solar radiation over ice under overcast sky
(  193)          zcldcor          ,  &   ! cloud correction
(  194)          zlsrise, zlsset  ,  &   ! sunrise and sunset
(  195)          zlmunoon         ,  &   ! local noon solar altitude
(  196)          zdlha            ,  &   ! length of the ninstr segments of the solar day
(  197)          zps              ,  &   ! sine of latitude per sine of solar decli.
(  198)          zpc                     ! cosine of latitude per cosine of solar decli. 
(  199) 
(  200)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  201)          zqlw_oce         ,  &   ! long-wave heat flux over ocean
(  202)          zqlw_ice         ,  &   ! long-wave heat flux over ice
(  203)          zqla_oce         ,  &   ! latent heat flux over ocean
(  204)          zqla_ice         ,  &   ! latent heat flux over ice
(  205)          zqsb_oce         ,  &   ! sensible heat flux over ocean
(  206)          zqsb_ice                ! sensible heat flux over ice
(  207)  
(  208)       REAL(wp), DIMENSION(jpi,jpj,jpintsr) ::    &
(  209)          zlha             ,  &   ! local hour angle
(  210)          zalbocs          ,  &   ! tempory var. of ocean albedo under clear sky
(  211)          zsqsro           ,  &   ! tempory var. of solar rad. over ocean 
(  212)          zsqsrics         ,  &   ! temp. var. of solar rad. over ice under clear sky
(  213)          zsqsrios                ! temp. var. of solar rad. over ice under overcast sky
(  214)       !!---------------------------------------------------------------------
(  215) 
(  216)       !---------------------
(  217)       !   Initilization    !
(  218)       !---------------------
(  219) #if ! defined key_ice_lim
(  220)       tn_ice(:,:) = psst(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 5

(  221) #endif
(  222) 
(  223)       !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
(  224)       !  and the correction factor for taking into account  the effect of clouds 
(  225)       !------------------------------------------------------
(  226)       IF( lbulk_init ) THEN
(  227)          DO jj = 1, jpj  
(  228)             DO ji = 1 , jpi
(  229)                zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
(  230)                zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
(  231)                indxb          = 1 + INT( zalat ) 
(  232)                !  correction factor to account for the effect of clouds 
(  233)                sbudyko(ji,jj) = budyko(indxb)  
(  234)                indxc          = 1 + INT( zclat )  
(  235)                zdl            = zclat - INT( zclat ) 
(  236)                zdr            = 1.0 - zdl
(  237)                stauc(ji,jj)   = zdr * tauco( indxc ) + zdl * tauco( indxc + 1 ) 
(  238)             END DO
(  239)          END DO
(  240)          IF( nleapy == 1 ) THEN
(  241)             yearday = 366.e0
(  242)          ELSE IF( nleapy == 0 ) THEN
(  243)             yearday = 365.e0
(  244)          ELSEIF( nleapy == 30) THEN
(  245)             yearday = 360.e0
(  246)          ENDIF
(  247)          lbulk_init = .FALSE.
(  248)       ENDIF
(  249) 
(  250)       zqlw_oce(:,:) = 0.e0
(  251)       zqla_oce(:,:) = 0.e0
(  252)       zqsb_oce(:,:) = 0.e0
(  253)       zqlw_ice(:,:) = 0.e0
(  254)       zqla_ice(:,:) = 0.e0
(  255)       zqsb_ice(:,:) = 0.e0
(  256) 
(  257)       zpis2       = rpi / 2.
(  258)       z2pi        = 2. * rpi
(  259) 
(  260)  !CDIR NOVERRCHK
(  261)       DO jj = 1, jpj
(  262)  !CDIR NOVERRCHK
(  263)          DO ji = 1, jpi
(  264) 
(  265)             ztatm (ji,jj) = 273.15 + tatm  (ji,jj)  !  air temperature in Kelvins 
(  266)             zcatm1(ji,jj) = 1.0    - catm  (ji,jj)  !  fractional cloud cover
(  267)             zfrld (ji,jj) = 1.0    - freeze(ji,jj)  !  fractional sea ice cover
(  268)             zpatm(ji,jj)  = 101000.               !  pressure 
(  269)       
(  270)             !  Computation of air density, obtained from the equation of state for dry air. 
(  271)             zrhoa(ji,jj) = zpatm(ji,jj) / ( 287.04 * ztatm(ji,jj) )
(  272)       
(  273)             !  zes : Saturation water vapour
(  274)             ztamr = ztatm(ji,jj) - rtt
(  275)             zmt1  = SIGN( 17.269, ztamr )
(  276)             zmt2  = SIGN( 21.875, ztamr )
(  277)             zmt3  = SIGN( 28.200, -ztamr )
(  278)             zes(ji,jj) = 611.0 * EXP (  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 6

(  279)                &                      / ( ztatm(ji,jj) - 35.86  + MAX( zzero, zmt3 ) ) )
(  280) 
(  281)             !  zev : vapour pressure  (hatm is relative humidity)  
(  282)             zev(ji,jj)   = hatm(ji,jj) * zes(ji,jj) 
(  283)             !  square-root of vapour pressure
(  284) !CDIR NOVERRCHK
(  285)             zevsqr(ji,jj) = SQRT( zev(ji,jj) * 0.01 )
(  286)             !  zqapb  : specific humidity 
(  287)             zqatm(ji,jj) = 0.622 * zev(ji,jj) / ( zpatm(ji,jj) - 0.378 * zev(ji,jj) )
(  288) 
(  289) 
(  290)             !----------------------------------------------------
(  291)             !   Computation of snow precipitation (Ledley, 1985) |
(  292)             !----------------------------------------------------
(  293) 
(  294)             zmt1  =   253.0 - ztatm(ji,jj)
(  295)             zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0 
(  296)             zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0
(  297)             zind1 = MAX( zzero, SIGN( zone, zmt1 ) )
(  298)             zind2 = MAX( zzero, SIGN( zone, zmt2 ) )
(  299)             zind3 = MAX( zzero, SIGN( zone, zmt3 ) )
(  300)             ! total precipitation
(  301)             tprecip(ji,jj) = watm(ji,jj)
(  302)             ! solid  (snow) precipitation
(  303)             sprecip(ji,jj) = tprecip(ji,jj) *       &
(  304)                &             (           zind1      &
(  305)                &               + ( 1.0 - zind1 ) * ( zind2 * ( 0.5 + zmt2 ) + ( 1.0 - zind2 ) *  zind3 * zmt3 ) ) 
(  306)          END DO
(  307)       END DO
(  308) 
(  309)       !----------------------------------------------------------
(  310)       !   Computation of albedo (need to calculates heat fluxes)|
(  311)       !-----------------------------------------------------------
(  312)       
(  313)       CALL flx_blk_albedo( zalbios, zalboos, zalbics, zalbomu )
(  314) 
(  315)       !-------------------------------------
(  316)       !   Computation of solar irradiance. |
(  317)       !----------------------------------------
(  318)       indaet   = 1  
(  319)       !  compution of the day of the year at which the fluxes have to be calculate 
(  320)       !--The date corresponds to the middle of the time step.
(  321)       zxday=nday_year + rdtbs2/rday
(  322) 
(  323)       iday   = INT( zxday )
(  324) 
(  325)       IF(ln_ctl) CALL prt_ctl_info('declin : iday ', ivar1=iday, clinfo2=' nfbulk= ', ivar2=nfbulk)
(  326) 
(  327)       !   computation of the solar declination, his sine and his cosine
(  328)       CALL flx_blk_declin( indaet, iday, zdecl )
(  329)       
(  330)       zdecl    = zdecl * rad
(  331)       zsdecl   = SIN( zdecl )
(  332)       zcdecl   = COS( zdecl )
(  333)       
(  334)       !  correction factor added for computation of shortwave flux to take into account the variation of
(  335)       !  the distance between the sun and the earth during the year (Oberhuber 1988)
(  336)       zdist    = zxday * z2pi / yearday






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 7

(  337)       zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )
(  338) 
(  339) !CDIR NOVERRCHK
(  340)       DO jj = 1, jpj
(  341) !CDIR NOVERRCHK
(  342)          DO ji = 1, jpi
(  343)             !  product of sine of latitude and sine of solar declination
(  344)             zps     (ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
(  345)             !  product of cosine of latitude and cosine of solar declination
(  346)             zpc     (ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
(  347)             !  computation of the both local time of sunrise and sunset
(  348)             zlsrise (ji,jj) = ACOS( - SIGN( zone, zps(ji,jj) ) * MIN( zone, SIGN( zone, zps(ji,jj) )  &
(  349)                &                     * ( zps(ji,jj) / zpc(ji,jj) ) ) ) 
(  350)             zlsset  (ji,jj) = - zlsrise(ji,jj)
(  351)             !  dividing the solar day into jpintsr segments of length zdlha
(  352)             zdlha   (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jpintsr )
(  353)             !  computation of the local noon solar altitude
(  354)             zlmunoon(ji,jj) = ASIN ( ( zps(ji,jj) + zpc(ji,jj) ) ) / rad
(  355)             
(  356)             !  cloud correction taken from Reed (1977) (imposed lower than 1)
(  357)             zcldcor (ji,jj) = MIN( zone, ( 1.e0 - 0.62 * catm(ji,jj) + 0.0019 * zlmunoon(ji,jj) ) )
(  358)          END DO
(  359)       END DO
(  360) 
(  361)          !  Computation of solar heat flux at each time of the day between sunrise and sunset. 
(  362)          !  We do this to a better optimisation of the code 
(  363)          !------------------------------------------------------       
(  364) 
(  365) !CDIR NOVERRCHK   
(  366)       DO jt = 1, jpintsr   
(  367)          zcoef = FLOAT( jt ) - 0.5
(  368) !CDIR NOVERRCHK     
(  369)          DO jj = 1, jpj
(  370) !CDIR NOVERRCHK
(  371)             DO ji = 1, jpi
(  372)                !  local hour angle
(  373)                zlha (ji,jj,jt) = COS ( zlsrise(ji,jj) - zcoef * zdlha(ji,jj) )
(  374) 
(  375)                ! cosine of local solar altitude
(  376)                zcmue              = MAX ( zzero ,   zps(ji,jj) + zpc(ji,jj) * zlha (ji,jj,jt)  )
(  377)                zcmue2             = 1368.0 * zcmue * zcmue
(  378)                zscmue             = SQRT ( zcmue )
(  379)                zpcmue             = zcmue**1.4
(  380)                ! computation of sea-water albedo (Payne, 1972)
(  381)                zalbocs(ji,jj,jt)  = 0.05 / ( 1.1 * zpcmue + 0.15 )
(  382)                zalbo              = zcatm1(ji,jj) * zalbocs(ji,jj,jt) + catm(ji,jj) * zalboos(ji,jj)
(  383)                ! solar heat flux absorbed at ocean surfaces (Zillman, 1972)
(  384)                zevo               = zev(ji,jj) * 1.0e-05
(  385)                zsqsro(ji,jj,jt)   =  ( 1.0 - zalbo ) * zdlha(ji,jj) * zcmue2                &
(  386)                                    / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
(  387)                !  solar heat flux absorbed at sea/ice surfaces 
(  388)                !  Formulation of Shine and Crane, 1984 adapted for high albedo surfaces 
(  389) 
(  390)                !  For clear sky        
(  391)                zevi               = zevo
(  392)                zalbi              = zalbics(ji,jj)
(  393)                zsqsrics(ji,jj,jt) =  ( 1.0 - zalbi ) * zdlha(ji,jj) * zcmue2                &
(  394)                   &                / ( ( 1.0 + zcmue ) * zevi + 1.2 * zcmue + 0.0455 )






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 8

(  395) 
(  396)                ! For overcast sky
(  397)                zalbi              = zalbios(ji,jj)
(  398)                zsqsrios(ji,jj,jt) = zdlha(ji,jj) *                                                           &
(  399)                   &                 ( ( 53.5 + 1274.5 * zcmue      ) *  zscmue * ( 1.0 - 0.996  * zalbi ) )  &
(  400)                   &                 / (  1.0 + 0.139  * stauc(ji,jj) *           ( 1.0 - 0.9435 * zalbi ) )
(  401)             END DO
(  402)          END DO
(  403)       END DO
(  404) 
(  405) 
(  406)       !  Computation of daily (between sunrise and sunset) solar heat flux absorbed 
(  407)       !  at the ocean and snow/ice surfaces.
(  408)       !--------------------------------------------------------------------
(  409) 
(  410)       zalbocsd(:,:) = 0.e0
(  411)       zqsro   (:,:) = 0.e0
(  412)       zqsrics (:,:) = 0.e0
(  413)       zqsrios (:,:) = 0.e0
(  414) 
(  415)       DO jt = 1, jpintsr 
(  416) #   if defined key_vectopt_loop && ! defined key_autotasking
(  417)          DO ji = 1, jpij  
(  418)             zalbocsd(ji,1) = zalbocsd(ji,1) + zdlha   (ji,1) * zalbocs(ji,1,jt)   &
(  419)                &                                             / MAX( 2.0 * zlsrise(ji,1) , zeps0 )
(  420)             zqsro   (ji,1) = zqsro   (ji,1) + zsqsro  (ji,1,jt)
(  421)             zqsrics (ji,1) = zqsrics (ji,1) + zsqsrics(ji,1,jt)
(  422)             zqsrios (ji,1) = zqsrios (ji,1) + zsqsrios(ji,1,jt)
(  423)          END DO
(  424) #  else
(  425)          DO jj = 1, jpj
(  426)             DO ji = 1, jpi  
(  427)                zalbocsd(ji,jj) = zalbocsd(ji,jj) + zdlha(ji,jj) * zalbocs(ji,jj,jt)   &
(  428)                   &                                              / MAX( 2.0 * zlsrise(ji,jj) , zeps0 )
(  429)                zqsro  (ji,jj)  = zqsro   (ji,jj) + zsqsro  (ji,jj,jt)
(  430)                zqsrics(ji,jj)  = zqsrics (ji,jj) + zsqsrics(ji,jj,jt)
(  431)                zqsrios(ji,jj)  = zqsrios (ji,jj) + zsqsrios(ji,jj,jt)
(  432)             END DO
(  433)          END DO
(  434) #  endif
(  435)       END DO
(  436) 
(  437)       DO jj = 1, jpj
(  438)          DO ji = 1, jpi 
(  439) 
(  440)             !------------------------------------------- 
(  441)             !  Computation of shortwave radiation.
(  442)             !-------------------------------------------
(  443) 
(  444)             ! the solar heat flux absorbed at ocean and snow/ice surfaces
(  445)             !------------------------------------------------------------
(  446) 
(  447)             ! For ocean
(  448)             qsr_oce(ji,jj) = srgamma * zcldcor(ji,jj) * zqsro(ji,jj) / z2pi
(  449)             zinda          = SIGN( zone , -( -0.5 - zfrld(ji,jj) ) )
(  450)             zinda          = 1.0 - MAX( zzero , zinda )
(  451)             qsr_oce(ji,jj) = ( 1.- zinda ) * qsr_oce(ji,jj)
(  452) 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 9

(  453)             ! For snow/ice 
(  454)             qsr_ice(ji,jj) = ( zcatm1(ji,jj) * zqsrics(ji,jj) + catm(ji,jj) * zqsrios(ji,jj) ) / z2pi
(  455) 
(  456) 
(  457)             ! Taking into account the ellipsity of the earth orbit
(  458)             !-----------------------------------------------------
(  459) 
(  460)             qsr_ice(ji,jj) = qsr_ice(ji,jj) * zdaycor
(  461)             qsr_oce(ji,jj) = qsr_oce(ji,jj) * zdaycor
(  462) 
(  463)             !  fraction of net shortwave radiation which is not absorbed in the 
(  464)             !  thin surface layer and penetrates inside the ice cover 
(  465)             !  ( Maykut and Untersteiner, 1971 ; Elbert anbd Curry, 1993 )
(  466)             !------------------------------------------------------------------
(  467) 
(  468)             fr1_i0(ji,jj) = 0.18  * zcatm1(ji,jj) + 0.35 * catm(ji,jj) 
(  469)             fr2_i0(ji,jj) = 0.82  * zcatm1(ji,jj) + 0.65 * catm(ji,jj)
(  470) 
(  471)             !---------------------------------------------------------------------------
(  472)             !   Computation of long-wave radiation  ( Berliand 1952 ; all latitudes )
(  473)             !---------------------------------------------------------------------------
(  474) 
(  475)             ! tempory variables
(  476)             ztatm3         = ztatm(ji,jj) * ztatm(ji,jj) * ztatm(ji,jj)
(  477)             ztatm4         = ztatm3 * ztatm(ji,jj)
(  478)             z4tatm3        = 4. * ztatm3
(  479)             zcldeff(ji,jj) = 1.0 - sbudyko(ji,jj) * catm(ji,jj) * catm(ji,jj)    
(  480)             ztaevbk        = ztatm4 * zcldeff(ji,jj) * ( 0.39 - 0.05 * zevsqr(ji,jj) ) 
(  481) 
(  482)             !  Long-Wave for Ice
(  483)             !----------------------
(  484)             zqlw_ice(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( tn_ice(ji,jj) - ztatm(ji,jj) ) ) 
(  485) 
(  486)             !  Long-Wave for Ocean
(  487)             !-----------------------
(  488)             zqlw_oce(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( psst  (ji,jj) - ztatm(ji,jj) ) ) 
(  489) 
(  490)          END DO
(  491)       END DO
(  492) 
(  493)       !----------------------------------------
(  494)       !  Computation of turbulent heat fluxes  ( Latent and sensible ) 
(  495)       !----------------------------------------        
(  496)       !CDIR NOVERRCHK
(  497)       DO jj = 2 , jpjm1
(  498) !ib   DO jj = 1 , jpj
(  499)          !CDIR NOVERRCHK
(  500)          DO  ji = 1, jpi
(  501) 
(  502)             !  Turbulent heat fluxes over water
(  503)             !----------------------------------
(  504) 
(  505)             ! zeso     : vapour pressure at saturation of ocean
(  506)             ! zqsato   : humidity close to the ocean surface (at saturation)
(  507)             zeso          =  611.0 * EXP ( 17.2693884 * ( psst(ji,jj) - rtt ) * tmask(ji,jj,1) / ( psst(ji,jj) - 35.86 ) )
(  508)             zqsato        = ( 0.622 * zeso ) / ( zpatm(ji,jj) - 0.378 * zeso )
(  509) 
(  510)             !  Drag coefficients from Large and Pond (1981,1982)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 10

(  511)             !---------------------------------------------------
(  512)     
(  513)             !  Stability parameters
(  514)             zdteta         = psst(ji,jj) - ztatm(ji,jj)
(  515)             zdeltaq        = zqatm(ji,jj) - zqsato
(  516)             ztvmoy         = ztatm(ji,jj) * ( 1. + 2.2e-3 * ztatm(ji,jj) * zqatm(ji,jj) )
(  517)             zdenum         = MAX( vatm(ji,jj) * vatm(ji,jj) * ztvmoy, zeps )
(  518) !i
(  519) !i          if( zdenum == 0.e0 ) then
(  520) !i               write(numout,*) 'flxblk  zdenum=0 ', ji,jj
(  521) !i               zdenum = zeps
(  522) !i          endif
(  523) !i
(  524)             zdtetar        = zdteta / zdenum
(  525)             ztvmoyr        = ztvmoy * ztvmoy * zdeltaq / zdenum
(  526)             
(  527)             ! For stable atmospheric conditions
(  528)             zobouks        = -70.0 * 10. * ( zdtetar + 3.2e-3 * ztvmoyr )
(  529)             zobouks        = MAX( zzero , zobouks )
(  530)             zpsims         = -7.0 * zobouks
(  531)             zpsihs         =  zpsims
(  532)             zpsils         =  zpsims
(  533) 
(  534)             !  For unstable atmospheric conditions
(  535)             zobouku        = -100.0 * 10.0 * ( zdtetar + 2.2e-3 * ztvmoyr )
(  536)             zobouku        = MIN( zzero , zobouku )
(  537)             zxins          = ( 1. - 16. * zobouku )**0.25
(  538)             zlxins         = LOG( ( 1. + zxins * zxins ) / 2. )
(  539)             zpsimu         = 2. * LOG( ( 1 + zxins ) / 2. )  + zlxins - 2. * ATAN( zxins ) + zpis2
(  540)             zpsihu         = 2. * zlxins
(  541)             zpsilu         = zpsihu
(  542) 
(  543)             ! computation of intermediate values
(  544)             zstab          = MAX( zzero , SIGN( zone , zdteta ) )
(  545)             zpsim          = zstab * zpsimu + (1.0 - zstab ) * zpsims
(  546)             zpsih          = zstab * zpsihu + (1.0 - zstab ) * zpsihs
(  547)             zpsil          = zpsih
(  548)             
(  549)             zvatmg         = MAX( 0.032 * 1.5e-3 * vatm(ji,jj) * vatm(ji,jj) / grav, zeps )
(  550) !i
(  551) !!          if( zvatmg == 0.e0 ) then
(  552) !!               write(numout,*) 'flxblk  zvatmg=0 ', ji,jj
(  553) !!               zvatmg = zeps
(  554) !!          endif
(  555) !i
(  556) 
(  557)             zcmn           = vkarmn / LOG ( 10. / zvatmg )
(  558)             zcmn2          = zcmn * zcmn
(  559)             zchn           = 0.0327 * zcmn
(  560)             zcln           = 0.0346 * zcmn
(  561)             zcmcmn         = 1 / ( 1 - zcmn * zpsim / vkarmn )
(  562)             zchcm          = zcmcmn / ( 1 - zchn * zpsih / ( vkarmn * zcmn ) )
(  563)             zclcm          = zchcm
(  564) 
(  565) 
(  566)             !  Transfer cofficient zcsho and zcleo over ocean according to Large and Pond (1981,1982)
(  567)             !-------------------------------------------------------------- 
(  568)             zcsho          = zchn * zchcm






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 11

(  569)             zcleo          = zcln * zclcm 
(  570) 
(  571) 
(  572)             !   Computation of sensible and latent fluxes over Ocean 
(  573)             !----------------------------------------------------------------
(  574) 
(  575)             !  computation of intermediate values
(  576)             zrhova         = zrhoa(ji,jj) * vatm(ji,jj)
(  577)             zrhovacsho     = zrhova * zcsho
(  578)             zrhovacleo     = zrhova * zcleo
(  579) 
(  580)             ! sensible heat flux
(  581)             zqsb_oce(ji,jj) = zrhovacsho * 1004.0  * ( psst(ji,jj) - ztatm(ji,jj) )  
(  582)          
(  583)             !  latent heat flux 
(  584)             zqla_oce(ji,jj) = MAX(0.e0, zrhovacleo * 2.5e+06 * ( zqsato      - zqatm(ji,jj) ) )
(  585)                
(  586)             !  Calculate evaporation over water. (kg/m2/s)
(  587)             !-------------------------------------------------
(  588)             evap(ji,jj)    = zqla_oce(ji,jj) / cevap
(  589)                
(  590)                
(  591)             !  Turbulent heat fluxes over snow/ice.
(  592)             !--------------------------------------------------
(  593)             
(  594)             !  zesi     : vapour pressure at saturation of ice
(  595)             !  zqsati   : humidity close to the ice surface (at saturation)
(  596)             zesi           =  611.0 * EXP ( 21.8745587 * tmask(ji,jj,1)   &   ! tmask needed to avoid overflow in the expone
(  597)                &                                       * ( tn_ice(ji,jj) - rtt ) / ( tn_ice(ji,jj) - 7.66 ) )
(  598)             zqsati         = ( 0.622 * zesi ) / ( zpatm(ji,jj) - 0.378 * zesi )
(  599)                
(  600)             !  computation of intermediate values
(  601)             zticemb        = ( tn_ice(ji,jj) - 7.66 )
(  602)             zticemb2       = zticemb * zticemb  
(  603)             ztice3         = tn_ice(ji,jj) * tn_ice(ji,jj) * tn_ice(ji,jj)
(  604)             zqsati2        = zqsati * zqsati
(  605)             zesi2          = zesi * zesi
(  606)             zdesidt        = zesi * ( 9.5 * LOG( 10.0 ) * ( rtt - 7.66 )  / zticemb2 )
(  607)                
(  608)             !  Transfer cofficient zcshi and zclei over ice. Assumed to be constant Parkinson 1979 ; Maykut 1982
(  609)             !--------------------------------------------------------------------
(  610)             zcshi          = 1.75e-03
(  611)             zclei          = zcshi
(  612)                
(  613)             !  Computation of sensible and latent fluxes over ice
(  614)             !----------------------------------------------------------------
(  615)                
(  616)             !  computation of intermediate values
(  617)             zrhova          = zrhoa(ji,jj) * vatm(ji,jj)
(  618)             zrhovacshi      = zrhova * zcshi * 2.834e+06
(  619)             zrhovaclei      = zrhova * zclei * 1004.0
(  620)             
(  621)             !  sensible heat flux
(  622)             zqsb_ice(ji,jj) = zrhovaclei * ( tn_ice(ji,jj) - ztatm(ji,jj) )
(  623)             
(  624)             !  latent heat flux 
(  625)             zqla_ice(ji,jj) = zrhovacshi * ( zqsati        - zqatm(ji,jj) )
(  626)             qla_ice (ji,jj) = MAX(0.e0, zqla_ice(ji,jj) )






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 12

(  627)               
(  628)             !  Computation of sensitivity of non solar fluxes (dQ/dT)
(  629)             !---------------------------------------------------------------
(  630)                
(  631)             !  computation of long-wave, sensible and latent flux sensitivity
(  632)             zdqlw_ice       = 4.0 * emic * stefan * ztice3
(  633)             zdqsb_ice       = zrhovaclei
(  634)             zdqla_ice       = zrhovacshi * ( zdesidt * ( zqsati2 / zesi2 ) * ( zpatm(ji,jj) / 0.622 ) )         
(  635)             
(  636)             !  total non solar sensitivity
(  637)             dqns_ice(ji,jj) = -( zdqlw_ice + zdqsb_ice + zdqla_ice ) 
(  638)             
(  639)             ! latent flux sensitivity
(  640)             dqla_ice(ji,jj) = zdqla_ice
(  641)             
(  642)          END DO
(  643)       END DO
(  644) 
(  645)       ! total non solar heat flux over ice
(  646)       qnsr_ice(:,:) = zqlw_ice(:,:) - zqsb_ice(:,:) - zqla_ice(:,:)
(  647)       ! total non solar heat flux over water 
(  648)       qnsr_oce(:,:) = zqlw_oce(:,:) - zqsb_oce(:,:) - zqla_oce(:,:)
(  649) 
(  650)       ! solid precipitations ( kg/m2/day -> kg/m2/s)
(  651)       tprecip(:,:) = tprecip  (:,:) / rday 
(  652)       ! snow  precipitations ( kg/m2/day -> kg/m2/s)
(  653)       sprecip(:,:) = sprecip  (:,:) / rday  
(  654) !i
(  655)       CALL lbc_lnk( qsr_oce (:,:) , 'T', 1. )
(  656)       CALL lbc_lnk( qnsr_oce(:,:) , 'T', 1. )
(  657)       CALL lbc_lnk( qsr_ice (:,:) , 'T', 1. )
(  658)       CALL lbc_lnk( qnsr_ice(:,:) , 'T', 1. )
(  659)       CALL lbc_lnk( qla_ice (:,:) , 'T', 1. )
(  660)       CALL lbc_lnk( dqns_ice(:,:) , 'T', 1. )
(  661)       CALL lbc_lnk( dqla_ice(:,:) , 'T', 1. )
(  662)       CALL lbc_lnk( fr1_i0  (:,:) , 'T', 1. )
(  663)       CALL lbc_lnk( fr2_i0  (:,:) , 'T', 1. )
(  664)       CALL lbc_lnk( tprecip (:,:) , 'T', 1. )
(  665)       CALL lbc_lnk( sprecip (:,:) , 'T', 1. )
(  666)       CALL lbc_lnk( evap    (:,:) , 'T', 1. )
(  667) !i
(  668) !i
(  669)       qsr_oce (:,:) = qsr_oce (:,:)*tmask(:,:,1)
(  670)       qnsr_oce(:,:) = qnsr_oce(:,:)*tmask(:,:,1)
(  671)       qsr_ice (:,:) = qsr_ice (:,:)*tmask(:,:,1)
(  672)       qnsr_ice(:,:) = qnsr_ice(:,:)*tmask(:,:,1)
(  673)       qla_ice (:,:) = qla_ice (:,:)*tmask(:,:,1)
(  674)       dqns_ice(:,:) = dqns_ice(:,:)*tmask(:,:,1)
(  675)       dqla_ice(:,:) = dqla_ice(:,:)*tmask(:,:,1)
(  676)       fr1_i0  (:,:) = fr1_i0  (:,:)*tmask(:,:,1)
(  677)       fr2_i0  (:,:) = fr2_i0  (:,:)*tmask(:,:,1)
(  678)       tprecip (:,:) = tprecip (:,:)*tmask(:,:,1)
(  679)       sprecip (:,:) = sprecip (:,:)*tmask(:,:,1)
(  680)       evap    (:,:) = evap    (:,:)*tmask(:,:,1)
(  681) !i
(  682) 
(  683)    END SUBROUTINE flx_blk
(  684) 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 13

(  685) 
(  686)    SUBROUTINE flx_blk_declin( ky, kday, pdecl )
(  687)       !!---------------------------------------------------------------------------
(  688)       !!               ***  ROUTINE flx_blk_declin  ***
(  689)       !!          
(  690)       !! ** Purpose :   Computation of the solar declination for the day
(  691)       !!	        kday ( in decimal degrees ).
(  692)       !!       
(  693)       !! ** Method  :
(  694)       !!
(  695)       !! History :
(  696)       !!         original    : 01-04 (LIM)
(  697)       !!         addition    : 02-08 (C. Ethe, G. Madec)
(  698)       !!---------------------------------------------------------------------
(  699)       !! * Arguments
(  700)       INTEGER, INTENT( in ) ::   &
(  701)          ky  ,        &  ! = -1, 0, 1 for odd, normal and leap years resp.
(  702)          kday            ! day of the year ( kday = 1 on january 1)
(  703)       REAL(wp), INTENT(out) ::  &
(  704)          pdecl            ! solar declination
(  705) 
(  706)       !! * Local variables
(  707)       REAL(wp) ::                & 
(  708)          zday                ,   &  ! corresponding day of type year (cf. ky)
(  709)          zp1, zp2, zp3, zp4         ! temporary scalars
(  710)       !!---------------------------------------------------------------------
(  711)       
(  712)       zday = FLOAT( kday ) 
(  713)       
(  714)       IF( ky == 1 )  THEN 
(  715)          zday = zday - 0.5
(  716)       ELSEIF( ky == 3 )  THEN
(  717)          zday = zday - 1.
(  718)       ELSE 
(  719)          zday = REAL( kday )
(  720)       ENDIF
(  721)       
(  722)       zp1 = rpi * ( 2.0 * zday - 367.0 ) / yearday
(  723)       zp2 = 2. * zp1
(  724)       zp3 = 3. * zp1
(  725)       zp4 = 4. * zp1
(  726)       
(  727)       pdecl  = a0                                                                      &
(  728)          &   + a1 * COS( zp1 ) + a2 * COS( zp2 ) + a3 * COS( zp3 ) + a4 * COS( zp4 )   &
(  729)          &   + b1 * SIN( zp1 ) + b2 * SIN( zp2 ) + b3 * SIN( zp3 ) + b4 * SIN( zp4 )
(  730)       
(  731)    END SUBROUTINE flx_blk_declin
(  732) 
(  733) #else
(  734)    !!----------------------------------------------------------------------
(  735)    !!   Default option :           Empty module                     NO bulk
(  736)    !!----------------------------------------------------------------------
(  737) CONTAINS
(  738)    SUBROUTINE flx_blk            ! Empty routine
(  739)    END SUBROUTINE flx_blk
(  740) #endif
(  741) 
(  742)    !!======================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 14

(  743) END MODULE flxblk































































PGF90 (Version     10.5)          09/23/2010  09:49:11      page 15

(    1) # 1 "flxblk.F90"
(    1) !!DB: original flxblk as of 2009.05.04-12:49
(    2) 
(    3) MODULE flxblk
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  flxblk  ***
(    6)    !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
(    7)    !!=====================================================================
(    8) # 9
(    9)    !!----------------------------------------------------------------------
(   10)    !!   'key_flx_bulk_monthly'   or                            MONTHLY bulk
(   11)    !!   'key_flx_bulk_daily'                                     DAILY bulk
(   12)    !!----------------------------------------------------------------------
(   13)    !!   flx_blk        : thermohaline fluxes from bulk
(   14)    !!   flx_blk_declin : solar declinaison
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE cpl_oce         ! ???
(   20)    USE phycst          ! physical constants
(   21)    USE daymod
(   22)    USE blk_oce         ! bulk variables
(   23)    USE flx_oce         ! forcings variables
(   24)    USE ocfzpt          ! ???
(   25)    USE in_out_manager
(   26)    USE lbclnk
(   27)    USE albedo
(   28)    USE prtctl          ! Print control
(   29) 
(   30)    IMPLICIT NONE
(   31)    PRIVATE
(   32) 
(   33)    !! * Accessibility
(   34)    PUBLIC flx_blk        ! routine called by flx.F90 
(   35) 
(   36)    !! * Module variables
(   37)    INTEGER, PARAMETER  ::   &
(   38)       jpintsr = 24          ! number of time step between sunrise and sunset
(   39)       !                     ! uses for heat flux computation
(   40)    LOGICAL ::   &
(   41)       lbulk_init = .TRUE.   ! flag, bulk initialization done or not)
(   42) 
(   43)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   44)       stauc            ,  &   ! cloud optical depth 
(   45)       sbudyko   
(   46) 
(   47)    !! * constants for bulk computation (flx_blk)
(   48)    REAL(wp), DIMENSION(19)  ::  &
(   49)       budyko                  ! BUDYKO's coefficient
(   50)    ! BUDYKO's coefficient (cloudiness effect on LW radiation):
(   51)    DATA budyko / 1.00, 0.98, 0.95, 0.92, 0.89, 0.86, 0.83, 0.80, 0.78, 0.75,  &
(   52)       &          0.72, 0.69, 0.67, 0.64, 0.61, 0.58, 0.56, 0.53, 0.50 /
(   53)    REAL(wp), DIMENSION(20)  :: &
(   54)       tauco                  ! cloud optical depth coefficient
(   55)    ! Cloud optical depth coefficient
(   56)    DATA tauco / 6.6, 6.6, 7.0, 7.2, 7.1, 6.8, 6.5, 6.6, 7.1, 7.6,   &
(   57)       &         6.6, 6.1, 5.6, 5.5, 5.8, 5.8, 5.6, 5.6, 5.6, 5.6 /






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 16

(   58)    REAL(wp)  ::            &  ! constant values
(   59)       zeps    = 1.e-20  ,  &
(   60)       zeps0   = 1.e-13  ,  &
(   61)       zeps1   = 1.e-06  ,  &
(   62)       zzero   = 0.e0    ,  &
(   63)       zone    = 1.0
(   64) 
(   65)    !! * constants for solar declinaison computation (flx_blk_declin)
(   66)    REAL(wp) ::                &
(   67)       a0  =  0.39507671   ,   &  ! coefficients
(   68)       a1  = 22.85684301   ,   &
(   69)       a2  = -0.38637317   ,   &
(   70)       a3  =  0.15096535   ,   &
(   71)       a4  = -0.00961411   ,   &
(   72)       b1  = -4.29692073   ,   &
(   73)       b2  =  0.05702074   ,   &
(   74)       b3  = -0.09028607   ,   &
(   75)       b4  =  0.00592797
(   76)    !!----------------------------------------------------------------------
(   77)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   78)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxblk.F90,v 1.8 2005/09/02 15:45:31 opalod Exp $ 
(   79)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   80)    !!----------------------------------------------------------------------
(   81) 
(   82) CONTAINS
(   83) 
(   84)    SUBROUTINE flx_blk( psst )
(   85)       !!---------------------------------------------------------------------------
(   86)       !!                     ***  ROUTINE flx_blk  ***
(   87)       !!                 
(   88)       !!  ** Purpose :   Computation of the heat fluxes at ocean and snow/ice
(   89)       !!       surface the solar heat at ocean and snow/ice surfaces and the 
(   90)       !!       sensitivity of total heat fluxes to the SST variations
(   91)       !!         
(   92)       !!  ** Method  :   The flux of heat at the ice and ocean surfaces are derived
(   93)       !!       from semi-empirical ( or bulk ) formulae which relate the flux to 
(   94)       !!       the properties of the surface and of the lower atmosphere. Here, we
(   95)       !!       follow the work of Oberhuber, 1988   
(   96)       !!
(   97)       !!  ** Action  :   call flx_blk_albedo to compute ocean and ice albedo 
(   98)       !!          computation of snow precipitation
(   99)       !!          computation of solar flux at the ocean and ice surfaces
(  100)       !!          computation of the long-wave radiation for the ocean and sea/ice
(  101)       !!          computation of turbulent heat fluxes over water and ice
(  102)       !!          computation of evaporation over water
(  103)       !!          computation of total heat fluxes sensitivity over ice (dQ/dT)
(  104)       !!          computation of latent heat flux sensitivity over ice (dQla/dT)
(  105)       !!
(  106)       !! History :
(  107)       !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
(  108)       !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
(  109)       !!----------------------------------------------------------------------
(  110)       !! * Arguments
(  111)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  112)          &                          psst      ! Sea Surface Temperature 
(  113) 
(  114)       !! * Local variables
(  115)       INTEGER  ::             &






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 17

(  116)          ji, jj, jt        ,  &  ! dummy loop indices
(  117)          indaet            ,  &  !  = -1, 0, 1 for odd, normal and leap years resp.
(  118)          iday              ,  &  ! integer part of day
(  119)          indxb             ,  &  ! index for budyko coefficient
(  120)          indxc                   ! index for cloud depth coefficient
(  121) 
(  122)       REAL(wp)  ::            & 
(  123)          zalat , zclat     ,  &  ! latitude in degrees 
(  124)          zmt1, zmt2, zmt3  ,  &  ! tempory air temperatures variables
(  125)          ztatm3, ztatm4    ,  &  ! power 3 and 4 of air temperature
(  126)          z4tatm3           ,  &  ! 4 * ztatm3
(  127)          zcmue             ,  &  ! cosine of local solar altitude
(  128)          zcmue2            ,  &  ! root of zcmue1 
(  129)          zscmue            ,  &  ! square-root of zcmue1 
(  130)          zpcmue            ,  &  ! zcmue1**1.4
(  131)          zdecl             ,  &  ! solar declination
(  132)          zsdecl , zcdecl   ,  &  ! sine and cosine of solar declination 
(  133)          zalbo             ,  &  ! albedo of sea-water
(  134)          zalbi             ,  &  ! albedo of ice
(  135)          ztamr             ,  &  ! air temperature minus triple point of water (rtt)
(  136)          ztaevbk           ,  &  ! part of net longwave radiation
(  137)          zevi , zevo       ,  &  ! vapour pressure of ice and ocean 
(  138)          zind1,zind2,zind3 ,  &  ! switch for testing the values of air temperature
(  139)          zinda             ,  &  ! switch for testing the values of sea ice cover
(  140)          zpis2             ,  &  ! pi / 2
(  141)          z2pi                    ! 2 * pi 
(  142) 
(  143)       REAL(wp)  ::            & 
(  144)          zxday             ,  &  ! day of year
(  145)          zdist             ,  &  ! distance between the sun and the earth during the year
(  146)          zdaycor           ,  &  ! corr. factor to take into account the variation of 
(  147)          !                       ! zday when calc. the solar rad.    
(  148)          zesi, zeso        ,  &  ! vapour pressure of ice and ocean at saturation
(  149)          zesi2             ,  &  ! root of zesi 
(  150)          zqsato            ,  &  ! humidity close to the ocean surface (at saturation)   
(  151)          zqsati            ,  &  ! humidity close to the ice surface (at saturation) 
(  152)          zqsati2           ,  &  ! root of  zqsati 
(  153)          zdesidt           ,  &  ! derivative of zesi, function of ice temperature
(  154)          zdteta            ,  &  ! diff. betw. sst and air temperature
(  155)          zdeltaq           ,  &  ! diff. betw. spec. hum. and hum. close to the surface
(  156)          ztvmoy, zobouks   ,  &  ! tempory scalars
(  157)          zpsims, zpsihs, zpsils, zobouku, zxins, zpsimu ,  & 
(  158)          zpsihu, zpsilu, zstab,zpsim, zpsih, zpsil      ,  & 
(  159)          zvatmg, zcmn, zchn, zcln, zcmcmn, zdenum       ,  & 
(  160)          zdtetar, ztvmoyr, zlxins, zcmn2, zchcm, zclcm , zcoef
(  161) 
(  162)       REAL(wp)  ::            & 
(  163)          zrhova            ,  &  ! air density per wind speed
(  164)          zcsho , zcleo     ,  &  ! transfer coefficient over ocean
(  165)          zcshi , zclei     ,  &  ! transfer coefficient over ice-free
(  166)          zrhovacleo        ,  &  ! air density per wind speed per transfer coef.
(  167)          zrhovacsho, zrhovaclei, zrhovacshi, & 
(  168)          ztice3            ,  &  ! power 3 of ice temperature
(  169)          zticemb, zticemb2 ,  &  ! tempory air temperatures variables
(  170)          zdqlw_ice         ,  &  ! sensitivity of long-wave flux over ice
(  171)          zdqsb_ice         ,  &  ! sensitivity of sensible heat flux over ice
(  172)          zdqla_ice         ,  &  ! sensitivity of latent heat flux over ice
(  173)          zdl, zdr                ! fractionnal part of latitude






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 18

(  174)       REAL(wp), DIMENSION(jpi,jpj) :: & 
(  175)          zpatm            ,  &   ! atmospheric pressure
(  176)          zqatm            ,  &   ! specific humidity
(  177)          zes              ,  &   ! vapour pressure at saturation
(  178)          zev, zevsqr      ,  &   ! vapour pressure and his square-root
(  179)          zrhoa            ,  &   ! air density
(  180)          ztatm            ,  &   ! air temperature in Kelvins
(  181)          zfrld            ,  &   ! fraction of sea ice cover 
(  182)          zcatm1           ,  &   ! fraction of cloud
(  183)          zcldeff                 ! correction factor to account cloud effect
(  184)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  185)          zalbocsd         ,  &   ! albedo of ocean
(  186)          zalboos          ,  &   ! albedo of ocean under overcast sky
(  187)          zalbics          ,  &   ! albedo of ice under clear sky
(  188)          zalbios          ,  &   ! albedo of ice under overcast sky
(  189)          zalbomu          ,  &   ! albedo of ocean when zcmue is 0.4
(  190)          zqsro            ,  &   ! solar radiation over ocean
(  191)          zqsrics          ,  &   ! solar radiation over ice under clear sky
(  192)          zqsrios          ,  &   ! solar radiation over ice under overcast sky
(  193)          zcldcor          ,  &   ! cloud correction
(  194)          zlsrise, zlsset  ,  &   ! sunrise and sunset
(  195)          zlmunoon         ,  &   ! local noon solar altitude
(  196)          zdlha            ,  &   ! length of the ninstr segments of the solar day
(  197)          zps              ,  &   ! sine of latitude per sine of solar decli.
(  198)          zpc                     ! cosine of latitude per cosine of solar decli. 
(  199) 
(  200)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  201)          zqlw_oce         ,  &   ! long-wave heat flux over ocean
(  202)          zqlw_ice         ,  &   ! long-wave heat flux over ice
(  203)          zqla_oce         ,  &   ! latent heat flux over ocean
(  204)          zqla_ice         ,  &   ! latent heat flux over ice
(  205)          zqsb_oce         ,  &   ! sensible heat flux over ocean
(  206)          zqsb_ice                ! sensible heat flux over ice
(  207)  
(  208)       REAL(wp), DIMENSION(jpi,jpj,jpintsr) ::    &
(  209)          zlha             ,  &   ! local hour angle
(  210)          zalbocs          ,  &   ! tempory var. of ocean albedo under clear sky
(  211)          zsqsro           ,  &   ! tempory var. of solar rad. over ocean 
(  212)          zsqsrics         ,  &   ! temp. var. of solar rad. over ice under clear sky
(  213)          zsqsrios                ! temp. var. of solar rad. over ice under overcast sky
(  214)       !!---------------------------------------------------------------------
(  215) 
(  216)       !---------------------
(  217)       !   Initilization    !
(  218)       !---------------------
(  219) 
(  220) # 223
(  223)       !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
(  224)       !  and the correction factor for taking into account  the effect of clouds 
(  225)       !------------------------------------------------------
(  226)       IF( lbulk_init ) THEN
(  227)          DO jj = 1, jpj  
(  228)             DO ji = 1 , jpi
(  229)                zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
(  230)                zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
(  231)                indxb          = 1 + INT( zalat ) 
(  232)                !  correction factor to account for the effect of clouds 
(  233)                sbudyko(ji,jj) = budyko(indxb)  






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 19

(  234)                indxc          = 1 + INT( zclat )  
(  235)                zdl            = zclat - INT( zclat ) 
(  236)                zdr            = 1.0 - zdl
(  237)                stauc(ji,jj)   = zdr * tauco( indxc ) + zdl * tauco( indxc + 1 ) 
(  238)             END DO
(  239)          END DO
(  240)          IF( nleapy == 1 ) THEN
(  241)             yearday = 366.e0
(  242)          ELSE IF( nleapy == 0 ) THEN
(  243)             yearday = 365.e0
(  244)          ELSEIF( nleapy == 30) THEN
(  245)             yearday = 360.e0
(  246)          ENDIF
(  247)          lbulk_init = .FALSE.
(  248)       ENDIF
(  249) 
(  250)       zqlw_oce(:,:) = 0.e0
(  251)       zqla_oce(:,:) = 0.e0
(  252)       zqsb_oce(:,:) = 0.e0
(  253)       zqlw_ice(:,:) = 0.e0
(  254)       zqla_ice(:,:) = 0.e0
(  255)       zqsb_ice(:,:) = 0.e0
(  256) 
(  257)       zpis2       = rpi / 2.
(  258)       z2pi        = 2. * rpi
(  259) 
(  260)  !CDIR NOVERRCHK
(  261)       DO jj = 1, jpj
(  262)  !CDIR NOVERRCHK
(  263)          DO ji = 1, jpi
(  264) 
(  265)             ztatm (ji,jj) = 273.15 + tatm  (ji,jj)  !  air temperature in Kelvins 
(  266)             zcatm1(ji,jj) = 1.0    - catm  (ji,jj)  !  fractional cloud cover
(  267)             zfrld (ji,jj) = 1.0    - freeze(ji,jj)  !  fractional sea ice cover
(  268)             zpatm(ji,jj)  = 101000.               !  pressure 
(  269)       
(  270)             !  Computation of air density, obtained from the equation of state for dry air. 
(  271)             zrhoa(ji,jj) = zpatm(ji,jj) / ( 287.04 * ztatm(ji,jj) )
(  272)       
(  273)             !  zes : Saturation water vapour
(  274)             ztamr = ztatm(ji,jj) - rtt
(  275)             zmt1  = SIGN( 17.269, ztamr )
(  276)             zmt2  = SIGN( 21.875, ztamr )
(  277)             zmt3  = SIGN( 28.200, -ztamr )
(  278)             zes(ji,jj) = 611.0 * EXP (  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &
(  279)                &                      / ( ztatm(ji,jj) - 35.86  + MAX( zzero, zmt3 ) ) )
(  280) 
(  281)             !  zev : vapour pressure  (hatm is relative humidity)  
(  282)             zev(ji,jj)   = hatm(ji,jj) * zes(ji,jj) 
(  283)             !  square-root of vapour pressure
(  284) !CDIR NOVERRCHK
(  285)             zevsqr(ji,jj) = SQRT( zev(ji,jj) * 0.01 )
(  286)             !  zqapb  : specific humidity 
(  287)             zqatm(ji,jj) = 0.622 * zev(ji,jj) / ( zpatm(ji,jj) - 0.378 * zev(ji,jj) )
(  288) 
(  289) 
(  290)             !----------------------------------------------------
(  291)             !   Computation of snow precipitation (Ledley, 1985) |






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 20

(  292)             !----------------------------------------------------
(  293) 
(  294)             zmt1  =   253.0 - ztatm(ji,jj)
(  295)             zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0 
(  296)             zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0
(  297)             zind1 = MAX( zzero, SIGN( zone, zmt1 ) )
(  298)             zind2 = MAX( zzero, SIGN( zone, zmt2 ) )
(  299)             zind3 = MAX( zzero, SIGN( zone, zmt3 ) )
(  300)             ! total precipitation
(  301)             tprecip(ji,jj) = watm(ji,jj)
(  302)             ! solid  (snow) precipitation
(  303)             sprecip(ji,jj) = tprecip(ji,jj) *       &
(  304)                &             (           zind1      &
(  305)                &               + ( 1.0 - zind1 ) * ( zind2 * ( 0.5 + zmt2 ) + ( 1.0 - zind2 ) *  zind3 * zmt3 ) ) 
(  306)          END DO
(  307)       END DO
(  308) 
(  309)       !----------------------------------------------------------
(  310)       !   Computation of albedo (need to calculates heat fluxes)|
(  311)       !-----------------------------------------------------------
(  312)       
(  313)       CALL flx_blk_albedo( zalbios, zalboos, zalbics, zalbomu )
(  314) 
(  315)       !-------------------------------------
(  316)       !   Computation of solar irradiance. |
(  317)       !----------------------------------------
(  318)       indaet   = 1  
(  319)       !  compution of the day of the year at which the fluxes have to be calculate 
(  320)       !--The date corresponds to the middle of the time step.
(  321)       zxday=nday_year + rdtbs2/rday
(  322) 
(  323)       iday   = INT( zxday )
(  324) 
(  325)       IF(ln_ctl) CALL prt_ctl_info('declin : iday ', ivar1=iday, clinfo2=' nfbulk= ', ivar2=nfbulk)
(  326) 
(  327)       !   computation of the solar declination, his sine and his cosine
(  328)       CALL flx_blk_declin( indaet, iday, zdecl )
(  329)       
(  330)       zdecl    = zdecl * rad
(  331)       zsdecl   = SIN( zdecl )
(  332)       zcdecl   = COS( zdecl )
(  333)       
(  334)       !  correction factor added for computation of shortwave flux to take into account the variation of
(  335)       !  the distance between the sun and the earth during the year (Oberhuber 1988)
(  336)       zdist    = zxday * z2pi / yearday
(  337)       zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )
(  338) 
(  339) !CDIR NOVERRCHK
(  340)       DO jj = 1, jpj
(  341) !CDIR NOVERRCHK
(  342)          DO ji = 1, jpi
(  343)             !  product of sine of latitude and sine of solar declination
(  344)             zps     (ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
(  345)             !  product of cosine of latitude and cosine of solar declination
(  346)             zpc     (ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
(  347)             !  computation of the both local time of sunrise and sunset
(  348)             zlsrise (ji,jj) = ACOS( - SIGN( zone, zps(ji,jj) ) * MIN( zone, SIGN( zone, zps(ji,jj) )  &
(  349)                &                     * ( zps(ji,jj) / zpc(ji,jj) ) ) ) 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 21

(  350)             zlsset  (ji,jj) = - zlsrise(ji,jj)
(  351)             !  dividing the solar day into jpintsr segments of length zdlha
(  352)             zdlha   (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jpintsr )
(  353)             !  computation of the local noon solar altitude
(  354)             zlmunoon(ji,jj) = ASIN ( ( zps(ji,jj) + zpc(ji,jj) ) ) / rad
(  355)             
(  356)             !  cloud correction taken from Reed (1977) (imposed lower than 1)
(  357)             zcldcor (ji,jj) = MIN( zone, ( 1.e0 - 0.62 * catm(ji,jj) + 0.0019 * zlmunoon(ji,jj) ) )
(  358)          END DO
(  359)       END DO
(  360) 
(  361)          !  Computation of solar heat flux at each time of the day between sunrise and sunset. 
(  362)          !  We do this to a better optimisation of the code 
(  363)          !------------------------------------------------------       
(  364) 
(  365) !CDIR NOVERRCHK   
(  366)       DO jt = 1, jpintsr   
(  367)          zcoef = FLOAT( jt ) - 0.5
(  368) !CDIR NOVERRCHK     
(  369)          DO jj = 1, jpj
(  370) !CDIR NOVERRCHK
(  371)             DO ji = 1, jpi
(  372)                !  local hour angle
(  373)                zlha (ji,jj,jt) = COS ( zlsrise(ji,jj) - zcoef * zdlha(ji,jj) )
(  374) 
(  375)                ! cosine of local solar altitude
(  376)                zcmue              = MAX ( zzero ,   zps(ji,jj) + zpc(ji,jj) * zlha (ji,jj,jt)  )
(  377)                zcmue2             = 1368.0 * zcmue * zcmue
(  378)                zscmue             = SQRT ( zcmue )
(  379)                zpcmue             = zcmue**1.4
(  380)                ! computation of sea-water albedo (Payne, 1972)
(  381)                zalbocs(ji,jj,jt)  = 0.05 / ( 1.1 * zpcmue + 0.15 )
(  382)                zalbo              = zcatm1(ji,jj) * zalbocs(ji,jj,jt) + catm(ji,jj) * zalboos(ji,jj)
(  383)                ! solar heat flux absorbed at ocean surfaces (Zillman, 1972)
(  384)                zevo               = zev(ji,jj) * 1.0e-05
(  385)                zsqsro(ji,jj,jt)   =  ( 1.0 - zalbo ) * zdlha(ji,jj) * zcmue2                &
(  386)                                    / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
(  387)                !  solar heat flux absorbed at sea/ice surfaces 
(  388)                !  Formulation of Shine and Crane, 1984 adapted for high albedo surfaces 
(  389) 
(  390)                !  For clear sky        
(  391)                zevi               = zevo
(  392)                zalbi              = zalbics(ji,jj)
(  393)                zsqsrics(ji,jj,jt) =  ( 1.0 - zalbi ) * zdlha(ji,jj) * zcmue2                &
(  394)                   &                / ( ( 1.0 + zcmue ) * zevi + 1.2 * zcmue + 0.0455 )
(  395) 
(  396)                ! For overcast sky
(  397)                zalbi              = zalbios(ji,jj)
(  398)                zsqsrios(ji,jj,jt) = zdlha(ji,jj) *                                                           &
(  399)                   &                 ( ( 53.5 + 1274.5 * zcmue      ) *  zscmue * ( 1.0 - 0.996  * zalbi ) )  &
(  400)                   &                 / (  1.0 + 0.139  * stauc(ji,jj) *           ( 1.0 - 0.9435 * zalbi ) )
(  401)             END DO
(  402)          END DO
(  403)       END DO
(  404) 
(  405) 
(  406)       !  Computation of daily (between sunrise and sunset) solar heat flux absorbed 
(  407)       !  at the ocean and snow/ice surfaces.






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 22

(  408)       !--------------------------------------------------------------------
(  409) 
(  410)       zalbocsd(:,:) = 0.e0
(  411)       zqsro   (:,:) = 0.e0
(  412)       zqsrics (:,:) = 0.e0
(  413)       zqsrios (:,:) = 0.e0
(  414) 
(  415)       DO jt = 1, jpintsr 
(  416) # 425
(  425)          DO jj = 1, jpj
(  426)             DO ji = 1, jpi  
(  427)                zalbocsd(ji,jj) = zalbocsd(ji,jj) + zdlha(ji,jj) * zalbocs(ji,jj,jt)   &
(  428)                   &                                              / MAX( 2.0 * zlsrise(ji,jj) , zeps0 )
(  429)                zqsro  (ji,jj)  = zqsro   (ji,jj) + zsqsro  (ji,jj,jt)
(  430)                zqsrics(ji,jj)  = zqsrics (ji,jj) + zsqsrics(ji,jj,jt)
(  431)                zqsrios(ji,jj)  = zqsrios (ji,jj) + zsqsrios(ji,jj,jt)
(  432)             END DO
(  433)          END DO
(  434) # 435
(  435)       END DO
(  436) 
(  437)       DO jj = 1, jpj
(  438)          DO ji = 1, jpi 
(  439) 
(  440)             !------------------------------------------- 
(  441)             !  Computation of shortwave radiation.
(  442)             !-------------------------------------------
(  443) 
(  444)             ! the solar heat flux absorbed at ocean and snow/ice surfaces
(  445)             !------------------------------------------------------------
(  446) 
(  447)             ! For ocean
(  448)             qsr_oce(ji,jj) = srgamma * zcldcor(ji,jj) * zqsro(ji,jj) / z2pi
(  449)             zinda          = SIGN( zone , -( -0.5 - zfrld(ji,jj) ) )
(  450)             zinda          = 1.0 - MAX( zzero , zinda )
(  451)             qsr_oce(ji,jj) = ( 1.- zinda ) * qsr_oce(ji,jj)
(  452) 
(  453)             ! For snow/ice 
(  454)             qsr_ice(ji,jj) = ( zcatm1(ji,jj) * zqsrics(ji,jj) + catm(ji,jj) * zqsrios(ji,jj) ) / z2pi
(  455) 
(  456) 
(  457)             ! Taking into account the ellipsity of the earth orbit
(  458)             !-----------------------------------------------------
(  459) 
(  460)             qsr_ice(ji,jj) = qsr_ice(ji,jj) * zdaycor
(  461)             qsr_oce(ji,jj) = qsr_oce(ji,jj) * zdaycor
(  462) 
(  463)             !  fraction of net shortwave radiation which is not absorbed in the 
(  464)             !  thin surface layer and penetrates inside the ice cover 
(  465)             !  ( Maykut and Untersteiner, 1971 ; Elbert anbd Curry, 1993 )
(  466)             !------------------------------------------------------------------
(  467) 
(  468)             fr1_i0(ji,jj) = 0.18  * zcatm1(ji,jj) + 0.35 * catm(ji,jj) 
(  469)             fr2_i0(ji,jj) = 0.82  * zcatm1(ji,jj) + 0.65 * catm(ji,jj)
(  470) 
(  471)             !---------------------------------------------------------------------------
(  472)             !   Computation of long-wave radiation  ( Berliand 1952 ; all latitudes )
(  473)             !---------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 23

(  474) 
(  475)             ! tempory variables
(  476)             ztatm3         = ztatm(ji,jj) * ztatm(ji,jj) * ztatm(ji,jj)
(  477)             ztatm4         = ztatm3 * ztatm(ji,jj)
(  478)             z4tatm3        = 4. * ztatm3
(  479)             zcldeff(ji,jj) = 1.0 - sbudyko(ji,jj) * catm(ji,jj) * catm(ji,jj)    
(  480)             ztaevbk        = ztatm4 * zcldeff(ji,jj) * ( 0.39 - 0.05 * zevsqr(ji,jj) ) 
(  481) 
(  482)             !  Long-Wave for Ice
(  483)             !----------------------
(  484)             zqlw_ice(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( tn_ice(ji,jj) - ztatm(ji,jj) ) ) 
(  485) 
(  486)             !  Long-Wave for Ocean
(  487)             !-----------------------
(  488)             zqlw_oce(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( psst  (ji,jj) - ztatm(ji,jj) ) ) 
(  489) 
(  490)          END DO
(  491)       END DO
(  492) 
(  493)       !----------------------------------------
(  494)       !  Computation of turbulent heat fluxes  ( Latent and sensible ) 
(  495)       !----------------------------------------        
(  496)       !CDIR NOVERRCHK
(  497)       DO jj = 2 , jpjm1
(  498) !ib   DO jj = 1 , jpj
(  499)          !CDIR NOVERRCHK
(  500)          DO  ji = 1, jpi
(  501) 
(  502)             !  Turbulent heat fluxes over water
(  503)             !----------------------------------
(  504) 
(  505)             ! zeso     : vapour pressure at saturation of ocean
(  506)             ! zqsato   : humidity close to the ocean surface (at saturation)
(  507)             zeso          =  611.0 * EXP ( 17.2693884 * ( psst(ji,jj) - rtt ) * tmask(ji,jj,1) / ( psst(ji,jj) - 35.86 ) )
(  508)             zqsato        = ( 0.622 * zeso ) / ( zpatm(ji,jj) - 0.378 * zeso )
(  509) 
(  510)             !  Drag coefficients from Large and Pond (1981,1982)
(  511)             !---------------------------------------------------
(  512)     
(  513)             !  Stability parameters
(  514)             zdteta         = psst(ji,jj) - ztatm(ji,jj)
(  515)             zdeltaq        = zqatm(ji,jj) - zqsato
(  516)             ztvmoy         = ztatm(ji,jj) * ( 1. + 2.2e-3 * ztatm(ji,jj) * zqatm(ji,jj) )
(  517)             zdenum         = MAX( vatm(ji,jj) * vatm(ji,jj) * ztvmoy, zeps )
(  518) !i
(  519) !i          if( zdenum == 0.e0 ) then
(  520) !i               write(numout,*) 'flxblk  zdenum=0 ', ji,jj
(  521) !i               zdenum = zeps
(  522) !i          endif
(  523) !i
(  524)             zdtetar        = zdteta / zdenum
(  525)             ztvmoyr        = ztvmoy * ztvmoy * zdeltaq / zdenum
(  526)             
(  527)             ! For stable atmospheric conditions
(  528)             zobouks        = -70.0 * 10. * ( zdtetar + 3.2e-3 * ztvmoyr )
(  529)             zobouks        = MAX( zzero , zobouks )
(  530)             zpsims         = -7.0 * zobouks
(  531)             zpsihs         =  zpsims






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 24

(  532)             zpsils         =  zpsims
(  533) 
(  534)             !  For unstable atmospheric conditions
(  535)             zobouku        = -100.0 * 10.0 * ( zdtetar + 2.2e-3 * ztvmoyr )
(  536)             zobouku        = MIN( zzero , zobouku )
(  537)             zxins          = ( 1. - 16. * zobouku )**0.25
(  538)             zlxins         = LOG( ( 1. + zxins * zxins ) / 2. )
(  539)             zpsimu         = 2. * LOG( ( 1 + zxins ) / 2. )  + zlxins - 2. * ATAN( zxins ) + zpis2
(  540)             zpsihu         = 2. * zlxins
(  541)             zpsilu         = zpsihu
(  542) 
(  543)             ! computation of intermediate values
(  544)             zstab          = MAX( zzero , SIGN( zone , zdteta ) )
(  545)             zpsim          = zstab * zpsimu + (1.0 - zstab ) * zpsims
(  546)             zpsih          = zstab * zpsihu + (1.0 - zstab ) * zpsihs
(  547)             zpsil          = zpsih
(  548)             
(  549)             zvatmg         = MAX( 0.032 * 1.5e-3 * vatm(ji,jj) * vatm(ji,jj) / grav, zeps )
(  550) !i
(  551) !!          if( zvatmg == 0.e0 ) then
(  552) !!               write(numout,*) 'flxblk  zvatmg=0 ', ji,jj
(  553) !!               zvatmg = zeps
(  554) !!          endif
(  555) !i
(  556) 
(  557)             zcmn           = vkarmn / LOG ( 10. / zvatmg )
(  558)             zcmn2          = zcmn * zcmn
(  559)             zchn           = 0.0327 * zcmn
(  560)             zcln           = 0.0346 * zcmn
(  561)             zcmcmn         = 1 / ( 1 - zcmn * zpsim / vkarmn )
(  562)             zchcm          = zcmcmn / ( 1 - zchn * zpsih / ( vkarmn * zcmn ) )
(  563)             zclcm          = zchcm
(  564) 
(  565) 
(  566)             !  Transfer cofficient zcsho and zcleo over ocean according to Large and Pond (1981,1982)
(  567)             !-------------------------------------------------------------- 
(  568)             zcsho          = zchn * zchcm
(  569)             zcleo          = zcln * zclcm 
(  570) 
(  571) 
(  572)             !   Computation of sensible and latent fluxes over Ocean 
(  573)             !----------------------------------------------------------------
(  574) 
(  575)             !  computation of intermediate values
(  576)             zrhova         = zrhoa(ji,jj) * vatm(ji,jj)
(  577)             zrhovacsho     = zrhova * zcsho
(  578)             zrhovacleo     = zrhova * zcleo
(  579) 
(  580)             ! sensible heat flux
(  581)             zqsb_oce(ji,jj) = zrhovacsho * 1004.0  * ( psst(ji,jj) - ztatm(ji,jj) )  
(  582)          
(  583)             !  latent heat flux 
(  584)             zqla_oce(ji,jj) = MAX(0.e0, zrhovacleo * 2.5e+06 * ( zqsato      - zqatm(ji,jj) ) )
(  585)                
(  586)             !  Calculate evaporation over water. (kg/m2/s)
(  587)             !-------------------------------------------------
(  588)             evap(ji,jj)    = zqla_oce(ji,jj) / cevap
(  589)                






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 25

(  590)                
(  591)             !  Turbulent heat fluxes over snow/ice.
(  592)             !--------------------------------------------------
(  593)             
(  594)             !  zesi     : vapour pressure at saturation of ice
(  595)             !  zqsati   : humidity close to the ice surface (at saturation)
(  596)             zesi           =  611.0 * EXP ( 21.8745587 * tmask(ji,jj,1)   &   ! tmask needed to avoid overflow in the exponential
(  597)                &                                       * ( tn_ice(ji,jj) - rtt ) / ( tn_ice(ji,jj) - 7.66 ) )
(  598)             zqsati         = ( 0.622 * zesi ) / ( zpatm(ji,jj) - 0.378 * zesi )
(  599)                
(  600)             !  computation of intermediate values
(  601)             zticemb        = ( tn_ice(ji,jj) - 7.66 )
(  602)             zticemb2       = zticemb * zticemb  
(  603)             ztice3         = tn_ice(ji,jj) * tn_ice(ji,jj) * tn_ice(ji,jj)
(  604)             zqsati2        = zqsati * zqsati
(  605)             zesi2          = zesi * zesi
(  606)             zdesidt        = zesi * ( 9.5 * LOG( 10.0 ) * ( rtt - 7.66 )  / zticemb2 )
(  607)                
(  608)             !  Transfer cofficient zcshi and zclei over ice. Assumed to be constant Parkinson 1979 ; Maykut 1982
(  609)             !--------------------------------------------------------------------
(  610)             zcshi          = 1.75e-03
(  611)             zclei          = zcshi
(  612)                
(  613)             !  Computation of sensible and latent fluxes over ice
(  614)             !----------------------------------------------------------------
(  615)                
(  616)             !  computation of intermediate values
(  617)             zrhova          = zrhoa(ji,jj) * vatm(ji,jj)
(  618)             zrhovacshi      = zrhova * zcshi * 2.834e+06
(  619)             zrhovaclei      = zrhova * zclei * 1004.0
(  620)             
(  621)             !  sensible heat flux
(  622)             zqsb_ice(ji,jj) = zrhovaclei * ( tn_ice(ji,jj) - ztatm(ji,jj) )
(  623)             
(  624)             !  latent heat flux 
(  625)             zqla_ice(ji,jj) = zrhovacshi * ( zqsati        - zqatm(ji,jj) )
(  626)             qla_ice (ji,jj) = MAX(0.e0, zqla_ice(ji,jj) )
(  627)               
(  628)             !  Computation of sensitivity of non solar fluxes (dQ/dT)
(  629)             !---------------------------------------------------------------
(  630)                
(  631)             !  computation of long-wave, sensible and latent flux sensitivity
(  632)             zdqlw_ice       = 4.0 * emic * stefan * ztice3
(  633)             zdqsb_ice       = zrhovaclei
(  634)             zdqla_ice       = zrhovacshi * ( zdesidt * ( zqsati2 / zesi2 ) * ( zpatm(ji,jj) / 0.622 ) )         
(  635)             
(  636)             !  total non solar sensitivity
(  637)             dqns_ice(ji,jj) = -( zdqlw_ice + zdqsb_ice + zdqla_ice ) 
(  638)             
(  639)             ! latent flux sensitivity
(  640)             dqla_ice(ji,jj) = zdqla_ice
(  641)             
(  642)          END DO
(  643)       END DO
(  644) 
(  645)       ! total non solar heat flux over ice
(  646)       qnsr_ice(:,:) = zqlw_ice(:,:) - zqsb_ice(:,:) - zqla_ice(:,:)
(  647)       ! total non solar heat flux over water 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 26

(  648)       qnsr_oce(:,:) = zqlw_oce(:,:) - zqsb_oce(:,:) - zqla_oce(:,:)
(  649) 
(  650)       ! solid precipitations ( kg/m2/day -> kg/m2/s)
(  651)       tprecip(:,:) = tprecip  (:,:) / rday 
(  652)       ! snow  precipitations ( kg/m2/day -> kg/m2/s)
(  653)       sprecip(:,:) = sprecip  (:,:) / rday  
(  654) !i
(  655)       CALL lbc_lnk( qsr_oce (:,:) , 'T', 1. )
(  656)       CALL lbc_lnk( qnsr_oce(:,:) , 'T', 1. )
(  657)       CALL lbc_lnk( qsr_ice (:,:) , 'T', 1. )
(  658)       CALL lbc_lnk( qnsr_ice(:,:) , 'T', 1. )
(  659)       CALL lbc_lnk( qla_ice (:,:) , 'T', 1. )
(  660)       CALL lbc_lnk( dqns_ice(:,:) , 'T', 1. )
(  661)       CALL lbc_lnk( dqla_ice(:,:) , 'T', 1. )
(  662)       CALL lbc_lnk( fr1_i0  (:,:) , 'T', 1. )
(  663)       CALL lbc_lnk( fr2_i0  (:,:) , 'T', 1. )
(  664)       CALL lbc_lnk( tprecip (:,:) , 'T', 1. )
(  665)       CALL lbc_lnk( sprecip (:,:) , 'T', 1. )
(  666)       CALL lbc_lnk( evap    (:,:) , 'T', 1. )
(  667) !i
(  668) !i
(  669)       qsr_oce (:,:) = qsr_oce (:,:)*tmask(:,:,1)
(  670)       qnsr_oce(:,:) = qnsr_oce(:,:)*tmask(:,:,1)
(  671)       qsr_ice (:,:) = qsr_ice (:,:)*tmask(:,:,1)
(  672)       qnsr_ice(:,:) = qnsr_ice(:,:)*tmask(:,:,1)
(  673)       qla_ice (:,:) = qla_ice (:,:)*tmask(:,:,1)
(  674)       dqns_ice(:,:) = dqns_ice(:,:)*tmask(:,:,1)
(  675)       dqla_ice(:,:) = dqla_ice(:,:)*tmask(:,:,1)
(  676)       fr1_i0  (:,:) = fr1_i0  (:,:)*tmask(:,:,1)
(  677)       fr2_i0  (:,:) = fr2_i0  (:,:)*tmask(:,:,1)
(  678)       tprecip (:,:) = tprecip (:,:)*tmask(:,:,1)
(  679)       sprecip (:,:) = sprecip (:,:)*tmask(:,:,1)
(  680)       evap    (:,:) = evap    (:,:)*tmask(:,:,1)
(  681) !i
(  682) 
(  683)    END SUBROUTINE flx_blk
(  684) 
(  685) 
(  686)    SUBROUTINE flx_blk_declin( ky, kday, pdecl )
(  687)       !!---------------------------------------------------------------------------
(  688)       !!               ***  ROUTINE flx_blk_declin  ***
(  689)       !!          
(  690)       !! ** Purpose :   Computation of the solar declination for the day
(  691)       !!	        kday ( in decimal degrees ).
(  692)       !!       
(  693)       !! ** Method  :
(  694)       !!
(  695)       !! History :
(  696)       !!         original    : 01-04 (LIM)
(  697)       !!         addition    : 02-08 (C. Ethe, G. Madec)
(  698)       !!---------------------------------------------------------------------
(  699)       !! * Arguments
(  700)       INTEGER, INTENT( in ) ::   &
(  701)          ky  ,        &  ! = -1, 0, 1 for odd, normal and leap years resp.
(  702)          kday            ! day of the year ( kday = 1 on january 1)
(  703)       REAL(wp), INTENT(out) ::  &
(  704)          pdecl            ! solar declination
(  705) 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 27

(  706)       !! * Local variables
(  707)       REAL(wp) ::                & 
(  708)          zday                ,   &  ! corresponding day of type year (cf. ky)
(  709)          zp1, zp2, zp3, zp4         ! temporary scalars
(  710)       !!---------------------------------------------------------------------
(  711)       
(  712)       zday = FLOAT( kday ) 
(  713)       
(  714)       IF( ky == 1 )  THEN 
(  715)          zday = zday - 0.5
(  716)       ELSEIF( ky == 3 )  THEN
(  717)          zday = zday - 1.
(  718)       ELSE 
(  719)          zday = REAL( kday )
(  720)       ENDIF
(  721)       
(  722)       zp1 = rpi * ( 2.0 * zday - 367.0 ) / yearday
(  723)       zp2 = 2. * zp1
(  724)       zp3 = 3. * zp1
(  725)       zp4 = 4. * zp1
(  726)       
(  727)       pdecl  = a0                                                                      &
(  728)          &   + a1 * COS( zp1 ) + a2 * COS( zp2 ) + a3 * COS( zp3 ) + a4 * COS( zp4 )   &
(  729)          &   + b1 * SIN( zp1 ) + b2 * SIN( zp2 ) + b3 * SIN( zp3 ) + b4 * SIN( zp4 )
(  730)       
(  731)    END SUBROUTINE flx_blk_declin
(  732) 
(  733) 
(  734) # 742
(  742)    !!======================================================================
(  743) END MODULE flxblk
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (flxblk.F90: 367)
