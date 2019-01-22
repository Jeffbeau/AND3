!!DB: original flxblk as of 2009.05.04-12:49

MODULE flxblk
   !!======================================================================
   !!                       ***  MODULE  flxblk  ***
   !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
   !!=====================================================================
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   !!----------------------------------------------------------------------
   !!   'key_flx_bulk_monthly'   or                            MONTHLY bulk
   !!   'key_flx_bulk_daily'                                     DAILY bulk
   !!----------------------------------------------------------------------
   !!   flx_blk        : thermohaline fluxes from bulk
   !!   flx_blk_declin : solar declinaison
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE cpl_oce         ! ???
   USE phycst          ! physical constants
   USE daymod
   USE blk_oce         ! bulk variables
   USE flx_oce         ! forcings variables
   USE ocfzpt          ! ???
   USE in_out_manager
   USE lbclnk
   USE albedo
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC flx_blk        ! routine called by flx.F90 

   !! * Module variables
   INTEGER, PARAMETER  ::   &
      jpintsr = 24          ! number of time step between sunrise and sunset
      !                     ! uses for heat flux computation
   LOGICAL ::   &
      lbulk_init = .TRUE.   ! flag, bulk initialization done or not)

   REAL(wp), DIMENSION(jpi,jpj) ::   &
      stauc            ,  &   ! cloud optical depth 
      sbudyko   

   !! * constants for bulk computation (flx_blk)
   REAL(wp), DIMENSION(19)  ::  &
      budyko                  ! BUDYKO's coefficient
   ! BUDYKO's coefficient (cloudiness effect on LW radiation):
   DATA budyko / 1.00, 0.98, 0.95, 0.92, 0.89, 0.86, 0.83, 0.80, 0.78, 0.75,  &
      &          0.72, 0.69, 0.67, 0.64, 0.61, 0.58, 0.56, 0.53, 0.50 /
   REAL(wp), DIMENSION(20)  :: &
      tauco                  ! cloud optical depth coefficient
   ! Cloud optical depth coefficient
   DATA tauco / 6.6, 6.6, 7.0, 7.2, 7.1, 6.8, 6.5, 6.6, 7.1, 7.6,   &
      &         6.6, 6.1, 5.6, 5.5, 5.8, 5.8, 5.6, 5.6, 5.6, 5.6 /
   REAL(wp)  ::            &  ! constant values
      zeps    = 1.e-20  ,  &
      zeps0   = 1.e-13  ,  &
      zeps1   = 1.e-06  ,  &
      zzero   = 0.e0    ,  &
      zone    = 1.0

   !! * constants for solar declinaison computation (flx_blk_declin)
   REAL(wp) ::                &
      a0  =  0.39507671   ,   &  ! coefficients
      a1  = 22.85684301   ,   &
      a2  = -0.38637317   ,   &
      a3  =  0.15096535   ,   &
      a4  = -0.00961411   ,   &
      b1  = -4.29692073   ,   &
      b2  =  0.05702074   ,   &
      b3  = -0.09028607   ,   &
      b4  =  0.00592797
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxblk.F90,v 1.8 2005/09/02 15:45:31 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flx_blk( psst )
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE flx_blk  ***
      !!                 
      !!  ** Purpose :   Computation of the heat fluxes at ocean and snow/ice
      !!       surface the solar heat at ocean and snow/ice surfaces and the 
      !!       sensitivity of total heat fluxes to the SST variations
      !!         
      !!  ** Method  :   The flux of heat at the ice and ocean surfaces are derived
      !!       from semi-empirical ( or bulk ) formulae which relate the flux to 
      !!       the properties of the surface and of the lower atmosphere. Here, we
      !!       follow the work of Oberhuber, 1988   
      !!
      !!  ** Action  :   call flx_blk_albedo to compute ocean and ice albedo 
      !!          computation of snow precipitation
      !!          computation of solar flux at the ocean and ice surfaces
      !!          computation of the long-wave radiation for the ocean and sea/ice
      !!          computation of turbulent heat fluxes over water and ice
      !!          computation of evaporation over water
      !!          computation of total heat fluxes sensitivity over ice (dQ/dT)
      !!          computation of latent heat flux sensitivity over ice (dQla/dT)
      !!
      !! History :
      !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
      !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
         &                          psst      ! Sea Surface Temperature 

      !! * Local variables
      INTEGER  ::             &
         ji, jj, jt        ,  &  ! dummy loop indices
         indaet            ,  &  !  = -1, 0, 1 for odd, normal and leap years resp.
         iday              ,  &  ! integer part of day
         indxb             ,  &  ! index for budyko coefficient
         indxc                   ! index for cloud depth coefficient

      REAL(wp)  ::            & 
         zalat , zclat     ,  &  ! latitude in degrees 
         zmt1, zmt2, zmt3  ,  &  ! tempory air temperatures variables
         ztatm3, ztatm4    ,  &  ! power 3 and 4 of air temperature
         z4tatm3           ,  &  ! 4 * ztatm3
         zcmue             ,  &  ! cosine of local solar altitude
         zcmue2            ,  &  ! root of zcmue1 
         zscmue            ,  &  ! square-root of zcmue1 
         zpcmue            ,  &  ! zcmue1**1.4
         zdecl             ,  &  ! solar declination
         zsdecl , zcdecl   ,  &  ! sine and cosine of solar declination 
         zalbo             ,  &  ! albedo of sea-water
         zalbi             ,  &  ! albedo of ice
         ztamr             ,  &  ! air temperature minus triple point of water (rtt)
         ztaevbk           ,  &  ! part of net longwave radiation
         zevi , zevo       ,  &  ! vapour pressure of ice and ocean 
         zind1,zind2,zind3 ,  &  ! switch for testing the values of air temperature
         zinda             ,  &  ! switch for testing the values of sea ice cover
         zpis2             ,  &  ! pi / 2
         z2pi                    ! 2 * pi 

      REAL(wp)  ::            & 
         zxday             ,  &  ! day of year
         zdist             ,  &  ! distance between the sun and the earth during the year
         zdaycor           ,  &  ! corr. factor to take into account the variation of 
         !                       ! zday when calc. the solar rad.    
         zesi, zeso        ,  &  ! vapour pressure of ice and ocean at saturation
         zesi2             ,  &  ! root of zesi 
         zqsato            ,  &  ! humidity close to the ocean surface (at saturation)   
         zqsati            ,  &  ! humidity close to the ice surface (at saturation) 
         zqsati2           ,  &  ! root of  zqsati 
         zdesidt           ,  &  ! derivative of zesi, function of ice temperature
         zdteta            ,  &  ! diff. betw. sst and air temperature
         zdeltaq           ,  &  ! diff. betw. spec. hum. and hum. close to the surface
         ztvmoy, zobouks   ,  &  ! tempory scalars
         zpsims, zpsihs, zpsils, zobouku, zxins, zpsimu ,  & 
         zpsihu, zpsilu, zstab,zpsim, zpsih, zpsil      ,  & 
         zvatmg, zcmn, zchn, zcln, zcmcmn, zdenum       ,  & 
         zdtetar, ztvmoyr, zlxins, zcmn2, zchcm, zclcm , zcoef

      REAL(wp)  ::            & 
         zrhova            ,  &  ! air density per wind speed
         zcsho , zcleo     ,  &  ! transfer coefficient over ocean
         zcshi , zclei     ,  &  ! transfer coefficient over ice-free
         zrhovacleo        ,  &  ! air density per wind speed per transfer coef.
         zrhovacsho, zrhovaclei, zrhovacshi, & 
         ztice3            ,  &  ! power 3 of ice temperature
         zticemb, zticemb2 ,  &  ! tempory air temperatures variables
         zdqlw_ice         ,  &  ! sensitivity of long-wave flux over ice
         zdqsb_ice         ,  &  ! sensitivity of sensible heat flux over ice
         zdqla_ice         ,  &  ! sensitivity of latent heat flux over ice
         zdl, zdr                ! fractionnal part of latitude
      REAL(wp), DIMENSION(jpi,jpj) :: & 
         zpatm            ,  &   ! atmospheric pressure
         zqatm            ,  &   ! specific humidity
         zes              ,  &   ! vapour pressure at saturation
         zev, zevsqr      ,  &   ! vapour pressure and his square-root
         zrhoa            ,  &   ! air density
         ztatm            ,  &   ! air temperature in Kelvins
         zfrld            ,  &   ! fraction of sea ice cover 
         zcatm1           ,  &   ! fraction of cloud
         zcldeff                 ! correction factor to account cloud effect
      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zalbocsd         ,  &   ! albedo of ocean
         zalboos          ,  &   ! albedo of ocean under overcast sky
         zalbics          ,  &   ! albedo of ice under clear sky
         zalbios          ,  &   ! albedo of ice under overcast sky
         zalbomu          ,  &   ! albedo of ocean when zcmue is 0.4
         zqsro            ,  &   ! solar radiation over ocean
         zqsrics          ,  &   ! solar radiation over ice under clear sky
         zqsrios          ,  &   ! solar radiation over ice under overcast sky
         zcldcor          ,  &   ! cloud correction
         zlsrise, zlsset  ,  &   ! sunrise and sunset
         zlmunoon         ,  &   ! local noon solar altitude
         zdlha            ,  &   ! length of the ninstr segments of the solar day
         zps              ,  &   ! sine of latitude per sine of solar decli.
         zpc                     ! cosine of latitude per cosine of solar decli. 

      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zqlw_oce         ,  &   ! long-wave heat flux over ocean
         zqlw_ice         ,  &   ! long-wave heat flux over ice
         zqla_oce         ,  &   ! latent heat flux over ocean
         zqla_ice         ,  &   ! latent heat flux over ice
         zqsb_oce         ,  &   ! sensible heat flux over ocean
         zqsb_ice                ! sensible heat flux over ice
 
      REAL(wp), DIMENSION(jpi,jpj,jpintsr) ::    &
         zlha             ,  &   ! local hour angle
         zalbocs          ,  &   ! tempory var. of ocean albedo under clear sky
         zsqsro           ,  &   ! tempory var. of solar rad. over ocean 
         zsqsrics         ,  &   ! temp. var. of solar rad. over ice under clear sky
         zsqsrios                ! temp. var. of solar rad. over ice under overcast sky
      !!---------------------------------------------------------------------

      !---------------------
      !   Initilization    !
      !---------------------
#if ! defined key_ice_lim
      tn_ice(:,:) = psst(:,:)
#endif

      !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
      !  and the correction factor for taking into account  the effect of clouds 
      !------------------------------------------------------
      IF( lbulk_init ) THEN
         DO jj = 1, jpj  
            DO ji = 1 , jpi
               zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
               zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
               indxb          = 1 + INT( zalat ) 
               !  correction factor to account for the effect of clouds 
               sbudyko(ji,jj) = budyko(indxb)  
               indxc          = 1 + INT( zclat )  
               zdl            = zclat - INT( zclat ) 
               zdr            = 1.0 - zdl
               stauc(ji,jj)   = zdr * tauco( indxc ) + zdl * tauco( indxc + 1 ) 
            END DO
         END DO
         IF( nleapy == 1 ) THEN
            yearday = 366.e0
         ELSE IF( nleapy == 0 ) THEN
            yearday = 365.e0
         ELSEIF( nleapy == 30) THEN
            yearday = 360.e0
         ENDIF
         lbulk_init = .FALSE.

         ! See test on line 575 
         if (lwp) print*, '    '
         if (lwp) print*, 'ATTENTION, Ch and Ce are multiplied by 1.2'  !DL
         if (lwp) print*, 'ATTENTION, albedo clear sky multiplied by 1.2'  !DL
      ENDIF

      zqlw_oce(:,:) = 0.e0
      zqla_oce(:,:) = 0.e0
      zqsb_oce(:,:) = 0.e0
      zqlw_ice(:,:) = 0.e0
      zqla_ice(:,:) = 0.e0
      zqsb_ice(:,:) = 0.e0

      zpis2       = rpi / 2.
      z2pi        = 2. * rpi

 !CDIR NOVERRCHK
      DO jj = 1, jpj
 !CDIR NOVERRCHK
         DO ji = 1, jpi

            ztatm (ji,jj) = 273.15 + tatm  (ji,jj)  !  air temperature in Kelvins 
            zcatm1(ji,jj) = 1.0    - catm  (ji,jj)  !  fractional cloud cover
            zfrld (ji,jj) = 1.0    - freeze(ji,jj)  !  fractional sea ice cover
            zpatm(ji,jj)  = 101000.               !  pressure 
      
            !  Computation of air density, obtained from the equation of state for dry air. 
            zrhoa(ji,jj) = zpatm(ji,jj) / ( 287.04 * ztatm(ji,jj) )
      
            !  zes : Saturation water vapour
            ztamr = ztatm(ji,jj) - rtt
            zmt1  = SIGN( 17.269, ztamr )
            zmt2  = SIGN( 21.875, ztamr )
            zmt3  = SIGN( 28.200, -ztamr )
            zes(ji,jj) = 611.0 * EXP (  ABS( ztamr ) * MIN ( zmt1, zmt2 )   &
               &                      / ( ztatm(ji,jj) - 35.86  + MAX( zzero, zmt3 ) ) )

            !  zev : vapour pressure  (hatm is relative humidity)  
            zev(ji,jj)   = hatm(ji,jj) * zes(ji,jj) 
            !  square-root of vapour pressure
!CDIR NOVERRCHK
            zevsqr(ji,jj) = SQRT( zev(ji,jj) * 0.01 )
            !  zqapb  : specific humidity 
            zqatm(ji,jj) = 0.622 * zev(ji,jj) / ( zpatm(ji,jj) - 0.378 * zev(ji,jj) )


            !----------------------------------------------------
            !   Computation of snow precipitation (Ledley, 1985) |
            !----------------------------------------------------

            zmt1  =   253.0 - ztatm(ji,jj)
            zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0 
            zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0
            zind1 = MAX( zzero, SIGN( zone, zmt1 ) )
            zind2 = MAX( zzero, SIGN( zone, zmt2 ) )
            zind3 = MAX( zzero, SIGN( zone, zmt3 ) )
            ! total precipitation
            tprecip(ji,jj) = watm(ji,jj)
            ! solid  (snow) precipitation
            sprecip(ji,jj) = tprecip(ji,jj) *       &
               &             (           zind1      &
               &               + ( 1.0 - zind1 ) * ( zind2 * ( 0.5 + zmt2 ) + ( 1.0 - zind2 ) *  zind3 * zmt3 ) ) 
         END DO
      END DO

      !----------------------------------------------------------
      !   Computation of albedo (need to calculates heat fluxes)|
      !-----------------------------------------------------------
      
      CALL flx_blk_albedo( zalbios, zalboos, zalbics, zalbomu )

      !-------------------------------------
      !   Computation of solar irradiance. |
      !----------------------------------------
      indaet   = 1  
      !  compution of the day of the year at which the fluxes have to be calculate 
      !--The date corresponds to the middle of the time step.
      zxday=nday_year + rdtbs2/rday

      iday   = INT( zxday )

      IF(ln_ctl) CALL prt_ctl_info('declin : iday ', ivar1=iday, clinfo2=' nfbulk= ', ivar2=nfbulk)

      !   computation of the solar declination, his sine and his cosine
      CALL flx_blk_declin( indaet, iday, zdecl )
      
      zdecl    = zdecl * rad
      zsdecl   = SIN( zdecl )
      zcdecl   = COS( zdecl )
      
      !  correction factor added for computation of shortwave flux to take into account the variation of
      !  the distance between the sun and the earth during the year (Oberhuber 1988)
      zdist    = zxday * z2pi / yearday
      zdaycor  = 1.0 + 0.0013 * SIN( zdist ) + 0.0342 * COS( zdist )

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            !  product of sine of latitude and sine of solar declination
            zps     (ji,jj) = SIN( gphit(ji,jj) * rad ) * zsdecl
            !  product of cosine of latitude and cosine of solar declination
            zpc     (ji,jj) = COS( gphit(ji,jj) * rad ) * zcdecl
            !  computation of the both local time of sunrise and sunset
            zlsrise (ji,jj) = ACOS( - SIGN( zone, zps(ji,jj) ) * MIN( zone, SIGN( zone, zps(ji,jj) )  &
               &                     * ( zps(ji,jj) / zpc(ji,jj) ) ) ) 
            zlsset  (ji,jj) = - zlsrise(ji,jj)
            !  dividing the solar day into jpintsr segments of length zdlha
            zdlha   (ji,jj) = ( zlsrise(ji,jj) - zlsset(ji,jj) ) / REAL( jpintsr )
            !  computation of the local noon solar altitude
            zlmunoon(ji,jj) = ASIN ( ( zps(ji,jj) + zpc(ji,jj) ) ) / rad
            
            !  cloud correction taken from Reed (1977) (imposed lower than 1)
            zcldcor (ji,jj) = MIN( zone, ( 1.e0 - 0.62 * catm(ji,jj) + 0.0019 * zlmunoon(ji,jj) ) )
         END DO
      END DO

         !  Computation of solar heat flux at each time of the day between sunrise and sunset. 
         !  We do this to a better optimisation of the code 
         !------------------------------------------------------       

!CDIR NOVERRCHK   
      DO jt = 1, jpintsr   
         zcoef = FLOAT( jt ) - 0.5
!CDIR NOVERRCHK     
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               !  local hour angle
               zlha (ji,jj,jt) = COS ( zlsrise(ji,jj) - zcoef * zdlha(ji,jj) )

               ! cosine of local solar altitude
               zcmue              = MAX ( zzero ,   zps(ji,jj) + zpc(ji,jj) * zlha (ji,jj,jt)  )
               zcmue2             = 1368.0 * zcmue * zcmue
               zscmue             = SQRT ( zcmue )
               zpcmue             = zcmue**1.4
               ! computation of sea-water albedo (Payne, 1972)
!DL               zalbocs(ji,jj,jt)  = 0.05 / ( 1.1 * zpcmue + 0.15 )
!DL test mars 2013, *1.2,  autre modif dans albedo.F90
               zalbocs(ji,jj,jt)  = 1.2 * 0.05 / ( 1.1 * zpcmue + 0.15 )
               zalbo              = zcatm1(ji,jj) * zalbocs(ji,jj,jt) + catm(ji,jj) * zalboos(ji,jj)
               ! solar heat flux absorbed at ocean surfaces (Zillman, 1972)
               zevo               = zev(ji,jj) * 1.0e-05
               zsqsro(ji,jj,jt)   =  ( 1.0 - zalbo ) * zdlha(ji,jj) * zcmue2                &
                                   / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
               !  solar heat flux absorbed at sea/ice surfaces 
               !  Formulation of Shine and Crane, 1984 adapted for high albedo surfaces 

               !  For clear sky        
               zevi               = zevo
               zalbi              = zalbics(ji,jj)
               zsqsrics(ji,jj,jt) =  ( 1.0 - zalbi ) * zdlha(ji,jj) * zcmue2                &
                  &                / ( ( 1.0 + zcmue ) * zevi + 1.2 * zcmue + 0.0455 )

               ! For overcast sky
               zalbi              = zalbios(ji,jj)
               zsqsrios(ji,jj,jt) = zdlha(ji,jj) *                                                           &
                  &                 ( ( 53.5 + 1274.5 * zcmue      ) *  zscmue * ( 1.0 - 0.996  * zalbi ) )  &
                  &                 / (  1.0 + 0.139  * stauc(ji,jj) *           ( 1.0 - 0.9435 * zalbi ) )
            END DO
         END DO
      END DO


      !  Computation of daily (between sunrise and sunset) solar heat flux absorbed 
      !  at the ocean and snow/ice surfaces.
      !--------------------------------------------------------------------

      zalbocsd(:,:) = 0.e0
      zqsro   (:,:) = 0.e0
      zqsrics (:,:) = 0.e0
      zqsrios (:,:) = 0.e0

      DO jt = 1, jpintsr 
#   if defined key_vectopt_loop && ! defined key_autotasking
         DO ji = 1, jpij  
            zalbocsd(ji,1) = zalbocsd(ji,1) + zdlha   (ji,1) * zalbocs(ji,1,jt)   &
               &                                             / MAX( 2.0 * zlsrise(ji,1) , zeps0 )
            zqsro   (ji,1) = zqsro   (ji,1) + zsqsro  (ji,1,jt)
            zqsrics (ji,1) = zqsrics (ji,1) + zsqsrics(ji,1,jt)
            zqsrios (ji,1) = zqsrios (ji,1) + zsqsrios(ji,1,jt)
         END DO
#  else
         DO jj = 1, jpj
            DO ji = 1, jpi  
               zalbocsd(ji,jj) = zalbocsd(ji,jj) + zdlha(ji,jj) * zalbocs(ji,jj,jt)   &
                  &                                              / MAX( 2.0 * zlsrise(ji,jj) , zeps0 )
               zqsro  (ji,jj)  = zqsro   (ji,jj) + zsqsro  (ji,jj,jt)
               zqsrics(ji,jj)  = zqsrics (ji,jj) + zsqsrics(ji,jj,jt)
               zqsrios(ji,jj)  = zqsrios (ji,jj) + zsqsrios(ji,jj,jt)
            END DO
         END DO
#  endif
      END DO

      DO jj = 1, jpj
         DO ji = 1, jpi 

            !------------------------------------------- 
            !  Computation of shortwave radiation.
            !-------------------------------------------

            ! the solar heat flux absorbed at ocean and snow/ice surfaces
            !------------------------------------------------------------

            ! For ocean
            qsr_oce(ji,jj) = srgamma * zcldcor(ji,jj) * zqsro(ji,jj) / z2pi
            zinda          = SIGN( zone , -( -0.5 - zfrld(ji,jj) ) )
            zinda          = 1.0 - MAX( zzero , zinda )
            qsr_oce(ji,jj) = ( 1.- zinda ) * qsr_oce(ji,jj)

            ! For snow/ice 
            qsr_ice(ji,jj) = ( zcatm1(ji,jj) * zqsrics(ji,jj) + catm(ji,jj) * zqsrios(ji,jj) ) / z2pi


            ! Taking into account the ellipsity of the earth orbit
            !-----------------------------------------------------

            qsr_ice(ji,jj) = qsr_ice(ji,jj) * zdaycor
            qsr_oce(ji,jj) = qsr_oce(ji,jj) * zdaycor

            !  fraction of net shortwave radiation which is not absorbed in the 
            !  thin surface layer and penetrates inside the ice cover 
            !  ( Maykut and Untersteiner, 1971 ; Elbert anbd Curry, 1993 )
            !------------------------------------------------------------------

            fr1_i0(ji,jj) = 0.18  * zcatm1(ji,jj) + 0.35 * catm(ji,jj) 
            fr2_i0(ji,jj) = 0.82  * zcatm1(ji,jj) + 0.65 * catm(ji,jj)

            !---------------------------------------------------------------------------
            !   Computation of long-wave radiation  ( Berliand 1952 ; all latitudes )
            !---------------------------------------------------------------------------

            ! tempory variables
            ztatm3         = ztatm(ji,jj) * ztatm(ji,jj) * ztatm(ji,jj)
            ztatm4         = ztatm3 * ztatm(ji,jj)
            z4tatm3        = 4. * ztatm3
            zcldeff(ji,jj) = 1.0 - sbudyko(ji,jj) * catm(ji,jj) * catm(ji,jj)    
            ztaevbk        = ztatm4 * zcldeff(ji,jj) * ( 0.39 - 0.05 * zevsqr(ji,jj) ) 

            !  Long-Wave for Ice
            !----------------------
            zqlw_ice(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( tn_ice(ji,jj) - ztatm(ji,jj) ) ) 

            !  Long-Wave for Ocean
            !-----------------------
            zqlw_oce(ji,jj) = - emic * stefan * ( ztaevbk + z4tatm3 * ( psst  (ji,jj) - ztatm(ji,jj) ) ) 

         END DO
      END DO

      !----------------------------------------
      !  Computation of turbulent heat fluxes  ( Latent and sensible ) 
      !----------------------------------------        

      !CDIR NOVERRCHK
      DO jj = 2 , jpjm1
!ib   DO jj = 1 , jpj
         !CDIR NOVERRCHK
         DO  ji = 1, jpi

            !  Turbulent heat fluxes over water
            !----------------------------------

            ! zeso     : vapour pressure at saturation of ocean
            ! zqsato   : humidity close to the ocean surface (at saturation)
            zeso          =  611.0 * EXP ( 17.2693884 * ( psst(ji,jj) - rtt ) * tmask(ji,jj,1) / ( psst(ji,jj) - 35.86 ) )
            zqsato        = ( 0.622 * zeso ) / ( zpatm(ji,jj) - 0.378 * zeso )

            !  Drag coefficients from Large and Pond (1981,1982)
            !---------------------------------------------------
    
            !  Stability parameters
            zdteta         = psst(ji,jj) - ztatm(ji,jj)
            zdeltaq        = zqatm(ji,jj) - zqsato
            ztvmoy         = ztatm(ji,jj) * ( 1. + 2.2e-3 * ztatm(ji,jj) * zqatm(ji,jj) )
            zdenum         = MAX( vatm(ji,jj) * vatm(ji,jj) * ztvmoy, zeps )
!i
!i          if( zdenum == 0.e0 ) then
!i               write(numout,*) 'flxblk  zdenum=0 ', ji,jj
!i               zdenum = zeps
!i          endif
!i
            zdtetar        = zdteta / zdenum
            ztvmoyr        = ztvmoy * ztvmoy * zdeltaq / zdenum
            
            ! For stable atmospheric conditions
            zobouks        = -70.0 * 10. * ( zdtetar + 3.2e-3 * ztvmoyr )
            zobouks        = MAX( zzero , zobouks )
            zpsims         = -7.0 * zobouks
            zpsihs         =  zpsims
            zpsils         =  zpsims

            !  For unstable atmospheric conditions
            zobouku        = -100.0 * 10.0 * ( zdtetar + 2.2e-3 * ztvmoyr )
            zobouku        = MIN( zzero , zobouku )
            zxins          = ( 1. - 16. * zobouku )**0.25
            zlxins         = LOG( ( 1. + zxins * zxins ) / 2. )
            zpsimu         = 2. * LOG( ( 1 + zxins ) / 2. )  + zlxins - 2. * ATAN( zxins ) + zpis2
            zpsihu         = 2. * zlxins
            zpsilu         = zpsihu

            ! computation of intermediate values
            zstab          = MAX( zzero , SIGN( zone , zdteta ) )
            zpsim          = zstab * zpsimu + (1.0 - zstab ) * zpsims
            zpsih          = zstab * zpsihu + (1.0 - zstab ) * zpsihs
            zpsil          = zpsih
            
            zvatmg         = MAX( 0.032 * 1.5e-3 * vatm(ji,jj) * vatm(ji,jj) / grav, zeps )
!i
!!          if( zvatmg == 0.e0 ) then
!!               write(numout,*) 'flxblk  zvatmg=0 ', ji,jj
!!               zvatmg = zeps
!!          endif
!i

            zcmn           = vkarmn / LOG ( 10. / zvatmg )
            zcmn2          = zcmn * zcmn
            zchn           = 0.0327 * zcmn
            zcln           = 0.0346 * zcmn
            zcmcmn         = 1 / ( 1 - zcmn * zpsim / vkarmn )
            zchcm          = zcmcmn / ( 1 - zchn * zpsih / ( vkarmn * zcmn ) )
            zclcm          = zchcm


            !  Transfer cofficient zcsho and zcleo over ocean according to Large and Pond (1981,1982)
            !-------------------------------------------------------------- 
            zcsho          = zchn * zchcm
            zcleo          = zcln * zclcm
 
!DL -------------------------------------
!Test to get more ice in GSL, jan 2013
!DL -------------------------------------
           zcsho = zcsho*1.2
           zcleo = zcsho*1.2
!DL -------------------------------------

            !   Computation of sensible and latent fluxes over Ocean 
            !----------------------------------------------------------------

            !  computation of intermediate values
            zrhova         = zrhoa(ji,jj) * vatm(ji,jj)
            zrhovacsho     = zrhova * zcsho
            zrhovacleo     = zrhova * zcleo

            ! sensible heat flux
            zqsb_oce(ji,jj) = zrhovacsho * 1004.0  * ( psst(ji,jj) - ztatm(ji,jj) )  
         
            !  latent heat flux 
            zqla_oce(ji,jj) = MAX(0.e0, zrhovacleo * 2.5e+06 * ( zqsato      - zqatm(ji,jj) ) )
               
            !  Calculate evaporation over water. (kg/m2/s)
            !-------------------------------------------------
            evap(ji,jj)    = zqla_oce(ji,jj) / cevap
               
               
            !  Turbulent heat fluxes over snow/ice.
            !--------------------------------------------------
            
            !  zesi     : vapour pressure at saturation of ice
            !  zqsati   : humidity close to the ice surface (at saturation)
            zesi           =  611.0 * EXP ( 21.8745587 * tmask(ji,jj,1)   &   ! tmask needed to avoid overflow in the exponential
               &                                       * ( tn_ice(ji,jj) - rtt ) / ( tn_ice(ji,jj) - 7.66 ) )
            zqsati         = ( 0.622 * zesi ) / ( zpatm(ji,jj) - 0.378 * zesi )
               
            !  computation of intermediate values
            zticemb        = ( tn_ice(ji,jj) - 7.66 )
            zticemb2       = zticemb * zticemb  
            ztice3         = tn_ice(ji,jj) * tn_ice(ji,jj) * tn_ice(ji,jj)
            zqsati2        = zqsati * zqsati
            zesi2          = zesi * zesi
            zdesidt        = zesi * ( 9.5 * LOG( 10.0 ) * ( rtt - 7.66 )  / zticemb2 )
               
            !  Transfer cofficient zcshi and zclei over ice. Assumed to be constant Parkinson 1979 ; Maykut 1982
            !--------------------------------------------------------------------
            zcshi          = 1.75e-03
            zclei          = zcshi
               
            !  Computation of sensible and latent fluxes over ice
            !----------------------------------------------------------------
               
            !  computation of intermediate values
            zrhova          = zrhoa(ji,jj) * vatm(ji,jj)
            zrhovacshi      = zrhova * zcshi * 2.834e+06
            zrhovaclei      = zrhova * zclei * 1004.0
            
            !  sensible heat flux
            zqsb_ice(ji,jj) = zrhovaclei * ( tn_ice(ji,jj) - ztatm(ji,jj) )
            
            !  latent heat flux 
            zqla_ice(ji,jj) = zrhovacshi * ( zqsati        - zqatm(ji,jj) )
            qla_ice (ji,jj) = MAX(0.e0, zqla_ice(ji,jj) )
              
            !  Computation of sensitivity of non solar fluxes (dQ/dT)
            !---------------------------------------------------------------
               
            !  computation of long-wave, sensible and latent flux sensitivity
            zdqlw_ice       = 4.0 * emic * stefan * ztice3
            zdqsb_ice       = zrhovaclei
            zdqla_ice       = zrhovacshi * ( zdesidt * ( zqsati2 / zesi2 ) * ( zpatm(ji,jj) / 0.622 ) )         
            
            !  total non solar sensitivity
            dqns_ice(ji,jj) = -( zdqlw_ice + zdqsb_ice + zdqla_ice ) 
            
            ! latent flux sensitivity
            dqla_ice(ji,jj) = zdqla_ice
            
         END DO
      END DO

      ! total non solar heat flux over ice
      qnsr_ice(:,:) = zqlw_ice(:,:) - zqsb_ice(:,:) - zqla_ice(:,:)
      ! total non solar heat flux over water 
      qnsr_oce(:,:) = zqlw_oce(:,:) - zqsb_oce(:,:) - zqla_oce(:,:)

      ! solid precipitations ( kg/m2/day -> kg/m2/s)
      tprecip(:,:) = tprecip  (:,:) / rday 
      ! snow  precipitations ( kg/m2/day -> kg/m2/s)
      sprecip(:,:) = sprecip  (:,:) / rday  
!i
      CALL lbc_lnk( qsr_oce (:,:) , 'T', 1. )
      CALL lbc_lnk( qnsr_oce(:,:) , 'T', 1. )
      CALL lbc_lnk( qsr_ice (:,:) , 'T', 1. )
      CALL lbc_lnk( qnsr_ice(:,:) , 'T', 1. )
      CALL lbc_lnk( qla_ice (:,:) , 'T', 1. )
      CALL lbc_lnk( dqns_ice(:,:) , 'T', 1. )
      CALL lbc_lnk( dqla_ice(:,:) , 'T', 1. )
      CALL lbc_lnk( fr1_i0  (:,:) , 'T', 1. )
      CALL lbc_lnk( fr2_i0  (:,:) , 'T', 1. )
      CALL lbc_lnk( tprecip (:,:) , 'T', 1. )
      CALL lbc_lnk( sprecip (:,:) , 'T', 1. )
      CALL lbc_lnk( evap    (:,:) , 'T', 1. )
!i
!i
      qsr_oce (:,:) = qsr_oce (:,:)*tmask(:,:,1)
      qnsr_oce(:,:) = qnsr_oce(:,:)*tmask(:,:,1)
      qsr_ice (:,:) = qsr_ice (:,:)*tmask(:,:,1)
      qnsr_ice(:,:) = qnsr_ice(:,:)*tmask(:,:,1)
      qla_ice (:,:) = qla_ice (:,:)*tmask(:,:,1)
      dqns_ice(:,:) = dqns_ice(:,:)*tmask(:,:,1)
      dqla_ice(:,:) = dqla_ice(:,:)*tmask(:,:,1)
      fr1_i0  (:,:) = fr1_i0  (:,:)*tmask(:,:,1)
      fr2_i0  (:,:) = fr2_i0  (:,:)*tmask(:,:,1)
      tprecip (:,:) = tprecip (:,:)*tmask(:,:,1)
      sprecip (:,:) = sprecip (:,:)*tmask(:,:,1)
      evap    (:,:) = evap    (:,:)*tmask(:,:,1)
!i

   END SUBROUTINE flx_blk


   SUBROUTINE flx_blk_declin( ky, kday, pdecl )
      !!---------------------------------------------------------------------------
      !!               ***  ROUTINE flx_blk_declin  ***
      !!          
      !! ** Purpose :   Computation of the solar declination for the day
      !!	        kday ( in decimal degrees ).
      !!       
      !! ** Method  :
      !!
      !! History :
      !!         original    : 01-04 (LIM)
      !!         addition    : 02-08 (C. Ethe, G. Madec)
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   &
         ky  ,        &  ! = -1, 0, 1 for odd, normal and leap years resp.
         kday            ! day of the year ( kday = 1 on january 1)
      REAL(wp), INTENT(out) ::  &
         pdecl            ! solar declination

      !! * Local variables
      REAL(wp) ::                & 
         zday                ,   &  ! corresponding day of type year (cf. ky)
         zp1, zp2, zp3, zp4         ! temporary scalars
      !!---------------------------------------------------------------------
      
      zday = FLOAT( kday ) 
      
      IF( ky == 1 )  THEN 
         zday = zday - 0.5
      ELSEIF( ky == 3 )  THEN
         zday = zday - 1.
      ELSE 
         zday = REAL( kday )
      ENDIF
      
      zp1 = rpi * ( 2.0 * zday - 367.0 ) / yearday
      zp2 = 2. * zp1
      zp3 = 3. * zp1
      zp4 = 4. * zp1
      
      pdecl  = a0                                                                      &
         &   + a1 * COS( zp1 ) + a2 * COS( zp2 ) + a3 * COS( zp3 ) + a4 * COS( zp4 )   &
         &   + b1 * SIN( zp1 ) + b2 * SIN( zp2 ) + b3 * SIN( zp3 ) + b4 * SIN( zp4 )
      
   END SUBROUTINE flx_blk_declin

#else
   !!----------------------------------------------------------------------
   !!   Default option :           Empty module                     NO bulk
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flx_blk            ! Empty routine
   END SUBROUTINE flx_blk
#endif

   !!======================================================================
END MODULE flxblk
