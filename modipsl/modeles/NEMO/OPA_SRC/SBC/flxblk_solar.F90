!!DB 2008.12.04
!!Modified routine to remove possible conflicts with lim-calls and
!!for easier use in BGCM_01
!!-1- changed qsr_oce? to ad_qsr_oce?
!!-2- the variable ad_qsr_oce() now contains a true daily ave
!!    (for BGCM_01 use)
!!-3- #ifdef key_lim_ice ==> ad_qsr_oce(:,:)=frld(:,:)*ad_qsr_oce(:,:)
!!    where frld is the fraction of open water. Note that this is not
!!    meant to be a perfect accounting of ice coverage, but rather
!!    something likely better than nothing. It will set qsr=0 for
!!    completely ice-covered cells. 
!!DB 2009.08.17 
!!Re -2- above: I adjusted the already-daily-averaged insolation value to be
!!daily-averaged which reduces the insolation. This turns out to be not
!!so bad as it looks like an average cloudiness correction ===> Keep as is
!!

MODULE flxblk_solar
   !!======================================================================
   !!                       ***  MODULE  flxblk_solar  ***
   !!
   !!=================== SUMMARY =========================================
   !!Purpose:         Stand alone Short wave flux calcuation
   !!
   !!Inputs:         None / CMC humidity cloud cover and air temp
   !!
   !!Outputs:             qsr_oce,     &        !: solar flux over ocean (daily average)
   !!                     qsr_oce1              !: solar flux over ocean (instantaneous) &
   !!                                                available  only with cmc
   !! 
   !! Usage:   -compilie with key_bulk_solar. No other keys required
   !!          -SUBROUTINE flx_blk_solar( kt) must be called at each timestep
   !!          -Access qsr_oce and qsr_oce1 with "USE flxblk_solar"
   !!================= DESCRIPTION =======================================
   !!
   !!    SUBROUTINE flx_blk_solar( kt) computes shortwave radiation based
   !!     on the Zillman (1972). It has been stripped out of flxblk.F90
   !!     which uses this formula for the open ocean and something else for
   !!     ice (and as a result is tied to runing opa with lim and bulk
   !!     forcing).
   !!   
   !!    The calculation from flxblk.F90 computes a daily average (kept
   !!    here as qsr_oce) and is independent of local time.  This module
   !!    also offers an instantanous flux based on the "real_time"
   !!    varialble read from the cmc forcing.  read_cmc_flx opens the cmc
   !!    netcdf forcing files (based on flx_bulk_cmc.h90) and reads the
   !!    nearest 3h forcast).  For this part, the nit000 is assumed to
   !!    coincide with the first fields in the cmc files. The correct
   !!    yearday is provided through ndate0 in namelist.  The reading of
   !!    the "real_time" in the cmc files is accomplished with
   !!    ncdf_readdate, a new subroutine in lib_ncdf, designed to read a
   !!    string of length 19 characters from an array of unlimited
   !!    length. The only time information from the real_time variable that
   !!    is actually used, is jhour which is the local gmt time. jhour is
   !!    used to calculate the appropriate time offset based on the
   !!    longitude of each grid point.
   !!   
   !!    If do_read_cmc logical parameter is set to false. The cmc files
   !!    are not read. The Zillman formula is computed from constant values
   !!    for humidity, cloud cover and air temperate. These are hardwired
   !!    in the code (search for do_read_cmc).  In this scenario qsr_oce1
   !!    is meaningless and only qsr_oce should be used.
   !!   
   !!    Last Modify: 28Sep2008 by A. Drozdowski
   !!   =====================================================================



#if defined key_bulk_solar || defined key_BGCM_01
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE cpl_oce         ! ???
   USE phycst          ! physical constants
   USE daymod

   USE ocfzpt          ! ???
   USE in_out_manager
   USE lbclnk
   USE albedo
   USE prtctl          ! Print control
#ifdef key_ice_lim
   USE ice, ONLY : frld
#endif



   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC flx_blk_solar 

   !! * Module variables
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
      ad_qsr_oce,     &        !: solar flux over ocean (daily average)
      ad_qsr_oce1,    &          !: solar flux over ocean (instantanous)    
      frac_day_length            !DB
   
   INTEGER, PARAMETER  ::   &
      jpintsr = 24          ! number of time step between sunrise and sunset
      !                     ! uses for heat flux computation
   LOGICAL ::   &
      lbulk_init = .TRUE.   ! flag, bulk initialization done or not)

   !! * constants for bulk computation (flx_blk)
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

!*** Control reading from cmc files
   LOGICAL do_read_cmc
!Joel:   parameter (do_read_cmc=.false.)
   parameter (do_read_cmc=.true.) !Joel: for cmc run
      
            
!!! import from blk_oce
   REAL(wp) ::        &
      yearday  ,      &  !: number of days per year
      rdtbs2             !: bulk time step divide by 2

   INTEGER ::         & !!: namdom : space/time domain (namlist)
      nfbulk =  5        !: bulk computation frequency 
      
            
   REAL(wp), PRIVATE, DIMENSION(jpi,jpj) ::   &  !:
!      watm     ,      &  !: precipitation
      tatm     ,      &  !: atmospheric temperature
      hatm     ,      &  !: relative humidity
!      vatm     ,      &  !: wind speed
      catm               !: percent of cloud cover


!import from flx_bulk_cmc: subroutine flx was change to read_cmc_flux
!!DB
   INTEGER :: nflx1, nflx2, num_recsf, ireadf 
      
!A.D: date string to be read from cmc file
   CHARACTER(LEN=19) :: cmc_date      
   REAL(wp) :: jhour  !julian hour of the day               
   
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxblk.F90,v 1.8 2005/09/02 15:45:31 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flx_blk_solar( kt)
      !!---------------------------------------------------------------------------
      !!                     ***  ROUTINE flx_blk  ***
      !!                 
      !!  ** Purpose :   Computation of the soloar heat flux over ocean: Zillman
      !!
      !! History :
      !!   AD: 11Sep08 Modified flx_blk
      !!----------------------------------------------------------------------
      !! * Arguments
      
      !! * Local variables
      INTEGER  ::             &
         kt,ji, jj, jt        ,  &  ! dummy loop indices
         indaet            ,  &  !  = -1, 0, 1 for odd, normal and leap years resp.
         iday              ,  &  ! integer part of day
         indxb             ,  &  ! index for budyko coefficient        
         indxc             , &     ! index for cloud depth coefficient
!A.D:         
         lhour         

      REAL(wp)  ::            & 
         zalat , zclat     ,  &  ! latitude in degrees 
         zmt1, zmt2, zmt3  ,  &  ! tempory air temperatures variables
         ztatm3, ztatm4    ,  &  ! power 3 and 4 of air temperature
         z4tatm3           ,  &  ! 4 * zftatm3
         zcmue             ,  &  ! cosine of local solar altitude
         zcmue2            ,  &  ! root of zcmue1 
         zpcmue            ,  &  ! zcmue1**1.4
         zdecl             ,  &  ! solar declination
         zsdecl , zcdecl   ,  &  ! sine and cosine of solar declination 
         zalbo             ,  &  ! albedo of sea-water
         zalbi             ,  &  ! albedo of ice
         ztamr             ,  &  ! air temperature minus triple point of water (rtt)
         ztaevbk           ,  &  ! part of net longwave radiation
         zevi , zevo       ,  &  ! vapour pressure of ice and ocean 
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
         zdl, zdr                ! fractionnal part of latitude
         
      REAL(wp), DIMENSION(jpi,jpj) :: & 
         zpatm            ,  &   ! atmospheric pressure
         zqatm            ,  &   ! specific humidity
         zes              ,  &   ! vapour pressure at saturation
         zev, zevsqr      ,  &   ! vapour pressure and his square-root
         zrhoa            ,  &   ! air density
         ztatm            ,  &   ! air temperature in Kelvins
         zcatm1                  ! fraction of cloud

      REAL(wp), DIMENSION(jpi,jpj) ::   & 
         zalboos          ,  &   ! albedo of ocean under overcast sky
         zalbics          ,  &   ! albedo of ice under clear sky
         zalbios          ,  &   ! albedo of ice under overcast sky
         zalbomu          ,  &   ! albedo of ocean when zcmue is 0.4
         zqsro            ,  &   ! solar radiation over ocean
         zcldcor          ,  &   ! cloud correction
         zlsrise, zlsset  ,  &   ! sunrise and sunset
         zlmunoon         ,  &   ! local noon solar altitude
         zdlha            ,  &   ! length of the ninstr segments of the solar day
         zps              ,  &   ! sine of latitude per sine of solar decli.
         zpc                     ! cosine of latitude per cosine of solar decli. 

      REAL(wp), DIMENSION(jpi,jpj,jpintsr) ::    &
         zlha             ,  &   ! local hour angle
         zalbocs          ,  &   ! tempory var. of ocean albedo under clear sky
         zsqsro                  ! tempory var. of solar rad. over ocean 

      !!---------------------------------------------------------------------

      !---------------------
      !   Initilization    !
      !---------------------

      !  Determine cloud optical depths as a function of latitude (Chou et al., 1981).
      !  and the correction factor for taking into account  the effect of clouds 
      !------------------------------------------------------
      IF( lbulk_init ) THEN
         DO jj = 1, jpj  
            DO ji = 1 , jpi
               zalat          = ( 90.e0 - ABS( gphit(ji,jj) ) ) /  5.e0
               zclat          = ( 95.e0 -      gphit(ji,jj)   ) / 10.e0
               indxb          = 1 + INT( zalat ) 
               indxc          = 1 + INT( zclat )  
               zdl            = zclat - INT( zclat ) 
               zdr            = 1.0 - zdl
            END DO
         END DO
         IF( nleapy == 1 ) THEN
            yearday = 366.e0
         ELSE IF( nleapy == 0 ) THEN
            yearday = 365.e0
         ELSEIF( nleapy == 30) THEN
            yearday = 360.e0
         ENDIF

!!! import from blk_oce                  
      ! computation of rdtbs2
        IF( nacc == 1 ) THEN
           rdtbs2 = nfbulk * rdtmin * 0.5
        ELSE
           rdtbs2 = nfbulk * rdt * 0.5
        ENDIF
                 
         
         lbulk_init = .FALSE.
      ENDIF


      

!#if defined key_flx_bulk_cmc
      IF (do_read_cmc)  then
!Joel:
       IF(nproc == 0) write(numout,*) PRINT*,'Warning CMC reading flag set to true'
           call read_cmc_flx( kt )
      ELSE
   !*** Setting Default values       
         tatm(:,:)=10.
         catm(:,:)=0.10      
         hatm(:,:)=0.4     
      ENDIF
!#endif      
        
            
      IF( MOD( kt - nit000 , nfbulk ) == 0 ) THEN
         zpis2       = rpi / 2.
         z2pi        = 2. * rpi
   
        
         
         DO jj = 1, jpj
            DO ji = 1, jpi
  
               ztatm (ji,jj) = 273.15 + tatm  (ji,jj)  !  air temperature in Kelvins 
               zcatm1(ji,jj) = 1.0    - catm  (ji,jj)  !  fractional cloud cover
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
               zevsqr(ji,jj) = SQRT( zev(ji,jj) * 0.01 )
               !  zqapb  : specific humidity 
               zqatm(ji,jj) = 0.622 * zev(ji,jj) / ( zpatm(ji,jj) - 0.378 * zev(ji,jj) )
   
   
               !----------------------------------------------------
               !   Computation of snow precipitation (Ledley, 1985) |
               !----------------------------------------------------
   
               zmt1  =   253.0 - ztatm(ji,jj)
               zmt2  = ( 272.0 - ztatm(ji,jj) ) / 38.0 
               zmt3  = ( 281.0 - ztatm(ji,jj) ) / 18.0
  
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
   
         DO jj = 1, jpj
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
   
         DO jt = 1, jpintsr   
            zcoef = FLOAT( jt ) - 0.5
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !  local hour angle
                  zlha (ji,jj,jt) = COS ( zlsrise(ji,jj) - zcoef * zdlha(ji,jj) )
   
                  ! cosine of local solar altitude
                  zcmue              = MAX ( zzero ,   zps(ji,jj) + zpc(ji,jj) * zlha (ji,jj,jt)  )
                  zcmue2             = 1368.0 * zcmue * zcmue
                  zpcmue             = zcmue**1.4
                  ! computation of sea-water albedo (Payne, 1972)
                  zalbocs(ji,jj,jt)  = 0.05 / ( 1.1 * zpcmue + 0.15 )
                  zalbo              = zcatm1(ji,jj) * zalbocs(ji,jj,jt) + catm(ji,jj) * zalboos(ji,jj)
                  
                  ! solar heat flux absorbed at ocean surfaces (Zillman, 1972)
                  zevo               = zev(ji,jj) * 1.0e-05
                  zsqsro(ji,jj,jt)   =  ( 1.0 - zalbo ) * zdlha(ji,jj) * zcmue2                &
                                       / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
!A.D Warning!!!! not appling albedo for testing purposes                  
!                  zsqsro(ji,jj,jt)   =  zdlha(ji,jj) * zcmue2                &
!                                      / ( ( zcmue + 2.7 ) * zevo + 1.085 * zcmue +  0.10 )
               
               END DO
            END DO
         END DO
   
  
         !  Computation of daily (between sunrise and sunset) solar heat flux absorbed 
         !  at the ocean and snow/ice surfaces.
         !--------------------------------------------------------------------
   
         zqsro   (:,:) = 0.e0
         DO jt = 1, jpintsr 
            DO jj = 1, jpj
               DO ji = 1, jpi  
                  zqsro  (ji,jj)  = zqsro   (ji,jj) + zsqsro  (ji,jj,jt)
               END DO
            END DO
         END DO

 
         
                 
         DO jj = 1, jpj
            DO ji = 1, jpi 
   
               !------------------------------------------- 
               !  Computation of shortwave radiation.
               !-------------------------------------------
               ! For ocean
               lhour=mod(int(jhour+glamt(ji,jj)/360.*24.),24)  !compute local hour  

               ad_qsr_oce(ji,jj)  = srgamma*zcldcor(ji,jj)*zqsro(ji,jj) / z2pi
               ad_qsr_oce1(ji,jj) = REAL( jpintsr )*srgamma*zcldcor(ji,jj)*zsqsro(ji,jj,lhour) / z2pi
!!DB -- NB: 
!               frac_day_length(ji,jj) = 24.0*(zlsrise(ji,jj)-zlsset(ji,jj))/z2pi
               frac_day_length(ji,jj) = (zlsrise(ji,jj)-zlsset(ji,jj))/z2pi
                              
!Diagnostic 1
!                if (jj .eq. 10 .AND. ji .eq. 10) then               
!                write(9090+narea,'(I7,F8.4,F8.1,2(F8.4))')kt,jhour,tatm(ji,jj),gphit(ji,jj),glamt(ji,jj)
!                endif
                                
               if ((jj -1 + njmpp) .eq. 80 .AND. (ji -1 + nimpp) .eq. 135) then
!                  write(9081,'(I7,30(F6.1))') kt,srgamma,zcldcor(ji,jj),z2pi,zone,zzero,zdaycor
!                  write(9080,'(I7,30(F6.1))') kt,ad_qsr_oce(ji,jj),zqsro(ji,jj),zsqsro(ji,jj,:)
!DB
                  write(9080,'(I7,I3,2(F9.4),30(F6.1,1x))') kt,lhour,gphit(ji,jj),glamt(ji,jj) &
                   ,ad_qsr_oce(ji,jj),ad_qsr_oce1(ji,jj), zlsrise(ji,jj),zlsset(ji,jj),frac_day_length(ji,jj), &
                   frac_day_length(ji,jj)*ad_qsr_oce(ji,jj)
               endif
!!DB: convert qsr to average over 24h instead of average over day_length
!!2009.08.17 --NB: This correction is in fact an error (see note at top of file) BUT
!!                 it serves as a sort-of cloudiness correction (reduces insolation)
!!                 ====> keep as is
               ad_qsr_oce(ji,jj)  =  ad_qsr_oce(ji,jj)*frac_day_length(ji,jj)
               
   
            END DO
         END DO


#ifdef key_ice_lim
!!DB: use frac open water to adjust insolation
!!Note the adjustment in case frld=0==> ice covered cell
!!This guards against zero insolation with a 1% error otherwise
         ad_qsr_oce(:,:) = ad_qsr_oce(:,:) * (frld(:,:)+0.01)
#endif
   
!!DBG: this causes some zero values which are not desired. 
!!===> do not call
!         CALL lbc_lnk( ad_qsr_oce (:,:) , 'T', 1. )

!!DB -- NB: DO NOT MASK
!         ad_qsr_oce (:,:) = ad_qsr_oce (:,:)*tmask(:,:,1)
!         write(3500+narea,*)-kt, minval(ad_qsr_oce(:,:)),minloc(ad_qsr_oce(:,:))

      ENDIF  
  
   END SUBROUTINE flx_blk_solar
  
   
   
  
   SUBROUTINE flx_blk_declin( ky, kday, pdecl )
      !!---------------------------------------------------------------------------
      !!               ***  ROUTINE flx_blk_declin  ***
      !!          
      !! ** Purpose :   Computation of the solar declination for the day
      !!         kday ( in decimal degrees ).
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
  
     
   SUBROUTINE read_cmc_flx( kt )
       !!---------------------------------------------------------------------
       !!                    ***  ROUTINE flx  ***
       !!        
       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
       !!      to the ocean at each time step.
       !!
       !!       Read several cmc fluxes file
       !!             temperature at 2m   tatm   (K)
       !!             relative humidite   hatm   (%)
 !      !!             wind speed          vatm   (m/s)
 !      !!             precip              watm   (kg/m2/day)
       !!             clouds              catm   (%)
       !!
       !! caution : now, in the opa global model, the net upward water flux is
       !! -------   with mm/day unit.
       !!
       !! History :
       !!        !  91-03  (O. Marti and Ph Dandin)  Original code
       !!        !  92-07  (M. Imbard)
       !!        !  96-11  (E. Guilyardi)  Daily AGCM input files
       !!        !  99-11  (M. Imbard)  NetCDF FORMAT with ioipsl
       !!        !  00-05  (K. Rodgers) Daily Netcdf
       !!   8.5  !  02-09  (C. Ethe and G. Madec)  F90: Free form and MODULE 
       !!----------------------------------------------------------------------
       !! * modules used
 !!DB
       USE oce, only :  perpetual_forcing, ramp
       USE lib_ncdf
 
       !! * arguments
       INTEGER, INTENT( in  ) ::   kt ! ocean time step
 
       !! * Local declarations      
       CHARACTER(len=45)  ::  &
          clname_n = 'tair_3h.nc',        &
          clname_c = 'hum_cloud_3h.nc',   &
          clname_x = 'rain_3h.nc',        &
          clname_w = 'wspd_3h.nc'
 !!DB
       REAL(wp) :: zxy, update_frac, rdt_3h
       INTEGER :: i, j, frame0, status, ndt_3h, ndt_1h
 
 
       !!---------------------------------------------------------------------
 
 
       ! Initialization
       ! --------------
       rdt_3h = 3.0*3600./rdt
       ndt_3h = nint(3.0*3600./rdt)
       ndt_1h = nint(1.0*3600./rdt) 
 !!DB  frame0 is a hardwired offset into the files. Use for the (likely restart) 
 !!    case where  ntau1 = (kt-nit000)/ndt_3h + 1 points to the wrong frame
       frame0 = 0
 
 
       ! 1. first call kt = nit000
       ! -------------------------
       
       IF( kt == nit000 ) THEN
 
 !!DB -- crucial assignments
 !!Use new routine to get #-records
          call ncdf_get_dim_size(clname_x, 'time_counter',num_recsf,status)
          if(status == 0) then
             if(lwp) write(numout,*)'read_cmc_flx: # records in ',clname_x, ' = ', num_recsf
          else
             if(lwp) write(numout,*)'read_cmc_flx: prob with ncdf_get_dim_size ====> STOP'
             stop
          endif
 
          nflx1 =  0
          if(nflx1 .gt. num_recsf) nflx1 = num_recsf
          nflx2 = nflx1 + 1
          if(nflx2 .gt. num_recsf) nflx2 = num_recsf
          ireadf = 0 ! next input timestep rel nit000 = 0 to start
 
 
     ENDIF
 
 
     ! 2. Read (3h) cmc data
     ! ----------------------------
 
     if(kt-nit000 == ireadf) then   !!orig DB statement
 
 
        nflx1 = nflx1 + 1
 !!DBG if using daily OMIP must turn off these checks
        if(nflx1 .gt. num_recsf) nflx1 = num_recsf
        nflx2 = nflx1 + 1
        if(nflx2 .gt. num_recsf) nflx2 = num_recsf
        ireadf = nint((nflx1)*rdt_3h) !! next input timestep rel nit000
        
 !!DB: perpetual_forcing from oce; 
        if(perpetual_forcing > 0) then
           nflx1 = 1
           nflx2 = nflx1 
 !!DBG
           if(lwp .AND. kt==nit000) then 
              write(numout,*)'read_cmc_flx: perpetual_forcing = ', perpetual_forcing
              write(numout,*)'read_cmc_flx: nflx2 = nflx1 = 1'
           endif
        endif
 
 !!DB
        if(lwp) write(numout,*)'read_cmc_flx:  -- reading field ',nflx1, 'at kt = ', kt 
        if(lwp) write(numout,*)'read_cmc_flx:  -- next input at ireadf = ', ireadf
        
        ! Read data
 !!Adjust using frame0 
        ! read temp at 2m (in K)
        call ncdf_read(clname_n, 'air', tatm, -(nflx1+frame0), status)
        if(status /= 0) then
           if(lwp) then
              write(numout,*)'read_cmc_flx: Prob reading tair file, status = ', status
              !               nstop = nstop + 1
           endif
        endif
        ! conversion of temperature Kelvin --> Celsius  [rt0=273.15]
        tatm(:,:) = ( tatm(:,:) - rt0 ) 
 
 !A.D: read date from tair file
        call ncdf_readdate(clname_n, 'real_time', cmc_date, (nflx1+frame0), status)        
        read(cmc_date(12:13),*) jhour                
           if(lwp) then
              write(numout,*)'read_cmc_flx: cmc_date = ', cmc_date,(nflx1+frame0)              
              write(numout,*)'read_cmc_flx: cmc_date (hour) = ', cmc_date(12:13)           
              write(numout,*)'read_cmc_flx: julian hour = ', jhour
           endif

        ! read humidity
        call ncdf_read(clname_c, 'socliohu', hatm, -(nflx1+frame0), status)
        if(status /= 0) then
           if(lwp) then
              write(numout,*)'read_cmc_flx: Prob reading hum_cloud file, status = ', status
              !               nstop = nstop + 1
           endif
        endif
        
        ! read clouds
        call ncdf_read(clname_c, 'socliocl', catm, -(nflx1+frame0), status)
        if(status /= 0) then
           if(lwp) then
              write(numout,*)'read_cmc_flx: Prob reading hum_cloud file, status = ', status
              !               nstop = nstop + 1
           endif
        endif
        
 
        call flush(numout)
          
     ENDIF
 
 
 !!DB: basic structure to do time interpolation every ndt_1h
 !!See tau_forced_cmc for more details
 !!NB: would have to have input the nflx2 fields as well (not done above)
 !!DB Use ndt_3h & ndt_1h
 !!Want a formula that also accounts for nature of ndt_* variables which may not
 !!be integers (depending on rdt), and still gives a consistent zxy upon restart
 !!NB: formula below is incorrect if perpetual_forcing is on so avoid this possibility 
 !      if(mod(kt-nit000,ndt_1h) == 0  .AND.  perpetual_forcing == 0) then
 !         zxy = (float(kt-nit000) - (nflx1-1)*rdt_3h)/rdt_3h
 !         zxy = min(zxy,1.0)  !to be safe
 !         taux(:,:) = (1.-zxy) * taux_dta(:,:,1) + zxy * taux_dta(:,:,2)
 !         tauy(:,:) = (1.-zxy) * tauy_dta(:,:,1) + zxy * tauy_dta(:,:,2)
 !      endif
 
 !       CALL blk(kt)
       
       CALL FLUSH(numout)
       
   END SUBROUTINE read_cmc_flx
      
   
#else
   !!----------------------------------------------------------------------
   !!   Default option :           Empty module                     NO bulk
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flx_blk_solar            ! Empty routine
   END SUBROUTINE flx_blk_solar
#endif

   !!======================================================================
END MODULE flxblk_solar
