MODULE phycst
   !!======================================================================
   !!                    ***  MODULE  phycst  ***
   !!     Definition of of both ocean and ice parameters used in the code
   !!=====================================================================
   !! * Modules used
   USE par_oce          ! ocean parameters
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC phy_cst          ! routine called by inipar.F90

   !! * Shared module variables
   INTEGER, PUBLIC, DIMENSION(12) ::   &  !:
      nbiss = (/ 31, 29, 31, 30, 31, 30,      &  !: number of days per month
         &       31, 31, 30, 31, 30, 31 /) ,  &  !  (leap-year)
      nobis = (/ 31, 28, 31, 30, 31, 30,      &  !: number of days per month
         &       31, 31, 30, 31, 30, 31 /)       !  (365 days a year)
   
   REAL(wp), PUBLIC ::                        &  !:
      rpi = 3.141592653589793_wp           ,  &  !: pi
      rad = 3.141592653589793_wp / 180._wp ,  &  !: conversion from degre into radian
      rsmall = 0.5 * EPSILON( 1. )               !: smallest real computer value
   
   REAL(wp), PUBLIC ::          & !:
      rday = 24.*60.*60.  ,     & !: day (s)
      rsiyea              ,     & !: sideral year (s)
      rsiday              ,     & !: sideral day (s)
      raajj = 365._wp     ,     & !: number of days in one year
      raamo =  12._wp     ,     & !: number of months in one year
      rjjhh =  24._wp     ,     & !: number of hours in one day
      rhhmm =  60._wp     ,     & !: number of minutes in one hour
      rmmss =  60._wp     ,     & !: number of seconds in one minute
      raass               ,     & !: number of seconds in one year
      rmoss               ,     & !: number of seconds in one month
      rjjss               ,     & !: number of seconds in one day
!!!   omega = 7.292115083046061e-5_wp ,  &  !: change the last digit!
      omega               ,    &  !: earth rotation parameter
      ra    = 6371229._wp ,    &  !: earth radius (meter)
      grav  = 9.80665_wp          !: gravity (m/s2)
   
   REAL(wp), PUBLIC ::         &  !:
      rtt      = 273.16_wp  ,  &  !: triple point of temperature (Kelvin)
      rt0      = 273.15_wp  ,  &  !: freezing point of water (Kelvin)
      rt0_snow = 273.15_wp  ,  &  !: melting point of snow  (Kelvin)
      rt0_ice  = 273.05_wp  ,  &  !: melting point of ice   (Kelvin)
      rau0     = 1020._wp   ,  &  !: volumic mass of reference (kg/m3)
      rauw     = 1000._wp   ,  &  !: density of pure water (kg/m3)
      rcp      =    4.e+3_wp,  &  !: ocean specific heat
      ro0cpr                      !: = 1. / ( rau0 * rcp )

   REAL(wp), PUBLIC ::            &  !:
      rcdsn   =   0.22_wp     ,   &  !: conductivity of the snow
      rcdic   =   2.034396_wp ,   &  !: conductivity of the ice
      rcpsn   =   6.9069e+5_wp,   &  !: density times specific heat for snow
      rcpic   =   1.8837e+6_wp,   &  !: volumetric latent heat fusion of sea ice
      xlsn    = 110.121e+6_wp ,   &  !: volumetric latent heat fusion of snow
      xlic    = 300.33e+6_wp  ,   &  !: volumetric latent heat fusion of ice
      xsn     =   2.8e+6      ,   &  !: latent heat of sublimation of snow
      rhoic   = 900._wp       ,   &  !: density of sea ice (kg/m3)
      rhosn   = 330._wp       ,   &  !: density of snow (kg/m3)
      emic    =   0.97_wp     ,   &  !: emissivity of snow or ice
      sice    =   6.0_wp      ,   &  !: salinity of ice (psu)
      soce    =  34.7_wp      ,   &  !: salinity of sea (psu)
      cevap   =   2.5e+6_wp   ,   &  !: latent heat of evaporation (water)
      srgamma =   0.9_wp      ,   &  !: correction factor for solar radiation (Oberhuber, 1974)
      vkarmn  =   0.4_wp      ,   &  !: von Karman constant
      stefan  =   5.67e-8_wp         !: Stefan-Boltzmann constant 
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/phycst.F90,v 1.4 2005/03/27 18:34:48 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!----------------------------------------------------------------------
   
CONTAINS
   
   SUBROUTINE phy_cst
      !!----------------------------------------------------------------------
      !!                       ***  ROUTINE phy_cst  ***
      !!
      !! ** Purpose :   Print model parameters and set and print the constants
      !!
      !! ** Method  :   no
      !!
      !! History :
      !!        !  90-10  (C. Levy - G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  91-12  (M. Imbard)
      !!   8.5  !  02-08  (G. Madec, C. Ethe)  F90, add ice constants 
      !!----------------------------------------------------------------------
      !! * Local variables
      CHARACTER (len=64) ::   cform = "(A9, 3(A13, I7) )" 
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' phy_cst : initialization of ocean parameters and constants'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~'

      ! Ocean Parameters
      ! ----------------
      IF(lwp) THEN
         WRITE(numout,*) '       parameter file'
         WRITE(numout,*)
         WRITE(numout,*) '          dimension of model'
         WRITE(numout,*) '              Local domain      Global domain       Data domain '
         WRITE(numout,cform) '         ','   jpi     : ', jpi, '   jpiglo  : ', jpiglo, '   jpidta  : ', jpidta
         WRITE(numout,cform) '         ','   jpj     : ', jpj, '   jpjglo  : ', jpjglo, '   jpjdta  : ', jpjdta
         WRITE(numout,cform) '         ','   jpk     : ', jpk, '   jpk     : ', jpk   , '   jpkdta  : ', jpkdta
         WRITE(numout,*)      '        ','   jpij    : ', jpij
         WRITE(numout,*)
         WRITE(numout,*) '          mpp local domain info (mpp)'
         WRITE(numout,*) '             jpni    : ', jpni, '   jpreci  : ', jpreci
         WRITE(numout,*) '             jpnj    : ', jpnj, '   jprecj  : ', jprecj
         WRITE(numout,*) '             jpnij   : ', jpnij

         WRITE(numout,*)
         WRITE(numout,*) '          lateral domain boundary condition type : jperio  = ', jperio
         WRITE(numout,*) '          domain island (use in rigid-lid case)  : jpisl   = ', jpisl 
         WRITE(numout,*) '                                                   jpnisl  = ', jpnisl
      ENDIF

      ! Define constants
      ! ----------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '       constants'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          mathematical constant                 rpi = ', rpi

      rsiyea = 365.25 * rday * 2. * rpi / 6.283076
      rsiday = rday / ( 1. + rday / rsiyea )
      omega  = 2. * rpi / rsiday 
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          day                                rday   = ', rday,   ' s'
      IF(lwp) WRITE(numout,*) '          sideral year                       rsiyea = ', rsiyea, ' s'
      IF(lwp) WRITE(numout,*) '          sideral day                        rsiday = ', rsiday, ' s'
      IF(lwp) WRITE(numout,*) '          omega                              omega  = ', omega,  ' s-1'

      rjjss = rjjhh * rhhmm * rmmss
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          nb of months per year               raamo = ', raamo, ' months'
      IF(lwp) WRITE(numout,*) '          nb of hours per day                 rjjhh = ', rjjhh, ' hours'
      IF(lwp) WRITE(numout,*) '          nb of minutes per hour              rhhmm = ', rhhmm, ' mn'
      IF(lwp) WRITE(numout,*) '          nb of seconds per minute            rmmss = ', rmmss, ' s'
      IF(lwp) WRITE(numout,*) '          nb of seconds per day               rjjss = ', rjjss, ' s'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          earth radius                         ra   = ', ra, ' m'
      IF(lwp) WRITE(numout,*) '          gravity                              grav = ', grav , ' m/s^2'

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          triple point of temperature      rtt      = ', rtt     , ' K'
      IF(lwp) WRITE(numout,*) '          freezing point of water          rt0      = ', rt0     , ' K'
      IF(lwp) WRITE(numout,*) '          melting point of snow            rt0_snow = ', rt0_snow, ' K'
      IF(lwp) WRITE(numout,*) '          melting point of ice             rt0_ice  = ', rt0_ice , ' K'

      ro0cpr = 1. / ( rau0 * rcp )
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          volumic mass of pure water         rauw   = ', rauw, ' kg/m^3'
      IF(lwp) WRITE(numout,*) '          volumic mass of reference          rau0   = ', rau0, ' kg/m^3'
      IF(lwp) WRITE(numout,*) '          ocean specific heat                rcp    = ', rcp
      IF(lwp) WRITE(numout,*) '                       1. / ( rau0 * rcp ) = ro0cpr = ', ro0cpr

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          thermal conductivity of the snow          = ', rcdsn   , ' J/s/m/K'
         WRITE(numout,*) '          thermal conductivity of the ice           = ', rcdic   , ' J/s/m/K'
         WRITE(numout,*) '          density times specific heat for snow      = ', rcpsn   , ' J/m^3/K' 
         WRITE(numout,*) '          density times specific heat for ice       = ', rcpic   , ' J/m^3/K'
         WRITE(numout,*) '          volumetric latent heat fusion of sea ice  = ', xlic    , ' J/m' 
         WRITE(numout,*) '          volumetric latent heat fusion of snow     = ', xlsn    , ' J/m' 
         WRITE(numout,*) '          latent heat of sublimation of snow        = ', xsn     , ' J/kg' 
         WRITE(numout,*) '          density of sea ice                        = ', rhoic   , ' kg/m^3'
         WRITE(numout,*) '          density of snow                           = ', rhosn   , ' kg/m^3'
         WRITE(numout,*) '          emissivity of snow or ice                 = ', emic  
         WRITE(numout,*) '          salinity of ice                           = ', sice    , ' psu'
         WRITE(numout,*) '          salinity of sea                           = ', soce    , ' psu'
         WRITE(numout,*) '          latent heat of evaporation (water)        = ', cevap   , ' J/m^3' 
         WRITE(numout,*) '          correction factor for solar radiation     = ', srgamma 
         WRITE(numout,*) '          von Karman constant                       = ', vkarmn 
         WRITE(numout,*) '          Stefan-Boltzmann constant                 = ', stefan  , ' J/s/m^2/K^4'

         WRITE(numout,*)
         WRITE(numout,*) '          conversion: degre ==> radian          rad = ', rad

         WRITE(numout,*)
         WRITE(numout,*) '          smallest real computer value       rsmall = ', rsmall
      ENDIF

   END SUBROUTINE phy_cst

   !!======================================================================
END MODULE phycst
