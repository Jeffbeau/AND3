


PGF90 (Version     10.5)          09/23/2010  09:49:00      page 1

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

Filename: phycst.F90

(    1) MODULE phycst
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  phycst  ***
(    4)    !!     Definition of of both ocean and ice parameters used in the code
(    5)    !!=====================================================================
(    6)    !! * Modules used
(    7)    USE par_oce          ! ocean parameters
(    8)    USE in_out_manager   ! I/O manager
(    9) 
(   10)    IMPLICIT NONE
(   11)    PRIVATE
(   12) 
(   13)    !! * Routine accessibility
(   14)    PUBLIC phy_cst          ! routine called by inipar.F90
(   15) 
(   16)    !! * Shared module variables
(   17)    INTEGER, PUBLIC, DIMENSION(12) ::   &  !:
(   18)       nbiss = (/ 31, 29, 31, 30, 31, 30,      &  !: number of days per month
(   19)          &       31, 31, 30, 31, 30, 31 /) ,  &  !  (leap-year)
(   20)       nobis = (/ 31, 28, 31, 30, 31, 30,      &  !: number of days per month
(   21)          &       31, 31, 30, 31, 30, 31 /)       !  (365 days a year)
(   22)    
(   23)    REAL(wp), PUBLIC ::                        &  !:
(   24)       rpi = 3.141592653589793_wp           ,  &  !: pi
(   25)       rad = 3.141592653589793_wp / 180._wp ,  &  !: conversion from degre into radian
(   26)       rsmall = 0.5 * EPSILON( 1. )               !: smallest real computer value
(   27)    
(   28)    REAL(wp), PUBLIC ::          & !:
(   29)       rday = 24.*60.*60.  ,     & !: day (s)
(   30)       rsiyea              ,     & !: sideral year (s)
(   31)       rsiday              ,     & !: sideral day (s)
(   32)       raajj = 365._wp     ,     & !: number of days in one year
(   33)       raamo =  12._wp     ,     & !: number of months in one year
(   34)       rjjhh =  24._wp     ,     & !: number of hours in one day
(   35)       rhhmm =  60._wp     ,     & !: number of minutes in one hour
(   36)       rmmss =  60._wp     ,     & !: number of seconds in one minute
(   37)       raass               ,     & !: number of seconds in one year
(   38)       rmoss               ,     & !: number of seconds in one month
(   39)       rjjss               ,     & !: number of seconds in one day
(   40) !!!   omega = 7.292115083046061e-5_wp ,  &  !: change the last digit!
(   41)       omega               ,    &  !: earth rotation parameter
(   42)       ra    = 6371229._wp ,    &  !: earth radius (meter)
(   43)       grav  = 9.80665_wp          !: gravity (m/s2)
(   44)    
(   45)    REAL(wp), PUBLIC ::         &  !:
(   46)       rtt      = 273.16_wp  ,  &  !: triple point of temperature (Kelvin)






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 2

(   47)       rt0      = 273.15_wp  ,  &  !: freezing point of water (Kelvin)
(   48)       rt0_snow = 273.15_wp  ,  &  !: melting point of snow  (Kelvin)
(   49)       rt0_ice  = 273.05_wp  ,  &  !: melting point of ice   (Kelvin)
(   50)       rau0     = 1020._wp   ,  &  !: volumic mass of reference (kg/m3)
(   51)       rauw     = 1000._wp   ,  &  !: density of pure water (kg/m3)
(   52)       rcp      =    4.e+3_wp,  &  !: ocean specific heat
(   53)       ro0cpr                      !: = 1. / ( rau0 * rcp )
(   54) 
(   55)    REAL(wp), PUBLIC ::            &  !:
(   56)       rcdsn   =   0.22_wp     ,   &  !: conductivity of the snow
(   57)       rcdic   =   2.034396_wp ,   &  !: conductivity of the ice
(   58)       rcpsn   =   6.9069e+5_wp,   &  !: density times specific heat for snow
(   59)       rcpic   =   1.8837e+6_wp,   &  !: volumetric latent heat fusion of sea ice
(   60)       xlsn    = 110.121e+6_wp ,   &  !: volumetric latent heat fusion of snow
(   61)       xlic    = 300.33e+6_wp  ,   &  !: volumetric latent heat fusion of ice
(   62)       xsn     =   2.8e+6      ,   &  !: latent heat of sublimation of snow
(   63)       rhoic   = 900._wp       ,   &  !: density of sea ice (kg/m3)
(   64)       rhosn   = 330._wp       ,   &  !: density of snow (kg/m3)
(   65)       emic    =   0.97_wp     ,   &  !: emissivity of snow or ice
(   66)       sice    =   6.0_wp      ,   &  !: salinity of ice (psu)
(   67)       soce    =  34.7_wp      ,   &  !: salinity of sea (psu)
(   68)       cevap   =   2.5e+6_wp   ,   &  !: latent heat of evaporation (water)
(   69)       srgamma =   0.9_wp      ,   &  !: correction factor for solar radiation (Oberhuber, 1974)
(   70)       vkarmn  =   0.4_wp      ,   &  !: von Karman constant
(   71)       stefan  =   5.67e-8_wp         !: Stefan-Boltzmann constant 
(   72)       !!----------------------------------------------------------------------
(   73)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   74)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/phycst.F90,v 1.4 2005/03/27 18:34:48 opalod Exp $ 
(   75)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   76)       !!----------------------------------------------------------------------
(   77)    
(   78) CONTAINS
(   79)    
(   80)    SUBROUTINE phy_cst
(   81)       !!----------------------------------------------------------------------
(   82)       !!                       ***  ROUTINE phy_cst  ***
(   83)       !!
(   84)       !! ** Purpose :   Print model parameters and set and print the constants
(   85)       !!
(   86)       !! ** Method  :   no
(   87)       !!
(   88)       !! History :
(   89)       !!        !  90-10  (C. Levy - G. Madec)  Original code
(   90)       !!        !  91-11  (G. Madec)
(   91)       !!        !  91-12  (M. Imbard)
(   92)       !!   8.5  !  02-08  (G. Madec, C. Ethe)  F90, add ice constants 
(   93)       !!----------------------------------------------------------------------
(   94)       !! * Local variables
(   95)       CHARACTER (len=64) ::   cform = "(A9, 3(A13, I7) )" 
(   96)       !!----------------------------------------------------------------------
(   97) 
(   98)       IF(lwp) WRITE(numout,*)
(   99)       IF(lwp) WRITE(numout,*) ' phy_cst : initialization of ocean parameters and constants'
(  100)       IF(lwp) WRITE(numout,*) ' ~~~~~~~'
(  101) 
(  102)       ! Ocean Parameters
(  103)       ! ----------------
(  104)       IF(lwp) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 3

(  105)          WRITE(numout,*) '       parameter file'
(  106)          WRITE(numout,*)
(  107)          WRITE(numout,*) '          dimension of model'
(  108)          WRITE(numout,*) '              Local domain      Global domain       Data domain '
(  109)          WRITE(numout,cform) '         ','   jpi     : ', jpi, '   jpiglo  : ', jpiglo, '   jpidta  : ', jpidta
(  110)          WRITE(numout,cform) '         ','   jpj     : ', jpj, '   jpjglo  : ', jpjglo, '   jpjdta  : ', jpjdta
(  111)          WRITE(numout,cform) '         ','   jpk     : ', jpk, '   jpk     : ', jpk   , '   jpkdta  : ', jpkdta
(  112)          WRITE(numout,*)      '        ','   jpij    : ', jpij
(  113)          WRITE(numout,*)
(  114)          WRITE(numout,*) '          mpp local domain info (mpp)'
(  115)          WRITE(numout,*) '             jpni    : ', jpni, '   jpreci  : ', jpreci
(  116)          WRITE(numout,*) '             jpnj    : ', jpnj, '   jprecj  : ', jprecj
(  117)          WRITE(numout,*) '             jpnij   : ', jpnij
(  118) 
(  119)          WRITE(numout,*)
(  120)          WRITE(numout,*) '          lateral domain boundary condition type : jperio  = ', jperio
(  121)          WRITE(numout,*) '          domain island (use in rigid-lid case)  : jpisl   = ', jpisl 
(  122)          WRITE(numout,*) '                                                   jpnisl  = ', jpnisl
(  123)       ENDIF
(  124) 
(  125)       ! Define constants
(  126)       ! ----------------
(  127)       IF(lwp) WRITE(numout,*)
(  128)       IF(lwp) WRITE(numout,*) '       constants'
(  129) 
(  130)       IF(lwp) WRITE(numout,*)
(  131)       IF(lwp) WRITE(numout,*) '          mathematical constant                 rpi = ', rpi
(  132) 
(  133)       rsiyea = 365.25 * rday * 2. * rpi / 6.283076
(  134)       rsiday = rday / ( 1. + rday / rsiyea )
(  135)       omega  = 2. * rpi / rsiday 
(  136)       IF(lwp) WRITE(numout,*)
(  137)       IF(lwp) WRITE(numout,*) '          day                                rday   = ', rday,   ' s'
(  138)       IF(lwp) WRITE(numout,*) '          sideral year                       rsiyea = ', rsiyea, ' s'
(  139)       IF(lwp) WRITE(numout,*) '          sideral day                        rsiday = ', rsiday, ' s'
(  140)       IF(lwp) WRITE(numout,*) '          omega                              omega  = ', omega,  ' s-1'
(  141) 
(  142)       rjjss = rjjhh * rhhmm * rmmss
(  143)       IF(lwp) WRITE(numout,*)
(  144)       IF(lwp) WRITE(numout,*) '          nb of months per year               raamo = ', raamo, ' months'
(  145)       IF(lwp) WRITE(numout,*) '          nb of hours per day                 rjjhh = ', rjjhh, ' hours'
(  146)       IF(lwp) WRITE(numout,*) '          nb of minutes per hour              rhhmm = ', rhhmm, ' mn'
(  147)       IF(lwp) WRITE(numout,*) '          nb of seconds per minute            rmmss = ', rmmss, ' s'
(  148)       IF(lwp) WRITE(numout,*) '          nb of seconds per day               rjjss = ', rjjss, ' s'
(  149) 
(  150)       IF(lwp) WRITE(numout,*)
(  151)       IF(lwp) WRITE(numout,*) '          earth radius                         ra   = ', ra, ' m'
(  152)       IF(lwp) WRITE(numout,*) '          gravity                              grav = ', grav , ' m/s^2'
(  153) 
(  154)       IF(lwp) WRITE(numout,*)
(  155)       IF(lwp) WRITE(numout,*) '          triple point of temperature      rtt      = ', rtt     , ' K'
(  156)       IF(lwp) WRITE(numout,*) '          freezing point of water          rt0      = ', rt0     , ' K'
(  157)       IF(lwp) WRITE(numout,*) '          melting point of snow            rt0_snow = ', rt0_snow, ' K'
(  158)       IF(lwp) WRITE(numout,*) '          melting point of ice             rt0_ice  = ', rt0_ice , ' K'
(  159) 
(  160)       ro0cpr = 1. / ( rau0 * rcp )
(  161)       IF(lwp) WRITE(numout,*)
(  162)       IF(lwp) WRITE(numout,*) '          volumic mass of pure water         rauw   = ', rauw, ' kg/m^3'






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 4

(  163)       IF(lwp) WRITE(numout,*) '          volumic mass of reference          rau0   = ', rau0, ' kg/m^3'
(  164)       IF(lwp) WRITE(numout,*) '          ocean specific heat                rcp    = ', rcp
(  165)       IF(lwp) WRITE(numout,*) '                       1. / ( rau0 * rcp ) = ro0cpr = ', ro0cpr
(  166) 
(  167)       IF(lwp) THEN
(  168)          WRITE(numout,*)
(  169)          WRITE(numout,*) '          thermal conductivity of the snow          = ', rcdsn   , ' J/s/m/K'
(  170)          WRITE(numout,*) '          thermal conductivity of the ice           = ', rcdic   , ' J/s/m/K'
(  171)          WRITE(numout,*) '          density times specific heat for snow      = ', rcpsn   , ' J/m^3/K' 
(  172)          WRITE(numout,*) '          density times specific heat for ice       = ', rcpic   , ' J/m^3/K'
(  173)          WRITE(numout,*) '          volumetric latent heat fusion of sea ice  = ', xlic    , ' J/m' 
(  174)          WRITE(numout,*) '          volumetric latent heat fusion of snow     = ', xlsn    , ' J/m' 
(  175)          WRITE(numout,*) '          latent heat of sublimation of snow        = ', xsn     , ' J/kg' 
(  176)          WRITE(numout,*) '          density of sea ice                        = ', rhoic   , ' kg/m^3'
(  177)          WRITE(numout,*) '          density of snow                           = ', rhosn   , ' kg/m^3'
(  178)          WRITE(numout,*) '          emissivity of snow or ice                 = ', emic  
(  179)          WRITE(numout,*) '          salinity of ice                           = ', sice    , ' psu'
(  180)          WRITE(numout,*) '          salinity of sea                           = ', soce    , ' psu'
(  181)          WRITE(numout,*) '          latent heat of evaporation (water)        = ', cevap   , ' J/m^3' 
(  182)          WRITE(numout,*) '          correction factor for solar radiation     = ', srgamma 
(  183)          WRITE(numout,*) '          von Karman constant                       = ', vkarmn 
(  184)          WRITE(numout,*) '          Stefan-Boltzmann constant                 = ', stefan  , ' J/s/m^2/K^4'
(  185) 
(  186)          WRITE(numout,*)
(  187)          WRITE(numout,*) '          conversion: degre ==> radian          rad = ', rad
(  188) 
(  189)          WRITE(numout,*)
(  190)          WRITE(numout,*) '          smallest real computer value       rsmall = ', rsmall
(  191)       ENDIF
(  192) 
(  193)    END SUBROUTINE phy_cst
(  194) 
(  195)    !!======================================================================
(  196) END MODULE phycst






























PGF90 (Version     10.5)          09/23/2010  09:49:00      page 5

(    1) # 1 "phycst.F90"
(    1) MODULE phycst
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  phycst  ***
(    4)    !!     Definition of of both ocean and ice parameters used in the code
(    5)    !!=====================================================================
(    6)    !! * Modules used
(    7)    USE par_oce          ! ocean parameters
(    8)    USE in_out_manager   ! I/O manager
(    9) 
(   10)    IMPLICIT NONE
(   11)    PRIVATE
(   12) 
(   13)    !! * Routine accessibility
(   14)    PUBLIC phy_cst          ! routine called by inipar.F90
(   15) 
(   16)    !! * Shared module variables
(   17)    INTEGER, PUBLIC, DIMENSION(12) ::   &  !:
(   18)       nbiss = (/ 31, 29, 31, 30, 31, 30,      &  !: number of days per month
(   19)          &       31, 31, 30, 31, 30, 31 /) ,  &  !  (leap-year)
(   20)       nobis = (/ 31, 28, 31, 30, 31, 30,      &  !: number of days per month
(   21)          &       31, 31, 30, 31, 30, 31 /)       !  (365 days a year)
(   22)    
(   23)    REAL(wp), PUBLIC ::                        &  !:
(   24)       rpi = 3.141592653589793_wp           ,  &  !: pi
(   25)       rad = 3.141592653589793_wp / 180._wp ,  &  !: conversion from degre into radian
(   26)       rsmall = 0.5 * EPSILON( 1. )               !: smallest real computer value
(   27)    
(   28)    REAL(wp), PUBLIC ::          & !:
(   29)       rday = 24.*60.*60.  ,     & !: day (s)
(   30)       rsiyea              ,     & !: sideral year (s)
(   31)       rsiday              ,     & !: sideral day (s)
(   32)       raajj = 365._wp     ,     & !: number of days in one year
(   33)       raamo =  12._wp     ,     & !: number of months in one year
(   34)       rjjhh =  24._wp     ,     & !: number of hours in one day
(   35)       rhhmm =  60._wp     ,     & !: number of minutes in one hour
(   36)       rmmss =  60._wp     ,     & !: number of seconds in one minute
(   37)       raass               ,     & !: number of seconds in one year
(   38)       rmoss               ,     & !: number of seconds in one month
(   39)       rjjss               ,     & !: number of seconds in one day
(   40) !!!   omega = 7.292115083046061e-5_wp ,  &  !: change the last digit!
(   41)       omega               ,    &  !: earth rotation parameter
(   42)       ra    = 6371229._wp ,    &  !: earth radius (meter)
(   43)       grav  = 9.80665_wp          !: gravity (m/s2)
(   44)    
(   45)    REAL(wp), PUBLIC ::         &  !:
(   46)       rtt      = 273.16_wp  ,  &  !: triple point of temperature (Kelvin)
(   47)       rt0      = 273.15_wp  ,  &  !: freezing point of water (Kelvin)
(   48)       rt0_snow = 273.15_wp  ,  &  !: melting point of snow  (Kelvin)
(   49)       rt0_ice  = 273.05_wp  ,  &  !: melting point of ice   (Kelvin)
(   50)       rau0     = 1020._wp   ,  &  !: volumic mass of reference (kg/m3)
(   51)       rauw     = 1000._wp   ,  &  !: density of pure water (kg/m3)
(   52)       rcp      =    4.e+3_wp,  &  !: ocean specific heat
(   53)       ro0cpr                      !: = 1. / ( rau0 * rcp )
(   54) 
(   55)    REAL(wp), PUBLIC ::            &  !:
(   56)       rcdsn   =   0.22_wp     ,   &  !: conductivity of the snow
(   57)       rcdic   =   2.034396_wp ,   &  !: conductivity of the ice






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 6

(   58)       rcpsn   =   6.9069e+5_wp,   &  !: density times specific heat for snow
(   59)       rcpic   =   1.8837e+6_wp,   &  !: volumetric latent heat fusion of sea ice
(   60)       xlsn    = 110.121e+6_wp ,   &  !: volumetric latent heat fusion of snow
(   61)       xlic    = 300.33e+6_wp  ,   &  !: volumetric latent heat fusion of ice
(   62)       xsn     =   2.8e+6      ,   &  !: latent heat of sublimation of snow
(   63)       rhoic   = 900._wp       ,   &  !: density of sea ice (kg/m3)
(   64)       rhosn   = 330._wp       ,   &  !: density of snow (kg/m3)
(   65)       emic    =   0.97_wp     ,   &  !: emissivity of snow or ice
(   66)       sice    =   6.0_wp      ,   &  !: salinity of ice (psu)
(   67)       soce    =  34.7_wp      ,   &  !: salinity of sea (psu)
(   68)       cevap   =   2.5e+6_wp   ,   &  !: latent heat of evaporation (water)
(   69)       srgamma =   0.9_wp      ,   &  !: correction factor for solar radiation (Oberhuber, 1974)
(   70)       vkarmn  =   0.4_wp      ,   &  !: von Karman constant
(   71)       stefan  =   5.67e-8_wp         !: Stefan-Boltzmann constant 
(   72)       !!----------------------------------------------------------------------
(   73)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   74)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/phycst.F90,v 1.4 2005/03/27 18:34:48 opalod Exp $ 
(   75)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   76)       !!----------------------------------------------------------------------
(   77)    
(   78) CONTAINS
(   79)    
(   80)    SUBROUTINE phy_cst
(   81)       !!----------------------------------------------------------------------
(   82)       !!                       ***  ROUTINE phy_cst  ***
(   83)       !!
(   84)       !! ** Purpose :   Print model parameters and set and print the constants
(   85)       !!
(   86)       !! ** Method  :   no
(   87)       !!
(   88)       !! History :
(   89)       !!        !  90-10  (C. Levy - G. Madec)  Original code
(   90)       !!        !  91-11  (G. Madec)
(   91)       !!        !  91-12  (M. Imbard)
(   92)       !!   8.5  !  02-08  (G. Madec, C. Ethe)  F90, add ice constants 
(   93)       !!----------------------------------------------------------------------
(   94)       !! * Local variables
(   95)       CHARACTER (len=64) ::   cform = "(A9, 3(A13, I7) )" 
(   96)       !!----------------------------------------------------------------------
(   97) 
(   98)       IF(lwp) WRITE(numout,*)
(   99)       IF(lwp) WRITE(numout,*) ' phy_cst : initialization of ocean parameters and constants'
(  100)       IF(lwp) WRITE(numout,*) ' ~~~~~~~'
(  101) 
(  102)       ! Ocean Parameters
(  103)       ! ----------------
(  104)       IF(lwp) THEN
(  105)          WRITE(numout,*) '       parameter file'
(  106)          WRITE(numout,*)
(  107)          WRITE(numout,*) '          dimension of model'
(  108)          WRITE(numout,*) '              Local domain      Global domain       Data domain '
(  109)          WRITE(numout,cform) '         ','   jpi     : ', jpi, '   jpiglo  : ', jpiglo, '   jpidta  : ', jpidta
(  110)          WRITE(numout,cform) '         ','   jpj     : ', jpj, '   jpjglo  : ', jpjglo, '   jpjdta  : ', jpjdta
(  111)          WRITE(numout,cform) '         ','   jpk     : ', jpk, '   jpk     : ', jpk   , '   jpkdta  : ', jpkdta
(  112)          WRITE(numout,*)      '        ','   jpij    : ', jpij
(  113)          WRITE(numout,*)
(  114)          WRITE(numout,*) '          mpp local domain info (mpp)'
(  115)          WRITE(numout,*) '             jpni    : ', jpni, '   jpreci  : ', jpreci






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 7

(  116)          WRITE(numout,*) '             jpnj    : ', jpnj, '   jprecj  : ', jprecj
(  117)          WRITE(numout,*) '             jpnij   : ', jpnij
(  118) 
(  119)          WRITE(numout,*)
(  120)          WRITE(numout,*) '          lateral domain boundary condition type : jperio  = ', jperio
(  121)          WRITE(numout,*) '          domain island (use in rigid-lid case)  : jpisl   = ', jpisl 
(  122)          WRITE(numout,*) '                                                   jpnisl  = ', jpnisl
(  123)       ENDIF
(  124) 
(  125)       ! Define constants
(  126)       ! ----------------
(  127)       IF(lwp) WRITE(numout,*)
(  128)       IF(lwp) WRITE(numout,*) '       constants'
(  129) 
(  130)       IF(lwp) WRITE(numout,*)
(  131)       IF(lwp) WRITE(numout,*) '          mathematical constant                 rpi = ', rpi
(  132) 
(  133)       rsiyea = 365.25 * rday * 2. * rpi / 6.283076
(  134)       rsiday = rday / ( 1. + rday / rsiyea )
(  135)       omega  = 2. * rpi / rsiday 
(  136)       IF(lwp) WRITE(numout,*)
(  137)       IF(lwp) WRITE(numout,*) '          day                                rday   = ', rday,   ' s'
(  138)       IF(lwp) WRITE(numout,*) '          sideral year                       rsiyea = ', rsiyea, ' s'
(  139)       IF(lwp) WRITE(numout,*) '          sideral day                        rsiday = ', rsiday, ' s'
(  140)       IF(lwp) WRITE(numout,*) '          omega                              omega  = ', omega,  ' s-1'
(  141) 
(  142)       rjjss = rjjhh * rhhmm * rmmss
(  143)       IF(lwp) WRITE(numout,*)
(  144)       IF(lwp) WRITE(numout,*) '          nb of months per year               raamo = ', raamo, ' months'
(  145)       IF(lwp) WRITE(numout,*) '          nb of hours per day                 rjjhh = ', rjjhh, ' hours'
(  146)       IF(lwp) WRITE(numout,*) '          nb of minutes per hour              rhhmm = ', rhhmm, ' mn'
(  147)       IF(lwp) WRITE(numout,*) '          nb of seconds per minute            rmmss = ', rmmss, ' s'
(  148)       IF(lwp) WRITE(numout,*) '          nb of seconds per day               rjjss = ', rjjss, ' s'
(  149) 
(  150)       IF(lwp) WRITE(numout,*)
(  151)       IF(lwp) WRITE(numout,*) '          earth radius                         ra   = ', ra, ' m'
(  152)       IF(lwp) WRITE(numout,*) '          gravity                              grav = ', grav , ' m/s^2'
(  153) 
(  154)       IF(lwp) WRITE(numout,*)
(  155)       IF(lwp) WRITE(numout,*) '          triple point of temperature      rtt      = ', rtt     , ' K'
(  156)       IF(lwp) WRITE(numout,*) '          freezing point of water          rt0      = ', rt0     , ' K'
(  157)       IF(lwp) WRITE(numout,*) '          melting point of snow            rt0_snow = ', rt0_snow, ' K'
(  158)       IF(lwp) WRITE(numout,*) '          melting point of ice             rt0_ice  = ', rt0_ice , ' K'
(  159) 
(  160)       ro0cpr = 1. / ( rau0 * rcp )
(  161)       IF(lwp) WRITE(numout,*)
(  162)       IF(lwp) WRITE(numout,*) '          volumic mass of pure water         rauw   = ', rauw, ' kg/m^3'
(  163)       IF(lwp) WRITE(numout,*) '          volumic mass of reference          rau0   = ', rau0, ' kg/m^3'
(  164)       IF(lwp) WRITE(numout,*) '          ocean specific heat                rcp    = ', rcp
(  165)       IF(lwp) WRITE(numout,*) '                       1. / ( rau0 * rcp ) = ro0cpr = ', ro0cpr
(  166) 
(  167)       IF(lwp) THEN
(  168)          WRITE(numout,*)
(  169)          WRITE(numout,*) '          thermal conductivity of the snow          = ', rcdsn   , ' J/s/m/K'
(  170)          WRITE(numout,*) '          thermal conductivity of the ice           = ', rcdic   , ' J/s/m/K'
(  171)          WRITE(numout,*) '          density times specific heat for snow      = ', rcpsn   , ' J/m^3/K' 
(  172)          WRITE(numout,*) '          density times specific heat for ice       = ', rcpic   , ' J/m^3/K'
(  173)          WRITE(numout,*) '          volumetric latent heat fusion of sea ice  = ', xlic    , ' J/m' 






PGF90 (Version     10.5)          09/23/2010  09:49:00      page 8

(  174)          WRITE(numout,*) '          volumetric latent heat fusion of snow     = ', xlsn    , ' J/m' 
(  175)          WRITE(numout,*) '          latent heat of sublimation of snow        = ', xsn     , ' J/kg' 
(  176)          WRITE(numout,*) '          density of sea ice                        = ', rhoic   , ' kg/m^3'
(  177)          WRITE(numout,*) '          density of snow                           = ', rhosn   , ' kg/m^3'
(  178)          WRITE(numout,*) '          emissivity of snow or ice                 = ', emic  
(  179)          WRITE(numout,*) '          salinity of ice                           = ', sice    , ' psu'
(  180)          WRITE(numout,*) '          salinity of sea                           = ', soce    , ' psu'
(  181)          WRITE(numout,*) '          latent heat of evaporation (water)        = ', cevap   , ' J/m^3' 
(  182)          WRITE(numout,*) '          correction factor for solar radiation     = ', srgamma 
(  183)          WRITE(numout,*) '          von Karman constant                       = ', vkarmn 
(  184)          WRITE(numout,*) '          Stefan-Boltzmann constant                 = ', stefan  , ' J/s/m^2/K^4'
(  185) 
(  186)          WRITE(numout,*)
(  187)          WRITE(numout,*) '          conversion: degre ==> radian          rad = ', rad
(  188) 
(  189)          WRITE(numout,*)
(  190)          WRITE(numout,*) '          smallest real computer value       rsmall = ', rsmall
(  191)       ENDIF
(  192) 
(  193)    END SUBROUTINE phy_cst
(  194) 
(  195)    !!======================================================================
(  196) END MODULE phycst
