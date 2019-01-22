


PGF90 (Version     10.5)          09/23/2010  09:49:01      page 1

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

Filename: flx_oce.F90

(    1) MODULE flx_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  flx_oce  ***
(    4)    !!        parameter and  variables defined in memory in forced mode
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flx_oce.F90,v 1.5 2005/03/27 18:35:13 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_oce          ! ocean parameters
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18)    
(   19) 
(   20)    !!----------------------------------------------------------------------
(   21)    !! fluxes common variables
(   22)    !!----------------------------------------------------------------------
(   23) #if defined key_flx_forced_daily
(   24)    !!----------------------------------------------------------------------
(   25)    !! 'key_flx_forced_daily'
(   26)    !!----------------------------------------------------------------------
(   27)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   28)       p_qt ,        &   !: total heat flux ( solar + non solar)
(   29)       p_qsr,        &   !: solar heat flux
(   30)       p_emp             !: evaporation minus precipitation            
(   31) #elif defined key_flx_forced_monthly   
(   32) !byoung
(   33)     !----------------------------------------------------------------
(   34)     !'key_flx_forced_monthly'
(   35)     !-------------------------------------------------------------
(   36)     REAL(wp), PUBLIC,DIMENSION(jpi,jpj):: &
(   37)     p_bqt,   &  !:total heat flux
(   38)     p_bqsr,  &  !:solar heat flux
(   39)     p_bemp      !:evaporation minus precipitation
(   40) 
(   41) #elif defined key_ice_lim || defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(   42)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
(   43)       qsr_ice  ,      &  !: solar flux over ice
(   44)       qsr_oce  ,      &  !: solar flux over ocean
(   45)       qnsr_oce ,      &  !: total non solar heat flux (Longwave downward radiation) over ocean 
(   46)       qnsr_ice ,      &  !: total non solar heat flux (Longwave downward radiation) over ice






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 2

(   47)       tprecip  ,      &  !: total precipitation ( or liquid precip minus evaporation in coupled mode)
(   48)       sprecip  ,      &  !: solid (snow) precipitation
(   49)       dqns_ice ,      &  !: total non solar sensibility over ice (LW+SEN+LA)
(   50)       tn_ice   ,      &  !: ice surface temperature
(   51)       evap     ,      &  !: evaporation over ocean
(   52)       fr1_i0   ,      &  !: 1st part of the fraction of sol. rad.  which penetrate inside the ice cover
(   53)       fr2_i0   ,      &  !: 2nd part of the fraction of sol. rad.  which penetrate inside the ice cover 
(   54) #if ! defined key_coupled
(   55)       qla_ice  ,      &  !: latent flux over ice  
(   56)       dqla_ice           !: latent sensibility over ice
(   57) #else
(   58)       rrunoff  ,      &  !: runoff
(   59)       calving  ,      &  !: calving
(   60)       alb_ice            !: albedo of ice      
(   61) #endif
(   62) 
(   63) #else
(   64)    !!----------------------------------------------------------------------
(   65)    !!   Default option                                         Empty module
(   66)    !!----------------------------------------------------------------------
(   67)    
(   68) #endif
(   69) 
(   70)    !!----------------------------------------------------------------------
(   71) END MODULE flx_oce







































PGF90 (Version     10.5)          09/23/2010  09:49:01      page 3

(    1) # 1 "flx_oce.F90"
(    1) MODULE flx_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  flx_oce  ***
(    4)    !!        parameter and  variables defined in memory in forced mode
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flx_oce.F90,v 1.5 2005/03/27 18:35:13 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_oce          ! ocean parameters
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18)    
(   19) 
(   20)    !!----------------------------------------------------------------------
(   21)    !! fluxes common variables
(   22)    !!----------------------------------------------------------------------
(   23) # 42
(   42)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
(   43)       qsr_ice  ,      &  !: solar flux over ice
(   44)       qsr_oce  ,      &  !: solar flux over ocean
(   45)       qnsr_oce ,      &  !: total non solar heat flux (Longwave downward radiation) over ocean 
(   46)       qnsr_ice ,      &  !: total non solar heat flux (Longwave downward radiation) over ice
(   47)       tprecip  ,      &  !: total precipitation ( or liquid precip minus evaporation in coupled mode)
(   48)       sprecip  ,      &  !: solid (snow) precipitation
(   49)       dqns_ice ,      &  !: total non solar sensibility over ice (LW+SEN+LA)
(   50)       tn_ice   ,      &  !: ice surface temperature
(   51)       evap     ,      &  !: evaporation over ocean
(   52)       fr1_i0   ,      &  !: 1st part of the fraction of sol. rad.  which penetrate inside the ice cover
(   53)       fr2_i0   ,      &  !: 2nd part of the fraction of sol. rad.  which penetrate inside the ice cover 
(   55)       qla_ice  ,      &  !: latent flux over ice  
(   56)       dqla_ice           !: latent sensibility over ice
(   57) 
(   58) 
(   59) # 70
(   70)    !!----------------------------------------------------------------------
(   71) END MODULE flx_oce
