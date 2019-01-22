


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

Filename: ice_oce.F90

(    1) MODULE ice_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  ice_oce  ***
(    4)    !! Ocean - ice  :  ice variables defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ice_oce.F90,v 1.3 2005/03/27 18:34:46 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) #if defined key_ice_lim
(   14)    !!----------------------------------------------------------------------
(   15)    !!   'key_ice_lim'   :                                     LIM ice model
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce         ! ocean parameters
(   19)    USE blk_oce         ! bulk parameters
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23)  
(   24)    !! Shared module variables
(   25)    LOGICAL, PUBLIC, PARAMETER ::   lk_ice_lim = .TRUE.    !: LIM ice model
(   26) 
(   27)    !!----------------------------------------------------------------------
(   28)    !! ice-ocean common variables
(   29)    !!----------------------------------------------------------------------
(   30) # if defined key_coupled
(   31)    REAL(wp), PUBLIC, DIMENSION(jpiglo,jpjglo) ::   &  !: cumulated fields
(   32)       fqsr_oce ,      &   !: Net short wave heat flux on free ocean 
(   33)       fqsr_ice ,      &   !: Net short wave het flux on sea ice 
(   34)       fqnsr_oce,      &   !: Net longwave heat flux on free ocean
(   35)       fqnsr_ice,      &   !: Net longwave heat flux on sea ice
(   36)       fdqns_ice,      &   !: Derivative of non solar heat flux on sea ice
(   37)       ftprecip ,      &   !: Water flux (liquid precipitation - evaporation) 
(   38)       fsprecip ,      &   !: Solid (snow) precipitation
(   39)       frunoff  ,      &   !: runoff
(   40)       fcalving            !: Iceberg calving 
(   41) # endif
(   42) 
(   43)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !: field exchanges with ice model to ocean
(   44)       sst_io, sss_io , &  !: sea surface temperature (C) and salinity (PSU)
(   45)       u_io  , v_io   , &  !: velocity at ice surface (m/s)
(   46)       fsolar, fnsolar, &  !: solar and non-solar heat fluxes (W/m2)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 2

(   47)       fsalt , fmass  , &  !: salt and freshwater fluxes
(   48)       ftaux , ftauy  , &  !: wind stresses
(   49)       gtaux , gtauy       !: wind stresses
(   50)    
(   51)    REAL(wp), PUBLIC ::   &  !:
(   52)       rdt_ice,           &  !: ice time step
(   53)       dtsd2                 !: ice time step divide by 2
(   54) 
(   55) #else
(   56)    !!----------------------------------------------------------------------
(   57)    !!   Default option                                 NO LIM sea-ice model
(   58)    !!----------------------------------------------------------------------
(   59)    LOGICAL, PUBLIC, PARAMETER ::   lk_ice_lim = .FALSE.  !: No LIM ice model
(   60) #endif
(   61) 
(   62)    INTEGER, PUBLIC ::   &  !: namdom : space/time domain (namlist)
(   63)       nfice =  5           !: coupling frequency OPA ICELLN  nfice 
(   64) 
(   65)    !!----------------------------------------------------------------------
(   66) END MODULE ice_oce












































PGF90 (Version     10.5)          09/23/2010  09:49:11      page 3

(    1) # 1 "ice_oce.F90"
(    1) MODULE ice_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  ice_oce  ***
(    4)    !! Ocean - ice  :  ice variables defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ice_oce.F90,v 1.3 2005/03/27 18:34:46 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) # 14
(   14)    !!----------------------------------------------------------------------
(   15)    !!   'key_ice_lim'   :                                     LIM ice model
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce         ! ocean parameters
(   19)    USE blk_oce         ! bulk parameters
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23)  
(   24)    !! Shared module variables
(   25)    LOGICAL, PUBLIC, PARAMETER ::   lk_ice_lim = .TRUE.    !: LIM ice model
(   26) 
(   27)    !!----------------------------------------------------------------------
(   28)    !! ice-ocean common variables
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) # 43
(   43)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !: field exchanges with ice model to ocean
(   44)       sst_io, sss_io , &  !: sea surface temperature (C) and salinity (PSU)
(   45)       u_io  , v_io   , &  !: velocity at ice surface (m/s)
(   46)       fsolar, fnsolar, &  !: solar and non-solar heat fluxes (W/m2)
(   47)       fsalt , fmass  , &  !: salt and freshwater fluxes
(   48)       ftaux , ftauy  , &  !: wind stresses
(   49)       gtaux , gtauy       !: wind stresses
(   50)    
(   51)    REAL(wp), PUBLIC ::   &  !:
(   52)       rdt_ice,           &  !: ice time step
(   53)       dtsd2                 !: ice time step divide by 2
(   54) 
(   55) 
(   56) # 62
(   62)    INTEGER, PUBLIC ::   &  !: namdom : space/time domain (namlist)
(   63)       nfice =  5           !: coupling frequency OPA ICELLN  nfice 
(   64) 
(   65)    !!----------------------------------------------------------------------
(   66) END MODULE ice_oce
