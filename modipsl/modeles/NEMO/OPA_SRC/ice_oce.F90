MODULE ice_oce
   !!======================================================================
   !!                 ***  MODULE  ice_oce  ***
   !! Ocean - ice  :  ice variables defined in memory 
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ice_oce.F90,v 1.3 2005/03/27 18:34:46 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'   :                                     LIM ice model
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters
   USE blk_oce         ! bulk parameters

   IMPLICIT NONE
   PRIVATE
 
   !! Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_ice_lim = .TRUE.    !: LIM ice model

   !!----------------------------------------------------------------------
   !! ice-ocean common variables
   !!----------------------------------------------------------------------
# if defined key_coupled
   REAL(wp), PUBLIC, DIMENSION(jpiglo,jpjglo) ::   &  !: cumulated fields
      fqsr_oce ,      &   !: Net short wave heat flux on free ocean 
      fqsr_ice ,      &   !: Net short wave het flux on sea ice 
      fqnsr_oce,      &   !: Net longwave heat flux on free ocean
      fqnsr_ice,      &   !: Net longwave heat flux on sea ice
      fdqns_ice,      &   !: Derivative of non solar heat flux on sea ice
      ftprecip ,      &   !: Water flux (liquid precipitation - evaporation) 
      fsprecip ,      &   !: Solid (snow) precipitation
      frunoff  ,      &   !: runoff
      fcalving            !: Iceberg calving 
# endif

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !: field exchanges with ice model to ocean
      sst_io, sss_io , &  !: sea surface temperature (C) and salinity (PSU)
      u_io  , v_io   , &  !: velocity at ice surface (m/s)
      fsolar, fnsolar, &  !: solar and non-solar heat fluxes (W/m2)
      fsalt , fmass  , &  !: salt and freshwater fluxes
      ftaux , ftauy  , &  !: wind stresses
      gtaux , gtauy       !: wind stresses
   
   REAL(wp), PUBLIC ::   &  !:
      rdt_ice,           &  !: ice time step
      dtsd2                 !: ice time step divide by 2

#else
   !!----------------------------------------------------------------------
   !!   Default option                                 NO LIM sea-ice model
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_ice_lim = .FALSE.  !: No LIM ice model
#endif

   INTEGER, PUBLIC ::   &  !: namdom : space/time domain (namlist)
      nfice =  5           !: coupling frequency OPA ICELLN  nfice 

   !!----------------------------------------------------------------------
END MODULE ice_oce
