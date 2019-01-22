MODULE flx_oce
   !!======================================================================
   !!                 ***  MODULE  flx_oce  ***
   !!        parameter and  variables defined in memory in forced mode
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (C. Ethe)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flx_oce.F90,v 1.5 2005/03/27 18:35:13 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce          ! ocean parameters

   IMPLICIT NONE
   PRIVATE
   

   !!----------------------------------------------------------------------
   !! fluxes common variables
   !!----------------------------------------------------------------------
#if defined key_flx_forced_daily
   !!----------------------------------------------------------------------
   !! 'key_flx_forced_daily'
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      p_qt ,        &   !: total heat flux ( solar + non solar)
      p_qsr,        &   !: solar heat flux
      p_emp             !: evaporation minus precipitation            
#elif defined key_flx_forced_monthly   
!byoung
    !----------------------------------------------------------------
    !'key_flx_forced_monthly'
    !-------------------------------------------------------------
    REAL(wp), PUBLIC,DIMENSION(jpi,jpj):: &
    p_bqt,   &  !:total heat flux
    p_bqsr,  &  !:solar heat flux
    p_bemp      !:evaporation minus precipitation

#elif defined key_ice_lim || defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
      qsr_ice  ,      &  !: solar flux over ice
      qsr_oce  ,      &  !: solar flux over ocean
      qnsr_oce ,      &  !: total non solar heat flux (Longwave downward radiation) over ocean 
      qnsr_ice ,      &  !: total non solar heat flux (Longwave downward radiation) over ice
      tprecip  ,      &  !: total precipitation ( or liquid precip minus evaporation in coupled mode)
      sprecip  ,      &  !: solid (snow) precipitation
      dqns_ice ,      &  !: total non solar sensibility over ice (LW+SEN+LA)
      tn_ice   ,      &  !: ice surface temperature
      evap     ,      &  !: evaporation over ocean
      fr1_i0   ,      &  !: 1st part of the fraction of sol. rad.  which penetrate inside the ice cover
      fr2_i0   ,      &  !: 2nd part of the fraction of sol. rad.  which penetrate inside the ice cover 
#if ! defined key_coupled
      qla_ice  ,      &  !: latent flux over ice  
      dqla_ice           !: latent sensibility over ice
#else
      rrunoff  ,      &  !: runoff
      calving  ,      &  !: calving
      alb_ice            !: albedo of ice      
#endif

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   
#endif

   !!----------------------------------------------------------------------
END MODULE flx_oce
