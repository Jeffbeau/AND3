MODULE ice
   !!======================================================================
   !!                        ***  MODULE ice  ***
   !! Sea Ice physics:  diagnostics variables of ice defined in memory
   !!=====================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   LIM sea-ice model
   !!----------------------------------------------------------------------
   !! History :
   !!   2.0  !  03-08  (C. Ethe)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/ice.F90,v 1.5 2006/03/21 08:38:39 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_ice          ! LIM sea-ice parameters

   IMPLICIT NONE
   PRIVATE

   !! * Share Module variables
   INTEGER , PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
      nbiter = 1      ,  &  !: number of sub-time steps for relaxation
      nbitdr = 250          !: maximum number of iterations for relaxation

   REAL(wp), PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
      epsd   = 1.0e-20,  &  !: tolerance parameter for dynamic
      alpha  = 0.5    ,  &  !: coefficient for semi-implicit coriolis
      dm     = 0.6e+03,  &  !: diffusion constant for dynamics
      om     = 0.5    ,  &  !: relaxation constant
      resl   = 5.0e-05,  &  !: maximum value for the residual of relaxation
      cw     = 5.0e-03,  &  !: drag coefficient for oceanic stress
      angvg  = 0.e0   ,  &  !: turning angle for oceanic stress
      pstar  = 1.0e+04,  &  !: first bulk-rheology parameter
      c_rhg  = 20.e0  ,  &  !: second bulk-rhelogy parameter
      etamn  = 0.e+07,   &  !: minimun value for viscosity
      creepl = 2.e-08,   &  !: creep limit
      ecc    = 2.e0   ,  &  !: eccentricity of the elliptical yield curve
      ahi0   = 350.e0       !: sea-ice hor. eddy diffusivity coeff. (m2/s)

   REAL(wp), PUBLIC ::   &  !:
      usecc2          ,  &  !:  = 1.0 / ( ecc * ecc )
      rhoco           ,  &  !: = rau0 * cw
      sangvg, cangvg  ,  &  !: sin and cos of the turning angle for ocean stress
      pstarh                !: pstar / 2.0

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
      u_oce, v_oce,      &  !: surface ocean velocity used in ice dynamics
      ahiu , ahiv ,      &  !: hor. diffusivity coeff. at ocean U- and V-points (m2/s)
      pahu , pahv ,      &  !: ice hor. eddy diffusivity coef. at ocean U- and V-points
      hsnm , hicm ,      &  !: mean snow and ice thicknesses
      ust2s                 !: friction velocity

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
        sst_ini,         &  !: sst read from a file for ice model initialization 
        sss_ini             !: sss read from a file for ice model initialization 

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      firic  ,   &  !: IR flux over the ice (only used for outputs)
      fcsic  ,   &  !: Sensible heat flux over the ice (only used for outputs)
      fleic  ,   &  !: Latent heat flux over the ice (only used for outputs)
      qlatic ,   &  !: latent flux
      rdvosif,   &  !: Variation of volume at surface (only used for outputs)
      rdvobif,   &  !: Variation of ice volume at the bottom ice (only used for outputs)
      fdvolif,   &  !: Total variation of ice volume (only used for outputs)
      rdvonif,   &  !: Lateral Variation of ice volume (only used for outputs)
      sist   ,   &  !: Sea-Ice Surface Temperature (Kelvin ??? degree ??? I don't know)
      tfu    ,   &  !: Melting point temperature of sea water
      hsnif  ,   &  !: Snow thickness
      hicif  ,   &  !: Ice thickness
      hicifp ,   &  !: Ice production/melting
      frld   ,   &  !: Leads fraction = 1-a/totalarea
      phicif ,   &  !: ice thickness  at previous time 
      pfrld  ,   &  !: Leads fraction at previous time  
      qstoif ,   &  !: Energy stored in the brine pockets
      fbif   ,   &  !: Heat flux at the ice base
      rdmsnif,   &  !: Variation of snow mass
      rdmicif,   &  !: Variation of ice mass
      qldif  ,   &  !: heat balance of the lead (or of the open ocean)
      qcmif  ,   &  !: Energy needed to bring the ocean surface layer until its freezing 
      fdtcn  ,   &  !: net downward heat flux from the ice to the ocean
      qdtcn  ,   &  !: energy from the ice to the ocean
      !             !  point (at a factor 2)
      thcm   ,   &  !: part of the solar energy used in the lead heat budget
      fstric ,   &  !: Solar flux transmitted trough the ice
      ffltbif,   &  !: Array linked with the max heat contained in brine pockets (?)
      fscmbq ,   &  !: Linked with the solar flux below the ice (?)
      fsbbq  ,   &  !: Also linked with the solar flux below the ice (?)
      qfvbq  ,   &  !: Array used to store energy in case of toral lateral ablation (?)
      dmgwi         !: Variation of the mass of snow ice

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      albege ,   &  !: Albedo of the snow or ice (only for outputs)
      albecn ,   &  !: Albedo of the ocean (only for outputs)
      tauc   ,   &  !: Cloud optical depth
      sdvt          !: u*^2/(Stress/density)


   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      u_ice, v_ice,   &  !: two components of the ice velocity (m/s)
      tio_u, tio_v       !: two components of the ice-ocean stress (N/m2)

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpsmax) ::   &  !:
      scal0              !: ???

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jplayersp1) ::   &  !:
      tbif          !: Temperature inside the ice/snow layer

   REAL(wp), DIMENSION(jpi,jpj,0:jpkmax+1) ::    &  !:
      reslum        !: Relative absorption of solar radiation in each ocean level

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
         sxice, syice, sxxice, syyice, sxyice,      &  !: moments for advection
         sxsn,  sysn,  sxxsn,  syysn,  sxysn,       &  !:
         sxa,   sya,   sxxa,   syya,   sxya,        &  !:
         sxc0,  syc0,  sxxc0,  syyc0,  sxyc0,       &  !:
         sxc1,  syc1,  sxxc1,  syyc1,  sxyc1,       &  !:
         sxc2,  syc2,  sxxc2,  syyc2,  sxyc2,       &  !:
         sxst,  syst,  sxxst,  syyst,  sxyst           !:

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module            NO LIM sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE ice
