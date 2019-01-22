MODULE trdicp_oce
   !!======================================================================
   !!                   ***  MODULE trdicp_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdicp_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_trdtra'   or                         tracer trends diagnostics
   !!   'key_trddyn'                            momentum trends diagnostics
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce                 ! ocean parameters

   IMPLICIT NONE
   PUBLIC

   !! Namelist parameters
   !!----------------------------------------------------------------------
   INTEGER  ::      & !!: namdia :  diagnostics on dynamics and/or tracer trends
      ntrd  = 10 ,  &  !: time step frequency dynamics and tracers trends
      nctls =  0       !: control surface type for trends vertical integration

   !! Tracers trends diagnostics parameters
   !!---------------------------------------------------------------------
   INTEGER, PARAMETER ::            &  !: trends index
      jpttdlad = 1,   &  !: tracer horizontal advection
      jpttdzad = 2,   &  !: tracer vertical advection
      jpttdldf = 3,   &  !: tracer horizontal diffusion
      jpttdzdf = 4,   &  !: tracer vertical diffusion
      jpttdnpc = 5,   &  !: tracer non penetrative convection
      jpttddoe = 6,   &  !: tracer D.amping O.r vertical E.iv
      jpttdqsr = 7,   &  !: tracer penetrative solar radiation
      jpttdnsr = 8       !: tracer non solar radiation

   !! Momentum trends diagnostics parameters
   !!---------------------------------------------------------------------
   INTEGER, PARAMETER ::            &  !: trends index
      jpdtdhpg =  1,   &  !: dynamic hydrostatic pressure gradient 
      jpdtdkeg =  2,   &  !: dynamic kinetic energy gradient
      jpdtdrvo =  3,   &  !: dynamic relative vorticity
      jpdtdpvo =  4,   &  !: dynamic planetary vorticity
      jpdtdldf =  5,   &  !: dynamic lateral diffusion
      jpdtdzad =  6,   &  !: dynamic vertical advection
      jpdtdzdf =  7,   &  !: dynamic vertical diffusion
      jpdtdspg =  8,   &  !: dynamic surface pressure gradient
      jpdtddat =  9,   &  !: dynamic damping term
      jpdtdswf = 10,   &  !: dynamic surface wind forcing
      jpdtdbfr = 11       !: dynamic bottom friction 

   REAL, DIMENSION(jpi,jpj) ::   &  !:
      tldfbbl, sldfbbl,          &  ! Temperature/salinity lateral diffusion trends
      !                             ! in the BBL  
      tladbbl, sladbbl              ! Temperature/salinity lateral advection trends 
      !                             ! in the BBL

   REAL, DIMENSION(jpi,jpj,jpk) ::   &  !:
      tladi, sladi,                  &  ! Temp./sal. MUSCL OR TVD advection fluxes 
      !                                 ! terms along i- 
      tladj, sladj                      ! Temp./sal. MUSCL OR TVD advection fluxes 
      !                                 ! terms along j- 
#if defined key_ldfslp
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  !:
      uldftrd, vldftrd     !: lateral diffusion trend in isopycnal case
#endif
#if   defined key_trdtra   ||   defined key_trddyn   ||   defined key_esopa

   !! Variables used for diagnostics
   !!---------------------------------------------------------------------
   REAL(wp) ::   &  !:
      tvolt,     &  !: volume of the whole ocean computed at t-points
      tvolu,     &  !: volume of the whole ocean computed at u-points
      tvolv         !: volume of the whole ocean computed at v-points

   !! Tracers trends diagnostics variables
   !!---------------------------------------------------------------------
   REAL(wp), DIMENSION(10) ::   &  !:
      tmo, smo         !: tracers trends average 
      !                !  tmo(1) : horizontal advection
      !                !  tmo(2) : vertical advection
      !                !  tmo(3) : horizontal diffusion
      !                !  tmo(4) : vertical diffusion
      !                !  tmo(5) : static instability
      !                !  tmo(6) : damping OR vertical EIV
      !                !  tmo(7) : penetrative solar radiation (T only)
   REAL(wp), DIMENSION(10) ::   &  !:
      t2, s2           !: tracers square trends average 
      !                !  t2(1) : horizontal advection
      !                !  t2(2) : vertical advection
      !                !  t2(3) : horizontal diffusion
      !                !  t2(4) : vertical diffusion
      !                !  t2(5) : static instability
      !                !  t2(6) : damping OR vertical EIV
      !                !  t2(7) : penetrative solar radiation (T only)
   
   !! Momentum trends diagnostics variables
   !!---------------------------------------------------------------------
   REAL(wp), DIMENSION(11) ::   &  !:
      umo, vmo         !: momentum trends average 
      !                !  umo(1) : hydrostatic pressure gradient
      !                !  umo(2) : kinetic energy
      !                !  umo(3) : lateral diffusion geo-pot
      !                !  umo(4) : 
      !                !  umo(5) : lateral diffusion
      !                !  umo(6) : vertical advection
      !                !  umo(7) : vertical diffusion
      !                !  umo(8) : surface pressure gradient
      !                !  umo(9) : 

   REAL(wp), DIMENSION(10) ::   &  !:
      hke              !: momentum square trends average 
      !                !  hke(1) : horizontal advection
      !                !  hke(2) : vertical advection

   REAL(wp) ::   &  !:
      rpktrd,    &  !: potential to kinetic energy conversion
      peke          !: conversion potential energy - kinetic energy trend

#endif

  !!======================================================================
END MODULE trdicp_oce
