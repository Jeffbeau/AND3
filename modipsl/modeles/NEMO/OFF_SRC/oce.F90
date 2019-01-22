MODULE oce
   !!======================================================================
   !!                      ***  MODULE  oce  ***
   !! Ocean        :  dynamics and active tracers defined in memory 
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/oce.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   !! dynamics and tracer fields
   !! --------------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      ! before !  now      !  after  !      ! the after trends becomes the fields
      ! fields !  fields   !  trends !      ! only in dyn(tra)_zdf and dyn(tra)_nxt
                  un       ,  ua     ,   &  !: i-horizontal velocity (m/s)
                  vn       ,  va     ,   &  !: j-horizontal velocity (m/s)
                  wn       ,             &  !: vertical velocity (m/s)
                  hdivn    ,             &  !: horizontal divergence (1/s)
                  tn       ,  ta     ,   &  !: potential temperature (celcius)
                  sn       ,  sa            !: salinity (psu)
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      rhd ,                              &  !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
      rhop,                              &  !: potential volumic mass (kg/m3)
      rn2                                   !: brunt-vaisala frequency (1/s2)

#if defined key_partial_steps     ||   defined key_esopa
   !! interpolated gradient
   !! ---------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      gtu, gsu, gru,      &  !: t-, s- and rd horizontal gradient at u- and 
      gtv, gsv, grv          !: v-points at bottom ocean level 
#else
   REAL(wp), PUBLIC ::   &   !:
      gtu, gsu, gru,      &  !: dummy scalars
      gtv, gsv, grv          !:
#endif

   !!----------------------------------------------------------------------
END MODULE oce
