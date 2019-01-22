MODULE oce
   !!======================================================================
   !!                      ***  MODULE  oce  ***
   !! Ocean        :  dynamics and active tracers defined in memory 
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/oce.F90,v 1.7 2005/12/21 10:46:29 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   !! Physics and algorithm flags
   !! ---------------------------
   LOGICAL, PUBLIC ::   ln_dynhpg_imp   = .FALSE.  !: semi-implicite hpg flag

   !! dynamics and tracer fields
   !! --------------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      ! before !  now      !  after  !      ! the after trends becomes the fields
      ! fields !  fields   !  trends !      ! only in dyn(tra)_zdf and dyn(tra)_nxt
      ub       ,  un       ,  ua     ,   &  !: i-horizontal velocity (m/s)
      vb       ,  vn       ,  va     ,   &  !: j-horizontal velocity (m/s)
                  wn       ,             &  !: vertical velocity (m/s)
      rotb     ,  rotn     ,             &  !: relative vorticity (1/s)
      hdivb    ,  hdivn    ,             &  !: horizontal divergence (1/s)
      tb       ,  tn       ,  ta     ,   &  !: potential temperature (celcius)
      sb       ,  sn       ,  sa            !: salinity (psu)
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      rhd ,                              &  !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
      rhop,                              &  !: potential volumic mass (kg/m3)
      rn2                                   !: brunt-vaisala frequency (1/s2)

      !! advection scheme choice
      !! -----------------------
      CHARACTER(len=3), PUBLIC  ::   l_adv   !: 'ce2' centre scheme used
              !                              !: 'tvd' TVD scheme used
              !                              !: 'mus' MUSCL scheme used
              !                              !: 'mu2' MUSCL2 scheme used

   !! surface pressure gradient
   !! -------------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      spgu, spgv             !: horizontal surface pressure gradient

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

   !! free surface
   !! ------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      sshb, sshn,         &  !: before, now sea surface height (meters)
      hu  , hv  ,         &  !: depth at u- and v-points (meters)
      hur , hvr              !: inverse of u and v-points ocean depth (1/m)
#if defined key_obc
   REAL(wp), PUBLIC ::    &  !:
      obcsurftot       !: Total lateral surface of open boundaries
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      obcumask, obcvmask     !: u-, v- Force filtering mask for the open 
      !                      !  boundary condition on grad D
#endif

#if defined key_dynspg_rl   ||   defined key_esopa
   !! rigid-lid formulation
   !! ---------------------
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      bsfb, bsfn,         &  !: before, now barotropic streamfunction (m3/s)
      bsfd                   !: now trend of barotropic streamfunction (m3/s2)
#endif


!!DB: 2007.12.05
!!DB control parameters -- default values
   INTEGER, PUBLIC ::   &
        perpetual_forcing =  0,   &  ! switch for perpetual forcing (1=on) -- controls call to day()
        M2_ave = 0,       &  ! > 0 ====> daily ave of last M2 dt's 
	ioutput_ave = 0         ! flag for TSUV... output, averaged over ioutput_ave timesteps
                                ! and output every mod(kt-nit000+1,ioutput_ave) 
   REAL(wp), PUBLIC :: ramp     ! ramp for forcing, used in obcdta,tau_forced_*.h90
!!DBG 2009.10.20 
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   sw_rad

   !!----------------------------------------------------------------------
END MODULE oce
