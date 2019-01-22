


PGF90 (Version     10.5)          09/23/2010  09:48:59      page 1

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

Filename: oce.F90

(    1) MODULE oce
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  oce  ***
(    4)    !! Ocean        :  dynamics and active tracers defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(    9)    !!----------------------------------------------------------------------
(   10)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/oce.F90,v 1.7 2005/12/21 10:46:29 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE par_oce      ! ocean parameters
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! Physics and algorithm flags
(   21)    !! ---------------------------
(   22)    LOGICAL, PUBLIC ::   ln_dynhpg_imp   = .FALSE.  !: semi-implicite hpg flag
(   23) 
(   24)    !! dynamics and tracer fields
(   25)    !! --------------------------
(   26)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   27)       ! before !  now      !  after  !      ! the after trends becomes the fields
(   28)       ! fields !  fields   !  trends !      ! only in dyn(tra)_zdf and dyn(tra)_nxt
(   29)       ub       ,  un       ,  ua     ,   &  !: i-horizontal velocity (m/s)
(   30)       vb       ,  vn       ,  va     ,   &  !: j-horizontal velocity (m/s)
(   31)                   wn       ,             &  !: vertical velocity (m/s)
(   32)       rotb     ,  rotn     ,             &  !: relative vorticity (1/s)
(   33)       hdivb    ,  hdivn    ,             &  !: horizontal divergence (1/s)
(   34)       tb       ,  tn       ,  ta     ,   &  !: potential temperature (celcius)
(   35)       sb       ,  sn       ,  sa            !: salinity (psu)
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   37)       rhd ,                              &  !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
(   38)       rhop,                              &  !: potential volumic mass (kg/m3)
(   39)       rn2                                   !: brunt-vaisala frequency (1/s2)
(   40) 
(   41)       !! advection scheme choice
(   42)       !! -----------------------
(   43)       CHARACTER(len=3), PUBLIC  ::   l_adv   !: 'ce2' centre scheme used
(   44)               !                              !: 'tvd' TVD scheme used
(   45)               !                              !: 'mus' MUSCL scheme used
(   46)               !                              !: 'mu2' MUSCL2 scheme used






PGF90 (Version     10.5)          09/23/2010  09:48:59      page 2

(   47) 
(   48)    !! surface pressure gradient
(   49)    !! -------------------------
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   51)       spgu, spgv             !: horizontal surface pressure gradient
(   52) 
(   53) #if defined key_partial_steps     ||   defined key_esopa
(   54)    !! interpolated gradient
(   55)    !! ---------------------
(   56)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   57)       gtu, gsu, gru,      &  !: t-, s- and rd horizontal gradient at u- and 
(   58)       gtv, gsv, grv          !: v-points at bottom ocean level 
(   59) #else
(   60)    REAL(wp), PUBLIC ::   &   !:
(   61)       gtu, gsu, gru,      &  !: dummy scalars
(   62)       gtv, gsv, grv          !:
(   63) #endif
(   64) 
(   65)    !! free surface
(   66)    !! ------------
(   67)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   68)       sshb, sshn,         &  !: before, now sea surface height (meters)
(   69)       hu  , hv  ,         &  !: depth at u- and v-points (meters)
(   70)       hur , hvr              !: inverse of u and v-points ocean depth (1/m)
(   71) #if defined key_obc
(   72)    REAL(wp), PUBLIC ::    &  !:
(   73)       obcsurftot       !: Total lateral surface of open boundaries
(   74)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   75)       obcumask, obcvmask     !: u-, v- Force filtering mask for the open 
(   76)       !                      !  boundary condition on grad D
(   77) #endif
(   78) 
(   79) #if defined key_dynspg_rl   ||   defined key_esopa
(   80)    !! rigid-lid formulation
(   81)    !! ---------------------
(   82)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   83)       bsfb, bsfn,         &  !: before, now barotropic streamfunction (m3/s)
(   84)       bsfd                   !: now trend of barotropic streamfunction (m3/s2)
(   85) #endif
(   86) 
(   87) 
(   88) !!DB: 2007.12.05
(   89) !!DB control parameters -- default values
(   90)    INTEGER, PUBLIC ::   &
(   91)         perpetual_forcing =  0,   &  ! switch for perpetual forcing (1=on) -- controls call to day()
(   92)         M2_ave = 0,       &  ! > 0 ====> daily ave of last M2 dt's 
(   93) 	ioutput_ave = 0         ! flag for TSUV... output, averaged over ioutput_ave timesteps
(   94)                                 ! and output every mod(kt-nit000+1,ioutput_ave) 
(   95)    REAL(wp), PUBLIC :: ramp     ! ramp for forcing, used in obcdta,tau_forced_*.h90
(   96) !!DBG 2009.10.20 
(   97)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   sw_rad
(   98) 
(   99)    !!----------------------------------------------------------------------
(  100) END MODULE oce










PGF90 (Version     10.5)          09/23/2010  09:48:59      page 3

(    1) # 1 "oce.F90"
(    1) MODULE oce
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  oce  ***
(    4)    !! Ocean        :  dynamics and active tracers defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(    9)    !!----------------------------------------------------------------------
(   10)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/oce.F90,v 1.7 2005/12/21 10:46:29 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE par_oce      ! ocean parameters
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! Physics and algorithm flags
(   21)    !! ---------------------------
(   22)    LOGICAL, PUBLIC ::   ln_dynhpg_imp   = .FALSE.  !: semi-implicite hpg flag
(   23) 
(   24)    !! dynamics and tracer fields
(   25)    !! --------------------------
(   26)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   27)       ! before !  now      !  after  !      ! the after trends becomes the fields
(   28)       ! fields !  fields   !  trends !      ! only in dyn(tra)_zdf and dyn(tra)_nxt
(   29)       ub       ,  un       ,  ua     ,   &  !: i-horizontal velocity (m/s)
(   30)       vb       ,  vn       ,  va     ,   &  !: j-horizontal velocity (m/s)
(   31)                   wn       ,             &  !: vertical velocity (m/s)
(   32)       rotb     ,  rotn     ,             &  !: relative vorticity (1/s)
(   33)       hdivb    ,  hdivn    ,             &  !: horizontal divergence (1/s)
(   34)       tb       ,  tn       ,  ta     ,   &  !: potential temperature (celcius)
(   35)       sb       ,  sn       ,  sa            !: salinity (psu)
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   37)       rhd ,                              &  !: in situ density anomalie rhd=(rho-rau0)/rau0 (no units)
(   38)       rhop,                              &  !: potential volumic mass (kg/m3)
(   39)       rn2                                   !: brunt-vaisala frequency (1/s2)
(   40) 
(   41)       !! advection scheme choice
(   42)       !! -----------------------
(   43)       CHARACTER(len=3), PUBLIC  ::   l_adv   !: 'ce2' centre scheme used
(   44)               !                              !: 'tvd' TVD scheme used
(   45)               !                              !: 'mus' MUSCL scheme used
(   46)               !                              !: 'mu2' MUSCL2 scheme used
(   47) 
(   48)    !! surface pressure gradient
(   49)    !! -------------------------
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   51)       spgu, spgv             !: horizontal surface pressure gradient
(   52) 
(   53) # 60
(   60)    REAL(wp), PUBLIC ::   &   !:
(   61)       gtu, gsu, gru,      &  !: dummy scalars
(   62)       gtv, gsv, grv          !:
(   63) 






PGF90 (Version     10.5)          09/23/2010  09:48:59      page 4

(   64) # 65
(   65)    !! free surface
(   66)    !! ------------
(   67)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   68)       sshb, sshn,         &  !: before, now sea surface height (meters)
(   69)       hu  , hv  ,         &  !: depth at u- and v-points (meters)
(   70)       hur , hvr              !: inverse of u and v-points ocean depth (1/m)
(   71) # 72
(   72)    REAL(wp), PUBLIC ::    &  !:
(   73)       obcsurftot       !: Total lateral surface of open boundaries
(   74)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   75)       obcumask, obcvmask     !: u-, v- Force filtering mask for the open 
(   76)       !                      !  boundary condition on grad D
(   77) 
(   78) 
(   79) 
(   80) # 88
(   88) !!DB: 2007.12.05
(   89) !!DB control parameters -- default values
(   90)    INTEGER, PUBLIC ::   &
(   91)         perpetual_forcing =  0,   &  ! switch for perpetual forcing (1=on) -- controls call to day()
(   92)         M2_ave = 0,       &  ! > 0 ====> daily ave of last M2 dt's 
(   93) 	ioutput_ave = 0         ! flag for TSUV... output, averaged over ioutput_ave timesteps
(   94)                                 ! and output every mod(kt-nit000+1,ioutput_ave) 
(   95)    REAL(wp), PUBLIC :: ramp     ! ramp for forcing, used in obcdta,tau_forced_*.h90
(   96) !!DBG 2009.10.20 
(   97)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   sw_rad
(   98) 
(   99)    !!----------------------------------------------------------------------
(  100) END MODULE oce
