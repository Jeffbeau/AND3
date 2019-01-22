


PGF90 (Version     10.5)          09/23/2010  09:49:02      page 1

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

Filename: ice.F90

(    1) MODULE ice
(    2)    !!======================================================================
(    3)    !!                        ***  MODULE ice  ***
(    4)    !! Sea Ice physics:  diagnostics variables of ice defined in memory
(    5)    !!=====================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !! History :
(   11)    !!   2.0  !  03-08  (C. Ethe)  F90: Free form and module
(   12)    !!----------------------------------------------------------------------
(   13)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/ice.F90,v 1.5 2006/03/21 08:38:39 opalod Exp $
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_ice          ! LIM sea-ice parameters
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Share Module variables
(   24)    INTEGER , PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
(   25)       nbiter = 1      ,  &  !: number of sub-time steps for relaxation
(   26)       nbitdr = 250          !: maximum number of iterations for relaxation
(   27) 
(   28)    REAL(wp), PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
(   29)       epsd   = 1.0e-20,  &  !: tolerance parameter for dynamic
(   30)       alpha  = 0.5    ,  &  !: coefficient for semi-implicit coriolis
(   31)       dm     = 0.6e+03,  &  !: diffusion constant for dynamics
(   32)       om     = 0.5    ,  &  !: relaxation constant
(   33)       resl   = 5.0e-05,  &  !: maximum value for the residual of relaxation
(   34)       cw     = 5.0e-03,  &  !: drag coefficient for oceanic stress
(   35)       angvg  = 0.e0   ,  &  !: turning angle for oceanic stress
(   36)       pstar  = 1.0e+04,  &  !: first bulk-rheology parameter
(   37)       c_rhg  = 20.e0  ,  &  !: second bulk-rhelogy parameter
(   38)       etamn  = 0.e+07,   &  !: minimun value for viscosity
(   39)       creepl = 2.e-08,   &  !: creep limit
(   40)       ecc    = 2.e0   ,  &  !: eccentricity of the elliptical yield curve
(   41)       ahi0   = 350.e0       !: sea-ice hor. eddy diffusivity coeff. (m2/s)
(   42) 
(   43)    REAL(wp), PUBLIC ::   &  !:
(   44)       usecc2          ,  &  !:  = 1.0 / ( ecc * ecc )
(   45)       rhoco           ,  &  !: = rau0 * cw
(   46)       sangvg, cangvg  ,  &  !: sin and cos of the turning angle for ocean stress






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 2

(   47)       pstarh                !: pstar / 2.0
(   48) 
(   49)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
(   50)       u_oce, v_oce,      &  !: surface ocean velocity used in ice dynamics
(   51)       ahiu , ahiv ,      &  !: hor. diffusivity coeff. at ocean U- and V-points (m2/s)
(   52)       pahu , pahv ,      &  !: ice hor. eddy diffusivity coef. at ocean U- and V-points
(   53)       hsnm , hicm ,      &  !: mean snow and ice thicknesses
(   54)       ust2s                 !: friction velocity
(   55) 
(   56)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
(   57)         sst_ini,         &  !: sst read from a file for ice model initialization 
(   58)         sss_ini             !: sss read from a file for ice model initialization 
(   59) 
(   60)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   61)       firic  ,   &  !: IR flux over the ice (only used for outputs)
(   62)       fcsic  ,   &  !: Sensible heat flux over the ice (only used for outputs)
(   63)       fleic  ,   &  !: Latent heat flux over the ice (only used for outputs)
(   64)       qlatic ,   &  !: latent flux
(   65)       rdvosif,   &  !: Variation of volume at surface (only used for outputs)
(   66)       rdvobif,   &  !: Variation of ice volume at the bottom ice (only used for outputs)
(   67)       fdvolif,   &  !: Total variation of ice volume (only used for outputs)
(   68)       rdvonif,   &  !: Lateral Variation of ice volume (only used for outputs)
(   69)       sist   ,   &  !: Sea-Ice Surface Temperature (Kelvin ??? degree ??? I don't know)
(   70)       tfu    ,   &  !: Melting point temperature of sea water
(   71)       hsnif  ,   &  !: Snow thickness
(   72)       hicif  ,   &  !: Ice thickness
(   73)       hicifp ,   &  !: Ice production/melting
(   74)       frld   ,   &  !: Leads fraction = 1-a/totalarea
(   75)       phicif ,   &  !: ice thickness  at previous time 
(   76)       pfrld  ,   &  !: Leads fraction at previous time  
(   77)       qstoif ,   &  !: Energy stored in the brine pockets
(   78)       fbif   ,   &  !: Heat flux at the ice base
(   79)       rdmsnif,   &  !: Variation of snow mass
(   80)       rdmicif,   &  !: Variation of ice mass
(   81)       qldif  ,   &  !: heat balance of the lead (or of the open ocean)
(   82)       qcmif  ,   &  !: Energy needed to bring the ocean surface layer until its freezing 
(   83)       fdtcn  ,   &  !: net downward heat flux from the ice to the ocean
(   84)       qdtcn  ,   &  !: energy from the ice to the ocean
(   85)       !             !  point (at a factor 2)
(   86)       thcm   ,   &  !: part of the solar energy used in the lead heat budget
(   87)       fstric ,   &  !: Solar flux transmitted trough the ice
(   88)       ffltbif,   &  !: Array linked with the max heat contained in brine pockets (?)
(   89)       fscmbq ,   &  !: Linked with the solar flux below the ice (?)
(   90)       fsbbq  ,   &  !: Also linked with the solar flux below the ice (?)
(   91)       qfvbq  ,   &  !: Array used to store energy in case of toral lateral ablation (?)
(   92)       dmgwi         !: Variation of the mass of snow ice
(   93) 
(   94)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   95)       albege ,   &  !: Albedo of the snow or ice (only for outputs)
(   96)       albecn ,   &  !: Albedo of the ocean (only for outputs)
(   97)       tauc   ,   &  !: Cloud optical depth
(   98)       sdvt          !: u*^2/(Stress/density)
(   99) 
(  100) 
(  101)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(  102)       u_ice, v_ice,   &  !: two components of the ice velocity (m/s)
(  103)       tio_u, tio_v       !: two components of the ice-ocean stress (N/m2)
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 3

(  105)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpsmax) ::   &  !:
(  106)       scal0              !: ???
(  107) 
(  108)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jplayersp1) ::   &  !:
(  109)       tbif          !: Temperature inside the ice/snow layer
(  110) 
(  111)    REAL(wp), DIMENSION(jpi,jpj,0:jpkmax+1) ::    &  !:
(  112)       reslum        !: Relative absorption of solar radiation in each ocean level
(  113) 
(  114)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(  115)          sxice, syice, sxxice, syyice, sxyice,      &  !: moments for advection
(  116)          sxsn,  sysn,  sxxsn,  syysn,  sxysn,       &  !:
(  117)          sxa,   sya,   sxxa,   syya,   sxya,        &  !:
(  118)          sxc0,  syc0,  sxxc0,  syyc0,  sxyc0,       &  !:
(  119)          sxc1,  syc1,  sxxc1,  syyc1,  sxyc1,       &  !:
(  120)          sxc2,  syc2,  sxxc2,  syyc2,  sxyc2,       &  !:
(  121)          sxst,  syst,  sxxst,  syyst,  sxyst           !:
(  122) 
(  123) #else
(  124)    !!----------------------------------------------------------------------
(  125)    !!   Default option         Empty module            NO LIM sea-ice model
(  126)    !!----------------------------------------------------------------------
(  127) #endif
(  128) 
(  129)    !!======================================================================
(  130) END MODULE ice






































PGF90 (Version     10.5)          09/23/2010  09:49:02      page 4

(    1) # 1 "ice.F90"
(    1) MODULE ice
(    2)    !!======================================================================
(    3)    !!                        ***  MODULE ice  ***
(    4)    !! Sea Ice physics:  diagnostics variables of ice defined in memory
(    5)    !!=====================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !! History :
(   11)    !!   2.0  !  03-08  (C. Ethe)  F90: Free form and module
(   12)    !!----------------------------------------------------------------------
(   13)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/ice.F90,v 1.5 2006/03/21 08:38:39 opalod Exp $
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_ice          ! LIM sea-ice parameters
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Share Module variables
(   24)    INTEGER , PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
(   25)       nbiter = 1      ,  &  !: number of sub-time steps for relaxation
(   26)       nbitdr = 250          !: maximum number of iterations for relaxation
(   27) 
(   28)    REAL(wp), PUBLIC ::   & !!: ** ice-dynamic namelist (namicedyn) **
(   29)       epsd   = 1.0e-20,  &  !: tolerance parameter for dynamic
(   30)       alpha  = 0.5    ,  &  !: coefficient for semi-implicit coriolis
(   31)       dm     = 0.6e+03,  &  !: diffusion constant for dynamics
(   32)       om     = 0.5    ,  &  !: relaxation constant
(   33)       resl   = 5.0e-05,  &  !: maximum value for the residual of relaxation
(   34)       cw     = 5.0e-03,  &  !: drag coefficient for oceanic stress
(   35)       angvg  = 0.e0   ,  &  !: turning angle for oceanic stress
(   36)       pstar  = 1.0e+04,  &  !: first bulk-rheology parameter
(   37)       c_rhg  = 20.e0  ,  &  !: second bulk-rhelogy parameter
(   38)       etamn  = 0.e+07,   &  !: minimun value for viscosity
(   39)       creepl = 2.e-08,   &  !: creep limit
(   40)       ecc    = 2.e0   ,  &  !: eccentricity of the elliptical yield curve
(   41)       ahi0   = 350.e0       !: sea-ice hor. eddy diffusivity coeff. (m2/s)
(   42) 
(   43)    REAL(wp), PUBLIC ::   &  !:
(   44)       usecc2          ,  &  !:  = 1.0 / ( ecc * ecc )
(   45)       rhoco           ,  &  !: = rau0 * cw
(   46)       sangvg, cangvg  ,  &  !: sin and cos of the turning angle for ocean stress
(   47)       pstarh                !: pstar / 2.0
(   48) 
(   49)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
(   50)       u_oce, v_oce,      &  !: surface ocean velocity used in ice dynamics
(   51)       ahiu , ahiv ,      &  !: hor. diffusivity coeff. at ocean U- and V-points (m2/s)
(   52)       pahu , pahv ,      &  !: ice hor. eddy diffusivity coef. at ocean U- and V-points
(   53)       hsnm , hicm ,      &  !: mean snow and ice thicknesses
(   54)       ust2s                 !: friction velocity
(   55) 
(   56)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::  &  !:
(   57)         sst_ini,         &  !: sst read from a file for ice model initialization 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 5

(   58)         sss_ini             !: sss read from a file for ice model initialization 
(   59) 
(   60)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   61)       firic  ,   &  !: IR flux over the ice (only used for outputs)
(   62)       fcsic  ,   &  !: Sensible heat flux over the ice (only used for outputs)
(   63)       fleic  ,   &  !: Latent heat flux over the ice (only used for outputs)
(   64)       qlatic ,   &  !: latent flux
(   65)       rdvosif,   &  !: Variation of volume at surface (only used for outputs)
(   66)       rdvobif,   &  !: Variation of ice volume at the bottom ice (only used for outputs)
(   67)       fdvolif,   &  !: Total variation of ice volume (only used for outputs)
(   68)       rdvonif,   &  !: Lateral Variation of ice volume (only used for outputs)
(   69)       sist   ,   &  !: Sea-Ice Surface Temperature (Kelvin ??? degree ??? I don't know)
(   70)       tfu    ,   &  !: Melting point temperature of sea water
(   71)       hsnif  ,   &  !: Snow thickness
(   72)       hicif  ,   &  !: Ice thickness
(   73)       hicifp ,   &  !: Ice production/melting
(   74)       frld   ,   &  !: Leads fraction = 1-a/totalarea
(   75)       phicif ,   &  !: ice thickness  at previous time 
(   76)       pfrld  ,   &  !: Leads fraction at previous time  
(   77)       qstoif ,   &  !: Energy stored in the brine pockets
(   78)       fbif   ,   &  !: Heat flux at the ice base
(   79)       rdmsnif,   &  !: Variation of snow mass
(   80)       rdmicif,   &  !: Variation of ice mass
(   81)       qldif  ,   &  !: heat balance of the lead (or of the open ocean)
(   82)       qcmif  ,   &  !: Energy needed to bring the ocean surface layer until its freezing 
(   83)       fdtcn  ,   &  !: net downward heat flux from the ice to the ocean
(   84)       qdtcn  ,   &  !: energy from the ice to the ocean
(   85)       !             !  point (at a factor 2)
(   86)       thcm   ,   &  !: part of the solar energy used in the lead heat budget
(   87)       fstric ,   &  !: Solar flux transmitted trough the ice
(   88)       ffltbif,   &  !: Array linked with the max heat contained in brine pockets (?)
(   89)       fscmbq ,   &  !: Linked with the solar flux below the ice (?)
(   90)       fsbbq  ,   &  !: Also linked with the solar flux below the ice (?)
(   91)       qfvbq  ,   &  !: Array used to store energy in case of toral lateral ablation (?)
(   92)       dmgwi         !: Variation of the mass of snow ice
(   93) 
(   94)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   95)       albege ,   &  !: Albedo of the snow or ice (only for outputs)
(   96)       albecn ,   &  !: Albedo of the ocean (only for outputs)
(   97)       tauc   ,   &  !: Cloud optical depth
(   98)       sdvt          !: u*^2/(Stress/density)
(   99) 
(  100) 
(  101)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(  102)       u_ice, v_ice,   &  !: two components of the ice velocity (m/s)
(  103)       tio_u, tio_v       !: two components of the ice-ocean stress (N/m2)
(  104) 
(  105)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpsmax) ::   &  !:
(  106)       scal0              !: ???
(  107) 
(  108)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jplayersp1) ::   &  !:
(  109)       tbif          !: Temperature inside the ice/snow layer
(  110) 
(  111)    REAL(wp), DIMENSION(jpi,jpj,0:jpkmax+1) ::    &  !:
(  112)       reslum        !: Relative absorption of solar radiation in each ocean level
(  113) 
(  114)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(  115)          sxice, syice, sxxice, syyice, sxyice,      &  !: moments for advection






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 6

(  116)          sxsn,  sysn,  sxxsn,  syysn,  sxysn,       &  !:
(  117)          sxa,   sya,   sxxa,   syya,   sxya,        &  !:
(  118)          sxc0,  syc0,  sxxc0,  syyc0,  sxyc0,       &  !:
(  119)          sxc1,  syc1,  sxxc1,  syyc1,  sxyc1,       &  !:
(  120)          sxc2,  syc2,  sxxc2,  syyc2,  sxyc2,       &  !:
(  121)          sxst,  syst,  sxxst,  syyst,  sxyst           !:
(  122) 
(  123) 
(  124) # 129
(  129)    !!======================================================================
(  130) END MODULE ice
