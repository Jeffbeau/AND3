


PGF90 (Version     10.5)          09/23/2010  09:49:08      page 1

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

Filename: trdicp_oce.F90

(    1) MODULE trdicp_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdicp_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdicp_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_trdtra'   or                         tracer trends diagnostics
(   14)    !!   'key_trddyn'                            momentum trends diagnostics
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce                 ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PUBLIC
(   21) 
(   22)    !! Namelist parameters
(   23)    !!----------------------------------------------------------------------
(   24)    INTEGER  ::      & !!: namdia :  diagnostics on dynamics and/or tracer trends
(   25)       ntrd  = 10 ,  &  !: time step frequency dynamics and tracers trends
(   26)       nctls =  0       !: control surface type for trends vertical integration
(   27) 
(   28)    !! Tracers trends diagnostics parameters
(   29)    !!---------------------------------------------------------------------
(   30)    INTEGER, PARAMETER ::            &  !: trends index
(   31)       jpttdlad = 1,   &  !: tracer horizontal advection
(   32)       jpttdzad = 2,   &  !: tracer vertical advection
(   33)       jpttdldf = 3,   &  !: tracer horizontal diffusion
(   34)       jpttdzdf = 4,   &  !: tracer vertical diffusion
(   35)       jpttdnpc = 5,   &  !: tracer non penetrative convection
(   36)       jpttddoe = 6,   &  !: tracer D.amping O.r vertical E.iv
(   37)       jpttdqsr = 7,   &  !: tracer penetrative solar radiation
(   38)       jpttdnsr = 8       !: tracer non solar radiation
(   39) 
(   40)    !! Momentum trends diagnostics parameters
(   41)    !!---------------------------------------------------------------------
(   42)    INTEGER, PARAMETER ::            &  !: trends index
(   43)       jpdtdhpg =  1,   &  !: dynamic hydrostatic pressure gradient 
(   44)       jpdtdkeg =  2,   &  !: dynamic kinetic energy gradient
(   45)       jpdtdrvo =  3,   &  !: dynamic relative vorticity
(   46)       jpdtdpvo =  4,   &  !: dynamic planetary vorticity






PGF90 (Version     10.5)          09/23/2010  09:49:08      page 2

(   47)       jpdtdldf =  5,   &  !: dynamic lateral diffusion
(   48)       jpdtdzad =  6,   &  !: dynamic vertical advection
(   49)       jpdtdzdf =  7,   &  !: dynamic vertical diffusion
(   50)       jpdtdspg =  8,   &  !: dynamic surface pressure gradient
(   51)       jpdtddat =  9,   &  !: dynamic damping term
(   52)       jpdtdswf = 10,   &  !: dynamic surface wind forcing
(   53)       jpdtdbfr = 11       !: dynamic bottom friction 
(   54) 
(   55)    REAL, DIMENSION(jpi,jpj) ::   &  !:
(   56)       tldfbbl, sldfbbl,          &  ! Temperature/salinity lateral diffusion trends
(   57)       !                             ! in the BBL  
(   58)       tladbbl, sladbbl              ! Temperature/salinity lateral advection trends 
(   59)       !                             ! in the BBL
(   60) 
(   61)    REAL, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   62)       tladi, sladi,                  &  ! Temp./sal. MUSCL OR TVD advection fluxes 
(   63)       !                                 ! terms along i- 
(   64)       tladj, sladj                      ! Temp./sal. MUSCL OR TVD advection fluxes 
(   65)       !                                 ! terms along j- 
(   66) #if defined key_ldfslp
(   67)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  !:
(   68)       uldftrd, vldftrd     !: lateral diffusion trend in isopycnal case
(   69) #endif
(   70) #if   defined key_trdtra   ||   defined key_trddyn   ||   defined key_esopa
(   71) 
(   72)    !! Variables used for diagnostics
(   73)    !!---------------------------------------------------------------------
(   74)    REAL(wp) ::   &  !:
(   75)       tvolt,     &  !: volume of the whole ocean computed at t-points
(   76)       tvolu,     &  !: volume of the whole ocean computed at u-points
(   77)       tvolv         !: volume of the whole ocean computed at v-points
(   78) 
(   79)    !! Tracers trends diagnostics variables
(   80)    !!---------------------------------------------------------------------
(   81)    REAL(wp), DIMENSION(10) ::   &  !:
(   82)       tmo, smo         !: tracers trends average 
(   83)       !                !  tmo(1) : horizontal advection
(   84)       !                !  tmo(2) : vertical advection
(   85)       !                !  tmo(3) : horizontal diffusion
(   86)       !                !  tmo(4) : vertical diffusion
(   87)       !                !  tmo(5) : static instability
(   88)       !                !  tmo(6) : damping OR vertical EIV
(   89)       !                !  tmo(7) : penetrative solar radiation (T only)
(   90)    REAL(wp), DIMENSION(10) ::   &  !:
(   91)       t2, s2           !: tracers square trends average 
(   92)       !                !  t2(1) : horizontal advection
(   93)       !                !  t2(2) : vertical advection
(   94)       !                !  t2(3) : horizontal diffusion
(   95)       !                !  t2(4) : vertical diffusion
(   96)       !                !  t2(5) : static instability
(   97)       !                !  t2(6) : damping OR vertical EIV
(   98)       !                !  t2(7) : penetrative solar radiation (T only)
(   99)    
(  100)    !! Momentum trends diagnostics variables
(  101)    !!---------------------------------------------------------------------
(  102)    REAL(wp), DIMENSION(11) ::   &  !:
(  103)       umo, vmo         !: momentum trends average 
(  104)       !                !  umo(1) : hydrostatic pressure gradient






PGF90 (Version     10.5)          09/23/2010  09:49:08      page 3

(  105)       !                !  umo(2) : kinetic energy
(  106)       !                !  umo(3) : lateral diffusion geo-pot
(  107)       !                !  umo(4) : 
(  108)       !                !  umo(5) : lateral diffusion
(  109)       !                !  umo(6) : vertical advection
(  110)       !                !  umo(7) : vertical diffusion
(  111)       !                !  umo(8) : surface pressure gradient
(  112)       !                !  umo(9) : 
(  113) 
(  114)    REAL(wp), DIMENSION(10) ::   &  !:
(  115)       hke              !: momentum square trends average 
(  116)       !                !  hke(1) : horizontal advection
(  117)       !                !  hke(2) : vertical advection
(  118) 
(  119)    REAL(wp) ::   &  !:
(  120)       rpktrd,    &  !: potential to kinetic energy conversion
(  121)       peke          !: conversion potential energy - kinetic energy trend
(  122) 
(  123) #endif
(  124) 
(  125)   !!======================================================================
(  126) END MODULE trdicp_oce










































PGF90 (Version     10.5)          09/23/2010  09:49:08      page 4

(    1) # 1 "trdicp_oce.F90"
(    1) MODULE trdicp_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdicp_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdicp_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12)    !!----------------------------------------------------------------------
(   13)    !!   'key_trdtra'   or                         tracer trends diagnostics
(   14)    !!   'key_trddyn'                            momentum trends diagnostics
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce                 ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PUBLIC
(   21) 
(   22)    !! Namelist parameters
(   23)    !!----------------------------------------------------------------------
(   24)    INTEGER  ::      & !!: namdia :  diagnostics on dynamics and/or tracer trends
(   25)       ntrd  = 10 ,  &  !: time step frequency dynamics and tracers trends
(   26)       nctls =  0       !: control surface type for trends vertical integration
(   27) 
(   28)    !! Tracers trends diagnostics parameters
(   29)    !!---------------------------------------------------------------------
(   30)    INTEGER, PARAMETER ::            &  !: trends index
(   31)       jpttdlad = 1,   &  !: tracer horizontal advection
(   32)       jpttdzad = 2,   &  !: tracer vertical advection
(   33)       jpttdldf = 3,   &  !: tracer horizontal diffusion
(   34)       jpttdzdf = 4,   &  !: tracer vertical diffusion
(   35)       jpttdnpc = 5,   &  !: tracer non penetrative convection
(   36)       jpttddoe = 6,   &  !: tracer D.amping O.r vertical E.iv
(   37)       jpttdqsr = 7,   &  !: tracer penetrative solar radiation
(   38)       jpttdnsr = 8       !: tracer non solar radiation
(   39) 
(   40)    !! Momentum trends diagnostics parameters
(   41)    !!---------------------------------------------------------------------
(   42)    INTEGER, PARAMETER ::            &  !: trends index
(   43)       jpdtdhpg =  1,   &  !: dynamic hydrostatic pressure gradient 
(   44)       jpdtdkeg =  2,   &  !: dynamic kinetic energy gradient
(   45)       jpdtdrvo =  3,   &  !: dynamic relative vorticity
(   46)       jpdtdpvo =  4,   &  !: dynamic planetary vorticity
(   47)       jpdtdldf =  5,   &  !: dynamic lateral diffusion
(   48)       jpdtdzad =  6,   &  !: dynamic vertical advection
(   49)       jpdtdzdf =  7,   &  !: dynamic vertical diffusion
(   50)       jpdtdspg =  8,   &  !: dynamic surface pressure gradient
(   51)       jpdtddat =  9,   &  !: dynamic damping term
(   52)       jpdtdswf = 10,   &  !: dynamic surface wind forcing
(   53)       jpdtdbfr = 11       !: dynamic bottom friction 
(   54) 
(   55)    REAL, DIMENSION(jpi,jpj) ::   &  !:
(   56)       tldfbbl, sldfbbl,          &  ! Temperature/salinity lateral diffusion trends
(   57)       !                             ! in the BBL  






PGF90 (Version     10.5)          09/23/2010  09:49:08      page 5

(   58)       tladbbl, sladbbl              ! Temperature/salinity lateral advection trends 
(   59)       !                             ! in the BBL
(   60) 
(   61)    REAL, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   62)       tladi, sladi,                  &  ! Temp./sal. MUSCL OR TVD advection fluxes 
(   63)       !                                 ! terms along i- 
(   64)       tladj, sladj                      ! Temp./sal. MUSCL OR TVD advection fluxes 
(   65)       !                                 ! terms along j- 
(   66) # 67
(   67)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  !:
(   68)       uldftrd, vldftrd     !: lateral diffusion trend in isopycnal case
(   69) 
(   70) # 125
(  125)   !!======================================================================
(  126) END MODULE trdicp_oce
