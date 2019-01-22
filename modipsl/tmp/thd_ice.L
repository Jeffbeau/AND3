


PGF90 (Version     10.5)          09/23/2010  09:49:55      page 1

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

Filename: thd_ice.F90

(    1) MODULE thd_ice
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE thd_ice  ***
(    4)    !! LIM sea-ice :   Ice thermodynamics in 1D
(    5)    !!=====================================================================
(    6)    !! History :
(    7)    !!   2.0  !  02-11  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/thd_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_ice
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Share Module variables
(   20)    REAL(wp) , PUBLIC ::   & !!! ** ice-thermo namelist (namicethd) **
(   21)       hmelt   = -0.15  ,  &  !: maximum melting at the bottom
(   22)       hicmin  = 0.2    ,  &  !: ice th. corr. to max. ener. in brine pocket
(   23)       hiclim  = 0.05   ,  &  !: minimum ice thickness
(   24)       amax    = 0.999  ,  &  !: maximum lead fraction
(   25)       swiqst  = 1.0    ,  &  !: energy stored in brine pocket (1) or not (0)
(   26)       sbeta   = 1.0    ,  &  !: numerical scheme for diffusion in ice 
(   27)       parlat  = 0.0    ,  &  !: percent. of energy used for lateral ablation
(   28)       hakspl  = 0.5    ,  &  !: slope of distr. for Hakkinen-Mellro's lat. melt
(   29)       hibspl  = 0.5    ,  &  !: slope of distribution for Hibler's lat. melt
(   30)       exld    = 2.0    ,  &  !: exponent for leads-closure rate
(   31)       hakdif  = 1.0    ,  &  !: coefficient for diffusions of ice and snow
(   32)       thth    = 0.2    ,  &  !: thick. for comp. of eq. thermal conduct
(   33)       hnzst   = 0.1    ,  &  !: thick. of the surf. layer in temp. comp.
(   34)       parsub  = 1.0    ,  &  !: switch for snow sublimation or not
(   35)       alphs   = 1.0          !: coef. for snow density when snow-ice formation
(   36) 
(   37)    REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:   
(   38)       hiccrit = (/0.3,0.3/)  !: ice th. for lateral accretion in the NH (SH) (m)
(   39) 
(   40)    REAL(wp) , PUBLIC ::   &  !:
(   41)       uscomi,             &  !: inverse of minimum lead fraction
(   42)       cnscg                  !: ratio  rcpsn/rcpic
(   43) 
(   44)    INTEGER , PUBLIC, DIMENSION(jpij) ::   &  !:
(   45)       npb     ,   &   !: number of points where computations has to be done
(   46)       npac            !: correspondance between the points






PGF90 (Version     10.5)          09/23/2010  09:49:55      page 2

(   47) 
(   48)    REAL(wp), PUBLIC, DIMENSION(jpij) ::   &  !: 
(   49)       qldif_1d    ,     &  !: corresponding to the 2D var  qldif
(   50)       qcmif_1d    ,     &  !: corresponding to the 2D var  qcmif
(   51)       thcm_1d     ,     &  !:    "                  "      thcm
(   52)       fstbif_1d   ,     &  !:    "                  "      fstric
(   53)       fltbif_1d   ,     &  !:    "                  "      ffltbif
(   54)       fscbq_1d    ,     &  !:    "                  "      fscmcbq
(   55)       qsr_ice_1d  ,     &  !:    "                  "      qsr_ice
(   56)       fr1_i0_1d   ,     &  !:    "                  "      fr1_i0
(   57)       fr2_i0_1d   ,     &  !:    "                  "      fr2_i0
(   58)       qnsr_ice_1d ,     &  !:    "                  "      qns_ice
(   59)       qfvbq_1d    ,     &  !:    "                  "      qfvbq
(   60)       sist_1d     ,     &  !:    "                  "      sist
(   61)       tfu_1d      ,     &  !:    "                  "      tfu
(   62)       sprecip_1d  ,     &  !:    "                  "      sprecip
(   63)       h_snow_1d   ,     &  !:    "                  "      h_snow
(   64)       h_ice_1d    ,     &  !:    "                  "      h_ice
(   65)       frld_1d     ,     &  !:    "                  "      frld
(   66)       qstbif_1d   ,     &  !:    "                  "      qstoif
(   67)       fbif_1d     ,     &  !:    "                  "      fbif
(   68)       rdmicif_1d  ,     &  !:    "                  "      rdmicif
(   69)       rdmsnif_1d  ,     &  !:    "                  "      rdmsnif
(   70)       qlbbq_1d    ,     &  !:    "                  "      qlbsbq
(   71)       dmgwi_1d    ,     &  !:    "                  "      dmgwi
(   72)       dvsbq_1d    ,     &  !:    "                  "      rdvosif
(   73)       dvbbq_1d    ,     &  !:    "                  "      rdvobif
(   74)       dvlbq_1d    ,     &  !:    "                  "      rdvolif
(   75)       dvnbq_1d    ,     &  !:    "                  "      rdvolif
(   76)       dqns_ice_1d ,     &  !:    "                  "      dqns_ice
(   77)       qla_ice_1d  ,     &  !:    "                  "      qla_ice
(   78)       dqla_ice_1d          !:    "                  "      dqla_ice
(   79) 
(   80)    REAL(wp), PUBLIC, DIMENSION(jpij,jplayersp1) ::   &  !:
(   81)       tbif_1d              !: corresponding to the 2D var  tbif
(   82) 
(   83)    !!======================================================================
(   84) END MODULE thd_ice


























PGF90 (Version     10.5)          09/23/2010  09:49:55      page 3

(    1) # 1 "thd_ice.F90"
(    1) MODULE thd_ice
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE thd_ice  ***
(    4)    !! LIM sea-ice :   Ice thermodynamics in 1D
(    5)    !!=====================================================================
(    6)    !! History :
(    7)    !!   2.0  !  02-11  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/thd_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_ice
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Share Module variables
(   20)    REAL(wp) , PUBLIC ::   & !!! ** ice-thermo namelist (namicethd) **
(   21)       hmelt   = -0.15  ,  &  !: maximum melting at the bottom
(   22)       hicmin  = 0.2    ,  &  !: ice th. corr. to max. ener. in brine pocket
(   23)       hiclim  = 0.05   ,  &  !: minimum ice thickness
(   24)       amax    = 0.999  ,  &  !: maximum lead fraction
(   25)       swiqst  = 1.0    ,  &  !: energy stored in brine pocket (1) or not (0)
(   26)       sbeta   = 1.0    ,  &  !: numerical scheme for diffusion in ice 
(   27)       parlat  = 0.0    ,  &  !: percent. of energy used for lateral ablation
(   28)       hakspl  = 0.5    ,  &  !: slope of distr. for Hakkinen-Mellro's lat. melt
(   29)       hibspl  = 0.5    ,  &  !: slope of distribution for Hibler's lat. melt
(   30)       exld    = 2.0    ,  &  !: exponent for leads-closure rate
(   31)       hakdif  = 1.0    ,  &  !: coefficient for diffusions of ice and snow
(   32)       thth    = 0.2    ,  &  !: thick. for comp. of eq. thermal conduct
(   33)       hnzst   = 0.1    ,  &  !: thick. of the surf. layer in temp. comp.
(   34)       parsub  = 1.0    ,  &  !: switch for snow sublimation or not
(   35)       alphs   = 1.0          !: coef. for snow density when snow-ice formation
(   36) 
(   37)    REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:   
(   38)       hiccrit = (/0.3,0.3/)  !: ice th. for lateral accretion in the NH (SH) (m)
(   39) 
(   40)    REAL(wp) , PUBLIC ::   &  !:
(   41)       uscomi,             &  !: inverse of minimum lead fraction
(   42)       cnscg                  !: ratio  rcpsn/rcpic
(   43) 
(   44)    INTEGER , PUBLIC, DIMENSION(jpij) ::   &  !:
(   45)       npb     ,   &   !: number of points where computations has to be done
(   46)       npac            !: correspondance between the points
(   47) 
(   48)    REAL(wp), PUBLIC, DIMENSION(jpij) ::   &  !: 
(   49)       qldif_1d    ,     &  !: corresponding to the 2D var  qldif
(   50)       qcmif_1d    ,     &  !: corresponding to the 2D var  qcmif
(   51)       thcm_1d     ,     &  !:    "                  "      thcm
(   52)       fstbif_1d   ,     &  !:    "                  "      fstric
(   53)       fltbif_1d   ,     &  !:    "                  "      ffltbif
(   54)       fscbq_1d    ,     &  !:    "                  "      fscmcbq
(   55)       qsr_ice_1d  ,     &  !:    "                  "      qsr_ice
(   56)       fr1_i0_1d   ,     &  !:    "                  "      fr1_i0
(   57)       fr2_i0_1d   ,     &  !:    "                  "      fr2_i0






PGF90 (Version     10.5)          09/23/2010  09:49:55      page 4

(   58)       qnsr_ice_1d ,     &  !:    "                  "      qns_ice
(   59)       qfvbq_1d    ,     &  !:    "                  "      qfvbq
(   60)       sist_1d     ,     &  !:    "                  "      sist
(   61)       tfu_1d      ,     &  !:    "                  "      tfu
(   62)       sprecip_1d  ,     &  !:    "                  "      sprecip
(   63)       h_snow_1d   ,     &  !:    "                  "      h_snow
(   64)       h_ice_1d    ,     &  !:    "                  "      h_ice
(   65)       frld_1d     ,     &  !:    "                  "      frld
(   66)       qstbif_1d   ,     &  !:    "                  "      qstoif
(   67)       fbif_1d     ,     &  !:    "                  "      fbif
(   68)       rdmicif_1d  ,     &  !:    "                  "      rdmicif
(   69)       rdmsnif_1d  ,     &  !:    "                  "      rdmsnif
(   70)       qlbbq_1d    ,     &  !:    "                  "      qlbsbq
(   71)       dmgwi_1d    ,     &  !:    "                  "      dmgwi
(   72)       dvsbq_1d    ,     &  !:    "                  "      rdvosif
(   73)       dvbbq_1d    ,     &  !:    "                  "      rdvobif
(   74)       dvlbq_1d    ,     &  !:    "                  "      rdvolif
(   75)       dvnbq_1d    ,     &  !:    "                  "      rdvolif
(   76)       dqns_ice_1d ,     &  !:    "                  "      dqns_ice
(   77)       qla_ice_1d  ,     &  !:    "                  "      qla_ice
(   78)       dqla_ice_1d          !:    "                  "      dqla_ice
(   79) 
(   80)    REAL(wp), PUBLIC, DIMENSION(jpij,jplayersp1) ::   &  !:
(   81)       tbif_1d              !: corresponding to the 2D var  tbif
(   82) 
(   83)    !!======================================================================
(   84) END MODULE thd_ice
