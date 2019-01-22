MODULE thd_ice
   !!======================================================================
   !!                       ***  MODULE thd_ice  ***
   !! LIM sea-ice :   Ice thermodynamics in 1D
   !!=====================================================================
   !! History :
   !!   2.0  !  02-11  (C. Ethe)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/thd_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_ice

   IMPLICIT NONE
   PRIVATE

   !! * Share Module variables
   REAL(wp) , PUBLIC ::   & !!! ** ice-thermo namelist (namicethd) **
      hmelt   = -0.15  ,  &  !: maximum melting at the bottom
      hicmin  = 0.2    ,  &  !: ice th. corr. to max. ener. in brine pocket
      hiclim  = 0.05   ,  &  !: minimum ice thickness
      amax    = 0.999  ,  &  !: maximum lead fraction
      swiqst  = 1.0    ,  &  !: energy stored in brine pocket (1) or not (0)
      sbeta   = 1.0    ,  &  !: numerical scheme for diffusion in ice 
      parlat  = 0.0    ,  &  !: percent. of energy used for lateral ablation
      hakspl  = 0.5    ,  &  !: slope of distr. for Hakkinen-Mellro's lat. melt
      hibspl  = 0.5    ,  &  !: slope of distribution for Hibler's lat. melt
      exld    = 2.0    ,  &  !: exponent for leads-closure rate
      hakdif  = 1.0    ,  &  !: coefficient for diffusions of ice and snow
      thth    = 0.2    ,  &  !: thick. for comp. of eq. thermal conduct
      hnzst   = 0.1    ,  &  !: thick. of the surf. layer in temp. comp.
      parsub  = 1.0    ,  &  !: switch for snow sublimation or not
      alphs   = 1.0          !: coef. for snow density when snow-ice formation

   REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:   
      hiccrit = (/0.3,0.3/)  !: ice th. for lateral accretion in the NH (SH) (m)

   REAL(wp) , PUBLIC ::   &  !:
      uscomi,             &  !: inverse of minimum lead fraction
      cnscg                  !: ratio  rcpsn/rcpic

   INTEGER , PUBLIC, DIMENSION(jpij) ::   &  !:
      npb     ,   &   !: number of points where computations has to be done
      npac            !: correspondance between the points

   REAL(wp), PUBLIC, DIMENSION(jpij) ::   &  !: 
      qldif_1d    ,     &  !: corresponding to the 2D var  qldif
      qcmif_1d    ,     &  !: corresponding to the 2D var  qcmif
      thcm_1d     ,     &  !:    "                  "      thcm
      fstbif_1d   ,     &  !:    "                  "      fstric
      fltbif_1d   ,     &  !:    "                  "      ffltbif
      fscbq_1d    ,     &  !:    "                  "      fscmcbq
      qsr_ice_1d  ,     &  !:    "                  "      qsr_ice
      fr1_i0_1d   ,     &  !:    "                  "      fr1_i0
      fr2_i0_1d   ,     &  !:    "                  "      fr2_i0
      qnsr_ice_1d ,     &  !:    "                  "      qns_ice
      qfvbq_1d    ,     &  !:    "                  "      qfvbq
      sist_1d     ,     &  !:    "                  "      sist
      tfu_1d      ,     &  !:    "                  "      tfu
      sprecip_1d  ,     &  !:    "                  "      sprecip
      h_snow_1d   ,     &  !:    "                  "      h_snow
      h_ice_1d    ,     &  !:    "                  "      h_ice
      frld_1d     ,     &  !:    "                  "      frld
      qstbif_1d   ,     &  !:    "                  "      qstoif
      fbif_1d     ,     &  !:    "                  "      fbif
      rdmicif_1d  ,     &  !:    "                  "      rdmicif
      rdmsnif_1d  ,     &  !:    "                  "      rdmsnif
      qlbbq_1d    ,     &  !:    "                  "      qlbsbq
      dmgwi_1d    ,     &  !:    "                  "      dmgwi
      dvsbq_1d    ,     &  !:    "                  "      rdvosif
      dvbbq_1d    ,     &  !:    "                  "      rdvobif
      dvlbq_1d    ,     &  !:    "                  "      rdvolif
      dvnbq_1d    ,     &  !:    "                  "      rdvolif
      dqns_ice_1d ,     &  !:    "                  "      dqns_ice
      qla_ice_1d  ,     &  !:    "                  "      qla_ice
      dqla_ice_1d          !:    "                  "      dqla_ice

   REAL(wp), PUBLIC, DIMENSION(jpij,jplayersp1) ::   &  !:
      tbif_1d              !: corresponding to the 2D var  tbif

   !!======================================================================
END MODULE thd_ice
