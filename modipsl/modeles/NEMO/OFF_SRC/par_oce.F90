MODULE par_oce
   !!======================================================================
   !!                        ***  par_oce  ***
   !! Ocean :   set the ocean parameters
   !!======================================================================
   !! History :
   !!   4.0  !  91     (Imbard, Levy, Madec)  Original code
   !!   9.0  !  04-01  (G. Madec, J.-M. Molines)  Free form and module
   !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_oce.F90,v 1.11 2006/03/10 10:55:34 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!   Domain decomposition
   !!----------------------------------------------------------------------
   !! * if we dont use massively parallel computer (parameters jpni=jpnj=1)
   !!      so jpiglo=jpi and jpjglo=jpj

   INTEGER, PUBLIC, PARAMETER ::    &  !:
      jpni   = 2,                   &  !: number of processors following i 
      jpnj   = 4,                   &  !: number of processors following j
      jpnij  = 8,                   &  !: nb of local domain = nb of processors 
      !                                !  ( <= jpni x jpnj )
      jpr2di = 0,                   &  !: number of columns for extra outer halo 
      jpr2dj = 0,                   &  !: number of rows    for extra outer halo 
      jpreci = 1,                   &  !: number of columns for overlap 
      jprecj = 1                       !: number of rows    for overlap 

   !! Ocean Domain sizes
   !! ------------------
   !!   data           domain   (jpidta,jpjdta)
   !!   global or zoom domain   (jpiglo,jpjglo)
   !!   local          domain   ( jpi  , jpj  )
   
#if   defined key_orca_r4
   !!---------------------------------------------------------------------
   !!   'key_orca_r4'   :                           global ocean : ORCA R4
   !!---------------------------------------------------------------------
#             include "par_ORCA_R4.h90"
#elif defined key_orca_r2
   !!---------------------------------------------------------------------
   !!   'key_orca_r2'   :                           global ocean : ORCA R4
   !!---------------------------------------------------------------------
#             include "par_ORCA_R2.h90"
#elif defined key_orca_r05
   !!---------------------------------------------------------------------
   !!   'key_orca_r05'  :                          global ocean : ORCA R05
   !!---------------------------------------------------------------------
#             include "par_ORCA_R05.h90"
#elif defined key_orca_r025
   !!---------------------------------------------------------------------
   !!   'key_orca_r025' :                         global ocean : ORCA R025
   !!---------------------------------------------------------------------
#             include "par_ORCA_R025.h90"
#elif defined key_eel_r2
   !!---------------------------------------------------------------------
   !!   'key_eel_r2'    :                                 channel : EEL R2
   !!---------------------------------------------------------------------
#             include "par_EEL_R2.h90"
#elif defined key_eel_r5
   !!---------------------------------------------------------------------
   !!   'key_eel_r5'    :                                 channel : EEL R5
   !!---------------------------------------------------------------------
#             include "par_EEL_R5.h90"
#elif defined key_eel_r6
   !!---------------------------------------------------------------------
   !!   'key_eel_r6'    :                                 channel : EEL R6
   !!---------------------------------------------------------------------
#             include "par_EEL_R6.h90"
#elif defined key_gyre
   !!---------------------------------------------------------------------
   !!   'key_gyre'      :                        mid-latitude basin : GYRE
   !!---------------------------------------------------------------------
#             include "par_GYRE.h90"
#elif defined key_ss_bio
   !!---------------------------------------------------------------------
   !!   'key_ss_bio'   :                       Scotian Shelf
   !!---------------------------------------------------------------------
#             include "par_SS_R008.h90"
   !byoung
#else
   !!---------------------------------------------------------------------
   !!   default option  :                               small closed basin
   !!---------------------------------------------------------------------
   CHARACTER(len=16), PUBLIC, PARAMETER ::   &  !:
      cp_cfg = "default"               !: name of the configuration
   INTEGER, PARAMETER ::            &  !:
      jp_cfg = 0  ,                 &  !: resolution of the configuration

      ! data size                     !!! * size of all input files *
      jpidta  = 10,                 &  !: 1st lateral dimension ( >= jpi )
      jpjdta  = 12,                 &  !: 2nd    "         "    ( >= jpj )
      jpkdta  = 31,                 &  !: number of levels      ( >= jpk )

      ! global or zoom domain size    !!! * computational domain *
      jpiglo  = jpidta,             &  !: 1st dimension of global domain --> i
      jpjglo  = jpjdta,             &  !: 2nd    "                  "    --> j
      jpk     = jpkdta,             &  !: number of vertical levels
      ! zoom starting position 
      jpizoom =   1   ,             &  !: left bottom (i,j) indices of the zoom
      jpjzoom =   1   ,             &  !: in data domain indices

      ! Domain characteristics
      jperio  =  0,                 &  !: lateral cond. type (between 0 and 6)
         !                             !  = 0 closed
         !                             !  = 1 cyclic East-West
         !                             !  = 2 equatorial symmetric
         !                             !  = 3 North fold T-point pivot
         !                             !  = 4 cyclic East-West AND North fold T-point pivot
         !                             !  = 5 North fold F-point pivot
         !                             !  = 6 cyclic East-West AND North fold F-point pivot
      jpisl   =  0,                 &  !: number of islands (rigid-lid only)
      jpnisl  =  0                     !: maximum number of points per island

      !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
      !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
      REAL(wp), PARAMETER ::   &  !:
         pp_not_used       = 999999._wp , &  !:
         pp_to_be_computed = 999999._wp      !:


   !! Horizontal grid parameters for domhgr
   !! =====================================

   INTEGER, PUBLIC, PARAMETER   ::   &  !:
      jphgr_msh = 0            !: type of horizontal mesh
      !                        !  = 0 curvilinear coordinate on the sphere
      !                        !      read in coordinate.nc file
      !                        !  = 1 geographical mesh on the sphere
      !                        !      with regular grid-spacing
      !                        !  = 2 f-plane with regular grid-spacing
      !                        !  = 3 beta-plane with regular grid-spacing
      !                        !  = 4 Mercator grid with T/U point at the equator  with
      !                        !      isotropic resolution (e1_deg)

   REAL(wp) , PUBLIC, PARAMETER ::   &   !:
      ppglam0  =    0.0_wp,   &  !: longitude of first raw and column T-point (jphgr_msh = 1)
      ppgphi0  =  -35.0_wp,   &  !: latitude  of first raw and column T-point (jphgr_msh = 1)
      !                          !  latitude for the Coriolis or Beta parameter (jphgr_msh = 2 or 3)
      ppe1_deg =    1.0_wp,   &  !: zonal      grid-spacing (degrees)
      ppe2_deg =    0.5_wp,   &  !: meridional grid-spacing (degrees)
      ppe1_m   = 5000.0_wp,   &  !: zonal      grid-spacing (degrees)
      ppe2_m   = 5000.0_wp       !: meridional grid-spacing (degrees)

   !! Vertical grid parameter for domzgr
   !! ==================================

   REAL(wp), PUBLIC, PARAMETER  ::   &  !:
      &     ppsur = -4762.96143546300_wp ,  &  !: ORCA r4, r2 and r05 coefficients
      &     ppa0  =   255.58049070440_wp ,  &  !: (default coefficients)
      &     ppa1  =   245.58132232490_wp ,  &  !:
      &     ppkth =    21.43336197938_wp ,  &  !:
      &     ppacr =     3.00000000000_wp       !:

   !!  If both ppa0 ppa1 and ppsur are specified to 0, then
   !!  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr

   REAL(wp), PUBLIC, PARAMETER ::   &  !:
      &     ppdzmin = 10._wp             ,  &  !: Minimum vertical spacing
      &     pphmax  = 5000._wp                 !: Maximum depth

   !!---------------------------------------------------------------------
#endif

   !!---------------------------------------------------------------------
   !! Domain Matrix size
   !!---------------------------------------------------------------------
   INTEGER  &  !:
#if !defined key_agrif
      ,PARAMETER  &
#endif
    :: &
      jpi = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci ,   &  !: first  dimension
      jpj = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj ,   &  !: second dimension
      jpim1 = jpi-1,                                             &  !: inner domain indices
      jpjm1 = jpj-1,                                             &  !:   "            "
      jpkm1 = jpk-1,                                             &  !:   "            "
      jpij  = jpi*jpj                                               !:  jpi x jpj

#if defined key_agrif
   !!---------------------------------------------------------------------
   !! Agrif variables
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER :: nbghostcells = 1
   INTEGER, PUBLIC :: nbcellsx = jpiglo - 2 - 2*nbghostcells
   INTEGER, PUBLIC :: nbcellsy = jpjglo - 2 - 2*nbghostcells
#endif
   !!---------------------------------------------------------------------
   !! Optimization/control flags
   !!---------------------------------------------------------------------
#if defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .TRUE.   !: flag to activate the all options
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options
#endif

#if defined key_vectopt_memory
   LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_mem  = .TRUE.   !: vector optimization flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_mem  = .FALSE.  !: vector optimization flag
#endif

#if defined key_vectopt_loop
   LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .TRUE.   !: vector optimization flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .FALSE.  !: vector optimization flag
#endif

#if defined key_autotasking
   LOGICAL, PUBLIC, PARAMETER ::   lk_jki = .TRUE.   !: j-k-i loop flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_jki = .FALSE.  !: k-j-i loop flag
#endif

   !!======================================================================
END MODULE par_oce
