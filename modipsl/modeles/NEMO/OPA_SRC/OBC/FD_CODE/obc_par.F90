MODULE obc_par
   !!==============================================================================
   !!                  ***  MODULE obc_par   ***
   !! Open Boundary Cond. :   define related parameters
   !!==============================================================================
#if defined key_obc
   !!----------------------------------------------------------------------
   !!   'key_obc' :                                 Open Boundary Condition
   !!----------------------------------------------------------------------
   !! history :
   !!  8.0   01/91   (CLIPPER)  Original code 
   !!  9.0   06/02   (C. Talandier)  modules
   !!        06/04   (F. Durand) ORCA_R2_ZIND config
   !!        06/04   (F. Durand) jptobc is defined as a parameter, 
   !!	     	     in order to allow time-dependent OBCs fields on input
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters

   IMPLICIT NONE
   PUBLIC
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_par.F90,v 1.5 2005/12/12 14:20:26 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_agrif
   LOGICAL, PUBLIC ::     &  !:
#else
   LOGICAL, PUBLIC, PARAMETER :: &  !:
#endif
      lk_obc = .TRUE.   !: Ocean Boundary Condition flag

# if defined key_eel_r5
   !!----------------------------------------------------------------------
   !!   'key_eel_r5' :                                 EEL R5 configuration
   !!----------------------------------------------------------------------
#    include "obc_par_EEL_R5.h90"

# elif defined key_na_bio
   !!---------------------------------------------------------------------
   !!   'key_na_bio'      :                    North Atlantic 0.25 degree
   !!---------------------------------------------------------------------
#    include "obc_par_NA025.h90"

# else
   !!---------------------------------------------------------------------
   !! open boundary parameter
   !!---------------------------------------------------------------------
   INTEGER, PARAMETER ::     &  !: time dimension of the BCS fields on input
      jptobc  =	      2 

   !! * EAST open boundary
#if defined key_agrif
   LOGICAL ::     &  !:
#else
   LOGICAL, PARAMETER ::     &  !:
#endif
      lp_obc_east = .TRUE.     !: to active or not the East open boundary
   INTEGER  &  !:
#if !defined key_agrif
      ,PARAMETER &
#endif
    :: &
      jpieob  = jpiglo-2,    &  !: i-localization of the East open boundary (must be ocean U-point)
      jpjed   =        2,    &  !: j-starting indice of the East open boundary (must be land T-point)
      jpjef   = jpjglo-1,    &  !: j-ending   indice of the East open boundary (must be land T-point)
      jpjedp1 =  jpjed+1,    &  !: first ocean point         "                 "
      jpjefm1 =  jpjef-1        !: last  ocean point         "                 "

   !! * WEST open boundary
#if defined key_agrif
   LOGICAL ::     &  !:
#else
   LOGICAL, PARAMETER ::     &  !:
#endif
      lp_obc_west = .FALSE.     !: to active or not the West open boundary
   INTEGER  &  !:
#if !defined key_agrif
      ,PARAMETER &
#endif
    :: &
      jpiwob  =	       2,    &  !: i-localization of the West open boundary (must be ocean U-point)
      jpjwd   =	       2,    &  !: j-starting indice of the West open boundary (must be land T-point)
      jpjwf   = jpjglo-1,    &  !: j-ending   indice of the West open boundary (must be land T-point)
      jpjwdp1 =  jpjwd+1,    &  !: first ocean point         "                 "
      jpjwfm1 =  jpjwf-1        !: last  ocean point         "                 "

   !! * NORTH open boundary
#if defined key_agrif
   LOGICAL ::     &  !:
#else
   LOGICAL, PARAMETER ::     &  !:
#endif
      lp_obc_north = .TRUE.    !: to active or not the North open boundary
   INTEGER  &  !:
#if !defined key_agrif
      ,PARAMETER &
#endif
    :: &
      jpjnob  = jpjglo-2,    &  !: j-localization of the North open boundary (must be ocean V-point)
      jpind   =        2,    &  !: i-starting indice of the North open boundary (must be land T-point)
      jpinf   = jpiglo-1,    &  !: i-ending   indice of the North open boundary (must be land T-point)
      jpindp1 =  jpind+1,    &  !: first ocean point         "                 "
      jpinfm1 =  jpinf-1        !: last  ocean point         "                 "

   !! * SOUTH open boundary
!dbg DW open the south boundary
#if defined key_agrif
   LOGICAL ::     &  !:
#else
   LOGICAL, PARAMETER ::     &  !:
#endif
      lp_obc_south = .TRUE.    !: to active or not the South open boundary
   INTEGER  &  !:
#if !defined key_agrif
      ,PARAMETER &
#endif
    :: &
      jpjsob  =        2,    &  !: j-localization of the South open boundary (must be ocean V-point)
      jpisd   =        2,    &  !: i-starting indice of the South open boundary (must be land T-point)
      jpisf   = jpiglo-1,    &  !: i-ending   indice of the South open boundary (must be land T-point)
      jpisdp1 =  jpisd+1,    &  !: first ocean point         "                 "
      jpisfm1 =  jpisf-1        !: last  ocean point         "                 "
   
   INTEGER, PARAMETER ::     &  !:
      jpnic = 2700              !: maximum number of isolated coastlines points 

# endif

#else
   !!----------------------------------------------------------------------
   !!   Default option :                         NO open boundary condition
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_obc = .FALSE.  !: Ocean Boundary Condition flag
#endif

   !!======================================================================
END MODULE obc_par
