MODULE trc
   !!======================================================================
   !!                      ***  MODULE  trc  ***
   !! Passive tracers   :  module for tracers defined
   !!======================================================================
   !! History :
   !!   8.2  !  96-01  (M. Levy)  Original code
   !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
   !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
   !!   9.0  !  04-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trc.F90,v 1.9 2006/04/11 13:48:58 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   'key_passivetrc'   :                               Passive tracer
   !!---------------------------------------------------------------------
   !! * Modules used
   USE par_oce
   USE par_trc
   IMPLICIT NONE

   PUBLIC


   !! passive tracers names and units (read in namelist)
   !! --------------------------------------------------
   CHARACTER(len=12), PUBLIC, DIMENSION(jptra) :: &
      ctrcnm  ,   &   !!: tracer name 
      ctrcun          !!: tracer unit

   CHARACTER(len=80), PUBLIC, DIMENSION(jptra) :: &
      ctrcnl          !!: tracer long name 
   
   
   !! parameters for the control of passive tracers
   !! --------------------------------------------------
   INTEGER, PUBLIC ::  &
      numnat          !!: the number of the passive tracer NAMELIST
   
   LOGICAL, PUBLIC, DIMENSION(jptra) ::   &
      lutini          !!:  initialisation from FILE or not (NAMELIST)

   INTEGER , PUBLIC, DIMENSION(jptra) :: &
      nutini          !!: FORTRAN LOGICAL UNIT for initialisation file

   !! passive tracers fields (before,now,after)
   !! --------------------------------------------------
   REAL(wp), PUBLIC, SAVE  ::  &
      trai    ,   &   !!: initial total tracer
      areatot         !!: total volume 

   REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk,jptra) :: &
      trn     ,   &   !!: traceur concentration for actual time step
      tra     ,   &   !!: traceur concentration for next time step
      trb             !!: traceur concentration for before time step


   !! numerical parameter (NAMELIST)
   !! --------------------------------------------------
   REAL(wp), PUBLIC  ::  &
      rsc     ,   &   !!: tuning coefficient for anti-diffusion
      rtrn            !!: value for truncation

   !! namelist parameters
   !! --------------------------------------------------
   INTEGER , PUBLIC  ::  & 
      ncortrc ,   &   !!: number of corrective phases
      ndttrc  ,   &   !!: frequency of step on passive tracers
      nittrc000       !!: first time step of passive tracers model  

   LOGICAL, PUBLIC  ::  & 
      crosster        !!: logical if true computes crossterms


   !! isopycnal scheme for passive tracers
   !! --------------------------------------------------  
   REAL(wp), PUBLIC  ::  &
      ahtrb0  ,   &   !!: background diffusivity coefficient for passive tracer (m2/s)
      trcrat  ,   &   !!: ratio between passive and active tracer coeff for diffusion
      ahtrc0  ,   &   !!: horizontal eddy diffusivity for passive tracers (m2/s)
      aeivtr0         !!: eddy induced velocity coefficient (m2/s)
   
   
   !! passive tracers restart (input and output)
   !! --------------------------------------------------  
   LOGICAL, PUBLIC  ::  &
      lrsttr          !!: boolean term for restart i/o for passive tracers (namelist)
   
   INTEGER , PUBLIC  ::  &
      nutwrs  ,   &   !!: output FILE for passive tracers restart
      nutrst  ,   &   !!: logical unit for restart FILE for passive tracers
      nrsttr          !!: control of the time step ( 0 or 1 ) for pass. tr.
   
#if defined key_partial_steps
   
   !! interpolated gradient
   !!--------------------------------------------------  
   REAL (wp), PUBLIC, DIMENSION (jpi,jpj,jptra) :: &
      gtru    ,   &   !!: horizontal gradient at u-points at bottom ocean level
      gtrv            !!: horizontal gradient at v-points at bottom ocean level
#else
   REAL (wp), PUBLIC :: &
      gtru    ,   &   !!: horizontal gradient at u-points at bottom ocean level
      gtrv            !!: horizontal gradient at v-points at bottom ocean level
   
#endif
   
#if defined key_trcldf_eiv && defined key_diaeiv
   !! The three component of the eddy induced velocity
   !! --------------------------------------------------
   REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk) :: &
      u_trc_eiv,  &   !!: u-eiv (m/s)
      v_trc_eiv,  &   !!: v-eiv (m/s)
      w_trc_eiv       !!: w-eiv (m/s)
#endif
   
   
   !! information for outputs
   !! --------------------------------------------------
   INTEGER , PUBLIC   ::  & 
      nwritetrc       !!: time step frequency for concentration outputs (namelist)
   
#if defined key_trc_diaadd
   !! additional 2D/3D outputs namelist
   !! --------------------------------------------------
   CHARACTER(len=8), PUBLIC, DIMENSION (jpdia2d) ::  & 
      ctrc2d  ,   &   !!: 2d output field name
      ctrc2u          !!: 2d output field unit
   
   CHARACTER(len=8), PUBLIC, DIMENSION (jpdia3d) ::  & 
      ctrc3d ,    &   !!: 3d output field name
      ctrc3u          !!: 3d output field unit
   
   CHARACTER(len=80), PUBLIC, DIMENSION (jpdia2d) ::  & 
      ctrc2l          !!: 2d output field long name
   
   CHARACTER(len=80), PUBLIC, DIMENSION (jpdia3d) ::  & 
      ctrc3l          !!: 3d output field long name
   
   REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpdia2d) ::  &  
      trc2d           !!:  additional 2d outputs  
   
   REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk,jpdia3d) ::  &  
      trc3d           !!:  additional 3d outputs  
   
   
   !! netcdf files and index common
   !! --------------------------------------------------
   INTEGER , PUBLIC :: &
      nwriteadd     !!: frequency of additional arrays outputs(namelist)
#endif
   
#if defined key_trc_diatrd
   
   !!  non conservative trends (biological, ...)
   !! --------------------------------------------------
   LOGICAL, PUBLIC, DIMENSION (jptra)  ::  &  
      luttrd          !!: large trends diagnostic to write or not (namelist)
   
   !!  dynamical trends
   !!	trtrd()   : trends of the tracer equations
   !!           1 : X advection
   !!           2 : Y advection
   !!           3 : Z advection
   !!           4 : X diffusion
   !!           5 : Y diffusion
   !!           6 : Z diffusion
   !!           7 : X gent velocity
   !!           8 : Y gent velocity
   !!           9 : Z gent velocity
   !! --------------------------------------------------
   
   
   REAL(wp), PUBLIC, DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE :: &
      trtrd           !!: trends of the tracer equations
   
   INTEGER, PUBLIC, DIMENSION(jptra), SAVE :: ikeep ! indice of tracer for which dyn trends are stored
   INTEGER, PUBLIC, SAVE                   :: nkeep ! number of tracers for which dyn trends are stored 
                                                    ! (used to allocate trtrd buffer)

   !! netcdf files and index common
   !! --------------------------------------------------
   INTEGER , PUBLIC :: &
      nwritetrd       !!: frequency of additional arrays outputs(namelist)
   
#endif 
   
   !! passive tracers data read and at given time_step
   !! --------------------------------------------------
#if defined key_trc_dta
   
   INTEGER , PUBLIC, DIMENSION(jptra) :: &
      numtr          !!: logical unit for passive tracers data
   
#endif

  !!  1D configuration
  !! --------------------------------------------------
#if defined key_cfg_1d
      LOGICAL, PARAMETER ::   lk_trccfg_1d   = .TRUE.   !: 1D pass. tracer configuration flag
#else   
      LOGICAL, PARAMETER ::   lk_trccfg_1d   = .FALSE.  !: 1D pass. tracer configuration flag
#endif


#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE trc
