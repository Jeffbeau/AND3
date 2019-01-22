MODULE flo_oce
   !!======================================================================
   !!                     ***  MODULE flo_oce  ***
   !!                
   !! ** Purpose : - Define in memory all floats parameters and variables
   !!
   !! History :
   !!   8.0  !  99-10  (CLIPPER projet)
   !!   9.0  !  02-11  (G. Madec, A. Bozec)  F90: Free form and module
   !!======================================================================
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flo_oce.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                        drifting floats
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters

   IMPLICIT NONE

   LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .TRUE.    !: float flag

!!DB
   INTEGER, PUBLIC :: jpnfl


   !! float parameters
   !! ----------------
   INTEGER, PARAMETER ::   &
!!DB -- jpnfl is now determined from the init_floats file
!      jpnfl     = 3 ,            &  ! total number of floats during the run
      jpnnewflo =  0 ,            &  ! number of floats added in a new run
!      jpnrstflo = jpnfl-jpnnewflo    ! number of floats for the restart
      jpnrstflo = 0    ! number of floats for the restart

   !! float variables
   !! ---------------
!!DB
   INTEGER, DIMENSION(:),ALLOCATABLE  ::    &
      nisobfl,    &  ! 0 for a isobar float
      !              ! 1 for a float following the w velocity
      ngrpfl         ! number to identify searcher group
   REAL(wp), DIMENSION(:),ALLOCATABLE ::    &
      flxx,       &  ! longitude of float (decimal degree)
      flyy,       &  ! latitude of float (decimal degree)
      flzz,       &  ! depth of float (m, positive)
      tpifl,      &  ! index of float position on zonal axe
      tpjfl,      &  ! index of float position on meridien axe
      tpkfl          ! index of float position on z axe
   
   REAL(wp), DIMENSION(jpi, jpj, jpk) ::    & 
      wb             ! vertical velocity at previous time step (m s-1).
   
   ! floats unit
   
   LOGICAL  ::                & !!! * namelist namflo *
      ln_rstflo = .FALSE. ,   &  ! T/F float restart 
      ln_argo   = .FALSE. ,   &  ! T/F argo type floats
      ln_flork4 = .FALSE.        ! T/F 4th order Runge-Kutta
   INTEGER  ::               & !!! * namelist namflo *
      nwritefl,              &  ! frequency of float output file 
      nstockfl                  ! frequency of float restart file

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                 NO drifting floats
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .FALSE.   !: float flag
#endif

   !!======================================================================
END MODULE flo_oce
