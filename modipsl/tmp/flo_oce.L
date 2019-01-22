


PGF90 (Version     10.5)          09/23/2010  09:49:46      page 1

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

Filename: flo_oce.F90

(    1) MODULE flo_oce
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE flo_oce  ***
(    4)    !!                
(    5)    !! ** Purpose : - Define in memory all floats parameters and variables
(    6)    !!
(    7)    !! History :
(    8)    !!   8.0  !  99-10  (CLIPPER projet)
(    9)    !!   9.0  !  02-11  (G. Madec, A. Bozec)  F90: Free form and module
(   10)    !!======================================================================
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flo_oce.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15) #if   defined key_floats   ||   defined key_esopa
(   16)    !!----------------------------------------------------------------------
(   17)    !!   'key_floats'                                        drifting floats
(   18)    !!----------------------------------------------------------------------
(   19)    !! * Modules used
(   20)    USE par_oce         ! ocean parameters
(   21) 
(   22)    IMPLICIT NONE
(   23) 
(   24)    LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .TRUE.    !: float flag
(   25) 
(   26) !!DB
(   27)    INTEGER, PUBLIC :: jpnfl
(   28) 
(   29) 
(   30)    !! float parameters
(   31)    !! ----------------
(   32)    INTEGER, PARAMETER ::   &
(   33) !!DB -- jpnfl is now determined from the init_floats file
(   34) !      jpnfl     = 3 ,            &  ! total number of floats during the run
(   35)       jpnnewflo =  0 ,            &  ! number of floats added in a new run
(   36) !      jpnrstflo = jpnfl-jpnnewflo    ! number of floats for the restart
(   37)       jpnrstflo = 0    ! number of floats for the restart
(   38) 
(   39)    !! float variables
(   40)    !! ---------------
(   41) !!DB
(   42)    INTEGER, DIMENSION(:),ALLOCATABLE  ::    &
(   43)       nisobfl,    &  ! 0 for a isobar float
(   44)       !              ! 1 for a float following the w velocity
(   45)       ngrpfl         ! number to identify searcher group
(   46)    REAL(wp), DIMENSION(:),ALLOCATABLE ::    &






PGF90 (Version     10.5)          09/23/2010  09:49:46      page 2

(   47)       flxx,       &  ! longitude of float (decimal degree)
(   48)       flyy,       &  ! latitude of float (decimal degree)
(   49)       flzz,       &  ! depth of float (m, positive)
(   50)       tpifl,      &  ! index of float position on zonal axe
(   51)       tpjfl,      &  ! index of float position on meridien axe
(   52)       tpkfl          ! index of float position on z axe
(   53)    
(   54)    REAL(wp), DIMENSION(jpi, jpj, jpk) ::    & 
(   55)       wb             ! vertical velocity at previous time step (m s-1).
(   56)    
(   57)    ! floats unit
(   58)    
(   59)    LOGICAL  ::                & !!! * namelist namflo *
(   60)       ln_rstflo = .FALSE. ,   &  ! T/F float restart 
(   61)       ln_argo   = .FALSE. ,   &  ! T/F argo type floats
(   62)       ln_flork4 = .FALSE.        ! T/F 4th order Runge-Kutta
(   63)    INTEGER  ::               & !!! * namelist namflo *
(   64)       nwritefl,              &  ! frequency of float output file 
(   65)       nstockfl                  ! frequency of float restart file
(   66) 
(   67) #else
(   68)    !!----------------------------------------------------------------------
(   69)    !!   Default option :                                 NO drifting floats
(   70)    !!----------------------------------------------------------------------
(   71)    LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .FALSE.   !: float flag
(   72) #endif
(   73) 
(   74)    !!======================================================================
(   75) END MODULE flo_oce



































PGF90 (Version     10.5)          09/23/2010  09:49:46      page 3

(    1) # 1 "flo_oce.F90"
(    1) MODULE flo_oce
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE flo_oce  ***
(    4)    !!                
(    5)    !! ** Purpose : - Define in memory all floats parameters and variables
(    6)    !!
(    7)    !! History :
(    8)    !!   8.0  !  99-10  (CLIPPER projet)
(    9)    !!   9.0  !  02-11  (G. Madec, A. Bozec)  F90: Free form and module
(   10)    !!======================================================================
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flo_oce.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15) # 68
(   68)    !!----------------------------------------------------------------------
(   69)    !!   Default option :                                 NO drifting floats
(   70)    !!----------------------------------------------------------------------
(   71)    LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .FALSE.   !: float flag
(   72) 
(   73) # 74
(   74)    !!======================================================================
(   75) END MODULE flo_oce
