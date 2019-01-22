


PGF90 (Version     10.5)          09/23/2010  09:49:03      page 1

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

Filename: ldfdyn_oce.F90

(    1) MODULE ldfdyn_oce
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE  ldfdyn_oce  ***
(    4)    !! Ocean physics:  lateral momentum mixing coefficient defined in memory 
(    5)    !!======================================================================
(    6)    !!
(    7)    !! ** Purpose :
(    8)    !!       - Define in memory lateral momentum mixing coefficients
(    9)    !!
(   10)    !! History :
(   11)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(   12)    !!----------------------------------------------------------------------
(   13)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_oce.F90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce      ! ocean parameters
(   19) 
(   20)    IMPLICIT NONE
(   21)    PUBLIC
(   22) 
(   23)    !!----------------------------------------------------------------------
(   24)    !! Lateral eddy viscosity coefficients (dynamics)
(   25)    !!----------------------------------------------------------------------
(   26) 
(   27)    LOGICAL  ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
(   28)       ln_dynldf_lap   = .TRUE.  ,   &  ! laplacian operator
(   29)       ln_dynldf_bilap = .FALSE. ,   &  ! bilaplacian operator
(   30)       ln_dynldf_level = .FALSE. ,   &  ! iso-level direction
(   31)       ln_dynldf_hor   = .TRUE.  ,   &  ! horizontal (geopotential) direction
(   32)       ln_dynldf_iso   = .FALSE.        ! iso-neutral direction
(   33) 
(   34)    REAL(wp) ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
(   35)       ahm0  = 40000._wp ,   &  ! lateral eddy viscosity (m2/s)
(   36)       ahmb0 =     0._wp        ! lateral background eddy viscosity (m2/s)
(   37) 
(   38)    LOGICAL  ::                      &  ! flag of the lateral diff. scheme used
(   39)       l_dynldf_lap         ,        &  ! iso-level laplacian operator
(   40)       l_dynldf_bilap       ,        &  ! iso-level bilaplacian operator
(   41)       l_dynldf_bilapg      ,        &  ! geopotential bilap. (s-coord)
(   42)       l_dynldf_iso         ,        &  ! iso-neutral laplacian or horizontal lapacian (s-coord)
(   43)       l_dynzdf_iso                     ! iso-neutral laplacian or horizontal lapacian (s-coord)
(   44)    
(   45) 
(   46) #if defined key_dynldf_c3d






PGF90 (Version     10.5)          09/23/2010  09:49:03      page 2

(   47)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  ! ** 3D coefficients **
(   48) #elif defined key_dynldf_smag
(   49)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  ! ** 3D coefficients **
(   50) #elif defined key_dynldf_c2d
(   51)    REAL(wp), DIMENSION(jpi,jpj)     ::   &  ! ** 2D coefficients **
(   52) #elif defined key_dynldf_c1d
(   53)    REAL(wp), DIMENSION(jpk)         ::   &  ! ** 2D coefficients **
(   54) #else
(   55)    REAL(wp)                         ::   &  ! ** 0D coefficients **
(   56) #endif
(   57)       ahm1, ahm2, ahm3, ahm4,tmph                ! ????
(   58) 
(   59)    !!----------------------------------------------------------------------
(   60) END MODULE ldfdyn_oce


















































PGF90 (Version     10.5)          09/23/2010  09:49:03      page 3

(    1) # 1 "ldfdyn_oce.F90"
(    1) MODULE ldfdyn_oce
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE  ldfdyn_oce  ***
(    4)    !! Ocean physics:  lateral momentum mixing coefficient defined in memory 
(    5)    !!======================================================================
(    6)    !!
(    7)    !! ** Purpose :
(    8)    !!       - Define in memory lateral momentum mixing coefficients
(    9)    !!
(   10)    !! History :
(   11)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(   12)    !!----------------------------------------------------------------------
(   13)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_oce.F90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce      ! ocean parameters
(   19) 
(   20)    IMPLICIT NONE
(   21)    PUBLIC
(   22) 
(   23)    !!----------------------------------------------------------------------
(   24)    !! Lateral eddy viscosity coefficients (dynamics)
(   25)    !!----------------------------------------------------------------------
(   26) 
(   27)    LOGICAL  ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
(   28)       ln_dynldf_lap   = .TRUE.  ,   &  ! laplacian operator
(   29)       ln_dynldf_bilap = .FALSE. ,   &  ! bilaplacian operator
(   30)       ln_dynldf_level = .FALSE. ,   &  ! iso-level direction
(   31)       ln_dynldf_hor   = .TRUE.  ,   &  ! horizontal (geopotential) direction
(   32)       ln_dynldf_iso   = .FALSE.        ! iso-neutral direction
(   33) 
(   34)    REAL(wp) ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
(   35)       ahm0  = 40000._wp ,   &  ! lateral eddy viscosity (m2/s)
(   36)       ahmb0 =     0._wp        ! lateral background eddy viscosity (m2/s)
(   37) 
(   38)    LOGICAL  ::                      &  ! flag of the lateral diff. scheme used
(   39)       l_dynldf_lap         ,        &  ! iso-level laplacian operator
(   40)       l_dynldf_bilap       ,        &  ! iso-level bilaplacian operator
(   41)       l_dynldf_bilapg      ,        &  ! geopotential bilap. (s-coord)
(   42)       l_dynldf_iso         ,        &  ! iso-neutral laplacian or horizontal lapacian (s-coord)
(   43)       l_dynzdf_iso                     ! iso-neutral laplacian or horizontal lapacian (s-coord)
(   44)    
(   45) 
(   46) # 49
(   49)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  ! ** 3D coefficients **
(   57)       ahm1, ahm2, ahm3, ahm4,tmph                ! ????
(   58) 
(   59)    !!----------------------------------------------------------------------
(   60) END MODULE ldfdyn_oce
