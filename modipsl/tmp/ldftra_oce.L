


PGF90 (Version     10.5)          09/23/2010  09:49:04      page 1

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

Filename: ldftra_oce.F90

(    1) MODULE ldftra_oce
(    2)    !!=====================================================================
(    3)    !!                      ***  MODULE  ldftra_oce  ***
(    4)    !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
(    5)    !!=====================================================================
(    6)    !!
(    7)    !! ** Purpose : - Define in memory lateral tracer mixing coefficients
(    8)    !!
(    9)    !! History :
(   10)    !!   9.0  !  02-11  (G. Madec)  Original code (from common.h)
(   11)    !!----------------------------------------------------------------------
(   12)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_oce.F90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !!----------------------------------------------------------------------
(   23)    !! Lateral eddy diffusivity coefficients (tracers)
(   24)    !!----------------------------------------------------------------------
(   25) 
(   26)    LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_traldf) **
(   27)       ln_traldf_lap   = .TRUE.  ,   &  !: laplacian operator
(   28)       ln_traldf_bilap = .FALSE. ,   &  !: bilaplacian operator
(   29)       ln_traldf_level = .FALSE. ,   &  !: iso-level direction
(   30)       ln_traldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
(   31)       ln_traldf_iso   = .TRUE.         !: iso-neutral direction
(   32) 
(   33)    REAL(wp), PUBLIC ::              & !!: ** lateral mixing namelist (namldf) **
(   34)       aht0  = 2000._wp     ,        &  !: lateral eddy diffusivity (m2/s)
(   35)       ahtb0 =    0._wp     ,        &  !: lateral background eddy diffusivity (m2/s)
(   36)       aeiv0 = 2000._wp                 !: eddy induced velocity coefficient (m2/s)
(   37) 
(   38)    LOGICAL , PUBLIC ::              &  !: flag of the lateral diff. scheme used 
(   39)       l_traldf_lap         ,        &  !: iso-level laplacian operator
(   40)       l_traldf_bilap       ,        &  !: iso-level bilaplacian operator
(   41)       l_traldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
(   42)       l_traldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
(   43)       l_trazdf_iso         ,        &  !: idem for the vertical component
(   44)       l_trazdf_iso_vo      ,        &  !: idem with vectopt_memory
(   45)       l_traldf_iso_zps                 !: iso-neutral laplacian (partial steps)
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 2

(   47) #if defined key_traldf_c3d
(   48)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
(   49) #elif defined key_traldf_smag
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
(   51) #elif defined key_traldf_c2d
(   52)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)     ::   &  !: ** 2D coefficients **
(   53) #elif defined key_traldf_c1d
(   54)    REAL(wp), PUBLIC, DIMENSION(jpk)         ::   &  !: ** 1D coefficients **
(   55) #else
(   56)    REAL(wp), PUBLIC                         ::   &  !: ** 0D coefficients **
(   57) #endif
(   58)       ahtt, ahtu, ahtv, ahtw                !: T-, U-, V-, W-points coefficients
(   59) 
(   60) 
(   61) #if defined key_traldf_eiv
(   62)    !!----------------------------------------------------------------------
(   63)    !!   'key_traldf_eiv'                              eddy induced velocity
(   64)    !!----------------------------------------------------------------------
(   65)    LOGICAL, PUBLIC, PARAMETER ::   lk_traldf_eiv   = .TRUE.   !: eddy induced velocity flag
(   66)       
(   67) # if defined key_traldf_c3d
(   68)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
(   69) # elif defined key_traldf_c2d
(   70)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)     ::   &  !: ** 2D coefficients **
(   71) # elif defined key_traldf_c1d
(   72)    REAL(wp), PUBLIC, DIMENSION(jpk)         ::   &  !: ** 1D coefficients **
(   73) # else
(   74)    REAL(wp), PUBLIC                         ::   &  !: ** 0D coefficients **
(   75) # endif
(   76)       aeiu, aeiv, aeiw                              !: U-, V-, W-points  induced velocity coef. (m2/s)
(   77) 
(   78) # if defined key_diaeiv
(   79)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::    &  !:
(   80)       u_eiv, v_eiv, w_eiv     !: The three component of the eddy induced velocity (m/s)
(   81) # endif
(   82) 
(   83) #else
(   84)    !!----------------------------------------------------------------------
(   85)    !!   Default option :                           NO eddy induced velocity
(   86)    !!----------------------------------------------------------------------
(   87)    LOGICAL , PUBLIC, PARAMETER ::   lk_traldf_eiv   = .FALSE.   !: eddy induced velocity flag
(   88)    REAL(wp), PUBLIC ::   aeiu, aeiv, aeiw
(   89) #endif
(   90) 
(   91)    !!----------------------------------------------------------------------
(   92) END MODULE ldftra_oce


















PGF90 (Version     10.5)          09/23/2010  09:49:04      page 3

(    1) # 1 "ldftra_oce.F90"
(    1) MODULE ldftra_oce
(    2)    !!=====================================================================
(    3)    !!                      ***  MODULE  ldftra_oce  ***
(    4)    !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
(    5)    !!=====================================================================
(    6)    !!
(    7)    !! ** Purpose : - Define in memory lateral tracer mixing coefficients
(    8)    !!
(    9)    !! History :
(   10)    !!   9.0  !  02-11  (G. Madec)  Original code (from common.h)
(   11)    !!----------------------------------------------------------------------
(   12)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_oce.F90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !!----------------------------------------------------------------------
(   23)    !! Lateral eddy diffusivity coefficients (tracers)
(   24)    !!----------------------------------------------------------------------
(   25) 
(   26)    LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_traldf) **
(   27)       ln_traldf_lap   = .TRUE.  ,   &  !: laplacian operator
(   28)       ln_traldf_bilap = .FALSE. ,   &  !: bilaplacian operator
(   29)       ln_traldf_level = .FALSE. ,   &  !: iso-level direction
(   30)       ln_traldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
(   31)       ln_traldf_iso   = .TRUE.         !: iso-neutral direction
(   32) 
(   33)    REAL(wp), PUBLIC ::              & !!: ** lateral mixing namelist (namldf) **
(   34)       aht0  = 2000._wp     ,        &  !: lateral eddy diffusivity (m2/s)
(   35)       ahtb0 =    0._wp     ,        &  !: lateral background eddy diffusivity (m2/s)
(   36)       aeiv0 = 2000._wp                 !: eddy induced velocity coefficient (m2/s)
(   37) 
(   38)    LOGICAL , PUBLIC ::              &  !: flag of the lateral diff. scheme used 
(   39)       l_traldf_lap         ,        &  !: iso-level laplacian operator
(   40)       l_traldf_bilap       ,        &  !: iso-level bilaplacian operator
(   41)       l_traldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
(   42)       l_traldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
(   43)       l_trazdf_iso         ,        &  !: idem for the vertical component
(   44)       l_trazdf_iso_vo      ,        &  !: idem with vectopt_memory
(   45)       l_traldf_iso_zps                 !: iso-neutral laplacian (partial steps)
(   46) 
(   47) # 50
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
(   58)       ahtt, ahtu, ahtv, ahtw                !: T-, U-, V-, W-points coefficients
(   59) 
(   60) 
(   61) # 84
(   84)    !!----------------------------------------------------------------------
(   85)    !!   Default option :                           NO eddy induced velocity
(   86)    !!----------------------------------------------------------------------
(   87)    LOGICAL , PUBLIC, PARAMETER ::   lk_traldf_eiv   = .FALSE.   !: eddy induced velocity flag
(   88)    REAL(wp), PUBLIC ::   aeiu, aeiv, aeiw






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 4

(   89) 
(   90) # 91
(   91)    !!----------------------------------------------------------------------
(   92) END MODULE ldftra_oce
