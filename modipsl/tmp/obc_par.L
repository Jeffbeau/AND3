


PGF90 (Version     10.5)          09/23/2010  09:49:05      page 1

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

Filename: obc_par.F90

(    1) MODULE obc_par
(    2)    !!==============================================================================
(    3)    !!                  ***  MODULE obc_par   ***
(    4)    !! Open Boundary Cond. :   define related parameters
(    5)    !!==============================================================================
(    6) #if defined key_obc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_obc' :                                 Open Boundary Condition
(    9)    !!----------------------------------------------------------------------
(   10)    !! history :
(   11)    !!  8.0   01/91   (CLIPPER)  Original code 
(   12)    !!  9.0   06/02   (C. Talandier)  modules
(   13)    !!        06/04   (F. Durand) jptobc is defined as a parameter, 
(   14)    !!	     	     in order to allow time-dependent OBCs fields on input
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PUBLIC
(   21)    !!----------------------------------------------------------------------
(   22)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_par.F90,v 1.5 2005/12/12 14:20:26 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26)    LOGICAL, PUBLIC, PARAMETER ::   lk_obc = .TRUE.   !: Ocean Boundary Condition flag
(   27) 
(   28)    !!---------------------------------------------------------------------
(   29)    !! open boundary parameter
(   30)    !!---------------------------------------------------------------------
(   31)    INTEGER, PARAMETER ::     &  !: time dimension of the BCS fields on input
(   32)       jptobc  =	      2 
(   33)    !! * EAST open boundary
(   34)    LOGICAL, PARAMETER ::     &  !:
(   35)       lp_obc_east = .TRUE.     !: to active or not the East open boundary !byoung
(   36)    INTEGER, PARAMETER ::     &  !:
(   37)       jpieob  = jpiglo-2,    &  !: i-localization of the East open boundary (must be ocean U-point)
(   38)       jpjed   =        2,    &  !: j-starting indice of the East open boundary (must be land T-point)
(   39)       jpjef   = jpjglo-1,    &  !: j-ending   indice of the East open boundary (must be land T-point)
(   40)       jpjedp1 =  jpjed+1,    &  !: first ocean point         "                 "
(   41)       jpjefm1 =  jpjef-1        !: last  ocean point         "                 "
(   42) 
(   43)    !! * WEST open boundary
(   44)    LOGICAL, PARAMETER ::     &  !:
(   45)       lp_obc_west = .TRUE.     !: to active or not the West open boundary
(   46)    INTEGER, PARAMETER ::     &  !:






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 2

(   47)       jpiwob  =	       2,    &  !: i-localization of the West open boundary (must be ocean U-point)
(   48)       jpjwd   =	       2,    &  !: j-starting indice of the West open boundary (must be land T-point)
(   49)       jpjwf   = jpjglo-1,    &  !: j-ending   indice of the West open boundary (must be land T-point)
(   50)       jpjwdp1 =  jpjwd+1,    &  !: first ocean point         "                 "
(   51)       jpjwfm1 =  jpjwf-1        !: last  ocean point         "                 "
(   52) 
(   53)    !! * NORTH open boundary
(   54)    LOGICAL, PARAMETER ::     &  !:
(   55)       lp_obc_north = .FALSE.    !: to active or not the North open boundary
(   56)    INTEGER, PARAMETER ::     &  !:
(   57)       jpjnob  = jpjglo-2,    &  !: j-localization of the North open boundary (must be ocean V-point)
(   58)       jpind   =        2,    &  !: i-starting indice of the North open boundary (must be land T-point)
(   59)       jpinf   = jpiglo-1,    &  !: i-ending   indice of the North open boundary (must be land T-point)
(   60)       jpindp1 =  jpind+1,    &  !: first ocean point         "                 "
(   61)       jpinfm1 =  jpinf-1        !: last  ocean point         "                 "
(   62) 
(   63)    !! * SOUTH open boundary
(   64)    LOGICAL, PARAMETER ::     &  !:
(   65)       lp_obc_south = .TRUE.    !: to active or not the South open boundary
(   66)    INTEGER, PARAMETER ::     &  !:
(   67)       jpjsob  =        2,    &  !: j-localization of the South open boundary (must be ocean V-point)
(   68)       jpisd   =        2,    &  !: i-starting indice of the South open boundary (must be land T-point)
(   69)       jpisf   = jpiglo-1,    &  !: i-ending   indice of the South open boundary (must be land T-point)
(   70)       jpisdp1 =  jpisd+1,    &  !: first ocean point         "                 "
(   71)       jpisfm1 =  jpisf-1        !: last  ocean point         "                 "
(   72)    
(   73)    INTEGER, PARAMETER ::     &  !:
(   74)       jpnic = 2700              !: maximum number of isolated coastlines points 
(   75) 
(   76) 
(   77) #else
(   78)    !!----------------------------------------------------------------------
(   79)    !!   Default option :                         NO open boundary condition
(   80)    !!----------------------------------------------------------------------
(   81)    LOGICAL, PUBLIC, PARAMETER ::   lk_obc = .FALSE.  !: Ocean Boundary Condition flag
(   82) #endif
(   83) 
(   84)    !!======================================================================
(   85) END MODULE obc_par

























PGF90 (Version     10.5)          09/23/2010  09:49:05      page 3

(    1) # 1 "obc_par.F90"
(    1) MODULE obc_par
(    2)    !!==============================================================================
(    3)    !!                  ***  MODULE obc_par   ***
(    4)    !! Open Boundary Cond. :   define related parameters
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_obc' :                                 Open Boundary Condition
(    9)    !!----------------------------------------------------------------------
(   10)    !! history :
(   11)    !!  8.0   01/91   (CLIPPER)  Original code 
(   12)    !!  9.0   06/02   (C. Talandier)  modules
(   13)    !!        06/04   (F. Durand) jptobc is defined as a parameter, 
(   14)    !!	     	     in order to allow time-dependent OBCs fields on input
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE par_oce         ! ocean parameters
(   18) 
(   19)    IMPLICIT NONE
(   20)    PUBLIC
(   21)    !!----------------------------------------------------------------------
(   22)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_par.F90,v 1.5 2005/12/12 14:20:26 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26)    LOGICAL, PUBLIC, PARAMETER ::   lk_obc = .TRUE.   !: Ocean Boundary Condition flag
(   27) 
(   28)    !!---------------------------------------------------------------------
(   29)    !! open boundary parameter
(   30)    !!---------------------------------------------------------------------
(   31)    INTEGER, PARAMETER ::     &  !: time dimension of the BCS fields on input
(   32)       jptobc  =	      2 
(   33)    !! * EAST open boundary
(   34)    LOGICAL, PARAMETER ::     &  !:
(   35)       lp_obc_east = .TRUE.     !: to active or not the East open boundary !byoung
(   36)    INTEGER, PARAMETER ::     &  !:
(   37)       jpieob  = jpiglo-2,    &  !: i-localization of the East open boundary (must be ocean U-point)
(   38)       jpjed   =        2,    &  !: j-starting indice of the East open boundary (must be land T-point)
(   39)       jpjef   = jpjglo-1,    &  !: j-ending   indice of the East open boundary (must be land T-point)
(   40)       jpjedp1 =  jpjed+1,    &  !: first ocean point         "                 "
(   41)       jpjefm1 =  jpjef-1        !: last  ocean point         "                 "
(   42) 
(   43)    !! * WEST open boundary
(   44)    LOGICAL, PARAMETER ::     &  !:
(   45)       lp_obc_west = .TRUE.     !: to active or not the West open boundary
(   46)    INTEGER, PARAMETER ::     &  !:
(   47)       jpiwob  =	       2,    &  !: i-localization of the West open boundary (must be ocean U-point)
(   48)       jpjwd   =	       2,    &  !: j-starting indice of the West open boundary (must be land T-point)
(   49)       jpjwf   = jpjglo-1,    &  !: j-ending   indice of the West open boundary (must be land T-point)
(   50)       jpjwdp1 =  jpjwd+1,    &  !: first ocean point         "                 "
(   51)       jpjwfm1 =  jpjwf-1        !: last  ocean point         "                 "
(   52) 
(   53)    !! * NORTH open boundary
(   54)    LOGICAL, PARAMETER ::     &  !:
(   55)       lp_obc_north = .FALSE.    !: to active or not the North open boundary
(   56)    INTEGER, PARAMETER ::     &  !:
(   57)       jpjnob  = jpjglo-2,    &  !: j-localization of the North open boundary (must be ocean V-point)






PGF90 (Version     10.5)          09/23/2010  09:49:05      page 4

(   58)       jpind   =        2,    &  !: i-starting indice of the North open boundary (must be land T-point)
(   59)       jpinf   = jpiglo-1,    &  !: i-ending   indice of the North open boundary (must be land T-point)
(   60)       jpindp1 =  jpind+1,    &  !: first ocean point         "                 "
(   61)       jpinfm1 =  jpinf-1        !: last  ocean point         "                 "
(   62) 
(   63)    !! * SOUTH open boundary
(   64)    LOGICAL, PARAMETER ::     &  !:
(   65)       lp_obc_south = .TRUE.    !: to active or not the South open boundary
(   66)    INTEGER, PARAMETER ::     &  !:
(   67)       jpjsob  =        2,    &  !: j-localization of the South open boundary (must be ocean V-point)
(   68)       jpisd   =        2,    &  !: i-starting indice of the South open boundary (must be land T-point)
(   69)       jpisf   = jpiglo-1,    &  !: i-ending   indice of the South open boundary (must be land T-point)
(   70)       jpisdp1 =  jpisd+1,    &  !: first ocean point         "                 "
(   71)       jpisfm1 =  jpisf-1        !: last  ocean point         "                 "
(   72)    
(   73)    INTEGER, PARAMETER ::     &  !:
(   74)       jpnic = 2700              !: maximum number of isolated coastlines points 
(   75) 
(   76) 
(   77) 
(   78) # 84
(   84)    !!======================================================================
(   85) END MODULE obc_par
