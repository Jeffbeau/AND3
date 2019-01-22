


PGF90 (Version     10.5)          09/23/2010  18:36:34      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: par_oce.F90

(    1) MODULE par_oce
(    2)    !!======================================================================
(    3)    !!                        ***  par_oce  ***
(    4)    !! Ocean :   set the ocean parameters
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   4.0  !  91     (Imbard, Levy, Madec)  Original code
(    8)    !!   9.0  !  04-01  (G. Madec, J.-M. Molines)  Free form and module
(    9)    !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   10)    !!----------------------------------------------------------------------
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_oce.F90,v 1.11 2006/03/10 10:55:34 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE par_kind          ! kind parameters
(   17) 
(   18)    IMPLICIT NONE
(   19)    PUBLIC
(   20) 
(   21)    !!----------------------------------------------------------------------
(   22)    !!   Domain decomposition
(   23)    !!----------------------------------------------------------------------
(   24)    !! * if we dont use massively parallel computer (parameters jpni=jpnj=1)
(   25)    !!      so jpiglo=jpi and jpjglo=jpj
(   26) 
(   27)    INTEGER, PUBLIC, PARAMETER ::    &  !:
(   28) !!DB: 2008.10 -- drakes standard config
(   29) !      jpni   = 3,                   &  !: number of processors following i 
(   30) !      jpnj   = 4,                   &  !: number of processors following j
(   31) !      jpnij  = 12,                   &  !: nb of local domain = nb of processors 
(   32) !
(   33) !!JC: sopa2x4
(   34)       jpni   = 2,                   &  !: number of processors following i 
(   35)       jpnj   = 4,                   &  !: number of processors following j
(   36)       jpnij  = 8,                   &  !: nb of local domain = nb of processors 
(   37) 
(   38) !!DB: sopa2x4
(   39) !      jpni   = 1,                   &  !: number of processors following i 
(   40) !      jpnj   = 8,                   &  !: number of processors following j
(   41) !      jpnij  = 8,                   &  !: nb of local domain = nb of processors 
(   42) 
(   43) !!DB: 2009.08 -- 20CPU HPC config
(   44) !      jpni   = 4,                   &  !: number of processors following i 
(   45) !      jpnj   = 5,                   &  !: number of processors following j
(   46) !      jpnij  = 20,                   &  !: nb of local domain = nb of processors 
(   47) 
(   48) 
(   49)       !                                !  ( <= jpni x jpnj )






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 2

(   50)       jpr2di = 0,                   &  !: number of columns for extra outer halo 
(   51)       jpr2dj = 0,                   &  !: number of rows    for extra outer halo 
(   52)       jpreci = 1,                   &  !: number of columns for overlap 
(   53)       jprecj = 1                       !: number of rows    for overlap 
(   54) 
(   55)    !! Ocean Domain sizes
(   56)    !! ------------------
(   57)    !!   data           domain   (jpidta,jpjdta)
(   58)    !!   global or zoom domain   (jpiglo,jpjglo)
(   59)    !!   local          domain   ( jpi  , jpj  )
(   60)    
(   61) !!DB: delete ORCA
(   62) #if defined key_ss_bio
(   63)    !!---------------------------------------------------------------------
(   64)    !!   'key_ss_bio'   :                       Scotian Shelf
(   65)    !!---------------------------------------------------------------------
(   66) #             include "par_SS_R008.h90"
(   67)    !byoung
(   68) #else
(   69)    !!---------------------------------------------------------------------
(   70)    !!   default option  :                               small closed basin
(   71)    !!---------------------------------------------------------------------
(   72)    CHARACTER(len=16), PUBLIC, PARAMETER ::   &  !:
(   73)       cp_cfg = "default"               !: name of the configuration
(   74)    INTEGER, PARAMETER ::            &  !:
(   75)       jp_cfg = 0  ,                 &  !: resolution of the configuration
(   76) 
(   77)       ! data size                     !!! * size of all input files *
(   78)       jpidta  = 10,                 &  !: 1st lateral dimension ( >= jpi )
(   79)       jpjdta  = 12,                 &  !: 2nd    "         "    ( >= jpj )
(   80)       jpkdta  = 31,                 &  !: number of levels      ( >= jpk )
(   81) 
(   82)       ! global or zoom domain size    !!! * computational domain *
(   83)       jpiglo  = jpidta,             &  !: 1st dimension of global domain --> i
(   84)       jpjglo  = jpjdta,             &  !: 2nd    "                  "    --> j
(   85)       jpk     = jpkdta,             &  !: number of vertical levels
(   86)       ! zoom starting position 
(   87)       jpizoom =   1   ,             &  !: left bottom (i,j) indices of the zoom
(   88)       jpjzoom =   1   ,             &  !: in data domain indices
(   89) 
(   90)       ! Domain characteristics
(   91)       jperio  =  0,                 &  !: lateral cond. type (between 0 and 6)
(   92)          !                             !  = 0 closed
(   93)          !                             !  = 1 cyclic East-West
(   94)          !                             !  = 2 equatorial symmetric
(   95)          !                             !  = 3 North fold T-point pivot
(   96)          !                             !  = 4 cyclic East-West AND North fold T-point pivot
(   97)          !                             !  = 5 North fold F-point pivot
(   98)          !                             !  = 6 cyclic East-West AND North fold F-point pivot
(   99)       jpisl   =  0,                 &  !: number of islands (rigid-lid only)
(  100)       jpnisl  =  0                     !: maximum number of points per island
(  101) 
(  102)       !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
(  103)       !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
(  104)       REAL(wp), PARAMETER ::   &  !:
(  105)          pp_not_used       = 999999._wp , &  !:
(  106)          pp_to_be_computed = 999999._wp      !:
(  107) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 3

(  108) 
(  109)    !! Horizontal grid parameters for domhgr
(  110)    !! =====================================
(  111) 
(  112)    INTEGER, PUBLIC, PARAMETER   ::   &  !:
(  113)       jphgr_msh = 0            !: type of horizontal mesh
(  114)       !                        !  = 0 curvilinear coordinate on the sphere
(  115)       !                        !      read in coordinate.nc file
(  116)       !                        !  = 1 geographical mesh on the sphere
(  117)       !                        !      with regular grid-spacing
(  118)       !                        !  = 2 f-plane with regular grid-spacing
(  119)       !                        !  = 3 beta-plane with regular grid-spacing
(  120)       !                        !  = 4 Mercator grid with T/U point at the equator  with
(  121)       !                        !      isotropic resolution (e1_deg)
(  122) 
(  123)    REAL(wp) , PUBLIC, PARAMETER ::   &   !:
(  124)       ppglam0  =    0.0_wp,   &  !: longitude of first raw and column T-point (jphgr_msh = 1)
(  125)       ppgphi0  =  -35.0_wp,   &  !: latitude  of first raw and column T-point (jphgr_msh = 1)
(  126)       !                          !  latitude for the Coriolis or Beta parameter (jphgr_msh = 2 or 3)
(  127)       ppe1_deg =    1.0_wp,   &  !: zonal      grid-spacing (degrees)
(  128)       ppe2_deg =    0.5_wp,   &  !: meridional grid-spacing (degrees)
(  129)       ppe1_m   = 5000.0_wp,   &  !: zonal      grid-spacing (degrees)
(  130)       ppe2_m   = 5000.0_wp       !: meridional grid-spacing (degrees)
(  131) 
(  132)    !! Vertical grid parameter for domzgr
(  133)    !! ==================================
(  134) 
(  135)    REAL(wp), PUBLIC, PARAMETER  ::   &  !:
(  136)       &     ppsur = -4762.96143546300_wp ,  &  !: ORCA r4, r2 and r05 coefficients
(  137)       &     ppa0  =   255.58049070440_wp ,  &  !: (default coefficients)
(  138)       &     ppa1  =   245.58132232490_wp ,  &  !:
(  139)       &     ppkth =    21.43336197938_wp ,  &  !:
(  140)       &     ppacr =     3.00000000000_wp       !:
(  141) 
(  142)    !!  If both ppa0 ppa1 and ppsur are specified to 0, then
(  143)    !!  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
(  144) 
(  145)    REAL(wp), PUBLIC, PARAMETER ::   &  !:
(  146)       &     ppdzmin = 10._wp             ,  &  !: Minimum vertical spacing
(  147)       &     pphmax  = 5000._wp                 !: Maximum depth
(  148) 
(  149)    !!---------------------------------------------------------------------
(  150) #endif
(  151) 
(  152)    !!---------------------------------------------------------------------
(  153)    !! Domain Matrix size
(  154)    !!---------------------------------------------------------------------
(  155)    INTEGER  &  !:
(  156) #if !defined key_agrif
(  157)       ,PARAMETER  &
(  158) #endif
(  159)     :: &
(  160)       jpi = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci ,   &  !: first  dimension
(  161)       jpj = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj ,   &  !: second dimension
(  162)       jpim1 = jpi-1,                                             &  !: inner domain indices
(  163)       jpjm1 = jpj-1,                                             &  !:   "            "
(  164)       jpkm1 = jpk-1,                                             &  !:   "            "
(  165)       jpij  = jpi*jpj                                               !:  jpi x jpj






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 4

(  166) 
(  167) #if defined key_agrif
(  168)    !!---------------------------------------------------------------------
(  169)    !! Agrif variables
(  170)    !!---------------------------------------------------------------------
(  171)    INTEGER, PUBLIC, PARAMETER :: nbghostcells = 1
(  172)    INTEGER, PUBLIC :: nbcellsx = jpiglo - 2 - 2*nbghostcells
(  173)    INTEGER, PUBLIC :: nbcellsy = jpjglo - 2 - 2*nbghostcells
(  174) #endif
(  175)    !!---------------------------------------------------------------------
(  176)    !! Optimization/control flags
(  177)    !!---------------------------------------------------------------------
(  178) #if defined key_esopa
(  179)    LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .TRUE.   !: flag to activate the all options
(  180) #else
(  181)    LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options
(  182) #endif
(  183) 
(  184) #if defined key_vectopt_memory
(  185)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_mem  = .TRUE.   !: vector optimization flag
(  186) #else
(  187)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_mem  = .FALSE.  !: vector optimization flag
(  188) #endif
(  189) 
(  190) #if defined key_vectopt_loop
(  191)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .TRUE.   !: vector optimization flag
(  192) #else
(  193)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .FALSE.  !: vector optimization flag
(  194) #endif
(  195) 
(  196) #if defined key_autotasking
(  197)    LOGICAL, PUBLIC, PARAMETER ::   lk_jki = .TRUE.   !: j-k-i loop flag
(  198) #else
(  199)    LOGICAL, PUBLIC, PARAMETER ::   lk_jki = .FALSE.  !: k-j-i loop flag
(  200) #endif
(  201) 
(  202)    !!======================================================================
(  203) END MODULE par_oce


























PGF90 (Version     10.5)          09/23/2010  18:36:34      page 5

(    1) # 1 "par_oce.F90"
(    1) MODULE par_oce
(    2)    !!======================================================================
(    3)    !!                        ***  par_oce  ***
(    4)    !! Ocean :   set the ocean parameters
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   4.0  !  91     (Imbard, Levy, Madec)  Original code
(    8)    !!   9.0  !  04-01  (G. Madec, J.-M. Molines)  Free form and module
(    9)    !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(   10)    !!----------------------------------------------------------------------
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_oce.F90,v 1.11 2006/03/10 10:55:34 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE par_kind          ! kind parameters
(   17) 
(   18)    IMPLICIT NONE
(   19)    PUBLIC
(   20) 
(   21)    !!----------------------------------------------------------------------
(   22)    !!   Domain decomposition
(   23)    !!----------------------------------------------------------------------
(   24)    !! * if we dont use massively parallel computer (parameters jpni=jpnj=1)
(   25)    !!      so jpiglo=jpi and jpjglo=jpj
(   26) 
(   27)    INTEGER, PUBLIC, PARAMETER ::    &  !:
(   28) !!DB: 2008.10 -- drakes standard config
(   29) !      jpni   = 3,                   &  !: number of processors following i 
(   30) !      jpnj   = 4,                   &  !: number of processors following j
(   31) !      jpnij  = 12,                   &  !: nb of local domain = nb of processors 
(   32) !
(   33) !!JC: sopa2x4
(   34)       jpni   = 2,                   &  !: number of processors following i 
(   35)       jpnj   = 4,                   &  !: number of processors following j
(   36)       jpnij  = 8,                   &  !: nb of local domain = nb of processors 
(   37) 
(   38) !!DB: sopa2x4
(   39) !      jpni   = 1,                   &  !: number of processors following i 
(   40) !      jpnj   = 8,                   &  !: number of processors following j
(   41) !      jpnij  = 8,                   &  !: nb of local domain = nb of processors 
(   42) 
(   43) !!DB: 2009.08 -- 20CPU HPC config
(   44) !      jpni   = 4,                   &  !: number of processors following i 
(   45) !      jpnj   = 5,                   &  !: number of processors following j
(   46) !      jpnij  = 20,                   &  !: nb of local domain = nb of processors 
(   47) 
(   48) 
(   49)       !                                !  ( <= jpni x jpnj )
(   50)       jpr2di = 0,                   &  !: number of columns for extra outer halo 
(   51)       jpr2dj = 0,                   &  !: number of rows    for extra outer halo 
(   52)       jpreci = 1,                   &  !: number of columns for overlap 
(   53)       jprecj = 1                       !: number of rows    for overlap 
(   54) 
(   55)    !! Ocean Domain sizes
(   56)    !! ------------------
(   57)    !!   data           domain   (jpidta,jpjdta)






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 6

(   58)    !!   global or zoom domain   (jpiglo,jpjglo)
(   59)    !!   local          domain   ( jpi  , jpj  )
(   60)    
(   61) !!DB: delete ORCA
(   62) # 63
(   63)    !!---------------------------------------------------------------------
(   64)    !!   'key_ss_bio'   :                       Scotian Shelf
(   65)    !!---------------------------------------------------------------------
(   66) # 1 "./par_SS_R008.h90"
(    1)    !!---------------------------------------------------------------------
(    2)    !!                     ***  par_SS_R008.h90  ***  
(    3)    !!   Ocean Domain : 0.08 degrees resolution Scotia Shelf
(    4)    !!                  (SS_R008 configuration)
(    5)    !!---------------------------------------------------------------------
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/par_ORCA_R008.h90,v 1.4 2005/03/30 10:28:50 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    CHARACTER (len=16), PARAMETER ::   &
(   12)       cp_cfg = "SS"           !: name of the configuration
(   13)    INTEGER, PARAMETER ::     &
(   14)       jp_cfg = 008  ,        &  !: resolution of the configuration (degrees)
(   15)       ! Original data size  !byoung
(   16)       jpidta  = 197,        &  !: first horizontal dimension > or = to jpi
(   17)       jpjdta  = 234,        &  !: second                     > or = to jpj
(   18)       jpkdta  = 46 ,       &  !: number of levels           > or = to jpk
(   19)       ! total domain matrix size
(   20)       jpiglo  = jpidta,      &  !: first  dimension of global domain --> i
(   21)       jpjglo  = jpjdta,      &  !: second dimension of global domain --> j
(   22)       jpk     = jpkdta,      &  !: number of vertical levels
(   23)       ! starting position of the zoom
(   24)       jpizoom =   1   ,      &  !: left bottom (i,j) indices of the zoom
(   25)       jpjzoom =   1   ,      &  !: in data indices
(   26)       ! Domain characteristics
(   27)       jperio  =    0  ,      &  !: lateral cond. type (between 0 and 6)
(   28)       jpisl   =  150  ,      &  !: number of islands
(   29)       jpnisl  = 3000            !: maximum number of points per island
(   30) 
(   31)    !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
(   32)    !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
(   33)    REAL,PARAMETER      ::  pp_not_used = 999999_wp , &
(   34)       &                    pp_to_be_computed = 0._wp
(   35)    !!
(   36)    !! Coefficients associated with the horizontal coordinate system (jphgr_msh /= 0 )
(   37)    !!
(   38)    INTEGER, PARAMETER ::     & !
(   39)       jphgr_msh = 0            !: type of horizontal mesh
(   40)       !                        !  = 0 curvilinear coordinate on the sphere
(   41)       !                        !      read in coordinate.nc file
(   42)       !                        !  = 1 geographical mesh on the sphere
(   43)       !                        !      with regular grid-spacing
(   44)       !                        !  = 2 f-plane with regular grid-spacing
(   45)       !                        !  = 3 beta-plane with regular grid-spacing
(   46)       !                        !  = 4 Mercator grid with T/U point at the equator  with
(   47)       !                        !      isotropic resolution (e1_deg)
(   48) 
(   49)       !   ppglam0 , ppgphi0: coordinates of the lower leftmost T point of the grid.






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 7

(   50)       !   The mercator grid starts only approximately at gphi0 because
(   51)       !   of the constraint that the equator be a T point.
(   52)    REAL(wp), PARAMETER ::       &  !
(   53)       ppglam0  = pp_not_used,   &  !: longitude of first raw and column T-point (jphgr_msh = 1)
(   54)       ppgphi0  = pp_not_used,   &  !: latitude  of first raw and column T-point (jphgr_msh = 1)
(   55)       !                            !  latitude for the Coriolis or Beta parameter (jphgr_msh = 2 or 3)
(   56)       ppe1_deg = pp_not_used,   &  !: zonal      grid-spacing (degrees)
(   57)       ppe2_deg = pp_not_used,   &  !: meridional grid-spacing (degrees)
(   58)       !
(   59)       ppe1_m   = pp_not_used,   &  !: zonal      grid-spacing (meters )
(   60)       ppe2_m   = pp_not_used       !: meridional grid-spacing (meters )
(   61) 
(   62)    !!  Coefficients associated with the vertical coordinate system
(   63) 
(   64)    REAL(wp), PARAMETER  ::       &
(   65)       &     ppsur = pp_to_be_computed ,  &  !: Computed in domzgr, set ppdzmin and pphmax below
(   66)       &     ppa0  = pp_to_be_computed ,  &  !:    "           "
(   67)       &     ppa1  = pp_to_be_computed ,  &  !:    "           "
(   68)       !
(   69)       &     ppkth =  23.563_wp        ,  &  !: (non dimensional): gives the approximate
(   70)       !                                     !: layer number above which  stretching will
(   71)       !                                     !: be maximum. Usually of order jpk/2.
(   72)       &     ppacr =    9.00000000000_wp     !: (non dimensional): stretching factor
(   73)       !                                     !: for the grid. The highest zacr, the smallest
(   74)       !                                     !: the stretching.
(   75) 
(   76)    !!
(   77)    !!  If both ppa0 ppa1 and ppsur are specified to 0, then
(   78)    !!  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
(   79)    !!
(   80)    REAL(wp), PARAMETER ::        &
(   81)       &     ppdzmin = 6._wp           ,  &  !: (meters) vertical thickness of the top layer
(   82)       &     pphmax  = 5750._wp              !: (meters) Maximum depth of the ocean gdepw(jpk)
(   83)    !!---------------------------------------------------------------------
(   84) # 67 "par_oce.F90"
(   67) # 67
(   67)    !byoung
(   68) 
(   69) # 152
(  152)    !!---------------------------------------------------------------------
(  153)    !! Domain Matrix size
(  154)    !!---------------------------------------------------------------------
(  155)    INTEGER  &  !:
(  157)       ,PARAMETER  &
(  159)     :: &
(  160)       jpi = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci ,   &  !: first  dimension
(  161)       jpj = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj ,   &  !: second dimension
(  162)       jpim1 = jpi-1,                                             &  !: inner domain indices
(  163)       jpjm1 = jpj-1,                                             &  !:   "            "
(  164)       jpkm1 = jpk-1,                                             &  !:   "            "
(  165)       jpij  = jpi*jpj                                               !:  jpi x jpj
(  166) 
(  167) # 175
(  175)    !!---------------------------------------------------------------------
(  176)    !! Optimization/control flags
(  177)    !!---------------------------------------------------------------------
(  178) # 181
(  181)    LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 8

(  182) 
(  183) # 187
(  187)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_mem  = .FALSE.  !: vector optimization flag
(  188) 
(  189) # 193
(  193)    LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .FALSE.  !: vector optimization flag
(  194) 
(  195) # 199
(  199)    LOGICAL, PUBLIC, PARAMETER ::   lk_jki = .FALSE.  !: k-j-i loop flag
(  200) 
(  201) # 202
(  202)    !!======================================================================
(  203) END MODULE par_oce
