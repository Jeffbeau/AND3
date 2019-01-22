


PGF90 (Version     10.5)          09/23/2010  09:49:16      page 1

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

Filename: trc.F90

(    1) MODULE trc
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  trc  ***
(    4)    !! Passive tracers   :  module for tracers defined
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.2  !  96-01  (M. Levy)  Original code
(    8)    !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
(    9)    !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
(   10)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  TOP 1.0,  LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trc.F90,v 1.9 2006/04/11 13:48:58 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16) #if defined key_passivetrc
(   17)    !!----------------------------------------------------------------------
(   18)    !!   'key_passivetrc'   :                               Passive tracer
(   19)    !!---------------------------------------------------------------------
(   20)    !! * Modules used
(   21)    USE par_oce
(   22)    USE par_trc
(   23)    IMPLICIT NONE
(   24) 
(   25)    PUBLIC
(   26) 
(   27) 
(   28)    !! passive tracers names and units (read in namelist)
(   29)    !! --------------------------------------------------
(   30)    CHARACTER(len=12), PUBLIC, DIMENSION(jptra) :: &
(   31)       ctrcnm  ,   &   !!: tracer name 
(   32)       ctrcun          !!: tracer unit
(   33) 
(   34)    CHARACTER(len=80), PUBLIC, DIMENSION(jptra) :: &
(   35)       ctrcnl          !!: tracer long name 
(   36)    
(   37)    
(   38)    !! parameters for the control of passive tracers
(   39)    !! --------------------------------------------------
(   40)    INTEGER, PUBLIC ::  &
(   41)       numnat          !!: the number of the passive tracer NAMELIST
(   42)    
(   43)    LOGICAL, PUBLIC, DIMENSION(jptra) ::   &
(   44)       lutini          !!:  initialisation from FILE or not (NAMELIST)
(   45) 
(   46)    INTEGER , PUBLIC, DIMENSION(jptra) :: &






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 2

(   47)       nutini          !!: FORTRAN LOGICAL UNIT for initialisation file
(   48) 
(   49)    !! passive tracers fields (before,now,after)
(   50)    !! --------------------------------------------------
(   51)    REAL(wp), PUBLIC, SAVE  ::  &
(   52)       trai    ,   &   !!: initial total tracer
(   53)       areatot         !!: total volume 
(   54) 
(   55) !!DB
(   56) !   REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk,jptra) :: &
(   57)    REAL(wp), PUBLIC, TARGET, DIMENSION (jpi,jpj,jpk,jptra) :: &
(   58)       trn     ,   &   !!: traceur concentration for actual time step
(   59)       tra     ,   &   !!: traceur concentration for next time step
(   60)       trb             !!: traceur concentration for before time step
(   61) 
(   62) 
(   63)    !! numerical parameter (NAMELIST)
(   64)    !! --------------------------------------------------
(   65)    REAL(wp), PUBLIC  ::  &
(   66)       rsc     ,   &   !!: tuning coefficient for anti-diffusion
(   67)       rtrn            !!: value for truncation
(   68) 
(   69)    !! namelist parameters
(   70)    !! --------------------------------------------------
(   71)    INTEGER , PUBLIC  ::  & 
(   72)       ncortrc ,   &   !!: number of corrective phases
(   73)       ndttrc  ,   &   !!: frequency of step on passive tracers
(   74)       nittrc000       !!: first time step of passive tracers model  
(   75) 
(   76)    LOGICAL, PUBLIC  ::  & 
(   77)       crosster        !!: logical if true computes crossterms
(   78) 
(   79) 
(   80)    !! isopycnal scheme for passive tracers
(   81)    !! --------------------------------------------------  
(   82)    REAL(wp), PUBLIC  ::  &
(   83)       ahtrb0  ,   &   !!: background diffusivity coefficient for passive tracer (m2/s)
(   84)       trcrat  ,   &   !!: ratio between passive and active tracer coeff for diffusion
(   85)       ahtrc0  ,   &   !!: horizontal eddy diffusivity for passive tracers (m2/s)
(   86)       aeivtr0         !!: eddy induced velocity coefficient (m2/s)
(   87)    
(   88)    
(   89)    !! passive tracers restart (input and output)
(   90)    !! --------------------------------------------------  
(   91)    LOGICAL, PUBLIC  ::  &
(   92)       lrsttr          !!: boolean term for restart i/o for passive tracers (namelist)
(   93)    
(   94)    INTEGER , PUBLIC  ::  &
(   95)       nutwrs  ,   &   !!: output FILE for passive tracers restart
(   96)       nutrst  ,   &   !!: logical unit for restart FILE for passive tracers
(   97)       nrsttr          !!: control of the time step ( 0 or 1 ) for pass. tr.
(   98)    
(   99) #if defined key_partial_steps
(  100)    
(  101)    !! interpolated gradient
(  102)    !!--------------------------------------------------  
(  103)    REAL (wp), PUBLIC, DIMENSION (jpi,jpj,jptra) :: &
(  104)       gtru    ,   &   !!: horizontal gradient at u-points at bottom ocean level






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 3

(  105)       gtrv            !!: horizontal gradient at v-points at bottom ocean level
(  106) #else
(  107)    REAL (wp), PUBLIC :: &
(  108)       gtru    ,   &   !!: horizontal gradient at u-points at bottom ocean level
(  109)       gtrv            !!: horizontal gradient at v-points at bottom ocean level
(  110)    
(  111) #endif
(  112)    
(  113) #if defined key_trcldf_eiv && defined key_diaeiv
(  114)    !! The three component of the eddy induced velocity
(  115)    !! --------------------------------------------------
(  116)    REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk) :: &
(  117)       u_trc_eiv,  &   !!: u-eiv (m/s)
(  118)       v_trc_eiv,  &   !!: v-eiv (m/s)
(  119)       w_trc_eiv       !!: w-eiv (m/s)
(  120) #endif
(  121)    
(  122)    
(  123)    !! information for outputs
(  124)    !! --------------------------------------------------
(  125)    INTEGER , PUBLIC   ::  & 
(  126)       nwritetrc       !!: time step frequency for concentration outputs (namelist)
(  127)    
(  128) #if defined key_trc_diaadd
(  129)    !! additional 2D/3D outputs namelist
(  130)    !! --------------------------------------------------
(  131)    CHARACTER(len=8), PUBLIC, DIMENSION (jpdia2d) ::  & 
(  132)       ctrc2d  ,   &   !!: 2d output field name
(  133)       ctrc2u          !!: 2d output field unit
(  134)    
(  135)    CHARACTER(len=8), PUBLIC, DIMENSION (jpdia3d) ::  & 
(  136)       ctrc3d ,    &   !!: 3d output field name
(  137)       ctrc3u          !!: 3d output field unit
(  138)    
(  139)    CHARACTER(len=80), PUBLIC, DIMENSION (jpdia2d) ::  & 
(  140)       ctrc2l          !!: 2d output field long name
(  141)    
(  142)    CHARACTER(len=80), PUBLIC, DIMENSION (jpdia3d) ::  & 
(  143)       ctrc3l          !!: 3d output field long name
(  144)    
(  145)    REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpdia2d) ::  &  
(  146)       trc2d           !!:  additional 2d outputs  
(  147)    
(  148)    REAL(wp), PUBLIC, DIMENSION (jpi,jpj,jpk,jpdia3d) ::  &  
(  149)       trc3d           !!:  additional 3d outputs  
(  150)    
(  151)    
(  152)    !! netcdf files and index common
(  153)    !! --------------------------------------------------
(  154)    INTEGER , PUBLIC :: &
(  155)       nwriteadd     !!: frequency of additional arrays outputs(namelist)
(  156) #endif
(  157)    
(  158) #if defined key_trc_diatrd
(  159)    
(  160)    !!  non conservative trends (biological, ...)
(  161)    !! --------------------------------------------------
(  162)    LOGICAL, PUBLIC, DIMENSION (jptra)  ::  &  






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 4

(  163)       luttrd          !!: large trends diagnostic to write or not (namelist)
(  164)    
(  165)    !!  dynamical trends
(  166)    !!	trtrd()   : trends of the tracer equations
(  167)    !!           1 : X advection
(  168)    !!           2 : Y advection
(  169)    !!           3 : Z advection
(  170)    !!           4 : X diffusion
(  171)    !!           5 : Y diffusion
(  172)    !!           6 : Z diffusion
(  173)    !!           7 : X gent velocity
(  174)    !!           8 : Y gent velocity
(  175)    !!           9 : Z gent velocity
(  176)    !! --------------------------------------------------
(  177)    
(  178)    
(  179)    REAL(wp), PUBLIC, DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE :: &
(  180)       trtrd           !!: trends of the tracer equations
(  181)    
(  182)    INTEGER, PUBLIC, DIMENSION(jptra), SAVE :: ikeep ! indice of tracer for which dyn trends are stored
(  183)    INTEGER, PUBLIC, SAVE                   :: nkeep ! number of tracers for which dyn trends are stored 
(  184)                                                     ! (used to allocate trtrd buffer)
(  185) 
(  186)    !! netcdf files and index common
(  187)    !! --------------------------------------------------
(  188)    INTEGER , PUBLIC :: &
(  189)       nwritetrd       !!: frequency of additional arrays outputs(namelist)
(  190)    
(  191) #endif 
(  192)    
(  193)    !! passive tracers data read and at given time_step
(  194)    !! --------------------------------------------------
(  195) #if defined key_trc_dta
(  196)    
(  197)    INTEGER , PUBLIC, DIMENSION(jptra) :: &
(  198)       numtr          !!: logical unit for passive tracers data
(  199)    
(  200) #endif
(  201) 
(  202)   !!  1D configuration
(  203)   !! --------------------------------------------------
(  204) #if defined key_cfg_1d
(  205)       LOGICAL, PARAMETER ::   lk_trccfg_1d   = .TRUE.   !: 1D pass. tracer configuration flag
(  206) #else   
(  207)       LOGICAL, PARAMETER ::   lk_trccfg_1d   = .FALSE.  !: 1D pass. tracer configuration flag
(  208) #endif
(  209) 
(  210) 
(  211) #else
(  212)    !!======================================================================
(  213)    !!  Empty module : No passive tracer 
(  214)    !!======================================================================
(  215) #endif
(  216) 
(  217) END MODULE trc









PGF90 (Version     10.5)          09/23/2010  09:49:16      page 5

(    1) # 1 "trc.F90"
(    1) MODULE trc
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  trc  ***
(    4)    !! Passive tracers   :  module for tracers defined
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.2  !  96-01  (M. Levy)  Original code
(    8)    !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
(    9)    !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
(   10)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  TOP 1.0,  LOCEAN-IPSL (2005) 
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trc.F90,v 1.9 2006/04/11 13:48:58 opalod Exp $ 
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   15)    !!----------------------------------------------------------------------
(   16) # 212
(  212)    !!======================================================================
(  213)    !!  Empty module : No passive tracer 
(  214)    !!======================================================================
(  215) 
(  216) # 217
(  217) END MODULE trc
