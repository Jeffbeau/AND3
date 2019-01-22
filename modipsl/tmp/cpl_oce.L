


PGF90 (Version     10.5)          09/23/2010  09:49:01      page 1

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

Filename: cpl_oce.F90

(    1) MODULE cpl_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  cpl_oce  ***
(    4)    !! Ocean coupling:  ocean-atmosphere-sea ice coupled exchanges 
(    5)    !!=====================================================================
(    6) #if defined key_coupled
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_coupled'                              Coupled Ocean/Atmosphere
(    9)    !!----------------------------------------------------------------------
(   10)    !! ** Purpose :   Atmosphere/Ice/Ocean Coupling
(   11)    !!      GASTON TEAM (CERFACS, Meteo-France, LMD, LSCE, IPSL, LODYC)
(   12)    !!
(   13)    !! history :
(   14)    !!  8.0   ! 08-98  (M.A. Foujols, M. Imbard)  Original code
(   15)    !!  8.5   ! 06/02  (G. Madec)  modules
(   16)    !!----------------------------------------------------------------------
(   17)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   18)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cpl_oce.F90,v 1.5 2005/03/27 18:34:46 opalod Exp $ 
(   19)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   20)    !!----------------------------------------------------------------------
(   21)    !! * Modules used
(   22)    USE par_oce          ! ocean parameters
(   23) 
(   24)    IMPLICIT NONE
(   25) 
(   26)    !! ---------------------------------------------------------------------
(   27)    !! Ocean/Ice/Atmosphere Coupling
(   28)    !! -----------------------------
(   29)    LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .TRUE.   !: coupled flag
(   30) 
(   31)    INTEGER, PARAMETER ::   &  !:  
(   32)       jpmaxfld = 40           !: Number of maximum fields exchange betwwen
(   33)       !                       ! the ocean and the coupler
(   34) 
(   35)    !!---------------------------------------------------------------------
(   36)    !! SIPC Method   (L. Terray, S. Valcke, CERFACS)
(   37)    !! -----------
(   38) 
(   39)    INTEGER, PARAMETER ::   &  !:
(   40)       jpbyteint = 4,       &  !: number of bytes per integer
(   41)       jpbyterea = 8,       &  !: number of bytes per real
(   42)       jpbytecha = 1           !: number of bytes per character  
(   43) 
(   44)    INTEGER, PARAMETER ::   &  !:
(   45)       jptest = 100            !: The models will test during 2*jptest 
(   46)       !                       ! seconds if the file DUMMY_SIPC has been






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 2

(   47)       !                       ! created by OASIS, signaling that the
(   48)       !                       ! SHM pools are opened. After, it aborts.
(   49) 
(   50)    !!---------------------------------------------------------------------
(   51)    !! PIPE Method   (L. Terray, CERFACS)
(   52)    !! -----------
(   53) 
(   54)    INTEGER, PARAMETER ::   &  !:
(   55)       jpread = 0,          &  !: 
(   56)       jpwrit = 1              !: 
(   57) 
(   58)    !!---------------------------------------------------------------------
(   59)    !! Messag Passing Method (CLIM)
(   60)    !! ----------------------------
(   61) !!!INCLUDE '../../CPL/include/clim.h90'
(   62) !!
(   63) !! -- clim.h   18-08-95   Version 2.0   Author: Laurent Terray
(   64) !!    ******
(   65) !!             26-10-99   Version 2.4   Jean Latour (F.S.E.) MPI-2 support
(   66) !!
(   67) !!   clim.h90  13-08-04  Change to F90 C. Levy
(   68) !!@
(   69) !!@  Contents : variables related to the CLIM library
(   70) !!@  --------
(   71) !!@ For complete definition, see the CLIM manual
(   72) !!@
(   73)       INTEGER (kind=4)	CLIM_MaxMod,    CLIM_MaxPort,  CLIM_MaxSegments, &
(   74)                CLIM_MaxTag, &
(   75)                CLIM_MaxLink, &
(   76)                CLIM_ParSize, & 
(   77)                CLIM_Clength, &
(   78)                CLIM_MaxCodes
(   79) !!
(   80)       INTEGER (kind=4) CLIM_Void
(   81) !!
(   82)       INTEGER (kind=4) CLIM_In,	CLIM_Out,	CLIM_InOut
(   83) !!
(   84)       INTEGER (kind=4) CLIM_Strategy,  CLIM_Segments,  &
(   85)                CLIM_Serial,    CLIM_Length,    CLIM_Orange, &
(   86)                CLIM_Apple,     CLIM_Offset, &
(   87)                CLIM_Box,	CLIM_SizeX,	CLIM_SizeY, &
(   88)                CLIM_LdX
(   89) !!
(   90)       INTEGER  (kind=4)CLIM_Integer,	CLIM_Real,	CLIM_Double
(   91) !!
(   92)       INTEGER  (kind=4)CLIM_StopPvm,   CLIM_ContPvm
(   93) !!
(   94)       INTEGER (kind=4)	CLIM_Ok
(   95)       INTEGER (kind=4) CLIM_FastExit, 	CLIM_BadName, 	CLIM_BadPort, &
(   96)                CLIM_BadType, 	CLIM_DoubleDef, CLIM_NotStep, &
(   97)                CLIM_IncStep, 	CLIM_IncSize, 	CLIM_NotClim, &
(   98)                CLIM_TimeOut, &
(   99)                CLIM_Pvm, 	CLIM_FirstCall, CLIM_PbRoute, &
(  100)                CLIM_Group, 	CLIM_BadTaskId, CLIM_NoTask, &
(  101)                CLIM_InitBuff, 	CLIM_Pack, 	CLIM_Unpack, &
(  102)                CLIM_Down, 	CLIM_PvmExit
(  103) !!
(  104)       INTEGER (kind=4) CLIM_jpmax, 	CLIM_jpmx8, 	CLIM_Mpi






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 3

(  105) !!
(  106) !!-----Parameter sizes
(  107) !!
(  108)       PARAMETER ( CLIM_Void    = 0  )
(  109)       PARAMETER ( CLIM_MaxMod  = 8 )
(  110)       PARAMETER ( CLIM_MaxPort = 50 )
(  111)       PARAMETER ( CLIM_MaxSegments = 50 )
(  112)       PARAMETER ( CLIM_MaxLink = CLIM_MaxMod * CLIM_MaxPort )
(  113)       PARAMETER ( CLIM_ParSize = 2*CLIM_MaxSegments+2 )
(  114)       PARAMETER ( CLIM_MaxTag  = 16777215 )
(  115)       PARAMETER ( CLIM_Clength = 32 )
(  116) !!
(  117) !!-----Dimension of buffer for packing / unpacking messages with MPI
(  118) !!     (must be equal to jpmax of Oasis)
(  119) !!
(  120)       PARAMETER ( CLIM_jpmax = 400000 )
(  121)       PARAMETER ( CLIM_jpmx8 = CLIM_jpmax*8 )
(  122) !!
(  123) !!-----Ports status
(  124) !!
(  125)       PARAMETER ( CLIM_In      = 1 )
(  126)       PARAMETER ( CLIM_Out     = 0 )
(  127)       PARAMETER ( CLIM_InOut   = 2 )
(  128) !!
(  129) !!-----Parallel distribution
(  130) !!
(  131)       PARAMETER ( CLIM_Strategy = 1 )
(  132)       PARAMETER ( CLIM_Segments = 2 )
(  133)       PARAMETER ( CLIM_Serial   = 0 )
(  134)       PARAMETER ( CLIM_Apple    = 1 )
(  135)       PARAMETER ( CLIM_Box      = 2 )
(  136)       PARAMETER ( CLIM_Orange   = 3 )
(  137)       PARAMETER ( CLIM_Offset   = 2 )
(  138)       PARAMETER ( CLIM_Length   = 3 )
(  139)       PARAMETER ( CLIM_SizeX    = 3 )
(  140)       PARAMETER ( CLIM_SizeY    = 4 )
(  141)       PARAMETER ( CLIM_LdX      = 5 )
(  142) !!
(  143) !!-----Datatypes
(  144) !!
(  145)       PARAMETER ( CLIM_Integer = 1 )
(  146)       PARAMETER ( CLIM_Real    = 4 ) 
(  147)       PARAMETER ( CLIM_Double  = 8 )
(  148) !!
(  149) !!-----Quit parameters
(  150) !!
(  151)       PARAMETER ( CLIM_ContPvm = 0 )
(  152)       PARAMETER ( CLIM_StopPvm = 1 )
(  153) !!
(  154) !!-----Error Codes
(  155) !!
(  156)       PARAMETER ( CLIM_MaxCodes  = -22 )
(  157) !!
(  158)       PARAMETER ( CLIM_Ok	 = 0 )
(  159)       PARAMETER ( CLIM_FastExit  = -1 )
(  160)       PARAMETER ( CLIM_BadName   = -2 )
(  161)       PARAMETER ( CLIM_BadPort   = -3 )
(  162)       PARAMETER ( CLIM_BadType   = -4 )






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 4

(  163)       PARAMETER ( CLIM_DoubleDef = -5 )
(  164)       PARAMETER ( CLIM_NotStep   = -6 )
(  165)       PARAMETER ( CLIM_IncStep   = -7 )
(  166)       PARAMETER ( CLIM_IncSize   = -8 )
(  167)       PARAMETER ( CLIM_NotClim   = -9 )
(  168)       PARAMETER ( CLIM_TimeOut   = -10 )
(  169)       PARAMETER ( CLIM_Pvm       = -11 )
(  170)       PARAMETER ( CLIM_FirstCall = -12 )
(  171)       PARAMETER ( CLIM_PbRoute   = -13 )
(  172)       PARAMETER	( CLIM_Group     = -14 )
(  173)       PARAMETER ( CLIM_BadTaskId = -15 )
(  174)       PARAMETER ( CLIM_NoTask    = -16 )
(  175)       PARAMETER ( CLIM_InitBuff  = -17 )
(  176)       PARAMETER ( CLIM_Pack      = -18 )
(  177)       PARAMETER ( CLIM_Unpack    = -19 )
(  178)       PARAMETER ( CLIM_Down      = -20 )
(  179)       PARAMETER ( CLIM_PvmExit   = -21 )
(  180)       PARAMETER ( CLIM_Mpi       = -22 )
(  181) 
(  182) !!
(  183) !     --- end of clim.h90
(  184) !!!END-----------------------------------------------------------------
(  185) 
(  186) !!!INCLUDE '../../CPL/include/mpiclim.h90'
(  187) !!
(  188) !! -- mpiclim.h  26-10-99   Version 2.4   Author: Jean Latour (F.S.E.)
(  189) !!    *********
(  190) !!    mpiclim.h90 13-08-04 change to F90 C. Levy
(  191) !!@
(  192) !!@  Contents : variables related to MPI-2 message passing
(  193) !!@  --------
(  194) !!@
(  195) !!@ -- mpi_totproc: number of processors on which to launch each model
(  196) !!@
(  197) !!@ -- mpi_nproc: number of processors involved in the coupling for
(  198) !!@               each model
(  199) !!@ -- cmpi_modnam: models name
(  200) !!     -----------------------------------------------------------------
(  201) !!
(  202)       INTEGER (kind=4) mpi_totproc(1:CLIM_MaxMod-1),mpi_nproc(0:CLIM_MaxMod-1)
(  203) !!
(  204)       CHARACTER (len=6) cmpi_modnam(1:CLIM_MaxMod-1)
(  205) !!
(  206)       common/CLIM_mpiclim/mpi_totproc, mpi_nproc, cmpi_modnam 
(  207) !!
(  208) !!!END-----------------------------------------------------------------
(  209) 
(  210) 
(  211)    !!----------------------------------------------------------------------
(  212)    !!  Atmosphere/Ice/Ocean Coupling
(  213)    !!---------------------------------
(  214) 
(  215)    REAL(wp), DIMENSION(jpi,jpj) ::   &   !: data from an atmospheric model
(  216)       qc  ,       &  !: total surf. total heat flux (wm-2)
(  217)       ec  ,       &  !: surface water flux (kg m-2s-1)
(  218)       qsrc           !: solar radiation (w m-2)
(  219) 
(  220) #  if defined key_ice_lim






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 5

(  221)    REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
(  222)       watm        ,    &  !:
(  223)       tatm        ,    &  !:
(  224)       hatm        ,    &  !:
(  225)       vatm        ,    &  !:
(  226)       catm                !:
(  227) #  endif
(  228) 
(  229)    !! Coupling
(  230) 
(  231)    INTEGER ::          &  !:
(  232)       npioc       ,    &  !: process-id of ocean PROGRAM
(  233)       nexco       ,    &  !: exchange frequency for fluxes
(  234)       nmodcpl     ,    &  !: coupling mode
(  235)       nflxc2o     ,    &  !: fluxes field number coupler to ocean
(  236)       ntauc2o     ,    &  !: stress field number coupler to ocean
(  237)       nfldo2c             !: surface field number ocean to coupler
(  238) 
(  239)    CHARACTER(len=4) ::   cchan       !: type of message passing (pipe or clim) 
(  240)    CHARACTER(len=6) ::   cplmodnam   !: model name
(  241)    CHARACTER(len=5) ::   cploasis    !: coupler name
(  242) 
(  243)    CHARACTER(len=8), DIMENSION(jpmaxfld) ::   &  !:
(  244)       cpl_f_readflx,   &  !: coupler to ocean file name for flx.coupled
(  245)       cpl_f_readtau,   &  !: coupler to ocean file name for tau.coupled
(  246)       cpl_f_writ   ,   &  !: ocean to coupler file name for cpl_stp
(  247)       cpl_readflx  ,   &  !: coupler to ocean field name for flx.coupled
(  248)       cpl_readtau  ,   &  !: coupler to ocean field name for tau.coupled
(  249)       cpl_writ            !: ocean to coupler field name for cpl_stp
(  250) 
(  251)   REAL(wp), DIMENSION(jpi,jpj) ::   &  !:
(  252)       sstoc,     &  !: work array to average sst
(  253)       sieoc,     &  !: work array to average Ice index
(  254)       alboc,     &  !: work array to average Ice Albedo
(  255)       ticoc         !: work array to average Ice temperature
(  256) 
(  257) 
(  258)    !! -- inc_sipc.h   97-08-11   Version 2.0   Author: S&A
(  259)    !!    **********
(  260)    !! variables describing pools formed of shared memory segments 
(  261) 
(  262)    INTEGER ::   &  !:
(  263)       mpoolinitr,  &  !: handles associated to model pools for passing
(  264)       mpoolinitw      !: initial info (r=read, w=write)
(  265) 
(  266)    INTEGER, DIMENSION(jpmaxfld) ::   &  !:
(  267)       mpoolwrit,   &  !: handles associated to pools used to pass fields
(  268)       !               !  exchanged from model to coupler 
(  269)       !               !  (see libsipc/SIPC_Write_Model.f)
(  270)       mpoolread       !: handles associated to pools used to pass fields
(  271)       !               !  exchanged from model to coupler
(  272)       !               !  (see libsipc/SIPC_Read_Model.f)
(  273) 
(  274) #else
(  275)    !!----------------------------------------------------------------------
(  276)    !!   Default case                                Forced Ocean/Atmosphere
(  277)    !!----------------------------------------------------------------------
(  278)    !!   Empty module






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 6

(  279)    LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .FALSE.   !: coupled flag
(  280) #endif
(  281) 
(  282)    !!----------------------------------------------------------------------
(  283) END MODULE cpl_oce



























































PGF90 (Version     10.5)          09/23/2010  09:49:01      page 7

(    1) # 1 "cpl_oce.F90"
(    1) MODULE cpl_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  cpl_oce  ***
(    4)    !! Ocean coupling:  ocean-atmosphere-sea ice coupled exchanges 
(    5)    !!=====================================================================
(    6) # 275
(  275)    !!----------------------------------------------------------------------
(  276)    !!   Default case                                Forced Ocean/Atmosphere
(  277)    !!----------------------------------------------------------------------
(  278)    !!   Empty module
(  279)    LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .FALSE.   !: coupled flag
(  280) 
(  281) # 282
(  282)    !!----------------------------------------------------------------------
(  283) END MODULE cpl_oce
