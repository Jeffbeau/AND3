


PGF90 (Version     10.5)          09/23/2010  18:36:35      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: in_out_manager.F90

(    1) MODULE in_out_manager   
(    2) 
(    3)   USE lib_print         ! formated print library
(    4)   USE par_kind
(    5)   USE par_oce
(    6) 
(    7)   PUBLIC
(    8)   !!----------------------------------------------------------------------
(    9)   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/in_out_manager.F90,v 1.12 2006/04/10 15:46:04 opalod Exp $ 
(   11)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)   !!----------------------------------------------------------------------
(   13) 
(   14)   !!----------------------------------------------------------------------
(   15)   !! namelist parameters
(   16)   !! -------------------------------------
(   17)   ! namrun:  parameters of the run
(   18)   CHARACTER (len=16) ::    &   !:
(   19)        cexper = "exp0"           !: experiment name used for output filename
(   20) 
(   21)   LOGICAL ::   &              !!: * namelist namrun *
(   22)        ln_rstart = .FALSE. ,  &  !: start from (F) rest or (T) a restart file
(   23)        ln_ctl    = .FALSE.       !: run control for debugging
(   24) 
(   25)   INTEGER ::                & !!: * namelist namrun *
(   26)        no     = 0        ,    &  !: job number
(   27)        nrstdt = 0        ,    &  !: control of the time step (0, 1 or 2)
(   28)        nit000 = 1        ,    &  !: index of the first time step
(   29)        nitend = 10       ,    &  !: index of the last time step
(   30)        ndate0 = 961115   ,    &  !: initial calendar date aammjj
(   31)        nleapy = 0        ,    &  !: Leap year calendar flag (0/1 or 30)
(   32)        ninist = 0        ,    &  !: initial state output flag (0/1)
(   33)        nbench = 0                !: benchmark parameter (0/1)
(   34)   !!----------------------------------------------------------------------
(   35)   !!                          Run control  
(   36)   !!----------------------------------------------------------------------
(   37) 
(   38)   INTEGER ::                &  !:
(   39)        nstop = 0 ,            &  !: e r r o r  flag (=number of reason for a
(   40)                                 !                         !                   prematurely stop the run)
(   41)        nwarn = 0                 !: w a r n i n g  flag (=number of warning
(   42)   !                         !                       found during the run)
(   43) 
(   44) 
(   45)   CHARACTER (len=64) ::        &                                                    !:
(   46)        cform_err="(/,' ===>>> : E R R O R',     /,'         ===========',/)"    ,   & !:
(   47)        cform_war="(/,' ===>>> : W A R N I N G', /,'         ===============',/)"      !:
(   48)   !!----------------------------------------------------------------------
(   49)   !! output monitoring






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 2

(   50)   !! -----------------------------------
(   51) 
(   52)   LOGICAL ::   &               !:
(   53)        lwp                ,   &  !: boolean : true on the 1st processor only
(   54)        lsp_area = .TRUE.         !: to make a control print over a specific area
(   55) 
(   56)   INTEGER ::                &  !:
(   57)        nstock =   10 ,        &  !: restart file frequency
(   58)        nprint =    0 ,        &  !: level of print (0 no print)
(   59)        nwrite =   10 ,        &  !: restart file frequency
(   60)        nictls =    0 ,        &  !: Start i indice for the SUM control
(   61)        nictle =    0 ,        &  !: End   i indice for the SUM control
(   62)        njctls =    0 ,        &  !: Start j indice for the SUM control
(   63)        njctle =    0 ,        &  !: End   j indice for the SUM control
(   64)        isplt  =    1 ,        &  !: number of processors following i
(   65)        jsplt  =    1 ,        &  !: number of processors following j
(   66)        ijsplt =    1             !: nb of local domain = nb of processors
(   67) 
(   68)   !!----------------------------------------------------------------------
(   69)   !! logical units
(   70)   !! ------------------------------
(   71)   INTEGER ::                &  !:
(   72)        numstp     =  1 ,      &  !: logical unit for time step
(   73)        numout     =  2 ,      &  !: logical unit for output print
(   74)        numnam     =  3 ,      &  !: logical unit for namelist
(   75)        numnam_ice =  4 ,      &  !: logical unit for ice namelist
(   76)        numevo_ice = 17 ,      &  !: logical unit for ice variables (temp. evolution)
(   77)        numice_dmp = 18 ,      &  !: logical unit for ice variables (damping)
(   78)        numsol     = 25 ,      &  !: logical unit for solver statistics
(   79)        numwri     = 40 ,      &  !: logical unit for output write
(   80)        numisp     = 41 ,      &  !: logical unit for island statistics
(   81)        numgap     = 45 ,      &  !: logical unit for differences diagnostic
(   82)        numwrs     = 46 ,      &  !: logical unit for output restart
(   83)        numtdt     = 62 ,      &  !: logical unit for data temperature
(   84)        numsdt     = 63 ,      &  !: logical unit for data salinity
(   85)        numrnf     = 64 ,      &  !: logical unit for runoff data
(   86)        numwso     = 71 ,      &  !: logical unit for 2d output write
(   87)        numwvo     = 72 ,      &  !: logical unit for 3d output write
(   88)        numsst     = 65 ,      &  !: logical unit for surface temperature data
(   89)        numsss     = 66 ,      &  !: logical unit for surface salinity data
(   90)        numbol     = 67 ,      &  !: logical unit for "bol" diagnostics
(   91)        numptr     = 68 ,      &  !: logical unit for Poleward TRansports
(   92)        numflo     = 69 ,      &  !: logical unit for drifting floats
(   93)                                 !                         !: * coupled units
(   94)        numlhf     = 71 ,      &  !: unit to transfer fluxes
(   95)        numlws     = 72 ,      &  !: unit to transfer stress
(   96)        numlts     = 73 ,      &  !: unit to transfer sst
(   97)        numlic     = 74           !: unit to transfer ice cover
(   98) 
(   99) !!DB: 2009.10.01
(  100)   INTEGER ::   numout2    = 9   ! logical unit for (new) run_info.out file
(  101) 
(  102)   !! Contral/debugging
(  103)   !! -----------------
(  104)   REAL(wp) ::               &  !:
(  105)        u_ctl, v_ctl,          &  !: sum of ua and va trend
(  106)        t_ctl, s_ctl              !: sum of ta and sa trend
(  107) 






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 3

(  108) END MODULE in_out_manager































































PGF90 (Version     10.5)          09/23/2010  18:36:35      page 4

(    1) # 1 "in_out_manager.F90"
(    1) MODULE in_out_manager   
(    2) 
(    3)   USE lib_print         ! formated print library
(    4)   USE par_kind
(    5)   USE par_oce
(    6) 
(    7)   PUBLIC
(    8)   !!----------------------------------------------------------------------
(    9)   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/in_out_manager.F90,v 1.12 2006/04/10 15:46:04 opalod Exp $ 
(   11)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)   !!----------------------------------------------------------------------
(   13) 
(   14)   !!----------------------------------------------------------------------
(   15)   !! namelist parameters
(   16)   !! -------------------------------------
(   17)   ! namrun:  parameters of the run
(   18)   CHARACTER (len=16) ::    &   !:
(   19)        cexper = "exp0"           !: experiment name used for output filename
(   20) 
(   21)   LOGICAL ::   &              !!: * namelist namrun *
(   22)        ln_rstart = .FALSE. ,  &  !: start from (F) rest or (T) a restart file
(   23)        ln_ctl    = .FALSE.       !: run control for debugging
(   24) 
(   25)   INTEGER ::                & !!: * namelist namrun *
(   26)        no     = 0        ,    &  !: job number
(   27)        nrstdt = 0        ,    &  !: control of the time step (0, 1 or 2)
(   28)        nit000 = 1        ,    &  !: index of the first time step
(   29)        nitend = 10       ,    &  !: index of the last time step
(   30)        ndate0 = 961115   ,    &  !: initial calendar date aammjj
(   31)        nleapy = 0        ,    &  !: Leap year calendar flag (0/1 or 30)
(   32)        ninist = 0        ,    &  !: initial state output flag (0/1)
(   33)        nbench = 0                !: benchmark parameter (0/1)
(   34)   !!----------------------------------------------------------------------
(   35)   !!                          Run control  
(   36)   !!----------------------------------------------------------------------
(   37) 
(   38)   INTEGER ::                &  !:
(   39)        nstop = 0 ,            &  !: e r r o r  flag (=number of reason for a
(   40)                                 !                         !                   prematurely stop the run)
(   41)        nwarn = 0                 !: w a r n i n g  flag (=number of warning
(   42)   !                         !                       found during the run)
(   43) 
(   44) 
(   45)   CHARACTER (len=64) ::        &                                                    !:
(   46)        cform_err="(/,' ===>>> : E R R O R',     /,'         ===========',/)"    ,   & !:
(   47)        cform_war="(/,' ===>>> : W A R N I N G', /,'         ===============',/)"      !:
(   48)   !!----------------------------------------------------------------------
(   49)   !! output monitoring
(   50)   !! -----------------------------------
(   51) 
(   52)   LOGICAL ::   &               !:
(   53)        lwp                ,   &  !: boolean : true on the 1st processor only
(   54)        lsp_area = .TRUE.         !: to make a control print over a specific area
(   55) 
(   56)   INTEGER ::                &  !:
(   57)        nstock =   10 ,        &  !: restart file frequency






PGF90 (Version     10.5)          09/23/2010  18:36:35      page 5

(   58)        nprint =    0 ,        &  !: level of print (0 no print)
(   59)        nwrite =   10 ,        &  !: restart file frequency
(   60)        nictls =    0 ,        &  !: Start i indice for the SUM control
(   61)        nictle =    0 ,        &  !: End   i indice for the SUM control
(   62)        njctls =    0 ,        &  !: Start j indice for the SUM control
(   63)        njctle =    0 ,        &  !: End   j indice for the SUM control
(   64)        isplt  =    1 ,        &  !: number of processors following i
(   65)        jsplt  =    1 ,        &  !: number of processors following j
(   66)        ijsplt =    1             !: nb of local domain = nb of processors
(   67) 
(   68)   !!----------------------------------------------------------------------
(   69)   !! logical units
(   70)   !! ------------------------------
(   71)   INTEGER ::                &  !:
(   72)        numstp     =  1 ,      &  !: logical unit for time step
(   73)        numout     =  2 ,      &  !: logical unit for output print
(   74)        numnam     =  3 ,      &  !: logical unit for namelist
(   75)        numnam_ice =  4 ,      &  !: logical unit for ice namelist
(   76)        numevo_ice = 17 ,      &  !: logical unit for ice variables (temp. evolution)
(   77)        numice_dmp = 18 ,      &  !: logical unit for ice variables (damping)
(   78)        numsol     = 25 ,      &  !: logical unit for solver statistics
(   79)        numwri     = 40 ,      &  !: logical unit for output write
(   80)        numisp     = 41 ,      &  !: logical unit for island statistics
(   81)        numgap     = 45 ,      &  !: logical unit for differences diagnostic
(   82)        numwrs     = 46 ,      &  !: logical unit for output restart
(   83)        numtdt     = 62 ,      &  !: logical unit for data temperature
(   84)        numsdt     = 63 ,      &  !: logical unit for data salinity
(   85)        numrnf     = 64 ,      &  !: logical unit for runoff data
(   86)        numwso     = 71 ,      &  !: logical unit for 2d output write
(   87)        numwvo     = 72 ,      &  !: logical unit for 3d output write
(   88)        numsst     = 65 ,      &  !: logical unit for surface temperature data
(   89)        numsss     = 66 ,      &  !: logical unit for surface salinity data
(   90)        numbol     = 67 ,      &  !: logical unit for "bol" diagnostics
(   91)        numptr     = 68 ,      &  !: logical unit for Poleward TRansports
(   92)        numflo     = 69 ,      &  !: logical unit for drifting floats
(   93)                                 !                         !: * coupled units
(   94)        numlhf     = 71 ,      &  !: unit to transfer fluxes
(   95)        numlws     = 72 ,      &  !: unit to transfer stress
(   96)        numlts     = 73 ,      &  !: unit to transfer sst
(   97)        numlic     = 74           !: unit to transfer ice cover
(   98) 
(   99) !!DB: 2009.10.01
(  100)   INTEGER ::   numout2    = 9   ! logical unit for (new) run_info.out file
(  101) 
(  102)   !! Contral/debugging
(  103)   !! -----------------
(  104)   REAL(wp) ::               &  !:
(  105)        u_ctl, v_ctl,          &  !: sum of ua and va trend
(  106)        t_ctl, s_ctl              !: sum of ta and sa trend
(  107) 
(  108) END MODULE in_out_manager
