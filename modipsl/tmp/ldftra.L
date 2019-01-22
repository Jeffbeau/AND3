


PGF90 (Version     10.5)          09/23/2010  09:50:04      page 1

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

Filename: ldftra.F90

(    1) MODULE ldftra
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldftra  ***
(    4)    !! Ocean physics:  lateral diffusivity coefficient 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   ldf_tra_init : initialization, namelist read, and parameters control
(    9)    !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
(   10)    !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
(   11)    !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE phycst          ! physical constants
(   17)    USE ldftra_oce      ! ocean tracer   lateral physics
(   18)    USE ldfslp          ! ???
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! *  Routine accessibility
(   27)    PUBLIC ldf_tra_init   ! called by opa.F90
(   28) !!DB: due to new code changes must make the relevant routine PUBLIC or step() cannot find it
(   29) #if defined key_traldf_c3d
(   30)    PUBLIC ldf_tra_c3d
(   31) #elif defined key_traldf_c2d
(   32)    PUBLIC ldf_tra_c2d
(   33) #elif defined key_traldf_smag
(   34)    PUBLIC ldf_tra_smag
(   35) #elif defined key_traldf_c1d
(   36)    PUBLIC ldf_tra_c1d
(   37) #endif
(   38) 
(   39)    !! * Substitutions
(   40) #  include "domzgr_substitute.h90"
(   41) #  include "vectopt_loop_substitute.h90"
(   42)    !!---------------------------------------------------------------------------------
(   43)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   44)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra.F90,v 1.4 2005/03/27 18:35:07 opalod Exp $ 
(   45)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   46)    !!---------------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 2

(   47) 
(   48) CONTAINS
(   49) 
(   50)    SUBROUTINE ldf_tra_init
(   51)       !!----------------------------------------------------------------------
(   52)       !!                  ***  ROUTINE ldf_tra_init  ***
(   53)       !! 
(   54)       !! ** Purpose :   initializations of the horizontal ocean tracer physics
(   55)       !!
(   56)       !! ** Method :
(   57)       !!      Direction of lateral diffusion (tracers and/or momentum)
(   58)       !!        ln_traldf_iso  = T : initialize the slope arrays to zero
(   59)       !!        ln_traldf_geop = T : initialise the slope arrays to the i- and
(   60)       !!                            j-slopes of s-surfaces
(   61)       !!      Eddy diffusivity and eddy induced velocity cefficients:
(   62)       !!         default option   : constant coef. aht0, aeiv0 (namelist)
(   63)       !!        'key_traldf_c1d': depth dependent coef. defined in 
(   64)       !!                            in ldf_tra_c1d routine
(   65)       !!        'key_traldf_c2d': latitude and longitude dependent coef.
(   66)       !!                            defined in ldf_tra_c2d routine
(   67)       !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
(   68)       !!                            defined in ldf_tra_c3d routine
(   69)       !!
(   70)       !!      N.B. User defined include files.  By default, 3d and 2d coef.
(   71)       !!      are set to a constant value given in the namelist and the 1d
(   72)       !!      coefficients are initialized to a hyperbolic tangent vertical
(   73)       !!      profile.
(   74)       !!
(   75)       !! Reference :
(   76)       !!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
(   77)       !!      the North Pole singularity, Climate Dynamics, 12, 381-388.
(   78)       !!
(   79)       !! History :
(   80)       !!        !  07-97  (G. Madec)  from inimix.F split in 2 routines
(   81)       !!        !  08-97  (G. Madec)  multi dimensional coefficients
(   82)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   83)       !!----------------------------------------------------------------------
(   84)       !! * Modules used
(   85)       USE ioipsl
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER ::   ioptio               ! ???
(   89)       LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
(   90)        
(   91)       NAMELIST/nam_traldf/ ln_traldf_lap  , ln_traldf_bilap,                &
(   92)          &                 ln_traldf_level, ln_traldf_hor, ln_traldf_iso,   &
(   93)          &                 aht0, ahtb0, aeiv0
(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       !  Define the lateral tracer physics parameters
(   97)       ! =============================================
(   98)     
(   99)       ! Read Namelist nam_traldf : Lateral physics on tracers
(  100)       REWIND( numnam )
(  101)       READ  ( numnam, nam_traldf )
(  102) 
(  103)       IF(lwp) THEN
(  104)          WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 3

(  105)          WRITE(numout,*) 'ldf_tra : lateral tracer physics'
(  106)          WRITE(numout,*) '~~~~~~~'
(  107)          WRITE(numout,*) '          Namelist nam_traldf : set lateral mixing parameters (type, direction, coefficients)'
(  108)          WRITE(numout,*) '             laplacian operator          ln_traldf_lap   = ', ln_traldf_lap
(  109)          WRITE(numout,*) '             bilaplacian operator        ln_traldf_bilap = ', ln_traldf_bilap
(  110)          WRITE(numout,*) '             iso-level                   ln_traldf_level = ', ln_traldf_level
(  111)          WRITE(numout,*) '             horizontal (geopotential)   ln_traldf_hor   = ', ln_traldf_hor
(  112)          WRITE(numout,*) '             iso-neutral                 ln_traldf_iso   = ', ln_traldf_iso
(  113)          WRITE(numout,*) '             lateral eddy diffusivity             aht0   = ', aht0
(  114)          WRITE(numout,*) '             background hor. diffusivity          ahtb0  = ', ahtb0
(  115)          WRITE(numout,*) '             eddy induced velocity coef.          aeiv0  = ', aeiv0
(  116)          WRITE(numout,*)
(  117)       ENDIF
(  118) 
(  119)       ! Parameter control
(  120) 
(  121)       ! control the input
(  122)       ioptio = 0
(  123)       IF( ln_traldf_lap   )   ioptio = ioptio + 1
(  124)       IF( ln_traldf_bilap )   ioptio = ioptio + 1
(  125)       IF( ioptio /= 1 )   THEN
(  126)           IF(lwp) WRITE(numout,cform_err)
(  127)           IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on tracer'
(  128)           nstop = nstop + 1
(  129)       ENDIF
(  130)       ioptio = 0
(  131)       IF( ln_traldf_level )   ioptio = ioptio + 1
(  132)       IF( ln_traldf_hor   )   ioptio = ioptio + 1
(  133)       IF( ln_traldf_iso   )   ioptio = ioptio + 1
(  134)       IF( ioptio /= 1 ) THEN
(  135)          IF(lwp) WRITE(numout,cform_err)
(  136)          IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
(  137)          nstop = nstop + 1
(  138)       ENDIF
(  139) 
(  140)       ! ... Choice of the lateral scheme used
(  141)       IF( lk_traldf_eiv ) THEN
(  142)          IF(lwp) WRITE(numout,*) '          eddy induced velocity on tracers'
(  143)             IF( .NOT.ln_traldf_iso .OR. ln_traldf_bilap ) THEN
(  144)             IF(lwp) WRITE(numout,cform_err)
(  145)             IF(lwp) WRITE(numout,*) ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion'
(  146)             nstop = nstop + 1
(  147)          ENDIF
(  148)       ENDIF
(  149) 
(  150)       IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal mixing
(  151)          IF( ( ln_traldf_iso .OR. ln_traldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
(  152)             IF(lwp) WRITE(numout,cform_err)
(  153)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  154)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  155)          ENDIF
(  156)       ELSE                       ! z-coordinates with/without partial step:
(  157)          ln_traldf_level = ln_traldf_level .OR. ln_traldf_hor      ! level diffusion = horizontal diffusion
(  158)          ln_traldf_hor   = .FALSE.
(  159)          IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_traldf_level = T'
(  160)          IF(lwp) WRITE(numout,*) '                                                  and    force ln_traldf_hor   = F'
(  161)          IF( ln_traldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
(  162)             IF(lwp) WRITE(numout,cform_err)






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 4

(  163)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  164)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  165)          ENDIF
(  166)       ENDIF
(  167) 
(  168)       l_traldf_lap     =       ln_traldf_lap   .AND. ln_traldf_level     ! iso-level   laplacian operator
(  169)       l_traldf_bilap   =       ln_traldf_bilap .AND. ln_traldf_level     ! iso-level bilaplacian operator
(  170)       l_traldf_bilapg  =       ln_traldf_bilap .AND. ln_traldf_hor       ! geopotential bilap. (s-coord)
(  171)       l_traldf_iso     =       ln_traldf_lap   .AND.                  &  ! laplacian operator
(  172)          &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (z-coord) or horizontal (s-coord)
(  173)          &                                     .AND. .NOT.lk_zps
(  174)       l_traldf_iso_zps =       ln_traldf_lap   .AND.                  &  ! laplacian operator
(  175)          &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (partial steps)
(  176)          &                                     .AND. lk_zps              ! or geopotential in mixed partial steps/s-coord
(  177)       l_trazdf_iso    = .FALSE.
(  178)       l_trazdf_iso_vo = .FALSE.
(  179)       IF( l_traldf_iso     )   l_trazdf_iso = .TRUE.
(  180)       IF( l_traldf_iso_zps )   l_trazdf_iso = .TRUE.
(  181) #if defined key_vectopt_memory
(  182)       IF( l_trazdf_iso ) THEN
(  183)          l_trazdf_iso    = .FALSE.
(  184)          l_trazdf_iso_vo = .TRUE.
(  185)       ENDIF
(  186) #endif
(  187) 
(  188)       ioptio = 0
(  189)       IF( l_traldf_lap     )   ioptio = ioptio + 1
(  190)       IF( l_traldf_bilap   )   ioptio = ioptio + 1
(  191)       IF( l_traldf_bilapg  )   ioptio = ioptio + 1
(  192)       IF( l_traldf_iso     )   ioptio = ioptio + 1
(  193)       IF( l_traldf_iso_zps )   ioptio = ioptio + 1
(  194)       IF( ioptio /= 1 ) THEN
(  195)          IF(lwp) WRITE(numout,cform_err)
(  196)          IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
(  197)          nstop = nstop + 1
(  198)       ENDIF
(  199)       IF( lk_esopa ) THEN
(  200)          l_traldf_lap = .TRUE.   ;   l_traldf_bilap   = .TRUE.   ;   l_traldf_bilapg  = .TRUE.
(  201)          l_traldf_iso = .TRUE.   ;   l_traldf_iso_zps = .TRUE.
(  202)          l_trazdf_iso = .TRUE.   ;   l_trazdf_iso_vo  = .TRUE.
(  203)          IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
(  204)       ENDIF
(  205) 
(  206) 
(  207)       ! ... Space variation of eddy coefficients
(  208)       ioptio = 0
(  209) #if defined key_traldf_c3d
(  210)       IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
(  211)       ioptio = ioptio + 1
(  212) #endif
(  213) #if defined key_traldf_smag
(  214)       IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
(  215)       ioptio = ioptio + 1
(  216) #endif
(  217) #if defined key_traldf_c2d
(  218)       IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude)'
(  219)       ioptio = ioptio + 1
(  220) #endif






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 5

(  221) #if defined key_traldf_c1d
(  222)       IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( depth )'
(  223)       ioptio = ioptio + 1
(  224)       IF( lk_sco ) THEN
(  225)          IF(lwp) WRITE(numout,cform_err)
(  226)          IF(lwp) WRITE(numout,*) '          key_traldf_c1d cannot be used in s-coordinate (key_s_coord)'
(  227)          nstop = nstop + 1
(  228)       ENDIF
(  229) #endif
(  230)       IF( ioptio == 0 ) THEN
(  231)           IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'
(  232)         ELSEIF( ioptio > 1 ) THEN
(  233)           IF(lwp) WRITE(numout,cform_err)
(  234)           IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
(  235)              &                    ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d'
(  236)           nstop = nstop + 1
(  237)       ENDIF
(  238) 
(  239)       IF( l_traldf_bilap .OR. l_traldf_bilapg ) THEN
(  240)          IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
(  241)          IF( aht0 > 0 .AND. .NOT. lk_esopa ) THEN
(  242)             IF(lwp) WRITE(numout,cform_err)
(  243)             IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be negative'
(  244)             nstop = nstop + 1
(  245)          ENDIF
(  246)       ELSE
(  247)          IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
(  248)          IF( aht0 < 0 .AND. .NOT. lk_esopa ) THEN
(  249)             IF(lwp) WRITE(numout,cform_err)
(  250)             IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be positive'
(  251)             nstop = nstop + 1
(  252)          ENDIF
(  253)       ENDIF
(  254) 
(  255) 
(  256)       !  Lateral eddy diffusivity and eddy induced velocity coefficients
(  257)       ! ================================================================
(  258) 
(  259) #if defined key_traldf_c3d
(  260)       CALL ldf_tra_c3d( ll_print )           ! aht = 3D coef. = F( longitude, latitude, depth )
(  261) #elif defined key_traldf_c2d
(  262)       CALL ldf_tra_c2d( ll_print )           ! aht = 2D coef. = F( longitude, latitude )
(  263) #elif defined key_traldf_smag
(  264)       CALL ldf_tra_smag( ll_print )           ! aht = 2D coef. = F( longitude, latitude )
(  265) #elif defined key_traldf_c1d
(  266)       CALL ldf_tra_c1d( ll_print )           ! aht = 1D coef. = F( depth )
(  267) #else
(  268)                                      ! Constant coefficients
(  269)       IF(lwp)WRITE(numout,*)
(  270)       IF(lwp)WRITE(numout,*) ' inildf: constant eddy diffusivity coef.'
(  271)       IF(lwp)WRITE(numout,*) ' ~~~~~~'
(  272)       IF(lwp)WRITE(numout,*) '        ahtu = ahtv = ahtw = aht0 = ', aht0
(  273)       IF( lk_traldf_eiv ) THEN
(  274)          IF(lwp)WRITE(numout,*)
(  275)          IF(lwp)WRITE(numout,*) ' inildf: constant eddy induced velocity coef.'
(  276)          IF(lwp)WRITE(numout,*) ' ~~~~~~  '
(  277)          IF(lwp)WRITE(numout,*) '         aeiu = aeiv = aeiw = aeiv0 = ', aeiv0
(  278)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 6

(  279) #endif
(  280) 
(  281)    END SUBROUTINE ldf_tra_init
(  282) 
(  283) #if defined key_traldf_c3d
(  284) #   include "ldftra_c3d.h90"
(  285) #elif defined key_traldf_c2d
(  286) #   include "ldftra_c2d.h90"
(  287) #elif defined key_traldf_smag
(  288) #   include "ldftra_smag.h90"
(  289) #elif defined key_traldf_c1d
(  290) #   include "ldftra_c1d.h90"
(  291) #endif
(  292) 
(  293)    !!======================================================================
(  294) END MODULE ldftra
















































PGF90 (Version     10.5)          09/23/2010  09:50:04      page 7

(    1) # 1 "ldftra.F90"
(    1) MODULE ldftra
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldftra  ***
(    4)    !! Ocean physics:  lateral diffusivity coefficient 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   ldf_tra_init : initialization, namelist read, and parameters control
(    9)    !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
(   10)    !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
(   11)    !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE phycst          ! physical constants
(   17)    USE ldftra_oce      ! ocean tracer   lateral physics
(   18)    USE ldfslp          ! ???
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! *  Routine accessibility
(   27)    PUBLIC ldf_tra_init   ! called by opa.F90
(   28) !!DB: due to new code changes must make the relevant routine PUBLIC or step() cannot find it
(   29) # 34
(   34)    PUBLIC ldf_tra_smag
(   35) 
(   36) # 39
(   39)    !! * Substitutions
(   40) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 41 "ldftra.F90"
(   41) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 8

(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 42 "ldftra.F90"
(   42) # 42
(   42)    !!---------------------------------------------------------------------------------
(   43)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   44)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra.F90,v 1.4 2005/03/27 18:35:07 opalod Exp $ 
(   45)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   46)    !!---------------------------------------------------------------------------------
(   47) 
(   48) CONTAINS
(   49) 
(   50)    SUBROUTINE ldf_tra_init
(   51)       !!----------------------------------------------------------------------
(   52)       !!                  ***  ROUTINE ldf_tra_init  ***
(   53)       !! 
(   54)       !! ** Purpose :   initializations of the horizontal ocean tracer physics
(   55)       !!
(   56)       !! ** Method :
(   57)       !!      Direction of lateral diffusion (tracers and/or momentum)
(   58)       !!        ln_traldf_iso  = T : initialize the slope arrays to zero
(   59)       !!        ln_traldf_geop = T : initialise the slope arrays to the i- and
(   60)       !!                            j-slopes of s-surfaces
(   61)       !!      Eddy diffusivity and eddy induced velocity cefficients:
(   62)       !!         default option   : constant coef. aht0, aeiv0 (namelist)
(   63)       !!        'key_traldf_c1d': depth dependent coef. defined in 
(   64)       !!                            in ldf_tra_c1d routine
(   65)       !!        'key_traldf_c2d': latitude and longitude dependent coef.
(   66)       !!                            defined in ldf_tra_c2d routine
(   67)       !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
(   68)       !!                            defined in ldf_tra_c3d routine
(   69)       !!
(   70)       !!      N.B. User defined include files.  By default, 3d and 2d coef.
(   71)       !!      are set to a constant value given in the namelist and the 1d
(   72)       !!      coefficients are initialized to a hyperbolic tangent vertical
(   73)       !!      profile.
(   74)       !!
(   75)       !! Reference :
(   76)       !!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
(   77)       !!      the North Pole singularity, Climate Dynamics, 12, 381-388.
(   78)       !!
(   79)       !! History :
(   80)       !!        !  07-97  (G. Madec)  from inimix.F split in 2 routines
(   81)       !!        !  08-97  (G. Madec)  multi dimensional coefficients
(   82)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   83)       !!----------------------------------------------------------------------
(   84)       !! * Modules used
(   85)       USE ioipsl
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER ::   ioptio               ! ???
(   89)       LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
(   90)        
(   91)       NAMELIST/nam_traldf/ ln_traldf_lap  , ln_traldf_bilap,                &
(   92)          &                 ln_traldf_level, ln_traldf_hor, ln_traldf_iso,   &
(   93)          &                 aht0, ahtb0, aeiv0






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 9

(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       !  Define the lateral tracer physics parameters
(   97)       ! =============================================
(   98)     
(   99)       ! Read Namelist nam_traldf : Lateral physics on tracers
(  100)       REWIND( numnam )
(  101)       READ  ( numnam, nam_traldf )
(  102) 
(  103)       IF(lwp) THEN
(  104)          WRITE(numout,*)
(  105)          WRITE(numout,*) 'ldf_tra : lateral tracer physics'
(  106)          WRITE(numout,*) '~~~~~~~'
(  107)          WRITE(numout,*) '          Namelist nam_traldf : set lateral mixing parameters (type, direction, coefficients)'
(  108)          WRITE(numout,*) '             laplacian operator          ln_traldf_lap   = ', ln_traldf_lap
(  109)          WRITE(numout,*) '             bilaplacian operator        ln_traldf_bilap = ', ln_traldf_bilap
(  110)          WRITE(numout,*) '             iso-level                   ln_traldf_level = ', ln_traldf_level
(  111)          WRITE(numout,*) '             horizontal (geopotential)   ln_traldf_hor   = ', ln_traldf_hor
(  112)          WRITE(numout,*) '             iso-neutral                 ln_traldf_iso   = ', ln_traldf_iso
(  113)          WRITE(numout,*) '             lateral eddy diffusivity             aht0   = ', aht0
(  114)          WRITE(numout,*) '             background hor. diffusivity          ahtb0  = ', ahtb0
(  115)          WRITE(numout,*) '             eddy induced velocity coef.          aeiv0  = ', aeiv0
(  116)          WRITE(numout,*)
(  117)       ENDIF
(  118) 
(  119)       ! Parameter control
(  120) 
(  121)       ! control the input
(  122)       ioptio = 0
(  123)       IF( ln_traldf_lap   )   ioptio = ioptio + 1
(  124)       IF( ln_traldf_bilap )   ioptio = ioptio + 1
(  125)       IF( ioptio /= 1 )   THEN
(  126)           IF(lwp) WRITE(numout,cform_err)
(  127)           IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on tracer'
(  128)           nstop = nstop + 1
(  129)       ENDIF
(  130)       ioptio = 0
(  131)       IF( ln_traldf_level )   ioptio = ioptio + 1
(  132)       IF( ln_traldf_hor   )   ioptio = ioptio + 1
(  133)       IF( ln_traldf_iso   )   ioptio = ioptio + 1
(  134)       IF( ioptio /= 1 ) THEN
(  135)          IF(lwp) WRITE(numout,cform_err)
(  136)          IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
(  137)          nstop = nstop + 1
(  138)       ENDIF
(  139) 
(  140)       ! ... Choice of the lateral scheme used
(  141)       IF( lk_traldf_eiv ) THEN
(  142)          IF(lwp) WRITE(numout,*) '          eddy induced velocity on tracers'
(  143)             IF( .NOT.ln_traldf_iso .OR. ln_traldf_bilap ) THEN
(  144)             IF(lwp) WRITE(numout,cform_err)
(  145)             IF(lwp) WRITE(numout,*) ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion'
(  146)             nstop = nstop + 1
(  147)          ENDIF
(  148)       ENDIF
(  149) 
(  150)       IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal mixing
(  151)          IF( ( ln_traldf_iso .OR. ln_traldf_hor ) .AND. .NOT.lk_ldfslp ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 10

(  152)             IF(lwp) WRITE(numout,cform_err)
(  153)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  154)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  155)          ENDIF
(  156)       ELSE                       ! z-coordinates with/without partial step:
(  157)          ln_traldf_level = ln_traldf_level .OR. ln_traldf_hor      ! level diffusion = horizontal diffusion
(  158)          ln_traldf_hor   = .FALSE.
(  159)          IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_traldf_level = T'
(  160)          IF(lwp) WRITE(numout,*) '                                                  and    force ln_traldf_hor   = F'
(  161)          IF( ln_traldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
(  162)             IF(lwp) WRITE(numout,cform_err)
(  163)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  164)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  165)          ENDIF
(  166)       ENDIF
(  167) 
(  168)       l_traldf_lap     =       ln_traldf_lap   .AND. ln_traldf_level     ! iso-level   laplacian operator
(  169)       l_traldf_bilap   =       ln_traldf_bilap .AND. ln_traldf_level     ! iso-level bilaplacian operator
(  170)       l_traldf_bilapg  =       ln_traldf_bilap .AND. ln_traldf_hor       ! geopotential bilap. (s-coord)
(  171)       l_traldf_iso     =       ln_traldf_lap   .AND.                  &  ! laplacian operator
(  172)          &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (z-coord) or horizontal (s-coord)
(  173)          &                                     .AND. .NOT.lk_zps
(  174)       l_traldf_iso_zps =       ln_traldf_lap   .AND.                  &  ! laplacian operator
(  175)          &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (partial steps)
(  176)          &                                     .AND. lk_zps              ! or geopotential in mixed partial steps/s-coord
(  177)       l_trazdf_iso    = .FALSE.
(  178)       l_trazdf_iso_vo = .FALSE.
(  179)       IF( l_traldf_iso     )   l_trazdf_iso = .TRUE.
(  180)       IF( l_traldf_iso_zps )   l_trazdf_iso = .TRUE.
(  181) 
(  182) # 188
(  188)       ioptio = 0
(  189)       IF( l_traldf_lap     )   ioptio = ioptio + 1
(  190)       IF( l_traldf_bilap   )   ioptio = ioptio + 1
(  191)       IF( l_traldf_bilapg  )   ioptio = ioptio + 1
(  192)       IF( l_traldf_iso     )   ioptio = ioptio + 1
(  193)       IF( l_traldf_iso_zps )   ioptio = ioptio + 1
(  194)       IF( ioptio /= 1 ) THEN
(  195)          IF(lwp) WRITE(numout,cform_err)
(  196)          IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
(  197)          nstop = nstop + 1
(  198)       ENDIF
(  199)       IF( lk_esopa ) THEN
(  200)          l_traldf_lap = .TRUE.   ;   l_traldf_bilap   = .TRUE.   ;   l_traldf_bilapg  = .TRUE.
(  201)          l_traldf_iso = .TRUE.   ;   l_traldf_iso_zps = .TRUE.
(  202)          l_trazdf_iso = .TRUE.   ;   l_trazdf_iso_vo  = .TRUE.
(  203)          IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
(  204)       ENDIF
(  205) 
(  206) 
(  207)       ! ... Space variation of eddy coefficients
(  208)       ioptio = 0
(  209) # 214
(  214)       IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
(  215)       ioptio = ioptio + 1
(  216) # 230
(  230)       IF( ioptio == 0 ) THEN
(  231)           IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 11

(  232)         ELSEIF( ioptio > 1 ) THEN
(  233)           IF(lwp) WRITE(numout,cform_err)
(  234)           IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
(  235)              &                    ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d'
(  236)           nstop = nstop + 1
(  237)       ENDIF
(  238) 
(  239)       IF( l_traldf_bilap .OR. l_traldf_bilapg ) THEN
(  240)          IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
(  241)          IF( aht0 > 0 .AND. .NOT. lk_esopa ) THEN
(  242)             IF(lwp) WRITE(numout,cform_err)
(  243)             IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be negative'
(  244)             nstop = nstop + 1
(  245)          ENDIF
(  246)       ELSE
(  247)          IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
(  248)          IF( aht0 < 0 .AND. .NOT. lk_esopa ) THEN
(  249)             IF(lwp) WRITE(numout,cform_err)
(  250)             IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be positive'
(  251)             nstop = nstop + 1
(  252)          ENDIF
(  253)       ENDIF
(  254) 
(  255) 
(  256)       !  Lateral eddy diffusivity and eddy induced velocity coefficients
(  257)       ! ================================================================
(  258) 
(  259) # 264
(  264)       CALL ldf_tra_smag( ll_print )           ! aht = 2D coef. = F( longitude, latitude )
(  265) 
(  266) # 281
(  281)    END SUBROUTINE ldf_tra_init
(  282) 
(  283) # 1 "./ldftra_smag.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                        ***  ldftra_smag.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4) 
(    5)    !!----------------------------------------------------------------------
(    6)    !!  OPA 9.0 , BIO (2005) 
(    7)    !!----------------------------------------------------------------------
(    8) 
(    9)    !!----------------------------------------------------------------------
(   10)    !!   'key_traldf_smag'             2D lateral eddy viscosity coefficients
(   11)    !!----------------------------------------------------------------------
(   12) 
(   13)    SUBROUTINE ldf_tra_smag( ld_print )
(   14)       !!----------------------------------------------------------------------
(   15)       !!                  ***  ROUTINE ldf_tra_smag  ***
(   16)       !!
(   17)       !! ** Purpose :   initializations of the horizontal eddy diffusivity coef.
(   18)       !!
(   19)       !! ** Method  :  SMAG
(   20)       !!       biharmonic operator  : ahtt (T)
(   21)       !!       harmonic operator    :
(   22)       !!                iso-model level   :ahtu,ahtv (U,V)
(   23)       !!                isopycnal         :ahtu,ahtv,ahtw (U,V,W)
(   24)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 12

(   25)       !!       eddy induced velocity
(   26)       !!         always harmonic   : aeiu, aeiv, aeiw defined at u-, v-, w-pts
(   27)       !!       History:
(   28)       !!              1.0 ! 07-05 Zeliang Wang ! Original Codes
(   29)       !!----------------------------------------------------------------------
(   30) 
(   31)       !!----------------------------------------------------------------------
(   32)       !! * Modules used
(   33)       !! * Arguments
(   34)       LOGICAL, INTENT (in) :: ld_print   ! If true, output arrays on numout
(   35) 
(   36)       !! * local variables
(   37)       INTEGER  ::   ji, jj, jk      ! dummy loop indices
(   38)       REAL(wp) ::   &
(   39)          alfa, bgvalue,dh_max,dh_min  ! SMAG coeff. 
(   40)       REAL(wp) ::   &
(   41)          dudx,dvdy,dvdx,dudy   ! SMAG coeff. 
(   42)       REAL(wp) ::   &
(   43)          dvdx1,dvdx2,dudy1,dudy2  ! SMAG coeff. 
(   44)       REAL(wp), DIMENSION(jpk) ::   zcoef       ! temporary workspace
(   45)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: tmph  ! temporary workspace
(   46)         alfa=0.1
(   47)       !!----------------------------------------------------------------------
(   48) 
(   49) 
(   50)       
(   51)       ! =================                       whatever its orientation is)
(   52)         DO ji=2,jpi-1
(   53)            DO jj=2,jpj-1
(   54)               DO jk=1,jpk
(   55)                  
(   56)                  dudx = (ub(ji,jj,jk)-ub(ji-1,jj,jk)) / e1t(ji,jj)*tmask(ji,jj,jk)
(   57)                  
(   58)                  dvdy = (vb(ji,jj,jk)-vb(ji,jj-1,jk)) / e2t(ji,jj)*tmask(ji,jj,jk)
(   59)                  
(   60)                  dudy1 = 0.5*((ub(ji,jj+1,jk)-ub(ji,jj,jk)) /  e2f(ji,jj) *fmask(ji,jj,jk)& 
(   61)                       &    + (ub(ji,jj,jk)-ub(ji,jj-1,jk)) / e2f(ji,jj-1)*fmask(ji,jj-1,jk))
(   62)                  
(   63)                  dudy2 = 0.5*((ub(ji-1,jj+1,jk)-ub(ji-1,jj,jk)) / e2f(ji-1,jj)*fmask(ji-1,jj,jk) &
(   64)                       &    + (ub(ji-1,jj,jk)-ub(ji-1,jj-1,jk)) / e2f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   65)                  
(   66)                  dvdx1 = 0.5*((vb(ji+1,jj,jk)-vb(ji,jj,jk)) / e1f(ji,jj)*fmask(ji,jj,jk) &
(   67)                       &    + (vb(ji,jj,jk)-vb(ji-1,jj,jk)) / e1f(ji-1,jj)*fmask(ji-1,jj,jk))
(   68)                  
(   69)                  dvdx2 = 0.5*((vb(ji+1,jj-1,jk)-vb(ji,jj-1,jk)) / e1f(ji,jj-1)*fmask(ji,jj-1,jk) &
(   70)                       &    + (vb(ji,jj-1,jk)-vb(ji-1,jj-1,jk)) / e1f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   71)                  
(   72)                  dudy=0.5*(dudy1+dudy2)
(   73)                  dvdx=0.5*(dvdx1+dvdx2)
(   74)                  
(   75)                  
(   76)                  tmph(ji,jj,jk)=alfa*e1t(ji,jj)*e2t(ji,jj)*SQRT(dudx**2+dvdy**2+0.5*(dvdx+dudy)**2)  
(   77)                  
(   78)               END DO
(   79)            END DO
(   80)         END DO
(   81)         
(   82)         CALL lbc_lnk(tmph,'T',1.)






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 13

(   83) 
(   84)         ahtw(:,:,:)=tmph(:,:,:)+aht0
(   85)         CALL lbc_lnk(ahtw,'W',1.)
(   86)         
(   87)         IF( ln_traldf_lap ) THEN
(   88)            
(   89)            DO ji=1,jpi-1
(   90)               DO jj=1,jpj          
(   91)                  DO jk=1,jpk
(   92)                     ahtu(ji,jj,jk)=(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji+1,jj,jk)*tmask(ji+1,jj,jk)+aht0)/  &
(   93)                          &         max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk),1.0)
(   94)                  END DO
(   95)               END DO
(   96)            END DO
(   97)            
(   98)            DO ji=1,jpi
(   99)               DO jj=1,jpj-1          
(  100)                  DO jk=1,jpk
(  101)                     ahtv(ji,jj,jk)=(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji,jj+1,jk)*tmask(ji,jj+1,jk)+aht0)/  &
(  102)                          &         max(tmask(ji,jj,jk)+tmask(ji,jj+1,jk),1.0)
(  103)                  END DO
(  104)               END DO
(  105)            END DO
(  106) 
(  107)            CALL lbc_lnk(ahtu,'U',1.)
(  108)            CALL lbc_lnk(ahtv,'V',1.)
(  109)            
(  110)            
(  111)            
(  112)         ENDIF
(  113)         
(  114)         
(  115)         
(  116)         IF( ln_traldf_bilap ) THEN
(  117)            
(  118)            DO jk=1,jpk
(  119)               ahtt(:,:,jk)=aht0-tmph(:,:,jk)*e1t(:,:)*e2t(:,:)/8.0
(  120)            END DO
(  121)            
(  122)            
(  123)         END IF
(  124)         
(  125)       END SUBROUTINE ldf_tra_smag
(  126) # 289 "ldftra.F90"
(  289) 
(  290) # 293
(  293)    !!======================================================================
(  294) END MODULE ldftra
