


PGF90 (Version     10.5)          09/23/2010  09:50:03      page 1

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

Filename: ldfdyn.F90

(    1) MODULE ldfdyn
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldfdyn  ***
(    4)    !! Ocean physics:  lateral viscosity coefficient 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   ldf_dyn_init : initialization, namelist read, and parameters control
(    9)    !!   ldf_dyn_c3d   : 3D eddy viscosity coefficient initialization
(   10)    !!   ldf_dyn_c2d   : 2D eddy viscosity coefficient initialization
(   11)    !!   ldf_dyn_c1d   : 1D eddy viscosity coefficient initialization
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers   
(   15)    USE dom_oce         ! ocean space and time domain 
(   16)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   17)    USE phycst          ! physical constants
(   18)    USE ldfslp          ! ???
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! *  Routine accessibility
(   27)    PUBLIC ldf_dyn_init   ! called by opa.F90
(   28) !!DB: due to new code changes must make the relevant routine PUBLIC or step() cannot find it
(   29) #if defined key_dynldf_smag || defined key_traldf_smag
(   30)    PUBLIC ldf_smag   ! possibly called by step.F90
(   31) #endif
(   32) #if defined key_dynldf_c3d
(   33)    PUBLIC ldf_dyn_c3d
(   34) #elif defined key_dynldf_c2d
(   35)    PUBLIC ldf_dyn_c2d
(   36) #elif defined key_dynldf_smag
(   37)    PUBLIC ldf_dyn_smag
(   38) #elif defined key_dynldf_c1d
(   39)    PUBLIC ldf_dyn_c1d
(   40) #else
(   41) !! do nothing
(   42) #endif
(   43) 
(   44) 
(   45) 
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 2

(   47)   INTERFACE ldf_zpf
(   48)      MODULE PROCEDURE ldf_zpf_1d, ldf_zpf_1d_3d, ldf_zpf_3d
(   49)   END INTERFACE
(   50) 
(   51)    !! * Substitutions
(   52) #  include "domzgr_substitute.h90"
(   53)    !!----------------------------------------------------------------------
(   54)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   55)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn.F90,v 1.5 2005/03/27 18:35:06 opalod Exp $ 
(   56)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   57)    !!----------------------------------------------------------------------
(   58) 
(   59) CONTAINS
(   60) 
(   61)    SUBROUTINE ldf_dyn_init
(   62)       !!----------------------------------------------------------------------
(   63)       !!                  ***  ROUTINE ldf_dyn_init  ***
(   64)       !!                   
(   65)       !! ** Purpose :   set the horizontal ocean dynamics physics
(   66)       !!
(   67)       !! ** Method  :  
(   68)       !!      Eddy viscosity coefficients:
(   69)       !!         default option   : constant coef. ahm0 (namelist)
(   70)       !!        'key_dynldf_c1d': depth dependent coef. defined in 
(   71)       !!                        in ldf_dyn_c1d routine
(   72)       !!        'key_dynldf_c2d': latitude and longitude dependent coef.
(   73)       !!                        defined in ldf_dyn_c2d routine
(   74)       !!        'key_dynldf_c3d': latitude, longitude, depth dependent coef.
(   75)       !!                        defined in ldf_dyn_c3d routine
(   76)       !!      N.B. User defined include files.  By default, 3d and 2d coef.
(   77)       !!      are set to a constant value given in the namelist and the 1d
(   78)       !!      coefficients are initialized to a hyperbolic tangent vertical
(   79)       !!      profile.
(   80)       !!
(   81)       !! Reference :
(   82)       !!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
(   83)       !!      the North Pole singularity, Climate Dynamics, 12, 381-388.
(   84)       !!
(   85)       !! History :
(   86)       !!        !  07-97  (G. Madec)  from inimix.F split in 2 routines
(   87)       !!        !  08-97  (G. Madec)  multi dimensional coefficients
(   88)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used
(   91)       USE ioipsl
(   92) 
(   93)       !! * Local declarations
(   94)       INTEGER ::   ioptio         ! ???
(   95)       LOGICAL :: ll_print = .FALSE.    ! Logical flag for printing viscosity coef.
(   96) 
(   97)        
(   98)       NAMELIST/nam_dynldf/ ln_dynldf_lap  , ln_dynldf_bilap,                &
(   99)          &                 ln_dynldf_level, ln_dynldf_hor, ln_dynldf_iso,   &
(  100)          &                 ahm0, ahmb0
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103) 
(  104)       ! Define the lateral physics parameters






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 3

(  105)       ! ======================================
(  106)     
(  107)       ! Read Namelist nam_dynldf : Lateral physics
(  108)       REWIND( numnam )
(  109)       READ  ( numnam, nam_dynldf )
(  110) 
(  111)       ! Parameter print
(  112)       IF(lwp) THEN
(  113)          WRITE(numout,*)
(  114)          WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
(  115)          WRITE(numout,*) '~~~~~~~'
(  116)          WRITE(numout,*) '          Namelist nam_dynldf : set lateral mixing parameters'
(  117)          WRITE(numout,*) '             laplacian operator          ln_dynldf_lap   = ', ln_dynldf_lap
(  118)          WRITE(numout,*) '             bilaplacian operator        ln_dynldf_bilap = ', ln_dynldf_bilap
(  119)          WRITE(numout,*) '             iso-level                   ln_dynldf_level = ', ln_dynldf_level
(  120)          WRITE(numout,*) '             horizontal (geopotential)   ln_dynldf_hor   = ', ln_dynldf_hor
(  121)          WRITE(numout,*) '             iso-neutral                 ln_dynldf_iso   = ', ln_dynldf_iso
(  122)          WRITE(numout,*) '             horizontal eddy viscosity            ahm0   = ', ahm0
(  123)          WRITE(numout,*) '             background viscosity                 ahmb0  = ', ahmb0
(  124)          WRITE(numout,*)
(  125)       ENDIF
(  126) 
(  127)       ! Parameter control
(  128) 
(  129)       ! control the input
(  130)       ioptio = 0
(  131)       IF( ln_dynldf_lap   )   ioptio = ioptio + 1
(  132)       IF( ln_dynldf_bilap )   ioptio = ioptio + 1
(  133)       IF( ioptio /= 1 )   THEN
(  134)           IF(lwp) WRITE(numout,cform_err)
(  135)           IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on momentum'
(  136)           nstop = nstop + 1
(  137)       ENDIF
(  138)       ioptio = 0
(  139)       IF( ln_dynldf_level )   ioptio = ioptio + 1
(  140)       IF( ln_dynldf_hor   )   ioptio = ioptio + 1
(  141)       IF( ln_dynldf_iso   )   ioptio = ioptio + 1
(  142)       IF( ioptio /= 1 ) THEN
(  143)          IF(lwp) WRITE(numout,cform_err)
(  144)          IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
(  145)          nstop = nstop + 1
(  146)       ENDIF
(  147) 
(  148)       IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal direction
(  149)          IF( ( ln_dynldf_iso .OR. ln_dynldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
(  150)             IF(lwp) WRITE(numout,cform_err)
(  151)             IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
(  152)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  153)          ENDIF
(  154)       ELSE                       ! z-coordinates with/without partial step:
(  155)          ln_dynldf_level = ln_dynldf_level .OR. ln_dynldf_hor      ! level mixing = horizontal mixing
(  156)          ln_dynldf_hor   = .FALSE.
(  157)          IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_dynldf_level = T'
(  158)          IF(lwp) WRITE(numout,*) '                                                  and    force ln_dynldf_hor   = F'
(  159)          IF( ln_dynldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
(  160)             IF(lwp) WRITE(numout,cform_err)
(  161)             IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
(  162)             IF( .NOT.lk_esopa )   nstop = nstop + 1






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 4

(  163)          ENDIF
(  164)       ENDIF
(  165) 
(  166)       l_dynldf_lap     =       ln_dynldf_lap   .AND. ln_dynldf_level     ! iso-level   laplacian operator
(  167)       l_dynldf_bilap   =       ln_dynldf_bilap .AND. ln_dynldf_level     ! iso-level bilaplacian operator
(  168)       l_dynldf_bilapg  =       ln_dynldf_bilap .AND. ln_dynldf_hor       ! geopotential bilap. (s-coord)
(  169)       l_dynldf_iso     =       ln_dynldf_lap   .AND.                  &  ! laplacian operator
(  170)          &                   ( ln_dynldf_iso   .OR.  ln_dynldf_hor )     ! iso-neutral (z-coord) or horizontal (s-coord)
(  171) 
(  172)       l_dynzdf_iso    = .FALSE.
(  173)       IF( l_dynldf_iso )   l_dynzdf_iso = .TRUE.
(  174) 
(  175)       ioptio = 0
(  176)       IF( l_dynldf_lap     )   ioptio = ioptio + 1
(  177)       IF( l_dynldf_bilap   )   ioptio = ioptio + 1
(  178)       IF( l_dynldf_bilapg  )   ioptio = ioptio + 1
(  179)       IF( l_dynldf_iso     )   ioptio = ioptio + 1
(  180)       IF( ioptio /= 1 ) THEN
(  181)          IF(lwp) WRITE(numout,cform_err)
(  182)          IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
(  183)          nstop = nstop + 1
(  184)       ENDIF
(  185)       IF( lk_esopa ) THEN
(  186)          l_dynldf_lap = .TRUE.   ;   l_dynldf_bilap   = .TRUE.   ;   l_dynldf_bilapg  = .TRUE.
(  187)          l_dynldf_iso = .TRUE.   ;   l_dynzdf_iso     = .TRUE.
(  188)          IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
(  189)       ENDIF
(  190) 
(  191)       ! control print
(  192)       IF( l_dynldf_lap    .AND. lwp ) WRITE(numout,*) '          iso-level laplacian momentum operator'
(  193)       IF( l_dynldf_bilap  .AND. lwp ) WRITE(numout,*) '          iso-level bilaplacian momentum operator'
(  194)       IF( l_dynldf_bilapg .AND. lwp ) WRITE(numout,*) '          geopotential bilaplacian momentum operator'
(  195)       IF( l_dynldf_iso    .AND. lwp ) WRITE(numout,*) '          iso-neutral laplacian momentum operator'
(  196) 
(  197)       ! ... Space variation of eddy coefficients
(  198)       ioptio = 0
(  199) #if defined key_dynldf_c3d
(  200)       IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth)'
(  201)       ioptio = ioptio+1
(  202) #endif
(  203) #if defined key_dynldf_smag
(  204)       IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth)'
(  205)       ioptio = ioptio+1
(  206) #endif
(  207) #if defined key_dynldf_c2d
(  208)       IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude)'
(  209)       ioptio = ioptio+1
(  210) #endif
(  211) #if defined key_dynldf_c1d
(  212)       IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( depth )'
(  213)       ioptio = ioptio+1
(  214)       IF( lk_sco ) THEN
(  215)          IF(lwp) WRITE(numout,cform_err)
(  216)          IF(lwp) WRITE(numout,*) '          key_dynldf_c1d cannot be used in s-coordinate (key_s_coord)'
(  217)          nstop = nstop + 1
(  218)       ENDIF
(  219) #endif
(  220)       IF( ioptio == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 5

(  221)           IF(lwp) WRITE(numout,*) '          momentum mixing coef. = constant  (default option)'
(  222)         ELSEIF( ioptio > 1 ) THEN
(  223)           IF(lwp) WRITE(numout,cform_err)
(  224)           IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
(  225)                                   ' key_dynldf_c3d, key_dynldf_c2d, key_dynldf_c1d'
(  226)           nstop = nstop + 1
(  227)       ENDIF
(  228) 
(  229) 
(  230)       IF( l_dynldf_bilap .OR. l_dynldf_bilapg ) THEN
(  231)          IF(lwp) WRITE(numout,*) '          biharmonic momentum diffusion'
(  232)          IF( ahm0 > 0 .AND. .NOT. lk_esopa ) THEN
(  233)             IF(lwp) WRITE(numout,cform_err)
(  234)             IF(lwp) WRITE(numout,*) 'The horizontal viscosity coef. ahm0 must be negative'
(  235)             nstop = nstop + 1
(  236)          ENDIF
(  237)       ELSE
(  238)          IF(lwp) WRITE(numout,*) '          harmonic momentum diff. (default)'
(  239)          IF( ahm0 < 0 .AND. .NOT. lk_esopa ) THEN
(  240)             IF(lwp) WRITE(numout,cform_err)
(  241)             IF(lwp) WRITE(numout,*) '          The horizontal viscosity coef. ahm0 must be positive'
(  242)             nstop = nstop + 1
(  243)          ENDIF
(  244)       ENDIF
(  245) 
(  246) 
(  247)       ! Lateral eddy viscosity
(  248)       ! ======================
(  249) 
(  250) #if defined key_dynldf_c3d
(  251)       CALL ldf_dyn_c3d( ll_print )   ! ahm = 3D coef. = F( longitude, latitude, depth )
(  252) #elif defined key_dynldf_c2d
(  253)       CALL ldf_dyn_c2d( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )
(  254) #elif defined key_dynldf_smag
(  255)       CALL ldf_dyn_smag( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )
(  256) #elif defined key_dynldf_c1d
(  257)       CALL ldf_dyn_c1d( ll_print )   ! ahm = 1D coef. = F( depth )
(  258) #else
(  259)                                      ! Constant coefficients
(  260)       IF(lwp) WRITE(numout,*)
(  261)       IF(lwp) WRITE(numout,*) 'inildf: constant eddy viscosity coef. '
(  262)       IF(lwp) WRITE(numout,*) '~~~~~~'
(  263)       IF(lwp) WRITE(numout,*) '        ahm1 = ahm2 = ahm0 =  ',ahm0
(  264) #endif
(  265) 
(  266)    END SUBROUTINE ldf_dyn_init
(  267) 
(  268) #if defined key_dynldf_c3d
(  269) #   include "ldfdyn_c3d.h90"
(  270) #elif defined key_dynldf_c2d
(  271) #   include "ldfdyn_c2d.h90"
(  272) #elif defined key_dynldf_smag
(  273) #   include "ldfdyn_smag.h90"
(  274) #elif defined key_dynldf_c1d
(  275) #   include "ldfdyn_c1d.h90"
(  276) #endif
(  277) 
(  278) !!DB: 2008.04.09 -- changes to smag procedures






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 6

(  279) #if defined key_dynldf_smag || defined key_traldf_smag
(  280) #   include "ldfsmag.h90"
(  281) #endif 
(  282) 
(  283)    SUBROUTINE ldf_zpf_1d( ld_print, pdam, pwam, pbot, pdep, pah )
(  284)       !!----------------------------------------------------------------------
(  285)       !!                  ***  ROUTINE ldf_zpf  ***
(  286)       !!                   
(  287)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  288)       !!
(  289)       !! ** Method  :   1D eddy viscosity coefficients ( depth )
(  290)       !!
(  291)       !!----------------------------------------------------------------------
(  292)       !! * Arguments
(  293)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  294)       REAL(wp), INTENT (in   ) ::   &
(  295)           pdam,     &  ! depth of the inflection point
(  296)           pwam,     &  ! width of inflection
(  297)           pbot         ! battom value (0<pbot<= 1)
(  298)       REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
(  299)           pdep         ! depth of the gridpoint (T, U, V, F)
(  300)       REAL(wp), INTENT (inout), DIMENSION(jpk) ::   &
(  301)           pah          ! adimensional vertical profile
(  302) 
(  303)       !! * Local variables
(  304)       INTEGER  ::   jk           ! dummy loop indices
(  305)       REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
(  306)       !!----------------------------------------------------------------------
(  307) 
(  308)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
(  309)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  310)       zmhs = zm00 / zm01
(  311)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  312) 
(  313)       DO jk = 1, jpk
(  314)          pah(jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
(  315)       END DO
(  316) 
(  317)       ! Control print
(  318)       IF(lwp .AND. ld_print ) THEN
(  319)          WRITE(numout,*)
(  320)          WRITE(numout,*) '         ahm profile : '
(  321)          WRITE(numout,*)
(  322)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  323)          DO jk = 1, jpk
(  324)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(jk), pdep(jk)
(  325)          END DO
(  326)       ENDIF
(  327) 
(  328)    END SUBROUTINE ldf_zpf_1d
(  329) 
(  330) 
(  331)    SUBROUTINE ldf_zpf_1d_3d( ld_print, pdam, pwam, pbot, pdep, pah )
(  332)       !!----------------------------------------------------------------------
(  333)       !!                  ***  ROUTINE ldf_zpf  ***
(  334)       !!
(  335)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  336)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 7

(  337)       !! ** Method  :   1D eddy viscosity coefficients ( depth )
(  338)       !!
(  339)       !!----------------------------------------------------------------------
(  340)       !! * Arguments
(  341)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  342)       REAL(wp), INTENT (in   ) ::   &
(  343)           pdam,     &  ! depth of the inflection point
(  344)           pwam,     &  ! width of inflection
(  345)           pbot         ! battom value (0<pbot<= 1)
(  346)       REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
(  347)           pdep         ! depth of the gridpoint (T, U, V, F)
(  348)       REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
(  349)           pah          ! adimensional vertical profile
(  350) 
(  351)       !! * Local variables
(  352)       INTEGER  ::   jk           ! dummy loop indices
(  353)       REAL(wp) ::   zm00, zm01, zmhb, zmhs, zcf  ! temporary scalars
(  354)       !!----------------------------------------------------------------------
(  355) 
(  356)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
(  357)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  358)       zmhs = zm00 / zm01
(  359)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  360) 
(  361)       DO jk = 1, jpk
(  362)          zcf = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
(  363)          pah(:,:,jk) = zcf
(  364)       END DO
(  365) 
(  366)       ! Control print
(  367)       IF(lwp .AND. ld_print ) THEN
(  368)          WRITE(numout,*)
(  369)          WRITE(numout,*) '         ahm profile : '
(  370)          WRITE(numout,*)
(  371)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  372)          DO jk = 1, jpk
(  373)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(jk)
(  374)          END DO
(  375)       ENDIF
(  376) 
(  377)    END SUBROUTINE ldf_zpf_1d_3d
(  378) 
(  379) 
(  380)    SUBROUTINE ldf_zpf_3d( ld_print, pdam, pwam, pbot, pdep, pah )
(  381)       !!----------------------------------------------------------------------
(  382)       !!                  ***  ROUTINE ldf_zpf  ***
(  383)       !!
(  384)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  385)       !!
(  386)       !! ** Method  :   3D for partial step or s-coordinate
(  387)       !!
(  388)       !!----------------------------------------------------------------------
(  389)       !! * Arguments
(  390)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  391)       REAL(wp), INTENT (in   ) ::   &
(  392)           pdam,     &  ! depth of the inflection point
(  393)           pwam,     &  ! width of inflection
(  394)           pbot         ! reduction factor (surface value / bottom value)






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 8

(  395)       REAL(wp), INTENT (in   ), DIMENSION(jpi,jpj,jpk) ::   &
(  396)           pdep         ! dep of the gridpoint (T, U, V, F)
(  397)       REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
(  398)           pah          ! adimensional vertical profile
(  399) 
(  400)       !! * Local variables
(  401)       INTEGER  ::   jk           ! dummy loop indices
(  402)       REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
(  403)       !!----------------------------------------------------------------------
(  404) 
(  405)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )   
(  406)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  407)       zmhs = zm00 / zm01
(  408)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  409) 
(  410)       DO jk = 1, jpk
(  411)          pah(:,:,jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(:,:,jk) ) / pwam )  )
(  412)       END DO
(  413) 
(  414)       ! Control print
(  415)       IF(lwp .AND. ld_print ) THEN
(  416)          WRITE(numout,*)
(  417)          WRITE(numout,*) '         ahm profile : '
(  418)          WRITE(numout,*)
(  419)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  420)          DO jk = 1, jpk
(  421)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(1,1,jk)
(  422)          END DO
(  423)       ENDIF
(  424) 
(  425)    END SUBROUTINE ldf_zpf_3d
(  426)    !!======================================================================
(  427) END MODULE ldfdyn































PGF90 (Version     10.5)          09/23/2010  09:50:03      page 9

(    1) # 1 "ldfdyn.F90"
(    1) MODULE ldfdyn
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldfdyn  ***
(    4)    !! Ocean physics:  lateral viscosity coefficient 
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   ldf_dyn_init : initialization, namelist read, and parameters control
(    9)    !!   ldf_dyn_c3d   : 3D eddy viscosity coefficient initialization
(   10)    !!   ldf_dyn_c2d   : 2D eddy viscosity coefficient initialization
(   11)    !!   ldf_dyn_c1d   : 1D eddy viscosity coefficient initialization
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers   
(   15)    USE dom_oce         ! ocean space and time domain 
(   16)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   17)    USE phycst          ! physical constants
(   18)    USE ldfslp          ! ???
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lib_mpp         ! distribued memory computing library
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! *  Routine accessibility
(   27)    PUBLIC ldf_dyn_init   ! called by opa.F90
(   28) !!DB: due to new code changes must make the relevant routine PUBLIC or step() cannot find it
(   29) # 30
(   30)    PUBLIC ldf_smag   ! possibly called by step.F90
(   31) # 37
(   37)    PUBLIC ldf_dyn_smag
(   38) 
(   39) 
(   40) 
(   41) 
(   42) # 47
(   47)   INTERFACE ldf_zpf
(   48)      MODULE PROCEDURE ldf_zpf_1d, ldf_zpf_1d_3d, ldf_zpf_3d
(   49)   END INTERFACE
(   50) 
(   51)    !! * Substitutions
(   52) # 1 "./domzgr_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 10

(   48) # 53 "ldfdyn.F90"
(   53) # 53
(   53)    !!----------------------------------------------------------------------
(   54)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   55)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn.F90,v 1.5 2005/03/27 18:35:06 opalod Exp $ 
(   56)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   57)    !!----------------------------------------------------------------------
(   58) 
(   59) CONTAINS
(   60) 
(   61)    SUBROUTINE ldf_dyn_init
(   62)       !!----------------------------------------------------------------------
(   63)       !!                  ***  ROUTINE ldf_dyn_init  ***
(   64)       !!                   
(   65)       !! ** Purpose :   set the horizontal ocean dynamics physics
(   66)       !!
(   67)       !! ** Method  :  
(   68)       !!      Eddy viscosity coefficients:
(   69)       !!         default option   : constant coef. ahm0 (namelist)
(   70)       !!        'key_dynldf_c1d': depth dependent coef. defined in 
(   71)       !!                        in ldf_dyn_c1d routine
(   72)       !!        'key_dynldf_c2d': latitude and longitude dependent coef.
(   73)       !!                        defined in ldf_dyn_c2d routine
(   74)       !!        'key_dynldf_c3d': latitude, longitude, depth dependent coef.
(   75)       !!                        defined in ldf_dyn_c3d routine
(   76)       !!      N.B. User defined include files.  By default, 3d and 2d coef.
(   77)       !!      are set to a constant value given in the namelist and the 1d
(   78)       !!      coefficients are initialized to a hyperbolic tangent vertical
(   79)       !!      profile.
(   80)       !!
(   81)       !! Reference :
(   82)       !!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
(   83)       !!      the North Pole singularity, Climate Dynamics, 12, 381-388.
(   84)       !!
(   85)       !! History :
(   86)       !!        !  07-97  (G. Madec)  from inimix.F split in 2 routines
(   87)       !!        !  08-97  (G. Madec)  multi dimensional coefficients
(   88)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used
(   91)       USE ioipsl
(   92) 
(   93)       !! * Local declarations
(   94)       INTEGER ::   ioptio         ! ???
(   95)       LOGICAL :: ll_print = .FALSE.    ! Logical flag for printing viscosity coef.
(   96) 
(   97)        
(   98)       NAMELIST/nam_dynldf/ ln_dynldf_lap  , ln_dynldf_bilap,                &
(   99)          &                 ln_dynldf_level, ln_dynldf_hor, ln_dynldf_iso,   &
(  100)          &                 ahm0, ahmb0
(  101)       !!----------------------------------------------------------------------
(  102) 
(  103) 
(  104)       ! Define the lateral physics parameters
(  105)       ! ======================================
(  106)     
(  107)       ! Read Namelist nam_dynldf : Lateral physics
(  108)       REWIND( numnam )






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 11

(  109)       READ  ( numnam, nam_dynldf )
(  110) 
(  111)       ! Parameter print
(  112)       IF(lwp) THEN
(  113)          WRITE(numout,*)
(  114)          WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
(  115)          WRITE(numout,*) '~~~~~~~'
(  116)          WRITE(numout,*) '          Namelist nam_dynldf : set lateral mixing parameters'
(  117)          WRITE(numout,*) '             laplacian operator          ln_dynldf_lap   = ', ln_dynldf_lap
(  118)          WRITE(numout,*) '             bilaplacian operator        ln_dynldf_bilap = ', ln_dynldf_bilap
(  119)          WRITE(numout,*) '             iso-level                   ln_dynldf_level = ', ln_dynldf_level
(  120)          WRITE(numout,*) '             horizontal (geopotential)   ln_dynldf_hor   = ', ln_dynldf_hor
(  121)          WRITE(numout,*) '             iso-neutral                 ln_dynldf_iso   = ', ln_dynldf_iso
(  122)          WRITE(numout,*) '             horizontal eddy viscosity            ahm0   = ', ahm0
(  123)          WRITE(numout,*) '             background viscosity                 ahmb0  = ', ahmb0
(  124)          WRITE(numout,*)
(  125)       ENDIF
(  126) 
(  127)       ! Parameter control
(  128) 
(  129)       ! control the input
(  130)       ioptio = 0
(  131)       IF( ln_dynldf_lap   )   ioptio = ioptio + 1
(  132)       IF( ln_dynldf_bilap )   ioptio = ioptio + 1
(  133)       IF( ioptio /= 1 )   THEN
(  134)           IF(lwp) WRITE(numout,cform_err)
(  135)           IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on momentum'
(  136)           nstop = nstop + 1
(  137)       ENDIF
(  138)       ioptio = 0
(  139)       IF( ln_dynldf_level )   ioptio = ioptio + 1
(  140)       IF( ln_dynldf_hor   )   ioptio = ioptio + 1
(  141)       IF( ln_dynldf_iso   )   ioptio = ioptio + 1
(  142)       IF( ioptio /= 1 ) THEN
(  143)          IF(lwp) WRITE(numout,cform_err)
(  144)          IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
(  145)          nstop = nstop + 1
(  146)       ENDIF
(  147) 
(  148)       IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal direction
(  149)          IF( ( ln_dynldf_iso .OR. ln_dynldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
(  150)             IF(lwp) WRITE(numout,cform_err)
(  151)             IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
(  152)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  153)          ENDIF
(  154)       ELSE                       ! z-coordinates with/without partial step:
(  155)          ln_dynldf_level = ln_dynldf_level .OR. ln_dynldf_hor      ! level mixing = horizontal mixing
(  156)          ln_dynldf_hor   = .FALSE.
(  157)          IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_dynldf_level = T'
(  158)          IF(lwp) WRITE(numout,*) '                                                  and    force ln_dynldf_hor   = F'
(  159)          IF( ln_dynldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
(  160)             IF(lwp) WRITE(numout,cform_err)
(  161)             IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
(  162)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  163)          ENDIF
(  164)       ENDIF
(  165) 
(  166)       l_dynldf_lap     =       ln_dynldf_lap   .AND. ln_dynldf_level     ! iso-level   laplacian operator






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 12

(  167)       l_dynldf_bilap   =       ln_dynldf_bilap .AND. ln_dynldf_level     ! iso-level bilaplacian operator
(  168)       l_dynldf_bilapg  =       ln_dynldf_bilap .AND. ln_dynldf_hor       ! geopotential bilap. (s-coord)
(  169)       l_dynldf_iso     =       ln_dynldf_lap   .AND.                  &  ! laplacian operator
(  170)          &                   ( ln_dynldf_iso   .OR.  ln_dynldf_hor )     ! iso-neutral (z-coord) or horizontal (s-coord)
(  171) 
(  172)       l_dynzdf_iso    = .FALSE.
(  173)       IF( l_dynldf_iso )   l_dynzdf_iso = .TRUE.
(  174) 
(  175)       ioptio = 0
(  176)       IF( l_dynldf_lap     )   ioptio = ioptio + 1
(  177)       IF( l_dynldf_bilap   )   ioptio = ioptio + 1
(  178)       IF( l_dynldf_bilapg  )   ioptio = ioptio + 1
(  179)       IF( l_dynldf_iso     )   ioptio = ioptio + 1
(  180)       IF( ioptio /= 1 ) THEN
(  181)          IF(lwp) WRITE(numout,cform_err)
(  182)          IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
(  183)          nstop = nstop + 1
(  184)       ENDIF
(  185)       IF( lk_esopa ) THEN
(  186)          l_dynldf_lap = .TRUE.   ;   l_dynldf_bilap   = .TRUE.   ;   l_dynldf_bilapg  = .TRUE.
(  187)          l_dynldf_iso = .TRUE.   ;   l_dynzdf_iso     = .TRUE.
(  188)          IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
(  189)       ENDIF
(  190) 
(  191)       ! control print
(  192)       IF( l_dynldf_lap    .AND. lwp ) WRITE(numout,*) '          iso-level laplacian momentum operator'
(  193)       IF( l_dynldf_bilap  .AND. lwp ) WRITE(numout,*) '          iso-level bilaplacian momentum operator'
(  194)       IF( l_dynldf_bilapg .AND. lwp ) WRITE(numout,*) '          geopotential bilaplacian momentum operator'
(  195)       IF( l_dynldf_iso    .AND. lwp ) WRITE(numout,*) '          iso-neutral laplacian momentum operator'
(  196) 
(  197)       ! ... Space variation of eddy coefficients
(  198)       ioptio = 0
(  199) # 204
(  204)       IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth)'
(  205)       ioptio = ioptio+1
(  206) # 220
(  220)       IF( ioptio == 0 ) THEN
(  221)           IF(lwp) WRITE(numout,*) '          momentum mixing coef. = constant  (default option)'
(  222)         ELSEIF( ioptio > 1 ) THEN
(  223)           IF(lwp) WRITE(numout,cform_err)
(  224)           IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
(  225)                                   ' key_dynldf_c3d, key_dynldf_c2d, key_dynldf_c1d'
(  226)           nstop = nstop + 1
(  227)       ENDIF
(  228) 
(  229) 
(  230)       IF( l_dynldf_bilap .OR. l_dynldf_bilapg ) THEN
(  231)          IF(lwp) WRITE(numout,*) '          biharmonic momentum diffusion'
(  232)          IF( ahm0 > 0 .AND. .NOT. lk_esopa ) THEN
(  233)             IF(lwp) WRITE(numout,cform_err)
(  234)             IF(lwp) WRITE(numout,*) 'The horizontal viscosity coef. ahm0 must be negative'
(  235)             nstop = nstop + 1
(  236)          ENDIF
(  237)       ELSE
(  238)          IF(lwp) WRITE(numout,*) '          harmonic momentum diff. (default)'
(  239)          IF( ahm0 < 0 .AND. .NOT. lk_esopa ) THEN
(  240)             IF(lwp) WRITE(numout,cform_err)
(  241)             IF(lwp) WRITE(numout,*) '          The horizontal viscosity coef. ahm0 must be positive'






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 13

(  242)             nstop = nstop + 1
(  243)          ENDIF
(  244)       ENDIF
(  245) 
(  246) 
(  247)       ! Lateral eddy viscosity
(  248)       ! ======================
(  249) 
(  250) # 255
(  255)       CALL ldf_dyn_smag( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )
(  256) 
(  257) # 266
(  266)    END SUBROUTINE ldf_dyn_init
(  267) 
(  268) # 1 "./ldfdyn_smag.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                        ***  ldfdyn_smag.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4) 
(    5)    !!----------------------------------------------------------------------
(    6)    !!  OPA 9.0 , BIO (2005) 
(    7)    !!----------------------------------------------------------------------
(    8) 
(    9)    !!----------------------------------------------------------------------
(   10)    !!   'key_dynldf_smag'             2D lateral eddy viscosity coefficients
(   11)    !!----------------------------------------------------------------------
(   12) 
(   13)    SUBROUTINE ldf_dyn_smag( ld_print )
(   14)       !!----------------------------------------------------------------------
(   15)       !!                  ***  ROUTINE ldf_dyn_smag  ***
(   16)       !!                   
(   17)       !! ** Purpose :   initializations of the horizontal ocean physics
(   18)       !!
(   19)       !! ** Method  :   3D eddy viscosity coef. ( longitude, latitude, depth )
(   20)       !!       laplacian operator   : ahm1, ahm2 defined at T- and F-points
(   21)       !!                              ahm2, ahm4 never used
(   22)       !!       bilaplacian operator : ahm1, ahm2 never used
(   23)       !!                           :  ahm3, ahm4 defined at U- and V-points
(   24)       !!       ??? explanation of the default is missing
(   25)       !!       History: 
(   26)       !!              1.0 ! 07-05 Zeliang Wang ! Original Codes 
(   27)       !!----------------------------------------------------------------------
(   28)       !! * Modules used
(   29)       !! * Arguments
(   30)       LOGICAL, INTENT (in) :: ld_print   ! If true, output arrays on numout
(   31) 
(   32)       !! * local variables
(   33)       INTEGER  ::   ji, jj, jk      ! dummy loop indices
(   34)       REAL(wp) ::   &
(   35)          alfa, bgvalue,dh_max,dh_min  ! SMAG coeff. 
(   36)       REAL(wp) ::   &
(   37)          dudx,dvdy,dvdx,dudy   ! SMAG coeff. 
(   38)       REAL(wp) ::   &
(   39)          dvdx1,dvdx2,dudy1,dudy2  ! SMAG coeff. 
(   40)       REAL(wp), DIMENSION(jpk) ::   zcoef       ! temporary workspace
(   41)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: tmph  ! temporary workspace
(   42)         alfa=0.1
(   43) 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 14

(   44)    !   IF(lwp) WRITE(numout,*)
(   45)    !   IF(lwp) WRITE(numout,*) 'ldf_dyn_smag : SMAG lateral eddy viscosity coefficient'
(   46)    !   IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(   47) 
(   48)       
(   49)       ! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operators
(   50)       ! =================                       whatever its orientation is)
(   51)          DO ji=2,jpi-1
(   52)             DO jj=2,jpj-1
(   53)                DO jk=1,jpk
(   54) 
(   55)            dudx = (ub(ji,jj,jk)-ub(ji-1,jj,jk)) / e1t(ji,jj)*tmask(ji,jj,jk)
(   56) 
(   57)            dvdy = (vb(ji,jj,jk)-vb(ji,jj-1,jk)) / e2t(ji,jj)*tmask(ji,jj,jk)
(   58) 
(   59)            dudy1 = 0.5*((ub(ji,jj+1,jk)-ub(ji,jj,jk)) /  e2f(ji,jj) *fmask(ji,jj,jk)& 
(   60)                  &    + (ub(ji,jj,jk)-ub(ji,jj-1,jk)) / e2f(ji,jj-1)*fmask(ji,jj-1,jk))
(   61) 
(   62)            dudy2 = 0.5*((ub(ji-1,jj+1,jk)-ub(ji-1,jj,jk)) / e2f(ji-1,jj)*fmask(ji-1,jj,jk) &
(   63)                  &    + (ub(ji-1,jj,jk)-ub(ji-1,jj-1,jk)) / e2f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   64) 
(   65)            dvdx1 = 0.5*((vb(ji+1,jj,jk)-vb(ji,jj,jk)) / e1f(ji,jj)*fmask(ji,jj,jk) &
(   66)                  &    + (vb(ji,jj,jk)-vb(ji-1,jj,jk)) / e1f(ji-1,jj)*fmask(ji-1,jj,jk))
(   67) 
(   68)            dvdx2 = 0.5*((vb(ji+1,jj-1,jk)-vb(ji,jj-1,jk)) / e1f(ji,jj-1)*fmask(ji,jj-1,jk) &
(   69)                  &    + (vb(ji,jj-1,jk)-vb(ji-1,jj-1,jk)) / e1f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   70) 
(   71)            dudy=0.5*(dudy1+dudy2)
(   72)            dvdx=0.5*(dvdx1+dvdx2)
(   73)             
(   74)   
(   75)            tmph(ji,jj,jk)=alfa*e1t(ji,jj)*e2t(ji,jj)*SQRT(dudx**2+dvdy**2+0.5*(dvdx+dudy)**2)  
(   76) 
(   77)               END DO
(   78)             END DO
(   79)          END DO
(   80) 
(   81)           CALL lbc_lnk(tmph,'T',1.)
(   82) 
(   83)       IF( ln_dynldf_lap ) THEN
(   84)          ! define ahm1 and ahm2 at the right grid point position
(   85)          ! (USER: modify ahm1 and ahm2 following your desiderata)
(   86) 
(   87) 
(   88)          IF(lwp) WRITE(numout,*) '              laplacian operator: ahm calculated by SMAG'
(   89)          IF(lwp) WRITE(numout,*) '              alfa can be changed!                      '    
(   90)          IF(lwp) WRITE(numout,*) '              test case!                                ' 
(   91) 
(   92) 
(   93) 
(   94)          IF( ln_dynldf_iso ) THEN
(   95)             IF(lwp) WRITE(numout,*) '              Caution, as implemented now, the isopycnal part of momentum'
(   96)             IF(lwp) WRITE(numout,*) '                 mixing use aht0 as eddy viscosity coefficient. Thus, it is'
(   97)             IF(lwp) WRITE(numout,*) '                 uniform and you must be sure that your ahm is greater than'
(   98)             IF(lwp) WRITE(numout,*) '                 aht0 everywhere in the model domain.'
(   99)          ENDIF
(  100) 
(  101)          ahm1=ahm0+tmph






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 15

(  102)          CALL lbc_lnk(ahm1,'T',1.)
(  103)          DO ji=1,jpi-1
(  104)             DO jj=1,jpj-1
(  105)                
(  106)                DO jk=1,jpk
(  107)                   ahm2(ji,jj,jk)=(ahm1(ji,jj,jk)*tmask(ji,jj,jk)+ahm1(ji+1,jj,jk)*tmask(ji+1,jj,jk)+      &
(  108)                        &             ahm1(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)+ahm1(ji,jj+1,jk)*tmask(ji,jj+1,jk))/   &
(  109)                        &             max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk)+tmask(ji+1,jj+1,jk)+tmask(ji,jj+1,jk),1.0)
(  110)                END DO
(  111)                
(  112)             END DO
(  113)          END DO
(  114) 
(  115)          CALL lbc_lnk(ahm2,'F',1.)
(  116)          
(  117)          
(  118)          dh_max=MAXVAL(ahm1(:,:,:))
(  119)          IF(lk_mpp) CALL mpp_max(dh_max)
(  120)          dh_min=MINVAL(ahm1(:,:,:))
(  121)          IF(lk_mpp) CALL mpp_min(dh_min)
(  122)          
(  123)          IF(lwp) WRITE(numout,*) '              laplacian operator: ahm calculated by SMAG'
(  124)          IF(lwp) WRITE(numout,*) ' MAX_DH',dh_max,'MIN_DH',dh_min
(  125)          
(  126)          
(  127)          
(  128)       ENDIF
(  129)       
(  130)       
(  131) 
(  132)       ! ahm3 and ahm4 at U- and V-points (used for bilaplacian operator
(  133)       ! ================================  whatever its orientation is)
(  134)       ! (USER: modify ahm3 and ahm4 following your desiderata)
(  135)       ! Here: ahm is proportional to the cube of the maximum of the gridspacing
(  136)       !       in the to horizontal direction
(  137)       !  Warning SHOUlD BE PUT AT U,V POINT!!!!!!!!
(  138)       IF( ln_dynldf_bilap ) THEN
(  139) 
(  140) 
(  141) !        IF(lwp) WRITE(numout,*) '              SAMG , but is going to be soon!! '
(  142)    !      IF(lwp) WRITE(numout,*) '              bi-laplacian operator: ahm proportional to e1**3 '
(  143)    !      IF(lwp) WRITE(numout,*) '              Caution, here we assume your mesh is isotropic ...'
(  144) 
(  145)         DO ji=1,jpi-1
(  146)            DO jj=1,jpj-1
(  147)               DO jk=1,jpk 
(  148)          ahm3(ji,jj,jk) =ahm0-(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji+1,jj,jk)*tmask(ji+1,jj,jk))* &
(  149)            &             min(e1u(ji,jj),e2u(ji,jj))**2/8. &
(  150)            &             /max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk),1.0) 
(  151)          ahm4(ji,jj,jk) =ahm0-(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji,jj+1,jk)*tmask(ji,jj+1,jk))* &
(  152)            &             min(e1v(ji,jj),e2v(ji,jj))**2/8. &
(  153)            &             /max(tmask(ji,jj,jk)+tmask(ji,jj+1,jk),1.0) 
(  154)               END DO
(  155)            END DO
(  156)         END DO 
(  157)    !    dh_max=MAXVAL(ABS(ahm3(:,:,:)))
(  158)    !    IF(lk_mpp) CALL mpp_max(dh_max)
(  159)    !    IF(lwp) WRITE(numout,*) 'Maximum Biharmonic Coef. is ', dh_max 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 16

(  160) 
(  161) 
(  162) 
(  163)           CALL lbc_lnk(ahm3,'U',1.)
(  164)           CALL lbc_lnk(ahm4,'V',1.)
(  165)       END IF
(  166) 
(  167)    END SUBROUTINE ldf_dyn_smag
(  168) # 274 "ldfdyn.F90"
(  274) 
(  275) # 278
(  278) !!DB: 2008.04.09 -- changes to smag procedures
(  279) # 1 "./ldfsmag.h90"
(    1) !!DB 2008.04.11
(    2) !!Combine ldfdyn_smag and ldftra_smag into 1 routine as much of the code is identical
(    3) !!This forces 1 call to _smag_ to handle dyn & tra cases
(    4) 
(    5)    !!----------------------------------------------------------------------
(    6)    !!                        ***  ldf_smag.h90  ***
(    7)    !!----------------------------------------------------------------------
(    8) 
(    9)    !!----------------------------------------------------------------------
(   10)    !!  OPA 9.0 , BIO (2005) 
(   11)    !!----------------------------------------------------------------------
(   12) 
(   13)    !!----------------------------------------------------------------------
(   14)    !!   'key_dynldf_smag || key_traldf_smag'      2D lateral eddy viscosity coefficients
(   15)    !!----------------------------------------------------------------------
(   16) 
(   17)    SUBROUTINE ldf_smag( kstp )
(   18)       !!----------------------------------------------------------------------
(   19)       !!                  ***  ROUTINE ldf_dyn_smag  ***
(   20)       !!                   
(   21)       !! ** Purpose :   initializations of the horizontal ocean physics
(   22)       !!
(   23)       !! ** Method  :   3D eddy viscosity coef. ( longitude, latitude, depth )
(   24)       !!       laplacian operator   : ahm1, ahm2 defined at T- and F-points
(   25)       !!                              ahm3, ahm4 never used
(   26)       !!       bilaplacian operator : ahm1, ahm2 never used
(   27)       !!                           :  ahm3, ahm4 defined at U- and V-points
(   28)       !!       ??? explanation of the default is missing
(   29)       !!       History: 
(   30)       !!              1.0 ! 07-05 Zeliang Wang ! Original Codes 
(   31)       !!----------------------------------------------------------------------
(   32)       !! * Modules used
(   33) # 34
(   34)      USE ldftra_oce
(   35) # 36
(   36)       !! * Arguments
(   37) !DB: NB -- not used
(   38)      INTEGER, INTENT(in) :: kstp  
(   39) 
(   40)       !! * local variables
(   41)       INTEGER  ::   ji, jj, jk      ! dummy loop indices
(   42)       REAL(wp) ::   &
(   43)          alfa, bgvalue,dh_max,dh_min  ! SMAG coeff. 
(   44)       REAL(wp) ::   &
(   45)          dudx,dvdy,dvdx,dudy   ! SMAG coeff. 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 17

(   46)       REAL(wp) ::   &
(   47)          dvdx1,dvdx2,dudy1,dudy2  ! SMAG coeff. 
(   48)       REAL(wp), DIMENSION(jpk) ::   zcoef       ! temporary workspace
(   49)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: tmph  ! temporary workspace
(   50) # 51
(   51) !DB: Prandtl number
(   52)       REAL(wp) ::   Pr, &
(   53)          dt_max,dt_min! SMAG coeff. 
(   54) 
(   55) 
(   56) !DB: take Pr = diffusivity(tracer/momentum) which may be wrong 
(   57) !    Note that this value was suggested by ZW
(   58)       Pr = 1.0/10.0
(   59) 
(   60) # 61
(   61) !!NB: using same alfa for tracers as momentum
(   62)       alfa=0.1
(   63)       !!----------------------------------------------------------------------
(   64)       
(   65)       ! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operators
(   66)       ! =================                       whatever its orientation is)
(   67)       DO ji=2,jpi-1
(   68)          DO jj=2,jpj-1
(   69)             DO jk=1,jpk
(   70)                
(   71)                dudx = (ub(ji,jj,jk)-ub(ji-1,jj,jk)) / e1t(ji,jj)*tmask(ji,jj,jk)
(   72)                
(   73)                dvdy = (vb(ji,jj,jk)-vb(ji,jj-1,jk)) / e2t(ji,jj)*tmask(ji,jj,jk)
(   74)                
(   75)                dudy1 = 0.5*((ub(ji,jj+1,jk)-ub(ji,jj,jk)) /  e2f(ji,jj) *fmask(ji,jj,jk)& 
(   76)                     &    + (ub(ji,jj,jk)-ub(ji,jj-1,jk)) / e2f(ji,jj-1)*fmask(ji,jj-1,jk))
(   77)                
(   78)                dudy2 = 0.5*((ub(ji-1,jj+1,jk)-ub(ji-1,jj,jk)) / e2f(ji-1,jj)*fmask(ji-1,jj,jk) &
(   79)                     &    + (ub(ji-1,jj,jk)-ub(ji-1,jj-1,jk)) / e2f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   80)                
(   81)                dvdx1 = 0.5*((vb(ji+1,jj,jk)-vb(ji,jj,jk)) / e1f(ji,jj)*fmask(ji,jj,jk) &
(   82)                     &    + (vb(ji,jj,jk)-vb(ji-1,jj,jk)) / e1f(ji-1,jj)*fmask(ji-1,jj,jk))
(   83)                
(   84)                dvdx2 = 0.5*((vb(ji+1,jj-1,jk)-vb(ji,jj-1,jk)) / e1f(ji,jj-1)*fmask(ji,jj-1,jk) &
(   85)                     &    + (vb(ji,jj-1,jk)-vb(ji-1,jj-1,jk)) / e1f(ji-1,jj-1)*fmask(ji-1,jj-1,jk))
(   86)                
(   87)                dudy=0.5*(dudy1+dudy2)
(   88)                dvdx=0.5*(dvdx1+dvdx2)
(   89)                
(   90)                
(   91)                tmph(ji,jj,jk) = alfa*e1t(ji,jj)*e2t(ji,jj)*SQRT(dudx**2+dvdy**2+0.5*(dvdx+dudy)**2)  
(   92)                
(   93)             END DO
(   94)          END DO
(   95)       END DO
(   96)       
(   97)       CALL lbc_lnk(tmph,'T',1.)
(   98)       
(   99)       IF( ln_dynldf_lap ) THEN
(  100)          ! define ahm1 and ahm2 at the right grid point position
(  101)          ! (USER: modify ahm1 and ahm2 following your desiderata)
(  102)          
(  103)          ahm1=ahm0+tmph






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 18

(  104)          CALL lbc_lnk(ahm1,'T',1.)
(  105)          DO ji=1,jpi-1
(  106)             DO jj=1,jpj-1
(  107)                
(  108)                DO jk=1,jpk
(  109)                   ahm2(ji,jj,jk)=(ahm1(ji,jj,jk)*tmask(ji,jj,jk)+ahm1(ji+1,jj,jk)*tmask(ji+1,jj,jk)+      &
(  110)                        &             ahm1(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)+ahm1(ji,jj+1,jk)*tmask(ji,jj+1,jk))/   &
(  111)                        &             max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk)+tmask(ji+1,jj+1,jk)+tmask(ji,jj+1,jk),1.0)
(  112)                END DO
(  113)                
(  114)             END DO
(  115)          END DO
(  116)          CALL lbc_lnk(ahm2,'F',1.)
(  117)          
(  118) !!DB: 2009.10.02 -- move to below and get rid of minval as it is always zero
(  119) !!                  Also, report once-per-day         
(  120) !         dh_max=MAXVAL(ahm1(:,:,:))     
(  121) !         IF(lk_mpp) CALL mpp_max(dh_max)
(  122) !         dh_min=MINVAL(ahm1(:,:,:))
(  123) !         IF(lk_mpp) CALL mpp_min(dh_min)
(  124) !         IF(lwp) WRITE(numout,*) 'laplacian operator: ahm calculated by SMAG'
(  125) !         IF(lwp) WRITE(numout,*) ' MAX_DH',dh_max,'MIN_DH',dh_min
(  126)          
(  127)          
(  128)       ENDIF
(  129)       
(  130)       
(  131)       
(  132)       ! ahm3 and ahm4 at U- and V-points (used for bilaplacian operator
(  133)       ! ================================  whatever its orientation is)
(  134)       ! (USER: modify ahm3 and ahm4 following your desiderata)
(  135)       ! Here: ahm is proportional to the cube of the maximum of the gridspacing
(  136)       !       in the to horizontal direction
(  137)       !  Warning SHOUlD BE PUT AT U,V POINT!!!!!!!!
(  138)       IF( ln_dynldf_bilap ) THEN
(  139)          
(  140)          
(  141)          DO ji=1,jpi-1
(  142)             DO jj=1,jpj-1
(  143)                DO jk=1,jpk 
(  144)                   ahm3(ji,jj,jk) =ahm0-(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji+1,jj,jk)*tmask(ji+1,jj,jk))* &
(  145)                        &             min(e1u(ji,jj),e2u(ji,jj))**2/8. &
(  146)                        &             /max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk),1.0) 
(  147)                   ahm4(ji,jj,jk) =ahm0-(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji,jj+1,jk)*tmask(ji,jj+1,jk))* &
(  148)                        &             min(e1v(ji,jj),e2v(ji,jj))**2/8. &
(  149)                        &             /max(tmask(ji,jj,jk)+tmask(ji,jj+1,jk),1.0) 
(  150)                END DO
(  151)             END DO
(  152)          END DO
(  153)          
(  154)          
(  155)          CALL lbc_lnk(ahm3,'U',1.)
(  156)          CALL lbc_lnk(ahm4,'V',1.)
(  157) 
(  158)       END IF
(  159) !!END of ZW momentum SMAG
(  160) 
(  161) 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 19

(  162) # 163
(  163)       tmph(:,:,:) = Pr*tmph(:,:,:)       !DB
(  164)       ahtw(:,:,:) = tmph(:,:,:)+aht0
(  165)       CALL lbc_lnk(ahtw,'W',1.)
(  166)       
(  167)       IF( ln_traldf_lap ) THEN
(  168)          
(  169)          DO ji=1,jpi-1
(  170)             DO jj=1,jpj          
(  171)                DO jk=1,jpk
(  172)                   ahtu(ji,jj,jk)=(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji+1,jj,jk)*tmask(ji+1,jj,jk)+aht0)/  &
(  173)                        &         max(tmask(ji,jj,jk)+tmask(ji+1,jj,jk),1.0)
(  174)                END DO
(  175)             END DO
(  176)          END DO
(  177)          
(  178)          DO ji=1,jpi
(  179)             DO jj=1,jpj-1          
(  180)                DO jk=1,jpk
(  181)                   ahtv(ji,jj,jk)=(tmph(ji,jj,jk)*tmask(ji,jj,jk)+tmph(ji,jj+1,jk)*tmask(ji,jj+1,jk)+aht0)/  &
(  182)                        &         max(tmask(ji,jj,jk)+tmask(ji,jj+1,jk),1.0)
(  183)                END DO
(  184)             END DO
(  185)          END DO
(  186) 
(  187)          CALL lbc_lnk(ahtu,'U',1.)
(  188)          CALL lbc_lnk(ahtv,'V',1.)
(  189)          
(  190)          
(  191)       ENDIF
(  192)         
(  193)       
(  194)       IF( ln_traldf_bilap ) THEN
(  195)          
(  196)          DO jk=1,jpk
(  197)             ahtt(:,:,jk)=aht0-tmph(:,:,jk)*e1t(:,:)*e2t(:,:)/8.0
(  198)          END DO
(  199)          
(  200)          
(  201)       END IF
(  202) 
(  203) 
(  204) !!DB: 2009.10.02 -- report values 4x-per-day         
(  205) !! Just do tracer 
(  206)       if(mod(kstp-nit000,int(0.25*rday/rdt)) == 0) then
(  207)          call mpp_maxloc(ahtw,tmask,dt_max,ji,jj,jk)
(  208)          IF(lwp) then
(  209)             write(numout2,'(a30,i10,2x,f10.0,2x,3(i4,1x))') 'kt; SMAG max Ah(i,j,k): ',kstp, dt_max,ji,jj,jk
(  210)          endif
(  211)       endif
(  212) 
(  213) 
(  214) 
(  215) # 216
(  216)       
(  217)     END SUBROUTINE ldf_smag
(  218)     
(  219) # 281 "ldfdyn.F90"






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 20

(  281) 
(  282) # 283
(  283)    SUBROUTINE ldf_zpf_1d( ld_print, pdam, pwam, pbot, pdep, pah )
(  284)       !!----------------------------------------------------------------------
(  285)       !!                  ***  ROUTINE ldf_zpf  ***
(  286)       !!                   
(  287)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  288)       !!
(  289)       !! ** Method  :   1D eddy viscosity coefficients ( depth )
(  290)       !!
(  291)       !!----------------------------------------------------------------------
(  292)       !! * Arguments
(  293)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  294)       REAL(wp), INTENT (in   ) ::   &
(  295)           pdam,     &  ! depth of the inflection point
(  296)           pwam,     &  ! width of inflection
(  297)           pbot         ! battom value (0<pbot<= 1)
(  298)       REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
(  299)           pdep         ! depth of the gridpoint (T, U, V, F)
(  300)       REAL(wp), INTENT (inout), DIMENSION(jpk) ::   &
(  301)           pah          ! adimensional vertical profile
(  302) 
(  303)       !! * Local variables
(  304)       INTEGER  ::   jk           ! dummy loop indices
(  305)       REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
(  306)       !!----------------------------------------------------------------------
(  307) 
(  308)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
(  309)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  310)       zmhs = zm00 / zm01
(  311)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  312) 
(  313)       DO jk = 1, jpk
(  314)          pah(jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
(  315)       END DO
(  316) 
(  317)       ! Control print
(  318)       IF(lwp .AND. ld_print ) THEN
(  319)          WRITE(numout,*)
(  320)          WRITE(numout,*) '         ahm profile : '
(  321)          WRITE(numout,*)
(  322)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  323)          DO jk = 1, jpk
(  324)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(jk), pdep(jk)
(  325)          END DO
(  326)       ENDIF
(  327) 
(  328)    END SUBROUTINE ldf_zpf_1d
(  329) 
(  330) 
(  331)    SUBROUTINE ldf_zpf_1d_3d( ld_print, pdam, pwam, pbot, pdep, pah )
(  332)       !!----------------------------------------------------------------------
(  333)       !!                  ***  ROUTINE ldf_zpf  ***
(  334)       !!
(  335)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  336)       !!
(  337)       !! ** Method  :   1D eddy viscosity coefficients ( depth )
(  338)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 21

(  339)       !!----------------------------------------------------------------------
(  340)       !! * Arguments
(  341)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  342)       REAL(wp), INTENT (in   ) ::   &
(  343)           pdam,     &  ! depth of the inflection point
(  344)           pwam,     &  ! width of inflection
(  345)           pbot         ! battom value (0<pbot<= 1)
(  346)       REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
(  347)           pdep         ! depth of the gridpoint (T, U, V, F)
(  348)       REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
(  349)           pah          ! adimensional vertical profile
(  350) 
(  351)       !! * Local variables
(  352)       INTEGER  ::   jk           ! dummy loop indices
(  353)       REAL(wp) ::   zm00, zm01, zmhb, zmhs, zcf  ! temporary scalars
(  354)       !!----------------------------------------------------------------------
(  355) 
(  356)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
(  357)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  358)       zmhs = zm00 / zm01
(  359)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  360) 
(  361)       DO jk = 1, jpk
(  362)          zcf = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
(  363)          pah(:,:,jk) = zcf
(  364)       END DO
(  365) 
(  366)       ! Control print
(  367)       IF(lwp .AND. ld_print ) THEN
(  368)          WRITE(numout,*)
(  369)          WRITE(numout,*) '         ahm profile : '
(  370)          WRITE(numout,*)
(  371)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  372)          DO jk = 1, jpk
(  373)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(jk)
(  374)          END DO
(  375)       ENDIF
(  376) 
(  377)    END SUBROUTINE ldf_zpf_1d_3d
(  378) 
(  379) 
(  380)    SUBROUTINE ldf_zpf_3d( ld_print, pdam, pwam, pbot, pdep, pah )
(  381)       !!----------------------------------------------------------------------
(  382)       !!                  ***  ROUTINE ldf_zpf  ***
(  383)       !!
(  384)       !! ** Purpose :   vertical adimensional profile for eddy coefficient
(  385)       !!
(  386)       !! ** Method  :   3D for partial step or s-coordinate
(  387)       !!
(  388)       !!----------------------------------------------------------------------
(  389)       !! * Arguments
(  390)       LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
(  391)       REAL(wp), INTENT (in   ) ::   &
(  392)           pdam,     &  ! depth of the inflection point
(  393)           pwam,     &  ! width of inflection
(  394)           pbot         ! reduction factor (surface value / bottom value)
(  395)       REAL(wp), INTENT (in   ), DIMENSION(jpi,jpj,jpk) ::   &
(  396)           pdep         ! dep of the gridpoint (T, U, V, F)






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 22

(  397)       REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
(  398)           pah          ! adimensional vertical profile
(  399) 
(  400)       !! * Local variables
(  401)       INTEGER  ::   jk           ! dummy loop indices
(  402)       REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
(  403)       !!----------------------------------------------------------------------
(  404) 
(  405)       zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )   
(  406)       zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
(  407)       zmhs = zm00 / zm01
(  408)       zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01
(  409) 
(  410)       DO jk = 1, jpk
(  411)          pah(:,:,jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(:,:,jk) ) / pwam )  )
(  412)       END DO
(  413) 
(  414)       ! Control print
(  415)       IF(lwp .AND. ld_print ) THEN
(  416)          WRITE(numout,*)
(  417)          WRITE(numout,*) '         ahm profile : '
(  418)          WRITE(numout,*)
(  419)          WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
(  420)          DO jk = 1, jpk
(  421)             WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(1,1,jk)
(  422)          END DO
(  423)       ENDIF
(  424) 
(  425)    END SUBROUTINE ldf_zpf_3d
(  426)    !!======================================================================
(  427) END MODULE ldfdyn
