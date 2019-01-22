


PGF90 (Version     10.5)          09/23/2010  09:50:23      page 1

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

Filename: trctrp_ctl.F90

(    1) MODULE trctrp_ctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trctrp_ctl  ***
(    4)    !! Ocean passive tracers:  transport option control
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_trp_ctl  : control the different options of transport
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce_trc             ! ocean dynamics and active tracers variables
(   12)    USE trc                 ! ocean passive tracers variables
(   13)    USE trctrp_lec          ! passive tracers transport
(   14) 
(   15)    IMPLICIT NONE
(   16)    PRIVATE
(   17) 
(   18)    !! * Accessibility
(   19)    PUBLIC trc_trp_ctl   
(   20) 
(   21)    !! * Module variable
(   22) #if defined key_trcldf_eiv
(   23)       LOGICAL, PARAMETER ::   lk_trcldf_eiv   = .TRUE.   !: eddy induced velocity flag
(   24) #else   
(   25)       LOGICAL, PARAMETER ::   lk_trcldf_eiv   = .FALSE.  !: eddy induced velocity flag
(   26) #endif
(   27) 
(   28)    !!----------------------------------------------------------------------
(   29)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp_ctl.F90,v 1.10 2006/04/11 13:49:00 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE trc_trp_ctl
(   37)       !!---------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE trc_trp_ctl  ***
(   39)       !!                
(   40)       !! ** Purpose :   Control the consistency between cpp options for 
(   41)       !!      tracer transport
(   42)       !!
(   43)       !! History :
(   44)       !!   9.0  !  04-0.  (C. Ethe) 
(   45)       !!----------------------------------------------------------------------
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(   49)       !!----------------------------------------------------------------------
(   50) 
(   51)       !! Control of Advection scheme options
(   52)       CALL trc_adv_ctl
(   53) 
(   54)       !! Control of Lateral diffusion scheme options
(   55)       CALL trc_ldf_ctl
(   56) 
(   57)       !! Control of Vertival diffusion scheme options
(   58)       CALL trc_zdf_ctl
(   59) 
(   60)       !! Control of Newtonian damping  options
(   61)       IF(lwp) THEN
(   62)          WRITE(numout,*) ' *** Tracer damping option'
(   63)          WRITE(numout,*)
(   64)       ENDIF
(   65) 
(   66) #if defined key_trcdmp
(   67)       IF(lwp) THEN 
(   68)          WRITE(numout,*)' key_trcdmp is defined'
(   69)          WRITE(numout,*)' Check trcdmp ROUTINE '
(   70)          WRITE(numout,*)'  '
(   71)       ENDIF 
(   72)       CALL trc_dmp_ctl
(   73) #else
(   74)       IF (lwp) WRITE(numout,*) ' No tracer damping'
(   75) #endif
(   76) 
(   77) 
(   78)    END SUBROUTINE trc_trp_ctl
(   79) 
(   80)    SUBROUTINE trc_adv_ctl
(   81)       !!---------------------------------------------------------------------
(   82)       !!                  ***  ROUTINE trc_adv_ctl  ***
(   83)       !!                
(   84)       !! ** Purpose :   Control the consistency between cpp options for 
(   85)       !!      tracer advection schemes
(   86)       !!
(   87)       !! History :
(   88)       !!   8.5  !  02-11  (G. Madec)  Original code
(   89)       !!   9.0  !  04-0.  (C. Ethe)  adapted for passive tracers
(   90)       !!----------------------------------------------------------------------
(   91) 
(   92)       !! * Local declarations
(   93)       INTEGER ::   ioptio
(   94) 
(   95) 
(   96)       !!----------------------------------------------------------------------
(   97)       !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(   98)       !!----------------------------------------------------------------------
(   99) 
(  100)       ! Control of Advection scheme options
(  101)       ! -----------------------------------
(  102)       ioptio = 0
(  103)       IF( ln_trcadv_cen2   )   ioptio = ioptio + 1
(  104)       IF( ln_trcadv_tvd    )   ioptio = ioptio + 1






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 3

(  105)       IF( ln_trcadv_muscl  )   ioptio = ioptio + 1
(  106)       IF( ln_trcadv_muscl2 )   ioptio = ioptio + 1
(  107)       IF( ln_trcadv_smolar )   ioptio = ioptio + 1
(  108) 
(  109)       IF( lk_esopa ) THEN
(  110)          IF(lwp) WRITE(numout,*) ' esopa control : the use of all scheme is forced'
(  111)          ln_trcadv_cen2   = .TRUE.
(  112)          ln_trcadv_tvd    = .TRUE.
(  113)          ln_trcadv_muscl  = .TRUE.
(  114)          ln_trcadv_muscl2 = .TRUE.
(  115)          ln_trcadv_smolar = .TRUE.
(  116)       ELSEIF( ioptio > 1 .OR. ioptio == 0 ) THEN
(  117)          IF(lwp) WRITE(numout,cform_err)
(  118)          IF(lwp) WRITE(numout,*) ' Choose one advection scheme in namelist nam_trcadv'
(  119)          IF(lwp) WRITE(numout,*) '        ***                              ***********'
(  120)          nstop = nstop + 1
(  121)       ENDIF
(  122) 
(  123)       IF( n_cla == 1 .AND. .NOT. ln_trcadv_cen2 ) THEN
(  124)          IF(lwp) WRITE(numout,cform_err)
(  125)          IF(lwp) WRITE(numout,*) '     cross-land advection only with 2nd order advection scheme'
(  126)          nstop = nstop + 1
(  127)       ENDIF
(  128) 
(  129)       IF( lk_trccfg_1d ) THEN
(  130)          ln_trcadv_cen2   = .FALSE.    ;  ln_trcadv_tvd    = .FALSE. ; ln_trcadv_muscl  = .FALSE.
(  131)          ln_trcadv_muscl2 = .FALSE.    ;  ln_trcadv_smolar = .FALSE.
(  132)          IF(lwp) WRITE(numout,*) ' *******  1D configuration : No advection on passive tracers *******'
(  133)          IF(lwp) WRITE(numout,*) ' *******                                                     *******'
(  134)       ENDIF
(  135) 
(  136)    END SUBROUTINE trc_adv_ctl
(  137) 
(  138)    SUBROUTINE trc_ldf_ctl
(  139)       !!----------------------------------------------------------------------
(  140)       !!                  ***  ROUTINE trc_ldf_ctl  ***
(  141)       !! 
(  142)       !! ** Purpose :   Control the consistency between cpp options for 
(  143)       !!      tracer lateral diffusion 
(  144)       !!
(  145)       !! History :
(  146)       !!   9.0  !  03-04  (C. Ethe) 
(  147)       !!----------------------------------------------------------------------
(  148)       !! * Local declarations
(  149)       INTEGER ::   ioptio               ! ???
(  150)       LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout      
(  151) 
(  152)       !!----------------------------------------------------------------------
(  153)       !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(  154)       !!----------------------------------------------------------------------
(  155) 
(  156)       ! Parameter control
(  157) 
(  158)       ! control the input
(  159)       ioptio = 0
(  160)       IF( ln_trcldf_lap   )   ioptio = ioptio + 1
(  161)       IF( ln_trcldf_bilap )   ioptio = ioptio + 1
(  162)       IF( ioptio /= 1 )   THEN






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 4

(  163)           IF(lwp) WRITE(numout,cform_err)
(  164)           IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on tracer'
(  165)           nstop = nstop + 1
(  166)       ENDIF
(  167)       ioptio = 0
(  168)       IF( ln_trcldf_level )   ioptio = ioptio + 1
(  169)       IF( ln_trcldf_hor   )   ioptio = ioptio + 1
(  170)       IF( ln_trcldf_iso   )   ioptio = ioptio + 1
(  171)       IF( ioptio /= 1 ) THEN
(  172)          IF(lwp) WRITE(numout,cform_err)
(  173)          IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
(  174)          nstop = nstop + 1
(  175)       ENDIF
(  176) 
(  177)       ! ... Choice of the lateral scheme used
(  178)       IF( lk_trcldf_eiv ) THEN
(  179)          IF(lwp) WRITE(numout,*) '          eddy induced velocity on tracers'
(  180)             IF( .NOT.ln_trcldf_iso .OR. ln_trcldf_bilap ) THEN
(  181)             IF(lwp) WRITE(numout,cform_err)
(  182)             IF(lwp) WRITE(numout,*) ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion'
(  183)             nstop = nstop + 1
(  184)          ENDIF
(  185)       ENDIF
(  186) 
(  187)       IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal mixing
(  188)          IF( ( ln_trcldf_iso .OR. ln_trcldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
(  189)             IF(lwp) WRITE(numout,cform_err)
(  190)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  191)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  192)          ENDIF
(  193)       ELSE                       ! z-coordinates with/without partial step:
(  194)          ln_trcldf_level = ln_trcldf_level .OR. ln_trcldf_hor      ! level diffusion = horizontal diffusion
(  195)          ln_trcldf_hor   = .FALSE.
(  196)          IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_trcldf_level = T'
(  197)          IF(lwp) WRITE(numout,*) '                                                  and    force ln_trcldf_hor   = F'
(  198)          IF( ln_trcldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
(  199)             IF(lwp) WRITE(numout,cform_err)
(  200)             IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
(  201)             IF( .NOT.lk_esopa )   nstop = nstop + 1
(  202)          ENDIF
(  203)       ENDIF
(  204) 
(  205)       l_trcldf_lap     =  ln_trcldf_lap   .AND. ln_trcldf_level     ! iso-level   laplacian operator
(  206)       l_trcldf_bilap   =  ln_trcldf_bilap .AND. ln_trcldf_level     ! iso-level bilaplacian operator
(  207)       l_trcldf_bilapg  =  ln_trcldf_bilap .AND. ln_trcldf_hor       ! geopotential bilap. (s-coord)
(  208)       l_trcldf_iso     =  ln_trcldf_lap   .AND.                  &  ! laplacian operator
(  209)          &                   ( ln_trcldf_iso   .OR.  ln_trcldf_hor )  &  ! iso-neutral (z-coord) or horizontal (s-coord)
(  210)          &                                     .AND. .NOT.lk_zps
(  211)       l_trcldf_iso_zps =       ln_trcldf_lap   .AND.                  &  ! laplacian operator
(  212)          &                   ( ln_trcldf_iso   .OR.  ln_trcldf_hor )  &  ! iso-neutral (partial steps)
(  213)          &                                     .AND. lk_zps              ! or geopotential in mixed partial steps/s-coord
(  214)       l_trczdf_iso    = .FALSE.
(  215)       l_trczdf_iso_vo = .FALSE.
(  216)       IF( l_trcldf_iso     )   l_trczdf_iso = .TRUE.
(  217)       IF( l_trcldf_iso_zps )   l_trczdf_iso = .TRUE.
(  218) #if defined key_vectopt_memory
(  219)       IF( l_trczdf_iso ) THEN
(  220)          l_trczdf_iso    = .FALSE.






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 5

(  221)          l_trczdf_iso_vo = .TRUE.
(  222)       ENDIF
(  223) #endif
(  224) 
(  225)  
(  226)       ioptio = 0
(  227)       IF( l_trcldf_lap     )   ioptio = ioptio + 1
(  228)       IF( l_trcldf_bilap   )   ioptio = ioptio + 1
(  229)       IF( l_trcldf_bilapg  )   ioptio = ioptio + 1
(  230)       IF( l_trcldf_iso     )   ioptio = ioptio + 1
(  231)       IF( l_trcldf_iso_zps )   ioptio = ioptio + 1
(  232)       IF( ioptio /= 1 ) THEN
(  233)          IF(lwp) WRITE(numout,cform_err)
(  234)          IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
(  235)          nstop = nstop + 1
(  236)       ENDIF
(  237) 
(  238)       IF( lk_esopa ) THEN
(  239)          l_trcldf_lap = .TRUE.   ;   l_trcldf_bilap   = .TRUE.   ;   l_trcldf_bilapg  = .TRUE.
(  240)          l_trcldf_iso = .TRUE.   ;   l_trcldf_iso_zps = .TRUE.
(  241)          l_trczdf_iso = .TRUE.   ;   l_trczdf_iso_vo  = .TRUE.
(  242)          IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
(  243)       ENDIF
(  244) 
(  245)       IF( .NOT. ln_trcldf_diff .OR. lk_trccfg_1d ) THEN
(  246)          l_trcldf_lap = .FALSE.   ;   l_trcldf_bilap   = .FALSE.   ;   l_trcldf_bilapg  = .FALSE.
(  247)          l_trcldf_iso = .FALSE.   ;   l_trcldf_iso_zps = .FALSE.
(  248)          l_trczdf_iso = .FALSE.   ;   l_trczdf_iso_vo  = .FALSE.
(  249)          IF(lwp ) WRITE(numout,*) '************* No lateral physics on passive tracers *****************'
(  250)          IF(lwp ) WRITE(numout,*) '*************                                       *****************'
(  251)       ELSE
(  252)          ! ... Space variation of eddy coefficients
(  253)          ioptio = 0
(  254) #if defined key_traldf_c3d
(  255)          IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( latitude, longitude, depth)'
(  256)          ioptio = ioptio + 1
(  257) #endif
(  258) #if defined key_traldf_c2d
(  259)          IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( latitude, longitude)'
(  260)          ioptio = ioptio + 1
(  261) #endif
(  262) #if defined key_traldf_c1d
(  263)          IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( depth )'
(  264)          ioptio = ioptio + 1
(  265)          IF( lk_sco ) THEN
(  266)             IF(lwp) WRITE(numout,cform_err)
(  267)             IF(lwp) WRITE(numout,*) '          key_traldf_c1d cannot be used in s-coordinate (key_s_coord)'
(  268)             nstop = nstop + 1
(  269)          ENDIF
(  270) #endif
(  271)          IF( ioptio == 0 ) THEN
(  272)             IF(lwp) WRITE(numout,*) ' tracer mixing coef. = constant (default option)'
(  273)          ELSEIF( ioptio > 1 ) THEN
(  274)             IF(lwp) WRITE(numout,cform_err)
(  275)             IF(lwp) WRITE(numout,*) ' use only one of the following keys:',   &
(  276)                &                    ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d'
(  277)             nstop = nstop + 1
(  278)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 6

(  279)          
(  280)          IF( l_trcldf_bilap .OR. l_trcldf_bilapg ) THEN
(  281)             IF(lwp) WRITE(numout,*) '  biharmonic tracer diffusion'
(  282)             IF( ahtrc0 > 0 .AND. .NOT. lk_esopa ) THEN
(  283)                IF(lwp) WRITE(numout,cform_err)
(  284)                IF(lwp) WRITE(numout,*) ' The horizontal diffusivity coef. aht0 must be negative'
(  285)                nstop = nstop + 1
(  286)             ENDIF
(  287)          ELSE
(  288)             IF(lwp) WRITE(numout,*) ' harmonic tracer diffusion (default)'
(  289)             IF( ahtrc0 < 0 .AND. .NOT. lk_esopa ) THEN
(  290)                IF(lwp) WRITE(numout,cform_err)
(  291)                IF(lwp) WRITE(numout,*) 'The horizontal diffusivity coef. aht0 must be positive'
(  292)                nstop = nstop + 1
(  293)             ENDIF
(  294)          ENDIF
(  295)       ENDIF
(  296) 
(  297)    END SUBROUTINE trc_ldf_ctl
(  298) 
(  299)    SUBROUTINE trc_zdf_ctl
(  300)       !!----------------------------------------------------------------------
(  301)       !!                  ***  ROUTINE trc_zdf_ctl  ***
(  302)       !! 
(  303)       !! ** Purpose :     Control the consistency between cpp options for 
(  304)       !!      tracer vertical diffusion
(  305)       !!
(  306)       !!   9.0  !  04-03  (C. Ethe)  
(  307)       !!----------------------------------------------------------------------
(  308)       !! * Local declarations
(  309) 
(  310)       !!----------------------------------------------------------------------
(  311)       !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(  312)       !!----------------------------------------------------------------------
(  313) 
(  314)       ! Parameter & key controls
(  315)       ! ------------------------
(  316)       ! ... vertical mixing
(  317)       ! time stepping scheme (N.B. TKE scheme => force the use of implicit scheme)
(  318) #if defined key_zdftke
(  319)       l_trczdf_exp = .FALSE.          ! use implicit scheme
(  320)       l_trczdf_imp = .TRUE. 
(  321) #else
(  322) !!DB: I think ln_zdfexp is incorrect, should be ln_trczdf_exp as in passivetrc.namelist file
(  323)       IF( ln_zdfexp  ) THEN  
(  324)          l_trczdf_exp = .TRUE.           ! use explicit scheme
(  325)          l_trczdf_imp = .FALSE.
(  326)       ELSE
(  327)          l_trczdf_exp = .FALSE.          ! use implicit scheme
(  328)          l_trczdf_imp = .TRUE. 
(  329)       ENDIF
(  330) #endif
(  331) 
(  332)       IF( l_trczdf_iso .OR. l_trczdf_iso_vo ) THEN  
(  333)          l_trczdf_exp = .FALSE.          ! iso-neutral diffusion : 
(  334)          l_trczdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
(  335)       ENDIF
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 7

(  337) #if defined key_esopa
(  338)       l_trczdf_exp = .TRUE.           ! esopa: use all options
(  339)       l_trczdf_imp = .TRUE.
(  340) #endif
(  341) 
(  342) 
(  343)    END SUBROUTINE trc_zdf_ctl
(  344) 
(  345)    SUBROUTINE trc_dmp_ctl
(  346)       !!----------------------------------------------------------------------
(  347)       !!                  ***  ROUTINE trc_dmp_ctl  ***
(  348)       !! 
(  349)       !! ** Purpose :    Control the consistency between cpp options for 
(  350)       !!      tracer newtonian damping 
(  351)       !!
(  352)       !!
(  353)       !! History :
(  354)       !!   9.0  !  04-03  (C. Ethe) 
(  355)       !!----------------------------------------------------------------------
(  356) #if defined key_trcdmp
(  357) 
(  358)       SELECT CASE ( ndmptr )
(  359) 
(  360)       CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
(  361)          IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'
(  362) 
(  363)       CASE ( 1:90 )             ! Damping poleward of 'ndmptr' degrees
(  364)          IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmptr, ' degrees'
(  365) 
(  366)       CASE DEFAULT
(  367)          IF(lwp) WRITE(numout,cform_err)
(  368)          IF(lwp) WRITE(numout,*) '          bad flag value for ndmptr = ', ndmptr
(  369)          nstop = nstop + 1
(  370) 
(  371)       END SELECT
(  372) 
(  373) 
(  374)       SELECT CASE ( nmldmptr )
(  375) 
(  376)       CASE ( 0 )                ! newtonian damping throughout the water column
(  377)          IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'
(  378) 
(  379)       CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
(  380)          IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'
(  381) 
(  382)       CASE ( 2 )                ! no damping in the mixed layer 
(  383)          IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'
(  384) 
(  385)       CASE DEFAULT
(  386)          IF(lwp) WRITE(numout,cform_err)
(  387)          IF(lwp) WRITE(numout,*) '          bad flag value for nmldmptr = ', nmldmptr
(  388)          nstop = nstop + 1
(  389) 
(  390)       END SELECT
(  391) #endif
(  392)  
(  393)    END SUBROUTINE trc_dmp_ctl
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 8

(  395) #else
(  396)    !!----------------------------------------------------------------------
(  397)    !!   Dummy module :                      NO passive tracer
(  398)    !!----------------------------------------------------------------------
(  399) CONTAINS
(  400)    SUBROUTINE trc_trp_ctl             ! Empty routine
(  401)    END SUBROUTINE trc_trp_ctl
(  402) #endif
(  403)    
(  404)   !!======================================================================
(  405) END MODULE trctrp_ctl





















































PGF90 (Version     10.5)          09/23/2010  09:50:23      page 9

(    1) # 1 "trctrp_ctl.F90"
(    1) MODULE trctrp_ctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trctrp_ctl  ***
(    4)    !! Ocean passive tracers:  transport option control
(    5)    !!==============================================================================
(    6) # 396
(  396)    !!----------------------------------------------------------------------
(  397)    !!   Dummy module :                      NO passive tracer
(  398)    !!----------------------------------------------------------------------
(  399) CONTAINS
(  400)    SUBROUTINE trc_trp_ctl             ! Empty routine
(  401)    END SUBROUTINE trc_trp_ctl
(  402) # 403
(  403)    
(  404)   !!======================================================================
(  405) END MODULE trctrp_ctl
