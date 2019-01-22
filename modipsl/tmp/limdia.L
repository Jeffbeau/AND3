


PGF90 (Version     10.5)          09/23/2010  09:49:57      page 1

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

Filename: limdia.F90

(    1) MODULE limdia
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limdia   ***
(    4)    !!                      diagnostics of ice model 
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_dia      : computation of the time evolution of keys var.
(   11)    !!   lim_dia_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst          ! 
(   15)    USE par_ice         ! ice parameters
(   16)    USE ice_oce         ! ice variables
(   17)    USE daymod          !
(   18)    USE dom_ice         !
(   19)    USE ice             !
(   20)    USE iceini          !
(   21)    USE limistate       !
(   22)    USE in_out_manager  ! I/O manager
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC lim_dia       ! called by ice_step
(   29) 
(   30)    !! * Shared module variables
(   31)    INTEGER, PUBLIC  ::  &  !:
(   32)       ntmoy   = 1 ,     &  !: instantaneous values of ice evolution or averaging ntmoy
(   33)       ninfo   = 1          !: frequency of ouputs on file ice_evolu in case of averaging
(   34) 
(   35)    !! * Module variables
(   36)    INTEGER, PARAMETER ::   &  ! Parameters for outputs to files "evolu"
(   37)       jpinfmx = 100         ,    &  ! maximum number of key variables
(   38)       jpchinf = 5           ,    &  ! ???
(   39)       jpchsep = jpchinf + 2         ! ???
(   40) 
(   41)    INTEGER ::   &
(   42)       nfrinf  = 4 ,     &  ! number of variables written in one line 
(   43)       nferme ,          &  ! last time step at which the var. are written on file
(   44)       nvinfo ,          &  ! number of total variables 
(   45)       nbvt   ,          &  ! number of time variables
(   46)       naveg                ! number of step for accumulation before averaging






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 2

(   47) 
(   48)    CHARACTER(len=8) ::   &
(   49)       fmtinf  = '1PE13.5 ' ! format of the output values  
(   50)    CHARACTER(len=30) ::   &
(   51)       fmtw  ,           &  ! formats
(   52)       fmtr  ,           &  ! ???
(   53)       fmtitr               ! ???
(   54)    CHARACTER(len=jpchsep), DIMENSION(jpinfmx) ::   &
(   55)       titvar               ! title of key variables
(   56)  
(   57)    REAL(wp) ::   &
(   58)       epsi06 = 1.e-06      ! ???
(   59)    REAL(wp), DIMENSION(jpinfmx) ::  &
(   60)       vinfom               ! temporary working space
(   61)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   62)       aire                 ! masked grid cell area
(   63) 
(   64)    !! * Substitutions
(   65) #  include "vectopt_loop_substitute.h90"
(   66)    !!----------------------------------------------------------------------
(   67)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   68)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdia.F90,v 1.5 2005/03/27 18:34:41 opalod Exp $ 
(   69)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   70)    !!----------------------------------------------------------------------
(   71) 
(   72) CONTAINS
(   73) 
(   74)    SUBROUTINE lim_dia
(   75)       !!--------------------------------------------------------------------
(   76)       !!                  ***  ROUTINE lim_dia  ***
(   77)       !!   
(   78)       !! ** Purpose : Computation and outputs on file ice.evolu 
(   79)       !!      the temporal evolution of some key variables
(   80)       !!
(   81)       !! History :
(   82)       !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
(   83)       !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
(   84)       !!-------------------------------------------------------------------
(   85)       !! * Local variables
(   86)        INTEGER  ::   jv,ji, jj   ! dummy loop indices
(   87)        INTEGER  ::   nv          ! indice of variable 
(   88)        REAL(wp), DIMENSION(jpinfmx) ::  & 
(   89)           vinfor           ! temporary working space 
(   90)        REAL(wp) ::    &
(   91)           zarea    ,    &  ! sea ice area
(   92)           zldarea  ,    &  ! leads area
(   93)           zextent15,    &  ! sea ice extent (15%)
(   94)           zextent85,    &  ! sea ice extent (85%)
(   95)           zicevol  ,    &  ! sea ice volume
(   96)           zsnwvol  ,    &  ! snow volume over sea ice
(   97)           zicespd          ! sea ice velocity
(   98)        !!-------------------------------------------------------------------
(   99) 
(  100)        IF( numit == nstart )   CALL lim_dia_init   ! initialisation of ice_evolu file      
(  101) 
(  102)        ! computation of key variables at each time step   
(  103) 
(  104)        nv = 1 






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 3

(  105)        vinfor(nv) = REAL( numit )
(  106)        nv = nv + 1
(  107)        vinfor(nv) = nyear
(  108)  
(  109)        DO jv = nbvt + 1, nvinfo
(  110)           vinfor(jv) = 0.e0
(  111)        END DO
(  112) 
(  113)        zextent15 = 0.e0
(  114)        zextent85 = 0.e0
(  115)        ! variables in northern Hemis
(  116)        DO jj = njeq, jpjm1
(  117)           DO ji = fs_2, fs_jpim1   ! vector opt.
(  118)              IF( tms(ji,jj) == 1 ) THEN
(  119)                 zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
(  120)                 IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
(  121)                 IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
(  122)                 zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
(  123)                 zicevol = zarea   * hicif(ji,jj)
(  124)                 zsnwvol = zarea   * hsnif(ji,jj)
(  125)                 zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
(  126)                    &                + v_ice(ji,jj) * v_ice(ji,jj) )
(  127)                 vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
(  128)                 vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
(  129)                 vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
(  130)                 vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
(  131)                 vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
(  132)                 vinfor(nv+11) = vinfor(nv+11) + zsnwvol
(  133)                 vinfor(nv+13) = vinfor(nv+13) + zicespd
(  134)              ENDIF
(  135)           END DO
(  136)        END DO
(  137)        vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )
(  138) 
(  139) 
(  140)       ! variables in southern Hemis
(  141)        nv = nv + 1
(  142)        DO jj = 2, njeqm1
(  143)           DO ji = fs_2, fs_jpim1   ! vector opt.
(  144)              IF( tms(ji,jj) == 1 ) THEN
(  145)                 zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
(  146)                 IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
(  147)                 IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
(  148)                 zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
(  149)                 zicevol = zarea   * hicif(ji,jj)
(  150)                 zsnwvol = zarea   * hsnif(ji,jj)
(  151)                 zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
(  152)                    &                + v_ice(ji,jj) * v_ice(ji,jj) )
(  153)                 vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
(  154)                 vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
(  155)                 vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
(  156)                 vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
(  157)                 vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
(  158)                 vinfor(nv+11) = vinfor(nv+11) + zsnwvol
(  159)                 vinfor(nv+13) = vinfor(nv+13) + zicespd
(  160)              ENDIF
(  161)           END DO
(  162)        END DO






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 4

(  163)        vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )    
(  164) 
(  165)        !  Accumulation before averaging 
(  166)        DO jv = 1, nvinfo
(  167)           vinfom(jv) = vinfom(jv) + vinfor(jv)
(  168)        END DO
(  169)        naveg = naveg + 1  
(  170)     
(  171)        ! oututs on file ice_evolu    
(  172)        IF( MOD( numit , ninfo ) == 0 ) THEN
(  173)           WRITE(numevo_ice,fmtw) ( titvar(jv), vinfom(jv)/naveg, jv = 1, nvinfo )
(  174)           naveg = 0
(  175)           DO jv = 1, nvinfo
(  176)              vinfom(jv) = 0.e0
(  177)           END DO
(  178)        ENDIF
(  179)   
(  180)     END SUBROUTINE lim_dia
(  181)  
(  182) 
(  183)     SUBROUTINE lim_dia_init
(  184)        !!-------------------------------------------------------------------
(  185)        !!                  ***  ROUTINE lim_dia_init  ***
(  186)        !!             
(  187)        !! ** Purpose : Preparation of the file ice_evolu for the output of
(  188)        !!      the temporal evolution of key variables
(  189)        !!
(  190)        !! ** input   : Namelist namicedia
(  191)        !!
(  192)        !! history :
(  193)        !!  8.5  ! 03-08 (C. Ethe) original code
(  194)        !!-------------------------------------------------------------------
(  195)        NAMELIST/namicedia/fmtinf, nfrinf, ninfo, ntmoy
(  196) 
(  197)        INTEGER  ::   jv   ,     &  ! dummy loop indice
(  198)           &          ntot ,     &
(  199)           &          ndeb ,     &
(  200)           &          irecl
(  201) 
(  202)        INTEGER  ::   nv            ! indice of variable 
(  203) 
(  204)        REAL(wp) ::   zxx0, zxx1    ! temporary scalars
(  205) 
(  206)        CHARACTER(len=jpchinf) ::   titinf
(  207)        !!-------------------------------------------------------------------
(  208) 
(  209)        ! Read Namelist namicedia
(  210)        REWIND ( numnam_ice )
(  211)        READ   ( numnam_ice  , namicedia )
(  212)        IF(lwp) THEN
(  213)           WRITE(numout,*)
(  214)           WRITE(numout,*) 'lim_dia_init : ice parameters for ice diagnostics '
(  215)           WRITE(numout,*) '~~~~~~~~~~~~'
(  216)           WRITE(numout,*) '   format of the output values                                 fmtinf = ', fmtinf
(  217)           WRITE(numout,*) '   number of variables written in one line                     nfrinf = ', nfrinf 
(  218)           WRITE(numout,*) '   Instantaneous values of ice evolution or averaging          ntmoy  = ', ntmoy
(  219)           WRITE(numout,*) '   frequency of ouputs on file ice_evolu in case of averaging  ninfo  = ', ninfo
(  220)        ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 5

(  221) 
(  222)        ! masked grid cell area
(  223)        aire(:,:) = area(:,:) * tms(:,:)
(  224) 
(  225)        ! Titles of ice key variables :
(  226)        nv = 1
(  227)        titvar(nv) = 'NoIt'  ! iteration number
(  228)        nv = nv + 1
(  229)        titvar(nv) = 'T yr'  ! time step in years
(  230)        nv = nv + 1
(  231) 
(  232)        nbvt = nv - 1
(  233) 
(  234)        titvar(nv) = 'AEFN' ! sea ice area in the northern Hemisp.(10^12 km2)
(  235)        nv = nv + 1
(  236)        titvar(nv) = 'AEFS' ! sea ice area in the southern Hemisp.(10^12 km2)
(  237)        nv = nv + 1
(  238)        titvar(nv) = 'A15N'  ! sea ice extent (15%) in the northern Hemisp.(10^12 km2)
(  239)        nv = nv + 1
(  240)        titvar(nv) = 'A15S'  ! sea ice extent (15%) in the southern Hemisp.(10^12 km2)
(  241)        nv = nv + 1
(  242)        titvar(nv) = 'A85N'  ! sea ice extent (85%) in the northern Hemisp.(10^12 km2)
(  243)        nv = nv + 1
(  244)        titvar(nv) = 'A85S'  ! sea ice extent (85%) in the southern Hemisp.(10^12 km2)
(  245)        nv = nv + 1
(  246)        titvar(nv) = 'ALEN'  ! leads area in the northern Hemisp.(10^12 km2)
(  247)        nv = nv + 1
(  248)        titvar(nv) = 'ALES'  ! leads area in the southern Hemisp.(10^12 km2)
(  249)        nv = nv + 1
(  250)        titvar(nv) = 'VOLN'  ! sea ice volume in the northern Hemisp.(10^3 km3)
(  251)        nv = nv + 1
(  252)        titvar(nv) = 'VOLS'  ! sea ice volume in the southern Hemisp.(10^3 km3)
(  253)        nv = nv + 1
(  254)        titvar(nv) = 'VONN'  ! snow volume over sea ice in the northern Hemisp.(10^3 km3)
(  255)        nv = nv + 1
(  256)        titvar(nv) = 'VONS'  ! snow volume over sea ice in the southern Hemisp.(10^3 km3)
(  257)        nv = nv + 1
(  258)        titvar(nv) = 'ECGN'  ! mean sea ice velocity in the northern Hemisp.(m/s)
(  259)        nv = nv + 1
(  260)        titvar(nv) = 'ECGS'  ! mean sea ice velocity in the southern Hemisp.(m/s)
(  261) 
(  262)        nvinfo = nv
(  263) 
(  264)        ! Definition et Ecriture de l'entete : nombre d'enregistrements 
(  265)        ndeb   = ( nstart - 1 ) / ninfo
(  266)        IF( nstart == 1 ) ndeb = -1
(  267) 
(  268)        nferme = ( nstart - 1 + nitrun) / ninfo
(  269)        ntot   = nferme - ndeb
(  270)        ndeb   = ninfo * ( 1 + ndeb )
(  271)        nferme = ninfo * nferme
(  272) 
(  273)        ! definition of formats 
(  274)        WRITE( fmtw  , '(A,I3,A2,I1,A)' )  '(', nfrinf, '(A', jpchsep, ','//fmtinf//'))'
(  275)        WRITE( fmtr  , '(A,I3,A,I1,A)'  )  '(', nfrinf, '(', jpchsep, 'X,'//fmtinf//'))'
(  276)        WRITE( fmtitr, '(A,I3,A,I1,A)'  )  '(', nvinfo, 'A', jpchinf, ')'
(  277) 
(  278)        ! opening  "ice_evolu" file






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 6

(  279)        irecl = ( jpchinf + 1 ) * nvinfo 
(  280)        OPEN( numevo_ice, file='ice.evolu', status='unknown', RECL = irecl)
(  281)        OPEN( numevo_ice, file='ice.evolu', status='unknown')
(  282) 
(  283)        !- ecriture de 2 lignes d''entete :
(  284)        WRITE(numevo_ice,1000) fmtr, fmtw, fmtitr, nvinfo, ntot, 0, nfrinf
(  285)        zxx0 = 0.001 * REAL( ninfo )
(  286)        zxx1 = 0.001 * REAL( ndeb  )
(  287)        WRITE(numevo_ice,1111) REAL(jpchinf), 0., zxx1, zxx0, 0., 0., 0
(  288) 
(  289)        !- ecriture de 2 lignes de titre :
(  290)        WRITE(numevo_ice,'(A,I8,A,I8,A,I5)')                                      &
(  291)           'Evolution chronologique - Experience '//cexper   &
(  292)           //'   de', ndeb, ' a', nferme, ' pas', ninfo
(  293)        WRITE(numevo_ice,fmtitr) ( titvar(jv), jv = 1, nvinfo )
(  294) 
(  295) 
(  296)        !--preparation de "titvar" pour l''ecriture parmi les valeurs numeriques :
(  297)        DO  jv = 2 , nvinfo
(  298)           titinf     = titvar(jv)(:jpchinf)
(  299)           titvar(jv) = '  '//titinf
(  300)        END DO
(  301) 
(  302)        !--Initialisation of the arrays for the accumulation
(  303)        DO  jv = 1, nvinfo
(  304)           vinfom(jv) = 0.
(  305)        END DO
(  306)        naveg = 0
(  307) 
(  308) 1000   FORMAT( 3(A20),4(1x,I6) )
(  309) 1111   FORMAT( 3(F7.1,1X,F7.3,1X),I3,A )  
(  310) 
(  311)     END SUBROUTINE lim_dia_init
(  312) 
(  313) #else
(  314)    !!----------------------------------------------------------------------
(  315)    !!   Default option :                               NO LIM sea-ice model
(  316)    !!----------------------------------------------------------------------
(  317) CONTAINS
(  318)    SUBROUTINE lim_dia         ! Empty routine
(  319)    END SUBROUTINE lim_dia
(  320) #endif
(  321) 
(  322)    !!======================================================================
(  323) END MODULE limdia



















PGF90 (Version     10.5)          09/23/2010  09:49:57      page 7

(    1) # 1 "limdia.F90"
(    1) MODULE limdia
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limdia   ***
(    4)    !!                      diagnostics of ice model 
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_dia      : computation of the time evolution of keys var.
(   11)    !!   lim_dia_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst          ! 
(   15)    USE par_ice         ! ice parameters
(   16)    USE ice_oce         ! ice variables
(   17)    USE daymod          !
(   18)    USE dom_ice         !
(   19)    USE ice             !
(   20)    USE iceini          !
(   21)    USE limistate       !
(   22)    USE in_out_manager  ! I/O manager
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC lim_dia       ! called by ice_step
(   29) 
(   30)    !! * Shared module variables
(   31)    INTEGER, PUBLIC  ::  &  !:
(   32)       ntmoy   = 1 ,     &  !: instantaneous values of ice evolution or averaging ntmoy
(   33)       ninfo   = 1          !: frequency of ouputs on file ice_evolu in case of averaging
(   34) 
(   35)    !! * Module variables
(   36)    INTEGER, PARAMETER ::   &  ! Parameters for outputs to files "evolu"
(   37)       jpinfmx = 100         ,    &  ! maximum number of key variables
(   38)       jpchinf = 5           ,    &  ! ???
(   39)       jpchsep = jpchinf + 2         ! ???
(   40) 
(   41)    INTEGER ::   &
(   42)       nfrinf  = 4 ,     &  ! number of variables written in one line 
(   43)       nferme ,          &  ! last time step at which the var. are written on file
(   44)       nvinfo ,          &  ! number of total variables 
(   45)       nbvt   ,          &  ! number of time variables
(   46)       naveg                ! number of step for accumulation before averaging
(   47) 
(   48)    CHARACTER(len=8) ::   &
(   49)       fmtinf  = '1PE13.5 ' ! format of the output values  
(   50)    CHARACTER(len=30) ::   &
(   51)       fmtw  ,           &  ! formats
(   52)       fmtr  ,           &  ! ???
(   53)       fmtitr               ! ???
(   54)    CHARACTER(len=jpchsep), DIMENSION(jpinfmx) ::   &
(   55)       titvar               ! title of key variables
(   56)  
(   57)    REAL(wp) ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 8

(   58)       epsi06 = 1.e-06      ! ???
(   59)    REAL(wp), DIMENSION(jpinfmx) ::  &
(   60)       vinfom               ! temporary working space
(   61)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   62)       aire                 ! masked grid cell area
(   63) 
(   64)    !! * Substitutions
(   65) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 66 "limdia.F90"
(   66) # 66
(   66)    !!----------------------------------------------------------------------
(   67)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   68)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdia.F90,v 1.5 2005/03/27 18:34:41 opalod Exp $ 
(   69)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   70)    !!----------------------------------------------------------------------
(   71) 
(   72) CONTAINS
(   73) 
(   74)    SUBROUTINE lim_dia
(   75)       !!--------------------------------------------------------------------
(   76)       !!                  ***  ROUTINE lim_dia  ***
(   77)       !!   
(   78)       !! ** Purpose : Computation and outputs on file ice.evolu 
(   79)       !!      the temporal evolution of some key variables
(   80)       !!
(   81)       !! History :
(   82)       !!   8.0  !  97-06  (Louvain-La-Neuve)  Original code
(   83)       !!   8.5  !  02-09  (C. Ethe , G. Madec )  F90: Free form and module
(   84)       !!-------------------------------------------------------------------
(   85)       !! * Local variables
(   86)        INTEGER  ::   jv,ji, jj   ! dummy loop indices
(   87)        INTEGER  ::   nv          ! indice of variable 
(   88)        REAL(wp), DIMENSION(jpinfmx) ::  & 
(   89)           vinfor           ! temporary working space 
(   90)        REAL(wp) ::    &
(   91)           zarea    ,    &  ! sea ice area
(   92)           zldarea  ,    &  ! leads area
(   93)           zextent15,    &  ! sea ice extent (15%)
(   94)           zextent85,    &  ! sea ice extent (85%)
(   95)           zicevol  ,    &  ! sea ice volume
(   96)           zsnwvol  ,    &  ! snow volume over sea ice
(   97)           zicespd          ! sea ice velocity
(   98)        !!-------------------------------------------------------------------
(   99) 
(  100)        IF( numit == nstart )   CALL lim_dia_init   ! initialisation of ice_evolu file      
(  101) 
(  102)        ! computation of key variables at each time step   






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 9

(  103) 
(  104)        nv = 1 
(  105)        vinfor(nv) = REAL( numit )
(  106)        nv = nv + 1
(  107)        vinfor(nv) = nyear
(  108)  
(  109)        DO jv = nbvt + 1, nvinfo
(  110)           vinfor(jv) = 0.e0
(  111)        END DO
(  112) 
(  113)        zextent15 = 0.e0
(  114)        zextent85 = 0.e0
(  115)        ! variables in northern Hemis
(  116)        DO jj = njeq, jpjm1
(  117)           DO ji = 2, jpim1   ! vector opt.
(  118)              IF( tms(ji,jj) == 1 ) THEN
(  119)                 zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
(  120)                 IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
(  121)                 IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
(  122)                 zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
(  123)                 zicevol = zarea   * hicif(ji,jj)
(  124)                 zsnwvol = zarea   * hsnif(ji,jj)
(  125)                 zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
(  126)                    &                + v_ice(ji,jj) * v_ice(ji,jj) )
(  127)                 vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
(  128)                 vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
(  129)                 vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
(  130)                 vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
(  131)                 vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
(  132)                 vinfor(nv+11) = vinfor(nv+11) + zsnwvol
(  133)                 vinfor(nv+13) = vinfor(nv+13) + zicespd
(  134)              ENDIF
(  135)           END DO
(  136)        END DO
(  137)        vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )
(  138) 
(  139) 
(  140)       ! variables in southern Hemis
(  141)        nv = nv + 1
(  142)        DO jj = 2, njeqm1
(  143)           DO ji = 2, jpim1   ! vector opt.
(  144)              IF( tms(ji,jj) == 1 ) THEN
(  145)                 zarea = ( 1.0 - frld(ji,jj) ) * aire(ji,jj)
(  146)                 IF (frld(ji,jj) <= 0.15 ) zextent15 = aire(ji,jj)    
(  147)                 IF (frld(ji,jj) <= 0.85 ) zextent85 = aire(ji,jj)   
(  148)                 zldarea = zarea   / MAX( ( 1 - frld(ji,jj) ) , epsi06 )
(  149)                 zicevol = zarea   * hicif(ji,jj)
(  150)                 zsnwvol = zarea   * hsnif(ji,jj)
(  151)                 zicespd = zicevol * ( u_ice(ji,jj) * u_ice(ji,jj)   &
(  152)                    &                + v_ice(ji,jj) * v_ice(ji,jj) )
(  153)                 vinfor(nv+ 1) = vinfor(nv+ 1) + zarea
(  154)                 vinfor(nv+ 3) = vinfor(nv+ 3) + zextent15
(  155)                 vinfor(nv+ 5) = vinfor(nv+ 5) + zextent85
(  156)                 vinfor(nv+ 7) = vinfor(nv+ 7) + zldarea
(  157)                 vinfor(nv+ 9) = vinfor(nv+ 9) + zicevol
(  158)                 vinfor(nv+11) = vinfor(nv+11) + zsnwvol
(  159)                 vinfor(nv+13) = vinfor(nv+13) + zicespd
(  160)              ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 10

(  161)           END DO
(  162)        END DO
(  163)        vinfor(nv+13) = SQRT( vinfor(nv+13) / MAX( vinfor(nv+9) , epsi06 ) )    
(  164) 
(  165)        !  Accumulation before averaging 
(  166)        DO jv = 1, nvinfo
(  167)           vinfom(jv) = vinfom(jv) + vinfor(jv)
(  168)        END DO
(  169)        naveg = naveg + 1  
(  170)     
(  171)        ! oututs on file ice_evolu    
(  172)        IF( MOD( numit , ninfo ) == 0 ) THEN
(  173)           WRITE(numevo_ice,fmtw) ( titvar(jv), vinfom(jv)/naveg, jv = 1, nvinfo )
(  174)           naveg = 0
(  175)           DO jv = 1, nvinfo
(  176)              vinfom(jv) = 0.e0
(  177)           END DO
(  178)        ENDIF
(  179)   
(  180)     END SUBROUTINE lim_dia
(  181)  
(  182) 
(  183)     SUBROUTINE lim_dia_init
(  184)        !!-------------------------------------------------------------------
(  185)        !!                  ***  ROUTINE lim_dia_init  ***
(  186)        !!             
(  187)        !! ** Purpose : Preparation of the file ice_evolu for the output of
(  188)        !!      the temporal evolution of key variables
(  189)        !!
(  190)        !! ** input   : Namelist namicedia
(  191)        !!
(  192)        !! history :
(  193)        !!  8.5  ! 03-08 (C. Ethe) original code
(  194)        !!-------------------------------------------------------------------
(  195)        NAMELIST/namicedia/fmtinf, nfrinf, ninfo, ntmoy
(  196) 
(  197)        INTEGER  ::   jv   ,     &  ! dummy loop indice
(  198)           &          ntot ,     &
(  199)           &          ndeb ,     &
(  200)           &          irecl
(  201) 
(  202)        INTEGER  ::   nv            ! indice of variable 
(  203) 
(  204)        REAL(wp) ::   zxx0, zxx1    ! temporary scalars
(  205) 
(  206)        CHARACTER(len=jpchinf) ::   titinf
(  207)        !!-------------------------------------------------------------------
(  208) 
(  209)        ! Read Namelist namicedia
(  210)        REWIND ( numnam_ice )
(  211)        READ   ( numnam_ice  , namicedia )
(  212)        IF(lwp) THEN
(  213)           WRITE(numout,*)
(  214)           WRITE(numout,*) 'lim_dia_init : ice parameters for ice diagnostics '
(  215)           WRITE(numout,*) '~~~~~~~~~~~~'
(  216)           WRITE(numout,*) '   format of the output values                                 fmtinf = ', fmtinf
(  217)           WRITE(numout,*) '   number of variables written in one line                     nfrinf = ', nfrinf 
(  218)           WRITE(numout,*) '   Instantaneous values of ice evolution or averaging          ntmoy  = ', ntmoy






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 11

(  219)           WRITE(numout,*) '   frequency of ouputs on file ice_evolu in case of averaging  ninfo  = ', ninfo
(  220)        ENDIF
(  221) 
(  222)        ! masked grid cell area
(  223)        aire(:,:) = area(:,:) * tms(:,:)
(  224) 
(  225)        ! Titles of ice key variables :
(  226)        nv = 1
(  227)        titvar(nv) = 'NoIt'  ! iteration number
(  228)        nv = nv + 1
(  229)        titvar(nv) = 'T yr'  ! time step in years
(  230)        nv = nv + 1
(  231) 
(  232)        nbvt = nv - 1
(  233) 
(  234)        titvar(nv) = 'AEFN' ! sea ice area in the northern Hemisp.(10^12 km2)
(  235)        nv = nv + 1
(  236)        titvar(nv) = 'AEFS' ! sea ice area in the southern Hemisp.(10^12 km2)
(  237)        nv = nv + 1
(  238)        titvar(nv) = 'A15N'  ! sea ice extent (15%) in the northern Hemisp.(10^12 km2)
(  239)        nv = nv + 1
(  240)        titvar(nv) = 'A15S'  ! sea ice extent (15%) in the southern Hemisp.(10^12 km2)
(  241)        nv = nv + 1
(  242)        titvar(nv) = 'A85N'  ! sea ice extent (85%) in the northern Hemisp.(10^12 km2)
(  243)        nv = nv + 1
(  244)        titvar(nv) = 'A85S'  ! sea ice extent (85%) in the southern Hemisp.(10^12 km2)
(  245)        nv = nv + 1
(  246)        titvar(nv) = 'ALEN'  ! leads area in the northern Hemisp.(10^12 km2)
(  247)        nv = nv + 1
(  248)        titvar(nv) = 'ALES'  ! leads area in the southern Hemisp.(10^12 km2)
(  249)        nv = nv + 1
(  250)        titvar(nv) = 'VOLN'  ! sea ice volume in the northern Hemisp.(10^3 km3)
(  251)        nv = nv + 1
(  252)        titvar(nv) = 'VOLS'  ! sea ice volume in the southern Hemisp.(10^3 km3)
(  253)        nv = nv + 1
(  254)        titvar(nv) = 'VONN'  ! snow volume over sea ice in the northern Hemisp.(10^3 km3)
(  255)        nv = nv + 1
(  256)        titvar(nv) = 'VONS'  ! snow volume over sea ice in the southern Hemisp.(10^3 km3)
(  257)        nv = nv + 1
(  258)        titvar(nv) = 'ECGN'  ! mean sea ice velocity in the northern Hemisp.(m/s)
(  259)        nv = nv + 1
(  260)        titvar(nv) = 'ECGS'  ! mean sea ice velocity in the southern Hemisp.(m/s)
(  261) 
(  262)        nvinfo = nv
(  263) 
(  264)        ! Definition et Ecriture de l'entete : nombre d'enregistrements 
(  265)        ndeb   = ( nstart - 1 ) / ninfo
(  266)        IF( nstart == 1 ) ndeb = -1
(  267) 
(  268)        nferme = ( nstart - 1 + nitrun) / ninfo
(  269)        ntot   = nferme - ndeb
(  270)        ndeb   = ninfo * ( 1 + ndeb )
(  271)        nferme = ninfo * nferme
(  272) 
(  273)        ! definition of formats 
(  274)        WRITE( fmtw  , '(A,I3,A2,I1,A)' )  '(', nfrinf, '(A', jpchsep, ','//fmtinf//'))'
(  275)        WRITE( fmtr  , '(A,I3,A,I1,A)'  )  '(', nfrinf, '(', jpchsep, 'X,'//fmtinf//'))'
(  276)        WRITE( fmtitr, '(A,I3,A,I1,A)'  )  '(', nvinfo, 'A', jpchinf, ')'






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 12

(  277) 
(  278)        ! opening  "ice_evolu" file
(  279)        irecl = ( jpchinf + 1 ) * nvinfo 
(  280)        OPEN( numevo_ice, file='ice.evolu', status='unknown', RECL = irecl)
(  281)        OPEN( numevo_ice, file='ice.evolu', status='unknown')
(  282) 
(  283)        !- ecriture de 2 lignes d''entete :
(  284)        WRITE(numevo_ice,1000) fmtr, fmtw, fmtitr, nvinfo, ntot, 0, nfrinf
(  285)        zxx0 = 0.001 * REAL( ninfo )
(  286)        zxx1 = 0.001 * REAL( ndeb  )
(  287)        WRITE(numevo_ice,1111) REAL(jpchinf), 0., zxx1, zxx0, 0., 0., 0
(  288) 
(  289)        !- ecriture de 2 lignes de titre :
(  290)        WRITE(numevo_ice,'(A,I8,A,I8,A,I5)')                                      &
(  291)           'Evolution chronologique - Experience '//cexper   &
(  292)           //'   de', ndeb, ' a', nferme, ' pas', ninfo
(  293)        WRITE(numevo_ice,fmtitr) ( titvar(jv), jv = 1, nvinfo )
(  294) 
(  295) 
(  296)        !--preparation de "titvar" pour l''ecriture parmi les valeurs numeriques :
(  297)        DO  jv = 2 , nvinfo
(  298)           titinf     = titvar(jv)(:jpchinf)
(  299)           titvar(jv) = '  '//titinf
(  300)        END DO
(  301) 
(  302)        !--Initialisation of the arrays for the accumulation
(  303)        DO  jv = 1, nvinfo
(  304)           vinfom(jv) = 0.
(  305)        END DO
(  306)        naveg = 0
(  307) 
(  308) 1000   FORMAT( 3(A20),4(1x,I6) )
(  309) 1111   FORMAT( 3(F7.1,1X,F7.3,1X),I3,A )  
(  310) 
(  311)     END SUBROUTINE lim_dia_init
(  312) 
(  313) 
(  314) # 322
(  322)    !!======================================================================
(  323) END MODULE limdia
