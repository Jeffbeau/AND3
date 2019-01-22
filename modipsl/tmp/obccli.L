


PGF90 (Version     10.5)          09/23/2010  09:49:24      page 1

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

Filename: obccli.F90

(    1) MODULE obccli
(    2)    !!===================================================================================
(    3)    !!                       ***  MODULE  obccli  ***
(    4)    !! Ocean dynamics:   Baroclinic componant of velocities on each open boundary
(    5)    !!===================================================================================
(    6) #if defined key_obc && defined key_dynspg_rl
(    7)    !!-----------------------------------------------------------------------------------
(    8)    !!   'key_obc'               and 
(    9)    !!   'key_dynspg_rl'
(   10)    !!-----------------------------------------------------------------------------------
(   11)    !!   obc_cli_dyn : Compute the baroclinic componant after the radiation phase
(   12)    !!   obc_cli_dta : Compute the baroclinic componant for the climatological velocities
(   13)    !!-----------------------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers   
(   16)    USE dom_oce         ! ocean space and time domain 
(   17)    USE phycst          ! physical constants
(   18)    USE obc_oce         ! ocean open boundary conditions
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC obc_cli    ! routine called in obcdyn.F90 and obcdta.F90 (rigid lid case)
(   25) 
(   26)    INTERFACE obc_cli
(   27)      MODULE PROCEDURE obc_cli_dyn, obc_cli_dta
(   28)    END INTERFACE
(   29) 
(   30)    !! * Substitutions
(   31) #  include "domzgr_substitute.h90"
(   32)    !!-----------------------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obccli.F90,v 1.3 2005/12/28 09:25:07 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!-----------------------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE obc_cli_dyn( obvel, velcli, obd, obf, obtyp, obl)
(   41)       !!--------------------------------------------------------------------------------
(   42)       !!                 ***  SUBROUTINE obc_cli_dyn  ***
(   43)       !!                   
(   44)       !! ** Purpose :   Compute the baroclinic velocities at the open boundaries.
(   45)       !!
(   46)       !! ** Method  :






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 2

(   47)       !!      - Compute de barotropic velocity along the considered Open Boundary 
(   48)       !!        and substract it to the total velocity to have baroclinic velotity.
(   49)       !!      - obtyp must be set to | 0 when traiting an East or West OB 
(   50)       !!                             | 1 when traiting a North or South OB.
(   51)       !!      - obl is the lenght of the OB (jpi or jpj) 
(   52)       !!
(   53)       !! History :
(   54)       !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(   55)       !!--------------------------------------------------------------------------------
(   56)       !! * Arguments
(   57)       INTEGER, INTENT( in ) ::   & ! OB localization:jpieob or jpiwob for East or West 
(   58)          obd, obf,               & !                 jpjnob or jpjsob for North or South 
(   59)          obl,                    & ! Lenght of the Open Boundary
(   60)          obtyp                     ! Type of Open Boundary: zonal or Meridional 
(   61)       REAL(wp), DIMENSION(:,:), INTENT( out) ::   &
(   62)          velcli                    ! Baroclinic velocity calculated
(   63)       REAL(wp), DIMENSION(:,:,:), INTENT( in ) ::   &
(   64)          obvel                     ! ua or va velocities from obcdyn.F90 routine
(   65) 
(   66)       !! * Local declarations
(   67)       INTEGER ::   &   
(   68)          ji, jj, jk, jle, jol         ! loop indices  
(   69)       REAL(wp) ::   zcbl              ! Temporary Baroclinic velocity 
(   70)       REAL(wp), DIMENSION(obl) ::   & 
(   71)          zvelbtpe,                  & ! Barotropic velocity 
(   72)          zhinv                        ! Invert of the local depth 1/H
(   73)       REAL(wp), DIMENSION(obl,jpk) ::   &
(   74)          zmskob,                    & ! Velocity mask
(   75)          zvel                         ! 2D Local velocity on OB
(   76) # if defined key_partial_steps
(   77)       REAL(wp), DIMENSION(obl,jpk) ::   &
(   78)          ze3ob                        ! Vertical scale factor
(   79) # else
(   80)       REAL(wp), DIMENSION(jpk) ::   &
(   81)          ze3ob                        ! Vertical scale factor
(   82) # endif
(   83)       !!--------------------------------------------------------------------------------
(   84) 
(   85)       ! 0. Array initialization
(   86)       ! -----------------------
(   87) 
(   88)       zhinv(:) = 0.e0
(   89)       zmskob(:,:) = 0.e0
(   90)       zvel(:,:) = 0.e0
(   91) # if defined key_partial_steps
(   92)       ze3ob(:,:) = 0.e0
(   93) # else
(   94)       ze3ob(:) = 0.e0
(   95) # endif
(   96) 
(   97)       IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
(   98)          DO ji = obd, obf
(   99)             zhinv(:) = hur(ji,:)
(  100)             zmskob(:,:) = umask(ji,:,:)
(  101)             zvel(:,:) = obvel(ji,:,:)
(  102) # if defined key_partial_steps
(  103)             ze3ob(:,:) = fse3u(ji,:,:)
(  104) # else






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 3

(  105)             ze3ob(:) = fse3u(:,:,:)
(  106) # endif
(  107)          END DO
(  108)       ELSE                             ! Zonal Open Boundary ( North or South OB )
(  109)          DO jj = obd, obf
(  110)             zhinv(:) = hvr(:,jj)
(  111)             zmskob(:,:) = vmask(:,jj,:)
(  112)             zvel(:,:) = obvel(:,jj,:)
(  113) # if defined key_partial_steps
(  114)             ze3ob(:,:) = fse3v(:,jj,:)
(  115) # else
(  116)             ze3ob(:) = fse3v(:,:,:)
(  117) # endif
(  118)          END DO
(  119)       END IF
(  120) 
(  121)       zvelbtpe(:) = 0.e0
(  122) 
(  123)       ! 1. vertical sum
(  124)       ! ----------------
(  125) # if defined key_vectopt_loop
(  126) !CDIR NOLOOPCHG
(  127) # endif
(  128)       DO jol = obd, obf ! Vector opt.
(  129)          DO jk = 1, jpkm1
(  130)             DO jle = 1, obl
(  131)                zvelbtpe(jle) = zvelbtpe(jle) + zvel(jle,jk)*zmskob(jle,jk) &
(  132) # if defined key_partial_steps
(  133)                                            * ze3ob(jol,jle,jk)
(  134) # else
(  135)                                            * ze3ob(jk)
(  136) # endif
(  137)             END DO
(  138)          END DO
(  139)       END DO
(  140) 
(  141)       ! 2. divide by the depth
(  142)       ! -----------------------
(  143)       DO jle = 1, obl
(  144)          zvelbtpe(jle) = zvelbtpe(jle) * zhinv(jle) * zmskob(jle,1) 
(  145)       END DO
(  146) 
(  147)       ! 3. substract zvelbtpe to the total velocity
(  148)       !    and save the baroclinic velocity in velcli()
(  149)       ! ------------------------------------------------
(  150)       DO jk = 1, jpkm1
(  151)          DO jle = 1, obl
(  152)             zcbl = zvel(jle,jk) - zvelbtpe(jle)*zmskob(jle,jk)
(  153)             velcli(jle,jk) = zcbl * zmskob(jle,jk)
(  154)          END DO
(  155)       END DO
(  156) 
(  157)    END SUBROUTINE obc_cli_dyn
(  158) 
(  159) 
(  160)    SUBROUTINE obc_cli_dta( obvel, velcli, obd, obf, obtyp, obl, mpp )
(  161)       !!--------------------------------------------------------------------------------
(  162)       !!                 ***  SUBROUTINE obc_cli_dta  ***






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 4

(  163)       !!                   
(  164)       !! ** Purpose :
(  165)       !!      Compute the baroclinic velocities for the climatological velocities.
(  166)       !!
(  167)       !! ** Method  :
(  168)       !!      - Compute de barotropic velocity along the considered Open Boundary 
(  169)       !!        and substract it to the total velocity to have baroclinic velotity.
(  170)       !!      - obtyp must be set to | 0 when traiting an East or West OB 
(  171)       !!                             | 1 when traiting a North or South OB.
(  172)       !!      - obl is the lenght of the OB (jpi or jpj) 
(  173)       !!
(  174)       !! History :
(  175)       !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  176)       !!--------------------------------------------------------------------------------
(  177)       !! * Arguments
(  178)       INTEGER, INTENT( in ) ::   & ! OB localization: jpieob or jpiwob for East or West 
(  179)          obd, obf,               & !                  jpjnob or jpjsob for North or South 
(  180)          obl,                    & ! Lenght of the Open Boundary
(  181)          mpp,                    & ! MPP index
(  182)          obtyp                     ! Type of Open Boundary: zonal or Meridional 
(  183)       REAL(wp), INTENT( out), DIMENSION(:,:) ::   &
(  184)          velcli                    ! Baroclinic velocity calculated
(  185)       REAL(wp), INTENT( inout ), DIMENSION(:,:,:) ::   &
(  186)          obvel                     ! uXdta or vXdta climatological velocities from 
(  187)                                    ! obcdta.F90 routine
(  188) 
(  189)       !! * Local declarations
(  190)       INTEGER ::   &
(  191)          ji, jj, jk, jle, jol, ij     ! loop indices  
(  192)       REAL(wp), DIMENSION(obl) ::   & 
(  193)          zvelbtpe,                  & ! Barotropic velocity 
(  194)          zhinv                        ! Invert of the local depth 1/H
(  195)       REAL(wp), DIMENSION(obl,jpk) ::   &
(  196)          zmskob                       ! Velocity mask
(  197) # if defined key_partial_steps
(  198)       REAL(wp), DIMENSION(obl,jpk) ::   &
(  199)          ze3ob                        ! Vertical scale factor
(  200) # else
(  201)       REAL(wp), DIMENSION(jpk) ::   &
(  202)          ze3ob                        ! Vertical scale factor
(  203) # endif
(  204)       !!--------------------------------------------------------------------------------
(  205) 
(  206)       ! 0. Array initialization
(  207)       ! -----------------------
(  208) 
(  209)       zhinv(:) = 0.e0
(  210)       zmskob(:,:) = 0.e0
(  211) # if defined key_partial_steps
(  212)       ze3ob(:,:) = 0.e0
(  213) # else
(  214)       ze3ob(:) = 0.e0
(  215) # endif
(  216) 
(  217)       IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
(  218)          DO ji = obd, obf
(  219)             zhinv(:) = hur(ji,:)
(  220)             zmskob(:,:) = umask(ji,:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 5

(  221) # if defined key_partial_steps
(  222)             ze3ob(:,:) = fse3u(ji,:,:)
(  223) # else
(  224)             ze3ob(:) = fse3u(:,:,:)
(  225) # endif
(  226)          END DO
(  227)       ELSE                             ! Zonal Open Boundary ( North or South OB )
(  228)          DO jj = obd, obf
(  229)             zhinv(:) = hvr(:,jj)
(  230)             zmskob(:,:) = vmask(:,jj,:)
(  231) # if defined key_partial_steps
(  232)             ze3ob(:,:) = fse3v(:,jj,:)
(  233) # else
(  234)             ze3ob(:) = fse3v(:,:,:)
(  235) # endif
(  236)          END DO
(  237)       END IF
(  238) 
(  239)       zvelbtpe(:) = 0.e0
(  240) 
(  241)       ! 1. vertical sum
(  242)       ! ----------------
(  243) # if defined key_vectopt_loop
(  244) !CDIR NOLOOPCHG
(  245) # endif
(  246)       DO jol = obd, obf ! Vector opt.
(  247)          DO jk = 1, jpkm1
(  248)             DO jle = 1, obl
(  249)                ij = jle -1 + mpp
(  250)                zvelbtpe(jle) = zvelbtpe(jle) + obvel(ij,jk,1)*zmskob(jle,jk) &
(  251) # if defined key_partial_steps
(  252)                                            * ze3ob(jol,jle,jk)
(  253) # else
(  254)                                            * ze3ob(jk)
(  255) # endif
(  256)             END DO
(  257)          END DO
(  258)       END DO
(  259) 
(  260)       ! 2. divide by the depth
(  261)       ! -----------------------
(  262)       DO jle = 1, obl
(  263)          zvelbtpe(jle) = zvelbtpe(jle) * zhinv(jle) * zmskob(jle,1) 
(  264)       END DO 
(  265) 
(  266)       ! 3. substract zvelbtpe to the total velocity
(  267)       !    and save the baroclinic velocity in velcli()
(  268)       ! ------------------------------------------------
(  269)       DO jk = 1, jpkm1
(  270)          DO jle = 1, obl
(  271)             ij = jle -1 + mpp
(  272)             obvel(ij,jk,1) = obvel(ij,jk,1) - zvelbtpe(jle)*zmskob(jle,jk)
(  273)             velcli(jle,jk) = obvel(ij,jk,1) * zmskob(jle,jk)
(  274)          END DO
(  275)       END DO
(  276) 
(  277)    END SUBROUTINE obc_cli_dta
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:24      page 6

(  279) #else
(  280)    !!----------------------------------------------------------------------------------
(  281)    !!   Default options :                                                  Empty module
(  282)    !!----------------------------------------------------------------------------------
(  283) CONTAINS
(  284)    SUBROUTINE obc_cli_dyn       ! Empty routine
(  285)    END SUBROUTINE obc_cli_dyn
(  286)    SUBROUTINE obc_cli_dta       ! Empty routine
(  287)    END SUBROUTINE obc_cli_dta
(  288) #endif
(  289) 
(  290)    !!==================================================================================
(  291) END MODULE obccli



















































PGF90 (Version     10.5)          09/23/2010  09:49:24      page 7

(    1) # 1 "obccli.F90"
(    1) MODULE obccli
(    2)    !!===================================================================================
(    3)    !!                       ***  MODULE  obccli  ***
(    4)    !! Ocean dynamics:   Baroclinic componant of velocities on each open boundary
(    5)    !!===================================================================================
(    6) # 280
(  280)    !!----------------------------------------------------------------------------------
(  281)    !!   Default options :                                                  Empty module
(  282)    !!----------------------------------------------------------------------------------
(  283) CONTAINS
(  284)    SUBROUTINE obc_cli_dyn       ! Empty routine
(  285)    END SUBROUTINE obc_cli_dyn
(  286)    SUBROUTINE obc_cli_dta       ! Empty routine
(  287)    END SUBROUTINE obc_cli_dta
(  288) 
(  289) # 290
(  290)    !!==================================================================================
(  291) END MODULE obccli
