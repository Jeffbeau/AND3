


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

Filename: ldfeiv.F90

(    1) MODULE ldfeiv
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ldfeiv  ***
(    4)    !! Ocean physics:  variable eddy induced velocity coefficients
(    5)    !!======================================================================
(    6) #if   defined key_traldf_eiv   &&   defined key_traldf_c2d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_traldf_eiv'      and                     eddy induced velocity
(    9)    !!   'key_traldf_c2d'                    2D tracer lateral  mixing coef.
(   10)    !!----------------------------------------------------------------------
(   11)    !!   ldf_eiv      : compute the eddy induced velocity coefficients
(   12)    !!                  Same results but not same routine if 'key_autotasking'
(   13)    !!                  is defined or not
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE ldftra_oce      ! ocean tracer   lateral physics
(   19)    USE phycst          ! physical constants
(   20)    USE ldfslp          ! iso-neutral slopes
(   21)    USE flxrnf          ! 
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28)    
(   29)    !! * Routine accessibility
(   30)    PUBLIC ldf_eiv               ! routine called by step.F90
(   31)    !!----------------------------------------------------------------------
(   32)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfeiv.F90,v 1.8 2006/04/28 12:24:20 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36)    !! * Substitutions
(   37) #  include "domzgr_substitute.h90"
(   38) #  include "vectopt_loop_substitute.h90"
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43) # if defined key_autotasking
(   44)    !!----------------------------------------------------------------------
(   45)    !!   'key_autotasking' :                            autotasking (j-slab)
(   46)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 2

(   47) 
(   48)    SUBROUTINE ldf_eiv( kt )
(   49)       !!----------------------------------------------------------------------
(   50)       !!                  ***  ROUTINE ldf_eiv  ***
(   51)       !!
(   52)       !! ** Purpose :   Compute the eddy induced velocity coefficient from the
(   53)       !!      growth rate of baroclinic instability.
(   54)       !!
(   55)       !! ** Method :
(   56)       !!
(   57)       !! ** Action :   uslp(),   : i- and j-slopes of neutral surfaces
(   58)       !!               vslp()      at u- and v-points, resp.
(   59)       !!               wslpi(),  : i- and j-slopes of neutral surfaces
(   60)       !!               wslpj()     at w-points. 
(   61)       !!
(   62)       !! History :
(   63)       !!   8.1  !  99-03  (G. Madec, A. Jouzeau)  Original code
(   64)       !!   8.5  !  02-06  (G. Madec)  Free form, F90
(   65)       !!----------------------------------------------------------------------
(   66)       !! * Arguments
(   67)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx
(   68)       
(   69)       !! * Local declarations
(   70)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(   71)       REAL(wp) ::   &
(   72)          zfw, ze3w, zn2, zf20,       &  ! temporary scalars
(   73)          zaht, zaht_min
(   74)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   75)          zn, zah, zhw, zross            ! workspace
(   76)       !!----------------------------------------------------------------------
(   77) 
(   78)       IF( kt == nit000 ) THEN
(   79)          IF(lwp) WRITE(numout,*)
(   80)          IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
(   81)          IF(lwp) WRITE(numout,*) '~~~~~~~   key_autotasking'
(   82)       ENDIF
(   83)       
(   84)       !                                                ! ===============
(   85)       DO jj = 2, jpjm1                                 !  Vertical slab
(   86)          !                                             ! ===============
(   87)          
(   88)          ! 0. Local initialization
(   89)          ! -----------------------
(   90)          zn   (:,jj) = 0.e0
(   91)          zhw  (:,jj) = 5.e0
(   92)          zah  (:,jj) = 0.e0
(   93)          zross(:,jj) = 0.e0
(   94)          
(   95)          ! 1. Compute lateral diffusive coefficient 
(   96)          ! ----------------------------------------
(   97) 
(   98) !CDIR NOVERRCHK 
(   99)          DO jk = 1, jpk
(  100) !CDIR NOVERRCHK 
(  101)             DO ji = 2, jpim1
(  102)                ! Take the max of N^2 and zero then take the vertical sum 
(  103)                ! of the square root of the resulting N^2 ( required to compute 
(  104)                ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 3

(  105)                zn2 = MAX( rn2(ji,jj,jk), 0.e0 )
(  106)                ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
(  107)                zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
(  108)                ! Compute elements required for the inverse time scale of baroclinic
(  109)                ! eddies using the isopycnal slopes calculated in ldfslp.F : 
(  110)                ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
(  111)                zah(ji,jj) = zah(ji,jj) + zn2   &
(  112)                               * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)    &
(  113)                                 + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) )   &
(  114)                               * ze3w
(  115)                zhw(ji,jj) = zhw(ji,jj) + ze3w
(  116)             END DO 
(  117)          END DO 
(  118)  
(  119) !CDIR NOVERRCHK 
(  120)          DO ji = 2, jpim1
(  121)             zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
(  122)             ! Rossby radius at w-point taken < 40km and  > 2km
(  123)             zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
(  124)             ! Compute aeiw by multiplying Ro^2 and T^-1
(  125)             aeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
(  126)             ! Take the minimum between aeiw and 1000m^2/s for depth levels
(  127)             ! lower than 20 (21 in w- point)
(  128)             IF( mbathy(ji,jj) <= 21. ) aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
(  129)          END DO
(  130) 
(  131)          ! Decrease the coefficient in the tropics (20N-20S) 
(  132)          zf20 = 2. * omega * sin( rad * 20. )
(  133)          DO ji = 2, jpim1
(  134)             aeiw(ji,jj) = MIN( 1., ABS( ff(ji,jj) / zf20 ) ) * aeiw(ji,jj)
(  135)          END DO
(  136)   
(  137)          ! ORCA R05: Take the minimum between aeiw  and 1000m2/s
(  138)          IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05
(  139)             DO ji = 2, jpim1
(  140)                aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
(  141)             END DO
(  142)          ENDIF
(  143)          !                                             ! ===============
(  144)       END DO                                           !   End of slab
(  145)       !                                                ! ===============
(  146) 
(  147)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  148) 
(  149)       ! lateral boundary condition on aeiw 
(  150)       CALL lbc_lnk( aeiw, 'W', 1. )
(  151) 
(  152)       ! Average the diffusive coefficient at u- v- points 
(  153)       DO jj = 2, jpjm1
(  154)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  155)             aeiu(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji+1,jj  ))
(  156)             aeiv(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji  ,jj+1))
(  157)          END DO 
(  158)       END DO 
(  159)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  160) 
(  161)       ! lateral boundary condition on aeiu, aeiv 
(  162)       CALL lbc_lnk( aeiu, 'U', 1. )






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 4

(  163)       CALL lbc_lnk( aeiv, 'V', 1. )
(  164) 
(  165)       IF(ln_ctl)   THEN
(  166)          CALL prt_ctl(tab2d_1=aeiu, clinfo1=' eiv  - u: ', ovlap=1)
(  167)          CALL prt_ctl(tab2d_1=aeiv, clinfo1=' eiv  - v: ', ovlap=1)
(  168)       ENDIF
(  169)       
(  170)       ! ORCA R05: add a space variation on aht (=aeiv except at the equator and river mouth)
(  171)       IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
(  172)          zf20     = 2. * omega * SIN( rad * 20. )
(  173)          zaht_min = 100.                              ! minimum value for aht
(  174)          DO jj = 1, jpj
(  175)             DO ji = 1, jpi
(  176)                zaht      = ( 1. -  MIN( 1., ABS( ff(ji,jj) / zf20 ) ) ) * ( aht0 - zaht_min )  &
(  177)                   &      + aht0 * upsrnfh(ji,jj)                          ! enhanced near river mouths
(  178)                ahtu(ji,jj) = MAX( MAX( zaht_min, aeiu(ji,jj) ) + zaht, aht0 )
(  179)                ahtv(ji,jj) = MAX( MAX( zaht_min, aeiv(ji,jj) ) + zaht, aht0 )
(  180)                ahtw(ji,jj) = MAX( MAX( zaht_min, aeiw(ji,jj) ) + zaht, aht0 )
(  181)             END DO
(  182)          END DO
(  183)          IF(ln_ctl) THEN
(  184)             CALL prt_ctl(tab2d_1=ahtu, clinfo1=' aht  - u: ', ovlap=1)
(  185)             CALL prt_ctl(tab2d_1=ahtv, clinfo1=' aht  - v: ', ovlap=1)
(  186)             CALL prt_ctl(tab2d_1=ahtw, clinfo1=' aht  - w: ', ovlap=1)
(  187)          ENDIF
(  188)       ENDIF
(  189) 
(  190)       IF( aeiv0 == 0.e0 ) THEN
(  191)          aeiu(:,:) = 0.e0
(  192)          aeiv(:,:) = 0.e0
(  193)          aeiw(:,:) = 0.e0
(  194)       ENDIF
(  195) 
(  196)    END SUBROUTINE ldf_eiv
(  197) 
(  198) # else
(  199)    !!----------------------------------------------------------------------
(  200)    !!   Default key                                             k-j-i loops
(  201)    !!----------------------------------------------------------------------
(  202) 
(  203)    SUBROUTINE ldf_eiv( kt )
(  204)       !!----------------------------------------------------------------------
(  205)       !!                  ***  ROUTINE ldf_eiv  ***
(  206)       !!
(  207)       !! ** Purpose :   Compute the eddy induced velocity coefficient from the
(  208)       !!      growth rate of baroclinic instability.
(  209)       !!
(  210)       !! ** Method :
(  211)       !!
(  212)       !! ** Action : - uslp(),  : i- and j-slopes of neutral surfaces
(  213)       !!             - vslp()      at u- and v-points, resp.
(  214)       !!             - wslpi(),  : i- and j-slopes of neutral surfaces
(  215)       !!             - wslpj()     at w-points. 
(  216)       !!
(  217)       !! History :
(  218)       !!   8.1  !  99-03  (G. Madec, A. Jouzeau)  Original code
(  219)       !!   8.5  !  02-06  (G. Madec)  Free form, F90
(  220)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 5

(  221)       !! * Arguments
(  222)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx
(  223)       
(  224)       !! * Local declarations
(  225)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  226)       REAL(wp) ::   &
(  227)          zfw, ze3w, zn2, zf20,       &  ! temporary scalars
(  228)          zaht, zaht_min
(  229)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  230)          zn, zah, zhw, zross            ! workspace
(  231)       !!----------------------------------------------------------------------
(  232)       
(  233)       IF( kt == nit000 ) THEN
(  234)          IF(lwp) WRITE(numout,*)
(  235)          IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
(  236)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  237)       ENDIF
(  238)       
(  239)       ! 0. Local initialization
(  240)       ! -----------------------
(  241)       zn   (:,:) = 0.e0
(  242)       zhw  (:,:) = 5.e0
(  243)       zah  (:,:) = 0.e0
(  244)       zross(:,:) = 0.e0
(  245) 
(  246) 
(  247)       ! 1. Compute lateral diffusive coefficient 
(  248)       ! ----------------------------------------
(  249) 
(  250)       DO jk = 1, jpk
(  251) #  if defined key_vectopt_loop  &&  ! defined key_autotasking
(  252) !CDIR NOVERRCHK 
(  253)          DO ji = 1, jpij   ! vector opt.
(  254)             ! Take the max of N^2 and zero then take the vertical sum
(  255)             ! of the square root of the resulting N^2 ( required to compute
(  256)             ! internal Rossby radius Ro = .5 * sum_jpk(N) / f
(  257)             zn2 = MAX( rn2(ji,1,jk), 0.e0 )
(  258)             zn(ji,1) = zn(ji,1) + SQRT( zn2 ) * fse3w(ji,1,jk)
(  259)             ! Compute elements required for the inverse time scale of baroclinic
(  260)             ! eddies using the isopycnal slopes calculated in ldfslp.F :
(  261)             ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
(  262)             ze3w = fse3w(ji,1,jk) * tmask(ji,1,jk)
(  263)                zah(ji,1) = zah(ji,1) + zn2   &
(  264)                               * ( wslpi(ji,1,jk) * wslpi(ji,1,jk)    &
(  265)                                 + wslpj(ji,1,jk) * wslpj(ji,1,jk) )   &
(  266)                               * ze3w
(  267)             zhw(ji,1) = zhw(ji,1) + ze3w
(  268)          END DO
(  269) #  else
(  270)          DO jj = 2, jpjm1
(  271) !CDIR NOVERRCHK 
(  272)             DO ji = 2, jpim1
(  273)                ! Take the max of N^2 and zero then take the vertical sum 
(  274)                ! of the square root of the resulting N^2 ( required to compute 
(  275)                ! internal Rossby radius Ro = .5 * sum_jpk(N) / f 
(  276)                zn2 = MAX( rn2(ji,jj,jk), 0.e0 )
(  277)                zn(ji,jj) = zn(ji,jj) + SQRT( zn2 ) * fse3w(ji,jj,jk)
(  278)                ! Compute elements required for the inverse time scale of baroclinic






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 6

(  279)                ! eddies using the isopycnal slopes calculated in ldfslp.F : 
(  280)                ! T^-1 = sqrt(m_jpk(N^2*(r1^2+r2^2)*e3w))
(  281)                ze3w = fse3w(ji,jj,jk) * tmask(ji,jj,jk)
(  282)                zah(ji,jj) = zah(ji,jj) + zn2   &
(  283)                               * ( wslpi(ji,jj,jk) * wslpi(ji,jj,jk)    &
(  284)                                 + wslpj(ji,jj,jk) * wslpj(ji,jj,jk) )  &
(  285)                               * ze3w
(  286)                zhw(ji,jj) = zhw(ji,jj) + ze3w
(  287)             END DO 
(  288)          END DO 
(  289) #  endif
(  290)       END DO 
(  291) 
(  292)       DO jj = 2, jpjm1
(  293) !CDIR NOVERRCHK 
(  294)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  295)             zfw = MAX( ABS( 2. * omega * SIN( rad * gphit(ji,jj) ) ) , 1.e-10 )
(  296)             ! Rossby radius at w-point taken < 40km and  > 2km
(  297)             zross(ji,jj) = MAX( MIN( .4 * zn(ji,jj) / zfw, 40.e3 ), 2.e3 )
(  298)             ! Compute aeiw by multiplying Ro^2 and T^-1
(  299)             aeiw(ji,jj) = zross(ji,jj) * zross(ji,jj) * SQRT( zah(ji,jj) / zhw(ji,jj) ) * tmask(ji,jj,1)
(  300)             ! Take the minimum between aeiw and 1000m^2/s for depth levels
(  301)             ! lower than 20 (21 in w- point)
(  302)             IF( mbathy(ji,jj) <= 21. ) aeiw(ji,jj) = MIN( aeiw(ji,jj), 1000. )
(  303)          END DO
(  304)       END DO
(  305) 
(  306)       ! Decrease the coefficient in the tropics (20N-20S) 
(  307)          zf20 = 2. * omega * sin( rad * 20. )
(  308)       DO jj = 2, jpjm1
(  309)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  310)             aeiw(ji,jj) = MIN( 1., ABS( ff(ji,jj) / zf20 ) ) * aeiw(ji,jj)
(  311)          END DO
(  312)       END DO
(  313) 
(  314)       ! ORCA R05: Take the minimum between aeiw  and 1000m2/s
(  315)       IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
(  316)          DO jj = 2, jpjm1
(  317)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  318)                aeiw(ji,jj) = MIN( aeiw(ji,jj), aeiv0 )
(  319)             END DO
(  320)          END DO
(  321)       ENDIF
(  322) 
(  323)       ! lateral boundary condition on aeiw 
(  324)       CALL lbc_lnk( aeiw, 'W', 1. )
(  325) 
(  326)       ! Average the diffusive coefficient at u- v- points 
(  327)       DO jj = 2, jpjm1
(  328)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  329)             aeiu(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji+1,jj  ) )
(  330)             aeiv(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji  ,jj+1) )
(  331)          END DO 
(  332)       END DO 
(  333) 
(  334)       ! lateral boundary condition on aeiu, aeiv
(  335)       CALL lbc_lnk( aeiu, 'U', 1. )
(  336)       CALL lbc_lnk( aeiv, 'V', 1. )






PGF90 (Version     10.5)          09/23/2010  09:50:04      page 7

(  337) 
(  338)       IF(ln_ctl)   THEN
(  339)          CALL prt_ctl(tab2d_1=aeiu, clinfo1=' eiv  - u: ', ovlap=1)
(  340)          CALL prt_ctl(tab2d_1=aeiv, clinfo1=' eiv  - v: ', ovlap=1)
(  341)       ENDIF
(  342) 
(  343)       ! ORCA R05: add a space variation on aht (=aeiv except at the equator and river mouth)
(  344)       IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN
(  345)          zf20     = 2. * omega * SIN( rad * 20. )
(  346)          zaht_min = 100.                              ! minimum value for aht
(  347)          DO jj = 1, jpj
(  348)             DO ji = 1, jpi
(  349)                zaht      = ( 1. -  MIN( 1., ABS( ff(ji,jj) / zf20 ) ) ) * ( aht0 - zaht_min )  &
(  350)                   &      + aht0 * upsrnfh(ji,jj)                          ! enhanced near river mouths
(  351)                ahtu(ji,jj) = MAX( MAX( zaht_min, aeiu(ji,jj) ) + zaht, aht0 )
(  352)                ahtv(ji,jj) = MAX( MAX( zaht_min, aeiv(ji,jj) ) + zaht, aht0 )
(  353)                ahtw(ji,jj) = MAX( MAX( zaht_min, aeiw(ji,jj) ) + zaht, aht0 )
(  354)             END DO
(  355)          END DO
(  356)          IF(ln_ctl) THEN
(  357)             CALL prt_ctl(tab2d_1=ahtu, clinfo1=' aht  - u: ', ovlap=1)
(  358)             CALL prt_ctl(tab2d_1=ahtv, clinfo1=' aht  - v: ', ovlap=1)
(  359)             CALL prt_ctl(tab2d_1=ahtw, clinfo1=' aht  - w: ', ovlap=1)
(  360)          ENDIF
(  361)       ENDIF
(  362)       
(  363)       IF( aeiv0 == 0.e0 ) THEN
(  364)          aeiu(:,:) = 0.e0
(  365)          aeiv(:,:) = 0.e0
(  366)          aeiw(:,:) = 0.e0
(  367)       ENDIF
(  368) 
(  369)    END SUBROUTINE ldf_eiv
(  370) 
(  371) # endif
(  372) 
(  373) #else
(  374)    !!----------------------------------------------------------------------
(  375)    !!   Default option                                         Dummy module
(  376)    !!----------------------------------------------------------------------
(  377)    USE in_out_manager
(  378) CONTAINS
(  379)    SUBROUTINE ldf_eiv( kt )       ! Empty routine
(  380)       if(lwp)WRITE(numout,*) 'ldf_eiv: You should not have seen this print! error?', kt
(  381)    END SUBROUTINE ldf_eiv
(  382) #endif
(  383) 
(  384)    !!======================================================================
(  385) END MODULE ldfeiv















PGF90 (Version     10.5)          09/23/2010  09:50:04      page 8

(    1) # 1 "ldfeiv.F90"
(    1) MODULE ldfeiv
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ldfeiv  ***
(    4)    !! Ocean physics:  variable eddy induced velocity coefficients
(    5)    !!======================================================================
(    6) # 374
(  374)    !!----------------------------------------------------------------------
(  375)    !!   Default option                                         Dummy module
(  376)    !!----------------------------------------------------------------------
(  377)    USE in_out_manager
(  378) CONTAINS
(  379)    SUBROUTINE ldf_eiv( kt )       ! Empty routine
(  380)       if(lwp)WRITE(numout,*) 'ldf_eiv: You should not have seen this print! error?', kt
(  381)    END SUBROUTINE ldf_eiv
(  382) 
(  383) # 384
(  384)    !!======================================================================
(  385) END MODULE ldfeiv
