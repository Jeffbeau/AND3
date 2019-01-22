


PGF90 (Version     10.5)          09/23/2010  09:50:08      page 1

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

Filename: obcrad.F90

(    1) MODULE obcrad 
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcrad  ***
(    4)    !! Ocean dynamic :   Phase velocities for each open boundary
(    5)    !!=================================================================================
(    6) #if defined key_obc
(    7)    !!---------------------------------------------------------------------------------
(    8)    !!   obc_rad        : call the subroutine for each open boundary
(    9)    !!   obc_rad_east   : compute the east phase velocities
(   10)    !!   obc_rad_west   : compute the west phase velocities
(   11)    !!   obc_rad_north  : compute the north phase velocities
(   12)    !!   obc_rad_south  : compute the south phase velocities
(   13)    !!---------------------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE phycst          ! physical constants
(   19)    USE obc_oce         ! ocean open boundary conditions
(   20)    USE lib_mpp         ! for mppobc
(   21)    USE in_out_manager  ! I/O units
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Accessibility
(   27)    PUBLIC obc_rad        ! routine called by step.F90
(   28) 
(   29)    !! * Module variables
(   30)    INTEGER ::   ji, jj, jk     ! dummy loop indices
(   31) 
(   32)    INTEGER ::      & ! ... boundary space indices 
(   33)       nib   = 1,   & ! nib   = boundary point
(   34)       nibm  = 2,   & ! nibm  = 1st interior point
(   35)       nibm2 = 3,   & ! nibm2 = 2nd interior point
(   36)                      ! ... boundary time indices 
(   37)       nit   = 1,   & ! nit    = now
(   38)       nitm  = 2,   & ! nitm   = before
(   39)       nitm2 = 3      ! nitm2  = before-before
(   40) 
(   41)    !! * Substitutions
(   42) #  include "obc_vectopt_loop_substitute.h90"
(   43)    !!---------------------------------------------------------------------------------
(   44)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcrad.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 2

(   47)    !!---------------------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE obc_rad ( kt )
(   52)       !!------------------------------------------------------------------------------
(   53)       !!                     SUBROUTINE obc_rad
(   54)       !!                    ********************
(   55)       !! ** Purpose :
(   56)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(   57)       !!      boundaries and calculate those phase speeds if the open boundaries are 
(   58)       !!      not fixed. In case of fixed open boundaries does nothing.
(   59)       !!
(   60)       !!     The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   61)       !!     and/or lp_obc_south allow the user to determine which boundary is an
(   62)       !!     open one (must be done in the param_obc.h90 file).
(   63)       !! 
(   64)       !! ** Reference : 
(   65)       !!     Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   66)       !!
(   67)       !!  History :
(   68)       !!    8.5  !  02-10  (C. Talandier, A-M. Treguier) Free surface, F90 from the 
(   69)       !!                                                 J. Molines and G. Madec version
(   70)       !!------------------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt
(   73)       !!----------------------------------------------------------------------
(   74) 
(   75)       IF( lp_obc_east  .AND. .NOT.lfbceast  )   CALL obc_rad_east ( kt )   ! East open boundary
(   76) 
(   77)       IF( lp_obc_west  .AND. .NOT.lfbcwest  )   CALL obc_rad_west ( kt )   ! West open boundary
(   78) 
(   79)       IF( lp_obc_north .AND. .NOT.lfbcnorth )   CALL obc_rad_north( kt )   ! North open boundary
(   80) 
(   81)       IF( lp_obc_south .AND. .NOT.lfbcsouth )   CALL obc_rad_south( kt )   ! South open boundary
(   82) 
(   83)    END SUBROUTINE obc_rad
(   84) 
(   85) 
(   86)    SUBROUTINE obc_rad_east ( kt )
(   87)       !!------------------------------------------------------------------------------
(   88)       !!                     ***  SUBROUTINE obc_rad_east  ***
(   89)       !!                   
(   90)       !! ** Purpose :
(   91)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(   92)       !!      east boundary and calculate those phase speeds if this OBC is not fixed.
(   93)       !!      In case of fixed OBC, this subrountine is not called.
(   94)       !!
(   95)       !!  History :
(   96)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(   97)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(   98)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(   99)       !!         ! 00-06 (J.-M. Molines) 
(  100)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  101)       !!------------------------------------------------------------------------------
(  102)       !! * Arguments
(  103)       INTEGER, INTENT( in ) ::   kt
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 3

(  105)       !! * Local declarations
(  106)       INTEGER  ::   ij
(  107)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  108)       REAL(wp) ::   zucb, zucbm, zucbm2
(  109)       !!------------------------------------------------------------------------------
(  110) 
(  111)       ! 1. Swap arrays before calculating radiative velocities
(  112)       ! ------------------------------------------------------
(  113) 
(  114)       ! 1.1  zonal velocity 
(  115)       ! -------------------
(  116) 
(  117)       IF( kt > nit000 ) THEN 
(  118) 
(  119)          ! ... advance in time (time filter, array swap) 
(  120)          DO jk = 1, jpkm1
(  121)             DO jj = 1, jpj
(  122)                uebnd(jj,jk,nib  ,nitm2) = uebnd(jj,jk,nib  ,nitm)*uemsk(jj,jk)
(  123)                uebnd(jj,jk,nibm ,nitm2) = uebnd(jj,jk,nibm ,nitm)*uemsk(jj,jk)
(  124)                uebnd(jj,jk,nibm2,nitm2) = uebnd(jj,jk,nibm2,nitm)*uemsk(jj,jk)
(  125)             END DO
(  126)          END DO
(  127)          ! ... fields nitm <== nit  plus time filter at the boundary 
(  128)          DO ji = fs_nie0, fs_nie1 ! Vector opt.
(  129)             DO jk = 1, jpkm1
(  130)                DO jj = 1, jpj
(  131)                   uebnd(jj,jk,nib  ,nitm) = uebnd(jj,jk,nib,  nit)*uemsk(jj,jk)
(  132)                   uebnd(jj,jk,nibm ,nitm) = uebnd(jj,jk,nibm ,nit)*uemsk(jj,jk)
(  133)                   uebnd(jj,jk,nibm2,nitm) = uebnd(jj,jk,nibm2,nit)*uemsk(jj,jk)
(  134)          ! ... fields nit <== now (kt+1) 
(  135)          ! ... Total or baroclinic velocity at b, bm and bm2
(  136) # if defined key_dynspg_rl
(  137)                   zucb   = uclie(jj,jk)
(  138) # else
(  139)                   zucb   = un(ji,jj,jk)
(  140) # endif
(  141) # if defined key_dynspg_rl
(  142)                   zucbm  = un(ji-1,jj,jk) + hur(ji-1,jj) / e2u(ji-1,jj) &
(  143)                            * ( bsfn(ji-1,jj) - bsfn(ji-1,jj-1) )
(  144) # else
(  145)                   zucbm  = un(ji-1,jj,jk)
(  146) # endif
(  147) # if defined key_dynspg_rl
(  148)                   zucbm2 = un(ji-2,jj,jk) + hur(ji-2,jj) / e2u(ji-2,jj) &
(  149)                            * ( bsfn(ji-2,jj) - bsfn(ji-2,jj-1) )
(  150) # else
(  151)                   zucbm2 = un(ji-2,jj,jk)
(  152) # endif
(  153)                   uebnd(jj,jk,nib  ,nit) = zucb   *uemsk(jj,jk)
(  154)                   uebnd(jj,jk,nibm ,nit) = zucbm  *uemsk(jj,jk) 
(  155)                   uebnd(jj,jk,nibm2,nit) = zucbm2 *uemsk(jj,jk) 
(  156)                END DO
(  157)             END DO
(  158)          END DO
(  159)          IF( lk_mpp )   CALL mppobc(uebnd,jpjed,jpjef,jpieob,jpk*3*3,2,jpj)
(  160) 
(  161)          ! ... extremeties nie0, nie1
(  162)          ij = jpjed +1 - njmpp






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 4

(  163)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  164)             DO jk = 1,jpkm1
(  165)                uebnd(ij,jk,nibm,nitm) = uebnd(ij+1 ,jk,nibm,nitm)
(  166)             END DO
(  167)          END IF
(  168)          ij = jpjef +1 - njmpp
(  169)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  170)             DO jk = 1,jpkm1
(  171)                uebnd(ij,jk,nibm,nitm) = uebnd(ij-1 ,jk,nibm,nitm)
(  172)             END DO
(  173)          END IF
(  174) 
(  175)          ! 1.2 tangential velocity
(  176)          ! -----------------------
(  177) 
(  178)          ! ... advance in time (time filter, array swap)
(  179)          DO jk = 1, jpkm1
(  180)             DO jj = 1, jpj
(  181)          ! ... fields nitm2 <== nitm
(  182)                vebnd(jj,jk,nib  ,nitm2) = vebnd(jj,jk,nib  ,nitm)*vemsk(jj,jk)
(  183)                vebnd(jj,jk,nibm ,nitm2) = vebnd(jj,jk,nibm ,nitm)*vemsk(jj,jk)
(  184)                vebnd(jj,jk,nibm2,nitm2) = vebnd(jj,jk,nibm2,nitm)*vemsk(jj,jk)
(  185)             END DO
(  186)          END DO
(  187) 
(  188)          DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
(  189)             DO jk = 1, jpkm1
(  190)                DO jj = 1, jpj
(  191)                   vebnd(jj,jk,nib  ,nitm) = vebnd(jj,jk,nib,  nit)*vemsk(jj,jk)
(  192)                   vebnd(jj,jk,nibm ,nitm) = vebnd(jj,jk,nibm ,nit)*vemsk(jj,jk)
(  193)                   vebnd(jj,jk,nibm2,nitm) = vebnd(jj,jk,nibm2,nit)*vemsk(jj,jk)
(  194)          ! ... fields nit <== now (kt+1)
(  195)                   vebnd(jj,jk,nib  ,nit) = vn(ji  ,jj,jk)*vemsk(jj,jk)
(  196)                   vebnd(jj,jk,nibm ,nit) = vn(ji-1,jj,jk)*vemsk(jj,jk)
(  197)                   vebnd(jj,jk,nibm2,nit) = vn(ji-2,jj,jk)*vemsk(jj,jk)
(  198)                END DO
(  199)             END DO
(  200)          END DO
(  201)          IF( lk_mpp )   CALL mppobc(vebnd,jpjed,jpjef,jpieob+1,jpk*3*3,2,jpj)
(  202) 
(  203)          !... extremeties nie0, nie1
(  204)          ij = jpjed +1 - njmpp
(  205)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  206)             DO jk = 1,jpkm1
(  207)                vebnd(ij,jk,nibm,nitm) = vebnd(ij+1 ,jk,nibm,nitm)
(  208)             END DO 
(  209)          END IF 
(  210)          ij = jpjef +1 - njmpp 
(  211)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  212)             DO jk = 1,jpkm1 
(  213)                vebnd(ij,jk,nibm,nitm) = vebnd(ij-1 ,jk,nibm,nitm)
(  214)             END DO 
(  215)          END IF 
(  216) 
(  217)          ! 1.3 Temperature and salinity
(  218)          ! ----------------------------
(  219) 
(  220)          ! ... advance in time (time filter, array swap)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 5

(  221)          DO jk = 1, jpkm1
(  222)             DO jj = 1, jpj
(  223)          ! ... fields nitm <== nit  plus time filter at the boundary
(  224)                tebnd(jj,jk,nib,nitm) = tebnd(jj,jk,nib,nit)*temsk(jj,jk)
(  225)                sebnd(jj,jk,nib,nitm) = sebnd(jj,jk,nib,nit)*temsk(jj,jk)
(  226)             END DO
(  227)          END DO
(  228) 
(  229)          DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
(  230)             DO jk = 1, jpkm1
(  231)                DO jj = 1, jpj
(  232)                   tebnd(jj,jk,nibm,nitm) = tebnd(jj,jk,nibm,nit)*temsk(jj,jk)
(  233)                   sebnd(jj,jk,nibm,nitm) = sebnd(jj,jk,nibm,nit)*temsk(jj,jk)
(  234)          ! ... fields nit <== now (kt+1)
(  235)                   tebnd(jj,jk,nib  ,nit) = tn(ji  ,jj,jk)*temsk(jj,jk)
(  236)                   tebnd(jj,jk,nibm ,nit) = tn(ji-1,jj,jk)*temsk(jj,jk)
(  237)                   sebnd(jj,jk,nib  ,nit) = sn(ji  ,jj,jk)*temsk(jj,jk)
(  238)                   sebnd(jj,jk,nibm ,nit) = sn(ji-1,jj,jk)*temsk(jj,jk)
(  239)                END DO
(  240)             END DO
(  241)          END DO
(  242)          IF( lk_mpp )   CALL mppobc(tebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  243)          IF( lk_mpp )   CALL mppobc(sebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  244) 
(  245)          ! ... extremeties nie0, nie1
(  246)          ij = jpjed +1 - njmpp
(  247)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  248)             DO jk = 1,jpkm1
(  249)                tebnd(ij,jk,nibm,nitm) = tebnd(ij+1 ,jk,nibm,nitm)
(  250)                sebnd(ij,jk,nibm,nitm) = sebnd(ij+1 ,jk,nibm,nitm)
(  251)             END DO
(  252)          END IF
(  253)          ij = jpjef +1 - njmpp
(  254)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  255)             DO jk = 1,jpkm1
(  256)                tebnd(ij,jk,nibm,nitm) = tebnd(ij-1 ,jk,nibm,nitm)
(  257)                sebnd(ij,jk,nibm,nitm) = sebnd(ij-1 ,jk,nibm,nitm)
(  258)             END DO
(  259)          END IF
(  260) 
(  261)       END IF     ! End of array swap
(  262) 
(  263)       ! 2 - Calculation of radiation velocities
(  264)       ! ---------------------------------------
(  265) 
(  266)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  267) 
(  268)          ! 2.1  Calculate the normal velocity U based on phase velocity u_cxebnd
(  269)          ! ---------------------------------------------------------------------
(  270)          !
(  271)          !          nibm2      nibm      nib
(  272)          !            |  nibm   |   nib   |///
(  273)          !            |    |    |    |    |///
(  274)          !  jj-line --f----v----f----v----f---
(  275)          !            |    |    |    |    |///
(  276)          !            |         |         |///
(  277)          !  jj-line   u    T    u    T    u///
(  278)          !            |         |         |///






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 6

(  279)          !            |    |    |    |    |///
(  280)          !          jpieob-2   jpieob-1   jpieob
(  281)          !                 |         |        
(  282)          !              jpieob-1    jpieob      
(  283)          !
(  284)          ! ... (jpjedp1, jpjefm1),jpieob
(  285)          DO ji = fs_nie0, fs_nie1 ! Vector opt.
(  286)             DO jk = 1, jpkm1
(  287)                DO jj = 2, jpjm1
(  288)          ! ... 2* gradi(u) (T-point i=nibm, time mean)
(  289)                   z2dx = ( uebnd(jj,jk,nibm ,nit) + uebnd(jj,jk,nibm ,nitm2) &
(  290)                            - 2.*uebnd(jj,jk,nibm2,nitm) ) / e1t(ji-1,jj)
(  291)          ! ... 2* gradj(u) (u-point i=nibm, time nitm)
(  292)                   z2dy = ( uebnd(jj+1,jk,nibm,nitm) - uebnd(jj-1,jk,nibm,nitm) ) / e2u(ji-1,jj)
(  293)          ! ... square of the norm of grad(u)
(  294)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  295)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  296)                   zdt = uebnd(jj,jk,nibm,nitm2) - uebnd(jj,jk,nibm,nit)
(  297)          ! ... i-phase speed ratio (bounded by 1)               
(  298)                   IF( z4nor2 == 0. ) THEN
(  299)                      z4nor2=.00001
(  300)                   END IF
(  301)                   z05cx = zdt * z2dx / z4nor2
(  302)                   u_cxebnd(jj,jk) = z05cx*uemsk(jj,jk)
(  303)                END DO
(  304)             END DO
(  305)          END DO
(  306) 
(  307)          ! 2.2  Calculate the tangential velocity based on phase velocity v_cxebnd
(  308)          ! -----------------------------------------------------------------------
(  309)          !
(  310)          !          nibm2      nibm      nib
(  311)          !            |   nibm  |   nib///|///
(  312)          !            |    |    |    |////|///
(  313)          !  jj-line --v----f----v----f----v---
(  314)          !            |    |    |    |////|///
(  315)          !            |    |    |    |////|///
(  316)          !            | jpieob-1| jpieob /|///
(  317)          !            |         |         |   
(  318)          !         jpieob-1    jpieob     jpieob+1
(  319)          !
(  320)          ! ... (jpjedp1, jpjefm1), jpieob+1
(  321)          DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
(  322)             DO jk = 1, jpkm1
(  323)                DO jj = 2, jpjm1
(  324)          ! ... 2* i-gradient of v (f-point i=nibm, time mean)
(  325)                   z2dx = ( vebnd(jj,jk,nibm ,nit) + vebnd(jj,jk,nibm ,nitm2) &
(  326)                           - 2.*vebnd(jj,jk,nibm2,nitm) ) / e1f(ji-2,jj)
(  327)          ! ... 2* j-gradient of v (v-point i=nibm, time nitm)
(  328)                   z2dy = ( vebnd(jj+1,jk,nibm,nitm) -  vebnd(jj-1,jk,nibm,nitm) ) / e2v(ji-1,jj)
(  329)          ! ... square of the norm of grad(v)
(  330)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  331)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  332)                   zdt = vebnd(jj,jk,nibm,nitm2) - vebnd(jj,jk,nibm,nit)
(  333)          ! ... i-phase speed ratio (bounded by 1) and save the unbounded phase
(  334)          !     velocity ratio no divided by e1f for the tracer radiation
(  335)                   IF( z4nor2 == 0. ) THEN
(  336)                      z4nor2=.000001






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 7

(  337)                   END IF
(  338)                   z05cx = zdt * z2dx / z4nor2
(  339)                   v_cxebnd(jj,jk) = z05cx*vemsk(jj,jk)
(  340)                END DO
(  341)             END DO
(  342)          END DO
(  343)          IF( lk_mpp )   CALL mppobc(v_cxebnd,jpjed,jpjef,jpieob+1,jpk,2,jpj)
(  344) 
(  345)          ! ... extremeties nie0, nie1
(  346)          ij = jpjed +1 - njmpp
(  347)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  348)             DO jk = 1,jpkm1
(  349)                v_cxebnd(ij,jk) = v_cxebnd(ij+1 ,jk)
(  350)             END DO
(  351)          END IF
(  352)          ij = jpjef +1 - njmpp
(  353)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  354)             DO jk = 1,jpkm1
(  355)                v_cxebnd(ij,jk) = v_cxebnd(ij-1 ,jk)
(  356)             END DO
(  357)          END IF
(  358) 
(  359)       END IF
(  360) 
(  361)    END SUBROUTINE obc_rad_east
(  362) 
(  363) 
(  364)    SUBROUTINE obc_rad_west ( kt )
(  365)       !!------------------------------------------------------------------------------
(  366)       !!                  ***  SUBROUTINE obc_rad_west  ***
(  367)       !!                    
(  368)       !! ** Purpose :
(  369)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  370)       !!      west boundary and calculate those phase speeds if this OBC is not fixed.
(  371)       !!      In case of fixed OBC, this subrountine is not called.
(  372)       !!
(  373)       !!  History :
(  374)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  375)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  376)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  377)       !!         ! 00-06 (J.-M. Molines) 
(  378)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  379)       !!------------------------------------------------------------------------------
(  380)       !! * Arguments
(  381)       INTEGER, INTENT( in ) ::   kt
(  382) 
(  383)       !! * Local declarations
(  384)       INTEGER ::   ij
(  385)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  386)       REAL(wp) ::   zucb, zucbm, zucbm2
(  387)       !!------------------------------------------------------------------------------
(  388) 
(  389)       ! 1. Swap arrays before calculating radiative velocities
(  390)       ! ------------------------------------------------------
(  391) 
(  392)       ! 1.1  zonal velocity 
(  393)       ! -------------------
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 8

(  395)       IF( kt > nit000 ) THEN
(  396) 
(  397)          ! ... advance in time (time filter, array swap) 
(  398)          DO jk = 1, jpkm1
(  399)             DO jj = 1, jpj 
(  400)                uwbnd(jj,jk,nib  ,nitm2) = uwbnd(jj,jk,nib  ,nitm)*uwmsk(jj,jk)
(  401)                uwbnd(jj,jk,nibm ,nitm2) = uwbnd(jj,jk,nibm ,nitm)*uwmsk(jj,jk)
(  402)                uwbnd(jj,jk,nibm2,nitm2) = uwbnd(jj,jk,nibm2,nitm)*uwmsk(jj,jk)
(  403)             END DO
(  404)          END DO
(  405) 
(  406)          ! ... fields nitm <== nit  plus time filter at the boundary 
(  407)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  408)             DO jk = 1, jpkm1
(  409)                DO jj = 1, jpj
(  410)                   uwbnd(jj,jk,nib  ,nitm) = uwbnd(jj,jk,nib  ,nit)*uwmsk(jj,jk)
(  411)                   uwbnd(jj,jk,nibm ,nitm) = uwbnd(jj,jk,nibm ,nit)*uwmsk(jj,jk)
(  412)                   uwbnd(jj,jk,nibm2,nitm) = uwbnd(jj,jk,nibm2,nit)*uwmsk(jj,jk)
(  413)          ! ... total or baroclinic velocity at b, bm and bm2
(  414) # if defined key_dynspg_rl
(  415)                   zucb   = ucliw(jj,jk) 
(  416) # else
(  417)                   zucb   = un (ji,jj,jk)
(  418) # endif
(  419) # if defined key_dynspg_rl
(  420)                   zucbm  = un (ji+1,jj,jk) + hur (ji+1,jj) / e2u (ji+1,jj) &
(  421)                            * ( bsfn(ji+1,jj) - bsfn(ji+1,jj-1) )
(  422) # else
(  423)                   zucbm  = un (ji+1,jj,jk)
(  424) # endif
(  425) # if defined key_dynspg_rl
(  426)                   zucbm2 = un (ji+2,jj,jk) + hur (ji+2,jj) / e2u (ji+2,jj) &
(  427)                            * ( bsfn(ji+2,jj) - bsfn(ji+2,jj-1) )
(  428) # else
(  429)                   zucbm2 = un (ji+2,jj,jk)
(  430) # endif
(  431) 
(  432)          ! ... fields nit <== now (kt+1) 
(  433)                   uwbnd(jj,jk,nib  ,nit) = zucb  *uwmsk(jj,jk)
(  434)                   uwbnd(jj,jk,nibm ,nit) = zucbm *uwmsk(jj,jk)
(  435)                   uwbnd(jj,jk,nibm2,nit) = zucbm2*uwmsk(jj,jk)
(  436)                END DO
(  437)             END DO
(  438)          END DO
(  439)          IF( lk_mpp )   CALL mppobc(uwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  440) 
(  441)          ! ... extremeties niw0, niw1
(  442)          ij = jpjwd +1 - njmpp
(  443)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  444)             DO jk = 1,jpkm1
(  445)                uwbnd(ij,jk,nibm,nitm) = uwbnd(ij+1 ,jk,nibm,nitm)
(  446)             END DO
(  447)          END IF
(  448)          ij = jpjwf +1 - njmpp
(  449)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  450)             DO jk = 1,jpkm1
(  451)                uwbnd(ij,jk,nibm,nitm) = uwbnd(ij-1 ,jk,nibm,nitm)
(  452)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 9

(  453)          END IF
(  454) 
(  455)          ! 1.2 tangential velocity
(  456)          ! -----------------------
(  457) 
(  458)          ! ... advance in time (time filter, array swap)
(  459)          DO jk = 1, jpkm1
(  460)             DO jj = 1, jpj 
(  461)          ! ... fields nitm2 <== nitm
(  462)                   vwbnd(jj,jk,nib  ,nitm2) = vwbnd(jj,jk,nib  ,nitm)*vwmsk(jj,jk)
(  463)                   vwbnd(jj,jk,nibm ,nitm2) = vwbnd(jj,jk,nibm ,nitm)*vwmsk(jj,jk)
(  464)                   vwbnd(jj,jk,nibm2,nitm2) = vwbnd(jj,jk,nibm2,nitm)*vwmsk(jj,jk)
(  465)             END DO
(  466)          END DO
(  467) 
(  468)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  469)             DO jk = 1, jpkm1
(  470)                DO jj = 1, jpj
(  471)                   vwbnd(jj,jk,nib  ,nitm) = vwbnd(jj,jk,nib,  nit)*vwmsk(jj,jk)
(  472)                   vwbnd(jj,jk,nibm ,nitm) = vwbnd(jj,jk,nibm ,nit)*vwmsk(jj,jk)
(  473)                   vwbnd(jj,jk,nibm2,nitm) = vwbnd(jj,jk,nibm2,nit)*vwmsk(jj,jk)
(  474)          ! ... fields nit <== now (kt+1)
(  475)                   vwbnd(jj,jk,nib  ,nit) = vn(ji  ,jj,jk)*vwmsk(jj,jk)
(  476)                   vwbnd(jj,jk,nibm ,nit) = vn(ji+1,jj,jk)*vwmsk(jj,jk)
(  477)                   vwbnd(jj,jk,nibm2,nit) = vn(ji+2,jj,jk)*vwmsk(jj,jk)
(  478)                END DO
(  479)             END DO
(  480)          END DO
(  481)          IF( lk_mpp )   CALL mppobc(vwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  482) 
(  483)          ! ... extremeties niw0, niw1 
(  484)          ij = jpjwd +1 - njmpp 
(  485)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  486)             DO jk = 1,jpkm1 
(  487)                vwbnd(ij,jk,nibm,nitm) = vwbnd(ij+1 ,jk,nibm,nitm)
(  488)             END DO 
(  489)          END IF
(  490)          ij = jpjwf +1 - njmpp 
(  491)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  492)             DO jk = 1,jpkm1 
(  493)                vwbnd(ij,jk,nibm,nitm) = vwbnd(ij-1 ,jk,nibm,nitm)
(  494)             END DO 
(  495)          END IF 
(  496)  
(  497)          ! 1.3 Temperature and salinity
(  498)          ! ----------------------------
(  499)  
(  500)          ! ... advance in time (time filter, array swap)
(  501)          DO jk = 1, jpkm1
(  502)             DO jj = 1, jpj
(  503)          ! ... fields nitm <== nit  plus time filter at the boundary
(  504)                twbnd(jj,jk,nib,nitm) = twbnd(jj,jk,nib,nit)*twmsk(jj,jk)
(  505)                swbnd(jj,jk,nib,nitm) = swbnd(jj,jk,nib,nit)*twmsk(jj,jk)
(  506)             END DO
(  507)          END DO
(  508)  
(  509)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  510)             DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 10

(  511)                DO jj = 1, jpj
(  512)                   twbnd(jj,jk,nibm ,nitm) = twbnd(jj,jk,nibm ,nit)*twmsk(jj,jk)
(  513)                   swbnd(jj,jk,nibm ,nitm) = swbnd(jj,jk,nibm ,nit)*twmsk(jj,jk)
(  514)          ! ... fields nit <== now (kt+1)
(  515)                   twbnd(jj,jk,nib  ,nit) = tn(ji   ,jj,jk)*twmsk(jj,jk)
(  516)                   twbnd(jj,jk,nibm ,nit) = tn(ji+1 ,jj,jk)*twmsk(jj,jk)
(  517)                   swbnd(jj,jk,nib  ,nit) = sn(ji   ,jj,jk)*twmsk(jj,jk)
(  518)                   swbnd(jj,jk,nibm ,nit) = sn(ji+1 ,jj,jk)*twmsk(jj,jk)
(  519)                END DO
(  520)             END DO
(  521)          END DO
(  522)          IF( lk_mpp )   CALL mppobc(twbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  523)          IF( lk_mpp )   CALL mppobc(swbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  524) 
(  525)          ! ... extremeties niw0, niw1
(  526)          ij = jpjwd +1 - njmpp
(  527)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  528)             DO jk = 1,jpkm1
(  529)                twbnd(ij,jk,nibm,nitm) = twbnd(ij+1 ,jk,nibm,nitm)
(  530)                swbnd(ij,jk,nibm,nitm) = swbnd(ij+1 ,jk,nibm,nitm)
(  531)             END DO
(  532)          END IF
(  533)          ij = jpjwf +1 - njmpp
(  534)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  535)             DO jk = 1,jpkm1
(  536)                twbnd(ij,jk,nibm,nitm) = twbnd(ij-1 ,jk,nibm,nitm)
(  537)                swbnd(ij,jk,nibm,nitm) = swbnd(ij-1 ,jk,nibm,nitm)
(  538)             END DO
(  539)          END IF
(  540)  
(  541)       END IF     ! End of array swap
(  542) 
(  543)       ! 2 - Calculation of radiation velocities
(  544)       ! ---------------------------------------
(  545)    
(  546)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  547)   
(  548)          ! 2.1  Calculate the normal velocity U based on phase velocity u_cxwbnd
(  549)          ! ----------------------------------------------------------------------
(  550)          !
(  551)          !          nib       nibm      nibm2
(  552)          !        ///|   nib   |   nibm  |
(  553)          !        ///|    |    |    |    |
(  554)          !        ---f----v----f----v----f-- jj-line
(  555)          !        ///|    |    |    |    |
(  556)          !        ///|         |         |
(  557)          !        ///u    T    u    T    u   jj-line
(  558)          !        ///|         |         |
(  559)          !        ///|    |    |    |    |
(  560)          !         jpiwob    jpiwob+1    jpiwob+2
(  561)          !                |         |        
(  562)          !              jpiwob+1    jpiwob+2     
(  563)          !
(  564)          ! ... If free surface formulation:
(  565)          ! ... radiative conditions on the total part + relaxation toward climatology
(  566)          ! ... (jpjwdp1, jpjwfm1), jpiwob
(  567)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  568)             DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 11

(  569)                DO jj = 2, jpjm1
(  570)          ! ... 2* gradi(u) (T-point i=nibm, time mean)
(  571)                   z2dx = ( - uwbnd(jj,jk,nibm ,nit) -  uwbnd(jj,jk,nibm ,nitm2) &
(  572)                            + 2.*uwbnd(jj,jk,nibm2,nitm) ) / e1t(ji+2,jj)
(  573)          ! ... 2* gradj(u) (u-point i=nibm, time nitm)
(  574)                   z2dy = ( uwbnd(jj+1,jk,nibm,nitm) - uwbnd(jj-1,jk,nibm,nitm) ) / e2u(ji+1,jj)
(  575)          ! ... square of the norm of grad(u)
(  576)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  577)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  578)                   zdt = uwbnd(jj,jk,nibm,nitm2) - uwbnd(jj,jk,nibm,nit)
(  579)          ! ... i-phase speed ratio (bounded by -1)
(  580)                   IF( z4nor2 == 0. ) THEN
(  581)                      z4nor2=0.00001
(  582)                   END IF
(  583)                   z05cx = zdt * z2dx / z4nor2
(  584)                   u_cxwbnd(jj,jk)=z05cx*uwmsk(jj,jk)
(  585)                END DO
(  586)             END DO
(  587)          END DO
(  588) 
(  589)          ! 2.2  Calculate the tangential velocity based on phase velocity v_cxwbnd
(  590)          ! -----------------------------------------------------------------------
(  591)          !
(  592)          !      nib       nibm     nibm2
(  593)          !    ///|///nib   |   nibm  |  nibm2
(  594)          !    ///|////|    |    |    |    |    |
(  595)          !    ---v----f----v----f----v----f----v-- jj-line
(  596)          !    ///|////|    |    |    |    |    |
(  597)          !    ///|////|    |    |    |    |    |
(  598)          !   jpiwob     jpiwob+1    jpiwob+2
(  599)          !            |         |         |   
(  600)          !          jpiwob   jpiwob+1   jpiwob+2    
(  601)          !
(  602)          ! ... radiative condition plus Raymond-Kuo
(  603)          ! ... (jpjwdp1, jpjwfm1),jpiwob
(  604)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  605)             DO jk = 1, jpkm1
(  606)                DO jj = 2, jpjm1
(  607)          ! ... 2* i-gradient of v (f-point i=nibm, time mean)
(  608)                   z2dx = ( - vwbnd(jj,jk,nibm ,nit) - vwbnd(jj,jk,nibm ,nitm2) &
(  609)                            + 2.*vwbnd(jj,jk,nibm2,nitm) ) / e1f(ji+1,jj)
(  610)          ! ... 2* j-gradient of v (v-point i=nibm, time nitm)
(  611)                   z2dy = ( vwbnd(jj+1,jk,nibm,nitm) - vwbnd(jj-1,jk,nibm,nitm) ) / e2v(ji+1,jj)
(  612)          ! ... square of the norm of grad(v)
(  613)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  614)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  615)                   zdt = vwbnd(jj,jk,nibm,nitm2) - vwbnd(jj,jk,nibm,nit)
(  616)          ! ... i-phase speed ratio (bounded by -1) and save the unbounded phase
(  617)          !     velocity ratio no divided by e1f for the tracer radiation
(  618)                   IF( z4nor2 == 0) THEN
(  619)                      z4nor2=0.000001
(  620)                   endif
(  621)                   z05cx = zdt * z2dx / z4nor2
(  622)                   v_cxwbnd(jj,jk) = z05cx*vwmsk(jj,jk)
(  623)                END DO
(  624)             END DO
(  625)          END DO
(  626)          IF( lk_mpp )   CALL mppobc(v_cxwbnd,jpjwd,jpjwf,jpiwob,jpk,2,jpj)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 12

(  627) 
(  628)          ! ... extremeties niw0, niw1
(  629)          ij = jpjwd +1 - njmpp
(  630)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  631)             DO jk = 1,jpkm1
(  632)                v_cxwbnd(ij,jk) = v_cxwbnd(ij+1 ,jk)
(  633)             END DO
(  634)          END IF
(  635)          ij = jpjwf +1 - njmpp
(  636)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  637)             DO jk = 1,jpkm1
(  638)                v_cxwbnd(ij,jk) = v_cxwbnd(ij-1 ,jk)
(  639)             END DO
(  640)          END IF
(  641) 
(  642)       END IF
(  643) 
(  644)    END SUBROUTINE obc_rad_west
(  645) 
(  646) 
(  647)    SUBROUTINE obc_rad_north ( kt )
(  648)       !!------------------------------------------------------------------------------
(  649)       !!                  ***  SUBROUTINE obc_rad_north  ***
(  650)       !!           
(  651)       !! ** Purpose :
(  652)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  653)       !!      north boundary and calculate those phase speeds if this OBC is not fixed.
(  654)       !!      In case of fixed OBC, this subrountine is not called.
(  655)       !!
(  656)       !!  History :
(  657)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  658)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  659)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  660)       !!         ! 00-06 (J.-M. Molines) 
(  661)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  662)       !!------------------------------------------------------------------------------
(  663)       !! * Arguments
(  664)       INTEGER, INTENT( in ) ::   kt
(  665) 
(  666)       !! * Local declarations
(  667)       INTEGER  ::   ii
(  668)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  669)       REAL(wp) ::   zvcb, zvcbm, zvcbm2
(  670)       !!------------------------------------------------------------------------------
(  671) 
(  672)       ! 1. Swap arrays before calculating radiative velocities
(  673)       ! ------------------------------------------------------
(  674) 
(  675)       ! 1.1  zonal velocity 
(  676)       ! -------------------
(  677) 
(  678)       IF( kt > nit000 ) THEN 
(  679) 
(  680)          ! ... advance in time (time filter, array swap)
(  681)          DO jk = 1, jpkm1
(  682)             DO ji = 1, jpi
(  683)          ! ... fields nitm2 <== nitm
(  684)                unbnd(ji,jk,nib  ,nitm2) = unbnd(ji,jk,nib  ,nitm)*unmsk(ji,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 13

(  685)                unbnd(ji,jk,nibm ,nitm2) = unbnd(ji,jk,nibm ,nitm)*unmsk(ji,jk)
(  686)                unbnd(ji,jk,nibm2,nitm2) = unbnd(ji,jk,nibm2,nitm)*unmsk(ji,jk)
(  687)             END DO
(  688)          END DO
(  689) 
(  690)          DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
(  691)             DO jk = 1, jpkm1
(  692)                DO ji = 1, jpi
(  693)                   unbnd(ji,jk,nib  ,nitm) = unbnd(ji,jk,nib,  nit)*unmsk(ji,jk)
(  694)                   unbnd(ji,jk,nibm ,nitm) = unbnd(ji,jk,nibm ,nit)*unmsk(ji,jk)
(  695)                   unbnd(ji,jk,nibm2,nitm) = unbnd(ji,jk,nibm2,nit)*unmsk(ji,jk)
(  696)          ! ... fields nit <== now (kt+1)
(  697)                   unbnd(ji,jk,nib  ,nit) = un(ji,jj,  jk)*unmsk(ji,jk)
(  698)                   unbnd(ji,jk,nibm ,nit) = un(ji,jj-1,jk)*unmsk(ji,jk)
(  699)                   unbnd(ji,jk,nibm2,nit) = un(ji,jj-2,jk)*unmsk(ji,jk)
(  700)                END DO
(  701)             END DO
(  702)          END DO
(  703)          IF( lk_mpp )   CALL mppobc(unbnd,jpind,jpinf,jpjnob+1,jpk*3*3,1,jpi)
(  704) 
(  705)          ! ... extremeties njn0,njn1 
(  706)          ii = jpind + 1 - nimpp 
(  707)          IF( ii >= 2 .AND. ii < jpim1 ) THEN 
(  708)             DO jk = 1, jpkm1
(  709)                 unbnd(ii,jk,nibm,nitm) = unbnd(ii+1,jk,nibm,nitm)
(  710)             END DO
(  711)          END IF 
(  712)          ii = jpinf + 1 - nimpp 
(  713)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  714)             DO jk = 1, jpkm1
(  715)                unbnd(ii,jk,nibm,nitm) = unbnd(ii-1,jk,nibm,nitm)
(  716)             END DO
(  717)          END IF
(  718)  
(  719)          ! 1.2. normal velocity 
(  720)          ! --------------------
(  721) 
(  722)          ! ... advance in time (time filter, array swap) 
(  723)          DO jk = 1, jpkm1
(  724)             DO ji = 1, jpi
(  725)          ! ... fields nitm2 <== nitm 
(  726)                vnbnd(ji,jk,nib  ,nitm2) = vnbnd(ji,jk,nib  ,nitm)*vnmsk(ji,jk)
(  727)                vnbnd(ji,jk,nibm ,nitm2) = vnbnd(ji,jk,nibm ,nitm)*vnmsk(ji,jk)
(  728)                vnbnd(ji,jk,nibm2,nitm2) = vnbnd(ji,jk,nibm2,nitm)*vnmsk(ji,jk)
(  729)             END DO
(  730)          END DO
(  731) 
(  732)          DO jj = fs_njn0, fs_njn1  ! Vector opt.
(  733)             DO jk = 1, jpkm1
(  734)                DO ji = 1, jpi
(  735)                   vnbnd(ji,jk,nib  ,nitm) = vnbnd(ji,jk,nib,  nit)*vnmsk(ji,jk)
(  736)                   vnbnd(ji,jk,nibm ,nitm) = vnbnd(ji,jk,nibm ,nit)*vnmsk(ji,jk)
(  737)                   vnbnd(ji,jk,nibm2,nitm) = vnbnd(ji,jk,nibm2,nit)*vnmsk(ji,jk)
(  738)          ! ... fields nit <== now (kt+1)
(  739)          ! ... total or baroclinic velocity at b, bm and bm2
(  740) # if defined key_dynspg_rl
(  741)                   zvcb   = vclin(ji,jk)
(  742) # else






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 14

(  743)                   zvcb   = vn (ji,jj,jk)
(  744) # endif
(  745) # if defined key_dynspg_rl
(  746)                   zvcbm  = vn (ji,jj-1,jk) - hvr (ji,jj-1) / e1v (ji,jj-1) &
(  747)                            * ( bsfn(ji,jj-1) - bsfn(ji-1,jj-1) )
(  748) # else
(  749)                   zvcbm  = vn (ji,jj-1,jk)
(  750) # endif
(  751) # if defined key_dynspg_rl
(  752)                   zvcbm2 = vn (ji,jj-2,jk) - hvr (ji,jj-2) / e1v (ji,jj-2) &
(  753)                            * ( bsfn(ji,jj-2) - bsfn(ji-1,jj-2) )
(  754) # else
(  755)                   zvcbm2 = vn (ji,jj-2,jk)
(  756) # endif
(  757)          ! ... fields nit <== now (kt+1) 
(  758)                   vnbnd(ji,jk,nib  ,nit) = zvcb  *vnmsk(ji,jk)
(  759)                   vnbnd(ji,jk,nibm ,nit) = zvcbm *vnmsk(ji,jk)
(  760)                   vnbnd(ji,jk,nibm2,nit) = zvcbm2*vnmsk(ji,jk)
(  761)                END DO
(  762)             END DO
(  763)          END DO
(  764)          IF( lk_mpp )   CALL mppobc(vnbnd,jpind,jpinf,jpjnob,jpk*3*3,1,jpi)
(  765) 
(  766)          ! ... extremeties njn0,njn1
(  767)          ii = jpind + 1 - nimpp
(  768)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  769)             DO jk = 1, jpkm1
(  770)                vnbnd(ii,jk,nibm,nitm) = vnbnd(ii+1,jk,nibm,nitm)
(  771)             END DO
(  772)          END IF
(  773)          ii = jpinf + 1 - nimpp
(  774)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  775)             DO jk = 1, jpkm1
(  776)                vnbnd(ii,jk,nibm,nitm) = vnbnd(ii-1,jk,nibm,nitm)
(  777)             END DO
(  778)          END IF
(  779) 
(  780)          ! 1.3 Temperature and salinity
(  781)          ! ----------------------------
(  782) 
(  783)          ! ... advance in time (time filter, array swap)
(  784)          DO jk = 1, jpkm1
(  785)             DO ji = 1, jpi
(  786)          ! ... fields nitm <== nit  plus time filter at the boundary
(  787)                tnbnd(ji,jk,nib ,nitm) = tnbnd(ji,jk,nib,nit)*tnmsk(ji,jk)
(  788)                snbnd(ji,jk,nib ,nitm) = snbnd(ji,jk,nib,nit)*tnmsk(ji,jk)
(  789)             END DO
(  790)          END DO
(  791) 
(  792)          DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
(  793)             DO jk = 1, jpkm1
(  794)                DO ji = 1, jpi
(  795)                   tnbnd(ji,jk,nibm ,nitm) = tnbnd(ji,jk,nibm ,nit)*tnmsk(ji,jk)
(  796)                   snbnd(ji,jk,nibm ,nitm) = snbnd(ji,jk,nibm ,nit)*tnmsk(ji,jk)
(  797)          ! ... fields nit <== now (kt+1)
(  798)                   tnbnd(ji,jk,nib  ,nit) = tn(ji,jj,  jk)*tnmsk(ji,jk)
(  799)                   tnbnd(ji,jk,nibm ,nit) = tn(ji,jj-1,jk)*tnmsk(ji,jk)
(  800)                   snbnd(ji,jk,nib  ,nit) = sn(ji,jj,  jk)*tnmsk(ji,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 15

(  801)                   snbnd(ji,jk,nibm ,nit) = sn(ji,jj-1,jk)*tnmsk(ji,jk)
(  802)                END DO
(  803)             END DO
(  804)          END DO
(  805)          IF( lk_mpp )   CALL mppobc(tnbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  806)          IF( lk_mpp )   CALL mppobc(snbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  807) 
(  808)          ! ... extremeties  njn0,njn1
(  809)          ii = jpind + 1 - nimpp
(  810)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  811)             DO jk = 1, jpkm1
(  812)                tnbnd(ii,jk,nibm,nitm) = tnbnd(ii+1,jk,nibm,nitm)
(  813)                snbnd(ii,jk,nibm,nitm) = snbnd(ii+1,jk,nibm,nitm)
(  814)             END DO
(  815)          END IF
(  816)          ii = jpinf + 1 - nimpp
(  817)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  818)             DO jk = 1, jpkm1
(  819)                tnbnd(ii,jk,nibm,nitm) = tnbnd(ii-1,jk,nibm,nitm)
(  820)                snbnd(ii,jk,nibm,nitm) = snbnd(ii-1,jk,nibm,nitm)
(  821)             END DO
(  822)          END IF
(  823) 
(  824)       END IF     ! End of array swap
(  825) 
(  826)       ! 2 - Calculation of radiation velocities
(  827)       ! ---------------------------------------
(  828) 
(  829)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  830) 
(  831)          ! 2.1  Calculate the normal velocity based on phase velocity u_cynbnd
(  832)          ! -------------------------------------------------------------------
(  833)          !
(  834)          !           ji-row
(  835)          !             |
(  836)          !     nib -///u//////  jpjnob + 1
(  837)          !        /////|//////
(  838)          !   nib  -----f-----   jpjnob
(  839)          !             |    
(  840)          !     nibm--  u   ---- jpjnob
(  841)          !             |        
(  842)          !  nibm  -----f-----   jpjnob-1
(  843)          !             |        
(  844)          !    nibm2--  u   ---- jpjnob-1
(  845)          !             |        
(  846)          !  nibm2 -----f-----   jpjnob-2
(  847)          !             |
(  848)          ! ... radiative condition
(  849)          ! ... jpjnob+1,(jpindp1, jpinfm1)
(  850)          DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
(  851)             DO jk = 1, jpkm1
(  852)                DO ji = 2, jpim1
(  853)          ! ... 2* j-gradient of u (f-point i=nibm, time mean)
(  854)                   z2dx = ( unbnd(ji,jk,nibm ,nit) + unbnd(ji,jk,nibm ,nitm2) &
(  855)                         - 2.*unbnd(ji,jk,nibm2,nitm)) / e2f(ji,jj-2)
(  856)          ! ... 2* i-gradient of u (u-point i=nibm, time nitm)
(  857)                   z2dy = ( unbnd(ji+1,jk,nibm,nitm) - unbnd(ji-1,jk,nibm,nitm) ) / e1u(ji,jj-1)
(  858)          ! ... square of the norm of grad(v)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 16

(  859)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  860)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  861)                   zdt = unbnd(ji,jk,nibm,nitm2) - unbnd(ji,jk,nibm,nit)
(  862)          ! ... i-phase speed ratio (bounded by 1) and save the unbounded phase
(  863)          !     velocity ratio no divided by e1f for the tracer radiation
(  864)                   IF( z4nor2 == 0.) THEN
(  865)                      z4nor2=.000001
(  866)                   END IF
(  867)                   z05cx = zdt * z2dx / z4nor2
(  868)                   u_cynbnd(ji,jk) = z05cx *unmsk(ji,jk)
(  869)                END DO
(  870)             END DO
(  871)          END DO
(  872)          IF( lk_mpp )   CALL mppobc(u_cynbnd,jpind,jpinf,jpjnob+1,jpk,1,jpi)
(  873) 
(  874)          ! ... extremeties  njn0,njn1
(  875)          ii = jpind + 1 - nimpp
(  876)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  877)             DO jk = 1, jpkm1
(  878)                u_cynbnd(ii,jk) = u_cynbnd(ii+1,jk)
(  879)             END DO
(  880)          END IF
(  881)          ii = jpinf + 1 - nimpp
(  882)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  883)             DO jk = 1, jpkm1
(  884)                u_cynbnd(ii,jk) = u_cynbnd(ii-1,jk)
(  885)             END DO
(  886)          END IF
(  887) 
(  888)          ! 2.2 Calculate the normal velocity based on phase velocity v_cynbnd 
(  889)          ! ------------------------------------------------------------------
(  890)          !
(  891)          !                ji-row  ji-row
(  892)          !                     |
(  893)          !        /////|/////////////////
(  894)          !   nib  -----f----v----f----  jpjnob
(  895)          !             |         |
(  896)          !     nib  -  u -- T -- u ---- jpjnob
(  897)          !             |         |
(  898)          !  nibm  -----f----v----f----  jpjnob-1
(  899)          !             |         |
(  900)          !    nibm --  u -- T -- u ---  jpjnob-1
(  901)          !             |         |
(  902)          !  nibm2 -----f----v----f----  jpjnob-2
(  903)          !             |         |
(  904)          ! ... If rigidlid formulation:
(  905)          ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  906)          ! ... If free surface formulation:
(  907)          ! ... radiative conditions on the total part + relaxation toward climatology 
(  908)          ! ... jpjnob,(jpindp1, jpinfm1)
(  909)          DO jj = fs_njn0, fs_njn1  ! Vector opt.
(  910)             DO jk = 1, jpkm1
(  911)                DO ji = 2, jpim1
(  912)          ! ... 2* gradj(v) (T-point i=nibm, time mean)
(  913)                   ii = ji -1 + nimpp
(  914)                   z2dx = ( vnbnd(ji,jk,nibm ,nit) + vnbnd(ji,jk,nibm ,nitm2) &
(  915)                           - 2.*vnbnd(ji,jk,nibm2,nitm)) / e2t(ji,jj-1)
(  916)          ! ... 2* gradi(v) (v-point i=nibm, time nitm)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 17

(  917)                   z2dy = ( vnbnd(ji+1,jk,nibm,nitm) - vnbnd(ji-1,jk,nibm,nitm) ) / e1v(ji,jj-1)
(  918)          ! ... square of the norm of grad(u)
(  919)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  920)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  921)                   zdt = vnbnd(ji,jk,nibm,nitm2) - vnbnd(ji,jk,nibm,nit)
(  922)          ! ... j-phase speed ratio (bounded by 1)
(  923)                   IF( z4nor2 == 0. ) THEN
(  924)                      z4nor2=.00001
(  925)                   END IF
(  926)                   z05cx = zdt * z2dx / z4nor2
(  927)                   v_cynbnd(ji,jk)=z05cx *vnmsk(ji,jk)
(  928)                END DO
(  929)             END DO
(  930)          END DO
(  931) 
(  932)       END IF
(  933) 
(  934)    END SUBROUTINE obc_rad_north
(  935) 
(  936) 
(  937)    SUBROUTINE obc_rad_south ( kt )
(  938)       !!------------------------------------------------------------------------------
(  939)       !!                  ***  SUBROUTINE obc_rad_south  ***
(  940)       !!           
(  941)       !! ** Purpose :
(  942)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  943)       !!      south boundary and calculate those phase speeds if this OBC is not fixed.
(  944)       !!      In case of fixed OBC, this subrountine is not called.
(  945)       !!
(  946)       !!  History :
(  947)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  948)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  949)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  950)       !!         ! 00-06 (J.-M. Molines) 
(  951)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  952)       !!------------------------------------------------------------------------------
(  953)       !! * Arguments
(  954)       INTEGER, INTENT( in ) ::   kt
(  955) 
(  956)       !! * Local declarations
(  957)       INTEGER ::   ii
(  958)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  959)       REAL(wp) ::   zvcb, zvcbm, zvcbm2
(  960)       !!------------------------------------------------------------------------------
(  961) 
(  962)       ! 1. Swap arrays before calculating radiative velocities
(  963)       ! ------------------------------------------------------
(  964) 
(  965)       ! 1.1  zonal velocity 
(  966)       ! --------------------
(  967)   
(  968)       IF( kt > nit000) THEN 
(  969) 
(  970)          ! ... advance in time (time filter, array swap)
(  971)          DO jk = 1, jpkm1
(  972)             DO ji = 1, jpi
(  973)          ! ... fields nitm2 <== nitm
(  974)                   usbnd(ji,jk,nib  ,nitm2) = usbnd(ji,jk,nib  ,nitm)*usmsk(ji,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 18

(  975)                   usbnd(ji,jk,nibm ,nitm2) = usbnd(ji,jk,nibm ,nitm)*usmsk(ji,jk)
(  976)                   usbnd(ji,jk,nibm2,nitm2) = usbnd(ji,jk,nibm2,nitm)*usmsk(ji,jk)
(  977)             END DO
(  978)          END DO
(  979)  
(  980)          DO jj = fs_njs0, fs_njs1  ! Vector opt.
(  981)             DO jk = 1, jpkm1
(  982)                DO ji = 1, jpi
(  983)                   usbnd(ji,jk,nib  ,nitm) = usbnd(ji,jk,nib,  nit)*usmsk(ji,jk)
(  984)                   usbnd(ji,jk,nibm ,nitm) = usbnd(ji,jk,nibm ,nit)*usmsk(ji,jk)
(  985)                   usbnd(ji,jk,nibm2,nitm) = usbnd(ji,jk,nibm2,nit)*usmsk(ji,jk)
(  986)          ! ... fields nit <== now (kt+1)
(  987)                   usbnd(ji,jk,nib  ,nit) = un(ji,jj  ,jk)*usmsk(ji,jk)
(  988)                   usbnd(ji,jk,nibm ,nit) = un(ji,jj+1,jk)*usmsk(ji,jk)
(  989)                   usbnd(ji,jk,nibm2,nit) = un(ji,jj+2,jk)*usmsk(ji,jk)
(  990)                END DO
(  991)             END DO
(  992)          END DO
(  993)          IF( lk_mpp )   CALL mppobc(usbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  994) 
(  995)          ! ... extremeties njs0,njs1
(  996)          ii = jpisd + 1 - nimpp
(  997)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  998)             DO jk = 1, jpkm1
(  999)                usbnd(ii,jk,nibm,nitm) = usbnd(ii+1,jk,nibm,nitm)
( 1000)             END DO
( 1001)          END IF
( 1002)          ii = jpisf + 1 - nimpp
( 1003)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1004)             DO jk = 1, jpkm1
( 1005)                usbnd(ii,jk,nibm,nitm) = usbnd(ii-1,jk,nibm,nitm)
( 1006)             END DO
( 1007)          END IF
( 1008)  
( 1009)          ! 1.2 normal velocity
( 1010)          ! -------------------
( 1011)  
( 1012)          !.. advance in time (time filter, array swap) 
( 1013)          DO jk = 1, jpkm1
( 1014)             DO ji = 1, jpi
( 1015)          ! ... fields nitm2 <== nitm 
( 1016)                vsbnd(ji,jk,nib  ,nitm2) = vsbnd(ji,jk,nib  ,nitm)*vsmsk(ji,jk)
( 1017)                vsbnd(ji,jk,nibm ,nitm2) = vsbnd(ji,jk,nibm ,nitm)*vsmsk(ji,jk)
( 1018)             END DO
( 1019)          END DO
( 1020) 
( 1021)          DO jj = fs_njs0, fs_njs1  ! Vector opt.
( 1022)             DO jk = 1, jpkm1
( 1023)                DO ji = 1, jpi
( 1024)                   vsbnd(ji,jk,nib  ,nitm) = vsbnd(ji,jk,nib,  nit)*vsmsk(ji,jk)
( 1025)                   vsbnd(ji,jk,nibm ,nitm) = vsbnd(ji,jk,nibm ,nit)*vsmsk(ji,jk)
( 1026)                   vsbnd(ji,jk,nibm2,nitm) = vsbnd(ji,jk,nibm2,nit)*vsmsk(ji,jk)
( 1027)          ! ... total or baroclinic velocity at b, bm and bm2
( 1028) # if defined key_dynspg_rl
( 1029)                   zvcb   = vclis(ji,jk)
( 1030) # else
( 1031)                   zvcb   = vn (ji,jj,jk)
( 1032) # endif






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 19

( 1033) # if defined key_dynspg_rl
( 1034)                   zvcbm  = vn (ji,jj+1,jk) - hvr (ji,jj+1) / e1v (ji,jj+1) & 
( 1035)                            * ( bsfn(ji,jj+1) - bsfn(ji-1,jj+1) )
( 1036) # else
( 1037)                   zvcbm  = vn (ji,jj+1,jk)
( 1038) # endif
( 1039) # if defined key_dynspg_rl  
( 1040)                   zvcbm2 = vn (ji,jj+2,jk) - hvr (ji,jj+2) / e1v (ji,jj+2) &
( 1041)                            * ( bsfn(ji,jj+2) - bsfn(ji-1,jj+2) )
( 1042) # else
( 1043)                   zvcbm2 = vn (ji,jj+2,jk)
( 1044) # endif
( 1045)          ! ... fields nit <== now (kt+1) 
( 1046)                   vsbnd(ji,jk,nib  ,nit) = zvcb   *vsmsk(ji,jk)
( 1047)                   vsbnd(ji,jk,nibm ,nit) = zvcbm  *vsmsk(ji,jk)
( 1048)                   vsbnd(ji,jk,nibm2,nit) = zvcbm2 *vsmsk(ji,jk)
( 1049)                END DO
( 1050)             END DO
( 1051)          END DO
( 1052)          IF( lk_mpp )   CALL mppobc(vsbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
( 1053) 
( 1054)          ! ... extremeties njs0,njs1
( 1055)          ii = jpisd + 1 - nimpp
( 1056)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1057)             DO jk = 1, jpkm1
( 1058)                vsbnd(ii,jk,nibm,nitm) = vsbnd(ii+1,jk,nibm,nitm)
( 1059)             END DO
( 1060)          END IF
( 1061)          ii = jpisf + 1 - nimpp
( 1062)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1063)             DO jk = 1, jpkm1
( 1064)                vsbnd(ii,jk,nibm,nitm) = vsbnd(ii-1,jk,nibm,nitm)
( 1065)             END DO
( 1066)          END IF
( 1067) 
( 1068)          ! 1.3 Temperature and salinity
( 1069)          ! ----------------------------
( 1070) 
( 1071)          ! ... advance in time (time filter, array swap)
( 1072)          DO jk = 1, jpkm1
( 1073)             DO ji = 1, jpi
( 1074)          ! ... fields nitm <== nit  plus time filter at the boundary
( 1075)                tsbnd(ji,jk,nib,nitm) = tsbnd(ji,jk,nib,nit)*tsmsk(ji,jk)
( 1076)                ssbnd(ji,jk,nib,nitm) = ssbnd(ji,jk,nib,nit)*tsmsk(ji,jk)
( 1077)             END DO
( 1078)          END DO
( 1079) 
( 1080)          DO jj = fs_njs0, fs_njs1  ! Vector opt.
( 1081)             DO jk = 1, jpkm1
( 1082)                DO ji = 1, jpi
( 1083)                   tsbnd(ji,jk,nibm ,nitm) = tsbnd(ji,jk,nibm ,nit)*tsmsk(ji,jk)
( 1084)                   ssbnd(ji,jk,nibm ,nitm) = ssbnd(ji,jk,nibm ,nit)*tsmsk(ji,jk)
( 1085)          ! ... fields nit <== now (kt+1)
( 1086)                   tsbnd(ji,jk,nib  ,nit) = tn(ji,jj   ,jk)*tsmsk(ji,jk)
( 1087)                   tsbnd(ji,jk,nibm ,nit) = tn(ji,jj+1 ,jk)*tsmsk(ji,jk)
( 1088)                   ssbnd(ji,jk,nib  ,nit) = sn(ji,jj   ,jk)*tsmsk(ji,jk)
( 1089)                   ssbnd(ji,jk,nibm ,nit) = sn(ji,jj+1 ,jk)*tsmsk(ji,jk)
( 1090)                END DO






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 20

( 1091)             END DO
( 1092)          END DO
( 1093)          IF( lk_mpp )   CALL mppobc(tsbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
( 1094)          IF( lk_mpp )   CALL mppobc(ssbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
( 1095) 
( 1096)          ! ... extremeties  njs0,njs1
( 1097)          ii = jpisd + 1 - nimpp
( 1098)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1099)             DO jk = 1, jpkm1
( 1100)                tsbnd(ii,jk,nibm,nitm) = tsbnd(ii+1,jk,nibm,nitm)
( 1101)                ssbnd(ii,jk,nibm,nitm) = ssbnd(ii+1,jk,nibm,nitm)
( 1102)             END DO
( 1103)          END IF
( 1104)          ii = jpisf + 1 - nimpp
( 1105)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1106)             DO jk = 1, jpkm1
( 1107)                tsbnd(ii,jk,nibm,nitm) = tsbnd(ii-1,jk,nibm,nitm)
( 1108)                ssbnd(ii,jk,nibm,nitm) = ssbnd(ii-1,jk,nibm,nitm)
( 1109)             END DO
( 1110)          END IF
( 1111) 
( 1112)       END IF     ! End of array swap
( 1113) 
( 1114)       ! 2 - Calculation of radiation velocities
( 1115)       ! ---------------------------------------
( 1116) 
( 1117)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
( 1118) 
( 1119)          ! 2.1  Calculate the normal velocity based on phase velocity u_cysbnd
( 1120)          ! -------------------------------------------------------------------
( 1121)          !
( 1122)          !          ji-row
( 1123)          !            |
( 1124)          ! nibm2 -----f-----   jpjsob +2
( 1125)          !            |    
( 1126)          !  nibm2 --  u  ----- jpjsob +2 
( 1127)          !            |        
( 1128)          !  nibm -----f-----   jpjsob +1
( 1129)          !            |        
( 1130)          !  nibm  --  u  ----- jpjsob +1
( 1131)          !            |        
( 1132)          !  nib  -----f-----   jpjsob
( 1133)          !       /////|//////
( 1134)          !  nib   ////u/////   jpjsob 
( 1135)          !
( 1136)          ! ... radiative condition plus Raymond-Kuo
( 1137)          ! ... jpjsob,(jpisdp1, jpisfm1)
( 1138)          DO jj = fs_njs0, fs_njs1  ! Vector opt.
( 1139)             DO jk = 1, jpkm1
( 1140)                DO ji = 2, jpim1
( 1141)          ! ... 2* j-gradient of u (f-point i=nibm, time mean)
( 1142)                   z2dx = (- usbnd(ji,jk,nibm ,nit) - usbnd(ji,jk,nibm ,nitm2) &
( 1143)                           + 2.*usbnd(ji,jk,nibm2,nitm) ) / e2f(ji,jj+1)
( 1144)          ! ... 2* i-gradient of u (u-point i=nibm, time nitm)
( 1145)                   z2dy = ( usbnd(ji+1,jk,nibm,nitm) - usbnd(ji-1,jk,nibm,nitm) ) / e1u(ji, jj+1)
( 1146)          ! ... square of the norm of grad(v)
( 1147)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
( 1148)                   IF( z4nor2 == 0.) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 21

( 1149)                      z4nor2 = 0.000001
( 1150)                   END IF
( 1151)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
( 1152)                   zdt = usbnd(ji,jk,nibm,nitm2) - usbnd(ji,jk,nibm,nit)
( 1153)          ! ... i-phase speed ratio (bounded by -1) and save the unbounded phase
( 1154)          !     velocity ratio no divided by e1f for the tracer radiation
( 1155)                   z05cx = zdt * z2dx / z4nor2
( 1156)                   u_cysbnd(ji,jk) = z05cx*usmsk(ji,jk)
( 1157)                END DO
( 1158)             END DO
( 1159)          END DO
( 1160)          IF( lk_mpp )   CALL mppobc(u_cysbnd,jpisd,jpisf,jpjsob,jpk,1,jpi)
( 1161) 
( 1162)          ! ... extremeties  njs0,njs1
( 1163)          ii = jpisd + 1 - nimpp
( 1164)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1165)             DO jk = 1, jpkm1
( 1166)                u_cysbnd(ii,jk) = u_cysbnd(ii+1,jk)
( 1167)             END DO
( 1168)          END IF
( 1169)          ii = jpisf + 1 - nimpp
( 1170)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1171)             DO jk = 1, jpkm1
( 1172)                u_cysbnd(ii,jk) = u_cysbnd(ii-1,jk)
( 1173)             END DO
( 1174)          END IF
( 1175) 
( 1176)          ! 2.2 Calculate the normal velocity based on phase velocity v_cysbnd 
( 1177)          ! -------------------------------------------------------------------
( 1178)          !
( 1179)          !               ji-row  ji-row
( 1180)          !            |         |
( 1181)          ! nibm2 -----f----v----f----  jpjsob+2
( 1182)          !            |         |
( 1183)          !  nibm   -  u -- T -- u ---- jpjsob+2
( 1184)          !            |         |
( 1185)          ! nibm  -----f----v----f----  jpjsob+1 
( 1186)          !            |         |
( 1187)          ! nib    --  u -- T -- u ---  jpjsob+1
( 1188)          !            |         |
( 1189)          ! nib   -----f----v----f----  jpjsob
( 1190)          !       /////////////////////
( 1191)          !
( 1192)          ! ... If rigidlid formulation:
( 1193)          ! ... radiative conditions on the baroclinic part only + relaxation toward climatology 
( 1194)          ! ... If free surface formulation:
( 1195)          ! ... radiative conditions on the total part + relaxation toward climatology
( 1196)          ! ... jpjsob,(jpisdp1,jpisfm1)
( 1197)          DO jj = fs_njs0, fs_njs1 ! Vector opt.
( 1198)             DO jk = 1, jpkm1
( 1199)                DO ji = 2, jpim1
( 1200)          ! ... 2* gradj(v) (T-point i=nibm, time mean)
( 1201)                   z2dx = ( - vsbnd(ji,jk,nibm ,nit) - vsbnd(ji,jk,nibm ,nitm2) &
( 1202)                            + 2.*vsbnd(ji,jk,nibm2,nitm) ) / e2t(ji,jj+1)
( 1203)          ! ... 2* gradi(v) (v-point i=nibm, time nitm)
( 1204)                   z2dy = ( vsbnd(ji+1,jk,nibm,nitm) - vsbnd(ji-1,jk,nibm,nitm) ) / e1v(ji,jj+1)
( 1205)          ! ... square of the norm of grad(u)
( 1206)                   z4nor2 = z2dx * z2dx + z2dy * z2dy






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 22

( 1207)                   IF( z4nor2 == 0.) THEN
( 1208)                      z4nor2 = 0.000001
( 1209)                   END IF
( 1210)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
( 1211)                   zdt = vsbnd(ji,jk,nibm,nitm2) - vsbnd(ji,jk,nibm,nit)
( 1212)          ! ... j-phase speed ratio (bounded by -1)
( 1213)                   z05cx = zdt * z2dx / z4nor2
( 1214)                   v_cysbnd(ji,jk)=z05cx*vsmsk(ji,jk)
( 1215)                END DO
( 1216)             END DO
( 1217)          END DO
( 1218) 
( 1219)       ENDIF
( 1220)  
( 1221)    END SUBROUTINE obc_rad_south
( 1222) 
( 1223) #else
( 1224)    !!=================================================================================
( 1225)    !!                       ***  MODULE  obcrad  ***
( 1226)    !! Ocean dynamic :   Phase velocities for each open boundary
( 1227)    !!=================================================================================
( 1228)    USE in_out_manager
( 1229) CONTAINS
( 1230)    SUBROUTINE obc_rad( kt )            ! No open boundaries ==> empty routine
( 1231)       INTEGER, INTENT(in) :: kt
( 1232)       if(lwp) WRITE(numout,*) 'obc_rad: You should not have seen this print! error?', kt
( 1233)    END SUBROUTINE obc_rad
( 1234) #endif
( 1235) 
( 1236) END MODULE obcrad


































PGF90 (Version     10.5)          09/23/2010  09:50:08      page 23

(    1) # 1 "obcrad.F90"
(    1) MODULE obcrad 
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcrad  ***
(    4)    !! Ocean dynamic :   Phase velocities for each open boundary
(    5)    !!=================================================================================
(    6) # 7
(    7)    !!---------------------------------------------------------------------------------
(    8)    !!   obc_rad        : call the subroutine for each open boundary
(    9)    !!   obc_rad_east   : compute the east phase velocities
(   10)    !!   obc_rad_west   : compute the west phase velocities
(   11)    !!   obc_rad_north  : compute the north phase velocities
(   12)    !!   obc_rad_south  : compute the south phase velocities
(   13)    !!---------------------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE phycst          ! physical constants
(   19)    USE obc_oce         ! ocean open boundary conditions
(   20)    USE lib_mpp         ! for mppobc
(   21)    USE in_out_manager  ! I/O units
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Accessibility
(   27)    PUBLIC obc_rad        ! routine called by step.F90
(   28) 
(   29)    !! * Module variables
(   30)    INTEGER ::   ji, jj, jk     ! dummy loop indices
(   31) 
(   32)    INTEGER ::      & ! ... boundary space indices 
(   33)       nib   = 1,   & ! nib   = boundary point
(   34)       nibm  = 2,   & ! nibm  = 1st interior point
(   35)       nibm2 = 3,   & ! nibm2 = 2nd interior point
(   36)                      ! ... boundary time indices 
(   37)       nit   = 1,   & ! nit    = now
(   38)       nitm  = 2,   & ! nitm   = before
(   39)       nitm2 = 3      ! nitm2  = before-before
(   40) 
(   41)    !! * Substitutions
(   42) # 1 "./obc_vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                 ***  obc_vectopt_loop_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_vectopt_loop_substitute.h90,v 1.3 2005/12/12 14:20:26 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 43 "obcrad.F90"
(   43) # 43
(   43)    !!---------------------------------------------------------------------------------
(   44)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 24

(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcrad.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   47)    !!---------------------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE obc_rad ( kt )
(   52)       !!------------------------------------------------------------------------------
(   53)       !!                     SUBROUTINE obc_rad
(   54)       !!                    ********************
(   55)       !! ** Purpose :
(   56)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(   57)       !!      boundaries and calculate those phase speeds if the open boundaries are 
(   58)       !!      not fixed. In case of fixed open boundaries does nothing.
(   59)       !!
(   60)       !!     The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   61)       !!     and/or lp_obc_south allow the user to determine which boundary is an
(   62)       !!     open one (must be done in the param_obc.h90 file).
(   63)       !! 
(   64)       !! ** Reference : 
(   65)       !!     Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   66)       !!
(   67)       !!  History :
(   68)       !!    8.5  !  02-10  (C. Talandier, A-M. Treguier) Free surface, F90 from the 
(   69)       !!                                                 J. Molines and G. Madec version
(   70)       !!------------------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt
(   73)       !!----------------------------------------------------------------------
(   74) 
(   75)       IF( lp_obc_east  .AND. .NOT.lfbceast  )   CALL obc_rad_east ( kt )   ! East open boundary
(   76) 
(   77)       IF( lp_obc_west  .AND. .NOT.lfbcwest  )   CALL obc_rad_west ( kt )   ! West open boundary
(   78) 
(   79)       IF( lp_obc_north .AND. .NOT.lfbcnorth )   CALL obc_rad_north( kt )   ! North open boundary
(   80) 
(   81)       IF( lp_obc_south .AND. .NOT.lfbcsouth )   CALL obc_rad_south( kt )   ! South open boundary
(   82) 
(   83)    END SUBROUTINE obc_rad
(   84) 
(   85) 
(   86)    SUBROUTINE obc_rad_east ( kt )
(   87)       !!------------------------------------------------------------------------------
(   88)       !!                     ***  SUBROUTINE obc_rad_east  ***
(   89)       !!                   
(   90)       !! ** Purpose :
(   91)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(   92)       !!      east boundary and calculate those phase speeds if this OBC is not fixed.
(   93)       !!      In case of fixed OBC, this subrountine is not called.
(   94)       !!
(   95)       !!  History :
(   96)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(   97)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(   98)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(   99)       !!         ! 00-06 (J.-M. Molines) 
(  100)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  101)       !!------------------------------------------------------------------------------
(  102)       !! * Arguments






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 25

(  103)       INTEGER, INTENT( in ) ::   kt
(  104) 
(  105)       !! * Local declarations
(  106)       INTEGER  ::   ij
(  107)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  108)       REAL(wp) ::   zucb, zucbm, zucbm2
(  109)       !!------------------------------------------------------------------------------
(  110) 
(  111)       ! 1. Swap arrays before calculating radiative velocities
(  112)       ! ------------------------------------------------------
(  113) 
(  114)       ! 1.1  zonal velocity 
(  115)       ! -------------------
(  116) 
(  117)       IF( kt > nit000 ) THEN 
(  118) 
(  119)          ! ... advance in time (time filter, array swap) 
(  120)          DO jk = 1, jpkm1
(  121)             DO jj = 1, jpj
(  122)                uebnd(jj,jk,nib  ,nitm2) = uebnd(jj,jk,nib  ,nitm)*uemsk(jj,jk)
(  123)                uebnd(jj,jk,nibm ,nitm2) = uebnd(jj,jk,nibm ,nitm)*uemsk(jj,jk)
(  124)                uebnd(jj,jk,nibm2,nitm2) = uebnd(jj,jk,nibm2,nitm)*uemsk(jj,jk)
(  125)             END DO
(  126)          END DO
(  127)          ! ... fields nitm <== nit  plus time filter at the boundary 
(  128)          DO ji = nie0 , nie1 ! Vector opt.
(  129)             DO jk = 1, jpkm1
(  130)                DO jj = 1, jpj
(  131)                   uebnd(jj,jk,nib  ,nitm) = uebnd(jj,jk,nib,  nit)*uemsk(jj,jk)
(  132)                   uebnd(jj,jk,nibm ,nitm) = uebnd(jj,jk,nibm ,nit)*uemsk(jj,jk)
(  133)                   uebnd(jj,jk,nibm2,nitm) = uebnd(jj,jk,nibm2,nit)*uemsk(jj,jk)
(  134)          ! ... fields nit <== now (kt+1) 
(  135)          ! ... Total or baroclinic velocity at b, bm and bm2
(  136) # 139
(  139)                   zucb   = un(ji,jj,jk)
(  140) # 145
(  145)                   zucbm  = un(ji-1,jj,jk)
(  146) # 151
(  151)                   zucbm2 = un(ji-2,jj,jk)
(  152) # 153
(  153)                   uebnd(jj,jk,nib  ,nit) = zucb   *uemsk(jj,jk)
(  154)                   uebnd(jj,jk,nibm ,nit) = zucbm  *uemsk(jj,jk) 
(  155)                   uebnd(jj,jk,nibm2,nit) = zucbm2 *uemsk(jj,jk) 
(  156)                END DO
(  157)             END DO
(  158)          END DO
(  159)          IF( lk_mpp )   CALL mppobc(uebnd,jpjed,jpjef,jpieob,jpk*3*3,2,jpj)
(  160) 
(  161)          ! ... extremeties nie0, nie1
(  162)          ij = jpjed +1 - njmpp
(  163)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  164)             DO jk = 1,jpkm1
(  165)                uebnd(ij,jk,nibm,nitm) = uebnd(ij+1 ,jk,nibm,nitm)
(  166)             END DO
(  167)          END IF
(  168)          ij = jpjef +1 - njmpp
(  169)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  170)             DO jk = 1,jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 26

(  171)                uebnd(ij,jk,nibm,nitm) = uebnd(ij-1 ,jk,nibm,nitm)
(  172)             END DO
(  173)          END IF
(  174) 
(  175)          ! 1.2 tangential velocity
(  176)          ! -----------------------
(  177) 
(  178)          ! ... advance in time (time filter, array swap)
(  179)          DO jk = 1, jpkm1
(  180)             DO jj = 1, jpj
(  181)          ! ... fields nitm2 <== nitm
(  182)                vebnd(jj,jk,nib  ,nitm2) = vebnd(jj,jk,nib  ,nitm)*vemsk(jj,jk)
(  183)                vebnd(jj,jk,nibm ,nitm2) = vebnd(jj,jk,nibm ,nitm)*vemsk(jj,jk)
(  184)                vebnd(jj,jk,nibm2,nitm2) = vebnd(jj,jk,nibm2,nitm)*vemsk(jj,jk)
(  185)             END DO
(  186)          END DO
(  187) 
(  188)          DO ji = nie0 +1, nie1+1 ! Vector opt.
(  189)             DO jk = 1, jpkm1
(  190)                DO jj = 1, jpj
(  191)                   vebnd(jj,jk,nib  ,nitm) = vebnd(jj,jk,nib,  nit)*vemsk(jj,jk)
(  192)                   vebnd(jj,jk,nibm ,nitm) = vebnd(jj,jk,nibm ,nit)*vemsk(jj,jk)
(  193)                   vebnd(jj,jk,nibm2,nitm) = vebnd(jj,jk,nibm2,nit)*vemsk(jj,jk)
(  194)          ! ... fields nit <== now (kt+1)
(  195)                   vebnd(jj,jk,nib  ,nit) = vn(ji  ,jj,jk)*vemsk(jj,jk)
(  196)                   vebnd(jj,jk,nibm ,nit) = vn(ji-1,jj,jk)*vemsk(jj,jk)
(  197)                   vebnd(jj,jk,nibm2,nit) = vn(ji-2,jj,jk)*vemsk(jj,jk)
(  198)                END DO
(  199)             END DO
(  200)          END DO
(  201)          IF( lk_mpp )   CALL mppobc(vebnd,jpjed,jpjef,jpieob+1,jpk*3*3,2,jpj)
(  202) 
(  203)          !... extremeties nie0, nie1
(  204)          ij = jpjed +1 - njmpp
(  205)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  206)             DO jk = 1,jpkm1
(  207)                vebnd(ij,jk,nibm,nitm) = vebnd(ij+1 ,jk,nibm,nitm)
(  208)             END DO 
(  209)          END IF 
(  210)          ij = jpjef +1 - njmpp 
(  211)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  212)             DO jk = 1,jpkm1 
(  213)                vebnd(ij,jk,nibm,nitm) = vebnd(ij-1 ,jk,nibm,nitm)
(  214)             END DO 
(  215)          END IF 
(  216) 
(  217)          ! 1.3 Temperature and salinity
(  218)          ! ----------------------------
(  219) 
(  220)          ! ... advance in time (time filter, array swap)
(  221)          DO jk = 1, jpkm1
(  222)             DO jj = 1, jpj
(  223)          ! ... fields nitm <== nit  plus time filter at the boundary
(  224)                tebnd(jj,jk,nib,nitm) = tebnd(jj,jk,nib,nit)*temsk(jj,jk)
(  225)                sebnd(jj,jk,nib,nitm) = sebnd(jj,jk,nib,nit)*temsk(jj,jk)
(  226)             END DO
(  227)          END DO
(  228) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 27

(  229)          DO ji = nie0 +1, nie1+1 ! Vector opt.
(  230)             DO jk = 1, jpkm1
(  231)                DO jj = 1, jpj
(  232)                   tebnd(jj,jk,nibm,nitm) = tebnd(jj,jk,nibm,nit)*temsk(jj,jk)
(  233)                   sebnd(jj,jk,nibm,nitm) = sebnd(jj,jk,nibm,nit)*temsk(jj,jk)
(  234)          ! ... fields nit <== now (kt+1)
(  235)                   tebnd(jj,jk,nib  ,nit) = tn(ji  ,jj,jk)*temsk(jj,jk)
(  236)                   tebnd(jj,jk,nibm ,nit) = tn(ji-1,jj,jk)*temsk(jj,jk)
(  237)                   sebnd(jj,jk,nib  ,nit) = sn(ji  ,jj,jk)*temsk(jj,jk)
(  238)                   sebnd(jj,jk,nibm ,nit) = sn(ji-1,jj,jk)*temsk(jj,jk)
(  239)                END DO
(  240)             END DO
(  241)          END DO
(  242)          IF( lk_mpp )   CALL mppobc(tebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  243)          IF( lk_mpp )   CALL mppobc(sebnd,jpjed,jpjef,jpieob+1,jpk*2*2,2,jpj)
(  244) 
(  245)          ! ... extremeties nie0, nie1
(  246)          ij = jpjed +1 - njmpp
(  247)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  248)             DO jk = 1,jpkm1
(  249)                tebnd(ij,jk,nibm,nitm) = tebnd(ij+1 ,jk,nibm,nitm)
(  250)                sebnd(ij,jk,nibm,nitm) = sebnd(ij+1 ,jk,nibm,nitm)
(  251)             END DO
(  252)          END IF
(  253)          ij = jpjef +1 - njmpp
(  254)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  255)             DO jk = 1,jpkm1
(  256)                tebnd(ij,jk,nibm,nitm) = tebnd(ij-1 ,jk,nibm,nitm)
(  257)                sebnd(ij,jk,nibm,nitm) = sebnd(ij-1 ,jk,nibm,nitm)
(  258)             END DO
(  259)          END IF
(  260) 
(  261)       END IF     ! End of array swap
(  262) 
(  263)       ! 2 - Calculation of radiation velocities
(  264)       ! ---------------------------------------
(  265) 
(  266)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  267) 
(  268)          ! 2.1  Calculate the normal velocity U based on phase velocity u_cxebnd
(  269)          ! ---------------------------------------------------------------------
(  270)          !
(  271)          !          nibm2      nibm      nib
(  272)          !            |  nibm   |   nib   |///
(  273)          !            |    |    |    |    |///
(  274)          !  jj-line --f----v----f----v----f---
(  275)          !            |    |    |    |    |///
(  276)          !            |         |         |///
(  277)          !  jj-line   u    T    u    T    u///
(  278)          !            |         |         |///
(  279)          !            |    |    |    |    |///
(  280)          !          jpieob-2   jpieob-1   jpieob
(  281)          !                 |         |        
(  282)          !              jpieob-1    jpieob      
(  283)          !
(  284)          ! ... (jpjedp1, jpjefm1),jpieob
(  285)          DO ji = nie0 , nie1 ! Vector opt.
(  286)             DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 28

(  287)                DO jj = 2, jpjm1
(  288)          ! ... 2* gradi(u) (T-point i=nibm, time mean)
(  289)                   z2dx = ( uebnd(jj,jk,nibm ,nit) + uebnd(jj,jk,nibm ,nitm2) &
(  290)                            - 2.*uebnd(jj,jk,nibm2,nitm) ) / e1t(ji-1,jj)
(  291)          ! ... 2* gradj(u) (u-point i=nibm, time nitm)
(  292)                   z2dy = ( uebnd(jj+1,jk,nibm,nitm) - uebnd(jj-1,jk,nibm,nitm) ) / e2u(ji-1,jj)
(  293)          ! ... square of the norm of grad(u)
(  294)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  295)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  296)                   zdt = uebnd(jj,jk,nibm,nitm2) - uebnd(jj,jk,nibm,nit)
(  297)          ! ... i-phase speed ratio (bounded by 1)               
(  298)                   IF( z4nor2 == 0. ) THEN
(  299)                      z4nor2=.00001
(  300)                   END IF
(  301)                   z05cx = zdt * z2dx / z4nor2
(  302)                   u_cxebnd(jj,jk) = z05cx*uemsk(jj,jk)
(  303)                END DO
(  304)             END DO
(  305)          END DO
(  306) 
(  307)          ! 2.2  Calculate the tangential velocity based on phase velocity v_cxebnd
(  308)          ! -----------------------------------------------------------------------
(  309)          !
(  310)          !          nibm2      nibm      nib
(  311)          !            |   nibm  |   nib///|///
(  312)          !            |    |    |    |////|///
(  313)          !  jj-line --v----f----v----f----v---
(  314)          !            |    |    |    |////|///
(  315)          !            |    |    |    |////|///
(  316)          !            | jpieob-1| jpieob /|///
(  317)          !            |         |         |   
(  318)          !         jpieob-1    jpieob     jpieob+1
(  319)          !
(  320)          ! ... (jpjedp1, jpjefm1), jpieob+1
(  321)          DO ji = nie0 +1, nie1+1 ! Vector opt.
(  322)             DO jk = 1, jpkm1
(  323)                DO jj = 2, jpjm1
(  324)          ! ... 2* i-gradient of v (f-point i=nibm, time mean)
(  325)                   z2dx = ( vebnd(jj,jk,nibm ,nit) + vebnd(jj,jk,nibm ,nitm2) &
(  326)                           - 2.*vebnd(jj,jk,nibm2,nitm) ) / e1f(ji-2,jj)
(  327)          ! ... 2* j-gradient of v (v-point i=nibm, time nitm)
(  328)                   z2dy = ( vebnd(jj+1,jk,nibm,nitm) -  vebnd(jj-1,jk,nibm,nitm) ) / e2v(ji-1,jj)
(  329)          ! ... square of the norm of grad(v)
(  330)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  331)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  332)                   zdt = vebnd(jj,jk,nibm,nitm2) - vebnd(jj,jk,nibm,nit)
(  333)          ! ... i-phase speed ratio (bounded by 1) and save the unbounded phase
(  334)          !     velocity ratio no divided by e1f for the tracer radiation
(  335)                   IF( z4nor2 == 0. ) THEN
(  336)                      z4nor2=.000001
(  337)                   END IF
(  338)                   z05cx = zdt * z2dx / z4nor2
(  339)                   v_cxebnd(jj,jk) = z05cx*vemsk(jj,jk)
(  340)                END DO
(  341)             END DO
(  342)          END DO
(  343)          IF( lk_mpp )   CALL mppobc(v_cxebnd,jpjed,jpjef,jpieob+1,jpk,2,jpj)
(  344) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 29

(  345)          ! ... extremeties nie0, nie1
(  346)          ij = jpjed +1 - njmpp
(  347)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  348)             DO jk = 1,jpkm1
(  349)                v_cxebnd(ij,jk) = v_cxebnd(ij+1 ,jk)
(  350)             END DO
(  351)          END IF
(  352)          ij = jpjef +1 - njmpp
(  353)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  354)             DO jk = 1,jpkm1
(  355)                v_cxebnd(ij,jk) = v_cxebnd(ij-1 ,jk)
(  356)             END DO
(  357)          END IF
(  358) 
(  359)       END IF
(  360) 
(  361)    END SUBROUTINE obc_rad_east
(  362) 
(  363) 
(  364)    SUBROUTINE obc_rad_west ( kt )
(  365)       !!------------------------------------------------------------------------------
(  366)       !!                  ***  SUBROUTINE obc_rad_west  ***
(  367)       !!                    
(  368)       !! ** Purpose :
(  369)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  370)       !!      west boundary and calculate those phase speeds if this OBC is not fixed.
(  371)       !!      In case of fixed OBC, this subrountine is not called.
(  372)       !!
(  373)       !!  History :
(  374)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  375)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  376)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  377)       !!         ! 00-06 (J.-M. Molines) 
(  378)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  379)       !!------------------------------------------------------------------------------
(  380)       !! * Arguments
(  381)       INTEGER, INTENT( in ) ::   kt
(  382) 
(  383)       !! * Local declarations
(  384)       INTEGER ::   ij
(  385)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  386)       REAL(wp) ::   zucb, zucbm, zucbm2
(  387)       !!------------------------------------------------------------------------------
(  388) 
(  389)       ! 1. Swap arrays before calculating radiative velocities
(  390)       ! ------------------------------------------------------
(  391) 
(  392)       ! 1.1  zonal velocity 
(  393)       ! -------------------
(  394) 
(  395)       IF( kt > nit000 ) THEN
(  396) 
(  397)          ! ... advance in time (time filter, array swap) 
(  398)          DO jk = 1, jpkm1
(  399)             DO jj = 1, jpj 
(  400)                uwbnd(jj,jk,nib  ,nitm2) = uwbnd(jj,jk,nib  ,nitm)*uwmsk(jj,jk)
(  401)                uwbnd(jj,jk,nibm ,nitm2) = uwbnd(jj,jk,nibm ,nitm)*uwmsk(jj,jk)
(  402)                uwbnd(jj,jk,nibm2,nitm2) = uwbnd(jj,jk,nibm2,nitm)*uwmsk(jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 30

(  403)             END DO
(  404)          END DO
(  405) 
(  406)          ! ... fields nitm <== nit  plus time filter at the boundary 
(  407)          DO ji = niw0 , niw1 ! Vector opt.
(  408)             DO jk = 1, jpkm1
(  409)                DO jj = 1, jpj
(  410)                   uwbnd(jj,jk,nib  ,nitm) = uwbnd(jj,jk,nib  ,nit)*uwmsk(jj,jk)
(  411)                   uwbnd(jj,jk,nibm ,nitm) = uwbnd(jj,jk,nibm ,nit)*uwmsk(jj,jk)
(  412)                   uwbnd(jj,jk,nibm2,nitm) = uwbnd(jj,jk,nibm2,nit)*uwmsk(jj,jk)
(  413)          ! ... total or baroclinic velocity at b, bm and bm2
(  414) # 417
(  417)                   zucb   = un (ji,jj,jk)
(  418) # 423
(  423)                   zucbm  = un (ji+1,jj,jk)
(  424) # 429
(  429)                   zucbm2 = un (ji+2,jj,jk)
(  430) 
(  431) # 432
(  432)          ! ... fields nit <== now (kt+1) 
(  433)                   uwbnd(jj,jk,nib  ,nit) = zucb  *uwmsk(jj,jk)
(  434)                   uwbnd(jj,jk,nibm ,nit) = zucbm *uwmsk(jj,jk)
(  435)                   uwbnd(jj,jk,nibm2,nit) = zucbm2*uwmsk(jj,jk)
(  436)                END DO
(  437)             END DO
(  438)          END DO
(  439)          IF( lk_mpp )   CALL mppobc(uwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  440) 
(  441)          ! ... extremeties niw0, niw1
(  442)          ij = jpjwd +1 - njmpp
(  443)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  444)             DO jk = 1,jpkm1
(  445)                uwbnd(ij,jk,nibm,nitm) = uwbnd(ij+1 ,jk,nibm,nitm)
(  446)             END DO
(  447)          END IF
(  448)          ij = jpjwf +1 - njmpp
(  449)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  450)             DO jk = 1,jpkm1
(  451)                uwbnd(ij,jk,nibm,nitm) = uwbnd(ij-1 ,jk,nibm,nitm)
(  452)             END DO
(  453)          END IF
(  454) 
(  455)          ! 1.2 tangential velocity
(  456)          ! -----------------------
(  457) 
(  458)          ! ... advance in time (time filter, array swap)
(  459)          DO jk = 1, jpkm1
(  460)             DO jj = 1, jpj 
(  461)          ! ... fields nitm2 <== nitm
(  462)                   vwbnd(jj,jk,nib  ,nitm2) = vwbnd(jj,jk,nib  ,nitm)*vwmsk(jj,jk)
(  463)                   vwbnd(jj,jk,nibm ,nitm2) = vwbnd(jj,jk,nibm ,nitm)*vwmsk(jj,jk)
(  464)                   vwbnd(jj,jk,nibm2,nitm2) = vwbnd(jj,jk,nibm2,nitm)*vwmsk(jj,jk)
(  465)             END DO
(  466)          END DO
(  467) 
(  468)          DO ji = niw0 , niw1 ! Vector opt.
(  469)             DO jk = 1, jpkm1
(  470)                DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 31

(  471)                   vwbnd(jj,jk,nib  ,nitm) = vwbnd(jj,jk,nib,  nit)*vwmsk(jj,jk)
(  472)                   vwbnd(jj,jk,nibm ,nitm) = vwbnd(jj,jk,nibm ,nit)*vwmsk(jj,jk)
(  473)                   vwbnd(jj,jk,nibm2,nitm) = vwbnd(jj,jk,nibm2,nit)*vwmsk(jj,jk)
(  474)          ! ... fields nit <== now (kt+1)
(  475)                   vwbnd(jj,jk,nib  ,nit) = vn(ji  ,jj,jk)*vwmsk(jj,jk)
(  476)                   vwbnd(jj,jk,nibm ,nit) = vn(ji+1,jj,jk)*vwmsk(jj,jk)
(  477)                   vwbnd(jj,jk,nibm2,nit) = vn(ji+2,jj,jk)*vwmsk(jj,jk)
(  478)                END DO
(  479)             END DO
(  480)          END DO
(  481)          IF( lk_mpp )   CALL mppobc(vwbnd,jpjwd,jpjwf,jpiwob,jpk*3*3,2,jpj)
(  482) 
(  483)          ! ... extremeties niw0, niw1 
(  484)          ij = jpjwd +1 - njmpp 
(  485)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  486)             DO jk = 1,jpkm1 
(  487)                vwbnd(ij,jk,nibm,nitm) = vwbnd(ij+1 ,jk,nibm,nitm)
(  488)             END DO 
(  489)          END IF
(  490)          ij = jpjwf +1 - njmpp 
(  491)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN 
(  492)             DO jk = 1,jpkm1 
(  493)                vwbnd(ij,jk,nibm,nitm) = vwbnd(ij-1 ,jk,nibm,nitm)
(  494)             END DO 
(  495)          END IF 
(  496)  
(  497)          ! 1.3 Temperature and salinity
(  498)          ! ----------------------------
(  499)  
(  500)          ! ... advance in time (time filter, array swap)
(  501)          DO jk = 1, jpkm1
(  502)             DO jj = 1, jpj
(  503)          ! ... fields nitm <== nit  plus time filter at the boundary
(  504)                twbnd(jj,jk,nib,nitm) = twbnd(jj,jk,nib,nit)*twmsk(jj,jk)
(  505)                swbnd(jj,jk,nib,nitm) = swbnd(jj,jk,nib,nit)*twmsk(jj,jk)
(  506)             END DO
(  507)          END DO
(  508)  
(  509)          DO ji = niw0 , niw1 ! Vector opt.
(  510)             DO jk = 1, jpkm1
(  511)                DO jj = 1, jpj
(  512)                   twbnd(jj,jk,nibm ,nitm) = twbnd(jj,jk,nibm ,nit)*twmsk(jj,jk)
(  513)                   swbnd(jj,jk,nibm ,nitm) = swbnd(jj,jk,nibm ,nit)*twmsk(jj,jk)
(  514)          ! ... fields nit <== now (kt+1)
(  515)                   twbnd(jj,jk,nib  ,nit) = tn(ji   ,jj,jk)*twmsk(jj,jk)
(  516)                   twbnd(jj,jk,nibm ,nit) = tn(ji+1 ,jj,jk)*twmsk(jj,jk)
(  517)                   swbnd(jj,jk,nib  ,nit) = sn(ji   ,jj,jk)*twmsk(jj,jk)
(  518)                   swbnd(jj,jk,nibm ,nit) = sn(ji+1 ,jj,jk)*twmsk(jj,jk)
(  519)                END DO
(  520)             END DO
(  521)          END DO
(  522)          IF( lk_mpp )   CALL mppobc(twbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  523)          IF( lk_mpp )   CALL mppobc(swbnd,jpjwd,jpjwf,jpiwob,jpk*2*2,2,jpj)
(  524) 
(  525)          ! ... extremeties niw0, niw1
(  526)          ij = jpjwd +1 - njmpp
(  527)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  528)             DO jk = 1,jpkm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 32

(  529)                twbnd(ij,jk,nibm,nitm) = twbnd(ij+1 ,jk,nibm,nitm)
(  530)                swbnd(ij,jk,nibm,nitm) = swbnd(ij+1 ,jk,nibm,nitm)
(  531)             END DO
(  532)          END IF
(  533)          ij = jpjwf +1 - njmpp
(  534)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  535)             DO jk = 1,jpkm1
(  536)                twbnd(ij,jk,nibm,nitm) = twbnd(ij-1 ,jk,nibm,nitm)
(  537)                swbnd(ij,jk,nibm,nitm) = swbnd(ij-1 ,jk,nibm,nitm)
(  538)             END DO
(  539)          END IF
(  540)  
(  541)       END IF     ! End of array swap
(  542) 
(  543)       ! 2 - Calculation of radiation velocities
(  544)       ! ---------------------------------------
(  545)    
(  546)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  547)   
(  548)          ! 2.1  Calculate the normal velocity U based on phase velocity u_cxwbnd
(  549)          ! ----------------------------------------------------------------------
(  550)          !
(  551)          !          nib       nibm      nibm2
(  552)          !        ///|   nib   |   nibm  |
(  553)          !        ///|    |    |    |    |
(  554)          !        ---f----v----f----v----f-- jj-line
(  555)          !        ///|    |    |    |    |
(  556)          !        ///|         |         |
(  557)          !        ///u    T    u    T    u   jj-line
(  558)          !        ///|         |         |
(  559)          !        ///|    |    |    |    |
(  560)          !         jpiwob    jpiwob+1    jpiwob+2
(  561)          !                |         |        
(  562)          !              jpiwob+1    jpiwob+2     
(  563)          !
(  564)          ! ... If free surface formulation:
(  565)          ! ... radiative conditions on the total part + relaxation toward climatology
(  566)          ! ... (jpjwdp1, jpjwfm1), jpiwob
(  567)          DO ji = niw0 , niw1 ! Vector opt.
(  568)             DO jk = 1, jpkm1
(  569)                DO jj = 2, jpjm1
(  570)          ! ... 2* gradi(u) (T-point i=nibm, time mean)
(  571)                   z2dx = ( - uwbnd(jj,jk,nibm ,nit) -  uwbnd(jj,jk,nibm ,nitm2) &
(  572)                            + 2.*uwbnd(jj,jk,nibm2,nitm) ) / e1t(ji+2,jj)
(  573)          ! ... 2* gradj(u) (u-point i=nibm, time nitm)
(  574)                   z2dy = ( uwbnd(jj+1,jk,nibm,nitm) - uwbnd(jj-1,jk,nibm,nitm) ) / e2u(ji+1,jj)
(  575)          ! ... square of the norm of grad(u)
(  576)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  577)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  578)                   zdt = uwbnd(jj,jk,nibm,nitm2) - uwbnd(jj,jk,nibm,nit)
(  579)          ! ... i-phase speed ratio (bounded by -1)
(  580)                   IF( z4nor2 == 0. ) THEN
(  581)                      z4nor2=0.00001
(  582)                   END IF
(  583)                   z05cx = zdt * z2dx / z4nor2
(  584)                   u_cxwbnd(jj,jk)=z05cx*uwmsk(jj,jk)
(  585)                END DO
(  586)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 33

(  587)          END DO
(  588) 
(  589)          ! 2.2  Calculate the tangential velocity based on phase velocity v_cxwbnd
(  590)          ! -----------------------------------------------------------------------
(  591)          !
(  592)          !      nib       nibm     nibm2
(  593)          !    ///|///nib   |   nibm  |  nibm2
(  594)          !    ///|////|    |    |    |    |    |
(  595)          !    ---v----f----v----f----v----f----v-- jj-line
(  596)          !    ///|////|    |    |    |    |    |
(  597)          !    ///|////|    |    |    |    |    |
(  598)          !   jpiwob     jpiwob+1    jpiwob+2
(  599)          !            |         |         |   
(  600)          !          jpiwob   jpiwob+1   jpiwob+2    
(  601)          !
(  602)          ! ... radiative condition plus Raymond-Kuo
(  603)          ! ... (jpjwdp1, jpjwfm1),jpiwob
(  604)          DO ji = niw0 , niw1 ! Vector opt.
(  605)             DO jk = 1, jpkm1
(  606)                DO jj = 2, jpjm1
(  607)          ! ... 2* i-gradient of v (f-point i=nibm, time mean)
(  608)                   z2dx = ( - vwbnd(jj,jk,nibm ,nit) - vwbnd(jj,jk,nibm ,nitm2) &
(  609)                            + 2.*vwbnd(jj,jk,nibm2,nitm) ) / e1f(ji+1,jj)
(  610)          ! ... 2* j-gradient of v (v-point i=nibm, time nitm)
(  611)                   z2dy = ( vwbnd(jj+1,jk,nibm,nitm) - vwbnd(jj-1,jk,nibm,nitm) ) / e2v(ji+1,jj)
(  612)          ! ... square of the norm of grad(v)
(  613)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  614)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  615)                   zdt = vwbnd(jj,jk,nibm,nitm2) - vwbnd(jj,jk,nibm,nit)
(  616)          ! ... i-phase speed ratio (bounded by -1) and save the unbounded phase
(  617)          !     velocity ratio no divided by e1f for the tracer radiation
(  618)                   IF( z4nor2 == 0) THEN
(  619)                      z4nor2=0.000001
(  620)                   endif
(  621)                   z05cx = zdt * z2dx / z4nor2
(  622)                   v_cxwbnd(jj,jk) = z05cx*vwmsk(jj,jk)
(  623)                END DO
(  624)             END DO
(  625)          END DO
(  626)          IF( lk_mpp )   CALL mppobc(v_cxwbnd,jpjwd,jpjwf,jpiwob,jpk,2,jpj)
(  627) 
(  628)          ! ... extremeties niw0, niw1
(  629)          ij = jpjwd +1 - njmpp
(  630)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  631)             DO jk = 1,jpkm1
(  632)                v_cxwbnd(ij,jk) = v_cxwbnd(ij+1 ,jk)
(  633)             END DO
(  634)          END IF
(  635)          ij = jpjwf +1 - njmpp
(  636)          IF( ij >= 2 .AND. ij < jpjm1 ) THEN
(  637)             DO jk = 1,jpkm1
(  638)                v_cxwbnd(ij,jk) = v_cxwbnd(ij-1 ,jk)
(  639)             END DO
(  640)          END IF
(  641) 
(  642)       END IF
(  643) 
(  644)    END SUBROUTINE obc_rad_west






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 34

(  645) 
(  646) 
(  647)    SUBROUTINE obc_rad_north ( kt )
(  648)       !!------------------------------------------------------------------------------
(  649)       !!                  ***  SUBROUTINE obc_rad_north  ***
(  650)       !!           
(  651)       !! ** Purpose :
(  652)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  653)       !!      north boundary and calculate those phase speeds if this OBC is not fixed.
(  654)       !!      In case of fixed OBC, this subrountine is not called.
(  655)       !!
(  656)       !!  History :
(  657)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  658)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  659)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  660)       !!         ! 00-06 (J.-M. Molines) 
(  661)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  662)       !!------------------------------------------------------------------------------
(  663)       !! * Arguments
(  664)       INTEGER, INTENT( in ) ::   kt
(  665) 
(  666)       !! * Local declarations
(  667)       INTEGER  ::   ii
(  668)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  669)       REAL(wp) ::   zvcb, zvcbm, zvcbm2
(  670)       !!------------------------------------------------------------------------------
(  671) 
(  672)       ! 1. Swap arrays before calculating radiative velocities
(  673)       ! ------------------------------------------------------
(  674) 
(  675)       ! 1.1  zonal velocity 
(  676)       ! -------------------
(  677) 
(  678)       IF( kt > nit000 ) THEN 
(  679) 
(  680)          ! ... advance in time (time filter, array swap)
(  681)          DO jk = 1, jpkm1
(  682)             DO ji = 1, jpi
(  683)          ! ... fields nitm2 <== nitm
(  684)                unbnd(ji,jk,nib  ,nitm2) = unbnd(ji,jk,nib  ,nitm)*unmsk(ji,jk)
(  685)                unbnd(ji,jk,nibm ,nitm2) = unbnd(ji,jk,nibm ,nitm)*unmsk(ji,jk)
(  686)                unbnd(ji,jk,nibm2,nitm2) = unbnd(ji,jk,nibm2,nitm)*unmsk(ji,jk)
(  687)             END DO
(  688)          END DO
(  689) 
(  690)          DO jj = njn0+1, njn1+1  ! Vector opt.
(  691)             DO jk = 1, jpkm1
(  692)                DO ji = 1, jpi
(  693)                   unbnd(ji,jk,nib  ,nitm) = unbnd(ji,jk,nib,  nit)*unmsk(ji,jk)
(  694)                   unbnd(ji,jk,nibm ,nitm) = unbnd(ji,jk,nibm ,nit)*unmsk(ji,jk)
(  695)                   unbnd(ji,jk,nibm2,nitm) = unbnd(ji,jk,nibm2,nit)*unmsk(ji,jk)
(  696)          ! ... fields nit <== now (kt+1)
(  697)                   unbnd(ji,jk,nib  ,nit) = un(ji,jj,  jk)*unmsk(ji,jk)
(  698)                   unbnd(ji,jk,nibm ,nit) = un(ji,jj-1,jk)*unmsk(ji,jk)
(  699)                   unbnd(ji,jk,nibm2,nit) = un(ji,jj-2,jk)*unmsk(ji,jk)
(  700)                END DO
(  701)             END DO
(  702)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 35

(  703)          IF( lk_mpp )   CALL mppobc(unbnd,jpind,jpinf,jpjnob+1,jpk*3*3,1,jpi)
(  704) 
(  705)          ! ... extremeties njn0,njn1 
(  706)          ii = jpind + 1 - nimpp 
(  707)          IF( ii >= 2 .AND. ii < jpim1 ) THEN 
(  708)             DO jk = 1, jpkm1
(  709)                 unbnd(ii,jk,nibm,nitm) = unbnd(ii+1,jk,nibm,nitm)
(  710)             END DO
(  711)          END IF 
(  712)          ii = jpinf + 1 - nimpp 
(  713)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  714)             DO jk = 1, jpkm1
(  715)                unbnd(ii,jk,nibm,nitm) = unbnd(ii-1,jk,nibm,nitm)
(  716)             END DO
(  717)          END IF
(  718)  
(  719)          ! 1.2. normal velocity 
(  720)          ! --------------------
(  721) 
(  722)          ! ... advance in time (time filter, array swap) 
(  723)          DO jk = 1, jpkm1
(  724)             DO ji = 1, jpi
(  725)          ! ... fields nitm2 <== nitm 
(  726)                vnbnd(ji,jk,nib  ,nitm2) = vnbnd(ji,jk,nib  ,nitm)*vnmsk(ji,jk)
(  727)                vnbnd(ji,jk,nibm ,nitm2) = vnbnd(ji,jk,nibm ,nitm)*vnmsk(ji,jk)
(  728)                vnbnd(ji,jk,nibm2,nitm2) = vnbnd(ji,jk,nibm2,nitm)*vnmsk(ji,jk)
(  729)             END DO
(  730)          END DO
(  731) 
(  732)          DO jj = njn0, njn1  ! Vector opt.
(  733)             DO jk = 1, jpkm1
(  734)                DO ji = 1, jpi
(  735)                   vnbnd(ji,jk,nib  ,nitm) = vnbnd(ji,jk,nib,  nit)*vnmsk(ji,jk)
(  736)                   vnbnd(ji,jk,nibm ,nitm) = vnbnd(ji,jk,nibm ,nit)*vnmsk(ji,jk)
(  737)                   vnbnd(ji,jk,nibm2,nitm) = vnbnd(ji,jk,nibm2,nit)*vnmsk(ji,jk)
(  738)          ! ... fields nit <== now (kt+1)
(  739)          ! ... total or baroclinic velocity at b, bm and bm2
(  740) # 743
(  743)                   zvcb   = vn (ji,jj,jk)
(  744) # 749
(  749)                   zvcbm  = vn (ji,jj-1,jk)
(  750) # 755
(  755)                   zvcbm2 = vn (ji,jj-2,jk)
(  756) # 757
(  757)          ! ... fields nit <== now (kt+1) 
(  758)                   vnbnd(ji,jk,nib  ,nit) = zvcb  *vnmsk(ji,jk)
(  759)                   vnbnd(ji,jk,nibm ,nit) = zvcbm *vnmsk(ji,jk)
(  760)                   vnbnd(ji,jk,nibm2,nit) = zvcbm2*vnmsk(ji,jk)
(  761)                END DO
(  762)             END DO
(  763)          END DO
(  764)          IF( lk_mpp )   CALL mppobc(vnbnd,jpind,jpinf,jpjnob,jpk*3*3,1,jpi)
(  765) 
(  766)          ! ... extremeties njn0,njn1
(  767)          ii = jpind + 1 - nimpp
(  768)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  769)             DO jk = 1, jpkm1
(  770)                vnbnd(ii,jk,nibm,nitm) = vnbnd(ii+1,jk,nibm,nitm)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 36

(  771)             END DO
(  772)          END IF
(  773)          ii = jpinf + 1 - nimpp
(  774)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  775)             DO jk = 1, jpkm1
(  776)                vnbnd(ii,jk,nibm,nitm) = vnbnd(ii-1,jk,nibm,nitm)
(  777)             END DO
(  778)          END IF
(  779) 
(  780)          ! 1.3 Temperature and salinity
(  781)          ! ----------------------------
(  782) 
(  783)          ! ... advance in time (time filter, array swap)
(  784)          DO jk = 1, jpkm1
(  785)             DO ji = 1, jpi
(  786)          ! ... fields nitm <== nit  plus time filter at the boundary
(  787)                tnbnd(ji,jk,nib ,nitm) = tnbnd(ji,jk,nib,nit)*tnmsk(ji,jk)
(  788)                snbnd(ji,jk,nib ,nitm) = snbnd(ji,jk,nib,nit)*tnmsk(ji,jk)
(  789)             END DO
(  790)          END DO
(  791) 
(  792)          DO jj = njn0+1, njn1+1  ! Vector opt.
(  793)             DO jk = 1, jpkm1
(  794)                DO ji = 1, jpi
(  795)                   tnbnd(ji,jk,nibm ,nitm) = tnbnd(ji,jk,nibm ,nit)*tnmsk(ji,jk)
(  796)                   snbnd(ji,jk,nibm ,nitm) = snbnd(ji,jk,nibm ,nit)*tnmsk(ji,jk)
(  797)          ! ... fields nit <== now (kt+1)
(  798)                   tnbnd(ji,jk,nib  ,nit) = tn(ji,jj,  jk)*tnmsk(ji,jk)
(  799)                   tnbnd(ji,jk,nibm ,nit) = tn(ji,jj-1,jk)*tnmsk(ji,jk)
(  800)                   snbnd(ji,jk,nib  ,nit) = sn(ji,jj,  jk)*tnmsk(ji,jk)
(  801)                   snbnd(ji,jk,nibm ,nit) = sn(ji,jj-1,jk)*tnmsk(ji,jk)
(  802)                END DO
(  803)             END DO
(  804)          END DO
(  805)          IF( lk_mpp )   CALL mppobc(tnbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  806)          IF( lk_mpp )   CALL mppobc(snbnd,jpind,jpinf,jpjnob+1,jpk*2*2,1,jpi)
(  807) 
(  808)          ! ... extremeties  njn0,njn1
(  809)          ii = jpind + 1 - nimpp
(  810)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  811)             DO jk = 1, jpkm1
(  812)                tnbnd(ii,jk,nibm,nitm) = tnbnd(ii+1,jk,nibm,nitm)
(  813)                snbnd(ii,jk,nibm,nitm) = snbnd(ii+1,jk,nibm,nitm)
(  814)             END DO
(  815)          END IF
(  816)          ii = jpinf + 1 - nimpp
(  817)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  818)             DO jk = 1, jpkm1
(  819)                tnbnd(ii,jk,nibm,nitm) = tnbnd(ii-1,jk,nibm,nitm)
(  820)                snbnd(ii,jk,nibm,nitm) = snbnd(ii-1,jk,nibm,nitm)
(  821)             END DO
(  822)          END IF
(  823) 
(  824)       END IF     ! End of array swap
(  825) 
(  826)       ! 2 - Calculation of radiation velocities
(  827)       ! ---------------------------------------
(  828) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 37

(  829)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
(  830) 
(  831)          ! 2.1  Calculate the normal velocity based on phase velocity u_cynbnd
(  832)          ! -------------------------------------------------------------------
(  833)          !
(  834)          !           ji-row
(  835)          !             |
(  836)          !     nib -///u//////  jpjnob + 1
(  837)          !        /////|//////
(  838)          !   nib  -----f-----   jpjnob
(  839)          !             |    
(  840)          !     nibm--  u   ---- jpjnob
(  841)          !             |        
(  842)          !  nibm  -----f-----   jpjnob-1
(  843)          !             |        
(  844)          !    nibm2--  u   ---- jpjnob-1
(  845)          !             |        
(  846)          !  nibm2 -----f-----   jpjnob-2
(  847)          !             |
(  848)          ! ... radiative condition
(  849)          ! ... jpjnob+1,(jpindp1, jpinfm1)
(  850)          DO jj = njn0+1, njn1+1  ! Vector opt.
(  851)             DO jk = 1, jpkm1
(  852)                DO ji = 2, jpim1
(  853)          ! ... 2* j-gradient of u (f-point i=nibm, time mean)
(  854)                   z2dx = ( unbnd(ji,jk,nibm ,nit) + unbnd(ji,jk,nibm ,nitm2) &
(  855)                         - 2.*unbnd(ji,jk,nibm2,nitm)) / e2f(ji,jj-2)
(  856)          ! ... 2* i-gradient of u (u-point i=nibm, time nitm)
(  857)                   z2dy = ( unbnd(ji+1,jk,nibm,nitm) - unbnd(ji-1,jk,nibm,nitm) ) / e1u(ji,jj-1)
(  858)          ! ... square of the norm of grad(v)
(  859)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  860)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  861)                   zdt = unbnd(ji,jk,nibm,nitm2) - unbnd(ji,jk,nibm,nit)
(  862)          ! ... i-phase speed ratio (bounded by 1) and save the unbounded phase
(  863)          !     velocity ratio no divided by e1f for the tracer radiation
(  864)                   IF( z4nor2 == 0.) THEN
(  865)                      z4nor2=.000001
(  866)                   END IF
(  867)                   z05cx = zdt * z2dx / z4nor2
(  868)                   u_cynbnd(ji,jk) = z05cx *unmsk(ji,jk)
(  869)                END DO
(  870)             END DO
(  871)          END DO
(  872)          IF( lk_mpp )   CALL mppobc(u_cynbnd,jpind,jpinf,jpjnob+1,jpk,1,jpi)
(  873) 
(  874)          ! ... extremeties  njn0,njn1
(  875)          ii = jpind + 1 - nimpp
(  876)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  877)             DO jk = 1, jpkm1
(  878)                u_cynbnd(ii,jk) = u_cynbnd(ii+1,jk)
(  879)             END DO
(  880)          END IF
(  881)          ii = jpinf + 1 - nimpp
(  882)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  883)             DO jk = 1, jpkm1
(  884)                u_cynbnd(ii,jk) = u_cynbnd(ii-1,jk)
(  885)             END DO
(  886)          END IF






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 38

(  887) 
(  888)          ! 2.2 Calculate the normal velocity based on phase velocity v_cynbnd 
(  889)          ! ------------------------------------------------------------------
(  890)          !
(  891)          !                ji-row  ji-row
(  892)          !                     |
(  893)          !        /////|/////////////////
(  894)          !   nib  -----f----v----f----  jpjnob
(  895)          !             |         |
(  896)          !     nib  -  u -- T -- u ---- jpjnob
(  897)          !             |         |
(  898)          !  nibm  -----f----v----f----  jpjnob-1
(  899)          !             |         |
(  900)          !    nibm --  u -- T -- u ---  jpjnob-1
(  901)          !             |         |
(  902)          !  nibm2 -----f----v----f----  jpjnob-2
(  903)          !             |         |
(  904)          ! ... If rigidlid formulation:
(  905)          ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  906)          ! ... If free surface formulation:
(  907)          ! ... radiative conditions on the total part + relaxation toward climatology 
(  908)          ! ... jpjnob,(jpindp1, jpinfm1)
(  909)          DO jj = njn0, njn1  ! Vector opt.
(  910)             DO jk = 1, jpkm1
(  911)                DO ji = 2, jpim1
(  912)          ! ... 2* gradj(v) (T-point i=nibm, time mean)
(  913)                   ii = ji -1 + nimpp
(  914)                   z2dx = ( vnbnd(ji,jk,nibm ,nit) + vnbnd(ji,jk,nibm ,nitm2) &
(  915)                           - 2.*vnbnd(ji,jk,nibm2,nitm)) / e2t(ji,jj-1)
(  916)          ! ... 2* gradi(v) (v-point i=nibm, time nitm)
(  917)                   z2dy = ( vnbnd(ji+1,jk,nibm,nitm) - vnbnd(ji-1,jk,nibm,nitm) ) / e1v(ji,jj-1)
(  918)          ! ... square of the norm of grad(u)
(  919)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
(  920)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  921)                   zdt = vnbnd(ji,jk,nibm,nitm2) - vnbnd(ji,jk,nibm,nit)
(  922)          ! ... j-phase speed ratio (bounded by 1)
(  923)                   IF( z4nor2 == 0. ) THEN
(  924)                      z4nor2=.00001
(  925)                   END IF
(  926)                   z05cx = zdt * z2dx / z4nor2
(  927)                   v_cynbnd(ji,jk)=z05cx *vnmsk(ji,jk)
(  928)                END DO
(  929)             END DO
(  930)          END DO
(  931) 
(  932)       END IF
(  933) 
(  934)    END SUBROUTINE obc_rad_north
(  935) 
(  936) 
(  937)    SUBROUTINE obc_rad_south ( kt )
(  938)       !!------------------------------------------------------------------------------
(  939)       !!                  ***  SUBROUTINE obc_rad_south  ***
(  940)       !!           
(  941)       !! ** Purpose :
(  942)       !!      Perform swap of arrays to calculate radiative phase speeds at the open 
(  943)       !!      south boundary and calculate those phase speeds if this OBC is not fixed.
(  944)       !!      In case of fixed OBC, this subrountine is not called.






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 39

(  945)       !!
(  946)       !!  History :
(  947)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  948)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  949)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  950)       !!         ! 00-06 (J.-M. Molines) 
(  951)       !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  952)       !!------------------------------------------------------------------------------
(  953)       !! * Arguments
(  954)       INTEGER, INTENT( in ) ::   kt
(  955) 
(  956)       !! * Local declarations
(  957)       INTEGER ::   ii
(  958)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  959)       REAL(wp) ::   zvcb, zvcbm, zvcbm2
(  960)       !!------------------------------------------------------------------------------
(  961) 
(  962)       ! 1. Swap arrays before calculating radiative velocities
(  963)       ! ------------------------------------------------------
(  964) 
(  965)       ! 1.1  zonal velocity 
(  966)       ! --------------------
(  967)   
(  968)       IF( kt > nit000) THEN 
(  969) 
(  970)          ! ... advance in time (time filter, array swap)
(  971)          DO jk = 1, jpkm1
(  972)             DO ji = 1, jpi
(  973)          ! ... fields nitm2 <== nitm
(  974)                   usbnd(ji,jk,nib  ,nitm2) = usbnd(ji,jk,nib  ,nitm)*usmsk(ji,jk)
(  975)                   usbnd(ji,jk,nibm ,nitm2) = usbnd(ji,jk,nibm ,nitm)*usmsk(ji,jk)
(  976)                   usbnd(ji,jk,nibm2,nitm2) = usbnd(ji,jk,nibm2,nitm)*usmsk(ji,jk)
(  977)             END DO
(  978)          END DO
(  979)  
(  980)          DO jj = njs0, njs1  ! Vector opt.
(  981)             DO jk = 1, jpkm1
(  982)                DO ji = 1, jpi
(  983)                   usbnd(ji,jk,nib  ,nitm) = usbnd(ji,jk,nib,  nit)*usmsk(ji,jk)
(  984)                   usbnd(ji,jk,nibm ,nitm) = usbnd(ji,jk,nibm ,nit)*usmsk(ji,jk)
(  985)                   usbnd(ji,jk,nibm2,nitm) = usbnd(ji,jk,nibm2,nit)*usmsk(ji,jk)
(  986)          ! ... fields nit <== now (kt+1)
(  987)                   usbnd(ji,jk,nib  ,nit) = un(ji,jj  ,jk)*usmsk(ji,jk)
(  988)                   usbnd(ji,jk,nibm ,nit) = un(ji,jj+1,jk)*usmsk(ji,jk)
(  989)                   usbnd(ji,jk,nibm2,nit) = un(ji,jj+2,jk)*usmsk(ji,jk)
(  990)                END DO
(  991)             END DO
(  992)          END DO
(  993)          IF( lk_mpp )   CALL mppobc(usbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
(  994) 
(  995)          ! ... extremeties njs0,njs1
(  996)          ii = jpisd + 1 - nimpp
(  997)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
(  998)             DO jk = 1, jpkm1
(  999)                usbnd(ii,jk,nibm,nitm) = usbnd(ii+1,jk,nibm,nitm)
( 1000)             END DO
( 1001)          END IF
( 1002)          ii = jpisf + 1 - nimpp






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 40

( 1003)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1004)             DO jk = 1, jpkm1
( 1005)                usbnd(ii,jk,nibm,nitm) = usbnd(ii-1,jk,nibm,nitm)
( 1006)             END DO
( 1007)          END IF
( 1008)  
( 1009)          ! 1.2 normal velocity
( 1010)          ! -------------------
( 1011)  
( 1012)          !.. advance in time (time filter, array swap) 
( 1013)          DO jk = 1, jpkm1
( 1014)             DO ji = 1, jpi
( 1015)          ! ... fields nitm2 <== nitm 
( 1016)                vsbnd(ji,jk,nib  ,nitm2) = vsbnd(ji,jk,nib  ,nitm)*vsmsk(ji,jk)
( 1017)                vsbnd(ji,jk,nibm ,nitm2) = vsbnd(ji,jk,nibm ,nitm)*vsmsk(ji,jk)
( 1018)             END DO
( 1019)          END DO
( 1020) 
( 1021)          DO jj = njs0, njs1  ! Vector opt.
( 1022)             DO jk = 1, jpkm1
( 1023)                DO ji = 1, jpi
( 1024)                   vsbnd(ji,jk,nib  ,nitm) = vsbnd(ji,jk,nib,  nit)*vsmsk(ji,jk)
( 1025)                   vsbnd(ji,jk,nibm ,nitm) = vsbnd(ji,jk,nibm ,nit)*vsmsk(ji,jk)
( 1026)                   vsbnd(ji,jk,nibm2,nitm) = vsbnd(ji,jk,nibm2,nit)*vsmsk(ji,jk)
( 1027)          ! ... total or baroclinic velocity at b, bm and bm2
( 1028) # 1031
( 1031)                   zvcb   = vn (ji,jj,jk)
( 1032) # 1037
( 1037)                   zvcbm  = vn (ji,jj+1,jk)
( 1038) # 1043
( 1043)                   zvcbm2 = vn (ji,jj+2,jk)
( 1044) # 1045
( 1045)          ! ... fields nit <== now (kt+1) 
( 1046)                   vsbnd(ji,jk,nib  ,nit) = zvcb   *vsmsk(ji,jk)
( 1047)                   vsbnd(ji,jk,nibm ,nit) = zvcbm  *vsmsk(ji,jk)
( 1048)                   vsbnd(ji,jk,nibm2,nit) = zvcbm2 *vsmsk(ji,jk)
( 1049)                END DO
( 1050)             END DO
( 1051)          END DO
( 1052)          IF( lk_mpp )   CALL mppobc(vsbnd,jpisd,jpisf,jpjsob,jpk*3*3,1,jpi)
( 1053) 
( 1054)          ! ... extremeties njs0,njs1
( 1055)          ii = jpisd + 1 - nimpp
( 1056)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1057)             DO jk = 1, jpkm1
( 1058)                vsbnd(ii,jk,nibm,nitm) = vsbnd(ii+1,jk,nibm,nitm)
( 1059)             END DO
( 1060)          END IF
( 1061)          ii = jpisf + 1 - nimpp
( 1062)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1063)             DO jk = 1, jpkm1
( 1064)                vsbnd(ii,jk,nibm,nitm) = vsbnd(ii-1,jk,nibm,nitm)
( 1065)             END DO
( 1066)          END IF
( 1067) 
( 1068)          ! 1.3 Temperature and salinity
( 1069)          ! ----------------------------
( 1070) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 41

( 1071)          ! ... advance in time (time filter, array swap)
( 1072)          DO jk = 1, jpkm1
( 1073)             DO ji = 1, jpi
( 1074)          ! ... fields nitm <== nit  plus time filter at the boundary
( 1075)                tsbnd(ji,jk,nib,nitm) = tsbnd(ji,jk,nib,nit)*tsmsk(ji,jk)
( 1076)                ssbnd(ji,jk,nib,nitm) = ssbnd(ji,jk,nib,nit)*tsmsk(ji,jk)
( 1077)             END DO
( 1078)          END DO
( 1079) 
( 1080)          DO jj = njs0, njs1  ! Vector opt.
( 1081)             DO jk = 1, jpkm1
( 1082)                DO ji = 1, jpi
( 1083)                   tsbnd(ji,jk,nibm ,nitm) = tsbnd(ji,jk,nibm ,nit)*tsmsk(ji,jk)
( 1084)                   ssbnd(ji,jk,nibm ,nitm) = ssbnd(ji,jk,nibm ,nit)*tsmsk(ji,jk)
( 1085)          ! ... fields nit <== now (kt+1)
( 1086)                   tsbnd(ji,jk,nib  ,nit) = tn(ji,jj   ,jk)*tsmsk(ji,jk)
( 1087)                   tsbnd(ji,jk,nibm ,nit) = tn(ji,jj+1 ,jk)*tsmsk(ji,jk)
( 1088)                   ssbnd(ji,jk,nib  ,nit) = sn(ji,jj   ,jk)*tsmsk(ji,jk)
( 1089)                   ssbnd(ji,jk,nibm ,nit) = sn(ji,jj+1 ,jk)*tsmsk(ji,jk)
( 1090)                END DO
( 1091)             END DO
( 1092)          END DO
( 1093)          IF( lk_mpp )   CALL mppobc(tsbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
( 1094)          IF( lk_mpp )   CALL mppobc(ssbnd,jpisd,jpisf,jpjsob,jpk*2*2,1,jpi)
( 1095) 
( 1096)          ! ... extremeties  njs0,njs1
( 1097)          ii = jpisd + 1 - nimpp
( 1098)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1099)             DO jk = 1, jpkm1
( 1100)                tsbnd(ii,jk,nibm,nitm) = tsbnd(ii+1,jk,nibm,nitm)
( 1101)                ssbnd(ii,jk,nibm,nitm) = ssbnd(ii+1,jk,nibm,nitm)
( 1102)             END DO
( 1103)          END IF
( 1104)          ii = jpisf + 1 - nimpp
( 1105)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1106)             DO jk = 1, jpkm1
( 1107)                tsbnd(ii,jk,nibm,nitm) = tsbnd(ii-1,jk,nibm,nitm)
( 1108)                ssbnd(ii,jk,nibm,nitm) = ssbnd(ii-1,jk,nibm,nitm)
( 1109)             END DO
( 1110)          END IF
( 1111) 
( 1112)       END IF     ! End of array swap
( 1113) 
( 1114)       ! 2 - Calculation of radiation velocities
( 1115)       ! ---------------------------------------
( 1116) 
( 1117)       IF( kt >= nit000 +3 .OR. ln_rstart ) THEN
( 1118) 
( 1119)          ! 2.1  Calculate the normal velocity based on phase velocity u_cysbnd
( 1120)          ! -------------------------------------------------------------------
( 1121)          !
( 1122)          !          ji-row
( 1123)          !            |
( 1124)          ! nibm2 -----f-----   jpjsob +2
( 1125)          !            |    
( 1126)          !  nibm2 --  u  ----- jpjsob +2 
( 1127)          !            |        
( 1128)          !  nibm -----f-----   jpjsob +1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 42

( 1129)          !            |        
( 1130)          !  nibm  --  u  ----- jpjsob +1
( 1131)          !            |        
( 1132)          !  nib  -----f-----   jpjsob
( 1133)          !       /////|//////
( 1134)          !  nib   ////u/////   jpjsob 
( 1135)          !
( 1136)          ! ... radiative condition plus Raymond-Kuo
( 1137)          ! ... jpjsob,(jpisdp1, jpisfm1)
( 1138)          DO jj = njs0, njs1  ! Vector opt.
( 1139)             DO jk = 1, jpkm1
( 1140)                DO ji = 2, jpim1
( 1141)          ! ... 2* j-gradient of u (f-point i=nibm, time mean)
( 1142)                   z2dx = (- usbnd(ji,jk,nibm ,nit) - usbnd(ji,jk,nibm ,nitm2) &
( 1143)                           + 2.*usbnd(ji,jk,nibm2,nitm) ) / e2f(ji,jj+1)
( 1144)          ! ... 2* i-gradient of u (u-point i=nibm, time nitm)
( 1145)                   z2dy = ( usbnd(ji+1,jk,nibm,nitm) - usbnd(ji-1,jk,nibm,nitm) ) / e1u(ji, jj+1)
( 1146)          ! ... square of the norm of grad(v)
( 1147)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
( 1148)                   IF( z4nor2 == 0.) THEN
( 1149)                      z4nor2 = 0.000001
( 1150)                   END IF
( 1151)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
( 1152)                   zdt = usbnd(ji,jk,nibm,nitm2) - usbnd(ji,jk,nibm,nit)
( 1153)          ! ... i-phase speed ratio (bounded by -1) and save the unbounded phase
( 1154)          !     velocity ratio no divided by e1f for the tracer radiation
( 1155)                   z05cx = zdt * z2dx / z4nor2
( 1156)                   u_cysbnd(ji,jk) = z05cx*usmsk(ji,jk)
( 1157)                END DO
( 1158)             END DO
( 1159)          END DO
( 1160)          IF( lk_mpp )   CALL mppobc(u_cysbnd,jpisd,jpisf,jpjsob,jpk,1,jpi)
( 1161) 
( 1162)          ! ... extremeties  njs0,njs1
( 1163)          ii = jpisd + 1 - nimpp
( 1164)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1165)             DO jk = 1, jpkm1
( 1166)                u_cysbnd(ii,jk) = u_cysbnd(ii+1,jk)
( 1167)             END DO
( 1168)          END IF
( 1169)          ii = jpisf + 1 - nimpp
( 1170)          IF( ii >= 2 .AND. ii < jpim1 ) THEN
( 1171)             DO jk = 1, jpkm1
( 1172)                u_cysbnd(ii,jk) = u_cysbnd(ii-1,jk)
( 1173)             END DO
( 1174)          END IF
( 1175) 
( 1176)          ! 2.2 Calculate the normal velocity based on phase velocity v_cysbnd 
( 1177)          ! -------------------------------------------------------------------
( 1178)          !
( 1179)          !               ji-row  ji-row
( 1180)          !            |         |
( 1181)          ! nibm2 -----f----v----f----  jpjsob+2
( 1182)          !            |         |
( 1183)          !  nibm   -  u -- T -- u ---- jpjsob+2
( 1184)          !            |         |
( 1185)          ! nibm  -----f----v----f----  jpjsob+1 
( 1186)          !            |         |






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 43

( 1187)          ! nib    --  u -- T -- u ---  jpjsob+1
( 1188)          !            |         |
( 1189)          ! nib   -----f----v----f----  jpjsob
( 1190)          !       /////////////////////
( 1191)          !
( 1192)          ! ... If rigidlid formulation:
( 1193)          ! ... radiative conditions on the baroclinic part only + relaxation toward climatology 
( 1194)          ! ... If free surface formulation:
( 1195)          ! ... radiative conditions on the total part + relaxation toward climatology
( 1196)          ! ... jpjsob,(jpisdp1,jpisfm1)
( 1197)          DO jj = njs0, njs1 ! Vector opt.
( 1198)             DO jk = 1, jpkm1
( 1199)                DO ji = 2, jpim1
( 1200)          ! ... 2* gradj(v) (T-point i=nibm, time mean)
( 1201)                   z2dx = ( - vsbnd(ji,jk,nibm ,nit) - vsbnd(ji,jk,nibm ,nitm2) &
( 1202)                            + 2.*vsbnd(ji,jk,nibm2,nitm) ) / e2t(ji,jj+1)
( 1203)          ! ... 2* gradi(v) (v-point i=nibm, time nitm)
( 1204)                   z2dy = ( vsbnd(ji+1,jk,nibm,nitm) - vsbnd(ji-1,jk,nibm,nitm) ) / e1v(ji,jj+1)
( 1205)          ! ... square of the norm of grad(u)
( 1206)                   z4nor2 = z2dx * z2dx + z2dy * z2dy
( 1207)                   IF( z4nor2 == 0.) THEN
( 1208)                      z4nor2 = 0.000001
( 1209)                   END IF
( 1210)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
( 1211)                   zdt = vsbnd(ji,jk,nibm,nitm2) - vsbnd(ji,jk,nibm,nit)
( 1212)          ! ... j-phase speed ratio (bounded by -1)
( 1213)                   z05cx = zdt * z2dx / z4nor2
( 1214)                   v_cysbnd(ji,jk)=z05cx*vsmsk(ji,jk)
( 1215)                END DO
( 1216)             END DO
( 1217)          END DO
( 1218) 
( 1219)       ENDIF
( 1220)  
( 1221)    END SUBROUTINE obc_rad_south
( 1222) 
( 1223) 
( 1224) # 1236
( 1236) END MODULE obcrad
