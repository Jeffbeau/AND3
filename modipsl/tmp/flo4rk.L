


PGF90 (Version     10.5)          09/23/2010  09:49:47      page 1

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

Filename: flo4rk.F90

(    1) MODULE flo4rk
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  flo4rk  ***
(    4)    !! Ocean floats :   trajectory computation using a 4th order Runge-Kutta
(    5)    !!======================================================================
(    6) #if   defined key_floats   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_floats'                                     float trajectories
(    9)    !!----------------------------------------------------------------------
(   10)    !!   flo_4rk        : Compute the geographical position of floats
(   11)    !!   flo_interp     : interpolation
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE flo_oce         ! ocean drifting floats
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE in_out_manager  ! I/O manager
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC flo_4rk      ! routine called by floats.F90
(   24) 
(   25)    !! * Module variables
(   26)    REAL(wp), DIMENSION (4) ::   &   ! RK4 and Lagrange interpolation
(   27)       tcoef1 = (/  1.0  ,  0.5  ,  0.5  ,  0.0  /) ,  &  ! coeffients  for 
(   28)       tcoef2 = (/  0.0  ,  0.5  ,  0.5  ,  1.0  /) ,  &  ! lagrangian interp.
(   29)       scoef2 = (/  1.0  ,  2.0  ,  2.0  ,  1.0  /) ,  &  ! RK4 coefficients
(   30)       rcoef  = (/-1./6. , 1./2. ,-1./2. , 1./6. /)       ! ???
(   31)    REAL(wp), DIMENSION (3) ::   &
(   32)       scoef1 = (/ .5, .5, 1. /)       ! compute position with interpolated
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flo4rk.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE flo_4rk( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE flo_4rk  ***
(   44)       !!
(   45)       !!  ** Purpose :   Compute the geographical position (lat,lon,depth)
(   46)       !!       of each float at each time step.






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 2

(   47)       !! 
(   48)       !!  ** Method  :   The position of a float is computed with a 4th order
(   49)       !!       Runge-Kutta scheme and and Lagrange interpolation.
(   50)       !!         We need to know the velocity field, the old positions of the
(   51)       !!       floats and the grid defined on the domain.
(   52)       !!
(   53)       !!----------------------------------------------------------------------
(   54)       !! * Arguments
(   55)       INTEGER, INTENT(in) ::   kt     ! ocean time-step index
(   56) 
(   57)       !! * Local declarations
(   58)       INTEGER ::  jfl, jind           ! dummy loop indices
(   59)       REAL(wp), DIMENSION ( jpnfl)  ::   &
(   60)          zgifl, zgjfl, zgkfl,   &     ! index RK  positions
(   61)          zufl, zvfl, zwfl             ! interpolated velocity at the
(   62)       ! float position 
(   63)       REAL(wp), DIMENSION ( jpnfl, 4 )  ::    &
(   64)          zrkxfl, zrkyfl, zrkzfl       ! RK coefficients
(   65)       !!---------------------------------------------------------------------
(   66)     
(   67)       IF( kt == nit000 ) THEN
(   68)          IF(lwp) WRITE(numout,*)
(   69)          IF(lwp) WRITE(numout,*) 'flo_4rk : compute Runge Kutta trajectories for floats '
(   70)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   71)       ENDIF
(   72) 
(   73)       ! Verification of the floats positions. If one of them leave the domain
(   74)       ! domain we replace the float near the border.
(   75)       DO jfl = 1, jpnfl
(   76)          ! i-direction
(   77)          IF( tpifl(jfl) <= 1.5 ) THEN
(   78)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(   79)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the WEST border.'
(   80)             tpifl(jfl) = tpifl(jfl) + 1.
(   81)             IF(lwp)WRITE(numout,*)'New initialisation for this float at i=',tpifl(jfl)
(   82)          ENDIF
(   83)           
(   84)          IF( tpifl(jfl) >= jpi-.5 ) THEN
(   85)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(   86)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the EAST border.'
(   87)             tpifl(jfl) = tpifl(jfl) - 1.
(   88)             IF(lwp)WRITE(numout,*)'New initialisation for this float at i=', tpifl(jfl)
(   89)          ENDIF
(   90)          ! j-direction
(   91)          IF( tpjfl(jfl) <= 1.5 ) THEN
(   92)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(   93)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the SOUTH border.'
(   94)             tpjfl(jfl) = tpjfl(jfl) + 1.
(   95)             IF(lwp)WRITE(numout,*)'New initialisation for this float at j=', tpjfl(jfl)
(   96)          ENDIF
(   97)            
(   98)          IF( tpjfl(jfl) >= jpj-.5 ) THEN
(   99)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(  100)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the NORTH border.'
(  101)             tpjfl(jfl) = tpjfl(jfl) - 1.
(  102)             IF(lwp)WRITE(numout,*)'New initialisation for this float at j=', tpjfl(jfl)
(  103)          ENDIF
(  104)          ! k-direction






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 3

(  105)          IF( tpkfl(jfl) <= .5 ) THEN
(  106)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(  107)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the TOP border.'
(  108)             tpkfl(jfl) = tpkfl(jfl) + 1.
(  109)             IF(lwp)WRITE(numout,*)'New initialisation for this float at k=', tpkfl(jfl)
(  110)          ENDIF
(  111)          
(  112)          IF( tpkfl(jfl) >= jpk-.5 )  THEN
(  113)             IF(lwp)WRITE(numout,*)'!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!'
(  114)             IF(lwp)WRITE(numout,*)'The float',jfl,'is out of the domain at the BOTTOM border.'
(  115)             tpkfl(jfl) = tpkfl(jfl) - 1.
(  116)             IF(lwp)WRITE(numout,*)'New initialisation for this float at k=', tpkfl(jfl)
(  117)          ENDIF
(  118)       END DO
(  119)       
(  120)       ! 4 steps of Runge-Kutta algorithme
(  121)       ! initialisation of the positions 
(  122)       
(  123)       DO jfl = 1, jpnfl
(  124)          zgifl(jfl) = tpifl(jfl)
(  125)          zgjfl(jfl) = tpjfl(jfl)
(  126)          zgkfl(jfl) = tpkfl(jfl)
(  127)       END DO
(  128)        
(  129)       DO  jind = 1, 4
(  130)          
(  131)          ! for each step we compute the compute the velocity with Lagrange interpolation
(  132)          CALL flo_interp(zgifl,zgjfl,zgkfl,zufl,zvfl,zwfl,jind)
(  133)          
(  134)          ! computation of Runge-Kutta factor
(  135)          
(  136)          DO jfl = 1, jpnfl
(  137)             zrkxfl(jfl,jind) = rdt*zufl(jfl)
(  138)             zrkyfl(jfl,jind) = rdt*zvfl(jfl)
(  139)             zrkzfl(jfl,jind) = rdt*zwfl(jfl)
(  140)          END DO
(  141)          IF( jind /= 4 ) THEN
(  142)             DO jfl = 1, jpnfl
(  143)                zgifl(jfl) = (tpifl(jfl)) + scoef1(jind)*zrkxfl(jfl,jind)
(  144)                zgjfl(jfl) = (tpjfl(jfl)) + scoef1(jind)*zrkyfl(jfl,jind)
(  145)                zgkfl(jfl) = (tpkfl(jfl)) + scoef1(jind)*zrkzfl(jfl,jind)
(  146)             END DO
(  147)          ENDIF
(  148)       END DO
(  149)       DO jind = 1, 4
(  150)          DO jfl = 1, jpnfl
(  151)             tpifl(jfl) = tpifl(jfl) + scoef2(jind)*zrkxfl(jfl,jind)/6.
(  152)             tpjfl(jfl) = tpjfl(jfl) + scoef2(jind)*zrkyfl(jfl,jind)/6.
(  153)             tpkfl(jfl) = tpkfl(jfl) + scoef2(jind)*zrkzfl(jfl,jind)/6.
(  154)          END DO
(  155)       END DO
(  156) 
(  157)    END SUBROUTINE flo_4rk
(  158) 
(  159) 
(  160)    SUBROUTINE flo_interp( pxt , pyt , pzt ,      &
(  161)       &                   pufl, pvfl, pwfl, kind )
(  162)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 4

(  163)       !!                ***  ROUTINE flointerp  ***
(  164)       !!
(  165)       !! ** Purpose :   Interpolation of the velocity on the float position
(  166)       !! 
(  167)       !! ** Method  :   Lagrange interpolation with the 64 neighboring
(  168)       !!      points. This routine is call 4 time at each time step to
(  169)       !!      compute velocity at the date and the position we need to
(  170)       !!      integrated with RK method.
(  171)       !!
(  172)       !!----------------------------------------------------------------------
(  173)       !! * Local declarations
(  174)       INTEGER ::   &
(  175)          kind,   &
(  176)          jfl, jind1, jind2, jind3,   &
(  177)          zsumu, zsumv, zsumw
(  178)       INTEGER , DIMENSION  ( jpnfl  ) ::   &
(  179)          iilu, ijlu, iklu,     &           ! nearest neighbour INDEX-u
(  180)          iilv, ijlv, iklv,     &           ! nearest neighbour INDEX-v
(  181)          iilw, ijlw, iklw                  ! nearest neighbour INDEX-w
(  182)       INTEGER , DIMENSION  ( jpnfl, 4  ) ::   &
(  183)          iidu, ijdu, ikdu,    &            ! 64 nearest neighbour INDEX-u
(  184)          iidv, ijdv, ikdv,    &            ! 64 nearest neighbour INDEX-v
(  185)          iidw, ijdw, ikdw                  ! 64 nearest neighbour INDEX-w
(  186)       REAL(wp) , DIMENSION  ( jpnfl  ) ::   &
(  187)          pxt , pyt , pzt,     &            ! position of the float
(  188)          pufl, pvfl, pwfl                  ! velocity at this position
(  189)       REAL(wp) , DIMENSION  ( jpnfl, 4, 4, 4 ) ::   &
(  190)          ztufl, ztvfl, ztwfl               ! velocity at choosen time step
(  191)       REAL(wp) , DIMENSION  ( jpnfl, 4 ) ::   &
(  192)          zlagxu, zlagyu, zlagzu,   &       ! Lagrange  coefficients
(  193)          zlagxv, zlagyv, zlagzv,   &
(  194)          zlagxw, zlagyw, zlagzw
(  195)       !!---------------------------------------------------------------------
(  196)       
(  197)       ! Interpolation of U velocity
(  198) 
(  199)       ! nearest neightboring point for computation of u       
(  200)       DO jfl = 1, jpnfl
(  201)          iilu(jfl) = INT(pxt(jfl)-.5)
(  202)          ijlu(jfl) = INT(pyt(jfl)-.5)
(  203)          iklu(jfl) = INT(pzt(jfl))
(  204)       END DO
(  205)       
(  206)       !  64 neightboring points for computation of u 
(  207)       DO jind1 = 1, 4
(  208)          DO jfl = 1, jpnfl
(  209)             !  i-direction
(  210)             IF( iilu(jfl) <= 2 ) THEN
(  211)                iidu(jfl,jind1) = jind1
(  212)             ELSE
(  213)                IF( iilu(jfl) >= jpi-1 ) THEN
(  214)                   iidu(jfl,jind1) = jpi + jind1 - 4
(  215)                ELSE
(  216)                   iidu(jfl,jind1) = iilu(jfl) + jind1 - 2
(  217)                ENDIF
(  218)             ENDIF
(  219)             !  j-direction
(  220)             IF( ijlu(jfl) <= 2 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 5

(  221)                ijdu(jfl,jind1) = jind1
(  222)             ELSE
(  223)                IF( ijlu(jfl) >= jpj-1 ) THEN
(  224)                   ijdu(jfl,jind1) = jpj + jind1 - 4
(  225)                ELSE
(  226)                   ijdu(jfl,jind1) = ijlu(jfl) + jind1 - 2
(  227)                ENDIF
(  228)             ENDIF
(  229)             ! k-direction
(  230)             IF( iklu(jfl) <= 2 ) THEN
(  231)                ikdu(jfl,jind1) = jind1
(  232)             ELSE
(  233)                IF( iklu(jfl) >= jpk-1 ) THEN
(  234)                   ikdu(jfl,jind1) = jpk + jind1 - 4
(  235)                ELSE
(  236)                   ikdu(jfl,jind1) = iklu(jfl) + jind1 - 2
(  237)                ENDIF
(  238)             ENDIF
(  239)          END DO
(  240)       END DO
(  241)       
(  242)       ! Lagrange coefficients
(  243)       
(  244)       DO jfl = 1, jpnfl
(  245)          DO jind1 = 1, 4
(  246)             zlagxu(jfl,jind1) = 1.
(  247)             zlagyu(jfl,jind1) = 1.
(  248)             zlagzu(jfl,jind1) = 1.
(  249)          END DO
(  250)       END DO
(  251)       
(  252)       DO jind1 = 1, 4
(  253)          DO jind2 = 1, 4
(  254)             DO jfl= 1, jpnfl
(  255)                IF( jind1 /= jind2 ) THEN
(  256)                   zlagxu(jfl,jind1) = zlagxu(jfl,jind1) * ( pxt(jfl)-(float(iidu(jfl,jind2))+.5) )
(  257)                   zlagyu(jfl,jind1) = zlagyu(jfl,jind1) * ( pyt(jfl)-(float(ijdu(jfl,jind2))) )
(  258)                   zlagzu(jfl,jind1) = zlagzu(jfl,jind1) * ( pzt(jfl)-(float(ikdu(jfl,jind2))) )
(  259)                ENDIF
(  260)             END DO
(  261)          END DO
(  262)       END DO
(  263)       
(  264)       ! velocity when we compute at middle time step
(  265)       
(  266)       DO jfl = 1, jpnfl
(  267)          DO jind1 = 1, 4
(  268)             DO jind2 = 1, 4
(  269)                DO jind3 = 1, 4
(  270)                   ztufl(jfl,jind1,jind2,jind3) =   &
(  271)                      &   (  tcoef1(kind) * ub(iidu(jfl,jind1),ijdu(jfl,jind2),ikdu(jfl,jind3)) +   &
(  272)                      &      tcoef2(kind) * un(iidu(jfl,jind1),ijdu(jfl,jind2),ikdu(jfl,jind3)) )   &
(  273)                      &      / e1u(iidu(jfl,jind1),ijdu(jfl,jind2)) 
(  274)                END DO
(  275)             END DO
(  276)          END DO
(  277)          
(  278)          zsumu = 0.






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 6

(  279)          DO jind1 = 1, 4
(  280)             DO jind2 = 1, 4
(  281)                DO jind3 = 1, 4
(  282)                   zsumu = zsumu + ztufl(jfl,jind1,jind2,jind3) * zlagxu(jfl,jind1) * zlagyu(jfl,jind2)   &
(  283)                      &  * zlagzu(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3) 
(  284)                END DO
(  285)             END DO
(  286)          END DO
(  287)          pufl(jfl) = zsumu
(  288)       END DO
(  289)       
(  290)       ! Interpolation of V velocity 
(  291) 
(  292)       ! nearest neightboring point for computation of v 
(  293)       DO jfl = 1, jpnfl
(  294)          iilv(jfl) = INT(pxt(jfl)-.5)
(  295)          ijlv(jfl) = INT(pyt(jfl)-.5)
(  296)          iklv(jfl) = INT(pzt(jfl))
(  297)       END DO
(  298)       
(  299)       ! 64 neightboring points for computation of v 
(  300)       DO jind1 = 1, 4
(  301)          DO jfl = 1, jpnfl
(  302)             ! i-direction
(  303)             IF( iilv(jfl) <= 2 ) THEN
(  304)                iidv(jfl,jind1) = jind1
(  305)             ELSE
(  306)                IF( iilv(jfl) >= jpi-1 ) THEN
(  307)                   iidv(jfl,jind1) = jpi + jind1 - 4
(  308)                ELSE
(  309)                   iidv(jfl,jind1) = iilv(jfl) + jind1 - 2
(  310)                ENDIF
(  311)             ENDIF
(  312)             ! j-direction
(  313)             IF( ijlv(jfl) <= 2 ) THEN
(  314)                ijdv(jfl,jind1) = jind1
(  315)             ELSE
(  316)                IF( ijlv(jfl) >= jpj-1 ) THEN
(  317)                   ijdv(jfl,jind1) = jpj + jind1 - 4
(  318)                ELSE
(  319)                   ijdv(jfl,jind1) = ijlv(jfl) + jind1 - 2
(  320)                ENDIF
(  321)             ENDIF
(  322)             ! k-direction
(  323)             IF( iklv(jfl) <= 2 ) THEN
(  324)                ikdv(jfl,jind1) = jind1
(  325)             ELSE
(  326)                IF( iklv(jfl) >= jpk-1 ) THEN
(  327)                   ikdv(jfl,jind1) = jpk + jind1 - 4
(  328)                ELSE
(  329)                   ikdv(jfl,jind1) = iklv(jfl) + jind1 - 2
(  330)                ENDIF
(  331)             ENDIF
(  332)          END DO
(  333)       END DO
(  334)       
(  335)       ! Lagrange coefficients
(  336)       






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 7

(  337)       DO jfl = 1, jpnfl
(  338)          DO jind1 = 1, 4
(  339)             zlagxv(jfl,jind1) = 1.
(  340)             zlagyv(jfl,jind1) = 1.
(  341)             zlagzv(jfl,jind1) = 1.
(  342)          END DO
(  343)       END DO
(  344)       
(  345)       DO jind1 = 1, 4
(  346)          DO jind2 = 1, 4
(  347)             DO jfl = 1, jpnfl
(  348)                IF( jind1 /= jind2 ) THEN
(  349)                   zlagxv(jfl,jind1)= zlagxv(jfl,jind1)*(pxt(jfl) - (float(iidv(jfl,jind2))) )
(  350)                   zlagyv(jfl,jind1)= zlagyv(jfl,jind1)*(pyt(jfl) - (float(ijdv(jfl,jind2))+.5) )
(  351)                   zlagzv(jfl,jind1)= zlagzv(jfl,jind1)*(pzt(jfl) - (float(ikdv(jfl,jind2))) )
(  352)                ENDIF
(  353)             END DO
(  354)          END DO
(  355)       END DO
(  356)       
(  357)       ! velocity when we compute at middle time step
(  358)       
(  359)       DO jfl = 1, jpnfl
(  360)          DO jind1 = 1, 4
(  361)             DO jind2 = 1, 4
(  362)                DO jind3 = 1 ,4
(  363)                   ztvfl(jfl,jind1,jind2,jind3)=   &
(  364)                      &   ( tcoef1(kind) * vb(iidv(jfl,jind1),ijdv(jfl,jind2),ikdv(jfl,jind3))  +   &
(  365)                      &     tcoef2(kind) * vn(iidv(jfl,jind1),ijdv(jfl,jind2),ikdv(jfl,jind3)) )    & 
(  366)                      &     / e2v(iidv(jfl,jind1),ijdv(jfl,jind2))
(  367)                END DO
(  368)             END DO
(  369)          END DO
(  370)          
(  371)          zsumv=0.
(  372)          DO jind1 = 1, 4
(  373)             DO jind2 = 1, 4
(  374)                DO jind3 = 1, 4
(  375)                   zsumv = zsumv + ztvfl(jfl,jind1,jind2,jind3) * zlagxv(jfl,jind1) * zlagyv(jfl,jind2)   &
(  376)                      &  * zlagzv(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3)
(  377)                END DO
(  378)             END DO
(  379)          END DO
(  380)          pvfl(jfl) = zsumv
(  381)       END DO
(  382)       
(  383)       ! Interpolation of W velocity
(  384) 
(  385)       ! nearest neightboring point for computation of w 
(  386)       DO jfl = 1, jpnfl
(  387)          iilw(jfl) = INT(pxt(jfl))
(  388)          ijlw(jfl) = INT(pyt(jfl))
(  389)          iklw(jfl) = INT(pzt(jfl)+.5)
(  390)       END DO
(  391)       
(  392)       ! 64 neightboring points for computation of w 
(  393)       DO jind1 = 1, 4
(  394)          DO jfl = 1, jpnfl






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 8

(  395)             ! i-direction
(  396)             IF( iilw(jfl) <= 2 ) THEN
(  397)                iidw(jfl,jind1) = jind1
(  398)             ELSE
(  399)                IF( iilw(jfl) >= jpi-1 ) THEN
(  400)                   iidw(jfl,jind1) = jpi + jind1 - 4
(  401)                ELSE
(  402)                   iidw(jfl,jind1) = iilw(jfl) + jind1 - 2
(  403)                ENDIF
(  404)             ENDIF
(  405)             ! j-direction
(  406)             IF( ijlw(jfl) <= 2 ) THEN
(  407)                ijdw(jfl,jind1) = jind1
(  408)             ELSE
(  409)                IF( ijlw(jfl) >= jpj-1 ) THEN
(  410)                   ijdw(jfl,jind1) = jpj + jind1 - 4
(  411)                ELSE
(  412)                   ijdw(jfl,jind1) = ijlw(jfl) + jind1 - 2
(  413)                ENDIF
(  414)             ENDIF
(  415)             ! k-direction
(  416)             IF( iklw(jfl) <= 2 ) THEN
(  417)                ikdw(jfl,jind1) = jind1
(  418)             ELSE
(  419)                IF( iklw(jfl) >= jpk-1 ) THEN
(  420)                   ikdw(jfl,jind1) = jpk + jind1 - 4
(  421)                ELSE
(  422)                   ikdw(jfl,jind1) = iklw(jfl) + jind1 - 2
(  423)                ENDIF
(  424)             ENDIF
(  425)          END DO
(  426)       END DO
(  427)       DO jind1 = 1, 4
(  428)          DO jfl = 1, jpnfl
(  429)             IF( iklw(jfl) <= 2 ) THEN
(  430)                ikdw(jfl,jind1) = jind1
(  431)             ELSE
(  432)                IF( iklw(jfl) >= jpk-1 ) THEN
(  433)                   ikdw(jfl,jind1) = jpk + jind1 - 4
(  434)                ELSE
(  435)                   ikdw(jfl,jind1) = iklw(jfl) + jind1 - 2
(  436)                ENDIF
(  437)             ENDIF
(  438)          END DO
(  439)       END DO
(  440)       
(  441)       ! Lagrange coefficients  for w interpolation
(  442)       
(  443)       DO jfl = 1, jpnfl
(  444)          DO jind1 = 1, 4
(  445)             zlagxw(jfl,jind1) = 1.
(  446)             zlagyw(jfl,jind1) = 1.
(  447)             zlagzw(jfl,jind1) = 1.
(  448)          END DO
(  449)       END DO
(  450)        
(  451)       DO jind1 = 1, 4
(  452)          DO jind2 = 1, 4






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 9

(  453)             DO jfl = 1, jpnfl
(  454)                IF( jind1 /= jind2 ) THEN
(  455)                   zlagxw(jfl,jind1) = zlagxw(jfl,jind1) * (pxt(jfl) - (float(iidw(jfl,jind2))) )
(  456)                   zlagyw(jfl,jind1) = zlagyw(jfl,jind1) * (pyt(jfl) - (float(ijdw(jfl,jind2))) )
(  457)                   zlagzw(jfl,jind1) = zlagzw(jfl,jind1) * (pzt(jfl) - (float(ikdw(jfl,jind2))-.5) )
(  458)                ENDIF
(  459)             END DO
(  460)          END DO
(  461)       END DO
(  462)       
(  463)       ! velocity w  when we compute at middle time step
(  464)       
(  465)       DO jfl = 1, jpnfl
(  466)          DO jind1 = 1, 4
(  467)             DO jind2 = 1, 4
(  468)                DO jind3 = 1, 4
(  469)                   ztwfl(jfl,jind1,jind2,jind3)=   &
(  470)                      &   ( tcoef1(kind) * wb(iidw(jfl,jind1),ijdw(jfl,jind2),ikdw(jfl,jind3))+   &
(  471)                      &     tcoef2(kind) * wn(iidw(jfl,jind1),ijdw(jfl,jind2),ikdw(jfl,jind3)) )  &
(  472) !!bug e3w instead of fse3 
(  473)                      &   / e3w(ikdw(jfl,jind3))
(  474)                END DO
(  475)             END DO
(  476)          END DO
(  477)          
(  478)          zsumw=0.
(  479)          DO jind1 = 1, 4
(  480)             DO jind2 = 1, 4
(  481)                DO jind3 = 1, 4
(  482)                   zsumw = zsumw + ztwfl(jfl,jind1,jind2,jind3) * zlagxw(jfl,jind1) * zlagyw(jfl,jind2)   &
(  483)                      &  * zlagzw(jfl,jind3) * rcoef(jind1)*rcoef(jind2)*rcoef(jind3)
(  484)                END DO
(  485)             END DO
(  486)          END DO
(  487)          pwfl(jfl) = zsumw
(  488)       END DO
(  489)       
(  490)    END SUBROUTINE flo_interp
(  491) 
(  492) #  else
(  493)    !!----------------------------------------------------------------------
(  494)    !!   Default option                                         Empty module
(  495)    !!----------------------------------------------------------------------
(  496) CONTAINS
(  497)    SUBROUTINE flo_4rk                 ! Empty routine
(  498)    END SUBROUTINE flo_4rk
(  499) #endif
(  500)    
(  501)    !!======================================================================
(  502) END MODULE flo4rk














PGF90 (Version     10.5)          09/23/2010  09:49:47      page 10

(    1) # 1 "flo4rk.F90"
(    1) MODULE flo4rk
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE  flo4rk  ***
(    4)    !! Ocean floats :   trajectory computation using a 4th order Runge-Kutta
(    5)    !!======================================================================
(    6) # 493
(  493)    !!----------------------------------------------------------------------
(  494)    !!   Default option                                         Empty module
(  495)    !!----------------------------------------------------------------------
(  496) CONTAINS
(  497)    SUBROUTINE flo_4rk                 ! Empty routine
(  498)    END SUBROUTINE flo_4rk
(  499) # 500
(  500)    
(  501)    !!======================================================================
(  502) END MODULE flo4rk
