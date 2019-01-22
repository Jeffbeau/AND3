


PGF90 (Version     10.5)          09/23/2010  09:49:25      page 1

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

Filename: obcdyn.F90

(    1) MODULE obcdyn
(    2) #if defined key_obc
(    3)   !!=================================================================================
(    4)   !!                       ***  MODULE  obcdyn  ***
(    5)   !! Ocean dynamics:   Radiation of velocities on each open boundary
(    6)   !!=================================================================================
(    7) 
(    8)   !!---------------------------------------------------------------------------------
(    9)   !!   obc_dyn        : call the subroutine for each open boundary
(   10)   !!   obc_dyn_east   : radiation of the east open boundary velocities
(   11)   !!   obc_dyn_west   : radiation of the west open boundary velocities
(   12)   !!   obc_dyn_north  : radiation of the north open boundary velocities
(   13)   !!   obc_dyn_south  : radiation of the south open boundary velocities
(   14)   !!----------------------------------------------------------------------------------
(   15) 
(   16)   !!----------------------------------------------------------------------------------
(   17)   !! * Modules used
(   18)   USE oce             ! ocean dynamics and tracers
(   19)   USE dom_oce         ! ocean space and time domain
(   20)   USE phycst          ! physical constants
(   21)   USE obc_oce         ! ocean open boundary conditions
(   22)   USE lbclnk          ! ???
(   23)   USE lib_mpp         ! ???
(   24)   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   25)   USE obccli          ! ocean open boundary conditions: climatology
(   26)   USE in_out_manager  ! I/O manager
(   27) 
(   28)   IMPLICIT NONE
(   29)   PRIVATE
(   30) 
(   31)   !! * Accessibility
(   32)   PUBLIC obc_dyn     ! routine called in dynspg_flt (free surface case)
(   33)                      ! routine called in dynnxt.F90 (rigid lid case)
(   34) 
(   35)   !! * Module variables
(   36)   INTEGER ::   ji, jj, jk     ! dummy loop indices
(   37) 
(   38)   INTEGER ::      & ! ... boundary space indices
(   39)      nib   = 1,   & ! nib   = boundary point
(   40)      nibm  = 2,   & ! nibm  = 1st interior point
(   41)      nibm2 = 3,   & ! nibm2 = 2nd interior point
(   42)                     ! ... boundary time indices
(   43)      nit   = 1,   & ! nit    = now
(   44)      nitm  = 2,   & ! nitm   = before
(   45)      nitm2 = 3      ! nitm2  = before-before
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 2

(   47)   REAL(wp) ::   rtaue  , rtauw  , rtaun  , rtaus  ,  &
(   48)                 rtauein, rtauwin, rtaunin, rtausin
(   49) 
(   50)   !!---------------------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)   SUBROUTINE obc_dyn ( kt )
(   55)      !!------------------------------------------------------------------------------
(   56)      !!                      SUBROUTINE obc_dyn
(   57)      !!                     ********************
(   58)      !! ** Purpose :
(   59)      !!      Compute  dynamics (u,v) at the open boundaries.
(   60)      !!      if defined key_dynspg_flt:
(   61)      !!                 this routine is called by dynspg_flt and updates
(   62)      !!                 ua, va which are the actual velocities (not trends)
(   63)      !!      else  (rigid lid case) ,
(   64)      !!                 this routine is called in dynnxt.F routine and updates ua, va.
(   65)      !!
(   66)      !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   67)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   68)      !!      open one (must be done in the param_obc.h90 file).
(   69)      !!
(   70)      !! ** Reference :
(   71)      !!      Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   72)      !!
(   73)      !! History :
(   74)      !!        !  95-03 (J.-M. Molines) Original, SPEM
(   75)      !!        !  97-07 (G. Madec, J.-M. Molines) addition
(   76)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(   77)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   78)      !!----------------------------------------------------------------------
(   79)      !! * Arguments
(   80)      INTEGER, INTENT( in ) ::   kt
(   81) 
(   82)      !!----------------------------------------------------------------------
(   83)      !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   84)      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdyn.F90,v 1.5 2005/12/28 09:25:07 opalod Exp $
(   85)      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   86)      !!----------------------------------------------------------------------
(   87) 
(   88)      ! 0. Local constant initialization
(   89)      ! --------------------------------
(   90) 
(   91)      IF( kt == nit000 .OR. ln_rstart) THEN
(   92)         ! ... Boundary restoring coefficient
(   93)         rtaue = 2. * rdt / rdpeob
(   94)         rtauw = 2. * rdt / rdpwob
(   95)         rtaun = 2. * rdt / rdpnob
(   96)         rtaus = 2. * rdt / rdpsob
(   97)         ! ... Boundary restoring coefficient for inflow ( all boundaries)
(   98)         rtauein = 2. * rdt / rdpein
(   99)         rtauwin = 2. * rdt / rdpwin
(  100)         rtaunin = 2. * rdt / rdpnin
(  101)         rtausin = 2. * rdt / rdpsin
(  102)      END IF
(  103) 
(  104)      IF( lp_obc_east  )   CALL obc_dyn_east ( kt )






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 3

(  105)      IF( lp_obc_west  )   CALL obc_dyn_west ( kt )
(  106)      IF( lp_obc_north )   CALL obc_dyn_north( kt )
(  107)      IF( lp_obc_south )   CALL obc_dyn_south( kt )
(  108) 
(  109)      IF( lk_mpp ) THEN
(  110)         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
(  111)            CALL lbc_lnk( ub, 'U', -1. )
(  112)            CALL lbc_lnk( vb, 'V', -1. )
(  113)         END IF
(  114)         CALL lbc_lnk( ua, 'U', -1. )
(  115)         CALL lbc_lnk( va, 'V', -1. )
(  116)      ENDIF
(  117) 
(  118)   END SUBROUTINE obc_dyn
(  119) 
(  120) 
(  121)   SUBROUTINE obc_dyn_east ( kt )
(  122)      !!------------------------------------------------------------------------------
(  123)      !!                  ***  SUBROUTINE obc_dyn_east  ***
(  124)      !!
(  125)      !! ** Purpose :
(  126)      !!      Apply the radiation algorithm on east OBC velocities ua, va using the
(  127)      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
(  128)      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
(  129)      !!
(  130)      !!  History :
(  131)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  132)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  133)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  134)      !!         ! 00-06 (J.-M. Molines)
(  135)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  136)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  137)      !!------------------------------------------------------------------------------
(  138)      !! * Arguments
(  139)      INTEGER, INTENT( in ) ::   kt
(  140) 
(  141)      !! * Local declaration
(  142)      REAL(wp) ::   z05cx, ztau, zin
(  143) !!DB
(  144)      INTEGER :: jji, jjj, jjk
(  145)      REAL(wp) :: ts1
(  146) 
(  147)      !!------------------------------------------------------------------------------
(  148) 
(  149)      ! 1. First three time steps and more if lfbceast is .TRUE.
(  150)      !    In that case open boundary conditions are FIXED.
(  151)      ! --------------------------------------------------------
(  152) 
(  153)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast .OR. lk_dynspg_exp ) THEN
(  154) 
(  155)         ! 1.1 U zonal velocity
(  156)         ! --------------------
(  157)         DO ji = nie0, nie1
(  158)            DO jk = 1, jpkm1
(  159)               DO jj = 1, jpj
(  160) # if defined key_dynspg_rl
(  161)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uemsk(jj,jk)) +                     &
(  162)                                 uemsk(jj,jk)*( ufoe(jj,jk) - hur (ji,jj) / e2u (ji,jj) &






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 4

(  163)                                 * ( bsfn(ji,jj) - bsfn(ji,jj-1) ) )
(  164) # else
(  165)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uemsk(jj,jk)) + &
(  166)                                 uemsk(jj,jk)*ufoe(jj,jk)
(  167) # endif
(  168)               END DO
(  169)            END DO
(  170)         END DO
(  171) 
(  172) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  173) !!
(  174)         zin = 2. * rdt / (1.0*rday)
(  175) !        do jji = nie0, nie1   !!isolate processor 
(  176) !DBG -- skip 
(  177)         do jji = 1, 0 
(  178)            do ji = nie0-4, nie0-1 !!5 cells - 1
(  179) !restore based on position -- but do not use
(  180)               ztau = zin * exp(float( -(nie0-ji)/4 )) 
(  181)               do jk = 1, jpkm1
(  182)                  do jj = 1, jpj
(  183)                     ua(ji,jj,jk) = ua(ji,jj,jk) + uemsk(jj,jk)*&
(  184)                          zin*(ufoe(jj,jk)-ua(ji,jj,jk)) * umask(ji,jj,jk)
(  185)                  enddo
(  186)               enddo
(  187)            enddo
(  188)         enddo
(  189) 
(  190) 
(  191) 
(  192)         ! 1.2 V meridional velocity
(  193)         ! -------------------------
(  194)         DO ji = nie0+1, nie1+1
(  195)            DO jk = 1, jpkm1
(  196)               DO jj = 1, jpj
(  197)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vemsk(jj,jk)) + &
(  198)                                 vfoe(jj,jk)*vemsk(jj,jk)
(  199)               END DO
(  200)            END DO
(  201)         END DO
(  202) 
(  203)      ELSE
(  204) 
(  205)      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
(  206)      ! -----------------------------------------------------
(  207) 
(  208)         ! 2.1. u-component of the velocity
(  209)         ! ---------------------------------
(  210)         !
(  211)         !          nibm2      nibm      nib
(  212)         !            |   nibm  |   nib   |///
(  213)         !            |    |    |    |    |///
(  214)         !  jj-line --f----v----f----v----f---
(  215)         !            |    |    |    |    |///
(  216)         !            |         |         |///
(  217)         !  jj-line   u    T    u    T    u///
(  218)         !            |         |         |///
(  219)         !            |    |    |    |    |///
(  220)         !          jpieob-2   jpieob-1   jpieob






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 5

(  221)         !                 |         |
(  222)         !              jpieob-1    jpieob
(  223)         !
(  224)         ! ... If free surface formulation:
(  225)         ! ... radiative conditions on the total part + relaxation toward climatology
(  226)         ! ... (jpjedp1, jpjefm1),jpieob
(  227)         DO ji = nie0, nie1
(  228)            DO jk = 1, jpkm1
(  229)               DO jj = 1, jpj
(  230)                  z05cx = u_cxebnd(jj,jk)
(  231)                  z05cx = z05cx / e1t(ji,jj)
(  232)                  z05cx = min( z05cx, 1. )
(  233)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  234)         !           > 0, outflow zin=1, ztau=rtaue
(  235)                  zin = sign( 1., z05cx )
(  236)                  zin = 0.5*( zin + abs(zin) )
(  237)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
(  238)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  239)                  z05cx = z05cx * zin
(  240)         ! ... update ua with radiative or climatological velocity
(  241)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uemsk(jj,jk) ) +          &
(  242)                                 uemsk(jj,jk) * (  ( 1. - z05cx - ztau )         &
(  243)                                 * uebnd(jj,jk,nib ,nitm) + 2.*z05cx               &
(  244)                                 * uebnd(jj,jk,nibm,nit ) + ztau * ufoe (jj,jk) )  &
(  245)                                 / (1. + z05cx)
(  246)               END DO
(  247)            END DO
(  248)         END DO
(  249) # if defined key_dynspg_rl
(  250)         ! ... ua must be a baroclinic velocity uclie()
(  251)         CALL obc_cli( ua, uclie, nie0, nie1, 0, jpj )
(  252) 
(  253)         ! ... add the correct barotropic radiative velocity (calculated from bsfn) to the
(  254)         !     baroclinc velocity uclie() to have the total velocity
(  255)         DO ji = nie0, nie1
(  256)            DO jk = 1, jpkm1
(  257)               DO jj = 1, jpj
(  258)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uemsk(jj,jk) ) +                     &
(  259)                                 uemsk(jj,jk) * ( uclie(jj,jk) -  hur (ji,jj) / e2u (ji,jj) &
(  260)                                 * ( bsfn(ji,jj) - bsfn(ji,jj-1) ) )
(  261)               END DO
(  262)            END DO
(  263)         END DO
(  264) # endif
(  265)         ! 2.2 v-component of the velocity
(  266)         ! -------------------------------
(  267)         !
(  268)         !          nibm2       nibm     nib
(  269)         !            |   nibm  |   nib///|///
(  270)         !            |    |    |    |////|///
(  271)         !  jj-line --v----f----v----f----v---
(  272)         !            |    |    |    |////|///
(  273)         !            |    |    |    |////|///
(  274)         !            | jpieob-1 |  jpieob /|///
(  275)         !            |         |         |
(  276)         !         jpieob-1    jpieob     jpieob+1
(  277)         !
(  278)         ! ... radiative condition






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 6

(  279)         ! ... (jpjedp1, jpjefm1), jpieob+1
(  280)         DO ji = nie0+1, nie1+1
(  281)            DO jk = 1, jpkm1
(  282)               DO jj = 1, jpj
(  283)                  z05cx = v_cxebnd(jj,jk)
(  284)                  z05cx = z05cx / e1f(ji-1,jj)
(  285)                  z05cx = min( z05cx, 1. )
(  286)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  287)         !           > 0, outflow zin=1, ztau=rtaue
(  288)                  zin = sign( 1., z05cx )
(  289)                  zin = 0.5*( zin + abs(zin) )
(  290)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
(  291)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  292)                  z05cx = z05cx * zin
(  293)         ! ... update va with radiative or climatological velocity
(  294)                  va(ji,jj,jk) = va(ji,jj,jk) * (1. - vemsk(jj,jk) ) +          &
(  295)                                 vemsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  296)                                 * vebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  297)                                 * vebnd(jj,jk,nibm,nit ) + ztau * vfoe(jj,jk) )  &
(  298)                                 / (1. + z05cx)
(  299)               END DO
(  300)            END DO
(  301)         END DO
(  302) 
(  303)      END IF
(  304) 
(  305)   END SUBROUTINE obc_dyn_east
(  306) 
(  307) 
(  308)   SUBROUTINE obc_dyn_west ( kt )
(  309)      !!------------------------------------------------------------------------------
(  310)      !!                  ***  SUBROUTINE obc_dyn_west  ***
(  311)      !!
(  312)      !! ** Purpose :
(  313)      !!      Apply the radiation algorithm on west OBC velocities ua, va using the
(  314)      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
(  315)      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
(  316)      !!
(  317)      !!  History :
(  318)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  319)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  320)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  321)      !!         ! 00-06 (J.-M. Molines)
(  322)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  323)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  324)      !!------------------------------------------------------------------------------
(  325)      !! * Arguments
(  326)      INTEGER, INTENT( in ) ::   kt
(  327) 
(  328)      !! * Local declaration
(  329)      REAL(wp) ::   z05cx, ztau, zin
(  330) !!DB
(  331)      INTEGER :: jji, jjj, jjk
(  332) 
(  333)      !!------------------------------------------------------------------------------
(  334) 
(  335)      ! 1. First three time steps and more if lfbcwest is .TRUE.
(  336)      !    In that case open boundary conditions are FIXED.






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 7

(  337)      ! --------------------------------------------------------
(  338) 
(  339)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest .OR. lk_dynspg_exp ) THEN
(  340) 
(  341)         ! 1.1 U zonal velocity
(  342)         ! ---------------------
(  343)         DO ji = niw0, niw1
(  344)            DO jk = 1, jpkm1
(  345)               DO jj = 1, jpj
(  346) # if defined key_dynspg_rl
(  347)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uwmsk(jj,jk)) +                     &
(  348)                                 uwmsk(jj,jk)*( ufow(jj,jk) - hur (ji,jj) / e2u (ji,jj) &
(  349)                                 * ( bsfn(ji,jj) - bsfn(ji,jj-1) ) )
(  350) # else
(  351)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uwmsk(jj,jk)) + &
(  352)                                 uwmsk(jj,jk)*ufow(jj,jk)
(  353) # endif
(  354)               END DO
(  355)            END DO
(  356)         END DO
(  357) 
(  358) 
(  359) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  360) !!
(  361)         zin = 2. * rdt / (1.0*rday)
(  362) !        do jji = niw0, niw1   !!isolate processor 
(  363)        do jji = 1, 0    !!isolate processor 
(  364)            do ji = niw0+1, niw0+4 !!5 cells - 1
(  365) !!restore based on position -- but do not use
(  366)               ztau = zin * exp(float( -(ji-niw0)/4 )) 
(  367)               do jk = 1, jpkm1
(  368)                  do jj = 1, jpj
(  369)                     ua(ji,jj,jk) = ua(ji,jj,jk) + uwmsk(jj,jk)*&
(  370)                          zin*(ufow(jj,jk)-ua(ji,jj,jk)) * umask(ji,jj,jk)
(  371)                  enddo
(  372)               enddo
(  373)            enddo
(  374)         enddo
(  375) 
(  376) 
(  377) 
(  378)         ! 1.2 V meridional velocity
(  379)         ! -------------------------
(  380)         DO ji = niw0, niw1
(  381)            DO jk = 1, jpkm1
(  382)               DO jj = 1, jpj
(  383)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vwmsk(jj,jk)) + &
(  384)                                 vfow(jj,jk)*vwmsk(jj,jk)
(  385)               END DO
(  386)            END DO
(  387)         END DO
(  388) 
(  389)      ELSE
(  390) 
(  391)      ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
(  392)      ! -----------------------------------------------------
(  393) 
(  394)         ! 2.1. u-component of the velocity






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 8

(  395)         ! ---------------------------------
(  396)         !
(  397)         !        nib       nibm     nibm2
(  398)         !      ///|   nib   |   nibm  |
(  399)         !      ///|    |    |    |    |
(  400)         !      ---f----v----f----v----f-- jj-line
(  401)         !      ///|    |    |    |    |
(  402)         !      ///|         |         |
(  403)         !      ///u    T    u    T    u   jj-line
(  404)         !      ///|         |         |
(  405)         !      ///|    |    |    |    |
(  406)         !       jpiwob    jpiwob+1    jpiwob+2
(  407)         !              |         |
(  408)         !            jpiwob+1    jpiwob+2
(  409)         !
(  410)         ! ... If free surface formulation:
(  411)         ! ... radiative conditions on the total part + relaxation toward climatology
(  412)         ! ... (jpjwdp1, jpjwfm1), jpiwob
(  413)         DO ji = niw0, niw1
(  414)            DO jk = 1, jpkm1
(  415)               DO jj = 1, jpj
(  416)                  z05cx = u_cxwbnd(jj,jk)
(  417)                  z05cx = z05cx / e1t(ji+1,jj)
(  418)                  z05cx = max( z05cx, -1. )
(  419)         ! ... z05c  > 0, inflow  zin=0, ztau=1
(  420)         !          =< 0, outflow zin=1, ztau=rtauw
(  421)                  zin = sign( 1., -1. * z05cx )
(  422)                  zin = 0.5*( zin + abs(zin) )
(  423)                  ztau = (1.-zin )* rtauwin + zin * rtauw
(  424)                  z05cx = z05cx * zin
(  425)         ! ... update un with radiative or climatological velocity
(  426)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uwmsk(jj,jk) ) +          &
(  427)                                 uwmsk(jj,jk) * ( ( 1. + z05cx - ztau )          &
(  428)                                 * uwbnd(jj,jk,nib ,nitm) - 2.*z05cx               &
(  429)                                 * uwbnd(jj,jk,nibm,nit ) + ztau  * ufow (jj,jk) ) &
(  430)                                 / (1. - z05cx)
(  431)               END DO
(  432)            END DO
(  433)         END DO
(  434) # if defined key_dynspg_rl
(  435)         ! ... ua must be a baroclinic velocity ucliw()
(  436)         CALL obc_cli( ua, ucliw, niw0, niw1, 0, jpj )
(  437) 
(  438)         ! ... add the correct barotropic radiative velocity (calculated from bsfn)
(  439)         !     to the baroclinc velocity ucliw() to have the total velocity
(  440)         DO ji = niw0, niw1
(  441)            DO jk = 1, jpkm1
(  442)               DO jj = 1, jpj
(  443)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uwmsk(jj,jk) ) +                    &
(  444)                                 uwmsk(jj,jk)*( ucliw(jj,jk) - hur (ji,jj) / e2u (ji,jj)   &
(  445)                                 * ( bsfn(ji,jj) - bsfn(ji,jj-1) ) )
(  446)               END DO
(  447)            END DO
(  448)         END DO
(  449) # endif
(  450) 
(  451)         ! 2.2 v-component of the velocity
(  452)         ! -------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 9

(  453)         !
(  454)         !    nib       nibm     nibm2
(  455)         !  ///|///nib   |   nibm  |  nibm2
(  456)         !  ///|////|    |    |    |    |    |
(  457)         !  ---v----f----v----f----v----f----v-- jj-line
(  458)         !  ///|////|    |    |    |    |    |
(  459)         !  ///|////|    |    |    |    |    |
(  460)         ! jpiwob     jpiwob+1    jpiwob+2
(  461)         !          |         |         |
(  462)         !        jpiwob   jpiwob+1   jpiwob+2
(  463)         !
(  464)         ! ... radiative condition plus Raymond-Kuo
(  465)         ! ... (jpjwdp1, jpjwfm1),jpiwob
(  466)         DO ji = niw0, niw1
(  467)            DO jk = 1, jpkm1
(  468)               DO jj = 1, jpj
(  469)                  z05cx = v_cxwbnd(jj,jk)
(  470)                  z05cx = z05cx / e1f(ji,jj)
(  471)                  z05cx = max( z05cx, -1. )
(  472)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  473)         !          =< 0, outflow zin=1, ztau=rtauw
(  474)                  zin = sign( 1., -1. * z05cx )
(  475)                  zin = 0.5*( zin + abs(zin) )
(  476)                  ztau = (1.-zin )*rtauwin + zin * rtauw
(  477)                  z05cx = z05cx * zin
(  478)         ! ... update va with radiative or climatological velocity
(  479)                  va(ji,jj,jk) = va(ji,jj,jk) * (1. - vwmsk(jj,jk) ) +          &
(  480)                                 vwmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  481)                                 * vwbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  482)                                 * vwbnd(jj,jk,nibm,nit ) + ztau * vfow (jj,jk) ) &
(  483)                                 / (1. - z05cx)
(  484)                END DO
(  485)             END DO
(  486)         END DO
(  487) 
(  488)      END IF
(  489) 
(  490)   END SUBROUTINE obc_dyn_west
(  491) 
(  492)   SUBROUTINE obc_dyn_north ( kt )
(  493)      !!------------------------------------------------------------------------------
(  494)      !!                     SUBROUTINE obc_dyn_north
(  495)      !!                    *************************
(  496)      !! ** Purpose :
(  497)      !!      Apply the radiation algorithm on north OBC velocities ua, va using the
(  498)      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
(  499)      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
(  500)      !!
(  501)      !!  History :
(  502)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  503)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  504)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  505)      !!         ! 00-06 (J.-M. Molines)
(  506)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  507)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  508)      !!------------------------------------------------------------------------------
(  509)      !! * Arguments
(  510)      INTEGER, INTENT( in ) ::   kt






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 10

(  511) 
(  512)      !! * Local declaration
(  513)      REAL(wp) ::   z05cx, ztau, zin
(  514)      !!------------------------------------------------------------------------------
(  515) 
(  516)      ! 1. First three time steps and more if lfbcnorth is .TRUE.
(  517)      !    In that case open boundary conditions are FIXED.
(  518)      ! ---------------------------------------------------------
(  519) 
(  520)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth  .OR. lk_dynspg_exp ) THEN
(  521) 
(  522)         ! 1.1 U zonal velocity
(  523)         ! --------------------
(  524)         DO jj = njn0+1, njn1+1
(  525)            DO jk = 1, jpkm1
(  526)               DO ji = 1, jpi
(  527)                  ua(ji,jj,jk)= ua(ji,jj,jk) * (1.-unmsk(ji,jk)) + &
(  528)                                ufon(ji,jk)*unmsk(ji,jk)
(  529)               END DO
(  530)            END DO
(  531)         END DO
(  532) 
(  533)         ! 1.2 V meridional velocity
(  534)         ! -------------------------
(  535)         DO jj = njn0, njn1
(  536)            DO jk = 1, jpkm1
(  537)               DO ji = 1, jpi
(  538) # if defined key_dynspg_rl
(  539)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vnmsk(ji,jk)) +                       &
(  540)                                vnmsk(ji,jk) * ( vfon(ji,jk) + hvr (ji,jj) / e1v (ji,jj) &
(  541)                                * ( bsfn(ji,jj) - bsfn(ji-1,jj) ) )
(  542) # else
(  543)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vnmsk(ji,jk)) + &
(  544)                                vfon(ji,jk)*vnmsk(ji,jk)
(  545) # endif
(  546)               END DO
(  547)            END DO
(  548)         END DO
(  549) 
(  550)      ELSE
(  551) 
(  552)      ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
(  553)      ! ------------------------------------------------------
(  554) 
(  555)         ! 2.1. u-component of the velocity
(  556)         ! --------------------------------
(  557)         !
(  558)         !            ji-row
(  559)         !              |
(  560)         !       nib ///u//////  jpjnob + 1
(  561)         !         /////|//////
(  562)         !     nib -----f-----   jpjnob
(  563)         !              |
(  564)         !      nibm--  u   ---- jpjnob
(  565)         !              |
(  566)         !    nibm -----f-----   jpjnob-1
(  567)         !              |
(  568)         !     nibm2--  u   ---- jpjnob-1






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 11

(  569)         !              |
(  570)         !   nibm2 -----f-----   jpjnob-2
(  571)         !              |
(  572)         !
(  573)         ! ... radiative condition
(  574)         ! ... jpjnob+1,(jpindp1, jpinfm1)
(  575)         DO jj = njn0+1, njn1+1
(  576)            DO jk = 1, jpkm1
(  577)               DO ji = 1, jpi
(  578)                  z05cx= u_cynbnd(ji,jk)
(  579)                  z05cx = z05cx / e2f(ji, jj-1)
(  580)                  z05cx = min( z05cx, 1. )
(  581)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  582)         !           > 0, outflow zin=1, ztau=rtaun
(  583)                  zin = sign( 1., z05cx )
(  584)                  zin = 0.5*( zin + abs(zin) )
(  585)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  586)                  ztau = (1.-zin ) * rtaunin  + zin * rtaun
(  587)         ! ... for u, when inflow, ufon is prescribed
(  588)                  z05cx = z05cx * zin
(  589)         ! ... update un with radiative or climatological velocity
(  590)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-unmsk(ji,jk)) +             &
(  591)                                 unmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  592)                                 * unbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  593)                                 * unbnd(ji,jk,nibm,nit ) + ztau * ufon (ji,jk) ) &
(  594)                                 / (1. + z05cx)
(  595)               END DO
(  596)            END DO
(  597)         END DO
(  598) 
(  599)         ! 2.2 v-component of the velocity
(  600)         ! -------------------------------
(  601)         !
(  602)         !                ji-row    ji-row
(  603)         !              |         |
(  604)         !         /////|/////////////////
(  605)         !    nib  -----f----v----f----  jpjnob
(  606)         !              |         |
(  607)         !      nib  -  u -- T -- u ---- jpjnob
(  608)         !              |         |
(  609)         !   nibm  -----f----v----f----  jpjnob-1
(  610)         !              |         |
(  611)         !     nibm --  u -- T -- u ---  jpjnob-1
(  612)         !              |         |
(  613)         !   nibm2 -----f----v----f----  jpjnob-2
(  614)         !              |         |
(  615)         !
(  616)         ! ... If rigidlid formulation:
(  617)         ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  618)         ! ... If free surface formulation:
(  619)         ! ... radiative conditions on the total part + relaxation toward climatology
(  620)         ! ... jpjnob,(jpindp1, jpinfm1)
(  621)         DO jj = njn0, njn1
(  622)            DO jk = 1, jpkm1
(  623)               DO ji = 1, jpi
(  624)         ! ... 2* gradj(v) (T-point i=nibm, time mean)
(  625)                  z05cx = v_cynbnd(ji,jk)
(  626)                  z05cx = z05cx / e2t(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 12

(  627)                  z05cx = min( z05cx, 1. )
(  628)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  629)         !           > 0, outflow zin=1, ztau=rtaun
(  630)                  zin = sign( 1., z05cx )
(  631)                  zin = 0.5*( zin + abs(zin) )
(  632)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  633)                  ztau = (1.-zin ) * rtaunin + zin * rtaun
(  634)                  z05cx = z05cx * zin
(  635)         ! ... update va with radiative or climatological velocity
(  636)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vnmsk(ji,jk)) +             &
(  637)                                 vnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  638)                                 * vnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  639)                                 * vnbnd(ji,jk,nibm,nit ) + ztau * vfon (ji,jk) ) &
(  640)                                 / (1. + z05cx)
(  641)               END DO
(  642)            END DO
(  643)         END DO
(  644) # if defined key_dynspg_rl
(  645)         ! ... va must be a baroclinic velocity vclin()
(  646)         CALL obc_cli( va, vclin, njn0, njn1, 1, jpi )
(  647) 
(  648)         ! ... add the correct barotropic radiative velocity (calculated from bsfn)
(  649)         !     to the baroclinc velocity vclin() to have the total velocity
(  650)         DO jj = njn0, njn1
(  651)            DO jk = 1, jpkm1
(  652)               DO ji = 1, jpi
(  653)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vnmsk(ji,jk)) +                         &
(  654)                                 vnmsk(ji,jk) * ( vclin(ji,jk) +  hvr (ji,jj) / e1v (ji,jj) &
(  655)                                 * ( bsfn(ji,jj) - bsfn(ji-1,jj) ) )
(  656)               END DO
(  657)            END DO
(  658)         END DO
(  659) # endif
(  660)      END IF
(  661) 
(  662)   END SUBROUTINE obc_dyn_north
(  663) 
(  664)   SUBROUTINE obc_dyn_south ( kt )
(  665)      !!------------------------------------------------------------------------------
(  666)      !!                     SUBROUTINE obc_dyn_south
(  667)      !!                    *************************
(  668)      !! ** Purpose :
(  669)      !!      Apply the radiation algorithm on south OBC velocities ua, va using the
(  670)      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
(  671)      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
(  672)      !!
(  673)      !!  History :
(  674)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  675)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  676)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  677)      !!         ! 00-06 (J.-M. Molines)
(  678)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  679)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  680)      !!------------------------------------------------------------------------------
(  681)      !! * Arguments
(  682)      INTEGER, INTENT( in ) ::   kt
(  683) 
(  684)      !! * Local declaration






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 13

(  685)      REAL(wp) ::   z05cx, ztau, zin
(  686) !!DB
(  687)      INTEGER :: jji, jjj, jjk
(  688) 
(  689) 
(  690)      !!------------------------------------------------------------------------------
(  691)      !!  OPA 8.5, LODYC-IPSL (2002)
(  692)      !!------------------------------------------------------------------------------
(  693) 
(  694)      ! 1. First three time steps and more if lfbcsouth is .TRUE.
(  695)      !    In that case open boundary conditions are FIXED.
(  696)      ! ---------------------------------------------------------
(  697) 
(  698)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth  .OR. lk_dynspg_exp ) THEN
(  699) 
(  700)         ! 1.1 U zonal velocity
(  701)         ! --------------------
(  702)         DO jj = njs0, njs1
(  703)            DO jk = 1, jpkm1
(  704)               DO ji = 1, jpi
(  705)                  ua(ji,jj,jk)= ua(ji,jj,jk) * (1.-usmsk(ji,jk)) + &
(  706)                                usmsk(ji,jk) * ufos(ji,jk)
(  707)               END DO
(  708)            END DO
(  709)         END DO
(  710) 
(  711) 
(  712) 
(  713)         ! 1.2 V meridional velocity
(  714)         ! -------------------------
(  715)         DO jj = njs0, njs1
(  716)            DO jk = 1, jpkm1
(  717)               DO ji = 1, jpi
(  718) # if defined key_dynspg_rl
(  719)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vsmsk(ji,jk)) +                      &
(  720)                                vsmsk(ji,jk) * (vfos(ji,jk) + hvr (ji,jj) / e1v (ji,jj) &
(  721)                                * ( bsfn(ji,jj) - bsfn(ji-1,jj) ) )
(  722) # else
(  723)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vsmsk(ji,jk)) + &
(  724)                                vsmsk(ji,jk) * vfos(ji,jk)
(  725) # endif
(  726)               END DO
(  727)            END DO
(  728)         END DO
(  729) 
(  730) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  731) !!
(  732)         zin = 2. * rdt / (1.0*rday)
(  733) !        do jji = njs0, njs1   !!isolate processor 
(  734)         do jji = 1, 0   !!isolate processor 
(  735)            do jj = njs0+1, njs0+4 !!5 cells - 1
(  736) !!restore based on position -- but do not use
(  737)               ztau = zin * exp(float( -(jj-njs0)/4 )) 
(  738)               do jk = 1, jpkm1
(  739)                  do ji = 1, jpi
(  740)                     va(ji,jj,jk) = va(ji,jj,jk) + vsmsk(ji,jk)*&
(  741)                          zin*(vfos(ji,jk)-va(ji,jj,jk)) * vmask(ji,jj,jk)
(  742)                  enddo






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 14

(  743)               enddo
(  744)            enddo
(  745)         enddo
(  746) 
(  747) 
(  748) 
(  749) 
(  750) 
(  751) 
(  752)      ELSE
(  753) 
(  754)      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
(  755)      ! ------------------------------------------------------
(  756) 
(  757)         ! 2.1. u-component of the velocity
(  758)         ! --------------------------------
(  759)         !
(  760)         !            ji-row
(  761)         !              |
(  762)         !   nibm2 -----f-----   jpjsob +2
(  763)         !              |
(  764)         !    nibm2 --  u   ---- jpjsob +2
(  765)         !              |
(  766)         !    nibm -----f-----   jpjsob +1
(  767)         !              |
(  768)         !    nibm  --  u   ---- jpjsob +1
(  769)         !              |
(  770)         !    nib  -----f-----   jpjsob
(  771)         !         /////|//////
(  772)         !    nib   ////u/////   jpjsob
(  773)         !
(  774)         ! ... radiative condition plus Raymond-Kuo
(  775)         ! ... jpjsob,(jpisdp1, jpisfm1)
(  776)         DO jj = njs0, njs1
(  777)            DO jk = 1, jpkm1
(  778)               DO ji = 1, jpi
(  779)                  z05cx= u_cysbnd(ji,jk)
(  780)                  z05cx = z05cx / e2f(ji, jj)
(  781)                  z05cx = max( z05cx, -1. )
(  782)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  783)         !          =< 0, outflow zin=1, ztau=rtaus
(  784)                  zin = sign( 1., -1. * z05cx )
(  785)                  zin = 0.5*( zin + abs(zin) )
(  786)                  ztau = (1.-zin ) * rtausin + zin * rtaus
(  787)                  z05cx = z05cx * zin
(  788)         ! ... update ua with radiative or climatological velocity
(  789)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-usmsk(ji,jk)) +              &
(  790)                                 usmsk(ji,jk) * (  ( 1. + z05cx - ztau )         &
(  791)                                 * usbnd(ji,jk,nib ,nitm) - 2.*z05cx               &
(  792)                                 * usbnd(ji,jk,nibm,nit ) + ztau * ufos (ji,jk) )  &
(  793)                                 / (1. - z05cx)
(  794)               END DO
(  795)            END DO
(  796)         END DO
(  797) 
(  798)         ! 2.2 v-component of the velocity
(  799)         ! -------------------------------
(  800)         !






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 15

(  801)         !                ji-row    ji-row
(  802)         !              |         |
(  803)         !  nibm2  -----f----v----f----  jpjsob+2
(  804)         !              |         |
(  805)         !    nibm   -  u -- T -- u ---- jpjsob+2
(  806)         !              |         |
(  807)         !   nibm  -----f----v----f----  jpjsob+1
(  808)         !              |         |
(  809)         !   nib    --  u -- T -- u ---  jpjsob+1
(  810)         !              |         |
(  811)         !   nib   -----f----v----f----  jpjsob
(  812)         !         /////////////////////
(  813)         !
(  814)         ! ... If rigidlid formulation:
(  815)         ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  816)         ! ... If free surface formulation:
(  817)         ! ... radiative conditions on the total part + relaxation toward climatology
(  818)         ! ... jpjsob,(jpisdp1,jpisfm1)
(  819)         DO jj = njs0, njs1
(  820)            DO jk = 1, jpkm1
(  821)               DO ji = 1, jpi
(  822)                  z05cx = v_cysbnd(ji,jk)
(  823)                  z05cx = z05cx / e2t(ji,jj+1)
(  824)                  z05cx = max( z05cx, -1. )
(  825)         ! ... z05c > 0, inflow  zin=0, ztau=1
(  826)         !         =< 0, outflow zin=1, ztau=rtaus
(  827)                  zin = sign( 1., -1. * z05cx )
(  828)                  zin = 0.5*( zin + abs(zin) )
(  829)                  ztau = (1.-zin )*rtausin + zin * rtaus
(  830)                  z05cx = z05cx * zin
(  831)         ! ... update va with radiative or climatological velocity
(  832)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vsmsk(ji,jk)) +             &
(  833)                                 vsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
(  834)                                 * vsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
(  835)                                 * vsbnd(ji,jk,nibm,nit ) + ztau * vfos (ji,jk) ) &
(  836)                                 / (1. - z05cx)
(  837)               END DO
(  838)            END DO
(  839)         END DO
(  840) # if defined key_dynspg_rl
(  841)         ! ... va must be a baroclinic velocity vclis()
(  842)         CALL obc_cli( va, vclis, njs0, njs1, 1, jpi )
(  843) 
(  844)         ! ... add the correct barotropic radiative velocity (calculated from bsfn)
(  845)         !     to the baroclinic velocity vclis() to have the total velocity
(  846)         DO jj = njs0, njs1
(  847)            DO jk = 1, jpkm1
(  848)               DO ji = 1, jpi
(  849)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vsmsk(ji,jk)) +                         &
(  850)                                 vsmsk(ji,jk) * ( vclis(ji,jk) +  hvr (ji,jj) / e1v (ji,jj) &
(  851)                                 * ( bsfn(ji,jj) - bsfn(ji-1,jj) ) )
(  852)               END DO
(  853)            END DO
(  854)         END DO
(  855) # endif
(  856)      END IF
(  857) 
(  858)   END SUBROUTINE obc_dyn_south






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 16

(  859) #else
(  860)   !!=================================================================================
(  861)   !!                       ***  MODULE  obcdyn  ***
(  862)   !! Ocean dynamics:   Radiation of velocities on each open boundary
(  863)   !!=================================================================================
(  864) CONTAINS
(  865) 
(  866)   SUBROUTINE obc_dyn
(  867)                              ! No open boundaries ==> empty routine
(  868)   END SUBROUTINE obc_dyn
(  869) #endif
(  870) 
(  871) END MODULE obcdyn



















































PGF90 (Version     10.5)          09/23/2010  09:49:25      page 17

(    1) # 1 "obcdyn.F90"
(    1) MODULE obcdyn
(    2) # 3
(    3)   !!=================================================================================
(    4)   !!                       ***  MODULE  obcdyn  ***
(    5)   !! Ocean dynamics:   Radiation of velocities on each open boundary
(    6)   !!=================================================================================
(    7) 
(    8)   !!---------------------------------------------------------------------------------
(    9)   !!   obc_dyn        : call the subroutine for each open boundary
(   10)   !!   obc_dyn_east   : radiation of the east open boundary velocities
(   11)   !!   obc_dyn_west   : radiation of the west open boundary velocities
(   12)   !!   obc_dyn_north  : radiation of the north open boundary velocities
(   13)   !!   obc_dyn_south  : radiation of the south open boundary velocities
(   14)   !!----------------------------------------------------------------------------------
(   15) 
(   16)   !!----------------------------------------------------------------------------------
(   17)   !! * Modules used
(   18)   USE oce             ! ocean dynamics and tracers
(   19)   USE dom_oce         ! ocean space and time domain
(   20)   USE phycst          ! physical constants
(   21)   USE obc_oce         ! ocean open boundary conditions
(   22)   USE lbclnk          ! ???
(   23)   USE lib_mpp         ! ???
(   24)   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   25)   USE obccli          ! ocean open boundary conditions: climatology
(   26)   USE in_out_manager  ! I/O manager
(   27) 
(   28)   IMPLICIT NONE
(   29)   PRIVATE
(   30) 
(   31)   !! * Accessibility
(   32)   PUBLIC obc_dyn     ! routine called in dynspg_flt (free surface case)
(   33)                      ! routine called in dynnxt.F90 (rigid lid case)
(   34) 
(   35)   !! * Module variables
(   36)   INTEGER ::   ji, jj, jk     ! dummy loop indices
(   37) 
(   38)   INTEGER ::      & ! ... boundary space indices
(   39)      nib   = 1,   & ! nib   = boundary point
(   40)      nibm  = 2,   & ! nibm  = 1st interior point
(   41)      nibm2 = 3,   & ! nibm2 = 2nd interior point
(   42)                     ! ... boundary time indices
(   43)      nit   = 1,   & ! nit    = now
(   44)      nitm  = 2,   & ! nitm   = before
(   45)      nitm2 = 3      ! nitm2  = before-before
(   46) 
(   47)   REAL(wp) ::   rtaue  , rtauw  , rtaun  , rtaus  ,  &
(   48)                 rtauein, rtauwin, rtaunin, rtausin
(   49) 
(   50)   !!---------------------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)   SUBROUTINE obc_dyn ( kt )
(   55)      !!------------------------------------------------------------------------------
(   56)      !!                      SUBROUTINE obc_dyn
(   57)      !!                     ********************






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 18

(   58)      !! ** Purpose :
(   59)      !!      Compute  dynamics (u,v) at the open boundaries.
(   60)      !!      if defined key_dynspg_flt:
(   61)      !!                 this routine is called by dynspg_flt and updates
(   62)      !!                 ua, va which are the actual velocities (not trends)
(   63)      !!      else  (rigid lid case) ,
(   64)      !!                 this routine is called in dynnxt.F routine and updates ua, va.
(   65)      !!
(   66)      !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   67)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   68)      !!      open one (must be done in the param_obc.h90 file).
(   69)      !!
(   70)      !! ** Reference :
(   71)      !!      Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   72)      !!
(   73)      !! History :
(   74)      !!        !  95-03 (J.-M. Molines) Original, SPEM
(   75)      !!        !  97-07 (G. Madec, J.-M. Molines) addition
(   76)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(   77)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   78)      !!----------------------------------------------------------------------
(   79)      !! * Arguments
(   80)      INTEGER, INTENT( in ) ::   kt
(   81) 
(   82)      !!----------------------------------------------------------------------
(   83)      !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   84)      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdyn.F90,v 1.5 2005/12/28 09:25:07 opalod Exp $
(   85)      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   86)      !!----------------------------------------------------------------------
(   87) 
(   88)      ! 0. Local constant initialization
(   89)      ! --------------------------------
(   90) 
(   91)      IF( kt == nit000 .OR. ln_rstart) THEN
(   92)         ! ... Boundary restoring coefficient
(   93)         rtaue = 2. * rdt / rdpeob
(   94)         rtauw = 2. * rdt / rdpwob
(   95)         rtaun = 2. * rdt / rdpnob
(   96)         rtaus = 2. * rdt / rdpsob
(   97)         ! ... Boundary restoring coefficient for inflow ( all boundaries)
(   98)         rtauein = 2. * rdt / rdpein
(   99)         rtauwin = 2. * rdt / rdpwin
(  100)         rtaunin = 2. * rdt / rdpnin
(  101)         rtausin = 2. * rdt / rdpsin
(  102)      END IF
(  103) 
(  104)      IF( lp_obc_east  )   CALL obc_dyn_east ( kt )
(  105)      IF( lp_obc_west  )   CALL obc_dyn_west ( kt )
(  106)      IF( lp_obc_north )   CALL obc_dyn_north( kt )
(  107)      IF( lp_obc_south )   CALL obc_dyn_south( kt )
(  108) 
(  109)      IF( lk_mpp ) THEN
(  110)         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
(  111)            CALL lbc_lnk( ub, 'U', -1. )
(  112)            CALL lbc_lnk( vb, 'V', -1. )
(  113)         END IF
(  114)         CALL lbc_lnk( ua, 'U', -1. )
(  115)         CALL lbc_lnk( va, 'V', -1. )






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 19

(  116)      ENDIF
(  117) 
(  118)   END SUBROUTINE obc_dyn
(  119) 
(  120) 
(  121)   SUBROUTINE obc_dyn_east ( kt )
(  122)      !!------------------------------------------------------------------------------
(  123)      !!                  ***  SUBROUTINE obc_dyn_east  ***
(  124)      !!
(  125)      !! ** Purpose :
(  126)      !!      Apply the radiation algorithm on east OBC velocities ua, va using the
(  127)      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
(  128)      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
(  129)      !!
(  130)      !!  History :
(  131)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  132)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  133)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  134)      !!         ! 00-06 (J.-M. Molines)
(  135)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  136)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  137)      !!------------------------------------------------------------------------------
(  138)      !! * Arguments
(  139)      INTEGER, INTENT( in ) ::   kt
(  140) 
(  141)      !! * Local declaration
(  142)      REAL(wp) ::   z05cx, ztau, zin
(  143) !!DB
(  144)      INTEGER :: jji, jjj, jjk
(  145)      REAL(wp) :: ts1
(  146) 
(  147)      !!------------------------------------------------------------------------------
(  148) 
(  149)      ! 1. First three time steps and more if lfbceast is .TRUE.
(  150)      !    In that case open boundary conditions are FIXED.
(  151)      ! --------------------------------------------------------
(  152) 
(  153)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast .OR. lk_dynspg_exp ) THEN
(  154) 
(  155)         ! 1.1 U zonal velocity
(  156)         ! --------------------
(  157)         DO ji = nie0, nie1
(  158)            DO jk = 1, jpkm1
(  159)               DO jj = 1, jpj
(  160) # 165
(  165)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uemsk(jj,jk)) + &
(  166)                                 uemsk(jj,jk)*ufoe(jj,jk)
(  167) # 168
(  168)               END DO
(  169)            END DO
(  170)         END DO
(  171) 
(  172) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  173) !!
(  174)         zin = 2. * rdt / (1.0*rday)
(  175) !        do jji = nie0, nie1   !!isolate processor 
(  176) !DBG -- skip 
(  177)         do jji = 1, 0 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 20

(  178)            do ji = nie0-4, nie0-1 !!5 cells - 1
(  179) !restore based on position -- but do not use
(  180)               ztau = zin * exp(float( -(nie0-ji)/4 )) 
(  181)               do jk = 1, jpkm1
(  182)                  do jj = 1, jpj
(  183)                     ua(ji,jj,jk) = ua(ji,jj,jk) + uemsk(jj,jk)*&
(  184)                          zin*(ufoe(jj,jk)-ua(ji,jj,jk)) * umask(ji,jj,jk)
(  185)                  enddo
(  186)               enddo
(  187)            enddo
(  188)         enddo
(  189) 
(  190) 
(  191) 
(  192)         ! 1.2 V meridional velocity
(  193)         ! -------------------------
(  194)         DO ji = nie0+1, nie1+1
(  195)            DO jk = 1, jpkm1
(  196)               DO jj = 1, jpj
(  197)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vemsk(jj,jk)) + &
(  198)                                 vfoe(jj,jk)*vemsk(jj,jk)
(  199)               END DO
(  200)            END DO
(  201)         END DO
(  202) 
(  203)      ELSE
(  204) 
(  205)      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
(  206)      ! -----------------------------------------------------
(  207) 
(  208)         ! 2.1. u-component of the velocity
(  209)         ! ---------------------------------
(  210)         !
(  211)         !          nibm2      nibm      nib
(  212)         !            |   nibm  |   nib   |///
(  213)         !            |    |    |    |    |///
(  214)         !  jj-line --f----v----f----v----f---
(  215)         !            |    |    |    |    |///
(  216)         !            |         |         |///
(  217)         !  jj-line   u    T    u    T    u///
(  218)         !            |         |         |///
(  219)         !            |    |    |    |    |///
(  220)         !          jpieob-2   jpieob-1   jpieob
(  221)         !                 |         |
(  222)         !              jpieob-1    jpieob
(  223)         !
(  224)         ! ... If free surface formulation:
(  225)         ! ... radiative conditions on the total part + relaxation toward climatology
(  226)         ! ... (jpjedp1, jpjefm1),jpieob
(  227)         DO ji = nie0, nie1
(  228)            DO jk = 1, jpkm1
(  229)               DO jj = 1, jpj
(  230)                  z05cx = u_cxebnd(jj,jk)
(  231)                  z05cx = z05cx / e1t(ji,jj)
(  232)                  z05cx = min( z05cx, 1. )
(  233)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  234)         !           > 0, outflow zin=1, ztau=rtaue
(  235)                  zin = sign( 1., z05cx )






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 21

(  236)                  zin = 0.5*( zin + abs(zin) )
(  237)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
(  238)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  239)                  z05cx = z05cx * zin
(  240)         ! ... update ua with radiative or climatological velocity
(  241)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uemsk(jj,jk) ) +          &
(  242)                                 uemsk(jj,jk) * (  ( 1. - z05cx - ztau )         &
(  243)                                 * uebnd(jj,jk,nib ,nitm) + 2.*z05cx               &
(  244)                                 * uebnd(jj,jk,nibm,nit ) + ztau * ufoe (jj,jk) )  &
(  245)                                 / (1. + z05cx)
(  246)               END DO
(  247)            END DO
(  248)         END DO
(  249) # 265
(  265)         ! 2.2 v-component of the velocity
(  266)         ! -------------------------------
(  267)         !
(  268)         !          nibm2       nibm     nib
(  269)         !            |   nibm  |   nib///|///
(  270)         !            |    |    |    |////|///
(  271)         !  jj-line --v----f----v----f----v---
(  272)         !            |    |    |    |////|///
(  273)         !            |    |    |    |////|///
(  274)         !            | jpieob-1 |  jpieob /|///
(  275)         !            |         |         |
(  276)         !         jpieob-1    jpieob     jpieob+1
(  277)         !
(  278)         ! ... radiative condition
(  279)         ! ... (jpjedp1, jpjefm1), jpieob+1
(  280)         DO ji = nie0+1, nie1+1
(  281)            DO jk = 1, jpkm1
(  282)               DO jj = 1, jpj
(  283)                  z05cx = v_cxebnd(jj,jk)
(  284)                  z05cx = z05cx / e1f(ji-1,jj)
(  285)                  z05cx = min( z05cx, 1. )
(  286)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  287)         !           > 0, outflow zin=1, ztau=rtaue
(  288)                  zin = sign( 1., z05cx )
(  289)                  zin = 0.5*( zin + abs(zin) )
(  290)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
(  291)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  292)                  z05cx = z05cx * zin
(  293)         ! ... update va with radiative or climatological velocity
(  294)                  va(ji,jj,jk) = va(ji,jj,jk) * (1. - vemsk(jj,jk) ) +          &
(  295)                                 vemsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  296)                                 * vebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  297)                                 * vebnd(jj,jk,nibm,nit ) + ztau * vfoe(jj,jk) )  &
(  298)                                 / (1. + z05cx)
(  299)               END DO
(  300)            END DO
(  301)         END DO
(  302) 
(  303)      END IF
(  304) 
(  305)   END SUBROUTINE obc_dyn_east
(  306) 
(  307) 
(  308)   SUBROUTINE obc_dyn_west ( kt )






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 22

(  309)      !!------------------------------------------------------------------------------
(  310)      !!                  ***  SUBROUTINE obc_dyn_west  ***
(  311)      !!
(  312)      !! ** Purpose :
(  313)      !!      Apply the radiation algorithm on west OBC velocities ua, va using the
(  314)      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
(  315)      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
(  316)      !!
(  317)      !!  History :
(  318)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  319)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  320)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  321)      !!         ! 00-06 (J.-M. Molines)
(  322)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  323)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  324)      !!------------------------------------------------------------------------------
(  325)      !! * Arguments
(  326)      INTEGER, INTENT( in ) ::   kt
(  327) 
(  328)      !! * Local declaration
(  329)      REAL(wp) ::   z05cx, ztau, zin
(  330) !!DB
(  331)      INTEGER :: jji, jjj, jjk
(  332) 
(  333)      !!------------------------------------------------------------------------------
(  334) 
(  335)      ! 1. First three time steps and more if lfbcwest is .TRUE.
(  336)      !    In that case open boundary conditions are FIXED.
(  337)      ! --------------------------------------------------------
(  338) 
(  339)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest .OR. lk_dynspg_exp ) THEN
(  340) 
(  341)         ! 1.1 U zonal velocity
(  342)         ! ---------------------
(  343)         DO ji = niw0, niw1
(  344)            DO jk = 1, jpkm1
(  345)               DO jj = 1, jpj
(  346) # 351
(  351)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-uwmsk(jj,jk)) + &
(  352)                                 uwmsk(jj,jk)*ufow(jj,jk)
(  353) # 354
(  354)               END DO
(  355)            END DO
(  356)         END DO
(  357) 
(  358) 
(  359) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  360) !!
(  361)         zin = 2. * rdt / (1.0*rday)
(  362) !        do jji = niw0, niw1   !!isolate processor 
(  363)        do jji = 1, 0    !!isolate processor 
(  364)            do ji = niw0+1, niw0+4 !!5 cells - 1
(  365) !!restore based on position -- but do not use
(  366)               ztau = zin * exp(float( -(ji-niw0)/4 )) 
(  367)               do jk = 1, jpkm1
(  368)                  do jj = 1, jpj
(  369)                     ua(ji,jj,jk) = ua(ji,jj,jk) + uwmsk(jj,jk)*&
(  370)                          zin*(ufow(jj,jk)-ua(ji,jj,jk)) * umask(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 23

(  371)                  enddo
(  372)               enddo
(  373)            enddo
(  374)         enddo
(  375) 
(  376) 
(  377) 
(  378)         ! 1.2 V meridional velocity
(  379)         ! -------------------------
(  380)         DO ji = niw0, niw1
(  381)            DO jk = 1, jpkm1
(  382)               DO jj = 1, jpj
(  383)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vwmsk(jj,jk)) + &
(  384)                                 vfow(jj,jk)*vwmsk(jj,jk)
(  385)               END DO
(  386)            END DO
(  387)         END DO
(  388) 
(  389)      ELSE
(  390) 
(  391)      ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
(  392)      ! -----------------------------------------------------
(  393) 
(  394)         ! 2.1. u-component of the velocity
(  395)         ! ---------------------------------
(  396)         !
(  397)         !        nib       nibm     nibm2
(  398)         !      ///|   nib   |   nibm  |
(  399)         !      ///|    |    |    |    |
(  400)         !      ---f----v----f----v----f-- jj-line
(  401)         !      ///|    |    |    |    |
(  402)         !      ///|         |         |
(  403)         !      ///u    T    u    T    u   jj-line
(  404)         !      ///|         |         |
(  405)         !      ///|    |    |    |    |
(  406)         !       jpiwob    jpiwob+1    jpiwob+2
(  407)         !              |         |
(  408)         !            jpiwob+1    jpiwob+2
(  409)         !
(  410)         ! ... If free surface formulation:
(  411)         ! ... radiative conditions on the total part + relaxation toward climatology
(  412)         ! ... (jpjwdp1, jpjwfm1), jpiwob
(  413)         DO ji = niw0, niw1
(  414)            DO jk = 1, jpkm1
(  415)               DO jj = 1, jpj
(  416)                  z05cx = u_cxwbnd(jj,jk)
(  417)                  z05cx = z05cx / e1t(ji+1,jj)
(  418)                  z05cx = max( z05cx, -1. )
(  419)         ! ... z05c  > 0, inflow  zin=0, ztau=1
(  420)         !          =< 0, outflow zin=1, ztau=rtauw
(  421)                  zin = sign( 1., -1. * z05cx )
(  422)                  zin = 0.5*( zin + abs(zin) )
(  423)                  ztau = (1.-zin )* rtauwin + zin * rtauw
(  424)                  z05cx = z05cx * zin
(  425)         ! ... update un with radiative or climatological velocity
(  426)                  ua(ji,jj,jk) = ua(ji,jj,jk) * ( 1. - uwmsk(jj,jk) ) +          &
(  427)                                 uwmsk(jj,jk) * ( ( 1. + z05cx - ztau )          &
(  428)                                 * uwbnd(jj,jk,nib ,nitm) - 2.*z05cx               &






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 24

(  429)                                 * uwbnd(jj,jk,nibm,nit ) + ztau  * ufow (jj,jk) ) &
(  430)                                 / (1. - z05cx)
(  431)               END DO
(  432)            END DO
(  433)         END DO
(  434) 
(  435) # 451
(  451)         ! 2.2 v-component of the velocity
(  452)         ! -------------------------------
(  453)         !
(  454)         !    nib       nibm     nibm2
(  455)         !  ///|///nib   |   nibm  |  nibm2
(  456)         !  ///|////|    |    |    |    |    |
(  457)         !  ---v----f----v----f----v----f----v-- jj-line
(  458)         !  ///|////|    |    |    |    |    |
(  459)         !  ///|////|    |    |    |    |    |
(  460)         ! jpiwob     jpiwob+1    jpiwob+2
(  461)         !          |         |         |
(  462)         !        jpiwob   jpiwob+1   jpiwob+2
(  463)         !
(  464)         ! ... radiative condition plus Raymond-Kuo
(  465)         ! ... (jpjwdp1, jpjwfm1),jpiwob
(  466)         DO ji = niw0, niw1
(  467)            DO jk = 1, jpkm1
(  468)               DO jj = 1, jpj
(  469)                  z05cx = v_cxwbnd(jj,jk)
(  470)                  z05cx = z05cx / e1f(ji,jj)
(  471)                  z05cx = max( z05cx, -1. )
(  472)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  473)         !          =< 0, outflow zin=1, ztau=rtauw
(  474)                  zin = sign( 1., -1. * z05cx )
(  475)                  zin = 0.5*( zin + abs(zin) )
(  476)                  ztau = (1.-zin )*rtauwin + zin * rtauw
(  477)                  z05cx = z05cx * zin
(  478)         ! ... update va with radiative or climatological velocity
(  479)                  va(ji,jj,jk) = va(ji,jj,jk) * (1. - vwmsk(jj,jk) ) +          &
(  480)                                 vwmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  481)                                 * vwbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  482)                                 * vwbnd(jj,jk,nibm,nit ) + ztau * vfow (jj,jk) ) &
(  483)                                 / (1. - z05cx)
(  484)                END DO
(  485)             END DO
(  486)         END DO
(  487) 
(  488)      END IF
(  489) 
(  490)   END SUBROUTINE obc_dyn_west
(  491) 
(  492)   SUBROUTINE obc_dyn_north ( kt )
(  493)      !!------------------------------------------------------------------------------
(  494)      !!                     SUBROUTINE obc_dyn_north
(  495)      !!                    *************************
(  496)      !! ** Purpose :
(  497)      !!      Apply the radiation algorithm on north OBC velocities ua, va using the
(  498)      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
(  499)      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
(  500)      !!
(  501)      !!  History :






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 25

(  502)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  503)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  504)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  505)      !!         ! 00-06 (J.-M. Molines)
(  506)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  507)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  508)      !!------------------------------------------------------------------------------
(  509)      !! * Arguments
(  510)      INTEGER, INTENT( in ) ::   kt
(  511) 
(  512)      !! * Local declaration
(  513)      REAL(wp) ::   z05cx, ztau, zin
(  514)      !!------------------------------------------------------------------------------
(  515) 
(  516)      ! 1. First three time steps and more if lfbcnorth is .TRUE.
(  517)      !    In that case open boundary conditions are FIXED.
(  518)      ! ---------------------------------------------------------
(  519) 
(  520)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth  .OR. lk_dynspg_exp ) THEN
(  521) 
(  522)         ! 1.1 U zonal velocity
(  523)         ! --------------------
(  524)         DO jj = njn0+1, njn1+1
(  525)            DO jk = 1, jpkm1
(  526)               DO ji = 1, jpi
(  527)                  ua(ji,jj,jk)= ua(ji,jj,jk) * (1.-unmsk(ji,jk)) + &
(  528)                                ufon(ji,jk)*unmsk(ji,jk)
(  529)               END DO
(  530)            END DO
(  531)         END DO
(  532) 
(  533)         ! 1.2 V meridional velocity
(  534)         ! -------------------------
(  535)         DO jj = njn0, njn1
(  536)            DO jk = 1, jpkm1
(  537)               DO ji = 1, jpi
(  538) # 543
(  543)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vnmsk(ji,jk)) + &
(  544)                                vfon(ji,jk)*vnmsk(ji,jk)
(  545) # 546
(  546)               END DO
(  547)            END DO
(  548)         END DO
(  549) 
(  550)      ELSE
(  551) 
(  552)      ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
(  553)      ! ------------------------------------------------------
(  554) 
(  555)         ! 2.1. u-component of the velocity
(  556)         ! --------------------------------
(  557)         !
(  558)         !            ji-row
(  559)         !              |
(  560)         !       nib ///u//////  jpjnob + 1
(  561)         !         /////|//////
(  562)         !     nib -----f-----   jpjnob
(  563)         !              |






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 26

(  564)         !      nibm--  u   ---- jpjnob
(  565)         !              |
(  566)         !    nibm -----f-----   jpjnob-1
(  567)         !              |
(  568)         !     nibm2--  u   ---- jpjnob-1
(  569)         !              |
(  570)         !   nibm2 -----f-----   jpjnob-2
(  571)         !              |
(  572)         !
(  573)         ! ... radiative condition
(  574)         ! ... jpjnob+1,(jpindp1, jpinfm1)
(  575)         DO jj = njn0+1, njn1+1
(  576)            DO jk = 1, jpkm1
(  577)               DO ji = 1, jpi
(  578)                  z05cx= u_cynbnd(ji,jk)
(  579)                  z05cx = z05cx / e2f(ji, jj-1)
(  580)                  z05cx = min( z05cx, 1. )
(  581)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  582)         !           > 0, outflow zin=1, ztau=rtaun
(  583)                  zin = sign( 1., z05cx )
(  584)                  zin = 0.5*( zin + abs(zin) )
(  585)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  586)                  ztau = (1.-zin ) * rtaunin  + zin * rtaun
(  587)         ! ... for u, when inflow, ufon is prescribed
(  588)                  z05cx = z05cx * zin
(  589)         ! ... update un with radiative or climatological velocity
(  590)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-unmsk(ji,jk)) +             &
(  591)                                 unmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  592)                                 * unbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  593)                                 * unbnd(ji,jk,nibm,nit ) + ztau * ufon (ji,jk) ) &
(  594)                                 / (1. + z05cx)
(  595)               END DO
(  596)            END DO
(  597)         END DO
(  598) 
(  599)         ! 2.2 v-component of the velocity
(  600)         ! -------------------------------
(  601)         !
(  602)         !                ji-row    ji-row
(  603)         !              |         |
(  604)         !         /////|/////////////////
(  605)         !    nib  -----f----v----f----  jpjnob
(  606)         !              |         |
(  607)         !      nib  -  u -- T -- u ---- jpjnob
(  608)         !              |         |
(  609)         !   nibm  -----f----v----f----  jpjnob-1
(  610)         !              |         |
(  611)         !     nibm --  u -- T -- u ---  jpjnob-1
(  612)         !              |         |
(  613)         !   nibm2 -----f----v----f----  jpjnob-2
(  614)         !              |         |
(  615)         !
(  616)         ! ... If rigidlid formulation:
(  617)         ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  618)         ! ... If free surface formulation:
(  619)         ! ... radiative conditions on the total part + relaxation toward climatology
(  620)         ! ... jpjnob,(jpindp1, jpinfm1)
(  621)         DO jj = njn0, njn1






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 27

(  622)            DO jk = 1, jpkm1
(  623)               DO ji = 1, jpi
(  624)         ! ... 2* gradj(v) (T-point i=nibm, time mean)
(  625)                  z05cx = v_cynbnd(ji,jk)
(  626)                  z05cx = z05cx / e2t(ji,jj)
(  627)                  z05cx = min( z05cx, 1. )
(  628)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  629)         !           > 0, outflow zin=1, ztau=rtaun
(  630)                  zin = sign( 1., z05cx )
(  631)                  zin = 0.5*( zin + abs(zin) )
(  632)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  633)                  ztau = (1.-zin ) * rtaunin + zin * rtaun
(  634)                  z05cx = z05cx * zin
(  635)         ! ... update va with radiative or climatological velocity
(  636)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vnmsk(ji,jk)) +             &
(  637)                                 vnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  638)                                 * vnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  639)                                 * vnbnd(ji,jk,nibm,nit ) + ztau * vfon (ji,jk) ) &
(  640)                                 / (1. + z05cx)
(  641)               END DO
(  642)            END DO
(  643)         END DO
(  644) # 660
(  660)      END IF
(  661) 
(  662)   END SUBROUTINE obc_dyn_north
(  663) 
(  664)   SUBROUTINE obc_dyn_south ( kt )
(  665)      !!------------------------------------------------------------------------------
(  666)      !!                     SUBROUTINE obc_dyn_south
(  667)      !!                    *************************
(  668)      !! ** Purpose :
(  669)      !!      Apply the radiation algorithm on south OBC velocities ua, va using the
(  670)      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
(  671)      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
(  672)      !!
(  673)      !!  History :
(  674)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  675)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  676)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  677)      !!         ! 00-06 (J.-M. Molines)
(  678)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  679)      !!    9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
(  680)      !!------------------------------------------------------------------------------
(  681)      !! * Arguments
(  682)      INTEGER, INTENT( in ) ::   kt
(  683) 
(  684)      !! * Local declaration
(  685)      REAL(wp) ::   z05cx, ztau, zin
(  686) !!DB
(  687)      INTEGER :: jji, jjj, jjk
(  688) 
(  689) 
(  690)      !!------------------------------------------------------------------------------
(  691)      !!  OPA 8.5, LODYC-IPSL (2002)
(  692)      !!------------------------------------------------------------------------------
(  693) 
(  694)      ! 1. First three time steps and more if lfbcsouth is .TRUE.






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 28

(  695)      !    In that case open boundary conditions are FIXED.
(  696)      ! ---------------------------------------------------------
(  697) 
(  698)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth  .OR. lk_dynspg_exp ) THEN
(  699) 
(  700)         ! 1.1 U zonal velocity
(  701)         ! --------------------
(  702)         DO jj = njs0, njs1
(  703)            DO jk = 1, jpkm1
(  704)               DO ji = 1, jpi
(  705)                  ua(ji,jj,jk)= ua(ji,jj,jk) * (1.-usmsk(ji,jk)) + &
(  706)                                usmsk(ji,jk) * ufos(ji,jk)
(  707)               END DO
(  708)            END DO
(  709)         END DO
(  710) 
(  711) 
(  712) 
(  713)         ! 1.2 V meridional velocity
(  714)         ! -------------------------
(  715)         DO jj = njs0, njs1
(  716)            DO jk = 1, jpkm1
(  717)               DO ji = 1, jpi
(  718) # 723
(  723)                  va(ji,jj,jk)= va(ji,jj,jk) * (1.-vsmsk(ji,jk)) + &
(  724)                                vsmsk(ji,jk) * vfos(ji,jk)
(  725) # 726
(  726)               END DO
(  727)            END DO
(  728)         END DO
(  729) 
(  730) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  731) !!
(  732)         zin = 2. * rdt / (1.0*rday)
(  733) !        do jji = njs0, njs1   !!isolate processor 
(  734)         do jji = 1, 0   !!isolate processor 
(  735)            do jj = njs0+1, njs0+4 !!5 cells - 1
(  736) !!restore based on position -- but do not use
(  737)               ztau = zin * exp(float( -(jj-njs0)/4 )) 
(  738)               do jk = 1, jpkm1
(  739)                  do ji = 1, jpi
(  740)                     va(ji,jj,jk) = va(ji,jj,jk) + vsmsk(ji,jk)*&
(  741)                          zin*(vfos(ji,jk)-va(ji,jj,jk)) * vmask(ji,jj,jk)
(  742)                  enddo
(  743)               enddo
(  744)            enddo
(  745)         enddo
(  746) 
(  747) 
(  748) 
(  749) 
(  750) 
(  751) 
(  752)      ELSE
(  753) 
(  754)      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
(  755)      ! ------------------------------------------------------
(  756) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 29

(  757)         ! 2.1. u-component of the velocity
(  758)         ! --------------------------------
(  759)         !
(  760)         !            ji-row
(  761)         !              |
(  762)         !   nibm2 -----f-----   jpjsob +2
(  763)         !              |
(  764)         !    nibm2 --  u   ---- jpjsob +2
(  765)         !              |
(  766)         !    nibm -----f-----   jpjsob +1
(  767)         !              |
(  768)         !    nibm  --  u   ---- jpjsob +1
(  769)         !              |
(  770)         !    nib  -----f-----   jpjsob
(  771)         !         /////|//////
(  772)         !    nib   ////u/////   jpjsob
(  773)         !
(  774)         ! ... radiative condition plus Raymond-Kuo
(  775)         ! ... jpjsob,(jpisdp1, jpisfm1)
(  776)         DO jj = njs0, njs1
(  777)            DO jk = 1, jpkm1
(  778)               DO ji = 1, jpi
(  779)                  z05cx= u_cysbnd(ji,jk)
(  780)                  z05cx = z05cx / e2f(ji, jj)
(  781)                  z05cx = max( z05cx, -1. )
(  782)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  783)         !          =< 0, outflow zin=1, ztau=rtaus
(  784)                  zin = sign( 1., -1. * z05cx )
(  785)                  zin = 0.5*( zin + abs(zin) )
(  786)                  ztau = (1.-zin ) * rtausin + zin * rtaus
(  787)                  z05cx = z05cx * zin
(  788)         ! ... update ua with radiative or climatological velocity
(  789)                  ua(ji,jj,jk) = ua(ji,jj,jk) * (1.-usmsk(ji,jk)) +              &
(  790)                                 usmsk(ji,jk) * (  ( 1. + z05cx - ztau )         &
(  791)                                 * usbnd(ji,jk,nib ,nitm) - 2.*z05cx               &
(  792)                                 * usbnd(ji,jk,nibm,nit ) + ztau * ufos (ji,jk) )  &
(  793)                                 / (1. - z05cx)
(  794)               END DO
(  795)            END DO
(  796)         END DO
(  797) 
(  798)         ! 2.2 v-component of the velocity
(  799)         ! -------------------------------
(  800)         !
(  801)         !                ji-row    ji-row
(  802)         !              |         |
(  803)         !  nibm2  -----f----v----f----  jpjsob+2
(  804)         !              |         |
(  805)         !    nibm   -  u -- T -- u ---- jpjsob+2
(  806)         !              |         |
(  807)         !   nibm  -----f----v----f----  jpjsob+1
(  808)         !              |         |
(  809)         !   nib    --  u -- T -- u ---  jpjsob+1
(  810)         !              |         |
(  811)         !   nib   -----f----v----f----  jpjsob
(  812)         !         /////////////////////
(  813)         !
(  814)         ! ... If rigidlid formulation:






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 30

(  815)         ! ... radiative conditions on the baroclinic part only + relaxation toward climatology
(  816)         ! ... If free surface formulation:
(  817)         ! ... radiative conditions on the total part + relaxation toward climatology
(  818)         ! ... jpjsob,(jpisdp1,jpisfm1)
(  819)         DO jj = njs0, njs1
(  820)            DO jk = 1, jpkm1
(  821)               DO ji = 1, jpi
(  822)                  z05cx = v_cysbnd(ji,jk)
(  823)                  z05cx = z05cx / e2t(ji,jj+1)
(  824)                  z05cx = max( z05cx, -1. )
(  825)         ! ... z05c > 0, inflow  zin=0, ztau=1
(  826)         !         =< 0, outflow zin=1, ztau=rtaus
(  827)                  zin = sign( 1., -1. * z05cx )
(  828)                  zin = 0.5*( zin + abs(zin) )
(  829)                  ztau = (1.-zin )*rtausin + zin * rtaus
(  830)                  z05cx = z05cx * zin
(  831)         ! ... update va with radiative or climatological velocity
(  832)                  va(ji,jj,jk) = va(ji,jj,jk) * (1.-vsmsk(ji,jk)) +             &
(  833)                                 vsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
(  834)                                 * vsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
(  835)                                 * vsbnd(ji,jk,nibm,nit ) + ztau * vfos (ji,jk) ) &
(  836)                                 / (1. - z05cx)
(  837)               END DO
(  838)            END DO
(  839)         END DO
(  840) # 856
(  856)      END IF
(  857) 
(  858)   END SUBROUTINE obc_dyn_south
(  859) 
(  860) # 871
(  871) END MODULE obcdyn
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (obcdyn.F90: 180)
