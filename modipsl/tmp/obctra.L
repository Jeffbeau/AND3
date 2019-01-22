


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

Filename: obctra.F90

(    1) MODULE obctra
(    2)   !!=================================================================================
(    3)   !!                       ***  MODULE  obctra  ***
(    4)   !! Ocean tracers:   Radiation of tracers on each open boundary
(    5)   !!=================================================================================
(    6) #if defined key_obc
(    7)   !!---------------------------------------------------------------------------------
(    8)   !!   'key_obc'      :                                      Open Boundary Conditions
(    9)   !!---------------------------------------------------------------------------------
(   10)   !!   obc_tra        : call the subroutine for each open boundary
(   11)   !!   obc_tra_east   : radiation of the east open boundary tracers
(   12)   !!   obc_tra_west   : radiation of the west open boundary tracers
(   13)   !!   obc_tra_north  : radiation of the north open boundary tracers
(   14)   !!   obc_tra_south  : radiation of the south open boundary tracers
(   15)   !!----------------------------------------------------------------------------------
(   16)   !! * Modules used
(   17)   USE oce             ! ocean dynamics and tracers variables
(   18)   USE dom_oce         ! ocean space and time domain variables
(   19)   USE phycst          ! physical constants
(   20)   USE obc_oce         ! ocean open boundary conditions
(   21)   USE lib_mpp         ! ???
(   22)   USE lbclnk          ! ???
(   23)   USE in_out_manager  ! I/O manager
(   24) 
(   25)   IMPLICIT NONE
(   26)   PRIVATE
(   27) 
(   28)   !! * Accessibility
(   29)   PUBLIC obc_tra     ! routine called in tranxt.F90
(   30) 
(   31)   !! * Module variables
(   32)   INTEGER ::      & ! ... boundary space indices
(   33)      nib   = 1,   & ! nib   = boundary point
(   34)      nibm  = 2,   & ! nibm  = 1st interior point
(   35)      nibm2 = 3,   & ! nibm2 = 2nd interior point
(   36)                     ! ... boundary time indices
(   37)      nit   = 1,   & ! nit    = now
(   38)      nitm  = 2,   & ! nitm   = before
(   39)      nitm2 = 3      ! nitm2  = before-before
(   40) 
(   41)   REAL(wp) ::     &
(   42)      rtaue  , rtauw  , rtaun  , rtaus  ,  &  ! Boundary restoring coefficient
(   43)      rtauein, rtauwin, rtaunin, rtausin      ! Boundary restoring coefficient for inflow
(   44) 
(   45)   !! * Substitutions
(   46) #  include "obc_vectopt_loop_substitute.h90"






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 2

(   47)   !!---------------------------------------------------------------------------------
(   48)   !!   OPA 9.0 , LOCEAN-IPSL (2005)
(   49)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obctra.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $
(   50)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   51)   !!---------------------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55)   SUBROUTINE obc_tra( kt )
(   56)      !!-------------------------------------------------------------------------------
(   57)      !!                 ***  SUBROUTINE obc_tra  ***
(   58)      !!
(   59)      !! ** Purpose :   Compute tracer fields (t,s) along the open boundaries.
(   60)      !!      This routine is called by the tranxt.F routine and updates ta,sa
(   61)      !!      which are the actual temperature and salinity fields.
(   62)      !!        The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   63)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   64)      !!      open one (must be done in the param_obc.h90 file).
(   65)      !!
(   66)      !! Reference :
(   67)      !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   68)      !!
(   69)      !!  History :
(   70)      !!        !  95-03 (J.-M. Molines) Original, SPEM
(   71)      !!        !  97-07 (G. Madec, J.-M. Molines) addition
(   72)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) F90
(   73)      !!----------------------------------------------------------------------
(   74)      !! * Arguments
(   75)      INTEGER, INTENT( in ) ::   kt
(   76)      !!----------------------------------------------------------------------
(   77) 
(   78)      ! 0. Local constant initialization
(   79) 
(   80)      IF( kt == nit000 .OR. ln_rstart) THEN
(   81)         ! ... Boundary restoring coefficient
(   82)         rtaue = 2. * rdt / rdpeob
(   83)         rtauw = 2. * rdt / rdpwob
(   84)         rtaun = 2. * rdt / rdpnob
(   85)         rtaus = 2. * rdt / rdpsob
(   86)         ! ... Boundary restoring coefficient for inflow ( all boundaries)
(   87)         rtauein = 2. * rdt / rdpein
(   88)         rtauwin = 2. * rdt / rdpwin
(   89)         rtaunin = 2. * rdt / rdpnin
(   90)         rtausin = 2. * rdt / rdpsin
(   91)      END IF
(   92) 
(   93)      IF( lp_obc_east  )   CALL obc_tra_east ( kt )    ! East open boundary
(   94) 
(   95)      IF( lp_obc_west  )   CALL obc_tra_west ( kt )    ! West open boundary
(   96) 
(   97)      IF( lp_obc_north )   CALL obc_tra_north( kt )    ! North open boundary
(   98) 
(   99)      IF( lp_obc_south )   CALL obc_tra_south( kt )    ! South open boundary
(  100) 
(  101)      IF( lk_mpp ) THEN                  !!bug ???
(  102)         IF( kt >= nit000+3 .AND. ln_rstart ) THEN
(  103)            CALL lbc_lnk( tb, 'T', 1. )
(  104)            CALL lbc_lnk( sb, 'T', 1. )






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 3

(  105)         END IF
(  106)         CALL lbc_lnk( ta, 'T', 1. )
(  107)         CALL lbc_lnk( sa, 'T', 1. )
(  108)      ENDIF
(  109) 
(  110)   END SUBROUTINE obc_tra
(  111) 
(  112) 
(  113)   SUBROUTINE obc_tra_east ( kt )
(  114)      !!------------------------------------------------------------------------------
(  115)      !!                ***  SUBROUTINE obc_tra_east  ***
(  116)      !!
(  117)      !! ** Purpose :
(  118)      !!      Apply the radiation algorithm on east OBC tracers ta, sa using the
(  119)      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
(  120)      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
(  121)      !!
(  122)      !!  History :
(  123)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  124)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  125)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  126)      !!         ! 00-06 (J.-M. Molines)
(  127)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  128)      !!------------------------------------------------------------------------------
(  129)      !! * Arguments
(  130)      INTEGER, INTENT( in ) ::   kt
(  131) 
(  132)      !! * Local declaration
(  133)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  134)      REAL(wp) ::   z05cx, ztau, zin
(  135) !!DB
(  136)      INTEGER :: jji, jjj, jjk
(  137) 
(  138)      !!------------------------------------------------------------------------------
(  139) 
(  140)      ! 1. First three time steps and more if lfbceast is .TRUE.
(  141)      !    In that case open boundary conditions are FIXED.
(  142)      ! --------------------------------------------------------
(  143) 
(  144)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast ) THEN
(  145)         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
(  146)            DO jk = 1, jpkm1
(  147)               DO jj = 1, jpj
(  148)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-temsk(jj,jk)) + &
(  149)                                temsk(jj,jk) * tfoe(jj,jk,1)
(  150)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-temsk(jj,jk)) + &
(  151)                                temsk(jj,jk) * sfoe(jj,jk,1)
(  152)               END DO
(  153)            END DO
(  154)         END DO
(  155) 
(  156) !!DB DBG TS sponge; 5 cells; timescale = 1day(?)
(  157) !!
(  158)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  159) !        zin = 2. * rdt / (20.0*rday)
(  160) !        do jji = nie0, nie1   !!isolate processor 
(  161) !DBG -- skip 
(  162)         do jji = 1, 0 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 4

(  163)            do ji = nie0-4, nie0-1 !!5 cells - 1
(  164) !restore based on position -- but do not use
(  165)               ztau = zin * exp(float( -(nie0-ji)/4 )) 
(  166)               do jk = 1, jpkm1
(  167)                  do jj = 1, jpj
(  168)                     ta(ji,jj,jk) = ta(ji,jj,jk) + temsk(jj,jk)*&
(  169)                          zin*(tfoe(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  170)                     sa(ji,jj,jk) = sa(ji,jj,jk) + temsk(jj,jk)*&
(  171)                          zin*(sfoe(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  172)                  enddo
(  173)               enddo
(  174)            enddo
(  175)         enddo
(  176) 
(  177) 
(  178) 
(  179) 
(  180) 
(  181) 
(  182) 
(  183)      ELSE
(  184) 
(  185)      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
(  186)      ! -----------------------------------------------------
(  187) 
(  188)         ! Temperature and salinity radiation
(  189)         ! ----------------------------------
(  190)         !
(  191)         !            nibm2      nibm      nib
(  192)         !              |   nibm  |   nib///|///
(  193)         !              |    |    |    |////|///
(  194)         !  jj   line --v----f----v----f----v---
(  195)         !              |    |    |    |////|///
(  196)         !                   |         |///   //
(  197)         !  jj   line   T    u    T    u/// T //
(  198)         !                   |         |///   //
(  199)         !              |    |    |    |////|///
(  200)         !  jj-1 line --v----f----v----f----v---
(  201)         !              |    |    |    |////|///
(  202)         !                jpieob-1    jpieob / ///
(  203)         !              |         |         |
(  204)         !           jpieob-1    jpieob     jpieob+1
(  205)         !
(  206)         ! ... radiative conditions + relaxation toward a climatology
(  207)         !     the phase velocity is taken as the phase velocity of the tangen-
(  208)         !     tial velocity (here vn), which have been saved in (u_cxebnd,v_cxebnd)
(  209)         ! ... (jpjedp1, jpjefm1), jpieob+1
(  210)         DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
(  211)            DO jk = 1, jpkm1
(  212)               DO jj = 2, jpjm1
(  213)         ! ... i-phase speed ratio (from averaged of v_cxebnd)
(  214)                  z05cx = ( 0.5 * ( v_cxebnd(jj,jk) + v_cxebnd(jj-1,jk) ) ) / e1t(ji-1,jj)
(  215)                  z05cx = min( z05cx, 1. )
(  216)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  217)         !           > 0, outflow zin=1, ztau=rtaue
(  218)                  zin = sign( 1., z05cx )
(  219)                  zin = 0.5*( zin + abs(zin) )
(  220)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 5

(  221)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  222)                  z05cx = z05cx * zin
(  223)         ! ... update ( ta, sa ) with radiative or climatological (t, s)
(  224)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
(  225)                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  226)                                 * tebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  227) !byoung                                 * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk) ) &
(  228)                                 * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk,1) ) &
(  229)                                 / (1. + z05cx)
(  230)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
(  231)                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  232)                                 * sebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  233) !byoung                                * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk) ) &
(  234)                                 * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk,1) ) &
(  235)                                 / (1. + z05cx)
(  236)               END DO
(  237)            END DO
(  238)         END DO
(  239) 
(  240)      END IF
(  241) 
(  242)   END SUBROUTINE obc_tra_east
(  243) 
(  244) 
(  245)   SUBROUTINE obc_tra_west ( kt )
(  246)      !!------------------------------------------------------------------------------
(  247)      !!                 ***  SUBROUTINE obc_tra_west  ***
(  248)      !!
(  249)      !! ** Purpose :
(  250)      !!      Apply the radiation algorithm on west OBC tracers ta, sa using the
(  251)      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
(  252)      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
(  253)      !!
(  254)      !!  History :
(  255)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  256)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  257)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  258)      !!         ! 00-06 (J.-M. Molines)
(  259)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  260)      !!------------------------------------------------------------------------------
(  261)      !! * Arguments
(  262)      INTEGER, INTENT( in ) ::   kt
(  263) 
(  264)      !! * Local declaration
(  265)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  266)      REAL(wp) ::   z05cx, ztau, zin
(  267) !!DB
(  268)      INTEGER :: jji, jjj, jjk
(  269) 
(  270)      !!------------------------------------------------------------------------------
(  271) 
(  272)      ! 1. First three time steps and more if lfbcwest is .TRUE.
(  273)      !    In that case open boundary conditions are FIXED.
(  274)      ! --------------------------------------------------------
(  275) 
(  276)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest ) THEN
(  277) 
(  278)         DO ji = fs_niw0, fs_niw1 ! Vector opt.






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 6

(  279)            DO jk = 1, jpkm1
(  280)               DO jj = 1, jpj
(  281)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
(  282)                                twmsk(jj,jk) * tfow(jj,jk,1)
(  283)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
(  284)                                twmsk(jj,jk) * sfow(jj,jk,1)
(  285)               END DO
(  286)            END DO
(  287)         END DO
(  288) !!DB DBG TS sponge; 5 cells; timescale = 1day(?)
(  289) !!
(  290)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  291) !        zin = 2. * rdt / (20.0*rday)
(  292) !        do jji = niw0, niw1   !!isolate processor 
(  293)         do jji = 1, 0    !!isolate processor 
(  294)            do ji = niw0+1, niw0+4 !!5 cells - 1
(  295) !!restore based on position -- but do not use
(  296)               ztau = zin * exp(float( -(ji-niw0)/4 )) 
(  297)               do jk = 1, jpkm1
(  298)                  do jj = 1, jpj
(  299)                     ta(ji,jj,jk) = ta(ji,jj,jk) + twmsk(jj,jk)*&
(  300)                          zin*(tfow(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  301)                     sa(ji,jj,jk) = sa(ji,jj,jk) + twmsk(jj,jk)*&
(  302)                          zin*(sfow(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  303)                  enddo
(  304)               enddo
(  305)            enddo
(  306)         enddo
(  307) 
(  308)      ELSE
(  309) 
(  310)      ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
(  311)      ! -----------------------------------------------------
(  312) 
(  313)         ! Temperature and salinity radiation
(  314)         ! ----------------------------------
(  315)         !
(  316)         !          nib       nibm     nibm2
(  317)         !     nib///|   nibm  |  nibm2  |
(  318)         !   ///|////|    |    |    |    |
(  319)         !   ---v----f----v----f----v----f-- jj   line
(  320)         !   ///|////|    |    |    |    |
(  321)         !   //   ///|         |         |
(  322)         !   // T ///u    T    u    T    u   jj   line
(  323)         !   //   ///|         |         |
(  324)         !   ///|////|    |    |    |    |
(  325)         !   ---v----f----v----f----v----f-- jj-1 line
(  326)         !   ///|////|    |    |    |    |
(  327)         !         jpiwob    jpiwob+1    jpiwob+2
(  328)         !      |         |         |
(  329)         !    jpiwob    jpiwob+1   jpiwob+2
(  330)         !
(  331)         ! ... radiative conditions + relaxation toward a climatology
(  332)         ! ... the phase velocity is taken as the phase velocity of the tangen-
(  333)         ! ... tial velocity (here vn), which have been saved in (v_cxwbnd)
(  334)         DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  335)            DO jk = 1, jpkm1
(  336)               DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 7

(  337)         ! ... i-phase speed ratio (from averaged of v_cxwbnd)
(  338)                  z05cx = (  0.5 * ( v_cxwbnd(jj,jk) + v_cxwbnd(jj-1,jk) ) ) / e1t(ji+1,jj)
(  339)                  z05cx = max( z05cx, -1. )
(  340)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  341)         !           < 0, outflow zin=1, ztau=rtauw
(  342)                  zin = sign( 1., -1.* z05cx )
(  343)                  zin = 0.5*( zin + abs(zin) )
(  344)                  ztau = (1.-zin )*rtauwin + zin * rtauw
(  345)                  z05cx = z05cx * zin
(  346)         ! ... update (ta,sa) with radiative or climatological (t, s)
(  347)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
(  348)                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  349)                                 * twbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  350) !byoung                                * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk) ) &
(  351)                                 * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk,1) ) &
(  352)                                 / (1. - z05cx)
(  353)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
(  354)                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  355)                                 * swbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  356) !byoung                                * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk) ) &
(  357)                                 * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk,1) ) &
(  358)                                 / (1. - z05cx)
(  359)               END DO
(  360)            END DO
(  361)         END DO
(  362) 
(  363)      END IF
(  364) 
(  365)   END SUBROUTINE obc_tra_west
(  366) 
(  367) 
(  368)   SUBROUTINE obc_tra_north ( kt )
(  369)      !!------------------------------------------------------------------------------
(  370)      !!                 ***  SUBROUTINE obc_tra_north  ***
(  371)      !!
(  372)      !! ** Purpose :
(  373)      !!      Apply the radiation algorithm on north OBC tracers ta, sa using the
(  374)      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
(  375)      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
(  376)      !!
(  377)      !!  History :
(  378)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  379)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  380)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  381)      !!         ! 00-06 (J.-M. Molines)
(  382)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  383)      !!------------------------------------------------------------------------------
(  384)      !! * Arguments
(  385)      INTEGER, INTENT( in ) ::   kt
(  386) 
(  387)      !! * Local declaration
(  388) !sujie      INTEGER ::   ji, jj, jk      ! dummy loop indices
(  389)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  390)      REAL(wp) ::   z05cx, ztau, zin
(  391)      !!------------------------------------------------------------------------------
(  392) 
(  393)      ! 1. First three time steps and more if lfbcnorth is .TRUE.
(  394)      !    In that case open boundary conditions are FIXED.






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 8

(  395)      ! --------------------------------------------------------
(  396) 
(  397)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth ) THEN
(  398) 
(  399) !sujie
(  400)         DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
(  401)            DO jk = 1, jpkm1
(  402)               DO ji = 1, jpi
(  403) !                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  404) !                                tnmsk(ji,jk) * tfon(ji,jk)
(  405) !                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  406) !                                tnmsk(ji,jk) * sfon(ji,jk)
(  407)              if(sfon(ji,jk,1).gt.0.) then
(  408)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  409)                                tnmsk(ji,jk) * tfon(ji,jk,1)
(  410)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  411)                                tnmsk(ji,jk) * sfon(ji,jk,1)
(  412)              endif
(  413) !             do j1=1,4
(  414) !             if(sfon(ji,jk,j1+1).gt.0.) then
(  415) !             ta(ji,jj-j1,jk)= ta(ji,jj-j1,jk)  - &
(  416) !                          (ta(ji,jj-j1,jk)  - tfon(ji,jk,j1+1)) &
(  417) !                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  418) !                         * tnmsk5(ji,jk,j1+1)
(  419) !             sa(ji,jj-j1,jk)= sa(ji,jj-j1,jk)  - &
(  420) !                          (sa(ji,jj-j1,jk)  - sfon(ji,jk,j1+1)) &
(  421) !                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  422) !                         * tnmsk5(ji,jk,j1+1)
(  423) !             endif
(  424) !             end do
(  425) 
(  426)               END DO
(  427)            END DO
(  428)         END DO
(  429) 
(  430) !ylu
(  431)         j1=0
(  432)         DO jj = mj0(jpjglo-2), mj1(jpjglo-5),-1
(  433)            j1=j1+1
(  434)            DO jk = 1, jpkm1
(  435)               DO ji = 2, jpi-1
(  436)             if(sfon(ji,jk,j1+1).gt.0.) then
(  437)             ta(ji,jj,jk)= ta(ji,jj,jk)  - &
(  438)                          (ta(ji,jj,jk)  - tfon(ji,jk,j1+1)) &
(  439)                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  440)                         * tnmsk5(ji,jk,j1+1)
(  441)             sa(ji,jj,jk)= sa(ji,jj,jk)  - &
(  442)                          (sa(ji,jj,jk)  - sfon(ji,jk,j1+1)) &
(  443)                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  444)                         * tnmsk5(ji,jk,j1+1)
(  445)             endif
(  446)               END DO
(  447)            END DO
(  448)         END DO
(  449) 
(  450)      ELSE
(  451) 
(  452)      ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 9

(  453)      ! -------------------------------------------------------
(  454) 
(  455)         ! Temperature and salinity radiation
(  456)         ! ----------------------------------
(  457)         !
(  458)         !           ji-1   ji   ji   ji +1
(  459)         !             |
(  460)         !    nib //// u // T // u // T //   jpjnob + 1
(  461)         !        /////|//////////////////
(  462)         !    nib  ----f----v----f----v---   jpjnob
(  463)         !             |         |
(  464)         !      nibm-- u -- T -- u -- T --   jpjnob
(  465)         !             |         |
(  466)         !   nibm  ----f----v----f----v---  jpjnob-1
(  467)         !             |         |
(  468)         !     nibm2-- u -- T -- T -- T --  jpjnob-1
(  469)         !             |         |
(  470)         !   nibm2 ----f----v----f----v---  jpjnob-2
(  471)         !             |         |
(  472)         !
(  473)         ! ... radiative conditions + relaxation toward a climatology
(  474)         ! ... the phase velocity is taken as the normal phase velocity of the tangen-
(  475)         ! ... tial velocity (here un), which has been saved in (u_cynbnd)
(  476)         ! ... jpjnob+1,(jpindp1, jpinfm1)
(  477)         DO jj = fs_njn0+1, fs_njn1+1 ! Vector opt.
(  478)            DO jk = 1, jpkm1
(  479)               DO ji = 2, jpim1
(  480)         ! ... j-phase speed ratio (from averaged of vtnbnd)
(  481)         !        (bounded by 1)
(  482)                  z05cx = ( 0.5 * ( u_cynbnd(ji,jk) + u_cynbnd(ji-1,jk) ) ) / e2t(ji,jj-1)
(  483)                  z05cx = min( z05cx, 1. )
(  484)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  485)         !           > 0, outflow zin=1, ztau=rtaun
(  486)                  zin = sign( 1., z05cx )
(  487)                  zin = 0.5*( zin + abs(zin) )
(  488)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  489)                  ztau = (1.-zin ) * rtaunin + zin * rtaun
(  490)                  z05cx = z05cx * zin
(  491)         ! ... update (ta,sa) with radiative or climatological (t, s)
(  492)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
(  493)                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  494)                                 * tnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  495) !sujie                                 * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk) ) &
(  496)                                 * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk,1) ) &
(  497)                                 / (1. + z05cx)
(  498)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
(  499)                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  500)                                 * snbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  501) !sujie                                 * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk) ) &
(  502)                                 * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk,1) ) &
(  503)                                 / (1. + z05cx)
(  504)               END DO
(  505)            END DO
(  506)         END DO
(  507) 
(  508)      END IF
(  509) 
(  510)   END SUBROUTINE obc_tra_north






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 10

(  511) 
(  512) 
(  513)   SUBROUTINE obc_tra_south ( kt )
(  514)      !!------------------------------------------------------------------------------
(  515)      !!                ***  SUBROUTINE obc_tra_south  ***
(  516)      !!
(  517)      !! ** Purpose :
(  518)      !!      Apply the radiation algorithm on south OBC tracers ta, sa using the
(  519)      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
(  520)      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
(  521)      !!
(  522)      !!  History :
(  523)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  524)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  525)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  526)      !!         ! 00-06 (J.-M. Molines)
(  527)      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  528)      !!------------------------------------------------------------------------------
(  529)      !! * Arguments
(  530)      INTEGER, INTENT( in ) ::   kt
(  531) 
(  532)      !! * Local declaration
(  533)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  534)      REAL(wp) ::   z05cx, ztau, zin
(  535) !!DB
(  536)      INTEGER :: jji, jjj, jjk
(  537) 
(  538)      !!------------------------------------------------------------------------------
(  539) 
(  540)      ! 1. First three time steps and more if lfbcsouth is .TRUE.
(  541)      !    In that case open boundary conditions are FIXED.
(  542)      ! --------------------------------------------------------
(  543) 
(  544)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth ) THEN
(  545) 
(  546) 
(  547)         DO jj = fs_njs0, fs_njs1  ! Vector opt.
(  548)            DO jk = 1, jpkm1
(  549)               DO ji = 1, jpi
(  550)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
(  551)                                tsmsk(ji,jk) * tfos(ji,jk,1)
(  552)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
(  553)                                tsmsk(ji,jk) * sfos(ji,jk,1)
(  554)               END DO
(  555)            END DO
(  556)         END DO
(  557) 
(  558) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  559) !!
(  560)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  561) !        zin = 2. * rdt / (20.0*rday)
(  562) !        do jji = njs0, njs1   !!isolate processor 
(  563)         do jji = 1, 0   !!isolate processor 
(  564)            do jj = njs0+1, njs0+4 !!5 cells - 1
(  565) !!restore based on position -- but do not use
(  566)               ztau = zin * exp(float( -(jj-njs0)/4 )) 
(  567) !!DBG
(  568) !              zin = ztau






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 11

(  569)               do jk = 1, jpkm1
(  570)                  do ji = 1, jpi
(  571)                     ta(ji,jj,jk) = ta(ji,jj,jk) + tsmsk(ji,jk)*&
(  572)                          zin*(tfos(ji,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  573)                     sa(ji,jj,jk) = sa(ji,jj,jk) + tsmsk(ji,jk)*&
(  574)                          zin*(sfos(ji,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  575)                  enddo
(  576)               enddo
(  577)            enddo
(  578)         enddo
(  579) 
(  580) 
(  581) 
(  582)      ELSE
(  583) 
(  584)      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
(  585)      ! -------------------------------------------------------
(  586) 
(  587)         ! Temperature and salinity radiation
(  588)         ! ----------------------------------
(  589)         !
(  590)         !           ji-1   ji   ji   ji +1
(  591)         !             |         |
(  592)         !   nibm2 ----f----v----f----v---   jpjsob+2
(  593)         !             |         |
(  594)         !   nibm2 --  u -- T -- u -- T --   jpjsob+2
(  595)         !             |         |
(  596)         !   nibm  ----f----v----f----v---   jpjsob+1
(  597)         !             |         |
(  598)         !    nibm --  u -- T -- T -- T --   jpjsob+1
(  599)         !             |         |
(  600)         !   nib  -----f----v----f----v---   jpjsob
(  601)         !       //////|/////////|////////
(  602)         !    nib //// u // T // u // T //   jpjsob
(  603)         !
(  604)         !... radiative conditions + relaxation toward a climatology
(  605)         !... the phase velocity is taken as the phase velocity of the tangen-
(  606)         !... tial velocity (here un), which has been saved in (u_cysbnd)
(  607)         !... jpjsob,(jpisdp1, jpisfm1)
(  608)         DO jj = fs_njs0, fs_njs1  ! Vector opt.
(  609)            DO jk = 1, jpkm1
(  610)               DO ji = 2, jpim1
(  611)         !... j-phase speed ratio (from averaged of u_cysbnd)
(  612)         !       (bounded by 1)
(  613)                  z05cx = ( 0.5 * ( u_cysbnd(ji,jk) + u_cysbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
(  614)                  z05cx = max( z05cx, -1. )
(  615)         !... z05cx > 0, inflow  zin=0, ztau=1
(  616)         !          < 0, outflow zin=1, ztau=rtaus
(  617)                  zin = sign( 1., -1.* z05cx )
(  618)                  zin = 0.5*( zin + abs(zin) )
(  619)                  ztau = (1.-zin ) + zin * rtaus
(  620)                  z05cx = z05cx * zin
(  621)         !... update (ta,sa) with radiative or climatological (t, s)
(  622)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
(  623)                                 tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
(  624)                                 * tsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
(  625) !sujie                                 * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk) ) &
(  626)                                 * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk,1) ) &






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 12

(  627)                                 / (1. - z05cx)
(  628)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
(  629)                                 tsmsk(ji,jk) * (  ( 1. + z05cx - ztau )        &
(  630)                                 * ssbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
(  631) !sujie                                 * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk) ) &
(  632)                                 * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk,1) ) &
(  633)                                 / (1. - z05cx)
(  634)               END DO
(  635)            END DO
(  636)         END DO
(  637) 
(  638)      END IF
(  639) 
(  640)   END SUBROUTINE obc_tra_south
(  641) 
(  642) #else
(  643)   !!---------------------------------------------------------------------------------
(  644)   !!   Default option                                                    Empty module
(  645)   !!---------------------------------------------------------------------------------
(  646) CONTAINS
(  647)   SUBROUTINE obc_tra      ! Empty routine
(  648)   END SUBROUTINE obc_tra
(  649) #endif
(  650) 
(  651)   !!=================================================================================
(  652) END MODULE obctra






































PGF90 (Version     10.5)          09/23/2010  09:50:08      page 13

(    1) # 1 "obctra.F90"
(    1) MODULE obctra
(    2)   !!=================================================================================
(    3)   !!                       ***  MODULE  obctra  ***
(    4)   !! Ocean tracers:   Radiation of tracers on each open boundary
(    5)   !!=================================================================================
(    6) # 7
(    7)   !!---------------------------------------------------------------------------------
(    8)   !!   'key_obc'      :                                      Open Boundary Conditions
(    9)   !!---------------------------------------------------------------------------------
(   10)   !!   obc_tra        : call the subroutine for each open boundary
(   11)   !!   obc_tra_east   : radiation of the east open boundary tracers
(   12)   !!   obc_tra_west   : radiation of the west open boundary tracers
(   13)   !!   obc_tra_north  : radiation of the north open boundary tracers
(   14)   !!   obc_tra_south  : radiation of the south open boundary tracers
(   15)   !!----------------------------------------------------------------------------------
(   16)   !! * Modules used
(   17)   USE oce             ! ocean dynamics and tracers variables
(   18)   USE dom_oce         ! ocean space and time domain variables
(   19)   USE phycst          ! physical constants
(   20)   USE obc_oce         ! ocean open boundary conditions
(   21)   USE lib_mpp         ! ???
(   22)   USE lbclnk          ! ???
(   23)   USE in_out_manager  ! I/O manager
(   24) 
(   25)   IMPLICIT NONE
(   26)   PRIVATE
(   27) 
(   28)   !! * Accessibility
(   29)   PUBLIC obc_tra     ! routine called in tranxt.F90
(   30) 
(   31)   !! * Module variables
(   32)   INTEGER ::      & ! ... boundary space indices
(   33)      nib   = 1,   & ! nib   = boundary point
(   34)      nibm  = 2,   & ! nibm  = 1st interior point
(   35)      nibm2 = 3,   & ! nibm2 = 2nd interior point
(   36)                     ! ... boundary time indices
(   37)      nit   = 1,   & ! nit    = now
(   38)      nitm  = 2,   & ! nitm   = before
(   39)      nitm2 = 3      ! nitm2  = before-before
(   40) 
(   41)   REAL(wp) ::     &
(   42)      rtaue  , rtauw  , rtaun  , rtaus  ,  &  ! Boundary restoring coefficient
(   43)      rtauein, rtauwin, rtaunin, rtausin      ! Boundary restoring coefficient for inflow
(   44) 
(   45)   !! * Substitutions
(   46) # 1 "./obc_vectopt_loop_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 14

(   12) # 47 "obctra.F90"
(   47) # 47
(   47)   !!---------------------------------------------------------------------------------
(   48)   !!   OPA 9.0 , LOCEAN-IPSL (2005)
(   49)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obctra.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $
(   50)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   51)   !!---------------------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55)   SUBROUTINE obc_tra( kt )
(   56)      !!-------------------------------------------------------------------------------
(   57)      !!                 ***  SUBROUTINE obc_tra  ***
(   58)      !!
(   59)      !! ** Purpose :   Compute tracer fields (t,s) along the open boundaries.
(   60)      !!      This routine is called by the tranxt.F routine and updates ta,sa
(   61)      !!      which are the actual temperature and salinity fields.
(   62)      !!        The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   63)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   64)      !!      open one (must be done in the param_obc.h90 file).
(   65)      !!
(   66)      !! Reference :
(   67)      !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   68)      !!
(   69)      !!  History :
(   70)      !!        !  95-03 (J.-M. Molines) Original, SPEM
(   71)      !!        !  97-07 (G. Madec, J.-M. Molines) addition
(   72)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) F90
(   73)      !!----------------------------------------------------------------------
(   74)      !! * Arguments
(   75)      INTEGER, INTENT( in ) ::   kt
(   76)      !!----------------------------------------------------------------------
(   77) 
(   78)      ! 0. Local constant initialization
(   79) 
(   80)      IF( kt == nit000 .OR. ln_rstart) THEN
(   81)         ! ... Boundary restoring coefficient
(   82)         rtaue = 2. * rdt / rdpeob
(   83)         rtauw = 2. * rdt / rdpwob
(   84)         rtaun = 2. * rdt / rdpnob
(   85)         rtaus = 2. * rdt / rdpsob
(   86)         ! ... Boundary restoring coefficient for inflow ( all boundaries)
(   87)         rtauein = 2. * rdt / rdpein
(   88)         rtauwin = 2. * rdt / rdpwin
(   89)         rtaunin = 2. * rdt / rdpnin
(   90)         rtausin = 2. * rdt / rdpsin
(   91)      END IF
(   92) 
(   93)      IF( lp_obc_east  )   CALL obc_tra_east ( kt )    ! East open boundary
(   94) 
(   95)      IF( lp_obc_west  )   CALL obc_tra_west ( kt )    ! West open boundary
(   96) 
(   97)      IF( lp_obc_north )   CALL obc_tra_north( kt )    ! North open boundary
(   98) 
(   99)      IF( lp_obc_south )   CALL obc_tra_south( kt )    ! South open boundary
(  100) 
(  101)      IF( lk_mpp ) THEN                  !!bug ???
(  102)         IF( kt >= nit000+3 .AND. ln_rstart ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 15

(  103)            CALL lbc_lnk( tb, 'T', 1. )
(  104)            CALL lbc_lnk( sb, 'T', 1. )
(  105)         END IF
(  106)         CALL lbc_lnk( ta, 'T', 1. )
(  107)         CALL lbc_lnk( sa, 'T', 1. )
(  108)      ENDIF
(  109) 
(  110)   END SUBROUTINE obc_tra
(  111) 
(  112) 
(  113)   SUBROUTINE obc_tra_east ( kt )
(  114)      !!------------------------------------------------------------------------------
(  115)      !!                ***  SUBROUTINE obc_tra_east  ***
(  116)      !!
(  117)      !! ** Purpose :
(  118)      !!      Apply the radiation algorithm on east OBC tracers ta, sa using the
(  119)      !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
(  120)      !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
(  121)      !!
(  122)      !!  History :
(  123)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  124)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  125)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  126)      !!         ! 00-06 (J.-M. Molines)
(  127)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  128)      !!------------------------------------------------------------------------------
(  129)      !! * Arguments
(  130)      INTEGER, INTENT( in ) ::   kt
(  131) 
(  132)      !! * Local declaration
(  133)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  134)      REAL(wp) ::   z05cx, ztau, zin
(  135) !!DB
(  136)      INTEGER :: jji, jjj, jjk
(  137) 
(  138)      !!------------------------------------------------------------------------------
(  139) 
(  140)      ! 1. First three time steps and more if lfbceast is .TRUE.
(  141)      !    In that case open boundary conditions are FIXED.
(  142)      ! --------------------------------------------------------
(  143) 
(  144)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast ) THEN
(  145)         DO ji = nie0 +1, nie1+1 ! Vector opt.
(  146)            DO jk = 1, jpkm1
(  147)               DO jj = 1, jpj
(  148)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-temsk(jj,jk)) + &
(  149)                                temsk(jj,jk) * tfoe(jj,jk,1)
(  150)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-temsk(jj,jk)) + &
(  151)                                temsk(jj,jk) * sfoe(jj,jk,1)
(  152)               END DO
(  153)            END DO
(  154)         END DO
(  155) 
(  156) !!DB DBG TS sponge; 5 cells; timescale = 1day(?)
(  157) !!
(  158)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  159) !        zin = 2. * rdt / (20.0*rday)
(  160) !        do jji = nie0, nie1   !!isolate processor 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 16

(  161) !DBG -- skip 
(  162)         do jji = 1, 0 
(  163)            do ji = nie0-4, nie0-1 !!5 cells - 1
(  164) !restore based on position -- but do not use
(  165)               ztau = zin * exp(float( -(nie0-ji)/4 )) 
(  166)               do jk = 1, jpkm1
(  167)                  do jj = 1, jpj
(  168)                     ta(ji,jj,jk) = ta(ji,jj,jk) + temsk(jj,jk)*&
(  169)                          zin*(tfoe(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  170)                     sa(ji,jj,jk) = sa(ji,jj,jk) + temsk(jj,jk)*&
(  171)                          zin*(sfoe(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  172)                  enddo
(  173)               enddo
(  174)            enddo
(  175)         enddo
(  176) 
(  177) 
(  178) 
(  179) 
(  180) 
(  181) 
(  182) 
(  183)      ELSE
(  184) 
(  185)      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
(  186)      ! -----------------------------------------------------
(  187) 
(  188)         ! Temperature and salinity radiation
(  189)         ! ----------------------------------
(  190)         !
(  191)         !            nibm2      nibm      nib
(  192)         !              |   nibm  |   nib///|///
(  193)         !              |    |    |    |////|///
(  194)         !  jj   line --v----f----v----f----v---
(  195)         !              |    |    |    |////|///
(  196)         !                   |         |///   //
(  197)         !  jj   line   T    u    T    u/// T //
(  198)         !                   |         |///   //
(  199)         !              |    |    |    |////|///
(  200)         !  jj-1 line --v----f----v----f----v---
(  201)         !              |    |    |    |////|///
(  202)         !                jpieob-1    jpieob / ///
(  203)         !              |         |         |
(  204)         !           jpieob-1    jpieob     jpieob+1
(  205)         !
(  206)         ! ... radiative conditions + relaxation toward a climatology
(  207)         !     the phase velocity is taken as the phase velocity of the tangen-
(  208)         !     tial velocity (here vn), which have been saved in (u_cxebnd,v_cxebnd)
(  209)         ! ... (jpjedp1, jpjefm1), jpieob+1
(  210)         DO ji = nie0 +1, nie1+1 ! Vector opt.
(  211)            DO jk = 1, jpkm1
(  212)               DO jj = 2, jpjm1
(  213)         ! ... i-phase speed ratio (from averaged of v_cxebnd)
(  214)                  z05cx = ( 0.5 * ( v_cxebnd(jj,jk) + v_cxebnd(jj-1,jk) ) ) / e1t(ji-1,jj)
(  215)                  z05cx = min( z05cx, 1. )
(  216)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  217)         !           > 0, outflow zin=1, ztau=rtaue
(  218)                  zin = sign( 1., z05cx )






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 17

(  219)                  zin = 0.5*( zin + abs(zin) )
(  220)         ! ... for inflow rtauein is used for relaxation coefficient else rtaue
(  221)                  ztau = (1.-zin ) * rtauein  + zin * rtaue
(  222)                  z05cx = z05cx * zin
(  223)         ! ... update ( ta, sa ) with radiative or climatological (t, s)
(  224)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
(  225)                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  226)                                 * tebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  227) !byoung                                 * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk) ) &
(  228)                                 * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk,1) ) &
(  229)                                 / (1. + z05cx)
(  230)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
(  231)                                 temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
(  232)                                 * sebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
(  233) !byoung                                * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk) ) &
(  234)                                 * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk,1) ) &
(  235)                                 / (1. + z05cx)
(  236)               END DO
(  237)            END DO
(  238)         END DO
(  239) 
(  240)      END IF
(  241) 
(  242)   END SUBROUTINE obc_tra_east
(  243) 
(  244) 
(  245)   SUBROUTINE obc_tra_west ( kt )
(  246)      !!------------------------------------------------------------------------------
(  247)      !!                 ***  SUBROUTINE obc_tra_west  ***
(  248)      !!
(  249)      !! ** Purpose :
(  250)      !!      Apply the radiation algorithm on west OBC tracers ta, sa using the
(  251)      !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
(  252)      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
(  253)      !!
(  254)      !!  History :
(  255)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  256)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  257)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  258)      !!         ! 00-06 (J.-M. Molines)
(  259)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  260)      !!------------------------------------------------------------------------------
(  261)      !! * Arguments
(  262)      INTEGER, INTENT( in ) ::   kt
(  263) 
(  264)      !! * Local declaration
(  265)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  266)      REAL(wp) ::   z05cx, ztau, zin
(  267) !!DB
(  268)      INTEGER :: jji, jjj, jjk
(  269) 
(  270)      !!------------------------------------------------------------------------------
(  271) 
(  272)      ! 1. First three time steps and more if lfbcwest is .TRUE.
(  273)      !    In that case open boundary conditions are FIXED.
(  274)      ! --------------------------------------------------------
(  275) 
(  276)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 18

(  277) 
(  278)         DO ji = niw0 , niw1 ! Vector opt.
(  279)            DO jk = 1, jpkm1
(  280)               DO jj = 1, jpj
(  281)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
(  282)                                twmsk(jj,jk) * tfow(jj,jk,1)
(  283)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
(  284)                                twmsk(jj,jk) * sfow(jj,jk,1)
(  285)               END DO
(  286)            END DO
(  287)         END DO
(  288) !!DB DBG TS sponge; 5 cells; timescale = 1day(?)
(  289) !!
(  290)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  291) !        zin = 2. * rdt / (20.0*rday)
(  292) !        do jji = niw0, niw1   !!isolate processor 
(  293)         do jji = 1, 0    !!isolate processor 
(  294)            do ji = niw0+1, niw0+4 !!5 cells - 1
(  295) !!restore based on position -- but do not use
(  296)               ztau = zin * exp(float( -(ji-niw0)/4 )) 
(  297)               do jk = 1, jpkm1
(  298)                  do jj = 1, jpj
(  299)                     ta(ji,jj,jk) = ta(ji,jj,jk) + twmsk(jj,jk)*&
(  300)                          zin*(tfow(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  301)                     sa(ji,jj,jk) = sa(ji,jj,jk) + twmsk(jj,jk)*&
(  302)                          zin*(sfow(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  303)                  enddo
(  304)               enddo
(  305)            enddo
(  306)         enddo
(  307) 
(  308)      ELSE
(  309) 
(  310)      ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
(  311)      ! -----------------------------------------------------
(  312) 
(  313)         ! Temperature and salinity radiation
(  314)         ! ----------------------------------
(  315)         !
(  316)         !          nib       nibm     nibm2
(  317)         !     nib///|   nibm  |  nibm2  |
(  318)         !   ///|////|    |    |    |    |
(  319)         !   ---v----f----v----f----v----f-- jj   line
(  320)         !   ///|////|    |    |    |    |
(  321)         !   //   ///|         |         |
(  322)         !   // T ///u    T    u    T    u   jj   line
(  323)         !   //   ///|         |         |
(  324)         !   ///|////|    |    |    |    |
(  325)         !   ---v----f----v----f----v----f-- jj-1 line
(  326)         !   ///|////|    |    |    |    |
(  327)         !         jpiwob    jpiwob+1    jpiwob+2
(  328)         !      |         |         |
(  329)         !    jpiwob    jpiwob+1   jpiwob+2
(  330)         !
(  331)         ! ... radiative conditions + relaxation toward a climatology
(  332)         ! ... the phase velocity is taken as the phase velocity of the tangen-
(  333)         ! ... tial velocity (here vn), which have been saved in (v_cxwbnd)
(  334)         DO ji = niw0 , niw1 ! Vector opt.






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 19

(  335)            DO jk = 1, jpkm1
(  336)               DO jj = 2, jpjm1
(  337)         ! ... i-phase speed ratio (from averaged of v_cxwbnd)
(  338)                  z05cx = (  0.5 * ( v_cxwbnd(jj,jk) + v_cxwbnd(jj-1,jk) ) ) / e1t(ji+1,jj)
(  339)                  z05cx = max( z05cx, -1. )
(  340)         ! ... z05cx > 0, inflow  zin=0, ztau=1
(  341)         !           < 0, outflow zin=1, ztau=rtauw
(  342)                  zin = sign( 1., -1.* z05cx )
(  343)                  zin = 0.5*( zin + abs(zin) )
(  344)                  ztau = (1.-zin )*rtauwin + zin * rtauw
(  345)                  z05cx = z05cx * zin
(  346)         ! ... update (ta,sa) with radiative or climatological (t, s)
(  347)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
(  348)                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  349)                                 * twbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  350) !byoung                                * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk) ) &
(  351)                                 * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk,1) ) &
(  352)                                 / (1. - z05cx)
(  353)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
(  354)                                 twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
(  355)                                 * swbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
(  356) !byoung                                * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk) ) &
(  357)                                 * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk,1) ) &
(  358)                                 / (1. - z05cx)
(  359)               END DO
(  360)            END DO
(  361)         END DO
(  362) 
(  363)      END IF
(  364) 
(  365)   END SUBROUTINE obc_tra_west
(  366) 
(  367) 
(  368)   SUBROUTINE obc_tra_north ( kt )
(  369)      !!------------------------------------------------------------------------------
(  370)      !!                 ***  SUBROUTINE obc_tra_north  ***
(  371)      !!
(  372)      !! ** Purpose :
(  373)      !!      Apply the radiation algorithm on north OBC tracers ta, sa using the
(  374)      !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
(  375)      !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
(  376)      !!
(  377)      !!  History :
(  378)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  379)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  380)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  381)      !!         ! 00-06 (J.-M. Molines)
(  382)      !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(  383)      !!------------------------------------------------------------------------------
(  384)      !! * Arguments
(  385)      INTEGER, INTENT( in ) ::   kt
(  386) 
(  387)      !! * Local declaration
(  388) !sujie      INTEGER ::   ji, jj, jk      ! dummy loop indices
(  389)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  390)      REAL(wp) ::   z05cx, ztau, zin
(  391)      !!------------------------------------------------------------------------------
(  392) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 20

(  393)      ! 1. First three time steps and more if lfbcnorth is .TRUE.
(  394)      !    In that case open boundary conditions are FIXED.
(  395)      ! --------------------------------------------------------
(  396) 
(  397)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth ) THEN
(  398) 
(  399) !sujie
(  400)         DO jj = njn0+1, njn1+1  ! Vector opt.
(  401)            DO jk = 1, jpkm1
(  402)               DO ji = 1, jpi
(  403) !                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  404) !                                tnmsk(ji,jk) * tfon(ji,jk)
(  405) !                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  406) !                                tnmsk(ji,jk) * sfon(ji,jk)
(  407)              if(sfon(ji,jk,1).gt.0.) then
(  408)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  409)                                tnmsk(ji,jk) * tfon(ji,jk,1)
(  410)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
(  411)                                tnmsk(ji,jk) * sfon(ji,jk,1)
(  412)              endif
(  413) !             do j1=1,4
(  414) !             if(sfon(ji,jk,j1+1).gt.0.) then
(  415) !             ta(ji,jj-j1,jk)= ta(ji,jj-j1,jk)  - &
(  416) !                          (ta(ji,jj-j1,jk)  - tfon(ji,jk,j1+1)) &
(  417) !                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  418) !                         * tnmsk5(ji,jk,j1+1)
(  419) !             sa(ji,jj-j1,jk)= sa(ji,jj-j1,jk)  - &
(  420) !                          (sa(ji,jj-j1,jk)  - sfon(ji,jk,j1+1)) &
(  421) !                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  422) !                         * tnmsk5(ji,jk,j1+1)
(  423) !             endif
(  424) !             end do
(  425) 
(  426)               END DO
(  427)            END DO
(  428)         END DO
(  429) 
(  430) !ylu
(  431)         j1=0
(  432)         DO jj = mj0(jpjglo-2), mj1(jpjglo-5),-1
(  433)            j1=j1+1
(  434)            DO jk = 1, jpkm1
(  435)               DO ji = 2, jpi-1
(  436)             if(sfon(ji,jk,j1+1).gt.0.) then
(  437)             ta(ji,jj,jk)= ta(ji,jj,jk)  - &
(  438)                          (ta(ji,jj,jk)  - tfon(ji,jk,j1+1)) &
(  439)                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  440)                         * tnmsk5(ji,jk,j1+1)
(  441)             sa(ji,jj,jk)= sa(ji,jj,jk)  - &
(  442)                          (sa(ji,jj,jk)  - sfon(ji,jk,j1+1)) &
(  443)                         * rdt/(86400.* exp(float(j1+1))/15.0) &
(  444)                         * tnmsk5(ji,jk,j1+1)
(  445)             endif
(  446)               END DO
(  447)            END DO
(  448)         END DO
(  449) 
(  450)      ELSE






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 21

(  451) 
(  452)      ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
(  453)      ! -------------------------------------------------------
(  454) 
(  455)         ! Temperature and salinity radiation
(  456)         ! ----------------------------------
(  457)         !
(  458)         !           ji-1   ji   ji   ji +1
(  459)         !             |
(  460)         !    nib //// u // T // u // T //   jpjnob + 1
(  461)         !        /////|//////////////////
(  462)         !    nib  ----f----v----f----v---   jpjnob
(  463)         !             |         |
(  464)         !      nibm-- u -- T -- u -- T --   jpjnob
(  465)         !             |         |
(  466)         !   nibm  ----f----v----f----v---  jpjnob-1
(  467)         !             |         |
(  468)         !     nibm2-- u -- T -- T -- T --  jpjnob-1
(  469)         !             |         |
(  470)         !   nibm2 ----f----v----f----v---  jpjnob-2
(  471)         !             |         |
(  472)         !
(  473)         ! ... radiative conditions + relaxation toward a climatology
(  474)         ! ... the phase velocity is taken as the normal phase velocity of the tangen-
(  475)         ! ... tial velocity (here un), which has been saved in (u_cynbnd)
(  476)         ! ... jpjnob+1,(jpindp1, jpinfm1)
(  477)         DO jj = njn0+1, njn1+1 ! Vector opt.
(  478)            DO jk = 1, jpkm1
(  479)               DO ji = 2, jpim1
(  480)         ! ... j-phase speed ratio (from averaged of vtnbnd)
(  481)         !        (bounded by 1)
(  482)                  z05cx = ( 0.5 * ( u_cynbnd(ji,jk) + u_cynbnd(ji-1,jk) ) ) / e2t(ji,jj-1)
(  483)                  z05cx = min( z05cx, 1. )
(  484)         ! ... z05cx=< 0, inflow  zin=0, ztau=1
(  485)         !           > 0, outflow zin=1, ztau=rtaun
(  486)                  zin = sign( 1., z05cx )
(  487)                  zin = 0.5*( zin + abs(zin) )
(  488)         ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
(  489)                  ztau = (1.-zin ) * rtaunin + zin * rtaun
(  490)                  z05cx = z05cx * zin
(  491)         ! ... update (ta,sa) with radiative or climatological (t, s)
(  492)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
(  493)                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  494)                                 * tnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  495) !sujie                                 * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk) ) &
(  496)                                 * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk,1) ) &
(  497)                                 / (1. + z05cx)
(  498)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
(  499)                                 tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
(  500)                                 * snbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
(  501) !sujie                                 * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk) ) &
(  502)                                 * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk,1) ) &
(  503)                                 / (1. + z05cx)
(  504)               END DO
(  505)            END DO
(  506)         END DO
(  507) 
(  508)      END IF






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 22

(  509) 
(  510)   END SUBROUTINE obc_tra_north
(  511) 
(  512) 
(  513)   SUBROUTINE obc_tra_south ( kt )
(  514)      !!------------------------------------------------------------------------------
(  515)      !!                ***  SUBROUTINE obc_tra_south  ***
(  516)      !!
(  517)      !! ** Purpose :
(  518)      !!      Apply the radiation algorithm on south OBC tracers ta, sa using the
(  519)      !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
(  520)      !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
(  521)      !!
(  522)      !!  History :
(  523)      !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  524)      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  525)      !!         ! 97-12 (M. Imbard) Mpp adaptation
(  526)      !!         ! 00-06 (J.-M. Molines)
(  527)      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  528)      !!------------------------------------------------------------------------------
(  529)      !! * Arguments
(  530)      INTEGER, INTENT( in ) ::   kt
(  531) 
(  532)      !! * Local declaration
(  533)      INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
(  534)      REAL(wp) ::   z05cx, ztau, zin
(  535) !!DB
(  536)      INTEGER :: jji, jjj, jjk
(  537) 
(  538)      !!------------------------------------------------------------------------------
(  539) 
(  540)      ! 1. First three time steps and more if lfbcsouth is .TRUE.
(  541)      !    In that case open boundary conditions are FIXED.
(  542)      ! --------------------------------------------------------
(  543) 
(  544)      IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth ) THEN
(  545) 
(  546) 
(  547)         DO jj = njs0, njs1  ! Vector opt.
(  548)            DO jk = 1, jpkm1
(  549)               DO ji = 1, jpi
(  550)                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
(  551)                                tsmsk(ji,jk) * tfos(ji,jk,1)
(  552)                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
(  553)                                tsmsk(ji,jk) * sfos(ji,jk,1)
(  554)               END DO
(  555)            END DO
(  556)         END DO
(  557) 
(  558) !!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
(  559) !!
(  560)         zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
(  561) !        zin = 2. * rdt / (20.0*rday)
(  562) !        do jji = njs0, njs1   !!isolate processor 
(  563)         do jji = 1, 0   !!isolate processor 
(  564)            do jj = njs0+1, njs0+4 !!5 cells - 1
(  565) !!restore based on position -- but do not use
(  566)               ztau = zin * exp(float( -(jj-njs0)/4 )) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 23

(  567) !!DBG
(  568) !              zin = ztau
(  569)               do jk = 1, jpkm1
(  570)                  do ji = 1, jpi
(  571)                     ta(ji,jj,jk) = ta(ji,jj,jk) + tsmsk(ji,jk)*&
(  572)                          zin*(tfos(ji,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
(  573)                     sa(ji,jj,jk) = sa(ji,jj,jk) + tsmsk(ji,jk)*&
(  574)                          zin*(sfos(ji,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
(  575)                  enddo
(  576)               enddo
(  577)            enddo
(  578)         enddo
(  579) 
(  580) 
(  581) 
(  582)      ELSE
(  583) 
(  584)      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
(  585)      ! -------------------------------------------------------
(  586) 
(  587)         ! Temperature and salinity radiation
(  588)         ! ----------------------------------
(  589)         !
(  590)         !           ji-1   ji   ji   ji +1
(  591)         !             |         |
(  592)         !   nibm2 ----f----v----f----v---   jpjsob+2
(  593)         !             |         |
(  594)         !   nibm2 --  u -- T -- u -- T --   jpjsob+2
(  595)         !             |         |
(  596)         !   nibm  ----f----v----f----v---   jpjsob+1
(  597)         !             |         |
(  598)         !    nibm --  u -- T -- T -- T --   jpjsob+1
(  599)         !             |         |
(  600)         !   nib  -----f----v----f----v---   jpjsob
(  601)         !       //////|/////////|////////
(  602)         !    nib //// u // T // u // T //   jpjsob
(  603)         !
(  604)         !... radiative conditions + relaxation toward a climatology
(  605)         !... the phase velocity is taken as the phase velocity of the tangen-
(  606)         !... tial velocity (here un), which has been saved in (u_cysbnd)
(  607)         !... jpjsob,(jpisdp1, jpisfm1)
(  608)         DO jj = njs0, njs1  ! Vector opt.
(  609)            DO jk = 1, jpkm1
(  610)               DO ji = 2, jpim1
(  611)         !... j-phase speed ratio (from averaged of u_cysbnd)
(  612)         !       (bounded by 1)
(  613)                  z05cx = ( 0.5 * ( u_cysbnd(ji,jk) + u_cysbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
(  614)                  z05cx = max( z05cx, -1. )
(  615)         !... z05cx > 0, inflow  zin=0, ztau=1
(  616)         !          < 0, outflow zin=1, ztau=rtaus
(  617)                  zin = sign( 1., -1.* z05cx )
(  618)                  zin = 0.5*( zin + abs(zin) )
(  619)                  ztau = (1.-zin ) + zin * rtaus
(  620)                  z05cx = z05cx * zin
(  621)         !... update (ta,sa) with radiative or climatological (t, s)
(  622)                  ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
(  623)                                 tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
(  624)                                 * tsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 24

(  625) !sujie                                 * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk) ) &
(  626)                                 * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk,1) ) &
(  627)                                 / (1. - z05cx)
(  628)                  sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
(  629)                                 tsmsk(ji,jk) * (  ( 1. + z05cx - ztau )        &
(  630)                                 * ssbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
(  631) !sujie                                 * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk) ) &
(  632)                                 * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk,1) ) &
(  633)                                 / (1. - z05cx)
(  634)               END DO
(  635)            END DO
(  636)         END DO
(  637) 
(  638)      END IF
(  639) 
(  640)   END SUBROUTINE obc_tra_south
(  641) 
(  642) 
(  643) # 651
(  651)   !!=================================================================================
(  652) END MODULE obctra
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (obctra.F90: 165)
