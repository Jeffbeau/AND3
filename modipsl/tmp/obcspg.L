


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

Filename: obcspg.F90

(    1) MODULE obcspg 
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  obcspg  ***
(    4)    !! Open Boundaries  :   Radiation of barotropic stream function on each
(    5)    !!                      open boundary
(    6)    !!======================================================================
(    7) #if   defined key_obc   &&   defined key_dynspg_rl
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_obc'    and                            Open Boundary Condition
(   10)    !!   'key_dynspg_rl'                                 Rigid-Lid
(   11)    !!----------------------------------------------------------------------
(   12)    !!   obc_spg       : call the subroutine for each open boundary
(   13)    !!   obc_spg_east  : radiation of the east open boundary streamfunction
(   14)    !!   obc_spg_west  : radiation of the west open boundary streamfunction
(   15)    !!   obc_spg_north : radiation of the north open boundary streamfunction
(   16)    !!   obc_spg_south : radiation of the south open boundary streamfunction
(   17)    !!----------------------------------------------------------------------
(   18)    !! * Modules used
(   19)    USE oce             ! ocean dynamics and tracers variables
(   20)    USE dom_oce         ! ocean space and time domain variables 
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE phycst          ! physical constants
(   23)    USE obc_oce         ! ocean open boundary conditions
(   24)    USE lib_mpp         ! for mppobc
(   25)    USE in_out_manager  ! I/O manager
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC obc_spg     ! routine called in step.F90 (rigid lid case)
(   32) 
(   33)    !! * Module variables
(   34)    INTEGER ::   ji, jj, jk, jnic   ! dummy loop indices
(   35) 
(   36)    INTEGER ::      & ! ... boundary space indices 
(   37)       nib   = 1,   & ! nib   = boundary point
(   38)       nibm  = 2,   & ! nibm  = 1st interior point
(   39)       nibm2 = 3,   & ! nibm2 = 2nd interior point
(   40)       !              ! ... boundary time indices 
(   41)       nit   = 1,   & ! nit    = now
(   42)       nitm  = 2,   & ! nitm   = before
(   43)       nitm2 = 3      ! nitm2  = before-before
(   44) 
(   45)    REAL(wp) ::   rtaue  , rtauw  , rtaun  , rtaus  ,   &  !
(   46)                  rtauein, rtauwin, rtaunin, rtausin       !






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 2

(   47) 
(   48)    !! * Substitutions
(   49) #  include "obc_vectopt_loop_substitute.h90"
(   50)    !!----------------------------------------------------------------------
(   51)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   52)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcspg.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $ 
(   53)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56) CONTAINS
(   57) 
(   58)    SUBROUTINE obc_spg ( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                    ***  ROUTINE obc_spg  ***
(   61)       !!
(   62)       !! **  Purpose :
(   63)       !!       Compute now barotropic stream function at the open boundaries.
(   64)       !!       (lp_obc_east, and/or lp_obc_west, and/or lp_obc_north, and/or lp_obc_south).
(   65)       !!       Deduce the correct bsf trend on the open boundaries and isolated 
(   66)       !!       coastlines previous to the call to the barotropic solver.
(   67)       !!
(   68)       !! ** Method :
(   69)       !!      In case of open boundaries, there can be a net barotropic flow
(   70)       !!      through the boundaries, hence the potential on the coastlines
(   71)       !!      on each side of the OBC is different.
(   72)       !!      This routine:
(   73)       !!           1. compute the contribution of the isolated coastlines to the
(   74)       !!              rhs of the barotropic equation
(   75)       !!           2. compute the contribution of the OBC to the rhs of the
(   76)       !!              barotropic equation using a radiation equation as explained
(   77)       !!              in the OBC routine.
(   78)       !!
(   79)       !! Reference : 
(   80)       !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
(   81)       !! History :
(   82)       !!        ! 95-03 (J.-M. Molines) Original, SPEM
(   83)       !!        ! 97-07 (G. Madec, J.-M. Molines) additions
(   84)       !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
(   85)       !!----------------------------------------------------------------------
(   86)       !! * Arguments
(   87)       INTEGER, INTENT( in ) ::   kt
(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       IF( kt == nit000 .OR. ln_rstart ) THEN      ! Initialization
(   91)          ! ... Boundary restoring coefficient
(   92)          rtaue = 2. * rdt / rdpeob
(   93)          rtauw = 2. * rdt / rdpwob
(   94)          rtaun = 2. * rdt / rdpnob
(   95)          rtaus = 2. * rdt / rdpsob
(   96)          ! ... Boundary restoring coefficient for inflow ( all boundaries)
(   97)          rtauein = 2. * rdt / rdpein 
(   98)          rtauwin = 2. * rdt / rdpwin
(   99)          rtaunin = 2. * rdt / rdpnin
(  100)          rtausin = 2. * rdt / rdpsin 
(  101)       ENDIF
(  102) 
(  103)       ! right hand side of the barotropic elliptic equation
(  104)       ! ---------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 3

(  105) 
(  106)       ! Isolated coastline contribution to the RHS of the barotropic Eq.
(  107)       gcbob(:,:) = 0.e0
(  108)       DO jnic = 1, nbobc-1
(  109)          gcbob(:,:) = gcbob(:,:) + gcfobc(:,:,jnic) * gcbic(jnic)
(  110)       END DO
(  111) 
(  112)       IF( lp_obc_east  )   CALL obc_spg_east ( kt )    ! East open boundary
(  113) 
(  114)       IF( lp_obc_west  )   CALL obc_spg_west ( kt )    ! West open boundary
(  115) 
(  116)       IF( lp_obc_north )   CALL obc_spg_north( kt )    ! North open boundary
(  117) 
(  118)       IF( lp_obc_south )   CALL obc_spg_south( kt )    ! South open boundary
(  119) 
(  120)       IF( lk_mpp )   CALL lbc_lnk( gcbob, 'G', 1. )
(  121)  
(  122)    END SUBROUTINE obc_spg
(  123) 
(  124) 
(  125)    SUBROUTINE obc_spg_east ( kt )
(  126)       !!------------------------------------------------------------------------------
(  127)       !!                ***  SUBROUTINE obc_spg_east  ***
(  128)       !!                 
(  129)       !! ** Purpose :   Apply the radiation algorithm on east OBC stream function.
(  130)       !!      If lfbceast=T , there is no radiation but only fixed OBC
(  131)       !!
(  132)       !!  History :
(  133)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  134)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  135)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  136)       !!         ! 00-06 (J.-M. Molines) 
(  137)       !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  138)       !!------------------------------------------------------------------------------
(  139)       !! * Arguments
(  140)       INTEGER, INTENT( in ) ::   kt
(  141) 
(  142)       !! * Local declarations
(  143)       INTEGER ::   ij
(  144)       REAL(wp) ::   z2dtr, ztau, zin
(  145)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  146)       !!------------------------------------------------------------------------------
(  147) 
(  148)       ! 1. First three time steps and more if lfbceast is .TRUE.
(  149)       !    In that case open boundary conditions are FIXED.
(  150)       ! --------------------------------------------------------
(  151) 
(  152)       IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart) .OR. lfbceast ) THEN
(  153) 
(  154)          ! 1.1 Fixed barotropic stream function
(  155)          ! ------------------------------------
(  156)          DO jj = nje0m1, nje1 
(  157)             ij = jj -1 + njmpp
(  158)             bsfeob(jj)=bfoe(ij)
(  159)          END DO
(  160) 
(  161)       ELSE
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 4

(  163)       ! 2. Beyond the fourth time step if lfbceast is .FALSE.
(  164)       ! -----------------------------------------------------
(  165) 
(  166)          ! 2.1. Barotropic stream function radiation
(  167)          ! ----------------------------------------
(  168)          !
(  169)          !          nibm2      nibm      nib
(  170)          !            |   nibm  |   nib   |///
(  171)          !            |    |    |    |    |///
(  172)          !  jj-line --f----v----f----v----f---
(  173)          !            |         |         |///
(  174)          !            |    |    |    |    |///
(  175)          !         jpieob-2   jpieob-1    jpieob
(  176)          !                 |         |        
(  177)          !              jpieob-1    jpieob      
(  178)          !
(  179)          ! ... radiative conditions plus restoring term toward climatology
(  180)          ! ... Initialize bsfeob to clim in any case, at the first step
(  181)          !     to ensure proper values at the ends of the open line.
(  182)          ! ... Take care that the j indices starts at nje0 (jpjed) and finish 
(  183)          !     at nje1 (jpjef) to be sure that jpjefm1 and jpjef are set OK.
(  184)          DO ji = fs_nie0-1, fs_nie1-1 ! Vector opt.
(  185)             DO jj = nje0p1, nje1m2 
(  186)                ij = jj -1 + njmpp
(  187)          ! ... 2* gradi(bsf) (v-point i=nibm, time mean)
(  188)                z2dx = ( bebnd(jj,nibm ,nit) + bebnd(jj,nibm ,nitm2) - 2.*bebnd(jj,nibm2,nitm) ) &
(  189)                       / e1v(ji,jj)
(  190)          ! ... 2* gradj(bsf) (f-point i=nibm, time nitm)
(  191)                z2dy = ( bebnd(jj+1,nibm,nitm) - bebnd(jj-1,nibm,nitm) ) / e2f(ji,jj)
(  192)          ! ... square of the norm of grad(bsf)
(  193)                z4nor2 = z2dx * z2dx + z2dy * z2dy
(  194)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  195)                zdt = bebnd(jj,nibm,nitm2) - bebnd(jj,nibm,nit)
(  196)          ! ... i-phase speed ratio (bounded by 1) and MASKED!
(  197)                IF( z4nor2 == 0 ) THEN
(  198)                   IF(lwp) WRITE(numout,*)' PB dans obc_spg_east au pt ',jj,' : z4nor=0'
(  199)                   z4nor2 = 0.001
(  200)                ENDIF
(  201)                z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
(  202)                z05cx = z05cx / e1v(ji+1,jj)
(  203)                z05cx = min( z05cx, 1. )
(  204)          ! ... z05cx < 0, inflow  zin=0, ztau=1  
(  205)          !          => 0, outflow zin=1, ztau=rtaue
(  206)                zin = sign( 1., z05cx )
(  207)                zin = 0.5*( zin + abs(zin) )
(  208)          ! ... Modification JM:  We maintain a restoring term toward
(  209)          !                   bsfb even in case of inflow 
(  210)          ! But restoring is stronger when in flow (10 days) (ztau in set in obcspg.F)
(  211)                ztau = (1.-zin ) * rtauein + zin * rtaue
(  212)                z05cx = z05cx * zin
(  213)          ! ... update bsfn with radiative or climatological bsf (not mask!)
(  214)                bsfeob(jj) = ( ( 1. - z05cx - ztau ) * bebnd(jj,nib ,nitm) + 2.*z05cx  &
(  215)                                * bebnd(jj,nibm,nit) + ztau * bfoe (ij) )              &
(  216)                             / (1. + z05cx)
(  217)             END DO
(  218)          END DO
(  219) 
(  220)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 5

(  221)       IF( lk_mpp )   CALL mppobc(bsfeob,jpjed,jpjef,jpieob-1,1,2,jpj)
(  222) 
(  223) 
(  224)       ! 3. right hand side of the barotropic elliptic equation
(  225)       ! ------------------------------------------------------
(  226)  
(  227)       IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
(  228)          z2dtr = 1.0 / rdt
(  229)       ELSE
(  230)          z2dtr = 0.5 / rdt
(  231)       ENDIF
(  232)       DO ji = fs_nie0-1, fs_nie1-1 ! Vector opt.
(  233)          DO jj = nje0m1, nje1 
(  234)             gcbob(ji,jj) = gcbob(ji,jj) - hvr(ji+1,jj) * e2v(ji+1,jj) / e1v(ji+1,jj) &
(  235)                            * ( bsfeob(jj) - bsfb(ji+1,jj) ) * z2dtr * bmask(ji,jj) 
(  236)          END DO
(  237)       END DO
(  238) 
(  239)    END SUBROUTINE obc_spg_east
(  240) 
(  241)    SUBROUTINE obc_spg_west ( kt )
(  242)       !!------------------------------------------------------------------------------
(  243)       !!                  ***  SUBROUTINE obc_spg_west  ***
(  244)       !!                    
(  245)       !! ** Purpose :
(  246)       !!      Apply the radiation algorithm on west OBC stream function.
(  247)       !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
(  248)       !!
(  249)       !!  History :
(  250)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  251)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  252)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  253)       !!         ! 00-06 (J.-M. Molines) 
(  254)       !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  255)       !!------------------------------------------------------------------------------
(  256)       !! * Arguments
(  257)       INTEGER, INTENT( in ) ::   kt
(  258) 
(  259)       !! * Local declarations
(  260)       INTEGER ::   ij
(  261) 
(  262)       REAL(wp) ::   z2dtr, ztau, zin
(  263)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  264) 
(  265)       !!------------------------------------------------------------------------------
(  266)       !!  OPA 8.5, LODYC-IPSL (2002)
(  267)       !!------------------------------------------------------------------------------
(  268) 
(  269)       ! 1. First three time steps and more if lfbcwest is .TRUE.
(  270)       !    In that case open boundary conditions are FIXED.
(  271)       ! --------------------------------------------------------
(  272) 
(  273)       IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcwest ) THEN
(  274) 
(  275)          ! 1.1 Fixed barotropic stream function
(  276)          ! ------------------------------------
(  277)          DO jj = njw0m1, njw1
(  278)             ij = jj -1 + njmpp






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 6

(  279)             bsfwob(jj)=bfow(ij)
(  280)          END DO
(  281) 
(  282)       ELSE
(  283) 
(  284)       ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
(  285)       ! -----------------------------------------------------
(  286) 
(  287)          ! 2.1. Barotropic stream function radiation
(  288)          ! ----------------------------------------
(  289)          !
(  290)          !         nib       nibm     nibm2
(  291)          !       ///|   nib   |   nibm  |
(  292)          !       ///|    |    |    |    |
(  293)          !       ---f----v----f----v----f-- jj-line
(  294)          !       ///|         |         |
(  295)          !       ///|    |    |    |    |
(  296)          !        jpiwob    jpiwob+1    jpiwob+2
(  297)          !               |         |        
(  298)          !             jpiwob+1    jpiwob+2     
(  299)          !
(  300)          ! ... radiative conditions plus restoring term toward climatology
(  301)          ! ... Initialize bsfwob to clim in any case, at the first step
(  302)          !     to ensure proper values at the ends of the open line.
(  303)          ! ... Take care that the j indices starts at njw0 (jpjwd) and finish
(  304)          ! ... at njw1 (jpjwf) to be sure that jpjwfm1 and jpjwf are set OK.
(  305)          DO ji = fs_niw0+1, fs_niw1+1 ! Vector opt.
(  306)             DO jj = njw0p1, njw1m2
(  307)                ij = jj -1 + njmpp
(  308)          ! ... 2* gradi(bsf) (v-point i=nibm, time mean)
(  309)                z2dx = ( -  bwbnd(jj,nibm ,nit ) - bwbnd(jj,nibm ,nitm2) + 2.*bwbnd(jj,nibm2,nitm ) ) &
(  310)                       / e1v(ji+1,jj)
(  311)          ! ... 2* gradj(bsf) (f-point i=nibm, time nitm)
(  312)                z2dy = ( bwbnd(jj+1,nibm,nitm) - bwbnd(jj-1,nibm,nitm) ) / e2f(ji,jj)
(  313)          ! ... square of the norm of grad(bsf)
(  314)                z4nor2 = z2dx * z2dx + z2dy * z2dy
(  315)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  316)                zdt = bwbnd(jj,nibm,nitm2) - bwbnd(jj,nibm,nit)
(  317)          ! ... i-phase speed ratio (bounded by 1) and MASKED!
(  318)                IF( z4nor2 == 0 ) THEN
(  319)                   IF(lwp) WRITE(numout,*)' PB dans obc_spg_west au pt ',jj,' : z4nor =0'
(  320)                   z4nor2=0.0001
(  321)                ENDIF
(  322)                z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
(  323)                z05cx = z05cx / e1v(ji,jj)
(  324)                z05cx = max( z05cx, -1. )
(  325)          ! ... z05cx => 0, inflow  zin=0, ztau=1  
(  326)          !           <  0, outflow zin=1, ztau=rtauw
(  327)                zin = sign( 1., -1. * z05cx )
(  328)                zin = 0.5*( zin + abs(zin) )
(  329)                ztau = (1.-zin )*rtauwin + zin * rtauw
(  330)                z05cx = z05cx * zin
(  331)          !  ... update bsfn with radiative or climatological bsf (not mask!)
(  332)                bsfwob(jj) = ( ( 1. + z05cx - ztau ) * bwbnd(jj,nib ,nitm) - 2.*z05cx &
(  333)                                * bwbnd(jj,nibm,nit) + ztau * bfow (ij) )             &
(  334)                             / (1. - z05cx)
(  335)             END DO
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 7

(  337) 
(  338)       ENDIF
(  339)       IF( lk_mpp )   CALL mppobc(bsfwob,jpjwd,jpjwf,jpiwob+1,1,2,jpj) 
(  340) 
(  341) 
(  342)       ! 3. right hand side of the barotropic elliptic equation
(  343)       ! -------------------------------------------------------
(  344) 
(  345)       IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
(  346)          z2dtr = 1.0 / rdt
(  347)       ELSE
(  348)          z2dtr = 0.5 / rdt
(  349)       ENDIF
(  350)       DO ji = fs_niw0+1, fs_niw1+1 ! Vector opt.
(  351)          DO jj = njw0m1, njw1
(  352)             gcbob(ji,jj) = gcbob(ji,jj) - hvr(ji,jj) * e2v(ji,jj) / e1v(ji,jj)     &
(  353)                            * ( bsfwob(jj) - bsfb(ji-1,jj) ) * z2dtr * bmask(ji,jj)
(  354)          END DO
(  355)       END DO
(  356) 
(  357)    END SUBROUTINE obc_spg_west
(  358) 
(  359)    SUBROUTINE obc_spg_north ( kt )
(  360)       !!------------------------------------------------------------------------------
(  361)       !!                 ***  SUBROUTINE obc_spg_north  ***
(  362)       !! 
(  363)       !! ** Purpose :   Apply the radiation algorithm on north OBC stream function.
(  364)       !!      If lfbcnorth=T, there is no radiation but only fixed OBC
(  365)       !!
(  366)       !!  History :
(  367)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  368)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  369)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  370)       !!         ! 00-06 (J.-M. Molines) 
(  371)       !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  372)       !!------------------------------------------------------------------------------
(  373)       !! * Arguments
(  374)       INTEGER, INTENT( in ) ::   kt
(  375) 
(  376)       !! * Local declarations
(  377)       INTEGER ::   ii
(  378)       REAL(wp) ::   z2dtr, ztau, zin
(  379)       REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
(  380)       !!------------------------------------------------------------------------------
(  381) 
(  382)       ! 1. First three time steps and more if lfbcnorth is .TRUE.
(  383)       !    In that case open boundary conditions are FIXED.
(  384)       ! --------------------------------------------------------
(  385) 
(  386)       IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth ) THEN
(  387) 
(  388)          ! 1.1 Fixed barotropic stream function
(  389)          ! ------------------------------------
(  390)          DO ji = nin0m1, nin1
(  391)             ii = ji -1 + nimpp
(  392)             bsfnob(ji)=bfon(ii)
(  393)          END DO
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 8

(  395)       ELSE      
(  396) 
(  397)       ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
(  398)       ! -----------------------------------------------------
(  399) 
(  400)          ! 2.1. Barotropic stream function radiation
(  401)          ! -----------------------------------------
(  402)          !
(  403)          !           ji-row
(  404)          !             |
(  405)          !        ////////////
(  406)          !        ////////////
(  407)          !   nib  -----f------  jpjnob
(  408)          !             |    
(  409)          !      nib--  u   ---- jpjnob
(  410)          !             |        
(  411)          !  nibm  -----f-----   jpjnob-1
(  412)          !             |        
(  413)          !     nibm--  u   ---- jpjnob-1
(  414)          !             |        
(  415)          !  nibm2 -----f-----   jpjnob-2
(  416)          !             |
(  417)          !
(  418)          ! ... radiative conditions plus restoring term toward climatology
(  419)          ! ... z05cx is always the cross boundary phase velocity
(  420)          ! ... Initialize bsfnob to clim in any case, at the first step
(  421)          !     to ensure proper values at the ends of the open line.
(  422)          ! ... Take care that the i indices starts at nin0 (jpind) and finish
(  423)          ! ... at nin1 (jpinf) to be sure that jpinfm1 and jpinf are set OK.
(  424)          DO jj = fs_njn0-1, fs_njn1-1 ! Vector opt.
(  425)             DO ji = nin0p1, nin1m2
(  426)                ii = ji -1 + nimpp
(  427)          ! ... 2* gradj(bsf) (u-point i=nibm, time mean)
(  428)                z2dx = ( bnbnd(ji,nibm ,nit) + bnbnd(ji,nibm ,nitm2) - 2.*bnbnd(ji,nibm2,nitm) ) &
(  429)                       / e2u(ji,jj)
(  430)          ! ... 2* gradi(bsf) (f-point i=nibm, time nitm)
(  431)                z2dy = ( bnbnd(ji+1,nibm,nitm) - bnbnd(ji-1,nibm,nitm) ) / e1f(ji,jj)
(  432)          ! ... square of the norm of grad(bsf)
(  433)                z4nor2 = z2dx * z2dx + z2dy * z2dy
(  434)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  435)                zdt = bnbnd(ji,nibm,nitm2) - bnbnd(ji,nibm,nit)
(  436)          ! ... j-phase speed ratio (bounded by 1) and MASKED!
(  437)                IF( z4nor2 == 0 ) THEN
(  438)                   IF(lwp) WRITE(numout,*)' PB dans obc_spg_north au pt',ji,' : z4nor =0'
(  439)                ENDIF
(  440)                z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
(  441)                z05cx = z05cx / e2u(ji,jj+1)
(  442)                z05cx = min( z05cx, 1. )
(  443)          ! ... z05cx < 0, inflow  zin=0, ztau=1 
(  444)          !          => 0, outflow zin=1, ztau=rtaun
(  445)                zin = sign( 1., z05cx )
(  446)                zin = 0.5*( zin + abs(zin) )
(  447)                ztau = (1.-zin ) * rtaunin + zin * rtaun
(  448)                z05cx = z05cx * zin
(  449)          ! ... update bsfn with radiative or climatological bsf (not mask!)
(  450)                bsfnob(ji) = ( ( 1. - z05cx - ztau ) * bnbnd(ji,nib ,nitm) + 2.*z05cx  &
(  451)                                * bnbnd(ji,nibm,nit) + ztau * bfon (ii) )              &
(  452)                             / (1. + z05cx)






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 9

(  453)             END DO
(  454)          END DO
(  455) 
(  456)       ENDIF
(  457)       IF( lk_mpp )   CALL mppobc(bsfnob,jpind,jpinf,jpjnob-1,1,1,jpi)
(  458) 
(  459) 
(  460)       ! 3. right hand side of the barotropic elliptic equation
(  461)       !-------------------------------------------------------
(  462) 
(  463)       IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
(  464)          z2dtr = 1.0 / rdt
(  465)       ELSE
(  466)          z2dtr = 0.5 / rdt
(  467)       ENDIF
(  468)       DO jj = fs_njn0-1, fs_njn1-1 ! Vector opt.
(  469)          DO ji = nin0m1, nin1
(  470)             gcbob(ji,jj) = gcbob(ji,jj) - hur(ji,jj+1) *  e1u(ji,jj+1) / e2u(ji,jj+1) &
(  471)                            * ( bsfnob(ji) - bsfb(ji,jj+1) ) * z2dtr * bmask(ji,jj)
(  472)          END DO
(  473)       END DO
(  474) 
(  475)    END SUBROUTINE obc_spg_north
(  476) 
(  477) 
(  478)    SUBROUTINE obc_spg_south ( kt )
(  479)       !!------------------------------------------------------------------------------
(  480)       !!                  ***  SUBROUTINE obc_spg_south  ***
(  481)       !!                
(  482)       !! ** Purpose :   Apply the radiation algorithm on south OBC stream function.
(  483)       !!      If lfbcsouth=T, there is no radiation but only fixed OBC
(  484)       !!
(  485)       !!  History :
(  486)       !!         ! 95-03 (J.-M. Molines) Original from SPEM
(  487)       !!         ! 97-07 (G. Madec, J.-M. Molines) additions
(  488)       !!         ! 97-12 (M. Imbard) Mpp adaptation
(  489)       !!         ! 00-06 (J.-M. Molines) 
(  490)       !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
(  491)       !!------------------------------------------------------------------------------
(  492)       !! * Arguments
(  493)       INTEGER, INTENT(in) ::   kt  ! ocean time-step index
(  494) 
(  495)       !! * Local declarations
(  496)       INTEGER  ::   ii             ! temporary integers
(  497)       REAL(wp) ::   &
(  498)          z2dtr, ztau, zin   ,   &  ! temporary scalars
(  499)          z05cx, zdt , z4nor2,   &  !    "         "
(  500)          z2dx , z2dy               !    "         "
(  501)       !!------------------------------------------------------------------------------
(  502) 
(  503)       ! 1. First three time steps and more if lfbcsouth is .TRUE.
(  504)       !    In that case open boundary conditions are FIXED.
(  505)       ! --------------------------------------------------------
(  506)       
(  507)       IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth ) THEN
(  508) 
(  509)          ! 1.1 Fixed barotropic stream function
(  510)          ! ------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 10

(  511)          DO ji = nis0m1, nis1 
(  512)             ii = ji -1 + nimpp 
(  513)             bsfsob(ji)=bfos(ii)
(  514)          END DO
(  515) 
(  516)       ELSE
(  517) 
(  518)       ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
(  519)       ! -----------------------------------------------------
(  520) 
(  521)          ! 2.1. Barotropic stream function radiation
(  522)          ! -----------------------------------------
(  523)          !
(  524)          !           ji-row
(  525)          !             |
(  526)          ! nibm2  -----f------  jpjsob + 2
(  527)          !             |    
(  528)          !   nibm  --  u  ----- jpjsob + 2
(  529)          !             |        
(  530)          !  nibm  -----f-----   jpjsob + 1
(  531)          !             |        
(  532)          !    nib  --  u  ----- jpjsob + 1
(  533)          !             |        
(  534)          !    nib -----f-----   jpjsob
(  535)          !        ///////////     
(  536)          !        ///////////
(  537)          !
(  538)          ! ... radiative conditions plus restoring term toward climatology
(  539)          ! ... z05cx is always the cross boundary phase velocity
(  540)          ! ... Initialize bsfsob to clim in any case, at the first step
(  541)          !     to ensure proper values at the ends of the open line.
(  542)          ! ... Take care that the i indices starts at nis0 (jpisd) and finish
(  543)          ! ... at nis1 (jpisf) to be sure that jpisfm1 and jpisf are set OK.
(  544)          DO jj = fs_njs0+1, fs_njs1+1 ! Vector opt.
(  545)             DO ji = nis0p1, nis1m2
(  546)                ii = ji -1 + nimpp
(  547)          ! ... 2* gradj(bsf) (u-point i=nibm, time mean)
(  548)                z2dx = ( - bsbnd(ji,nibm ,nit) - bsbnd(ji,nibm ,nitm2) + 2.*bsbnd(ji,nibm2,nitm) ) &
(  549)                       / e2u(ji,jj+1)
(  550)          ! ... 2* gradi(bsf) (f-point i=nibm, time nitm)
(  551)                z2dy = ( bsbnd(ji+1,nibm,nitm) - bsbnd(ji-1,nibm,nitm) ) / e1f(ji,jj)
(  552)          ! ... square of the norm of grad(bsf)
(  553)                z4nor2 = z2dx * z2dx + z2dy * z2dy
(  554)          ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
(  555)                zdt = bsbnd(ji,nibm,nitm2) - bsbnd(ji,nibm,nit)
(  556)          ! ... j-phase speed ratio (bounded by -1) and MASKED!
(  557)                IF( z4nor2 == 0 ) THEN
(  558)                   IF(lwp) WRITE(numout,*)' PB dans obc_spg_south au pt ',ji,' : z4nor =0'
(  559)                ENDIF
(  560)                z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
(  561)                z05cx = z05cx / e2u(ji,jj)
(  562)                z05cx = max( z05cx, -1. )
(  563)          ! ... z05cx => 0, inflow  zin=0, ztau=1
(  564)          !           <  0, outflow zin=1, ztau=rtaus
(  565)                zin = sign( 1., -1. * z05cx )
(  566)                zin = 0.5*( zin + abs(zin) )
(  567)                ztau = (1.-zin ) *rtausin  + zin * rtaus
(  568)                z05cx = z05cx * zin






PGF90 (Version     10.5)          09/23/2010  09:50:08      page 11

(  569)          ! ... update bsfn with radiative or climatological bsf (not mask!)
(  570)                bsfsob(ji) = ( ( 1. + z05cx - ztau ) * bsbnd(ji,nib ,nitm) - 2.*z05cx  & 
(  571)                                * bsbnd(ji,nibm,nit) + ztau * bfos (ii) )              &
(  572)                             / (1. - z05cx)
(  573)             END DO
(  574)          END DO
(  575) 
(  576)       ENDIF
(  577)       IF( lk_mpp )   CALL mppobc(bsfsob,jpisd,jpisf,jpjsob+1,1,1,jpi)
(  578) 
(  579)  
(  580)       ! 3. right hand side of the barotropic elliptic equation
(  581)       ! -------------------------------------------------------
(  582) 
(  583)       IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
(  584)          z2dtr = 1.0 / rdt
(  585)       ELSE
(  586)          z2dtr = 0.5 / rdt
(  587)       ENDIF
(  588)       DO jj = fs_njs0+1, fs_njs1+1 ! Vector opt.
(  589)          DO ji = nis0m1, nis1 
(  590)             gcbob(ji,jj) = gcbob(ji,jj) - hur(ji,jj) * e1u(ji,jj) / e2u(ji,jj) &
(  591)                            * ( bsfsob(ji) - bsfb(ji,jj-1) ) * z2dtr * bmask(ji,jj)
(  592)          END DO
(  593)       END DO
(  594) 
(  595)    END SUBROUTINE obc_spg_south
(  596) 
(  597) #else
(  598)    !!----------------------------------------------------------------------
(  599)    !!   Default case :                                         Empty module
(  600)    !!----------------------------------------------------------------------
(  601) CONTAINS
(  602)    SUBROUTINE obc_spg( kt )        ! Empty routine
(  603)       INTEGER, INTENT( in ) :: kt
(  604) !      WRITE(*,*) 'obc_spg: You should not have seen this print! error?', kt
(  605)    END SUBROUTINE obc_spg
(  606) #endif
(  607) 
(  608)    !!======================================================================
(  609) END MODULE obcspg























PGF90 (Version     10.5)          09/23/2010  09:50:08      page 12

(    1) # 1 "obcspg.F90"
(    1) MODULE obcspg 
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  obcspg  ***
(    4)    !! Open Boundaries  :   Radiation of barotropic stream function on each
(    5)    !!                      open boundary
(    6)    !!======================================================================
(    7) # 598
(  598)    !!----------------------------------------------------------------------
(  599)    !!   Default case :                                         Empty module
(  600)    !!----------------------------------------------------------------------
(  601) CONTAINS
(  602)    SUBROUTINE obc_spg( kt )        ! Empty routine
(  603)       INTEGER, INTENT( in ) :: kt
(  604) !      WRITE(*,*) 'obc_spg: You should not have seen this print! error?', kt
(  605)    END SUBROUTINE obc_spg
(  606) 
(  607) # 608
(  608)    !!======================================================================
(  609) END MODULE obcspg
