


PGF90 (Version     10.5)          09/23/2010  09:49:29      page 1

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

Filename: diaspr.F90

(    1) MODULE diaspr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diaspr  ***
(    4)    !! Ocean diagnostics:  surface pressure (rigid-lid case) 
(    5)    !!=====================================================================
(    6) #if   defined key_diaspr   &&   defined key_dynspg_rl
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_diaspr'        and                surface pressure diagnostics
(    9)    !!   'key_dynspg_rl'                                      rigid-lid case
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dia_spr      : update momentum and tracer Kz from a tke scheme
(   12)    !!   sprmat       : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE phycst          ! physical constants
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE sol_oce         ! ocean elliptic solver
(   20)    USE solpcg          ! preconditionned conjugate gradient solver
(   21)    USE solsor          ! Successive Over-relaxation solver
(   22)    USE solfet          ! FETI solver
(   23)    USE lib_mpp         ! distributed memory computing library
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC dia_spr   ! routine called by step.F90
(   30) 
(   31)    !! * Shared module variables
(   32)    LOGICAL, PUBLIC, PARAMETER ::   lk_diaspr = .TRUE.    !: surface pressure diag. flag
(   33)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   gps         !: surface pressure
(   34) 
(   35)    !! * Module variables
(   36)    INTEGER ::                 &
(   37)       nmoyps,                 &  ! time step for average
(   38)       nindic,                 &  ! indicator of convergence of the solver
(   39)       !                          ! namspr  surface pressure diagnostic
(   40)       nmaxp ,                 &  ! maximum of iterations for the solver
(   41)       niterp                     ! number of iteration done by the solver
(   42) 
(   43)    REAL(wp) ::     &
(   44)       ! namspr  surface pressure diagnostic
(   45)       epsp                       ! absolute precision of the solver
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 2

(   47)       !! * Namelist
(   48)       NAMELIST/namspr/ nmaxp, epsp, niterp
(   49) 
(   50)    REAL(wp) ::     &
(   51)       e1e2t                      ! ???
(   52) 
(   53)    REAL(wp), PUBLIC DIMENSION(jpi,jpj) ::   &
(   54)       spgum, spgvm,           &  ! average value of the surface pressure gradients
(   55)       gpsuu, gpsvv,           &  ! surface pressure gradients computed from comp. PS
(   56)       gcdpsc,                 &  ! inverse diagonal preconditioning matrix
(   57)       gcsmat,                 &  ! diagonal preconditioning matrix
(   58)       spmsk                      ! surface pressure Mask
(   59) 
(   60)    REAL(wp), DIMENSION(jpi,jpj,4) ::   &
(   61)       gcps                       ! extra-diagonal elements of SPG matrix
(   62)    !!----------------------------------------------------------------------
(   63)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   64)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diaspr.F90,v 1.4 2005/12/21 10:46:34 opalod Exp $ 
(   65)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68) CONTAINS
(   69) 
(   70)    SUBROUTINE dia_spr( kt )
(   71)       !!---------------------------------------------------------------------
(   72)       !!                  ***  ROUTINE dia_spr  ***
(   73)       !!
(   74)       !! ** Purpose :   compute the surface pressure from its gradient
(   75)       !!
(   76)       !! ** Method  :   rigid-lid appromimation: the surface pressure 
(   77)       !!      gradient is given by:
(   78)       !!           spgu = 1/rau0 d/dx(ps) = Mu + 1/(hu e2u) dj-1(bsfd)
(   79)       !!           spgv = 1/rau0 d/dy(ps) = Mv - 1/(hv e1v) di-1(bsfd)
(   80)       !!
(   81)       !!      where (Mu,Mv) is the vertically averaged momentum trend, i.e.
(   82)       !!      the vertical ponderated sum of the general momentum trend.
(   83)       !!      where bsfd is the trend of the barotropic stream function.
(   84)       !!
(   85)       !!       taking the divergence of the surface pressure gradient provides
(   86)       !!      an elliptic equation for ps which is solved using either a
(   87)       !!      diagonal preconditioned conjugate gradient method (solpcg.f) or
(   88)       !!      an successive-over-relaxation method (solsor.f) or FETI method
(   89)       !!      (solfet.F).
(   90)       !!
(   91)       !!      n.b. this resolution is valid with topography, cyclic east-west
(   92)       !!      boundary conditions and islands.
(   93)       !!
(   94)       !! History :
(   95)       !!        !  98-01  (G. Madec & M. Ioualalen)  Original code
(   96)       !!        !  98-02  (M. Guyon)  FETI method
(   97)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   98)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   99)       !!----------------------------------------------------------------------
(  100)       !! * Arguments
(  101)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  102) 
(  103)       !! * Local declarations
(  104)       INTEGER  ::  ji, jj






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 3

(  105)       INTEGER  ::   imax, ijt, iju
(  106)       REAL(wp) ::   zpsmea, zeps, zmoyr
(  107)       REAL(wp) ::   ztab(jpi,jpj,8)
(  108)       REAL(wp) ::   zemin1, zemax1, zemin2, zemax2, zgwgt
(  109)       REAL(wp) ::   z1, z2, zcompt,z3,z4
(  110)       REAL(wp) ::   zdif1, zdif2, zvar1, zvar2
(  111)       !!----------------------------------------------------------------------
(  112) 
(  113) 
(  114)       ! 0. initialisation (the first time step)
(  115)       ! ---------------------------------------
(  116)       
(  117)       IF( kt == nit000 ) THEN
(  118) 
(  119)          ! Namelist namspr : surface pressure
(  120) 
(  121)          nmaxp  = 2000
(  122)          epsp   = 1.e-6
(  123)          niterp = 16
(  124) 
(  125)          ! Read Namelist namspr : surface pressure diagnostics
(  126)          REWIND ( numnam )
(  127)          READ(numnam,namspr)
(  128) 
(  129)          IF(lwp) THEN
(  130)             WRITE(numout,*) 'dia_spr : surface pressure diagnostic (rigid-lid case)'
(  131)             WRITE(numout,*) '~~~~~~~'
(  132)             WRITE(numout,*)
(  133)             WRITE(numout,*) '          Namelist namspr : set solver parameters'
(  134)             WRITE(numout,*)
(  135)             WRITE(numout,*) '             maximum iterations for solver  nmaxp  = ', nmaxp
(  136)             WRITE(numout,*) '             absolute precision of solver   epsp   = ', epsp
(  137)             WRITE(numout,*) '             number of solver iterations    niterp = ', niterp
(  138)             WRITE(numout,*) '             frequeny of averaged output    nwrite = ', nwrite
(  139)             WRITE(numout,*)
(  140)          ENDIF
(  141) 
(  142)          ! control
(  143) # if ! defined key_dynspg_rl
(  144)       IF(lwp) WRITE(numout,cform_err)
(  145)       IF(lwp) WRITE(numout,*) '          surface pressure already explicitly computed !!'
(  146)       nstop = nstop + 1
(  147) # endif
(  148) 
(  149)          ! compute the ocean surface
(  150)          e1e2t = 0.e0
(  151)          DO jj = 2, jpjm1
(  152)             DO ji = 2, jpim1
(  153)                e1e2t = e1e2t + e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
(  154)             END DO
(  155)          END DO
(  156)          IF( lk_mpp )   CALL  mpp_sum( e1e2t )   ! sum over the global domain
(  157)          
(  158)          ! build the matrix for the surface pressure
(  159)          CALL sprmat
(  160)          
(  161)          ! set to zero the mean surface pressure gradient
(  162)          nmoyps = 0






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 4

(  163)          spgum(:,:) = 0.e0
(  164)          spgvm(:,:) = 0.e0
(  165) 
(  166)       ENDIF
(  167) 
(  168)       ! 1. cumulate the surface pressure gradient (at each time step)
(  169)       ! -----------------------------------------
(  170) 
(  171)       nmoyps = nmoyps + 1
(  172)       spgum(:,:) = spgum(:,:) + spgu(:,:)
(  173)       spgvm(:,:) = spgvm(:,:) + spgv(:,:)
(  174)       
(  175) 
(  176)       ! 2. ps computation each nwrite time step
(  177)       ! ---------------------------------------
(  178)       
(  179)       ! RETURN IF not the right time to compute ps
(  180)       IF ( MOD(kt-nit000+1,nwrite) /= 0 ) RETURN
(  181)       
(  182)       
(  183)       ! mean surface pressure gradient
(  184)       !   averaging and mask
(  185)       zmoyr = 1./float(nmoyps)
(  186)       DO jj = 2, jpjm1
(  187)          DO ji = 2, jpim1
(  188)             spgum(ji,jj) = spgum(ji,jj) * zmoyr * umask(ji,jj,1)
(  189)             spgvm(ji,jj) = spgvm(ji,jj) * zmoyr * vmask(ji,jj,1)
(  190)          END DO
(  191)       END DO
(  192) 
(  193)       CALL  lbc_lnk(spgum, 'U', -1. )
(  194)       CALL  lbc_lnk(spgvm, 'V', -1. )
(  195) 
(  196)       
(  197)       ! SAVE in local arrays and variables of solver informations
(  198)       zeps   = eps
(  199)       imax  = nmax 
(  200)       ztab(:,:,1) = gcp   (:,:,1)
(  201)       ztab(:,:,2) = gcp   (:,:,2)
(  202)       ztab(:,:,3) = gcp   (:,:,3)
(  203)       ztab(:,:,4) = gcp   (:,:,4)
(  204)       ztab(:,:,5) = gcdprc(:,:  )
(  205)       ztab(:,:,6) = gcdmat(:,:  )
(  206)       ztab(:,:,7) = gcx   (:,:  )
(  207)       ztab(:,:,8) = bmask (:,:  )
(  208) 
(  209)       ! replace bsf solver informations by ps solver one
(  210)       eps    = epsp
(  211)       nmax   = nmaxp
(  212)       gcp   (:,:,1) = gcps  (:,:,1)
(  213)       gcp   (:,:,2) = gcps  (:,:,2)
(  214)       gcp   (:,:,3) = gcps  (:,:,3)
(  215)       gcp   (:,:,4) = gcps  (:,:,4)
(  216)       gcdprc(:,:  ) = gcdpsc(:,:  )
(  217)       gcdmat(:,:  ) = gcsmat(:,:  )
(  218)       bmask (:,:  ) = spmsk (:,:  )
(  219)       !    first guess: ps
(  220)       gcx   (:,:  ) = gps   (:,:  )






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 5

(  221) 
(  222)       !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
(  223) 
(  224)       ! right hand side: 2d div. of the surface pressure gradient
(  225)       DO jj = 2, jpjm1
(  226)          DO ji = 2, jpim1
(  227)             gcb(ji,jj) = -gcdpsc(ji,jj)*   &
(  228)                (  e2u(ji,jj)*spgum(ji,jj) - e2u(ji-1,jj)*spgum(ji-1,jj)   &
(  229)                + e1v(ji,jj)*spgvm(ji,jj) - e1v(ji,jj-1)*spgvm(ji,jj-1) )
(  230)          END DO
(  231)       END DO
(  232)       
(  233)       !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
(  234)       
(  235)       ! relative PRECISION
(  236)       rnorme = 0.
(  237)       DO jj = 1, jpj
(  238)          DO ji = 1, jpi
(  239)             rnorme = rnorme + gcb(ji,jj) * gcsmat(ji,jj) * gcb(ji,jj)
(  240)          END DO
(  241)       END DO
(  242)       IF( lk_mpp )   CALL  mpp_sum( rnorme )   ! sum over the global domain
(  243) 
(  244)       epsr=eps*eps*rnorme
(  245)       ncut=0
(  246)       !   IF the second member is 0 the solution is 0, solpcg isn't called
(  247)       IF ( rnorme == 0.e0 ) THEN
(  248)          gps(:,:) = 0.e0
(  249)          res   = 0.e0
(  250)          niter = 0
(  251)          ncut  = 999
(  252)       ENDIF
(  253)       
(  254)       !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
(  255)       
(  256)       nindic = 0
(  257) 
(  258)       ! iterarive solver of the spg system (except IF sol.=0)
(  259)       !     (OUTPUT in gcx with boundary conditions applied)
(  260)       IF ( ncut == 0 ) THEN
(  261)          IF ( nsolv == 1 ) THEN
(  262)             CALL sol_pcg( nindic )         !   diagonal preconditioned conjuguate gradient
(  263)          ELSE IF ( nsolv == 2 ) THEN
(  264)             CALL sol_sor( nindic )     !   successive-over-relaxation
(  265)          ELSE IF(nsolv == 3) THEN
(  266)             CALL sol_fet( nindic )         !   FETI solver
(  267)          ELSE
(  268)             !   e r r o r  in nsolv namelist PARAMETER
(  269)             IF(lwp) THEN
(  270)                WRITE(numout,*) ' dia_spr: e r r o r, nsolv = 1 or 2'
(  271)                WRITE(numout,*) ' ******               not = ',nsolv
(  272)             ENDIF
(  273)             STOP 'dia_spr'
(  274)          ENDIF
(  275)       ENDIF
(  276)       
(  277)       !,,,,,,,,,,,,,,,,,,,,,,,,synchro IF macrotasking,,,,,,,,,,,,,,,,,,,,,,,
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 6

(  279)       
(  280)       ! sp solver statistics  (i.e. problem for the solver) 
(  281)       IF ( epsr < 0.) THEN
(  282)          IF(lwp) THEN 
(  283)             WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  284)             IF(lwp)WRITE(numout,*)'dia_spr-1:',epsr
(  285)             IF(lwp)WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  286)          ENDIF
(  287)       ENDIF
(  288)       IF(lwp)WRITE(numout,9300) kt, niter, res, SQRT(epsr)/eps
(  289)       IF (nindic < 0) THEN 
(  290)          IF(lwp) THEN 
(  291)             WRITE(numout,9100)
(  292)             WRITE(numout,*) ' dia_spr : the surface pressure solver DO not converge'
(  293)             WRITE(numout,*) ' ====== ' 
(  294)             WRITE(numout,*) 
(  295)          ENDIF
(  296)       ENDIF
(  297) 9100  FORMAT( /,' ===>>>> : w a r n i n g',/,'          ===============',/ )
(  298) 9300  FORMAT(' it :', i8, ' niter :', i4, ' res :',e20.10,' b :', e20.10)
(  299) 
(  300)       ! recover bsf solver informations and SAVE ps for next computation
(  301)       eps    = zeps
(  302)       nmax   = imax 
(  303)       gps   (:,:  ) = gcx (:,:)
(  304)       gcp   (:,:,1) = ztab(:,:,1)
(  305)       gcp   (:,:,2) = ztab(:,:,2)
(  306)       gcp   (:,:,3) = ztab(:,:,3)
(  307)       gcp   (:,:,4) = ztab(:,:,4)
(  308)       gcdprc(:,:  ) = ztab(:,:,5)
(  309)       gcdmat(:,:  ) = ztab(:,:,6)
(  310)       gcx   (:,:  ) = ztab(:,:,7)
(  311)       bmask (:,:  ) = ztab(:,:,8)
(  312)       
(  313)       ! compute and substract the mean value
(  314)       
(  315)       zpsmea = 0.e0
(  316)       DO jj=2,jpjm1
(  317)          DO ji=2,jpim1
(  318)             zpsmea = zpsmea + gps(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask(ji,jj,1)
(  319)          END DO
(  320)       END DO
(  321)       IF( lk_mpp )   CALL  mpp_sum( zpsmea )   ! sum over the global domain
(  322) 
(  323)       zpsmea = zpsmea / e1e2t
(  324)       gps(:,:) = ( gps(:,:) - zpsmea ) * tmask(:,:,1)
(  325)  
(  326)       IF(lwp)WRITE(numout,*) ' mean value of ps = ',zpsmea,' is substracted'
(  327)       ! ----------------------------------------
(  328)       ! i. compute the surface pressure gradient
(  329)       !    from the computed surface pressure
(  330)       ! ----------------------------------------
(  331) 
(  332)       DO jj=2,jpjm1
(  333)          DO ji=2,jpim1
(  334)             gpsuu(ji,jj)=(gps(ji+1,jj)-gps(ji,jj))/e1u(ji,jj) * umask(ji,jj,1)
(  335)             gpsvv(ji,jj)=(gps(ji,jj+1)-gps(ji,jj))/e2v(ji,jj) * vmask(ji,jj,1)
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 7

(  337)       END DO
(  338)       
(  339)       ! compute the max and min error
(  340)       
(  341)       zemax1 = 0.e0
(  342)       zemin1 = 0.e0
(  343)       zemax2 = 0.e0
(  344)       zemin2 = 0.e0
(  345)       DO jj = 2,jpj-1
(  346)          DO ji = 2,jpi-1
(  347)             z1 = ABS( spgum(ji,jj)-gpsuu(ji,jj) )*umask(ji,jj,1)
(  348)             z2 = ABS( spgvm(ji,jj)-gpsvv(ji,jj) )*vmask(ji,jj,1)
(  349)             z3 = MAX ( ABS( spgum(ji,jj) ), ABS( spgvm(ji,jj) ) )
(  350)             z4 = MAX ( ABS( gpsuu(ji,jj) ), ABS( gpsvv(ji,jj) ) )
(  351)             zemax1 = MAX(z1,zemax1)
(  352)             zemax2 = MAX(z2,zemax2)
(  353)             zemin1 = MAX(z3,zemin1)
(  354)             zemin2 = MAX(z4,zemin2)
(  355)          END DO
(  356)       END DO
(  357)       IF( lk_mpp )   CALL  mpp_sum( zemax1 )   ! sum over the global domain
(  358)       IF( lk_mpp )   CALL  mpp_sum( zemax2 )   ! sum over the global domain
(  359)       IF( lk_mpp )   CALL  mpp_sum( zemin1 )   ! sum over the global domain
(  360)       IF( lk_mpp )   CALL  mpp_sum( zemin2 )   ! sum over the global domain
(  361) 
(  362)       IF(lwp) THEN
(  363)          WRITE(numout,*)
(  364)          WRITE(numout,*) 'pserro : time step = ', kt
(  365)          WRITE(numout,*) '******** ------------------'
(  366)          WRITE(numout,*)
(  367)          WRITE(numout,*) '         gpsx error  max=',zemax1
(  368)          WRITE(numout,*) '         gpsy error  max=',zemax2
(  369)          WRITE(numout,*) '         gps max =',zemin1
(  370)          WRITE(numout,*) '         gpsc max =',zemin2
(  371)          WRITE(numout,*)
(  372)       ENDIF
(  373) 
(  374)       ! compute the norme and variance of this error
(  375) 
(  376)       zcompt = 0.e0
(  377)       zdif1  = 0.e0
(  378)       zdif2  = 0.e0
(  379)       zvar1  = 0.e0
(  380)       zvar2  = 0.e0
(  381)       DO jj = 2, jpj-1
(  382)          DO ji = 2, jpi-1
(  383)             z1 = ( spgum(ji,jj)-gpsuu(ji,jj) ) * umask(ji,jj,1)
(  384)             z2 = ( spgvm(ji,jj)-gpsvv(ji,jj) ) * vmask(ji,jj,1)
(  385)             zcompt=zcompt+tmask(ji,jj,1)
(  386)             zdif1=zdif1+z1
(  387)             zdif2=zdif2+z2
(  388)             zvar1=zvar1+z1*z1
(  389)             zvar2=zvar2+z2*z2
(  390)          END DO
(  391)       END DO
(  392)       IF( lk_mpp )   CALL  mpp_sum( zcompt )   ! sum over the global domain
(  393)       IF( lk_mpp )   CALL  mpp_sum( zdif1  )   ! sum over the global domain
(  394)       IF( lk_mpp )   CALL  mpp_sum( zdif2  )   ! sum over the global domain






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 8

(  395)       IF( lk_mpp )   CALL  mpp_sum( zvar1  )   ! sum over the global domain
(  396)       IF( lk_mpp )   CALL  mpp_sum( zvar2  )   ! sum over the global domain
(  397) 
(  398)       IF(lwp) WRITE(numout,*) '        zcompt = ',zcompt
(  399)       zdif1=zdif1/zcompt
(  400)       zdif2=zdif2/zcompt
(  401)       IF( zvar1 < 0.) THEN 
(  402)          IF(lwp) THEN
(  403)             WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  404)             WRITE(numout,*)'dia_spr-2:',zvar1
(  405)             WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  406)          ENDIF
(  407)       ENDIF
(  408)       zvar1 = SQRT(zvar1)/zcompt
(  409)       IF( zvar2 < 0. ) THEN 
(  410)          IF(lwp)THEN
(  411)             WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  412)             WRITE(numout,*)'dia_spr-3:',zvar2
(  413)             WRITE(numout,*)'rrrrrrrrrrrrrrr'
(  414)          ENDIF
(  415)       ENDIF
(  416)       zvar2 = SQRT(zvar2)/zcompt
(  417)       
(  418)       IF(lwp) THEN 
(  419)          WRITE(numout,*)
(  420)          WRITE(numout,*) '         gpsx mean error = ',zdif1
(  421)          WRITE(numout,*) '         gpsy mean error = ',zdif2
(  422)          WRITE(numout,*)
(  423)          WRITE(numout,*) '         gpsx var. error = ',zvar1
(  424)          WRITE(numout,*) '         gpsy var. error = ',zvar2
(  425)          WRITE(numout,*)
(  426)          WRITE(numout,*)
(  427)       ENDIF
(  428)       
(  429)       ! reset to zero nmoyps and the mean surface pressure gradient
(  430)       nmoyps = 0
(  431)       spgum(:,:) = 0.e0
(  432)       spgvm(:,:) = 0.e0
(  433)       
(  434)    END SUBROUTINE dia_spr
(  435) 
(  436) 
(  437)    SUBROUTINE sprmat
(  438)       !!---------------------------------------------------------------------
(  439)       !!                  ***  ROUTINE sprmat  ***
(  440)       !!               
(  441)       !! ** Purpose :   construction of the matrix of the surface pressure
(  442)       !!      system and the diagonal preconditioning matrix.
(  443)       !!
(  444)       !! ** Method :
(  445)       !!
(  446)       !! History :
(  447)       !!        !  98-01  (G. Madec & M. Ioualalen)  Original code
(  448)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(  449)       !!----------------------------------------------------------------------
(  450)       !! * Local declarations
(  451)       INTEGER ::   ji, jj, jl                ! dummy loop indices
(  452)       REAL(wp) ::   zcoefs, zcoefw, zcoefe, zcoefn






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 9

(  453)       !!----------------------------------------------------------------------
(  454)       
(  455)       
(  456)       ! 0. ocean/land mask at ps-point (computed from surface tmask): spmsk
(  457)       ! --------------------------------------------------------------------
(  458)       
(  459)       ! computation
(  460)       spmsk(:,:) = tmask(:,:,1)
(  461)       
(  462)       ! boundary conditions
(  463)       ! south symmetry: psmsk must be set to 0. on 1
(  464)       IF( nperio == 2 ) THEN
(  465)          spmsk(:, 1 ) = 0.e0
(  466)       ENDIF
(  467)       
(  468)       ! east-west cyclic: spmsk must be set to 0. on 1 and jpi
(  469)       IF( nperio == 1 .OR. nperio == 4 .OR.nperio == 6) THEN 
(  470)          spmsk( 1 ,:) = 0.e0
(  471)          spmsk(jpi,:) = 0.e0
(  472)       ENDIF
(  473)       
(  474)       ! north fold: spmsk must be set to 0. on ligne jpj and on half
(  475)       !                   ligne jpj-1
(  476)       ! T-point pivot
(  477)       IF( nperio == 3 .OR. nperio == 4 ) THEN
(  478)          spmsk(:,jpj) = 0.e0
(  479)          DO ji = jpi/2+1, jpi
(  480)             spmsk(ji,jpjm1) = 0.e0
(  481)          END DO
(  482)       ENDIF
(  483)       ! F-point pivot
(  484)       IF( nperio == 5 .OR. nperio == 6 ) THEN
(  485)          spmsk(:,jpj) = 0.e0
(  486)       ENDIF
(  487)       
(  488)       ! mpp boundary cond.: spmsk is initialized at zero on the overlap
(  489)       ! region for both the preconjugate gradient and the sor algorithms
(  490) 
(  491)       IF( nbondi /= -1 .AND. nbondi /= 2 ) THEN
(  492)          DO jl = 1, jpreci
(  493)             spmsk(jl,:) = 0.e0
(  494)          END DO
(  495)       ENDIF
(  496)       IF( nbondi /= 1 .AND. nbondi /= 2 ) THEN
(  497)          DO ji = nlci, jpi
(  498)             spmsk(ji,:) = 0.e0
(  499)          END DO
(  500)       ENDIF
(  501)       IF( nbondj /= -1 .AND. nbondj /= 2 ) THEN
(  502)          DO jl=1,jprecj
(  503)             spmsk(:,jl) = 0.e0
(  504)          END DO
(  505)       ENDIF
(  506)       IF( nbondj /= 1 .AND. nbondj /= 2 ) THEN
(  507)          DO jj = nlcj, jpj
(  508)             spmsk(:,jj) = 0.e0
(  509)          END DO
(  510)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 10

(  511)       
(  512)       ! 1. construction of the matrix
(  513)       ! -----------------------------
(  514)       
(  515)       DO jj = 1, jpj
(  516)          DO ji = 1, jpi
(  517)             
(  518)             IF( spmsk(ji,jj) == 0. ) THEN
(  519)                ! land points
(  520)                gcps  (ji,jj,1) = 0.e0
(  521)                gcps  (ji,jj,2) = 0.e0
(  522)                gcps  (ji,jj,3) = 0.e0
(  523)                gcps  (ji,jj,4) = 0.e0
(  524)                gcdpsc(ji,jj  ) = 0.e0
(  525)                gcsmat(ji,jj  ) = 0.e0
(  526)             ELSE
(  527)                ! south coefficient
(  528)                zcoefs = -e1v(ji,jj-1) / e2v(ji,jj-1) * vmask(ji,jj-1,1)
(  529)                gcps(ji,jj,1) = zcoefs
(  530)                ! west coefficient
(  531)                zcoefw = -e2u(ji-1,jj) / e1u(ji-1,jj) * umask(ji-1,jj,1)
(  532)                gcps(ji,jj,2) = zcoefw
(  533)                ! east coefficient
(  534)                zcoefe = -e2u(ji  ,jj) / e1u(ji  ,jj) * umask(ji  ,jj,1)
(  535)                gcps(ji,jj,3) = zcoefe
(  536)                ! north coefficient
(  537)                zcoefn = -e1v(ji,jj  ) / e2v(ji,jj  ) * vmask(ji,jj  ,1)
(  538)                gcps(ji,jj,4) = zcoefn
(  539)                
(  540)                ! diagonal coefficient
(  541)                gcsmat(ji,jj) = -zcoefs-zcoefw-zcoefe-zcoefn
(  542)             ENDIF
(  543)          END DO
(  544)       END DO
(  545)       
(  546) 
(  547)       ! 2. boundary conditions 
(  548)       ! ----------------------
(  549)       
(  550)       ! cyclic east-west boundary conditions
(  551)       ! ji=2 is the column east of ji=jpim1 and reciprocally,
(  552)       ! ji=jpim1 is the column west of ji=2
(  553)       ! all the coef are already set to zero as spmask is initialized to
(  554)       ! zero for ji=1 and ji=jpj.
(  555)       
(  556)       ! symetrical conditions
(  557)       ! the diagonal coefficient of the southern grid points must be modify to
(  558)       ! account for the existence of the south symmetric bassin.
(  559)       IF( nperio == 2 ) THEN
(  560)          DO ji = 1, jpi
(  561)             IF( spmsk(ji,2) /= 0 ) THEN
(  562)                zcoefs = e1v(ji,1) / e2v(ji,1)
(  563)                gcsmat(ji,2) = gcsmat(ji,2) - zcoefs
(  564)             ENDIF
(  565)          END DO
(  566)       ENDIF
(  567)       
(  568)       ! North fold boundary condition






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 11

(  569)       ! all the coef are already set to zero as bmask is initialized to
(  570)       ! zero on duplicated lignes and portion of lignes
(  571)       
(  572)       
(  573)       ! 3. preconditioned matrix
(  574)       ! ------------------------
(  575)       
(  576)       DO jj = 1, jpj
(  577)          DO ji = 1, jpi
(  578)             IF( spmsk(ji,jj) /= 0. ) gcdpsc(ji,jj) = 1.e0 / gcsmat(ji,jj)
(  579)          END DO
(  580)       END DO
(  581)       
(  582)       gcps(:,:,1) = gcps(:,:,1) * gcdpsc(:,:)
(  583)       gcps(:,:,2) = gcps(:,:,2) * gcdpsc(:,:)
(  584)       gcps(:,:,3) = gcps(:,:,3) * gcdpsc(:,:)
(  585)       gcps(:,:,4) = gcps(:,:,4) * gcdpsc(:,:)
(  586)       
(  587)       
(  588)       ! 3. initialization the arrays used in sp solver
(  589)       ! ----------------------------------------------
(  590)       
(  591)       gps  (:,:) = 0.e0
(  592)       gpsuu(:,:) = 0.e0
(  593)       gpsvv(:,:) = 0.e0
(  594)       
(  595)    END SUBROUTINE sprmat
(  596) 
(  597) #else
(  598)    !!----------------------------------------------------------------------
(  599)    !!   Default option :                    NO surface pressure diagnostics
(  600)    !!----------------------------------------------------------------------
(  601)    USE in_out_manager  
(  602)    LOGICAL, PUBLIC, PARAMETER ::   lk_diaspr = .FALSE.   !: surface pressure diag. flag
(  603) CONTAINS
(  604)    SUBROUTINE dia_spr( kt )      ! Empty routine
(  605)       if(lwp) WRITE(numout,*) 'dia_spr: You should not have seen this print! error?', kt
(  606)    END SUBROUTINE dia_spr
(  607) #endif
(  608) 
(  609)    !!======================================================================
(  610) END MODULE diaspr






















PGF90 (Version     10.5)          09/23/2010  09:49:29      page 12

(    1) # 1 "diaspr.F90"
(    1) MODULE diaspr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diaspr  ***
(    4)    !! Ocean diagnostics:  surface pressure (rigid-lid case) 
(    5)    !!=====================================================================
(    6) # 598
(  598)    !!----------------------------------------------------------------------
(  599)    !!   Default option :                    NO surface pressure diagnostics
(  600)    !!----------------------------------------------------------------------
(  601)    USE in_out_manager  
(  602)    LOGICAL, PUBLIC, PARAMETER ::   lk_diaspr = .FALSE.   !: surface pressure diag. flag
(  603) CONTAINS
(  604)    SUBROUTINE dia_spr( kt )      ! Empty routine
(  605)       if(lwp) WRITE(numout,*) 'dia_spr: You should not have seen this print! error?', kt
(  606)    END SUBROUTINE dia_spr
(  607) 
(  608) # 609
(  609)    !!======================================================================
(  610) END MODULE diaspr
