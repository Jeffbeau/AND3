


PGF90 (Version     10.5)          09/23/2010  09:49:10      page 1

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

Filename: trdmld.F90

(    1) MODULE trdmld
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdmld  ***
(    4)    !! Ocean diagnostics:  mixed layer T-S trends 
(    5)    !!=====================================================================
(    6) #if   defined key_trdmld   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_trdmld'                          mixed layer trend diagnostics
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trd_mld          : T and S cumulated trends averaged over the mixed layer
(   11)    !!   trd_mld_zint     : T and S trends vertical integration
(   12)    !!   trd_mld_init     : initialization step
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE ldftra_oce      ! ocean active tracers lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE phycst          ! Define parameters for the routines
(   22)    USE daymod          ! calendar
(   23)    USE dianam          ! build the name of file (routine)
(   24)    USE ldfslp          ! iso-neutral slopes 
(   25)    USE zdfmxl          ! mixed layer depth
(   26)    USE zdfddm          ! ocean vertical physics: double diffusion
(   27)    USE ioipsl          ! NetCDF library
(   28)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   29) !   USE diadimg         ! dimg direct access file format output
(   30) 
(   31)    IMPLICIT NONE
(   32)    PRIVATE
(   33) 
(   34)    !! * Accessibility
(   35)    PUBLIC trd_mld        ! routine called by step.F90
(   36)    PUBLIC trd_mld_init   ! routine called by opa.F90
(   37)    PUBLIC trd_mld_zint   ! routine called by tracers routines
(   38) 
(   39)    !! * Shared module variables
(   40)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .TRUE.    !: momentum trend flag
(   41) 
(   42)    !! * Module variables
(   43)    INTEGER ::   &
(   44)       nh_t, nmoymltrd,             &  ! ???
(   45)       nidtrd,                      &
(   46)       ndextrd1(jpi*jpj),           &






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 2

(   47)       ndimtrd1
(   48)    INTEGER, SAVE ::   &
(   49)       ionce, icount,               &
(   50)       idebug                          ! (0/1) set it to 1 in case of problem to have more print
(   51) 
(   52)    INTEGER, DIMENSION(jpi,jpj) ::   &
(   53)       nmld,                         & ! mixed layer depth
(   54)       nbol                
(   55) 
(   56)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   57)       rmld   ,          &  ! mld depth (m) corresponding to nmld
(   58)       tml    , sml  ,   &  ! average T and S over mixed layer
(   59)       tmlb   , smlb ,   &  ! before tml and sml (kt-1)
(   60)       tmlbb  , smlbb,   &  ! tml and sml at begining of the nwrite-1 
(   61)       !                    ! timestep averaging period
(   62)       tmlbn  , smlbn,   &  ! after tml and sml at time step after the
(   63)       !                    ! begining of the NWRITE-1 timesteps
(   64)       tmltrdm, smltrdm     !
(   65) 
(   66)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   67)       tmltrd ,          &  ! total cumulative trends of temperature and 
(   68)       smltrd ,          &  ! salinity over nwrite-1 time steps
(   69)       wkx
(   70) 
(   71)    CHARACTER(LEN=80) :: clname
(   72) 
(   73)    !! * Substitutions
(   74) #  include "domzgr_substitute.h90"
(   75) #  include "ldftra_substitute.h90"
(   76) #  include "zdfddm_substitute.h90"
(   77)    !!----------------------------------------------------------------------
(   78)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   79)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmld.F90,v 1.7 2005/12/12 14:18:08 opalod Exp $ 
(   80)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   81)    !!----------------------------------------------------------------------
(   82) 
(   83) CONTAINS
(   84) 
(   85) SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
(   86)       !!----------------------------------------------------------------------
(   87)       !!                  ***  ROUTINE trd_mld_zint  ***
(   88)       !! 
(   89)       !! ** Purpose :   computation of vertically integrated T and S budgets
(   90)       !!                from ocean surface down to control surface 
(   91)       !!
(   92)       !! ** Method/usage :
(   93)       !!      integration done over nwrite-1 time steps 
(   94)       !!      Control surface can be either a mixed layer depth (time varying)
(   95)       !!      or a fixed surface (jk level or bowl). 
(   96)       !!      Choose control surface with nctls in namelist NAMDIA.
(   97)       !!      nctls = 0  : use mixed layer with density criterion 
(   98)       !!      nctls = 1  : read index from file 'ctlsurf_idx'
(   99)       !!      nctls > 1  : use fixed level surface jk = nctls
(  100)       !!      Note: in the remainder of the routine, the volume between the 
(  101)       !!            surface and the control surface is called "mixed-layer"
(  102)       !!      Method check : if the control surface is fixed, the residual dh/dt
(  103)       !!                     entrainment should be zero
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 3

(  105)       !! ** Action :
(  106)       !!            /commld/   : rmld         mld depth corresponding to nmld
(  107)       !!                         tml          average T over mixed layer
(  108)       !!                         tmlb         tml at kt-1
(  109)       !!                         tmlbb        tml at begining of the NWRITE-1 
(  110)       !!                                      time steps averaging period
(  111)       !!                         tmlbn        tml at time step after the 
(  112)       !!                                      begining of the NWRITE-1 time
(  113)       !!                                      steps averaging period
(  114)       !!
(  115)       !!                  mixed layer trends :
(  116)       !!
(  117)       !!                  tmltrd (,,1) = zonal advection
(  118)       !!                  tmltrd (,,2) = meridional advection
(  119)       !!                  tmltrd (,,3) = vertical advection
(  120)       !!                  tmltrd (,,4) = lateral diffusion (horiz. component+Beckman)
(  121)       !!                  tmltrd (,,5) = forcing
(  122)       !!                  tmltrd (,,6) = entrainment due to vertical diffusion (TKE)
(  123)       !!          if iso  tmltrd (,,7) = lateral diffusion (vertical component)
(  124)       !!                  tmltrd (,,8) = eddy induced zonal advection
(  125)       !!                  tmltrd (,,9) = eddy induced meridional advection
(  126)       !!                  tmltrd (,,10) = eddy induced vertical advection
(  127)       !!
(  128)       !!                  tmltrdm(,) : total cumulative trends over nwrite-1 time steps
(  129)       !!                  ztmltot(,) : dT/dt over the NWRITE-1 time steps 
(  130)       !!                               averaging period (including Asselin 
(  131)       !!                               terms)
(  132)       !!                  ztmlres(,) : residual = dh/dt entrainment
(  133)       !!
(  134)       !!      trends output in netCDF format using ioipsl
(  135)       !!
(  136)       !! History :
(  137)       !!        !  95-04  (J. Vialard)  Original code
(  138)       !!        !  97-02  (E. Guilyardi)  Adaptation global + base cmo
(  139)       !!        !  99-09  (E. Guilyardi)  Re-writing + netCDF output
(  140)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  141)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  142)       !!----------------------------------------------------------------------
(  143)       !! * Arguments
(  144)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  145) 
(  146)       CHARACTER(len=2), INTENT( in ) ::   &
(  147)          ctype                                ! surface/bottom (2D arrays) or
(  148)                                               ! interior (3D arrays) physics
(  149) 
(  150)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  151)          pttrdmld,                         &  ! Temperature trend 
(  152)          pstrdmld                             ! Salinity    trend
(  153) 
(  154)       !! * Local declarations
(  155)       INTEGER ::   ji, jj, jk, isum
(  156) # if defined key_trabbl_dif
(  157)       INTEGER ::   ikb
(  158) # endif
(  159) 
(  160)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  161)          zvlmsk
(  162)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 4

(  163) 
(  164)       IF( icount == 1 ) THEN        
(  165) 
(  166)          zvlmsk(:,:)   = 0.e0
(  167)          tmltrd(:,:,:) = 0.e0
(  168)          smltrd(:,:,:) = 0.e0
(  169)          
(  170)          ! This computation should be done only once per time step
(  171) 
(  172)          !  ========================================================
(  173)          !   I. definition of control surface and associated fields
(  174)          !  ========================================================
(  175) 
(  176)          !    I.1 set nmld(ji,jj) = index of first T point below control surface
(  177)          !    -------------------                       or outside mixed-layer
(  178) 
(  179)          IF( nctls == 0 ) THEN
(  180)             ! control surface = mixed-layer with density criterion 
(  181)             ! (array nmln computed in zdfmxl.F90)
(  182)             nmld(:,:) = nmln(:,:)
(  183)          ELSE IF( nctls == 1 ) THEN
(  184)             ! control surface = read index from file 
(  185)             nmld(:,:) = nbol(:,:)
(  186)          ELSE IF( nctls >= 2 ) THEN
(  187)             ! control surface = model level
(  188)             nctls = MIN( nctls, jpktrd - 1 )
(  189)             nmld(:,:) = nctls + 1
(  190)          ENDIF
(  191) 
(  192)          IF( ionce == 1 ) THEN  ! compute ndextrd1 and ndimtrd1 only once
(  193)             ! Check of validity : nmld(ji,jj) =< jpktrd
(  194)             isum = 0
(  195) 
(  196)             IF( jpktrd < jpk ) THEN 
(  197)                DO jj = 1, jpj
(  198)                   DO ji = 1, jpi
(  199)                      IF( nmld(ji,jj) <= jpktrd ) THEN
(  200)                         zvlmsk(ji,jj) = tmask(ji,jj,1)
(  201)                      ELSE
(  202)                         isum = isum + 1
(  203)                         zvlmsk(ji,jj) = 0.
(  204)                      ENDIF
(  205)                   END DO
(  206)                END DO
(  207)             ENDIF
(  208) 
(  209)             ! Index of ocean points (2D only)
(  210)             IF( isum > 0 ) THEN
(  211)                if(lwp) WRITE(numout,*)' Number of invalid points nmld > jpktrd', isum 
(  212)                CALL wheneq( jpi*jpj, zvlmsk(:,:) , 1, 1., ndextrd1, ndimtrd1 )    ! surface
(  213)             ELSE 
(  214)                CALL wheneq( jpi*jpj, tmask(:,:,1), 1, 1., ndextrd1, ndimtrd1 )    ! surface
(  215)             ENDIF                                
(  216) 
(  217)             ! no more pass here
(  218)             ionce = 0
(  219) 
(  220)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 5

(  221)          
(  222)          IF( idebug /= 0 ) THEN
(  223)             ! CALL prihre (zvlmsk,jpi,jpj,1,jpi,2,1,jpj,2,3,numout)
(  224)             WRITE(numout,*) ' debuging trd_mld_zint: I.1 done '  
(  225)             CALL FLUSH(numout)
(  226)          ENDIF
(  227) 
(  228) 
(  229)          ! I.2 probability density function of presence in mixed-layer
(  230)          ! --------------------------------
(  231)          ! (i.e. weight of each grid point in vertical integration : wkx(ji,jj,jk)
(  232) 
(  233) 
(  234)          ! initialize wkx with vertical scale factor in mixed-layer
(  235) 
(  236)          wkx(:,:,:) = 0.e0
(  237)          DO jk = 1, jpktrd
(  238)             DO jj = 1,jpj
(  239)                DO ji = 1,jpi
(  240)                   IF( jk - nmld(ji,jj) < 0. )   wkx(ji,jj,jk) = fse3t(ji,jj,jk) * tmask(ji,jj,jk)
(  241)                END DO
(  242)             END DO
(  243)          END DO
(  244)          
(  245)          ! compute mixed-layer depth : rmld
(  246)          
(  247)          rmld(:,:) = 0.
(  248)          DO jk = 1, jpktrd
(  249)             rmld(:,:) = rmld(:,:) + wkx(:,:,jk)
(  250)          END DO
(  251)          
(  252)          ! compute PDF
(  253) 
(  254)          DO jk = 1, jpktrd
(  255)             wkx(:,:,jk) = wkx(:,:,jk) / MAX( 1., rmld(:,:) )
(  256)          END DO
(  257) 
(  258)          IF( idebug /= 0 ) THEN
(  259)              if(lwp) WRITE(numout,*) ' debuging trd_mld_zint: I.2 done '  
(  260)             CALL FLUSH(numout)
(  261)          ENDIF
(  262) 
(  263)          ! Set counter icount to 0 to avoid this part at each time step
(  264)          icount = 0
(  265) 
(  266)       ENDIF
(  267) 
(  268) 
(  269)       !  ====================================================
(  270)       !   II. vertical integration of trends in mixed-layer
(  271)       !  ====================================================
(  272) 
(  273)       ! II.1 vertical integration of 3D and 2D trends
(  274)       ! ---------------------------------------------
(  275) 
(  276)       SELECT CASE (ctype)
(  277) 
(  278)       CASE ('3D')       ! 3D treatment






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 6

(  279) 
(  280)          ! trends terms in the mixed-layer
(  281)          DO jk = 1, jpktrd
(  282)             ! Temperature
(  283)             tmltrd(:,:,ktrd) = tmltrd(:,:,ktrd) + pttrdmld(:,:,jk) * wkx(:,:,jk)   
(  284) 
(  285)             ! Salinity
(  286)             smltrd(:,:,ktrd) = smltrd(:,:,ktrd) + pstrdmld(:,:,jk) * wkx(:,:,jk)   
(  287)          ENDDO
(  288) 
(  289)       CASE ('2D')       ! 2D treatment
(  290) 
(  291)          SELECT CASE (ktrd) 
(  292) 
(  293)          CASE (jpmldldf)
(  294) 
(  295) # if defined key_trabbl_dif
(  296)                ! trends terms from Beckman over-flow parameterization
(  297)                DO jj = 1,jpj
(  298)                   DO ji = 1,jpi
(  299)                      ikb = MAX( mbathy(ji,jj)-1, 1 )
(  300)                      ! beckmann component -> horiz. part of lateral diffusion
(  301)                      tmltrd(ji,jj,ktrd) = tmltrd(ji,jj,ktrd) + pttrdmld(ji,jj,1) * wkx(ji,jj,ikb)
(  302)                      smltrd(ji,jj,ktrd) = smltrd(ji,jj,ktrd) + pstrdmld(ji,jj,1) * wkx(ji,jj,ikb)
(  303)                   END DO
(  304)                END DO
(  305) # endif
(  306) 
(  307)          CASE DEFAULT
(  308) 
(  309)             ! trends terms at upper boundary of mixed-layer
(  310) 
(  311)             ! forcing term (non penetrative)
(  312)             ! Temperature
(  313)             tmltrd(:,:,ktrd) = tmltrd(:,:,ktrd) + pttrdmld(:,:,1) * wkx(:,:,1)   
(  314) 
(  315)             ! forcing term
(  316)             ! Salinity
(  317)             smltrd(:,:,ktrd) = smltrd(:,:,ktrd) + pstrdmld(:,:,1) * wkx(:,:,1)   
(  318) 
(  319)          END SELECT
(  320) 
(  321)       END SELECT
(  322) 
(  323)       IF( idebug /= 0 ) THEN
(  324)          IF(lwp) WRITE(numout,*) ' debuging trd_mld_zint: II.1 done'  
(  325)          CALL FLUSH(numout)
(  326)       ENDIF
(  327) 
(  328)    END SUBROUTINE trd_mld_zint
(  329) 
(  330) 
(  331) 
(  332)    SUBROUTINE trd_mld( kt )
(  333)       !!----------------------------------------------------------------------
(  334)       !!                  ***  ROUTINE trd_mld  ***
(  335)       !! 
(  336)       !! ** Purpose :  computation of cumulated trends over analysis period






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 7

(  337)       !!               and make outputs (NetCDF or DIMG format)
(  338)       !!
(  339)       !! ** Method/usage :
(  340)       !!
(  341)       !! History :
(  342)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  343)       !!----------------------------------------------------------------------
(  344)       !! * Arguments
(  345)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(  346) 
(  347)       !! * Local declarations
(  348)       INTEGER :: ji, jj, jk, jl, ik, it
(  349) 
(  350)       REAL(wp) :: zmean, zavt
(  351) 
(  352)       REAL(wp) ,DIMENSION(jpi,jpj) ::   &
(  353)          ztmltot, ztmlres,              &
(  354)          zsmltot, zsmlres,              & 
(  355)          z2d
(  356) 
(  357) !#if defined key_dimgout
(  358) !      INTEGER ::  iyear,imon,iday
(  359) !      CHARACTER(LEN=80) :: cltext, clmode
(  360) !#endif
(  361)       !!----------------------------------------------------------------------
(  362) 
(  363)       ! I. trends terms at lower boundary of mixed-layer
(  364)       ! ------------------------------------------------
(  365) 
(  366)       DO jj = 1,jpj
(  367)          DO ji = 1,jpi
(  368)             
(  369)             ik = nmld(ji,jj)
(  370)             
(  371)             ! Temperature
(  372)             ! entrainment due to vertical diffusion
(  373)             !       - due to vertical mixing scheme (TKE)
(  374)             zavt = avt(ji,jj,ik)
(  375)             tmltrd(ji,jj,jpmldevd) = - 1. * zavt / fse3w(ji,jj,ik) * tmask(ji,jj,ik)   &
(  376)                &                   * ( tn(ji,jj,ik-1) - tn(ji,jj,ik) )   &
(  377)                &                   / MAX( 1., rmld(ji,jj) ) * tmask(ji,jj,1)
(  378)             ! Salinity
(  379)             ! entrainment due to vertical diffusion
(  380)             !       - due to vertical mixing scheme (TKE)
(  381)             zavt = fsavs(ji,jj,ik)
(  382)             smltrd(ji,jj,jpmldevd) = -1. * zavt / fse3w(ji,jj,ik) * tmask(ji,jj,ik)   &
(  383)                &                  * ( sn(ji,jj,ik-1) - sn(ji,jj,ik) )   &
(  384)                &                  / MAX( 1., rmld(ji,jj) ) * tmask(ji,jj,1)
(  385)          END DO
(  386)       END DO
(  387) 
(  388)       IF( ln_traldf_iso ) THEN
(  389)          ! We substract to the TOTAL vertical diffusion tmltrd(:,:,jpmldzdf) 
(  390)          ! computed in subroutines trazdf_iso.F90 or trazdf_imp.F90
(  391)          ! the vertical part du to the Kz in order to keep only the vertical
(  392)          ! isopycnal diffusion (i.e the isopycnal diffusion componant on the vertical):
(  393)          tmltrd(:,:,jpmldzdf) = tmltrd(:,:,jpmldzdf) - tmltrd(:,:,jpmldevd)   ! - due to isopycnal mixing scheme (implicit p
(  394)          smltrd(:,:,jpmldzdf) = smltrd(:,:,jpmldzdf) - smltrd(:,:,jpmldevd)   ! - due to isopycnal mixing scheme (implicit p






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 8

(  395)       ENDIF
(  396) 
(  397)       ! Boundary conditions
(  398)       CALL lbc_lnk( tmltrd, 'T', 1. )
(  399)       CALL lbc_lnk( smltrd, 'T', 1. )
(  400) 
(  401)       IF( idebug /= 0 ) THEN
(  402)           if(lwp) WRITE(numout,*) ' debuging trd_mld: I. done'  
(  403)          CALL FLUSH(numout)
(  404)       ENDIF
(  405) 
(  406)       !  =================================
(  407)       !   II. Cumulated trends
(  408)       !  =================================
(  409) 
(  410)       ! II.1 set before values of vertically average T and S 
(  411)       ! ---------------------------------------------------
(  412) 
(  413)       IF( kt > nit000 ) THEN
(  414)          tmlb(:,:) = tml(:,:)
(  415)          smlb(:,:) = sml(:,:)
(  416)       ENDIF
(  417) 
(  418)       ! II.2 vertically integrated T and S
(  419)       ! ---------------------------------
(  420) 
(  421)       tml(:,:) = 0.
(  422)       sml(:,:) = 0.
(  423) 
(  424)       DO jk = 1, jpktrd - 1
(  425)          tml(:,:) = tml(:,:) + wkx(:,:,jk) * tn(:,:,jk)
(  426)          sml(:,:) = sml(:,:) + wkx(:,:,jk) * sn(:,:,jk) 
(  427)       END DO
(  428) 
(  429)       IF(idebug /= 0) THEN
(  430)          if(lwp) WRITE(numout,*) ' debuging trd_mld: II.2 done'  
(  431)          CALL FLUSH(numout)
(  432)       ENDIF
(  433) 
(  434)       ! II.3 set `before' mixed layer values for kt = nit000+1
(  435)       ! --------------------------------------------------------
(  436) 
(  437)       IF( kt == nit000+1 ) THEN
(  438)          tmlbb(:,:) = tmlb(:,:)
(  439)          tmlbn(:,:) = tml (:,:)
(  440)          smlbb(:,:) = smlb(:,:)
(  441)          smlbn(:,:) = sml (:,:)
(  442)       ENDIF
(  443) 
(  444)       IF( idebug /= 0 ) THEN
(  445)          if(lwp) WRITE(numout,*) ' debuging trd_mld: II.3 done'  
(  446)          CALL FLUSH(numout)
(  447)       ENDIF
(  448) 
(  449)       ! II.4 cumulated trends over analysis period (kt=2 to nwrite)
(  450)       ! -----------------------------------------------------------
(  451) 
(  452)       ! trends cumulated over nwrite-2 time steps






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 9

(  453) 
(  454)       IF( kt >= nit000+2 ) THEN
(  455)          nmoymltrd = nmoymltrd + 1
(  456)          DO jl = 1, jpltrd
(  457)             tmltrdm(:,:) = tmltrdm(:,:) + tmltrd(:,:,jl)
(  458)             smltrdm(:,:) = smltrdm(:,:) + smltrd(:,:,jl)
(  459)          END DO
(  460)       ENDIF
(  461) 
(  462)       IF( idebug /= 0 ) THEN
(  463)          if(lwp) WRITE(numout,*) ' debuging trd_mld: II.4 done'  
(  464)          CALL FLUSH(numout)
(  465)       ENDIF
(  466) 
(  467)       !  =============================================
(  468)       !   III. Output in netCDF + residual computation
(  469)       !  =============================================
(  470) 
(  471)       ztmltot(:,:) = 0.
(  472)       zsmltot(:,:) = 0.
(  473)       ztmlres(:,:) = 0.
(  474)       zsmlres(:,:) = 0.
(  475) 
(  476)       IF( MOD( kt - nit000+1, nwrite ) == 0 ) THEN
(  477) 
(  478)          ! III.1 compute total trend 
(  479)          ! ------------------------
(  480) 
(  481)          zmean = float(nmoymltrd)
(  482)          
(  483)          ztmltot(:,:) = ( tml(:,:) - tmlbn(:,:) + tmlb(:,:) - tmlbb(:,:) ) /  (zmean * 2. * rdt)
(  484)          zsmltot(:,:) = ( sml(:,:) - smlbn(:,:) + smlb(:,:) - smlbb(:,:) ) /  (zmean * 2. * rdt)
(  485) 
(  486)          IF(idebug /= 0) THEN
(  487)             if(lwp) WRITE(numout,*) ' zmean = ',zmean  
(  488)             if(lwp) WRITE(numout,*) ' debuging trd_mld: III.1 done'  
(  489)             CALL FLUSH(numout)
(  490)          ENDIF
(  491)           
(  492) 
(  493)          ! III.2 compute residual 
(  494)          ! ---------------------
(  495) 
(  496)          ztmlres(:,:) = ztmltot(:,:) - tmltrdm(:,:) / zmean
(  497)          zsmlres(:,:) = zsmltot(:,:) - smltrdm(:,:) / zmean
(  498) 
(  499) 
(  500)          ! Boundary conditions
(  501) 
(  502)          CALL lbc_lnk( ztmltot, 'T', 1. )
(  503)          CALL lbc_lnk( ztmlres, 'T', 1. )
(  504)          CALL lbc_lnk( zsmltot, 'T', 1. )
(  505)          CALL lbc_lnk( zsmlres, 'T', 1. )
(  506) 
(  507)          IF( idebug /= 0 ) THEN
(  508)             if(lwp) WRITE(numout,*) ' debuging trd_mld: III.2 done'  
(  509)             CALL FLUSH(numout)
(  510)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 10

(  511) 
(  512) 
(  513)          ! III.3 time evolution array swap
(  514)          ! ------------------------------
(  515) 
(  516)          tmlbb(:,:) = tmlb(:,:)
(  517)          tmlbn(:,:) = tml (:,:)
(  518)          smlbb(:,:) = smlb(:,:)
(  519)          smlbn(:,:) = sml (:,:)
(  520) 
(  521)          IF( idebug /= 0 ) THEN
(  522)             if(lwp) WRITE(numout,*) ' debuging trd_mld: III.3 done'  
(  523)             CALL FLUSH(numout)
(  524)          ENDIF
(  525) 
(  526) 
(  527)          ! III.4 zero cumulative array
(  528)          ! ---------------------------
(  529) 
(  530)           nmoymltrd = 0
(  531) 
(  532)           tmltrdm(:,:) = 0.
(  533)           smltrdm(:,:) = 0.
(  534) 
(  535)           IF(idebug /= 0) THEN
(  536)              if(lwp) WRITE(numout,*) ' debuging trd_mld: III.4 done'  
(  537)              CALL FLUSH(numout)
(  538)           ENDIF
(  539)           
(  540)       ENDIF
(  541) 
(  542)       ! III.5 write trends to output
(  543)       ! ---------------------------
(  544) 
(  545)       IF( kt >=  nit000+1 ) THEN
(  546) 
(  547)          ! define time axis
(  548)          it= kt-nit000+1
(  549)          IF( lwp .AND. MOD( kt, nwrite ) == 0 ) THEN
(  550)             WRITE(numout,*) '     trd_mld : write NetCDF fields'
(  551)          ENDIF
(  552)          
(  553)          CALL histwrite( nidtrd,"somlttml",it,rmld          ,ndimtrd1,ndextrd1) ! Mixed-layer depth
(  554)          
(  555)          ! Temperature trends
(  556)          ! ------------------
(  557)          CALL histwrite( nidtrd,"somltemp",it,tml           ,ndimtrd1,ndextrd1) ! Mixed-layer temperature
(  558)          CALL histwrite( nidtrd,"somlttto",it,ztmltot       ,ndimtrd1,ndextrd1) ! total 
(  559)          CALL histwrite( nidtrd,"somlttax",it,tmltrd(:,:, 1),ndimtrd1,ndextrd1) ! i- adv.
(  560)          CALL histwrite( nidtrd,"somlttay",it,tmltrd(:,:, 2),ndimtrd1,ndextrd1) ! j- adv.
(  561)          CALL histwrite( nidtrd,"somlttaz",it,tmltrd(:,:, 3),ndimtrd1,ndextrd1) ! vertical adv.
(  562)          CALL histwrite( nidtrd,"somlttdh",it,tmltrd(:,:, 4),ndimtrd1,ndextrd1) ! hor. lateral diff.
(  563)          CALL histwrite( nidtrd,"somlttfo",it,tmltrd(:,:, 5),ndimtrd1,ndextrd1) ! forcing
(  564) 
(  565)          CALL histwrite( nidtrd,"somlbtdz",it,tmltrd(:,:, 6),ndimtrd1,ndextrd1) ! vert. diffusion 
(  566)          CALL histwrite( nidtrd,"somlbtdt",it,ztmlres       ,ndimtrd1,ndextrd1) ! dh/dt entrainment (residual)
(  567)          IF( ln_traldf_iso ) THEN
(  568)             CALL histwrite( nidtrd,"somlbtdv",it,tmltrd(:,:, 7),ndimtrd1,ndextrd1) ! vert. lateral diff.






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 11

(  569)          ENDIF
(  570) #if defined key_traldf_eiv
(  571)          CALL histwrite( nidtrd,"somlgtax",it,tmltrd(:,:, 8),ndimtrd1,ndextrd1) ! i- adv. (eiv)
(  572)          CALL histwrite( nidtrd,"somlgtay",it,tmltrd(:,:, 9),ndimtrd1,ndextrd1) ! j- adv. (eiv)
(  573)          CALL histwrite( nidtrd,"somlgtaz",it,tmltrd(:,:,10),ndimtrd1,ndextrd1) ! vert. adv. (eiv)
(  574)          z2d(:,:) = tmltrd(:,:,8) + tmltrd(:,:,9) + tmltrd(:,:,10)
(  575)          CALL histwrite( nidtrd,"somlgtat",it,z2d           ,ndimtrd1,ndextrd1) ! total adv. (eiv)
(  576) #endif   
(  577) 
(  578)          ! Salinity trends
(  579)          ! ---------------
(  580)          CALL histwrite( nidtrd,"somlsalt",it,sml           ,ndimtrd1,ndextrd1) ! Mixed-layer salinity
(  581)          CALL histwrite( nidtrd,"somltsto",it,zsmltot       ,ndimtrd1,ndextrd1) ! total 
(  582)          CALL histwrite( nidtrd,"somltsax",it,smltrd(:,:, 1),ndimtrd1,ndextrd1) ! i- adv.
(  583)          CALL histwrite( nidtrd,"somltsay",it,smltrd(:,:, 2),ndimtrd1,ndextrd1) ! j- adv.
(  584)          CALL histwrite( nidtrd,"somltsaz",it,smltrd(:,:, 3),ndimtrd1,ndextrd1) ! vert. adv.
(  585)          CALL histwrite( nidtrd,"somltsdh",it,smltrd(:,:, 4),ndimtrd1,ndextrd1) ! hor. lateral diff.
(  586)          CALL histwrite( nidtrd,"somltsfo",it,smltrd(:,:, 5),ndimtrd1,ndextrd1) ! forcing
(  587)          CALL histwrite( nidtrd,"somlbsdz",it,smltrd(:,:, 6),ndimtrd1,ndextrd1) ! vert. diff.
(  588)          CALL histwrite( nidtrd,"somlbsdt",it,zsmlres       ,ndimtrd1,ndextrd1) ! dh/dt entrainment (residual)
(  589)          IF( ln_traldf_iso ) THEN
(  590)             CALL histwrite( nidtrd,"somlbsdv",it,smltrd(:,:, 7),ndimtrd1,ndextrd1) ! vert. lateral diff;
(  591)          ENDIF
(  592) #if defined key_traldf_eiv
(  593)          CALL histwrite( nidtrd,"somlgsax",it,smltrd(:,:, 8),ndimtrd1,ndextrd1) ! i-adv. (eiv)
(  594)          CALL histwrite( nidtrd,"somlgsay",it,smltrd(:,:, 9),ndimtrd1,ndextrd1) ! j-adv. (eiv)
(  595)          CALL histwrite( nidtrd,"somlgsaz",it,smltrd(:,:,10),ndimtrd1,ndextrd1) ! vert. adv. (eiv)
(  596)          z2d(:,:) = smltrd(:,:,8) + smltrd(:,:,9) + smltrd(:,:,10)
(  597)          CALL histwrite( nidtrd,"somlgsat",it,z2d           ,ndimtrd1,ndextrd1) ! total adv. (eiv)
(  598) #endif
(  599) 
(  600)          IF( idebug /= 0 ) THEN
(  601)             if(lwp) WRITE(numout,*) ' debuging trd_mld: III.5 done'  
(  602)             CALL FLUSH(numout)
(  603)          ENDIF
(  604) 
(  605)          ! set counter icount to one to allow the calculation
(  606)          ! of the surface control in the next time step in the trd_mld_zint subroutine
(  607)          icount = 1
(  608) 
(  609)       ENDIF
(  610) 
(  611)       ! At the end of the 1st time step, set icount to 1 to be
(  612)       ! able to compute the surface control at the beginning of
(  613)       ! the second time step
(  614)       IF( kt == nit000 )   icount = 1
(  615) 
(  616)       IF( kt == nitend )   CALL histclo( nidtrd )
(  617) 
(  618)    END SUBROUTINE trd_mld
(  619) 
(  620) 
(  621) 
(  622)    SUBROUTINE trd_mld_init
(  623)       !!----------------------------------------------------------------------
(  624)       !!                  ***  ROUTINE trd_mld_init  ***
(  625)       !! 
(  626)       !! ** Purpose :   computation of vertically integrated T and S budgets






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 12

(  627)       !!      from ocean surface down to control surface (NetCDF output)
(  628)       !!
(  629)       !! ** Method/usage :
(  630)       !!
(  631)       !! History :
(  632)       !!        !  95-04  (J. Vialard)  Original code
(  633)       !!        !  97-02  (E. Guilyardi)  Adaptation global + base cmo
(  634)       !!        !  99-09  (E. Guilyardi)  Re-writing + netCDF output
(  635)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  636)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  637)       !!----------------------------------------------------------------------
(  638)       !! * Local declarations
(  639)       INTEGER :: ilseq
(  640) 
(  641)       REAL(wp) ::   zjulian, zsto, zout
(  642) 
(  643)       CHARACTER (LEN=21) ::   &
(  644)          clold ='OLD'        , & ! open specifier (direct access files)
(  645)          clunf ='UNFORMATTED', & ! open specifier (direct access files)
(  646)          clseq ='SEQUENTIAL'     ! open specifier (direct access files)
(  647)       CHARACTER (LEN=40) ::   clhstnam
(  648)       CHARACTER (LEN=40) ::   clop
(  649)       CHARACTER (LEN=12) ::   clmxl
(  650) 
(  651)       NAMELIST/namtrd/ ntrd, nctls
(  652)       !!----------------------------------------------------------------------
(  653) 
(  654)       !  ===================
(  655)       !   I. initialization
(  656)       !  ===================
(  657) 
(  658)       ! Open specifier
(  659)       ilseq  = 1
(  660)       idebug = 0      ! set it to 1 in case of problem to have more print
(  661)       icount = 1      
(  662)       ionce  = 1
(  663) 
(  664)       ! namelist namtrd : trend diagnostic
(  665)       REWIND( numnam )
(  666)       READ  ( numnam, namtrd )
(  667) 
(  668)       IF(lwp) THEN
(  669)          WRITE(numout,*) ' '
(  670)          WRITE(numout,*) 'trd_mld_init: mixed layer heat & freshwater budget trends'
(  671)          WRITE(numout,*) '~~~~~~~~~~~~~'
(  672)          WRITE(numout,*) ' '
(  673)          WRITE(numout,*) '          Namelist namtrd : '
(  674)          WRITE(numout,*) '             control surface for trends      nctls = ',nctls
(  675)          WRITE(numout,*) ' '
(  676)       ENDIF
(  677) 
(  678)       ! cumulated trends array init
(  679)       nmoymltrd = 0
(  680)       tmltrdm(:,:) = 0.e0
(  681)       smltrdm(:,:) = 0.e0
(  682) 
(  683)       !  read control surface from file ctlsurf_idx
(  684) 






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 13

(  685)       IF( nctls == 1 ) THEN
(  686)          clname ='ctlsurf_idx'
(  687)          CALL ctlopn(numbol,clname,clold,clunf,clseq,   &
(  688)               ilseq,numout,lwp,1)
(  689)          REWIND (numbol)
(  690)          READ(numbol) nbol
(  691)       ENDIF
(  692) 
(  693) 
(  694)       IF( idebug /= 0 ) THEN
(  695)          if(lwp) WRITE(numout,*) ' debuging trd_mld_init: 0. done '  
(  696)          CALL FLUSH(numout)
(  697)       ENDIF
(  698) 
(  699)       !  ===================================
(  700)       !   II. netCDF output initialization
(  701)       !  ===================================
(  702) 
(  703)       !     clmxl = legend root for netCDF output
(  704)       IF( nctls == 0 ) THEN
(  705)          ! control surface = mixed-layer with density criterion 
(  706)          ! (array nmln computed in zdfmxl.F90)
(  707)          clmxl = 'Mixed Layer '
(  708)       ELSE IF( nctls == 1 ) THEN
(  709)          ! control surface = read index from file 
(  710)          clmxl = '      Bowl '
(  711)       ELSE IF( nctls >= 2 ) THEN
(  712)          ! control surface = model level
(  713)          WRITE(clmxl,'(A9,I2,1X)') 'Levels 1-', nctls
(  714)       ENDIF
(  715) 
(  716)       !-----------------------------------------
(  717)       ! II.1 Define frequency of output and means
(  718)       ! -----------------------------------------
(  719) 
(  720) #if defined key_diainstant
(  721)       zsto = nwrite*rdt
(  722)       clop ="inst(x)"
(  723) #else
(  724)       zsto = rdt
(  725)       clop ="ave(x)"
(  726) #endif
(  727)       zout = nwrite*rdt
(  728) 
(  729)       IF(lwp) WRITE (numout,*) ' trdmld_ncinit: netCDF initialization'
(  730) 
(  731)       ! II.2 Compute julian date from starting date of the run
(  732)       ! ------------------------
(  733) 
(  734)       CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )
(  735)       IF (lwp) WRITE(numout,*)' '  
(  736)       IF (lwp) WRITE(numout,*)' Date 0 used :',nit000   &
(  737)            ,' YEAR ', nyear,' MONTH ', nmonth,' DAY ', nday   &
(  738)            ,'Julian day : ', zjulian
(  739) 
(  740) 
(  741)       ! II.3 Define the T grid trend file (nidtrd)
(  742)       ! ---------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 14

(  743) 
(  744)       CALL dia_nam( clhstnam, nwrite, 'trends' )                  ! filename
(  745)       IF(lwp) WRITE(numout,*) ' Name of NETCDF file ', clhstnam
(  746)       CALL histbeg( clhstnam, jpi, glamt, jpj, gphit,1, jpi,   &  ! Horizontal grid : glamt and gphit
(  747)          &          1, jpj, 0, zjulian, rdt, nh_t, nidtrd, domain_id=nidom )
(  748) 
(  749)       ! Declare output fields as netCDF variables
(  750) 
(  751)       ! Mixed layer Depth
(  752)       CALL histdef( nidtrd, "somlttml", clmxl//"Depth"              , "m"   ,   &  ! hmlp
(  753)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  754) 
(  755)       ! Temperature
(  756)       CALL histdef( nidtrd, "somltemp", clmxl//"Temperature"        , "C"   ,   &  ! ???
(  757)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  758)       ! Temperature trends
(  759)       CALL histdef( nidtrd, "somlttto", clmxl//"T Total"             , "C/s",   &  ! total
(  760)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zout, zout )
(  761)       CALL histdef( nidtrd, "somlttax", clmxl//"T Zonal Advection", "C/s",       & ! i-adv.
(  762)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  763)       CALL histdef( nidtrd, "somlttay", clmxl//"T Meridional Advection", "C/s",   & ! j-adv.
(  764)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  765)       CALL histdef( nidtrd, "somlttaz", clmxl//"T Vertical Advection", "C/s",   & ! vert. adv.
(  766)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  767)       CALL histdef( nidtrd, "somlttdh", clmxl//"T Horizontal Diffusion ", "C/s",   & ! hor. lateral diff.
(  768)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  769)       CALL histdef( nidtrd, "somlttfo", clmxl//"T Forcing", "C/s",   & ! forcing
(  770)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  771)       CALL histdef( nidtrd, "somlbtdz", clmxl//"T Vertical Diffusion", "C/s",   & ! vert. diff.
(  772)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  773)       CALL histdef( nidtrd, "somlbtdt", clmxl//"T dh/dt Entrainment (Residual)", "C/s",   & ! T * dh/dt 
(  774)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zout, zout )
(  775)       IF( ln_traldf_iso ) THEN
(  776)       CALL histdef( nidtrd, "somlbtdv", clmxl//"T Vert. lateral Diffusion","C/s",   & ! vertical diffusion entrainment (ISO)
(  777)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  778)       ENDIF
(  779) #if defined key_traldf_eiv
(  780)       CALL histdef( nidtrd, "somlgtax", clmxl//"T Zonal EIV Advection", "C/s",   & ! i-adv. (eiv)
(  781)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  782)       CALL histdef( nidtrd, "somlgtay", clmxl//"T Meridional EIV Advection", "C/s",   & ! j-adv. (eiv)
(  783)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  784)       CALL histdef( nidtrd, "somlgtaz", clmxl//"T Vertical EIV Advection", "C/s",   & ! vert. adv. (eiv)
(  785)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  786)       CALL histdef( nidtrd, "somlgtat", clmxl//"T Total EIV Advection", "C/s",   & ! total advection (eiv)
(  787)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  788) #endif
(  789)       ! Salinity
(  790)       CALL histdef( nidtrd, "somlsalt", clmxl//"Salinity", "PSU",   & ! Mixed-layer salinity
(  791)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  792)       ! Salinity trends
(  793)       CALL histdef( nidtrd, "somltsto", clmxl//"S Total", "PSU/s",   & ! total 
(  794)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  795)       CALL histdef( nidtrd, "somltsax", clmxl//"S Zonal Advection", "PSU/s",   & ! i-advection
(  796)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  797)       CALL histdef( nidtrd, "somltsay", clmxl//"S Meridional Advection", "PSU/s",   & ! j-advection
(  798)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  799)       CALL histdef( nidtrd, "somltsaz", clmxl//"S Vertical Advection", "PSU/s",   & ! vertical advection
(  800)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 15

(  801)       CALL histdef( nidtrd, "somltsdh", clmxl//"S Horizontal Diffusion ", "PSU/s",   & ! hor. lat. diff.
(  802)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  803)       CALL histdef( nidtrd, "somltsfo", clmxl//"S Forcing", "PSU/s",   & ! forcing
(  804)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  805) 
(  806)       CALL histdef( nidtrd, "somlbsdz", clmxl//"S Vertical Diffusion", "PSU/s",   & ! vert. diff.
(  807)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  808)       CALL histdef( nidtrd, "somlbsdt", clmxl//"S dh/dt Entrainment (Residual)", "PSU/s",   & ! S * dh/dt 
(  809)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  810)       IF( ln_traldf_iso ) THEN
(  811)       ! vertical diffusion entrainment (ISO)
(  812)       CALL histdef( nidtrd, "somlbsdv", clmxl//"S Vertical lateral Diffusion", "PSU/s",   & ! vert. lat. diff.
(  813)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  814)       ENDIF
(  815) #if defined key_traldf_eiv
(  816)       CALL histdef( nidtrd, "somlgsax", clmxl//"S Zonal EIV Advection", "PSU/s",   & ! i-advection (eiv)
(  817)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  818)       CALL histdef( nidtrd, "somlgsay", clmxl//"S Meridional EIV Advection", "PSU/s",   & ! j-advection (eiv)
(  819)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  820)       CALL histdef( nidtrd, "somlgsaz", clmxl//"S Vertical EIV Advection", "PSU/s",   & ! vert. adv. (eiv)
(  821)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  822)       CALL histdef( nidtrd, "somlgsat", clmxl//"S Total EIV Advection", "PSU/s",   & ! total adv. (eiv)
(  823)          &          jpi, jpj, nh_t, 1  , 1, 1  , -99 , 32, clop, zsto, zout )
(  824) #endif
(  825)       CALL histend( nidtrd )
(  826) 
(  827)       IF( idebug /= 0 ) THEN
(  828)          if(lwp) WRITE(numout,*) ' debuging trd_mld_init: II. done'  
(  829)          CALL FLUSH(numout)
(  830)       ENDIF
(  831) 
(  832) 
(  833)       END SUBROUTINE trd_mld_init
(  834) 
(  835) #else
(  836)    !!----------------------------------------------------------------------
(  837)    !!   Default option :                                       Empty module
(  838)    !!----------------------------------------------------------------------
(  839)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .FALSE.   !: momentum trend flag
(  840) CONTAINS
(  841)    SUBROUTINE trd_mld( kt )             ! Empty routine
(  842)       INTEGER, INTENT( in) ::   kt
(  843) !      WRITE(*,*) 'trd_mld: You should not have seen this print! error?', kt
(  844)    END SUBROUTINE trd_mld
(  845)    SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
(  846)       REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
(  847)          pttrdmld, pstrdmld                   ! Temperature and Salinity trends
(  848)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  849)       CHARACTER(len=2), INTENT( in ) ::   &  
(  850)          ctype                                ! surface/bottom (2D arrays) or
(  851)          !                                    ! interior (3D arrays) physics
(  852) !      WRITE(*,*) 'trd_mld_zint: You should not have seen this print! error?', pttrdmld(1,1,1)
(  853) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', pstrdmld(1,1,1)
(  854) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ctype
(  855) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ktrd
(  856)    END SUBROUTINE trd_mld_zint
(  857)    SUBROUTINE trd_mld_init              ! Empty routine
(  858) !      WRITE(*,*) 'trd_mld_init: You should not have seen this print! error?'






PGF90 (Version     10.5)          09/23/2010  09:49:10      page 16

(  859)    END SUBROUTINE trd_mld_init
(  860) #endif
(  861) 
(  862)    !!======================================================================
(  863) END MODULE trdmld



























































PGF90 (Version     10.5)          09/23/2010  09:49:10      page 17

(    1) # 1 "trdmld.F90"
(    1) MODULE trdmld
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  trdmld  ***
(    4)    !! Ocean diagnostics:  mixed layer T-S trends 
(    5)    !!=====================================================================
(    6) # 836
(  836)    !!----------------------------------------------------------------------
(  837)    !!   Default option :                                       Empty module
(  838)    !!----------------------------------------------------------------------
(  839)    LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .FALSE.   !: momentum trend flag
(  840) CONTAINS
(  841)    SUBROUTINE trd_mld( kt )             ! Empty routine
(  842)       INTEGER, INTENT( in) ::   kt
(  843) !      WRITE(*,*) 'trd_mld: You should not have seen this print! error?', kt
(  844)    END SUBROUTINE trd_mld
(  845)    SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
(  846)       REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
(  847)          pttrdmld, pstrdmld                   ! Temperature and Salinity trends
(  848)       INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
(  849)       CHARACTER(len=2), INTENT( in ) ::   &  
(  850)          ctype                                ! surface/bottom (2D arrays) or
(  851)          !                                    ! interior (3D arrays) physics
(  852) !      WRITE(*,*) 'trd_mld_zint: You should not have seen this print! error?', pttrdmld(1,1,1)
(  853) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', pstrdmld(1,1,1)
(  854) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ctype
(  855) !      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ktrd
(  856)    END SUBROUTINE trd_mld_zint
(  857)    SUBROUTINE trd_mld_init              ! Empty routine
(  858) !      WRITE(*,*) 'trd_mld_init: You should not have seen this print! error?'
(  859)    END SUBROUTINE trd_mld_init
(  860) 
(  861) # 862
(  862)    !!======================================================================
(  863) END MODULE trdmld
