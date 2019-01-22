


PGF90 (Version     10.5)          09/23/2010  09:49:14      page 1

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

Filename: ocesbc.F90

(    1) !!DB: 2009.08.31 -- Eliminated GYRE config
(    2) MODULE ocesbc
(    3)    !!======================================================================
(    4)    !!                     ***  MODULE  ocesbc  ***
(    5)    !!                     Ocean surface boundary conditions
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   oce_sbc     : ???
(   10)    !!   oce_sbc_dmp : ???
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce            ! dynamics and tracers variables
(   14)    USE dom_oce        ! ocean space domain variables
(   15)    USE cpl_oce        ! coupled ocean-atmosphere variables
(   16)    USE ice_oce        ! sea-ice variable
(   17)    USE blk_oce        ! bulk variables
(   18)    USE flx_oce        ! sea-ice/ocean forcings variables
(   19)    USE phycst         ! Define parameters for the routines
(   20)    USE taumod         ! surface stress forcing
(   21)    USE flxmod         ! thermohaline fluxes
(   22)    USE flxrnf         ! runoffs forcing
(   23)    USE tradmp         ! damping salinity trend
(   24)    USE dtatem         ! ocean temperature data
(   25)    USE dtasal         ! ocean salinity data
(   26)    USE ocfzpt         ! surface ocean freezing point
(   27)    USE lbclnk         ! ocean lateral boundary condition
(   28)    USE lib_mpp        ! distribued memory computing library
(   29)    USE in_out_manager ! I/O manager
(   30)    USE prtctl         ! Print control
(   31) 
(   32)    IMPLICIT NONE
(   33)    PRIVATE
(   34) 
(   35)    !! * Accessibility
(   36)    PUBLIC oce_sbc    ! routine called by step
(   37) 
(   38)    !! * Shared module variables
(   39)    REAL(wp), PUBLIC ::   &  !:
(   40)       aplus, aminus,     &  !:
(   41)       empold = 0.e0         !: current year freshwater budget correction
(   42)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   43)       qt  ,         &  !: total surface heat flux (w/m2)
(   44)       qsr ,         &  !: solar radiation (w/m2)
(   45)       emp ,         &  !: evaporation minus precipitation (kg/m2/s = mm/s)
(   46)       emps,         &  !: evaporation - precipitation (free surface)






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 2

(   47)       qrp ,         &  !: heat flux damping (w/m2)
(   48)       erp              !: evaporation damping (kg/m2/s = mm/s)
(   49) #if ! defined key_dynspg_rl
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   51)       dmp              !: internal dampind term
(   52) #endif
(   53) 
(   54) #  include "domzgr_substitute.h90"
(   55) #  include "vectopt_loop_substitute.h90"
(   56) 
(   57)    !!----------------------------------------------------------------------
(   58)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   59)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/ocesbc.F90,v 1.16 2006/04/19 14:43:16 opalod Exp $ 
(   60)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   61)    !!----------------------------------------------------------------------
(   62) CONTAINS
(   63) 
(   64) #if defined key_ice_lim
(   65)    !!----------------------------------------------------------------------
(   66)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   67)    !!----------------------------------------------------------------------
(   68) # if defined key_coupled
(   69)       !!----------------------------------------------------------------------
(   70)       !!   'key_coupled' :                            Coupled Ocean/Atmosphere
(   71)       !!----------------------------------------------------------------------
(   72) 
(   73)    SUBROUTINE oce_sbc( kt )
(   74)       !!---------------------------------------------------------------------
(   75)       !!                 ***  ROUTINE oce_sbc  ***
(   76)       !!                    
(   77)       !! ** Purpose :   Ocean surface boundaries conditions with 
(   78)       !!        Louvain la Neuve Sea Ice Model in coupled mode
(   79)       !!
(   80)       !! History :
(   81)       !!   1.0  !  00-10  (O. Marti)  Original code
(   82)       !!   2.0  !  02-12  (G. Madec)  F90: Free form and module
(   83)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   84)       !!----------------------------------------------------------------------
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER ::   ji, jj                   ! dummy loop indices
(   90)       REAL(wp) ::   ztx, ztaux, zty, ztauy
(   91)       REAL(wp) ::   ztdta, ztgel, zqrp
(   92)       !!----------------------------------------------------------------------
(   93)  
(   94)       ! 1. initialization to zero at kt = nit000
(   95)       ! ---------------------------------------
(   96)       
(   97)       IF( kt == nit000 ) THEN     
(   98)          qsr   (:,:) = 0.e0
(   99)          freeze(:,:) = 0.e0
(  100)          qt    (:,:) = 0.e0
(  101)          qrp   (:,:) = 0.e0
(  102)          emp   (:,:) = 0.e0
(  103)          emps  (:,:) = 0.e0
(  104)          erp   (:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 3

(  105) #if ! defined key_dynspg_rl 
(  106)          dmp   (:,:) = 0.e0
(  107) #endif
(  108)       ENDIF
(  109) 
(  110)       IF( MOD( kt-1, nfice ) == 0 ) THEN 
(  111) 
(  112)          CALL oce_sbc_dmp   ! Computation of internal and evaporation damping terms       
(  113) 
(  114)          ! Surface heat flux (W/m2)
(  115)          ! -----------------------
(  116) 
(  117)          ! restoring heat flux
(  118)          DO jj = 1, jpj
(  119)             DO ji = 1, jpi
(  120)                ztgel = fzptn(ji,jj)
(  121) #if defined key_dtasst
(  122)                ztdta = MAX( sst(ji,jj),    ztgel )
(  123) #else
(  124)                ztdta = MAX( t_dta(ji,jj,1), ztgel )
(  125) #endif
(  126)                zqrp = dqdt0 * ( tb(ji,jj,1) - ztdta )
(  127) 
(  128)                qrp(ji,jj) = (1.0-freeze(ji,jj) ) * zqrp
(  129)             END DO
(  130)          END DO
(  131) 
(  132)          ! non solar heat flux + solar flux + restoring
(  133)          qt (:,:) = fnsolar(:,:) + fsolar(:,:) + qrp(:,:)
(  134) 
(  135)          ! solar flux
(  136)          qsr(:,:) = fsolar(:,:)
(  137) 
(  138) #if ! defined key_dynspg_rl  
(  139)          ! total concentration/dilution effect (use on SSS)
(  140)          emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)
(  141) 
(  142)          ! total volume flux (use on sea-surface height)
(  143)          emp (:,:) = fmass(:,:)  -  dmp(:,:) + runoff(:,:) + erp(:,:)
(  144) #else
(  145)          ! Rigid-lid (emp=emps=E-P-R+Erp) 
(  146)          ! freshwater flux
(  147)          emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)
(  148)          emp (:,:) = emps(:,:)
(  149) #endif
(  150) 
(  151)          DO jj = 1, jpjm1
(  152)             DO ji = 1, fs_jpim1   ! vertor opt.
(  153)                ztx   = 0.5 * ( freeze(ji+1,jj) + freeze(ji+1,jj+1) )
(  154)                ztaux = 0.5 * ( ftaux (ji+1,jj) + ftaux (ji+1,jj+1) )
(  155)                taux(ji,jj) = (1.0-ztx) * taux(ji,jj) + ztx * ztaux
(  156) 
(  157)                zty   = 0.5 * ( freeze(ji,jj+1) + freeze(ji+1,jj+1) )
(  158)                ztauy = 0.5 * ( ftauy (ji,jj+1) + ftauy (ji+1,jj+1) )
(  159)                tauy(ji,jj) = (1.0-zty) * tauy(ji,jj) + zty * ztauy
(  160)             END DO
(  161)          END DO
(  162)          CALL lbc_lnk( taux, 'U', -1. )






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 4

(  163)          CALL lbc_lnk( tauy, 'V', -1. )    
(  164) 
(  165)          ! Re-initialization of fluxes
(  166)          sst_io(:,:) = 0.e0
(  167)          sss_io(:,:) = 0.e0
(  168)          u_io  (:,:) = 0.e0
(  169)          v_io  (:,:) = 0.e0
(  170)          gtaux (:,:) = 0.e0
(  171)          gtauy (:,:) = 0.e0
(  172) 
(  173)       ENDIF
(  174) 
(  175)    END SUBROUTINE oce_sbc
(  176) 
(  177) # elif defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  178)       !!----------------------------------------------------------------------
(  179)       !!   'key_ice_lim'                              with  LIM sea-ice model
(  180)       !!----------------------------------------------------------------------
(  181) 
(  182)    SUBROUTINE oce_sbc( kt )
(  183)       !!---------------------------------------------------------------------
(  184)       !!                   ***  ROUTINE oce_sbc  ***
(  185)       !!                    
(  186)       !! ** Purpose : - Ocean surface boundary conditions with LIM sea-ice
(  187)       !!        model in forced mode using bulk formulea
(  188)       !!
(  189)       !! History :
(  190)       !!   1.0  !  99-11  (M. Imbard)  Original code
(  191)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(  192)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(  193)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  194)       !!----------------------------------------------------------------------
(  195)       !! * arguments
(  196)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  197) 
(  198)       !! * Local declarations
(  199)       INTEGER  ::   ji, jj                   ! dummy loop indices
(  200)       REAL(wp) ::   ztx, ztaux, zty, ztauy
(  201)       !!----------------------------------------------------------------------
(  202) 
(  203)       ! 1. initialization to zero at kt = nit000
(  204)       ! ---------------------------------------
(  205)       
(  206)       IF( kt == nit000 ) THEN     
(  207)          qsr    (:,:) = 0.e0
(  208)          qt     (:,:) = 0.e0
(  209)          qrp    (:,:) = 0.e0
(  210)          emp    (:,:) = 0.e0
(  211)          emps   (:,:) = 0.e0
(  212)          erp    (:,:) = 0.e0
(  213) #if ! defined key_dynspg_rl 
(  214)          dmp    (:,:) = 0.e0
(  215) #endif
(  216)       ENDIF
(  217) 
(  218)       IF( MOD( kt-1, nfice ) == 0 ) THEN
(  219) 
(  220)          CALL oce_sbc_dmp       ! Computation of internal and evaporation damping terms       






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 5

(  221) 
(  222)          ! Surface Ocean fluxes
(  223)          ! ====================
(  224) 
(  225)          ! Surface heat flux (W/m2)
(  226)          ! -----------------
(  227) 
(  228)          qt  (:,:) = fnsolar(:,:) + fsolar(:,:)     ! non solar heat flux + solar flux
(  229)          qsr (:,:) = fsolar(:,:)                     ! solar flux
(  230) 
(  231) #if ! defined key_dynspg_rl     
(  232)          ! total concentration/dilution effect (use on SSS)
(  233)          emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:) + empold
(  234) 
(  235)          ! total volume flux (use on sea-surface height)
(  236)          emp (:,:) = fmass(:,:) -   dmp(:,:) + runoff(:,:) + erp(:,:) + empold      
(  237) #else
(  238)          ! Rigid-lid (emp=emps=E-P-R+Erp)
(  239)          emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)     ! freshwater flux
(  240)          emp (:,:) = emps(:,:)
(  241) 
(  242) #endif
(  243) 
(  244)          ! Surface stress
(  245)          ! --------------
(  246) 
(  247)          ! update the stress beloww sea-ice area
(  248)          DO jj = 1, jpjm1
(  249)             DO ji = 1, fs_jpim1   ! vertor opt.
(  250)                ztx         = MAX( freezn(ji,jj), freezn(ji,jj+1) )   ! ice/ocean indicator at U- and V-points
(  251)                zty         = MAX( freezn(ji,jj), freezn(ji+1,jj) )
(  252)                ztaux       = 0.5 *( ftaux(ji+1,jj) + ftaux(ji+1,jj+1) ) ! ice-ocean stress at U- and V-points
(  253)                ztauy       = 0.5 *( ftauy(ji,jj+1) + ftauy(ji+1,jj+1) )
(  254)                taux(ji,jj) = (1.-ztx) * taux(ji,jj) + ztx * ztaux    ! stress at the ocean surface
(  255)                tauy(ji,jj) = (1.-zty) * tauy(ji,jj) + zty * ztauy
(  256)             END DO
(  257)          END DO
(  258) 
(  259)          ! boundary condition on the stress (taux,tauy)
(  260)          CALL lbc_lnk( taux, 'U', -1. )
(  261)          CALL lbc_lnk( tauy, 'V', -1. )
(  262) 
(  263)          ! Re-initialization of fluxes
(  264)          sst_io(:,:) = 0.e0
(  265)          sss_io(:,:) = 0.e0
(  266)          u_io  (:,:) = 0.e0
(  267)          v_io  (:,:) = 0.e0
(  268) 
(  269)       ENDIF
(  270) 
(  271)    END SUBROUTINE oce_sbc
(  272) 
(  273) # else
(  274)       !!----------------------------------------------------------------------
(  275)       !!   Error option               LIM sea-ice model requires bulk formulea
(  276)       !!----------------------------------------------------------------------
(  277)       This line forced a compilation error
(  278) # endif






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 6

(  279) 
(  280) #else
(  281)    !!----------------------------------------------------------------------
(  282)    !!   Default option                                 NO LIM sea-ice model
(  283)    !!----------------------------------------------------------------------
(  284) # if defined key_coupled
(  285)       !!----------------------------------------------------------------------
(  286)       !!   'key_coupled' :                            Coupled Ocean/Atmosphere
(  287)       !!----------------------------------------------------------------------
(  288) 
(  289)    SUBROUTINE oce_sbc( kt )
(  290)       !!---------------------------------------------------------------------
(  291)       !!                 ***  ROUTINE oce_sbc  ***
(  292)       !!                    
(  293)       !! ** Purpose :   Ocean surface boundaries conditions in 
(  294)       !!                coupled ocean/atmosphere case without sea-ice
(  295)       !!
(  296)       !! History :
(  297)       !!   1.0  !  00-10  (O. Marti)  Original code
(  298)       !!   2.0  !  02-12  (G. Madec)  F90: Free form and module
(  299)       !!----------------------------------------------------------------------
(  300)       !! * Modules used
(  301)       USE cpl_oce       ! coupled ocean-atmosphere variables
(  302) 
(  303)       !! * Arguments
(  304)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step index
(  305) 
(  306)       !! * Local declarations
(  307)       INTEGER  ::   ji, jj, jf         ! dummy loop indices
(  308)       REAL(wp) ::   ztgel,          &  ! temporary scalars
(  309)          zice, zhemis, zqrp, zqri,  &  !    "         "
(  310)          zq, zqi, zerp, ze, zei, zro   !    "         "
(  311)       !!----------------------------------------------------------------------
(  312) 
(  313)       ! Compute fluxes
(  314)       ! --------------
(  315) 
(  316)       DO jj = 1, jpj
(  317)          DO ji = 1, jpi
(  318) 
(  319)             ztgel = fzptn(ji,jj)   ! local freezing temperature
(  320) 
(  321)             ! opa model ice freeze()
(  322) 
(  323)             zice = tmask(ji,jj,1)
(  324)             IF( tn(ji,jj,1) >=  ztgel ) zice = 0.
(  325)             freeze(ji,jj) = zice
(  326) 
(  327)             ! hemisphere indicator (=1 north, =-1 south)
(  328)             
(  329)             zhemis = float(isign(1, mjg(jj)-(jpjglo/2+1)))
(  330)             
(  331)             ! a) net downward radiative flux qsr()
(  332)             ! - AGCM qsrc if no ice
(  333)             ! - zero under ice (zice=1)
(  334) 
(  335)             qsr(ji,jj) = (1.-zice)*qsrc(ji,jj)*tmask(ji,jj,1)
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 7

(  337)             ! b) heat flux damping term qrp()
(  338)             ! - no damping if no  ice      (zice=0) 
(  339)             ! - gamma*min(0,t-tgel) if ice (zice=1)
(  340) 
(  341)             zqrp = 0.
(  342)             zqri = dqdt0*MIN( 0., tb(ji,jj,1)-ztgel )
(  343)             qrp(ji,jj) = ( ( 1. - zice ) * zqrp + zice * zqri ) * tmask(ji,jj,1)
(  344) 
(  345) 
(  346)             ! c) net downward heat flux q() = q0 + qrp()
(  347)             ! for q0
(  348)             ! - AGCM qc if no  ice (zice=0)
(  349)             ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if ice (zice=1)
(  350)             zq  = qc(ji,jj)
(  351)             zqi = -3. + zhemis
(  352)             qt(ji,jj) = ( (1.-zice) * zq + zice * zqi ) * tmask(ji,jj,1) + qrp(ji,jj)
(  353)             
(  354)             ! d) water flux damping term erp()
(  355)             ! - no damping
(  356)             zerp = 0.
(  357)             erp(ji,jj) = zerp
(  358)             
(  359)             ! e) net upward water flux e() = eo + runoff() + erp()
(  360)             ! for e0
(  361)             ! - AGCM if no ice (zice=0)
(  362)             ! - 1.mm/day if climatological and opa ice (zice=1)
(  363)             ze  = ec(ji,jj)
(  364)             zei = 1./rday
(  365)             zro = runoff(ji,jj)
(  366)             emp(ji,jj) = ( ( 1. - zice ) *  ze + zice * zei + zro ) * tmask(ji,jj,1) + erp(ji,jj)
(  367)             
(  368)             ! f) net upward water flux for the salinity surface 
(  369)             !    boundary condition
(  370)             emps(:,:) = emp(:,:)
(  371) 
(  372)          END DO
(  373)       END DO
(  374) 
(  375)    END SUBROUTINE oce_sbc
(  376) 
(  377) # elif defined key_flx_bulk_monthly || defined key_flx_bulk_daily || defined key_flx_forced_daily
(  378)       !!-------------------------------------------------------------------------
(  379)       !!   'key_flx_bulk_monthly' or 'key_flx_bulk_daily' or        bulk formulea
(  380)       !!   'key_flx_forced_daily'                                or no bulk case 
(  381)       !!-------------------------------------------------------------------------
(  382) 
(  383)    SUBROUTINE oce_sbc( kt )
(  384)       !!---------------------------------------------------------------------
(  385)       !!                   ***  ROUTINE oce_sbc  ***
(  386)       !!                    
(  387)       !! ** Purpose :   Ocean surface boundary conditions in forced mode
(  388)       !!      using either flux or bulk formulation.
(  389)       !!
(  390)       !! History :
(  391)       !!   1.0  !  99-11  (M. Imbard)  Original code
(  392)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(  393)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(  394)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 8

(  395)       !!----------------------------------------------------------------------
(  396)       !! * Modules used
(  397)       USE daymod                 ! calendar
(  398) #if ! defined key_dtasst
(  399)       USE dtasst, ONLY : rclice  ! sea surface temperature data
(  400) #endif
(  401) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  402)       USE blk_oce                ! bulk variables
(  403) #endif
(  404) #if defined key_flx_forced_daily
(  405)       USE flx_oce                ! sea-ice/ocean forcings variables
(  406) #endif
(  407) 
(  408)       !! * arguments
(  409)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  410) 
(  411)       !! * local declarations
(  412)       INTEGER ::   ji, jj        ! dummy loop arguments
(  413)       INTEGER ::   i15, ifreq             !      
(  414)       REAL(wp) ::  zxy
(  415)       REAL(wp) ::  zsice, zqri, zqrp, ztdta, zqrj
(  416)       REAL(wp) ::  zq, zqi, zhemis
(  417)       REAL(wp), DIMENSION(jpi,jpj) :: zeri, zerps, ziclim
(  418)       REAL(wp), DIMENSION(jpi,jpj) :: zqt, zqsr, zemp  
(  419)       !!----------------------------------------------------------------------
(  420)  
(  421)       ! 1. initialization to zero at kt = nit000
(  422)       ! ---------------------------------------
(  423)       
(  424)       IF( kt == nit000 ) THEN     
(  425)          qsr    (:,:) = 0.e0
(  426)          freeze (:,:) = 0.e0
(  427)          qt     (:,:) = 0.e0
(  428)          qrp    (:,:) = 0.e0
(  429)          emp    (:,:) = 0.e0
(  430)          emps   (:,:) = 0.e0
(  431)          erp    (:,:) = 0.e0
(  432) #if ! defined key_dynspg_rl 
(  433)          dmp    (:,:) = 0.e0
(  434) #endif
(  435)       ENDIF
(  436) 
(  437) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(  438)       ifreq      = nfbulk
(  439)       zqt (:,:)  = qsr_oce(:,:) + qnsr_oce(:,:)
(  440)       zqsr(:,:)  = qsr_oce(:,:)
(  441)       zemp(:,:)  = evap(:,:) - tprecip(:,:)
(  442) #endif 
(  443) 
(  444) #if defined key_flx_forced_daily
(  445)       ifreq      = 1
(  446)       zqt (:,:)  = p_qt (:,:)
(  447)       zqsr(:,:)  = p_qsr(:,:)
(  448)       zemp(:,:)  = p_emp(:,:)
(  449) #endif 
(  450) 
(  451)       IF( MOD( kt-1, ifreq) == 0 ) THEN
(  452)          ! Computation of internal and evaporation damping terms       






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 9

(  453)          CALL oce_sbc_dmp
(  454) 
(  455)          zsice = - 0.04 / 0.8    ! ratio of isohaline compressibility over isotherme compressibility 
(  456)                                  ! ( d rho / dt ) / ( d rho / ds )      ( s = 34, t = -1.8 )
(  457)          ! Flux computation
(  458)          DO jj = 1, jpj
(  459)             DO ji = 1, jpi      
(  460)                ! climatological ice 
(  461)                ziclim(ji,jj) = FLOAT( NINT( rclice(ji,jj,1) ) )
(  462) 
(  463)                ! avoid surfreezing point            
(  464)                tn(ji,jj,1) = MAX( tn(ji,jj,1), fzptn(ji,jj) )
(  465) 
(  466)                ! hemisphere indicator (=1 north, =-1 south)           
(  467)                zhemis = FLOAT( isign(1, mjg(jj) - (jpjdta/2+1) ) )
(  468) 
(  469)                ! restoring temperature (ztdta >= to local freezing temperature)            
(  470) #if defined key_dtasst
(  471)                ztdta = MAX( sst(ji,jj),    fzptn(ji,jj) )
(  472) #else
(  473)                ztdta = MAX( t_dta(ji,jj,1), fzptn(ji,jj) )
(  474) #endif
(  475) 
(  476)                ! a) net downward radiative flux qsr()           
(  477)                qsr(ji,jj) = (1.-ziclim(ji,jj)) * zqsr(ji,jj) * tmask(ji,jj,1)
(  478) 
(  479)                ! b) heat flux damping term qrp()
(  480)                ! - gamma*(t-tlevitus) if no  climatological ice (ziclim=0)
(  481)                ! - gamma*(t-(tgel-1.))  if climatological ice and no opa ice   (ziclim=1 zicopa=0)
(  482)                ! - gamma*min(0,t-tgel) if climatological and opa ice (ziclim=1 zicopa=1)
(  483) 
(  484)                zqri = dqdt0 * ( tb(ji,jj,1) - ( fzptn(ji,jj) - 1.) )
(  485)                zqrj = dqdt0 * MIN( 0., tb(ji,jj,1) - fzptn(ji,jj) )
(  486) 
(  487)                qrp(ji,jj) =  ( ziclim(ji,jj) * ( (1 - freeze(ji,jj)) * zqri    &
(  488)                  &                                  + freeze(ji,jj)  * zqrj ) ) * tmask(ji,jj,1)
(  489) 
(  490) #if ! defined key_flx_bulk_monthly || ! defined key_flx_bulk_daily
(  491)                zqrp = dqdt0 * ( tb(ji,jj,1) - ztdta )
(  492)                qrp(ji,jj) = qrp(ji,jj) + (1. - ziclim(ji,jj)) * zqrp
(  493) # endif
(  494) 
(  495)                ! c) net downward heat flux q() = q0 + qrp()
(  496)                ! for q0
(  497)                ! - ECMWF fluxes if no climatological ice      (ziclim=0)
(  498)                ! - qrp if climatological ice and no opa ice   (ziclim=1 zicopa=0)
(  499)                ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if climatological and opa ice 
(  500)                !                                              (ziclim=1 zicopa=1)
(  501)                zq  = zqt(ji,jj)
(  502)                zqi = -3. + zhemis
(  503)                qt (ji,jj) = ( (1.-ziclim(ji,jj)) * zq   &
(  504)                   +ziclim(ji,jj)  * freeze(ji,jj) * zqi )   &
(  505)                   * tmask(ji,jj,1)   &
(  506)                   + qrp(ji,jj)
(  507) 
(  508)             END DO
(  509)          END DO
(  510) 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 10

(  511) #if ! defined key_dynspg_rl 
(  512)          ! Free-surface
(  513) 
(  514)          ! Water flux for zero buoyancy flux if no opa ice and ice clim
(  515)          zeri(:,:) = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
(  516)          zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
(  517) 
(  518)          ! Contribution to sea level:
(  519)          ! net upward water flux emp() = e-p + runoff() + erp() + dmp + empold
(  520)          emp (:,:) = zemp(:,:)     &   ! e-p data
(  521)             &      + runoff(:,:)   &   ! runoff data
(  522)             &      + erp(:,:)      &   ! restoring term to SSS data
(  523)             &      + dmp(:,:)      &   ! freshwater flux associated with internal damping
(  524)             &      + empold            ! domain averaged annual mean correction
(  525) 
(  526)          ! Contribution to salinity:
(  527)          ! net upward water flux emps() = e-p + runoff() + erp() + zerps + empold
(  528)          emps(:,:) = zemp(:,:)     &
(  529)             &      + runoff(:,:)   &
(  530)             &      + erp(:,:)      &
(  531)             &      + zerps(:,:)    &
(  532)             &      + empold 
(  533) #else
(  534)          ! Rigid-lid (emp=emps=E-P-R+Erp)
(  535)          ! freshwater flux
(  536)          zeri(:,:)  = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
(  537)          zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
(  538)          emps (:,:) = zemp(:,:)     &
(  539)             &       + runoff(:,:)   &
(  540)             &       + erp(:,:)      &
(  541)             &       + zerps(:,:)
(  542)          emp (:,:) = emps(:,:)
(  543) #endif  
(  544) 
(  545)          ! Boundary condition on emp for free surface option
(  546)          ! -------------------------------------------------
(  547)          CALL lbc_lnk( emp, 'T', 1. )
(  548)       
(  549)       ENDIF
(  550) 
(  551)    END SUBROUTINE oce_sbc
(  552) 
(  553) !----------------------------------------------------------------!byoung
(  554) # elif defined key_flx_forced_monthly
(  555)       !!-------------------------------------------------------------------------
(  556)       !!   'key_flx_bulk_monthly' or 'key_flx_bulk_daily' or        bulk formulea
(  557)       !!   'key_flx_forced_daily'                                or no bulk case
(  558)       !!-------------------------------------------------------------------------
(  559) 
(  560)    SUBROUTINE oce_sbc( kt )
(  561)       !!---------------------------------------------------------------------
(  562)       !!                   ***  ROUTINE oce_sbc  ***
(  563)       !!
(  564)       !! ** Purpose :   Ocean surface boundary conditions in forced mode
(  565)       !!      using either flux or bulk formulation.
(  566)       !!
(  567)       !! History :
(  568)       !!   1.0  !  99-11  (M. Imbard)  Original code






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 11

(  569)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(  570)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(  571)       !!----------------------------------------------------------------------
(  572)       !! * Modules used
(  573)       USE daymod                 ! calendar
(  574) #if ! defined key_dtasst
(  575)       USE dtasst, ONLY : rclice  ! sea surface temperature data
(  576) #endif
(  577) #if defined key_flx_forced_monthly
(  578)       USE flx_oce                ! sea-ice/ocean forcings variables
(  579) #endif
(  580) 
(  581)       !! * arguments
(  582)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  583) 
(  584)       !! * local declarations
(  585)       INTEGER ::   ji, jj        ! dummy loop arguments
(  586)       INTEGER ::   i15, ifreq             !
(  587)       REAL(wp) ::  zxy
(  588)       REAL(wp) ::  zsice, zqri, zqrp, ztdta, zqrj
(  589)       REAL(wp) ::  zq, zqi, zhemis, ztrp
(  590)       REAL(wp), DIMENSION(jpi,jpj) :: zeri, zerps, ziclim
(  591)       REAL(wp), DIMENSION(jpi,jpj) :: zqt, zqsr, zemp
(  592)       !!----------------------------------------------------------------------
(  593) 
(  594)       ! 1. initialization to zero at kt = nit000
(  595)       ! ---------------------------------------
(  596) 
(  597)       IF( kt == nit000 ) THEN
(  598)          qsr    (:,:) = 0.e0
(  599)          freeze (:,:) = 0.e0
(  600)          qt     (:,:) = 0.e0
(  601)          qrp    (:,:) = 0.e0
(  602)          emp    (:,:) = 0.e0
(  603)          emps   (:,:) = 0.e0
(  604)          erp    (:,:) = 0.e0
(  605) #if ! defined key_dynspg_rl
(  606)          dmp    (:,:) = 0.e0
(  607) #endif
(  608)       ENDIF
(  609) 
(  610) 
(  611) #if defined key_flx_forced_monthly
(  612)       ifreq      = 1
(  613)       zqt (:,:)  = p_bqt (:,:)
(  614)       zqsr(:,:)  = p_bqsr(:,:)
(  615)       zemp(:,:)  = p_bemp(:,:)
(  616) #endif
(  617) 
(  618)       IF( MOD( kt-1, ifreq) == 0 ) THEN
(  619)          ! Computation of internal and evaporation damping terms
(  620)          CALL oce_sbc_dmp
(  621) 
(  622)          ztrp = -40.             ! restoring terme for temperature (w/m2/k)
(  623)          zsice = - 0.04 / 0.8    ! ratio of isohaline compressibility over isotherme compressibility
(  624)                                  ! ( d rho / dt ) / ( d rho / ds )      ( s = 34, t = -1.8 )
(  625)          ! Flux computation
(  626)          DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 12

(  627)             DO ji = 1, jpi
(  628)                ! climatological ice
(  629)                ziclim(ji,jj) = FLOAT( NINT( rclice(ji,jj,1) ) )
(  630) 
(  631)                ! avoid surfreezing point
(  632)                tn(ji,jj,1) = MAX( tn(ji,jj,1), fzptn(ji,jj) )
(  633) 
(  634)                ! hemisphere indicator (=1 north, =-1 south)
(  635)                zhemis = FLOAT( isign(1, mjg(jj) - (jpjdta/2+1) ) )
(  636) 
(  637)                ! restoring temperature (ztdta >= to local freezing temperature)
(  638) #if defined key_dtasst
(  639)                ztdta = MAX( sst(ji,jj),    fzptn(ji,jj) )
(  640) #else
(  641)                ztdta = MAX( t_dta(ji,jj,1), fzptn(ji,jj) )
(  642) #endif
(  643) 
(  644)                ! a) net downward radiative flux qsr()
(  645)                qsr(ji,jj) = (1.-ziclim(ji,jj)) * zqsr(ji,jj) * tmask(ji,jj,1)
(  646) 
(  647)                ! b) heat flux damping term qrp()
(  648)                ! - gamma*(t-tlevitus) if no  climatological ice (ziclim=0)
(  649)                ! - gamma*(t-(tgel-1.))  if climatological ice and no opa ice   (ziclim=1 zicopa=0)
(  650)                ! - gamma*min(0,t-tgel) if climatological and opa ice (ziclim=1 zicopa=1)
(  651) 
(  652)                zqri = ztrp * ( tb(ji,jj,1) - ( fzptn(ji,jj) - 1.) )
(  653)                zqrj = ztrp * MIN( 0., tb(ji,jj,1) - fzptn(ji,jj) )
(  654) 
(  655)                qrp(ji,jj) =  ( ziclim(ji,jj) * ( (1 - freeze(ji,jj)) * zqri    &
(  656)                  &                                  + freeze(ji,jj)  * zqrj ) ) * tmask(ji,jj,1)
(  657) 
(  658) 
(  659)                ! c) net downward heat flux q() = q0 + qrp()
(  660)                ! for q0
(  661)                ! - ECMWF fluxes if no climatological ice      (ziclim=0)
(  662)                ! - qrp if climatological ice and no opa ice   (ziclim=1 zicopa=0)
(  663)                ! - -2 watt/m2 (arctic) or -4 watt/m2 (antarctic) if climatological and opa ice
(  664)                !                                              (ziclim=1 zicopa=1)
(  665)                zq  = zqt(ji,jj)
(  666)                zqi = -3. + zhemis
(  667)                qt (ji,jj) = ( (1.-ziclim(ji,jj)) * zq   &
(  668)                   +ziclim(ji,jj)  * freeze(ji,jj) * zqi )   &
(  669)                   * tmask(ji,jj,1)   &
(  670)                   + qrp(ji,jj)
(  671) 
(  672)             END DO
(  673)          END DO
(  674) 
(  675) #if ! defined key_dynspg_rl
(  676)          ! Free-surface
(  677) 
(  678)          ! Water flux for zero buoyancy flux if no opa ice and ice clim
(  679)          zeri(:,:) = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
(  680)          zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
(  681) 
(  682)          ! Contribution to sea level:
(  683)          ! net upward water flux emp() = e-p + runoff() + erp() + dmp + empold
(  684)          emp (:,:) = zemp(:,:)     &   ! e-p data






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 13

(  685)             &      + runoff(:,:)   &   ! runoff data
(  686)             &      + erp(:,:)      &   ! restoring term to SSS data
(  687)             &      + dmp(:,:)      &   ! freshwater flux associated with internal damping
(  688)             &      + empold            ! domain averaged annual mean correction
(  689) 
(  690)          ! Contribution to salinity:
(  691)          ! net upward water flux emps() = e-p + runoff() + erp() + zerps + empold
(  692)          emps(:,:) = zemp(:,:)     &
(  693)             &      + runoff(:,:)   &
(  694)             &      + erp(:,:)      &
(  695)             &      + zerps(:,:)    &
(  696)             &      + empold
(  697) #else
(  698)          ! Rigid-lid (emp=emps=E-P-R+Erp)
(  699)          ! freshwater flux
(  700)          zeri(:,:)  = -zsice * qrp(:,:) * ro0cpr * rauw / 34.0
(  701)          zerps(:,:) = ziclim(:,:) * ( (1-freeze(:,:)) * zeri(:,:) )
(  702)          emps (:,:) = zemp(:,:)     &
(  703)             &       + runoff(:,:)   &
(  704)             &       + erp(:,:)      &
(  705)             &       + zerps(:,:)
(  706)          emp (:,:) = emps(:,:)
(  707) #endif
(  708) 
(  709)          ! Boundary condition on emp for free surface option
(  710)          ! -------------------------------------------------
(  711)          CALL lbc_lnk( emp, 'T', 1. )
(  712) 
(  713)       ENDIF
(  714) 
(  715)    END SUBROUTINE oce_sbc
(  716) 
(  717) 
(  718) 
(  719) 
(  720) 
(  721) !--------------------------------BIO WZL END
(  722) 
(  723) # else
(  724)       !!----------------------------------------------------------------------
(  725)       !!   Default option :                                 Analytical forcing
(  726)       !!----------------------------------------------------------------------
(  727) 
(  728)    SUBROUTINE oce_sbc( kt )
(  729)       !!---------------------------------------------------------------------
(  730)       !!                    ***  ROUTINE oce_sbc  ***
(  731)       !!              
(  732)       !! ** Purpose :   provide the thermohaline fluxes (heat and freshwater)
(  733)       !!                to the ocean at each time step.
(  734)       !!
(  735)       !! ** Method  :   Constant surface fluxes (read in namelist (namflx))
(  736)       !!
(  737)       !! ** Action  : - qt, qsr, emp, emps, qrp, erp
(  738)       !!
(  739)       !! History :
(  740)       !!        !  91-03  ()  Original code
(  741)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  742)       !!   9.0  !  04-05  (A. Koch-Larrouy) Add Gyre configuration 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 14

(  743)       !!----------------------------------------------------------------------
(  744)       !! * Modules used
(  745)       USE flxrnf                       ! ocean runoffs
(  746)       USE daymod, ONLY : nyear         ! calendar
(  747)       USE dtasss                       ! sea surface salinity data
(  748) 
(  749)       !! * arguments
(  750)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  751) 
(  752)       !! * local declarations
(  753)       REAL(wp) ::                   & !!! surface fluxes namelist (namflx)
(  754)          q0   = 0.e0,               &  ! net heat flux
(  755)          qsr0 = 0.e0,               &  ! solar heat flux
(  756)          emp0 = 0.e0                   ! net freshwater flux
(  757)       REAL(wp) ::   zsrp,           &
(  758)          zemp_S, zemp_N, zemp_sais, &
(  759)          zTstar, zcos_sais1, zconv, &
(  760)          zcos_sais2
(  761)       REAL(wp) ::           &
(  762)          zsumemp,           &          ! tampon used for the emp sum
(  763)          zsurf,             &          ! tampon used for the domain sum
(  764)          ztime,             &          ! time in hour
(  765)          ztimemax1, ztimemin1, &       ! 21th june,   and 21th december if date0 = 1st january
(  766)          ztimemax2, ztimemin2          ! 21th august, and 21th february if date0 = 1st january
(  767)       REAL(wp), DIMENSION(jpi,jpj) :: t_star
(  768)       INTEGER  ::   ji, jj             ! dummy loop indices
(  769) 
(  770)       INTEGER  ::           &
(  771)          zyear0,            &          ! initial year
(  772)          zmonth0,           &          ! initial month
(  773)          zday0,             &          ! initial day
(  774)          zday_year0                    ! initial day since january 1st
(  775) 
(  776)       NAMELIST/namflx/ q0, qsr0, emp0
(  777)       !!---------------------------------------------------------------------
(  778) 
(  779)       !same temperature, E-P as in HAZELEGER 2000
(  780) 
(  781)          ! Constant surface fluxes
(  782) 
(  783)          IF( kt == nit000 ) THEN
(  784)             IF(lwp) THEN
(  785)                WRITE(numout,*)' '
(  786)                WRITE(numout,*)' ocesbc  : Constant surface fluxes read in namelist'
(  787)                WRITE(numout,*)' ~~~~~~~ '
(  788)                WRITE(numout,*)'           Namelist namflx: set the constant flux values'
(  789)                WRITE(numout,*)'              net heat flux          q0   = ', q0  , ' W/m2'
(  790)                WRITE(numout,*)'              solar heat flux        qsr0 = ', qsr0, ' W/m2'
(  791)                WRITE(numout,*)'              net heat flux          emp0 = ', emp0, ' W/m2'
(  792)             ENDIF
(  793) 
(  794)             qt    (:,:) = q0
(  795)             qsr   (:,:) = qsr0
(  796)             emp   (:,:) = emp0
(  797)             emps  (:,:) = emp0
(  798)             qrp   (:,:) = 0.e0
(  799)             erp   (:,:) = 0.e0
(  800)    






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 15

(  801)             runoff(:,:) = 0.e0
(  802)          ENDIF
(  803) 
(  804) 
(  805)    END SUBROUTINE oce_sbc
(  806) 
(  807) # endif
(  808) #endif
(  809) 
(  810) #if defined key_dtasal
(  811)    !!----------------------------------------------------------------------
(  812)    !!   'key_dtasal'                                          salinity data
(  813)    !!----------------------------------------------------------------------
(  814)    SUBROUTINE oce_sbc_dmp
(  815)       !!---------------------------------------------------------------------
(  816)       !!                   ***  ROUTINE oce_sbc_dmp  ***
(  817)       !!                    
(  818)       !! ** Purpose : Computation of internal and evaporation damping terms 
(  819)       !!        for ocean surface boundary conditions 
(  820)       !!
(  821)       !! History :
(  822)       !!   9.0  !  04-01  (G. Madec, C. Ethe)  Original code
(  823)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  824)       !!----------------------------------------------------------------------
(  825)       !! * Local declarations
(  826)       INTEGER ::   ji, jj                   ! dummy loop indices
(  827)       REAL(wp), DIMENSION(jpi,jpj)  :: zsss, zfreeze
(  828)       REAL(wp) ::   zerp, zsrp
(  829)       CHARACTER (len=71) :: charout
(  830) #if ! defined key_dynspg_rl
(  831)       REAL(wp) ::   zwei
(  832)       REAL(wp) ::   zerpplus(jpi,jpj), zerpminus(jpi,jpj)
(  833)       REAL(wp) ::   zplus, zminus, zadefi
(  834) # if defined key_tradmp
(  835)       INTEGER jk
(  836)       REAL(wp), DIMENSION(jpi,jpj) ::   zstrdmp
(  837) # endif
(  838) #endif
(  839)       !!----------------------------------------------------------------------
(  840) 
(  841) #if defined key_ice_lim
(  842)       ! sea ice indicator (1 or 0)
(  843)       DO jj = 1, jpj
(  844)          DO ji = 1, jpi
(  845)             freezn(ji,jj) = MAX(0., SIGN(1., freeze(ji,jj)-rsmall) )
(  846)          END DO
(  847)       END DO
(  848)       zsss   (:,:) = sss_io(:,:)
(  849)       zfreeze(:,:) = freezn(:,:)
(  850) #else
(  851)       zsss   (:,:) = sb    (:,:,1)
(  852)       zfreeze(:,:) = freeze(:,:)
(  853) #endif
(  854) 
(  855)       ! Initialisation
(  856)       ! --------------
(  857)       ! Restoring coefficients on SST and SSS   
(  858)       zsrp = dqdt0 * ro0cpr * rauw   ! (Kg/m2/s) 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 16

(  859) 
(  860) #if ! defined key_dynspg_rl 
(  861)       ! Free-surface
(  862)          
(  863)       ! Internal damping
(  864) # if defined key_tradmp
(  865)       ! Vertical mean of dampind trend (computed in tradmp module)
(  866)       zstrdmp(:,:) = 0.e0
(  867)       DO jk = 1, jpk
(  868)          zstrdmp(:,:) = zstrdmp(:,:) + strdmp(:,:,jk) * fse3t(:,:,jk)
(  869)       END DO
(  870)       ! volume flux associated to internal damping to climatology
(  871)       dmp(:,:) = zstrdmp(:,:) * rauw / ( zsss(:,:) + 1.e-20 )
(  872) # else
(  873)       dmp(:,:) = 0.e0            ! No internal damping
(  874) # endif
(  875)       
(  876)       !   evaporation damping term ( Surface restoring )
(  877)       zerpplus (:,:) = 0.e0
(  878)       zerpminus(:,:) = 0.e0
(  879)       zplus  =  15. / rday
(  880)       zminus = -15. / rday
(  881)       
(  882)       DO jj = 1, jpj
(  883)          DO ji = 1, jpi
(  884)             zerp = ( 1. - 2.*upsrnfh(ji,jj) ) * zsrp   &
(  885)                & * ( zsss(ji,jj) - s_dta(ji,jj,1) )     &
(  886)                & / ( zsss(ji,jj) + 1.e-20        )
(  887)             
(  888)             zerp = MIN( zerp, zplus  )
(  889)             zerp = MAX( zerp, zminus )
(  890)             erp(ji,jj) = zerp
(  891)             zerpplus (ji,jj) = MAX( erp(ji,jj), 0.e0 )
(  892)             zerpminus(ji,jj) = MIN( erp(ji,jj), 0.e0 )
(  893)          END DO
(  894)       END DO
(  895) 
(  896)       aplus  = 0.e0
(  897)       aminus = 0.e0
(  898)       DO jj = 1, jpj
(  899)          DO ji = 1, jpi
(  900)             zwei   = e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
(  901)             aplus  = aplus  + zerpplus (ji,jj) * zwei
(  902)             aminus = aminus - zerpminus(ji,jj) * zwei
(  903)          END DO
(  904)       END DO
(  905)       IF( lk_mpp )   CALL mpp_sum( aplus  )   ! sums over the global domain
(  906)       IF( lk_mpp )   CALL mpp_sum( aminus )
(  907) 
(  908)       IF(ln_ctl)   THEN
(  909)          WRITE(charout,FMT="('oce_sbc_dmp : a+ = ',D23.16, ' a- = ',D23.16)") aplus, aminus
(  910)          CALL prt_ctl_info(charout)
(  911)       ENDIF
(  912) 
(  913)       zadefi = MIN( aplus, aminus )
(  914)       IF( zadefi == 0.e0 ) THEN 
(  915)          erp(:,:) = 0.e0
(  916)       ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 17

(  917)          erp(:,:) = zadefi * ( zerpplus(:,:) / aplus + zerpminus(:,:) / aminus )
(  918)       ENDIF
(  919)       erp(:,:)=0 !byoung for 008
(  920) #else
(  921)       ! Rigid-lid (emp=emps=E-P-R+Erp)
(  922)       
(  923)       erp(:,:) = ( 1. - zfreeze(:,:) ) * zsrp    &   ! surface restoring term
(  924)          &     * ( zsss(:,:) - s_dta(:,:,1) )     &
(  925)          &     / ( zsss(:,:) + 1.e-20      )
(  926) #endif
(  927) 
(  928)    END SUBROUTINE oce_sbc_dmp
(  929) 
(  930) #else
(  931)    !!----------------------------------------------------------------------
(  932)    !!   Dummy routine                                      NO salinity data
(  933)    !!----------------------------------------------------------------------
(  934)    SUBROUTINE oce_sbc_dmp         ! Dummy routine
(  935) !      WRITE(*,*) 'oce_sbc_dmp: you should not have seen that print! error?'
(  936)    END SUBROUTINE oce_sbc_dmp
(  937) #endif
(  938) 
(  939)    !!======================================================================
(  940) END MODULE ocesbc








































PGF90 (Version     10.5)          09/23/2010  09:49:14      page 18

(    1) # 1 "ocesbc.F90"
(    1) !!DB: 2009.08.31 -- Eliminated GYRE config
(    2) MODULE ocesbc
(    3)    !!======================================================================
(    4)    !!                     ***  MODULE  ocesbc  ***
(    5)    !!                     Ocean surface boundary conditions
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   oce_sbc     : ???
(   10)    !!   oce_sbc_dmp : ???
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce            ! dynamics and tracers variables
(   14)    USE dom_oce        ! ocean space domain variables
(   15)    USE cpl_oce        ! coupled ocean-atmosphere variables
(   16)    USE ice_oce        ! sea-ice variable
(   17)    USE blk_oce        ! bulk variables
(   18)    USE flx_oce        ! sea-ice/ocean forcings variables
(   19)    USE phycst         ! Define parameters for the routines
(   20)    USE taumod         ! surface stress forcing
(   21)    USE flxmod         ! thermohaline fluxes
(   22)    USE flxrnf         ! runoffs forcing
(   23)    USE tradmp         ! damping salinity trend
(   24)    USE dtatem         ! ocean temperature data
(   25)    USE dtasal         ! ocean salinity data
(   26)    USE ocfzpt         ! surface ocean freezing point
(   27)    USE lbclnk         ! ocean lateral boundary condition
(   28)    USE lib_mpp        ! distribued memory computing library
(   29)    USE in_out_manager ! I/O manager
(   30)    USE prtctl         ! Print control
(   31) 
(   32)    IMPLICIT NONE
(   33)    PRIVATE
(   34) 
(   35)    !! * Accessibility
(   36)    PUBLIC oce_sbc    ! routine called by step
(   37) 
(   38)    !! * Shared module variables
(   39)    REAL(wp), PUBLIC ::   &  !:
(   40)       aplus, aminus,     &  !:
(   41)       empold = 0.e0         !: current year freshwater budget correction
(   42)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   43)       qt  ,         &  !: total surface heat flux (w/m2)
(   44)       qsr ,         &  !: solar radiation (w/m2)
(   45)       emp ,         &  !: evaporation minus precipitation (kg/m2/s = mm/s)
(   46)       emps,         &  !: evaporation - precipitation (free surface)
(   47)       qrp ,         &  !: heat flux damping (w/m2)
(   48)       erp              !: evaporation damping (kg/m2/s = mm/s)
(   49) # 50
(   50)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   51)       dmp              !: internal dampind term
(   52) 
(   53) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 19

(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 55 "ocesbc.F90"
(   55) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 56 "ocesbc.F90"
(   56) 
(   57) # 57
(   57)    !!----------------------------------------------------------------------
(   58)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   59)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/ocesbc.F90,v 1.16 2006/04/19 14:43:16 opalod Exp $ 
(   60)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   61)    !!----------------------------------------------------------------------
(   62) CONTAINS
(   63) 
(   64) # 65
(   65)    !!----------------------------------------------------------------------
(   66)    !!   'key_ice_lim' :                                   LIM sea-ice model
(   67)    !!----------------------------------------------------------------------
(   68) # 178
(  178)       !!----------------------------------------------------------------------
(  179)       !!   'key_ice_lim'                              with  LIM sea-ice model
(  180)       !!----------------------------------------------------------------------
(  181) 
(  182)    SUBROUTINE oce_sbc( kt )
(  183)       !!---------------------------------------------------------------------
(  184)       !!                   ***  ROUTINE oce_sbc  ***
(  185)       !!                    
(  186)       !! ** Purpose : - Ocean surface boundary conditions with LIM sea-ice
(  187)       !!        model in forced mode using bulk formulea
(  188)       !!
(  189)       !! History :
(  190)       !!   1.0  !  99-11  (M. Imbard)  Original code
(  191)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(  192)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(  193)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  194)       !!----------------------------------------------------------------------
(  195)       !! * arguments
(  196)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  197) 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 20

(  198)       !! * Local declarations
(  199)       INTEGER  ::   ji, jj                   ! dummy loop indices
(  200)       REAL(wp) ::   ztx, ztaux, zty, ztauy
(  201)       !!----------------------------------------------------------------------
(  202) 
(  203)       ! 1. initialization to zero at kt = nit000
(  204)       ! ---------------------------------------
(  205)       
(  206)       IF( kt == nit000 ) THEN     
(  207)          qsr    (:,:) = 0.e0
(  208)          qt     (:,:) = 0.e0
(  209)          qrp    (:,:) = 0.e0
(  210)          emp    (:,:) = 0.e0
(  211)          emps   (:,:) = 0.e0
(  212)          erp    (:,:) = 0.e0
(  213) # 214
(  214)          dmp    (:,:) = 0.e0
(  215) # 216
(  216)       ENDIF
(  217) 
(  218)       IF( MOD( kt-1, nfice ) == 0 ) THEN
(  219) 
(  220)          CALL oce_sbc_dmp       ! Computation of internal and evaporation damping terms       
(  221) 
(  222)          ! Surface Ocean fluxes
(  223)          ! ====================
(  224) 
(  225)          ! Surface heat flux (W/m2)
(  226)          ! -----------------
(  227) 
(  228)          qt  (:,:) = fnsolar(:,:) + fsolar(:,:)     ! non solar heat flux + solar flux
(  229)          qsr (:,:) = fsolar(:,:)                     ! solar flux
(  230) 
(  231) # 232
(  232)          ! total concentration/dilution effect (use on SSS)
(  233)          emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:) + empold
(  234) 
(  235)          ! total volume flux (use on sea-surface height)
(  236)          emp (:,:) = fmass(:,:) -   dmp(:,:) + runoff(:,:) + erp(:,:) + empold      
(  237) 
(  238) # 244
(  244)          ! Surface stress
(  245)          ! --------------
(  246) 
(  247)          ! update the stress beloww sea-ice area
(  248)          DO jj = 1, jpjm1
(  249)             DO ji = 1, jpim1   ! vertor opt.
(  250)                ztx         = MAX( freezn(ji,jj), freezn(ji,jj+1) )   ! ice/ocean indicator at U- and V-points
(  251)                zty         = MAX( freezn(ji,jj), freezn(ji+1,jj) )
(  252)                ztaux       = 0.5 *( ftaux(ji+1,jj) + ftaux(ji+1,jj+1) ) ! ice-ocean stress at U- and V-points
(  253)                ztauy       = 0.5 *( ftauy(ji,jj+1) + ftauy(ji+1,jj+1) )
(  254)                taux(ji,jj) = (1.-ztx) * taux(ji,jj) + ztx * ztaux    ! stress at the ocean surface
(  255)                tauy(ji,jj) = (1.-zty) * tauy(ji,jj) + zty * ztauy
(  256)             END DO
(  257)          END DO
(  258) 
(  259)          ! boundary condition on the stress (taux,tauy)
(  260)          CALL lbc_lnk( taux, 'U', -1. )






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 21

(  261)          CALL lbc_lnk( tauy, 'V', -1. )
(  262) 
(  263)          ! Re-initialization of fluxes
(  264)          sst_io(:,:) = 0.e0
(  265)          sss_io(:,:) = 0.e0
(  266)          u_io  (:,:) = 0.e0
(  267)          v_io  (:,:) = 0.e0
(  268) 
(  269)       ENDIF
(  270) 
(  271)    END SUBROUTINE oce_sbc
(  272) 
(  273) 
(  274) 
(  275) # 811
(  811)    !!----------------------------------------------------------------------
(  812)    !!   'key_dtasal'                                          salinity data
(  813)    !!----------------------------------------------------------------------
(  814)    SUBROUTINE oce_sbc_dmp
(  815)       !!---------------------------------------------------------------------
(  816)       !!                   ***  ROUTINE oce_sbc_dmp  ***
(  817)       !!                    
(  818)       !! ** Purpose : Computation of internal and evaporation damping terms 
(  819)       !!        for ocean surface boundary conditions 
(  820)       !!
(  821)       !! History :
(  822)       !!   9.0  !  04-01  (G. Madec, C. Ethe)  Original code
(  823)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  824)       !!----------------------------------------------------------------------
(  825)       !! * Local declarations
(  826)       INTEGER ::   ji, jj                   ! dummy loop indices
(  827)       REAL(wp), DIMENSION(jpi,jpj)  :: zsss, zfreeze
(  828)       REAL(wp) ::   zerp, zsrp
(  829)       CHARACTER (len=71) :: charout
(  830) # 831
(  831)       REAL(wp) ::   zwei
(  832)       REAL(wp) ::   zerpplus(jpi,jpj), zerpminus(jpi,jpj)
(  833)       REAL(wp) ::   zplus, zminus, zadefi
(  834) # 835
(  835)       INTEGER jk
(  836)       REAL(wp), DIMENSION(jpi,jpj) ::   zstrdmp
(  837) # 839
(  839)       !!----------------------------------------------------------------------
(  840) 
(  841) # 842
(  842)       ! sea ice indicator (1 or 0)
(  843)       DO jj = 1, jpj
(  844)          DO ji = 1, jpi
(  845)             freezn(ji,jj) = MAX(0., SIGN(1., freeze(ji,jj)-rsmall) )
(  846)          END DO
(  847)       END DO
(  848)       zsss   (:,:) = sss_io(:,:)
(  849)       zfreeze(:,:) = freezn(:,:)
(  850) 
(  851) # 855
(  855)       ! Initialisation
(  856)       ! --------------
(  857)       ! Restoring coefficients on SST and SSS   






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 22

(  858)       zsrp = dqdt0 * ro0cpr * rauw   ! (Kg/m2/s) 
(  859) 
(  860) # 861
(  861)       ! Free-surface
(  862)          
(  863)       ! Internal damping
(  864) # 865
(  865)       ! Vertical mean of dampind trend (computed in tradmp module)
(  866)       zstrdmp(:,:) = 0.e0
(  867)       DO jk = 1, jpk
(  868)          zstrdmp(:,:) = zstrdmp(:,:) + strdmp(:,:,jk) * e3t(jk)
(  869)       END DO
(  870)       ! volume flux associated to internal damping to climatology
(  871)       dmp(:,:) = zstrdmp(:,:) * rauw / ( zsss(:,:) + 1.e-20 )
(  872) # 875
(  875)       
(  876)       !   evaporation damping term ( Surface restoring )
(  877)       zerpplus (:,:) = 0.e0
(  878)       zerpminus(:,:) = 0.e0
(  879)       zplus  =  15. / rday
(  880)       zminus = -15. / rday
(  881)       
(  882)       DO jj = 1, jpj
(  883)          DO ji = 1, jpi
(  884)             zerp = ( 1. - 2.*upsrnfh(ji,jj) ) * zsrp   &
(  885)                & * ( zsss(ji,jj) - s_dta(ji,jj,1) )     &
(  886)                & / ( zsss(ji,jj) + 1.e-20        )
(  887)             
(  888)             zerp = MIN( zerp, zplus  )
(  889)             zerp = MAX( zerp, zminus )
(  890)             erp(ji,jj) = zerp
(  891)             zerpplus (ji,jj) = MAX( erp(ji,jj), 0.e0 )
(  892)             zerpminus(ji,jj) = MIN( erp(ji,jj), 0.e0 )
(  893)          END DO
(  894)       END DO
(  895) 
(  896)       aplus  = 0.e0
(  897)       aminus = 0.e0
(  898)       DO jj = 1, jpj
(  899)          DO ji = 1, jpi
(  900)             zwei   = e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
(  901)             aplus  = aplus  + zerpplus (ji,jj) * zwei
(  902)             aminus = aminus - zerpminus(ji,jj) * zwei
(  903)          END DO
(  904)       END DO
(  905)       IF( lk_mpp )   CALL mpp_sum( aplus  )   ! sums over the global domain
(  906)       IF( lk_mpp )   CALL mpp_sum( aminus )
(  907) 
(  908)       IF(ln_ctl)   THEN
(  909)          WRITE(charout,FMT="('oce_sbc_dmp : a+ = ',D23.16, ' a- = ',D23.16)") aplus, aminus
(  910)          CALL prt_ctl_info(charout)
(  911)       ENDIF
(  912) 
(  913)       zadefi = MIN( aplus, aminus )
(  914)       IF( zadefi == 0.e0 ) THEN 
(  915)          erp(:,:) = 0.e0
(  916)       ELSE
(  917)          erp(:,:) = zadefi * ( zerpplus(:,:) / aplus + zerpminus(:,:) / aminus )






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 23

(  918)       ENDIF
(  919)       erp(:,:)=0 !byoung for 008
(  920) 
(  921) # 928
(  928)    END SUBROUTINE oce_sbc_dmp
(  929) 
(  930) 
(  931) # 939
(  939)    !!======================================================================
(  940) END MODULE ocesbc
