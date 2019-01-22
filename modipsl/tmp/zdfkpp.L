


PGF90 (Version     10.5)          09/23/2010  09:49:30      page 1

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

Filename: zdfkpp.F90

(    1) MODULE zdfkpp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfkpp  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the KPP 
(    5)    !!                 turbulent closure parameterization
(    6)    !!=====================================================================
(    7) #if defined key_zdfkpp   ||   defined key_esopa
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_zdfkpp'                                             KPP scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !!   zdf_kpp      : update momentum and tracer Kz from a kpp scheme
(   12)    !!   zdf_kpp_init : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and active tracers 
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE zdf_oce         ! ocean vertical physics
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE phycst          ! physical constants
(   21)    USE taumod          ! surface stress
(   22)    USE eosbn2          ! equation of state
(   23)    USE ocesbc          ! thermohaline fluxes
(   24)    USE zdfddm          ! double diffusion mixing
(   25)    USE prtctl          ! Print control
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC zdf_kpp   ! routine called by step.F90
(   32) 
(   33)    !! * Share Module variables
(   34)    LOGICAL, PUBLIC, PARAMETER ::  &
(   35)       lk_zdfkpp = .TRUE.          !: KPP vertical mixing flag
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   & 
(   37)       ghats                       ! non-local scalar mixing term (gamma/<ws>o)
(   38)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   & 
(   39)       wt0                   ,  &  ! surface temperature flux for non local flux
(   40)       ws0                   ,  &  ! surface salinity flux for non local flux
(   41)       hkpp                        ! boundary layer depht
(   42)    !! * Module variables
(   43)    INTEGER ::                  & !!! ** kpp namelist (namkpp) **
(   44)       nave  =  1            ,  &  ! = 0/1 flag for horizontal average on avt, avmu, avmv
(   45)       navb  =  0                  ! = 0/1 flag for constant or profile background avt
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 2

(   47)    REAL(wp) ::                 & !!! ** Interior Mixing 
(   48)       difmiw  = 1.2e-04_wp  ,  &  ! constant internal wave viscosity (m2/s)
(   49)       difsiw  = 1.2e-05_wp  ,  &  ! constant internal wave diffusivity (m2/s)
(   50)       Riinfty = 0.8_wp      ,  &  ! local Richardson Number limit for shear instability
(   51)       difri   = 5.e-03_wp   ,  &  ! maximum shear mixing at Rig = 0    (m2/s)
(   52)       bvsqcon = -1.e-09_wp  ,  &  ! Brunt-Vaisala squared (1/s**2) for maximum convection
(   53)       difcon  = 1._wp             ! maximum mixing in interior convection (m2/s)
(   54) #if defined key_zdfddm
(   55)    REAL(wp) ::                 & !!! ** Double diffusion Mixing
(   56)       difssf  = 1.e-03_wp   ,  &  ! maximum salt fingering mixing 
(   57)       Rrho0   = 1.9_wp      ,  &  ! limit for salt  fingering mixing 
(   58)       difsdc  = 1.5e-06_wp       ! maximum diffusive convection mixing
(   59) #endif
(   60)    LOGICAL  ::                 &
(   61)       ln_kpprimix  = .TRUE.       ! Shear instability mixing 
(   62) 
(   63)    REAL(wp) ::                 & !!! ** General constants  **
(   64)       epsln   = 1.0e-20_wp   , &  ! a small positive number    
(   65)       pthird  = 1._wp/3._wp  , &  ! 1/3
(   66)       pfourth = 1._wp/4._wp       ! 1/4
(   67) 
(   68)    REAL(wp) ::                 & !!! ** Boundary Layer Turbulence Parameters  **
(   69)       vonk     = 0.4_wp     ,  &  ! von Karman's constant
(   70)       epsilon  = 0.1_wp     ,  &  ! nondimensional extent of the surface layer
(   71)       rconc1   = 5.0_wp     ,  &  ! standard flux profile function parmaeters
(   72)       rconc2   = 16.0_wp    ,  &  !         "        "
(   73)       rconcm   = 8.38_wp    ,  &  ! momentum flux profile fit
(   74)       rconam   = 1.26_wp    ,  &  !         "       "
(   75)       rzetam   = -.20_wp    ,  &  !         "       "       
(   76)       rconcs   = 98.96_wp   ,  &  !  scalar  flux profile fit
(   77)       rconas   = -28.86_wp  ,  &  !         "       "
(   78)       rzetas   = -1.0_wp          !         "       "  
(   79)    REAL(wp) ::                 & !!! ** Boundary Layer Depth Diagnostic  **
(   80)       Ricr     = 0.3_wp     ,  &  ! critical bulk Richardson Number
(   81)       rcekman  = 0.7_wp     ,  &  ! coefficient for ekman depth  
(   82)       rcmonob  = 1.0_wp     ,  &  ! coefficient for Monin-Obukhov depth 
(   83)       rconcv   = 1.7_wp     ,  &  ! ratio of interior buoyancy frequency to buoyancy frequency at entrainment depth
(   84)       hbf      = 1.0_wp     ,  &  ! fraction of bound. layer depth to which absorbed solar 
(   85)       !                           ! rad. and contributes to surf. buo. forcing
(   86)       Vtc                         ! function of rconcv,rconcs,epsilon,vonk,Ricr
(   87)    REAL(wp) ::                 & !!! ** Nonlocal Boundary Layer Mixing **
(   88)       rcstar   = 5.0_wp     ,  &  ! coefficient for convective nonlocal transport
(   89)       rcs      = 1.0e-3_wp  ,  &  ! conversion: mm/s ==> m/s   
(   90)       rcg                         ! non-dimensional coefficient for nonlocal transport
(   91) 
(   92) #if ! defined key_kppcustom
(   93)    REAL(wp), DIMENSION(jpk,jpk) ::   & 
(   94)       del                         ! array for reference mean values of vertical integration 
(   95) #endif
(   96) 
(   97) #if defined key_kpplktb
(   98)    INTEGER, PARAMETER ::       & !!! ** Parameters for lookup table for turbulent velocity scales ** 
(   99)       nilktb   = 892        ,  &  ! number of values for zehat in KPP lookup table
(  100)       njlktb   = 482        ,  &  ! number of values for ustar in KPP lookup table
(  101)       nilktbm1 = nilktb - 1 ,  &  !
(  102)       njlktbm1 = njlktb - 1       !
(  103) 
(  104)    REAL(wp), DIMENSION(nilktb,njlktb) ::  &






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 3

(  105)       wmlktb                ,  &  ! lookup table for the turbulent vertical velocity scale for momentum
(  106)       wslktb                       ! lookup table for the turbulent vertical velocity scale for tracers
(  107) 
(  108)    REAL(wp) ::                 &
(  109)       dehatmin = -4.e-7_wp  ,  &  ! minimum limit for zhat in lookup table (m3/s3) 
(  110)       dehatmax = 0._wp      ,  &  ! maximum limit for zhat in lookup table (m3/s3)
(  111)       ustmin   = 0._wp      ,  &  ! minimum limit for ustar in lookup table (m/s)
(  112)       ustmax   = 0.04_wp    ,  &  ! maximum limit for ustar in lookup table (m/s)    
(  113)       dezehat               ,  &  ! delta zhat in lookup table
(  114)       deustar                     ! delta ustar in lookup table
(  115) #endif
(  116)    REAL(wp), DIMENSION(jpk) :: &  !!! attenuation coef  
(  117)       ratt        
(  118)    !! already defines in module traqsr, but only if the solar radiation penetration is considered
(  119)    REAL(wp) ::                 & !!! * penetrative solar radiation coefficient *
(  120)       rabs = 0.58_wp        ,  &  ! fraction associated with xsi1
(  121)       xsi1 = 0.35_wp        ,  &  ! first depth of extinction 
(  122)       xsi2 = 23.0_wp              ! second depth of extinction 
(  123)       !                           ! (default values: water type Ib) 
(  124) 
(  125)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  126)       etmean                ,  &  ! coefficient used for horizontal smoothing
(  127)       eumean                ,  &  ! at t-, u- and v-points
(  128)       evmean  
(  129)  
(  130) #if defined key_cfg_1d
(  131)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &
(  132)       rig                   ,  &  ! gradient Richardson number
(  133)       rib                   ,  &  ! bulk Richardson number
(  134)       buof                  ,  &  ! buoyancy forcing
(  135)       mols                        ! moning-Obukhov length scale 
(  136)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
(  137)       ekdp                        ! Ekman depth
(  138) #endif
(  139) 
(  140)    INTEGER  ::  &                 !
(  141)       jip = 62 , jjp = 111
(  142) 
(  143)    !! * Substitutions
(  144) #  include "domzgr_substitute.h90"
(  145) #  include "vectopt_loop_substitute.h90"
(  146)    !!----------------------------------------------------------------------
(  147)    !!   OPA 9.0 , LOCEAN-IPSL   (2005)
(  148)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfkpp.F90,v 1.2 2005/09/02 15:45:43 opalod Exp $ 
(  149)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(  150)    !!----------------------------------------------------------------------
(  151) 
(  152) CONTAINS
(  153) 
(  154) 
(  155)    SUBROUTINE zdf_kpp ( kt )
(  156)       !!----------------------------------------------------------------------
(  157)       !!                   ***  ROUTINE zdf_kpp  ***
(  158)       !!
(  159)       !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
(  160)       !!      coefficients and non local mixing using K-profile parameterization
(  161)       !!
(  162)       !! ** Method :   The boundary layer depth hkpp is diagnosed at tracer points






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 4

(  163)       !!      from profiles of buoyancy, and shear, and the surface forcing.
(  164)       !!      Above hbl (sigma=-z/hbl <1) the mixing coefficients are computed from 
(  165)       !!
(  166)       !!                      Kx =  hkpp  Wx(sigma) G(sigma)  
(  167)       !!
(  168)       !!             and the non local term ghat = Cs / Ws(sigma) / hkpp
(  169)       !!      Below hkpp  the coefficients are the sum of mixing due to internal waves
(  170)       !!      shear instability and double diffusion.
(  171)       !!
(  172)       !!      -1- Compute the now interior vertical mixing coefficients at all depths. 
(  173)       !!      -2- Diagnose the boundary layer depth. 
(  174)       !!      -3- Compute the now boundary layer vertical mixing coefficients. 
(  175)       !!      -4- Compute the now vertical eddy vicosity and diffusivity.
(  176)       !!      -5- Smoothing
(  177)       !!
(  178)       !!        N.B. The computation is done from jk=2 to jpkm1 
(  179)       !!             Surface value of avt avmu avmv are set once a time to zero
(  180)       !!             in routine zdf_kpp_init.
(  181)       !!
(  182)       !! ** Action  :   update the non-local terms ghats
(  183)       !!                update avt, avmu, avmv (before vertical eddy coef.)
(  184)       !!
(  185)       !!   References :       
(  186)       !!         Large W.G., Mc Williams J.C. and Doney S.C.              
(  187)       !!         Reviews of Geophysics, 32, 4, November 1994
(  188)       !!         Comments in the code refer to this paper, particularly 
(  189)       !!         the equation number. (LMD94, here after)
(  190)       !!
(  191)       !!   Modifications :
(  192)       !!   --------------
(  193)       !!     original  :  00-03 (LARGE W.G.)
(  194)       !!     additions :  00-04 (CHANUT J.)
(  195)       !!               :  02-06 (MOLINES J.M. for real case CLIPPER)
(  196)       !!               :  03-10 (CHANUT J.)
(  197)       !!
(  198)       !! History :
(  199)       !!   8.1  !  00-03  (J. Chanut)  Original code
(  200)       !!   8.1  !  00-04  (J.M. Molines) for real case CLIPPER  
(  201)       !!   9.0  !  05-01  (C. Ethe) Free form, F90
(  202)       !!----------------------------------------------------------------------
(  203)       !! * Modules used
(  204) #if defined  key_zdfddm
(  205)       USE oce     , zviscos => ua,      &  ! temp. array for viscosities use ua as workspace
(  206)          &          zdiffut => ta,      &  ! temp. array for diffusivities use sa as workspace
(  207)          &          zdiffus => sa          ! temp. array for diffusivities use sa as workspace
(  208) #else
(  209)       USE oce     , zviscos => ua,      &  ! temp. array for viscosities use ua as workspace
(  210)          &          zdiffut => ta          ! temp. array for diffusivities use sa as workspace
(  211) #endif
(  212) 
(  213) 
(  214)       !! * arguments
(  215)       INTEGER, INTENT( in  ) ::         &
(  216)          kt                                ! ocean time step
(  217) 
(  218)       !! * local declarations
(  219)       INTEGER ::                        &
(  220)          ji, jj, jk                        ! dummy loop indices






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 5

(  221)       INTEGER ::                        &  !
(  222)          ikbot, jkmax, jkm1, jkp2          !
(  223) 
(  224)       REAL(wp), DIMENSION(jpi,jpj) ::   & !!! Surface buoyancy forcing, friction velocity
(  225)          zBo, zBosol, zustar              !
(  226)                       !
(  227)       REAL(wp) ::                       &  !
(  228)          ztx, zty, ztau, zflageos,      &  !
(  229)          zstabl, zbuofdep,zucube,       &  !
(  230)          zrhos, zalbet, zbeta,          &  !
(  231)          zthermal, zhalin, zatt1           !
(  232)      
(  233)       REAL(wp) ::                       & !!! Bulk richardson number
(  234)          zref, zt, zs, zh,              &  !
(  235)          zu, zv, zrh,                   &  !
(  236)          zrib, zrinum,                  &  !
(  237)          zdVsq, zVtsq                      !
(  238)       
(  239)       REAL(wp) ::                       & !!! Velocity scales
(  240)          zehat, zeta, zhrib, zsig,      &  !
(  241)          zscale, zwst, zws, zwm
(  242) 
(  243) #if defined key_kpplktb
(  244)       INTEGER ::                        & !!! Lookup table or Analytical functions 
(  245)          il, jl                            !
(  246)       REAL(wp) ::                       &  !
(  247)          ud, zfrac, ufrac,              &  !
(  248)          zwam, zwbm, zwas, zwbs            !
(  249) #else
(  250)      REAL(wp) ::                        &  !
(  251)         zwsun, zwmun,                   & 
(  252)         zcons, zconm, zwcons, zwconm      !
(  253) #endif
(  254)  
(  255)      REAL(wp) ::                       & !!! In situ density
(  256)          zsr, zbw, ze,                  &  !
(  257)          zb, zd, zc, zaw, za,           &  !
(  258)          zb1, za1, zkw, zk0,            &  !
(  259)          zcomp , zrhd, zrhdr,zbvzed       !
(  260) 
(  261) #if ! defined key_kppcustom     
(  262)      !! * local declarations
(  263)       INTEGER ::                        &
(  264)          jm                                ! dummy loop indices
(  265)       REAL(wp) ::                       & !!! Compression terms
(  266)          zr1, zr2, zr3, zr4,            &  !
(  267)          zrhop                             !
(  268) #endif
(  269)  
(  270)       REAL(wp) ::                       &  !
(  271)          zflag, ztemp, zrn2,            &  !
(  272)          zdep21, zdep32, zdep43
(  273) 
(  274)       REAL(wp) ::                       & !!! Interior richardson mixing
(  275)          zdku2, zdkv2, ze3sqr,          &  !
(  276)          zsh2, zri, zfri                   !
(  277) 
(  278)       REAL(wp), DIMENSION(jpi,0:2) ::  &  !!! Moning-Obukov limitation






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 6

(  279)          zmoek
(  280)       REAL(wp), DIMENSION(jpi)     ::  &
(  281)          zmoa, zekman                
(  282)       REAL(wp)                     ::  &
(  283)          zmob, zek
(  284) 
(  285)       REAL(wp), DIMENSION(jpi,4) ::     &  !!! The pipe 
(  286)          zdepw, zdift, zvisc
(  287)       REAL(wp), DIMENSION(jpi,3) ::     & 
(  288)          zdept
(  289)       REAL(wp), DIMENSION(jpi,2) ::     &  
(  290)          zriblk
(  291)       REAL(wp), DIMENSION(jpi,jpk) ::   &  !
(  292)          zmask                          
(  293)       REAL(wp), DIMENSION(jpi) ::       &  ! 
(  294)          zhmax, zria, zhbl 
(  295)       REAL(wp) ::                       &  !
(  296)          zflagri, zflagek,              &  !
(  297)          zflagmo, zflagh, zflagkb          !
(  298)       REAL(wp), DIMENSION(jpi)     ::   & !!! Shape function (G)
(  299)          za2m, za3m, zkmpm,             &
(  300)          za2t, za3t, zkmpt
(  301)       REAL(wp) ::                       &  !
(  302)          zdelta, zdelta2,               &  !
(  303)          zdzup, zdzdn, zdzh,            &  !
(  304)          zvath, zgat1, zdat1,           &  !
(  305)          zkm1m, zkm1t
(  306)       REAL(wp), DIMENSION(jpi,jpk) ::   & !!! Boundary layer diffusivities/viscosities
(  307)          zblcm, zblct                          
(  308) #if defined key_zdfddm
(  309)       REAL(wp) ::                       & !!! double diffusion mixing
(  310)          zrrau, zds,                    &
(  311)          zavdds, zavddt,zinr 
(  312)       REAL(wp), DIMENSION(jpi,4) ::     &  
(  313)         zdifs
(  314)       REAL(wp), DIMENSION(jpi)     ::   &
(  315)          za2s, za3s, zkmps
(  316)       REAL(wp) ::                       & 
(  317)          zkm1s
(  318)       REAL(wp), DIMENSION(jpi,jpk) ::   & 
(  319)          zblcs                     
(  320) #endif
(  321)       !!--------------------------------------------------------------------
(  322) 
(  323) 
(  324)       ! Initialization (first time-step only)
(  325)       ! --------------
(  326)       IF( kt == nit000  )   CALL zdf_kpp_init
(  327)      
(  328)       zviscos(:,:,:) = 0.
(  329)       zblcm  (:,:  ) = 0. 
(  330)       zdiffut(:,:,:) = 0.
(  331)       zblct  (:,:  ) = 0. 
(  332) #if defined key_zdfddm
(  333)       zdiffus(:,:,:) = 0.
(  334)       zblcs  (:,:  ) = 0. 
(  335) #endif
(  336)       ghats(:,:,:) = 0.






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 7

(  337)      
(  338)       zBo   (:,:) = 0.
(  339)       zBosol(:,:) = 0.
(  340)       zustar(:,:) = 0.
(  341) 
(  342) 
(  343)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  344)       ! I. Interior diffusivity and viscosity at w points ( T interfaces)
(  345)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
(  346)       DO jk = 2, jpkm1
(  347)          DO jj = 2, jpjm1
(  348)             DO ji = fs_2, fs_jpim1 
(  349)                ! Mixing due to internal waves breaking
(  350)                ! -------------------------------------
(  351)                avmu(ji,jj,jk)  = difmiw 
(  352)                avt (ji,jj,jk)  = difsiw             
(  353)                ! Mixing due to vertical shear instability
(  354)                ! -------------------------------------               
(  355)                IF( ln_kpprimix ) THEN          
(  356)                   ! Compute the gradient Richardson  number at interfaces (zri):
(  357)                   ! LMD94, eq. 27 (is vertical smoothing needed : Rig=N^2 / (dz(u))^2 + (dz(v))^2
(  358)                   zdku2 =   ( un(ji - 1,jj,jk - 1) - un(ji - 1,jj,jk) ) &
(  359)                      &    * ( un(ji - 1,jj,jk - 1) - un(ji - 1,jj,jk) ) &
(  360)                      &    + ( un(ji    ,jj,jk - 1) - un(ji    ,jj,jk) ) &
(  361)                      &    * ( un(ji    ,jj,jk - 1) - un(ji    ,jj,jk) )  
(  362)                   
(  363)                   zdkv2 =   ( vn(ji,jj - 1,jk - 1) - vn(ji,jj - 1,jk) ) &
(  364)                      &    * ( vn(ji,jj - 1,jk - 1) - vn(ji,jj - 1,jk) ) &
(  365)                      &    + ( vn(ji,    jj,jk - 1) - vn(ji,    jj,jk) ) &
(  366)                      &    * ( vn(ji,    jj,jk - 1) - vn(ji,    jj,jk) )  
(  367) 
(  368)                   ze3sqr = 1. / ( fse3w(ji,jj,jk) * fse3w(ji,jj,jk) )
(  369)                   ! Square of vertical shear  at interfaces
(  370)                   zsh2   = 0.5 * 0.5 * ( zdku2 + zdkv2 ) * ze3sqr
(  371)                   zri    = MAX( rn2(ji,jj,jk), 0. ) / ( zsh2 + epsln ) 
(  372) #if defined key_cfg_1d
(  373)                   ! save the gradient richardson number
(  374)                   rig(ji,jj,jk) = zri * tmask(ji,jj,jk)
(  375) #endif                  
(  376)                   ! Evaluate f of Ri (zri) for shear instability store in zfri
(  377)                   ! LMD94, eq. 28a,b,c, figure 3 ; Rem: p1 is 3, hard coded
(  378)                   zfri  = MAX( zri , 0. )
(  379)                   zfri  = MIN( zfri / Riinfty , 1.0 )
(  380)                   zfri  = ( 1.0 - zfri * zfri )
(  381)                   zfri  = zfri * zfri  * zfri
(  382)                   ! add shear contribution to mixing coef. 
(  383)                   avmu(ji,jj,jk) =  avmu(ji,jj,jk) + difri * zfri   
(  384)                   avt (ji,jj,jk) =  avt (ji,jj,jk) + difri * zfri    
(  385)                ENDIF
(  386) #if defined key_zdfddm 
(  387)                avs (ji,jj,jk) =  avt (ji,jj,jk)              
(  388)                !  Double diffusion mixing ; NOT IN ROUTINE ZDFDDM.F90
(  389)                ! ------------------------------------------------------------------
(  390)                ! only retains positive value of rrau
(  391)                zrrau = MAX( rrau(ji,jj,jk), epsln )
(  392)                zds   = sn(ji,jj,jk-1) - sn(ji,jj,jk)
(  393)                IF( zrrau > 1. .AND. zds > 0.) THEN
(  394)                   !






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 8

(  395)                   ! Salt fingering case.
(  396)                   !---------------------
(  397)                   ! Compute interior diffusivity for double diffusive mixing of
(  398)                   ! salinity. Upper bound "zrrau" by "Rrho0"; (Rrho0=1.9, difcoefnuf=0.001).
(  399)                   ! After that set interior diffusivity for double diffusive mixing
(  400)                   ! of temperature
(  401)                   zavdds = MIN( zrrau, Rrho0 )
(  402)                   zavdds = ( zavdds - 1.0 ) / ( Rrho0 - 1.0 )
(  403)                   zavdds = 1.0 - zavdds * zavdds 
(  404)                   zavdds = zavdds * zavdds * zavdds 
(  405)                   zavdds = difssf * zavdds 
(  406)                   zavddt = 0.7 * zavdds
(  407)                ELSEIF( zrrau < 1. .AND. zrrau > 0. .AND. zds < 0.) THEN
(  408)                   !
(  409)                   ! Diffusive convection case.
(  410)                   !---------------------------
(  411)                   ! Compute interior diffusivity for double diffusive mixing of
(  412)                   ! temperature (Marmorino and Caldwell, 1976); 
(  413)                   ! Compute interior diffusivity for double diffusive mixing of salinity 
(  414)                   zinr   = 1. / zrrau
(  415)                   zavddt = 0.909 * EXP( 4.6 * EXP( -0.54* ( zinr - 1. ) ) ) 
(  416)                   zavddt = difsdc * zavddt
(  417)                   IF( zrrau < 0.5) THEN
(  418)                      zavdds = zavddt * 0.15 * zrrau
(  419)                   ELSE
(  420)                      zavdds = zavddt * (1.85 * zrrau - 0.85 ) 
(  421)                   ENDIF
(  422)                ELSE
(  423)                   zavddt = 0.
(  424)                   zavdds = 0.
(  425)                ENDIF 
(  426)                ! Add double diffusion contribution to temperature and salinity  mixing coefficients.
(  427)                avt (ji,jj,jk) =  avt (ji,jj,jk) +  zavddt 
(  428)                avs (ji,jj,jk) =  avs (ji,jj,jk) +  zavdds         
(  429) #endif                      
(  430)             END DO
(  431)          END DO
(  432)       END DO
(  433) 
(  434) 
(  435)       ! Radiative (zBosol) and non radiative (zBo) surface buoyancy
(  436)       !JMM at the time zdfkpp is called, q still holds the sum q + qsr
(  437)       !---------------------------------------------------------------------
(  438)       DO jj = 2, jpjm1
(  439)          DO ji = fs_2, fs_jpim1     
(  440)             IF( neos < 1) THEN   
(  441)                zt     = tn(ji,jj,1)
(  442)                zs     = sn(ji,jj,1) - 35.0
(  443)                zh     = fsdept(ji,jj,1)
(  444)                !  potential volumic mass
(  445)                zrhos  = rhop(ji,jj,1)
(  446)                zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
(  447)                   &                               - 0.203814e-03 ) * zt   &
(  448)                   &                               + 0.170907e-01 ) * zt   &
(  449)                   &   + 0.665157e-01                                      &
(  450)                   &   +     ( - 0.678662e-05 * zs                         &
(  451)                   &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
(  452)                   &   +   ( ( - 0.302285e-13 * zh                         &






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 9

(  453)                   &           - 0.251520e-11 * zs                         &
(  454)                   &           + 0.512857e-12 * zt * zt           ) * zh   &
(  455)                   &           - 0.164759e-06 * zs                         &
(  456)                   &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
(  457)                   &                               + 0.380374e-04 ) * zh
(  458) 
(  459)                zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
(  460)                   &                            - 0.301985e-05 ) * zt      &
(  461)                   &   + 0.785567e-03                                      &
(  462)                   &   + (     0.515032e-08 * zs                           &
(  463)                   &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
(  464)                   &   +(  (   0.121551e-17 * zh                           &
(  465)                   &         - 0.602281e-15 * zs                           &
(  466)                   &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
(  467)                   &                             + 0.408195e-10   * zs     &
(  468)                   &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
(  469)                   &                             - 0.121555e-07 ) * zh
(  470) 
(  471)                zthermal = zbeta * zalbet / ( rcp * zrhos + epsln )
(  472)                zhalin   = zbeta * sn(ji,jj,1) * rcs
(  473)             ELSE
(  474)                zrhos    = rhop(ji,jj,1) + rau0 * ( 1. - tmask(ji,jj,1) )
(  475)                zthermal = ralpha / ( rcp * zrhos + epsln )
(  476)                zhalin   = rbeta * sn(ji,jj,1) * rcs
(  477)             ENDIF
(  478)             ! Radiative surface buoyancy force
(  479)             zBosol(ji,jj) = grav * zthermal * qsr(ji,jj)
(  480)             ! Non radiative surface buoyancy force
(  481)             zBo   (ji,jj) = grav * zthermal * ( qt(ji,jj) - qsr(ji,jj) ) -  grav * zhalin * emp(ji,jj)
(  482)             ! Surface Temperature flux for non-local term
(  483)             wt0(ji,jj) = - qt(ji,jj) * ro0cpr * tmask(ji,jj,1)
(  484)             ! Surface salinity flux for non-local term
(  485)             ws0(ji,jj) = - ( emp(ji,jj) * sn(ji,jj,1) * rcs ) * tmask(ji,jj,1)
(  486)          ENDDO
(  487)       ENDDO
(  488) 
(  489)       zflageos = 0.5 + SIGN( 0.5, neos - 1. ) 
(  490)       !  Compute surface buoyancy forcing, Monin Obukhov and Ekman depths  
(  491)       !------------------------------------------------------------------    
(  492)       DO jj = 2, jpjm1
(  493)          DO ji = fs_2, fs_jpim1
(  494)             !  Reference surface density = density at first T point level   
(  495)             zrhos         = rhop(ji,jj,1) + zflageos * rau0 * ( 1. - tmask(ji,jj,1) )  
(  496)             ! Friction velocity (zustar), at T-point : LMD94 eq. 2
(  497)             ztx           = 0.5 * ( taux(ji,jj) + taux(ji - 1, jj    ) )
(  498)             zty           = 0.5 * ( tauy(ji,jj) + tauy(ji    , jj - 1) )
(  499)             ztau          = SQRT( ztx * ztx + zty * zty )
(  500)             zustar(ji,jj) = SQRT( ztau / ( zrhos +  epsln ) )
(  501)          ENDDO
(  502)       ENDDO
(  503) 
(  504) !CDIR NOVERRCHK  
(  505)       !                                               ! ===============
(  506)       DO jj = 2, jpjm1                                 !  Vertical slab
(  507)          !                                             ! ===============
(  508)          
(  509)          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  510)          ! II Compute Boundary layer mixing coef. and diagnose the new boundary layer depth






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 10

(  511)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  512)          
(  513)          ! Initialization
(  514)          jkmax       = 0
(  515)          zdept (:,:) = 0.
(  516)          zdepw (:,:) = 0.
(  517)          zriblk(:,:) = 0.
(  518)          zmoek (:,:) = 0.
(  519)          zvisc (:,:) = 0.
(  520)          zdift (:,:) = 0.
(  521) #if defined key_zdfddm
(  522)          zdifs (:,:) = 0.
(  523) #endif
(  524)          zmask (:,:) = 0.
(  525)          DO ji = fs_2, fs_jpim1
(  526)             zria(ji ) = 0.
(  527)             ! Maximum boundary layer depth
(  528)             ikbot     = mbathy(ji,jj) - 1 ! ikbot is the last T point in the water
(  529)             zhmax(ji) = fsdept(ji,jj,ikbot) - 0.001      
(  530)             ! Compute Monin obukhov length scale at the surface and Ekman depth:
(  531)             zbuofdep   = zBo(ji,jj) + zBosol(ji,jj) * ratt(1)
(  532)             zekman(ji) = rcekman * zustar(ji,jj) / ( ABS( ff(ji,jj) ) + epsln )
(  533)             zucube     = zustar(ji,jj) * zustar(ji,jj) * zustar(ji,jj) 
(  534)             zmoa(ji)   = zucube / ( vonk * ( zbuofdep + epsln ) )   
(  535) #if defined key_cfg_1d
(  536)             ! store the surface buoyancy forcing
(  537)             zstabl        = 0.5 + SIGN( 0.5, zbuofdep )
(  538)             buof(ji,jj,1) = zbuofdep * tmask(ji,jj,1)
(  539)             ! store the moning-oboukov length scale at surface
(  540)             zmob          = zstabl * zmoa(ji) + ( 1.0 - zstabl ) * fsdept(ji,jj,1)
(  541)             mols(ji,jj,1) = MIN( zmob , zhmax(ji) ) * tmask(ji,jj,1)
(  542)             ! store Ekman depth
(  543)             zek           = zstabl * zekman(ji) + ( 1.0 - zstabl ) * fsdept(ji,jj,1)  
(  544)             ekdp(ji,jj )  = MIN( zek , zhmax(ji) ) * tmask(ji,jj,1)  
(  545) #endif 
(  546)          END DO     
(  547)          ! Compute the pipe
(  548)          ! ---------------------
(  549)          DO jk = 2, jpkm1
(  550)             DO ji = fs_2, fs_jpim1
(  551)                ! Compute bfsfc = Bo + radiative contribution down to hbf*depht
(  552)                zbuofdep = zBo(ji,jj) + zBosol(ji,jj) * ratt(jk)
(  553)                ! Flag (zstabl  = 1) if positive forcing
(  554)                zstabl   =  0.5 + SIGN(  0.5, zbuofdep)
(  555) 
(  556)                !   Compute bulk richardson number zrib at depht 
(  557)                !-------------------------------------------------------
(  558)                !                           [Br - B(d)] * d         zrinum
(  559)                !             Rib(z) = ----------------------- = -------------
(  560)                !                       |Vr - V(d)|^2 + Vt(d)^2   zdVsq + zVtsq
(  561)                !
(  562)                ! First compute zt,zs,zu,zv = means in the surface layer < epsilon*depht  
(  563)                ! Else surface values are taken at the first T level.
(  564)                ! For stability, resolved vertical shear is computed with "before velocities".
(  565)                zref = epsilon * fsdept(ji,jj,jk)
(  566) #if defined key_kppcustom
(  567)                ! zref = gdept(1)
(  568)                zref = fsdept(ji,jj,1)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 11

(  569)                zt   = tn(ji,jj,1)
(  570)                zs   = sn(ji,jj,1)
(  571)                zrh  = rhop(ji,jj,1)
(  572)                zu   = ( ub(ji,jj,1) + ub(ji - 1,jj    ,1) ) / MAX( 1. , umask(ji,jj,1) + umask(ji - 1,jj   ,1) )
(  573)                zv   = ( vb(ji,jj,1) + vb(ji    ,jj - 1,1) ) / MAX( 1. , vmask(ji,jj,1) + vmask(ji   ,jj - 1,1) )
(  574) #else
(  575)                zt   = 0.
(  576)                zs   = 0.
(  577)                zu   = 0.
(  578)                zv   = 0.
(  579)                zrh  = 0.
(  580)                ! vertically integration over the upper epsilon*gdept(jk) ; del () array is computed once in zdf_kpp_init
(  581)                DO jm = 1, jpkm1
(  582)                   zt   = zt  + del(jk,jm) * tn(ji,jj,jm)
(  583)                   zs   = zs  + del(jk,jm) * sn(ji,jj,jm)
(  584)                   zu   = zu  + 0.5 * del(jk,jm) &
(  585)                      &            * ( ub(ji,jj,jm) + ub(ji - 1,jj,jm) ) &
(  586)                      &            / MAX( 1. , umask(ji,jj,jm) + umask(ji - 1,jj,jm) )
(  587)                   zv   = zv  + 0.5 * del(jk,jm) &
(  588)                      &            * ( vb(ji,jj,jm) + vb(ji,jj - 1,jm) ) &
(  589)                      &            / MAX( 1. , vmask(ji,jj,jm) + vmask(ji,jj - 1,jm) )
(  590)                   zrh  = zrh + del(jk,jm) * rhop(ji,jj,jm)
(  591)                END DO
(  592) #endif
(  593)                zsr = SQRT( ABS( sn(ji,jj,jk) ) )
(  594)                ! depth
(  595)                zh = fsdept(ji,jj,jk)
(  596)                ! compute compression terms on density
(  597)                ze  = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  598)                zbw = (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  599)                zb  = zbw + ze * zs
(  600)                
(  601)                zd  = -2.042967e-2
(  602)                zc  =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  603)                zaw = ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
(  604)                za  = ( zd*zsr + zc ) *zs + zaw
(  605)                
(  606)                zb1 =   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  607)                za1 = ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  608)                zkw = ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
(  609)                zk0 = ( zb1*zsr + za1 )*zs + zkw
(  610)                zcomp =   1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )
(  611)                
(  612) #if defined key_kppcustom
(  613)                ! potential density of water(zrh = zt,zs at level jk):
(  614)                zrhdr = zrh / zcomp
(  615) #else
(  616)                ! potential density of water(ztref,zsref at level jk):
(  617)                ! compute volumic mass pure water at atm pressure
(  618)                IF ( neos < 1 ) THEN
(  619)                   zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  620)                      &       -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  621)                   ! seawater volumic mass atm pressure
(  622)                   zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
(  623)                      &   -4.0899e-3 ) *zt+0.824493
(  624)                   zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  625)                   zr4= 4.8314e-4              
(  626)                   ! potential volumic mass (reference to the surface)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 12

(  627)                   zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1                 
(  628)                   zrhdr = zrhop / zcomp
(  629)                ELSE
(  630)                   zrhdr = zrh / zcomp
(  631)                ENDIF
(  632) #endif
(  633)                
(  634)                ! potential density of ambiant water at level jk :
(  635)                zrhd   = ( rhd(ji,jj,jk) * rau0 + rau0 )  
(  636)                
(  637)                ! And now the Rib number numerator .
(  638)                zrinum = grav * ( zrhd - zrhdr ) / rau0
(  639)                zrinum = zrinum * ( fsdept(ji,jj,jk) - zref ) * tmask(ji,jj,jk)
(  640)            
(  641)                ! Resolved shear contribution to Rib at depth T-point (zdVsq)
(  642)                ztx    =   ( ub( ji , jj ,jk)   +  ub(ji - 1, jj ,jk) ) &
(  643)                   &     / MAX( 1. , umask( ji , jj ,jk) + umask(ji - 1, jj ,jk) )   
(  644)                zty    =   ( vb( ji , jj ,jk)   +  vb(ji  ,jj - 1,jk) ) &
(  645)                   &     / MAX( 1., vmask( ji , jj ,jk) + vmask(ji  ,jj - 1,jk) ) 
(  646)                
(  647)                zdVsq  = ( zu - ztx ) * ( zu - ztx ) + ( zv - zty ) * ( zv - zty )
(  648)                
(  649)                ! Scalar turbulent velocity scale zws for hbl=gdept
(  650)                zscale = zstabl + ( 1.0 - zstabl ) * epsilon
(  651)                zehat  = vonk * zscale * fsdept(ji,jj,jk) * zbuofdep
(  652)                zucube = zustar(ji,jj) * zustar(ji,jj) * zustar(ji,jj)              
(  653)                zeta   = zehat / ( zucube + epsln )
(  654)                
(  655)                IF( zehat > 0. ) THEN
(  656)                   ! Stable case
(  657)                   zws  = vonk * zustar(ji,jj) / ( 1.0 + rconc1 * zeta )
(  658)                ELSE
(  659)                   ! Unstable case
(  660) #if defined key_kpplktb
(  661)                   ! use lookup table
(  662)                   zd     = zehat - dehatmin
(  663)                   il     = INT( zd / dezehat )
(  664)                   il     = MIN( il, nilktbm1 )
(  665)                   il     = MAX( il, 1 )
(  666)                   
(  667)                   ud     = zustar(ji,jj) - ustmin
(  668)                   jl     = INT( ud / deustar )
(  669)                   jl     = MIN( jl, njlktbm1 )
(  670)                   jl     = MAX( jl, 1 )
(  671)                   
(  672)                   zfrac  = zd / dezehat - FLOAT( il )  
(  673)                   ufrac  = ud / deustar - FLOAT( jl )
(  674)                   zwas   = ( 1. - zfrac ) * wslktb(il,jl+1) + zfrac * wslktb(il+1,jl+1)
(  675)                   zwbs   = ( 1. - zfrac ) * wslktb(il,jl  ) + zfrac * wslktb(il+1,jl  )
(  676)                   !
(  677)                   zws    = ( 1. - ufrac ) * zwbs + ufrac * zwas
(  678) #else
(  679)                   ! use analytical functions:
(  680)                   zcons  = 0.5 + SIGN( 0.5 , ( rzetas - zeta ) )
(  681)                   zwcons = vonk * zustar(ji,jj) * ( ( ABS( rconas - rconcs * zeta ) )**pthird ) 
(  682)                   zwsun  = vonk * zustar(ji,jj) * SQRT( ABS ( 1.0 - rconc2 * zeta ) )
(  683)                   !
(  684)                   zws    = zcons * zwcons +  ( 1.0 - zcons) * zwsun






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 13

(  685) #endif
(  686)                ENDIF
(  687)                
(  688)                ! Turbulent shear contribution to Rib (zVtsq) bv frequency at levels  ( ie T-point jk)
(  689)                zrn2   = 0.5 * ( rn2(ji,jj,jk) + rn2(ji,jj,jk+1) )   
(  690)                zbvzed = SQRT( ABS( zrn2 ) ) 
(  691)                zVtsq  = fsdept(ji,jj,jk) * zws * zbvzed  * Vtc
(  692)                
(  693)                ! Finally, the bulk Richardson number at depth fsdept(i,j,k) 
(  694)                zrib  = zrinum   / ( zdVsq + zVtsq + epsln )
(  695)  
(  696)                ! Find subscripts around the boundary layer depth, build the pipe
(  697)                ! ----------------------------------------------------------------
(  698) 
(  699)                ! Flag (zflagri = 1) if zrib < Ricr  
(  700)                zflagri = 0.5 + SIGN( 0.5, ( Ricr - zrib ) )
(  701)                !  Flag (zflagh  = 1) if still within overall boundary layer
(  702)                zflagh  = 0.5 + SIGN( 0.5, ( fsdept(ji,jj,1) - zdept(ji,2) ) )
(  703)                
(  704)                ! Ekman layer depth
(  705)                zek     = zstabl * zekman(ji) + ( 1.0 - zstabl ) * zhmax(ji)
(  706)                zflag   = 0.5 + SIGN( 0.5, ( zek - fsdept(ji,jj,jk-1) ) )
(  707)                zek     = zflag * zek + ( 1.0 - zflag ) * zhmax(ji)
(  708)                zflagek = 0.5 + SIGN( 0.5, ( zek - fsdept(ji,jj,jk) ) )
(  709)                ! Flag (zflagmo = 1) if still within stable Monin-Obukhov and in water
(  710)                zmob    = zucube / ( vonk * ( zbuofdep + epsln ) )  
(  711)                ztemp   = zstabl * zmob + ( 1.0 - zstabl) * zhmax(ji)
(  712)                ztemp   = MIN( ztemp , zhmax(ji) ) 
(  713)                zflagmo = 0.5 + SIGN( 0.5, ( ztemp - fsdept(ji,jj,jk) ) )             
(  714) 
(  715)                ! No limitation by Monin Obukhov or Ekman depths:
(  716) !               zflagek = 1.0
(  717) !               zflagmo = 0.5 + SIGN( 0.5, ( zhmax(ji) - fsdept(ji,jj,jk) ) )
(  718) 
(  719)                ! Load  pipe via zflagkb  for later calculations
(  720)                ! Flag (zflagkb = 1) if zflagh = 1 and (zflagri = 0 or zflagek = 0 or zflagmo = 0)
(  721)                zflagkb = zflagh * ( 1.0 - ( zflagri * zflagek * zflagmo ) )
(  722) 
(  723)                zmask(ji,jk) = zflagh
(  724)                jkp2         = MIN( jk+2 , ikbot )
(  725)                jkm1         = MAX( jk-1 , 2 )
(  726)                jkmax        = MAX( jkmax, jk * INT( REAL( zflagh+epsln ) ) ) 
(  727) 
(  728)                zdept(ji,1)  = zdept(ji,1) + zflagkb * fsdept(ji,jj,jk-1) 
(  729)                zdept(ji,2)  = zdept(ji,2) + zflagkb * fsdept(ji,jj,jk  ) 
(  730)                zdept(ji,3)  = zdept(ji,3) + zflagkb * fsdept(ji,jj,jk+1) 
(  731) 
(  732)                zdepw(ji,1)  = zdepw(ji,1) + zflagkb * fsdepw(ji,jj,jk-1) 
(  733)                zdepw(ji,2)  = zdepw(ji,2) + zflagkb * fsdepw(ji,jj,jk  ) 
(  734)                zdepw(ji,3)  = zdepw(ji,3) + zflagkb * fsdepw(ji,jj,jk+1)
(  735)                zdepw(ji,4)  = zdepw(ji,4) + zflagkb * fsdepw(ji,jj,jkp2)  
(  736) 
(  737)                zriblk(ji,1) = zriblk(ji,1) + zflagkb * zria(ji)
(  738)                zriblk(ji,2) = zriblk(ji,2) + zflagkb * zrib
(  739) 
(  740)                zmoek (ji,0) = zmoek (ji,0) + zflagkb * zek
(  741)                zmoek (ji,1) = zmoek (ji,1) + zflagkb * zmoa(ji)
(  742)                zmoek (ji,2) = zmoek (ji,2) + zflagkb * ztemp  






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 14

(  743)                ! Save Monin Obukhov depth
(  744)                zmoa  (ji)   = zmob
(  745)            
(  746)                zvisc(ji,1) = zvisc(ji,1) + zflagkb * avmu(ji,jj,jkm1)
(  747)                zvisc(ji,2) = zvisc(ji,2) + zflagkb * avmu(ji,jj,jk  )
(  748)                zvisc(ji,3) = zvisc(ji,3) + zflagkb * avmu(ji,jj,jk+1)
(  749)                zvisc(ji,4) = zvisc(ji,4) + zflagkb * avmu(ji,jj,jkp2)
(  750)                
(  751)                zdift(ji,1) = zdift(ji,1) + zflagkb * avt (ji,jj,jkm1)
(  752)                zdift(ji,2) = zdift(ji,2) + zflagkb * avt (ji,jj,jk  )
(  753)                zdift(ji,3) = zdift(ji,3) + zflagkb * avt (ji,jj,jk+1)
(  754)                zdift(ji,4) = zdift(ji,4) + zflagkb * avt (ji,jj,jkp2)
(  755) 
(  756) #if defined key_zdfddm
(  757)                zdifs(ji,1) = zdifs(ji,1) + zflagkb * avs (ji,jj,jkm1)
(  758)                zdifs(ji,2) = zdifs(ji,2) + zflagkb * avs (ji,jj,jk  )
(  759)                zdifs(ji,3) = zdifs(ji,3) + zflagkb * avs (ji,jj,jk+1)
(  760)                zdifs(ji,4) = zdifs(ji,4) + zflagkb * avs (ji,jj,jkp2)
(  761) #endif               
(  762)                ! Save the Richardson number 
(  763)                zria  (ji)   = zrib  
(  764) #if defined key_cfg_1d
(  765)                ! store buoyancy length scale
(  766)                buof(ji,jj,jk) = zbuofdep * tmask(ji,jj,jk) 
(  767)                ! store Monin Obukhov
(  768)                zmob           = zstabl * zmob + ( 1.0 - zstabl) * fsdept(ji,jj,1)
(  769)                mols(ji,jj,jk) = MIN( zmob , zhmax(ji) ) * tmask(ji,jj,jk) 
(  770)                ! Bulk Richardson number
(  771)                rib(ji,jj,jk)  = zrib * tmask(ji,jj,jk)             
(  772) #endif               
(  773)             END DO
(  774)          END DO
(  775)          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  776)          ! III PROCESS THE PIPE
(  777)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  778)          
(  779)          DO ji = fs_2, fs_jpim1 
(  780)             
(  781)             ! Find the boundary layer depth zhbl
(  782)             ! ----------------------------------------
(  783)             
(  784)             ! Interpolate monin Obukhov and critical Ri mumber depths   
(  785)             ztemp = zdept(ji,2) - zdept(ji,1)
(  786)             zflag = ( Ricr - zriblk(ji,1) ) / ( zriblk(ji,2) - zriblk(ji,1)  + epsln )
(  787)             zhrib = zdept(ji,1) + zflag * ztemp      
(  788) 
(  789)             IF( zriblk(ji,2) < Ricr ) zhrib = zhmax(ji) 
(  790)          
(  791)             IF( zmoek(ji,2) < zdept(ji,2) ) THEN
(  792)                IF ( zmoek(ji,1) < 0. ) THEN
(  793)                   zmob = zdept(ji,2) - epsln
(  794)                ELSE
(  795)                   zmob = ztemp + zmoek(ji,1) - zmoek(ji,2)
(  796)                   zmob = ( zmoek(ji,1) * zdept(ji,2) - zmoek(ji,2) * zdept(ji,1) ) / zmob
(  797)                   zmob = MAX( zmob , zdept(ji,1) + epsln )                
(  798)                ENDIF
(  799)             ELSE           
(  800)                zmob = zhmax(ji) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 15

(  801)             ENDIF
(  802)             ztemp   = MIN( zmob , zmoek(ji,0) )
(  803)                          
(  804)             ! Finally, the boundary layer depth, zhbl 
(  805)             zhbl(ji) = MAX( fsdept(ji,jj,1) + epsln, MIN( zhrib , ztemp ) )
(  806)             
(  807)             ! Save hkpp for further diagnostics (optional)
(  808)             hkpp(ji,jj) = zhbl(ji) * tmask(ji,jj,1) 
(  809)           
(  810)             ! Correct mask if zhbl < fsdepw(ji,jj,2) for no viscosity/diffusivity enhancement at fsdepw(ji,jj,2)
(  811)             !     zflag = 1 if zhbl(ji) > fsdepw(ji,jj,2)
(  812)             IF( zhbl(ji) < fsdepw(ji,jj,2) ) zmask(ji,2) = 0.
(  813)           
(  814)             
(  815)             !  Velocity scales at depth zhbl
(  816)             ! -----------------------------------
(  817)             
(  818)             !  Compute bouyancy forcing down to zhbl
(  819)             ztemp    = -hbf * zhbl(ji)
(  820)             zatt1    = 1.0 - ( rabs * EXP( ztemp / xsi1 ) + ( 1.0 - rabs ) * EXP( ztemp / xsi2 ) )
(  821)             zbuofdep = zBo(ji,jj) + zBosol(ji,jj) * zatt1
(  822)             zstabl   = 0.5 + SIGN( 0.5 , zbuofdep ) 
(  823) 
(  824)             zbuofdep = zbuofdep + zstabl * epsln
(  825) 
(  826)             zscale = zstabl + ( 1.0 - zstabl ) * epsilon          
(  827)             zehat  = vonk * zscale * zhbl(ji) * zbuofdep
(  828)             zucube = zustar(ji,jj) * zustar(ji,jj) * zustar(ji,jj)              
(  829)             zeta   = zehat / ( zucube + epsln )
(  830)             
(  831)             IF( zehat > 0. ) THEN
(  832)                ! Stable case
(  833)                zws  = vonk * zustar(ji,jj) / ( 1.0 + rconc1 * zeta )
(  834)                zwm  = zws
(  835)             ELSE
(  836)                ! Unstable case
(  837) #if defined key_kpplktb
(  838)                ! use lookup table
(  839)                zd     = zehat - dehatmin
(  840)                il     = INT( zd / dezehat )
(  841)                il     = MIN( il, nilktbm1 )
(  842)                il     = MAX( il, 1 )
(  843)                
(  844)                ud     = zustar(ji,jj) - ustmin
(  845)                jl     = INT( ud / deustar )
(  846)                jl     = MIN( jl, njlktbm1 )
(  847)                jl     = MAX( jl, 1 )
(  848)                
(  849)                zfrac  = zd / dezehat - FLOAT( il )  
(  850)                ufrac  = ud / deustar - FLOAT( jl )
(  851)                zwas   = ( 1. - zfrac ) * wslktb(il,jl+1) + zfrac * wslktb(il+1,jl+1)
(  852)                zwbs   = ( 1. - zfrac ) * wslktb(il,jl  ) + zfrac * wslktb(il+1,jl  )
(  853)                zwam   = ( 1. - zfrac ) * wmlktb(il,jl+1) + zfrac * wmlktb(il+1,jl+1)
(  854)                zwbm   = ( 1. - zfrac ) * wmlktb(il,jl  ) + zfrac * wmlktb(il+1,jl  )
(  855)                !
(  856)                zws    = ( 1. - ufrac ) * zwbs + ufrac * zwas
(  857)                zwm    = ( 1. - ufrac ) * zwbm + ufrac * zwam
(  858) #else






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 16

(  859)                ! use analytical functions
(  860)                zconm  = 0.5 + SIGN( 0.5, ( rzetam - zeta) )
(  861)                zcons  = 0.5 + SIGN( 0.5, ( rzetas - zeta) )
(  862)                
(  863)                ! Momentum : zeta < rzetam (zconm = 1)
(  864)                ! Scalars  : zeta < rzetas (zcons = 1) 
(  865)                zwconm = zustar(ji,jj) * vonk * ( ( ABS( rconam - rconcm * zeta) )**pthird )
(  866)                zwcons = zustar(ji,jj) * vonk * ( ( ABS( rconas - rconcs * zeta) )**pthird )
(  867)                
(  868)                ! Momentum : rzetam <= zeta < 0 (zconm = 0)
(  869)                ! Scalars  : rzetas <= zeta < 0 (zcons = 0)	
(  870)                zwmun  = SQRT( ABS( 1.0 - rconc2 * zeta ) )
(  871)                zwsun  = vonk * zustar(ji,jj) * zwmun
(  872)                zwmun  = vonk * zustar(ji,jj) * SQRT(zwmun)
(  873)                !
(  874)                zwm    = zconm * zwconm + ( 1.0 - zconm ) * zwmun
(  875)                zws    = zcons * zwcons + ( 1.0 - zcons ) * zwsun
(  876)                
(  877) #endif
(  878)             ENDIF
(  879)             
(  880)             
(  881)             ! Viscosity, diffusivity values and derivatives at h
(  882)             ! --------------------------------------------------------
(  883)             
(  884)             ! check between at which interfaces is located zhbl(ji)
(  885)             ! ztemp = 1, zdepw(ji,2) < zhbl <  zdepw(ji,3)
(  886)             ! ztemp = 0, zdepw(ji,1) < zhbl <  zdepw(ji,2)
(  887)             ztemp  =  0.5 + SIGN( 0.5, ( zhbl(ji) - zdepw(ji,2) ) )  
(  888)             zdep21 =   zdepw(ji,2) - zdepw(ji,1) + epsln
(  889)             zdep32 =   zdepw(ji,3) - zdepw(ji,2) + epsln
(  890)             zdep43 =   zdepw(ji,4) - zdepw(ji,3) + epsln  
(  891)             
(  892)             ! Compute R as in LMD94, eq D5b
(  893)             zdelta =  ( zhbl(ji) - zdepw(ji,2) ) *         ztemp   / zdep32   &
(  894)                &    + ( zhbl(ji) - zdepw(ji,1) ) * ( 1.0 - ztemp ) / zdep21 
(  895)             
(  896)             ! Compute the vertical derivative of viscosities (zdzh) at z=zhbl(ji)
(  897)             zdzup  =  ( zvisc(ji,2) - zvisc(ji,3) ) *         ztemp   / zdep32 &
(  898)                &    + ( zvisc(ji,1) - zvisc(ji,2) ) * ( 1.0 - ztemp ) / zdep21
(  899)             
(  900)             zdzdn  =  ( zvisc(ji,3) - zvisc(ji,4) ) *         ztemp   / zdep43 &
(  901)                &    + ( zvisc(ji,2) - zvisc(ji,3) ) * ( 1.0 - ztemp ) / zdep32
(  902)             
(  903)             ! LMD94, eq D5b :          
(  904)             zdzh   = ( 1.0 - zdelta ) * zdzup + zdelta * zdzdn
(  905)             zdzh   = MAX( zdzh , 0. )
(  906)             
(  907)             ! Compute viscosities (zvath) at z=zhbl(ji), LMD94 eq D5a
(  908)             zvath  =          ztemp   * ( zvisc(ji,3) + zdzh * ( zdepw(ji,3) - zhbl(ji) ) ) &
(  909)                &    + ( 1.0 - ztemp ) * ( zvisc(ji,2) + zdzh * ( zdepw(ji,2) - zhbl(ji) ) )
(  910)             
(  911)             ! Compute G (zgat1) and its derivative (zdat1) at z=hbl(ji), LMD94 eq 18
(  912)             
(  913)             ! Vertical derivative of velocity scale divided by velocity scale squared at z=hbl(ji) 
(  914)             ! (non zero only in stable conditions)
(  915)             zflag  =  -zstabl * rconc1 * zbuofdep / ( zucube * zustar(ji,jj) + epsln )
(  916)             






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 17

(  917)             ! G at its derivative at z=hbl:
(  918)             zgat1  = zvath  / ( zhbl(ji) * ( zwm + epsln )  )
(  919)             zdat1  = -zdzh  / ( zwm + epsln ) -  zflag * zvath / zhbl(ji)
(  920)             
(  921)             ! G coefficients, LMD94 eq 17
(  922)             za2m(ji) = -2.0 + 3.0 * zgat1 - zdat1
(  923)             za3m(ji) =  1.0 - 2.0 * zgat1 + zdat1
(  924) 
(  925)             
(  926)             ! Compute the vertical derivative of temperature diffusivities (zdzh) at z=zhbl(ji)
(  927)             zdzup  =  ( zdift(ji,2) - zdift(ji,3) ) *         ztemp   / zdep32 &
(  928)                &    + ( zdift(ji,1) - zdift(ji,2) ) * ( 1.0 - ztemp ) / zdep21
(  929)             
(  930)             zdzdn  =  ( zdift(ji,3) - zdift(ji,4) ) *         ztemp   / zdep43 &
(  931)                &    + ( zdift(ji,2) - zdift(ji,3) ) * ( 1.0 - ztemp ) / zdep32
(  932)             
(  933)             ! LMD94, eq D5b :          
(  934)             zdzh   = ( 1.0 - zdelta ) * zdzup + zdelta * zdzdn
(  935)             zdzh   = MAX( zdzh , 0. )
(  936)             
(  937)             
(  938)             ! Compute diffusivities (zvath) at z=zhbl(ji), LMD94 eq D5a
(  939)             zvath  =          ztemp   * ( zdift(ji,3) + zdzh * ( zdepw(ji,3) - zhbl(ji) ) ) &
(  940)                &    + ( 1.0 - ztemp ) * ( zdift(ji,2) + zdzh * ( zdepw(ji,2) - zhbl(ji) ) )
(  941)                         
(  942)             ! G at its derivative at z=hbl:
(  943)             zgat1  = zvath  / ( zhbl(ji) * ( zws + epsln )  )
(  944)             zdat1  = -zdzh  / ( zws + epsln ) -  zflag * zvath / zhbl(ji)
(  945)             
(  946)             ! G coefficients, LMD94 eq 17
(  947)             za2t(ji) = -2.0 + 3.0 * zgat1 - zdat1
(  948)             za3t(ji) =  1.0 - 2.0 * zgat1 + zdat1
(  949) 
(  950) #if defined key_zdfddm
(  951)             ! Compute the vertical derivative of salinities diffusivities (zdzh) at z=zhbl(ji)
(  952)             zdzup  =  ( zdifs(ji,2) - zdifs(ji,3) ) *         ztemp   / zdep32 &
(  953)                &    + ( zdifs(ji,1) - zdifs(ji,2) ) * ( 1.0 - ztemp ) / zdep21
(  954)             
(  955)             zdzdn  =  ( zdifs(ji,3) - zdifs(ji,4) ) *         ztemp   / zdep43 &
(  956)                &    + ( zdifs(ji,2) - zdifs(ji,3) ) * ( 1.0 - ztemp ) / zdep32
(  957)             
(  958)             ! LMD94, eq D5b :          
(  959)             zdzh   = ( 1.0 - zdelta ) * zdzup + zdelta * zdzdn
(  960)             zdzh   = MAX( zdzh , 0. )           
(  961)             
(  962)             ! Compute diffusivities (zvath) at z=zhbl(ji), LMD94 eq D5a
(  963)             zvath  =          ztemp   * ( zdifs(ji,3) + zdzh * ( zdepw(ji,3) - zhbl(ji) ) ) &
(  964)                &    + ( 1.0 - ztemp ) * ( zdifs(ji,2) + zdzh * ( zdepw(ji,2) - zhbl(ji) ) )
(  965)                         
(  966)             ! G at its derivative at z=hbl:
(  967)             zgat1  = zvath  / ( zhbl(ji) * ( zws + epsln )  )
(  968)             zdat1  = -zdzh  / ( zws + epsln ) -  zflag * zvath / zhbl(ji)
(  969)             
(  970)             ! G coefficients, LMD94 eq 17
(  971)             za2s(ji) = -2.0 + 3.0 * zgat1 - zdat1
(  972)             za3s(ji) =  1.0 - 2.0 * zgat1 + zdat1
(  973) #endif
(  974) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 18

(  975)             !-------------------turn off interior matching here------
(  976)             !          za2(ji,1) = -2.0
(  977)             !          za3(ji,1) =  1.0
(  978)             !          za2(ji,2) = -2.0
(  979)             !          za3(ji,2) =  1.0
(  980)             !--------------------------------------------------------
(  981)             
(  982)             !  Compute Enhanced Mixing Coefficients (LMD94,eq D6)
(  983)             ! ---------------------------------------------------------------
(  984)             
(  985)             ! Delta 
(  986)             zdelta  = ( zhbl(ji)  - zdept(ji,1) ) / ( zdept(ji,2) - zdept(ji,1) + epsln )
(  987)             zdelta2 = zdelta * zdelta
(  988)             
(  989)             !  Mixing coefficients at first level above h (zdept(ji,1)) 
(  990)             ! and at first interface in the pipe (zdepw(ji,2))
(  991)             
(  992)             ! At first T level above h (zdept(ji,1)) (always in the boundary layer)
(  993)             zsig    = zdept(ji,1) / zhbl(ji)
(  994)             ztemp   = zstabl * zsig  + ( 1.0 - zstabl ) * MIN( zsig , epsilon )
(  995)             zehat   = vonk * ztemp * zhbl(ji) * zbuofdep
(  996)             zeta    = zehat / ( zucube + epsln)
(  997)             zwst    = vonk * zustar(ji,jj) / ( ABS( 1.0 + rconc1 * zeta ) + epsln)
(  998)             zwm     = zstabl * zwst + ( 1.0 - zstabl ) * zwm
(  999)             zws     = zstabl * zwst + ( 1.0 - zstabl ) * zws
( 1000) 
( 1001)             zkm1m  = zhbl(ji) * zwm * zsig * ( 1.0 + zsig * ( za2m(ji) + zsig * za3m(ji) ) )
( 1002)             zkm1t  = zhbl(ji) * zws * zsig * ( 1.0 + zsig * ( za2t(ji) + zsig * za3t(ji) ) )
( 1003) #if defined key_zdfddm
( 1004)             zkm1s  = zhbl(ji) * zws * zsig * ( 1.0 + zsig * ( za2s(ji) + zsig * za3s(ji) ) )
( 1005) #endif                        
( 1006)             ! At first W level in the pipe (zdepw(ji,2)) (not always in the boundary layer ):
( 1007)             zsig    = MIN( zdepw(ji,2) / zhbl(ji) , 1.0 )
( 1008)             ztemp   = zstabl * zsig + ( 1.0 - zstabl ) * MIN( zsig , epsilon )
( 1009)             zehat   = vonk * ztemp * zhbl(ji) * zbuofdep
( 1010)             zeta    = zehat / ( zucube + epsln )
( 1011)             zwst    = vonk * zustar(ji,jj) / ( ABS( 1.0 + rconc1 * zeta ) + epsln)
( 1012)             zws     = zstabl * zws + ( 1.0 - zstabl ) * zws
( 1013)             zwm     = zstabl * zws + ( 1.0 - zstabl ) * zwm
( 1014) 
( 1015)             zkmpm(ji) = zhbl(ji) * zwm * zsig * ( 1.0 + zsig * ( za2m(ji) + zsig * za3m(ji) ) )
( 1016)             zkmpt(ji) = zhbl(ji) * zws * zsig * ( 1.0 + zsig * ( za2t(ji) + zsig * za3t(ji) ) )
( 1017) #if defined key_zdfddm
( 1018)             zkmps(ji) = zhbl(ji) * zws * zsig * ( 1.0 + zsig * ( za2s(ji) + zsig * za3s(ji) ) )
( 1019) #endif  
( 1020)       
( 1021)             ! check if this point is in the boundary layer,else take interior viscosity/diffusivity:
( 1022)             zflag       = 0.5 + SIGN( 0.5, ( zhbl(ji) - zdepw(ji,2) ) )
( 1023)             zkmpm(ji) = zkmpm(ji) * zflag + ( 1.0 - zflag ) * zvisc(ji,2)
( 1024)             zkmpt(ji) = zkmpt(ji) * zflag + ( 1.0 - zflag ) * zdift(ji,2)
( 1025) #if defined key_zdfddm
( 1026)             zkmps(ji) = zkmps(ji) * zflag + ( 1.0 - zflag ) * zdifs(ji,2)
( 1027) #endif
( 1028) 
( 1029)             ! Enhanced viscosity/diffusivity at zdepw(ji,2)
( 1030)             ztemp     = ( 1.0 - 2.0 * zdelta + zdelta2 ) * zkm1m + zdelta2 * zkmpm(ji)
( 1031)             zkmpm(ji) = ( 1.0 - zdelta ) * zvisc(ji,2) + zdelta * ztemp
( 1032)             ztemp     = ( 1.0 - 2.0 * zdelta + zdelta2 ) * zkm1t + zdelta2 * zkmpt(ji)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 19

( 1033)             zkmpt(ji) = ( 1.0 - zdelta ) * zdift(ji,2) + zdelta * ztemp
( 1034) #if defined key_zdfddm
( 1035)             ztemp     = ( 1.0 - 2.0 * zdelta + zdelta2 ) * zkm1s + zdelta2 * zkmps(ji)
( 1036)             zkmps(ji) = ( 1.0 - zdelta ) * zdifs(ji,2) + zdelta * ztemp
( 1037) #endif            
( 1038) 
( 1039)          END DO
( 1040)          !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
( 1041)          ! IV. Compute vertical eddy viscosity and diffusivity coefficients
( 1042)          !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
( 1043)          
( 1044)          DO jk  = 2, jkmax
( 1045)             
( 1046)             ! Compute turbulent velocity scales on the interfaces
( 1047)             ! --------------------------------------------------------
( 1048)             DO  ji = fs_2, fs_jpim1
( 1049)                zbuofdep = zBo(ji,jj) + zBosol(ji,jj) * zatt1
( 1050)                zstabl   = 0.5 + SIGN( 0.5 , zbuofdep ) 
( 1051)                zbuofdep = zbuofdep + zstabl * epsln          
( 1052)                zsig    = fsdepw(ji,jj,jk) / zhbl(ji)
( 1053)                ztemp   = zstabl * zsig + ( 1. - zstabl ) * MIN( zsig , epsilon )
( 1054)                zehat   = vonk * ztemp * zhbl(ji) * zbuofdep
( 1055)                zucube  = zustar(ji,jj) * zustar(ji,jj) * zustar(ji,jj)
( 1056)                zeta    = zehat / ( zucube + epsln )
( 1057) 
( 1058)                IF( zehat > 0. ) THEN
( 1059)                   ! Stable case
( 1060)                   zws  = vonk * zustar(ji,jj) / ( 1.0 + rconc1 * zeta )
( 1061)                   zwm  = zws
( 1062)                ELSE
( 1063)                   ! Unstable case
( 1064) #if defined key_kpplktb
( 1065)                   ! use lookup table
( 1066)                   zd     = zehat - dehatmin
( 1067)                   il     = INT( zd / dezehat )
( 1068)                   il     = MIN( il, nilktbm1 )
( 1069)                   il     = MAX( il, 1 )
( 1070)                   
( 1071)                   ud     = zustar(ji,jj) - ustmin
( 1072)                   jl     = INT( ud / deustar )
( 1073)                   jl     = MIN( jl, njlktbm1 )
( 1074)                   jl     = MAX( jl, 1 )
( 1075)                   
( 1076)                   zfrac  = zd / dezehat - FLOAT( il )  
( 1077)                   ufrac  = ud / deustar - FLOAT( jl )
( 1078)                   zwas   = ( 1. - zfrac ) * wslktb(il,jl+1) + zfrac * wslktb(il+1,jl+1)
( 1079)                   zwbs   = ( 1. - zfrac ) * wslktb(il,jl  ) + zfrac * wslktb(il+1,jl  )
( 1080)                   zwam   = ( 1. - zfrac ) * wmlktb(il,jl+1) + zfrac * wmlktb(il+1,jl+1)
( 1081)                   zwbm   = ( 1. - zfrac ) * wmlktb(il,jl  ) + zfrac * wmlktb(il+1,jl  )
( 1082)                   !
( 1083)                   zws    = ( 1. - ufrac ) * zwbs + ufrac * zwas
( 1084)                   zwm    = ( 1. - ufrac ) * zwbm + ufrac * zwam
( 1085) #else
( 1086)                   ! use analytical functions
( 1087)                   zconm  = 0.5 + SIGN( 0.5, ( rzetam - zeta) )
( 1088)                   zcons  = 0.5 + SIGN( 0.5, ( rzetas - zeta) )
( 1089)                   
( 1090)                   ! Momentum : zeta < rzetam (zconm = 1)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 20

( 1091)                   ! Scalars  : zeta < rzetas (zcons = 1) 
( 1092)                   zwconm = zustar(ji,jj) * vonk * ( ( ABS( rconam - rconcm * zeta) )**pthird )
( 1093)                   zwcons = zustar(ji,jj) * vonk * ( ( ABS( rconas - rconcs * zeta) )**pthird )
( 1094)                   
( 1095)                   ! Momentum : rzetam <= zeta < 0 (zconm = 0)
( 1096)                   ! Scalars  : rzetas <= zeta < 0 (zcons = 0)	
( 1097)                   zwmun  = SQRT( ABS( 1.0 - rconc2 * zeta ) )
( 1098)                   zwsun  = vonk * zustar(ji,jj) * zwmun
( 1099)                   zwmun  = vonk * zustar(ji,jj) * SQRT(zwmun)
( 1100)                   !
( 1101)                   zwm    = zconm * zwconm + ( 1.0 - zconm ) * zwmun
( 1102)                   zws    = zcons * zwcons + ( 1.0 - zcons ) * zwsun
( 1103)                   
( 1104) #endif
( 1105)                ENDIF
( 1106)                
( 1107)                zblcm(ji,jk) = zhbl(ji) * zwm * zsig  * ( 1.0 + zsig * ( za2m(ji) + zsig * za3m(ji) ) )
( 1108)                zblct(ji,jk) = zhbl(ji) * zws * zsig  * ( 1.0 + zsig * ( za2t(ji) + zsig * za3t(ji) ) )
( 1109) #if defined key_zdfddm
( 1110)                zblcs(ji,jk) = zhbl(ji) * zws * zsig  * ( 1.0 + zsig * ( za2s(ji) + zsig * za3s(ji) ) )
( 1111) #endif              
( 1112)                !  Compute Nonlocal transport term = ghats * <ws>o
( 1113)                ! ----------------------------------------------------
( 1114)                ghats(ji,jj,jk-1) = ( 1. - zstabl ) * rcg / ( zws * zhbl(ji) + epsln ) * tmask(ji,jj,jk)
( 1115) 
( 1116)             END DO
( 1117)          END DO     
( 1118)          !  Combine interior and boundary layer coefficients and nonlocal term
( 1119)          ! -----------------------------------------------------------------------
( 1120)          DO jk = 2, jpkm1   
( 1121)             DO ji = fs_2, fs_jpim1
( 1122)                zflag = zmask(ji,jk) * zmask(ji,jk+1)
( 1123)                zviscos(ji,jj,jk) = ( 1.0 - zmask(ji,jk) )         * avmu (ji,jj,jk) & ! interior viscosities
( 1124)                   &              +                        zflag   * zblcm(ji,jk    ) & ! boundary layer viscosities
( 1125)                   &              + zmask(ji,jk) * ( 1.0 - zflag ) * zkmpm(ji       )   ! viscosity enhancement at W_level ne
( 1126)                
( 1127)                zviscos(ji,jj,jk) = zviscos(ji,jj,jk) * tmask(ji,jj,jk)   
( 1128) 
( 1129)             
( 1130)                zdiffut(ji,jj,jk) = ( 1.0 - zmask(ji,jk) )          * avt (ji,jj,jk) & ! interior diffusivities 
( 1131)                   &              +                        zflag   * zblct(ji,jk   ) & ! boundary layer diffusivities
( 1132)                   &              + zmask(ji,jk) * ( 1.0 - zflag ) * zkmpt(ji      )   ! diffusivity enhancement at W_level n
( 1133)                        
( 1134)                zdiffut(ji,jj,jk) = zdiffut(ji,jj,jk) * tmask(ji,jj,jk) 
( 1135) #if defined key_zdfddm
( 1136)                zdiffus(ji,jj,jk) = ( 1.0 - zmask(ji,jk) )          * avs (ji,jj,jk) & ! interior diffusivities 
( 1137)                   &              +                        zflag   * zblcs(ji,jk   ) & ! boundary layer diffusivities
( 1138)                   &              + zmask(ji,jk) * ( 1.0 - zflag ) * zkmps(ji      )   ! diffusivity enhancement at W_level n
( 1139)                        
( 1140)                zdiffus(ji,jj,jk) = zdiffus(ji,jj,jk) * tmask(ji,jj,jk) 
( 1141) #endif               
( 1142)                ! Non local flux in the boundary layer only
( 1143)                ghats(ji,jj,jk-1) = zmask(ji,jk) * ghats(ji,jj,jk-1)
( 1144) 
( 1145)             ENDDO
( 1146)          END DO
( 1147)          !                                             ! ===============
( 1148)       END DO                                           !   End of slab






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 21

( 1149)       !                                                ! ===============
( 1150) 
( 1151)       ! Lateral boundary conditions on zvicos and zdiffus  (sign unchanged)
( 1152)       CALL lbc_lnk( zviscos(:,:,:), 'U', 1. )  ; CALL lbc_lnk( zdiffut(:,:,:), 'W', 1. )  
( 1153) #if defined key_zdfddm  
( 1154)       CALL lbc_lnk( zdiffus(:,:,:), 'W', 1. ) 
( 1155) #endif
( 1156) 
( 1157)       SELECT CASE ( nave )
( 1158)          !
( 1159)          CASE ( 0 )             ! no viscosity and diffusivity smoothing
( 1160) 
( 1161)             DO jk = 2, jpkm1
( 1162)                DO jj = 2, jpjm1
( 1163)                   DO ji = fs_2, fs_jpim1
( 1164)                      avmu(ji,jj,jk) = ( zviscos(ji,jj,jk) + zviscos(ji+1,jj,jk) ) &
( 1165)                         &  / MAX( 1., tmask(ji,jj,jk) + tmask (ji + 1,jj,jk) ) * umask(ji,jj,jk)
( 1166)                      
( 1167)                      avmv(ji,jj,jk) = ( zviscos(ji,jj,jk) + zviscos(ji,jj+1,jk) ) &
( 1168)                         &  / MAX( 1., tmask(ji,jj,jk) + tmask (ji,jj+1,jk) ) * vmask(ji,jj,jk)
( 1169)                      
( 1170)                      avt (ji,jj,jk) =  zdiffut(ji,jj,jk) * tmask(ji,jj,jk)  
( 1171) #if defined key_zdfddm     
( 1172)                      avs (ji,jj,jk) =  zdiffus(ji,jj,jk) * tmask(ji,jj,jk)  
( 1173) #endif
( 1174)                   END DO
( 1175)                END DO
( 1176)             END DO
( 1177)             
( 1178)          CASE ( 1 )                ! viscosity and diffusivity smoothing
( 1179)             !                      
( 1180)             !           ( 1/2  1  1/2 )              ( 1/2  1/2 )             ( 1/2  1  1/2 )
( 1181)             ! avt = 1/8 ( 1    2  1   )   avmu = 1/4 ( 1    1   )   avmv= 1/4 ( 1/2  1  1/2 )
( 1182)             !           ( 1/2  1  1/2 )              ( 1/2  1/2 )
( 1183)   
( 1184)             DO jk = 2, jpkm1
( 1185)                DO jj = 2, jpjm1
( 1186)                   DO ji = fs_2, fs_jpim1
( 1187) 
( 1188)                      avmu(ji,jj,jk) = (      zviscos(ji  ,jj  ,jk) + zviscos(ji+1,jj  ,jk)   &
( 1189)                         &              +.5*( zviscos(ji  ,jj-1,jk) + zviscos(ji+1,jj-1,jk)   &
( 1190)                         &                   +zviscos(ji  ,jj+1,jk) + zviscos(ji+1,jj+1,jk) ) ) * eumean(ji,jj,jk)
( 1191)                      
( 1192)                      avmv(ji,jj,jk) = (      zviscos(ji  ,jj  ,jk) + zviscos(ji  ,jj+1,jk)   &
( 1193)                         &              +.5*( zviscos(ji-1,jj  ,jk) + zviscos(ji-1,jj+1,jk)   &
( 1194)                         &                   +zviscos(ji+1,jj  ,jk) + zviscos(ji+1,jj+1,jk) ) ) * evmean(ji,jj,jk)
( 1195)  
( 1196)                      avt (ji,jj,jk) = ( .5*( zdiffut(ji-1,jj+1,jk) + zdiffut(ji-1,jj-1,jk)    &
( 1197)                         &                   +zdiffut(ji+1,jj+1,jk) + zdiffut(ji+1,jj-1,jk) )  &
( 1198)                         &              +1.*( zdiffut(ji-1,jj  ,jk) + zdiffut(ji  ,jj+1,jk)    &
( 1199)                         &                   +zdiffut(ji  ,jj-1,jk) + zdiffut(ji+1,jj  ,jk) )  &
( 1200)                         &              +2.*  zdiffut(ji  ,jj  ,jk)                          ) * etmean(ji,jj,jk)
( 1201) #if defined key_zdfddm   
( 1202)                      avs (ji,jj,jk) = ( .5*( zdiffus(ji-1,jj+1,jk) + zdiffus(ji-1,jj-1,jk)    &
( 1203)                         &                   +zdiffus(ji+1,jj+1,jk) + zdiffus(ji+1,jj-1,jk) )  &
( 1204)                         &              +1.*( zdiffus(ji-1,jj  ,jk) + zdiffus(ji  ,jj+1,jk)    &
( 1205)                         &                   +zdiffus(ji  ,jj-1,jk) + zdiffus(ji+1,jj  ,jk) )  &
( 1206)                         &              +2.*  zdiffus(ji  ,jj  ,jk)                          ) * etmean(ji,jj,jk)  






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 22

( 1207) #endif                
( 1208)                   END DO
( 1209)                END DO
( 1210)             END DO
( 1211)          
( 1212)          END SELECT
( 1213) 
( 1214)          DO jk = 2, jpkm1                       ! vertical slab
( 1215)             !
( 1216)             !  Minimum value on the eddy diffusivity
( 1217)             ! ----------------------------------------
( 1218)             DO jj = 2, jpjm1
( 1219)                DO ji = fs_2, fs_jpim1   ! vector opt.
( 1220)                   avt(ji,jj,jk) = MAX( avt(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
( 1221) #if defined key_zdfddm  
( 1222)                   avs(ji,jj,jk) = MAX( avs(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
( 1223) #endif
( 1224)                END DO
( 1225)             END DO
( 1226) 
( 1227)             !
( 1228)             ! Minimum value on the eddy viscosity
( 1229)             ! ----------------------------------------
( 1230)             DO jj = 1, jpj
( 1231)                DO ji = 1, jpi
( 1232)                   avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk), avmb(jk) ) * umask(ji,jj,jk)
( 1233)                   avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk), avmb(jk) ) * vmask(ji,jj,jk)
( 1234)                END DO
( 1235)             END DO
( 1236)             !
( 1237)          END DO
( 1238) 
( 1239)          ! Lateral boundary conditions on avt  (sign unchanged)
( 1240)          CALL lbc_lnk( hkpp(:,:), 'T', 1. )
( 1241) 
( 1242)          ! Lateral boundary conditions on avt  (sign unchanged)
( 1243)          CALL lbc_lnk( avt(:,:,:), 'W', 1. )
( 1244) #if defined key_zdfddm  
( 1245)          CALL lbc_lnk( avs(:,:,:), 'W', 1. )  
( 1246) #endif
( 1247)          ! Lateral boundary conditions (avmu,avmv) (U- and V- points, sign unchanged)
( 1248)          CALL lbc_lnk( avmu(:,:,:), 'U', 1. )   ;    CALL lbc_lnk( avmv(:,:,:), 'V', 1. )  
( 1249)  
( 1250)          IF(ln_ctl) THEN
( 1251) #if defined key_zdfddm
( 1252)             CALL prt_ctl(tab3d_1=avt  , clinfo1=' kpp  - t: ', tab3d_2=avs , clinfo2=' s: ', &
( 1253)                &         ovlap=1, kdim=jpk)
( 1254) #else
( 1255)             CALL prt_ctl(tab3d_1=avt  , clinfo1=' kpp  - t: ', ovlap=1, kdim=jpk)
( 1256) #endif
( 1257)             CALL prt_ctl(tab3d_1=avmu  , clinfo1='       u: ', tab3d_2=avmv , clinfo2=' v: ', &
( 1258)                &         ovlap=1, kdim=jpk)
( 1259)          ENDIF
( 1260) 
( 1261)    END SUBROUTINE zdf_kpp
( 1262) 
( 1263) 
( 1264) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 23

( 1265)    SUBROUTINE zdf_kpp_init
( 1266)       !!----------------------------------------------------------------------
( 1267)       !!                  ***  ROUTINE zdf_kpp_init  ***
( 1268)       !!                     
( 1269)       !! ** Purpose :   Initialization of the vertical eddy diffivity and 
( 1270)       !!      viscosity when using a kpp turbulent closure scheme
( 1271)       !!
( 1272)       !! ** Method  :   Read the namkpp namelist and check the parameters
( 1273)       !!      called at the first timestep (nit000)
( 1274)       !!
( 1275)       !! ** input   :   Namlist namkpp
( 1276)       !!
( 1277)       !!
( 1278)       !! history :
( 1279)       !!     8.1  ! 00-02 (J. Chanut) KPP Mixing
( 1280)       !!     9.0  ! 05-01 (C. Ethe) F90 : free form
( 1281)       !!----------------------------------------------------------------------
( 1282)       !! * local declarations
( 1283) 
( 1284)       INTEGER    ::   &
( 1285)          ji, jj, jk             ! dummy loop indices
( 1286)       
( 1287) #if ! defined key_kppcustom
( 1288)       INTEGER    ::   &
( 1289)          jm                       ! dummy loop indices     
( 1290)       REAL(wp)   ::              & !!! tempory scalars
( 1291)          zref, zdist
( 1292) #endif
( 1293) 
( 1294) #if defined key_kpplktb
( 1295)       REAL(wp)   ::              & !!! tempory scalars
( 1296)          zustar,    &
( 1297)          zucube, zustvk,         & 
( 1298)          zeta, zehat
( 1299) #endif
( 1300)       REAL(wp)   ::             & !!! tempory scalars
( 1301)          zhbf
( 1302)       LOGICAL ::                &
( 1303)          ll_kppcustom,          &  ! 1st ocean level taken as surface layer
( 1304)          ll_kpplktb                ! Lookup table for turbul. velocity scales 
( 1305) 
( 1306)       NAMELIST/namkpp/ ln_kpprimix, difmiw, difsiw, Riinfty, difri, bvsqcon, difcon, nave, navb
( 1307) 
( 1308)       !!----------------------------------------------------------------------
( 1309) 
( 1310)       ! Read Namelist namtke : Turbulente Kinetic Energy
( 1311)       ! --------------------
( 1312)       REWIND ( numnam )
( 1313)       READ   ( numnam, namkpp )
( 1314) 
( 1315)       ! Parameter control and print
( 1316)       ! ---------------------------
( 1317)       REWIND( numnam )
( 1318)       READ  ( numnam, namkpp )
( 1319) 
( 1320)       ! Control print
( 1321)       IF(lwp) THEN
( 1322)          WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 24

( 1323)          WRITE(numout,*) 'zdf_kpp_init : kpp turbulent closure scheme'
( 1324)          WRITE(numout,*) '~~~~~~~~~~~~'
( 1325)          WRITE(numout,*) '   Namelist namkpp : set tke mixing parameters'
( 1326)          WRITE(numout,*) '     Shear instability mixing                 ln_kpprimix  = ', ln_kpprimix
( 1327)          WRITE(numout,*) '     max. internal wave viscosity                   difmiw = ', difmiw
( 1328)          WRITE(numout,*) '     max. internal wave diffusivity                 difsiw = ', difsiw
( 1329)          WRITE(numout,*) '     Richardson Number limit for shear instability Riinfty = ', Riinfty
( 1330)          WRITE(numout,*) '     max. shear mixing at Rig = 0                  difri   = ', difri
( 1331)          WRITE(numout,*) '     Brunt-Vaisala squared for max. convection     bvsqcon = ', bvsqcon
( 1332)          WRITE(numout,*) '     max. mix. in interior convec.                 difcon  = ', difcon
( 1333)          WRITE(numout,*) '     horizontal average flag                        nave   = ', nave
( 1334)          WRITE(numout,*) '     constant background or profile                 navb   = ', navb
( 1335)       ENDIF
( 1336) 
( 1337)       ll_kppcustom = .FALSE.
( 1338)       ll_kpplktb   = .FALSE.
( 1339) 
( 1340) #if defined key_kppcustom
( 1341)       ll_kppcustom = .TRUE.
( 1342) #endif
( 1343) #if defined key_kpplktb
( 1344)       ll_kpplktb   = .TRUE.
( 1345) #endif
( 1346)       IF(lwp) THEN
( 1347)          WRITE(numout,*) '     Lookup table for turbul. velocity scales ll_kpplktb   = ', ll_kpplktb
( 1348)          WRITE(numout,*) '     1st ocean level taken as surface layer   ll_kppcustom = ', ll_kppcustom
( 1349)          WRITE(numout,*) ' '
( 1350)       ENDIF
( 1351)       
( 1352)       IF( lk_zdfddm) THEN
( 1353)          IF(lwp) THEN
( 1354)             WRITE(numout,*)
( 1355)             WRITE(numout,*) '    Double diffusion mixing on temperature and salinity '
( 1356)             WRITE(numout,*) '    CAUTION : done in routine zdfkpp, not in routine zdfddm '
( 1357)             WRITE(numout,*) ' '
( 1358)          ENDIF
( 1359)       ENDIF
( 1360)       
( 1361)          
( 1362) 
( 1363)       !set constants not in namelist
( 1364)       !-----------------------------
( 1365)       Vtc  = rconcv * SQRT( 0.2 / ( rconcs * epsilon ) ) / ( vonk * vonk * Ricr )
( 1366)       rcg  = rcstar * vonk * ( rconcs * vonk * epsilon )**pthird
( 1367) 
( 1368)       IF(lwp) THEN
( 1369)          WRITE(numout,*) '     Constant value for unreso. turbul. velocity shear Vtc = ', Vtc
( 1370)          WRITE(numout,*) '     Non-dimensional coef. for nonlocal transport      rcg = ', rcg
( 1371)          WRITE(numout,*) ' '
( 1372)        ENDIF
( 1373) 
( 1374)       ! ratt is the attenuation coefficient for solar flux
( 1375)       ! Should be different is s_coordinate
( 1376)       DO jk = 1, jpk
( 1377)          zhbf     = - fsdept(1,1,jk) * hbf
( 1378)          ratt(jk) = 1.0 - ( rabs * EXP( zhbf / xsi1 ) + ( 1.0 - rabs ) * EXP( zhbf / xsi2 ) )        
( 1379)       ENDDO
( 1380) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 25

( 1381)       ! Horizontal average : initialization of weighting arrays 
( 1382)       ! -------------------
( 1383)       
( 1384)       SELECT CASE ( nave )
( 1385) 
( 1386)       CASE ( 0 )                ! no horizontal average
( 1387)          IF(lwp) WRITE(numout,*) '          no horizontal average on avt, avmu, avmv'
( 1388)          IF(lwp) WRITE(numout,*) '          only in very high horizontal resolution !'
( 1389)          ! weighting mean arrays etmean, eumean and evmean
( 1390)          !           ( 1  1 )                                          ( 1 )
( 1391)          ! avt = 1/4 ( 1  1 )     avmu = 1/2 ( 1  1 )       avmv=  1/2 ( 1 )
( 1392)          !                          
( 1393)          etmean(:,:,:) = 0.e0
( 1394)          eumean(:,:,:) = 0.e0
( 1395)          evmean(:,:,:) = 0.e0
( 1396)          
( 1397)          DO jk = 1, jpkm1
( 1398)             DO jj = 2, jpjm1
( 1399)                DO ji = 2, jpim1   ! vector opt.
( 1400)                   etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
( 1401)                   &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
( 1402)                   &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
( 1403)                   
( 1404)                   eumean(ji,jj,jk) = umask(ji,jj,jk)                     &
( 1405)                   &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk)  )
( 1406) 
( 1407)                   evmean(ji,jj,jk) = vmask(ji,jj,jk)                     &
( 1408)                   &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk)  )
( 1409)                END DO
( 1410)             END DO
( 1411)          END DO
( 1412) 
( 1413)       CASE ( 1 )                ! horizontal average 
( 1414)          IF(lwp) WRITE(numout,*) '          horizontal average on avt, avmu, avmv'
( 1415)          ! weighting mean arrays etmean, eumean and evmean
( 1416)          !           ( 1/2  1  1/2 )              ( 1/2  1/2 )             ( 1/2  1  1/2 )
( 1417)          ! avt = 1/8 ( 1    2  1   )   avmu = 1/4 ( 1    1   )   avmv= 1/4 ( 1/2  1  1/2 )
( 1418)          !           ( 1/2  1  1/2 )              ( 1/2  1/2 )
( 1419)          etmean(:,:,:) = 0.e0
( 1420)          eumean(:,:,:) = 0.e0
( 1421)          evmean(:,:,:) = 0.e0
( 1422)          
( 1423)          DO jk = 1, jpkm1
( 1424)             DO jj = 2, jpjm1
( 1425)                DO ji = fs_2, fs_jpim1   ! vector opt.
( 1426)                   etmean(ji,jj,jk) = tmask(ji, jj,jk)                           &
( 1427)                      & / MAX( 1., 2.* tmask(ji,jj,jk)                           &
( 1428)                      &      +.5 * ( tmask(ji-1,jj+1,jk) + tmask(ji-1,jj-1,jk)   &
( 1429)                      &             +tmask(ji+1,jj+1,jk) + tmask(ji+1,jj-1,jk) ) &
( 1430)                      &      +1. * ( tmask(ji-1,jj  ,jk) + tmask(ji  ,jj+1,jk)   &
( 1431)                      &             +tmask(ji  ,jj-1,jk) + tmask(ji+1,jj  ,jk) ) )
( 1432)                   
( 1433)                   eumean(ji,jj,jk) = umask(ji,jj,jk)                        &
( 1434)                      &  / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
( 1435)                      &       +.5 * ( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &
( 1436)                      &              +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )
( 1437)                   
( 1438)                   evmean(ji,jj,jk) = vmask(ji,jj,jk)                        &






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 26

( 1439)                      &  / MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
( 1440)                      &       +.5 * ( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
( 1441)                      &              +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
( 1442)                END DO
( 1443)             END DO
( 1444)          END DO
( 1445) 
( 1446)       CASE DEFAULT
( 1447)          IF(lwp) WRITE(numout,cform_err)
( 1448)          IF(lwp) WRITE(numout,*) '          bad flag value for nave = ', nave
( 1449)          nstop = nstop + 1
( 1450) 
( 1451)       END SELECT
( 1452) 
( 1453) 
( 1454)  
( 1455)       ! Background eddy viscosity and diffusivity profil
( 1456)       ! ------------------------------------------------
( 1457)       IF( navb == 0 ) THEN
( 1458)          !   Define avmb, avtb from namelist parameter
( 1459)          avmb(:) = avm0
( 1460)          avtb(:) = avt0
( 1461)       ELSE
( 1462)          !   Background profile of avt (fit a theoretical/observational profile (Krauss 1990) 
( 1463)          avmb(:) = avm0
( 1464)          avtb(:) = 1.e-5 + 2.8e-8 * gdepw(:)   ! m2/s
( 1465)       ENDIF
( 1466) 
( 1467)       !   Increase the background in the surface layers
( 1468) !!      avmb(1) = 1. * avmb(1)      ;      avtb(1) = 1. * avtb(1)
( 1469) !!      avmb(2) = 1. * avmb(2)      ;      avtb(2) = 1. * avtb(2)
( 1470) !!      avmb(3) = 1. * avmb(3)      ;      avtb(3) = 1. * avtb(3)
( 1471) !!      avmb(4) = 1. * avmb(4)      ;      avtb(4) = 1. * avtb(4)
( 1472) 
( 1473) 
( 1474)       ! Initialization of vertical eddy coef. to the background value
( 1475)       ! -------------------------------------------------------------
( 1476)       DO jk = 1, jpk
( 1477)          avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
( 1478)          avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
( 1479)          avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
( 1480)       END DO
( 1481) 
( 1482)       ! zero the surface flux for non local term and kpp mixed layer depth
( 1483)       ! ------------------------------------------------------------------
( 1484)       ghats(:,:,:) = 0.
( 1485)       wt0  (:,:  ) = 0.
( 1486)       ws0  (:,:  ) = 0.
( 1487)       hkpp (:,:  ) = 0. ! just a diagnostic (not essential)
( 1488) 
( 1489) #if ! defined key_kppcustom
( 1490)       ! compute arrays (del, wz) for reference mean values 
( 1491)       ! (increase speed for vectorization key_kppcustom not defined)
( 1492)       del(1:jpk, 1:jpk) = 0.
( 1493)       DO jk = 1, jpk
( 1494)          zref = epsilon * fsdept(1,1,jk)    
( 1495)          DO jm = 1 , jpk
( 1496)             zdist = zref - fsdepw(1,1,jm)   






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 27

( 1497)             IF( zdist > 0.  ) THEN
( 1498)                del(jk,jm) = MIN( zdist, fse3t(1,1,jm) ) / zref   
( 1499)             ELSE
( 1500)                del(jk,jm) = 0.
( 1501)             ENDIF
( 1502)          ENDDO
( 1503)       ENDDO
( 1504) #endif
( 1505) 
( 1506) #if defined key_kpplktb
( 1507)       ! build lookup table for turbulent velocity scales
( 1508)       dezehat = ( dehatmax - dehatmin ) / nilktbm1
( 1509)       deustar = ( ustmax   - ustmin   ) / njlktbm1
( 1510)  
( 1511)       DO jj = 1, njlktb
( 1512)          zustar = ( jj - 1) * deustar + ustmin
( 1513)          zustvk = vonk * zustar 
( 1514)          zucube = zustar * zustar * zustar 
( 1515)          DO ji = 1 , nilktb
( 1516)             zehat = ( ji - 1 ) * dezehat + dehatmin
( 1517)             zeta   = zehat / ( zucube + epsln )
( 1518)             IF( zehat >= 0 ) THEN             ! Stable case
( 1519)                wmlktb(ji,jj) = zustvk / ABS( 1.0 + rconc1 * zeta + epsln )                        
( 1520)                wslktb(ji,jj) = wmlktb(ji,jj)
( 1521)             ELSE                                ! Unstable case 
( 1522)                IF( zeta > rzetam ) THEN
( 1523)                   wmlktb(ji,jj) = zustvk * ABS( 1.0    - rconc2 * zeta )**pfourth
( 1524)                ELSE
( 1525)                   wmlktb(ji,jj) = zustvk * ABS( rconam - rconcm * zeta )**pthird
( 1526)                ENDIF
( 1527)                
( 1528)                IF( zeta > rzetas ) THEN
( 1529)                   wslktb(ji,jj) = zustvk * SQRT( ABS( 1.0 - rconc2 * zeta ) )
( 1530)                ELSE
( 1531)                   wslktb(ji,jj) = zustvk * ABS( rconas - rconcs * zeta )**pthird
( 1532)                ENDIF
( 1533)             ENDIF
( 1534)          END DO
( 1535)       END DO
( 1536) #endif
( 1537) 
( 1538)    END SUBROUTINE zdf_kpp_init
( 1539) 
( 1540) #else
( 1541)    !!----------------------------------------------------------------------
( 1542)    !!   Dummy module :                                        NO KPP scheme
( 1543)    !!----------------------------------------------------------------------
( 1544)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfkpp = .FALSE.   !: KPP flag
( 1545) CONTAINS
( 1546)    SUBROUTINE zdf_kpp( kt )          ! Empty routine
( 1547) !      WRITE(*,*) 'zdf_kpp: You should not have seen this print! error?', kt
( 1548)    END SUBROUTINE zdf_kpp
( 1549) #endif
( 1550) 
( 1551)    !!======================================================================
( 1552) END MODULE zdfkpp








PGF90 (Version     10.5)          09/23/2010  09:49:30      page 28

(    1) # 1 "zdfkpp.F90"
(    1) MODULE zdfkpp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfkpp  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the KPP 
(    5)    !!                 turbulent closure parameterization
(    6)    !!=====================================================================
(    7) # 1541
( 1541)    !!----------------------------------------------------------------------
( 1542)    !!   Dummy module :                                        NO KPP scheme
( 1543)    !!----------------------------------------------------------------------
( 1544)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfkpp = .FALSE.   !: KPP flag
( 1545) CONTAINS
( 1546)    SUBROUTINE zdf_kpp( kt )          ! Empty routine
( 1547) !      WRITE(*,*) 'zdf_kpp: You should not have seen this print! error?', kt
( 1548)    END SUBROUTINE zdf_kpp
( 1549) 
( 1550) # 1551
( 1551)    !!======================================================================
( 1552) END MODULE zdfkpp
