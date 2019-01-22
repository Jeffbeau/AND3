


PGF90 (Version     10.5)          09/23/2010  09:49:31      page 1

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

Filename: zdftke.F90

(    1) MODULE zdftke
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdftke  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the tke 
(    5)    !!                 turbulent closure parameterization
(    6)    !!=====================================================================
(    7) #if defined key_zdftke   ||   defined key_esopa
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_zdftke'                                             TKE scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !!   zdf_tke      : update momentum and tracer Kz from a tke scheme
(   12)    !!   zdf_tke_init : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and active tracers 
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE zdf_oce         ! ocean vertical physics
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE phycst          ! physical constants
(   21)    USE taumod          ! surface stress
(   22)    USE prtctl          ! Print control
(   23)    USE traadv_ctl      ! advection scheme control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC zdf_tke   ! routine called by step.F90
(   30) 
(   31)    !! * Share Module variables
(   32)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .TRUE.    !: TKE vertical mixing flag
(   33)    LOGICAL, PUBLIC ::         & !!: ** tke namelist (namtke) **
(   34)      ln_rstke = .FALSE.          !: =T restart with tke from a run without tke with 
(   35)      !                           !  a none zero initial value for en
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   37)       en                         !: now turbulent kinetic energy
(   38) 
(   39)    !! * Module variables
(   40)    INTEGER ::                 & !!! ** tke namelist (namtke) **
(   41)       nitke = 50 ,            &  ! number of restart iterative loops
(   42)       nmxl  =  2 ,            &  ! = 0/1/2/3 flag for the type of mixing length used
(   43)       npdl  =  1 ,            &  ! = 0/1/2 flag on prandtl number on vert. eddy coeff.
(   44)       nave  =  1 ,            &  ! = 0/1 flag for horizontal average on avt, avmu, avmv
(   45)       navb  =  0                 ! = 0/1 flag for constant or profile background avt
(   46)    REAL(wp) ::                & !!! ** tke namlist (namtke) **






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 2

(   47)       ediff = 0.1_wp       ,  &  ! coeff. for vertical eddy coef.; avt=ediff*mxl*sqrt(e)
(   48)       ediss = 0.7_wp       ,  &  ! coef. of the Kolmogoroff dissipation 
(   49)       ebb   = 3.75_wp      ,  &  ! coef. of the surface input of tke
(   50)       efave = 1._wp        ,  &  ! coef. for the tke vert. diff. coeff.; avtke=efave*avm
(   51)       emin  = 0.7071e-6_wp ,  &  ! minimum value of tke (m2/s2)
(   52)       emin0 = 1.e-4_wp     ,  &  ! surface minimum value of tke (m2/s2)
(   53)       ri_c  = 2._wp / 9._wp      ! critic Richardson number
(   54)    REAL(wp) ::   &  
(   55)       eboost        ! multiplicative coeff of the shear product.
(   56) 
(   57)    !! caution vectopt_memory change the solution (last digit of the solver stat)
(   58) # if defined key_vectopt_memory
(   59)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   60)       etmean,    &  ! coefficient used for horizontal smoothing
(   61)       eumean,    &  ! at t-, u- and v-points
(   62)       evmean        !
(   63) # endif
(   64) 
(   65) # if defined key_cfg_1d
(   66)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   
(   67)       e_dis,    &   ! dissipation turbulent lengh scale
(   68)       e_mix,    &   ! mixing turbulent lengh scale
(   69)       e_pdl,    &   ! prandl number
(   70)       e_ric         ! local Richardson number
(   71) #endif
(   72) 
(   73)    !! * Substitutions
(   74) #  include "domzgr_substitute.h90"
(   75) #  include "vectopt_loop_substitute.h90"
(   76)    !!----------------------------------------------------------------------
(   77)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   78)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdftke.F90,v 1.9 2006/03/21 09:16:38 opalod Exp $ 
(   79)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   80)    !!----------------------------------------------------------------------
(   81) 
(   82) CONTAINS
(   83) 
(   84) # if defined key_autotasking
(   85)    !!----------------------------------------------------------------------
(   86)    !!   'key_autotasking' :                             j-k-i loop (j-slab)
(   87)    !!----------------------------------------------------------------------
(   88) #  include "zdftke_atsk.h90"
(   89) 
(   90) # else
(   91)    !!----------------------------------------------------------------------
(   92)    !!   default option :                                         k-j-i loop
(   93)    !!----------------------------------------------------------------------
(   94) 
(   95)    SUBROUTINE zdf_tke ( kt )
(   96)       !!----------------------------------------------------------------------
(   97)       !!                   ***  ROUTINE zdf_tke  ***
(   98)       !!
(   99)       !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
(  100)       !!      coefficients using a 1.5 turbulent closure scheme.
(  101)       !!
(  102)       !! ** Method  :   The time evolution of the turbulent kinetic energy 
(  103)       !!      (tke) is computed from a prognostic equation :
(  104)       !!         d(en)/dt = eboost eav (d(u)/dz)**2       ! shear production






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 3

(  105)       !!                  + d( efave eav d(en)/dz )/dz    ! diffusion of tke
(  106)       !!                  + grav/rau0 pdl eav d(rau)/dz   ! stratif. destruc.
(  107)       !!                  - ediss / emxl en**(2/3)        ! dissipation
(  108)       !!      with the boundary conditions:
(  109)       !!         surface: en = max( emin0,ebb sqrt(taux^2 + tauy^2) )
(  110)       !!         bottom : en = emin
(  111)       !!      -1- The dissipation and mixing turbulent lengh scales are computed
(  112)       !!      from the usual diagnostic buoyancy length scale:  
(  113)       !!         mxl= 1/(sqrt(en)/N)  WHERE N is the brunt-vaisala frequency
(  114)       !!      Four cases : 
(  115)       !!         nmxl=0 : mxl bounded by the distance to surface and bottom.
(  116)       !!                  zmxld = zmxlm = mxl
(  117)       !!         nmxl=1 : mxl bounded by the vertical scale factor.
(  118)       !!                  zmxld = zmxlm = mxl
(  119)       !!         nmxl=2 : mxl bounded such that the vertical derivative of mxl 
(  120)       !!                  is less than 1 (|d/dz(xml)|<1). 
(  121)       !!                  zmxld = zmxlm = mxl
(  122)       !!         nmxl=3 : lup = mxl bounded using |d/dz(xml)|<1 from the surface 
(  123)       !!                        to the bottom
(  124)       !!                  ldown = mxl bounded using |d/dz(xml)|<1 from the bottom 
(  125)       !!                        to the surface
(  126)       !!                  zmxld = sqrt (lup*ldown) ; zmxlm = min(lup,ldown)
(  127)       !!      -2- Compute the now Turbulent kinetic energy. The time differencing
(  128)       !!      is implicit for vertical diffusion term, linearized for kolmo-
(  129)       !!      goroff dissipation term, and explicit forward for both buoyancy
(  130)       !!      and dynamic production terms. Thus a tridiagonal linear system is
(  131)       !!      solved.
(  132)       !!         Note that - the shear production is multiplied by eboost in order
(  133)       !!      to set the critic richardson number to ri_c (namelist parameter)
(  134)       !!                   - the destruction by stratification term is multiplied
(  135)       !!      by the Prandtl number (defined by an empirical funtion of the local 
(  136)       !!      Richardson number) if npdl=1 (namelist parameter)
(  137)       !!      coefficient (zesh2):
(  138)       !!      -3- Compute the now vertical eddy vicosity and diffusivity
(  139)       !!      coefficients from en (before the time stepping) and zmxlm:
(  140)       !!              avm = max( avtb, ediff*zmxlm*en^1/2 )
(  141)       !!              avt = max( avmb, pdl*avm )  (pdl=1 if npdl=0)
(  142)       !!              eav = max( avmb, avm )
(  143)       !!      avt and avm are horizontally averaged to avoid numerical insta-
(  144)       !!      bilities.
(  145)       !!        N.B. The computation is done from jk=2 to jpkm1 except for
(  146)       !!      en. Surface value of avt avmu avmv are set once a time to
(  147)       !!      their background value in routine zdf_tke_init.
(  148)       !!
(  149)       !! ** Action  :   compute en (now turbulent kinetic energy)
(  150)       !!                update avt, avmu, avmv (before vertical eddy coef.)
(  151)       !!
(  152)       !! References :
(  153)       !!      Gaspar et al., jgr, 95, 1990,
(  154)       !!      Blanke and Delecluse, jpo, 1991
(  155)       !! History :
(  156)       !!   6.0  !  91-03  (b. blanke)  Original code
(  157)       !!   7.0  !  91-11  (G. Madec)   bug fix
(  158)       !!   7.1  !  92-10  (G. Madec)   new mixing length and eav
(  159)       !!   7.2  !  93-03  (M. Guyon)   symetrical conditions
(  160)       !!   7.3  !  94-08  (G. Madec, M. Imbard)   npdl flag
(  161)       !!   7.5  !  96-01  (G. Madec)   s-coordinates
(  162)       !!   8.0  !  97-07  (G. Madec)   lbc






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 4

(  163)       !!   8.1  !  99-01  (E. Stretta) new option for the mixing length
(  164)       !!   8.5  !  02-08  (G. Madec)  ri_c and Free form, F90
(  165)       !!   9.0  !  04-10  (C. Ethe )  1D configuration
(  166)       !!----------------------------------------------------------------------
(  167)       !! * Modules used
(  168)       USE oce     , zwd   => ua,  &  ! use ua as workspace
(  169)          &          zmxlm => ta,  &  ! use ta as workspace
(  170)          &          zmxld => sa      ! use sa as workspace
(  171) 
(  172)       !! * arguments
(  173)       INTEGER, INTENT( in  ) ::   kt ! ocean time step
(  174) 
(  175)       !! * local declarations
(  176)       INTEGER ::   ji, jj, jk        ! dummy loop arguments
(  177)       REAL(wp) ::   &
(  178)          zmlmin, zbbrau,          &  ! temporary scalars
(  179)          zfact1, zfact2, zfact3,  &  !
(  180)          zrn2, zesurf,            &  !
(  181)          ztx2, zty2, zav,         &  !
(  182)          zcoef, zcof, zsh2,       &  !
(  183)          zdku, zdkv, zpdl, zri,   &  !
(  184)          zsqen, zesh2,            &  !
(  185)          zemxl, zemlm, zemlp
(  186)       !!--------------------------------------------------------------------
(  187) 
(  188)       ! Initialization (first time-step only)
(  189)       ! --------------
(  190)       IF( kt == nit000  )   CALL zdf_tke_init
(  191) 
(  192)       ! Local constant initialization
(  193)       zmlmin = 1.e-8
(  194)       zbbrau =  .5 * ebb / rau0
(  195)       zfact1 = -.5 * rdt * efave
(  196)       zfact2 = 1.5 * rdt * ediss
(  197)       zfact3 = 0.5 * rdt * ediss
(  198) 
(  199) 
(  200)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  201)       ! I.  Mixing length
(  202)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  203) 
(  204) 
(  205)       ! Buoyancy length scale: l=sqrt(2*e/n**2)
(  206)       ! ---------------------
(  207)       zmxlm(:,:, 1 ) = zmlmin   ! surface set to the minimum value
(  208)       zmxlm(:,:,jpk) = zmlmin   ! bottom  set to the minimum value
(  209) !CDIR NOVERRCHK
(  210)       DO jk = 2, jpkm1
(  211) !CDIR NOVERRCHK
(  212)          DO jj = 2, jpjm1
(  213) !CDIR NOVERRCHK
(  214)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  215)                zrn2 = MAX( rn2(ji,jj,jk), rsmall )
(  216)                zmxlm(ji,jj,jk) = MAX( SQRT( 2. * en(ji,jj,jk) / zrn2 ), zmlmin  )
(  217)             END DO
(  218)          END DO
(  219)       END DO
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 5

(  221) 
(  222)       ! Physical limits for the mixing length
(  223)       ! -------------------------------------
(  224)       zmxld(:,:, 1 ) = zmlmin   ! surface set to the minimum value
(  225)       zmxld(:,:,jpk) = zmlmin   ! bottom  set to the minimum value
(  226) 
(  227)       SELECT CASE ( nmxl )
(  228) 
(  229)       CASE ( 0 )           ! bounded by the distance to surface and bottom
(  230)          DO jk = 2, jpkm1
(  231)             DO jj = 2, jpjm1
(  232)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  233)                   zemxl = MIN( fsdepw(ji,jj,jk), zmxlm(ji,jj,jk),   &
(  234)                   &            fsdepw(ji,jj,mbathy(ji,jj)) - fsdepw(ji,jj,jk) )
(  235)                   zmxlm(ji,jj,jk) = zemxl
(  236)                   zmxld(ji,jj,jk) = zemxl
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240) 
(  241)       CASE ( 1 )           ! bounded by the vertical scale factor
(  242)          DO jk = 2, jpkm1
(  243)             DO jj = 2, jpjm1
(  244)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  245)                   zemxl = MIN( fse3w(ji,jj,jk), zmxlm(ji,jj,jk) )
(  246)                   zmxlm(ji,jj,jk) = zemxl
(  247)                   zmxld(ji,jj,jk) = zemxl
(  248)                END DO
(  249)             END DO
(  250)          END DO
(  251) 
(  252)       CASE ( 2 )           ! |dk[xml]| bounded by e3t :
(  253)          DO jk = 2, jpkm1         ! from the surface to the bottom :
(  254)             DO jj = 2, jpjm1
(  255)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  256)                   zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk-1) + fse3t(ji,jj,jk-1), zmxlm(ji,jj,jk) )
(  257)                END DO
(  258)             END DO
(  259)          END DO
(  260)          DO jk = jpkm1, 2, -1     ! from the bottom to the surface :
(  261)             DO jj = 2, jpjm1
(  262)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  263)                   zemxl = MIN( zmxlm(ji,jj,jk+1) + fse3t(ji,jj,jk+1), zmxlm(ji,jj,jk) )
(  264)                   zmxlm(ji,jj,jk) = zemxl
(  265)                   zmxld(ji,jj,jk) = zemxl
(  266)                END DO
(  267)             END DO
(  268)          END DO
(  269) 
(  270)       CASE ( 3 )           ! lup and ldown, |dk[xml]| bounded by e3t :
(  271)          DO jk = 2, jpkm1         ! from the surface to the bottom : lup
(  272)             DO jj = 2, jpjm1
(  273)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  274)                   zmxld(ji,jj,jk) = MIN( zmxld(ji,jj,jk-1) + fse3t(ji,jj,jk-1), zmxlm(ji,jj,jk) )
(  275)                END DO
(  276)             END DO
(  277)          END DO
(  278)          DO jk = jpkm1, 2, -1     ! from the bottom to the surface : ldown






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 6

(  279)             DO jj = 2, jpjm1
(  280)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  281)                   zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk+1) + fse3t(ji,jj,jk+1), zmxlm(ji,jj,jk) )
(  282)                END DO
(  283)             END DO
(  284)          END DO
(  285) !CDIR NOVERRCHK
(  286)          DO jk = 2, jpkm1
(  287) !CDIR NOVERRCHK
(  288)             DO jj = 2, jpjm1
(  289) !CDIR NOVERRCHK
(  290)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  291)                   zemlm = MIN ( zmxld(ji,jj,jk),  zmxlm(ji,jj,jk) )
(  292)                   zemlp = SQRT( zmxld(ji,jj,jk) * zmxlm(ji,jj,jk) )
(  293)                   zmxlm(ji,jj,jk) = zemlm
(  294)                   zmxld(ji,jj,jk) = zemlp
(  295)                END DO
(  296)             END DO
(  297)          END DO
(  298) 
(  299)       END SELECT
(  300) 
(  301) # if defined key_cfg_1d
(  302)       ! save mixing and dissipation turbulent length scales
(  303)       e_dis(:,:,:) = zmxld(:,:,:)
(  304)       e_mix(:,:,:) = zmxlm(:,:,:)
(  305) # endif
(  306) 
(  307) 
(  308)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  309)       ! II  Tubulent kinetic energy time stepping
(  310)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  311) 
(  312) 
(  313)       ! 1. Vertical eddy viscosity on tke (put in zmxlm) and first estimate of avt
(  314)       ! ---------------------------------------------------------------------
(  315) !CDIR NOVERRCHK
(  316)       DO jk = 2, jpkm1
(  317) !CDIR NOVERRCHK
(  318)          DO jj = 2, jpjm1
(  319) !CDIR NOVERRCHK
(  320)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  321)                zsqen = SQRT( en(ji,jj,jk) )
(  322)                zav   = ediff * zmxlm(ji,jj,jk) * zsqen
(  323)                avt  (ji,jj,jk) = MAX( zav, avtb(jk) ) * tmask(ji,jj,jk)
(  324)                zmxlm(ji,jj,jk) = MAX( zav, avmb(jk) ) * tmask(ji,jj,jk)
(  325)                zmxld(ji,jj,jk) = zsqen / zmxld(ji,jj,jk)
(  326)             END DO
(  327)          END DO
(  328)       END DO
(  329) 
(  330)       ! 2. Surface boundary condition on tke and its eddy viscosity (zmxlm)
(  331)       ! -------------------------------------------------
(  332)       ! en(1)   = ebb sqrt(taux^2+tauy^2) / rau0  (min value emin0)
(  333)       ! zmxlm(1) = avmb(1) and zmxlm(jpk) = 0.
(  334) !CDIR NOVERRCHK
(  335)       DO jj = 2, jpjm1
(  336) !CDIR NOVERRCHK






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 7

(  337)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  338)             ztx2 = taux(ji-1,jj  ) + taux(ji,jj)
(  339)             zty2 = tauy(ji  ,jj-1) + tauy(ji,jj)
(  340)             zesurf = zbbrau * SQRT( ztx2 * ztx2 + zty2 * zty2 )
(  341)             en (ji,jj,1) = MAX( zesurf, emin0 ) * tmask(ji,jj,1)
(  342)             zmxlm(ji,jj,1  ) = avmb(1) * tmask(ji,jj,1)
(  343)             zmxlm(ji,jj,jpk) = 0.e0
(  344)          END DO
(  345)       END DO
(  346) 
(  347)       ! 3. Now Turbulent kinetic energy (output in en)
(  348)       ! -------------------------------
(  349)       ! Resolution of a tridiagonal linear system by a "methode de chasse"
(  350)       ! computation from level 2 to jpkm1  (e(1) already computed and
(  351)       ! e(jpk)=0 ). 
(  352) 
(  353)       SELECT CASE ( npdl )
(  354) 
(  355)       CASE ( 0 )           ! No Prandtl number
(  356)          DO jk = 2, jpkm1
(  357)             DO jj = 2, jpjm1
(  358)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  359)                   ! zesh2 = eboost * (du/dz)^2 - N^2
(  360)                   zcoef = 0.5 / fse3w(ji,jj,jk)
(  361)                   ! shear
(  362)                   zdku = zcoef * (   ub(ji-1, jj ,jk-1) + ub(ji,jj,jk-1)   &
(  363)                   &                - ub(ji-1, jj ,jk  ) - ub(ji,jj,jk  )  )
(  364)                   zdkv = zcoef * (   vb( ji ,jj-1,jk-1) + vb(ji,jj,jk-1)   &
(  365)                   &                - vb( ji ,jj-1,jk  ) - vb(ji,jj,jk  )  )
(  366)                   ! coefficient (zesh2)
(  367)                   zesh2 =  eboost * ( zdku*zdku + zdkv*zdkv ) - rn2(ji,jj,jk)
(  368) 
(  369)                   ! Matrix
(  370)                   zcof = zfact1 * tmask(ji,jj,jk)
(  371)                   ! lower diagonal
(  372)                   avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
(  373)                   &                    / ( fse3t(ji,jj,jk-1) * fse3w(ji,jj,jk  ) )
(  374)                   ! upper diagonal
(  375)                   avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
(  376)                   &                    / ( fse3t(ji,jj,jk  ) * fse3w(ji,jj,jk) )
(  377)                   ! diagonal
(  378)                   zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
(  379)                   ! right hand side in en 
(  380)                   en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
(  381)                   &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
(  382)                END DO
(  383)             END DO
(  384)          END DO
(  385) 
(  386)       CASE ( 1 )           ! Prandtl number
(  387)          DO jk = 2, jpkm1
(  388)             DO jj = 2, jpjm1
(  389)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  390)                   ! zesh2 =  eboost * (du/dz)^2 - pdl * N^2
(  391)                   zcoef = 0.5 / fse3w(ji,jj,jk)
(  392)                   ! shear
(  393)                   zdku = zcoef * (   ub(ji-1,jj  ,jk-1) + ub(ji,jj,jk-1) &
(  394)                   &                - ub(ji-1,jj  ,jk  ) - ub(ji,jj,jk  )   )






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 8

(  395)                   zdkv = zcoef * (   vb(ji  ,jj-1,jk-1) + vb(ji,jj,jk-1) &
(  396)                   &                - vb(ji  ,jj-1,jk  ) - vb(ji,jj,jk  )   )
(  397)                   ! square of vertical shear
(  398)                   zsh2 = zdku * zdku + zdkv * zdkv
(  399)                   ! local Richardson number
(  400)                   zri  = MAX( rn2(ji,jj,jk), 0. ) / ( zsh2 + 1.e-20 )
(  401) # if defined key_cfg_1d
(  402)                   ! save masked local Richardson number in zmxlm array
(  403)                   e_ric(ji,jj,jk) = zri * tmask(ji,jj,jk)
(  404) # endif
(  405)                   ! Prandtl number
(  406)                   zpdl = 1.0
(  407)                   IF( zri >= 0.2 ) zpdl = 0.2 / zri
(  408)                   zpdl = MAX( 0.1, zpdl )
(  409)                   ! coefficient (esh2)
(  410)                   zesh2 = eboost * zsh2 - zpdl * rn2(ji,jj,jk)
(  411) 
(  412)                   ! Matrix
(  413)                   zcof = zfact1 * tmask(ji,jj,jk)
(  414)                   ! lower diagonal
(  415)                   avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
(  416)                   &                     / ( fse3t(ji,jj,jk-1) * fse3w(ji,jj,jk  ) )
(  417)                   ! upper diagonal
(  418)                   avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
(  419)                   &                     / ( fse3t(ji,jj,jk  ) * fse3w(ji,jj,jk) )
(  420)                   ! diagonal
(  421)                   zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
(  422)                   ! right hand side in en 
(  423)                   en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
(  424)                   &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
(  425)                   ! save masked Prandlt number in zmxlm array
(  426)                   zmxld(ji,jj,jk) = zpdl * tmask(ji,jj,jk)
(  427)                END DO
(  428)             END DO
(  429)          END DO
(  430) 
(  431)       END SELECT
(  432) 
(  433) # if defined key_cfg_1d
(  434)       !  save masked Prandlt number
(  435)       e_pdl(:,:,2:jpkm1) = zmxld(:,:,2:jpkm1)
(  436)       e_pdl(:,:,      1) = e_pdl(:,:,      2)
(  437)       e_pdl(:,:,    jpk) = e_pdl(:,:,  jpkm1)      
(  438) # endif
(  439) 
(  440)       ! 4. Matrix inversion from level 2 (tke prescribed at level 1)
(  441)       !!--------------------------------
(  442)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
(  443)       DO jk = 3, jpkm1
(  444)          DO jj = 2, jpjm1
(  445)             DO ji = fs_2, fs_jpim1    ! vector opt.
(  446)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - avmv(ji,jj,jk) * avmu(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  447)             END DO
(  448)          END DO
(  449)       END DO
(  450) 
(  451)       ! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
(  452)       DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 9

(  453)          DO ji = fs_2, fs_jpim1    ! vector opt.
(  454)             avmv(ji,jj,2) = en(ji,jj,2) - avmv(ji,jj,2) * en(ji,jj,1)    ! Surface boudary conditions on tke
(  455)          END DO
(  456)       END DO
(  457)       DO jk = 3, jpkm1
(  458)          DO jj = 2, jpjm1
(  459)             DO ji = fs_2, fs_jpim1    ! vector opt.
(  460)                avmv(ji,jj,jk) = en(ji,jj,jk) - avmv(ji,jj,jk) / zwd(ji,jj,jk-1) *avmv(ji,jj,jk-1)
(  461)             END DO
(  462)          END DO
(  463)       END DO
(  464) 
(  465)       ! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk
(  466)       DO jj = 2, jpjm1
(  467)          DO ji = fs_2, fs_jpim1    ! vector opt.
(  468)             en(ji,jj,jpkm1) = avmv(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  469)          END DO
(  470)       END DO
(  471)       DO jk = jpk-2, 2, -1
(  472)          DO jj = 2, jpjm1
(  473)             DO ji = fs_2, fs_jpim1    ! vector opt.
(  474)                en(ji,jj,jk) = ( avmv(ji,jj,jk) - avmu(ji,jj,jk) * en(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  475)             END DO
(  476)          END DO
(  477)       END DO
(  478) 
(  479)       ! Save the result in en and set minimum value of tke : emin
(  480)       DO jk = 2, jpkm1
(  481)          DO jj = 2, jpjm1
(  482)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  483)                en(ji,jj,jk) = MAX( en(ji,jj,jk), emin ) * tmask(ji,jj,jk)
(  484)             END DO
(  485)          END DO
(  486)       END DO
(  487) 
(  488)       ! Lateral boundary conditions on ( avt, en )   (sign unchanged)
(  489)       CALL lbc_lnk( en , 'W', 1. )   ;   CALL lbc_lnk( avt, 'W', 1. )
(  490) 
(  491) 
(  492)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  493)       ! III.  Before vertical eddy vicosity and diffusivity coefficients
(  494)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  495) 
(  496)       SELECT CASE ( nave )
(  497)          
(  498)       CASE ( 0 )                ! no horizontal average
(  499) 
(  500)          ! Vertical eddy viscosity
(  501) 
(  502)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  503)             DO jj = 2, jpjm1
(  504)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  505)                   avmu(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji+1,jj  ,jk) ) * umask(ji,jj,jk)   &
(  506)                   &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk) )
(  507)                   avmv(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji  ,jj+1,jk) ) * vmask(ji,jj,jk)   &
(  508)                   &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk) )
(  509)                END DO
(  510)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 10

(  511)          END DO
(  512) 
(  513) !DB: jack up viscosity in shallow water 
(  514) !         do jj = 2, jpjm1
(  515) !            do ji = fs_2, fs_jpim1   ! vector opt.
(  516) !               if(mbathy(ji,jj) <= 5) then
(  517) !                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
(  518) !                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
(  519) !               endif
(  520) !            enddo
(  521) !         enddo
(  522) !!DB: END
(  523) 
(  524)          ! Lateral boundary conditions (avmu,avmv) (U- and V- points, sign unchanged)
(  525)          CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )
(  526)          
(  527)       CASE ( 1 )                ! horizontal average
(  528) 
(  529)          !                                                ( 1/2  1/2 )
(  530)          ! Eddy viscosity: horizontal average: avmu = 1/4 ( 1    1   )
(  531)          !                      ( 1/2  1 1/2 )            ( 1/2  1/2 )
(  532)          !           avmv = 1/4 ( 1/2  1 1/2 )      
(  533)          
(  534) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  535) #  if defined key_vectopt_memory
(  536)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  537)             DO jj = 2, jpjm1
(  538)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  539)                   avmu(ji,jj,jk) = (      avt(ji,jj  ,jk) + avt(ji+1,jj  ,jk)   &
(  540)                   &                 +.5*( avt(ji,jj-1,jk) + avt(ji+1,jj-1,jk)   &
(  541)                   &                      +avt(ji,jj+1,jk) + avt(ji+1,jj+1,jk) ) ) * eumean(ji,jj,jk)
(  542) 
(  543)                   avmv(ji,jj,jk) = (      avt(ji  ,jj,jk) + avt(ji  ,jj+1,jk)   &
(  544)                   &                 +.5*( avt(ji-1,jj,jk) + avt(ji-1,jj+1,jk)   &
(  545)                   &                      +avt(ji+1,jj,jk) + avt(ji+1,jj+1,jk) ) ) * evmean(ji,jj,jk)
(  546)                END DO
(  547)             END DO
(  548)          END DO
(  549) !DB: jack up viscosity in shallow water
(  550) !         do jj = 2, jpjm1
(  551) !            do ji = fs_2, fs_jpim1   ! vector opt.
(  552) !               if(mbathy(ji,jj) <= 5) then
(  553) !                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
(  554) !                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
(  555) !               endif
(  556) !            enddo
(  557) !         enddo
(  558) !!DB: END
(  559) 
(  560) #  else
(  561)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  562)             DO jj = 2, jpjm1
(  563)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  564)                   avmu(ji,jj,jk) = (   avt  (ji,jj  ,jk) + avt  (ji+1,jj  ,jk)   &
(  565)                   &              +.5*( avt  (ji,jj-1,jk) + avt  (ji+1,jj-1,jk)   &
(  566)                   &                   +avt  (ji,jj+1,jk) + avt  (ji+1,jj+1,jk) ) ) * umask(ji,jj,jk)  &
(  567)                   &       / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
(  568)                   &              +.5*( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 11

(  569)                   &                   +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )
(  570) 
(  571)                   avmv(ji,jj,jk) = (   avt  (ji  ,jj,jk) + avt  (ji  ,jj+1,jk)   &
(  572)                   &              +.5*( avt  (ji-1,jj,jk) + avt  (ji-1,jj+1,jk)   &
(  573)                   &                   +avt  (ji+1,jj,jk) + avt  (ji+1,jj+1,jk) ) ) * vmask(ji,jj,jk)  &
(  574)                   &      /  MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
(  575)                   &              +.5*( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
(  576)                   &                   +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
(  577)                END DO
(  578)             END DO
(  579)          END DO
(  580) 
(  581) !DB: jack up viscosity in shallow water
(  582) !         do jj = 2, jpjm1
(  583) !            do ji = fs_2, fs_jpim1   ! vector opt.
(  584) !               if(mbathy(ji,jj) <= 5) then
(  585) !                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
(  586) !                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
(  587) !               endif
(  588) !            enddo
(  589) !         enddo
(  590) !!DB: END
(  591) 
(  592) #  endif
(  593) 
(  594)          ! Lateral boundary conditions (avmu,avmv) (sign unchanged)
(  595)          CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )
(  596) 
(  597)          ! Vertical eddy diffusivity
(  598)          ! ------------------------------
(  599)          !                                (1 2 1)
(  600)          ! horizontal average  avt = 1/16 (2 4 2)
(  601)          !                                (1 2 1)
(  602)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  603) #  if defined key_vectopt_memory
(  604)             DO jj = 2, jpjm1
(  605)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  606)                   avt(ji,jj,jk) = ( avmu(ji,jj,jk) + avmu(ji-1,jj  ,jk)    &
(  607)                   &               + avmv(ji,jj,jk) + avmv(ji  ,jj-1,jk)  ) * etmean(ji,jj,jk)
(  608)                END DO
(  609)             END DO
(  610) #  else
(  611)             DO jj = 2, jpjm1
(  612)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  613)                   avt(ji,jj,jk) = ( avmu (ji,jj,jk) + avmu (ji-1,jj  ,jk)   &
(  614)                   &               + avmv (ji,jj,jk) + avmv (ji  ,jj-1,jk)  ) * tmask(ji,jj,jk)   &
(  615)                   &     / MAX( 1.,  umask(ji,jj,jk) + umask(ji-1,jj  ,jk)   &
(  616)                   &               + vmask(ji,jj,jk) + vmask(ji  ,jj-1,jk)  )
(  617)                END DO
(  618)             END DO
(  619) #  endif
(  620)          END DO
(  621) 
(  622)       END SELECT
(  623) 
(  624)       ! multiplied by the Prandtl number (npdl>1)
(  625)       ! ----------------------------------------
(  626)       IF( npdl == 1 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 12

(  627)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  628)             DO jj = 2, jpjm1
(  629)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  630)                   zpdl = zmxld(ji,jj,jk)
(  631)                   avt(ji,jj,jk) = MAX( zpdl * avt(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
(  632)                END DO
(  633)             END DO
(  634)          END DO
(  635)       ENDIF
(  636) 
(  637)       ! Minimum value on the eddy viscosity
(  638)       ! ----------------------------------------
(  639)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  640)          DO jj = 1, jpj
(  641)             DO ji = 1, jpi
(  642)                avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk), avmb(jk) ) * umask(ji,jj,jk)
(  643)                avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk), avmb(jk) ) * vmask(ji,jj,jk)
(  644)             END DO
(  645)          END DO
(  646)       END DO
(  647) 
(  648)       ! Lateral boundary conditions on avt  (sign unchanged)
(  649)       ! ------------------------------=====
(  650)       CALL lbc_lnk( avt, 'W', 1. )
(  651) 
(  652)       IF(ln_ctl) THEN
(  653)          CALL prt_ctl(tab3d_1=en  , clinfo1=' tke  - e: ', tab3d_2=avt , clinfo2=' t: ', ovlap=1, kdim=jpk)
(  654)          CALL prt_ctl(tab3d_1=avmu, clinfo1=' tke  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
(  655)       ENDIF
(  656) 
(  657)    END SUBROUTINE zdf_tke
(  658) 
(  659) # endif
(  660) 
(  661)    SUBROUTINE zdf_tke_init
(  662)       !!----------------------------------------------------------------------
(  663)       !!                  ***  ROUTINE zdf_tke_init  ***
(  664)       !!                     
(  665)       !! ** Purpose :   Initialization of the vertical eddy diffivity and 
(  666)       !!      viscosity when using a tke turbulent closure scheme
(  667)       !!
(  668)       !! ** Method  :   Read the namtke namelist and check the parameters
(  669)       !!      called at the first timestep (nit000)
(  670)       !!
(  671)       !! ** input   :   Namlist namtke
(  672)       !!
(  673)       !! ** Action  :   Increase by 1 the nstop flag is setting problem encounter
(  674)       !!
(  675)       !! history :
(  676)       !!  8.5  ! 02-06 (G. Madec) original code
(  677)       !!----------------------------------------------------------------------
(  678)       !! * Module used
(  679)       USE dynzdf_exp
(  680)       USE trazdf_exp
(  681) 
(  682)       !! * local declarations
(  683)       !! caution vectopt_memory change the solution (last digit of the solver stat)
(  684) # if defined key_vectopt_memory






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 13

(  685)       INTEGER ::   ji, jj, jk, jit   ! dummy loop indices
(  686) # else
(  687)       INTEGER ::           jk, jit   ! dummy loop indices
(  688) # endif
(  689) 
(  690)       NAMELIST/namtke/ ln_rstke, ediff, ediss, ebb, efave, emin, emin0,   &
(  691)          ri_c, nitke, nmxl, npdl, nave, navb
(  692)       !!----------------------------------------------------------------------
(  693) 
(  694)       ! Read Namelist namtke : Turbulente Kinetic Energy
(  695)       ! --------------------
(  696)       REWIND ( numnam )
(  697)       READ   ( numnam, namtke )
(  698) 
(  699)       ! Compute boost associated with the Richardson critic
(  700)       !     (control values: ri_c = 0.3   ==> eboost=1.25 for npdl=1 or 2)
(  701)       !     (                ri_c = 0.222 ==> eboost=1.                  )
(  702)       eboost = ri_c * ( 2. + ediss / ediff ) / 2.
(  703) 
(  704) 
(  705)       ! Parameter control and print
(  706)       ! ---------------------------
(  707)       ! Control print
(  708)       IF(lwp) THEN
(  709)          WRITE(numout,*)
(  710)          WRITE(numout,*) 'zdf_tke_init : tke turbulent closure scheme'
(  711)          WRITE(numout,*) '~~~~~~~~~~~~'
(  712)          WRITE(numout,*) '          Namelist namtke : set tke mixing parameters'
(  713)          WRITE(numout,*) '             restart with tke from no tke ln_rstke = ', ln_rstke
(  714)          WRITE(numout,*) '             coef. to compute avt           ediff  = ', ediff
(  715)          WRITE(numout,*) '             Kolmogoroff dissipation coef.  ediss  = ', ediss
(  716)          WRITE(numout,*) '             tke surface input coef.        ebb    = ', ebb
(  717)          WRITE(numout,*) '             tke diffusion coef.            efave  = ', efave
(  718)          WRITE(numout,*) '             minimum value of tke           emin   = ', emin
(  719)          WRITE(numout,*) '             surface minimum value of tke   emin0  = ', emin0
(  720)          WRITE(numout,*) '             number of restart iter loops   nitke  = ', nitke
(  721)          WRITE(numout,*) '             mixing length type             nmxl   = ', nmxl
(  722)          WRITE(numout,*) '             prandl number flag             npdl   = ', npdl
(  723)          WRITE(numout,*) '             horizontal average flag        nave   = ', nave
(  724)          WRITE(numout,*) '             critic Richardson nb           ri_c   = ', ri_c
(  725)          WRITE(numout,*) '                and its associated coeff.   eboost = ', eboost
(  726)          WRITE(numout,*) '             constant background or profile navb   = ', navb
(  727)          WRITE(numout,*)
(  728)       ENDIF
(  729) 
(  730)       ! Check nmxl and npdl values
(  731)       IF( nmxl < 0 .OR. nmxl > 3 ) THEN
(  732)          IF(lwp) WRITE(numout,cform_err)
(  733)          IF(lwp) WRITE(numout,*) '          bad flag: nmxl is < 0 or > 3 '
(  734)          nstop = nstop + 1
(  735)       ENDIF
(  736)       IF ( npdl < 0 .OR. npdl > 1 ) THEN
(  737)          IF(lwp) WRITE(numout,cform_err)
(  738)          IF(lwp) WRITE(numout,*) '          bad flag: npdl is < 0 or > 1 '
(  739)          nstop = nstop + 1
(  740)       ENDIF
(  741) 
(  742) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 14

(  743)       ! Horizontal average : initialization of weighting arrays 
(  744)       ! -------------------
(  745)       
(  746)       SELECT CASE ( nave )
(  747) 
(  748)       CASE ( 0 )                ! no horizontal average
(  749)          IF(lwp) WRITE(numout,*) '          no horizontal average on avt, avmu, avmv'
(  750)          IF(lwp) WRITE(numout,*) '          only in very high horizontal resolution !'
(  751) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  752) # if defined key_vectopt_memory
(  753)          ! weighting mean arrays etmean, eumean and evmean
(  754)          !           ( 1  1 )                                          ( 1 )
(  755)          ! avt = 1/4 ( 1  1 )     avmu = 1/2 ( 1  1 )       avmv=  1/2 ( 1 )
(  756)          !                          
(  757)          etmean(:,:,:) = 0.e0
(  758)          eumean(:,:,:) = 0.e0
(  759)          evmean(:,:,:) = 0.e0
(  760)          
(  761)          DO jk = 1, jpkm1
(  762)             DO jj = 2, jpjm1
(  763)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  764)                   etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
(  765)                   &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
(  766)                   &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
(  767)                   
(  768)                   eumean(ji,jj,jk) = umask(ji,jj,jk)                     &
(  769)                   &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk)  )
(  770) 
(  771)                   evmean(ji,jj,jk) = vmask(ji,jj,jk)                     &
(  772)                   &  / MAX( 1.,  tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk)  )
(  773)                END DO
(  774)             END DO
(  775)          END DO
(  776) # endif
(  777) 
(  778)       CASE ( 1 )                ! horizontal average 
(  779)          IF(lwp) WRITE(numout,*) '          horizontal average on avt, avmu, avmv'
(  780) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  781) # if defined key_vectopt_memory
(  782)          ! weighting mean arrays etmean, eumean and evmean
(  783)          !           ( 1  1 )              ( 1/2  1/2 )             ( 1/2  1  1/2 )
(  784)          ! avt = 1/4 ( 1  1 )   avmu = 1/4 ( 1    1   )   avmv= 1/4 ( 1/2  1  1/2 )
(  785)          !                                 ( 1/2  1/2 )
(  786)          etmean(:,:,:) = 0.e0
(  787)          eumean(:,:,:) = 0.e0
(  788)          evmean(:,:,:) = 0.e0
(  789)          
(  790)          DO jk = 1, jpkm1
(  791)             DO jj = 2, jpjm1
(  792)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  793)                   etmean(ji,jj,jk) = tmask(ji,jj,jk)                     &
(  794)                   &  / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
(  795)                   &            + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
(  796)                   
(  797)                   eumean(ji,jj,jk) = umask(ji,jj,jk)                        &
(  798)                   &  / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
(  799)                   &       +.5 * ( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &
(  800)                   &              +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 15

(  801) 
(  802)                   evmean(ji,jj,jk) = vmask(ji,jj,jk)                        &
(  803)                   &  / MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
(  804)                   &       +.5 * ( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
(  805)                   &              +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
(  806)                END DO
(  807)             END DO
(  808)          END DO
(  809) # endif
(  810) 
(  811)       CASE DEFAULT
(  812)          IF(lwp) WRITE(numout,cform_err)
(  813)          IF(lwp) WRITE(numout,*) '          bad flag value for nave = ', nave
(  814)          nstop = nstop + 1
(  815) 
(  816)       END SELECT
(  817) 
(  818) 
(  819)       ! Background eddy viscosity and diffusivity profil
(  820)       ! ------------------------------------------------
(  821)       IF( navb == 0 ) THEN
(  822)          !   Define avmb, avtb from namelist parameter
(  823)          avmb(:) = avm0
(  824)          avtb(:) = avt0
(  825)       ELSE
(  826)          !   Background profile of avt (fit a theoretical/observational profile (Krauss 1990) 
(  827)          avmb(:) = avm0
(  828)          avtb(:) = avt0 + ( 3.0e-4 - 2 * avt0 ) * 1.0e-4 * gdepw(:)   ! m2/s
(  829)       ENDIF
(  830) 
(  831) !!DB: delete ORCA
(  832) !      IF( cp_cfg == "orca"  .AND. jp_cfg == 2 .AND. ln_traadv_cen2 )  THEN
(  833) 
(  834)       ! Initialization of vertical eddy coef. to the background value
(  835)       ! -------------------------------------------------------------
(  836)       DO jk = 1, jpk
(  837)          avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
(  838)          avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
(  839)          avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
(  840)       END DO
(  841) 
(  842) 
(  843)       ! Initialization of turbulent kinetic energy ( en )
(  844)       ! -------------------------------------------------
(  845)       IF( ln_rstart ) THEN
(  846)          ! no en field in the restart file, en set by iterative loop
(  847)          IF( ln_rstke ) THEN
(  848)             en (:,:,:) = emin * tmask(:,:,:)
(  849)             DO jit = 2, nitke+1
(  850)                CALL zdf_tke( jit )
(  851)             END DO
(  852)          ENDIF
(  853)          ! otherwise en is already read in dtrlec called by inidtr
(  854)       ELSE
(  855)          ! no restart: en set to emin
(  856)          en(:,:,:) = emin * tmask(:,:,:)
(  857)       ENDIF
(  858) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 16

(  859)    END SUBROUTINE zdf_tke_init
(  860) 
(  861) #else
(  862)    !!----------------------------------------------------------------------
(  863)    !!   Dummy module :                                        NO TKE scheme
(  864)    !!----------------------------------------------------------------------
(  865)    USE in_out_manager
(  866)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .FALSE.   !: TKE flag
(  867) CONTAINS
(  868)    SUBROUTINE zdf_tke( kt )          ! Empty routine
(  869)       if(lwp) WRITE(numout,*) 'zdf_tke: You should not have seen this print! error?', kt
(  870)    END SUBROUTINE zdf_tke
(  871) #endif
(  872) 
(  873)    !!======================================================================
(  874) END MODULE zdftke
















































PGF90 (Version     10.5)          09/23/2010  09:49:31      page 17

(    1) # 1 "zdftke.F90"
(    1) MODULE zdftke
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdftke  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the tke 
(    5)    !!                 turbulent closure parameterization
(    6)    !!=====================================================================
(    7) # 8
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_zdftke'                                             TKE scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !!   zdf_tke      : update momentum and tracer Kz from a tke scheme
(   12)    !!   zdf_tke_init : initialization, namelist read, and parameters control
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and active tracers 
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE zdf_oce         ! ocean vertical physics
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE phycst          ! physical constants
(   21)    USE taumod          ! surface stress
(   22)    USE prtctl          ! Print control
(   23)    USE traadv_ctl      ! advection scheme control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC zdf_tke   ! routine called by step.F90
(   30) 
(   31)    !! * Share Module variables
(   32)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdftke = .TRUE.    !: TKE vertical mixing flag
(   33)    LOGICAL, PUBLIC ::         & !!: ** tke namelist (namtke) **
(   34)      ln_rstke = .FALSE.          !: =T restart with tke from a run without tke with 
(   35)      !                           !  a none zero initial value for en
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   37)       en                         !: now turbulent kinetic energy
(   38) 
(   39)    !! * Module variables
(   40)    INTEGER ::                 & !!! ** tke namelist (namtke) **
(   41)       nitke = 50 ,            &  ! number of restart iterative loops
(   42)       nmxl  =  2 ,            &  ! = 0/1/2/3 flag for the type of mixing length used
(   43)       npdl  =  1 ,            &  ! = 0/1/2 flag on prandtl number on vert. eddy coeff.
(   44)       nave  =  1 ,            &  ! = 0/1 flag for horizontal average on avt, avmu, avmv
(   45)       navb  =  0                 ! = 0/1 flag for constant or profile background avt
(   46)    REAL(wp) ::                & !!! ** tke namlist (namtke) **
(   47)       ediff = 0.1_wp       ,  &  ! coeff. for vertical eddy coef.; avt=ediff*mxl*sqrt(e)
(   48)       ediss = 0.7_wp       ,  &  ! coef. of the Kolmogoroff dissipation 
(   49)       ebb   = 3.75_wp      ,  &  ! coef. of the surface input of tke
(   50)       efave = 1._wp        ,  &  ! coef. for the tke vert. diff. coeff.; avtke=efave*avm
(   51)       emin  = 0.7071e-6_wp ,  &  ! minimum value of tke (m2/s2)
(   52)       emin0 = 1.e-4_wp     ,  &  ! surface minimum value of tke (m2/s2)
(   53)       ri_c  = 2._wp / 9._wp      ! critic Richardson number
(   54)    REAL(wp) ::   &  
(   55)       eboost        ! multiplicative coeff of the shear product.
(   56) 
(   57)    !! caution vectopt_memory change the solution (last digit of the solver stat)






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 18

(   58) 
(   59) 
(   60) # 73
(   73)    !! * Substitutions
(   74) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
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
(   48) # 75 "zdftke.F90"
(   75) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 76 "zdftke.F90"
(   76) # 76
(   76)    !!----------------------------------------------------------------------
(   77)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   78)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdftke.F90,v 1.9 2006/03/21 09:16:38 opalod Exp $ 
(   79)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   80)    !!----------------------------------------------------------------------
(   81) 
(   82) CONTAINS
(   83) 
(   84) # 91
(   91)    !!----------------------------------------------------------------------
(   92)    !!   default option :                                         k-j-i loop
(   93)    !!----------------------------------------------------------------------
(   94) 
(   95)    SUBROUTINE zdf_tke ( kt )
(   96)       !!----------------------------------------------------------------------
(   97)       !!                   ***  ROUTINE zdf_tke  ***
(   98)       !!
(   99)       !! ** Purpose :   Compute the vertical eddy viscosity and diffusivity
(  100)       !!      coefficients using a 1.5 turbulent closure scheme.
(  101)       !!
(  102)       !! ** Method  :   The time evolution of the turbulent kinetic energy 
(  103)       !!      (tke) is computed from a prognostic equation :
(  104)       !!         d(en)/dt = eboost eav (d(u)/dz)**2       ! shear production
(  105)       !!                  + d( efave eav d(en)/dz )/dz    ! diffusion of tke






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 19

(  106)       !!                  + grav/rau0 pdl eav d(rau)/dz   ! stratif. destruc.
(  107)       !!                  - ediss / emxl en**(2/3)        ! dissipation
(  108)       !!      with the boundary conditions:
(  109)       !!         surface: en = max( emin0,ebb sqrt(taux^2 + tauy^2) )
(  110)       !!         bottom : en = emin
(  111)       !!      -1- The dissipation and mixing turbulent lengh scales are computed
(  112)       !!      from the usual diagnostic buoyancy length scale:  
(  113)       !!         mxl= 1/(sqrt(en)/N)  WHERE N is the brunt-vaisala frequency
(  114)       !!      Four cases : 
(  115)       !!         nmxl=0 : mxl bounded by the distance to surface and bottom.
(  116)       !!                  zmxld = zmxlm = mxl
(  117)       !!         nmxl=1 : mxl bounded by the vertical scale factor.
(  118)       !!                  zmxld = zmxlm = mxl
(  119)       !!         nmxl=2 : mxl bounded such that the vertical derivative of mxl 
(  120)       !!                  is less than 1 (|d/dz(xml)|<1). 
(  121)       !!                  zmxld = zmxlm = mxl
(  122)       !!         nmxl=3 : lup = mxl bounded using |d/dz(xml)|<1 from the surface 
(  123)       !!                        to the bottom
(  124)       !!                  ldown = mxl bounded using |d/dz(xml)|<1 from the bottom 
(  125)       !!                        to the surface
(  126)       !!                  zmxld = sqrt (lup*ldown) ; zmxlm = min(lup,ldown)
(  127)       !!      -2- Compute the now Turbulent kinetic energy. The time differencing
(  128)       !!      is implicit for vertical diffusion term, linearized for kolmo-
(  129)       !!      goroff dissipation term, and explicit forward for both buoyancy
(  130)       !!      and dynamic production terms. Thus a tridiagonal linear system is
(  131)       !!      solved.
(  132)       !!         Note that - the shear production is multiplied by eboost in order
(  133)       !!      to set the critic richardson number to ri_c (namelist parameter)
(  134)       !!                   - the destruction by stratification term is multiplied
(  135)       !!      by the Prandtl number (defined by an empirical funtion of the local 
(  136)       !!      Richardson number) if npdl=1 (namelist parameter)
(  137)       !!      coefficient (zesh2):
(  138)       !!      -3- Compute the now vertical eddy vicosity and diffusivity
(  139)       !!      coefficients from en (before the time stepping) and zmxlm:
(  140)       !!              avm = max( avtb, ediff*zmxlm*en^1/2 )
(  141)       !!              avt = max( avmb, pdl*avm )  (pdl=1 if npdl=0)
(  142)       !!              eav = max( avmb, avm )
(  143)       !!      avt and avm are horizontally averaged to avoid numerical insta-
(  144)       !!      bilities.
(  145)       !!        N.B. The computation is done from jk=2 to jpkm1 except for
(  146)       !!      en. Surface value of avt avmu avmv are set once a time to
(  147)       !!      their background value in routine zdf_tke_init.
(  148)       !!
(  149)       !! ** Action  :   compute en (now turbulent kinetic energy)
(  150)       !!                update avt, avmu, avmv (before vertical eddy coef.)
(  151)       !!
(  152)       !! References :
(  153)       !!      Gaspar et al., jgr, 95, 1990,
(  154)       !!      Blanke and Delecluse, jpo, 1991
(  155)       !! History :
(  156)       !!   6.0  !  91-03  (b. blanke)  Original code
(  157)       !!   7.0  !  91-11  (G. Madec)   bug fix
(  158)       !!   7.1  !  92-10  (G. Madec)   new mixing length and eav
(  159)       !!   7.2  !  93-03  (M. Guyon)   symetrical conditions
(  160)       !!   7.3  !  94-08  (G. Madec, M. Imbard)   npdl flag
(  161)       !!   7.5  !  96-01  (G. Madec)   s-coordinates
(  162)       !!   8.0  !  97-07  (G. Madec)   lbc
(  163)       !!   8.1  !  99-01  (E. Stretta) new option for the mixing length






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 20

(  164)       !!   8.5  !  02-08  (G. Madec)  ri_c and Free form, F90
(  165)       !!   9.0  !  04-10  (C. Ethe )  1D configuration
(  166)       !!----------------------------------------------------------------------
(  167)       !! * Modules used
(  168)       USE oce     , zwd   => ua,  &  ! use ua as workspace
(  169)          &          zmxlm => ta,  &  ! use ta as workspace
(  170)          &          zmxld => sa      ! use sa as workspace
(  171) 
(  172)       !! * arguments
(  173)       INTEGER, INTENT( in  ) ::   kt ! ocean time step
(  174) 
(  175)       !! * local declarations
(  176)       INTEGER ::   ji, jj, jk        ! dummy loop arguments
(  177)       REAL(wp) ::   &
(  178)          zmlmin, zbbrau,          &  ! temporary scalars
(  179)          zfact1, zfact2, zfact3,  &  !
(  180)          zrn2, zesurf,            &  !
(  181)          ztx2, zty2, zav,         &  !
(  182)          zcoef, zcof, zsh2,       &  !
(  183)          zdku, zdkv, zpdl, zri,   &  !
(  184)          zsqen, zesh2,            &  !
(  185)          zemxl, zemlm, zemlp
(  186)       !!--------------------------------------------------------------------
(  187) 
(  188)       ! Initialization (first time-step only)
(  189)       ! --------------
(  190)       IF( kt == nit000  )   CALL zdf_tke_init
(  191) 
(  192)       ! Local constant initialization
(  193)       zmlmin = 1.e-8
(  194)       zbbrau =  .5 * ebb / rau0
(  195)       zfact1 = -.5 * rdt * efave
(  196)       zfact2 = 1.5 * rdt * ediss
(  197)       zfact3 = 0.5 * rdt * ediss
(  198) 
(  199) 
(  200)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  201)       ! I.  Mixing length
(  202)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  203) 
(  204) 
(  205)       ! Buoyancy length scale: l=sqrt(2*e/n**2)
(  206)       ! ---------------------
(  207)       zmxlm(:,:, 1 ) = zmlmin   ! surface set to the minimum value
(  208)       zmxlm(:,:,jpk) = zmlmin   ! bottom  set to the minimum value
(  209) !CDIR NOVERRCHK
(  210)       DO jk = 2, jpkm1
(  211) !CDIR NOVERRCHK
(  212)          DO jj = 2, jpjm1
(  213) !CDIR NOVERRCHK
(  214)             DO ji = 2, jpim1   ! vector opt.
(  215)                zrn2 = MAX( rn2(ji,jj,jk), rsmall )
(  216)                zmxlm(ji,jj,jk) = MAX( SQRT( 2. * en(ji,jj,jk) / zrn2 ), zmlmin  )
(  217)             END DO
(  218)          END DO
(  219)       END DO
(  220) 
(  221) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 21

(  222)       ! Physical limits for the mixing length
(  223)       ! -------------------------------------
(  224)       zmxld(:,:, 1 ) = zmlmin   ! surface set to the minimum value
(  225)       zmxld(:,:,jpk) = zmlmin   ! bottom  set to the minimum value
(  226) 
(  227)       SELECT CASE ( nmxl )
(  228) 
(  229)       CASE ( 0 )           ! bounded by the distance to surface and bottom
(  230)          DO jk = 2, jpkm1
(  231)             DO jj = 2, jpjm1
(  232)                DO ji = 2, jpim1   ! vector opt.
(  233)                   zemxl = MIN( gdepw(jk), zmxlm(ji,jj,jk),   &
(  234)                   &            gdepw(mbathy(ji,jj)) - gdepw(jk) )
(  235)                   zmxlm(ji,jj,jk) = zemxl
(  236)                   zmxld(ji,jj,jk) = zemxl
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240) 
(  241)       CASE ( 1 )           ! bounded by the vertical scale factor
(  242)          DO jk = 2, jpkm1
(  243)             DO jj = 2, jpjm1
(  244)                DO ji = 2, jpim1   ! vector opt.
(  245)                   zemxl = MIN( e3w(jk), zmxlm(ji,jj,jk) )
(  246)                   zmxlm(ji,jj,jk) = zemxl
(  247)                   zmxld(ji,jj,jk) = zemxl
(  248)                END DO
(  249)             END DO
(  250)          END DO
(  251) 
(  252)       CASE ( 2 )           ! |dk[xml]| bounded by e3t :
(  253)          DO jk = 2, jpkm1         ! from the surface to the bottom :
(  254)             DO jj = 2, jpjm1
(  255)                DO ji = 2, jpim1   ! vector opt.
(  256)                   zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk-1) + e3t(jk-1), zmxlm(ji,jj,jk) )
(  257)                END DO
(  258)             END DO
(  259)          END DO
(  260)          DO jk = jpkm1, 2, -1     ! from the bottom to the surface :
(  261)             DO jj = 2, jpjm1
(  262)                DO ji = 2, jpim1   ! vector opt.
(  263)                   zemxl = MIN( zmxlm(ji,jj,jk+1) + e3t(jk+1), zmxlm(ji,jj,jk) )
(  264)                   zmxlm(ji,jj,jk) = zemxl
(  265)                   zmxld(ji,jj,jk) = zemxl
(  266)                END DO
(  267)             END DO
(  268)          END DO
(  269) 
(  270)       CASE ( 3 )           ! lup and ldown, |dk[xml]| bounded by e3t :
(  271)          DO jk = 2, jpkm1         ! from the surface to the bottom : lup
(  272)             DO jj = 2, jpjm1
(  273)                DO ji = 2, jpim1   ! vector opt.
(  274)                   zmxld(ji,jj,jk) = MIN( zmxld(ji,jj,jk-1) + e3t(jk-1), zmxlm(ji,jj,jk) )
(  275)                END DO
(  276)             END DO
(  277)          END DO
(  278)          DO jk = jpkm1, 2, -1     ! from the bottom to the surface : ldown
(  279)             DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 22

(  280)                DO ji = 2, jpim1   ! vector opt.
(  281)                   zmxlm(ji,jj,jk) = MIN( zmxlm(ji,jj,jk+1) + e3t(jk+1), zmxlm(ji,jj,jk) )
(  282)                END DO
(  283)             END DO
(  284)          END DO
(  285) !CDIR NOVERRCHK
(  286)          DO jk = 2, jpkm1
(  287) !CDIR NOVERRCHK
(  288)             DO jj = 2, jpjm1
(  289) !CDIR NOVERRCHK
(  290)                DO ji = 2, jpim1   ! vector opt.
(  291)                   zemlm = MIN ( zmxld(ji,jj,jk),  zmxlm(ji,jj,jk) )
(  292)                   zemlp = SQRT( zmxld(ji,jj,jk) * zmxlm(ji,jj,jk) )
(  293)                   zmxlm(ji,jj,jk) = zemlm
(  294)                   zmxld(ji,jj,jk) = zemlp
(  295)                END DO
(  296)             END DO
(  297)          END DO
(  298) 
(  299)       END SELECT
(  300) 
(  301) 
(  302) 
(  303) # 308
(  308)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  309)       ! II  Tubulent kinetic energy time stepping
(  310)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  311) 
(  312) 
(  313)       ! 1. Vertical eddy viscosity on tke (put in zmxlm) and first estimate of avt
(  314)       ! ---------------------------------------------------------------------
(  315) !CDIR NOVERRCHK
(  316)       DO jk = 2, jpkm1
(  317) !CDIR NOVERRCHK
(  318)          DO jj = 2, jpjm1
(  319) !CDIR NOVERRCHK
(  320)             DO ji = 2, jpim1   ! vector opt.
(  321)                zsqen = SQRT( en(ji,jj,jk) )
(  322)                zav   = ediff * zmxlm(ji,jj,jk) * zsqen
(  323)                avt  (ji,jj,jk) = MAX( zav, avtb(jk) ) * tmask(ji,jj,jk)
(  324)                zmxlm(ji,jj,jk) = MAX( zav, avmb(jk) ) * tmask(ji,jj,jk)
(  325)                zmxld(ji,jj,jk) = zsqen / zmxld(ji,jj,jk)
(  326)             END DO
(  327)          END DO
(  328)       END DO
(  329) 
(  330)       ! 2. Surface boundary condition on tke and its eddy viscosity (zmxlm)
(  331)       ! -------------------------------------------------
(  332)       ! en(1)   = ebb sqrt(taux^2+tauy^2) / rau0  (min value emin0)
(  333)       ! zmxlm(1) = avmb(1) and zmxlm(jpk) = 0.
(  334) !CDIR NOVERRCHK
(  335)       DO jj = 2, jpjm1
(  336) !CDIR NOVERRCHK
(  337)          DO ji = 2, jpim1   ! vector opt.
(  338)             ztx2 = taux(ji-1,jj  ) + taux(ji,jj)
(  339)             zty2 = tauy(ji  ,jj-1) + tauy(ji,jj)
(  340)             zesurf = zbbrau * SQRT( ztx2 * ztx2 + zty2 * zty2 )
(  341)             en (ji,jj,1) = MAX( zesurf, emin0 ) * tmask(ji,jj,1)






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 23

(  342)             zmxlm(ji,jj,1  ) = avmb(1) * tmask(ji,jj,1)
(  343)             zmxlm(ji,jj,jpk) = 0.e0
(  344)          END DO
(  345)       END DO
(  346) 
(  347)       ! 3. Now Turbulent kinetic energy (output in en)
(  348)       ! -------------------------------
(  349)       ! Resolution of a tridiagonal linear system by a "methode de chasse"
(  350)       ! computation from level 2 to jpkm1  (e(1) already computed and
(  351)       ! e(jpk)=0 ). 
(  352) 
(  353)       SELECT CASE ( npdl )
(  354) 
(  355)       CASE ( 0 )           ! No Prandtl number
(  356)          DO jk = 2, jpkm1
(  357)             DO jj = 2, jpjm1
(  358)                DO ji = 2, jpim1   ! vector opt.
(  359)                   ! zesh2 = eboost * (du/dz)^2 - N^2
(  360)                   zcoef = 0.5 / e3w(jk)
(  361)                   ! shear
(  362)                   zdku = zcoef * (   ub(ji-1, jj ,jk-1) + ub(ji,jj,jk-1)   &
(  363)                   &                - ub(ji-1, jj ,jk  ) - ub(ji,jj,jk  )  )
(  364)                   zdkv = zcoef * (   vb( ji ,jj-1,jk-1) + vb(ji,jj,jk-1)   &
(  365)                   &                - vb( ji ,jj-1,jk  ) - vb(ji,jj,jk  )  )
(  366)                   ! coefficient (zesh2)
(  367)                   zesh2 =  eboost * ( zdku*zdku + zdkv*zdkv ) - rn2(ji,jj,jk)
(  368) 
(  369)                   ! Matrix
(  370)                   zcof = zfact1 * tmask(ji,jj,jk)
(  371)                   ! lower diagonal
(  372)                   avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
(  373)                   &                    / ( e3t(jk-1) * e3w(jk  ) )
(  374)                   ! upper diagonal
(  375)                   avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
(  376)                   &                    / ( e3t(jk  ) * e3w(jk) )
(  377)                   ! diagonal
(  378)                   zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
(  379)                   ! right hand side in en 
(  380)                   en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
(  381)                   &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
(  382)                END DO
(  383)             END DO
(  384)          END DO
(  385) 
(  386)       CASE ( 1 )           ! Prandtl number
(  387)          DO jk = 2, jpkm1
(  388)             DO jj = 2, jpjm1
(  389)                DO ji = 2, jpim1   ! vector opt.
(  390)                   ! zesh2 =  eboost * (du/dz)^2 - pdl * N^2
(  391)                   zcoef = 0.5 / e3w(jk)
(  392)                   ! shear
(  393)                   zdku = zcoef * (   ub(ji-1,jj  ,jk-1) + ub(ji,jj,jk-1) &
(  394)                   &                - ub(ji-1,jj  ,jk  ) - ub(ji,jj,jk  )   )
(  395)                   zdkv = zcoef * (   vb(ji  ,jj-1,jk-1) + vb(ji,jj,jk-1) &
(  396)                   &                - vb(ji  ,jj-1,jk  ) - vb(ji,jj,jk  )   )
(  397)                   ! square of vertical shear
(  398)                   zsh2 = zdku * zdku + zdkv * zdkv
(  399)                   ! local Richardson number






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 24

(  400)                   zri  = MAX( rn2(ji,jj,jk), 0. ) / ( zsh2 + 1.e-20 )
(  401) # 405
(  405)                   ! Prandtl number
(  406)                   zpdl = 1.0
(  407)                   IF( zri >= 0.2 ) zpdl = 0.2 / zri
(  408)                   zpdl = MAX( 0.1, zpdl )
(  409)                   ! coefficient (esh2)
(  410)                   zesh2 = eboost * zsh2 - zpdl * rn2(ji,jj,jk)
(  411) 
(  412)                   ! Matrix
(  413)                   zcof = zfact1 * tmask(ji,jj,jk)
(  414)                   ! lower diagonal
(  415)                   avmv(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk  ) + zmxlm(ji,jj,jk-1) )   &
(  416)                   &                     / ( e3t(jk-1) * e3w(jk  ) )
(  417)                   ! upper diagonal
(  418)                   avmu(ji,jj,jk) = zcof * ( zmxlm(ji,jj,jk+1) + zmxlm(ji,jj,jk  ) )   &
(  419)                   &                     / ( e3t(jk  ) * e3w(jk) )
(  420)                   ! diagonal
(  421)                   zwd(ji,jj,jk) = 1. - avmv(ji,jj,jk) - avmu(ji,jj,jk) + zfact2 * zmxld(ji,jj,jk)
(  422)                   ! right hand side in en 
(  423)                   en(ji,jj,jk) = en(ji,jj,jk) + zfact3 * zmxld(ji,jj,jk) * en   (ji,jj,jk)   &
(  424)                   &                           +   rdt  * zmxlm(ji,jj,jk) * zesh2
(  425)                   ! save masked Prandlt number in zmxlm array
(  426)                   zmxld(ji,jj,jk) = zpdl * tmask(ji,jj,jk)
(  427)                END DO
(  428)             END DO
(  429)          END DO
(  430) 
(  431)       END SELECT
(  432) 
(  433) 
(  434) # 440
(  440)       ! 4. Matrix inversion from level 2 (tke prescribed at level 1)
(  441)       !!--------------------------------
(  442)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1
(  443)       DO jk = 3, jpkm1
(  444)          DO jj = 2, jpjm1
(  445)             DO ji = 2, jpim1    ! vector opt.
(  446)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - avmv(ji,jj,jk) * avmu(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  447)             END DO
(  448)          END DO
(  449)       END DO
(  450) 
(  451)       ! Second recurrence : Lk = RHSk - Lk / Dk-1 * Lk-1
(  452)       DO jj = 2, jpjm1
(  453)          DO ji = 2, jpim1    ! vector opt.
(  454)             avmv(ji,jj,2) = en(ji,jj,2) - avmv(ji,jj,2) * en(ji,jj,1)    ! Surface boudary conditions on tke
(  455)          END DO
(  456)       END DO
(  457)       DO jk = 3, jpkm1
(  458)          DO jj = 2, jpjm1
(  459)             DO ji = 2, jpim1    ! vector opt.
(  460)                avmv(ji,jj,jk) = en(ji,jj,jk) - avmv(ji,jj,jk) / zwd(ji,jj,jk-1) *avmv(ji,jj,jk-1)
(  461)             END DO
(  462)          END DO
(  463)       END DO
(  464) 
(  465)       ! thrid recurrence : Ek = ( Lk - Uk * Ek+1 ) / Dk






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 25

(  466)       DO jj = 2, jpjm1
(  467)          DO ji = 2, jpim1    ! vector opt.
(  468)             en(ji,jj,jpkm1) = avmv(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  469)          END DO
(  470)       END DO
(  471)       DO jk = jpk-2, 2, -1
(  472)          DO jj = 2, jpjm1
(  473)             DO ji = 2, jpim1    ! vector opt.
(  474)                en(ji,jj,jk) = ( avmv(ji,jj,jk) - avmu(ji,jj,jk) * en(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  475)             END DO
(  476)          END DO
(  477)       END DO
(  478) 
(  479)       ! Save the result in en and set minimum value of tke : emin
(  480)       DO jk = 2, jpkm1
(  481)          DO jj = 2, jpjm1
(  482)             DO ji = 2, jpim1   ! vector opt.
(  483)                en(ji,jj,jk) = MAX( en(ji,jj,jk), emin ) * tmask(ji,jj,jk)
(  484)             END DO
(  485)          END DO
(  486)       END DO
(  487) 
(  488)       ! Lateral boundary conditions on ( avt, en )   (sign unchanged)
(  489)       CALL lbc_lnk( en , 'W', 1. )   ;   CALL lbc_lnk( avt, 'W', 1. )
(  490) 
(  491) 
(  492)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  493)       ! III.  Before vertical eddy vicosity and diffusivity coefficients
(  494)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  495) 
(  496)       SELECT CASE ( nave )
(  497)          
(  498)       CASE ( 0 )                ! no horizontal average
(  499) 
(  500)          ! Vertical eddy viscosity
(  501) 
(  502)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  503)             DO jj = 2, jpjm1
(  504)                DO ji = 2, jpim1   ! vector opt.
(  505)                   avmu(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji+1,jj  ,jk) ) * umask(ji,jj,jk)   &
(  506)                   &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji+1,jj  ,jk) )
(  507)                   avmv(ji,jj,jk) = ( avt  (ji,jj,jk) + avt  (ji  ,jj+1,jk) ) * vmask(ji,jj,jk)   &
(  508)                   &     / MAX( 1.,   tmask(ji,jj,jk) + tmask(ji  ,jj+1,jk) )
(  509)                END DO
(  510)             END DO
(  511)          END DO
(  512) 
(  513) !DB: jack up viscosity in shallow water 
(  514) !         do jj = 2, jpjm1
(  515) !            do ji = 2, jpim1   ! vector opt.
(  516) !               if(mbathy(ji,jj) <= 5) then
(  517) !                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
(  518) !                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
(  519) !               endif
(  520) !            enddo
(  521) !         enddo
(  522) !!DB: END
(  523) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 26

(  524)          ! Lateral boundary conditions (avmu,avmv) (U- and V- points, sign unchanged)
(  525)          CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )
(  526)          
(  527)       CASE ( 1 )                ! horizontal average
(  528) 
(  529)          !                                                ( 1/2  1/2 )
(  530)          ! Eddy viscosity: horizontal average: avmu = 1/4 ( 1    1   )
(  531)          !                      ( 1/2  1 1/2 )            ( 1/2  1/2 )
(  532)          !           avmv = 1/4 ( 1/2  1 1/2 )      
(  533)          
(  534) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  535) # 561
(  561)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  562)             DO jj = 2, jpjm1
(  563)                DO ji = 2, jpim1   ! vector opt.
(  564)                   avmu(ji,jj,jk) = (   avt  (ji,jj  ,jk) + avt  (ji+1,jj  ,jk)   &
(  565)                   &              +.5*( avt  (ji,jj-1,jk) + avt  (ji+1,jj-1,jk)   &
(  566)                   &                   +avt  (ji,jj+1,jk) + avt  (ji+1,jj+1,jk) ) ) * umask(ji,jj,jk)  &
(  567)                   &       / MAX( 1.,   tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk)   &
(  568)                   &              +.5*( tmask(ji,jj-1,jk) + tmask(ji+1,jj-1,jk)   &
(  569)                   &                   +tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) )  )
(  570) 
(  571)                   avmv(ji,jj,jk) = (   avt  (ji  ,jj,jk) + avt  (ji  ,jj+1,jk)   &
(  572)                   &              +.5*( avt  (ji-1,jj,jk) + avt  (ji-1,jj+1,jk)   &
(  573)                   &                   +avt  (ji+1,jj,jk) + avt  (ji+1,jj+1,jk) ) ) * vmask(ji,jj,jk)  &
(  574)                   &      /  MAX( 1.,   tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
(  575)                   &              +.5*( tmask(ji-1,jj,jk) + tmask(ji-1,jj+1,jk)   &
(  576)                   &                   +tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk) )  )
(  577)                END DO
(  578)             END DO
(  579)          END DO
(  580) 
(  581) !DB: jack up viscosity in shallow water
(  582) !         do jj = 2, jpjm1
(  583) !            do ji = 2, jpim1   ! vector opt.
(  584) !               if(mbathy(ji,jj) <= 5) then
(  585) !                  avmu(ji,jj,:) = 10. * avmu(ji,jj,:)
(  586) !                  avmv(ji,jj,:) = 10. * avmv(ji,jj,:)
(  587) !               endif
(  588) !            enddo
(  589) !         enddo
(  590) !!DB: END
(  591) 
(  592) 
(  593) # 594
(  594)          ! Lateral boundary conditions (avmu,avmv) (sign unchanged)
(  595)          CALL lbc_lnk( avmu, 'U', 1. )   ;    CALL lbc_lnk( avmv, 'V', 1. )
(  596) 
(  597)          ! Vertical eddy diffusivity
(  598)          ! ------------------------------
(  599)          !                                (1 2 1)
(  600)          ! horizontal average  avt = 1/16 (2 4 2)
(  601)          !                                (1 2 1)
(  602)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  603) # 611
(  611)             DO jj = 2, jpjm1
(  612)                DO ji = 2, jpim1   ! vector opt.
(  613)                   avt(ji,jj,jk) = ( avmu (ji,jj,jk) + avmu (ji-1,jj  ,jk)   &






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 27

(  614)                   &               + avmv (ji,jj,jk) + avmv (ji  ,jj-1,jk)  ) * tmask(ji,jj,jk)   &
(  615)                   &     / MAX( 1.,  umask(ji,jj,jk) + umask(ji-1,jj  ,jk)   &
(  616)                   &               + vmask(ji,jj,jk) + vmask(ji  ,jj-1,jk)  )
(  617)                END DO
(  618)             END DO
(  619) # 620
(  620)          END DO
(  621) 
(  622)       END SELECT
(  623) 
(  624)       ! multiplied by the Prandtl number (npdl>1)
(  625)       ! ----------------------------------------
(  626)       IF( npdl == 1 ) THEN
(  627)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  628)             DO jj = 2, jpjm1
(  629)                DO ji = 2, jpim1   ! vector opt.
(  630)                   zpdl = zmxld(ji,jj,jk)
(  631)                   avt(ji,jj,jk) = MAX( zpdl * avt(ji,jj,jk), avtb(jk) ) * tmask(ji,jj,jk)
(  632)                END DO
(  633)             END DO
(  634)          END DO
(  635)       ENDIF
(  636) 
(  637)       ! Minimum value on the eddy viscosity
(  638)       ! ----------------------------------------
(  639)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  640)          DO jj = 1, jpj
(  641)             DO ji = 1, jpi
(  642)                avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk), avmb(jk) ) * umask(ji,jj,jk)
(  643)                avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk), avmb(jk) ) * vmask(ji,jj,jk)
(  644)             END DO
(  645)          END DO
(  646)       END DO
(  647) 
(  648)       ! Lateral boundary conditions on avt  (sign unchanged)
(  649)       ! ------------------------------=====
(  650)       CALL lbc_lnk( avt, 'W', 1. )
(  651) 
(  652)       IF(ln_ctl) THEN
(  653)          CALL prt_ctl(tab3d_1=en  , clinfo1=' tke  - e: ', tab3d_2=avt , clinfo2=' t: ', ovlap=1, kdim=jpk)
(  654)          CALL prt_ctl(tab3d_1=avmu, clinfo1=' tke  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
(  655)       ENDIF
(  656) 
(  657)    END SUBROUTINE zdf_tke
(  658) 
(  659) 
(  660) # 661
(  661)    SUBROUTINE zdf_tke_init
(  662)       !!----------------------------------------------------------------------
(  663)       !!                  ***  ROUTINE zdf_tke_init  ***
(  664)       !!                     
(  665)       !! ** Purpose :   Initialization of the vertical eddy diffivity and 
(  666)       !!      viscosity when using a tke turbulent closure scheme
(  667)       !!
(  668)       !! ** Method  :   Read the namtke namelist and check the parameters
(  669)       !!      called at the first timestep (nit000)
(  670)       !!
(  671)       !! ** input   :   Namlist namtke






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 28

(  672)       !!
(  673)       !! ** Action  :   Increase by 1 the nstop flag is setting problem encounter
(  674)       !!
(  675)       !! history :
(  676)       !!  8.5  ! 02-06 (G. Madec) original code
(  677)       !!----------------------------------------------------------------------
(  678)       !! * Module used
(  679)       USE dynzdf_exp
(  680)       USE trazdf_exp
(  681) 
(  682)       !! * local declarations
(  683)       !! caution vectopt_memory change the solution (last digit of the solver stat)
(  684) # 687
(  687)       INTEGER ::           jk, jit   ! dummy loop indices
(  688) 
(  689) # 690
(  690)       NAMELIST/namtke/ ln_rstke, ediff, ediss, ebb, efave, emin, emin0,   &
(  691)          ri_c, nitke, nmxl, npdl, nave, navb
(  692)       !!----------------------------------------------------------------------
(  693) 
(  694)       ! Read Namelist namtke : Turbulente Kinetic Energy
(  695)       ! --------------------
(  696)       REWIND ( numnam )
(  697)       READ   ( numnam, namtke )
(  698) 
(  699)       ! Compute boost associated with the Richardson critic
(  700)       !     (control values: ri_c = 0.3   ==> eboost=1.25 for npdl=1 or 2)
(  701)       !     (                ri_c = 0.222 ==> eboost=1.                  )
(  702)       eboost = ri_c * ( 2. + ediss / ediff ) / 2.
(  703) 
(  704) 
(  705)       ! Parameter control and print
(  706)       ! ---------------------------
(  707)       ! Control print
(  708)       IF(lwp) THEN
(  709)          WRITE(numout,*)
(  710)          WRITE(numout,*) 'zdf_tke_init : tke turbulent closure scheme'
(  711)          WRITE(numout,*) '~~~~~~~~~~~~'
(  712)          WRITE(numout,*) '          Namelist namtke : set tke mixing parameters'
(  713)          WRITE(numout,*) '             restart with tke from no tke ln_rstke = ', ln_rstke
(  714)          WRITE(numout,*) '             coef. to compute avt           ediff  = ', ediff
(  715)          WRITE(numout,*) '             Kolmogoroff dissipation coef.  ediss  = ', ediss
(  716)          WRITE(numout,*) '             tke surface input coef.        ebb    = ', ebb
(  717)          WRITE(numout,*) '             tke diffusion coef.            efave  = ', efave
(  718)          WRITE(numout,*) '             minimum value of tke           emin   = ', emin
(  719)          WRITE(numout,*) '             surface minimum value of tke   emin0  = ', emin0
(  720)          WRITE(numout,*) '             number of restart iter loops   nitke  = ', nitke
(  721)          WRITE(numout,*) '             mixing length type             nmxl   = ', nmxl
(  722)          WRITE(numout,*) '             prandl number flag             npdl   = ', npdl
(  723)          WRITE(numout,*) '             horizontal average flag        nave   = ', nave
(  724)          WRITE(numout,*) '             critic Richardson nb           ri_c   = ', ri_c
(  725)          WRITE(numout,*) '                and its associated coeff.   eboost = ', eboost
(  726)          WRITE(numout,*) '             constant background or profile navb   = ', navb
(  727)          WRITE(numout,*)
(  728)       ENDIF
(  729) 
(  730)       ! Check nmxl and npdl values
(  731)       IF( nmxl < 0 .OR. nmxl > 3 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 29

(  732)          IF(lwp) WRITE(numout,cform_err)
(  733)          IF(lwp) WRITE(numout,*) '          bad flag: nmxl is < 0 or > 3 '
(  734)          nstop = nstop + 1
(  735)       ENDIF
(  736)       IF ( npdl < 0 .OR. npdl > 1 ) THEN
(  737)          IF(lwp) WRITE(numout,cform_err)
(  738)          IF(lwp) WRITE(numout,*) '          bad flag: npdl is < 0 or > 1 '
(  739)          nstop = nstop + 1
(  740)       ENDIF
(  741) 
(  742) 
(  743)       ! Horizontal average : initialization of weighting arrays 
(  744)       ! -------------------
(  745)       
(  746)       SELECT CASE ( nave )
(  747) 
(  748)       CASE ( 0 )                ! no horizontal average
(  749)          IF(lwp) WRITE(numout,*) '          no horizontal average on avt, avmu, avmv'
(  750)          IF(lwp) WRITE(numout,*) '          only in very high horizontal resolution !'
(  751) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  752) 
(  753) # 778
(  778)       CASE ( 1 )                ! horizontal average 
(  779)          IF(lwp) WRITE(numout,*) '          horizontal average on avt, avmu, avmv'
(  780) !! caution vectopt_memory change the solution (last digit of the solver stat)
(  781) 
(  782) # 811
(  811)       CASE DEFAULT
(  812)          IF(lwp) WRITE(numout,cform_err)
(  813)          IF(lwp) WRITE(numout,*) '          bad flag value for nave = ', nave
(  814)          nstop = nstop + 1
(  815) 
(  816)       END SELECT
(  817) 
(  818) 
(  819)       ! Background eddy viscosity and diffusivity profil
(  820)       ! ------------------------------------------------
(  821)       IF( navb == 0 ) THEN
(  822)          !   Define avmb, avtb from namelist parameter
(  823)          avmb(:) = avm0
(  824)          avtb(:) = avt0
(  825)       ELSE
(  826)          !   Background profile of avt (fit a theoretical/observational profile (Krauss 1990) 
(  827)          avmb(:) = avm0
(  828)          avtb(:) = avt0 + ( 3.0e-4 - 2 * avt0 ) * 1.0e-4 * gdepw(:)   ! m2/s
(  829)       ENDIF
(  830) 
(  831) !!DB: delete ORCA
(  832) !      IF( cp_cfg == "orca"  .AND. jp_cfg == 2 .AND. ln_traadv_cen2 )  THEN
(  833) 
(  834)       ! Initialization of vertical eddy coef. to the background value
(  835)       ! -------------------------------------------------------------
(  836)       DO jk = 1, jpk
(  837)          avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
(  838)          avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
(  839)          avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
(  840)       END DO
(  841) 






PGF90 (Version     10.5)          09/23/2010  09:49:31      page 30

(  842) 
(  843)       ! Initialization of turbulent kinetic energy ( en )
(  844)       ! -------------------------------------------------
(  845)       IF( ln_rstart ) THEN
(  846)          ! no en field in the restart file, en set by iterative loop
(  847)          IF( ln_rstke ) THEN
(  848)             en (:,:,:) = emin * tmask(:,:,:)
(  849)             DO jit = 2, nitke+1
(  850)                CALL zdf_tke( jit )
(  851)             END DO
(  852)          ENDIF
(  853)          ! otherwise en is already read in dtrlec called by inidtr
(  854)       ELSE
(  855)          ! no restart: en set to emin
(  856)          en(:,:,:) = emin * tmask(:,:,:)
(  857)       ENDIF
(  858) 
(  859)    END SUBROUTINE zdf_tke_init
(  860) 
(  861) 
(  862) # 873
(  873)    !!======================================================================
(  874) END MODULE zdftke
