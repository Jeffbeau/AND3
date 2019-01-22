


PGF90 (Version     10.5)          09/23/2010  09:49:09      page 1

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

Filename: trabbl.F90

(    1) MODULE trabbl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trabbl  ***
(    4)    !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
(    5)    !!==============================================================================
(    6) #if   defined key_trabbl_dif   ||   defined key_trabbl_adv   || defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_trabbl_dif'   or            diffusive bottom boundary layer
(    9)    !!   'key_trabbl_adv'                 advective bottom boundary layer
(   10)    !!----------------------------------------------------------------------
(   11)    !!   tra_bbl_dif  : update the active tracer trends due to the bottom
(   12)    !!                  boundary layer (diffusive only)
(   13)    !!   tra_bbl_adv  : update the active tracer trends due to the bottom
(   14)    !!                  boundary layer (advective and/or diffusive)
(   15)    !!   tra_bbl_init : initialization, namlist read, parameters control
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce                  ! ocean dynamics and active tracers
(   19)    USE dom_oce              ! ocean space and time domain
(   20)    USE trdmod_oce           ! ocean variables trends
(   21)    USE in_out_manager       ! I/O manager
(   22)    USE prtctl               ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC tra_bbl_dif    ! routine called by step.F90
(   29)    PUBLIC tra_bbl_adv    ! routine called by step.F90
(   30) 
(   31)    !! * Shared module variables
(   32)    REAL(wp), PUBLIC ::            &  !!: * bbl namelist *
(   33)       atrbbl = 1.e+3                  !: lateral coeff. for BBL scheme (m2/s)
(   34) #if defined key_trabbl_dif
(   35)    LOGICAL, PUBLIC, PARAMETER ::   &  !:
(   36)       lk_trabbl_dif = .TRUE.          !: diffusive bottom boundary layer flag
(   37) #else
(   38)    LOGICAL, PUBLIC, PARAMETER ::   &  !:
(   39)       lk_trabbl_dif = .FALSE.         !: diffusive bottom boundary layer flag
(   40) #endif
(   41) 
(   42) # if defined key_trabbl_adv
(   43)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   44)       lk_trabbl_adv = .TRUE.   !: advective bottom boundary layer flag
(   45)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   46)        u_bbl, v_bbl,  &  !: velocity involved in exhanges in the advective BBL






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 2

(   47)        w_bbl             !: vertical increment of velocity due to advective BBL
(   48)        !                 !  only affect tracer vertical advection
(   49) # else
(   50)    LOGICAL, PUBLIC, PARAMETER ::    &  !:
(   51)       lk_trabbl_adv = .FALSE.  !: advective bottom boundary layer flag
(   52) # endif
(   53) 
(   54)    !! * Module variables
(   55)    INTEGER, DIMENSION(jpi,jpj) ::   &  !:
(   56)       mbkt,           &   ! vertical index of the bottom ocean T-level
(   57)       mbku, mbkv          ! vertical index of the bottom ocean U/V-level
(   58) 
(   59)    !! * Substitutions
(   60) #  include "domzgr_substitute.h90"
(   61) #  include "vectopt_loop_substitute.h90"
(   62)    !!----------------------------------------------------------------------
(   63)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   64)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trabbl.F90,v 1.10 2006/03/27 09:55:49 opalod Exp $ 
(   65)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68) CONTAINS
(   69) 
(   70)    SUBROUTINE tra_bbl_dif( kt )
(   71)       !!----------------------------------------------------------------------
(   72)       !!                   ***  ROUTINE tra_bbl_dif  ***
(   73)       !!
(   74)       !! ** Purpose :   Compute the before tracer (t & s) trend associated 
(   75)       !!      with the bottom boundary layer and add it to the general trend 
(   76)       !!      of tracer equations. The bottom boundary layer is supposed to be
(   77)       !!      a purely diffusive bottom boundary layer.
(   78)       !!
(   79)       !! ** Method  :   When the product grad( rho) * grad(h) < 0 (where grad 
(   80)       !!      is an along bottom slope gradient) an additional lateral diffu-
(   81)       !!      sive trend along the bottom slope is added to the general tracer
(   82)       !!      trend, otherwise nothing is done.
(   83)       !!      Second order operator (laplacian type) with variable coefficient
(   84)       !!      computed as follow for temperature (idem on s): 
(   85)       !!         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
(   86)       !!                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
(   87)       !!      where ztb is a 2D array: the bottom ocean temperature and ahtb
(   88)       !!      is a time and space varying diffusive coefficient defined by:
(   89)       !!         ahbt = zahbp    if grad(rho).grad(h) < 0
(   90)       !!              = 0.       otherwise.
(   91)       !!      Note that grad(.) is the along bottom slope gradient. grad(rho)
(   92)       !!      is evaluated using the local density (i.e. referenced at the
(   93)       !!      local depth). Typical value of ahbt is 2000 m2/s (equivalent to
(   94)       !!      a downslope velocity of 20 cm/s if the condition for slope
(   95)       !!      convection is satified)
(   96)       !!      Add this before trend to the general trend (ta,sa) of the 
(   97)       !!      botton ocean tracer point:
(   98)       !!         ta = ta + difft
(   99)       !!
(  100)       !! ** Action  : - update (ta,sa) at the bottom level with the bottom
(  101)       !!                boundary layer trend
(  102)       !!              - save the trends in tldfbbl/sldfbbl ('key_trdtra')
(  103)       !!
(  104)       !! References :






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 3

(  105)       !!     Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
(  106)       !!
(  107)       !! History :
(  108)       !!   8.0  !  96-06  (L. Mortier)  Original code
(  109)       !!   8.0  !  97-11  (G. Madec)  Optimization
(  110)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(  111)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  112)       !!----------------------------------------------------------------------
(  113)       !! * Modules used     
(  114)       USE oce, ONLY :    ztdta => ua,     & ! use ua as 3D workspace   
(  115)                          ztdsa => va        ! use va as 3D workspace   
(  116)       USE eosbn2, ONLY : neos ! type of equation of state
(  117) 
(  118)       !! * Arguments 
(  119)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(  120) 
(  121)       !! * Local declarations
(  122)       INTEGER ::   ji, jj                   ! dummy loop indices
(  123)       INTEGER ::   ik
(  124)       INTEGER ::   ii0, ii1, ij0, ij1       ! temporary integers
(  125) #  if defined key_partial_steps
(  126)       INTEGER  ::   iku1, iku2, ikv1,ikv2   ! temporary intergers
(  127)       REAL(wp) ::   ze3u, ze3v              ! temporary scalars
(  128) #  else
(  129)       INTEGER ::   iku, ikv
(  130) #  endif
(  131)       REAL(wp) ::   &
(  132)          zsign, zt, zs, zh, zalbet,      &  ! temporary scalars
(  133)          zgdrho, zbtr, zta, zsa
(  134)       REAL(wp), DIMENSION(jpi,jpj) ::    &
(  135)         zki, zkj, zkw, zkx, zky, zkz,    &  ! temporary workspace arrays
(  136)         ztnb, zsnb, zdep,                &
(  137)         ztbb, zsbb, zahu, zahv
(  138)       REAL(wp) ::   &
(  139)          fsalbt, pft, pfs, pfh              ! statement function
(  140)       !!----------------------------------------------------------------------
(  141)       ! ratio alpha/beta
(  142)       ! ================
(  143)       !  fsalbt: ratio of thermal over saline expension coefficients
(  144)       !       pft :  potential temperature in degrees celcius
(  145)       !       pfs :  salinity anomaly (s-35) in psu
(  146)       !       pfh :  depth in meters
(  147) 
(  148)       fsalbt( pft, pfs, pfh ) =                                              &
(  149)          ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft                    &
(  150)                                    - 0.203814e-03 ) * pft                    &
(  151)                                    + 0.170907e-01 ) * pft                    &
(  152)                                    + 0.665157e-01                            &
(  153)          +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs   &
(  154)          +  ( ( - 0.302285e-13 * pfh                                         &
(  155)                 - 0.251520e-11 * pfs                                         &
(  156)                 + 0.512857e-12 * pft * pft          ) * pfh                  &
(  157)                                      - 0.164759e-06   * pfs                  &
(  158)              +(   0.791325e-08 * pft - 0.933746e-06 ) * pft                  &
(  159)                                      + 0.380374e-04 ) * pfh   
(  160)       !!----------------------------------------------------------------------
(  161) 
(  162)       IF( kt == nit000 )   CALL tra_bbl_init






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 4

(  163) 
(  164)       ! Save ta and sa trends
(  165)       IF( l_trdtra )   THEN
(  166)          ztdta(:,:,:) = ta(:,:,:) 
(  167)          ztdsa(:,:,:) = sa(:,:,:) 
(  168)       ENDIF
(  169) 
(  170)       ! 0. 2D fields of bottom temperature and salinity, and bottom slope
(  171)       ! -----------------------------------------------------------------
(  172)       ! mbathy= number of w-level, minimum value=1 (cf dommsk.F)
(  173) 
(  174) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  175)       jj = 1
(  176)       DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  177) #  else
(  178)       DO jj = 1, jpj
(  179)          DO ji = 1, jpi
(  180) #  endif
(  181)             ik = mbkt(ji,jj)                              ! index of the bottom ocean T-level
(  182)             ztnb(ji,jj) = tn(ji,jj,ik) * tmask(ji,jj,1)   ! masked now T and S at ocean bottom 
(  183)             zsnb(ji,jj) = sn(ji,jj,ik) * tmask(ji,jj,1)
(  184)             ztbb(ji,jj) = tb(ji,jj,ik) * tmask(ji,jj,1)   ! masked before T and S at ocean bottom 
(  185)             zsbb(ji,jj) = sb(ji,jj,ik) * tmask(ji,jj,1)
(  186)             zdep(ji,jj) = fsdept(ji,jj,ik)                ! depth of the ocean bottom T-level
(  187) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  188)          END DO
(  189) #  endif
(  190)       END DO
(  191) 
(  192) #  if defined key_partial_steps
(  193)       ! partial steps correction 
(  194) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(  195)       jj = 1
(  196)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  197) #   else
(  198)       DO jj = 1, jpjm1
(  199)          DO ji = 1, jpim1
(  200) #   endif
(  201)             iku1 = MAX( mbathy(ji+1,jj  )-1, 1 )
(  202)             iku2 = MAX( mbathy(ji  ,jj  )-1, 1 )
(  203)             ikv1 = MAX( mbathy(ji  ,jj+1)-1, 1 )
(  204)             ikv2 = MAX( mbathy(ji  ,jj  )-1, 1 )
(  205)             ze3u = MIN( fse3u(ji,jj,iku1), fse3u(ji,jj,iku2) ) 
(  206)             ze3v = MIN( fse3v(ji,jj,ikv1), fse3v(ji,jj,ikv2) ) 
(  207)             zahu(ji,jj) = atrbbl * e2u(ji,jj) * ze3u / e1u(ji,jj) * umask(ji,jj,1)
(  208)             zahv(ji,jj) = atrbbl * e1v(ji,jj) * ze3v / e2v(ji,jj) * vmask(ji,jj,1)
(  209) #   if ! defined key_vectopt_loop   ||   defined key_autotasking
(  210)          END DO
(  211) #   endif
(  212)       END DO
(  213) #  else
(  214) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(  215)       jj = 1
(  216)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  217) #   else
(  218)       DO jj = 1, jpjm1
(  219)          DO ji = 1, jpim1
(  220) #   endif






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 5

(  221)             iku = mbku(ji,jj)
(  222)             ikv = mbkv(ji,jj)
(  223)             zahu(ji,jj) = atrbbl * e2u(ji,jj) * fse3u(ji,jj,iku) / e1u(ji,jj) * umask(ji,jj,1)
(  224)             zahv(ji,jj) = atrbbl * e1v(ji,jj) * fse3v(ji,jj,ikv) / e2v(ji,jj) * vmask(ji,jj,1)
(  225) #   if ! defined key_vectopt_loop   ||   defined key_autotasking
(  226)          END DO
(  227) #   endif
(  228)       END DO
(  229) #  endif
(  230) 
(  231)       ! 1. Criteria of additional bottom diffusivity: grad(rho).grad(h)<0
(  232)       ! --------------------------------------------
(  233)       ! Sign of the local density gradient along the i- and j-slopes
(  234)       ! multiplied by the slope of the ocean bottom
(  235) 
(  236)       SELECT CASE ( neos )
(  237) 
(  238)       CASE ( 0   )               ! 0 :Jackett and McDougall (1994) formulation
(  239) 
(  240) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  241)       jj = 1
(  242)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  243) #  else
(  244)       DO jj = 1, jpjm1
(  245)          DO ji = 1, jpim1
(  246) #  endif
(  247)             ! temperature, salinity anomalie and depth
(  248)             zt = 0.5 * ( ztnb(ji,jj) + ztnb(ji+1,jj) )
(  249)             zs = 0.5 * ( zsnb(ji,jj) + zsnb(ji+1,jj) ) - 35.0
(  250)             zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
(  251)             ! masked ratio alpha/beta
(  252)             zalbet = fsalbt( zt, zs, zh )*umask(ji,jj,1)
(  253)             ! local density gradient along i-bathymetric slope
(  254)             zgdrho = zalbet * ( ztnb(ji+1,jj) - ztnb(ji,jj) )   &
(  255)                    -          ( zsnb(ji+1,jj) - zsnb(ji,jj) )
(  256)             ! sign of local i-gradient of density multiplied by the i-slope
(  257)             zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
(  258)             zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
(  259) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  260)          END DO
(  261) #  endif
(  262)       END DO
(  263) 
(  264) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  265)       jj = 1
(  266)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  267) #  else
(  268)       DO jj = 1, jpjm1
(  269)          DO ji = 1, jpim1
(  270) #  endif
(  271)             ! temperature, salinity anomalie and depth
(  272)             zt = 0.5 * ( ztnb(ji,jj+1) + ztnb(ji,jj) )
(  273)             zs = 0.5 * ( zsnb(ji,jj+1) + zsnb(ji,jj) ) - 35.0
(  274)             zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
(  275)             ! masked ratio alpha/beta
(  276)             zalbet = fsalbt( zt, zs, zh )*vmask(ji,jj,1)
(  277)             ! local density gradient along j-bathymetric slope
(  278)             zgdrho = zalbet * ( ztnb(ji,jj+1) - ztnb(ji,jj) )   &






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 6

(  279)                    -          ( zsnb(ji,jj+1) - zsnb(ji,jj) )
(  280)             ! sign of local j-gradient of density multiplied by the j-slope
(  281)             zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
(  282)             zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
(  283) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  284)          END DO
(  285) #  endif
(  286)       END DO
(  287) 
(  288)       CASE ( 1 )               ! Linear formulation function of temperature only
(  289)                                ! 
(  290) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  291)       jj = 1
(  292)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  293) #  else
(  294)       DO jj = 1, jpjm1
(  295)          DO ji = 1, jpim1
(  296) #  endif
(  297)             ! local 'density/temperature' gradient along i-bathymetric slope
(  298)             zgdrho =  ztnb(ji+1,jj) - ztnb(ji,jj) 
(  299)             ! sign of local i-gradient of density multiplied by the i-slope
(  300)             zsign = SIGN( 0.5, - zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
(  301)             zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
(  302) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  303)          END DO
(  304) #  endif
(  305)       END DO
(  306) 
(  307) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  308)       jj = 1
(  309)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  310) #  else
(  311)       DO jj = 1, jpjm1
(  312)          DO ji = 1, jpim1
(  313) #  endif
(  314)             ! local density gradient along j-bathymetric slope
(  315)             zgdrho =  ztnb(ji,jj+1) - ztnb(ji,jj) 
(  316)             ! sign of local j-gradient of density multiplied by the j-slope
(  317)             zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
(  318)             zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
(  319) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  320)          END DO
(  321) #  endif
(  322)       END DO
(  323) 
(  324)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  325) 
(  326)          IF(lwp) WRITE(numout,cform_err)
(  327)          IF(lwp) WRITE(numout,*) '          use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'
(  328)          IF(lwp) WRITE(numout,*) '          bbl not implented: easy to do it '
(  329)          nstop = nstop + 1
(  330) 
(  331)       CASE DEFAULT
(  332) 
(  333)          IF(lwp) WRITE(numout,cform_err)
(  334)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  335)          nstop = nstop + 1
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 7

(  337)       END SELECT
(  338) 
(  339)       ! 2. Additional second order diffusive trends
(  340)       ! -------------------------------------------
(  341) 
(  342)       ! first derivative (gradient)
(  343) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  344)       jj = 1
(  345)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  346) #  else
(  347)       DO jj = 1, jpjm1
(  348)          DO ji = 1, jpim1
(  349) #  endif
(  350)             zkx(ji,jj) = zki(ji,jj) * ( ztbb(ji+1,jj) - ztbb(ji,jj) )
(  351)             zkz(ji,jj) = zki(ji,jj) * ( zsbb(ji+1,jj) - zsbb(ji,jj) )
(  352) 
(  353)             zky(ji,jj) = zkj(ji,jj) * ( ztbb(ji,jj+1) - ztbb(ji,jj) )
(  354)             zkw(ji,jj) = zkj(ji,jj) * ( zsbb(ji,jj+1) - zsbb(ji,jj) )
(  355) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  356)          END DO
(  357) #  endif
(  358)       END DO
(  359) 
(  360) !!DB: deleted ORCA
(  361) !      IF( cp_cfg == "orca" ) THEN
(  362)       ! second derivative (divergence) and add to the general tracer trend
(  363) 
(  364) 
(  365) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  366)       jj = 1
(  367)       DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  368) #  else
(  369)       DO jj = 2, jpjm1
(  370)          DO ji = 2, jpim1
(  371) #  endif
(  372)             ik = max( mbathy(ji,jj)-1, 1 )
(  373)             zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,ik) )
(  374)             zta = (  zkx(ji,jj) - zkx(ji-1,jj  )    &
(  375)                    + zky(ji,jj) - zky(ji  ,jj-1)  ) * zbtr
(  376)             zsa = (  zkz(ji,jj) - zkz(ji-1,jj  )    &
(  377)                    + zkw(ji,jj) - zkw(ji  ,jj-1)  ) * zbtr
(  378)             ta(ji,jj,ik) = ta(ji,jj,ik) + zta
(  379)             sa(ji,jj,ik) = sa(ji,jj,ik) + zsa
(  380) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  381)          END DO
(  382) #  endif
(  383)       END DO
(  384) 
(  385)       ! save the trends for diagnostic
(  386)       ! BBL lateral diffusion tracers trends
(  387)       IF( l_trdtra )   THEN
(  388) #  if defined key_vectopt_loop   &&   ! defined key_autotasking
(  389)          jj = 1
(  390)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  391) #  else
(  392)          DO jj = 2, jpjm1
(  393)             DO ji = 2, jpim1
(  394) #  endif






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 8

(  395)             ik = max( mbathy(ji,jj)-1, 1 )
(  396)             tldfbbl(ji,jj) = ta(ji,jj,ik) - ztdta(ji,jj,ik)
(  397)             sldfbbl(ji,jj) = sa(ji,jj,ik) - ztdsa(ji,jj,ik)
(  398) #  if ! defined key_vectopt_loop   ||   defined key_autotasking
(  399)             END DO
(  400) #  endif
(  401)          END DO
(  402) 
(  403)       ENDIF
(  404) 
(  405)       IF(ln_ctl) THEN
(  406)          CALL prt_ctl(tab3d_1=ta, clinfo1=' bbl  - Ta: ', mask1=tmask, &
(  407)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  408)       ENDIF
(  409) 
(  410)    END SUBROUTINE tra_bbl_dif
(  411) 
(  412) # if defined key_trabbl_adv
(  413)    !!----------------------------------------------------------------------
(  414)    !!   'key_trabbl_adv'                    advective bottom boundary layer
(  415)    !!----------------------------------------------------------------------
(  416) #  include "trabbl_adv.h90"
(  417) # else
(  418)    !!----------------------------------------------------------------------
(  419)    !!   Default option :                 NO advective bottom boundary layer
(  420)    !!----------------------------------------------------------------------
(  421)    SUBROUTINE tra_bbl_adv (kt )              ! Empty routine
(  422)       INTEGER, INTENT(in) :: kt
(  423) !      WRITE(*,*) 'tra_bbl_adv: You should not have seen this print! error?', kt
(  424)    END SUBROUTINE tra_bbl_adv
(  425) # endif
(  426) 
(  427)    SUBROUTINE tra_bbl_init
(  428)       !!----------------------------------------------------------------------
(  429)       !!                  ***  ROUTINE tra_bbl_init  ***
(  430)       !!
(  431)       !! ** Purpose :   Initialization for the bottom boundary layer scheme.
(  432)       !!
(  433)       !! ** Method  :   Read the nambbl namelist and check the parameters
(  434)       !!      called by tra_bbl at the first timestep (nit000)
(  435)       !!
(  436)       !! History :
(  437)       !!    8.5  !  02-08  (G. Madec)  Original code
(  438)       !!----------------------------------------------------------------------
(  439)       !! * Local declarations
(  440)       INTEGER ::   ji, jj      ! dummy loop indices
(  441)       NAMELIST/nambbl/ atrbbl
(  442)       !!----------------------------------------------------------------------
(  443) 
(  444)       ! Read Namelist nambbl : bottom boundary layer scheme
(  445)       ! --------------------
(  446)       REWIND ( numnam )
(  447)       READ   ( numnam, nambbl )
(  448) 
(  449) 
(  450)       ! Parameter control and print
(  451)       ! ---------------------------
(  452)       IF(lwp) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:09      page 9

(  453)          WRITE(numout,*)
(  454)          WRITE(numout,*) 'tra_bbl_init : '
(  455)          WRITE(numout,*) '~~~~~~~~~~~~'
(  456)          IF (lk_trabbl_dif ) THEN
(  457)             WRITE(numout,*) '               * Diffusive Bottom Boundary Layer'
(  458)          ENDIF 
(  459)          IF( lk_trabbl_adv ) THEN
(  460)             WRITE(numout,*) '               * Advective Bottom Boundary Layer'
(  461)          ENDIF
(  462)          WRITE(numout,*) '          Namelist nambbl : set bbl parameters'
(  463)          WRITE(numout,*)
(  464)          WRITE(numout,*) '          bottom boundary layer coef.    atrbbl = ', atrbbl
(  465)          WRITE(numout,*)
(  466)       ENDIF
(  467)  
(  468)       DO jj = 1, jpj
(  469)          DO ji = 1, jpi
(  470)             mbkt(ji,jj) = MAX( mbathy(ji,jj) - 1, 1 )   ! vertical index of the bottom ocean T-level
(  471)          END DO
(  472)       END DO
(  473)       DO jj = 1, jpjm1
(  474)          DO ji = 1, jpim1
(  475)             mbku(ji,jj) = MAX( MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) ) - 1, 1 )
(  476)             mbkv(ji,jj) = MAX( MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) ) - 1, 1 )
(  477)          END DO
(  478)       END DO
(  479) !!bug ???
(  480) !!bug Caution : define the vakue of mbku & mbkv everywhere!!! but lbc mpp lnk : pb when closed (0)
(  481) 
(  482) # if defined key_trabbl_adv
(  483)       ! initialisation of w_bbl to zero
(  484)       w_bbl(:,:,:) = 0.e0    
(  485) # endif
(  486) 
(  487)    END SUBROUTINE tra_bbl_init
(  488) 
(  489) #else
(  490)    !!----------------------------------------------------------------------
(  491)    !!   Dummy module :                      No bottom boundary layer scheme
(  492)    !!----------------------------------------------------------------------
(  493)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_dif = .FALSE.   !: diff bbl flag
(  494)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_adv = .FALSE.   !: adv  bbl flag
(  495) CONTAINS
(  496)    SUBROUTINE tra_bbl_dif (kt )              ! Empty routine
(  497)       INTEGER, INTENT(in) :: kt
(  498) !      WRITE(*,*) 'tra_bbl_dif: You should not have seen this print! error?', kt
(  499)    END SUBROUTINE tra_bbl_dif
(  500)    SUBROUTINE tra_bbl_adv (kt )              ! Empty routine
(  501)       INTEGER, INTENT(in) :: kt
(  502) !      WRITE(*,*) 'tra_bbl_adv: You should not have seen this print! error?', kt
(  503)    END SUBROUTINE tra_bbl_adv
(  504) #endif
(  505) 
(  506)    !!======================================================================
(  507) END MODULE trabbl









PGF90 (Version     10.5)          09/23/2010  09:49:09      page 10

(    1) # 1 "trabbl.F90"
(    1) MODULE trabbl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trabbl  ***
(    4)    !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
(    5)    !!==============================================================================
(    6) # 490
(  490)    !!----------------------------------------------------------------------
(  491)    !!   Dummy module :                      No bottom boundary layer scheme
(  492)    !!----------------------------------------------------------------------
(  493)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_dif = .FALSE.   !: diff bbl flag
(  494)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbl_adv = .FALSE.   !: adv  bbl flag
(  495) CONTAINS
(  496)    SUBROUTINE tra_bbl_dif (kt )              ! Empty routine
(  497)       INTEGER, INTENT(in) :: kt
(  498) !      WRITE(*,*) 'tra_bbl_dif: You should not have seen this print! error?', kt
(  499)    END SUBROUTINE tra_bbl_dif
(  500)    SUBROUTINE tra_bbl_adv (kt )              ! Empty routine
(  501)       INTEGER, INTENT(in) :: kt
(  502) !      WRITE(*,*) 'tra_bbl_adv: You should not have seen this print! error?', kt
(  503)    END SUBROUTINE tra_bbl_adv
(  504) 
(  505) # 506
(  506)    !!======================================================================
(  507) END MODULE trabbl
