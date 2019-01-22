


PGF90 (Version     10.5)          09/23/2010  09:49:37      page 1

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

Filename: dynldf_bilapg.F90

(    1) MODULE dynldf_bilapg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynldf_bilapg  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                              Rotation of mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_ldf_bilapg : update the momentum trend with the horizontal part
(   11)    !!                    of the horizontal s-coord. bilaplacian diffusion
(   12)    !!   ldfguv         : 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   18)    USE zdf_oce         ! ocean vertical physics
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE trdmod          ! ocean dynamics trends 
(   21)    USE trdmod_oce      ! ocean variables trends
(   22)    USE ldfslp          ! iso-neutral slopes available
(   23)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_ldf_bilapg ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "ldfdyn_substitute.h90"
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_bilapg.F90,v 1.8 2005/09/02 15:45:23 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE dyn_ldf_bilapg( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                   ***  ROUTINE dyn_ldf_bilapg  ***
(   46)       !!                      






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 2

(   47)       !! ** Purpose :   Compute the before trend of the horizontal momentum
(   48)       !!      diffusion and add it to the general trend of momentum equation.
(   49)       !!
(   50)       !! ** Method  :   The lateral momentum diffusive trends is provided by a 
(   51)       !!      a 4th order operator rotated along geopotential surfaces. It is 
(   52)       !!      computed using before fields (forward in time) and geopotential
(   53)       !!      slopes computed in routine inildf.
(   54)       !!         -1- compute the geopotential harmonic operator applied to
(   55)       !!      (ub,vb) and multiply it by the eddy diffusivity coefficient
(   56)       !!      (done by a call to ldfgpu and ldfgpv routines) The result is in
(   57)       !!      (zwk1,zwk2) arrays. Applied the domain lateral boundary conditions
(   58)       !!      by call to lbc_lnk.
(   59)       !!         -2- applied to (zwk1,zwk2) the geopotential harmonic operator
(   60)       !!      by a second call to ldfgpu and ldfgpv routines respectively. The
(   61)       !!      result is in (zwk3,zwk4) arrays.
(   62)       !!         -3- Add this trend to the general trend (ta,sa):
(   63)       !!            (ua,va) = (ua,va) + (zwk3,zwk4)
(   64)       !!      'key_trddyn' defined: the trend is saved for diagnostics.
(   65)       !!
(   66)       !! ** Action  : - Update (ua,va) arrays with the before geopotential
(   67)       !!                biharmonic mixing trend.
(   68)       !!              - save the trend in (zwk3,zwk4) ('key_trddyn')
(   69)       !!
(   70)       !! History :
(   71)       !!   8.0  !  97-07  (G. Madec)  Original code
(   72)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   73)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   74)       !!----------------------------------------------------------------------
(   75)       !! * Modules used     
(   76)       USE oce, ONLY :    zwk3 => ta,   & ! use ta as 3D workspace   
(   77)                          zwk4 => sa      ! use sa as 3D workspace   
(   78) 
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(   84)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   85)          zwk1, zwk2                ! work array used for rotated biharmonic
(   86)          !                         ! operator on tracers and/or momentum
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89)       IF( kt == nit000 ) THEN
(   90)          IF(lwp) WRITE(numout,*)
(   91)          IF(lwp) WRITE(numout,*) 'dyn_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
(   92)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(   93)          zwk1(:,:,:) = 0.e0   ;   zwk3(:,:,:) = 0.e0
(   94)          zwk2(:,:,:) = 0.e0   ;   zwk4(:,:,:) = 0.e0
(   95)       ENDIF
(   96) 
(   97)       ! Laplacian of (ub,vb) multiplied by ahm
(   98)       ! --------------------------------------  
(   99)       ! rotated harmonic operator applied to (ub,vb)
(  100)       !     and multiply by ahmu, ahmv (output in (zwk1,zwk2) )
(  101) 
(  102)       CALL ldfguv ( ub, vb, zwk1, zwk2, 1 )
(  103) 
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 3

(  105)       ! Lateral boundary conditions on (zwk1,zwk2)
(  106)       CALL lbc_lnk( zwk1, 'U', -1. )
(  107)       CALL lbc_lnk( zwk2, 'V', -1. )
(  108) 
(  109) 
(  110)       ! Bilaplacian of (ub,vb)
(  111)       ! ---------------------- 
(  112)       ! rotated harmonic operator applied to (zwk1,zwk2) (output in (zwk3,zwk4) )
(  113) 
(  114)       CALL ldfguv ( zwk1, zwk2, zwk3, zwk4, 2 )
(  115) 
(  116) 
(  117)       ! Update the momentum trends			   (j-slab :   2, jpj-1)
(  118)       ! --------------------------
(  119)       !                                                ! ===============
(  120)       DO jj = 2, jpjm1                                 !  Vertical slab
(  121)          !                                             ! ===============
(  122)          DO jk = 1, jpkm1
(  123)             DO ji = 2, jpim1
(  124)                ! add the diffusive trend to the general momentum trends
(  125)                ua(ji,jj,jk) = ua(ji,jj,jk) + zwk3(ji,jj,jk)
(  126)                va(ji,jj,jk) = va(ji,jj,jk) + zwk4(ji,jj,jk)
(  127)             END DO
(  128)          END DO
(  129)          !                                             ! ===============
(  130)       END DO                                           !   End of slab
(  131)       !                                                ! ===============
(  132) 
(  133)       ! save the lateral diffusion trends for diagnostic
(  134)       ! momentum trends
(  135)       IF( l_trddyn )   THEN
(  136) 
(  137)          CALL trd_mod(zwk3, zwk4, jpdtdldf, 'DYN', kt)
(  138)       ENDIF
(  139) 
(  140)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  141)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  142)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  143)       ENDIF
(  144) 
(  145)    END SUBROUTINE dyn_ldf_bilapg
(  146) 
(  147) 
(  148)    SUBROUTINE ldfguv( pu, pv, plu, plv, kahm )
(  149)       !!----------------------------------------------------------------------
(  150)       !!                  ***  ROUTINE ldfguv  ***
(  151)       !!      
(  152)       !! ** Purpose :   Apply a geopotential harmonic operator to (pu,pv)
(  153)       !!      (defined at u- and v-points) and multiply it by the eddy
(  154)       !!      viscosity coefficient (if kahm=1).
(  155)       !!
(  156)       !! ** Method  :   The harmonic operator rotated along geopotential 
(  157)       !!      surfaces is applied to (pu,pv) using the slopes of geopotential
(  158)       !!      surfaces computed in inildf routine. The result is provided in
(  159)       !!      (plu,plv) arrays. It is computed in 2 stepv:
(  160)       !!
(  161)       !!      First step: horizontal part of the operator. It is computed on
(  162)       !!      ==========  pu as follows (idem on pv)






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 4

(  163)       !!      horizontal fluxes :
(  164)       !!         zftu = e2u*e3u/e1u di[ pu ] - e2u*uslp dk[ mi(mk(pu)) ]
(  165)       !!         zftv = e1v*e3v/e2v dj[ pu ] - e1v*vslp dk[ mj(mk(pu)) ]
(  166)       !!      take the horizontal divergence of the fluxes (no divided by
(  167)       !!      the volume element :
(  168)       !!         plu  = di-1[ zftu ] +  dj-1[ zftv ]
(  169)       !!
(  170)       !!      Second step: vertical part of the operator. It is computed on
(  171)       !!      ===========  pu as follows (idem on pv)
(  172)       !!      vertical fluxes :
(  173)       !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pu ]
(  174)       !!              -     e2t     *       wslpi        di[ mi(mk(pu)) ]
(  175)       !!              -     e1t     *       wslpj        dj[ mj(mk(pu)) ]
(  176)       !!      take the vertical divergence of the fluxes add it to the hori-
(  177)       !!      zontal component, divide the result by the volume element and
(  178)       !!      if kahm=1, multiply by the eddy diffusivity coefficient:
(  179)       !!         plu  = aht / (e1t*e2t*e3t) { plu + dk[ zftw ] }
(  180)       !!      else:
(  181)       !!         plu  =  1  / (e1t*e2t*e3t) { plu + dk[ zftw ] }
(  182)       !!
(  183)       !! ** Action :
(  184)       !!        plu, plv        : partial harmonic operator applied to
(  185)       !!                          pu and pv (all the components except
(  186)       !!                          second order vertical derivative term)
(  187)       !!      'key_trddyn' defined: the trend is saved for diagnostics.
(  188)       !!
(  189)       !! History :
(  190)       !!   8.0  !  97-07  (G. Madec)  Original code
(  191)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  192)       !!----------------------------------------------------------------------
(  193)       !! * Arguments
(  194)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  195)          pu, pv     ! momentum fields (before u and v for the 1st call, and
(  196)       !             ! laplacian of these fields multiplied by ahm for the 2nd
(  197)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  198)          plu, plv   ! partial harmonic operator applied to
(  199)       !             ! pu and pv (all the components except
(  200)       !             ! second order vertical derivative term)
(  201)       INTEGER, INTENT( in ) ::   &
(  202)          kahm       ! =1 the laplacian is multiplied by the eddy diffusivity coef.
(  203)       !             ! =2 no multiplication
(  204) 
(  205)       !! * Local declarations
(  206)       INTEGER  ::   ji, jj, jk       ! dummy loop indices
(  207)       REAL(wp) ::   &
(  208)          zabe1, zabe2, zcof1, zcof2,    &  ! temporary scalars
(  209)          zcoef0, zcoef3, zcoef4
(  210)       REAL(wp) ::   &
(  211)          zbur, zbvr, zmkt, zmkf, zuav, zvav,    &
(  212)          zuwslpi, zuwslpj, zvwslpi, zvwslpj
(  213)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  214)          ziut, zjuf , zjvt, zivf,       &  ! workspace
(  215)          zdku, zdk1u, zdkv, zdk1v
(  216)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  217)          zfuw, zfvw, zdiu, zdiv,        &  ! workspace
(  218)          zdju, zdj1u, zdjv, zdj1v 
(  219)       !!----------------------------------------------------------------------
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 5

(  221)       !                               ! ********** !   ! ===============
(  222)       DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
(  223)          !                            ! ********** !   ! ===============
(  224) 
(  225)          ! I.1 Vertical gradient of pu and pv at level jk and jk+1
(  226)          ! -------------------------------------------------------
(  227)          ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
(  228)          !                             zdkv(jk=1)=zdkv(jk=2)
(  229) 
(  230)          zdk1u(:,:) = ( pu(:,:,jk) - pu(:,:,jk+1) ) * umask(:,:,jk+1)
(  231)          zdk1v(:,:) = ( pv(:,:,jk) - pv(:,:,jk+1) ) * vmask(:,:,jk+1)
(  232) 
(  233)          IF( jk == 1 ) THEN
(  234)             zdku(:,:) = zdk1u(:,:)
(  235)             zdkv(:,:) = zdk1v(:,:)
(  236)          ELSE
(  237)             zdku(:,:) = ( pu(:,:,jk-1) - pu(:,:,jk) ) * umask(:,:,jk)
(  238)             zdkv(:,:) = ( pv(:,:,jk-1) - pv(:,:,jk) ) * vmask(:,:,jk)
(  239)          ENDIF
(  240) 
(  241)          !                                -----f-----
(  242)          ! I.2 Horizontal fluxes on U          |
(  243)          ! ------------------------===     t   u   t
(  244)          !                                     |
(  245)          ! i-flux at t-point              -----f-----
(  246)          DO jj = 1, jpjm1
(  247)             DO ji = 2, jpi
(  248)                zabe1 = e2t(ji,jj) * fse3t(ji,jj,jk) / e1t(ji,jj)
(  249) 
(  250)                zmkt  = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
(  251)                               + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )
(  252) 
(  253)                zcof1 = -e2t(ji,jj) * zmkt   &
(  254)                      * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )
(  255) 
(  256)                ziut(ji,jj) = tmask(ji,jj,jk) *   &
(  257)                            (  zabe1 * ( pu(ji,jj,jk) - pu(ji-1,jj,jk) )   &
(  258)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
(  259)                                        +zdk1u(ji,jj) + zdku (ji-1,jj) )  )
(  260)             END DO
(  261)          END DO
(  262) 
(  263)          ! j-flux at f-point
(  264)          DO jj = 1, jpjm1
(  265)             DO ji = 1, jpim1
(  266)                zabe2 = e1f(ji,jj) * fse3f(ji,jj,jk) / e2f(ji,jj)
(  267) 
(  268)                zmkf  = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
(  269)                               + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )
(  270) 
(  271)                zcof2 = -e1f(ji,jj) * zmkf   &
(  272)                      * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )
(  273) 
(  274)                zjuf(ji,jj) = fmask(ji,jj,jk) *   &
(  275)                            (  zabe2 * ( pu(ji,jj+1,jk) - pu(ji,jj,jk) )   &
(  276)                             + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
(  277)                                        +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
(  278)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 6

(  279)          END DO
(  280) 
(  281)          !                                 |   t   |
(  282)          ! I.3 Horizontal fluxes on V      |       |
(  283)          ! ------------------------===     f---v---f
(  284)          !                                 |       |
(  285)          ! i-flux at f-point               |   t   |
(  286)          DO jj = 1, jpjm1
(  287)             DO ji = 1, jpim1
(  288)                zabe1 = e2f(ji,jj) * fse3f(ji,jj,jk) / e1f(ji,jj)
(  289) 
(  290)                zmkf  = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
(  291)                               + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )
(  292) 
(  293)                zcof1 = -e2f(ji,jj) * zmkf   &
(  294)                      * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )
(  295) 
(  296)                zivf(ji,jj) = fmask(ji,jj,jk) *   &
(  297)                            (  zabe1 * ( pu(ji+1,jj,jk) - pu(ji,jj,jk) )   &
(  298)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji+1,jj)     &
(  299)                                        +zdk1u(ji,jj) + zdku (ji+1,jj) )  )
(  300)             END DO
(  301)          END DO
(  302) 
(  303)          ! j-flux at t-point
(  304)          DO jj = 2, jpj
(  305)             DO ji = 1, jpim1
(  306)                zabe2 = e1t(ji,jj) * fse3t(ji,jj,jk) / e2t(ji,jj)
(  307) 
(  308)                zmkt  = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
(  309)                               + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )
(  310) 
(  311)                zcof2 = -e1t(ji,jj) * zmkt   &
(  312)                      * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )
(  313) 
(  314)                zjvt(ji,jj) = tmask(ji,jj,jk) *   &
(  315)                            (  zabe2 * ( pu(ji,jj,jk) - pu(ji,jj-1,jk) )   &
(  316)                             + zcof2 * ( zdku (ji,jj-1) + zdk1u(ji,jj)     &
(  317)                                        +zdk1u(ji,jj-1) + zdku (ji,jj) )  )
(  318)             END DO
(  319)          END DO
(  320) 
(  321) 
(  322)          ! I.4 Second derivative (divergence) (not divided by the volume)
(  323)          ! ---------------------
(  324) 
(  325)          DO jj = 2, jpjm1
(  326)             DO ji = 2, jpim1
(  327)                plu(ji,jj,jk) = ziut (ji+1,jj) - ziut (ji,jj  )   &
(  328)                              + zjuf (ji  ,jj) - zjuf (ji,jj-1)
(  329)                plv(ji,jj,jk) = zivf (ji,jj  ) - zivf (ji-1,jj)   &
(  330)                              + zjvt (ji,jj+1) - zjvt (ji,jj  ) 
(  331)             END DO
(  332)          END DO
(  333) 
(  334)          !                                             ! ===============
(  335)       END DO                                           !   End of slab
(  336)       !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 7

(  337) 
(  338)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  339) 
(  340)       !                             ! ************ !   ! ===============
(  341)       DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
(  342)          !                          ! ************ !   ! ===============
(  343) 
(  344)          ! II.1 horizontal (pu,pv) gradients
(  345)          ! ---------------------------------
(  346) 
(  347)          DO jk = 1, jpk
(  348)             DO ji = 2, jpi
(  349)                ! i-gradient of u at jj
(  350)                zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( pu(ji,jj  ,jk) - pu(ji-1,jj  ,jk) )
(  351)                ! j-gradient of u and v at jj
(  352)                zdju (ji,jk) = fmask(ji,jj  ,jk) * ( pu(ji,jj+1,jk) - pu(ji  ,jj  ,jk) )
(  353)                zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( pv(ji,jj  ,jk) - pv(ji  ,jj-1,jk) )
(  354)                ! j-gradient of u and v at jj+1
(  355)                zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( pu(ji,jj  ,jk) - pu(ji  ,jj-1,jk) )
(  356)                zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( pv(ji,jj+1,jk) - pv(ji  ,jj  ,jk) )
(  357)             END DO
(  358)          END DO
(  359)          DO jk = 1, jpk
(  360)             DO ji = 1, jpim1
(  361)                ! i-gradient of v at jj
(  362)                zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( pv(ji+1,jj,jk) - pv(ji  ,jj  ,jk) )
(  363)             END DO
(  364)          END DO
(  365) 
(  366) 
(  367)          ! II.2 Vertical fluxes
(  368)          ! --------------------
(  369) 
(  370)          ! Surface and bottom vertical fluxes set to zero
(  371) 
(  372)          zfuw(:, 1 ) = 0.e0
(  373)          zfvw(:, 1 ) = 0.e0
(  374)          zfuw(:,jpk) = 0.e0
(  375)          zfvw(:,jpk) = 0.e0
(  376) 
(  377)          ! interior (2=<jk=<jpk-1) on pu field
(  378) 
(  379)          DO jk = 2, jpkm1
(  380)             DO ji = 2, jpim1
(  381)                ! i- and j-slopes at uw-point
(  382)                zuwslpi = 0.5 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
(  383)                zuwslpj = 0.5 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )
(  384)                ! coef. for the vertical dirative
(  385)                zcoef0 = e1u(ji,jj) * e2u(ji,jj) / fse3u(ji,jj,jk)   &
(  386)                       * ( zuwslpi * zuwslpi + zuwslpj * zuwslpj )
(  387)                ! weights for the i-k, j-k averaging at t- and f-points, resp.
(  388)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
(  389)                              + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
(  390)                zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
(  391)                              + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )
(  392)                ! coef. for the horitontal derivative
(  393)                zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
(  394)                zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 8

(  395)                ! vertical flux on u field
(  396)                zfuw(ji,jk) = umask(ji,jj,jk) *   &
(  397)                            (  zcoef0 * ( pu  (ji,jj,jk-1) - pu  (ji,jj,jk) )   &
(  398)                             + zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
(  399)                                         +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
(  400)                             + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
(  401)                                         +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) ) )
(  402)             END DO
(  403)          END DO
(  404) 
(  405)          ! interior (2=<jk=<jpk-1) on pv field
(  406) 
(  407)          DO jk = 2, jpkm1
(  408)             DO ji = 2, jpim1
(  409)                ! i- and j-slopes at vw-point
(  410)                zvwslpi = 0.5 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
(  411)                zvwslpj = 0.5 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )
(  412)                ! coef. for the vertical derivative
(  413)                zcoef0 = e1v(ji,jj) * e2v(ji,jj) / fse3v(ji,jj,jk)   &
(  414)                       * ( zvwslpi * zvwslpi + zvwslpj * zvwslpj )
(  415)                ! weights for the i-k, j-k averaging at f- and t-points, resp.
(  416)                zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
(  417)                              + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
(  418)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
(  419)                              + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )
(  420)                ! coef. for the horizontal derivatives
(  421)                zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
(  422)                zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
(  423)                ! vertical flux on pv field
(  424)                zfvw(ji,jk) = vmask(ji,jj,jk) *   &
(  425)                            (  zcoef0 * ( pv  (ji,jj,jk-1) - pv  (ji,jj,jk) )   &
(  426)                             + zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
(  427)                                         +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
(  428)                             + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
(  429)                                         +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )  )
(  430)             END DO
(  431)          END DO
(  432) 
(  433) 
(  434)          ! II.3 Divergence of vertical fluxes added to the horizontal divergence
(  435)          ! ---------------------------------------------------------------------
(  436) 
(  437)          IF( kahm == 1 ) THEN
(  438)             ! multiply the laplacian by the eddy viscosity coefficient
(  439)             DO jk = 1, jpkm1
(  440)                DO ji = 2, jpim1
(  441)                   ! eddy coef. divided by the volume element
(  442)                   zbur = fsahmu(ji,jj,jk) / ( e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) )
(  443)                   zbvr = fsahmv(ji,jj,jk) / ( e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) )
(  444)                   ! vertical divergence
(  445)                   zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
(  446)                   zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
(  447)                   ! harmonic operator applied to (pu,pv) and multiply by ahm
(  448)                   plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
(  449)                   plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
(  450)                END DO
(  451)             END DO
(  452)          ELSEIF( kahm == 2 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 9

(  453)             ! second call, no multiplication
(  454)             DO jk = 1, jpkm1
(  455)                DO ji = 2, jpim1
(  456)                   ! inverse of the volume element
(  457)                   zbur = 1. / ( e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) )
(  458)                   zbvr = 1. / ( e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) )
(  459)                   ! vertical divergence
(  460)                   zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
(  461)                   zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
(  462)                   ! harmonic operator applied to (pu,pv) 
(  463)                   plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
(  464)                   plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
(  465)                END DO
(  466)             END DO
(  467)          ELSE
(  468)             IF(lwp)WRITE(numout,*) ' ldfguv: kahm= 1 or 2, here =', kahm
(  469)             IF(lwp)WRITE(numout,*) '         We stop'
(  470)             STOP 'ldfguv'
(  471)          ENDIF
(  472)          !                                             ! ===============
(  473)       END DO                                           !   End of slab
(  474)       !                                                ! ===============
(  475)    END SUBROUTINE ldfguv
(  476) 
(  477) #else
(  478)    !!----------------------------------------------------------------------
(  479)    !!   Dummy module :                         NO rotation of mixing tensor
(  480)    !!----------------------------------------------------------------------
(  481)    USE in_out_manager
(  482) CONTAINS
(  483)    SUBROUTINE dyn_ldf_bilapg( kt )               ! Dummy routine
(  484)       if(lwp) WRITE(numout,*) 'dyn_ldf_bilapg: You should not have seen this print! error?', kt
(  485)    END SUBROUTINE dyn_ldf_bilapg
(  486) #endif
(  487) 
(  488)    !!======================================================================
(  489) END MODULE dynldf_bilapg



























PGF90 (Version     10.5)          09/23/2010  09:49:37      page 10

(    1) # 1 "dynldf_bilapg.F90"
(    1) MODULE dynldf_bilapg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynldf_bilapg  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                              Rotation of mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dyn_ldf_bilapg : update the momentum trend with the horizontal part
(   11)    !!                    of the horizontal s-coord. bilaplacian diffusion
(   12)    !!   ldfguv         : 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldfdyn_oce      ! ocean dynamics lateral physics
(   18)    USE zdf_oce         ! ocean vertical physics
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE trdmod          ! ocean dynamics trends 
(   21)    USE trdmod_oce      ! ocean variables trends
(   22)    USE ldfslp          ! iso-neutral slopes available
(   23)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_ldf_bilapg ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) # 1 "./domzgr_substitute.h90"
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
(   48) # 34 "dynldf_bilapg.F90"
(   34) # 1 "./ldfdyn_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  ldfdyn_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
(    5)    !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 11

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_substitute.h90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   ' key_dynldf_c3d' :                  3D coefficient
(   20) # 35 "dynldf_bilapg.F90"
(   35) # 35
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_bilapg.F90,v 1.8 2005/09/02 15:45:23 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE dyn_ldf_bilapg( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                   ***  ROUTINE dyn_ldf_bilapg  ***
(   46)       !!                      
(   47)       !! ** Purpose :   Compute the before trend of the horizontal momentum
(   48)       !!      diffusion and add it to the general trend of momentum equation.
(   49)       !!
(   50)       !! ** Method  :   The lateral momentum diffusive trends is provided by a 
(   51)       !!      a 4th order operator rotated along geopotential surfaces. It is 
(   52)       !!      computed using before fields (forward in time) and geopotential
(   53)       !!      slopes computed in routine inildf.
(   54)       !!         -1- compute the geopotential harmonic operator applied to
(   55)       !!      (ub,vb) and multiply it by the eddy diffusivity coefficient
(   56)       !!      (done by a call to ldfgpu and ldfgpv routines) The result is in
(   57)       !!      (zwk1,zwk2) arrays. Applied the domain lateral boundary conditions
(   58)       !!      by call to lbc_lnk.
(   59)       !!         -2- applied to (zwk1,zwk2) the geopotential harmonic operator
(   60)       !!      by a second call to ldfgpu and ldfgpv routines respectively. The
(   61)       !!      result is in (zwk3,zwk4) arrays.
(   62)       !!         -3- Add this trend to the general trend (ta,sa):
(   63)       !!            (ua,va) = (ua,va) + (zwk3,zwk4)
(   64)       !!      'key_trddyn' defined: the trend is saved for diagnostics.
(   65)       !!
(   66)       !! ** Action  : - Update (ua,va) arrays with the before geopotential
(   67)       !!                biharmonic mixing trend.
(   68)       !!              - save the trend in (zwk3,zwk4) ('key_trddyn')
(   69)       !!
(   70)       !! History :
(   71)       !!   8.0  !  97-07  (G. Madec)  Original code
(   72)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   73)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   74)       !!----------------------------------------------------------------------
(   75)       !! * Modules used     
(   76)       USE oce, ONLY :    zwk3 => ta,   & ! use ta as 3D workspace   
(   77)                          zwk4 => sa      ! use sa as 3D workspace   
(   78) 
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(   84)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   85)          zwk1, zwk2                ! work array used for rotated biharmonic






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 12

(   86)          !                         ! operator on tracers and/or momentum
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89)       IF( kt == nit000 ) THEN
(   90)          IF(lwp) WRITE(numout,*)
(   91)          IF(lwp) WRITE(numout,*) 'dyn_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
(   92)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(   93)          zwk1(:,:,:) = 0.e0   ;   zwk3(:,:,:) = 0.e0
(   94)          zwk2(:,:,:) = 0.e0   ;   zwk4(:,:,:) = 0.e0
(   95)       ENDIF
(   96) 
(   97)       ! Laplacian of (ub,vb) multiplied by ahm
(   98)       ! --------------------------------------  
(   99)       ! rotated harmonic operator applied to (ub,vb)
(  100)       !     and multiply by ahmu, ahmv (output in (zwk1,zwk2) )
(  101) 
(  102)       CALL ldfguv ( ub, vb, zwk1, zwk2, 1 )
(  103) 
(  104) 
(  105)       ! Lateral boundary conditions on (zwk1,zwk2)
(  106)       CALL lbc_lnk( zwk1, 'U', -1. )
(  107)       CALL lbc_lnk( zwk2, 'V', -1. )
(  108) 
(  109) 
(  110)       ! Bilaplacian of (ub,vb)
(  111)       ! ---------------------- 
(  112)       ! rotated harmonic operator applied to (zwk1,zwk2) (output in (zwk3,zwk4) )
(  113) 
(  114)       CALL ldfguv ( zwk1, zwk2, zwk3, zwk4, 2 )
(  115) 
(  116) 
(  117)       ! Update the momentum trends			   (j-slab :   2, jpj-1)
(  118)       ! --------------------------
(  119)       !                                                ! ===============
(  120)       DO jj = 2, jpjm1                                 !  Vertical slab
(  121)          !                                             ! ===============
(  122)          DO jk = 1, jpkm1
(  123)             DO ji = 2, jpim1
(  124)                ! add the diffusive trend to the general momentum trends
(  125)                ua(ji,jj,jk) = ua(ji,jj,jk) + zwk3(ji,jj,jk)
(  126)                va(ji,jj,jk) = va(ji,jj,jk) + zwk4(ji,jj,jk)
(  127)             END DO
(  128)          END DO
(  129)          !                                             ! ===============
(  130)       END DO                                           !   End of slab
(  131)       !                                                ! ===============
(  132) 
(  133)       ! save the lateral diffusion trends for diagnostic
(  134)       ! momentum trends
(  135)       IF( l_trddyn )   THEN
(  136) 
(  137)          CALL trd_mod(zwk3, zwk4, jpdtdldf, 'DYN', kt)
(  138)       ENDIF
(  139) 
(  140)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  141)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  142)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  143)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 13

(  144) 
(  145)    END SUBROUTINE dyn_ldf_bilapg
(  146) 
(  147) 
(  148)    SUBROUTINE ldfguv( pu, pv, plu, plv, kahm )
(  149)       !!----------------------------------------------------------------------
(  150)       !!                  ***  ROUTINE ldfguv  ***
(  151)       !!      
(  152)       !! ** Purpose :   Apply a geopotential harmonic operator to (pu,pv)
(  153)       !!      (defined at u- and v-points) and multiply it by the eddy
(  154)       !!      viscosity coefficient (if kahm=1).
(  155)       !!
(  156)       !! ** Method  :   The harmonic operator rotated along geopotential 
(  157)       !!      surfaces is applied to (pu,pv) using the slopes of geopotential
(  158)       !!      surfaces computed in inildf routine. The result is provided in
(  159)       !!      (plu,plv) arrays. It is computed in 2 stepv:
(  160)       !!
(  161)       !!      First step: horizontal part of the operator. It is computed on
(  162)       !!      ==========  pu as follows (idem on pv)
(  163)       !!      horizontal fluxes :
(  164)       !!         zftu = e2u*e3u/e1u di[ pu ] - e2u*uslp dk[ mi(mk(pu)) ]
(  165)       !!         zftv = e1v*e3v/e2v dj[ pu ] - e1v*vslp dk[ mj(mk(pu)) ]
(  166)       !!      take the horizontal divergence of the fluxes (no divided by
(  167)       !!      the volume element :
(  168)       !!         plu  = di-1[ zftu ] +  dj-1[ zftv ]
(  169)       !!
(  170)       !!      Second step: vertical part of the operator. It is computed on
(  171)       !!      ===========  pu as follows (idem on pv)
(  172)       !!      vertical fluxes :
(  173)       !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pu ]
(  174)       !!              -     e2t     *       wslpi        di[ mi(mk(pu)) ]
(  175)       !!              -     e1t     *       wslpj        dj[ mj(mk(pu)) ]
(  176)       !!      take the vertical divergence of the fluxes add it to the hori-
(  177)       !!      zontal component, divide the result by the volume element and
(  178)       !!      if kahm=1, multiply by the eddy diffusivity coefficient:
(  179)       !!         plu  = aht / (e1t*e2t*e3t) { plu + dk[ zftw ] }
(  180)       !!      else:
(  181)       !!         plu  =  1  / (e1t*e2t*e3t) { plu + dk[ zftw ] }
(  182)       !!
(  183)       !! ** Action :
(  184)       !!        plu, plv        : partial harmonic operator applied to
(  185)       !!                          pu and pv (all the components except
(  186)       !!                          second order vertical derivative term)
(  187)       !!      'key_trddyn' defined: the trend is saved for diagnostics.
(  188)       !!
(  189)       !! History :
(  190)       !!   8.0  !  97-07  (G. Madec)  Original code
(  191)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  192)       !!----------------------------------------------------------------------
(  193)       !! * Arguments
(  194)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  195)          pu, pv     ! momentum fields (before u and v for the 1st call, and
(  196)       !             ! laplacian of these fields multiplied by ahm for the 2nd
(  197)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  198)          plu, plv   ! partial harmonic operator applied to
(  199)       !             ! pu and pv (all the components except
(  200)       !             ! second order vertical derivative term)
(  201)       INTEGER, INTENT( in ) ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 14

(  202)          kahm       ! =1 the laplacian is multiplied by the eddy diffusivity coef.
(  203)       !             ! =2 no multiplication
(  204) 
(  205)       !! * Local declarations
(  206)       INTEGER  ::   ji, jj, jk       ! dummy loop indices
(  207)       REAL(wp) ::   &
(  208)          zabe1, zabe2, zcof1, zcof2,    &  ! temporary scalars
(  209)          zcoef0, zcoef3, zcoef4
(  210)       REAL(wp) ::   &
(  211)          zbur, zbvr, zmkt, zmkf, zuav, zvav,    &
(  212)          zuwslpi, zuwslpj, zvwslpi, zvwslpj
(  213)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  214)          ziut, zjuf , zjvt, zivf,       &  ! workspace
(  215)          zdku, zdk1u, zdkv, zdk1v
(  216)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  217)          zfuw, zfvw, zdiu, zdiv,        &  ! workspace
(  218)          zdju, zdj1u, zdjv, zdj1v 
(  219)       !!----------------------------------------------------------------------
(  220) 
(  221)       !                               ! ********** !   ! ===============
(  222)       DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
(  223)          !                            ! ********** !   ! ===============
(  224) 
(  225)          ! I.1 Vertical gradient of pu and pv at level jk and jk+1
(  226)          ! -------------------------------------------------------
(  227)          ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
(  228)          !                             zdkv(jk=1)=zdkv(jk=2)
(  229) 
(  230)          zdk1u(:,:) = ( pu(:,:,jk) - pu(:,:,jk+1) ) * umask(:,:,jk+1)
(  231)          zdk1v(:,:) = ( pv(:,:,jk) - pv(:,:,jk+1) ) * vmask(:,:,jk+1)
(  232) 
(  233)          IF( jk == 1 ) THEN
(  234)             zdku(:,:) = zdk1u(:,:)
(  235)             zdkv(:,:) = zdk1v(:,:)
(  236)          ELSE
(  237)             zdku(:,:) = ( pu(:,:,jk-1) - pu(:,:,jk) ) * umask(:,:,jk)
(  238)             zdkv(:,:) = ( pv(:,:,jk-1) - pv(:,:,jk) ) * vmask(:,:,jk)
(  239)          ENDIF
(  240) 
(  241)          !                                -----f-----
(  242)          ! I.2 Horizontal fluxes on U          |
(  243)          ! ------------------------===     t   u   t
(  244)          !                                     |
(  245)          ! i-flux at t-point              -----f-----
(  246)          DO jj = 1, jpjm1
(  247)             DO ji = 2, jpi
(  248)                zabe1 = e2t(ji,jj) * e3t(jk) / e1t(ji,jj)
(  249) 
(  250)                zmkt  = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
(  251)                               + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )
(  252) 
(  253)                zcof1 = -e2t(ji,jj) * zmkt   &
(  254)                      * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )
(  255) 
(  256)                ziut(ji,jj) = tmask(ji,jj,jk) *   &
(  257)                            (  zabe1 * ( pu(ji,jj,jk) - pu(ji-1,jj,jk) )   &
(  258)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
(  259)                                        +zdk1u(ji,jj) + zdku (ji-1,jj) )  )






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 15

(  260)             END DO
(  261)          END DO
(  262) 
(  263)          ! j-flux at f-point
(  264)          DO jj = 1, jpjm1
(  265)             DO ji = 1, jpim1
(  266)                zabe2 = e1f(ji,jj) * e3t(jk) / e2f(ji,jj)
(  267) 
(  268)                zmkf  = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
(  269)                               + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )
(  270) 
(  271)                zcof2 = -e1f(ji,jj) * zmkf   &
(  272)                      * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )
(  273) 
(  274)                zjuf(ji,jj) = fmask(ji,jj,jk) *   &
(  275)                            (  zabe2 * ( pu(ji,jj+1,jk) - pu(ji,jj,jk) )   &
(  276)                             + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
(  277)                                        +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
(  278)             END DO
(  279)          END DO
(  280) 
(  281)          !                                 |   t   |
(  282)          ! I.3 Horizontal fluxes on V      |       |
(  283)          ! ------------------------===     f---v---f
(  284)          !                                 |       |
(  285)          ! i-flux at f-point               |   t   |
(  286)          DO jj = 1, jpjm1
(  287)             DO ji = 1, jpim1
(  288)                zabe1 = e2f(ji,jj) * e3t(jk) / e1f(ji,jj)
(  289) 
(  290)                zmkf  = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
(  291)                               + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )
(  292) 
(  293)                zcof1 = -e2f(ji,jj) * zmkf   &
(  294)                      * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )
(  295) 
(  296)                zivf(ji,jj) = fmask(ji,jj,jk) *   &
(  297)                            (  zabe1 * ( pu(ji+1,jj,jk) - pu(ji,jj,jk) )   &
(  298)                             + zcof1 * ( zdku (ji,jj) + zdk1u(ji+1,jj)     &
(  299)                                        +zdk1u(ji,jj) + zdku (ji+1,jj) )  )
(  300)             END DO
(  301)          END DO
(  302) 
(  303)          ! j-flux at t-point
(  304)          DO jj = 2, jpj
(  305)             DO ji = 1, jpim1
(  306)                zabe2 = e1t(ji,jj) * e3t(jk) / e2t(ji,jj)
(  307) 
(  308)                zmkt  = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
(  309)                               + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )
(  310) 
(  311)                zcof2 = -e1t(ji,jj) * zmkt   &
(  312)                      * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )
(  313) 
(  314)                zjvt(ji,jj) = tmask(ji,jj,jk) *   &
(  315)                            (  zabe2 * ( pu(ji,jj,jk) - pu(ji,jj-1,jk) )   &
(  316)                             + zcof2 * ( zdku (ji,jj-1) + zdk1u(ji,jj)     &
(  317)                                        +zdk1u(ji,jj-1) + zdku (ji,jj) )  )






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 16

(  318)             END DO
(  319)          END DO
(  320) 
(  321) 
(  322)          ! I.4 Second derivative (divergence) (not divided by the volume)
(  323)          ! ---------------------
(  324) 
(  325)          DO jj = 2, jpjm1
(  326)             DO ji = 2, jpim1
(  327)                plu(ji,jj,jk) = ziut (ji+1,jj) - ziut (ji,jj  )   &
(  328)                              + zjuf (ji  ,jj) - zjuf (ji,jj-1)
(  329)                plv(ji,jj,jk) = zivf (ji,jj  ) - zivf (ji-1,jj)   &
(  330)                              + zjvt (ji,jj+1) - zjvt (ji,jj  ) 
(  331)             END DO
(  332)          END DO
(  333) 
(  334)          !                                             ! ===============
(  335)       END DO                                           !   End of slab
(  336)       !                                                ! ===============
(  337) 
(  338)       !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
(  339) 
(  340)       !                             ! ************ !   ! ===============
(  341)       DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
(  342)          !                          ! ************ !   ! ===============
(  343) 
(  344)          ! II.1 horizontal (pu,pv) gradients
(  345)          ! ---------------------------------
(  346) 
(  347)          DO jk = 1, jpk
(  348)             DO ji = 2, jpi
(  349)                ! i-gradient of u at jj
(  350)                zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( pu(ji,jj  ,jk) - pu(ji-1,jj  ,jk) )
(  351)                ! j-gradient of u and v at jj
(  352)                zdju (ji,jk) = fmask(ji,jj  ,jk) * ( pu(ji,jj+1,jk) - pu(ji  ,jj  ,jk) )
(  353)                zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( pv(ji,jj  ,jk) - pv(ji  ,jj-1,jk) )
(  354)                ! j-gradient of u and v at jj+1
(  355)                zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( pu(ji,jj  ,jk) - pu(ji  ,jj-1,jk) )
(  356)                zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( pv(ji,jj+1,jk) - pv(ji  ,jj  ,jk) )
(  357)             END DO
(  358)          END DO
(  359)          DO jk = 1, jpk
(  360)             DO ji = 1, jpim1
(  361)                ! i-gradient of v at jj
(  362)                zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( pv(ji+1,jj,jk) - pv(ji  ,jj  ,jk) )
(  363)             END DO
(  364)          END DO
(  365) 
(  366) 
(  367)          ! II.2 Vertical fluxes
(  368)          ! --------------------
(  369) 
(  370)          ! Surface and bottom vertical fluxes set to zero
(  371) 
(  372)          zfuw(:, 1 ) = 0.e0
(  373)          zfvw(:, 1 ) = 0.e0
(  374)          zfuw(:,jpk) = 0.e0
(  375)          zfvw(:,jpk) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 17

(  376) 
(  377)          ! interior (2=<jk=<jpk-1) on pu field
(  378) 
(  379)          DO jk = 2, jpkm1
(  380)             DO ji = 2, jpim1
(  381)                ! i- and j-slopes at uw-point
(  382)                zuwslpi = 0.5 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
(  383)                zuwslpj = 0.5 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )
(  384)                ! coef. for the vertical dirative
(  385)                zcoef0 = e1u(ji,jj) * e2u(ji,jj) / e3t(jk)   &
(  386)                       * ( zuwslpi * zuwslpi + zuwslpj * zuwslpj )
(  387)                ! weights for the i-k, j-k averaging at t- and f-points, resp.
(  388)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
(  389)                              + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
(  390)                zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
(  391)                              + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )
(  392)                ! coef. for the horitontal derivative
(  393)                zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
(  394)                zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj
(  395)                ! vertical flux on u field
(  396)                zfuw(ji,jk) = umask(ji,jj,jk) *   &
(  397)                            (  zcoef0 * ( pu  (ji,jj,jk-1) - pu  (ji,jj,jk) )   &
(  398)                             + zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
(  399)                                         +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
(  400)                             + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
(  401)                                         +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) ) )
(  402)             END DO
(  403)          END DO
(  404) 
(  405)          ! interior (2=<jk=<jpk-1) on pv field
(  406) 
(  407)          DO jk = 2, jpkm1
(  408)             DO ji = 2, jpim1
(  409)                ! i- and j-slopes at vw-point
(  410)                zvwslpi = 0.5 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
(  411)                zvwslpj = 0.5 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )
(  412)                ! coef. for the vertical derivative
(  413)                zcoef0 = e1v(ji,jj) * e2v(ji,jj) / e3t(jk)   &
(  414)                       * ( zvwslpi * zvwslpi + zvwslpj * zvwslpj )
(  415)                ! weights for the i-k, j-k averaging at f- and t-points, resp.
(  416)                zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
(  417)                              + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
(  418)                zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
(  419)                              + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )
(  420)                ! coef. for the horizontal derivatives
(  421)                zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
(  422)                zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
(  423)                ! vertical flux on pv field
(  424)                zfvw(ji,jk) = vmask(ji,jj,jk) *   &
(  425)                            (  zcoef0 * ( pv  (ji,jj,jk-1) - pv  (ji,jj,jk) )   &
(  426)                             + zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
(  427)                                         +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
(  428)                             + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
(  429)                                         +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )  )
(  430)             END DO
(  431)          END DO
(  432) 
(  433) 






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 18

(  434)          ! II.3 Divergence of vertical fluxes added to the horizontal divergence
(  435)          ! ---------------------------------------------------------------------
(  436) 
(  437)          IF( kahm == 1 ) THEN
(  438)             ! multiply the laplacian by the eddy viscosity coefficient
(  439)             DO jk = 1, jpkm1
(  440)                DO ji = 2, jpim1
(  441)                   ! eddy coef. divided by the volume element
(  442)                   zbur = ahm3(ji,jj,jk) / ( e1u(ji,jj)*e2u(ji,jj)*e3t(jk) )
(  443)                   zbvr = ahm4(ji,jj,jk) / ( e1v(ji,jj)*e2v(ji,jj)*e3t(jk) )
(  444)                   ! vertical divergence
(  445)                   zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
(  446)                   zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
(  447)                   ! harmonic operator applied to (pu,pv) and multiply by ahm
(  448)                   plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
(  449)                   plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
(  450)                END DO
(  451)             END DO
(  452)          ELSEIF( kahm == 2 ) THEN
(  453)             ! second call, no multiplication
(  454)             DO jk = 1, jpkm1
(  455)                DO ji = 2, jpim1
(  456)                   ! inverse of the volume element
(  457)                   zbur = 1. / ( e1u(ji,jj)*e2u(ji,jj)*e3t(jk) )
(  458)                   zbvr = 1. / ( e1v(ji,jj)*e2v(ji,jj)*e3t(jk) )
(  459)                   ! vertical divergence
(  460)                   zuav = zfuw(ji,jk) - zfuw(ji,jk+1)
(  461)                   zvav = zfvw(ji,jk) - zfvw(ji,jk+1)
(  462)                   ! harmonic operator applied to (pu,pv) 
(  463)                   plu(ji,jj,jk) = ( plu(ji,jj,jk) + zuav ) * zbur
(  464)                   plv(ji,jj,jk) = ( plv(ji,jj,jk) + zvav ) * zbvr
(  465)                END DO
(  466)             END DO
(  467)          ELSE
(  468)             IF(lwp)WRITE(numout,*) ' ldfguv: kahm= 1 or 2, here =', kahm
(  469)             IF(lwp)WRITE(numout,*) '         We stop'
(  470)             STOP 'ldfguv'
(  471)          ENDIF
(  472)          !                                             ! ===============
(  473)       END DO                                           !   End of slab
(  474)       !                                                ! ===============
(  475)    END SUBROUTINE ldfguv
(  476) 
(  477) 
(  478) # 488
(  488)    !!======================================================================
(  489) END MODULE dynldf_bilapg
