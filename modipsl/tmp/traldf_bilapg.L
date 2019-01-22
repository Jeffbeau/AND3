


PGF90 (Version     10.5)          09/23/2010  09:50:13      page 1

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

Filename: traldf_bilapg.F90

(    1) MODULE traldf_bilapg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traldf_bilapg  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_ldf_bilapg : update the tracer trend with the horizontal diffusion
(   11)    !!                    using an horizontal biharmonic operator in s-coordinate 
(   12)    !!   ldfght         :  ???
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE ldfslp          ! iso-neutral slopes available
(   22)    USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(   23)    USE diaptr          ! poleward transport diagnostics 
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC tra_ldf_bilapg    ! routine called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "ldftra_substitute.h90"
(   35) #  include "ldfeiv_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_bilapg.F90,v 1.7 2005/09/02 15:45:33 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41)    
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE tra_ldf_bilapg( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                 ***  ROUTINE tra_ldf_bilapg  ***






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 2

(   47)       !!                    
(   48)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   49)       !!      trend and add it to the general trend of tracer equation.
(   50)       !!
(   51)       !! ** Method  :   The lateral diffusive trends is provided by a 4th order
(   52)       !!      operator rotated along geopotential surfaces. It is computed
(   53)       !!      using before fields (forward in time) and geopotential slopes
(   54)       !!      computed in routine inildf.
(   55)       !!         -1- compute the geopotential harmonic operator applied to
(   56)       !!      (tb,sb) and multiply it by the eddy diffusivity coefficient
(   57)       !!      (done by a call to ldfght routine, result in (wk1,wk2) arrays).
(   58)       !!      Applied the domain lateral boundary conditions by call to lbc_lnk
(   59)       !!         -2- compute the geopotential harmonic operator applied to
(   60)       !!      (wk1,wk2) by a second call to ldfght routine (result in (wk3,wk4)
(   61)       !!      arrays).
(   62)       !!         -3- Add this trend to the general trend (ta,sa):
(   63)       !!            (ta,sa) = (ta,sa) + (wk3,wk4)
(   64)       !!
(   65)       !! ** Action : - Update (ta,sa) arrays with the before geopotential 
(   66)       !!               biharmonic mixing trend.
(   67)       !!             - Save the trends  in (ztdta,ztdsa) ('key_trdtra')
(   68)       !!
(   69)       !! History :
(   70)       !!   8.0  !  97-07  (G. Madec)  Original code
(   71)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   72)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   73)       !!----------------------------------------------------------------------
(   74)       !! * Modules used
(   75)       USE oce           , wk1 => ua,  &  ! use ua as workspace
(   76)          &                wk2 => va      ! use va as workspace
(   77) 
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(   83)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   84)          wk3, wk4                ! work array used for rotated biharmonic
(   85)          !                       ! operator on tracers and/or momentum
(   86)       !!----------------------------------------------------------------------
(   87) 
(   88)       IF( kt == nit000 ) THEN
(   89)          IF(lwp) WRITE(numout,*)
(   90)          IF(lwp) WRITE(numout,*) 'tra_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
(   91)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(   92)       ENDIF
(   93) 
(   94)       ! 1. Laplacian of (tb,sb) * aht
(   95)       ! ----------------------------- 
(   96)       ! rotated harmonic operator applied to (tb,sb)
(   97)       ! and multiply by aht (output in (wk1,wk2) )
(   98) 
(   99)       CALL ldfght ( kt, tb, sb, wk1, wk2, 1 )
(  100) 
(  101) 
(  102)       ! Lateral boundary conditions on (wk1,wk2)   (unchanged sign)
(  103)       CALL lbc_lnk( wk1, 'T', 1. )   ;   CALL lbc_lnk( wk2, 'T', 1. )
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 3

(  105)       ! 2. Bilaplacian of (tb,sb)
(  106)       ! -------------------------
(  107)       ! rotated harmonic operator applied to (wk1,wk2)
(  108)       ! (output in (wk3,wk4) )
(  109) 
(  110)       CALL ldfght ( kt, wk1, wk2, wk3, wk4, 2 )
(  111) 
(  112) 
(  113)       ! 3. Update the tracer trends                    (j-slab :   2, jpj-1)
(  114)       ! ---------------------------
(  115)       !                                                ! ===============
(  116)       DO jj = 2, jpjm1                                 !  Vertical slab
(  117)          !                                             ! ===============
(  118)          DO jk = 1, jpkm1
(  119)             DO ji = 2, jpim1
(  120)                ! add it to the general tracer trends
(  121)                ta(ji,jj,jk) = ta(ji,jj,jk) + wk3(ji,jj,jk)
(  122)                sa(ji,jj,jk) = sa(ji,jj,jk) + wk4(ji,jj,jk)
(  123)             END DO
(  124)          END DO
(  125)          !                                             ! ===============
(  126)       END DO                                           !   End of slab
(  127)       !                                                ! ===============
(  128) 
(  129)       ! save the trends for diagnostic
(  130)       ! save the horizontal diffusive trends
(  131)       IF( l_trdtra )   THEN
(  132) 
(  133)          CALL trd_mod(wk3, wk4, jpttdldf, 'TRA', kt)
(  134)       ENDIF
(  135) 
(  136)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  137)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  138)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  139)       ENDIF
(  140) 
(  141)    END SUBROUTINE tra_ldf_bilapg
(  142) 
(  143) 
(  144)    SUBROUTINE ldfght ( kt, pt, ps, plt, pls, kaht )
(  145)       !!----------------------------------------------------------------------
(  146)       !!                  ***  ROUTINE ldfght  ***
(  147)       !!          
(  148)       !! ** Purpose :   Apply a geopotential harmonic operator to (pt,ps) and 
(  149)       !!      multiply it by the eddy diffusivity coefficient (if kaht=1).
(  150)       !!      Routine only used in s-coordinates (l_sco=T) with bilaplacian
(  151)       !!      operator (ln_traldf_bilap=T) acting along geopotential surfaces
(  152)       !!      (ln_traldf_hor).
(  153)       !!
(  154)       !! ** Method  :   The harmonic operator rotated along geopotential 
(  155)       !!      surfaces is applied to (pt,ps) using the slopes of geopotential
(  156)       !!      surfaces computed in inildf routine. The result is provided in
(  157)       !!      (plt,pls) arrays. It is computed in 2 steps:
(  158)       !!
(  159)       !!      First step: horizontal part of the operator. It is computed on
(  160)       !!      ==========  pt as follows (idem on ps)
(  161)       !!      horizontal fluxes :
(  162)       !!         zftu = e2u*e3u/e1u di[ pt ] - e2u*uslp dk[ mi(mk(pt)) ]






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 4

(  163)       !!         zftv = e1v*e3v/e2v dj[ pt ] - e1v*vslp dk[ mj(mk(pt)) ]
(  164)       !!      take the horizontal divergence of the fluxes (no divided by
(  165)       !!      the volume element :
(  166)       !!         plt  = di-1[ zftu ] +  dj-1[ zftv ]
(  167)       !!
(  168)       !!      Second step: vertical part of the operator. It is computed on
(  169)       !!      ===========  pt as follows (idem on ps)
(  170)       !!      vertical fluxes :
(  171)       !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pt ]
(  172)       !!              -     e2t     *       wslpi        di[ mi(mk(pt)) ]
(  173)       !!              -     e1t     *       wslpj        dj[ mj(mk(pt)) ]
(  174)       !!      take the vertical divergence of the fluxes add it to the hori-
(  175)       !!      zontal component, divide the result by the volume element and
(  176)       !!      if kaht=1, multiply by the eddy diffusivity coefficient:
(  177)       !!         plt  = aht / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  178)       !!      else:
(  179)       !!         plt  =  1  / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  180)       !!
(  181)       !! * Action :
(  182)       !!      'key_trdtra' defined: the trend is saved for diagnostics.
(  183)       !!
(  184)       !! History :
(  185)       !!   8.0  !  97-07  (G. Madec)  Original code
(  186)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  187)       !!----------------------------------------------------------------------
(  188)       !! * Arguments
(  189)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  190)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in  ) ::   &
(  191)          pt, ps           ! tracer fields (before t and s for 1st call
(  192)       !                   ! and laplacian of these fields for 2nd call.
(  193)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  194)          plt, pls         ! partial harmonic operator applied to
(  195)       !                   ! pt & ps components except
(  196)       !                   ! second order vertical derivative term)
(  197)       INTEGER, INTENT( in ) ::   &
(  198)          kaht             ! =1 multiply the laplacian by the eddy diffusivity coeff.
(  199)       !                   ! =2 no multiplication
(  200) 
(  201)       !! * Local declarations
(  202)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(  203)       REAL(wp) ::   &
(  204)          zabe1, zabe2, zmku, zmkv,     &  ! temporary scalars
(  205)          zbtr, ztah, zsah, ztav, zsav, &
(  206)          zcof0, zcof1, zcof2,          &
(  207)          zcof3, zcof4
(  208)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(  209)          zftu, zfsu,                   &  ! workspace
(  210)          zdkt, zdk1t,                  &
(  211)          zdks, zdk1s
(  212)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  213)          zftv, zfsv                       ! workspace (only v components for ptr)
(  214)       REAL(wp), DIMENSION(jpi,jpk) ::  &
(  215)          zftw, zfsw,                   &  ! workspace
(  216)          zdit, zdjt, zdj1t,            &
(  217)          zdis, zdjs, zdj1s
(  218)       !!----------------------------------------------------------------------
(  219) 
(  220)       !                               ! ********** !   ! ===============






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 5

(  221)       DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
(  222)          !                            ! ********** !   ! ===============
(  223) 
(  224)          ! I.1 Vertical gradient of pt and ps at level jk and jk+1
(  225)          ! -------------------------------------------------------
(  226)          !     surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  227) 
(  228)          zdk1t(:,:) = ( pt(:,:,jk) - pt(:,:,jk+1) ) * tmask(:,:,jk+1)
(  229)          zdk1s(:,:) = ( ps(:,:,jk) - ps(:,:,jk+1) ) * tmask(:,:,jk+1)
(  230) 
(  231)          IF( jk == 1 ) THEN
(  232)             zdkt(:,:) = zdk1t(:,:)
(  233)             zdks(:,:) = zdk1s(:,:)
(  234)          ELSE
(  235)             zdkt(:,:) = ( pt(:,:,jk-1) - pt(:,:,jk) ) * tmask(:,:,jk)
(  236)             zdks(:,:) = ( ps(:,:,jk-1) - ps(:,:,jk) ) * tmask(:,:,jk)
(  237)          ENDIF
(  238) 
(  239) 
(  240)          ! I.2 Horizontal fluxes
(  241)          ! ---------------------
(  242) 
(  243)          DO jj = 1, jpjm1
(  244)             DO ji = 1, jpim1
(  245)                zabe1 = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  246)                zabe2 = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  247) 
(  248)                zmku=1./MAX( tmask(ji+1,jj,jk  )+tmask(ji,jj,jk+1)   &
(  249)                            +tmask(ji+1,jj,jk+1)+tmask(ji,jj,jk  ),1. )
(  250)                zmkv=1./MAX( tmask(ji,jj+1,jk  )+tmask(ji,jj,jk+1)   &
(  251)                            +tmask(ji,jj+1,jk+1)+tmask(ji,jj,jk  ),1. )
(  252) 
(  253)                zcof1= -e2u(ji,jj) * uslp(ji,jj,jk) * zmku
(  254)                zcof2= -e1v(ji,jj) * vslp(ji,jj,jk) * zmkv
(  255) 
(  256)                zftu(ji,jj)= umask(ji,jj,jk) *   &
(  257)                   (  zabe1 *( pt(ji+1,jj,jk) - pt(ji,jj,jk) )   &
(  258)                    + zcof1 *( zdkt (ji+1,jj) + zdk1t(ji,jj)     &
(  259)                              +zdk1t(ji+1,jj) + zdkt (ji,jj) )  )
(  260) 
(  261)                zftv(ji,jj,jk)= vmask(ji,jj,jk) *   &
(  262)                   (  zabe2 *( pt(ji,jj+1,jk) - pt(ji,jj,jk) )   &
(  263)                    + zcof2 *( zdkt (ji,jj+1) + zdk1t(ji,jj)     &
(  264)                              +zdk1t(ji,jj+1) + zdkt (ji,jj) )  )
(  265) 
(  266)                zfsu(ji,jj)= umask(ji,jj,jk) *   &
(  267)                   (  zabe1 *( ps(ji+1,jj,jk) - ps(ji,jj,jk) )   &
(  268)                    + zcof1 *( zdks (ji+1,jj) + zdk1s(ji,jj)     &
(  269)                              +zdk1s(ji+1,jj) + zdks (ji,jj) )  )
(  270) 
(  271)                zfsv(ji,jj,jk)= vmask(ji,jj,jk) *   &
(  272)                   (  zabe2 *( ps(ji,jj+1,jk) - ps(ji,jj,jk) )   &
(  273)                    + zcof2 *( zdks (ji,jj+1) + zdk1s(ji,jj)     &
(  274)                              +zdk1s(ji,jj+1) + zdks (ji,jj) )  )
(  275)             END DO
(  276)          END DO
(  277) 
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 6

(  279)          ! I.3 Second derivative (divergence) (not divided by the volume)
(  280)          ! ---------------------
(  281) 
(  282)          DO jj = 2 , jpjm1
(  283)             DO ji = 2 , jpim1
(  284)                ztah = zftu(ji,jj) - zftu(ji-1,jj) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)
(  285)                zsah = zfsu(ji,jj) - zfsu(ji-1,jj) + zfsv(ji,jj,jk) - zfsv(ji,jj-1,jk)
(  286)                plt(ji,jj,jk) = ztah
(  287)                pls(ji,jj,jk) = zsah
(  288)             END DO
(  289)          END DO
(  290)          !                                             ! ===============
(  291)       END DO                                           !   End of slab
(  292)       !                                                ! ===============
(  293) 
(  294)       ! "zonal" mean diffusive heat and salt transport
(  295)       IF( ln_diaptr .AND. ( kaht == 2 ) .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  296)          pht_ldf(:) = ptr_vj( zftv(:,:,:) )
(  297)          pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
(  298)       ENDIF
(  299) 
(  300)       !                             ! ************ !   ! ===============
(  301)       DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
(  302)          !                          ! ************ !   ! ===============
(  303) 
(  304)          ! II.1 horizontal tracer gradient
(  305)          ! -------------------------------
(  306) 
(  307)          DO jk = 1, jpk
(  308)             DO ji = 1, jpim1
(  309)                zdit (ji,jk) = ( pt(ji+1,jj  ,jk) - pt(ji,jj  ,jk) ) * umask(ji,jj  ,jk)
(  310)                zdis (ji,jk) = ( ps(ji+1,jj  ,jk) - ps(ji,jj  ,jk) ) * umask(ji,jj  ,jk)
(  311)                zdjt (ji,jk) = ( pt(ji  ,jj+1,jk) - pt(ji,jj  ,jk) ) * vmask(ji,jj  ,jk)
(  312)                zdjs (ji,jk) = ( ps(ji  ,jj+1,jk) - ps(ji,jj  ,jk) ) * vmask(ji,jj  ,jk)
(  313)                zdj1t(ji,jk) = ( pt(ji  ,jj  ,jk) - pt(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  314)                zdj1s(ji,jk) = ( ps(ji  ,jj  ,jk) - ps(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  315)             END DO
(  316)          END DO
(  317) 
(  318) 
(  319)          ! II.2 Vertical fluxes
(  320)          ! --------------------
(  321) 
(  322)          ! Surface and bottom vertical fluxes set to zero
(  323)          zftw(:, 1 ) = 0.e0
(  324)          zfsw(:, 1 ) = 0.e0
(  325)          zftw(:,jpk) = 0.e0
(  326)          zfsw(:,jpk) = 0.e0
(  327) 
(  328)          ! interior (2=<jk=<jpk-1)
(  329)          DO jk = 2, jpkm1
(  330)             DO ji = 2, jpim1
(  331)                zcof0 = e1t(ji,jj) * e2t(ji,jj) / fse3w(ji,jj,jk)   &
(  332)                      * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)        &
(  333)                         + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  334) 
(  335)                zmku =1./MAX(  umask(ji  ,jj,jk-1)+umask(ji-1,jj,jk)   &
(  336)                              +umask(ji-1,jj,jk-1)+umask(ji  ,jj,jk), 1. )






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 7

(  337) 
(  338)                zmkv =1./MAX(  vmask(ji,jj  ,jk-1)+vmask(ji,jj-1,jk)   &
(  339)                              +vmask(ji,jj-1,jk-1)+vmask(ji,jj  ,jk), 1. )
(  340) 
(  341)                zcof3 = - e2t(ji,jj) * wslpi (ji,jj,jk) * zmku
(  342)                zcof4 = - e1t(ji,jj) * wslpj (ji,jj,jk) * zmkv
(  343) 
(  344)                zftw(ji,jk) = tmask(ji,jj,jk) *   &
(  345)                   (  zcof0 * ( pt  (ji,jj,jk-1) - pt  (ji,jj,jk) )   &
(  346)                    + zcof3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
(  347)                               +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
(  348)                    + zcof4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
(  349)                               +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )  )
(  350) 
(  351)                zfsw(ji,jk) = tmask(ji,jj,jk) *   &
(  352)                   (  zcof0 * ( ps  (ji,jj,jk-1) - ps  (ji,jj,jk) )   &
(  353)                    + zcof3 * ( zdis (ji  ,jk-1) + zdis (ji-1,jk)     &
(  354)                               +zdis (ji-1,jk-1) + zdis (ji  ,jk) )   &
(  355)                    + zcof4 * ( zdjs (ji  ,jk-1) + zdj1s(ji  ,jk)     &
(  356)                               +zdj1s(ji  ,jk-1) + zdjs (ji  ,jk) )  )
(  357)             END DO
(  358)          END DO
(  359) 
(  360) 
(  361)          ! II.3 Divergence of vertical fluxes added to the horizontal divergence
(  362)          ! ---------------------------------------------------------------------
(  363) 
(  364)          IF( kaht == 1 ) THEN
(  365)             ! multiply the laplacian by the eddy diffusivity coefficient
(  366)             DO jk = 1, jpkm1
(  367)                DO ji = 2, jpim1
(  368)                   ! eddy coef. divided by the volume element
(  369)                   zbtr = fsahtt(ji,jj,jk) / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  370)                   ! vertical divergence
(  371)                   ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  372)                   zsav = zfsw(ji,jk) - zfsw(ji,jk+1)
(  373)                   ! harmonic operator applied to (pt,ps) and multiply by aht
(  374)                   plt(ji,jj,jk) = ( plt(ji,jj,jk) + ztav ) * zbtr
(  375)                   pls(ji,jj,jk) = ( pls(ji,jj,jk) + zsav ) * zbtr
(  376)                END DO
(  377)             END DO
(  378)          ELSEIF( kaht == 2 ) THEN
(  379)             ! second call, no multiplication
(  380)             DO jk = 1, jpkm1
(  381)                DO ji = 2, jpim1
(  382)                   ! inverse of the volume element
(  383)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  384)                   ! vertical divergence
(  385)                   ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  386)                   zsav = zfsw(ji,jk) - zfsw(ji,jk+1)
(  387)                   ! harmonic operator applied to (pt,ps) 
(  388)                   plt(ji,jj,jk) = ( plt(ji,jj,jk) + ztav ) * zbtr
(  389)                   pls(ji,jj,jk) = ( pls(ji,jj,jk) + zsav ) * zbtr
(  390)                END DO
(  391)             END DO
(  392)          ELSE
(  393)             IF(lwp) WRITE(numout,*) ' ldfght: kaht= 1 or 2, here =', kaht
(  394)             IF(lwp) WRITE(numout,*) '         We stop'






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 8

(  395)             STOP 'ldfght'
(  396)          ENDIF
(  397)          !                                             ! ===============
(  398)       END DO                                           !   End of slab
(  399)       !                                                ! ===============
(  400)    END SUBROUTINE ldfght
(  401) 
(  402) #else  
(  403)    !!----------------------------------------------------------------------
(  404)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  405)    !!----------------------------------------------------------------------
(  406) CONTAINS
(  407)    SUBROUTINE tra_ldf_bilapg( kt )               ! Dummy routine
(  408) !      WRITE(*,*) 'tra_ldf_bilapg: You should not have seen this print! error?', kt
(  409)    END SUBROUTINE tra_ldf_bilapg
(  410) #endif
(  411) 
(  412)    !!==============================================================================
(  413) END MODULE traldf_bilapg













































PGF90 (Version     10.5)          09/23/2010  09:50:13      page 9

(    1) # 1 "traldf_bilapg.F90"
(    1) MODULE traldf_bilapg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traldf_bilapg  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_ldf_bilapg : update the tracer trend with the horizontal diffusion
(   11)    !!                    using an horizontal biharmonic operator in s-coordinate 
(   12)    !!   ldfght         :  ???
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables
(   17)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE ldfslp          ! iso-neutral slopes available
(   22)    USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(   23)    USE diaptr          ! poleward transport diagnostics 
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC tra_ldf_bilapg    ! routine called by step.F90
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
(   48) # 34 "traldf_bilapg.F90"
(   34) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 10

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 35 "traldf_bilapg.F90"
(   35) # 1 "./ldfeiv_substitute.h90"
(    1) # 40 "./ldfeiv_substitute.h90"
(   40) !   Default option :                         eiv: dummy variables
(   41) # 36 "traldf_bilapg.F90"
(   36) # 36
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_bilapg.F90,v 1.7 2005/09/02 15:45:33 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41)    
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE tra_ldf_bilapg( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                 ***  ROUTINE tra_ldf_bilapg  ***
(   47)       !!                    
(   48)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   49)       !!      trend and add it to the general trend of tracer equation.
(   50)       !!
(   51)       !! ** Method  :   The lateral diffusive trends is provided by a 4th order
(   52)       !!      operator rotated along geopotential surfaces. It is computed
(   53)       !!      using before fields (forward in time) and geopotential slopes
(   54)       !!      computed in routine inildf.
(   55)       !!         -1- compute the geopotential harmonic operator applied to
(   56)       !!      (tb,sb) and multiply it by the eddy diffusivity coefficient
(   57)       !!      (done by a call to ldfght routine, result in (wk1,wk2) arrays).
(   58)       !!      Applied the domain lateral boundary conditions by call to lbc_lnk
(   59)       !!         -2- compute the geopotential harmonic operator applied to
(   60)       !!      (wk1,wk2) by a second call to ldfght routine (result in (wk3,wk4)
(   61)       !!      arrays).
(   62)       !!         -3- Add this trend to the general trend (ta,sa):
(   63)       !!            (ta,sa) = (ta,sa) + (wk3,wk4)
(   64)       !!
(   65)       !! ** Action : - Update (ta,sa) arrays with the before geopotential 
(   66)       !!               biharmonic mixing trend.
(   67)       !!             - Save the trends  in (ztdta,ztdsa) ('key_trdtra')
(   68)       !!
(   69)       !! History :
(   70)       !!   8.0  !  97-07  (G. Madec)  Original code
(   71)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   72)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   73)       !!----------------------------------------------------------------------
(   74)       !! * Modules used
(   75)       USE oce           , wk1 => ua,  &  ! use ua as workspace
(   76)          &                wk2 => va      ! use va as workspace
(   77) 
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER ::   ji, jj, jk                 ! dummy loop indices






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 11

(   83)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   84)          wk3, wk4                ! work array used for rotated biharmonic
(   85)          !                       ! operator on tracers and/or momentum
(   86)       !!----------------------------------------------------------------------
(   87) 
(   88)       IF( kt == nit000 ) THEN
(   89)          IF(lwp) WRITE(numout,*)
(   90)          IF(lwp) WRITE(numout,*) 'tra_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
(   91)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(   92)       ENDIF
(   93) 
(   94)       ! 1. Laplacian of (tb,sb) * aht
(   95)       ! ----------------------------- 
(   96)       ! rotated harmonic operator applied to (tb,sb)
(   97)       ! and multiply by aht (output in (wk1,wk2) )
(   98) 
(   99)       CALL ldfght ( kt, tb, sb, wk1, wk2, 1 )
(  100) 
(  101) 
(  102)       ! Lateral boundary conditions on (wk1,wk2)   (unchanged sign)
(  103)       CALL lbc_lnk( wk1, 'T', 1. )   ;   CALL lbc_lnk( wk2, 'T', 1. )
(  104) 
(  105)       ! 2. Bilaplacian of (tb,sb)
(  106)       ! -------------------------
(  107)       ! rotated harmonic operator applied to (wk1,wk2)
(  108)       ! (output in (wk3,wk4) )
(  109) 
(  110)       CALL ldfght ( kt, wk1, wk2, wk3, wk4, 2 )
(  111) 
(  112) 
(  113)       ! 3. Update the tracer trends                    (j-slab :   2, jpj-1)
(  114)       ! ---------------------------
(  115)       !                                                ! ===============
(  116)       DO jj = 2, jpjm1                                 !  Vertical slab
(  117)          !                                             ! ===============
(  118)          DO jk = 1, jpkm1
(  119)             DO ji = 2, jpim1
(  120)                ! add it to the general tracer trends
(  121)                ta(ji,jj,jk) = ta(ji,jj,jk) + wk3(ji,jj,jk)
(  122)                sa(ji,jj,jk) = sa(ji,jj,jk) + wk4(ji,jj,jk)
(  123)             END DO
(  124)          END DO
(  125)          !                                             ! ===============
(  126)       END DO                                           !   End of slab
(  127)       !                                                ! ===============
(  128) 
(  129)       ! save the trends for diagnostic
(  130)       ! save the horizontal diffusive trends
(  131)       IF( l_trdtra )   THEN
(  132) 
(  133)          CALL trd_mod(wk3, wk4, jpttdldf, 'TRA', kt)
(  134)       ENDIF
(  135) 
(  136)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  137)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  138)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  139)       ENDIF
(  140) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 12

(  141)    END SUBROUTINE tra_ldf_bilapg
(  142) 
(  143) 
(  144)    SUBROUTINE ldfght ( kt, pt, ps, plt, pls, kaht )
(  145)       !!----------------------------------------------------------------------
(  146)       !!                  ***  ROUTINE ldfght  ***
(  147)       !!          
(  148)       !! ** Purpose :   Apply a geopotential harmonic operator to (pt,ps) and 
(  149)       !!      multiply it by the eddy diffusivity coefficient (if kaht=1).
(  150)       !!      Routine only used in s-coordinates (l_sco=T) with bilaplacian
(  151)       !!      operator (ln_traldf_bilap=T) acting along geopotential surfaces
(  152)       !!      (ln_traldf_hor).
(  153)       !!
(  154)       !! ** Method  :   The harmonic operator rotated along geopotential 
(  155)       !!      surfaces is applied to (pt,ps) using the slopes of geopotential
(  156)       !!      surfaces computed in inildf routine. The result is provided in
(  157)       !!      (plt,pls) arrays. It is computed in 2 steps:
(  158)       !!
(  159)       !!      First step: horizontal part of the operator. It is computed on
(  160)       !!      ==========  pt as follows (idem on ps)
(  161)       !!      horizontal fluxes :
(  162)       !!         zftu = e2u*e3u/e1u di[ pt ] - e2u*uslp dk[ mi(mk(pt)) ]
(  163)       !!         zftv = e1v*e3v/e2v dj[ pt ] - e1v*vslp dk[ mj(mk(pt)) ]
(  164)       !!      take the horizontal divergence of the fluxes (no divided by
(  165)       !!      the volume element :
(  166)       !!         plt  = di-1[ zftu ] +  dj-1[ zftv ]
(  167)       !!
(  168)       !!      Second step: vertical part of the operator. It is computed on
(  169)       !!      ===========  pt as follows (idem on ps)
(  170)       !!      vertical fluxes :
(  171)       !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pt ]
(  172)       !!              -     e2t     *       wslpi        di[ mi(mk(pt)) ]
(  173)       !!              -     e1t     *       wslpj        dj[ mj(mk(pt)) ]
(  174)       !!      take the vertical divergence of the fluxes add it to the hori-
(  175)       !!      zontal component, divide the result by the volume element and
(  176)       !!      if kaht=1, multiply by the eddy diffusivity coefficient:
(  177)       !!         plt  = aht / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  178)       !!      else:
(  179)       !!         plt  =  1  / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  180)       !!
(  181)       !! * Action :
(  182)       !!      'key_trdtra' defined: the trend is saved for diagnostics.
(  183)       !!
(  184)       !! History :
(  185)       !!   8.0  !  97-07  (G. Madec)  Original code
(  186)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  187)       !!----------------------------------------------------------------------
(  188)       !! * Arguments
(  189)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  190)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in  ) ::   &
(  191)          pt, ps           ! tracer fields (before t and s for 1st call
(  192)       !                   ! and laplacian of these fields for 2nd call.
(  193)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  194)          plt, pls         ! partial harmonic operator applied to
(  195)       !                   ! pt & ps components except
(  196)       !                   ! second order vertical derivative term)
(  197)       INTEGER, INTENT( in ) ::   &
(  198)          kaht             ! =1 multiply the laplacian by the eddy diffusivity coeff.






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 13

(  199)       !                   ! =2 no multiplication
(  200) 
(  201)       !! * Local declarations
(  202)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(  203)       REAL(wp) ::   &
(  204)          zabe1, zabe2, zmku, zmkv,     &  ! temporary scalars
(  205)          zbtr, ztah, zsah, ztav, zsav, &
(  206)          zcof0, zcof1, zcof2,          &
(  207)          zcof3, zcof4
(  208)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(  209)          zftu, zfsu,                   &  ! workspace
(  210)          zdkt, zdk1t,                  &
(  211)          zdks, zdk1s
(  212)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  213)          zftv, zfsv                       ! workspace (only v components for ptr)
(  214)       REAL(wp), DIMENSION(jpi,jpk) ::  &
(  215)          zftw, zfsw,                   &  ! workspace
(  216)          zdit, zdjt, zdj1t,            &
(  217)          zdis, zdjs, zdj1s
(  218)       !!----------------------------------------------------------------------
(  219) 
(  220)       !                               ! ********** !   ! ===============
(  221)       DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
(  222)          !                            ! ********** !   ! ===============
(  223) 
(  224)          ! I.1 Vertical gradient of pt and ps at level jk and jk+1
(  225)          ! -------------------------------------------------------
(  226)          !     surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  227) 
(  228)          zdk1t(:,:) = ( pt(:,:,jk) - pt(:,:,jk+1) ) * tmask(:,:,jk+1)
(  229)          zdk1s(:,:) = ( ps(:,:,jk) - ps(:,:,jk+1) ) * tmask(:,:,jk+1)
(  230) 
(  231)          IF( jk == 1 ) THEN
(  232)             zdkt(:,:) = zdk1t(:,:)
(  233)             zdks(:,:) = zdk1s(:,:)
(  234)          ELSE
(  235)             zdkt(:,:) = ( pt(:,:,jk-1) - pt(:,:,jk) ) * tmask(:,:,jk)
(  236)             zdks(:,:) = ( ps(:,:,jk-1) - ps(:,:,jk) ) * tmask(:,:,jk)
(  237)          ENDIF
(  238) 
(  239) 
(  240)          ! I.2 Horizontal fluxes
(  241)          ! ---------------------
(  242) 
(  243)          DO jj = 1, jpjm1
(  244)             DO ji = 1, jpim1
(  245)                zabe1 = e2u(ji,jj) * e3t(jk) / e1u(ji,jj)
(  246)                zabe2 = e1v(ji,jj) * e3t(jk) / e2v(ji,jj)
(  247) 
(  248)                zmku=1./MAX( tmask(ji+1,jj,jk  )+tmask(ji,jj,jk+1)   &
(  249)                            +tmask(ji+1,jj,jk+1)+tmask(ji,jj,jk  ),1. )
(  250)                zmkv=1./MAX( tmask(ji,jj+1,jk  )+tmask(ji,jj,jk+1)   &
(  251)                            +tmask(ji,jj+1,jk+1)+tmask(ji,jj,jk  ),1. )
(  252) 
(  253)                zcof1= -e2u(ji,jj) * uslp(ji,jj,jk) * zmku
(  254)                zcof2= -e1v(ji,jj) * vslp(ji,jj,jk) * zmkv
(  255) 
(  256)                zftu(ji,jj)= umask(ji,jj,jk) *   &






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 14

(  257)                   (  zabe1 *( pt(ji+1,jj,jk) - pt(ji,jj,jk) )   &
(  258)                    + zcof1 *( zdkt (ji+1,jj) + zdk1t(ji,jj)     &
(  259)                              +zdk1t(ji+1,jj) + zdkt (ji,jj) )  )
(  260) 
(  261)                zftv(ji,jj,jk)= vmask(ji,jj,jk) *   &
(  262)                   (  zabe2 *( pt(ji,jj+1,jk) - pt(ji,jj,jk) )   &
(  263)                    + zcof2 *( zdkt (ji,jj+1) + zdk1t(ji,jj)     &
(  264)                              +zdk1t(ji,jj+1) + zdkt (ji,jj) )  )
(  265) 
(  266)                zfsu(ji,jj)= umask(ji,jj,jk) *   &
(  267)                   (  zabe1 *( ps(ji+1,jj,jk) - ps(ji,jj,jk) )   &
(  268)                    + zcof1 *( zdks (ji+1,jj) + zdk1s(ji,jj)     &
(  269)                              +zdk1s(ji+1,jj) + zdks (ji,jj) )  )
(  270) 
(  271)                zfsv(ji,jj,jk)= vmask(ji,jj,jk) *   &
(  272)                   (  zabe2 *( ps(ji,jj+1,jk) - ps(ji,jj,jk) )   &
(  273)                    + zcof2 *( zdks (ji,jj+1) + zdk1s(ji,jj)     &
(  274)                              +zdk1s(ji,jj+1) + zdks (ji,jj) )  )
(  275)             END DO
(  276)          END DO
(  277) 
(  278) 
(  279)          ! I.3 Second derivative (divergence) (not divided by the volume)
(  280)          ! ---------------------
(  281) 
(  282)          DO jj = 2 , jpjm1
(  283)             DO ji = 2 , jpim1
(  284)                ztah = zftu(ji,jj) - zftu(ji-1,jj) + zftv(ji,jj,jk) - zftv(ji,jj-1,jk)
(  285)                zsah = zfsu(ji,jj) - zfsu(ji-1,jj) + zfsv(ji,jj,jk) - zfsv(ji,jj-1,jk)
(  286)                plt(ji,jj,jk) = ztah
(  287)                pls(ji,jj,jk) = zsah
(  288)             END DO
(  289)          END DO
(  290)          !                                             ! ===============
(  291)       END DO                                           !   End of slab
(  292)       !                                                ! ===============
(  293) 
(  294)       ! "zonal" mean diffusive heat and salt transport
(  295)       IF( ln_diaptr .AND. ( kaht == 2 ) .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  296)          pht_ldf(:) = ptr_vj( zftv(:,:,:) )
(  297)          pst_ldf(:) = ptr_vj( zfsv(:,:,:) )
(  298)       ENDIF
(  299) 
(  300)       !                             ! ************ !   ! ===============
(  301)       DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
(  302)          !                          ! ************ !   ! ===============
(  303) 
(  304)          ! II.1 horizontal tracer gradient
(  305)          ! -------------------------------
(  306) 
(  307)          DO jk = 1, jpk
(  308)             DO ji = 1, jpim1
(  309)                zdit (ji,jk) = ( pt(ji+1,jj  ,jk) - pt(ji,jj  ,jk) ) * umask(ji,jj  ,jk)
(  310)                zdis (ji,jk) = ( ps(ji+1,jj  ,jk) - ps(ji,jj  ,jk) ) * umask(ji,jj  ,jk)
(  311)                zdjt (ji,jk) = ( pt(ji  ,jj+1,jk) - pt(ji,jj  ,jk) ) * vmask(ji,jj  ,jk)
(  312)                zdjs (ji,jk) = ( ps(ji  ,jj+1,jk) - ps(ji,jj  ,jk) ) * vmask(ji,jj  ,jk)
(  313)                zdj1t(ji,jk) = ( pt(ji  ,jj  ,jk) - pt(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  314)                zdj1s(ji,jk) = ( ps(ji  ,jj  ,jk) - ps(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 15

(  315)             END DO
(  316)          END DO
(  317) 
(  318) 
(  319)          ! II.2 Vertical fluxes
(  320)          ! --------------------
(  321) 
(  322)          ! Surface and bottom vertical fluxes set to zero
(  323)          zftw(:, 1 ) = 0.e0
(  324)          zfsw(:, 1 ) = 0.e0
(  325)          zftw(:,jpk) = 0.e0
(  326)          zfsw(:,jpk) = 0.e0
(  327) 
(  328)          ! interior (2=<jk=<jpk-1)
(  329)          DO jk = 2, jpkm1
(  330)             DO ji = 2, jpim1
(  331)                zcof0 = e1t(ji,jj) * e2t(ji,jj) / e3w(jk)   &
(  332)                      * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)        &
(  333)                         + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  334) 
(  335)                zmku =1./MAX(  umask(ji  ,jj,jk-1)+umask(ji-1,jj,jk)   &
(  336)                              +umask(ji-1,jj,jk-1)+umask(ji  ,jj,jk), 1. )
(  337) 
(  338)                zmkv =1./MAX(  vmask(ji,jj  ,jk-1)+vmask(ji,jj-1,jk)   &
(  339)                              +vmask(ji,jj-1,jk-1)+vmask(ji,jj  ,jk), 1. )
(  340) 
(  341)                zcof3 = - e2t(ji,jj) * wslpi (ji,jj,jk) * zmku
(  342)                zcof4 = - e1t(ji,jj) * wslpj (ji,jj,jk) * zmkv
(  343) 
(  344)                zftw(ji,jk) = tmask(ji,jj,jk) *   &
(  345)                   (  zcof0 * ( pt  (ji,jj,jk-1) - pt  (ji,jj,jk) )   &
(  346)                    + zcof3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
(  347)                               +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
(  348)                    + zcof4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
(  349)                               +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )  )
(  350) 
(  351)                zfsw(ji,jk) = tmask(ji,jj,jk) *   &
(  352)                   (  zcof0 * ( ps  (ji,jj,jk-1) - ps  (ji,jj,jk) )   &
(  353)                    + zcof3 * ( zdis (ji  ,jk-1) + zdis (ji-1,jk)     &
(  354)                               +zdis (ji-1,jk-1) + zdis (ji  ,jk) )   &
(  355)                    + zcof4 * ( zdjs (ji  ,jk-1) + zdj1s(ji  ,jk)     &
(  356)                               +zdj1s(ji  ,jk-1) + zdjs (ji  ,jk) )  )
(  357)             END DO
(  358)          END DO
(  359) 
(  360) 
(  361)          ! II.3 Divergence of vertical fluxes added to the horizontal divergence
(  362)          ! ---------------------------------------------------------------------
(  363) 
(  364)          IF( kaht == 1 ) THEN
(  365)             ! multiply the laplacian by the eddy diffusivity coefficient
(  366)             DO jk = 1, jpkm1
(  367)                DO ji = 2, jpim1
(  368)                   ! eddy coef. divided by the volume element
(  369)                   zbtr = ahtt(ji,jj,jk) / ( e1t(ji,jj)*e2t(ji,jj)*e3t(jk) )
(  370)                   ! vertical divergence
(  371)                   ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  372)                   zsav = zfsw(ji,jk) - zfsw(ji,jk+1)






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 16

(  373)                   ! harmonic operator applied to (pt,ps) and multiply by aht
(  374)                   plt(ji,jj,jk) = ( plt(ji,jj,jk) + ztav ) * zbtr
(  375)                   pls(ji,jj,jk) = ( pls(ji,jj,jk) + zsav ) * zbtr
(  376)                END DO
(  377)             END DO
(  378)          ELSEIF( kaht == 2 ) THEN
(  379)             ! second call, no multiplication
(  380)             DO jk = 1, jpkm1
(  381)                DO ji = 2, jpim1
(  382)                   ! inverse of the volume element
(  383)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*e3t(jk) )
(  384)                   ! vertical divergence
(  385)                   ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  386)                   zsav = zfsw(ji,jk) - zfsw(ji,jk+1)
(  387)                   ! harmonic operator applied to (pt,ps) 
(  388)                   plt(ji,jj,jk) = ( plt(ji,jj,jk) + ztav ) * zbtr
(  389)                   pls(ji,jj,jk) = ( pls(ji,jj,jk) + zsav ) * zbtr
(  390)                END DO
(  391)             END DO
(  392)          ELSE
(  393)             IF(lwp) WRITE(numout,*) ' ldfght: kaht= 1 or 2, here =', kaht
(  394)             IF(lwp) WRITE(numout,*) '         We stop'
(  395)             STOP 'ldfght'
(  396)          ENDIF
(  397)          !                                             ! ===============
(  398)       END DO                                           !   End of slab
(  399)       !                                                ! ===============
(  400)    END SUBROUTINE ldfght
(  401) 
(  402) 
(  403) # 412
(  412)    !!==============================================================================
(  413) END MODULE traldf_bilapg
