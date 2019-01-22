


PGF90 (Version     10.5)          09/23/2010  09:49:19      page 1

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

Filename: trcldf_bilapg.F90

(    1) MODULE trcldf_bilapg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcldf_bilapg  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if key_passivetrc && defined key_ldfslp
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_ldf_bilapg : update the tracer trend with the horizontal diffusion
(   11)    !!                    using an horizontal biharmonic operator in s-coordinate 
(   12)    !!   ldfght         :  ???
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce_trc             ! ocean dynamics and tracers variables
(   16)    USE trc                 ! ocean passive tracers variables
(   17)    USE lbclnk              ! ocean lateral boundary condition (or mpp link)
(   18)    USE prtctl_trc          ! Print control for debbuging
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Routine accessibility
(   24)    PUBLIC trc_ldf_bilapg    ! routine called by step.F90
(   25) 
(   26)    !! * Substitutions
(   27) #  include "passivetrc_substitute.h90"
(   28)    !!----------------------------------------------------------------------
(   29)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_bilapg.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33)    
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE trc_ldf_bilapg( kt )
(   37)       !!----------------------------------------------------------------------
(   38)       !!                 ***  ROUTINE trc_ldf_bilapg  ***
(   39)       !!                    
(   40)       !! ** Purpose :   Compute the before horizontal passive tracer diffusive 
(   41)       !!      trend and add it to the general trend of tracer equation.
(   42)       !!
(   43)       !! ** Method  :   The lateral diffusive trends is provided by a 4th order
(   44)       !!      operator rotated along geopotential surfaces. It is computed
(   45)       !!      using before fields (forward in time) and geopotential slopes
(   46)       !!      computed in routine inildf.






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 2

(   47)       !!         -1- compute the geopotential harmonic operator applied to
(   48)       !!      trb and multiply it by the eddy diffusivity coefficient
(   49)       !!      (done by a call to ldfght routine, result in wk1 array).
(   50)       !!      Applied the domain lateral boundary conditions by call to lbc_lnk
(   51)       !!         -2- compute the geopotential harmonic operator applied to
(   52)       !!      wk1 by a second call to ldfght routine (result in wk3 array).
(   53)       !!         -3- Add this trend to the general trend (ta,sa):
(   54)       !!            tra = tra + wk3
(   55)       !!
(   56)       !! ** Action : - Update tra arrays with the before geopotential 
(   57)       !!               biharmonic mixing trend.
(   58)       !!             - Save the trends  in trtrd ('key_trc_diatrd')
(   59)       !!
(   60)       !! History :
(   61)       !!   8.0  !  97-07  (G. Madec)  Original code
(   62)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   63)       !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
(   64)       !!----------------------------------------------------------------------
(   65)       !! * Arguments
(   66)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   67) 
(   68)       !! * Local declarations
(   69)       INTEGER ::   ji, jj, jk,jn              ! dummy loop indices
(   70)       REAL(wp) ::  ztra                       ! workspace   
(   71)       REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   &
(   72)          wk1, wk2               ! work array used for rotated biharmonic
(   73)                                 ! operator on tracers and/or momentum
(   74)       CHARACTER (len=22) :: charout
(   75)       !!----------------------------------------------------------------------
(   76) 
(   77)       IF( kt == nittrc000 ) THEN
(   78)          IF(lwp) WRITE(numout,*)
(   79)          IF(lwp) WRITE(numout,*) 'trc_ldf_bilapg : horizontal biharmonic operator in s-coordinate'
(   80)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~'
(   81)       ENDIF
(   82) 
(   83)       ! 1. Laplacian of passive tracers trb * aht
(   84)       ! ----------------------------- 
(   85)       ! rotated harmonic operator applied to trb 
(   86)       ! and multiply by aht (output in wk1 )
(   87) 
(   88)       CALL ldfght ( trb, wk1, 1 )
(   89) 
(   90)       DO jn = 1, jptra
(   91)       ! Lateral boundary conditions on wk1   (unchanged sign)
(   92)          CALL lbc_lnk( wk1(:,:,:,jn) , 'T', 1. )
(   93)       END DO
(   94) 
(   95)       ! 2. Bilaplacian of trb
(   96)       ! -------------------------
(   97)       ! rotated harmonic operator applied to wk1 (output in wk2 )
(   98) 
(   99)       CALL ldfght ( wk1, wk2, 2 )
(  100) 
(  101) 
(  102)       DO jn = 1, jptra
(  103)          ! 3. Update the tracer trends                    (j-slab :   2, jpj-1)
(  104)          ! ---------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 3

(  105)          !                                                ! ===============
(  106)          DO jj = 2, jpjm1                                 !  Vertical slab
(  107)             !                                             ! ===============
(  108)             DO jk = 1, jpkm1
(  109)                DO ji = 2, jpim1
(  110)                   ! add it to the general tracer trends
(  111)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + wk2(ji,jj,jk,jn)
(  112) #if defined key_trc_diatrd
(  113)                   ! save the horizontal diffusive trends
(  114)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = wk2(ji,jj,jk,jn)
(  115) #endif
(  116)                END DO
(  117)             END DO
(  118)             !                                             ! ===============
(  119)          END DO                                           !   End of slab
(  120)          !                                                ! ===============
(  121) 
(  122)       END DO
(  123) 
(  124)      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  125)          WRITE(charout, FMT="('ldf - bilapg')")
(  126)          CALL prt_ctl_trc_info(charout)
(  127)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  128)       ENDIF
(  129) 
(  130)    END SUBROUTINE trc_ldf_bilapg
(  131) 
(  132) 
(  133)    SUBROUTINE ldfght ( pt, plt, kaht )
(  134)       !!----------------------------------------------------------------------
(  135)       !!                  ***  ROUTINE ldfght  ***
(  136)       !!          
(  137)       !! ** Purpose :   Apply a geopotential harmonic operator to pt and 
(  138)       !!      multiply it by the eddy diffusivity coefficient (if kaht=1).
(  139)       !!      Routine only used in s-coordinates (l_sco=T) with bilaplacian
(  140)       !!      operator (ln_traldf_bilap=T) acting along geopotential surfaces
(  141)       !!      (ln_traldf_hor).
(  142)       !!
(  143)       !! ** Method  :   The harmonic operator rotated along geopotential 
(  144)       !!      surfaces is applied to pt using the slopes of geopotential
(  145)       !!      surfaces computed in inildf routine. The result is provided in
(  146)       !!      plt arrays. It is computed in 2 steps:
(  147)       !!
(  148)       !!      First step: horizontal part of the operator. It is computed on
(  149)       !!      ==========  pt as follows (idem on ps)
(  150)       !!      horizontal fluxes :
(  151)       !!         zftu = e2u*e3u/e1u di[ pt ] - e2u*uslp dk[ mi(mk(pt)) ]
(  152)       !!         zftv = e1v*e3v/e2v dj[ pt ] - e1v*vslp dk[ mj(mk(pt)) ]
(  153)       !!      take the horizontal divergence of the fluxes (no divided by
(  154)       !!      the volume element :
(  155)       !!         plt  = di-1[ zftu ] +  dj-1[ zftv ]
(  156)       !!
(  157)       !!      Second step: vertical part of the operator. It is computed on
(  158)       !!      ===========  pt as follows (idem on ps)
(  159)       !!      vertical fluxes :
(  160)       !!         zftw = e1t*e2t/e3w * (wslpi^2+wslpj^2)  dk-1[ pt ]
(  161)       !!              -     e2t     *       wslpi        di[ mi(mk(pt)) ]
(  162)       !!              -     e1t     *       wslpj        dj[ mj(mk(pt)) ]






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 4

(  163)       !!      take the vertical divergence of the fluxes add it to the hori-
(  164)       !!      zontal component, divide the result by the volume element and
(  165)       !!      if kaht=1, multiply by the eddy diffusivity coefficient:
(  166)       !!         plt  = aht / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  167)       !!      else:
(  168)       !!         plt  =  1  / (e1t*e2t*e3t) { plt + dk[ zftw ] }
(  169)       !!
(  170)       !! * Action :
(  171)       !!      'key_trdtra' defined: the trend is saved for diagnostics.
(  172)       !!
(  173)       !! History :
(  174)       !!   8.0  !  97-07  (G. Madec)  Original code
(  175)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  176)       !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
(  177)       !!----------------------------------------------------------------------
(  178)       !! * Arguments
(  179)       REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( in  ) ::   &
(  180)          pt               ! tracer fields before for 1st call
(  181)       !                   ! and laplacian of these fields for 2nd call.
(  182)       REAL(wp), DIMENSION(jpi,jpj,jpk,jptra), INTENT( out ) ::   &
(  183)          plt              ! partial harmonic operator applied to
(  184)       !                   ! pt components except
(  185)       !                   ! second order vertical derivative term)
(  186)       INTEGER, INTENT( in ) ::   &
(  187)          kaht             ! =1 multiply the laplacian by the eddy diffusivity coeff.
(  188)       !                   ! =2 no multiplication
(  189) 
(  190)       !! * Local declarations
(  191)       INTEGER ::   ji, jj, jk,jn             ! dummy loop indices
(  192)       REAL(wp) ::   &
(  193)          zabe1, zabe2, zmku, zmkv,     &  ! temporary scalars
(  194)          zbtr, ztah, ztav, &
(  195)          zcof0, zcof1, zcof2,          &
(  196)          zcof3, zcof4
(  197) 
(  198)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  199)          zftu, zftv ,                   &  ! workspace
(  200)          zdkt, zdk1t
(  201) 
(  202)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  203)          zftw,                          &  ! workspace
(  204)          zdit, zdjt, zdj1t
(  205) 
(  206)       !!----------------------------------------------------------------------
(  207) 
(  208) 
(  209)       DO jn = 1, jptra
(  210)          !                               ! ********** !   ! ===============
(  211)          DO jk = 1, jpkm1                ! First step !   ! Horizontal slab
(  212)             !                            ! ********** !   ! ===============
(  213) 
(  214)             ! I.1 Vertical gradient of pt and ps at level jk and jk+1
(  215)             ! -------------------------------------------------------
(  216)             !     surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  217) 
(  218)             zdk1t(:,:) = ( pt(:,:,jk,jn) - pt(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
(  219) 
(  220)             IF( jk == 1 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 5

(  221)                zdkt(:,:) = zdk1t(:,:)
(  222)             ELSE
(  223)                zdkt(:,:) = ( pt(:,:,jk-1,jn) - pt(:,:,jk,jn) ) * tmask(:,:,jk)
(  224)             ENDIF
(  225) 
(  226) 
(  227)             ! I.2 Horizontal fluxes
(  228)             ! ---------------------
(  229) 
(  230)             DO jj = 1, jpjm1
(  231)                DO ji = 1, jpim1
(  232)                   zabe1 = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  233)                   zabe2 = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  234) 
(  235)                   zmku=1./MAX( tmask(ji+1,jj,jk  )+tmask(ji,jj,jk+1)   &
(  236)                               +tmask(ji+1,jj,jk+1)+tmask(ji,jj,jk  ),1. )
(  237)                   zmkv=1./MAX( tmask(ji,jj+1,jk  )+tmask(ji,jj,jk+1)   &
(  238)                               +tmask(ji,jj+1,jk+1)+tmask(ji,jj,jk  ),1. )
(  239) 
(  240)                   zcof1= -e2u(ji,jj) * uslp(ji,jj,jk) * zmku
(  241)                   zcof2= -e1v(ji,jj) * vslp(ji,jj,jk) * zmkv
(  242) 
(  243)                   zftu(ji,jj)= umask(ji,jj,jk) *   &
(  244)                      (  zabe1 *( pt(ji+1,jj,jk,jn) - pt(ji,jj,jk,jn) )   &
(  245)                       + zcof1 *( zdkt (ji+1,jj) + zdk1t(ji,jj)           &
(  246)                                 +zdk1t(ji+1,jj) + zdkt (ji,jj) )  )
(  247) 
(  248)                   zftv(ji,jj)= vmask(ji,jj,jk) *   &
(  249)                      (  zabe2 *( pt(ji,jj+1,jk,jn) - pt(ji,jj,jk,jn) )   &
(  250)                       + zcof2 *( zdkt (ji,jj+1) + zdk1t(ji,jj)           &
(  251)                                 +zdk1t(ji,jj+1) + zdkt (ji,jj) )  )
(  252)                END DO
(  253)             END DO
(  254) 
(  255) 
(  256)             ! I.3 Second derivative (divergence) (not divided by the volume)
(  257)             ! ---------------------
(  258) 
(  259)             DO jj = 2 , jpjm1
(  260)                DO ji = 2 , jpim1
(  261)                   ztah = zftu(ji,jj) - zftu(ji-1,jj) + zftv(ji,jj) - zftv(ji,jj-1)
(  262)                   plt(ji,jj,jk,jn) = ztah
(  263)                END DO
(  264)             END DO
(  265)             !                                             ! ===============
(  266)          END DO                                           !   End of slab
(  267)          !                                                ! ===============
(  268) 
(  269) 
(  270)          !                             ! ************ !   ! ===============
(  271)          DO jj = 2, jpjm1              !  Second step !   ! Horizontal slab
(  272)             !                          ! ************ !   ! ===============
(  273) 
(  274)             ! II.1 horizontal tracer gradient
(  275)             ! -------------------------------
(  276) 
(  277)             DO jk = 1, jpk
(  278)                DO ji = 1, jpim1






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 6

(  279)                   zdit (ji,jk) = ( pt(ji+1,jj  ,jk,jn) - pt(ji,jj  ,jk,jn) ) * umask(ji,jj  ,jk)
(  280)                   zdjt (ji,jk) = ( pt(ji  ,jj+1,jk,jn) - pt(ji,jj  ,jk,jn) ) * vmask(ji,jj  ,jk)
(  281)                   zdj1t(ji,jk) = ( pt(ji  ,jj  ,jk,jn) - pt(ji,jj-1,jk,jn) ) * vmask(ji,jj-1,jk)
(  282)                END DO
(  283)             END DO
(  284) 
(  285) 
(  286)             ! II.2 Vertical fluxes
(  287)             ! --------------------
(  288) 
(  289)             ! Surface and bottom vertical fluxes set to zero
(  290)             zftw(:, 1 ) = 0.e0
(  291)             zftw(:,jpk) = 0.e0
(  292) 
(  293)             ! interior (2=<jk=<jpk-1)
(  294)             DO jk = 2, jpkm1
(  295)                DO ji = 2, jpim1
(  296)                   zcof0 = e1t(ji,jj) * e2t(ji,jj) / fse3w(ji,jj,jk)   &
(  297)                      * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)        &
(  298)                         + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  299) 
(  300)                   zmku =1./MAX(  umask(ji  ,jj,jk-1)+umask(ji-1,jj,jk)   &
(  301)                                 +umask(ji-1,jj,jk-1)+umask(ji  ,jj,jk), 1. )
(  302) 
(  303)                   zmkv =1./MAX(  vmask(ji,jj  ,jk-1)+vmask(ji,jj-1,jk)   &
(  304)                                 +vmask(ji,jj-1,jk-1)+vmask(ji,jj  ,jk), 1. )
(  305) 
(  306)                   zcof3 = - e2t(ji,jj) * wslpi (ji,jj,jk) * zmku
(  307)                   zcof4 = - e1t(ji,jj) * wslpj (ji,jj,jk) * zmkv
(  308) 
(  309)                   zftw(ji,jk) = tmask(ji,jj,jk) *   &
(  310)                      (  zcof0 * ( pt  (ji,jj,jk-1,jn) - pt  (ji,jj,jk,jn) )   &
(  311)                       + zcof3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)           &
(  312)                                  +zdit (ji-1,jk-1) + zdit (ji  ,jk) )         &
(  313)                       + zcof4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)           &
(  314)                                  +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )  )
(  315)                END DO
(  316)             END DO
(  317) 
(  318) 
(  319)             ! II.3 Divergence of vertical fluxes added to the horizontal divergence
(  320)             ! ---------------------------------------------------------------------
(  321) 
(  322)             IF( kaht == 1 ) THEN
(  323)                ! multiply the laplacian by the eddy diffusivity coefficient
(  324)                DO jk = 1, jpkm1
(  325)                   DO ji = 2, jpim1
(  326)                      ! eddy coef. divided by the volume element
(  327)                      zbtr = fsahtrt(ji,jj,jk) / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  328)                      ! vertical divergence
(  329)                      ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  330)                      ! harmonic operator applied to (pt,ps) and multiply by aht
(  331)                      plt(ji,jj,jk,jn) = ( plt(ji,jj,jk,jn) + ztav ) * zbtr
(  332)                   END DO
(  333)                END DO
(  334)             ELSEIF( kaht == 2 ) THEN
(  335)                ! second call, no multiplication
(  336)                DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 7

(  337)                   DO ji = 2, jpim1
(  338)                      ! inverse of the volume element
(  339)                      zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  340)                      ! vertical divergence
(  341)                      ztav = zftw(ji,jk) - zftw(ji,jk+1)
(  342)                      ! harmonic operator applied to (pt,ps) 
(  343)                      plt(ji,jj,jk,jn) = ( plt(ji,jj,jk,jn) + ztav ) * zbtr
(  344)                   END DO
(  345)                END DO
(  346)             ELSE
(  347)                IF(lwp) WRITE(numout,*) ' ldfght: kaht= 1 or 2, here =', kaht
(  348)                IF(lwp) WRITE(numout,*) '         We stop'
(  349)                STOP 'ldfght'
(  350)             ENDIF
(  351)             !                                             ! ===============
(  352)          END DO                                           !   End of slab
(  353)          !                                                ! ===============
(  354) 
(  355)       END DO
(  356) 
(  357)    END SUBROUTINE ldfght
(  358) 
(  359) #else  
(  360)    !!----------------------------------------------------------------------
(  361)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  362)    !!----------------------------------------------------------------------
(  363) CONTAINS
(  364)    SUBROUTINE trc_ldf_bilapg( kt )               ! Dummy routine
(  365)       INTEGER, INTENT(in) :: kt
(  366) !      WRITE(*,*) 'trc_ldf_bilapg: You should not have seen this print! error?', kt
(  367)    END SUBROUTINE trc_ldf_bilapg
(  368) #endif
(  369) 
(  370)    !!==============================================================================
(  371) END MODULE trcldf_bilapg





























PGF90 (Version     10.5)          09/23/2010  09:49:19      page 8

(    1) # 1 "trcldf_bilapg.F90"
(    1) MODULE trcldf_bilapg
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcldf_bilapg  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 360
(  360)    !!----------------------------------------------------------------------
(  361)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  362)    !!----------------------------------------------------------------------
(  363) CONTAINS
(  364)    SUBROUTINE trc_ldf_bilapg( kt )               ! Dummy routine
(  365)       INTEGER, INTENT(in) :: kt
(  366) !      WRITE(*,*) 'trc_ldf_bilapg: You should not have seen this print! error?', kt
(  367)    END SUBROUTINE trc_ldf_bilapg
(  368) 
(  369) # 370
(  370)    !!==============================================================================
(  371) END MODULE trcldf_bilapg
