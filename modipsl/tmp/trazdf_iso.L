


PGF90 (Version     10.5)          09/23/2010  09:50:17      page 1

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

Filename: trazdf_iso.F90

(    1) MODULE trazdf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_iso  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_zdf_iso  : update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables 
(   17)    USE ldfslp          ! Make iso-neutral slopes available
(   18)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics 
(   20)    USE zdfddm          ! ocean vertical physics: double diffusion
(   21)    USE trdmod          ! ocean active tracers trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   25)    USE zdfkpp          ! KPP parameterisation
(   26)    USE prtctl          ! Print control
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Accessibility
(   32)    PUBLIC tra_zdf_iso    ! routine called by step.F90
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36) #  include "ldftra_substitute.h90"
(   37) #  include "ldfeiv_substitute.h90"
(   38) #  include "zdfddm_substitute.h90"
(   39)    !!----------------------------------------------------------------------
(   40)    !!----------------------------------------------------------------------
(   41)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   42)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso.F90,v 1.13 2005/09/22 15:43:10 opalod Exp $ 
(   43)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   44)    !!----------------------------------------------------------------------
(   45) CONTAINS
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 2

(   47)    SUBROUTINE tra_zdf_iso( kt )
(   48)       !!----------------------------------------------------------------------
(   49)       !!                  ***  ROUTINE tra_zdf_iso  ***
(   50)       !!
(   51)       !! ** Purpose :
(   52)       !!     Compute the trend due to the vertical tracer diffusion inclu-
(   53)       !!     ding the vertical component of lateral mixing (only for second
(   54)       !!     order operator, for fourth order it is already computed and
(   55)       !!     add to the general trend in traldf.F) and add it to the general
(   56)       !!     trend of the tracer equations.
(   57)       !!
(   58)       !! ** Method :
(   59)       !!         The vertical component of the lateral diffusive trends is
(   60)       !!      provided by a 2nd order operator rotated along neural or geopo-
(   61)       !!      tential surfaces to which an eddy induced advection can be added
(   62)       !!      It is computed using before fields (forward in time) and isopyc-
(   63)       !!      nal or geopotential slopes computed in routine ldfslp.
(   64)       !!
(   65)       !!      First part: vertical trends associated with the lateral mixing
(   66)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(   67)       !!      vertical fluxes associated with the rotated lateral mixing:
(   68)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
(   69)       !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
(   70)       !!      save avt coef. resulting from vertical physics alone in zavt:
(   71)       !!         zavt = avt
(   72)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(   73)       !!         avt = avt + wslpi^2+wslj^2
(   74)       !!      add vertical Eddy Induced advective fluxes ('lk_traldf_eiv=T):
(   75)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(   76)       !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
(   77)       !!      take the horizontal divergence of the fluxes:
(   78)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(   79)       !!      Add this trend to the general trend (ta,sa):
(   80)       !!         ta = ta + difft
(   81)       !!
(   82)       !!      Second part: vertical trend associated with the vertical physics
(   83)       !!      ===========  (including the vertical flux proportional to dk[t]
(   84)       !!                  associated with the lateral mixing, through the
(   85)       !!                  update of avt)
(   86)       !!      The vertical diffusion of tracers (t & s) is given by:
(   87)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(   88)       !!      It is computed using a backward time scheme, t=ta.
(   89)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   90)       !!      both tracers (bottom, applied through the masked field avt).
(   91)       !!      Add this trend to the general trend ta,sa :
(   92)       !!         ta = ta + dz( avt dz(t) )
(   93)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
(   94)       !!
(   95)       !!      Third part: recover avt resulting from the vertical physics
(   96)       !!      ==========  alone, for further diagnostics (for example to
(   97)       !!                  compute the turbocline depth in zdfmxl.F90).
(   98)       !!         avt = zavt
(   99)       !!         (avs = zavs if lk_zdfddm=T )
(  100)       !!
(  101)       !!      'key_trdtra' defined: trend saved for futher diagnostics.
(  102)       !!
(  103)       !!      macro-tasked on vertical slab (jj-loop)
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 3

(  105)       !! ** Action :
(  106)       !!         Update (ta,sa) arrays with the before vertical diffusion trend
(  107)       !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
(  108)       !!
(  109)       !! History :
(  110)       !!   7.0  !  91-11  (G. Madec)  Original code
(  111)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(  112)       !!        !  96-01  (G. Madec)  statement function for e3
(  113)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(  114)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(  115)       !!        !  00-08  (G. Madec)  double diffusive mixing
(  116)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  117)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  118)       !!   9.0  !  05-06  (C. Ethe )  non-local flux in KPP vertical mixing scheme
(  119)       !!---------------------------------------------------------------------
(  120)       !! * Modules used
(  121)       USE oce                ,   &
(  122) # if defined key_zdfddm
(  123)          zavs => va,  &  ! use va as workspace
(  124) # endif
(  125)          zavt => ua      ! use ua as workspace
(  126) 
(  127)       !! * Arguments
(  128)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  129) 
(  130)       !! * Local save
(  131)       REAL(wp), DIMENSION(jpk), SAVE ::  &
(  132)          z2dt
(  133) 
(  134)       !! * Local declarations
(  135)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(  136)       INTEGER ::   ikst, ikenm2, ikstp1       ! temporary integers
(  137) #if defined key_partial_steps
(  138)       INTEGER ::   iku, ikv, ikv1             ! temporary integers
(  139) #endif
(  140)       REAL(wp) ::   &
(  141)          zcoef0, zcoef3,        &  ! ???
(  142)          zcoef4, zavi,          &  ! ???
(  143)          zbtr, zmku, zmkv,      &  !
(  144)          ztav, zsav
(  145)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  146)          zwd, zws, zwi,         &  ! ???
(  147)          zwx, zwy, zwz, zwt        ! ???
(  148)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  149)          ztfw, zdit, zdjt, zdj1t,   &
(  150)          zsfw, zdis, zdjs, zdj1s
(  151)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  152)          ztavg, zsavg,                      & ! workspace arrays
(  153)          ztdta, ztdsa                         ! workspace arrays
(  154) #if defined key_traldf_eiv   ||   defined key_esopa
(  155)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  156)          ztfwg, zsfwg 
(  157)       REAL(wp) ::                       &
(  158)          zcoeg3,                        &
(  159)          zuwk, zvwk,                    &
(  160)          zuwki, zvwki
(  161) #endif
(  162)       !!---------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 4

(  163)       !!  OPA 8.5, LODYC-IPSL (2002)
(  164)       !!---------------------------------------------------------------------
(  165) 
(  166)       IF( kt == nit000 ) THEN
(  167)          IF(lwp) WRITE(numout,*)
(  168)          IF(lwp) WRITE(numout,*) 'tra_zdf_iso : vertical mixing (including isopycnal component)'
(  169)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(  170) #if defined key_diaeiv
(  171)          w_eiv(:,:,:) = 0.e0
(  172) #endif
(  173)       ENDIF
(  174) 
(  175)       ! 0. Local constant initialization
(  176)       ! --------------------------------
(  177)       ztavg(:,:,:) = 0.e0
(  178)       zsavg(:,:,:) = 0.e0
(  179) 
(  180)       ! time step = 2 rdttra ex
(  181)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  182)          z2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(  183)       ELSEIF( kt <= nit000 + 1) THEN
(  184)          z2dt(:) = 2. * rdttra(:)          ! leapfrog
(  185)       ENDIF
(  186) 
(  187)       ! Save ta and sa trends
(  188)       IF( l_trdtra )   THEN
(  189)          ztdta(:,:,:) = ta(:,:,:) 
(  190)          ztdsa(:,:,:) = sa(:,:,:) 
(  191)       ENDIF
(  192) 
(  193)       !                                                ! ===============
(  194)       DO jj = 2, jpjm1                                 !  Vertical slab
(  195)          !                                             ! ===============
(  196) 
(  197)          ! I. vertical trends associated with the lateral mixing
(  198)          ! =====================================================
(  199)          !  (excluding the vertical flux proportional to dk[t]
(  200) 
(  201) 
(  202)          ! I.1 horizontal tracer gradient
(  203)          ! ------------------------------
(  204) 
(  205)          DO jk = 1, jpkm1
(  206)             DO ji = 1, jpim1
(  207)                ! i-gradient of T and S at jj
(  208)                zdit (ji,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
(  209)                zdis (ji,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
(  210)                ! j-gradient of T and S at jj
(  211)                zdjt (ji,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  212)                zdjs (ji,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  213)                ! j-gradient of T and S at jj+1
(  214)                zdj1t(ji,jk) = ( tb(ji,jj,jk)-tb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  215)                zdj1s(ji,jk) = ( sb(ji,jj,jk)-sb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  216)             END DO
(  217)          END DO
(  218) #  if defined key_partial_steps
(  219)          ! partial steps correction at the bottom ocean level 
(  220)          DO ji = 1, jpim1






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 5

(  221)             ! last ocean level
(  222)             iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  223)             ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  224)             ikv1 = MIN( mbathy(ji,jj), mbathy(ji  ,jj-1) ) - 1
(  225)             ! i-gradient of T and S at jj
(  226)             zdit (ji,iku) = gtu(ji,jj)
(  227)             zdis (ji,iku) = gsu(ji,jj)
(  228)             ! j-gradient of T and S at jj
(  229)             zdjt (ji,ikv) = gtv(ji,jj) 
(  230)             zdjs (ji,ikv) = gsv(ji,jj) 
(  231)             ! j-gradient of T and S at jj+1
(  232)             zdj1t(ji,ikv1)= gtv(ji,jj-1)
(  233)             zdj1s(ji,ikv1)= gsv(ji,jj-1)
(  234)          END DO
(  235) #endif
(  236) 
(  237) 
(  238)          ! I.2 Vertical fluxes
(  239)          ! -------------------
(  240) 
(  241)          ! Surface and bottom vertical fluxes set to zero
(  242)          ztfw(:, 1 ) = 0.e0
(  243)          zsfw(:, 1 ) = 0.e0
(  244)          ztfw(:,jpk) = 0.e0
(  245)          zsfw(:,jpk) = 0.e0
(  246) #if defined key_traldf_eiv
(  247)          ztfwg(:,:, 1 ) = 0.e0
(  248)          zsfwg(:,:, 1 ) = 0.e0
(  249)          ztfwg(:,:,jpk) = 0.e0
(  250)          zsfwg(:,:,jpk) = 0.e0
(  251) #endif
(  252) 
(  253)          ! interior (2=<jk=<jpk-1)
(  254)          DO jk = 2, jpkm1
(  255)             DO ji = 2, jpim1
(  256)                zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  257) 
(  258)                zmku = 1./MAX( umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)   &
(  259)                   &          +umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1. )
(  260) 
(  261)                zmkv = 1./MAX( vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)   &
(  262)                   &          +vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1. )
(  263) 
(  264)                zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  265)                zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  266) 
(  267)                ztfw(ji,jk) = zcoef3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
(  268)                   &                    +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
(  269)                   &        + zcoef4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
(  270)                   &                    +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )
(  271) 
(  272)                zsfw(ji,jk) = zcoef3 * ( zdis (ji  ,jk-1) + zdis (ji-1,jk)     &
(  273)                   &                    +zdis (ji-1,jk-1) + zdis (ji  ,jk) )   &
(  274)                   &        + zcoef4 * ( zdjs (ji  ,jk-1) + zdj1s(ji  ,jk)     &
(  275)                   &                    +zdj1s(ji  ,jk-1) + zdjs (ji  ,jk) )
(  276)             END DO
(  277)          END DO
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 6

(  279) 
(  280)          ! I.3  update and save of avt (and avs if double diffusive mixing)
(  281)          ! ---------------------------
(  282) 
(  283)          DO jk = 2, jpkm1
(  284)             DO ji = 2, jpim1
(  285) 
(  286)                zavi = fsahtw(ji,jj,jk)*( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)   &
(  287)                   &                     +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )
(  288) 
(  289)                ! save avt in zavt to recover avt for mixed layer depth diag.
(  290)                zavt(ji,jj,jk) = avt(ji,jj,jk)
(  291)                ! add isopycnal vertical coeff. to avt
(  292)                avt(ji,jj,jk) = avt(ji,jj,jk) + zavi
(  293)                ! same procedure on avs if necessary
(  294) #if defined key_zdfddm
(  295)                ! save avs in zavs to recover avs in output files
(  296)                zavs(ji,jj,jk) = fsavs(ji,jj,jk)
(  297)                ! add isopycnal vertical coeff. to avs
(  298)                fsavs(ji,jj,jk) = fsavs(ji,jj,jk) + zavi
(  299) #endif
(  300)             END DO
(  301)          END DO
(  302) 
(  303) #if defined key_traldf_eiv
(  304)          !                              ! ---------------------------------------!
(  305)          !                              ! Eddy induced vertical advective fluxes !
(  306)          !                              ! ---------------------------------------!
(  307) #if defined key_traldf_c2d || defined key_traldf_c3d
(  308)             DO jk = 2, jpkm1
(  309)                DO ji = 2, jpim1
(  310)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  311)                      &  * fsaeiu(ji-1,jj,jk) * e2u(ji-1,jj)*umask(ji-1,jj,jk)
(  312)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  313)                      &  * fsaeiu(ji  ,jj,jk) * e2u(ji  ,jj)*umask(ji  ,jj,jk)
(  314)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  315)                      &  * fsaeiv(ji,jj-1,jk) * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
(  316)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  317)                      &  * fsaeiv(ji,jj  ,jk) * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
(  318) 
(  319)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
(  320)    
(  321)                   ztfwg(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
(  322)                   zsfwg(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )
(  323)    
(  324)                   ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jj,jk)
(  325)                   zsfw(ji,jk) = zsfw(ji,jk) + zsfwg(ji,jj,jk)
(  326) # if defined key_diaeiv
(  327)                   w_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  328) # endif
(  329)                END DO
(  330)             END DO
(  331) 
(  332) #else
(  333)             DO jk = 2, jpkm1
(  334)                DO ji = 2, jpim1
(  335)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  336)                      &  * e2u(ji-1,jj)*umask(ji-1,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 7

(  337)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  338)                      &  * e2u(ji  ,jj)*umask(ji  ,jj,jk)
(  339)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  340)                      &  * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
(  341)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  342)                      &  * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
(  343)    
(  344)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
(  345)                      &            * ( zuwk - zuwki + zvwk - zvwki )
(  346) 
(  347)                   ztfwg(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
(  348)                   zsfwg(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )
(  349) 
(  350)                   ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jj,jk)
(  351)                   zsfw(ji,jk) = zsfw(ji,jk) + zsfwg(ji,jj,jk)
(  352) # if defined key_diaeiv
(  353)                   w_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  354) # endif
(  355)                END DO
(  356)             END DO
(  357) #endif
(  358) 
(  359) #endif
(  360) 
(  361)          ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  362)          ! -------------------------------------------------------------------
(  363) 
(  364)          DO jk = 1, jpkm1
(  365)             DO ji = 2, jpim1
(  366)                zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  367)                ztav = (  ztfw(ji,jk) - ztfw(ji,jk+1)  ) * zbtr
(  368)                zsav = (  zsfw(ji,jk) - zsfw(ji,jk+1)  ) * zbtr
(  369)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
(  370)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
(  371)             END DO
(  372)          END DO
(  373)          !                                             ! ===============
(  374)       END DO                                           !   End of slab
(  375)       !                                                ! ===============
(  376) 
(  377)       ! save the trends for diagnostic
(  378)       !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
(  379)       IF( l_trdtra )   THEN
(  380) #   if defined key_traldf_eiv
(  381)          ! Compute the vertical Gent velocity trend
(  382)          !                                                ! ===============
(  383)          DO jj = 2, jpjm1                                 !  Vertical slab
(  384)             !                                             ! ===============
(  385)             DO jk = 1, jpkm1
(  386)                DO ji = 2, jpim1
(  387)                   zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  388)                   ztavg(ji,jj,jk) = ( ztfwg(ji,jj,jk) - ztfwg(ji,jj,jk+1) ) * zbtr
(  389)                   zsavg(ji,jj,jk) = ( zsfwg(ji,jj,jk) - zsfwg(ji,jj,jk+1) ) * zbtr
(  390)                END DO
(  391)             END DO
(  392)             !                                             ! ===============
(  393)          END DO                                           !   End of slab
(  394)          !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 8

(  395) 
(  396)          CALL trd_mod(ztavg, zsavg, jpttddoe, 'TRA', kt)
(  397) #   endif
(  398)          ! Recompute the divergence of vertical fluxes ztav & zsav trends 
(  399)          ! computed at step 1.5 above in making the difference between the new 
(  400)          ! trend ta()/sa() and the previous one ztdta()/ztdsa() and substract 
(  401)          ! the vertical Gent velocity trend ztavg()/zsavg() (zero if not used)
(  402)          ztavg(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztavg(:,:,:) 
(  403)          zsavg(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsavg(:,:,:) 
(  404) 
(  405)          ! Save the new ta and sa trends
(  406)          ztdta(:,:,:) = ta(:,:,:) 
(  407)          ztdsa(:,:,:) = sa(:,:,:) 
(  408)       ENDIF
(  409) 
(  410)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  411)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
(  412)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  413)       ENDIF
(  414) 
(  415)       !                                                ! ===============
(  416)       DO jj = 2, jpjm1                                 !  Vertical slab
(  417)          !                                             ! ===============
(  418) 
(  419)          ! II. Vertical trend associated with the vertical physics
(  420)          ! =======================================================
(  421)          !     (including the vertical flux proportional to dk[t] associated
(  422)          !      with the lateral mixing, through the avt update)
(  423)          !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  424) 
(  425) 
(  426)          ! II.0 Matrix construction
(  427)          ! ------------------------
(  428) 
(  429)          ! Diagonal, inferior, superior
(  430)          ! (including the bottom boundary condition via avt masked)
(  431)          DO jk = 1, jpkm1
(  432)             DO ji = 2, jpim1
(  433)                zwi(ji,jk) = - z2dt(jk) * avt(ji,jj,jk  )   &
(  434)                                        / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  435)                zws(ji,jk) = - z2dt(jk) * avt(ji,jj,jk+1)   &
(  436)                                        / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  437)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  438)             END DO
(  439)          END DO
(  440) 
(  441)          ! Surface boudary conditions
(  442)          DO ji = 2, jpim1
(  443)             zwi(ji,1) = 0.e0
(  444)             zwd(ji,1) = 1. - zws(ji,1)
(  445)          END DO
(  446) 
(  447) 
(  448)          ! II.1. Vertical diffusion on t
(  449)          ! ---------------------------
(  450) 
(  451)          ! Second member construction
(  452) #if defined key_zdfkpp






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 9

(  453)          ! add non-local temperature flux ( in convective case only)
(  454)          DO jk = 1, jpkm1
(  455)             DO ji = 2, jpim1
(  456)                ! zrhs=right hand side 
(  457)                zwy(ji,jk) = tb(ji,jj,jk) + z2dt(jk) * ta(ji,jj,jk) &
(  458)                   &  - z2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
(  459)                   &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
(  460)             END DO
(  461)          END DO
(  462) #else
(  463)          DO jk = 1, jpkm1
(  464)             DO ji = 2, jpim1
(  465)                zwy(ji,jk) = tb(ji,jj,jk) + z2dt(jk) * ta(ji,jj,jk)
(  466)             END DO
(  467)          END DO
(  468) #endif
(  469) 
(  470)          ! Matrix inversion from the first level
(  471)          ikst = 1
(  472) #   include "zdf.matrixsolver.h90"
(  473) 
(  474)          ! Save the masked temperature after in ta
(  475)          ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
(  476)          !                  it will not be done in tranxt)
(  477)          DO jk = 1, jpkm1
(  478)             DO ji = 2, jpim1
(  479)                ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
(  480)             END DO
(  481)          END DO
(  482) 
(  483) 
(  484)          ! II.2 Vertical diffusion on s
(  485)          ! ---------------------------
(  486) 
(  487) #if defined key_zdfddm
(  488)          ! Rebuild the Matrix as avt /= avs
(  489) 
(  490)          ! Diagonal, inferior, superior
(  491)          ! (including the bottom boundary condition via avs masked)
(  492)          DO jk = 1, jpkm1
(  493)             DO ji = 2, jpim1
(  494)                zwi(ji,jk) = - z2dt(jk) * fsavs(ji,jj,jk  )   &
(  495)                                        /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  496)                zws(ji,jk) = - z2dt(jk) * fsavs(ji,jj,jk+1)   &
(  497)                                        /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  498)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  499)             END DO
(  500)          END DO
(  501) 
(  502)          ! Surface boudary conditions
(  503)          DO ji = 2, jpim1
(  504)             zwi(ji,1) = 0.e0
(  505)             zwd(ji,1) = 1. - zws(ji,1)
(  506)          END DO
(  507) #endif
(  508)          ! Second member construction
(  509) 
(  510) #if defined key_zdfkpp






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 10

(  511)          ! add non-local temperature flux ( in convective case only)
(  512)          DO jk = 1, jpkm1
(  513)             DO ji = 2, jpim1
(  514)                ! zrhs=right hand side 
(  515)                zwy(ji,jk) = sb(ji,jj,jk) + z2dt(jk) * sa(ji,jj,jk) &
(  516)                   &  - z2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
(  517)                   &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
(  518)             END DO
(  519)          END DO
(  520) #else
(  521)          DO jk = 1, jpkm1
(  522)             DO ji = 2, jpim1
(  523)                zwy(ji,jk) = sb(ji,jj,jk) + z2dt(jk) * sa(ji,jj,jk)
(  524)             END DO
(  525)          END DO
(  526) #endif
(  527) 
(  528)          ! Matrix inversion from the first level
(  529)          ikst = 1
(  530) #   include "zdf.matrixsolver.h90"
(  531) 
(  532)          ! Save the masked salinity after in sa
(  533)          ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
(  534)          !                  it will not be done in tranxt)
(  535)          DO jk = 1, jpkm1
(  536)             DO ji = 2, jpim1
(  537)                sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
(  538)             END DO
(  539)          END DO
(  540) 
(  541) 
(  542)          ! III. recover the avt (avs) resulting from vertical physics only
(  543)          ! ===============================================================
(  544) 
(  545)          DO jk = 2, jpkm1
(  546)             DO ji = 2, jpim1
(  547)                avt(ji,jj,jk) = zavt(ji,jj,jk)
(  548) #if defined key_zdfddm
(  549)                fsavs(ji,jj,jk) = zavs(ji,jj,jk)
(  550) #endif
(  551)             END DO
(  552)          END DO
(  553) 
(  554)          !                                             ! ===============
(  555)       END DO                                           !   End of slab
(  556)       !                                                ! ===============
(  557) 
(  558)       ! save the trends for diagnostic
(  559)       ! compute the vertical diffusive trends in substracting the previous 
(  560)       ! trends ztdta()/ztdsa() to the new one computed via dT/dt or dS/dt 
(  561)       ! i.e. with the new temperature/salinity ta/sa computed above
(  562)       IF( l_trdtra )   THEN
(  563)          IF( ln_traldf_iso)   THEN
(  564)             DO jk = 1, jpkm1
(  565)                ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk) + ztavg(:,:,jk) 
(  566)                ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk) + zsavg(:,:,jk) 
(  567)             END DO
(  568)          ELSE






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 11

(  569)             DO jk = 1, jpkm1
(  570)                ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk)                             
(  571)                ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk)                             
(  572)             END DO
(  573)          ENDIF
(  574) 
(  575)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  576)       ENDIF
(  577) 
(  578)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  579)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
(  580)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)
(  581)       ENDIF
(  582) 
(  583)    END SUBROUTINE tra_zdf_iso
(  584) 
(  585) #else
(  586)    !!----------------------------------------------------------------------
(  587)    !!   Dummy module               NO rotation of the lateral mixing tensor
(  588)    !!----------------------------------------------------------------------
(  589) CONTAINS
(  590)    SUBROUTINE tra_zdf_iso( kt )              ! empty routine
(  591) !      WRITE(*,*) 'tra_zdf_iso: You should not have seen this print! error?', kt
(  592)    END SUBROUTINE tra_zdf_iso
(  593) #endif
(  594) 
(  595)    !!==============================================================================
(  596) END MODULE trazdf_iso




































PGF90 (Version     10.5)          09/23/2010  09:50:17      page 12

(    1) # 1 "trazdf_iso.F90"
(    1) MODULE trazdf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_iso  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_zdf_iso  : update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers variables
(   16)    USE dom_oce         ! ocean space and time domain variables 
(   17)    USE ldfslp          ! Make iso-neutral slopes available
(   18)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   19)    USE zdf_oce         ! ocean vertical physics 
(   20)    USE zdfddm          ! ocean vertical physics: double diffusion
(   21)    USE trdmod          ! ocean active tracers trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   25)    USE zdfkpp          ! KPP parameterisation
(   26)    USE prtctl          ! Print control
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Accessibility
(   32)    PUBLIC tra_zdf_iso    ! routine called by step.F90
(   33) 
(   34)    !! * Substitutions
(   35) # 1 "./domzgr_substitute.h90"
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
(   48) # 36 "trazdf_iso.F90"
(   36) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 13

(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 37 "trazdf_iso.F90"
(   37) # 1 "./ldfeiv_substitute.h90"
(    1) # 40 "./ldfeiv_substitute.h90"
(   40) !   Default option :                         eiv: dummy variables
(   41) # 38 "trazdf_iso.F90"
(   38) # 1 "./zdfddm_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** zdfddm_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm_substitute.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 16
(   16) !   Defautl option :                     avs = avt
(   17) # 39 "trazdf_iso.F90"
(   39) # 39
(   39)    !!----------------------------------------------------------------------
(   40)    !!----------------------------------------------------------------------
(   41)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   42)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso.F90,v 1.13 2005/09/22 15:43:10 opalod Exp $ 
(   43)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   44)    !!----------------------------------------------------------------------
(   45) CONTAINS
(   46) 
(   47)    SUBROUTINE tra_zdf_iso( kt )
(   48)       !!----------------------------------------------------------------------
(   49)       !!                  ***  ROUTINE tra_zdf_iso  ***
(   50)       !!
(   51)       !! ** Purpose :
(   52)       !!     Compute the trend due to the vertical tracer diffusion inclu-
(   53)       !!     ding the vertical component of lateral mixing (only for second
(   54)       !!     order operator, for fourth order it is already computed and
(   55)       !!     add to the general trend in traldf.F) and add it to the general
(   56)       !!     trend of the tracer equations.
(   57)       !!
(   58)       !! ** Method :
(   59)       !!         The vertical component of the lateral diffusive trends is
(   60)       !!      provided by a 2nd order operator rotated along neural or geopo-
(   61)       !!      tential surfaces to which an eddy induced advection can be added
(   62)       !!      It is computed using before fields (forward in time) and isopyc-
(   63)       !!      nal or geopotential slopes computed in routine ldfslp.
(   64)       !!
(   65)       !!      First part: vertical trends associated with the lateral mixing
(   66)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(   67)       !!      vertical fluxes associated with the rotated lateral mixing:
(   68)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 14

(   69)       !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
(   70)       !!      save avt coef. resulting from vertical physics alone in zavt:
(   71)       !!         zavt = avt
(   72)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(   73)       !!         avt = avt + wslpi^2+wslj^2
(   74)       !!      add vertical Eddy Induced advective fluxes ('lk_traldf_eiv=T):
(   75)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(   76)       !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
(   77)       !!      take the horizontal divergence of the fluxes:
(   78)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(   79)       !!      Add this trend to the general trend (ta,sa):
(   80)       !!         ta = ta + difft
(   81)       !!
(   82)       !!      Second part: vertical trend associated with the vertical physics
(   83)       !!      ===========  (including the vertical flux proportional to dk[t]
(   84)       !!                  associated with the lateral mixing, through the
(   85)       !!                  update of avt)
(   86)       !!      The vertical diffusion of tracers (t & s) is given by:
(   87)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(   88)       !!      It is computed using a backward time scheme, t=ta.
(   89)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   90)       !!      both tracers (bottom, applied through the masked field avt).
(   91)       !!      Add this trend to the general trend ta,sa :
(   92)       !!         ta = ta + dz( avt dz(t) )
(   93)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
(   94)       !!
(   95)       !!      Third part: recover avt resulting from the vertical physics
(   96)       !!      ==========  alone, for further diagnostics (for example to
(   97)       !!                  compute the turbocline depth in zdfmxl.F90).
(   98)       !!         avt = zavt
(   99)       !!         (avs = zavs if lk_zdfddm=T )
(  100)       !!
(  101)       !!      'key_trdtra' defined: trend saved for futher diagnostics.
(  102)       !!
(  103)       !!      macro-tasked on vertical slab (jj-loop)
(  104)       !!
(  105)       !! ** Action :
(  106)       !!         Update (ta,sa) arrays with the before vertical diffusion trend
(  107)       !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
(  108)       !!
(  109)       !! History :
(  110)       !!   7.0  !  91-11  (G. Madec)  Original code
(  111)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(  112)       !!        !  96-01  (G. Madec)  statement function for e3
(  113)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(  114)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(  115)       !!        !  00-08  (G. Madec)  double diffusive mixing
(  116)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  117)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  118)       !!   9.0  !  05-06  (C. Ethe )  non-local flux in KPP vertical mixing scheme
(  119)       !!---------------------------------------------------------------------
(  120)       !! * Modules used
(  121)       USE oce                ,   &
(  125)          zavt => ua      ! use ua as workspace
(  126) 
(  127)       !! * Arguments
(  128)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  129) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 15

(  130)       !! * Local save
(  131)       REAL(wp), DIMENSION(jpk), SAVE ::  &
(  132)          z2dt
(  133) 
(  134)       !! * Local declarations
(  135)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(  136)       INTEGER ::   ikst, ikenm2, ikstp1       ! temporary integers
(  137) # 140
(  140)       REAL(wp) ::   &
(  141)          zcoef0, zcoef3,        &  ! ???
(  142)          zcoef4, zavi,          &  ! ???
(  143)          zbtr, zmku, zmkv,      &  !
(  144)          ztav, zsav
(  145)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  146)          zwd, zws, zwi,         &  ! ???
(  147)          zwx, zwy, zwz, zwt        ! ???
(  148)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  149)          ztfw, zdit, zdjt, zdj1t,   &
(  150)          zsfw, zdis, zdjs, zdj1s
(  151)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  152)          ztavg, zsavg,                      & ! workspace arrays
(  153)          ztdta, ztdsa                         ! workspace arrays
(  154) # 162
(  162)       !!---------------------------------------------------------------------
(  163)       !!  OPA 8.5, LODYC-IPSL (2002)
(  164)       !!---------------------------------------------------------------------
(  165) 
(  166)       IF( kt == nit000 ) THEN
(  167)          IF(lwp) WRITE(numout,*)
(  168)          IF(lwp) WRITE(numout,*) 'tra_zdf_iso : vertical mixing (including isopycnal component)'
(  169)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(  170) # 173
(  173)       ENDIF
(  174) 
(  175)       ! 0. Local constant initialization
(  176)       ! --------------------------------
(  177)       ztavg(:,:,:) = 0.e0
(  178)       zsavg(:,:,:) = 0.e0
(  179) 
(  180)       ! time step = 2 rdttra ex
(  181)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  182)          z2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(  183)       ELSEIF( kt <= nit000 + 1) THEN
(  184)          z2dt(:) = 2. * rdttra(:)          ! leapfrog
(  185)       ENDIF
(  186) 
(  187)       ! Save ta and sa trends
(  188)       IF( l_trdtra )   THEN
(  189)          ztdta(:,:,:) = ta(:,:,:) 
(  190)          ztdsa(:,:,:) = sa(:,:,:) 
(  191)       ENDIF
(  192) 
(  193)       !                                                ! ===============
(  194)       DO jj = 2, jpjm1                                 !  Vertical slab
(  195)          !                                             ! ===============
(  196) 
(  197)          ! I. vertical trends associated with the lateral mixing
(  198)          ! =====================================================






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 16

(  199)          !  (excluding the vertical flux proportional to dk[t]
(  200) 
(  201) 
(  202)          ! I.1 horizontal tracer gradient
(  203)          ! ------------------------------
(  204) 
(  205)          DO jk = 1, jpkm1
(  206)             DO ji = 1, jpim1
(  207)                ! i-gradient of T and S at jj
(  208)                zdit (ji,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
(  209)                zdis (ji,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
(  210)                ! j-gradient of T and S at jj
(  211)                zdjt (ji,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  212)                zdjs (ji,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  213)                ! j-gradient of T and S at jj+1
(  214)                zdj1t(ji,jk) = ( tb(ji,jj,jk)-tb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  215)                zdj1s(ji,jk) = ( sb(ji,jj,jk)-sb(ji,jj-1,jk) ) * vmask(ji,jj-1,jk)
(  216)             END DO
(  217)          END DO
(  218) 
(  219) 
(  220) # 238
(  238)          ! I.2 Vertical fluxes
(  239)          ! -------------------
(  240) 
(  241)          ! Surface and bottom vertical fluxes set to zero
(  242)          ztfw(:, 1 ) = 0.e0
(  243)          zsfw(:, 1 ) = 0.e0
(  244)          ztfw(:,jpk) = 0.e0
(  245)          zsfw(:,jpk) = 0.e0
(  246) 
(  247) # 253
(  253)          ! interior (2=<jk=<jpk-1)
(  254)          DO jk = 2, jpkm1
(  255)             DO ji = 2, jpim1
(  256)                zcoef0 = - ahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  257) 
(  258)                zmku = 1./MAX( umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)   &
(  259)                   &          +umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1. )
(  260) 
(  261)                zmkv = 1./MAX( vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)   &
(  262)                   &          +vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1. )
(  263) 
(  264)                zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  265)                zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  266) 
(  267)                ztfw(ji,jk) = zcoef3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
(  268)                   &                    +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
(  269)                   &        + zcoef4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
(  270)                   &                    +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )
(  271) 
(  272)                zsfw(ji,jk) = zcoef3 * ( zdis (ji  ,jk-1) + zdis (ji-1,jk)     &
(  273)                   &                    +zdis (ji-1,jk-1) + zdis (ji  ,jk) )   &
(  274)                   &        + zcoef4 * ( zdjs (ji  ,jk-1) + zdj1s(ji  ,jk)     &
(  275)                   &                    +zdj1s(ji  ,jk-1) + zdjs (ji  ,jk) )
(  276)             END DO
(  277)          END DO
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 17

(  279) 
(  280)          ! I.3  update and save of avt (and avs if double diffusive mixing)
(  281)          ! ---------------------------
(  282) 
(  283)          DO jk = 2, jpkm1
(  284)             DO ji = 2, jpim1
(  285) 
(  286)                zavi = ahtw(ji,jj,jk)*( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)   &
(  287)                   &                     +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )
(  288) 
(  289)                ! save avt in zavt to recover avt for mixed layer depth diag.
(  290)                zavt(ji,jj,jk) = avt(ji,jj,jk)
(  291)                ! add isopycnal vertical coeff. to avt
(  292)                avt(ji,jj,jk) = avt(ji,jj,jk) + zavi
(  293)                ! same procedure on avs if necessary
(  294) # 300
(  300)             END DO
(  301)          END DO
(  302) 
(  303) 
(  304) # 361
(  361)          ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  362)          ! -------------------------------------------------------------------
(  363) 
(  364)          DO jk = 1, jpkm1
(  365)             DO ji = 2, jpim1
(  366)                zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*e3t(jk) )
(  367)                ztav = (  ztfw(ji,jk) - ztfw(ji,jk+1)  ) * zbtr
(  368)                zsav = (  zsfw(ji,jk) - zsfw(ji,jk+1)  ) * zbtr
(  369)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
(  370)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
(  371)             END DO
(  372)          END DO
(  373)          !                                             ! ===============
(  374)       END DO                                           !   End of slab
(  375)       !                                                ! ===============
(  376) 
(  377)       ! save the trends for diagnostic
(  378)       !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
(  379)       IF( l_trdtra )   THEN
(  380) # 398
(  398)          ! Recompute the divergence of vertical fluxes ztav & zsav trends 
(  399)          ! computed at step 1.5 above in making the difference between the new 
(  400)          ! trend ta()/sa() and the previous one ztdta()/ztdsa() and substract 
(  401)          ! the vertical Gent velocity trend ztavg()/zsavg() (zero if not used)
(  402)          ztavg(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - ztavg(:,:,:) 
(  403)          zsavg(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - zsavg(:,:,:) 
(  404) 
(  405)          ! Save the new ta and sa trends
(  406)          ztdta(:,:,:) = ta(:,:,:) 
(  407)          ztdsa(:,:,:) = sa(:,:,:) 
(  408)       ENDIF
(  409) 
(  410)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  411)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
(  412)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  413)       ENDIF
(  414) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 18

(  415)       !                                                ! ===============
(  416)       DO jj = 2, jpjm1                                 !  Vertical slab
(  417)          !                                             ! ===============
(  418) 
(  419)          ! II. Vertical trend associated with the vertical physics
(  420)          ! =======================================================
(  421)          !     (including the vertical flux proportional to dk[t] associated
(  422)          !      with the lateral mixing, through the avt update)
(  423)          !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  424) 
(  425) 
(  426)          ! II.0 Matrix construction
(  427)          ! ------------------------
(  428) 
(  429)          ! Diagonal, inferior, superior
(  430)          ! (including the bottom boundary condition via avt masked)
(  431)          DO jk = 1, jpkm1
(  432)             DO ji = 2, jpim1
(  433)                zwi(ji,jk) = - z2dt(jk) * avt(ji,jj,jk  )   &
(  434)                                        / ( e3t(jk) * e3w(jk  ) )
(  435)                zws(ji,jk) = - z2dt(jk) * avt(ji,jj,jk+1)   &
(  436)                                        / ( e3t(jk) * e3w(jk+1) )
(  437)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  438)             END DO
(  439)          END DO
(  440) 
(  441)          ! Surface boudary conditions
(  442)          DO ji = 2, jpim1
(  443)             zwi(ji,1) = 0.e0
(  444)             zwd(ji,1) = 1. - zws(ji,1)
(  445)          END DO
(  446) 
(  447) 
(  448)          ! II.1. Vertical diffusion on t
(  449)          ! ---------------------------
(  450) 
(  451)          ! Second member construction
(  452) # 463
(  463)          DO jk = 1, jpkm1
(  464)             DO ji = 2, jpim1
(  465)                zwy(ji,jk) = tb(ji,jj,jk) + z2dt(jk) * ta(ji,jj,jk)
(  466)             END DO
(  467)          END DO
(  468) 
(  469) # 470
(  470)          ! Matrix inversion from the first level
(  471)          ikst = 1
(  472) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 19

(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !
(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------
(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO
(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 473 "trazdf_iso.F90"
(  473) 
(  474) # 474
(  474)          ! Save the masked temperature after in ta
(  475)          ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
(  476)          !                  it will not be done in tranxt)
(  477)          DO jk = 1, jpkm1
(  478)             DO ji = 2, jpim1
(  479)                ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
(  480)             END DO
(  481)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 20

(  482) 
(  483) 
(  484)          ! II.2 Vertical diffusion on s
(  485)          ! ---------------------------
(  486) 
(  487) # 508
(  508)          ! Second member construction
(  509) 
(  510) # 521
(  521)          DO jk = 1, jpkm1
(  522)             DO ji = 2, jpim1
(  523)                zwy(ji,jk) = sb(ji,jj,jk) + z2dt(jk) * sa(ji,jj,jk)
(  524)             END DO
(  525)          END DO
(  526) 
(  527) # 528
(  528)          ! Matrix inversion from the first level
(  529)          ikst = 1
(  530) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !
(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------
(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 21

(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 531 "trazdf_iso.F90"
(  531) 
(  532) # 532
(  532)          ! Save the masked salinity after in sa
(  533)          ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
(  534)          !                  it will not be done in tranxt)
(  535)          DO jk = 1, jpkm1
(  536)             DO ji = 2, jpim1
(  537)                sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
(  538)             END DO
(  539)          END DO
(  540) 
(  541) 
(  542)          ! III. recover the avt (avs) resulting from vertical physics only
(  543)          ! ===============================================================
(  544) 
(  545)          DO jk = 2, jpkm1
(  546)             DO ji = 2, jpim1
(  547)                avt(ji,jj,jk) = zavt(ji,jj,jk)
(  548) # 551
(  551)             END DO
(  552)          END DO
(  553) 
(  554)          !                                             ! ===============
(  555)       END DO                                           !   End of slab
(  556)       !                                                ! ===============
(  557) 
(  558)       ! save the trends for diagnostic
(  559)       ! compute the vertical diffusive trends in substracting the previous 
(  560)       ! trends ztdta()/ztdsa() to the new one computed via dT/dt or dS/dt 
(  561)       ! i.e. with the new temperature/salinity ta/sa computed above
(  562)       IF( l_trdtra )   THEN
(  563)          IF( ln_traldf_iso)   THEN
(  564)             DO jk = 1, jpkm1
(  565)                ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk) + ztavg(:,:,jk) 
(  566)                ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk) + zsavg(:,:,jk) 
(  567)             END DO
(  568)          ELSE
(  569)             DO jk = 1, jpkm1
(  570)                ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / z2dt(jk) ) - ztdta(:,:,jk)                             






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 22

(  571)                ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / z2dt(jk) ) - ztdsa(:,:,jk)                             
(  572)             END DO
(  573)          ENDIF
(  574) 
(  575)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  576)       ENDIF
(  577) 
(  578)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  579)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
(  580)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)
(  581)       ENDIF
(  582) 
(  583)    END SUBROUTINE tra_zdf_iso
(  584) 
(  585) 
(  586) # 595
(  595)    !!==============================================================================
(  596) END MODULE trazdf_iso
