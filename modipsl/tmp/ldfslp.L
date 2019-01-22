


PGF90 (Version     10.5)          09/23/2010  09:49:04      page 1

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

Filename: ldfslp.F90

(    1) MODULE ldfslp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldfslp  ***
(    4)    !! Ocean physics: slopes of neutral surfaces
(    5)    !!======================================================================
(    6) #if   defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                      Rotation of lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ldf_slp      : compute the slopes of neutral surface
(   11)    !!   ldf_slp_mxl  : compute the slopes of iso-neutral surface
(   12)    !!   ldf_slp_init : initialization of the slopes computation
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldftra_oce
(   18)    USE ldfdyn_oce
(   19)    USE phycst          ! physical constants
(   20)    USE zdfmxl          ! mixed layer depth
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC ldf_slp        ! routine called by step.F90
(   30) 
(   31)    !! * Share module variables
(   32)    LOGICAL , PUBLIC, PARAMETER ::   lk_ldfslp = .TRUE.     !: slopes flag
(   33)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   34)       uslp, wslpi,         &  !: i_slope at U- and W-points
(   35)       vslp, wslpj             !: j-slope at V- and W-points
(   36)    
(   37)    !! * Module variables
(   38)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   39)       omlmask                 ! mask of the surface mixed layer at T-pt
(   40)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   41)       uslpml, wslpiml,     &  ! i_slope at U- and W-points just below
(   42)    !                          ! the surface mixed layer
(   43)       vslpml, wslpjml         ! j_slope at V- and W-points just below
(   44)    !                          ! the surface mixed layer
(   45) 
(   46)    !! * Substitutions






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 2

(   47) #  include "domzgr_substitute.h90"
(   48) #  include "vectopt_loop_substitute.h90"
(   49)    !!----------------------------------------------------------------------
(   50)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfslp.F90,v 1.6 2005/09/02 15:45:30 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55) CONTAINS
(   56) 
(   57)    SUBROUTINE ldf_slp( kt, prd, pn2 )
(   58)       !!----------------------------------------------------------------------
(   59)       !!                 ***  ROUTINE ldf_slp  ***
(   60)       !! 
(   61)       !! ** Purpose :   Compute the slopes of neutral surface (slope of iso-
(   62)       !!      pycnal surfaces referenced locally) ('key_traldfiso').
(   63)       !!
(   64)       !! ** Method  :   The slope in the i-direction is computed at U- and 
(   65)       !!      W-points (uslp, wslpi) and the slope in the j-direction is 
(   66)       !!      computed at V- and W-points (vslp, wslpj).
(   67)       !!      They are bounded by 1/100 over the whole ocean, and within the
(   68)       !!      surface layer they are bounded by the distance to the surface
(   69)       !!      ( slope<= depth/l  where l is the length scale of horizontal
(   70)       !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
(   71)       !!      of 10cm/s)
(   72)       !!        A horizontal shapiro filter is applied to the slopes
(   73)       !!        'key_s_coord' defined: add to the previously computed slopes
(   74)       !!      the slope of the model level surface.
(   75)       !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
(   76)       !!      [slopes already set to zero at level 1, and to zero or the ocean
(   77)       !!      bottom slope ('key_s_coord' defined) at level jpk in inildf]
(   78)       !!
(   79)       !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes 
(   80)       !!               of now neutral surfaces at u-, w- and v- w-points, resp.
(   81)       !!
(   82)       !! History :
(   83)       !!   7.0  !  94-12  (G. Madec, M. Imbard)  Original code
(   84)       !!   8.0  !  97-06  (G. Madec)  optimization, lbc
(   85)       !!   8.1  !  99-10  (A. Jouzeau)  NEW profile
(   86)       !!   8.5  !  99-10  (G. Madec)  Free form, F90
(   87)       !!----------------------------------------------------------------------
(   88)       !! * Modules used
(   89)       USE oce            , zgru  => ua,  &  ! use ua as workspace
(   90)                            zgrv  => va,  &  ! use va as workspace
(   91)                            zwy   => ta,  &  ! use ta as workspace
(   92)                            zwz   => sa      ! use sa as workspace
(   93)       !! * Arguments
(   94)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   95)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(   96)          prd,                            &  ! in situ density
(   97)          pn2                                ! Brunt-Vaisala frequency (locally ref.)
(   98) 
(   99)       !! * Local declarations
(  100)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(  101)       INTEGER  ::   ii0, ii1, ij0, ij1      ! temporary integer
(  102) #if defined key_partial_steps
(  103)       INTEGER  ::   iku, ikv  ! temporary integers
(  104) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 3

(  105)       REAL(wp) ::   &
(  106)          zeps, zmg, zm05g, zcoef1, zcoef2,   &  ! temporary scalars
(  107)          zau, zbu, zav, zbv,                 &
(  108)          zai, zbi, zaj, zbj,                &
(  109)          zcofu, zcofv, zcofw,                &
(  110)          z1u, z1v, z1wu, z1wv,               &
(  111)          zalpha
(  112)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zww
(  113)       !!----------------------------------------------------------------------
(  114) 
(  115)       
(  116)       ! 0. Initialization (first time-step only)
(  117)       !    --------------
(  118)       
(  119)       IF( kt == nit000 ) CALL ldf_slp_init
(  120)       
(  121) 
(  122)       ! 0. Local constant initialization
(  123)       ! --------------------------------
(  124)       
(  125)       zeps  =  1.e-20
(  126)       zmg   = -1.0 / grav
(  127)       zm05g = -0.5 / grav
(  128) 
(  129)       zww(:,:,:) = 0.e0
(  130)       zwz(:,:,:) = 0.e0
(  131) 
(  132)       ! horizontal density gradient computation
(  133)       DO jk = 1, jpk
(  134)          DO jj = 1, jpjm1
(  135)             DO ji = 1, fs_jpim1   ! vector opt.
(  136)                zgru(ji,jj,jk) = umask(ji,jj,jk) * ( prd(ji+1,jj  ,jk) - prd(ji,jj,jk) ) 
(  137)                zgrv(ji,jj,jk) = vmask(ji,jj,jk) * ( prd(ji  ,jj+1,jk) - prd(ji,jj,jk) ) 
(  138)             END DO
(  139)          END DO
(  140)       END DO
(  141) 
(  142) #if defined key_partial_steps
(  143)       ! partial steps correction at the bottom ocean level (zps_hde routine)
(  144) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  145)       jj = 1
(  146)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  147) # else
(  148)       DO jj = 1, jpjm1
(  149)          DO ji = 1, jpim1
(  150) # endif
(  151)             ! last ocean level
(  152)             iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj) ) - 1
(  153)             ikv = MIN ( mbathy(ji,jj), mbathy(ji,jj+1) ) - 1
(  154)             zgru(ji,jj,iku) = gru(ji,jj) 
(  155)             zgrv(ji,jj,ikv) = grv(ji,jj)               
(  156) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  157)          END DO
(  158) # endif
(  159)       END DO
(  160) #endif
(  161) 
(  162)       ! Slopes of isopycnal surfaces just below the mixed layer






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 4

(  163)       ! -------------------------------------------------------
(  164)       
(  165)       CALL ldf_slp_mxl( prd, pn2 )
(  166)       
(  167)       !-------------------synchro---------------------------------------------
(  168) 
(  169)       !                                                ! ===============
(  170)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  171)          !                                             ! ===============
(  172) 
(  173)          ! I.  slopes at u and v point
(  174)          ! ===========================
(  175)          
(  176)          
(  177)          ! I.1. Slopes of isopycnal surfaces
(  178)          ! ---------------------------------
(  179)          ! uslp = d/di( prd ) / d/dz( prd )
(  180)          ! vslp = d/dj( prd ) / d/dz( prd )
(  181)          
(  182)          ! Local vertical density gradient evaluated from N^2
(  183)          ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  184) 
(  185)          DO jj = 1, jpj
(  186)             DO ji = 1, jpi
(  187)                zwy(ji,jj,jk) = zmg * ( prd(ji,jj,jk) + 1. )                &
(  188)                   &          * ( pn2(ji,jj,jk) + pn2(ji,jj,jk+1) )       &
(  189)                   &          / MAX( tmask(ji,jj,jk) + tmask (ji,jj,jk+1), 1. )
(  190)             END DO
(  191)          END DO
(  192)          
(  193)          ! Slope at u and v points
(  194)          DO jj = 2, jpjm1
(  195)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  196)                ! horizontal and vertical density gradient at u- and v-points
(  197)                zau = 1. / e1u(ji,jj) * zgru(ji,jj,jk)
(  198)                zav = 1. / e2v(ji,jj) * zgrv(ji,jj,jk)
(  199)                zbu = 0.5 * ( zwy(ji,jj,jk) + zwy(ji+1,jj  ,jk) )
(  200)                zbv = 0.5 * ( zwy(ji,jj,jk) + zwy(ji  ,jj+1,jk) )
(  201)                ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  202)                !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  203)                zbu = MIN( zbu, -100.*ABS( zau ), -7.e+3/fse3u(ji,jj,jk)*ABS( zau ) )
(  204)                zbv = MIN( zbv, -100.*ABS( zav ), -7.e+3/fse3v(ji,jj,jk)*ABS( zav ) )
(  205)                ! uslp and vslp output in zwz and zww, resp.
(  206)                zalpha = MAX( omlmask(ji,jj,jk), omlmask(ji+1,jj,jk) )
(  207) #if defined key_s_coord 
(  208)                zwz (ji,jj,jk) = ( zau / ( zbu - zeps ) * ( 1. - zalpha)   &
(  209)                   &        + zalpha * uslpml(ji,jj)   &
(  210)                   &        * ( fsdepu(ji,jj,jk) - .5*fse3u(ji,jj,1) )   &
(  211)                   &        / MAX( hmlpt(ji,jj), hmlpt(ji+1,jj), 5. ) )   &
(  212)                   &        * umask(ji,jj,jk)
(  213)                zalpha = MAX( omlmask(ji,jj,jk), omlmask(ji,jj+1,jk) )
(  214)                zww (ji,jj,jk) = ( zav / ( zbv - zeps ) * ( 1. - zalpha)   &
(  215)                   &        + zalpha * vslpml(ji,jj)   &
(  216)                   &        * ( fsdepv(ji,jj,jk) - .5*fse3v(ji,jj,1) )   &
(  217)                   &         / MAX( hmlpt(ji,jj), hmlpt(ji,jj+1), 5. ) )   &
(  218)                   &        * vmask(ji,jj,jk)
(  219) #else
(  220)                ! z-coord and partial steps slope computed in the same way






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 5

(  221)                zwz (ji,jj,jk) = ( zau / ( zbu - zeps ) * ( 1. - zalpha)    &
(  222)                   &        + zalpha * uslpml(ji,jj)    &
(  223)                   &        * ( fsdept(ji,jj,jk) - .5*fse3u(ji,jj,1))    &
(  224)                   &        / MAX (hmlpt(ji,jj),hmlpt(ji+1,jj),5.) )    &
(  225)                   &        * umask (ji,jj,jk)
(  226)                zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj+1,jk))
(  227)                zww (ji,jj,jk) = ( zav / ( zbv - zeps ) * ( 1. - zalpha)    &
(  228)                   &        + zalpha * vslpml(ji,jj)    &
(  229)                   &        * ( fsdept(ji,jj,jk) - .5*fse3v(ji,jj,1))    &
(  230)                   &        / MAX(hmlpt(ji,jj),hmlpt(ji,jj+1),5.) )    &
(  231)                   &        * vmask (ji,jj,jk)
(  232) #endif
(  233)             END DO
(  234)          END DO
(  235)          !                                             ! ===============
(  236)       END DO                                           !   end of slab
(  237)       !                                                ! ===============
(  238) 
(  239) 
(  240)          ! lateral boundary conditions on zww and zwz
(  241)          CALL lbc_lnk( zwz, 'U', -1. )
(  242)          CALL lbc_lnk( zww, 'V', -1. )
(  243) 
(  244)       !                                                ! ===============
(  245)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  246)          !                                             ! ===============
(  247)          
(  248)          ! Shapiro filter applied in the horizontal direction
(  249)          zcofu = 1. / 16.
(  250)          zcofv = 1. / 16.
(  251)          DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1 only
(  252)             DO ji = 2, jpim1  
(  253)                !uslop
(  254)                uslp(ji,jj,jk) = zcofu * (       zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
(  255)                   &                       +     zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
(  256)                   &                       + 2.*(zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
(  257)                   &                       +     zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
(  258)                   &                       + 4.* zwz(ji  ,jj  ,jk)                       )
(  259)                ! vslop
(  260)                vslp(ji,jj,jk) = zcofv * (       zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
(  261)                   &                       +     zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
(  262)                   &                       + 2.*(zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
(  263)                   &                       +     zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
(  264)                   &                       + 4.* zww(ji,jj    ,jk)                       )
(  265)             END DO
(  266)          END DO
(  267) 
(  268)          DO jj = 3, jpj-2
(  269)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  270)                ! uslop
(  271)                uslp(ji,jj,jk) = zcofu * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
(  272)                   &                       +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
(  273)                   &                       + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
(  274)                   &                       +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
(  275)                   &                       + 4.*  zwz(ji  ,jj  ,jk)                       )
(  276)                ! vslop
(  277)                vslp(ji,jj,jk) = zcofv * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
(  278)                   &                       +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 6

(  279)                   &                       + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
(  280)                   &                       +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
(  281)                   &                       + 4.*  zww(ji,jj    ,jk)                       )
(  282)             END DO
(  283)          END DO
(  284) 
(  285)          ! decrease along coastal boundaries
(  286)          DO jj = 2, jpjm1
(  287)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  288)                z1u  = ( umask(ji,jj+1,jk) + umask(ji,jj-1,jk) )*.5
(  289)                z1v  = ( vmask(ji+1,jj,jk) + vmask(ji-1,jj,jk) )*.5
(  290)                z1wu = ( umask(ji,jj,jk)   + umask(ji,jj,jk+1) )*.5
(  291)                z1wv = ( vmask(ji,jj,jk)   + vmask(ji,jj,jk+1) )*.5
(  292)                uslp(ji,jj,jk) = uslp(ji,jj,jk) * z1u * z1wu
(  293)                vslp(ji,jj,jk) = vslp(ji,jj,jk) * z1v * z1wv
(  294)             END DO
(  295)          END DO
(  296) 
(  297) 
(  298)          IF( lk_sco ) THEN
(  299)             ! Add the slope of level surfaces
(  300)             ! -----------------------------------
(  301)             ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
(  302)             ! in inildf, ldfslp never called
(  303)             ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
(  304)             ! is added to the slope of isopycnal surfaces.
(  305)             ! c a u t i o n : minus sign as fsdep has positive value 
(  306)          
(  307)             DO jj = 2, jpjm1
(  308)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  309)                   uslp(ji,jj,jk) = uslp(ji,jj,jk) - 1. / e1u(ji,jj)   &
(  310)                      &           * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )
(  311)                   vslp(ji,jj,jk) = vslp(ji,jj,jk) - 1. / e2v(ji,jj)   &
(  312)                      &           * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
(  313)                END DO
(  314)             END DO
(  315)          ENDIF
(  316) 
(  317) 
(  318)          ! II. Computation of slopes at w point
(  319)          ! ====================================
(  320)          
(  321)          
(  322)          ! II.1 Slopes of isopycnal surfaces
(  323)          ! ---------------------------------
(  324)          ! wslpi = mij( d/di( prd ) / d/dz( prd )
(  325)          ! wslpj = mij( d/dj( prd ) / d/dz( prd )
(  326)          
(  327)          
(  328)          ! Local vertical density gradient evaluated from N^2
(  329)          !     zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
(  330)          DO jj = 1, jpj
(  331)             DO ji = 1, jpi
(  332)                zwy (ji,jj,jk) = zm05g * pn2 (ji,jj,jk) *   &
(  333)                   &                     ( prd (ji,jj,jk) + prd (ji,jj,jk-1) + 2. )
(  334)             END DO
(  335)          END DO
(  336)          






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 7

(  337)          ! Slope at w point
(  338)          DO jj = 2, jpjm1
(  339)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  340)                ! horizontal density i-gradient at w-points
(  341)                zcoef1 = MAX( zeps, umask(ji-1,jj,jk  )+umask(ji,jj,jk  )    &
(  342)                   &               +umask(ji-1,jj,jk-1)+umask(ji,jj,jk-1) )
(  343)                zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
(  344)                zai = zcoef1 * (  zgru(ji  ,jj,jk  ) + zgru(ji  ,jj,jk-1)   &
(  345)                   &            + zgru(ji-1,jj,jk-1) + zgru(ji-1,jj,jk  ) ) * tmask (ji,jj,jk)
(  346)                ! horizontal density j-gradient at w-points
(  347)                zcoef2 = MAX( zeps, vmask(ji,jj-1,jk  )+vmask(ji,jj,jk-1)   &
(  348)                   &               +vmask(ji,jj-1,jk-1)+vmask(ji,jj,jk  ) )
(  349)                zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
(  350)                zaj = zcoef2 * (  zgrv(ji,jj  ,jk  ) + zgrv(ji,jj  ,jk-1)   &
(  351)                   &            + zgrv(ji,jj-1,jk-1) + zgrv(ji,jj-1,jk  ) ) * tmask (ji,jj,jk)
(  352)                ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
(  353)                !                   static instability:
(  354)                !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  355)                zbi = MIN( zwy (ji,jj,jk),- 100.*ABS(zai), -7.e+3/fse3w(ji,jj,jk)*ABS(zai) )
(  356)                zbj = MIN( zwy (ji,jj,jk), -100.*ABS(zaj), -7.e+3/fse3w(ji,jj,jk)*ABS(zaj) )
(  357)                ! wslpi and wslpj output in zwz and zww, resp.
(  358)                zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj,jk-1))
(  359)                zwz(ji,jj,jk) = ( zai / ( zbi - zeps) * ( 1. - zalpha )   &
(  360)                   &            + zalpha * wslpiml(ji,jj)   &
(  361)                   &            * fsdepw(ji,jj,jk) / MAX( hmlp(ji,jj),10. ) )   &
(  362)                   &            * tmask (ji,jj,jk)
(  363)                zww(ji,jj,jk) = ( zaj / ( zbj - zeps) * ( 1. - zalpha )   &
(  364)                   &            + zalpha * wslpjml(ji,jj)   &
(  365)                   &            * fsdepw(ji,jj,jk) / MAX( hmlp(ji,jj),10. ) )   &
(  366)                   &            * tmask (ji,jj,jk)
(  367)             END DO
(  368)          END DO
(  369)          !                                             ! ===============
(  370)       END DO                                           !   end of slab
(  371)       !                                                ! ===============
(  372) 
(  373) 
(  374)       ! lateral boundary conditions on zwz and zww
(  375)       CALL lbc_lnk( zwz, 'T', -1. )
(  376)       CALL lbc_lnk( zww, 'T', -1. )
(  377) 
(  378)       !                                                ! ===============
(  379)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  380)          !                                             ! ===============
(  381) 
(  382)          ! Shapiro filter applied in the horizontal direction
(  383) 
(  384)          DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1
(  385)             DO ji = 2, jpim1
(  386)                zcofw = tmask(ji,jj,jk)/16.
(  387)                wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
(  388)                   &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
(  389)                   &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
(  390)                   &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
(  391)                   &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw
(  392) 
(  393)                wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
(  394)                   &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 8

(  395)                   &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
(  396)                   &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
(  397)                   &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
(  398)             END DO
(  399)          END DO
(  400)          
(  401)          DO jj = 3, jpj-2
(  402)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  403)                zcofw = tmask(ji,jj,jk)/16.
(  404)                wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
(  405)                   &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
(  406)                   &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
(  407)                   &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
(  408)                   &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw
(  409) 
(  410)                wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
(  411)                   &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
(  412)                   &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
(  413)                   &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
(  414)                   &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
(  415)             END DO
(  416)          END DO
(  417)          
(  418)          ! decrease the slope along the boundaries
(  419)          DO jj = 2, jpjm1
(  420)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  421)                z1u = ( umask(ji,jj,jk) + umask(ji-1,jj,jk) ) *.5
(  422)                z1v = ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) *.5
(  423)                wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * z1u * z1v
(  424)                wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * z1u * z1v
(  425)             END DO
(  426)          END DO
(  427)          
(  428)          IF( lk_sco ) THEN
(  429)          
(  430)             ! Slope of level surfaces
(  431)             ! -----------------------
(  432)             ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
(  433)             ! in inildf, ldfslp never called
(  434)             ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
(  435)             ! is added to the slope of isopycnal surfaces.
(  436)          
(  437)             DO jj = 2, jpjm1
(  438)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  439)                   wslpi(ji,jj,jk) = wslpi(ji,jj,jk) - 1. / e1t(ji,jj)   &
(  440)                      &                                   * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
(  441)                   wslpj(ji,jj,jk) = wslpj(ji,jj,jk) - 1. / e2t(ji,jj)   &
(  442)                      &                                   * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
(  443)                END DO
(  444)             END DO
(  445)          ENDIF
(  446)          
(  447)          ! III. Specific grid points
(  448)          ! -------------------------
(  449) 
(  450)          IF( cp_cfg == "orca" .AND. jp_cfg == 4 ) THEN
(  451)             !                                        ! =======================
(  452)             ! Horizontal diffusion in                !  ORCA_R4 configuration 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 9

(  453)             ! specific area                          ! =======================
(  454)             !
(  455)             !                                             ! Gibraltar Strait
(  456)             ij0 =  50   ;   ij1 =  53
(  457)             ii0 =  69   ;   ii1 =  71   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  458)             ij0 =  51   ;   ij1 =  53
(  459)             ii0 =  68   ;   ii1 =  71   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  460)             ii0 =  69   ;   ii1 =  71   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  461)             ii0 =  69   ;   ii1 =  71   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  462) 
(  463)             !                                             ! Mediterrannean Sea
(  464)             ij0 =  49   ;   ij1 =  56
(  465)             ii0 =  71   ;   ii1 =  90   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  466)             ij0 =  50   ;   ij1 =  56
(  467)             ii0 =  70   ;   ii1 =  90   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  468)             ii0 =  71   ;   ii1 =  90   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  469)             ii0 =  71   ;   ii1 =  90   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  470)          ENDIF
(  471)          !                                             ! ===============
(  472)       END DO                                           !   end of slab
(  473)       !                                                ! ===============       
(  474) 
(  475)       
(  476)       ! III Lateral boundary conditions on all slopes (uslp , vslp, 
(  477)       ! -------------------------------                wslpi, wslpj )
(  478)       CALL lbc_lnk( uslp , 'U', -1. )
(  479)       CALL lbc_lnk( vslp , 'V', -1. )
(  480)       CALL lbc_lnk( wslpi, 'W', -1. )
(  481)       CALL lbc_lnk( wslpj, 'W', -1. )
(  482) 
(  483)       IF(ln_ctl) THEN
(  484)          CALL prt_ctl(tab3d_1=uslp , clinfo1=' slp  - u : ', tab3d_2=vslp,  clinfo2=' v : ', kdim=jpk)
(  485)          CALL prt_ctl(tab3d_1=wslpi, clinfo1=' slp  - wi: ', tab3d_2=wslpj, clinfo2=' wj: ', kdim=jpk)
(  486)       ENDIF
(  487) 
(  488)    END SUBROUTINE ldf_slp
(  489) 
(  490) 
(  491)    SUBROUTINE ldf_slp_mxl( prd, pn2 )
(  492)       !!----------------------------------------------------------------------
(  493)       !!                  ***  ROUTINE ldf_slp_mxl  ***
(  494)       !! ** Purpose :
(  495)       !!     Compute the slopes of iso-neutral surface (slope of isopycnal
(  496)       !!   surfaces referenced locally) just above the mixed layer.
(  497)       !!
(  498)       !! ** Method :
(  499)       !!      The slope in the i-direction is computed at u- and w-points
(  500)       !!   (uslp, wslpi) and the slope in the j-direction is computed at
(  501)       !!   v- and w-points (vslp, wslpj).
(  502)       !!   They are bounded by 1/100 over the whole ocean, and within the
(  503)       !!   surface layer they are bounded by the distance to the surface
(  504)       !!   ( slope<= depth/l  where l is the length scale of horizontal
(  505)       !!   diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
(  506)       !!   of 10cm/s)
(  507)       !!
(  508)       !! ** Action :
(  509)       !!      Compute uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
(  510)       !!   of now neutral surfaces at u-, w- and v- w-points, resp.






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 10

(  511)       !!
(  512)       !! History :
(  513)       !!   8.1  !  99-10  (A. Jouzeau)  Original code
(  514)       !!   8.5  !  99-10  (G. Madec)  Free form, F90
(  515)       !!----------------------------------------------------------------------
(  516)       !! * Modules used
(  517)       USE oce           , zgru  => ua,  &  ! ua, va used as workspace and set to hor. 
(  518)                           zgrv  => va      ! density gradient in ldf_slp
(  519) 
(  520)       !! * Arguments
(  521)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  522)          prd,                           &  ! in situ density
(  523)          pn2                              ! Brunt-Vaisala frequency (locally ref.)
(  524) 
(  525)       !! * Local declarations
(  526)       INTEGER  ::   ji, jj, jk             ! dummy loop indices
(  527)       INTEGER  ::   ik, ikm1               ! temporary integers
(  528)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  529)          zwy                               ! temporary workspace
(  530)       REAL(wp) ::   &
(  531)          zeps, zmg, zm05g,              &  ! temporary scalars
(  532)          zcoef1, zcoef2,                &  !    "         "
(  533)          zau, zbu, zav, zbv,            &  !    "         "
(  534)          zai, zbi, zaj, zbj                !    "         "
(  535)       !!----------------------------------------------------------------------
(  536) 
(  537) 
(  538)       ! 0. Local constant initialization
(  539)       ! --------------------------------
(  540) 
(  541)       zeps  =  1.e-20
(  542)       zmg   = -1.0 / grav
(  543)       zm05g = -0.5 / grav
(  544) 
(  545) 
(  546)       uslpml (1,:) = 0.e0      ;      uslpml (jpi,:) = 0.e0
(  547)       vslpml (1,:) = 0.e0      ;      vslpml (jpi,:) = 0.e0
(  548)       wslpiml(1,:) = 0.e0      ;      wslpiml(jpi,:) = 0.e0
(  549)       wslpjml(1,:) = 0.e0      ;      wslpjml(jpi,:) = 0.e0
(  550) 
(  551)       ! surface mixed layer mask
(  552) 
(  553)       ! mask for mixed layer
(  554)       DO jk = 1, jpk
(  555) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  556)          jj = 1
(  557)          DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  558) # else
(  559)          DO jj = 1, jpj
(  560)             DO ji = 1, jpi
(  561) # endif
(  562)                ! mixed layer interior (mask = 1) and exterior (mask = 0)
(  563)                ik = nmln(ji,jj) - 1
(  564)                IF( jk <= ik ) THEN
(  565)                   omlmask(ji,jj,jk) = 1.e0
(  566)                ELSE
(  567)                   omlmask(ji,jj,jk) = 0.e0
(  568)                ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 11

(  569) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  570)             END DO
(  571) # endif
(  572)          END DO
(  573)       END DO
(  574) 
(  575) 
(  576)       ! Slopes of isopycnal surfaces just before bottom of mixed layer
(  577)       ! --------------------------------------------------------------
(  578)       ! uslpml = d/di( prd ) / d/dz( prd )
(  579)       ! vslpml = d/dj( prd ) / d/dz( prd )
(  580) 
(  581)       ! Local vertical density gradient evaluated from N^2
(  582)       ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  583) 
(  584)       !-----------------------------------------------------------------------
(  585)       zwy(:,jpj) = 0.e0
(  586)       zwy(jpi,:) = 0.e0
(  587) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  588)       jj = 1
(  589)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  590) # else
(  591)       DO jj = 1, jpjm1
(  592)          DO ji = 1, jpim1
(  593) # endif
(  594)             ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
(  595)             ! if ik = jpk take jpkm1 values
(  596)             ik = MIN( ik,jpkm1 )
(  597)             zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
(  598)                &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
(  599)                &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
(  600) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  601)          END DO
(  602) # endif
(  603)       END DO
(  604)       ! lateral boundary conditions on zwy
(  605)       CALL lbc_lnk( zwy, 'U', -1. )
(  606) 
(  607)       ! Slope at u points
(  608) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  609)       jj = 1
(  610)       DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  611) # else
(  612)       DO jj = 2, jpjm1
(  613)          DO ji = 2, jpim1
(  614) # endif
(  615)             ! horizontal and vertical density gradient at u-points
(  616)             ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
(  617)             ik = MIN( ik,jpkm1 )
(  618)             zau = 1./ e1u(ji,jj) * zgru(ji,jj,ik)
(  619)             zbu = 0.5*( zwy(ji,jj) + zwy(ji+1,jj) )
(  620)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  621)             !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  622)             zbu = MIN( zbu, -100.*ABS(zau), -7.e+3/fse3u(ji,jj,ik)*ABS(zau) )
(  623)             ! uslpml
(  624)             uslpml (ji,jj) = zau / ( zbu - zeps ) * umask (ji,jj,ik)
(  625) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  626)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 12

(  627) # endif
(  628)       END DO
(  629) 
(  630)       ! lateral boundary conditions on uslpml
(  631)       CALL lbc_lnk( uslpml, 'U', -1. )
(  632) 
(  633)       ! Local vertical density gradient evaluated from N^2
(  634)       !     zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  635)       zwy ( :, jpj) = 0.e0
(  636)       zwy ( jpi, :) = 0.e0
(  637) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  638)       jj = 1
(  639)       DO ji = 1, jpij-jpi   ! vector opt. (forced unrolling)
(  640) # else
(  641)       DO jj = 1, jpjm1
(  642)          DO ji = 1, jpim1
(  643) # endif
(  644)             ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
(  645)             ik = MIN( ik,jpkm1 )
(  646)             zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
(  647)                &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
(  648)                &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
(  649) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  650)          END DO
(  651) # endif
(  652)       END DO
(  653) 
(  654)       ! lateral boundary conditions on zwy
(  655)       CALL lbc_lnk( zwy, 'V', -1. )
(  656) 
(  657)       ! Slope at v points
(  658) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  659)       jj = 1
(  660)       DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  661) # else
(  662)       DO jj = 2, jpjm1
(  663)          DO ji = 2, jpim1
(  664) # endif
(  665)             ! horizontal and vertical density gradient at v-points
(  666)             ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
(  667)             ik = MIN( ik,jpkm1 )
(  668)             zav = 1./ e2v(ji,jj) * zgrv(ji,jj,ik)
(  669)             zbv = 0.5*( zwy(ji,jj) + zwy(ji,jj+1) )
(  670)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  671)             !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  672)             zbv = MIN( zbv, -100.*ABS(zav), -7.e+3/fse3v(ji,jj,ik)*ABS( zav ) )
(  673)             ! vslpml
(  674)             vslpml (ji,jj) = zav / ( zbv - zeps ) * vmask (ji,jj,ik)
(  675) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  676)          END DO
(  677) # endif
(  678)       END DO
(  679) 
(  680)       ! lateral boundary conditions on vslpml
(  681)       CALL lbc_lnk( vslpml, 'V', -1. )
(  682) 
(  683)       ! wslpiml = mij( d/di( prd ) / d/dz( prd )
(  684)       ! wslpjml = mij( d/dj( prd ) / d/dz( prd )






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 13

(  685) 
(  686) 
(  687)       ! Local vertical density gradient evaluated from N^2
(  688)       ! zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
(  689) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  690)       jj = 1
(  691)       DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  692) # else
(  693)       DO jj = 1, jpj
(  694)          DO ji = 1, jpi
(  695) # endif
(  696)             ik = nmln(ji,jj)+1
(  697)             ik = MIN( ik,jpk )
(  698)             ikm1 = MAX ( 1, ik-1)
(  699)             zwy (ji,jj) = zm05g * pn2 (ji,jj,ik) *     &
(  700)                &             ( prd (ji,jj,ik) + prd (ji,jj,ikm1) + 2. )
(  701) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  702)          END DO
(  703) # endif
(  704)       END DO
(  705) 
(  706)       ! Slope at w point
(  707) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  708)       jj = 1
(  709)       DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  710) # else
(  711)       DO jj = 2, jpjm1
(  712)          DO ji = 2, jpim1
(  713) # endif
(  714)             ik = nmln(ji,jj)+1
(  715)             ik = MIN( ik,jpk )
(  716)             ikm1 = MAX ( 1, ik-1)
(  717)             ! horizontal density i-gradient at w-points
(  718)             zcoef1 = MAX( zeps, umask(ji-1,jj,ik  )+umask(ji,jj,ik  )   &
(  719)                &               +umask(ji-1,jj,ikm1)+umask(ji,jj,ikm1) )
(  720)             zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
(  721)             zai = zcoef1 * (  zgru(ji  ,jj,ik  ) + zgru(ji  ,jj,ikm1)   &
(  722)                &            + zgru(ji-1,jj,ikm1) + zgru(ji-1,jj,ik  ) ) * tmask (ji,jj,ik)
(  723)             ! horizontal density j-gradient at w-points
(  724)             zcoef2 = MAX( zeps, vmask(ji,jj-1,ik  )+vmask(ji,jj,ikm1)    &
(  725)                &               +vmask(ji,jj-1,ikm1)+vmask(ji,jj,ik  ) )
(  726)             zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
(  727)             zaj = zcoef2 * (  zgrv(ji,jj  ,ik  ) + zgrv(ji,jj  ,ikm1)   &
(  728)                &            + zgrv(ji,jj-1,ikm1) + zgrv(ji,jj-1,ik  ) ) * tmask (ji,jj,ik)
(  729)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
(  730)             !                   static instability:
(  731)             !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  732)             zbi = MIN ( zwy (ji,jj),- 100.*ABS(zai), -7.e+3/fse3w(ji,jj,ik)*ABS(zai) )
(  733)             zbj = MIN ( zwy (ji,jj), -100.*ABS(zaj), -7.e+3/fse3w(ji,jj,ik)*ABS(zaj) )
(  734)             ! wslpiml and wslpjml
(  735)             wslpiml (ji,jj) = zai / ( zbi - zeps) * tmask (ji,jj,ik)
(  736)             wslpjml (ji,jj) = zaj / ( zbj - zeps) * tmask (ji,jj,ik)
(  737) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  738)          END DO
(  739) # endif
(  740)       END DO
(  741) 
(  742)       ! lateral boundary conditions on wslpiml and wslpjml






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 14

(  743)       CALL lbc_lnk( wslpiml, 'W', -1. )
(  744)       CALL lbc_lnk( wslpjml, 'W', -1. )
(  745) 
(  746)    END SUBROUTINE ldf_slp_mxl
(  747) 
(  748) 
(  749)    SUBROUTINE ldf_slp_init
(  750)       !!----------------------------------------------------------------------
(  751)       !!                  ***  ROUTINE ldf_slp_init  ***
(  752)       !!
(  753)       !! ** Purpose :   Initialization for the isopycnal slopes computation
(  754)       !!
(  755)       !! ** Method  :   read the nammbf namelist and check the parameter 
(  756)       !!      values called by tra_dmp at the first timestep (nit000)
(  757)       !!
(  758)       !! History :
(  759)       !!  8.5  ! 02-06 (G. Madec) original code
(  760)       !!----------------------------------------------------------------------
(  761)       !! * local declarations
(  762)       INTEGER ::   ji, jj, jk   ! dummy loop indices
(  763)       !!----------------------------------------------------------------------
(  764)       
(  765)       
(  766)       ! Parameter control and print
(  767)       ! ---------------------------
(  768)       IF(lwp) THEN
(  769)          WRITE(numout,*)
(  770)          WRITE(numout,*) 'ldf_slp : direction of lateral mixing'
(  771)          WRITE(numout,*) '~~~~~~~'
(  772)       ENDIF
(  773) 
(  774)       ! Direction of lateral diffusion (tracers and/or momentum)
(  775)       ! ------------------------------
(  776)       ! set the slope to zero (even in s-coordinates)
(  777) 
(  778)       uslp (:,:,:) = 0.e0
(  779)       vslp (:,:,:) = 0.e0
(  780)       wslpi(:,:,:) = 0.e0
(  781)       wslpj(:,:,:) = 0.e0
(  782) 
(  783)       uslpml (:,:) = 0.e0
(  784)       vslpml (:,:) = 0.e0
(  785)       wslpiml(:,:) = 0.e0
(  786)       wslpjml(:,:) = 0.e0
(  787) 
(  788)       IF( ln_traldf_hor .OR. ln_dynldf_hor ) THEN
(  789) 
(  790)          ! geopotential diffusion in s-coordinates on tracers and/or momentum
(  791)          ! The slopes of s-surfaces are computed once (no call to ldfslp in step)
(  792)          ! The slopes for momentum diffusion are i- or j- averaged of those on tracers
(  793) 
(  794)          ! set the slope of diffusion to the slope of s-surfaces
(  795)          !      ( c a u t i o n : minus sign as fsdep has positive value )
(  796)          DO jk = 1, jpk
(  797)             DO jj = 2, jpjm1
(  798)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  799)                   uslp (ji,jj,jk) = -1. / e1u(ji,jj) * umask(ji,jj,jk)   &
(  800)                      &                               * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 15

(  801)                   vslp (ji,jj,jk) = -1. / e2v(ji,jj) * vmask(ji,jj,jk)   &
(  802)                      &                               * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
(  803)                   wslpi(ji,jj,jk) = -1. / e1t(ji,jj) * tmask(ji,jj,jk)   &
(  804)                      &                               * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
(  805)                   wslpj(ji,jj,jk) = -1. / e2t(ji,jj) * tmask(ji,jj,jk)   &
(  806)                      &                               * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
(  807)                END DO
(  808)             END DO
(  809)          END DO
(  810) 
(  811)          ! Lateral boundary conditions on the slopes
(  812)          CALL lbc_lnk( uslp , 'U', -1. )
(  813)          CALL lbc_lnk( vslp , 'V', -1. )
(  814)          CALL lbc_lnk( wslpi, 'W', -1. )
(  815)          CALL lbc_lnk( wslpj, 'W', -1. )
(  816)       ENDIF
(  817) 
(  818)    END SUBROUTINE ldf_slp_init
(  819) 
(  820) #else
(  821)    !!------------------------------------------------------------------------
(  822)    !!   Dummy module :                 NO Rotation of lateral mixing tensor
(  823)    !!------------------------------------------------------------------------
(  824)    USE in_out_manager
(  825)    LOGICAL, PUBLIC, PARAMETER ::   lk_ldfslp = .FALSE.    !: slopes flag
(  826) CONTAINS
(  827)    SUBROUTINE ldf_slp( kt, prd, pn2 )        ! Dummy routine
(  828)       INTEGER, INTENT(in) :: kt 
(  829)       REAL,DIMENSION(:,:,:), INTENT(in) :: prd, pn2
(  830)       if(lwp)WRITE(numout,*) 'ldf_slp: You should not have seen this print! error?', kt, prd(1,1,1), pn2(1,1,1)
(  831)    END SUBROUTINE ldf_slp
(  832) #endif
(  833) 
(  834)    !!======================================================================
(  835) END MODULE ldfslp





























PGF90 (Version     10.5)          09/23/2010  09:49:04      page 16

(    1) # 1 "ldfslp.F90"
(    1) MODULE ldfslp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ldfslp  ***
(    4)    !! Ocean physics: slopes of neutral surfaces
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                      Rotation of lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ldf_slp      : compute the slopes of neutral surface
(   11)    !!   ldf_slp_mxl  : compute the slopes of iso-neutral surface
(   12)    !!   ldf_slp_init : initialization of the slopes computation
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce             ! ocean dynamics and tracers
(   16)    USE dom_oce         ! ocean space and time domain
(   17)    USE ldftra_oce
(   18)    USE ldfdyn_oce
(   19)    USE phycst          ! physical constants
(   20)    USE zdfmxl          ! mixed layer depth
(   21)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC ldf_slp        ! routine called by step.F90
(   30) 
(   31)    !! * Share module variables
(   32)    LOGICAL , PUBLIC, PARAMETER ::   lk_ldfslp = .TRUE.     !: slopes flag
(   33)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
(   34)       uslp, wslpi,         &  !: i_slope at U- and W-points
(   35)       vslp, wslpj             !: j-slope at V- and W-points
(   36)    
(   37)    !! * Module variables
(   38)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   39)       omlmask                 ! mask of the surface mixed layer at T-pt
(   40)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   41)       uslpml, wslpiml,     &  ! i_slope at U- and W-points just below
(   42)    !                          ! the surface mixed layer
(   43)       vslpml, wslpjml         ! j_slope at V- and W-points just below
(   44)    !                          ! the surface mixed layer
(   45) 
(   46)    !! * Substitutions
(   47) # 1 "./domzgr_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 17

(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 48 "ldfslp.F90"
(   48) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 49 "ldfslp.F90"
(   49) # 49
(   49)    !!----------------------------------------------------------------------
(   50)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfslp.F90,v 1.6 2005/09/02 15:45:30 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55) CONTAINS
(   56) 
(   57)    SUBROUTINE ldf_slp( kt, prd, pn2 )
(   58)       !!----------------------------------------------------------------------
(   59)       !!                 ***  ROUTINE ldf_slp  ***
(   60)       !! 
(   61)       !! ** Purpose :   Compute the slopes of neutral surface (slope of iso-
(   62)       !!      pycnal surfaces referenced locally) ('key_traldfiso').
(   63)       !!
(   64)       !! ** Method  :   The slope in the i-direction is computed at U- and 
(   65)       !!      W-points (uslp, wslpi) and the slope in the j-direction is 
(   66)       !!      computed at V- and W-points (vslp, wslpj).
(   67)       !!      They are bounded by 1/100 over the whole ocean, and within the
(   68)       !!      surface layer they are bounded by the distance to the surface
(   69)       !!      ( slope<= depth/l  where l is the length scale of horizontal
(   70)       !!      diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
(   71)       !!      of 10cm/s)
(   72)       !!        A horizontal shapiro filter is applied to the slopes
(   73)       !!        'key_s_coord' defined: add to the previously computed slopes
(   74)       !!      the slope of the model level surface.
(   75)       !!        macro-tasked on horizontal slab (jk-loop)  (2, jpk-1)
(   76)       !!      [slopes already set to zero at level 1, and to zero or the ocean
(   77)       !!      bottom slope ('key_s_coord' defined) at level jpk in inildf]
(   78)       !!
(   79)       !! ** Action : - uslp, wslpi, and vslp, wslpj, the i- and  j-slopes 
(   80)       !!               of now neutral surfaces at u-, w- and v- w-points, resp.
(   81)       !!
(   82)       !! History :
(   83)       !!   7.0  !  94-12  (G. Madec, M. Imbard)  Original code
(   84)       !!   8.0  !  97-06  (G. Madec)  optimization, lbc
(   85)       !!   8.1  !  99-10  (A. Jouzeau)  NEW profile
(   86)       !!   8.5  !  99-10  (G. Madec)  Free form, F90
(   87)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 18

(   88)       !! * Modules used
(   89)       USE oce            , zgru  => ua,  &  ! use ua as workspace
(   90)                            zgrv  => va,  &  ! use va as workspace
(   91)                            zwy   => ta,  &  ! use ta as workspace
(   92)                            zwz   => sa      ! use sa as workspace
(   93)       !! * Arguments
(   94)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   95)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(   96)          prd,                            &  ! in situ density
(   97)          pn2                                ! Brunt-Vaisala frequency (locally ref.)
(   98) 
(   99)       !! * Local declarations
(  100)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(  101)       INTEGER  ::   ii0, ii1, ij0, ij1      ! temporary integer
(  102) # 105
(  105)       REAL(wp) ::   &
(  106)          zeps, zmg, zm05g, zcoef1, zcoef2,   &  ! temporary scalars
(  107)          zau, zbu, zav, zbv,                 &
(  108)          zai, zbi, zaj, zbj,                &
(  109)          zcofu, zcofv, zcofw,                &
(  110)          z1u, z1v, z1wu, z1wv,               &
(  111)          zalpha
(  112)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zww
(  113)       !!----------------------------------------------------------------------
(  114) 
(  115)       
(  116)       ! 0. Initialization (first time-step only)
(  117)       !    --------------
(  118)       
(  119)       IF( kt == nit000 ) CALL ldf_slp_init
(  120)       
(  121) 
(  122)       ! 0. Local constant initialization
(  123)       ! --------------------------------
(  124)       
(  125)       zeps  =  1.e-20
(  126)       zmg   = -1.0 / grav
(  127)       zm05g = -0.5 / grav
(  128) 
(  129)       zww(:,:,:) = 0.e0
(  130)       zwz(:,:,:) = 0.e0
(  131) 
(  132)       ! horizontal density gradient computation
(  133)       DO jk = 1, jpk
(  134)          DO jj = 1, jpjm1
(  135)             DO ji = 1, jpim1   ! vector opt.
(  136)                zgru(ji,jj,jk) = umask(ji,jj,jk) * ( prd(ji+1,jj  ,jk) - prd(ji,jj,jk) ) 
(  137)                zgrv(ji,jj,jk) = vmask(ji,jj,jk) * ( prd(ji  ,jj+1,jk) - prd(ji,jj,jk) ) 
(  138)             END DO
(  139)          END DO
(  140)       END DO
(  141) 
(  142) 
(  143) # 162
(  162)       ! Slopes of isopycnal surfaces just below the mixed layer
(  163)       ! -------------------------------------------------------
(  164)       
(  165)       CALL ldf_slp_mxl( prd, pn2 )






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 19

(  166)       
(  167)       !-------------------synchro---------------------------------------------
(  168) 
(  169)       !                                                ! ===============
(  170)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  171)          !                                             ! ===============
(  172) 
(  173)          ! I.  slopes at u and v point
(  174)          ! ===========================
(  175)          
(  176)          
(  177)          ! I.1. Slopes of isopycnal surfaces
(  178)          ! ---------------------------------
(  179)          ! uslp = d/di( prd ) / d/dz( prd )
(  180)          ! vslp = d/dj( prd ) / d/dz( prd )
(  181)          
(  182)          ! Local vertical density gradient evaluated from N^2
(  183)          ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  184) 
(  185)          DO jj = 1, jpj
(  186)             DO ji = 1, jpi
(  187)                zwy(ji,jj,jk) = zmg * ( prd(ji,jj,jk) + 1. )                &
(  188)                   &          * ( pn2(ji,jj,jk) + pn2(ji,jj,jk+1) )       &
(  189)                   &          / MAX( tmask(ji,jj,jk) + tmask (ji,jj,jk+1), 1. )
(  190)             END DO
(  191)          END DO
(  192)          
(  193)          ! Slope at u and v points
(  194)          DO jj = 2, jpjm1
(  195)             DO ji = 2, jpim1   ! vector opt.
(  196)                ! horizontal and vertical density gradient at u- and v-points
(  197)                zau = 1. / e1u(ji,jj) * zgru(ji,jj,jk)
(  198)                zav = 1. / e2v(ji,jj) * zgrv(ji,jj,jk)
(  199)                zbu = 0.5 * ( zwy(ji,jj,jk) + zwy(ji+1,jj  ,jk) )
(  200)                zbv = 0.5 * ( zwy(ji,jj,jk) + zwy(ji  ,jj+1,jk) )
(  201)                ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  202)                !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  203)                zbu = MIN( zbu, -100.*ABS( zau ), -7.e+3/e3t(jk)*ABS( zau ) )
(  204)                zbv = MIN( zbv, -100.*ABS( zav ), -7.e+3/e3t(jk)*ABS( zav ) )
(  205)                ! uslp and vslp output in zwz and zww, resp.
(  206)                zalpha = MAX( omlmask(ji,jj,jk), omlmask(ji+1,jj,jk) )
(  207) # 220
(  220)                ! z-coord and partial steps slope computed in the same way
(  221)                zwz (ji,jj,jk) = ( zau / ( zbu - zeps ) * ( 1. - zalpha)    &
(  222)                   &        + zalpha * uslpml(ji,jj)    &
(  223)                   &        * ( gdept(jk) - .5*e3t(1))    &
(  224)                   &        / MAX (hmlpt(ji,jj),hmlpt(ji+1,jj),5.) )    &
(  225)                   &        * umask (ji,jj,jk)
(  226)                zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj+1,jk))
(  227)                zww (ji,jj,jk) = ( zav / ( zbv - zeps ) * ( 1. - zalpha)    &
(  228)                   &        + zalpha * vslpml(ji,jj)    &
(  229)                   &        * ( gdept(jk) - .5*e3t(1))    &
(  230)                   &        / MAX(hmlpt(ji,jj),hmlpt(ji,jj+1),5.) )    &
(  231)                   &        * vmask (ji,jj,jk)
(  232) # 233
(  233)             END DO
(  234)          END DO
(  235)          !                                             ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 20

(  236)       END DO                                           !   end of slab
(  237)       !                                                ! ===============
(  238) 
(  239) 
(  240)          ! lateral boundary conditions on zww and zwz
(  241)          CALL lbc_lnk( zwz, 'U', -1. )
(  242)          CALL lbc_lnk( zww, 'V', -1. )
(  243) 
(  244)       !                                                ! ===============
(  245)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  246)          !                                             ! ===============
(  247)          
(  248)          ! Shapiro filter applied in the horizontal direction
(  249)          zcofu = 1. / 16.
(  250)          zcofv = 1. / 16.
(  251)          DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1 only
(  252)             DO ji = 2, jpim1  
(  253)                !uslop
(  254)                uslp(ji,jj,jk) = zcofu * (       zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
(  255)                   &                       +     zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
(  256)                   &                       + 2.*(zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
(  257)                   &                       +     zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
(  258)                   &                       + 4.* zwz(ji  ,jj  ,jk)                       )
(  259)                ! vslop
(  260)                vslp(ji,jj,jk) = zcofv * (       zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
(  261)                   &                       +     zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
(  262)                   &                       + 2.*(zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
(  263)                   &                       +     zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
(  264)                   &                       + 4.* zww(ji,jj    ,jk)                       )
(  265)             END DO
(  266)          END DO
(  267) 
(  268)          DO jj = 3, jpj-2
(  269)             DO ji = 2, jpim1   ! vector opt.
(  270)                ! uslop
(  271)                uslp(ji,jj,jk) = zcofu * (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)      &
(  272)                   &                       +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)      &
(  273)                   &                       + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)      &
(  274)                   &                       +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )    &
(  275)                   &                       + 4.*  zwz(ji  ,jj  ,jk)                       )
(  276)                ! vslop
(  277)                vslp(ji,jj,jk) = zcofv * (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)      &
(  278)                   &                       +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)      &
(  279)                   &                       + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)      &
(  280)                   &                       +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )    &
(  281)                   &                       + 4.*  zww(ji,jj    ,jk)                       )
(  282)             END DO
(  283)          END DO
(  284) 
(  285)          ! decrease along coastal boundaries
(  286)          DO jj = 2, jpjm1
(  287)             DO ji = 2, jpim1   ! vector opt.
(  288)                z1u  = ( umask(ji,jj+1,jk) + umask(ji,jj-1,jk) )*.5
(  289)                z1v  = ( vmask(ji+1,jj,jk) + vmask(ji-1,jj,jk) )*.5
(  290)                z1wu = ( umask(ji,jj,jk)   + umask(ji,jj,jk+1) )*.5
(  291)                z1wv = ( vmask(ji,jj,jk)   + vmask(ji,jj,jk+1) )*.5
(  292)                uslp(ji,jj,jk) = uslp(ji,jj,jk) * z1u * z1wu
(  293)                vslp(ji,jj,jk) = vslp(ji,jj,jk) * z1v * z1wv






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 21

(  294)             END DO
(  295)          END DO
(  296) 
(  297) 
(  298)          IF( lk_sco ) THEN
(  299)             ! Add the slope of level surfaces
(  300)             ! -----------------------------------
(  301)             ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
(  302)             ! in inildf, ldfslp never called
(  303)             ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
(  304)             ! is added to the slope of isopycnal surfaces.
(  305)             ! c a u t i o n : minus sign as fsdep has positive value 
(  306)          
(  307)             DO jj = 2, jpjm1
(  308)                DO ji = 2, jpim1   ! vector opt.
(  309)                   uslp(ji,jj,jk) = uslp(ji,jj,jk) - 1. / e1u(ji,jj)   &
(  310)                      &           * ( gdept(jk) - gdept(jk) )
(  311)                   vslp(ji,jj,jk) = vslp(ji,jj,jk) - 1. / e2v(ji,jj)   &
(  312)                      &           * ( gdept(jk) - gdept(jk) )
(  313)                END DO
(  314)             END DO
(  315)          ENDIF
(  316) 
(  317) 
(  318)          ! II. Computation of slopes at w point
(  319)          ! ====================================
(  320)          
(  321)          
(  322)          ! II.1 Slopes of isopycnal surfaces
(  323)          ! ---------------------------------
(  324)          ! wslpi = mij( d/di( prd ) / d/dz( prd )
(  325)          ! wslpj = mij( d/dj( prd ) / d/dz( prd )
(  326)          
(  327)          
(  328)          ! Local vertical density gradient evaluated from N^2
(  329)          !     zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
(  330)          DO jj = 1, jpj
(  331)             DO ji = 1, jpi
(  332)                zwy (ji,jj,jk) = zm05g * pn2 (ji,jj,jk) *   &
(  333)                   &                     ( prd (ji,jj,jk) + prd (ji,jj,jk-1) + 2. )
(  334)             END DO
(  335)          END DO
(  336)          
(  337)          ! Slope at w point
(  338)          DO jj = 2, jpjm1
(  339)             DO ji = 2, jpim1   ! vector opt.
(  340)                ! horizontal density i-gradient at w-points
(  341)                zcoef1 = MAX( zeps, umask(ji-1,jj,jk  )+umask(ji,jj,jk  )    &
(  342)                   &               +umask(ji-1,jj,jk-1)+umask(ji,jj,jk-1) )
(  343)                zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )
(  344)                zai = zcoef1 * (  zgru(ji  ,jj,jk  ) + zgru(ji  ,jj,jk-1)   &
(  345)                   &            + zgru(ji-1,jj,jk-1) + zgru(ji-1,jj,jk  ) ) * tmask (ji,jj,jk)
(  346)                ! horizontal density j-gradient at w-points
(  347)                zcoef2 = MAX( zeps, vmask(ji,jj-1,jk  )+vmask(ji,jj,jk-1)   &
(  348)                   &               +vmask(ji,jj-1,jk-1)+vmask(ji,jj,jk  ) )
(  349)                zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
(  350)                zaj = zcoef2 * (  zgrv(ji,jj  ,jk  ) + zgrv(ji,jj  ,jk-1)   &
(  351)                   &            + zgrv(ji,jj-1,jk-1) + zgrv(ji,jj-1,jk  ) ) * tmask (ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 22

(  352)                ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
(  353)                !                   static instability:
(  354)                !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  355)                zbi = MIN( zwy (ji,jj,jk),- 100.*ABS(zai), -7.e+3/e3w(jk)*ABS(zai) )
(  356)                zbj = MIN( zwy (ji,jj,jk), -100.*ABS(zaj), -7.e+3/e3w(jk)*ABS(zaj) )
(  357)                ! wslpi and wslpj output in zwz and zww, resp.
(  358)                zalpha = MAX(omlmask(ji,jj,jk),omlmask(ji,jj,jk-1))
(  359)                zwz(ji,jj,jk) = ( zai / ( zbi - zeps) * ( 1. - zalpha )   &
(  360)                   &            + zalpha * wslpiml(ji,jj)   &
(  361)                   &            * gdepw(jk) / MAX( hmlp(ji,jj),10. ) )   &
(  362)                   &            * tmask (ji,jj,jk)
(  363)                zww(ji,jj,jk) = ( zaj / ( zbj - zeps) * ( 1. - zalpha )   &
(  364)                   &            + zalpha * wslpjml(ji,jj)   &
(  365)                   &            * gdepw(jk) / MAX( hmlp(ji,jj),10. ) )   &
(  366)                   &            * tmask (ji,jj,jk)
(  367)             END DO
(  368)          END DO
(  369)          !                                             ! ===============
(  370)       END DO                                           !   end of slab
(  371)       !                                                ! ===============
(  372) 
(  373) 
(  374)       ! lateral boundary conditions on zwz and zww
(  375)       CALL lbc_lnk( zwz, 'T', -1. )
(  376)       CALL lbc_lnk( zww, 'T', -1. )
(  377) 
(  378)       !                                                ! ===============
(  379)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  380)          !                                             ! ===============
(  381) 
(  382)          ! Shapiro filter applied in the horizontal direction
(  383) 
(  384)          DO jj = 2, jpjm1, jpj-3   ! row jj=2 and =jpjm1
(  385)             DO ji = 2, jpim1
(  386)                zcofw = tmask(ji,jj,jk)/16.
(  387)                wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
(  388)                   &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
(  389)                   &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
(  390)                   &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
(  391)                   &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw
(  392) 
(  393)                wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
(  394)                   &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
(  395)                   &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
(  396)                   &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
(  397)                   &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
(  398)             END DO
(  399)          END DO
(  400)          
(  401)          DO jj = 3, jpj-2
(  402)             DO ji = 2, jpim1   ! vector opt.
(  403)                zcofw = tmask(ji,jj,jk)/16.
(  404)                wslpi(ji,jj,jk) = (        zwz(ji-1,jj-1,jk) + zwz(ji+1,jj-1,jk)     &
(  405)                   &                +      zwz(ji-1,jj+1,jk) + zwz(ji+1,jj+1,jk)     &
(  406)                   &                + 2.*( zwz(ji  ,jj-1,jk) + zwz(ji-1,jj  ,jk)     &
(  407)                   &                +      zwz(ji+1,jj  ,jk) + zwz(ji  ,jj+1,jk) )   &
(  408)                   &                + 4.*  zwz(ji  ,jj  ,jk)                        ) * zcofw
(  409) 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 23

(  410)                wslpj(ji,jj,jk) = (        zww(ji-1,jj-1,jk) + zww(ji+1,jj-1,jk)     &
(  411)                   &                +      zww(ji-1,jj+1,jk) + zww(ji+1,jj+1,jk)     &
(  412)                   &                + 2.*( zww(ji  ,jj-1,jk) + zww(ji-1,jj  ,jk)     &
(  413)                   &                +      zww(ji+1,jj  ,jk) + zww(ji  ,jj+1,jk) )   &
(  414)                   &                + 4.*  zww(ji  ,jj  ,jk)                        ) * zcofw
(  415)             END DO
(  416)          END DO
(  417)          
(  418)          ! decrease the slope along the boundaries
(  419)          DO jj = 2, jpjm1
(  420)             DO ji = 2, jpim1   ! vector opt.
(  421)                z1u = ( umask(ji,jj,jk) + umask(ji-1,jj,jk) ) *.5
(  422)                z1v = ( vmask(ji,jj,jk) + vmask(ji,jj-1,jk) ) *.5
(  423)                wslpi(ji,jj,jk) = wslpi(ji,jj,jk) * z1u * z1v
(  424)                wslpj(ji,jj,jk) = wslpj(ji,jj,jk) * z1u * z1v
(  425)             END DO
(  426)          END DO
(  427)          
(  428)          IF( lk_sco ) THEN
(  429)          
(  430)             ! Slope of level surfaces
(  431)             ! -----------------------
(  432)             ! 'key_s_coord' defined but not 'key_traldfiso' the computation is done
(  433)             ! in inildf, ldfslp never called
(  434)             ! 'key_s_coord' and 'key_traldfiso' defined, the slope of level surfaces
(  435)             ! is added to the slope of isopycnal surfaces.
(  436)          
(  437)             DO jj = 2, jpjm1
(  438)                DO ji = 2, jpim1   ! vector opt.
(  439)                   wslpi(ji,jj,jk) = wslpi(ji,jj,jk) - 1. / e1t(ji,jj)   &
(  440)                      &                                   * ( gdepw(jk) - gdepw(jk) )
(  441)                   wslpj(ji,jj,jk) = wslpj(ji,jj,jk) - 1. / e2t(ji,jj)   &
(  442)                      &                                   * ( gdepw(jk) - gdepw(jk) )
(  443)                END DO
(  444)             END DO
(  445)          ENDIF
(  446)          
(  447)          ! III. Specific grid points
(  448)          ! -------------------------
(  449) 
(  450)          IF( cp_cfg == "orca" .AND. jp_cfg == 4 ) THEN
(  451)             !                                        ! =======================
(  452)             ! Horizontal diffusion in                !  ORCA_R4 configuration 
(  453)             ! specific area                          ! =======================
(  454)             !
(  455)             !                                             ! Gibraltar Strait
(  456)             ij0 =  50   ;   ij1 =  53
(  457)             ii0 =  69   ;   ii1 =  71   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  458)             ij0 =  51   ;   ij1 =  53
(  459)             ii0 =  68   ;   ii1 =  71   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  460)             ii0 =  69   ;   ii1 =  71   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  461)             ii0 =  69   ;   ii1 =  71   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  462) 
(  463)             !                                             ! Mediterrannean Sea
(  464)             ij0 =  49   ;   ij1 =  56
(  465)             ii0 =  71   ;   ii1 =  90   ;   uslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  466)             ij0 =  50   ;   ij1 =  56
(  467)             ii0 =  70   ;   ii1 =  90   ;   vslp ( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 24

(  468)             ii0 =  71   ;   ii1 =  90   ;   wslpi( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  469)             ii0 =  71   ;   ii1 =  90   ;   wslpj( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , jk ) = 0.e0
(  470)          ENDIF
(  471)          !                                             ! ===============
(  472)       END DO                                           !   end of slab
(  473)       !                                                ! ===============       
(  474) 
(  475)       
(  476)       ! III Lateral boundary conditions on all slopes (uslp , vslp, 
(  477)       ! -------------------------------                wslpi, wslpj )
(  478)       CALL lbc_lnk( uslp , 'U', -1. )
(  479)       CALL lbc_lnk( vslp , 'V', -1. )
(  480)       CALL lbc_lnk( wslpi, 'W', -1. )
(  481)       CALL lbc_lnk( wslpj, 'W', -1. )
(  482) 
(  483)       IF(ln_ctl) THEN
(  484)          CALL prt_ctl(tab3d_1=uslp , clinfo1=' slp  - u : ', tab3d_2=vslp,  clinfo2=' v : ', kdim=jpk)
(  485)          CALL prt_ctl(tab3d_1=wslpi, clinfo1=' slp  - wi: ', tab3d_2=wslpj, clinfo2=' wj: ', kdim=jpk)
(  486)       ENDIF
(  487) 
(  488)    END SUBROUTINE ldf_slp
(  489) 
(  490) 
(  491)    SUBROUTINE ldf_slp_mxl( prd, pn2 )
(  492)       !!----------------------------------------------------------------------
(  493)       !!                  ***  ROUTINE ldf_slp_mxl  ***
(  494)       !! ** Purpose :
(  495)       !!     Compute the slopes of iso-neutral surface (slope of isopycnal
(  496)       !!   surfaces referenced locally) just above the mixed layer.
(  497)       !!
(  498)       !! ** Method :
(  499)       !!      The slope in the i-direction is computed at u- and w-points
(  500)       !!   (uslp, wslpi) and the slope in the j-direction is computed at
(  501)       !!   v- and w-points (vslp, wslpj).
(  502)       !!   They are bounded by 1/100 over the whole ocean, and within the
(  503)       !!   surface layer they are bounded by the distance to the surface
(  504)       !!   ( slope<= depth/l  where l is the length scale of horizontal
(  505)       !!   diffusion (here, aht=2000m2/s ==> l=20km with a typical velocity
(  506)       !!   of 10cm/s)
(  507)       !!
(  508)       !! ** Action :
(  509)       !!      Compute uslp, wslpi, and vslp, wslpj, the i- and  j-slopes
(  510)       !!   of now neutral surfaces at u-, w- and v- w-points, resp.
(  511)       !!
(  512)       !! History :
(  513)       !!   8.1  !  99-10  (A. Jouzeau)  Original code
(  514)       !!   8.5  !  99-10  (G. Madec)  Free form, F90
(  515)       !!----------------------------------------------------------------------
(  516)       !! * Modules used
(  517)       USE oce           , zgru  => ua,  &  ! ua, va used as workspace and set to hor. 
(  518)                           zgrv  => va      ! density gradient in ldf_slp
(  519) 
(  520)       !! * Arguments
(  521)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  522)          prd,                           &  ! in situ density
(  523)          pn2                              ! Brunt-Vaisala frequency (locally ref.)
(  524) 
(  525)       !! * Local declarations






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 25

(  526)       INTEGER  ::   ji, jj, jk             ! dummy loop indices
(  527)       INTEGER  ::   ik, ikm1               ! temporary integers
(  528)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  529)          zwy                               ! temporary workspace
(  530)       REAL(wp) ::   &
(  531)          zeps, zmg, zm05g,              &  ! temporary scalars
(  532)          zcoef1, zcoef2,                &  !    "         "
(  533)          zau, zbu, zav, zbv,            &  !    "         "
(  534)          zai, zbi, zaj, zbj                !    "         "
(  535)       !!----------------------------------------------------------------------
(  536) 
(  537) 
(  538)       ! 0. Local constant initialization
(  539)       ! --------------------------------
(  540) 
(  541)       zeps  =  1.e-20
(  542)       zmg   = -1.0 / grav
(  543)       zm05g = -0.5 / grav
(  544) 
(  545) 
(  546)       uslpml (1,:) = 0.e0      ;      uslpml (jpi,:) = 0.e0
(  547)       vslpml (1,:) = 0.e0      ;      vslpml (jpi,:) = 0.e0
(  548)       wslpiml(1,:) = 0.e0      ;      wslpiml(jpi,:) = 0.e0
(  549)       wslpjml(1,:) = 0.e0      ;      wslpjml(jpi,:) = 0.e0
(  550) 
(  551)       ! surface mixed layer mask
(  552) 
(  553)       ! mask for mixed layer
(  554)       DO jk = 1, jpk
(  555) # 559
(  559)          DO jj = 1, jpj
(  560)             DO ji = 1, jpi
(  561) # 562
(  562)                ! mixed layer interior (mask = 1) and exterior (mask = 0)
(  563)                ik = nmln(ji,jj) - 1
(  564)                IF( jk <= ik ) THEN
(  565)                   omlmask(ji,jj,jk) = 1.e0
(  566)                ELSE
(  567)                   omlmask(ji,jj,jk) = 0.e0
(  568)                ENDIF
(  569) # 570
(  570)             END DO
(  571) # 572
(  572)          END DO
(  573)       END DO
(  574) 
(  575) 
(  576)       ! Slopes of isopycnal surfaces just before bottom of mixed layer
(  577)       ! --------------------------------------------------------------
(  578)       ! uslpml = d/di( prd ) / d/dz( prd )
(  579)       ! vslpml = d/dj( prd ) / d/dz( prd )
(  580) 
(  581)       ! Local vertical density gradient evaluated from N^2
(  582)       ! zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  583) 
(  584)       !-----------------------------------------------------------------------
(  585)       zwy(:,jpj) = 0.e0
(  586)       zwy(jpi,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 26

(  587) # 591
(  591)       DO jj = 1, jpjm1
(  592)          DO ji = 1, jpim1
(  593) # 594
(  594)             ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
(  595)             ! if ik = jpk take jpkm1 values
(  596)             ik = MIN( ik,jpkm1 )
(  597)             zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
(  598)                &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
(  599)                &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
(  600) # 601
(  601)          END DO
(  602) # 603
(  603)       END DO
(  604)       ! lateral boundary conditions on zwy
(  605)       CALL lbc_lnk( zwy, 'U', -1. )
(  606) 
(  607)       ! Slope at u points
(  608) # 612
(  612)       DO jj = 2, jpjm1
(  613)          DO ji = 2, jpim1
(  614) # 615
(  615)             ! horizontal and vertical density gradient at u-points
(  616)             ik = MAX( 1, nmln(ji,jj) , nmln(ji+1,jj) )
(  617)             ik = MIN( ik,jpkm1 )
(  618)             zau = 1./ e1u(ji,jj) * zgru(ji,jj,ik)
(  619)             zbu = 0.5*( zwy(ji,jj) + zwy(ji+1,jj) )
(  620)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  621)             !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  622)             zbu = MIN( zbu, -100.*ABS(zau), -7.e+3/e3t(ik)*ABS(zau) )
(  623)             ! uslpml
(  624)             uslpml (ji,jj) = zau / ( zbu - zeps ) * umask (ji,jj,ik)
(  625) # 626
(  626)          END DO
(  627) # 628
(  628)       END DO
(  629) 
(  630)       ! lateral boundary conditions on uslpml
(  631)       CALL lbc_lnk( uslpml, 'U', -1. )
(  632) 
(  633)       ! Local vertical density gradient evaluated from N^2
(  634)       !     zwy = d/dz(prd)= - ( prd ) / grav * mk(pn2) -- at t point
(  635)       zwy ( :, jpj) = 0.e0
(  636)       zwy ( jpi, :) = 0.e0
(  637) # 641
(  641)       DO jj = 1, jpjm1
(  642)          DO ji = 1, jpim1
(  643) # 644
(  644)             ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
(  645)             ik = MIN( ik,jpkm1 )
(  646)             zwy(ji,jj) = zmg * ( prd(ji,jj,ik) + 1. )   &
(  647)                &             * ( pn2(ji,jj,ik) + pn2(ji,jj,ik+1) )   &
(  648)                &             / MAX( tmask(ji,jj,ik) + tmask (ji,jj,ik+1), 1. )
(  649) # 650
(  650)          END DO
(  651) # 652
(  652)       END DO
(  653) 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 27

(  654)       ! lateral boundary conditions on zwy
(  655)       CALL lbc_lnk( zwy, 'V', -1. )
(  656) 
(  657)       ! Slope at v points
(  658) # 662
(  662)       DO jj = 2, jpjm1
(  663)          DO ji = 2, jpim1
(  664) # 665
(  665)             ! horizontal and vertical density gradient at v-points
(  666)             ik = MAX( 1, nmln(ji,jj) , nmln(ji,jj+1) )
(  667)             ik = MIN( ik,jpkm1 )
(  668)             zav = 1./ e2v(ji,jj) * zgrv(ji,jj,ik)
(  669)             zbv = 0.5*( zwy(ji,jj) + zwy(ji,jj+1) )
(  670)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0
(  671)             !                         kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  672)             zbv = MIN( zbv, -100.*ABS(zav), -7.e+3/e3t(ik)*ABS( zav ) )
(  673)             ! vslpml
(  674)             vslpml (ji,jj) = zav / ( zbv - zeps ) * vmask (ji,jj,ik)
(  675) # 676
(  676)          END DO
(  677) # 678
(  678)       END DO
(  679) 
(  680)       ! lateral boundary conditions on vslpml
(  681)       CALL lbc_lnk( vslpml, 'V', -1. )
(  682) 
(  683)       ! wslpiml = mij( d/di( prd ) / d/dz( prd )
(  684)       ! wslpjml = mij( d/dj( prd ) / d/dz( prd )
(  685) 
(  686) 
(  687)       ! Local vertical density gradient evaluated from N^2
(  688)       ! zwy = d/dz(prd)= - mk ( prd ) / grav * pn2 -- at w point
(  689) # 693
(  693)       DO jj = 1, jpj
(  694)          DO ji = 1, jpi
(  695) # 696
(  696)             ik = nmln(ji,jj)+1
(  697)             ik = MIN( ik,jpk )
(  698)             ikm1 = MAX ( 1, ik-1)
(  699)             zwy (ji,jj) = zm05g * pn2 (ji,jj,ik) *     &
(  700)                &             ( prd (ji,jj,ik) + prd (ji,jj,ikm1) + 2. )
(  701) # 702
(  702)          END DO
(  703) # 704
(  704)       END DO
(  705) 
(  706)       ! Slope at w point
(  707) # 711
(  711)       DO jj = 2, jpjm1
(  712)          DO ji = 2, jpim1
(  713) # 714
(  714)             ik = nmln(ji,jj)+1
(  715)             ik = MIN( ik,jpk )
(  716)             ikm1 = MAX ( 1, ik-1)
(  717)             ! horizontal density i-gradient at w-points
(  718)             zcoef1 = MAX( zeps, umask(ji-1,jj,ik  )+umask(ji,jj,ik  )   &
(  719)                &               +umask(ji-1,jj,ikm1)+umask(ji,jj,ikm1) )
(  720)             zcoef1 = 1. / ( zcoef1 * e1t (ji,jj) )






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 28

(  721)             zai = zcoef1 * (  zgru(ji  ,jj,ik  ) + zgru(ji  ,jj,ikm1)   &
(  722)                &            + zgru(ji-1,jj,ikm1) + zgru(ji-1,jj,ik  ) ) * tmask (ji,jj,ik)
(  723)             ! horizontal density j-gradient at w-points
(  724)             zcoef2 = MAX( zeps, vmask(ji,jj-1,ik  )+vmask(ji,jj,ikm1)    &
(  725)                &               +vmask(ji,jj-1,ikm1)+vmask(ji,jj,ik  ) )
(  726)             zcoef2 = 1.0 / ( zcoef2 *  e2t (ji,jj) )
(  727)             zaj = zcoef2 * (  zgrv(ji,jj  ,ik  ) + zgrv(ji,jj  ,ikm1)   &
(  728)                &            + zgrv(ji,jj-1,ikm1) + zgrv(ji,jj-1,ik  ) ) * tmask (ji,jj,ik)
(  729)             ! bound the slopes: abs(zw.)<= 1/100 and zb..<0.
(  730)             !                   static instability:
(  731)             !                   kxz max= ah slope max =< e1 e3 /(pi**2 2 dt)
(  732)             zbi = MIN ( zwy (ji,jj),- 100.*ABS(zai), -7.e+3/e3w(ik)*ABS(zai) )
(  733)             zbj = MIN ( zwy (ji,jj), -100.*ABS(zaj), -7.e+3/e3w(ik)*ABS(zaj) )
(  734)             ! wslpiml and wslpjml
(  735)             wslpiml (ji,jj) = zai / ( zbi - zeps) * tmask (ji,jj,ik)
(  736)             wslpjml (ji,jj) = zaj / ( zbj - zeps) * tmask (ji,jj,ik)
(  737) # 738
(  738)          END DO
(  739) # 740
(  740)       END DO
(  741) 
(  742)       ! lateral boundary conditions on wslpiml and wslpjml
(  743)       CALL lbc_lnk( wslpiml, 'W', -1. )
(  744)       CALL lbc_lnk( wslpjml, 'W', -1. )
(  745) 
(  746)    END SUBROUTINE ldf_slp_mxl
(  747) 
(  748) 
(  749)    SUBROUTINE ldf_slp_init
(  750)       !!----------------------------------------------------------------------
(  751)       !!                  ***  ROUTINE ldf_slp_init  ***
(  752)       !!
(  753)       !! ** Purpose :   Initialization for the isopycnal slopes computation
(  754)       !!
(  755)       !! ** Method  :   read the nammbf namelist and check the parameter 
(  756)       !!      values called by tra_dmp at the first timestep (nit000)
(  757)       !!
(  758)       !! History :
(  759)       !!  8.5  ! 02-06 (G. Madec) original code
(  760)       !!----------------------------------------------------------------------
(  761)       !! * local declarations
(  762)       INTEGER ::   ji, jj, jk   ! dummy loop indices
(  763)       !!----------------------------------------------------------------------
(  764)       
(  765)       
(  766)       ! Parameter control and print
(  767)       ! ---------------------------
(  768)       IF(lwp) THEN
(  769)          WRITE(numout,*)
(  770)          WRITE(numout,*) 'ldf_slp : direction of lateral mixing'
(  771)          WRITE(numout,*) '~~~~~~~'
(  772)       ENDIF
(  773) 
(  774)       ! Direction of lateral diffusion (tracers and/or momentum)
(  775)       ! ------------------------------
(  776)       ! set the slope to zero (even in s-coordinates)
(  777) 
(  778)       uslp (:,:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 29

(  779)       vslp (:,:,:) = 0.e0
(  780)       wslpi(:,:,:) = 0.e0
(  781)       wslpj(:,:,:) = 0.e0
(  782) 
(  783)       uslpml (:,:) = 0.e0
(  784)       vslpml (:,:) = 0.e0
(  785)       wslpiml(:,:) = 0.e0
(  786)       wslpjml(:,:) = 0.e0
(  787) 
(  788)       IF( ln_traldf_hor .OR. ln_dynldf_hor ) THEN
(  789) 
(  790)          ! geopotential diffusion in s-coordinates on tracers and/or momentum
(  791)          ! The slopes of s-surfaces are computed once (no call to ldfslp in step)
(  792)          ! The slopes for momentum diffusion are i- or j- averaged of those on tracers
(  793) 
(  794)          ! set the slope of diffusion to the slope of s-surfaces
(  795)          !      ( c a u t i o n : minus sign as fsdep has positive value )
(  796)          DO jk = 1, jpk
(  797)             DO jj = 2, jpjm1
(  798)                DO ji = 2, jpim1   ! vector opt.
(  799)                   uslp (ji,jj,jk) = -1. / e1u(ji,jj) * umask(ji,jj,jk)   &
(  800)                      &                               * ( gdept(jk) - gdept(jk) )
(  801)                   vslp (ji,jj,jk) = -1. / e2v(ji,jj) * vmask(ji,jj,jk)   &
(  802)                      &                               * ( gdept(jk) - gdept(jk) )
(  803)                   wslpi(ji,jj,jk) = -1. / e1t(ji,jj) * tmask(ji,jj,jk)   &
(  804)                      &                               * ( gdepw(jk) - gdepw(jk) )
(  805)                   wslpj(ji,jj,jk) = -1. / e2t(ji,jj) * tmask(ji,jj,jk)   &
(  806)                      &                               * ( gdepw(jk) - gdepw(jk) )
(  807)                END DO
(  808)             END DO
(  809)          END DO
(  810) 
(  811)          ! Lateral boundary conditions on the slopes
(  812)          CALL lbc_lnk( uslp , 'U', -1. )
(  813)          CALL lbc_lnk( vslp , 'V', -1. )
(  814)          CALL lbc_lnk( wslpi, 'W', -1. )
(  815)          CALL lbc_lnk( wslpj, 'W', -1. )
(  816)       ENDIF
(  817) 
(  818)    END SUBROUTINE ldf_slp_init
(  819) 
(  820) 
(  821) # 834
(  834)    !!======================================================================
(  835) END MODULE ldfslp
