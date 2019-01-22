


PGF90 (Version     10.5)          09/23/2010  09:50:10      page 1

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

Filename: traadv_cen2.F90

(    1) MODULE traadv_cen2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_cen2  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6)  
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_cen2 : update the tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a 2nd order 
(   10)    !!                  centered finite difference scheme
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and active tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE flxrnf          !
(   18)    USE trabbl          ! advective term in the BBL
(   19)    USE ocfzpt          !
(   20)    USE lib_mpp
(   21)    USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE diaptr          ! poleward transport diagnostics
(   24)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   25)    USE prtctl          ! Print control
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC tra_adv_cen2    ! routine called by step.F90
(   32) 
(   33)    !! * Substitutions
(   34) #  include "domzgr_substitute.h90"
(   35) #  include "vectopt_loop_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_cen2.F90,v 1.11 2005/12/28 09:25:10 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44) #if defined key_autotasking
(   45)    !!----------------------------------------------------------------------
(   46)    !!   'key_autotasking' :      2nd order centered scheme (k- and j-slabs)






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 2

(   47)    !!----------------------------------------------------------------------
(   48) #  include "traadv_cen2_atsk.h90"
(   49) 
(   50) #else
(   51)    !!----------------------------------------------------------------------
(   52)    !!   Default option :             2nd order centered scheme (k-j-i loop)
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55)    SUBROUTINE tra_adv_cen2( kt )
(   56)       !!----------------------------------------------------------------------
(   57)       !!                  ***  ROUTINE tra_adv_cen2  ***
(   58)       !!                 
(   59)       !! ** Purpose :   Compute the now trend due to the advection of tracers
(   60)       !!      and add it to the general trend of passive tracer equations.
(   61)       !!
(   62)       !! ** Method  :   The advection is evaluated by a second order centered
(   63)       !!      scheme using now fields (leap-frog scheme). In specific areas
(   64)       !!      (vicinity of major river mouths, some straits, or where tn is
(   65)       !!      is approaching the freezing point) it is mixed with an upstream
(   66)       !!      scheme for stability reasons.
(   67)       !!        Part 0 : compute the upstream / centered flag
(   68)       !!                 (3D array, zind, defined at T-point (0<zind<1))
(   69)       !!        Part I : horizontal advection
(   70)       !!      * centered flux:
(   71)       !!         * s-coordinate (lk_sco=T) or
(   72)       !!         * z-coordinate with partial steps (lk_zps=T),
(   73)       !!        the vertical scale factors e3. are inside the derivatives:
(   74)       !!               zcenu = e2u*e3u  un  mi(tn)
(   75)       !!               zcenv = e1v*e3v  vn  mj(tn)
(   76)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   77)       !!               zcenu = e2u  un  mi(tn)
(   78)       !!               zcenv = e1v  vn  mj(tn)
(   79)       !!      * upstream flux:
(   80)       !!         * s-coordinate (lk_sco=T) or
(   81)       !!         * z-coordinate with partial steps (lk_zps=T)
(   82)       !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
(   83)       !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
(   84)       !!         * z-coordinate (default key)
(   85)       !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
(   86)       !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
(   87)       !!      * mixed upstream / centered horizontal advection scheme
(   88)       !!               zcofi = max(zind(i+1), zind(i))
(   89)       !!               zcofj = max(zind(j+1), zind(j))
(   90)       !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
(   91)       !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
(   92)       !!      * horizontal advective trend (divergence of the fluxes)
(   93)       !!         * s-coordinate (lk_sco=T) or
(   94)       !!         * z-coordinate with partial steps (lk_zps=T)
(   95)       !!               zta = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
(   96)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   97)       !!               zta = 1/(e1t*e2t) { di-1[zwx] + dj-1[zwy] }
(   98)       !!      * Add this trend now to the general trend of tracer (ta,sa):
(   99)       !!              (ta,sa) = (ta,sa) + ( zta , zsa )
(  100)       !!      * trend diagnostic ('key_trdtra'): the trend is saved
(  101)       !!      for diagnostics. The trends saved is expressed as
(  102)       !!      Uh.gradh(T), (save trend = zta + tn divn).
(  103)       !!         In addition, the advective trend in the two horizontal direc-
(  104)       !!      tion is also re-computed as Uh gradh(T). Indeed hadt+tn divn is






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 3

(  105)       !!      equal to (in s-coordinates, and similarly in z-coord.):
(  106)       !!         zta+tn*divn=1/(e1t*e2t*e3t) { mi-1( e2u*e3u  un  di[tn] )
(  107)       !!                                      +mj-1( e1v*e3v  vn  mj[tn] )  }
(  108)       !!         C A U T I O N : the trend saved is the centered trend only.
(  109)       !!      It doesn't take into account the upstream part of the scheme.
(  110)       !!
(  111)       !!         Part II : vertical advection
(  112)       !!      For temperature (idem for salinity) the advective trend is com-
(  113)       !!      puted as follows :
(  114)       !!            zta = 1/e3t dk+1[ zwz ]
(  115)       !!      where the vertical advective flux, zwz, is given by :
(  116)       !!            zwz = zcofk * zupst + (1-zcofk) * zcent
(  117)       !!      with 
(  118)       !!        zupsv = upstream flux = wn * (tb(k) or tb(k-1) ) [wn>0 or <0]
(  119)       !!        zcenu = centered flux = wn * mk(tn)
(  120)       !!         The surface boundary condition is : 
(  121)       !!      rigid-lid (default option) : zero advective flux
(  122)       !!      free-surf ("key_fresurf_cstvol") : wn(:,:,1) * tn(:,:,1)
(  123)       !!         Add this trend now to the general trend of tracer (ta,sa):
(  124)       !!            (ta,sa) = (ta,sa) + ( zta , zsa )
(  125)       !!         Trend diagnostic ('key_trdtra'): the trend is saved for
(  126)       !!      diagnostics. The trends saved is expressed as :
(  127)       !!             save trend =  w.gradz(T) = zta - tn divn.
(  128)       !!
(  129)       !! ** Action : - update (ta,sa) with the now advective tracer trends
(  130)       !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
(  131)       !!
(  132)       !! History :
(  133)       !!   8.2  !  01-08  (G. Madec, E. Durand)  trahad+trazad = traadv 
(  134)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  135)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  136)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  137)       !!----------------------------------------------------------------------
(  138)       !! * Modules used
(  139)       USE oce                , zwx => ua,  &  ! use ua as workspace
(  140)          &                     zwy => va      ! use va as workspace
(  141) #if defined key_trabbl_adv
(  142)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
(  143)          &         zun, zvn, zwn
(  144) #else
(  145)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(  146)          &                     zvn => vn,  &  ! When no bbl, zvn == vn
(  147)          &                     zwn => wn      ! When no bbl, zwn == wn
(  148) #endif
(  149) 
(  150)  
(  151)       !! * Arguments
(  152)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  153)  
(  154)       !! * Local save
(  155)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(  156)          zbtr2
(  157)  
(  158)       !! * Local declarations
(  159)       INTEGER  ::   ji, jj, jk                 ! dummy loop indices
(  160)       REAL(wp) ::                           &
(  161)          zbtr, zta, zsa, zfui, zfvj,        &  ! temporary scalars
(  162)          zhw, ze3tr, zcofi, zcofj,          &  !    "         "






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 4

(  163)          zupsut, zupsvt, zupsus, zupsvs,    &  !    "         "
(  164)          zfp_ui, zfp_vj, zfm_ui, zfm_vj,    &  !    "         "
(  165)          zcofk, zupst, zupss, zcent,        &  !    "         "
(  166)          zcens, zfp_w, zfm_w,               &  !    "         "
(  167)          zcenut, zcenvt, zcenus, zcenvs        !    "         "
(  168)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  169)          zwz, zww, zind,                    &  ! temporary workspace arrays
(  170)          ztdta, ztdsa                          !    "         "
(  171)       !!----------------------------------------------------------------------
(  172) 
(  173)       IF( kt == nit000 ) THEN
(  174)          IF(lwp) WRITE(numout,*)
(  175)          IF(lwp) WRITE(numout,*) 'tra_adv_cen2 : 2nd order centered advection scheme'
(  176)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
(  177)          IF(lwp) WRITE(numout,*)
(  178)    
(  179)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(  180)       ENDIF
(  181) 
(  182)       ! Save ta and sa trends
(  183)       IF( l_trdtra )   THEN
(  184)          ztdta(:,:,:) = ta(:,:,:) 
(  185)          ztdsa(:,:,:) = sa(:,:,:) 
(  186)          l_adv = 'ce2'
(  187)       ENDIF
(  188) 
(  189) #if defined key_trabbl_adv
(  190)       ! Advective bottom boundary layer 
(  191)       ! -------------------------------
(  192)       zun(:,:,:) = un(:,:,:) - u_bbl(:,:,:)
(  193)       zvn(:,:,:) = vn(:,:,:) - v_bbl(:,:,:)
(  194)       zwn(:,:,:) = wn(:,:,:) + w_bbl(:,:,:)
(  195) #endif
(  196) 
(  197)       ! Upstream / centered scheme indicator
(  198)       ! ------------------------------------
(  199)  
(  200)       DO jk = 1, jpk
(  201)          DO jj = 1, jpj
(  202)             DO ji = 1, jpi
(  203)                zind(ji,jj,jk) =  MAX ( upsrnfh(ji,jj) * upsrnfz(jk),     &  ! changing advection scheme near runoff
(  204)                   &                    upsadv(ji,jj)                     &  ! in the vicinity of some straits
(  205) #if defined key_ice_lim
(  206)                   &                  , tmask(ji,jj,jk)                   &  ! half upstream tracer fluxes
(  207)                   &                  * MAX( 0., SIGN( 1., fzptn(ji,jj)   &  ! if tn < ("freezing"+0.1 )
(  208)                   &                                +0.1-tn(ji,jj,jk) ) ) &
(  209) #endif
(  210)                   &                  )
(  211)             END DO
(  212)          END DO
(  213)       END DO
(  214) 
(  215) 
(  216)       ! I. Horizontal advective fluxes
(  217)       ! ------------------------------
(  218) 
(  219)       ! 1. Second order centered tracer flux at u and v-points
(  220)       ! -------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 5

(  221) 
(  222)       !                                                ! ===============
(  223)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  224)          !                                             ! ===============
(  225)          DO jj = 1, jpjm1
(  226)             DO ji = 1, fs_jpim1   ! vector opt.
(  227)                ! upstream indicator
(  228)                zcofi = MAX( zind(ji+1,jj,jk), zind(ji,jj,jk) )
(  229)                zcofj = MAX( zind(ji,jj+1,jk), zind(ji,jj,jk) )
(  230)                ! volume fluxes * 1/2
(  231) #if defined key_s_coord || defined key_partial_steps
(  232)                zfui = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  233)                zfvj = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  234) #else
(  235)                zfui = 0.5 * e2u(ji,jj) * zun(ji,jj,jk)
(  236)                zfvj = 0.5 * e1v(ji,jj) * zvn(ji,jj,jk)
(  237) #endif
(  238)                ! upstream scheme
(  239)                zfp_ui = zfui + ABS( zfui )
(  240)                zfp_vj = zfvj + ABS( zfvj )
(  241)                zfm_ui = zfui - ABS( zfui )
(  242)                zfm_vj = zfvj - ABS( zfvj )
(  243)                zupsut = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
(  244)                zupsvt = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
(  245)                zupsus = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
(  246)                zupsvs = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
(  247)                ! centered scheme
(  248)                zcenut = zfui * ( tn(ji,jj,jk) + tn(ji+1,jj  ,jk) )
(  249)                zcenvt = zfvj * ( tn(ji,jj,jk) + tn(ji  ,jj+1,jk) )
(  250)                zcenus = zfui * ( sn(ji,jj,jk) + sn(ji+1,jj  ,jk) )
(  251)                zcenvs = zfvj * ( sn(ji,jj,jk) + sn(ji  ,jj+1,jk) )
(  252)                ! mixed centered / upstream scheme
(  253)                zwx(ji,jj,jk) = zcofi * zupsut + (1.-zcofi) * zcenut
(  254)                zwy(ji,jj,jk) = zcofj * zupsvt + (1.-zcofj) * zcenvt
(  255)                zww(ji,jj,jk) = zcofi * zupsus + (1.-zcofi) * zcenus
(  256)                zwz(ji,jj,jk) = zcofj * zupsvs + (1.-zcofj) * zcenvs
(  257)             END DO
(  258)          END DO
(  259) 
(  260)          ! 2. Tracer flux divergence at t-point added to the general trend
(  261)          ! ---------------------------------------------------------------
(  262) 
(  263)          DO jj = 2, jpjm1
(  264)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  265) #if defined key_s_coord || defined key_partial_steps
(  266)                zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  267) #else
(  268)                zbtr = zbtr2(ji,jj) 
(  269) #endif
(  270)                ! horizontal advective trends
(  271)                zta = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
(  272)                   &            + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
(  273)                zsa = - zbtr * (  zww(ji,jj,jk) - zww(ji-1,jj  ,jk)   &
(  274)                   &            + zwz(ji,jj,jk) - zwz(ji  ,jj-1,jk)  )
(  275)                ! add it to the general tracer trends
(  276)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  277)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  278)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 6

(  279)          END DO
(  280)          !                                             ! ===============
(  281)       END DO                                           !   End of slab
(  282)       !                                                ! ===============
(  283) 
(  284)       ! 3. Save the horizontal advective trends for diagnostic
(  285)       ! ------------------------------------------------------
(  286)       IF( l_trdtra )   THEN
(  287)          ! Recompute the hoizontal advection zta & zsa trends computed 
(  288)          ! at the step 2. above in making the difference between the new 
(  289)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  290)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  291)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  292)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  293) 
(  294)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  295) 
(  296)          ! Save the new ta and sa trends
(  297)          ztdta(:,:,:) = ta(:,:,:) 
(  298)          ztdsa(:,:,:) = sa(:,:,:) 
(  299) 
(  300)       ENDIF
(  301) 
(  302)       IF(ln_ctl)   THEN
(  303)           CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 had  - Ta: ', mask1=tmask, &
(  304)              &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  305)       ENDIF
(  306) 
(  307)       ! "zonal" mean advective heat and salt transport 
(  308)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  309) # if defined key_s_coord || defined key_partial_steps
(  310)          pht_adv(:) = ptr_vj( zwy(:,:,:) )
(  311)          pst_adv(:) = ptr_vj( zwz(:,:,:) )
(  312) # else
(  313)          DO jk = 1, jpkm1
(  314)             DO jj = 2, jpjm1
(  315)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  316)                   zwy(ji,jj,jk) = zwy(ji,jj,jk) * fse3v(ji,jj,jk)
(  317)                   zwz(ji,jj,jk) = zwz(ji,jj,jk) * fse3v(ji,jj,jk)
(  318)                END DO
(  319)             END DO
(  320)          END DO
(  321)          pht_adv(:) = ptr_vj( zwy(:,:,:) )
(  322)          pst_adv(:) = ptr_vj( zwz(:,:,:) )
(  323) # endif
(  324)       ENDIF
(  325) 
(  326) 
(  327)       ! II. Vertical advection
(  328)       ! ----------------------
(  329) 
(  330)       ! Bottom value : flux set to zero
(  331)       zwx(:,:,jpk) = 0.e0     ;    zwy(:,:,jpk) = 0.e0
(  332) 
(  333)       ! Surface value
(  334)       IF( lk_dynspg_rl ) THEN
(  335)          ! rigid lid : flux set to zero
(  336)          zwx(:,:, 1 ) = 0.e0    ;    zwy(:,:, 1 ) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 7

(  337)       ELSE
(  338)          ! free surface
(  339)          zwx(:,:, 1 ) = zwn(:,:,1) * tn(:,:,1)
(  340)          zwy(:,:, 1 ) = zwn(:,:,1) * sn(:,:,1)
(  341)       ENDIF
(  342) 
(  343)       ! 1. Vertical advective fluxes
(  344)       ! ----------------------------
(  345) 
(  346)       ! Second order centered tracer flux at w-point
(  347) 
(  348)       DO jk = 2, jpk
(  349)          DO jj = 2, jpjm1
(  350)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  351)                ! upstream indicator
(  352)                zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) )
(  353)                ! velocity * 1/2
(  354)                zhw = 0.5 * zwn(ji,jj,jk)
(  355)                ! upstream scheme
(  356)                zfp_w = zhw + ABS( zhw )
(  357)                zfm_w = zhw - ABS( zhw )
(  358)                zupst = zfp_w * tb(ji,jj,jk) + zfm_w * tb(ji,jj,jk-1)
(  359)                zupss = zfp_w * sb(ji,jj,jk) + zfm_w * sb(ji,jj,jk-1)
(  360)                ! centered scheme
(  361)                zcent = zhw * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) )
(  362)                zcens = zhw * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) )
(  363)                ! mixed centered / upstream scheme
(  364)                zwx(ji,jj,jk) = zcofk * zupst + (1.-zcofk) * zcent
(  365)                zwy(ji,jj,jk) = zcofk * zupss + (1.-zcofk) * zcens
(  366)             END DO
(  367)          END DO
(  368)       END DO
(  369) 
(  370) 
(  371)       ! 2. Tracer flux divergence at t-point added to the general trend
(  372)       ! -------------------------
(  373)   
(  374)       DO jk = 1, jpkm1
(  375)          DO jj = 2, jpjm1
(  376)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  377)                ze3tr = 1. / fse3t(ji,jj,jk)
(  378)                ! vertical advective trends
(  379)                zta = - ze3tr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
(  380)                zsa = - ze3tr * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) )
(  381)                ! add it to the general tracer trends
(  382)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
(  383)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
(  384)             END DO
(  385)          END DO
(  386)       END DO
(  387) 
(  388)       ! 3. Save the vertical advective trends for diagnostic
(  389)       ! ----------------------------------------------------
(  390)       IF( l_trdtra )   THEN
(  391)          ! Recompute the vertical advection zta & zsa trends computed 
(  392)          ! at the step 2. above in making the difference between the new 
(  393)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  394)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 8

(  395)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  396)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  397) 
(  398)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
(  399)       ENDIF
(  400) 
(  401)       IF(ln_ctl)   THEN
(  402)           CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 zad  - Ta: ', mask1=tmask, &
(  403)              &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  404)       ENDIF
(  405) 
(  406)    END SUBROUTINE tra_adv_cen2
(  407) 
(  408) #endif
(  409) 
(  410)    !!======================================================================
(  411) END MODULE traadv_cen2















































PGF90 (Version     10.5)          09/23/2010  09:50:10      page 9

(    1) # 1 "traadv_cen2.F90"
(    1) MODULE traadv_cen2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  traadv_cen2  ***
(    4)    !! Ocean active tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6)  
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_cen2 : update the tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a 2nd order 
(   10)    !!                  centered finite difference scheme
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and active tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE flxrnf          !
(   18)    USE trabbl          ! advective term in the BBL
(   19)    USE ocfzpt          !
(   20)    USE lib_mpp
(   21)    USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE diaptr          ! poleward transport diagnostics
(   24)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   25)    USE prtctl          ! Print control
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC tra_adv_cen2    ! routine called by step.F90
(   32) 
(   33)    !! * Substitutions
(   34) # 1 "./domzgr_substitute.h90"
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
(   48) # 35 "traadv_cen2.F90"
(   35) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 10

(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 36 "traadv_cen2.F90"
(   36) # 36
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_cen2.F90,v 1.11 2005/12/28 09:25:10 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44) # 51
(   51)    !!----------------------------------------------------------------------
(   52)    !!   Default option :             2nd order centered scheme (k-j-i loop)
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55)    SUBROUTINE tra_adv_cen2( kt )
(   56)       !!----------------------------------------------------------------------
(   57)       !!                  ***  ROUTINE tra_adv_cen2  ***
(   58)       !!                 
(   59)       !! ** Purpose :   Compute the now trend due to the advection of tracers
(   60)       !!      and add it to the general trend of passive tracer equations.
(   61)       !!
(   62)       !! ** Method  :   The advection is evaluated by a second order centered
(   63)       !!      scheme using now fields (leap-frog scheme). In specific areas
(   64)       !!      (vicinity of major river mouths, some straits, or where tn is
(   65)       !!      is approaching the freezing point) it is mixed with an upstream
(   66)       !!      scheme for stability reasons.
(   67)       !!        Part 0 : compute the upstream / centered flag
(   68)       !!                 (3D array, zind, defined at T-point (0<zind<1))
(   69)       !!        Part I : horizontal advection
(   70)       !!      * centered flux:
(   71)       !!         * s-coordinate (lk_sco=T) or
(   72)       !!         * z-coordinate with partial steps (lk_zps=T),
(   73)       !!        the vertical scale factors e3. are inside the derivatives:
(   74)       !!               zcenu = e2u*e3u  un  mi(tn)
(   75)       !!               zcenv = e1v*e3v  vn  mj(tn)
(   76)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   77)       !!               zcenu = e2u  un  mi(tn)
(   78)       !!               zcenv = e1v  vn  mj(tn)
(   79)       !!      * upstream flux:
(   80)       !!         * s-coordinate (lk_sco=T) or
(   81)       !!         * z-coordinate with partial steps (lk_zps=T)
(   82)       !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
(   83)       !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
(   84)       !!         * z-coordinate (default key)
(   85)       !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
(   86)       !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
(   87)       !!      * mixed upstream / centered horizontal advection scheme
(   88)       !!               zcofi = max(zind(i+1), zind(i))
(   89)       !!               zcofj = max(zind(j+1), zind(j))
(   90)       !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
(   91)       !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
(   92)       !!      * horizontal advective trend (divergence of the fluxes)
(   93)       !!         * s-coordinate (lk_sco=T) or






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 11

(   94)       !!         * z-coordinate with partial steps (lk_zps=T)
(   95)       !!               zta = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
(   96)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   97)       !!               zta = 1/(e1t*e2t) { di-1[zwx] + dj-1[zwy] }
(   98)       !!      * Add this trend now to the general trend of tracer (ta,sa):
(   99)       !!              (ta,sa) = (ta,sa) + ( zta , zsa )
(  100)       !!      * trend diagnostic ('key_trdtra'): the trend is saved
(  101)       !!      for diagnostics. The trends saved is expressed as
(  102)       !!      Uh.gradh(T), (save trend = zta + tn divn).
(  103)       !!         In addition, the advective trend in the two horizontal direc-
(  104)       !!      tion is also re-computed as Uh gradh(T). Indeed hadt+tn divn is
(  105)       !!      equal to (in s-coordinates, and similarly in z-coord.):
(  106)       !!         zta+tn*divn=1/(e1t*e2t*e3t) { mi-1( e2u*e3u  un  di[tn] )
(  107)       !!                                      +mj-1( e1v*e3v  vn  mj[tn] )  }
(  108)       !!         C A U T I O N : the trend saved is the centered trend only.
(  109)       !!      It doesn't take into account the upstream part of the scheme.
(  110)       !!
(  111)       !!         Part II : vertical advection
(  112)       !!      For temperature (idem for salinity) the advective trend is com-
(  113)       !!      puted as follows :
(  114)       !!            zta = 1/e3t dk+1[ zwz ]
(  115)       !!      where the vertical advective flux, zwz, is given by :
(  116)       !!            zwz = zcofk * zupst + (1-zcofk) * zcent
(  117)       !!      with 
(  118)       !!        zupsv = upstream flux = wn * (tb(k) or tb(k-1) ) [wn>0 or <0]
(  119)       !!        zcenu = centered flux = wn * mk(tn)
(  120)       !!         The surface boundary condition is : 
(  121)       !!      rigid-lid (default option) : zero advective flux
(  122)       !!      free-surf ("key_fresurf_cstvol") : wn(:,:,1) * tn(:,:,1)
(  123)       !!         Add this trend now to the general trend of tracer (ta,sa):
(  124)       !!            (ta,sa) = (ta,sa) + ( zta , zsa )
(  125)       !!         Trend diagnostic ('key_trdtra'): the trend is saved for
(  126)       !!      diagnostics. The trends saved is expressed as :
(  127)       !!             save trend =  w.gradz(T) = zta - tn divn.
(  128)       !!
(  129)       !! ** Action : - update (ta,sa) with the now advective tracer trends
(  130)       !!             - save the trends in (ttrdh,strdh) ('key_trdtra')
(  131)       !!
(  132)       !! History :
(  133)       !!   8.2  !  01-08  (G. Madec, E. Durand)  trahad+trazad = traadv 
(  134)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  135)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  136)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  137)       !!----------------------------------------------------------------------
(  138)       !! * Modules used
(  139)       USE oce                , zwx => ua,  &  ! use ua as workspace
(  140)          &                     zwy => va      ! use va as workspace
(  141) # 145
(  145)       USE oce                , zun => un,  &  ! When no bbl, zun == un
(  146)          &                     zvn => vn,  &  ! When no bbl, zvn == vn
(  147)          &                     zwn => wn      ! When no bbl, zwn == wn
(  148) 
(  149) # 150
(  150)  
(  151)       !! * Arguments
(  152)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  153)  
(  154)       !! * Local save






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 12

(  155)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(  156)          zbtr2
(  157)  
(  158)       !! * Local declarations
(  159)       INTEGER  ::   ji, jj, jk                 ! dummy loop indices
(  160)       REAL(wp) ::                           &
(  161)          zbtr, zta, zsa, zfui, zfvj,        &  ! temporary scalars
(  162)          zhw, ze3tr, zcofi, zcofj,          &  !    "         "
(  163)          zupsut, zupsvt, zupsus, zupsvs,    &  !    "         "
(  164)          zfp_ui, zfp_vj, zfm_ui, zfm_vj,    &  !    "         "
(  165)          zcofk, zupst, zupss, zcent,        &  !    "         "
(  166)          zcens, zfp_w, zfm_w,               &  !    "         "
(  167)          zcenut, zcenvt, zcenus, zcenvs        !    "         "
(  168)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  169)          zwz, zww, zind,                    &  ! temporary workspace arrays
(  170)          ztdta, ztdsa                          !    "         "
(  171)       !!----------------------------------------------------------------------
(  172) 
(  173)       IF( kt == nit000 ) THEN
(  174)          IF(lwp) WRITE(numout,*)
(  175)          IF(lwp) WRITE(numout,*) 'tra_adv_cen2 : 2nd order centered advection scheme'
(  176)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
(  177)          IF(lwp) WRITE(numout,*)
(  178)    
(  179)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(  180)       ENDIF
(  181) 
(  182)       ! Save ta and sa trends
(  183)       IF( l_trdtra )   THEN
(  184)          ztdta(:,:,:) = ta(:,:,:) 
(  185)          ztdsa(:,:,:) = sa(:,:,:) 
(  186)          l_adv = 'ce2'
(  187)       ENDIF
(  188) 
(  189) 
(  190) # 197
(  197)       ! Upstream / centered scheme indicator
(  198)       ! ------------------------------------
(  199)  
(  200)       DO jk = 1, jpk
(  201)          DO jj = 1, jpj
(  202)             DO ji = 1, jpi
(  203)                zind(ji,jj,jk) =  MAX ( upsrnfh(ji,jj) * upsrnfz(jk),     &  ! changing advection scheme near runoff
(  204)                   &                    upsadv(ji,jj)                     &  ! in the vicinity of some straits
(  206)                   &                  , tmask(ji,jj,jk)                   &  ! half upstream tracer fluxes
(  207)                   &                  * MAX( 0., SIGN( 1., fzptn(ji,jj)   &  ! if tn < ("freezing"+0.1 )
(  208)                   &                                +0.1-tn(ji,jj,jk) ) ) &
(  210)                   &                  )
(  211)             END DO
(  212)          END DO
(  213)       END DO
(  214) 
(  215) 
(  216)       ! I. Horizontal advective fluxes
(  217)       ! ------------------------------
(  218) 
(  219)       ! 1. Second order centered tracer flux at u and v-points
(  220)       ! -------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 13

(  221) 
(  222)       !                                                ! ===============
(  223)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  224)          !                                             ! ===============
(  225)          DO jj = 1, jpjm1
(  226)             DO ji = 1, jpim1   ! vector opt.
(  227)                ! upstream indicator
(  228)                zcofi = MAX( zind(ji+1,jj,jk), zind(ji,jj,jk) )
(  229)                zcofj = MAX( zind(ji,jj+1,jk), zind(ji,jj,jk) )
(  230)                ! volume fluxes * 1/2
(  231) # 235
(  235)                zfui = 0.5 * e2u(ji,jj) * zun(ji,jj,jk)
(  236)                zfvj = 0.5 * e1v(ji,jj) * zvn(ji,jj,jk)
(  237) # 238
(  238)                ! upstream scheme
(  239)                zfp_ui = zfui + ABS( zfui )
(  240)                zfp_vj = zfvj + ABS( zfvj )
(  241)                zfm_ui = zfui - ABS( zfui )
(  242)                zfm_vj = zfvj - ABS( zfvj )
(  243)                zupsut = zfp_ui * tb(ji,jj,jk) + zfm_ui * tb(ji+1,jj  ,jk)
(  244)                zupsvt = zfp_vj * tb(ji,jj,jk) + zfm_vj * tb(ji  ,jj+1,jk)
(  245)                zupsus = zfp_ui * sb(ji,jj,jk) + zfm_ui * sb(ji+1,jj  ,jk)
(  246)                zupsvs = zfp_vj * sb(ji,jj,jk) + zfm_vj * sb(ji  ,jj+1,jk)
(  247)                ! centered scheme
(  248)                zcenut = zfui * ( tn(ji,jj,jk) + tn(ji+1,jj  ,jk) )
(  249)                zcenvt = zfvj * ( tn(ji,jj,jk) + tn(ji  ,jj+1,jk) )
(  250)                zcenus = zfui * ( sn(ji,jj,jk) + sn(ji+1,jj  ,jk) )
(  251)                zcenvs = zfvj * ( sn(ji,jj,jk) + sn(ji  ,jj+1,jk) )
(  252)                ! mixed centered / upstream scheme
(  253)                zwx(ji,jj,jk) = zcofi * zupsut + (1.-zcofi) * zcenut
(  254)                zwy(ji,jj,jk) = zcofj * zupsvt + (1.-zcofj) * zcenvt
(  255)                zww(ji,jj,jk) = zcofi * zupsus + (1.-zcofi) * zcenus
(  256)                zwz(ji,jj,jk) = zcofj * zupsvs + (1.-zcofj) * zcenvs
(  257)             END DO
(  258)          END DO
(  259) 
(  260)          ! 2. Tracer flux divergence at t-point added to the general trend
(  261)          ! ---------------------------------------------------------------
(  262) 
(  263)          DO jj = 2, jpjm1
(  264)             DO ji = 2, jpim1   ! vector opt.
(  265) # 268
(  268)                zbtr = zbtr2(ji,jj) 
(  269) # 270
(  270)                ! horizontal advective trends
(  271)                zta = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
(  272)                   &            + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
(  273)                zsa = - zbtr * (  zww(ji,jj,jk) - zww(ji-1,jj  ,jk)   &
(  274)                   &            + zwz(ji,jj,jk) - zwz(ji  ,jj-1,jk)  )
(  275)                ! add it to the general tracer trends
(  276)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  277)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  278)             END DO
(  279)          END DO
(  280)          !                                             ! ===============
(  281)       END DO                                           !   End of slab
(  282)       !                                                ! ===============
(  283) 






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 14

(  284)       ! 3. Save the horizontal advective trends for diagnostic
(  285)       ! ------------------------------------------------------
(  286)       IF( l_trdtra )   THEN
(  287)          ! Recompute the hoizontal advection zta & zsa trends computed 
(  288)          ! at the step 2. above in making the difference between the new 
(  289)          ! trends and the previous one ta()/sa - ztdta()/ztdsa() and add
(  290)          ! the term tn()/sn()*hdivn() to recover the Uh gradh(T/S) trends
(  291)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) + tn(:,:,:) * hdivn(:,:,:) 
(  292)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) + sn(:,:,:) * hdivn(:,:,:)
(  293) 
(  294)          CALL trd_mod(ztdta, ztdsa, jpttdlad, 'TRA', kt)
(  295) 
(  296)          ! Save the new ta and sa trends
(  297)          ztdta(:,:,:) = ta(:,:,:) 
(  298)          ztdsa(:,:,:) = sa(:,:,:) 
(  299) 
(  300)       ENDIF
(  301) 
(  302)       IF(ln_ctl)   THEN
(  303)           CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 had  - Ta: ', mask1=tmask, &
(  304)              &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  305)       ENDIF
(  306) 
(  307)       ! "zonal" mean advective heat and salt transport 
(  308)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  309) # 313
(  313)          DO jk = 1, jpkm1
(  314)             DO jj = 2, jpjm1
(  315)                DO ji = 2, jpim1   ! vector opt.
(  316)                   zwy(ji,jj,jk) = zwy(ji,jj,jk) * e3t(jk)
(  317)                   zwz(ji,jj,jk) = zwz(ji,jj,jk) * e3t(jk)
(  318)                END DO
(  319)             END DO
(  320)          END DO
(  321)          pht_adv(:) = ptr_vj( zwy(:,:,:) )
(  322)          pst_adv(:) = ptr_vj( zwz(:,:,:) )
(  323) # 324
(  324)       ENDIF
(  325) 
(  326) 
(  327)       ! II. Vertical advection
(  328)       ! ----------------------
(  329) 
(  330)       ! Bottom value : flux set to zero
(  331)       zwx(:,:,jpk) = 0.e0     ;    zwy(:,:,jpk) = 0.e0
(  332) 
(  333)       ! Surface value
(  334)       IF( lk_dynspg_rl ) THEN
(  335)          ! rigid lid : flux set to zero
(  336)          zwx(:,:, 1 ) = 0.e0    ;    zwy(:,:, 1 ) = 0.e0
(  337)       ELSE
(  338)          ! free surface
(  339)          zwx(:,:, 1 ) = zwn(:,:,1) * tn(:,:,1)
(  340)          zwy(:,:, 1 ) = zwn(:,:,1) * sn(:,:,1)
(  341)       ENDIF
(  342) 
(  343)       ! 1. Vertical advective fluxes
(  344)       ! ----------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 15

(  345) 
(  346)       ! Second order centered tracer flux at w-point
(  347) 
(  348)       DO jk = 2, jpk
(  349)          DO jj = 2, jpjm1
(  350)             DO ji = 2, jpim1   ! vector opt.
(  351)                ! upstream indicator
(  352)                zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) )
(  353)                ! velocity * 1/2
(  354)                zhw = 0.5 * zwn(ji,jj,jk)
(  355)                ! upstream scheme
(  356)                zfp_w = zhw + ABS( zhw )
(  357)                zfm_w = zhw - ABS( zhw )
(  358)                zupst = zfp_w * tb(ji,jj,jk) + zfm_w * tb(ji,jj,jk-1)
(  359)                zupss = zfp_w * sb(ji,jj,jk) + zfm_w * sb(ji,jj,jk-1)
(  360)                ! centered scheme
(  361)                zcent = zhw * ( tn(ji,jj,jk) + tn(ji,jj,jk-1) )
(  362)                zcens = zhw * ( sn(ji,jj,jk) + sn(ji,jj,jk-1) )
(  363)                ! mixed centered / upstream scheme
(  364)                zwx(ji,jj,jk) = zcofk * zupst + (1.-zcofk) * zcent
(  365)                zwy(ji,jj,jk) = zcofk * zupss + (1.-zcofk) * zcens
(  366)             END DO
(  367)          END DO
(  368)       END DO
(  369) 
(  370) 
(  371)       ! 2. Tracer flux divergence at t-point added to the general trend
(  372)       ! -------------------------
(  373)   
(  374)       DO jk = 1, jpkm1
(  375)          DO jj = 2, jpjm1
(  376)             DO ji = 2, jpim1   ! vector opt.
(  377)                ze3tr = 1. / e3t(jk)
(  378)                ! vertical advective trends
(  379)                zta = - ze3tr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
(  380)                zsa = - ze3tr * ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) )
(  381)                ! add it to the general tracer trends
(  382)                ta(ji,jj,jk) =  ta(ji,jj,jk) + zta
(  383)                sa(ji,jj,jk) =  sa(ji,jj,jk) + zsa
(  384)             END DO
(  385)          END DO
(  386)       END DO
(  387) 
(  388)       ! 3. Save the vertical advective trends for diagnostic
(  389)       ! ----------------------------------------------------
(  390)       IF( l_trdtra )   THEN
(  391)          ! Recompute the vertical advection zta & zsa trends computed 
(  392)          ! at the step 2. above in making the difference between the new 
(  393)          ! trends and the previous one: ta()/sa - ztdta()/ztdsa() and substract
(  394)          ! the term tn()/sn()*hdivn() to recover the W gradz(T/S) trends
(  395)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:) - tn(:,:,:) * hdivn(:,:,:) 
(  396)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:) - sn(:,:,:) * hdivn(:,:,:)
(  397) 
(  398)          CALL trd_mod(ztdta, ztdsa, jpttdzad, 'TRA', kt)
(  399)       ENDIF
(  400) 
(  401)       IF(ln_ctl)   THEN
(  402)           CALL prt_ctl(tab3d_1=ta, clinfo1=' centered2 zad  - Ta: ', mask1=tmask, &






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 16

(  403)              &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  404)       ENDIF
(  405) 
(  406)    END SUBROUTINE tra_adv_cen2
(  407) 
(  408) 
(  409) # 410
(  410)    !!======================================================================
(  411) END MODULE traadv_cen2
