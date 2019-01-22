


PGF90 (Version     10.5)          09/23/2010  09:49:18      page 1

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

Filename: trcadv_cen2.F90

(    1) MODULE trcadv_cen2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_cen2  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective tracer trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_adv_cen2 : update the tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a 2nd order 
(   10)    !!                  centered finite difference scheme
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce_trc             ! ocean dynamics and active tracers variables
(   14)    USE trc                 ! ocean passive tracers variables
(   15)    USE trcbbl              ! advective passive tracers in the BBL
(   16)    USE prtctl_trc
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC trc_adv_cen2    ! routine called by trcstp.F90
(   23) 
(   24)    !! * Substitutions
(   25) #  include "passivetrc_substitute.h90"
(   26)    !!----------------------------------------------------------------------
(   27)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_cen2.F90,v 1.11 2006/04/10 15:38:53 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    !!----------------------------------------------------------------------
(   35)    !!   Default option :             2nd order centered scheme (k-j-i loop)
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38)    SUBROUTINE trc_adv_cen2( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                  ***  ROUTINE trc_adv_cen2  ***
(   41)       !!                 
(   42)       !! ** Purpose :   Compute the now trend due to the advection of tracers
(   43)       !!      and add it to the general trend of passive tracer equations.
(   44)       !!
(   45)       !! ** Method  :   The advection is evaluated by a second order centered
(   46)       !!      scheme using now fields (leap-frog scheme). In specific areas






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 2

(   47)       !!      (vicinity of major river mouths, some straits, or where tn is
(   48)       !!      is approaching the freezing point) it is mixed with an upstream
(   49)       !!      scheme for stability reasons.
(   50)       !!        Part 0 : compute the upstream / centered flag
(   51)       !!                 (3D array, zind, defined at T-point (0<zind<1))
(   52)       !!        Part I : horizontal advection
(   53)       !!      * centered flux:
(   54)       !!         * s-coordinate (lk_sco=T) or
(   55)       !!         * z-coordinate with partial steps (lk_zps=T),
(   56)       !!        the vertical scale factors e3. are inside the derivatives:
(   57)       !!               zcenu = e2u*e3u  un  mi(tn)
(   58)       !!               zcenv = e1v*e3v  vn  mj(tn)
(   59)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   60)       !!               zcenu = e2u  un  mi(tn)
(   61)       !!               zcenv = e1v  vn  mj(tn)
(   62)       !!      * horizontal advective trend (divergence of the fluxes)
(   63)       !!         * s-coordinate (lk_sco=T) or
(   64)       !!         * z-coordinate with partial steps (lk_zps=T)
(   65)       !!               ztra = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
(   66)       !!         * z-coordinate (default key), e3t=e3u=e3v:
(   67)       !!               ztra = 1/(e1t*e2t) { di-1[zwx] + dj-1[zwy] }
(   68)       !!      * Add this trend now to the general trend of tracer tra:
(   69)       !!              tra = tra + ztra
(   70)       !!      * trend diagnostic ('key_trc_diatrd'): the trend is saved
(   71)       !!      for diagnostics. The trends saved is expressed as
(   72)       !!      Uh.gradh(T)
(   73)       !!
(   74)       !!         Part II : vertical advection
(   75)       !!      For any tracer  the advective trend is computed as follows :
(   76)       !!            ztra = 1/e3t dk+1[ zwz ]
(   77)       !!      where the vertical advective flux, zwz, is given by :
(   78)       !!            zwz = zcofk * zupst + (1-zcofk) * zcent
(   79)       !!      with 
(   80)       !!        zupsv = upstream flux = wn * (trb(k) or trb(k-1) ) [wn>0 or <0]
(   81)       !!        zcenu = centered flux = wn * mk(trn)
(   82)       !!         The surface boundary condition is : 
(   83)       !!      rigid-lid (default option) : zero advective flux
(   84)       !!      free-surf ("key_fresurf_cstvol") : wn(:,:,1) * trn(:,:,1)
(   85)       !!         Add this trend now to the general trend of tracer tra :
(   86)       !!            tra = tra + ztra
(   87)       !!         Trend diagnostic ('key_trc_diatrd'): the trend is saved for
(   88)       !!      diagnostics. The trends saved is expressed as :
(   89)       !!             save trend =  w.gradz(T) = ztra - trn divn.
(   90)       !!
(   91)       !! ** Action : - update tra with the now advective tracer trends
(   92)       !!             - save the trends in trtrd ('key_trc_diatrd')
(   93)       !!
(   94)       !! History :
(   95)       !!   8.2  !  01-08  (M-A Filiberti, and M.Levy)  trahad+trazad = traadv 
(   96)       !!   8.5  !  02-06  (G. Madec, C. Ethe)  F90: Free form and module
(   97)       !!----------------------------------------------------------------------
(   98)       !! * Modules used
(   99)       USE oce_trc            , zwx => ua,  &  ! use ua as workspace
(  100)          &                     zwy => va      ! use va as workspace
(  101) #if defined key_trcbbl_adv
(  102)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &  ! temporary arrays
(  103)          &         zun, zvn, zwn
(  104) #else






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 3

(  105)       USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
(  106)          &                     zvn => vn,  &  ! When no bbl, zvn == vn
(  107)          &                     zwn => wn      ! When no bbl, zwn == wn
(  108) #endif
(  109) !!DBG
(  110)       USE obc_oce
(  111) 
(  112)  
(  113)       !! * Arguments
(  114)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  115)  
(  116)       !! * Local save
(  117)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(  118)          zbtr2
(  119)  
(  120)       !! * Local declarations
(  121)       INTEGER  ::   ji, jj, jk, jn             ! dummy loop indices
(  122)       REAL(wp) ::                           &
(  123)          zbtr, ztra, zfui, zfvj,            &  ! temporary scalars
(  124)          zhw, ze3tr, zcofi, zcofj,          &  !    "         "
(  125)          zupsut, zupsvt,                    &  !    "         "
(  126)          zfp_ui, zfp_vj, zfm_ui, zfm_vj,    &  !    "         "
(  127)          zcofk, zupst, zcent,               &  !    "         "
(  128)          zfp_w, zfm_w,                      &  !    "         "
(  129)          zcenut, zcenvt                        ! 
(  130) 
(  131)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  132)          zind                              ! temporary workspace arrays
(  133) #if defined key_trc_diatrd
(  134)       REAL(wp) ::                           &
(  135)          ztai, ztaj,                        &  ! temporary scalars
(  136)          zfui1, zfvj1                          !    "         "
(  137) #endif
(  138)       CHARACTER (len=22) :: charout
(  139) !!DB
(  140)       INTEGER :: ii, ij, ik, ilast, ifirst
(  141)       !!----------------------------------------------------------------------
(  142) 
(  143)       IF( kt == nittrc000 ) THEN
(  144)          IF(lwp) WRITE(numout,*)
(  145)          IF(lwp) WRITE(numout,*) 'trc_adv_cen2 : 2nd order centered advection scheme'
(  146)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
(  147)          IF(lwp) WRITE(numout,*)
(  148)    
(  149)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(  150)       ENDIF
(  151) 
(  152) !!DBG
(  153)       zwx(:,:,:) = 0.0; zwy(:,:,:) = 0.0
(  154) 
(  155) 
(  156) #if defined key_trcbbl_adv
(  157) 
(  158)       ! Advective bottom boundary layer 
(  159)       ! -------------------------------
(  160)       zun(:,:,:) = un(:,:,:) - u_trc_bbl(:,:,:)
(  161)       zvn(:,:,:) = vn(:,:,:) - v_trc_bbl(:,:,:)
(  162)       zwn(:,:,:) = wn(:,:,:) + w_trc_bbl(:,:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 4

(  163) #endif
(  164) 
(  165) 
(  166)       ! Upstream / centered scheme indicator
(  167)       ! ------------------------------------
(  168)  
(  169)       DO jk = 1, jpk
(  170)          DO jj = 1, jpj
(  171)             DO ji = 1, jpi
(  172)                zind(ji,jj,jk) =  MAX ( upsrnfh(ji,jj) * upsrnfz(jk),     &  ! changing advection scheme near runoff
(  173)                   &                    upsadv(ji,jj)                     &  ! in the vicinity of some straits
(  174) #if defined key_ice_lim
(  175)                   &                  , tmask(ji,jj,jk)                   &  ! half upstream tracer fluxes
(  176)                   &                  * MAX( 0., SIGN( 1., fzptn(ji,jj)   &  ! if tn < ("freezing"+0.1 )
(  177)                   &                                +0.1-tn(ji,jj,jk) ) ) &
(  178) #endif
(  179)                   &                  )
(  180)             END DO
(  181)          END DO
(  182)       END DO
(  183) 
(  184) !!DBG: change indices on east & west boundaries
(  185)       ilast = fs_jpim1
(  186) !      do ji = nie0, nie1
(  187) !         if(fs_jpim1 > nie0) ilast = nie0
(  188) !      enddo
(  189)       ifirst = 1
(  190) !      do ji = niw0, niw1
(  191) !         ifirst = 2
(  192) !      enddo
(  193) 
(  194) 
(  195)       DO jn = 1, jptra
(  196)          ! I. Horizontal advective fluxes
(  197)          ! ------------------------------
(  198)          
(  199)          ! Second order centered tracer flux at u and v-points
(  200)          
(  201)          !                                                ! ===============
(  202)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  203)             !                                             ! ===============
(  204)             DO jj = 1, jpjm1
(  205) !               DO ji = 1, ilast   ! vector opt.
(  206) !!DBG
(  207)                DO ji = ifirst, ilast   ! vector opt.
(  208)                ! upstream indicator
(  209)                   zcofi = MAX( zind(ji+1,jj,jk), zind(ji,jj,jk) )
(  210)                   zcofj = MAX( zind(ji,jj+1,jk), zind(ji,jj,jk) )
(  211)                   ! volume fluxes * 1/2
(  212) #if defined key_s_coord || defined key_partial_steps
(  213)                   zfui = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  214)                   zfvj = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  215) #else
(  216)                   zfui = 0.5 * e2u(ji,jj) * zun(ji,jj,jk)
(  217)                   zfvj = 0.5 * e1v(ji,jj) * zvn(ji,jj,jk)
(  218) #endif
(  219)                ! upstream scheme
(  220)                   zfp_ui = zfui + ABS( zfui )






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 5

(  221)                   zfp_vj = zfvj + ABS( zfvj )
(  222)                   zfm_ui = zfui - ABS( zfui )
(  223)                   zfm_vj = zfvj - ABS( zfvj )
(  224)                   zupsut = zfp_ui * trb(ji,jj,jk,jn) + zfm_ui * trb(ji+1,jj  ,jk,jn)
(  225)                   zupsvt = zfp_vj * trb(ji,jj,jk,jn) + zfm_vj * trb(ji  ,jj+1,jk,jn)
(  226)                   ! centered scheme
(  227)                   zcenut = zfui * ( trn(ji,jj,jk,jn) + trn(ji+1,jj  ,jk,jn) )
(  228)                   zcenvt = zfvj * ( trn(ji,jj,jk,jn) + trn(ji  ,jj+1,jk,jn) )
(  229)                   ! mixed centered / upstream scheme
(  230)                   zwx(ji,jj,jk) = zcofi * zupsut + (1.-zcofi) * zcenut
(  231)                   zwy(ji,jj,jk) = zcofj * zupsvt + (1.-zcofj) * zcenvt                
(  232)                END DO
(  233)             END DO
(  234) 
(  235) 
(  236)             ! 2. Tracer flux divergence at t-point added to the general trend
(  237)             ! -------------------------
(  238) 
(  239)             DO jj = 2, jpjm1
(  240)                DO ji = fs_2, ilast   ! vector opt.
(  241) #if defined key_s_coord || defined key_partial_steps
(  242)                   zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  243) #else
(  244)                   zbtr = zbtr2(ji,jj) 
(  245) #endif
(  246)                   ! horizontal advective trends
(  247)                   ztra = - zbtr * (  zwx(ji,jj,jk) - zwx(ji-1,jj  ,jk)   &
(  248)                      &             + zwy(ji,jj,jk) - zwy(ji  ,jj-1,jk)  )
(  249) 
(  250)                   ! add it to the general tracer trends
(  251)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  252) 
(  253) #if defined key_trc_diatrd 
(  254)                   ! recompute the trends in i- and j-direction as Uh gradh(T)
(  255) # if defined key_s_coord || defined key_partial_steps
(  256)                   zfui = 0.5 * e2u(ji  ,jj) * fse3u(ji,  jj,jk) * zun(ji,  jj,jk)
(  257)                   zfui1= 0.5 * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * zun(ji-1,jj,jk)
(  258)                   zfvj = 0.5 * e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * zvn(ji,jj  ,jk)
(  259)                   zfvj1= 0.5 * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * zvn(ji,jj-1,jk)
(  260) # else
(  261)                   zfui = 0.5 * e2u(ji  ,jj) * zun(ji,  jj,jk)
(  262)                   zfui1= 0.5 * e2u(ji-1,jj) * zun(ji-1,jj,jk)
(  263)                   zfvj = 0.5 * e1v(ji,jj  ) * zvn(ji,jj  ,jk)
(  264)                   zfvj1= 0.5 * e1v(ji,jj-1) * zvn(ji,jj-1,jk)
(  265) # endif
(  266)                   ztai = - zbtr * ( zfui  * ( trn(ji+1,jj  ,jk,jn) - trn(ji,  jj,jk,jn) )   &
(  267)                      &                + zfui1 * ( trn(ji,  jj,  jk,jn) - trn(ji-1,jj,jk,jn) ) )
(  268)                   ztaj = - zbtr * ( zfvj  * ( trn(ji  ,jj+1,jk,jn) - trn(ji,jj  ,jk,jn) )    &
(  269)                      &                + zfvj1 * ( trn(ji  ,jj  ,jk,jn) - trn(ji,jj-1,jk,jn) ) )
(  270)                   ! save i- and j- advective trends computed as Uh gradh(T)
(  271)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = ztai
(  272)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = ztaj
(  273) #endif
(  274)                END DO
(  275)             END DO
(  276)             !                                             ! ===============
(  277)          END DO                                           !   End of slab
(  278)          !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 6

(  279)       ENDDO
(  280) 
(  281)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  282)          WRITE(charout, FMT="('centered2 - had')")
(  283)          CALL prt_ctl_trc_info(charout)
(  284)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  285)       ENDIF
(  286)       
(  287)       ! II. Vertical advection
(  288)       ! ----------------------
(  289)       DO jn = 1, jptra
(  290) 
(  291)          ! Bottom value : flux set to zero
(  292)          zwx(:,:,jpk) = 0.e0 
(  293) 
(  294)          ! Surface value
(  295)          IF ( lk_dynspg_rl ) THEN       ! rigid lid : flux set to zero
(  296)             zwx(:,:, 1 ) = 0.e0  
(  297)          ELSE                           ! free surface-constant volume
(  298)             zwx(:,:, 1 ) = zwn(:,:,1) * trn(:,:,1,jn)
(  299)          ENDIF
(  300) 
(  301)          ! 1. Vertical advective fluxes
(  302)          ! ----------------------------
(  303) 
(  304)          ! Second order centered tracer flux at w-point
(  305) 
(  306)          DO jk = 2, jpk
(  307)             DO jj = 2, jpjm1
(  308)                DO ji = fs_2, ilast   ! vector opt.
(  309)                   ! upstream indicator
(  310)                   zcofk = MAX( zind(ji,jj,jk-1), zind(ji,jj,jk) )
(  311)                   ! velocity * 1/2
(  312)                   zhw = 0.5 * zwn(ji,jj,jk)
(  313)                   ! upstream scheme
(  314)                   zfp_w = zhw + ABS( zhw )
(  315)                   zfm_w = zhw - ABS( zhw )
(  316)                   zupst = zfp_w * trb(ji,jj,jk,jn) + zfm_w * trb(ji,jj,jk-1,jn)
(  317)                   ! centered scheme
(  318)                   zcent = zhw * ( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) )
(  319)                   ! centered scheme
(  320)                   zwx(ji,jj,jk) = zcofk * zupst + (1.-zcofk) * zcent
(  321)                END DO
(  322)             END DO
(  323)          END DO
(  324) 
(  325) 
(  326)          ! 2. Tracer flux divergence at t-point added to the general trend
(  327)          ! -------------------------
(  328) 
(  329)          DO jk = 1, jpkm1
(  330)             DO jj = 2, jpjm1
(  331)                DO ji = fs_2, ilast   ! vector opt.
(  332)                   ze3tr = 1. / fse3t(ji,jj,jk)
(  333)                   ! vertical advective trends
(  334)                   ztra = - ze3tr * ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) )
(  335)                   ! add it to the general tracer trends
(  336)                   tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 7

(  337) #if defined key_trc_diatrd 
(  338)                   ! save the vertical advective trends computed as w gradz(T)
(  339)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = ztra - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
(  340) #endif
(  341)                END DO
(  342)             END DO
(  343)          END DO
(  344) 
(  345)       END DO
(  346) 
(  347)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  348)          WRITE(charout, FMT="('centered - zad')")
(  349)          CALL prt_ctl_trc_info(charout)
(  350)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  351)       ENDIF
(  352) 
(  353)    END SUBROUTINE trc_adv_cen2
(  354) #else
(  355) 
(  356)    !!----------------------------------------------------------------------
(  357)    !!   Default option                                         Empty module
(  358)    !!----------------------------------------------------------------------
(  359) CONTAINS
(  360)    SUBROUTINE trc_adv_cen2( kt )  
(  361)       INTEGER, INTENT(in) :: kt
(  362) !      WRITE(*,*) 'trc_adv_cen2: You should not have seen this print! error?', kt
(  363)    END SUBROUTINE trc_adv_cen2
(  364) #endif
(  365)    !!======================================================================
(  366) END MODULE trcadv_cen2


































PGF90 (Version     10.5)          09/23/2010  09:49:18      page 8

(    1) # 1 "trcadv_cen2.F90"
(    1) MODULE trcadv_cen2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_cen2  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective tracer trend
(    5)    !!==============================================================================
(    6) 
(    7) # 356
(  356)    !!----------------------------------------------------------------------
(  357)    !!   Default option                                         Empty module
(  358)    !!----------------------------------------------------------------------
(  359) CONTAINS
(  360)    SUBROUTINE trc_adv_cen2( kt )  
(  361)       INTEGER, INTENT(in) :: kt
(  362) !      WRITE(*,*) 'trc_adv_cen2: You should not have seen this print! error?', kt
(  363)    END SUBROUTINE trc_adv_cen2
(  364) # 365
(  365)    !!======================================================================
(  366) END MODULE trcadv_cen2
