


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

Filename: trcadv_tvd.F90

(    1) MODULE trcadv_tvd
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_tvd  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_adv_tvd  : update the passive tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a TVD scheme
(   10)    !!   nonosc       : compute monotonic tracer fluxes by a nonoscillatory
(   11)    !!                  algorithm 
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce_trc             ! ocean dynamics and active tracers variables
(   15)    USE trc                 ! ocean passive tracers variables
(   16)    USE lbclnk              ! ocean lateral boundary conditions (or mpp link)
(   17)    USE trcbbl              ! advective passive tracers in the BBL
(   18)    USE prtctl_trc      ! Print control for debbuging
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC trc_adv_tvd    ! routine called by trcstp.F90
(   25) 
(   26)    !! * Substitutions
(   27) #  include "passivetrc_substitute.h90"
(   28)    !!----------------------------------------------------------------------
(   29)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_tvd.F90,v 1.12 2006/04/10 15:38:54 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36)    SUBROUTINE trc_adv_tvd( kt )
(   37)       !!----------------------------------------------------------------------
(   38)       !!                  ***  ROUTINE trc_adv_tvd  ***
(   39)       !! 
(   40)       !! **  Purpose :   Compute the now trend due to total advection of 
(   41)       !!       tracers and add it to the general trend of tracer equations
(   42)       !!
(   43)       !! **  Method  :   TVD scheme, i.e. 2nd order centered scheme with
(   44)       !!       corrected flux (monotonic correction)
(   45)       !!       note: - this advection scheme needs a leap-frog time scheme
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 2

(   47)       !! ** Action : - update tra with the now advective tracer trends
(   48)       !!             - save the trends in trtrd ('key_trc_diatrd)
(   49)       !!
(   50)       !! History :
(   51)       !!        !  95-12  (L. Mortier)  Original code
(   52)       !!        !  00-01  (H. Loukos)  adapted to ORCA 
(   53)       !!        !  00-10  (MA Foujols E.Kestenare)  include file not routine
(   54)       !!        !  00-12  (E. Kestenare M. Levy)  fix bug in trtrd indexes
(   55)       !!        !  01-07  (E. Durand G. Madec)  adaptation to ORCA config
(   56)       !!   9.0  !  02-06  (C. Ethe, G. Madec)  F90: Free form and module
(   57)       !!----------------------------------------------------------------------
(   58)       !! * Modules used
(   59) #if defined key_trcbbl_adv
(   60)       USE oce_trc            , zun => ua,  &  ! use ua as workspace
(   61)          &                     zvn => va      ! use va as workspace
(   62)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   63) #else
(   64)       USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
(   65)                                zvn => vn,  &  !             zvn == vn
(   66)                                zwn => wn      !             zwn == wn
(   67) #endif
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   70) 
(   71)       !! * Local declarations
(   72)       INTEGER  ::   ji, jj, jk,jn           ! dummy loop indices
(   73)       REAL(wp) ::   ztra                    ! temporary scalar
(   74) 
(   75)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   76)          zti, ztu, ztv, ztw                ! temporary workspace
(   77) 
(   78)       REAL(wp) ::   &
(   79)          z2dtt, zbtr, zeu, zev, zew, z2, &  ! temporary scalar
(   80)          zfp_ui, zfp_vj, zfp_wk,         &  !    "         "
(   81)          zfm_ui, zfm_vj, zfm_wk             !    "         "
(   82) 
(   83) #if defined key_trc_diatrd
(   84)        REAL(wp) :: &
(   85)           zgm, zgz
(   86) #endif
(   87) 
(   88)       CHARACTER (len=22) :: charout
(   89)       !!----------------------------------------------------------------------
(   90) 
(   91)       zti(:,:,:) = 0.e0
(   92) 
(   93)       IF( kt == nittrc000  .AND. lwp ) THEN
(   94)          WRITE(numout,*)
(   95)          WRITE(numout,*) 'trc_adv_tvd : TVD advection scheme'
(   96)          WRITE(numout,*) '~~~~~~~~~~~'
(   97)       ENDIF
(   98) 
(   99)       IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(  100)          z2=1.
(  101)       ELSE
(  102)          z2=2.
(  103)       ENDIF
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 3

(  105) #if defined key_trcbbl_adv
(  106)       ! Advective Bottom boundary layer: add the velocity
(  107)       ! -------------------------------------------------
(  108)       zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
(  109)       zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
(  110)       zwn(:,:,:) = wn (:,:,:) + w_trc_bbl(:,:,:)
(  111) #endif
(  112) 
(  113)       DO jn = 1, jptra
(  114) 
(  115)          ! 1. Bottom value : flux set to zero
(  116)          ! ---------------
(  117)          ztu(:,:,jpk) = 0.e0
(  118)          ztv(:,:,jpk) = 0.e0
(  119)          ztw(:,:,jpk) = 0.e0
(  120)          zti(:,:,jpk) = 0.e0
(  121) 
(  122) 
(  123)          ! 2. upstream advection with initial mass fluxes & intermediate update
(  124)          ! --------------------------------------------------------------------
(  125)          ! upstream tracer flux in the i and j direction
(  126)          DO jk = 1, jpkm1
(  127)             DO jj = 1, jpjm1
(  128)                DO ji = 1, fs_jpim1   ! vector opt.
(  129)                   zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  130)                   zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  131)                   ! upstream scheme
(  132)                   zfp_ui = zeu + ABS( zeu )
(  133)                   zfm_ui = zeu - ABS( zeu )
(  134)                   zfp_vj = zev + ABS( zev )
(  135)                   zfm_vj = zev - ABS( zev )
(  136)                   ztu(ji,jj,jk) = zfp_ui * trb(ji,jj,jk,jn) + zfm_ui * trb(ji+1,jj  ,jk,jn)
(  137)                   ztv(ji,jj,jk) = zfp_vj * trb(ji,jj,jk,jn) + zfm_vj * trb(ji  ,jj+1,jk,jn)
(  138)                END DO
(  139)             END DO
(  140)          END DO
(  141) 
(  142)          ! upstream tracer flux in the k direction
(  143)          ! Surface value
(  144)          IF( lk_dynspg_rl ) THEN   ! rigid lid : flux set to zero
(  145)             ztw(:,:,1) = 0.e0
(  146)          ELSE                      ! free surface 
(  147)             DO jj = 1, jpj
(  148)                DO ji = 1, jpi
(  149)                   zew = e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,1)
(  150)                   ztw(ji,jj,1) = zew * trb(ji,jj,1,jn)
(  151)                END DO
(  152)             END DO
(  153)          ENDIF
(  154) 
(  155)          ! Interior value
(  156)          DO jk = 2, jpkm1
(  157)             DO jj = 1, jpj
(  158)                DO ji = 1, jpi
(  159)                   zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  160)                   zfp_wk = zew + ABS( zew )
(  161)                   zfm_wk = zew - ABS( zew )
(  162)                   ztw(ji,jj,jk) = zfp_wk * trb(ji,jj,jk,jn) + zfm_wk * trb(ji,jj,jk-1,jn)






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 4

(  163)                END DO
(  164)             END DO
(  165)          END DO
(  166) 
(  167)          ! total advective trend
(  168)          DO jk = 1, jpkm1
(  169)             DO jj = 2, jpjm1
(  170)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  171)                   zbtr = 1./ ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  172)                   zti(ji,jj,jk) = - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )   &
(  173)                      &              + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )   &
(  174)                      &              + ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
(  175) 
(  176) #if defined key_trc_diatrd
(  177)                   IF ( luttrd(jn) ) &
(  178)                      trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
(  179)                         &                          zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )                     
(  180)                   IF ( luttrd(jn) ) &
(  181)                      trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
(  182)                         &                          zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )                     
(  183)                   IF ( luttrd(jn) ) &
(  184)                      trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
(  185)                         &                          zbtr * ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) )
(  186) #endif
(  187)                END DO
(  188)             END DO
(  189)          END DO
(  190) 
(  191) 
(  192)          ! update and guess with monotonic sheme
(  193)          DO jk = 1, jpkm1
(  194)             z2dtt = z2 * rdttra(jk) * FLOAT(ndttrc)
(  195)             DO jj = 2, jpjm1
(  196)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  197)                   tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + zti(ji,jj,jk)
(  198)                   zti (ji,jj,jk) = ( trb(ji,jj,jk,jn) + z2dtt * zti(ji,jj,jk) ) * tmask(ji,jj,jk)
(  199)                END DO
(  200)             END DO
(  201)          END DO
(  202) 
(  203)          ! Lateral boundary conditions on zti, zsi   (unchanged sign)
(  204)          CALL lbc_lnk( zti, 'T', 1. )
(  205) 
(  206)          ! 3. antidiffusive flux : high order minus low order
(  207)          ! --------------------------------------------------
(  208)          ! antidiffusive flux on i and j
(  209)          DO jk = 1, jpkm1
(  210)             DO jj = 1, jpjm1
(  211)                DO ji = 1, fs_jpim1   ! vector opt.
(  212)                   zeu = 0.5 * e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  213)                   zev = 0.5 * e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  214)                   ztu(ji,jj,jk) = zeu * ( trn(ji,jj,jk,jn) + trn(ji+1,jj,jk,jn) ) - ztu(ji,jj,jk)
(  215)                   ztv(ji,jj,jk) = zev * ( trn(ji,jj,jk,jn) + trn(ji,jj+1,jk,jn) ) - ztv(ji,jj,jk)
(  216)                END DO
(  217)             END DO
(  218)          END DO
(  219) 
(  220)          ! antidiffusive flux on k






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 5

(  221)          ! Surface value
(  222)          ztw(:,:,1) = 0.
(  223) 
(  224)          ! Interior value
(  225)          DO jk = 2, jpkm1
(  226)             DO jj = 1, jpj
(  227)                DO ji = 1, jpi
(  228)                   zew = 0.5 * e1t(ji,jj) * e2t(ji,jj) * zwn(ji,jj,jk)
(  229)                   ztw(ji,jj,jk) = zew * ( trn(ji,jj,jk,jn) + trn(ji,jj,jk-1,jn) ) - ztw(ji,jj,jk)
(  230)                END DO
(  231)             END DO
(  232)          END DO
(  233) 
(  234)          ! Lateral bondary conditions
(  235)          CALL lbc_lnk( ztu, 'U', -1. )
(  236)          CALL lbc_lnk( ztv, 'V', -1. )
(  237)          CALL lbc_lnk( ztw, 'W',  1. )
(  238) 
(  239)          ! 4. monotonicity algorithm
(  240)          ! -------------------------
(  241)          CALL nonosc( trb(:,:,:,jn), ztu, ztv, ztw, zti, z2 )
(  242) 
(  243) 
(  244)          ! 5. final trend with corrected fluxes
(  245)          ! ------------------------------------
(  246)          DO jk = 1, jpkm1
(  247)             DO jj = 2, jpjm1
(  248)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  249)                   zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  250) #if defined key_trc_diatrd
(  251)                   IF ( luttrd(jn) ) &
(  252)                      trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
(  253)                         &                          zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )                     
(  254)                   IF ( luttrd(jn) ) &
(  255)                      trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
(  256)                         &                          zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )                     
(  257)                   IF ( luttrd(jn) ) &
(  258)                      trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
(  259)                         &                          zbtr * ( ztw(ji,jj,jk) - ztw(ji,jj,jk+1) )
(  260) #endif
(  261)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn)   &
(  262)                      &         - ( ztu(ji,jj,jk) - ztu(ji-1,jj  ,jk  )   &
(  263)                      &           + ztv(ji,jj,jk) - ztv(ji  ,jj-1,jk  )   &
(  264)                      &           + ztw(ji,jj,jk) - ztw(ji  ,jj  ,jk+1) ) * zbtr
(  265)                END DO
(  266)             END DO
(  267)          END DO
(  268)          ! 6.0 convert the transport trend into advection trend
(  269)          ! ----------------------------------------------------
(  270)          
(  271) #if defined key_trc_diatrd
(  272)          DO jk = 1,jpk
(  273)             DO jj = 2,jpjm1
(  274)                DO  ji = 2,jpim1
(  275)                   zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
(  276)                   zgm = zbtr * trn(ji,jj,jk,jn) *     &
(  277)                      &         (  zun(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)    &
(  278)                      &          - zun(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 6

(  279)                   
(  280)                   zgz = zbtr * trn(ji,jj,jk,jn) *     &
(  281)                      &         (  zvn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)    &
(  282)                      &          - zvn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk) )
(  283)                   
(  284)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) + zgm
(  285)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) + zgz
(  286)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3)    &
(  287)                      &            - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
(  288)                END DO
(  289)             END DO
(  290)          END DO
(  291)          
(  292)          ! Lateral boundary conditions on trtrd:
(  293)          
(  294)          IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),1), 'T', 1. )
(  295)          IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),2), 'T', 1. )
(  296)          IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),3), 'T', 1. )
(  297) #endif
(  298) 
(  299)       END DO
(  300) 
(  301)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  302)          WRITE(charout, FMT="('tvd - adv')")
(  303)          CALL prt_ctl_trc_info(charout)
(  304)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  305)       ENDIF
(  306) 
(  307)    END SUBROUTINE trc_adv_tvd
(  308) 
(  309) 
(  310)    SUBROUTINE nonosc( pbef, paa, pbb, pcc, paft, prdt )
(  311)       !!---------------------------------------------------------------------
(  312)       !!                    ***  ROUTINE nonosc  ***
(  313)       !!     
(  314)       !! **  Purpose :   compute monotonic tracer fluxes from the upstream 
(  315)       !!       scheme and the before field by a nonoscillatory algorithm 
(  316)       !!
(  317)       !! **  Method  :   ... ???
(  318)       !!       warning : pbef and paft must be masked, but the boundaries
(  319)       !!       conditions on the fluxes are not necessary zalezak (1979)
(  320)       !!       drange (1995) multi-dimensional forward-in-time and upstream-
(  321)       !!       in-space based differencing for fluid
(  322)       !!
(  323)       !! History :
(  324)       !!        !  97-04  (L. Mortier) Original code
(  325)       !!        !  00-02  (H. Loukos)  rewritting for opa8
(  326)       !!        !  00-10  (M.A Foujols, E. Kestenare)  lateral b.c.
(  327)       !!        !  01-03  (E. Kestenare)  add key_passivetrc
(  328)       !!        !  01-07  (E. Durand G. Madec)  adapted for T & S
(  329)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  330)       !!----------------------------------------------------------------------
(  331)       !! * Arguments
(  332)       REAL(wp), INTENT( in ) ::   &
(  333)          prdt                               ! ???
(  334)       REAL(wp), DIMENSION (jpi,jpj,jpk), INTENT( inout ) ::   &
(  335)          pbef,                            & ! before field
(  336)          paft,                            & ! after field






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 7

(  337)          paa,                             & ! monotonic flux in the i direction
(  338)          pbb,                             & ! monotonic flux in the j direction
(  339)          pcc                                ! monotonic flux in the k direction
(  340) 
(  341)       !! * Local declarations
(  342)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(  343)       INTEGER ::   ikm1
(  344)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   zbetup, zbetdo
(  345)       REAL(wp) ::   zpos, zneg, zbt, za, zb, zc, zbig, zrtrn, z2dtt
(  346)       !!----------------------------------------------------------------------
(  347) 
(  348)       zbig = 1.e+40
(  349)       zrtrn = 1.e-15
(  350)       zbetup(:,:,:) = 0.e0   ;   zbetdo(:,:,:) = 0.e0
(  351) 
(  352)       ! Search local extrema
(  353)       ! --------------------
(  354)       ! large negative value (-zbig) inside land
(  355)       ! large negative value (-zbig) inside land
(  356)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  357)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) - zbig * ( 1.e0 - tmask(:,:,:) )
(  358)       ! search maximum in neighbourhood
(  359)       DO jk = 1, jpkm1
(  360)          ikm1 = MAX(jk-1,1)
(  361)          DO jj = 2, jpjm1
(  362)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  363)                zbetup(ji,jj,jk) = MAX(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  364)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  365)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
(  366)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  367)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  368)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  369)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  370)             END DO
(  371)          END DO
(  372)       END DO
(  373)       ! large positive value (+zbig) inside land
(  374)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  375)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:) + zbig * ( 1.e0 - tmask(:,:,:) )
(  376)       ! search minimum in neighbourhood
(  377)       DO jk = 1, jpkm1
(  378)          ikm1 = MAX(jk-1,1)
(  379)          DO jj = 2, jpjm1
(  380)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  381)                zbetdo(ji,jj,jk) = MIN(  pbef(ji  ,jj  ,jk  ), paft(ji  ,jj  ,jk  ),   &
(  382)                   &                     pbef(ji-1,jj  ,jk  ), pbef(ji+1,jj  ,jk  ),   &
(  383)                   &                     paft(ji-1,jj  ,jk  ), paft(ji+1,jj  ,jk  ),   &
(  384)                   &                     pbef(ji  ,jj-1,jk  ), pbef(ji  ,jj+1,jk  ),   &
(  385)                   &                     paft(ji  ,jj-1,jk  ), paft(ji  ,jj+1,jk  ),   &
(  386)                   &                     pbef(ji  ,jj  ,ikm1), pbef(ji  ,jj  ,jk+1),   &
(  387)                   &                     paft(ji  ,jj  ,ikm1), paft(ji  ,jj  ,jk+1)  )
(  388)             END DO
(  389)          END DO
(  390)       END DO
(  391) 
(  392)       ! restore masked values to zero
(  393)       pbef(:,:,:) = pbef(:,:,:) * tmask(:,:,:)
(  394)       paft(:,:,:) = paft(:,:,:) * tmask(:,:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 8

(  395)  
(  396) 
(  397)       ! 2. Positive and negative part of fluxes and beta terms
(  398)       ! ------------------------------------------------------
(  399) 
(  400)       DO jk = 1, jpkm1
(  401)          z2dtt = prdt * rdttra(jk) * FLOAT(ndttrc)
(  402)          DO jj = 2, jpjm1
(  403)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  404)                ! positive & negative part of the flux
(  405)                zpos = MAX( 0., paa(ji-1,jj  ,jk  ) ) - MIN( 0., paa(ji  ,jj  ,jk  ) )   &
(  406)                   & + MAX( 0., pbb(ji  ,jj-1,jk  ) ) - MIN( 0., pbb(ji  ,jj  ,jk  ) )   &
(  407)                   & + MAX( 0., pcc(ji  ,jj  ,jk+1) ) - MIN( 0., pcc(ji  ,jj  ,jk  ) )
(  408)                zneg = MAX( 0., paa(ji  ,jj  ,jk  ) ) - MIN( 0., paa(ji-1,jj  ,jk  ) )   &
(  409)                   & + MAX( 0., pbb(ji  ,jj  ,jk  ) ) - MIN( 0., pbb(ji  ,jj-1,jk  ) )   &
(  410)                   & + MAX( 0., pcc(ji  ,jj  ,jk  ) ) - MIN( 0., pcc(ji  ,jj  ,jk+1) )
(  411)                ! up & down beta terms
(  412)                zbt = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) / z2dtt
(  413)                zbetup(ji,jj,jk) = ( zbetup(ji,jj,jk) - paft(ji,jj,jk) ) / (zpos+zrtrn) * zbt
(  414)                zbetdo(ji,jj,jk) = ( paft(ji,jj,jk) - zbetdo(ji,jj,jk) ) / (zneg+zrtrn) * zbt
(  415)             END DO
(  416)          END DO
(  417)       END DO
(  418) 
(  419)       ! lateral boundary condition on zbetup & zbetdo   (unchanged sign)
(  420)       CALL lbc_lnk( zbetup, 'T', 1. )
(  421)       CALL lbc_lnk( zbetdo, 'T', 1. )
(  422) 
(  423) 
(  424)       ! 3. monotonic flux in the i direction, i.e. paa
(  425)       ! ----------------------------------------------
(  426)       DO jk = 1, jpkm1
(  427)          DO jj = 2, jpjm1
(  428)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  429)                zc = paa(ji,jj,jk)
(  430)                IF( zc >= 0. ) THEN
(  431)                   za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji+1,jj,jk) )
(  432)                   paa(ji,jj,jk) = za * zc
(  433)                ELSE
(  434)                   zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji+1,jj,jk) )
(  435)                   paa(ji,jj,jk) = zb * zc
(  436)                ENDIF
(  437)             END DO
(  438)          END DO
(  439)       END DO
(  440) 
(  441)       ! lateral boundary condition on paa   (changed sign)
(  442)       CALL lbc_lnk( paa, 'U', -1. )
(  443) 
(  444) 
(  445)       ! 4. monotonic flux in the j direction, i.e. pbb
(  446)       ! ----------------------------------------------
(  447)       DO jk = 1, jpkm1
(  448)          DO jj = 2, jpjm1
(  449)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  450)                zc = pbb(ji,jj,jk)
(  451)                IF( zc >= 0. ) THEN
(  452)                   za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj+1,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 9

(  453)                   pbb(ji,jj,jk) = za * zc
(  454)                ELSE
(  455)                   zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj+1,jk) )
(  456)                   pbb(ji,jj,jk) = zb * zc
(  457)                ENDIF
(  458)             END DO
(  459)          END DO
(  460)       END DO
(  461) 
(  462)       ! lateral boundary condition on pbb   (changed sign)
(  463)       CALL lbc_lnk( pbb, 'V', -1. )
(  464) 
(  465) 
(  466)       ! monotonic flux in the k direction, i.e. pcc
(  467)       ! -------------------------------------------
(  468)       DO jk = 2, jpkm1
(  469)          DO jj = 2, jpjm1
(  470)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  471)                zc = pcc(ji,jj,jk)
(  472)                IF( zc >= 0. ) THEN
(  473)                   za = MIN( 1., zbetdo(ji,jj,jk), zbetup(ji,jj,jk-1) )
(  474)                   pcc(ji,jj,jk) = za * zc
(  475)                ELSE
(  476)                   zb = MIN( 1., zbetup(ji,jj,jk), zbetdo(ji,jj,jk-1) )
(  477)                   pcc(ji,jj,jk) = zb * zc
(  478)                ENDIF
(  479)             END DO
(  480)          END DO
(  481)       END DO
(  482) 
(  483)       ! lateral boundary condition on pcc   (unchanged sign)
(  484)       CALL lbc_lnk( pcc, 'W', 1. )
(  485) 
(  486)    END SUBROUTINE nonosc
(  487) 
(  488) #else
(  489)    !!----------------------------------------------------------------------
(  490)    !!   Default option                                         Empty module
(  491)    !!----------------------------------------------------------------------
(  492) CONTAINS
(  493)    SUBROUTINE trc_adv_tvd( kt )  
(  494)       INTEGER, INTENT(in) :: kt
(  495) !      WRITE(*,*) 'trc_adv_tvd: You should not have seen this print! error?', kt
(  496)    END SUBROUTINE trc_adv_tvd
(  497) #endif
(  498) 
(  499)    !!======================================================================
(  500) END MODULE trcadv_tvd
















PGF90 (Version     10.5)          09/23/2010  09:49:19      page 10

(    1) # 1 "trcadv_tvd.F90"
(    1) MODULE trcadv_tvd
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_tvd  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) # 489
(  489)    !!----------------------------------------------------------------------
(  490)    !!   Default option                                         Empty module
(  491)    !!----------------------------------------------------------------------
(  492) CONTAINS
(  493)    SUBROUTINE trc_adv_tvd( kt )  
(  494)       INTEGER, INTENT(in) :: kt
(  495) !      WRITE(*,*) 'trc_adv_tvd: You should not have seen this print! error?', kt
(  496)    END SUBROUTINE trc_adv_tvd
(  497) 
(  498) # 499
(  499)    !!======================================================================
(  500) END MODULE trcadv_tvd
