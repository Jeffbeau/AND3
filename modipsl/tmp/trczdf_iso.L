


PGF90 (Version     10.5)          09/23/2010  09:49:23      page 1

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

Filename: trczdf_iso.F90

(    1) MODULE trczdf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_iso  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc && ( defined key_ldfslp   ||   defined key_esopa )
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_zdf_iso  : update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce_trc          ! ocean dynamics and tracers variables
(   16)    USE trc              ! ocean passive tracers variables 
(   17)    USE lbclnk           ! ocean lateral boundary conditions (or mpp link)
(   18)    USE trctrp_lec       ! passive tracers transport
(   19)    USE prtctl_trc          ! Print control for debbuging
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC trc_zdf_iso    ! routine called by step.F90
(   26) 
(   27)    !! * Module variable
(   28)    REAL(wp), DIMENSION(jpk) ::   &
(   29)       rdttrc                     ! vertical profile of 2 x tracer time-step
(   30) 
(   31)    !! * Substitutions
(   32) #  include "passivetrc_substitute.h90"
(   33)    !!----------------------------------------------------------------------
(   34)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_iso.F90,v 1.10 2006/04/10 15:38:55 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE trc_zdf_iso( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE trc_zdf_iso  ***
(   44)       !!
(   45)       !! ** Purpose :
(   46)       !!     Compute the trend due to the vertical tracer diffusion inclu-






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47)       !!     ding the vertical component of lateral mixing (only for second
(   48)       !!     order operator, for fourth order it is already computed and
(   49)       !!     add to the general trend in trcldf.F) and add it to the general
(   50)       !!     trend of the tracer equations.
(   51)       !!
(   52)       !! ** Method :
(   53)       !!         The vertical component of the lateral diffusive trends is
(   54)       !!      provided by a 2nd order operator rotated along neural or geopo-
(   55)       !!      tential surfaces to which an eddy induced advection can be added
(   56)       !!      It is computed using before fields (forward in time) and isopyc-
(   57)       !!      nal or geopotential slopes computed in routine ldfslp.
(   58)       !!
(   59)       !!      First part: vertical trends associated with the lateral mixing
(   60)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(   61)       !!      vertical fluxes associated with the rotated lateral mixing:
(   62)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(trb)) ]
(   63)       !!                     + e1t*wslpj dj[ mj(mk(trb)) ]  }
(   64)       !!      save avt coef. resulting from vertical physics alone in zavt:
(   65)       !!         zavt = avt
(   66)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(   67)       !!         avt = avt + wslpi^2+wslj^2
(   68)       !!      add vertical Eddy Induced advective fluxes ('lk_trcldf_eiv=T):
(   69)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(   70)       !!                    +dj[aht e1v mj(wslpj)] } mk(trb)
(   71)       !!      take the horizontal divergence of the fluxes:
(   72)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(   73)       !!      Add this trend to the general trend tra :
(   74)       !!         tra = tra + difft
(   75)       !!
(   76)       !!      Second part: vertical trend associated with the vertical physics
(   77)       !!      ===========  (including the vertical flux proportional to dk[t]
(   78)       !!                  associated with the lateral mixing, through the
(   79)       !!                  update of avt)
(   80)       !!      The vertical diffusion of tracers tra  is given by:
(   81)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(   82)       !!      It is computed using a backward time scheme, t=ta.
(   83)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   84)       !!      both tracers (bottom, applied through the masked field avt).
(   85)       !!      Add this trend to the general trend tra  :
(   86)       !!         tra = tra + dz( avt dz(t) )
(   87)       !!         (tra = tra + dz( avs dz(t) ) if lk_trc_zdfddm=T )
(   88)       !!
(   89)       !!      Third part: recover avt resulting from the vertical physics
(   90)       !!      ==========  alone, for further diagnostics (for example to
(   91)       !!                  compute the turbocline depth in diamld).
(   92)       !!         avt = zavt
(   93)       !!         (avs = zavs if lk_trc_zdfddm=T )
(   94)       !!
(   95)       !!      'key_trc_diatrd' defined: trend saved for futher diagnostics.
(   96)       !!
(   97)       !!      macro-tasked on vertical slab (jj-loop)
(   98)       !!
(   99)       !! ** Action :
(  100)       !!         Update tra arrays with the before vertical diffusion trend
(  101)       !!         Save in trtrd arrays the trends if 'key_trc_diatrd' defined
(  102)       !!
(  103)       !! History :
(  104)       !!   7.0  !  91-11  (G. Madec)  Original code






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(  106)       !!        !  96-01  (G. Madec)  statement function for e3
(  107)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(  108)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(  109)       !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
(  110)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(  111)       !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
(  112)       !!        !                     avt multiple correction
(  113)       !!        !  00-08  (G. Madec)  double diffusive mixing
(  114)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  115)       !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
(  116)       !!---------------------------------------------------------------------
(  117)       !! * Modules used
(  118)       USE oce_trc               ,   &
(  119)          zavs => va
(  120) 
(  121)       !! * Arguments
(  122)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(  123) 
(  124)       !! * Local declarations
(  125)       INTEGER ::   ji, jj, jk, jn                 ! dummy loop indices
(  126)       INTEGER ::   ikst, ikenm2, ikstp1       ! temporary integers
(  127) #if defined key_partial_steps
(  128)       INTEGER ::   iku, ikv, ikv1             ! temporary integers
(  129) #endif
(  130)       REAL(wp) ::   ztra
(  131)       REAL(wp) ::   &
(  132)          ztavg,                 &  ! ???
(  133)          zcoef0, zcoef3,        &  ! ???
(  134)          zcoef4, zavi,          &  ! ???
(  135)          zbtr, zmku, zmkv,      &  !
(  136)          ztav
(  137)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  138)          zwd, zws, zwi,         &  ! ???
(  139)          zwx, zwy, zwz, zwt        ! ???
(  140)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  141)          ztfw, zdit, zdjt, zdj1t
(  142) #if defined key_trcldf_eiv   ||   defined key_esopa
(  143)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(  144)          ztfwg
(  145) 
(  146)       REAL(wp) ::         &
(  147)          zcoeg3,          &
(  148)          zuwk, zvwk,      &
(  149)          zuwki, zvwki
(  150) #endif
(  151)       CHARACTER (len=22) :: charout
(  152)       !!---------------------------------------------------------------------
(  153) 
(  154)       IF( kt == nittrc000 ) THEN
(  155)          IF(lwp) WRITE(numout,*)
(  156)          IF(lwp) WRITE(numout,*) 'trc_zdf_iso : vertical mixing (including isopycnal component)'
(  157)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
(  158) #if defined key_trcldf_eiv && defined key_diaeiv
(  159)          w_trc_eiv(:,:,:) = 0.e0
(  160) #endif
(  161)       ENDIF
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163)       ! 0.0  Local constant initialization
(  164)       ! --------------------------------
(  165)       IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
(  166)          ! time step = 2 rdttra with Arakawa or TVD advection scheme
(  167)          IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(  168)             rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
(  169)          ELSEIF( kt <= nittrc000 + 1 ) THEN
(  170)             rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
(  171)          ENDIF
(  172)       ELSE
(  173)          rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
(  174)       ENDIF
(  175) 
(  176) 
(  177) 
(  178)       ! 0.1 Save avs in zavs to recover avs in output files
(  179)       !---------------------------------------------------
(  180)       zavs(:,:,:) = fstravs(:,:,:)
(  181) 
(  182) 
(  183) 
(  184)       DO jn = 1, jptra
(  185) 
(  186)          ztavg = 0.e0
(  187) 
(  188)          !                                                ! ===============
(  189)          DO jj = 2, jpjm1                                 !  Vertical slab
(  190)             !                                             ! ===============
(  191) 
(  192)             ! I. vertical trends associated with the lateral mixing
(  193)             ! =====================================================
(  194)             !  (excluding the vertical flux proportional to dk[t]
(  195) 
(  196) 
(  197)             ! I.1 horizontal tracer gradient
(  198)             ! ------------------------------
(  199) 
(  200)             DO jk = 1, jpkm1
(  201)                DO ji = 1, jpim1
(  202)                   ! i-gradient of passive tracer at jj
(  203)                   zdit (ji,jk) = ( trb(ji+1,jj,jk,jn)-trb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
(  204)                   ! j-gradient of passive tracer at jj
(  205)                   zdjt (ji,jk) = ( trb(ji,jj+1,jk,jn)-trb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
(  206)                   ! j-gradient of passive tracer at jj+1
(  207)                   zdj1t(ji,jk) = ( trb(ji,jj,jk,jn)-trb(ji,jj-1,jk,jn) ) * vmask(ji,jj-1,jk)
(  208)                END DO
(  209)             END DO
(  210) #  if defined key_partial_steps
(  211)             ! partial steps correction at the bottom ocean level 
(  212)             DO ji = 1, jpim1
(  213)                ! last ocean level
(  214)                iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  215)                ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  216)                ikv1 = MIN( mbathy(ji,jj), mbathy(ji  ,jj-1) ) - 1
(  217)                ! i-gradient of of passive tracer at jj
(  218)                zdit (ji,iku) = gtru(ji,jj,jn)
(  219)                ! j-gradient of of passive tracer at jj
(  220)                zdjt (ji,ikv) = gtrv(ji,jj,jn) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(  221)                ! j-gradient of of passive tracer at jj+1
(  222)                zdj1t(ji,ikv1)= gtrv(ji,jj-1,jn)
(  223)             END DO
(  224) #endif
(  225) 
(  226) 
(  227)             ! I.2 Vertical fluxes
(  228)             ! -------------------
(  229) 
(  230)             ! Surface and bottom vertical fluxes set to zero
(  231)             ztfw(:, 1 ) = 0.e0
(  232)             ztfw(:,jpk) = 0.e0
(  233) 
(  234) #if defined key_trcldf_eiv
(  235)             ztfwg(:, 1 ) = 0.e0
(  236)             ztfwg(:,jpk) = 0.e0
(  237) #endif
(  238) 
(  239)             ! interior (2=<jk=<jpk-1)
(  240)             DO jk = 2, jpkm1
(  241)                DO ji = 2, jpim1
(  242)                   zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  243) 
(  244)                   zmku = 1./MAX( umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)   &
(  245)                      &          +umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1. )
(  246) 
(  247)                   zmkv = 1./MAX( vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)   &
(  248)                      &          +vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1. )
(  249) 
(  250)                   zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  251)                   zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  252) 
(  253)                   ztfw(ji,jk) = zcoef3 * ( zdit (ji  ,jk-1) + zdit (ji-1,jk)     &
(  254)                      &                    +zdit (ji-1,jk-1) + zdit (ji  ,jk) )   &
(  255)                      &        + zcoef4 * ( zdjt (ji  ,jk-1) + zdj1t(ji  ,jk)     &
(  256)                      &                    +zdj1t(ji  ,jk-1) + zdjt (ji  ,jk) )
(  257) 
(  258)                END DO
(  259)             END DO
(  260) 
(  261) 
(  262)             ! I.3  update and save of avt (and avs if double diffusive mixing)
(  263)             ! ---------------------------
(  264) 
(  265)             DO jk = 2, jpkm1
(  266)                DO ji = 2, jpim1
(  267) 
(  268)                   zavi = fsahtw(ji,jj,jk)*( wslpi(ji,jj,jk)*wslpi(ji,jj,jk)   &
(  269)                      &                     +wslpj(ji,jj,jk)*wslpj(ji,jj,jk) )
(  270) 
(  271)                   ! add isopycnal vertical coeff. to avs
(  272)                   fstravs(ji,jj,jk) = fstravs(ji,jj,jk) + zavi
(  273) 
(  274)                END DO
(  275)             END DO
(  276) 
(  277) #if defined key_trcldf_eiv
(  278)             !                              ! ---------------------------------------!






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 6

(  279)             !                              ! Eddy induced vertical advective fluxes !
(  280)             !                              ! ---------------------------------------!
(  281) #if defined key_traldf_c2d || defined key_traldf_c3d
(  282)             DO jk = 2, jpkm1
(  283)                DO ji = 2, jpim1
(  284)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  285)                      &  * fsaeitru(ji-1,jj,jk) * e2u(ji-1,jj)*umask(ji-1,jj,jk)
(  286)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  287)                      &  * fsaeitru(ji  ,jj,jk) * e2u(ji  ,jj)*umask(ji  ,jj,jk)
(  288)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  289)                      &  * fsaeitrv(ji,jj-1,jk) * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
(  290)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  291)                      &  * fsaeitrv(ji,jj  ,jk) * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
(  292) 
(  293)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
(  294) 
(  295)                   ztfwg(ji,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) ) 
(  296)                   ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jk)
(  297) 
(  298) # if defined key_diaeiv
(  299)                   w_trc_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  300) # endif
(  301)                END DO
(  302)             END DO
(  303) 
(  304) #else
(  305)             DO jk = 2, jpkm1
(  306)                DO ji = 2, jpim1
(  307)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  308)                      &  * e2u(ji-1,jj)*umask(ji-1,jj,jk)
(  309)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  310)                      &  * e2u(ji  ,jj)*umask(ji  ,jj,jk)
(  311)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  312)                      &  * e1v(ji,jj-1)*vmask(ji,jj-1,jk)
(  313)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  314)                      &  * e1v(ji  ,jj)*vmask(ji  ,jj,jk)
(  315) 
(  316)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeitrw(ji,jj,jk)   &
(  317)                      &            * ( zuwk - zuwki + zvwk - zvwki )
(  318) 
(  319)                   ztfwg(ji,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) )
(  320)                   ztfw(ji,jk) = ztfw(ji,jk) + ztfwg(ji,jk)
(  321) 
(  322) # if defined key_diaeiv
(  323)                   w_trc_eiv(ji,jj,jk) = -2. *  zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  324) # endif
(  325)                END DO
(  326)             END DO
(  327) #endif
(  328) 
(  329) #endif
(  330) 
(  331) 
(  332)             ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  333)             ! -------------------------------------------------------------------
(  334) 
(  335)             DO jk = 1, jpkm1
(  336)                DO ji = 2, jpim1






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 7

(  337)                   zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  338)                   ztav = (  ztfw(ji,jk) - ztfw(ji,jk+1)  ) * zbtr
(  339)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztav
(  340) 
(  341) #if defined key_trc_diatrd
(  342) #   if defined key_trcldf_eiv
(  343)                   ztavg = ( ztfwg(ji,jk) - ztfwg(ji,jk+1) ) * zbtr
(  344)                   !  WARNING trtrd(ji,jj,jk,6) used for vertical gent velocity trend
(  345)                   !                           not for damping !!!
(  346)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztavg
(  347) #   endif
(  348)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztav - ztavg
(  349) #endif
(  350)                END DO
(  351)             END DO
(  352)             !                                             ! ===============
(  353)          END DO                                           !   End of slab
(  354)          !                                                ! ===============
(  355) 
(  356)       END DO
(  357) 
(  358)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  359)          WRITE(charout, FMT="('zdf - 1')")
(  360)          CALL prt_ctl_trc_info(charout)
(  361)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  362)       ENDIF
(  363) 
(  364)       DO jn = 1, jptra
(  365)          !                                                ! ===============
(  366)          DO jj = 2, jpjm1                                 !  Vertical slab
(  367)             !                                             ! ===============
(  368) 
(  369)             ! II. Vertical trend associated with the vertical physics
(  370)             ! =======================================================
(  371)             !     (including the vertical flux proportional to dk[t] associated
(  372)             !      with the lateral mixing, through the avt update)
(  373)             !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  374) 
(  375) 
(  376)             ! Diagonal, inferior, superior
(  377)             ! (including the bottom boundary condition via avs masked)
(  378)             DO jk = 1, jpkm1
(  379)                DO ji = 2, jpim1
(  380)                   zwi(ji,jk) = - rdttrc(jk) * fstravs(ji,jj,jk  )   &
(  381)                      /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  382)                   zws(ji,jk) = - rdttrc(jk) * fstravs(ji,jj,jk+1)   &
(  383)                      /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  384)                   zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  385)                END DO
(  386)             END DO
(  387) 
(  388)             ! Surface boudary conditions
(  389)             DO ji = 2, jpim1
(  390)                zwi(ji,1) = 0.e0
(  391)                zwd(ji,1) = 1. - zws(ji,1)
(  392)             END DO
(  393) 
(  394)             ! Second member construction






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 8

(  395)             DO jk = 1, jpkm1
(  396)                DO ji = 2, jpim1
(  397)                   zwy(ji,jk) = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)
(  398)                END DO
(  399)             END DO
(  400) 
(  401)             ! Matrix inversion from the first level
(  402)             ikst = 1
(  403) #   include "zdf.matrixsolver.h90"
(  404) #if defined key_trc_diatrd
(  405)             ! Compute and save the vertical diffusive of tracers trends
(  406) #  if defined key_trc_ldfiso
(  407)             DO jk = 1, jpkm1
(  408)                DO ji = 2, jpim1
(  409)                   ztra = ( zwx(ji,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  410)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn) + trtrd(ji,jj,jk,ikeep(jn),6)
(  411)                END DO
(  412)             END DO
(  413) #  else
(  414)             DO jk = 1, jpkm1
(  415)                DO ji = 2, jpim1
(  416)                   ztra = ( zwx(ji,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  417)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn)
(  418)                END DO
(  419)             END DO
(  420) #  endif
(  421) #endif
(  422)             ! Save the masked passive tracer after in tra
(  423)             ! (c a u t i o n:  tracer not its trend, Leap-frog scheme done
(  424)             !                  it will not be done in tranxt)
(  425)             DO jk = 1, jpkm1
(  426)                DO ji = 2, jpim1
(  427)                   tra(ji,jj,jk,jn) = zwx(ji,jk)  * tmask(ji,jj,jk)
(  428)                END DO
(  429)             END DO
(  430)             !                                             ! ===============
(  431)          END DO                                           !   End of slab
(  432)          !                                                ! ===============
(  433) 
(  434)       END DO
(  435) 
(  436) 
(  437) 
(  438)       ! III. recover the avt (avs) resulting from vertical physics only
(  439)       !---------------------------------------------------------------
(  440)       fstravs(:,:,:) = zavs(:,:,:)
(  441) 
(  442)       
(  443)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  444)          WRITE(charout, FMT="('zdf - 2')")
(  445)          CALL prt_ctl_trc_info(charout)
(  446)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  447)       ENDIF
(  448) 
(  449)    END SUBROUTINE trc_zdf_iso
(  450) 
(  451) #else
(  452)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 9

(  453)    !!   Dummy module               NO rotation of the lateral mixing tensor
(  454)    !!----------------------------------------------------------------------
(  455) CONTAINS
(  456)    SUBROUTINE trc_zdf_iso( kt )              ! empty routine
(  457) !      WRITE(*,*) 'trc_zdf_iso: You should not have seen this print! error?', kt
(  458)    END SUBROUTINE trc_zdf_iso
(  459) #endif
(  460) 
(  461)    !!==============================================================================
(  462) END MODULE trczdf_iso






















































PGF90 (Version     10.5)          09/23/2010  09:49:23      page 10

(    1) # 1 "trczdf_iso.F90"
(    1) MODULE trczdf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_iso  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) # 452
(  452)    !!----------------------------------------------------------------------
(  453)    !!   Dummy module               NO rotation of the lateral mixing tensor
(  454)    !!----------------------------------------------------------------------
(  455) CONTAINS
(  456)    SUBROUTINE trc_zdf_iso( kt )              ! empty routine
(  457) !      WRITE(*,*) 'trc_zdf_iso: You should not have seen this print! error?', kt
(  458)    END SUBROUTINE trc_zdf_iso
(  459) 
(  460) # 461
(  461)    !!==============================================================================
(  462) END MODULE trczdf_iso
