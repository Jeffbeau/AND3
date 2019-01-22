


PGF90 (Version     10.5)          09/23/2010  09:49:36      page 1

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

Filename: dynhpg.F90

(    1) MODULE dynhpg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynhpg  ***
(    4)    !! Ocean dynamics:  hydrostatic pressure gradient trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_hpg      : update the momentum trend with the horizontal
(    9)    !!                  gradient of the hydrostatic pressure
(   10)    !!
(   11)    !!   default case : use of 3D work arrays (vector opt. available)
(   12)    !!   key_s_coord       : s-coordinate
(   13)    !!   key_partial_steps : z-coordinate with partial steps
(   14)    !!   default key       : z-coordinate
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE phycst          ! physical constants
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE trdmod          ! ocean dynamics trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC dyn_hpg                ! routine called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynhpg.F90,v 1.10 2005/12/21 10:46:37 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42) #if defined key_s_coord
(   43)    !!----------------------------------------------------------------------
(   44)    !!   'key_s_coord' :                                        s-coordinate
(   45)    !!----------------------------------------------------------------------   
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 2

(   47)    SUBROUTINE dyn_hpg( kt )
(   48)       !!---------------------------------------------------------------------
(   49)       !!                  ***  ROUTINE dyn_hpg  ***
(   50)       !!
(   51)       !! ** Purpose :   Compute the now momentum trend due to the hor. gradient
(   52)       !!      of the hydrostatic pressure. Add it to the general momentum trend.
(   53)       !!
(   54)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(   55)       !!      jk is computed by taking the vertical integral of the in-situ 
(   56)       !!      density gradient along the model level from the suface to that 
(   57)       !!      level. s-coordinates ('key_s_coord'): a corrective term is added
(   58)       !!      to the horizontal pressure gradient :
(   59)       !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
(   60)       !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
(   61)       !!      add it to the general momentum trend (ua,va).
(   62)       !!         ua = ua - 1/e1u * zhpi
(   63)       !!         va = va - 1/e2v * zhpj
(   64)       !!
(   65)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(   66)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(   67)       !!
(   68)       !! History :
(   69)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(   70)       !!        !  91-11  (G. Madec)
(   71)       !!        !  96-01  (G. Madec)  s-coordinates
(   72)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(   73)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module, vector opt.
(   74)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   75)       !!----------------------------------------------------------------------
(   76)       !! * modules used
(   77)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(   78)          &              zhpj => sa      ! use sa as 3D workspace
(   79) 
(   80)       !! * Arguments
(   81)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(   82)       
(   83)       !! * Local declarations
(   84)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(   85)       REAL(wp) ::   &
(   86)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(   87)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   88)          ztdua, ztdva                   ! temporary scalars
(   89)       !!----------------------------------------------------------------------
(   90) 
(   91)       IF( kt == nit000 ) THEN
(   92)          IF(lwp) WRITE(numout,*)
(   93)          IF(lwp) WRITE(numout,*) 'dyn_hpg : hydrostatic pressure gradient trend'
(   94)          IF(lwp) WRITE(numout,*) '~~~~~~~   s-coordinate case, vector opt. case'
(   95)       ENDIF
(   96) 
(   97)       ! Save ua and va trends
(   98)       IF( l_trddyn )   THEN
(   99)          ztdua(:,:,:) = ua(:,:,:) 
(  100)          ztdva(:,:,:) = va(:,:,:) 
(  101)       ENDIF
(  102) 
(  103)       ! 0. Local constant initialization
(  104)       ! --------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 3

(  105)       zcoef0 = - grav * 0.5
(  106)       zuap   = 0.e0
(  107)       zvap   = 0.e0
(  108) 
(  109)       ! 1. Surface value
(  110)       ! ----------------
(  111)       DO jj = 2, jpjm1
(  112)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  113)             ! hydrostatic pressure gradient along s-surfaces
(  114)             zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj)   &
(  115)                        * ( fse3w(ji+1,jj,1) * rhd(ji+1,jj,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
(  116)             zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj)   &
(  117)                        * ( fse3w(ji,jj+1,1) * rhd(ji,jj+1,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
(  118)             ! s-coordinate pressure gradient correction
(  119)             zuap = -zcoef0 * ( rhd(ji+1,jj,1) + rhd(ji,jj,1) )   &
(  120)                  * ( fsde3w(ji+1,jj,1) - fsde3w(ji,jj,1) ) / e1u(ji,jj)
(  121)             zvap = -zcoef0 * ( rhd(ji,jj+1,1) + rhd(ji,jj,1) )   &
(  122)                  * ( fsde3w(ji,jj+1,1) - fsde3w(ji,jj,1) ) / e2v(ji,jj)
(  123)             ! add to the general momentum trend
(  124)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
(  125)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
(  126)          END DO  
(  127)       END DO  
(  128) 
(  129)       ! 2. interior value (2=<jk=<jpkm1)
(  130)       ! -----------------
(  131)       DO jk = 2, jpkm1
(  132)          DO jj = 2, jpjm1 
(  133)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  134)                ! hydrostatic pressure gradient along s-surfaces
(  135)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
(  136)                   &           * ( fse3w(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
(  137)                   &              -fse3w(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  )
(  138)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
(  139)                   &           * ( fse3w(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
(  140)                   &              -fse3w(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  )
(  141)                ! s-coordinate pressure gradient correction 
(  142)                zuap = -zcoef0 * ( rhd(ji+1,jj  ,jk) + rhd(ji,jj,jk) )   &
(  143)                     * ( fsde3w(ji+1,jj,jk) - fsde3w(ji,jj,jk) ) / e1u(ji,jj)
(  144)                zvap = -zcoef0 * ( rhd(ji  ,jj+1,jk) + rhd(ji,jj,jk) )   &
(  145)                     * ( fsde3w(ji,jj+1,jk) - fsde3w(ji,jj,jk) ) / e2v(ji,jj)
(  146)                ! add to the general momentum trend
(  147)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
(  148)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
(  149)             END DO
(  150)          END DO
(  151)       END DO
(  152) 
(  153)       ! save the hydrostatic pressure gradient trends for diagnostic
(  154)       ! momentum trends
(  155)       IF( l_trddyn )   THEN
(  156)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  157)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  158)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
(  159)       ENDIF
(  160) 
(  161)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  162)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 4

(  163)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  164)       ENDIF
(  165) 
(  166)    END SUBROUTINE dyn_hpg
(  167) 
(  168) #elif defined key_partial_steps
(  169)    !!---------------------------------------------------------------------
(  170)    !!   'key_partial_steps'                     z-coordinate partial steps
(  171)    !!---------------------------------------------------------------------
(  172) 
(  173)    SUBROUTINE dyn_hpg( kt )
(  174)       !!---------------------------------------------------------------------
(  175)       !!                 ***  ROUTINE dyn_hpg  ***
(  176)       !!                    
(  177)       !! ** Purpose :   Compute the now momentum trend due to the horizontal 
(  178)       !!      gradient of the hydrostatic pressure. Add it to the general
(  179)       !!      momentum trend.
(  180)       !!
(  181)       !! ** Method  :   The now hydrostatic pressure gradient at a given level 
(  182)       !!      jk is computed by taking the vertical integral of the in-situ 
(  183)       !!      density gradient along the model level from the suface to that
(  184)       !!      level:   zhpi = grav .....
(  185)       !!               zhpj = grav .....
(  186)       !!      add it to the general momentum trend (ua,va).
(  187)       !!            ua = ua - 1/e1u * zhpi
(  188)       !!            va = va - 1/e2v * zhpj
(  189)       !!
(  190)       !! ** Action  : - Update (ua,va) with the now hydrastatic pressure trend
(  191)       !!              - Save the trend in (utrd,vtrd) ('key_trddyn')
(  192)       !!
(  193)       !! History :
(  194)       !!   8.5  !  02-08  (A. Bozec)  Original code
(  195)       !!----------------------------------------------------------------------
(  196)       !! * modules used
(  197)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  198)          &              zhpj => sa      ! use sa as 3D workspace
(  199) 
(  200)       !! * Arguments
(  201)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  202) 
(  203)       !! * local declarations
(  204)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  205)       INTEGER ::   iku, ikv             ! temporary integers
(  206)       REAL(wp) ::   &
(  207)          zcoef0, zcoef1, zuap,       &  ! temporary scalars
(  208)          zcoef2, zcoef3, zvap           !    "         "
(  209)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  210)          ztdua, ztdva                   ! temporary scalars
(  211)       !!----------------------------------------------------------------------
(  212) 
(  213)       IF( kt == nit000 ) THEN
(  214)          IF(lwp) WRITE(numout,*)
(  215)          IF(lwp) WRITE(numout,*) 'dyn_hpg : hydrostatic pressure gradient trend'
(  216)          IF(lwp) WRITE(numout,*) '~~~~~~~   z-coordinate with partial steps'
(  217)          IF(lwp) WRITE(numout,*) '          vector optimization, no autotasking'
(  218)       ENDIF
(  219) 
(  220)       ! Save ua and va trends






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 5

(  221)       IF( l_trddyn )   THEN
(  222)          ztdua(:,:,:) = ua(:,:,:) 
(  223)          ztdva(:,:,:) = va(:,:,:) 
(  224)       ENDIF
(  225) 
(  226)       ! 0. Local constant initialization
(  227)       ! --------------------------------
(  228)       zcoef0 = - grav * 0.5
(  229)       zuap   = 0.e0
(  230)       zvap   = 0.e0
(  231) 
(  232)       ! 1. Surface value
(  233)       ! ----------------
(  234)       DO jj = 2, jpjm1
(  235)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  236)             zcoef1 = zcoef0 * fse3w(ji,jj,1)
(  237)             ! hydrostatic pressure gradient
(  238)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  239)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  240)             ! add to the general momentum trend
(  241)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  242)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
(  243)          END DO
(  244)       END DO
(  245) 
(  246)       ! 2. interior value (2=<jk=<jpkm1)
(  247)       ! -----------------
(  248)       DO jk = 2, jpkm1
(  249)          DO jj = 2, jpjm1
(  250)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  251)                zcoef1 = zcoef0 * fse3w(ji,jj,jk)
(  252)                ! hydrostatic pressure gradient
(  253)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  254)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  255)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  256) 
(  257)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  258)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  259)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  260)                ! add to the general momentum trend
(  261)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  262)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  263)             END DO 
(  264)          END DO
(  265)       END DO
(  266) 
(  267)       ! partial steps correction at the last level  (new gradient with  intgrd.F)
(  268) # if defined key_vectopt_loop
(  269)          jj = 1
(  270)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  271) # else
(  272)       DO jj = 2, jpjm1
(  273)          DO ji = 2, jpim1
(  274) # endif
(  275)             iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj) ) - 1
(  276)             ikv = MIN ( mbathy(ji,jj), mbathy(ji,jj+1) ) - 1
(  277)             zcoef2 = zcoef0 * MIN( fse3w(ji,jj,iku), fse3w(ji+1,jj  ,iku) )
(  278)             zcoef3 = zcoef0 * MIN( fse3w(ji,jj,ikv), fse3w(ji  ,jj+1,ikv) )






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 6

(  279)             ! on i-direction
(  280)             IF ( iku > 2 ) THEN
(  281)                ! subtract old value  
(  282)                ua(ji,jj,iku) = ua(ji,jj,iku) - zhpi(ji,jj,iku)
(  283)                ! compute the new one   
(  284)                zhpi (ji,jj,iku) = zhpi(ji,jj,iku-1)   &
(  285)                   + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + gru(ji,jj) ) / e1u(ji,jj)
(  286)                ! add the new one to the general momentum trend
(  287)                ua(ji,jj,iku) = ua(ji,jj,iku) + zhpi(ji,jj,iku)
(  288)             ENDIF
(  289)             ! on j-direction
(  290)             IF ( ikv > 2 ) THEN
(  291)                ! subtract old value  
(  292)                va(ji,jj,ikv) = va(ji,jj,ikv) - zhpj(ji,jj,ikv)
(  293)                ! compute the new one   
(  294)                zhpj (ji,jj,ikv) = zhpj(ji,jj,ikv-1)   &
(  295)                   + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + grv(ji,jj) ) / e2v(ji,jj)
(  296)                ! add the new one to the general momentum trend
(  297)                va(ji,jj,ikv) = va(ji,jj,ikv) + zhpj(ji,jj,ikv)
(  298)             ENDIF
(  299) # if ! defined key_vectopt_loop
(  300)          END DO
(  301) # endif
(  302)       END DO
(  303) 
(  304)       ! save the hydrostatic pressure gradient trends for diagnostic
(  305)       ! momentum trends
(  306)       IF( l_trddyn )   THEN
(  307)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  308)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  309)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
(  310)       ENDIF
(  311) 
(  312)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  313)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
(  314)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  315)       ENDIF
(  316) 
(  317)    END SUBROUTINE dyn_hpg
(  318) 
(  319) #else
(  320)    !!---------------------------------------------------------------------
(  321)    !!   Default case :                                        z-coordinate
(  322)    !!---------------------------------------------------------------------
(  323) 
(  324)    SUBROUTINE dyn_hpg( kt )
(  325)       !!---------------------------------------------------------------------
(  326)       !!                  ***  ROUTINE dyn_hpg  ***
(  327)       !!
(  328)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(  329)       !!      gradient of the hydrostatic pressure. Add it to the general 
(  330)       !!      momentum trend.
(  331)       !!
(  332)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(  333)       !!      jk is computed by taking the vertical integral of the in-situ
(  334)       !!      density gradient along the model level from the suface to that
(  335)       !!      level:    zhpi = grav .....
(  336)       !!                zhpj = grav .....






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 7

(  337)       !!      add it to the general momentum trend (ua,va).
(  338)       !!            ua = ua - 1/e1u * zhpi
(  339)       !!            va = va - 1/e2v * zhpj
(  340)       !!
(  341)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(  342)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(  343)       !!
(  344)       !! History :
(  345)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(  346)       !!        !  91-11  (G. Madec)
(  347)       !!        !  96-01  (G. Madec)  s-coordinates
(  348)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(  349)       !!   8.5  !  02-07  (G. Madec)  F90: Free form and module
(  350)       !!----------------------------------------------------------------------
(  351)       !! * modules used
(  352)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  353)          &              zhpj => sa      ! use sa as 3D workspace
(  354) 
(  355)       !! * Arguments
(  356)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  357) 
(  358)       !! * local declarations
(  359)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  360)       REAL(wp) ::   &
(  361)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(  362)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  363)          ztdua, ztdva                   ! temporary scalars
(  364)       !!----------------------------------------------------------------------
(  365) 
(  366)       IF( kt == nit000 ) THEN
(  367)          IF(lwp) WRITE(numout,*)
(  368)          IF(lwp) WRITE(numout,*) 'dyn_hpg : hydrostatic pressure gradient trend'
(  369)          IF(lwp) WRITE(numout,*) '~~~~~~~   z-coordinate case '
(  370)       ENDIF
(  371) 
(  372)       ! Save ua and va trends
(  373)       IF( l_trddyn )   THEN
(  374)          ztdua(:,:,:) = ua(:,:,:) 
(  375)          ztdva(:,:,:) = va(:,:,:) 
(  376)       ENDIF
(  377) 
(  378)       ! 0. Local constant initialization
(  379)       ! --------------------------------
(  380)       zcoef0 = - grav * 0.5
(  381)       zuap   = 0.e0
(  382)       zvap   = 0.e0
(  383) 
(  384)       ! 1. Surface value
(  385)       ! ----------------
(  386)       DO jj = 2, jpjm1
(  387)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  388)             zcoef1 = zcoef0 * fse3w(ji,jj,1)
(  389)             ! hydrostatic pressure gradient
(  390)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  391)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  392)             ! add to the general momentum trend
(  393)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  394)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 8

(  395)          END DO
(  396)       END DO
(  397) 
(  398)       ! 2. interior value (2=<jk=<jpkm1)
(  399)       ! -----------------
(  400)       DO jk = 2, jpkm1
(  401)          DO jj = 2, jpjm1
(  402)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  403)                zcoef1 = zcoef0 * fse3w(ji,jj,jk)
(  404)                ! hydrostatic pressure gradient
(  405)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  406)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  407)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  408) 
(  409)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  410)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  411)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  412)                ! add to the general momentum trend
(  413)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  414)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  415)             END DO 
(  416)          END DO
(  417)       END DO
(  418) 
(  419)       ! save the hydrostatic pressure ggradient trends for diagnostic
(  420)       ! momentum trends
(  421)       IF( l_trddyn )   THEN
(  422)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  423)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  424) 
(  425)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
(  426)       ENDIF
(  427) 
(  428)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  429)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
(  430)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  431)       ENDIF
(  432) 
(  433)    END SUBROUTINE dyn_hpg
(  434) 
(  435) #endif
(  436) 
(  437)    !!======================================================================
(  438) END MODULE dynhpg




















PGF90 (Version     10.5)          09/23/2010  09:49:36      page 9

(    1) # 1 "dynhpg.F90"
(    1) MODULE dynhpg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynhpg  ***
(    4)    !! Ocean dynamics:  hydrostatic pressure gradient trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_hpg      : update the momentum trend with the horizontal
(    9)    !!                  gradient of the hydrostatic pressure
(   10)    !!
(   11)    !!   default case : use of 3D work arrays (vector opt. available)
(   12)    !!   key_s_coord       : s-coordinate
(   13)    !!   key_partial_steps : z-coordinate with partial steps
(   14)    !!   default key       : z-coordinate
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17)    USE oce             ! ocean dynamics and tracers
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE phycst          ! physical constants
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE trdmod          ! ocean dynamics trends 
(   22)    USE trdmod_oce      ! ocean variables trends
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Accessibility
(   29)    PUBLIC dyn_hpg                ! routine called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) # 1 "./domzgr_substitute.h90"
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
(   48) # 33 "dynhpg.F90"
(   33) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 10

(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 34 "dynhpg.F90"
(   34) # 34
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynhpg.F90,v 1.10 2005/12/21 10:46:37 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42) # 320
(  320)    !!---------------------------------------------------------------------
(  321)    !!   Default case :                                        z-coordinate
(  322)    !!---------------------------------------------------------------------
(  323) 
(  324)    SUBROUTINE dyn_hpg( kt )
(  325)       !!---------------------------------------------------------------------
(  326)       !!                  ***  ROUTINE dyn_hpg  ***
(  327)       !!
(  328)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(  329)       !!      gradient of the hydrostatic pressure. Add it to the general 
(  330)       !!      momentum trend.
(  331)       !!
(  332)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(  333)       !!      jk is computed by taking the vertical integral of the in-situ
(  334)       !!      density gradient along the model level from the suface to that
(  335)       !!      level:    zhpi = grav .....
(  336)       !!                zhpj = grav .....
(  337)       !!      add it to the general momentum trend (ua,va).
(  338)       !!            ua = ua - 1/e1u * zhpi
(  339)       !!            va = va - 1/e2v * zhpj
(  340)       !!
(  341)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(  342)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(  343)       !!
(  344)       !! History :
(  345)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(  346)       !!        !  91-11  (G. Madec)
(  347)       !!        !  96-01  (G. Madec)  s-coordinates
(  348)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(  349)       !!   8.5  !  02-07  (G. Madec)  F90: Free form and module
(  350)       !!----------------------------------------------------------------------
(  351)       !! * modules used
(  352)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  353)          &              zhpj => sa      ! use sa as 3D workspace
(  354) 
(  355)       !! * Arguments
(  356)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  357) 
(  358)       !! * local declarations
(  359)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  360)       REAL(wp) ::   &
(  361)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(  362)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  363)          ztdua, ztdva                   ! temporary scalars
(  364)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 11

(  365) 
(  366)       IF( kt == nit000 ) THEN
(  367)          IF(lwp) WRITE(numout,*)
(  368)          IF(lwp) WRITE(numout,*) 'dyn_hpg : hydrostatic pressure gradient trend'
(  369)          IF(lwp) WRITE(numout,*) '~~~~~~~   z-coordinate case '
(  370)       ENDIF
(  371) 
(  372)       ! Save ua and va trends
(  373)       IF( l_trddyn )   THEN
(  374)          ztdua(:,:,:) = ua(:,:,:) 
(  375)          ztdva(:,:,:) = va(:,:,:) 
(  376)       ENDIF
(  377) 
(  378)       ! 0. Local constant initialization
(  379)       ! --------------------------------
(  380)       zcoef0 = - grav * 0.5
(  381)       zuap   = 0.e0
(  382)       zvap   = 0.e0
(  383) 
(  384)       ! 1. Surface value
(  385)       ! ----------------
(  386)       DO jj = 2, jpjm1
(  387)          DO ji = 2, jpim1   ! vector opt.
(  388)             zcoef1 = zcoef0 * e3w(1)
(  389)             ! hydrostatic pressure gradient
(  390)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  391)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  392)             ! add to the general momentum trend
(  393)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  394)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
(  395)          END DO
(  396)       END DO
(  397) 
(  398)       ! 2. interior value (2=<jk=<jpkm1)
(  399)       ! -----------------
(  400)       DO jk = 2, jpkm1
(  401)          DO jj = 2, jpjm1
(  402)             DO ji = 2, jpim1   ! vector opt.
(  403)                zcoef1 = zcoef0 * e3w(jk)
(  404)                ! hydrostatic pressure gradient
(  405)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  406)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  407)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  408) 
(  409)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  410)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  411)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  412)                ! add to the general momentum trend
(  413)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  414)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  415)             END DO 
(  416)          END DO
(  417)       END DO
(  418) 
(  419)       ! save the hydrostatic pressure ggradient trends for diagnostic
(  420)       ! momentum trends
(  421)       IF( l_trddyn )   THEN
(  422)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 12

(  423)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  424) 
(  425)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
(  426)       ENDIF
(  427) 
(  428)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  429)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
(  430)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  431)       ENDIF
(  432) 
(  433)    END SUBROUTINE dyn_hpg
(  434) 
(  435) 
(  436) # 437
(  437)    !!======================================================================
(  438) END MODULE dynhpg
