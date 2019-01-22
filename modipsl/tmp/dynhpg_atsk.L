


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

Filename: dynhpg_atsk.F90

(    1) MODULE dynhpg_atsk
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynhpg_atsk  ***
(    4)    !! Ocean dynamics:  hydrostatic pressure gradient trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_hpg_atsk : update the momentum trend with the horizontal
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
(   29)    PUBLIC dyn_hpg_atsk ! routine called by step.F90
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynhpg_atsk.F90,v 1.9 2005/09/02 15:45:23 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42) #if defined key_s_coord
(   43)    !!---------------------------------------------------------------------
(   44)    !!                  ***  dynhpg_atsk.h90  ***
(   45)    !!---------------------------------------------------------------------
(   46)    !!   'key_s_coord'                                         s-coordinate






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 2

(   47)    !!---------------------------------------------------------------------
(   48) 
(   49)    SUBROUTINE dyn_hpg_atsk( kt )
(   50)       !!---------------------------------------------------------------------
(   51)       !!                  ***  ROUTINE dyn_hpg_atsk  ***
(   52)       !!        
(   53)       !! ** Purpose :   Compute the now momentum trend due to the horizontal 
(   54)       !!     gradient of the hydrostatic pressure. Add it to the general
(   55)       !!     momentum trend.
(   56)       !!
(   57)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(   58)       !!      jk is computed by taking the vertical integral of the in-situ 
(   59)       !!      density gradient along the model level from the suface to that
(   60)       !!      level. s-coordinate case ('key_s_coord'): a corrective term is
(   61)       !!      added to the horizontal pressure gradient :
(   62)       !!         zhpi = grav .....   + 1/e1u mi(rhd) di[ grav dep3w ]
(   63)       !!         zhpj = grav .....   + 1/e2v mj(rhd) dj[ grav dep3w ]
(   64)       !!      add it to the general momentum trend (ua,va).
(   65)       !!         ua = ua - 1/e1u * zhpi
(   66)       !!         va = va - 1/e2v * zhpj
(   67)       !!      j-k-i loop (j-slab) ('key_autotasking')
(   68)       !!
(   69)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(   70)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(   71)       !!
(   72)       !! History :
(   73)       !!   7.0  !  96-01  (G. Madec)  s-coordinates
(   74)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(   75)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   76)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   77)       !!----------------------------------------------------------------------
(   78)       !! * modules used
(   79)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(   80)          &              zhpj => sa      ! use sa as 3D workspace
(   81) 
(   82)       !! * Arguments
(   83)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(   84)       
(   85)       !! * Local declarations
(   86)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(   87)       REAL(wp) ::   &
(   88)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(   89)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   90)          ztdua, ztdva                   ! temporary scalars
(   91)       !!----------------------------------------------------------------------
(   92) 
(   93)       IF( kt == nit000 ) THEN
(   94)          IF(lwp) WRITE(numout,*)
(   95)          IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : s-coordinate hydrostatic pressure gradient trend'
(   96)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   autotasking case (j-k-i loop)'
(   97)       ENDIF
(   98) 
(   99)       ! Save ua and va trends
(  100)       IF( l_trddyn )   THEN
(  101)          ztdua(:,:,:) = ua(:,:,:) 
(  102)          ztdva(:,:,:) = va(:,:,:) 
(  103)       ENDIF
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 3

(  105)       ! 0. Local constant initialization
(  106)       ! --------------------------------
(  107)       zcoef0 = - grav * 0.5
(  108)       zuap   = 0.e0
(  109)       zvap   = 0.e0
(  110) 
(  111)       !                                                ! ===============
(  112)       DO jj = 2, jpjm1                                 !  Vertical slab
(  113)          !                                             ! ===============
(  114)          ! 1. Surface value
(  115)          ! ----------------
(  116)          DO ji = 2, jpim1
(  117)             ! hydrostatic pressure gradient along s-surfaces
(  118)             zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj)   &
(  119)                        * ( fse3w(ji+1,jj,1) * rhd(ji+1,jj,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
(  120)             zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj)   &
(  121)                        * ( fse3w(ji,jj+1,1) * rhd(ji,jj+1,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
(  122)             ! s-coordinate pressure gradient correction
(  123)             zuap = -zcoef0 * ( rhd(ji+1,jj,1) + rhd(ji,jj,1) )   &
(  124)                  * ( fsde3w(ji+1,jj,1) - fsde3w(ji,jj,1) ) / e1u(ji,jj)
(  125)             zvap = -zcoef0 * ( rhd(ji,jj+1,1) + rhd(ji,jj,1) )   &
(  126)                  * ( fsde3w(ji,jj+1,1) - fsde3w(ji,jj,1) ) / e2v(ji,jj)
(  127)             ! add to the general momentum trend
(  128)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
(  129)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
(  130)          END DO  
(  131) 
(  132)          ! 2. interior value (2=<jk=<jpkm1)
(  133)          ! -----------------
(  134)          DO jk = 2, jpkm1
(  135)             DO ji = 2, jpim1
(  136)                ! hydrostatic pressure gradient along s-surfaces
(  137)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
(  138)                   &           * ( fse3w(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
(  139)                   &              -fse3w(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  )
(  140)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
(  141)                   &           * ( fse3w(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
(  142)                   &              -fse3w(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  )
(  143)                ! s-coordinate pressure gradient correction 
(  144)                zuap = -zcoef0 * ( rhd(ji+1,jj  ,jk) + rhd(ji,jj,jk) )   &
(  145)                     * ( fsde3w(ji+1,jj,jk) - fsde3w(ji,jj,jk) ) / e1u(ji,jj)
(  146)                zvap = -zcoef0 * ( rhd(ji  ,jj+1,jk) + rhd(ji,jj,jk) )   &
(  147)                     * ( fsde3w(ji,jj+1,jk) - fsde3w(ji,jj,jk) ) / e2v(ji,jj)
(  148)                ! add to the general momentum trend
(  149)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
(  150)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
(  151)             END DO
(  152)          END DO
(  153)          !                                             ! ===============
(  154)       END DO                                           !   End of slab
(  155)       !                                                ! ===============
(  156) 
(  157)       ! save the hydrostatic pressure gradient trends for diagnostic
(  158)       ! momentum trends
(  159)       IF( l_trddyn )   THEN
(  160)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  161)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  162)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 4

(  163)       ENDIF
(  164) 
(  165)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  166)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
(  167)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  168)       ENDIF
(  169) 
(  170)    END SUBROUTINE dyn_hpg_atsk
(  171) 
(  172) #elif defined key_partial_steps
(  173)    !!---------------------------------------------------------------------
(  174)    !!   'key_partial_steps'                     z-coordinate partial steps
(  175)    !!---------------------------------------------------------------------
(  176) 
(  177)    SUBROUTINE dyn_hpg_atsk( kt )
(  178)       !!---------------------------------------------------------------------
(  179)       !!                  ***  ROUTINE dyn_hpg_atsk  ***
(  180)       !!  
(  181)       !! ** Purpose :   Compute the now momentum trend due to the hor. gradient
(  182)       !!      of the hydrostatic pressure. Add it to the general momentum trend.
(  183)       !!
(  184)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(  185)       !!      jk is computed by taking the vertical integral of the in-situ
(  186)       !!      density gradient along the model level from the suface to that 
(  187)       !!      level:    zhpi = grav .....
(  188)       !!                zhpj = grav .....
(  189)       !!      add it to the general momentum trend (ua,va).
(  190)       !!            ua = ua - 1/e1u * zhpi
(  191)       !!            va = va - 1/e2v * zhpj
(  192)       !!      j-k-i loop (j-slab) ('key_autotasking')
(  193)       !!
(  194)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(  195)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(  196)       !!
(  197)       !! History :
(  198)       !!   8.5  !  02-08  (A. Bozec)  Original code
(  199)       !!----------------------------------------------------------------------
(  200)       !! * modules used
(  201)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  202)          &              zhpj => sa      ! use sa as 3D workspace
(  203) 
(  204)       !! * Arguments
(  205)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  206) 
(  207)       !! * local declarations
(  208)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  209)       INTEGER ::   iku, ikv             ! temporary integers
(  210)       REAL(wp) ::   &
(  211)          zcoef0, zcoef1, zuap,       &  ! temporary scalars
(  212)          zcoef2, zcoef3, zvap           !    "         "
(  213)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  214)          ztdua, ztdva                   ! temporary scalars
(  215)       !!----------------------------------------------------------------------
(  216) 
(  217)       IF( kt == nit000 ) THEN
(  218)          IF(lwp) WRITE(numout,*)
(  219)          IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : z-coord. partial steps hydrostatic pressure gradient trend'
(  220)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   autotasking case (j-k-i loop)'






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 5

(  221)       ENDIF
(  222) 
(  223)       ! Save ua and va trends
(  224)       IF( l_trddyn )   THEN
(  225)          ztdua(:,:,:) = ua(:,:,:) 
(  226)          ztdva(:,:,:) = va(:,:,:) 
(  227)       ENDIF
(  228) 
(  229)       ! 0. Local constant initialization
(  230)       ! --------------------------------
(  231)       zcoef0 = - grav * 0.5
(  232)       zuap   = 0.e0
(  233)       zvap   = 0.e0
(  234)       !                                                ! ===============
(  235)       DO jj = 2, jpjm1                                 !  Vertical slab
(  236)          !                                             ! ===============
(  237)          ! 1. Surface value
(  238)          ! ----------------
(  239)          DO ji = 2, jpim1
(  240)             zcoef1 = zcoef0 * fse3w(ji,jj,1)
(  241)             ! hydrostatic pressure gradient
(  242)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  243)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  244)             ! add to the general momentum trend
(  245)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  246)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
(  247)          END DO
(  248) 
(  249)          ! 2. interior value (2=<jk=<jpkm1)
(  250)          ! -----------------
(  251)          DO jk = 2, jpkm1
(  252)             DO ji = 2, jpim1
(  253)                zcoef1 = zcoef0 * fse3w(ji,jj,jk)
(  254)                ! hydrostatic pressure gradient
(  255)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  256)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  257)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  258) 
(  259)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  260)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  261)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  262)                ! add to the general momentum trend
(  263)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  264)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  265)             END DO 
(  266)          END DO
(  267) 
(  268)          ! partial steps correction at the last level  (new gradient with  intgrd.F)
(  269)          DO ji = 2, jpim1
(  270)             iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj) ) - 1
(  271)             ikv = MIN ( mbathy(ji,jj), mbathy(ji,jj+1) ) - 1
(  272)             zcoef2 = zcoef0 * MIN( fse3w(ji,jj,iku), fse3w(ji+1,jj  ,iku) )
(  273)             zcoef3 = zcoef0 * MIN( fse3w(ji,jj,ikv), fse3w(ji  ,jj+1,ikv) )
(  274)             ! on i-direction
(  275)             IF ( iku > 2 ) THEN
(  276)                ! subtract old value  
(  277)                ua(ji,jj,iku) = ua(ji,jj,iku) - zhpi(ji,jj,iku)
(  278)                ! compute the new one   






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 6

(  279)                zhpi (ji,jj,iku) = zhpi(ji,jj,iku-1)   &
(  280)                   + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + gru(ji,jj) ) / e1u(ji,jj)
(  281)                ! add the new one to the general momentum trend
(  282)                ua(ji,jj,iku) = ua(ji,jj,iku) + zhpi(ji,jj,iku)
(  283)             ENDIF
(  284)             ! on j-direction
(  285)             IF ( ikv > 2 ) THEN
(  286)                ! subtract old value  
(  287)                va(ji,jj,ikv) = va(ji,jj,ikv) - zhpj(ji,jj,ikv)
(  288)                ! compute the new one   
(  289)                zhpj (ji,jj,ikv) = zhpj(ji,jj,ikv-1)   &
(  290)                   + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + grv(ji,jj) ) / e2v(ji,jj)
(  291)                ! add the new one to the general momentum trend
(  292)                va(ji,jj,ikv) = va(ji,jj,ikv) + zhpj(ji,jj,ikv)
(  293)             ENDIF
(  294)          END DO
(  295)          !                                             ! ===============
(  296)       END DO                                           !   End of slab
(  297)       !                                                ! ===============
(  298) 
(  299)       ! save the hydrostatic pressure gradient trends for diagnostic
(  300)       ! momentum trends
(  301)       IF( l_trddyn )   THEN
(  302)          zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  303)          zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  304)          CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
(  305)       ENDIF
(  306) 
(  307)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  308)          CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
(  309)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  310)       ENDIF   
(  311) 
(  312)    END SUBROUTINE dyn_hpg_atsk
(  313) 
(  314) #else
(  315)    !!---------------------------------------------------------------------
(  316)    !!   Default case :                                        z-coordinate
(  317)    !!---------------------------------------------------------------------
(  318) 
(  319)    SUBROUTINE dyn_hpg_atsk( kt )
(  320)       !!---------------------------------------------------------------------
(  321)       !!                  ***  ROUTINE dyn_hpg_atsk  ***
(  322)       !!   
(  323)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(  324)       !!      gradient of the hydrostatic pressure. Add it to the general
(  325)       !!      momentum trend.
(  326)       !!
(  327)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(  328)       !!      jk is computed by taking the vertical integral of the in-situ  
(  329)       !!      density gradient along the model level from the suface to that
(  330)       !!      level:    zhpi = grav .....
(  331)       !!                zhpj = grav .....
(  332)       !!      add it to the general momentum trend (ua,va).
(  333)       !!            ua = ua - 1/e1u * zhpi
(  334)       !!            va = va - 1/e2v * zhpj
(  335)       !!      j-k-i loop (j-slab) ('key_autotasking')
(  336)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 7

(  337)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(  338)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(  339)       !!
(  340)       !! History :
(  341)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(  342)       !!        !  91-11  (G. Madec)
(  343)       !!        !  96-01  (G. Madec)  s-coordinates
(  344)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(  345)       !!   8.5  !  02-07  (G. Madec)  F90: Free form and module
(  346)       !!----------------------------------------------------------------------
(  347)       !! * modules used
(  348)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  349)          &              zhpj => sa      ! use sa as 3D workspace
(  350) 
(  351)       !! * Arguments
(  352)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  353) 
(  354)       !! * local declarations
(  355)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  356)       REAL(wp) ::   &
(  357)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(  358)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  359)          ztdua, ztdva                   ! temporary scalars
(  360)       !!----------------------------------------------------------------------
(  361) 
(  362)       IF( kt == nit000 ) THEN
(  363)          IF(lwp) WRITE(numout,*)
(  364)          IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : z-coordinate hydrostatic pressure gradient trend'
(  365)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   auto-tasking (j-k-i loop)'
(  366)       ENDIF
(  367) 
(  368)       ! Save ua and va trends
(  369)       IF( l_trddyn )   THEN
(  370)          ztdua(:,:,:) = ua(:,:,:) 
(  371)          ztdva(:,:,:) = va(:,:,:) 
(  372)       ENDIF
(  373) 
(  374)       ! 0. Local constant initialization
(  375)       ! --------------------------------
(  376)       zcoef0 = - grav * 0.5
(  377)       zuap   = 0.e0
(  378)       zvap   = 0.e0
(  379) 
(  380)       !                                                ! ===============
(  381)       DO jj = 2, jpjm1                                 !  Vertical slab
(  382)          !                                             ! ===============
(  383)          ! 1. Surface value
(  384)          ! ----------------
(  385)          
(  386) 
(  387)          DO ji = 2, jpim1
(  388)             zcoef1 = zcoef0 * fse3w(ji,jj,1)
(  389)             ! hydrostatic pressure gradient
(  390)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  391)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  392)             ! add to the general momentum trend
(  393)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  394)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 8

(  395)          END DO
(  396) 
(  397)          ! 2. interior value (2=<jk=<jpkm1)
(  398)          ! -----------------
(  399)          DO jk = 2, jpkm1
(  400)             DO ji = 2, jpim1
(  401)                zcoef1 = zcoef0 * fse3w(ji,jj,jk)
(  402)                ! hydrostatic pressure gradient
(  403)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  404)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  405)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  406) 
(  407)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  408)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  409)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  410)                ! add to the general momentum trend
(  411)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  412)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  413)             END DO 
(  414)          END DO
(  415)          !                                             ! ===============
(  416)       END DO                                           !   End of slab
(  417)       !                                                ! ===============
(  418) 
(  419)       ! save the hydrostatic pressure gradient trends for diagnostic
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
(  433)    END SUBROUTINE dyn_hpg_atsk
(  434) 
(  435) #endif
(  436) 
(  437)    !!======================================================================
(  438) END MODULE dynhpg_atsk




















PGF90 (Version     10.5)          09/23/2010  09:49:36      page 9

(    1) # 1 "dynhpg_atsk.F90"
(    1) MODULE dynhpg_atsk
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynhpg_atsk  ***
(    4)    !! Ocean dynamics:  hydrostatic pressure gradient trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_hpg_atsk : update the momentum trend with the horizontal
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
(   29)    PUBLIC dyn_hpg_atsk ! routine called by step.F90
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
(   48) # 33 "dynhpg_atsk.F90"
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
(   12) # 34 "dynhpg_atsk.F90"
(   34) # 34
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynhpg_atsk.F90,v 1.9 2005/09/02 15:45:23 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42) # 315
(  315)    !!---------------------------------------------------------------------
(  316)    !!   Default case :                                        z-coordinate
(  317)    !!---------------------------------------------------------------------
(  318) 
(  319)    SUBROUTINE dyn_hpg_atsk( kt )
(  320)       !!---------------------------------------------------------------------
(  321)       !!                  ***  ROUTINE dyn_hpg_atsk  ***
(  322)       !!   
(  323)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(  324)       !!      gradient of the hydrostatic pressure. Add it to the general
(  325)       !!      momentum trend.
(  326)       !!
(  327)       !! ** Method  :   The now hydrostatic pressure gradient at a given level
(  328)       !!      jk is computed by taking the vertical integral of the in-situ  
(  329)       !!      density gradient along the model level from the suface to that
(  330)       !!      level:    zhpi = grav .....
(  331)       !!                zhpj = grav .....
(  332)       !!      add it to the general momentum trend (ua,va).
(  333)       !!            ua = ua - 1/e1u * zhpi
(  334)       !!            va = va - 1/e2v * zhpj
(  335)       !!      j-k-i loop (j-slab) ('key_autotasking')
(  336)       !!
(  337)       !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
(  338)       !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
(  339)       !!
(  340)       !! History :
(  341)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(  342)       !!        !  91-11  (G. Madec)
(  343)       !!        !  96-01  (G. Madec)  s-coordinates
(  344)       !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
(  345)       !!   8.5  !  02-07  (G. Madec)  F90: Free form and module
(  346)       !!----------------------------------------------------------------------
(  347)       !! * modules used
(  348)       USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
(  349)          &              zhpj => sa      ! use sa as 3D workspace
(  350) 
(  351)       !! * Arguments
(  352)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(  353) 
(  354)       !! * local declarations
(  355)       INTEGER ::   ji, jj, jk           ! dummy loop indices
(  356)       REAL(wp) ::   &
(  357)          zcoef0, zcoef1, zuap, zvap     ! temporary scalars
(  358)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  359)          ztdua, ztdva                   ! temporary scalars






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 11

(  360)       !!----------------------------------------------------------------------
(  361) 
(  362)       IF( kt == nit000 ) THEN
(  363)          IF(lwp) WRITE(numout,*)
(  364)          IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : z-coordinate hydrostatic pressure gradient trend'
(  365)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   auto-tasking (j-k-i loop)'
(  366)       ENDIF
(  367) 
(  368)       ! Save ua and va trends
(  369)       IF( l_trddyn )   THEN
(  370)          ztdua(:,:,:) = ua(:,:,:) 
(  371)          ztdva(:,:,:) = va(:,:,:) 
(  372)       ENDIF
(  373) 
(  374)       ! 0. Local constant initialization
(  375)       ! --------------------------------
(  376)       zcoef0 = - grav * 0.5
(  377)       zuap   = 0.e0
(  378)       zvap   = 0.e0
(  379) 
(  380)       !                                                ! ===============
(  381)       DO jj = 2, jpjm1                                 !  Vertical slab
(  382)          !                                             ! ===============
(  383)          ! 1. Surface value
(  384)          ! ----------------
(  385)          
(  386) 
(  387)          DO ji = 2, jpim1
(  388)             zcoef1 = zcoef0 * e3w(1)
(  389)             ! hydrostatic pressure gradient
(  390)             zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
(  391)             zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
(  392)             ! add to the general momentum trend
(  393)             ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
(  394)             va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
(  395)          END DO
(  396) 
(  397)          ! 2. interior value (2=<jk=<jpkm1)
(  398)          ! -----------------
(  399)          DO jk = 2, jpkm1
(  400)             DO ji = 2, jpim1
(  401)                zcoef1 = zcoef0 * e3w(jk)
(  402)                ! hydrostatic pressure gradient
(  403)                zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
(  404)                   &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
(  405)                   &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)
(  406) 
(  407)                zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
(  408)                   &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
(  409)                   &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
(  410)                ! add to the general momentum trend
(  411)                ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
(  412)                va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
(  413)             END DO 
(  414)          END DO
(  415)          !                                             ! ===============
(  416)       END DO                                           !   End of slab
(  417)       !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 12

(  418) 
(  419)       ! save the hydrostatic pressure gradient trends for diagnostic
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
(  433)    END SUBROUTINE dyn_hpg_atsk
(  434) 
(  435) 
(  436) # 437
(  437)    !!======================================================================
(  438) END MODULE dynhpg_atsk
