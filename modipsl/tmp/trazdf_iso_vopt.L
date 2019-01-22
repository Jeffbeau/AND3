


PGF90 (Version     10.5)          09/23/2010  09:50:18      page 1

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

Filename: trazdf_iso_vopt.F90

(    1) MODULE trazdf_iso_vopt
(    2)    !!==============================================================================
(    3)    !!                 ***  MODULE  trazdf_iso_vopt  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_ldfslp   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_zdf_iso_vopt : Update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion. vector optimization, use
(   13)    !!                  k-j-i loops.
(   14)    !!   tra_zdf_iso  :
(   15)    !!   tra_zdf_zdf  :
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce             ! ocean dynamics and tracers variables
(   19)    USE dom_oce         ! ocean space and time domain variables 
(   20)    USE zdf_oce         ! ocean vertical physics variables
(   21)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   22)    USE trdmod          ! ocean active tracers trends 
(   23)    USE trdmod_oce      ! ocean variables trends
(   24)    USE ldfslp          ! iso-neutral slopes 
(   25)    USE zdfddm          ! ocean vertical physics: double diffusion
(   26)    USE in_out_manager  ! I/O manager
(   27)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   28)    USE zdfkpp          ! KPP parameterisation
(   29)    USE prtctl          ! Print control
(   30) 
(   31)    IMPLICIT NONE
(   32)    PRIVATE
(   33) 
(   34)    !! * Routine accessibility
(   35)    PUBLIC tra_zdf_iso_vopt   !  routine called by step.F90
(   36) 
(   37)    !! * Module variables
(   38)    REAL(wp), DIMENSION(jpk) ::  &
(   39)       r2dt                          ! vertical profile of 2 x time-step
(   40)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
(   41)       tavg, savg,                     & ! workspace arrays
(   42)       tdta, tdsa                        ! workspace arrays
(   43) 
(   44)    !! * Substitutions
(   45) #  include "domzgr_substitute.h90"
(   46) #  include "ldftra_substitute.h90"






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 2

(   47) #  include "ldfeiv_substitute.h90"
(   48) #  include "zdfddm_substitute.h90"
(   49) #  include "vectopt_loop_substitute.h90"
(   50)    !!----------------------------------------------------------------------
(   51)    !!----------------------------------------------------------------------
(   52)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   53)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso_vopt.F90,v 1.10 2005/09/22 10:29:06 opalod Exp $ 
(   54)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   55)    !!----------------------------------------------------------------------
(   56) CONTAINS
(   57)    
(   58)    SUBROUTINE tra_zdf_iso_vopt( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                  ***  ROUTINE tra_zdf_iso_vopt  ***
(   61)       !!
(   62)       !! ** Purpose :
(   63)       !! ** Method  :
(   64)       !! ** Action  :
(   65)       !!
(   66)       !! History :
(   67)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   68)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   69)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   70)       !!---------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   73) 
(   74)       IF( kt == nit000 ) THEN
(   75)          IF(lwp)WRITE(numout,*)
(   76)          IF(lwp)WRITE(numout,*) 'tra_zdf_iso_vopt : vertical mixing computation'
(   77)          IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~~   is  iso-neutral diffusion : implicit vertical time stepping'
(   78) #if defined key_diaeiv 
(   79)          w_eiv(:,:,:) = 0.e0
(   80) #endif
(   81)       ENDIF
(   82) 
(   83)       ! initialization step
(   84)       tavg(:,:,:) = 0.e0
(   85)       savg(:,:,:) = 0.e0
(   86) 
(   87)       ! I. vertical extra-diagonal part of the rotated tensor
(   88)       ! -----------------------------------------------------
(   89) 
(   90)       CALL tra_zdf_iso( kt )
(   91) 
(   92)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(   93)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
(   94)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(   95)       ENDIF
(   96) 
(   97)       ! II. vertical diffusion (including the vertical diagonal part of the rotated tensor)
(   98)       ! ----------------------
(   99) 
(  100)       CALL tra_zdf_zdf( kt )
(  101) 
(  102)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  103)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
(  104)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 3

(  105)       ENDIF
(  106) 
(  107)    END SUBROUTINE tra_zdf_iso_vopt
(  108) 
(  109) 
(  110)    SUBROUTINE tra_zdf_zdf( kt )
(  111)       !!----------------------------------------------------------------------
(  112)       !!                 ***  ROUTINE tra_zdf_zdf  ***
(  113)       !!                    
(  114)       !! ** Purpose :   Compute the trend due to the vertical tracer diffusion
(  115)       !!     including the vertical component of lateral mixing (only for 2nd
(  116)       !!     order operator, for fourth order it is already computed and add
(  117)       !!     to the general trend in traldf.F) and add it to the general trend
(  118)       !!     of the tracer equations.
(  119)       !!
(  120)       !! ** Method  :   The vertical component of the lateral diffusive trends
(  121)       !!      is provided by a 2nd order operator rotated along neural or geo-
(  122)       !!      potential surfaces to which an eddy induced advection can be 
(  123)       !!      added. It is computed using before fields (forward in time) and 
(  124)       !!      isopycnal or geopotential slopes computed in routine ldfslp.
(  125)       !!
(  126)       !!      Second part: vertical trend associated with the vertical physics
(  127)       !!      ===========  (including the vertical flux proportional to dk[t]
(  128)       !!                  associated with the lateral mixing, through the
(  129)       !!                  update of avt)
(  130)       !!      The vertical diffusion of tracers (t & s) is given by:
(  131)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(  132)       !!      It is computed using a backward time scheme (t=ta).
(  133)       !!      Surface and bottom boundary conditions: no diffusive flux on
(  134)       !!      both tracers (bottom, applied through the masked field avt).
(  135)       !!      Add this trend to the general trend ta,sa :
(  136)       !!         ta = ta + dz( avt dz(t) )
(  137)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
(  138)       !!
(  139)       !!      Third part: recover avt resulting from the vertical physics
(  140)       !!      ==========  alone, for further diagnostics (for example to
(  141)       !!                  compute the turbocline depth in zdfmxl.F90).
(  142)       !!         avt = zavt
(  143)       !!         (avs = zavs if lk_zdfddm=T )
(  144)       !!
(  145)       !!      'key_trdtra' defined: trend saved for futher diagnostics.
(  146)       !!
(  147)       !!      macro-tasked on vertical slab (jj-loop)
(  148)       !!
(  149)       !! ** Action  : - Update (ta,sa) with before vertical diffusion trend
(  150)       !!              - Save the trend in (ztdta,ztdsa) ('key_trdtra')
(  151)       !!
(  152)       !! History :
(  153)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  154)       !!   7.0  !  91-11 (G. Madec)
(  155)       !!        !  92-06 (M. Imbard) correction on tracer trend loops
(  156)       !!        !  96-01 (G. Madec) statement function for e3
(  157)       !!        !  97-05 (G. Madec) vertical component of isopycnal
(  158)       !!        !  97-07 (G. Madec) geopotential diffusion in s-coord
(  159)       !!        !  00-08 (G. Madec) double diffusive mixing
(  160)       !!   8.5  !  02-08 (G. Madec)  F90: Free form and module
(  161)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(  162)       !!   9.0  !  05-06 (C. Ethe )  non-local flux in KPP vertical mixing scheme






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 4

(  163)       !!---------------------------------------------------------------------
(  164)       !! * Modules used
(  165)       USE oce    , ONLY :   zwd   => ua,  &  ! ua, va used as
(  166)                             zws   => va      ! workspace
(  167)       !! * Arguments
(  168)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  169) 
(  170)       !! * Local declarations
(  171)       INTEGER ::   ji, jj, jk                ! dummy loop indices
(  172)       REAL(wp) ::   &
(  173)          zavi, zrhs                          ! temporary scalars
(  174)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  175)          zwi, zwt, zavsi                     ! temporary workspace arrays
(  176) 
(  177) 
(  178) 
(  179)       ! I. Local constant initialization
(  180)       ! --------------------------------
(  181)       ! ... time step = 2 rdttra ex
(  182)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  183)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(  184)       ELSEIF( kt <= nit000 + 1) THEN
(  185)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  186)       ENDIF
(  187) 
(  188)       ! Save ta and sa trends
(  189)       IF( l_trdtra )   THEN
(  190)          tdta(:,:,:) = ta(:,:,:) 
(  191)          tdsa(:,:,:) = sa(:,:,:) 
(  192)       ENDIF
(  193) 
(  194)       zwd  (1,:, : ) = 0.e0     ;     zwd  (jpi,:,:) = 0.e0
(  195)       zws  (1,:, : ) = 0.e0     ;     zws  (jpi,:,:) = 0.e0
(  196)       zwi  (1,:, : ) = 0.e0     ;     zwi  (jpi,:,:) = 0.e0
(  197)       zwt  (1,:, : ) = 0.e0     ;     zwt  (jpi,:,:) = 0.e0
(  198)       zavsi(1,:, : ) = 0.e0     ;     zavsi(jpi,:,:) = 0.e0
(  199)       zwt  (:,:,jpk) = 0.e0     ;     zwt  ( : ,:,1) = 0.e0
(  200)       zavsi(:,:,jpk) = 0.e0     ;     zavsi( : ,:,1) = 0.e0
(  201) 
(  202)       ! II. Vertical trend associated with the vertical physics
(  203)       ! =======================================================
(  204)       !     (including the vertical flux proportional to dk[t] associated
(  205)       !      with the lateral mixing, through the avt update)
(  206)       !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  207) 
(  208) 
(  209)       ! II.0 Matrix construction
(  210)       ! ------------------------
(  211) 
(  212)       ! update and save of avt (and avs if double diffusive mixing)
(  213)       DO jk = 2, jpkm1
(  214)          DO jj = 2, jpjm1
(  215)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  216)                zavi = fsahtw(ji,jj,jk) * (                 &   ! vertical mixing coef. due to lateral mixing
(  217)                     wslpi(ji,jj,jk) * wslpi(ji,jj,jk)      &
(  218)                   + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  219)                zwt(ji,jj,jk) = avt(ji,jj,jk) + zavi            ! zwt=avt+zavi (total vertical mixing coef. on temperature)
(  220) #if defined key_zdfddm






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 5

(  221)                zavsi(ji,jj,jk) = fsavs(ji,jj,jk) + zavi        ! dd mixing: zavsi = total vertical mixing coef. on salinity
(  222) #endif
(  223)             END DO
(  224)          END DO
(  225)       END DO
(  226) 
(  227)       ! Diagonal, inferior, superior  (including the bottom boundary condition via avt masked)
(  228)       DO jk = 1, jpkm1
(  229)          DO jj = 2, jpjm1
(  230)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  231)                zwi(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  232)                zws(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  233)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
(  234)             END DO
(  235)          END DO
(  236)       END DO
(  237) 
(  238)       ! Surface boudary conditions
(  239)       DO jj = 2, jpjm1
(  240)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  241)             zwi(ji,jj,1) = 0.e0
(  242)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  243)          END DO
(  244)       END DO
(  245) 
(  246) 
(  247)       ! II.1. Vertical diffusion on t
(  248)       ! ---------------------------
(  249) 
(  250)       !! Matrix inversion from the first level
(  251)       !!----------------------------------------------------------------------
(  252)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  253)       !
(  254)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  255)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  256)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  257)       !        (        ...               )( ...  ) ( ...  )
(  258)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  259)       !
(  260)       !   m is decomposed in the product of an upper and lower triangular matrix
(  261)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  262)       !   The second member is in 2d array zwy
(  263)       !   The solution is in 2d array zwx
(  264)       !   The 3d arry zwt is a work space array
(  265)       !   zwy is used and then used as a work space array : its value is modified!
(  266) 
(  267)       ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
(  268)       DO jj = 2, jpjm1
(  269)          DO ji = fs_2, fs_jpim1
(  270)             zwt(ji,jj,1) = zwd(ji,jj,1)
(  271)          END DO
(  272)       END DO
(  273)       DO jk = 2, jpkm1
(  274)          DO jj = 2, jpjm1
(  275)             DO ji = fs_2, fs_jpim1
(  276)                zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
(  277)             END DO
(  278)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 6

(  279)       END DO
(  280) 
(  281)       ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
(  282) #if defined key_zdfkpp
(  283)          ! add non-local temperature flux ( in convective case only)
(  284)       DO jj = 2, jpjm1
(  285)          DO ji = fs_2, fs_jpim1
(  286)             ta(ji,jj,1) = tb(ji,jj,1) + r2dt(1) * ta(ji,jj,1) &
(  287)                   &  - r2dt(1) * ( ghats(ji,jj,1) * avt(ji,jj,1) - ghats(ji,jj,2) * avt(ji,jj,2) ) &
(  288)                   &               * wt0(ji,jj) / fse3t(ji,jj,2) 
(  289)          END DO
(  290)       END DO
(  291) 
(  292)       DO jk = 2, jpkm1
(  293)          DO jj = 2, jpjm1
(  294)             DO ji = fs_2, fs_jpim1
(  295)                ! zrhs=right hand side 
(  296)                zrhs = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)  &
(  297)                   &  - r2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
(  298)                   &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
(  299)                ta(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *ta(ji,jj,jk-1)
(  300)             END DO
(  301)          END DO
(  302)       END DO
(  303) #else
(  304)       DO jj = 2, jpjm1
(  305)          DO ji = fs_2, fs_jpim1
(  306)             ta(ji,jj,1) = tb(ji,jj,1) + r2dt(1) * ta(ji,jj,1)
(  307)          END DO
(  308)       END DO
(  309)       DO jk = 2, jpkm1
(  310)          DO jj = 2, jpjm1
(  311)             DO ji = fs_2, fs_jpim1
(  312)                zrhs = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)   ! zrhs=right hand side 
(  313)                ta(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *ta(ji,jj,jk-1)
(  314)             END DO
(  315)          END DO
(  316)       END DO
(  317) #endif
(  318) 
(  319)       ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
(  320)       ! Save the masked temperature after in ta
(  321)       ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
(  322)       DO jj = 2, jpjm1
(  323)          DO ji = fs_2, fs_jpim1
(  324)             ta(ji,jj,jpkm1) = ta(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
(  325)          END DO
(  326)       END DO
(  327)       DO jk = jpk-2, 1, -1
(  328)          DO jj = 2, jpjm1
(  329)             DO ji = fs_2, fs_jpim1
(  330)                ta(ji,jj,jk) = ( ta(ji,jj,jk) - zws(ji,jj,jk) * ta(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
(  331)             END DO
(  332)          END DO
(  333)       END DO
(  334) 
(  335)       ! II.2 Vertical diffusion on salinity
(  336)       ! -----------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 7

(  337) 
(  338) #if defined key_zdfddm
(  339)       ! Rebuild the Matrix as avt /= avs
(  340) 
(  341)       ! Diagonal, inferior, superior  (including the bottom boundary condition via avs masked)
(  342)       DO jk = 1, jpkm1
(  343)          DO jj = 2, jpjm1
(  344)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  345)                zwi(ji,jj,jk) = - r2dt(jk) * zavsi(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  346)                zws(ji,jj,jk) = - r2dt(jk) * zavsi(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  347)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
(  348)             END DO
(  349)          END DO
(  350)       END DO
(  351) 
(  352)       ! Surface boudary conditions
(  353)       DO jj = 2, jpjm1
(  354)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  355)             zwi(ji,jj,1) = 0.e0
(  356)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  357)          END DO
(  358)       END DO
(  359) #endif
(  360) 
(  361) 
(  362)       !! Matrix inversion from the first level
(  363)       !!----------------------------------------------------------------------
(  364)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  365)       !
(  366)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  367)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  368)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  369)       !        (        ...               )( ...  ) ( ...  )
(  370)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  371)       !
(  372)       !   m is decomposed in the product of an upper and lower triangular
(  373)       !   matrix
(  374)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  375)       !   The second member is in 2d array zwy
(  376)       !   The solution is in 2d array zwx
(  377)       !   The 3d arry zwt is a work space array
(  378)       !   zwy is used and then used as a work space array : its value is modified!
(  379) 
(  380)       ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
(  381)       DO jj = 2, jpjm1
(  382)          DO ji = fs_2, fs_jpim1
(  383)             zwt(ji,jj,1) = zwd(ji,jj,1)
(  384)          END DO
(  385)       END DO
(  386)       DO jk = 2, jpkm1
(  387)          DO jj = 2, jpjm1
(  388)             DO ji = fs_2, fs_jpim1
(  389)                zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
(  390)             END DO
(  391)          END DO
(  392)       END DO
(  393) 
(  394)       ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 8

(  395) 
(  396) #if defined key_zdfkpp
(  397)          ! add non-local temperature flux ( in convective case only)
(  398)       DO jj = 2, jpjm1
(  399)          DO ji = fs_2, fs_jpim1
(  400)             sa(ji,jj,1) = sb(ji,jj,1) + r2dt(1) * sa(ji,jj,1) &
(  401)                   &  - r2dt(1) * ( ghats(ji,jj,1) * fsavs(ji,jj,1) - ghats(ji,jj,2) * fsavs(ji,jj,2) ) &
(  402)                   &               * ws0(ji,jj) / fse3t(ji,jj,2) 
(  403)          END DO
(  404)       END DO
(  405) 
(  406)       DO jk = 2, jpkm1
(  407)          DO jj = 2, jpjm1
(  408)             DO ji = fs_2, fs_jpim1
(  409)                ! zrhs=right hand side 
(  410)                zrhs = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)  &
(  411)                   &  - r2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
(  412)                   &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
(  413)                sa(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *sa(ji,jj,jk-1)
(  414)             END DO
(  415)          END DO
(  416)       END DO
(  417) #else
(  418)       DO jj = 2, jpjm1
(  419)          DO ji = fs_2, fs_jpim1
(  420)             sa(ji,jj,1) = sb(ji,jj,1) + r2dt(1) * sa(ji,jj,1)
(  421)          END DO
(  422)       END DO
(  423)       DO jk = 2, jpkm1
(  424)          DO jj = 2, jpjm1
(  425)             DO ji = fs_2, fs_jpim1
(  426)                zrhs = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)   ! zrhs=right hand side
(  427)                sa(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *sa(ji,jj,jk-1)
(  428)             END DO
(  429)          END DO
(  430)       END DO
(  431) #endif
(  432) 
(  433)       ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
(  434)       ! Save the masked temperature after in ta
(  435)       ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
(  436)       DO jj = 2, jpjm1
(  437)          DO ji = fs_2, fs_jpim1
(  438)             sa(ji,jj,jpkm1) = sa(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
(  439)          END DO
(  440)       END DO
(  441)       DO jk = jpk-2, 1, -1
(  442)          DO jj = 2, jpjm1
(  443)             DO ji = fs_2, fs_jpim1
(  444)                sa(ji,jj,jk) = ( sa(ji,jj,jk) - zws(ji,jj,jk) * sa(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
(  445)             END DO
(  446)          END DO
(  447)       END DO
(  448) 
(  449) 
(  450)       ! save the trends for diagnostic
(  451)       ! compute the vertical diffusive trends in substracting the previous 
(  452)       ! trends tdta()/tdsa() to the new one computed via dT/dt or dS/dt 






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 9

(  453)       ! i.e. with the new temperature/salinity ta/sa computed above
(  454)       IF( l_trdtra )   THEN
(  455)          IF( ln_traldf_iso)   THEN
(  456)             DO jk = 1, jpkm1
(  457)                tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk) + tavg(:,:,jk) 
(  458)                tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk) + savg(:,:,jk) 
(  459)             END DO
(  460)          ELSE
(  461)             DO jk = 1, jpkm1
(  462)                tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk)                             
(  463)                tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk)                             
(  464)             END DO
(  465)          ENDIF
(  466) 
(  467)          CALL trd_mod(tdta, tdsa, jpttdzdf, 'TRA', kt)
(  468)       ENDIF
(  469) 
(  470)    END SUBROUTINE tra_zdf_zdf
(  471) 
(  472) 
(  473)    SUBROUTINE tra_zdf_iso( kt )
(  474)       !!----------------------------------------------------------------------
(  475)       !!                  ***  ROUTINE tra_zdf_iso  ***
(  476)       !!
(  477)       !! ** Purpose :
(  478)       !!     Compute the trend due to the vertical tracer diffusion inclu-
(  479)       !!     ding the vertical component of lateral mixing (only for second
(  480)       !!     order operator, for fourth order it is already computed and
(  481)       !!     add to the general trend in traldf.F) and add it to the general
(  482)       !!     trend of the tracer equations.
(  483)       !!
(  484)       !! ** Method :
(  485)       !!         The vertical component of the lateral diffusive trends is
(  486)       !!      provided by a 2nd order operator rotated along neural or geopo-
(  487)       !!      tential surfaces to which an eddy induced advection can be added
(  488)       !!      It is computed using before fields (forward in time) and isopyc-
(  489)       !!      nal or geopotential slopes computed in routine ldfslp.
(  490)       !!
(  491)       !!      First part: vertical trends associated with the lateral mixing
(  492)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(  493)       !!      vertical fluxes associated with the rotated lateral mixing:
(  494)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
(  495)       !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
(  496)       !!      save avt coef. resulting from vertical physics alone in zavt:
(  497)       !!         zavt = avt
(  498)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(  499)       !!         avt = avt + wslpi^2+wslj^2
(  500)       !!      add vertical Eddy Induced advective fluxes (lk_traldf_eiv=T):
(  501)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(  502)       !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
(  503)       !!      take the horizontal divergence of the fluxes:
(  504)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(  505)       !!      Add this trend to the general trend (ta,sa):
(  506)       !!         ta = ta + difft
(  507)       !!
(  508)       !! ** Action :
(  509)       !!         Update (ta,sa) arrays with the before vertical diffusion trend
(  510)       !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 10

(  511)       !!
(  512)       !! History :
(  513)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  514)       !!   7.0  !  91-11  (G. Madec)
(  515)       !!        !  92-06  (M. Imbard) correction on tracer trend loops
(  516)       !!        !  96-01  (G. Madec) statement function for e3
(  517)       !!        !  97-05  (G. Madec) vertical component of isopycnal
(  518)       !!        !  97-07  (G. Madec) geopotential diffusion in s-coord
(  519)       !!        !  00-08  (G. Madec) double diffusive mixing
(  520)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  521)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  522)       !!---------------------------------------------------------------------
(  523)       !! * Modules used
(  524)       USE oce    , ONLY :   zwx => ua,  &  ! use ua, va as
(  525)                             zwy => va      ! workspace arrays
(  526) 
(  527)       !! * Arguments
(  528)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  529) 
(  530)       !! * Local declarations
(  531)       INTEGER ::   ji, jj, jk       ! dummy loop indices
(  532) #if defined key_partial_steps
(  533)       INTEGER ::   iku, ikv
(  534) #endif
(  535)       REAL(wp) ::   &
(  536)          zcoef0, zcoef3,         &  !    "         "
(  537)          zcoef4,                 &  !    "         "
(  538)          zbtr, zmku, zmkv,       &  !    "         "
(  539) #if defined key_traldf_eiv
(  540)          zcoeg3,                 &  !    "         "
(  541)          zuwki, zvwki,           &  !    "         "
(  542)          zuwk, zvwk,             &  !    "         "
(  543) #endif
(  544)          ztav, zsav
(  545)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  546)          zwz, zwt, ztfw, zsfw       ! temporary workspace arrays
(  547) 
(  548) 
(  549)       ! 0. Local constant initialization
(  550)       ! --------------------------------
(  551)       zwx (1,:,:) = 0.e0     ;     zwx (jpi,:,:) = 0.e0
(  552)       zwy (1,:,:) = 0.e0     ;     zwy (jpi,:,:) = 0.e0
(  553)       zwz (1,:,:) = 0.e0     ;     zwz (jpi,:,:) = 0.e0
(  554)       zwt (1,:,:) = 0.e0     ;     zwt (jpi,:,:) = 0.e0
(  555)       ztfw(1,:,:) = 0.e0     ;     ztfw(jpi,:,:) = 0.e0
(  556)       zsfw(1,:,:) = 0.e0     ;     zsfw(jpi,:,:) = 0.e0
(  557) 
(  558)       ! Save ta and sa trends
(  559)       IF( l_trdtra )   THEN
(  560)          tdta(:,:,:) = ta(:,:,:) 
(  561)          tdsa(:,:,:) = sa(:,:,:) 
(  562)       ENDIF
(  563) 
(  564)       ! I. Vertical trends associated with lateral mixing
(  565)       ! -------------------------------------------------
(  566)       !    (excluding the vertical flux proportional to dk[t] )
(  567) 
(  568) 






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 11

(  569)       ! I.1 horizontal tracer gradient
(  570)       ! ------------------------------
(  571) 
(  572)       DO jk = 1, jpkm1
(  573)          DO jj = 1, jpjm1
(  574)             DO ji = 1, fs_jpim1   ! vector opt.
(  575)                ! i-gradient of T and S at jj
(  576)                zwx (ji,jj,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
(  577)                zwy (ji,jj,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
(  578)                ! j-gradient of T and S at jj
(  579)                zwz (ji,jj,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  580)                zwt (ji,jj,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  581)             END DO
(  582)          END DO
(  583)       END DO
(  584) #  if defined key_partial_steps
(  585)       ! partial steps correction at the bottom ocean level 
(  586)       DO jj = 1, jpjm1
(  587)          DO ji = 1, fs_jpim1   ! vector opt.
(  588)             ! last ocean level
(  589)             iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  590)             ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  591)             ! i-gradient of T and S
(  592)             zwx (ji,jj,iku) = gtu(ji,jj)
(  593)             zwy (ji,jj,iku) = gsu(ji,jj)
(  594)             ! j-gradient of T and S
(  595)             zwz (ji,jj,ikv) = gtv(ji,jj) 
(  596)             zwt (ji,jj,ikv) = gsv(ji,jj) 
(  597)          END DO
(  598)          END DO
(  599) #endif
(  600) 
(  601) 
(  602)       ! I.2 Vertical fluxes
(  603)       ! -------------------
(  604) 
(  605)       ! Surface and bottom vertical fluxes set to zero
(  606)       ztfw(:,:, 1 ) = 0.e0
(  607)       zsfw(:,:, 1 ) = 0.e0
(  608)       ztfw(:,:,jpk) = 0.e0
(  609)       zsfw(:,:,jpk) = 0.e0
(  610) 
(  611)       ! interior (2=<jk=<jpk-1)
(  612)       DO jk = 2, jpkm1
(  613)          DO jj = 2, jpjm1
(  614)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  615)                zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  616) 
(  617)                zmku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
(  618)                               + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )
(  619) 
(  620)                zmkv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
(  621)                               + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )
(  622) 
(  623)                zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  624)                zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  625) 
(  626)                ztfw(ji,jj,jk) = zcoef3 * (   zwx(ji  ,jj  ,jk-1) + zwx(ji-1,jj  ,jk)      &






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 12

(  627)                                            + zwx(ji-1,jj  ,jk-1) + zwx(ji  ,jj  ,jk)  )   &
(  628)                               + zcoef4 * (   zwz(ji  ,jj  ,jk-1) + zwz(ji  ,jj-1,jk)      &
(  629)                                            + zwz(ji  ,jj-1,jk-1) + zwz(ji  ,jj  ,jk)  )
(  630) 
(  631)                zsfw(ji,jj,jk) = zcoef3 * (   zwy(ji  ,jj  ,jk-1) + zwy(ji-1,jj  ,jk)      &
(  632)                                            + zwy(ji-1,jj  ,jk-1) + zwy(ji  ,jj  ,jk)  )   &
(  633)                               + zcoef4 * (   zwt(ji  ,jj  ,jk-1) + zwt(ji  ,jj-1,jk)      &
(  634)                                            + zwt(ji  ,jj-1,jk-1) + zwt(ji  ,jj  ,jk)  )
(  635)             END DO
(  636)          END DO
(  637)       END DO
(  638) 
(  639) #if defined key_traldf_eiv
(  640)       !                              ! ---------------------------------------!
(  641)       !                              ! Eddy induced vertical advective fluxes !
(  642)       !                              ! ---------------------------------------!
(  643)          zwx(:,:, 1 ) = 0.e0
(  644)          zwy(:,:, 1 ) = 0.e0
(  645)          zwx(:,:,jpk) = 0.e0
(  646)          zwy(:,:,jpk) = 0.e0
(  647) 
(  648)          DO jk = 2, jpkm1
(  649)             DO jj = 2, jpjm1
(  650)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  651) #   if defined key_traldf_c2d || defined key_traldf_c3d
(  652)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  653)                      &  * fsaeiu(ji-1,jj,jk) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
(  654)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  655)                      &  * fsaeiu(ji  ,jj,jk) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
(  656)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  657)                      &  * fsaeiv(ji,jj-1,jk) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
(  658)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  659)                      &  * fsaeiv(ji,jj  ,jk) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
(  660) 
(  661)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
(  662) #   else
(  663)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  664)                      &  * e2u(ji-1,jj) * umask(ji-1,jj,jk)
(  665)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  666)                      &  * e2u(ji  ,jj) * umask(ji  ,jj,jk)
(  667)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  668)                      &  * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
(  669)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  670)                      &  * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
(  671) 
(  672)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
(  673)                      &            * ( zuwk - zuwki + zvwk - zvwki )
(  674) #   endif
(  675)                   zwx(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
(  676)                   zwy(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )
(  677) 
(  678)                   ztfw(ji,jj,jk) = ztfw(ji,jj,jk) + zwx(ji,jj,jk)
(  679)                   zsfw(ji,jj,jk) = zsfw(ji,jj,jk) + zwy(ji,jj,jk)
(  680) #   if defined key_diaeiv
(  681)                   w_eiv(ji,jj,jk) = -2. * zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  682) #   endif
(  683)                END DO
(  684)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 13

(  685)          END DO
(  686) #endif
(  687) 
(  688)       ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  689)       ! -------------------------------------------------------------------
(  690) 
(  691)       DO jk = 1, jpkm1
(  692)          DO jj = 2, jpjm1
(  693)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  694)                zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  695)                ztav = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
(  696)                zsav = (  zsfw(ji,jj,jk) - zsfw(ji,jj,jk+1)  ) * zbtr
(  697)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
(  698)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
(  699)             END DO
(  700)          END DO
(  701)       END DO
(  702) 
(  703)       ! save the trends for diagnostic
(  704)       !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
(  705)       IF( l_trdtra )   THEN
(  706) #   if defined key_traldf_eiv
(  707)          ! Compute the vertical Gent velocity trend
(  708)          DO jk = 1, jpkm1
(  709)             DO jj = 2, jpjm1
(  710)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  711)                   zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  712)                   tavg(ji,jj,jk) = ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) ) * zbtr
(  713)                   savg(ji,jj,jk) = ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) ) * zbtr
(  714)                END DO
(  715)             END DO
(  716)          END DO
(  717) 
(  718)          CALL trd_mod(tavg, savg, jpttddoe, 'TRA', kt)
(  719) #   endif
(  720)          ! Recompute the divergence of vertical fluxes ztav & zsav trends 
(  721)          ! computed at step 1.5 above in making the difference between the new 
(  722)          ! trend ta()/sa() and the previous one tdta()/tdsa() and substract 
(  723)          ! the vertical Gent velocity trend tavg()/savg() (zero if not used)
(  724)          tavg(:,:,:) = ta(:,:,:) - tdta(:,:,:) - tavg(:,:,:) 
(  725)          savg(:,:,:) = sa(:,:,:) - tdsa(:,:,:) - savg(:,:,:) 
(  726)       ENDIF
(  727) 
(  728)    END SUBROUTINE tra_zdf_iso
(  729) 
(  730) #else
(  731)    !!----------------------------------------------------------------------
(  732)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  733)    !!----------------------------------------------------------------------
(  734) CONTAINS
(  735)    SUBROUTINE tra_zdf_iso_vopt( kt )              ! empty routine
(  736) !      WRITE(*,*) 'tra_zdf_iso_vopt: You should not have seen this print! error?', kt
(  737)    END SUBROUTINE tra_zdf_iso_vopt
(  738) #endif
(  739) 
(  740)    !!==============================================================================
(  741) END MODULE trazdf_iso_vopt







PGF90 (Version     10.5)          09/23/2010  09:50:18      page 14

(    1) # 1 "trazdf_iso_vopt.F90"
(    1) MODULE trazdf_iso_vopt
(    2)    !!==============================================================================
(    3)    !!                 ***  MODULE  trazdf_iso_vopt  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_zdf_iso_vopt : Update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion. vector optimization, use
(   13)    !!                  k-j-i loops.
(   14)    !!   tra_zdf_iso  :
(   15)    !!   tra_zdf_zdf  :
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce             ! ocean dynamics and tracers variables
(   19)    USE dom_oce         ! ocean space and time domain variables 
(   20)    USE zdf_oce         ! ocean vertical physics variables
(   21)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   22)    USE trdmod          ! ocean active tracers trends 
(   23)    USE trdmod_oce      ! ocean variables trends
(   24)    USE ldfslp          ! iso-neutral slopes 
(   25)    USE zdfddm          ! ocean vertical physics: double diffusion
(   26)    USE in_out_manager  ! I/O manager
(   27)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   28)    USE zdfkpp          ! KPP parameterisation
(   29)    USE prtctl          ! Print control
(   30) 
(   31)    IMPLICIT NONE
(   32)    PRIVATE
(   33) 
(   34)    !! * Routine accessibility
(   35)    PUBLIC tra_zdf_iso_vopt   !  routine called by step.F90
(   36) 
(   37)    !! * Module variables
(   38)    REAL(wp), DIMENSION(jpk) ::  &
(   39)       r2dt                          ! vertical profile of 2 x time-step
(   40)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
(   41)       tavg, savg,                     & ! workspace arrays
(   42)       tdta, tdsa                        ! workspace arrays
(   43) 
(   44)    !! * Substitutions
(   45) # 1 "./domzgr_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 15

(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 46 "trazdf_iso_vopt.F90"
(   46) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19
(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 47 "trazdf_iso_vopt.F90"
(   47) # 1 "./ldfeiv_substitute.h90"
(    1) # 40 "./ldfeiv_substitute.h90"
(   40) !   Default option :                         eiv: dummy variables
(   41) # 48 "trazdf_iso_vopt.F90"
(   48) # 1 "./zdfddm_substitute.h90"
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
(   17) # 49 "trazdf_iso_vopt.F90"
(   49) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 50 "trazdf_iso_vopt.F90"
(   50) # 50
(   50)    !!----------------------------------------------------------------------
(   51)    !!----------------------------------------------------------------------
(   52)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   53)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso_vopt.F90,v 1.10 2005/09/22 10:29:06 opalod Exp $ 
(   54)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   55)    !!----------------------------------------------------------------------
(   56) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 16

(   57)    
(   58)    SUBROUTINE tra_zdf_iso_vopt( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                  ***  ROUTINE tra_zdf_iso_vopt  ***
(   61)       !!
(   62)       !! ** Purpose :
(   63)       !! ** Method  :
(   64)       !! ** Action  :
(   65)       !!
(   66)       !! History :
(   67)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   68)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   69)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   70)       !!---------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   73) 
(   74)       IF( kt == nit000 ) THEN
(   75)          IF(lwp)WRITE(numout,*)
(   76)          IF(lwp)WRITE(numout,*) 'tra_zdf_iso_vopt : vertical mixing computation'
(   77)          IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~~   is  iso-neutral diffusion : implicit vertical time stepping'
(   78) # 81
(   81)       ENDIF
(   82) 
(   83)       ! initialization step
(   84)       tavg(:,:,:) = 0.e0
(   85)       savg(:,:,:) = 0.e0
(   86) 
(   87)       ! I. vertical extra-diagonal part of the rotated tensor
(   88)       ! -----------------------------------------------------
(   89) 
(   90)       CALL tra_zdf_iso( kt )
(   91) 
(   92)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(   93)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
(   94)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(   95)       ENDIF
(   96) 
(   97)       ! II. vertical diffusion (including the vertical diagonal part of the rotated tensor)
(   98)       ! ----------------------
(   99) 
(  100)       CALL tra_zdf_zdf( kt )
(  101) 
(  102)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  103)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
(  104)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)
(  105)       ENDIF
(  106) 
(  107)    END SUBROUTINE tra_zdf_iso_vopt
(  108) 
(  109) 
(  110)    SUBROUTINE tra_zdf_zdf( kt )
(  111)       !!----------------------------------------------------------------------
(  112)       !!                 ***  ROUTINE tra_zdf_zdf  ***
(  113)       !!                    
(  114)       !! ** Purpose :   Compute the trend due to the vertical tracer diffusion
(  115)       !!     including the vertical component of lateral mixing (only for 2nd
(  116)       !!     order operator, for fourth order it is already computed and add






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 17

(  117)       !!     to the general trend in traldf.F) and add it to the general trend
(  118)       !!     of the tracer equations.
(  119)       !!
(  120)       !! ** Method  :   The vertical component of the lateral diffusive trends
(  121)       !!      is provided by a 2nd order operator rotated along neural or geo-
(  122)       !!      potential surfaces to which an eddy induced advection can be 
(  123)       !!      added. It is computed using before fields (forward in time) and 
(  124)       !!      isopycnal or geopotential slopes computed in routine ldfslp.
(  125)       !!
(  126)       !!      Second part: vertical trend associated with the vertical physics
(  127)       !!      ===========  (including the vertical flux proportional to dk[t]
(  128)       !!                  associated with the lateral mixing, through the
(  129)       !!                  update of avt)
(  130)       !!      The vertical diffusion of tracers (t & s) is given by:
(  131)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(  132)       !!      It is computed using a backward time scheme (t=ta).
(  133)       !!      Surface and bottom boundary conditions: no diffusive flux on
(  134)       !!      both tracers (bottom, applied through the masked field avt).
(  135)       !!      Add this trend to the general trend ta,sa :
(  136)       !!         ta = ta + dz( avt dz(t) )
(  137)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
(  138)       !!
(  139)       !!      Third part: recover avt resulting from the vertical physics
(  140)       !!      ==========  alone, for further diagnostics (for example to
(  141)       !!                  compute the turbocline depth in zdfmxl.F90).
(  142)       !!         avt = zavt
(  143)       !!         (avs = zavs if lk_zdfddm=T )
(  144)       !!
(  145)       !!      'key_trdtra' defined: trend saved for futher diagnostics.
(  146)       !!
(  147)       !!      macro-tasked on vertical slab (jj-loop)
(  148)       !!
(  149)       !! ** Action  : - Update (ta,sa) with before vertical diffusion trend
(  150)       !!              - Save the trend in (ztdta,ztdsa) ('key_trdtra')
(  151)       !!
(  152)       !! History :
(  153)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  154)       !!   7.0  !  91-11 (G. Madec)
(  155)       !!        !  92-06 (M. Imbard) correction on tracer trend loops
(  156)       !!        !  96-01 (G. Madec) statement function for e3
(  157)       !!        !  97-05 (G. Madec) vertical component of isopycnal
(  158)       !!        !  97-07 (G. Madec) geopotential diffusion in s-coord
(  159)       !!        !  00-08 (G. Madec) double diffusive mixing
(  160)       !!   8.5  !  02-08 (G. Madec)  F90: Free form and module
(  161)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(  162)       !!   9.0  !  05-06 (C. Ethe )  non-local flux in KPP vertical mixing scheme
(  163)       !!---------------------------------------------------------------------
(  164)       !! * Modules used
(  165)       USE oce    , ONLY :   zwd   => ua,  &  ! ua, va used as
(  166)                             zws   => va      ! workspace
(  167)       !! * Arguments
(  168)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  169) 
(  170)       !! * Local declarations
(  171)       INTEGER ::   ji, jj, jk                ! dummy loop indices
(  172)       REAL(wp) ::   &
(  173)          zavi, zrhs                          ! temporary scalars
(  174)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 18

(  175)          zwi, zwt, zavsi                     ! temporary workspace arrays
(  176) 
(  177) 
(  178) 
(  179)       ! I. Local constant initialization
(  180)       ! --------------------------------
(  181)       ! ... time step = 2 rdttra ex
(  182)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(  183)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(  184)       ELSEIF( kt <= nit000 + 1) THEN
(  185)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  186)       ENDIF
(  187) 
(  188)       ! Save ta and sa trends
(  189)       IF( l_trdtra )   THEN
(  190)          tdta(:,:,:) = ta(:,:,:) 
(  191)          tdsa(:,:,:) = sa(:,:,:) 
(  192)       ENDIF
(  193) 
(  194)       zwd  (1,:, : ) = 0.e0     ;     zwd  (jpi,:,:) = 0.e0
(  195)       zws  (1,:, : ) = 0.e0     ;     zws  (jpi,:,:) = 0.e0
(  196)       zwi  (1,:, : ) = 0.e0     ;     zwi  (jpi,:,:) = 0.e0
(  197)       zwt  (1,:, : ) = 0.e0     ;     zwt  (jpi,:,:) = 0.e0
(  198)       zavsi(1,:, : ) = 0.e0     ;     zavsi(jpi,:,:) = 0.e0
(  199)       zwt  (:,:,jpk) = 0.e0     ;     zwt  ( : ,:,1) = 0.e0
(  200)       zavsi(:,:,jpk) = 0.e0     ;     zavsi( : ,:,1) = 0.e0
(  201) 
(  202)       ! II. Vertical trend associated with the vertical physics
(  203)       ! =======================================================
(  204)       !     (including the vertical flux proportional to dk[t] associated
(  205)       !      with the lateral mixing, through the avt update)
(  206)       !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  207) 
(  208) 
(  209)       ! II.0 Matrix construction
(  210)       ! ------------------------
(  211) 
(  212)       ! update and save of avt (and avs if double diffusive mixing)
(  213)       DO jk = 2, jpkm1
(  214)          DO jj = 2, jpjm1
(  215)             DO ji = 2, jpim1   ! vector opt.
(  216)                zavi = ahtw(ji,jj,jk) * (                 &   ! vertical mixing coef. due to lateral mixing
(  217)                     wslpi(ji,jj,jk) * wslpi(ji,jj,jk)      &
(  218)                   + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  219)                zwt(ji,jj,jk) = avt(ji,jj,jk) + zavi            ! zwt=avt+zavi (total vertical mixing coef. on temperature)
(  220) # 223
(  223)             END DO
(  224)          END DO
(  225)       END DO
(  226) 
(  227)       ! Diagonal, inferior, superior  (including the bottom boundary condition via avt masked)
(  228)       DO jk = 1, jpkm1
(  229)          DO jj = 2, jpjm1
(  230)             DO ji = 2, jpim1   ! vector opt.
(  231)                zwi(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk  ) / ( e3t(jk) * e3w(jk  ) )
(  232)                zws(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk+1) / ( e3t(jk) * e3w(jk+1) )
(  233)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
(  234)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 19

(  235)          END DO
(  236)       END DO
(  237) 
(  238)       ! Surface boudary conditions
(  239)       DO jj = 2, jpjm1
(  240)          DO ji = 2, jpim1   ! vector opt.
(  241)             zwi(ji,jj,1) = 0.e0
(  242)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  243)          END DO
(  244)       END DO
(  245) 
(  246) 
(  247)       ! II.1. Vertical diffusion on t
(  248)       ! ---------------------------
(  249) 
(  250)       !! Matrix inversion from the first level
(  251)       !!----------------------------------------------------------------------
(  252)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  253)       !
(  254)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  255)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  256)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  257)       !        (        ...               )( ...  ) ( ...  )
(  258)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  259)       !
(  260)       !   m is decomposed in the product of an upper and lower triangular matrix
(  261)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  262)       !   The second member is in 2d array zwy
(  263)       !   The solution is in 2d array zwx
(  264)       !   The 3d arry zwt is a work space array
(  265)       !   zwy is used and then used as a work space array : its value is modified!
(  266) 
(  267)       ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
(  268)       DO jj = 2, jpjm1
(  269)          DO ji = 2, jpim1
(  270)             zwt(ji,jj,1) = zwd(ji,jj,1)
(  271)          END DO
(  272)       END DO
(  273)       DO jk = 2, jpkm1
(  274)          DO jj = 2, jpjm1
(  275)             DO ji = 2, jpim1
(  276)                zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
(  277)             END DO
(  278)          END DO
(  279)       END DO
(  280) 
(  281)       ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
(  282) # 304
(  304)       DO jj = 2, jpjm1
(  305)          DO ji = 2, jpim1
(  306)             ta(ji,jj,1) = tb(ji,jj,1) + r2dt(1) * ta(ji,jj,1)
(  307)          END DO
(  308)       END DO
(  309)       DO jk = 2, jpkm1
(  310)          DO jj = 2, jpjm1
(  311)             DO ji = 2, jpim1
(  312)                zrhs = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)   ! zrhs=right hand side 
(  313)                ta(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *ta(ji,jj,jk-1)






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 20

(  314)             END DO
(  315)          END DO
(  316)       END DO
(  317) 
(  318) # 319
(  319)       ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
(  320)       ! Save the masked temperature after in ta
(  321)       ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
(  322)       DO jj = 2, jpjm1
(  323)          DO ji = 2, jpim1
(  324)             ta(ji,jj,jpkm1) = ta(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
(  325)          END DO
(  326)       END DO
(  327)       DO jk = jpk-2, 1, -1
(  328)          DO jj = 2, jpjm1
(  329)             DO ji = 2, jpim1
(  330)                ta(ji,jj,jk) = ( ta(ji,jj,jk) - zws(ji,jj,jk) * ta(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
(  331)             END DO
(  332)          END DO
(  333)       END DO
(  334) 
(  335)       ! II.2 Vertical diffusion on salinity
(  336)       ! -----------------------------------
(  337) 
(  338) 
(  339) 
(  340) # 362
(  362)       !! Matrix inversion from the first level
(  363)       !!----------------------------------------------------------------------
(  364)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  365)       !
(  366)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  367)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  368)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  369)       !        (        ...               )( ...  ) ( ...  )
(  370)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  371)       !
(  372)       !   m is decomposed in the product of an upper and lower triangular
(  373)       !   matrix
(  374)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  375)       !   The second member is in 2d array zwy
(  376)       !   The solution is in 2d array zwx
(  377)       !   The 3d arry zwt is a work space array
(  378)       !   zwy is used and then used as a work space array : its value is modified!
(  379) 
(  380)       ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
(  381)       DO jj = 2, jpjm1
(  382)          DO ji = 2, jpim1
(  383)             zwt(ji,jj,1) = zwd(ji,jj,1)
(  384)          END DO
(  385)       END DO
(  386)       DO jk = 2, jpkm1
(  387)          DO jj = 2, jpjm1
(  388)             DO ji = 2, jpim1
(  389)                zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
(  390)             END DO
(  391)          END DO
(  392)       END DO






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 21

(  393) 
(  394)       ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
(  395) 
(  396) # 418
(  418)       DO jj = 2, jpjm1
(  419)          DO ji = 2, jpim1
(  420)             sa(ji,jj,1) = sb(ji,jj,1) + r2dt(1) * sa(ji,jj,1)
(  421)          END DO
(  422)       END DO
(  423)       DO jk = 2, jpkm1
(  424)          DO jj = 2, jpjm1
(  425)             DO ji = 2, jpim1
(  426)                zrhs = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)   ! zrhs=right hand side
(  427)                sa(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *sa(ji,jj,jk-1)
(  428)             END DO
(  429)          END DO
(  430)       END DO
(  431) 
(  432) # 433
(  433)       ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
(  434)       ! Save the masked temperature after in ta
(  435)       ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
(  436)       DO jj = 2, jpjm1
(  437)          DO ji = 2, jpim1
(  438)             sa(ji,jj,jpkm1) = sa(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
(  439)          END DO
(  440)       END DO
(  441)       DO jk = jpk-2, 1, -1
(  442)          DO jj = 2, jpjm1
(  443)             DO ji = 2, jpim1
(  444)                sa(ji,jj,jk) = ( sa(ji,jj,jk) - zws(ji,jj,jk) * sa(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
(  445)             END DO
(  446)          END DO
(  447)       END DO
(  448) 
(  449) 
(  450)       ! save the trends for diagnostic
(  451)       ! compute the vertical diffusive trends in substracting the previous 
(  452)       ! trends tdta()/tdsa() to the new one computed via dT/dt or dS/dt 
(  453)       ! i.e. with the new temperature/salinity ta/sa computed above
(  454)       IF( l_trdtra )   THEN
(  455)          IF( ln_traldf_iso)   THEN
(  456)             DO jk = 1, jpkm1
(  457)                tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk) + tavg(:,:,jk) 
(  458)                tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk) + savg(:,:,jk) 
(  459)             END DO
(  460)          ELSE
(  461)             DO jk = 1, jpkm1
(  462)                tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk)                             
(  463)                tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk)                             
(  464)             END DO
(  465)          ENDIF
(  466) 
(  467)          CALL trd_mod(tdta, tdsa, jpttdzdf, 'TRA', kt)
(  468)       ENDIF
(  469) 
(  470)    END SUBROUTINE tra_zdf_zdf
(  471) 






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 22

(  472) 
(  473)    SUBROUTINE tra_zdf_iso( kt )
(  474)       !!----------------------------------------------------------------------
(  475)       !!                  ***  ROUTINE tra_zdf_iso  ***
(  476)       !!
(  477)       !! ** Purpose :
(  478)       !!     Compute the trend due to the vertical tracer diffusion inclu-
(  479)       !!     ding the vertical component of lateral mixing (only for second
(  480)       !!     order operator, for fourth order it is already computed and
(  481)       !!     add to the general trend in traldf.F) and add it to the general
(  482)       !!     trend of the tracer equations.
(  483)       !!
(  484)       !! ** Method :
(  485)       !!         The vertical component of the lateral diffusive trends is
(  486)       !!      provided by a 2nd order operator rotated along neural or geopo-
(  487)       !!      tential surfaces to which an eddy induced advection can be added
(  488)       !!      It is computed using before fields (forward in time) and isopyc-
(  489)       !!      nal or geopotential slopes computed in routine ldfslp.
(  490)       !!
(  491)       !!      First part: vertical trends associated with the lateral mixing
(  492)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(  493)       !!      vertical fluxes associated with the rotated lateral mixing:
(  494)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
(  495)       !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
(  496)       !!      save avt coef. resulting from vertical physics alone in zavt:
(  497)       !!         zavt = avt
(  498)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(  499)       !!         avt = avt + wslpi^2+wslj^2
(  500)       !!      add vertical Eddy Induced advective fluxes (lk_traldf_eiv=T):
(  501)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(  502)       !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
(  503)       !!      take the horizontal divergence of the fluxes:
(  504)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(  505)       !!      Add this trend to the general trend (ta,sa):
(  506)       !!         ta = ta + difft
(  507)       !!
(  508)       !! ** Action :
(  509)       !!         Update (ta,sa) arrays with the before vertical diffusion trend
(  510)       !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
(  511)       !!
(  512)       !! History :
(  513)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  514)       !!   7.0  !  91-11  (G. Madec)
(  515)       !!        !  92-06  (M. Imbard) correction on tracer trend loops
(  516)       !!        !  96-01  (G. Madec) statement function for e3
(  517)       !!        !  97-05  (G. Madec) vertical component of isopycnal
(  518)       !!        !  97-07  (G. Madec) geopotential diffusion in s-coord
(  519)       !!        !  00-08  (G. Madec) double diffusive mixing
(  520)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  521)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(  522)       !!---------------------------------------------------------------------
(  523)       !! * Modules used
(  524)       USE oce    , ONLY :   zwx => ua,  &  ! use ua, va as
(  525)                             zwy => va      ! workspace arrays
(  526) 
(  527)       !! * Arguments
(  528)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  529) 






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 23

(  530)       !! * Local declarations
(  531)       INTEGER ::   ji, jj, jk       ! dummy loop indices
(  532) # 535
(  535)       REAL(wp) ::   &
(  536)          zcoef0, zcoef3,         &  !    "         "
(  537)          zcoef4,                 &  !    "         "
(  538)          zbtr, zmku, zmkv,       &  !    "         "
(  544)          ztav, zsav
(  545)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  546)          zwz, zwt, ztfw, zsfw       ! temporary workspace arrays
(  547) 
(  548) 
(  549)       ! 0. Local constant initialization
(  550)       ! --------------------------------
(  551)       zwx (1,:,:) = 0.e0     ;     zwx (jpi,:,:) = 0.e0
(  552)       zwy (1,:,:) = 0.e0     ;     zwy (jpi,:,:) = 0.e0
(  553)       zwz (1,:,:) = 0.e0     ;     zwz (jpi,:,:) = 0.e0
(  554)       zwt (1,:,:) = 0.e0     ;     zwt (jpi,:,:) = 0.e0
(  555)       ztfw(1,:,:) = 0.e0     ;     ztfw(jpi,:,:) = 0.e0
(  556)       zsfw(1,:,:) = 0.e0     ;     zsfw(jpi,:,:) = 0.e0
(  557) 
(  558)       ! Save ta and sa trends
(  559)       IF( l_trdtra )   THEN
(  560)          tdta(:,:,:) = ta(:,:,:) 
(  561)          tdsa(:,:,:) = sa(:,:,:) 
(  562)       ENDIF
(  563) 
(  564)       ! I. Vertical trends associated with lateral mixing
(  565)       ! -------------------------------------------------
(  566)       !    (excluding the vertical flux proportional to dk[t] )
(  567) 
(  568) 
(  569)       ! I.1 horizontal tracer gradient
(  570)       ! ------------------------------
(  571) 
(  572)       DO jk = 1, jpkm1
(  573)          DO jj = 1, jpjm1
(  574)             DO ji = 1, jpim1   ! vector opt.
(  575)                ! i-gradient of T and S at jj
(  576)                zwx (ji,jj,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
(  577)                zwy (ji,jj,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
(  578)                ! j-gradient of T and S at jj
(  579)                zwz (ji,jj,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  580)                zwt (ji,jj,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
(  581)             END DO
(  582)          END DO
(  583)       END DO
(  584) 
(  585) 
(  586) # 602
(  602)       ! I.2 Vertical fluxes
(  603)       ! -------------------
(  604) 
(  605)       ! Surface and bottom vertical fluxes set to zero
(  606)       ztfw(:,:, 1 ) = 0.e0
(  607)       zsfw(:,:, 1 ) = 0.e0
(  608)       ztfw(:,:,jpk) = 0.e0
(  609)       zsfw(:,:,jpk) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 24

(  610) 
(  611)       ! interior (2=<jk=<jpk-1)
(  612)       DO jk = 2, jpkm1
(  613)          DO jj = 2, jpjm1
(  614)             DO ji = 2, jpim1   ! vector opt.
(  615)                zcoef0 = - ahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  616) 
(  617)                zmku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
(  618)                               + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )
(  619) 
(  620)                zmkv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
(  621)                               + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )
(  622) 
(  623)                zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  624)                zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  625) 
(  626)                ztfw(ji,jj,jk) = zcoef3 * (   zwx(ji  ,jj  ,jk-1) + zwx(ji-1,jj  ,jk)      &
(  627)                                            + zwx(ji-1,jj  ,jk-1) + zwx(ji  ,jj  ,jk)  )   &
(  628)                               + zcoef4 * (   zwz(ji  ,jj  ,jk-1) + zwz(ji  ,jj-1,jk)      &
(  629)                                            + zwz(ji  ,jj-1,jk-1) + zwz(ji  ,jj  ,jk)  )
(  630) 
(  631)                zsfw(ji,jj,jk) = zcoef3 * (   zwy(ji  ,jj  ,jk-1) + zwy(ji-1,jj  ,jk)      &
(  632)                                            + zwy(ji-1,jj  ,jk-1) + zwy(ji  ,jj  ,jk)  )   &
(  633)                               + zcoef4 * (   zwt(ji  ,jj  ,jk-1) + zwt(ji  ,jj-1,jk)      &
(  634)                                            + zwt(ji  ,jj-1,jk-1) + zwt(ji  ,jj  ,jk)  )
(  635)             END DO
(  636)          END DO
(  637)       END DO
(  638) 
(  639) 
(  640) # 688
(  688)       ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  689)       ! -------------------------------------------------------------------
(  690) 
(  691)       DO jk = 1, jpkm1
(  692)          DO jj = 2, jpjm1
(  693)             DO ji = 2, jpim1   ! vector opt.
(  694)                zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*e3t(jk) )
(  695)                ztav = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
(  696)                zsav = (  zsfw(ji,jj,jk) - zsfw(ji,jj,jk+1)  ) * zbtr
(  697)                ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
(  698)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
(  699)             END DO
(  700)          END DO
(  701)       END DO
(  702) 
(  703)       ! save the trends for diagnostic
(  704)       !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
(  705)       IF( l_trdtra )   THEN
(  706) # 720
(  720)          ! Recompute the divergence of vertical fluxes ztav & zsav trends 
(  721)          ! computed at step 1.5 above in making the difference between the new 
(  722)          ! trend ta()/sa() and the previous one tdta()/tdsa() and substract 
(  723)          ! the vertical Gent velocity trend tavg()/savg() (zero if not used)
(  724)          tavg(:,:,:) = ta(:,:,:) - tdta(:,:,:) - tavg(:,:,:) 
(  725)          savg(:,:,:) = sa(:,:,:) - tdsa(:,:,:) - savg(:,:,:) 
(  726)       ENDIF
(  727) 






PGF90 (Version     10.5)          09/23/2010  09:50:18      page 25

(  728)    END SUBROUTINE tra_zdf_iso
(  729) 
(  730) 
(  731) # 740
(  740)    !!==============================================================================
(  741) END MODULE trazdf_iso_vopt
