


PGF90 (Version     10.5)          09/23/2010  09:49:16      page 1

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

Filename: traqsr.F90

(    1) MODULE traqsr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  traqsr  ***
(    4)    !! Ocean physics: solar radiation penetration in the top ocean levels
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_qsr      : trend due to the solar radiation penetration
(    9)    !!   tra_qsr_init : solar radiation penetration initialization
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE trc_oce         ! share SMS/Ocean variables
(   18)    USE ocesbc          ! thermohaline fluxes
(   19)    USE phycst          ! physical constants
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_qsr      ! routine called by step.F90 (ln_traqsr=T)
(   27)    PUBLIC tra_qsr_init ! routine called by opa.F90
(   28) 
(   29)    !! * Shared module variables
(   30)    LOGICAL, PUBLIC ::   ln_traqsr = .TRUE.   !: qsr flag (Default=T)
(   31) 
(   32)    !! * Module variables
(   33)    REAL(wp), PUBLIC ::  & !!! * penetrative solar radiation namelist *
(   34)       rabs = 0.58_wp,   &  ! fraction associated with xsi1
(   35)       xsi1 = 0.35_wp,   &  ! first depth of extinction 
(   36)       xsi2 = 23.0_wp       ! second depth of extinction 
(   37)       !                    ! (default values: water type Ib)
(   38)    LOGICAL ::           & 
(   39)       ln_qsr_sms = .false. ! flag to use or not the biological 
(   40)       !                    ! fluxes for light 
(   41)    
(   42)    INTEGER ::    &
(   43)       nksr                 ! number of levels
(   44)    REAL(wp), DIMENSION(jpk) ::   &
(   45)       gdsr                 ! profile of the solar flux penetration
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 2

(   47)    !! * Substitutions
(   48) #  include "domzgr_substitute.h90"
(   49) #  include "vectopt_loop_substitute.h90"
(   50)    !!----------------------------------------------------------------------
(   51)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   52)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traqsr.F90,v 1.10 2005/09/22 10:55:46 opalod Exp $ 
(   53)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56) CONTAINS
(   57) 
(   58)    SUBROUTINE tra_qsr( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                  ***  ROUTINE tra_qsr  ***
(   61)       !!
(   62)       !! ** Purpose :   Compute the temperature trend due to the solar radiation
(   63)       !!      penetration and add it to the general temperature trend.
(   64)       !!
(   65)       !! ** Method  : The profile of the solar radiation within the ocean is
(   66)       !!      defined through two penetration length scale (xsr1,xsr2) and a 
(   67)       !!      ratio (rabs) as :
(   68)       !!         I(k) = Qsr*( rabs*EXP(z(k)/xsr1) + (1.-rabs)*EXP(z(k)/xsr2) )
(   69)       !!         The temperature trend associated with the solar radiation
(   70)       !!      penetration is given by :
(   71)       !!            zta = 1/e3t dk[ I ] / (rau0*Cp)
(   72)       !!         At the bottom, boudary condition for the radiation is no flux :
(   73)       !!      all heat which has not been absorbed in the above levels is put
(   74)       !!      in the last ocean level.
(   75)       !!         In z-coordinate case, the computation is only done down to the
(   76)       !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients 
(   77)       !!      used for the computation are calculated one for once as they
(   78)       !!      depends on k only.
(   79)       !!
(   80)       !! ** Action  : - update ta with the penetrative solar radiation trend
(   81)       !!              - save the trend in ttrd ('key_trdtra')
(   82)       !!
(   83)       !! History :
(   84)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   85)       !!   7.0  !  91-11  (G. Madec)
(   86)       !!        !  96-01  (G. Madec)  s-coordinates
(   87)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   88)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used     
(   91)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   92)                          ztdsa => va      ! use va as 3D workspace   
(   93) 
(   94)       !! * Arguments
(   95)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step
(   96) 
(   97)       !! * Local declarations
(   98)       INTEGER ::    ji, jj, jk            ! dummy loop indexes
(   99)       REAL(wp) ::   zc0, zta              ! temporary scalars
(  100)       REAL(wp) ::   zc1 , zc2 ,        &  ! temporary scalars
(  101)                     zdp1, zdp2            !
(  102)       !!----------------------------------------------------------------------
(  103) 
(  104)       IF( kt == nit000 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 3

(  105)          IF ( lwp )  WRITE(numout,*)
(  106)          IF ( lwp )  WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
(  107)          IF ( lwp )  WRITE(numout,*) '~~~~~~~'
(  108)       ENDIF
(  109) 
(  110)       ! Save ta and sa trends
(  111)       IF( l_trdtra )   THEN
(  112)          ztdta(:,:,:) = ta(:,:,:) 
(  113)          ztdsa(:,:,:) = 0.e0
(  114)       ENDIF
(  115) 
(  116)       IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN    !  Biological fluxes  !
(  117)          !                                      ! =================== !
(  118)          !
(  119)          !                                                ! ===============
(  120)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  121)             !                                             ! ===============
(  122)             DO jj = 2, jpjm1
(  123)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  124)                   
(  125)                   zc0 = ro0cpr  / fse3t(ji,jj,jk)      ! compute the qsr trend
(  126)                   zta = zc0 * ( etot3(ji,jj,jk  ) * tmask(ji,jj,jk)     &
(  127)                      &        - etot3(ji,jj,jk+1) * tmask(ji,jj,jk+1) )
(  128)                   
(  129)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  130)                   
(  131)                END DO
(  132)             END DO
(  133)             !                                             ! ===============
(  134)          END DO                                           !   End of slab
(  135)          !                                                ! ===============
(  136)          ! save the trends for diagnostic
(  137)          ! qsr tracers trends
(  138)          IF( l_trdtra )   THEN
(  139)             ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  140)             CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  141)          ENDIF
(  142) 
(  143)       ELSE
(  144)          !                                                ! =================== !
(  145)          IF( lk_sco ) THEN                                !    s-coordinate     !
(  146)          !                                                ! =================== !
(  147)          !
(  148)          !                                                   ! ===============
(  149)             DO jk = 1, jpkm1                                 ! Horizontal slab
(  150)             !                                                ! ===============
(  151)                DO jj = 2, jpjm1
(  152)                   DO ji = fs_2, fs_jpim1   ! vector opt.
(  153) 
(  154)                      zdp1 = -fsdepw(ji,jj,jk  )              ! compute the qsr trend
(  155)                      zdp2 = -fsdepw(ji,jj,jk+1)
(  156)                      zc0  = qsr(ji,jj) * ro0cpr / fse3t(ji,jj,jk)
(  157)                      zc1  =   (  rabs * EXP(zdp1/xsi1) + (1.-rabs) * EXP(zdp1/xsi2)  )
(  158)                      zc2  = - (  rabs * EXP(zdp2/xsi1) + (1.-rabs) * EXP(zdp2/xsi2)  )
(  159)                      zta  = zc0 * (  zc1 * tmask(ji,jj,jk) + zc2 * tmask(ji,jj,jk+1)  )
(  160)                      
(  161)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  162)                      






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 4

(  163)                   END DO
(  164)                END DO
(  165)                !                                                ! ===============
(  166)             END DO                                              !   End of slab
(  167)             !                                                   ! ===============
(  168)             ! save the trends for diagnostic
(  169)             ! qsr tracers trends
(  170)             IF( l_trdtra )   THEN
(  171)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  172)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  173)             ENDIF
(  174)             !
(  175)          ENDIF
(  176)          !                                                ! =================== !
(  177)          IF( lk_zps ) THEN                                !    partial steps    !
(  178)             !                                             ! =================== !
(  179)             !
(  180)             !                                                ! ===============
(  181)             DO jk = 1, nksr                                  ! Horizontal slab
(  182)                !                                             ! ===============
(  183)                DO jj = 2, jpjm1
(  184)                   DO ji = fs_2, fs_jpim1   ! vector opt.
(  185)                      
(  186)                      zc0 = qsr(ji,jj) / fse3t(ji,jj,jk)      ! compute the qsr trend
(  187)                      zta = zc0 * ( gdsr(jk) * tmask(ji,jj,jk) - gdsr(jk+1) * tmask(ji,jj,jk+1) )
(  188)                      
(  189)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  190)                      
(  191)                   END DO
(  192)                END DO
(  193)                !                                             ! ===============
(  194)             END DO                                           !   End of slab
(  195)             !                                                ! ===============
(  196)             ! save the trends for diagnostic
(  197)             ! qsr tracers trends
(  198)             IF( l_trdtra )   THEN
(  199)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  200)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  201)             ENDIF
(  202)             !
(  203)          ENDIF
(  204)          !                                                ! =================== !
(  205)          IF( lk_zco ) THEN                                !     z-coordinate    !
(  206)             !                                             ! =================== !
(  207)             !
(  208)             !                                                ! ===============
(  209)             DO jk = 1, nksr                                  ! Horizontal slab
(  210)                !                                             ! ===============
(  211)                zc0 = 1. / fse3t(1,1,jk)
(  212)                DO jj = 2, jpjm1
(  213)                   DO ji = fs_2, fs_jpim1   ! vector opt.
(  214)                      !                                       ! compute qsr forcing trend
(  215)                      zta = qsr(ji,jj) * zc0 * ( gdsr(jk)*tmask(ji,jj,jk) - gdsr(jk+1)*tmask(ji,jj,jk+1) )
(  216)                      
(  217)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  218)                      
(  219)                   END DO
(  220)                END DO






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 5

(  221)                !                                             ! ===============
(  222)             END DO                                           !   End of slab
(  223)             !                                                ! ===============
(  224)             ! save the trends for diagnostic
(  225)             ! qsr tracers trends
(  226)             IF( l_trdtra )   THEN
(  227)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  228)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  229)             ENDIF
(  230)             !
(  231)          ENDIF
(  232)          !
(  233)       ENDIF
(  234) 
(  235) 
(  236)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  237)          CALL prt_ctl(tab3d_1=ta, clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta')
(  238)       ENDIF
(  239) 
(  240)    END SUBROUTINE tra_qsr
(  241) 
(  242) 
(  243)    SUBROUTINE tra_qsr_init
(  244)       !!----------------------------------------------------------------------
(  245)       !!                  ***  ROUTINE tra_qsr_init  ***
(  246)       !!
(  247)       !! ** Purpose :   Initialization for the penetrative solar radiation
(  248)       !!
(  249)       !! ** Method  :   The profile of solar radiation within the ocean is set
(  250)       !!      from two length scale of penetration (xsr1,xsr2) and a ratio
(  251)       !!      (rabs). These parameters are read in the namqsr namelist. The
(  252)       !!      default values correspond to clear water (type I in Jerlov' 
(  253)       !!      (1968) classification.
(  254)       !!         called by tra_qsr at the first timestep (nit000)
(  255)       !!
(  256)       !! ** Action  : - initialize xsr1, xsr2 and rabs
(  257)       !!
(  258)       !! Reference :
(  259)       !!   Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
(  260)       !!
(  261)       !! History :
(  262)       !!   8.5  !  02-06  (G. Madec) Original code
(  263)       !!----------------------------------------------------------------------
(  264)       !! * Local declarations
(  265)       INTEGER ::    ji,jj,jk, &  ! dummy loop index
(  266)                     indic        ! temporary integer
(  267)       REAL(wp) ::   zdp1         ! temporary scalar
(  268) 
(  269)       NAMELIST/namqsr/ ln_traqsr, rabs, xsi1, xsi2, ln_qsr_sms
(  270)       !!----------------------------------------------------------------------
(  271) 
(  272)       ! Read Namelist namqsr : ratio and length of penetration
(  273)       ! --------------------
(  274)       REWIND ( numnam )
(  275)       READ   ( numnam, namqsr )
(  276) 
(  277)       ! Parameter control and print
(  278)       ! ---------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 6

(  279)       IF( ln_traqsr  ) THEN
(  280)         IF ( lwp ) THEN
(  281)          WRITE(numout,*)
(  282)          WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
(  283)          WRITE(numout,*) '~~~~~~~~~~~~'
(  284)          WRITE(numout,*) '    Namelist namqsr : set the parameter of penetration'
(  285)          WRITE(numout,*) '        fraction associated with xsi     rabs        = ',rabs
(  286)          WRITE(numout,*) '        first depth of extinction        xsi1        = ',xsi1
(  287)          WRITE(numout,*) '        second depth of extinction       xsi2        = ',xsi2
(  288)          IF( lk_qsr_sms ) THEN
(  289)             WRITE(numout,*) '     Biological fluxes for light(Y/N) ln_qsr_sms  = ',ln_qsr_sms
(  290)          ENDIF
(  291)          WRITE(numout,*) ' '
(  292)         END IF
(  293)       ELSE
(  294)         IF ( lwp ) THEN
(  295)          WRITE(numout,*)
(  296)          WRITE(numout,*) 'tra_qsr_init : NO solar flux penetration'
(  297)          WRITE(numout,*) '~~~~~~~~~~~~'
(  298)         END IF
(  299)       ENDIF
(  300) 
(  301)       IF( rabs > 1.e0 .OR. rabs < 0.e0 .OR. xsi1 < 0.e0 .OR. xsi2 < 0.e0 ) THEN
(  302)          IF(lwp) WRITE(numout,cform_err)
(  303)          IF(lwp) WRITE(numout,*) '             0<rabs<1, 0<xsi1, or 0<xsi2 not satisfied'
(  304)          nstop = nstop + 1
(  305)       ENDIF
(  306) 
(  307) 
(  308)       ! Initialization
(  309)       ! --------------
(  310)       IF( .NOT. lk_sco ) THEN
(  311)          ! z-coordinate with or without partial step : same before last ocean w-level everywhere
(  312)          gdsr(:) = 0.e0
(  313)          DO jk = 1, jpk
(  314)             zdp1 = -fsdepw(1,1,jk)
(  315)             gdsr(jk) = ro0cpr * (  rabs  * EXP( zdp1/xsi1 ) + (1.-rabs) * EXP( zdp1/xsi2 )  )
(  316)             IF ( gdsr(jk) <= 1.e-10 ) EXIT
(  317)          END DO
(  318)          indic = 0
(  319)          DO jk = 1, jpk
(  320)             IF( gdsr(jk) <= 1.e-15 .AND. indic == 0 ) THEN
(  321)                gdsr(jk) = 0.e0
(  322)                nksr = jk
(  323)                !!bug Edmee chg res   nksr = jk - 1
(  324)                indic = 1
(  325)             ENDIF
(  326)          END DO
(  327)          nksr = MIN( nksr, jpkm1 )
(  328)          IF(lwp) THEN
(  329)             WRITE(numout,*)
(  330)             WRITE(numout,*) '        - z-coordinate, level max of computation =', nksr
(  331)             WRITE(numout,*) '             profile of coef. of penetration:'
(  332)             WRITE(numout,"('              ',7e11.2)") ( gdsr(jk), jk = 1, nksr )
(  333)             WRITE(numout,*)
(  334)          ENDIF
(  335)          ! Initialisation of Biological fluxes for light here because
(  336)          ! the optical biological model is call after the dynamical one






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 7

(  337)          IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN
(  338)             DO jk = 1, jpkm1
(  339)                DO jj = 1, jpj
(  340)                   DO ji = 1, jpi
(  341)                      etot3(ji,jj,jk) = qsr(ji,jj) * gdsr(jk) * tmask(ji,jj,jk) / ro0cpr
(  342)                   END DO
(  343)                END DO
(  344)             END DO
(  345)          ENDIF
(  346) 
(  347)       ENDIF
(  348) 
(  349)    END SUBROUTINE tra_qsr_init
(  350) 
(  351)    !!======================================================================
(  352) END MODULE traqsr
















































PGF90 (Version     10.5)          09/23/2010  09:49:16      page 8

(    1) # 1 "traqsr.F90"
(    1) MODULE traqsr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  traqsr  ***
(    4)    !! Ocean physics: solar radiation penetration in the top ocean levels
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_qsr      : trend due to the solar radiation penetration
(    9)    !!   tra_qsr_init : solar radiation penetration initialization
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracers trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE trc_oce         ! share SMS/Ocean variables
(   18)    USE ocesbc          ! thermohaline fluxes
(   19)    USE phycst          ! physical constants
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_qsr      ! routine called by step.F90 (ln_traqsr=T)
(   27)    PUBLIC tra_qsr_init ! routine called by opa.F90
(   28) 
(   29)    !! * Shared module variables
(   30)    LOGICAL, PUBLIC ::   ln_traqsr = .TRUE.   !: qsr flag (Default=T)
(   31) 
(   32)    !! * Module variables
(   33)    REAL(wp), PUBLIC ::  & !!! * penetrative solar radiation namelist *
(   34)       rabs = 0.58_wp,   &  ! fraction associated with xsi1
(   35)       xsi1 = 0.35_wp,   &  ! first depth of extinction 
(   36)       xsi2 = 23.0_wp       ! second depth of extinction 
(   37)       !                    ! (default values: water type Ib)
(   38)    LOGICAL ::           & 
(   39)       ln_qsr_sms = .false. ! flag to use or not the biological 
(   40)       !                    ! fluxes for light 
(   41)    
(   42)    INTEGER ::    &
(   43)       nksr                 ! number of levels
(   44)    REAL(wp), DIMENSION(jpk) ::   &
(   45)       gdsr                 ! profile of the solar flux penetration
(   46) 
(   47)    !! * Substitutions
(   48) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 9

(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 49 "traqsr.F90"
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
(   12) # 50 "traqsr.F90"
(   50) # 50
(   50)    !!----------------------------------------------------------------------
(   51)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   52)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traqsr.F90,v 1.10 2005/09/22 10:55:46 opalod Exp $ 
(   53)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   54)    !!----------------------------------------------------------------------
(   55) 
(   56) CONTAINS
(   57) 
(   58)    SUBROUTINE tra_qsr( kt )
(   59)       !!----------------------------------------------------------------------
(   60)       !!                  ***  ROUTINE tra_qsr  ***
(   61)       !!
(   62)       !! ** Purpose :   Compute the temperature trend due to the solar radiation
(   63)       !!      penetration and add it to the general temperature trend.
(   64)       !!
(   65)       !! ** Method  : The profile of the solar radiation within the ocean is
(   66)       !!      defined through two penetration length scale (xsr1,xsr2) and a 
(   67)       !!      ratio (rabs) as :
(   68)       !!         I(k) = Qsr*( rabs*EXP(z(k)/xsr1) + (1.-rabs)*EXP(z(k)/xsr2) )
(   69)       !!         The temperature trend associated with the solar radiation
(   70)       !!      penetration is given by :
(   71)       !!            zta = 1/e3t dk[ I ] / (rau0*Cp)
(   72)       !!         At the bottom, boudary condition for the radiation is no flux :
(   73)       !!      all heat which has not been absorbed in the above levels is put
(   74)       !!      in the last ocean level.
(   75)       !!         In z-coordinate case, the computation is only done down to the
(   76)       !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients 
(   77)       !!      used for the computation are calculated one for once as they
(   78)       !!      depends on k only.
(   79)       !!
(   80)       !! ** Action  : - update ta with the penetrative solar radiation trend
(   81)       !!              - save the trend in ttrd ('key_trdtra')
(   82)       !!
(   83)       !! History :
(   84)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   85)       !!   7.0  !  91-11  (G. Madec)
(   86)       !!        !  96-01  (G. Madec)  s-coordinates
(   87)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 10

(   88)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   89)       !!----------------------------------------------------------------------
(   90)       !! * Modules used     
(   91)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   92)                          ztdsa => va      ! use va as 3D workspace   
(   93) 
(   94)       !! * Arguments
(   95)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step
(   96) 
(   97)       !! * Local declarations
(   98)       INTEGER ::    ji, jj, jk            ! dummy loop indexes
(   99)       REAL(wp) ::   zc0, zta              ! temporary scalars
(  100)       REAL(wp) ::   zc1 , zc2 ,        &  ! temporary scalars
(  101)                     zdp1, zdp2            !
(  102)       !!----------------------------------------------------------------------
(  103) 
(  104)       IF( kt == nit000 ) THEN
(  105)          IF ( lwp )  WRITE(numout,*)
(  106)          IF ( lwp )  WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
(  107)          IF ( lwp )  WRITE(numout,*) '~~~~~~~'
(  108)       ENDIF
(  109) 
(  110)       ! Save ta and sa trends
(  111)       IF( l_trdtra )   THEN
(  112)          ztdta(:,:,:) = ta(:,:,:) 
(  113)          ztdsa(:,:,:) = 0.e0
(  114)       ENDIF
(  115) 
(  116)       IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN    !  Biological fluxes  !
(  117)          !                                      ! =================== !
(  118)          !
(  119)          !                                                ! ===============
(  120)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  121)             !                                             ! ===============
(  122)             DO jj = 2, jpjm1
(  123)                DO ji = 2, jpim1   ! vector opt.
(  124)                   
(  125)                   zc0 = ro0cpr  / e3t(jk)      ! compute the qsr trend
(  126)                   zta = zc0 * ( etot3(ji,jj,jk  ) * tmask(ji,jj,jk)     &
(  127)                      &        - etot3(ji,jj,jk+1) * tmask(ji,jj,jk+1) )
(  128)                   
(  129)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  130)                   
(  131)                END DO
(  132)             END DO
(  133)             !                                             ! ===============
(  134)          END DO                                           !   End of slab
(  135)          !                                                ! ===============
(  136)          ! save the trends for diagnostic
(  137)          ! qsr tracers trends
(  138)          IF( l_trdtra )   THEN
(  139)             ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  140)             CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  141)          ENDIF
(  142) 
(  143)       ELSE
(  144)          !                                                ! =================== !
(  145)          IF( lk_sco ) THEN                                !    s-coordinate     !






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 11

(  146)          !                                                ! =================== !
(  147)          !
(  148)          !                                                   ! ===============
(  149)             DO jk = 1, jpkm1                                 ! Horizontal slab
(  150)             !                                                ! ===============
(  151)                DO jj = 2, jpjm1
(  152)                   DO ji = 2, jpim1   ! vector opt.
(  153) 
(  154)                      zdp1 = -gdepw(jk  )              ! compute the qsr trend
(  155)                      zdp2 = -gdepw(jk+1)
(  156)                      zc0  = qsr(ji,jj) * ro0cpr / e3t(jk)
(  157)                      zc1  =   (  rabs * EXP(zdp1/xsi1) + (1.-rabs) * EXP(zdp1/xsi2)  )
(  158)                      zc2  = - (  rabs * EXP(zdp2/xsi1) + (1.-rabs) * EXP(zdp2/xsi2)  )
(  159)                      zta  = zc0 * (  zc1 * tmask(ji,jj,jk) + zc2 * tmask(ji,jj,jk+1)  )
(  160)                      
(  161)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  162)                      
(  163)                   END DO
(  164)                END DO
(  165)                !                                                ! ===============
(  166)             END DO                                              !   End of slab
(  167)             !                                                   ! ===============
(  168)             ! save the trends for diagnostic
(  169)             ! qsr tracers trends
(  170)             IF( l_trdtra )   THEN
(  171)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  172)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  173)             ENDIF
(  174)             !
(  175)          ENDIF
(  176)          !                                                ! =================== !
(  177)          IF( lk_zps ) THEN                                !    partial steps    !
(  178)             !                                             ! =================== !
(  179)             !
(  180)             !                                                ! ===============
(  181)             DO jk = 1, nksr                                  ! Horizontal slab
(  182)                !                                             ! ===============
(  183)                DO jj = 2, jpjm1
(  184)                   DO ji = 2, jpim1   ! vector opt.
(  185)                      
(  186)                      zc0 = qsr(ji,jj) / e3t(jk)      ! compute the qsr trend
(  187)                      zta = zc0 * ( gdsr(jk) * tmask(ji,jj,jk) - gdsr(jk+1) * tmask(ji,jj,jk+1) )
(  188)                      
(  189)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  190)                      
(  191)                   END DO
(  192)                END DO
(  193)                !                                             ! ===============
(  194)             END DO                                           !   End of slab
(  195)             !                                                ! ===============
(  196)             ! save the trends for diagnostic
(  197)             ! qsr tracers trends
(  198)             IF( l_trdtra )   THEN
(  199)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  200)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  201)             ENDIF
(  202)             !
(  203)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 12

(  204)          !                                                ! =================== !
(  205)          IF( lk_zco ) THEN                                !     z-coordinate    !
(  206)             !                                             ! =================== !
(  207)             !
(  208)             !                                                ! ===============
(  209)             DO jk = 1, nksr                                  ! Horizontal slab
(  210)                !                                             ! ===============
(  211)                zc0 = 1. / e3t(jk)
(  212)                DO jj = 2, jpjm1
(  213)                   DO ji = 2, jpim1   ! vector opt.
(  214)                      !                                       ! compute qsr forcing trend
(  215)                      zta = qsr(ji,jj) * zc0 * ( gdsr(jk)*tmask(ji,jj,jk) - gdsr(jk+1)*tmask(ji,jj,jk+1) )
(  216)                      
(  217)                      ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
(  218)                      
(  219)                   END DO
(  220)                END DO
(  221)                !                                             ! ===============
(  222)             END DO                                           !   End of slab
(  223)             !                                                ! ===============
(  224)             ! save the trends for diagnostic
(  225)             ! qsr tracers trends
(  226)             IF( l_trdtra )   THEN
(  227)                ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  228)                CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
(  229)             ENDIF
(  230)             !
(  231)          ENDIF
(  232)          !
(  233)       ENDIF
(  234) 
(  235) 
(  236)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  237)          CALL prt_ctl(tab3d_1=ta, clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta')
(  238)       ENDIF
(  239) 
(  240)    END SUBROUTINE tra_qsr
(  241) 
(  242) 
(  243)    SUBROUTINE tra_qsr_init
(  244)       !!----------------------------------------------------------------------
(  245)       !!                  ***  ROUTINE tra_qsr_init  ***
(  246)       !!
(  247)       !! ** Purpose :   Initialization for the penetrative solar radiation
(  248)       !!
(  249)       !! ** Method  :   The profile of solar radiation within the ocean is set
(  250)       !!      from two length scale of penetration (xsr1,xsr2) and a ratio
(  251)       !!      (rabs). These parameters are read in the namqsr namelist. The
(  252)       !!      default values correspond to clear water (type I in Jerlov' 
(  253)       !!      (1968) classification.
(  254)       !!         called by tra_qsr at the first timestep (nit000)
(  255)       !!
(  256)       !! ** Action  : - initialize xsr1, xsr2 and rabs
(  257)       !!
(  258)       !! Reference :
(  259)       !!   Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
(  260)       !!
(  261)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 13

(  262)       !!   8.5  !  02-06  (G. Madec) Original code
(  263)       !!----------------------------------------------------------------------
(  264)       !! * Local declarations
(  265)       INTEGER ::    ji,jj,jk, &  ! dummy loop index
(  266)                     indic        ! temporary integer
(  267)       REAL(wp) ::   zdp1         ! temporary scalar
(  268) 
(  269)       NAMELIST/namqsr/ ln_traqsr, rabs, xsi1, xsi2, ln_qsr_sms
(  270)       !!----------------------------------------------------------------------
(  271) 
(  272)       ! Read Namelist namqsr : ratio and length of penetration
(  273)       ! --------------------
(  274)       REWIND ( numnam )
(  275)       READ   ( numnam, namqsr )
(  276) 
(  277)       ! Parameter control and print
(  278)       ! ---------------------------
(  279)       IF( ln_traqsr  ) THEN
(  280)         IF ( lwp ) THEN
(  281)          WRITE(numout,*)
(  282)          WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
(  283)          WRITE(numout,*) '~~~~~~~~~~~~'
(  284)          WRITE(numout,*) '    Namelist namqsr : set the parameter of penetration'
(  285)          WRITE(numout,*) '        fraction associated with xsi     rabs        = ',rabs
(  286)          WRITE(numout,*) '        first depth of extinction        xsi1        = ',xsi1
(  287)          WRITE(numout,*) '        second depth of extinction       xsi2        = ',xsi2
(  288)          IF( lk_qsr_sms ) THEN
(  289)             WRITE(numout,*) '     Biological fluxes for light(Y/N) ln_qsr_sms  = ',ln_qsr_sms
(  290)          ENDIF
(  291)          WRITE(numout,*) ' '
(  292)         END IF
(  293)       ELSE
(  294)         IF ( lwp ) THEN
(  295)          WRITE(numout,*)
(  296)          WRITE(numout,*) 'tra_qsr_init : NO solar flux penetration'
(  297)          WRITE(numout,*) '~~~~~~~~~~~~'
(  298)         END IF
(  299)       ENDIF
(  300) 
(  301)       IF( rabs > 1.e0 .OR. rabs < 0.e0 .OR. xsi1 < 0.e0 .OR. xsi2 < 0.e0 ) THEN
(  302)          IF(lwp) WRITE(numout,cform_err)
(  303)          IF(lwp) WRITE(numout,*) '             0<rabs<1, 0<xsi1, or 0<xsi2 not satisfied'
(  304)          nstop = nstop + 1
(  305)       ENDIF
(  306) 
(  307) 
(  308)       ! Initialization
(  309)       ! --------------
(  310)       IF( .NOT. lk_sco ) THEN
(  311)          ! z-coordinate with or without partial step : same before last ocean w-level everywhere
(  312)          gdsr(:) = 0.e0
(  313)          DO jk = 1, jpk
(  314)             zdp1 = -gdepw(jk)
(  315)             gdsr(jk) = ro0cpr * (  rabs  * EXP( zdp1/xsi1 ) + (1.-rabs) * EXP( zdp1/xsi2 )  )
(  316)             IF ( gdsr(jk) <= 1.e-10 ) EXIT
(  317)          END DO
(  318)          indic = 0
(  319)          DO jk = 1, jpk






PGF90 (Version     10.5)          09/23/2010  09:49:16      page 14

(  320)             IF( gdsr(jk) <= 1.e-15 .AND. indic == 0 ) THEN
(  321)                gdsr(jk) = 0.e0
(  322)                nksr = jk
(  323)                !!bug Edmee chg res   nksr = jk - 1
(  324)                indic = 1
(  325)             ENDIF
(  326)          END DO
(  327)          nksr = MIN( nksr, jpkm1 )
(  328)          IF(lwp) THEN
(  329)             WRITE(numout,*)
(  330)             WRITE(numout,*) '        - z-coordinate, level max of computation =', nksr
(  331)             WRITE(numout,*) '             profile of coef. of penetration:'
(  332)             WRITE(numout,"('              ',7e11.2)") ( gdsr(jk), jk = 1, nksr )
(  333)             WRITE(numout,*)
(  334)          ENDIF
(  335)          ! Initialisation of Biological fluxes for light here because
(  336)          ! the optical biological model is call after the dynamical one
(  337)          IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN
(  338)             DO jk = 1, jpkm1
(  339)                DO jj = 1, jpj
(  340)                   DO ji = 1, jpi
(  341)                      etot3(ji,jj,jk) = qsr(ji,jj) * gdsr(jk) * tmask(ji,jj,jk) / ro0cpr
(  342)                   END DO
(  343)                END DO
(  344)             END DO
(  345)          ENDIF
(  346) 
(  347)       ENDIF
(  348) 
(  349)    END SUBROUTINE tra_qsr_init
(  350) 
(  351)    !!======================================================================
(  352) END MODULE traqsr
