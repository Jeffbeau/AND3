


PGF90 (Version     10.5)          09/23/2010  09:50:00      page 1

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

Filename: icestp1d.F90

(    1) MODULE icestp1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE icestp1d   ***
(    4)    !!   Sea-Ice model : 1D LIM Sea ice model time-stepping
(    5)    !!======================================================================
(    6) #if defined key_cfg_1d && defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'                                       1D Configuration
(    9)    !!   'key_ice_lim' :                                   Lim sea-ice model
(   10)    !!----------------------------------------------------------------------
(   11)    !!   ice_stp_1d       : sea-ice model time-stepping
(   12)    !!----------------------------------------------------------------------
(   13)    USE dom_oce
(   14)    USE oce  ! dynamics and tracers variables
(   15)    USE in_out_manager
(   16)    USE ice_oce         ! ice variables
(   17)    USE flx_oce         ! forcings variables
(   18)    USE dom_ice
(   19)    USE cpl_oce
(   20)    USE blk_oce
(   21)    USE daymod
(   22)    USE phycst          ! Define parameters for the routines
(   23)    USE taumod
(   24)    USE ice
(   25)    USE iceini
(   26)    USE lbclnk
(   27)    USE limdyn
(   28)    USE limtrp
(   29)    USE limthd
(   30)    USE limflx
(   31)    USE limdia
(   32)    USE limwri
(   33)    USE limrst
(   34) 
(   35)    USE ocesbc
(   36)    USE flxmod
(   37)    USE flxrnf
(   38)    USE tradmp         ! damping salinity trend
(   39)    USE dtatem
(   40)    USE dtasal
(   41)    USE ocfzpt
(   42)    USE prtctl          ! Print control
(   43) 
(   44) 
(   45)    IMPLICIT NONE
(   46)    PRIVATE






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 2

(   47) 
(   48)    !! * Routine accessibility
(   49)    PUBLIC ice_stp_1d  ! called by step.F90
(   50) 
(   51)    !! * Substitutions
(   52) #  include "domzgr_substitute.h90"
(   53) #  include "vectopt_loop_substitute.h90"
(   54)    !!-----------------------------------------------------
(   55)    !!   LIM 2.0 , UCL-LOCEAN-IPSL (2005) 
(   56)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/icestp1d.F90,v 1.6 2006/04/19 14:43:13 opalod Exp $ 
(   57)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   58)    !!-----------------------------------------------------
(   59) 
(   60) CONTAINS
(   61) 
(   62)    SUBROUTINE ice_stp_1d ( kt )
(   63)       !!---------------------------------------------------------------------
(   64)       !!                  ***  ROUTINE ice_stp_1d  ***
(   65)       !!                   
(   66)       !! ** Purpose :   Louvain la Neuve Sea Ice Model time stepping 
(   67)       !!
(   68)       !! ** Action  : - call the ice dynamics routine 
(   69)       !!              - call the ice advection/diffusion routine 
(   70)       !!              - call the ice thermodynamics routine 
(   71)       !!              - call the routine that computes mass and 
(   72)       !!                heat fluxes at the ice/ocean interface
(   73)       !!              - save the outputs 
(   74)       !!              - save the outputs for restart when necessary
(   75)       !!
(   76)       !! History :
(   77)       !!   1.0  !  99-11  (M. Imbard)  Original code
(   78)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(   79)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(   80)       !!   9.0  !  04-10  (C. Ethe) 1D configuration
(   81)       !!----------------------------------------------------------------------
(   82)       !! * Arguments
(   83)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   84) 
(   85)       !! * Local declarations
(   86)       INTEGER   ::   ji, jj   ! dummy loop indices
(   87) 
(   88)       REAL(wp) , DIMENSION(jpi,jpj)    :: &
(   89)          zsss_io, zsss2_io, zsss3_io          ! tempory workspaces
(   90)       REAL(wp)  :: ztair2
(   91)       !!----------------------------------------------------------------------
(   92) 
(   93)       IF( kt == nit000 ) THEN
(   94)          IF(lwp) WRITE(numout,*)
(   95)          IF(lwp) WRITE(numout,*) 'ice_stp_1d : Louvain la Neuve Ice Model (LIM)' 
(   96)          IF(lwp) WRITE(numout,*) '~~~~~~~   forced case using bulk formulea'
(   97)          !  Initialize fluxes fields
(   98)          gtaux(:,:) = 0.e0
(   99)          gtauy(:,:) = 0.e0
(  100)       ENDIF
(  101) 
(  102)       ! Temperature , salinity and horizonta wind
(  103)       ! sst_io  and sss_io, u_io and v_io  are initialized at nit000 in limistate.F90 (or limrst.F90) with :
(  104)       !              sst_io = sst_io + (nfice - 1) * (tn(:,:,1)+rt0 )






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 3

(  105)       !              sss_io = sss_io + (nfice - 1) * sn(:,:,1)
(  106)       !              u_io   = u_io   + (nfice - 1) * un(:,:,1)
(  107)       !              v_io   = v_io   + (nfice - 1) * vn(:,:,1)
(  108)       !    cumulate fields
(  109)       !
(  110)       sst_io(:,:) = sst_io(:,:) + tn(:,:,1) + rt0
(  111)       sss_io(:,:) = sss_io(:,:) + sn(:,:,1)
(  112)       u_io  (:,:) = u_io  (:,:) + un(:,:,1)
(  113)       v_io  (:,:) = v_io  (:,:) + vn(:,:,1)
(  114)  
(  115) 
(  116)       
(  117)       IF( MOD( kt-1, nfice ) == 0 ) THEN
(  118)          
(  119)          ! The LIM model is going to be call
(  120)          sst_io(:,:) = sst_io(:,:) / FLOAT( nfice ) * tmask(:,:,1)
(  121)          sss_io(:,:) = sss_io(:,:) / FLOAT( nfice )
(  122)          u_io  (:,:) = u_io  (:,:) / FLOAT( nfice )
(  123)          v_io  (:,:) = v_io  (:,:) / FLOAT( nfice )
(  124)          gtaux (:,:) = taux  (:,:)
(  125)          gtauy (:,:) = tauy  (:,:)
(  126) 
(  127)          zsss_io (:,:) = SQRT( sss_io(:,:) ) 
(  128)          zsss2_io(:,:) =  sss_io(:,:) *  sss_io(:,:)
(  129)          zsss3_io(:,:) = zsss_io(:,:) * zsss_io(:,:) * zsss_io(:,:)
(  130) 
(  131)          DO jj = 1, jpj
(  132)             DO ji = 1, jpi
(  133)                tfu(ji,jj)  = ABS ( rt0 - 0.0575       *   sss_io(ji,jj)   &
(  134)                   &                    + 1.710523e-03 * zsss3_io(ji,jj)   &
(  135)                   &                    - 2.154996e-04 * zsss2_io(ji,jj) ) * tms(ji,jj)
(  136)             END DO
(  137)          END DO
(  138)          
(  139)          
(  140)          IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  141)             CALL prt_ctl_info('Ice Forcings ')
(  142)             CALL prt_ctl(tab2d_1=qsr_oce ,clinfo1=' qsr_oce  : ', tab2d_2=qsr_ice , clinfo2=' qsr_ice   : ')
(  143)             CALL prt_ctl(tab2d_1=qnsr_oce,clinfo1=' qnsr_oce : ', tab2d_2=qnsr_ice, clinfo2=' qnsr_ice  : ')
(  144)             CALL prt_ctl(tab2d_1=evap    ,clinfo1=' evap     : ')
(  145)             CALL prt_ctl(tab2d_1=tprecip ,clinfo1=' precip   : ', tab2d_2=sprecip , clinfo2=' Snow      : ')
(  146)             CALL prt_ctl(tab2d_1=gtaux   ,clinfo1=' u-stress : ', tab2d_2=gtauy   , clinfo2=' v-stress  : ')
(  147)             CALL prt_ctl(tab2d_1=sst_io  ,clinfo1=' sst      : ', tab2d_2=sss_io  , clinfo2=' sss       : ')
(  148)             CALL prt_ctl(tab2d_1=u_io    ,clinfo1=' u_io     : ', tab2d_2=v_io    , clinfo2=' v_io      : ')
(  149)             CALL prt_ctl(tab2d_1=hsnif   ,clinfo1=' hsnif  1 : ', tab2d_2=hicif   , clinfo2=' hicif     : ')
(  150)             CALL prt_ctl(tab2d_1=frld    ,clinfo1=' frld   1 : ', tab2d_2=sist    , clinfo2=' sist      : ')
(  151)          ENDIF
(  152) 
(  153)          ! Ice model call
(  154)          numit = numit + nfice                                       ! Friction velocity
(  155)                                                                       
(  156)          DO jj = 1, jpj
(  157)             DO ji = 1, jpi
(  158)                tio_u(ji,jj) = - gtaux(ji,jj) / rau0
(  159)                tio_v(ji,jj) = - gtauy(ji,jj) / rau0               
(  160)                ztair2       = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)           
(  161)                ust2s(ji,jj) = ( SQRT( ztair2  )  / rau0 )  * tms(ji,jj)  
(  162)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 4

(  163)          END DO
(  164) 
(  165)          !                                                           !--------------------!
(  166)          CALL lim_thd                                                ! Ice thermodynamics !
(  167)          !                                                           !--------------------!
(  168)          IF(ln_ctl) THEN
(  169)             CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  2 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  170)             CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   2 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  171)             CALL prt_ctl(tab2d_1=u_io  ,clinfo1=' u_io   4 : ', tab2d_2=v_io  , clinfo2=' v_io      : ')
(  172)             CALL prt_ctl(tab2d_1=tio_u  ,clinfo1=' tio_u  4 : ', tab2d_2=tio_v  , clinfo2=' tio_v     : ')
(  173)          ENDIF
(  174) 
(  175) 
(  176) 
(  177)          ! Mass and heat fluxes from ice to ocean
(  178)          !                                                           !------------------------------!
(  179)          CALL lim_flx                                                ! Ice/Ocean Mass & Heat fluxes !
(  180)          !                                                           !------------------------------!
(  181) 
(  182)          IF(ln_ctl) THEN
(  183)             CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  7 : ', tab2d_2=hicif , clinfo2=' hicif   : ')
(  184)             CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   7 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  185)             CALL prt_ctl(tab2d_1=tio_u  ,clinfo1=' tio_u  7 : ', tab2d_2=tio_v  , clinfo2=' tio_v     : ')
(  186)          ENDIF
(  187)          !                                                           !-------------!
(  188)          CALL lim_wri                                                ! Ice outputs !
(  189)          !                                                           !-------------!
(  190) 
(  191)          IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
(  192)             !                                                        !------------------!
(  193)             CALL lim_rst_write( numit )                              ! Ice restart file !
(  194)             !                                                        !------------------!
(  195)          ENDIF
(  196) 
(  197)          ! Re-initialization of forcings
(  198)          qsr_oce (:,:) = 0.e0
(  199)          qsr_ice (:,:) = 0.e0
(  200)          qnsr_oce(:,:) = 0.e0
(  201)          qnsr_ice(:,:) = 0.e0 
(  202)          dqns_ice(:,:) = 0.e0 
(  203)          tprecip (:,:) = 0.e0 
(  204)          sprecip (:,:) = 0.e0
(  205)          qla_ice (:,:) = 0.e0
(  206)          dqla_ice(:,:) = 0.e0
(  207)          fr1_i0  (:,:) = 0.e0
(  208)          fr2_i0  (:,:) = 0.e0
(  209)          evap    (:,:) = 0.e0
(  210) 
(  211)         CALL oce_sbc_1d ( kt )
(  212) 
(  213)       ENDIF
(  214) 
(  215)    END SUBROUTINE ice_stp_1d
(  216)    
(  217)    
(  218)    SUBROUTINE oce_sbc_1d( kt )
(  219)       !!---------------------------------------------------------------------
(  220)       !!                   ***  ROUTINE oce_sbc_1d  ***






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 5

(  221)       !!                    
(  222)       !! ** Purpose : - Ocean surface boundary conditions with LIM sea-ice
(  223)       !!        model in forced mode using bulk formulea
(  224)       !!
(  225)       !! History :
(  226)       !!   1.0  !  99-11  (M. Imbard)  Original code
(  227)       !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(  228)       !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(  229)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  230)       !!----------------------------------------------------------------------
(  231)       !! * arguments
(  232)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(  233) 
(  234)       !! * Local declarations
(  235)       INTEGER  ::   ji, jj                   ! dummy loop indices
(  236)       REAL(wp) ::   ztxy
(  237)       !!----------------------------------------------------------------------
(  238) 
(  239)       ! 1. initialization to zero at kt = nit000
(  240)       ! ---------------------------------------
(  241)       
(  242)       IF( kt == nit000 ) THEN     
(  243)          qsr    (:,:) = 0.e0
(  244)          qt     (:,:) = 0.e0
(  245)          qrp    (:,:) = 0.e0
(  246)          emp    (:,:) = 0.e0
(  247)          emps   (:,:) = 0.e0
(  248)          erp    (:,:) = 0.e0
(  249) #if ! defined key_dynspg_rl 
(  250)          dmp    (:,:) = 0.e0
(  251) #endif
(  252)       ENDIF
(  253) 
(  254)       CALL oce_sbc_dmp       ! Computation of internal and evaporation damping terms       
(  255) 
(  256)       ! Surface Ocean fluxes
(  257)       ! ====================
(  258)       
(  259)       ! Surface heat flux (W/m2)
(  260)       ! -----------------
(  261)       
(  262)       qt  (:,:) = fnsolar(:,:) + fsolar(:,:)     ! non solar heat flux + solar flux
(  263)       qsr (:,:) = fsolar(:,:)                     ! solar flux
(  264)       
(  265) #if ! defined key_dynspg_rl     
(  266)       ! total concentration/dilution effect (use on SSS)
(  267)       emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:) + empold
(  268)       
(  269)       ! total volume flux (use on sea-surface height)
(  270)       emp (:,:) = fmass(:,:) -   dmp(:,:) + runoff(:,:) + erp(:,:) + empold      
(  271) #else
(  272)       ! Rigid-lid (emp=emps=E-P-R+Erp)
(  273)       emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)     ! freshwater flux
(  274)       emp (:,:) = emps(:,:)
(  275)       
(  276) #endif
(  277)       
(  278)       ! Surface stress






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 6

(  279)       ! --------------
(  280)       
(  281)       ! update the stress beloww sea-ice area
(  282)       DO jj = 1, jpjm1
(  283)          DO ji = 1, fs_jpim1   ! vertor opt.
(  284)             ztxy        = freezn(ji,jj)             ! ice/ocean indicator at T-points
(  285)             taux(ji,jj) = (1.-ztxy) * taux(ji,jj) + ztxy * ftaux(ji,jj)    ! stress at the ocean surface
(  286)             tauy(ji,jj) = (1.-ztxy) * tauy(ji,jj) + ztxy * ftauy(ji,jj)
(  287)          END DO
(  288)       END DO
(  289)       
(  290)       ! boundary condition on the stress (taux,tauy)
(  291)       CALL lbc_lnk( taux, 'U', -1. )
(  292)       CALL lbc_lnk( tauy, 'V', -1. )
(  293)       
(  294)       ! Re-initialization of fluxes
(  295)       sst_io(:,:) = 0.e0
(  296)       sss_io(:,:) = 0.e0
(  297)       u_io  (:,:) = 0.e0
(  298)       v_io  (:,:) = 0.e0
(  299)       
(  300)       
(  301)    END SUBROUTINE oce_sbc_1d
(  302)    
(  303) #if defined key_dtasal
(  304)    !!----------------------------------------------------------------------
(  305)    !!   'key_dtasal'                                          salinity data
(  306)    !!----------------------------------------------------------------------
(  307)    SUBROUTINE oce_sbc_dmp
(  308)       !!---------------------------------------------------------------------
(  309)       !!                   ***  ROUTINE oce_sbc_dmp  ***
(  310)       !!                    
(  311)       !! ** Purpose : Computation of internal and evaporation damping terms 
(  312)       !!        for ocean surface boundary conditions 
(  313)       !!
(  314)       !! History :
(  315)       !!   9.0  !  04-01  (G. Madec, C. Ethe)  Original code
(  316)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  317)       !!----------------------------------------------------------------------
(  318)       !! * Local declarations
(  319)       INTEGER ::   ji, jj                   ! dummy loop indices
(  320)       REAL(wp) ::   zerp, zsrp
(  321) #if ! defined key_dynspg_rl
(  322)       REAL(wp) ::   zwei
(  323)       REAL(wp) ::   zerpplus(jpi,jpj), zerpminus(jpi,jpj)
(  324)       REAL(wp) ::   zplus, zminus, zadefi
(  325) # if defined key_tradmp
(  326)       INTEGER jk
(  327)       REAL(wp), DIMENSION(jpi,jpj) ::   zstrdmp
(  328) # endif
(  329) #endif
(  330)       !!----------------------------------------------------------------------
(  331) 
(  332) 
(  333)       ! sea ice indicator (1 or 0)
(  334)       DO jj = 1, jpj
(  335)          DO ji = 1, jpi
(  336)             freezn(ji,jj) = MAX(0., SIGN(1., freeze(ji,jj)-rsmall) )






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 7

(  337)          END DO
(  338)       END DO
(  339) 
(  340)       ! Initialisation
(  341)       ! --------------
(  342)       ! Restoring coefficients on SST and SSS   
(  343)       zsrp = dqdt0 * ro0cpr * rauw   ! (Kg/m2/s) 
(  344) 
(  345) #if ! defined key_dynspg_rl 
(  346)       ! Free-surface
(  347)          
(  348)       ! Internal damping
(  349) # if defined key_tradmp
(  350)       ! Vertical mean of dampind trend (computed in tradmp module)
(  351)       zstrdmp(:,:) = 0.e0
(  352)       DO jk = 1, jpk
(  353)          zstrdmp(:,:) = zstrdmp(:,:) + strdmp(:,:,jk) * fse3t(:,:,jk)
(  354)       END DO
(  355)       ! volume flux associated to internal damping to climatology
(  356)       dmp(:,:) = zstrdmp(:,:) * rauw / ( sss_io(:,:) + 1.e-20 )
(  357) # else
(  358)       dmp(:,:) = 0.e0            ! No internal damping
(  359) # endif
(  360)       
(  361)       !   evaporation damping term ( Surface restoring )
(  362)       zerpplus (:,:) = 0.e0
(  363)       zerpminus(:,:) = 0.e0
(  364)       zplus          =  15. / rday
(  365)       zminus         = -15. / rday
(  366)       
(  367)       DO jj = 1, jpj
(  368)          DO ji = 1, jpi
(  369)             zerp = ( 1. - 2.*upsrnfh(ji,jj) ) * zsrp   &
(  370)                & * ( sss_io(ji,jj) - s_dta(ji,jj,1) )     &
(  371)                & / ( sss_io(ji,jj) + 1.e-20        )
(  372)             erp(ji,jj) = zerp
(  373)             zerpplus (ji,jj) = MAX( erp(ji,jj), 0.e0 )
(  374)             zerpminus(ji,jj) = MIN( erp(ji,jj), 0.e0 )
(  375)          END DO
(  376)       END DO
(  377) 
(  378)       aplus  = 0.e0
(  379)       aminus = 0.e0
(  380)       DO jj = 1, jpj
(  381)          DO ji = 1, jpi
(  382)             zwei   = e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
(  383)             aplus  = aplus  + zerpplus (ji,jj) * zwei
(  384)             aminus = aminus - zerpminus(ji,jj) * zwei
(  385)          END DO
(  386)       END DO
(  387) 
(  388)       IF(ln_ctl .AND. lwp)   WRITE(numout,*) ' oce_sbc_dmp : a+ = ', aplus, ' a- = ', aminus
(  389) #else
(  390)       ! Rigid-lid (emp=emps=E-P-R+Erp)
(  391)       
(  392)       erp(:,:) = ( 1. - freezn(:,:) ) * zsrp    &   ! surface restoring term
(  393)          &     * ( sss_io(:,:) - s_dta(:,:,1) )     &
(  394)          &     / ( sss_io(:,:) + 1.e-20      )






PGF90 (Version     10.5)          09/23/2010  09:50:00      page 8

(  395) #endif
(  396) 
(  397)    END SUBROUTINE oce_sbc_dmp
(  398)    
(  399) #else
(  400)    !!----------------------------------------------------------------------
(  401)    !!   Dummy routine                                      NO salinity data
(  402)    !!----------------------------------------------------------------------
(  403)    USE in_out_manager
(  404)    SUBROUTINE oce_sbc_dmp         ! Dummy routine
(  405)       if(lwp) WRITE(numout,*) 'oce_sbc_dmp: you should not have seen that print! error?'
(  406)    END SUBROUTINE oce_sbc_dmp
(  407) #endif
(  408) #else
(  409)    !!----------------------------------------------------------------------
(  410)    !!   Default option           Dummy module   NO 1D && NO LIM sea-ice model
(  411)    !!----------------------------------------------------------------------
(  412)    USE in_out_manager
(  413) CONTAINS
(  414)    SUBROUTINE ice_stp_1d ( kt )     ! Dummy routine
(  415)       if(lwp) WRITE(numout,*) 'ice_stp_1d: You should not have seen this print! error?', kt
(  416)    END SUBROUTINE ice_stp_1d
(  417)    SUBROUTINE oce_sbc_1d ( kt )     ! Dummy routine
(  418)       if(lwp)WRITE(numout,*) 'oce_sbc_1d: You should not have seen this print! error?', kt
(  419)    END SUBROUTINE oce_sbc_1d
(  420) #endif   
(  421)    !!======================================================================
(  422) END MODULE icestp1d




































PGF90 (Version     10.5)          09/23/2010  09:50:00      page 9

(    1) # 1 "icestp1d.F90"
(    1) MODULE icestp1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE icestp1d   ***
(    4)    !!   Sea-Ice model : 1D LIM Sea ice model time-stepping
(    5)    !!======================================================================
(    6) # 409
(  409)    !!----------------------------------------------------------------------
(  410)    !!   Default option           Dummy module   NO 1D && NO LIM sea-ice model
(  411)    !!----------------------------------------------------------------------
(  412)    USE in_out_manager
(  413) CONTAINS
(  414)    SUBROUTINE ice_stp_1d ( kt )     ! Dummy routine
(  415)       if(lwp) WRITE(numout,*) 'ice_stp_1d: You should not have seen this print! error?', kt
(  416)    END SUBROUTINE ice_stp_1d
(  417)    SUBROUTINE oce_sbc_1d ( kt )     ! Dummy routine
(  418)       if(lwp)WRITE(numout,*) 'oce_sbc_1d: You should not have seen this print! error?', kt
(  419)    END SUBROUTINE oce_sbc_1d
(  420) # 421
(  421)    !!======================================================================
(  422) END MODULE icestp1d
