


PGF90 (Version     10.5)          09/23/2010  09:50:02      page 1

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

Filename: icestp.F90

(    1) MODULE icestp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE icestp   ***
(    4)    !!   Sea-Ice model : LIM Sea ice model time-stepping
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   Lim sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ice_stp       : sea-ice model time-stepping
(   11)    !!----------------------------------------------------------------------
(   12)    USE dom_oce
(   13)    USE oce  ! dynamics and tracers variables
(   14)    USE in_out_manager
(   15)    USE ice_oce         ! ice variables
(   16)    USE flx_oce         ! forcings variables
(   17)    USE dom_ice
(   18)    USE cpl_oce
(   19)    USE daymod
(   20)    USE phycst          ! Define parameters for the routines
(   21)    USE taumod
(   22)    USE ice
(   23)    USE iceini
(   24)    USE ocesbc
(   25)    USE lbclnk
(   26)    USE limdyn
(   27)    USE limtrp
(   28)    USE limthd
(   29)    USE limflx
(   30)    USE limdia
(   31)    USE limwri
(   32)    USE limrst
(   33)    USE limdmp          ! Ice damping
(   34)    USE prtctl          ! Print control
(   35) 
(   36) !!DB
(   37)    USE lib_ncdf        ! netCDF I/O library
(   38)    USE restart
(   39) !DBG
(   40) #if defined key_BGCM_01
(   41) !   USE bgcm_01_model, ONLY : sw_rad
(   42) #endif
(   43) 
(   44) 
(   45)    IMPLICIT NONE
(   46)    PRIVATE






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 2

(   47) 
(   48)    !! * Routine accessibility
(   49)    PUBLIC ice_stp  ! called by step.F90
(   50) 
(   51)    !! * Substitutions
(   52) #  include "domzgr_substitute.h90"
(   53) #  include "vectopt_loop_substitute.h90"
(   54)    !!-----------------------------------------------------
(   55)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   56)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/icestp.F90,v 1.7 2006/03/21 08:42:22 opalod Exp $ 
(   57)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   58)    !!-----------------------------------------------------
(   59) 
(   60) CONTAINS
(   61) 
(   62)   SUBROUTINE ice_stp ( kt )
(   63)     !!---------------------------------------------------------------------
(   64)     !!                  ***  ROUTINE ice_stp  ***
(   65)     !!                   
(   66)     !! ** Purpose :   Louvain la Neuve Sea Ice Model time stepping 
(   67)     !!
(   68)     !! ** Action  : - call the ice dynamics routine 
(   69)     !!              - call the ice advection/diffusion routine 
(   70)     !!              - call the ice thermodynamics routine 
(   71)     !!              - call the routine that computes mass and 
(   72)     !!                heat fluxes at the ice/ocean interface
(   73)     !!              - save the outputs 
(   74)     !!              - save the outputs for restart when necessary
(   75)     !!
(   76)     !! History :
(   77)     !!   1.0  !  99-11  (M. Imbard)  Original code
(   78)     !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(   79)     !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(   80)     !!----------------------------------------------------------------------
(   81)     !! * Arguments
(   82)     INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   83) 
(   84)     !! * Local declarations
(   85)     INTEGER   ::   ji, jj   ! dummy loop indices
(   86) 
(   87)     REAL(wp) , DIMENSION(jpi,jpj)    :: &
(   88)          zsss_io, zsss2_io, zsss3_io          ! tempory workspaces
(   89)     !!----------------------------------------------------------------------
(   90) 
(   91)     IF( kt == nit000 ) THEN
(   92)        IF( lk_cpl ) THEN
(   93)           IF(lwp) WRITE(numout,*)
(   94)           IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)'
(   95)           IF(lwp) WRITE(numout,*) '~~~~~~~   coupled case'
(   96)        ELSE
(   97)           IF(lwp) WRITE(numout,*)
(   98)           IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)' 
(   99)           IF(lwp) WRITE(numout,*) '~~~~~~~   forced case using bulk formulea'
(  100)        ENDIF
(  101)        !  Initialize fluxes fields
(  102)        gtaux(:,:) = 0.e0
(  103)        gtauy(:,:) = 0.e0
(  104)     ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 3

(  105) 
(  106)     ! Temperature , salinity and horizonta wind
(  107)     ! sst_io  and sss_io, u_io and v_io  are initialized at nit000 in limistate.F90 (or limrst.F90) with :
(  108)     !              sst_io = sst_io + (nfice - 1) * (tn(:,:,1)+rt0 )
(  109)     !              sss_io = sss_io + (nfice - 1) * sn(:,:,1)
(  110)     !              u_io  = u_io  + (nfice - 1) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  111)     !              v_io  = v_io  + (nfice - 1) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  112)     !    cumulate fields
(  113)     !
(  114)     sst_io(:,:) = sst_io(:,:) + tn(:,:,1) + rt0
(  115)     sss_io(:,:) = sss_io(:,:) + sn(:,:,1)
(  116) 
(  117) 
(  118)     ! vectors at F-point
(  119)     DO jj = 2, jpj
(  120)        DO ji = fs_2, jpi   ! vector opt.
(  121)           u_io(ji,jj) = u_io(ji,jj) + 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  122)           v_io(ji,jj) = v_io(ji,jj) + 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  123)        END DO
(  124)     END DO
(  125) 
(  126) !!DB: 2009.09.30
(  127) !    IF( MOD( kt-1, nfice ) == 0 ) THEN
(  128)     IF( MOD( kt-nit000, nfice ) == 0 ) THEN
(  129) 
(  130)        ! The LIM model is going to be call
(  131)        sst_io(:,:) = sst_io(:,:) / FLOAT( nfice ) * tmask(:,:,1)
(  132)        sss_io(:,:) = sss_io(:,:) / FLOAT( nfice )
(  133) 
(  134)        ! stress from ocean U- and V-points to ice U,V point
(  135)        DO jj = 2, jpj
(  136)           DO ji = fs_2, jpi   ! vector opt.
(  137)              gtaux(ji,jj) = 0.5 * ( taux(ji-1,jj  ) + taux(ji-1,jj-1) )
(  138)              gtauy(ji,jj) = 0.5 * ( tauy(ji  ,jj-1) + tauy(ji-1,jj-1) )
(  139)              u_io  (ji,jj) = u_io(ji,jj) / FLOAT( nfice )
(  140)              v_io  (ji,jj) = v_io(ji,jj) / FLOAT( nfice )
(  141)           END DO
(  142)        END DO
(  143) 
(  144)        ! lateral boundary condition
(  145)        CALL lbc_lnk( gtaux(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  146)        CALL lbc_lnk( gtauy(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  147)        CALL lbc_lnk( u_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  148)        CALL lbc_lnk( v_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  149) 
(  150)        !!gmbug  in the ocean freezing point computed as :
(  151)        !!gm           fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
(  152)        !!gm                                     - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
(  153)        !!gm           !!                        - 7.53e-4 * pressure
(  154)        !!gm
(  155)        !!gm!bug this is much more accurate and efficient computation
(  156)        !!gm       **************************************************
(  157)        !!gm freezing point from unesco:
(  158)        !!gm     real function tf(s,p)
(  159)        !!gm   function to compute the freezing point of seawater
(  160)        !!gm
(  161)        !!gm   reference: unesco tech. papers in the marine science no. 28. 1978
(  162)        !!gm   eighth report jpots






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 4

(  163)        !!gm   annex 6 freezing point of seawater f.j. millero pp.29-35.
(  164)        !!gm
(  165)        !!gm  units:
(  166)        !!gm         pressure      p          decibars
(  167)        !!gm         salinity      s          pss-78
(  168)        !!gm         temperature   tf         degrees celsius
(  169)        !!gm         freezing pt.
(  170)        !!gm************************************************************
(  171)        !!gm  checkvalue: tf= -2.588567 deg. c for s=40.0, p=500. decibars
(  172)        !!gm     tf=(-.0575+1.710523e-3*sqrt(abs(s))-2.154996e-4*s)*s-7.53e-4*p
(  173)        !!gm     return
(  174)        !!gm     end
(  175)        !!gm!bug
(  176) 
(  177) 
(  178)        !!gm      DO jj = 1, jpj
(  179)        !!gm         DO ji = 1, jpi
(  180)        !!gm            tfu(ji,jj)  = (  rt0 + ( - 0.0575                              &
(  181)        !!gm               &                     + 1.710523e-3 * SQRT( sss_io(ji,jj) )   &
(  182)        !!gm               &                     - 2.154996e-4 *       sss_io(ji,jj)   ) * sss_io(ji,jj)  ) * tms(ji,jj)
(  183)        !!gm         END DO
(  184)        !!gm      END DO
(  185)        !!gm
(  186)        zsss_io (:,:) = SQRT( sss_io(:,:) ) 
(  187)        zsss2_io(:,:) =  sss_io(:,:) *  sss_io(:,:)
(  188)        zsss3_io(:,:) = zsss_io(:,:) * zsss_io(:,:) * zsss_io(:,:)
(  189) 
(  190)        DO jj = 1, jpj
(  191)           DO ji = 1, jpi
(  192)              tfu(ji,jj)  = ABS ( rt0 - 0.0575       *   sss_io(ji,jj)   &
(  193)                   &                    + 1.710523e-03 * zsss3_io(ji,jj)   &
(  194)                   &                    - 2.154996e-04 * zsss2_io(ji,jj) ) * tms(ji,jj)
(  195)           END DO
(  196)        END DO
(  197) 
(  198) 
(  199) 
(  200)        IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  201)           CALL prt_ctl_info('Ice Forcings ')
(  202)           CALL prt_ctl(tab2d_1=qsr_oce ,clinfo1=' qsr_oce  : ', tab2d_2=qsr_ice , clinfo2=' qsr_ice   : ')
(  203)           CALL prt_ctl(tab2d_1=qnsr_oce,clinfo1=' qnsr_oce : ', tab2d_2=qnsr_ice, clinfo2=' qnsr_ice  : ')
(  204)           CALL prt_ctl(tab2d_1=evap    ,clinfo1=' evap     : ')
(  205)           CALL prt_ctl(tab2d_1=tprecip ,clinfo1=' precip   : ', tab2d_2=sprecip , clinfo2=' Snow      : ')
(  206)           CALL prt_ctl(tab2d_1=gtaux   ,clinfo1=' u-stress : ', tab2d_2=gtauy   , clinfo2=' v-stress  : ')
(  207)           CALL prt_ctl(tab2d_1=sst_io  ,clinfo1=' sst      : ', tab2d_2=sss_io  , clinfo2=' sss       : ')
(  208)           CALL prt_ctl(tab2d_1=u_io    ,clinfo1=' u_io     : ', tab2d_2=v_io    , clinfo2=' v_io      : ')
(  209)           CALL prt_ctl(tab2d_1=hsnif   ,clinfo1=' hsnif  1 : ', tab2d_2=hicif   , clinfo2=' hicif     : ')
(  210)           CALL prt_ctl(tab2d_1=frld    ,clinfo1=' frld   1 : ', tab2d_2=sist    , clinfo2=' sist      : ')
(  211)        ENDIF
(  212) 
(  213)        ! Ice model call
(  214)        numit = numit + nfice 
(  215) 
(  216)        !                                                           !--------------!
(  217)        CALL lim_dyn                                                ! Ice dynamics !   ( rheology/dynamics )
(  218)        !                                                           !--------------!
(  219)        IF(ln_ctl) THEN
(  220)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  2 : ', tab2d_2=hicif , clinfo2=' hicif     : ')






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 5

(  221)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   2 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  222)        ENDIF
(  223) 
(  224) 
(  225)        !                                                           !---------------!
(  226)        CALL lim_trp                                                ! Ice transport !  ( Advection/diffusion )
(  227)        !                                                           !---------------!
(  228)        IF(ln_ctl) THEN
(  229)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  3 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  230)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   3 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  231)        ENDIF
(  232) 
(  233)        !                                                           !-------------!
(  234)        IF( ln_limdmp ) CALL lim_dmp(kt)                            ! Ice damping !
(  235)        !                                                           !-------------!
(  236) 
(  237)        !                                                           !--------------------!
(  238)        CALL lim_thd                                                ! Ice thermodynamics !
(  239)        !                                                           !--------------------!
(  240)        IF(ln_ctl) THEN
(  241)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  4 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  242)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   4 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  243)        ENDIF
(  244) 
(  245) 
(  246)        ! Mass and heat fluxes from ice to ocean
(  247)        !                                                           !------------------------------!
(  248)        CALL lim_flx                                                ! Ice/Ocean Mass & Heat fluxes !
(  249)        !                                                           !------------------------------!
(  250) 
(  251)        IF( MOD( numit, ninfo ) == 0 .OR. ntmoy == 1 )  THEN        !-----------------!
(  252)           CALL lim_dia                                             ! Ice Diagnostics !
(  253)        ENDIF                                                       !-----------------!
(  254) 
(  255)        !                                                           !-------------!
(  256)        CALL lim_wri                                                ! Ice outputs !
(  257)        !                                                           !-------------!
(  258) !!DB: 2009.09.30: Moved outside of if-block
(  259) !       IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
(  260) !          !                                                        !------------------!
(  261) !          !!DB: old IOIPSL call
(  262) !          !            CALL lim_rst_write( numit )                              ! Ice restart file !
(  263) !          !                                                        !------------------!
(  264) !          !!DB
(  265) !          call rst_ice_write(numit)
(  266) !       ENDIF
(  267) 
(  268) 
(  269) #if defined key_passivetrc || defined key_BGCM_01
(  270)        !!DB 2009.08.17 -- need sw_rad variable for BGCM. For most situations this is 
(  271)        !! the variable qsr_oce, which is zeroed below ====> save a copy of it
(  272) !!DBG
(  273) !       sw_rad(:,:) = qsr_oce(:,:)
(  274) #endif
(  275) 
(  276)        ! Re-initialization of forcings
(  277)        qsr_oce (:,:) = 0.e0
(  278)        qsr_ice (:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 6

(  279)        qnsr_oce(:,:) = 0.e0
(  280)        qnsr_ice(:,:) = 0.e0 
(  281)        dqns_ice(:,:) = 0.e0 
(  282)        tprecip (:,:) = 0.e0 
(  283)        sprecip (:,:) = 0.e0
(  284)        fr1_i0  (:,:) = 0.e0
(  285)        fr2_i0  (:,:) = 0.e0
(  286)        evap    (:,:) = 0.e0
(  287) #if defined key_coupled 
(  288)        rrunoff (:,:) = 0.e0
(  289)        calving (:,:) = 0.e0
(  290) #else
(  291)        qla_ice (:,:) = 0.e0
(  292)        dqla_ice(:,:) = 0.e0
(  293) #endif
(  294) 
(  295)     ENDIF   !!  IF( MOD( kt-1, nfice ) == 0 ) 
(  296) 
(  297) !!DB: 2009.09.30: Force ice restarts at same kt as restart()
(  298)     if( mod( kt, nstock ) == 0 .OR. kt == nitend ) THEN
(  299)        call rst_ice_write(kt)
(  300)     endif
(  301) 
(  302) 
(  303) 
(  304)   END SUBROUTINE ice_stp
(  305) 
(  306) #else
(  307)    !!----------------------------------------------------------------------
(  308)    !!   Default option           Dummy module          NO LIM sea-ice model
(  309)    !!----------------------------------------------------------------------
(  310)    USE in_out_manager
(  311) CONTAINS
(  312)    SUBROUTINE ice_stp ( kt )     ! Dummy routine
(  313)       if(lpw)WRITE(numout,*) 'ice_stp: You should not have seen this print! error?', kt
(  314)    END SUBROUTINE ice_stp
(  315) #endif
(  316) 
(  317)    !!======================================================================
(  318) END MODULE icestp
























PGF90 (Version     10.5)          09/23/2010  09:50:02      page 7

(    1) # 1 "icestp.F90"
(    1) MODULE icestp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE icestp   ***
(    4)    !!   Sea-Ice model : LIM Sea ice model time-stepping
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   Lim sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ice_stp       : sea-ice model time-stepping
(   11)    !!----------------------------------------------------------------------
(   12)    USE dom_oce
(   13)    USE oce  ! dynamics and tracers variables
(   14)    USE in_out_manager
(   15)    USE ice_oce         ! ice variables
(   16)    USE flx_oce         ! forcings variables
(   17)    USE dom_ice
(   18)    USE cpl_oce
(   19)    USE daymod
(   20)    USE phycst          ! Define parameters for the routines
(   21)    USE taumod
(   22)    USE ice
(   23)    USE iceini
(   24)    USE ocesbc
(   25)    USE lbclnk
(   26)    USE limdyn
(   27)    USE limtrp
(   28)    USE limthd
(   29)    USE limflx
(   30)    USE limdia
(   31)    USE limwri
(   32)    USE limrst
(   33)    USE limdmp          ! Ice damping
(   34)    USE prtctl          ! Print control
(   35) 
(   36) !!DB
(   37)    USE lib_ncdf        ! netCDF I/O library
(   38)    USE restart
(   39) !DBG
(   40) 
(   41) 
(   42) # 45
(   45)    IMPLICIT NONE
(   46)    PRIVATE
(   47) 
(   48)    !! * Routine accessibility
(   49)    PUBLIC ice_stp  ! called by step.F90
(   50) 
(   51)    !! * Substitutions
(   52) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 8

(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 53 "icestp.F90"
(   53) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 54 "icestp.F90"
(   54) # 54
(   54)    !!-----------------------------------------------------
(   55)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   56)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/icestp.F90,v 1.7 2006/03/21 08:42:22 opalod Exp $ 
(   57)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   58)    !!-----------------------------------------------------
(   59) 
(   60) CONTAINS
(   61) 
(   62)   SUBROUTINE ice_stp ( kt )
(   63)     !!---------------------------------------------------------------------
(   64)     !!                  ***  ROUTINE ice_stp  ***
(   65)     !!                   
(   66)     !! ** Purpose :   Louvain la Neuve Sea Ice Model time stepping 
(   67)     !!
(   68)     !! ** Action  : - call the ice dynamics routine 
(   69)     !!              - call the ice advection/diffusion routine 
(   70)     !!              - call the ice thermodynamics routine 
(   71)     !!              - call the routine that computes mass and 
(   72)     !!                heat fluxes at the ice/ocean interface
(   73)     !!              - save the outputs 
(   74)     !!              - save the outputs for restart when necessary
(   75)     !!
(   76)     !! History :
(   77)     !!   1.0  !  99-11  (M. Imbard)  Original code
(   78)     !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
(   79)     !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
(   80)     !!----------------------------------------------------------------------
(   81)     !! * Arguments
(   82)     INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   83) 
(   84)     !! * Local declarations
(   85)     INTEGER   ::   ji, jj   ! dummy loop indices
(   86) 
(   87)     REAL(wp) , DIMENSION(jpi,jpj)    :: &
(   88)          zsss_io, zsss2_io, zsss3_io          ! tempory workspaces
(   89)     !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 9

(   90) 
(   91)     IF( kt == nit000 ) THEN
(   92)        IF( lk_cpl ) THEN
(   93)           IF(lwp) WRITE(numout,*)
(   94)           IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)'
(   95)           IF(lwp) WRITE(numout,*) '~~~~~~~   coupled case'
(   96)        ELSE
(   97)           IF(lwp) WRITE(numout,*)
(   98)           IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)' 
(   99)           IF(lwp) WRITE(numout,*) '~~~~~~~   forced case using bulk formulea'
(  100)        ENDIF
(  101)        !  Initialize fluxes fields
(  102)        gtaux(:,:) = 0.e0
(  103)        gtauy(:,:) = 0.e0
(  104)     ENDIF
(  105) 
(  106)     ! Temperature , salinity and horizonta wind
(  107)     ! sst_io  and sss_io, u_io and v_io  are initialized at nit000 in limistate.F90 (or limrst.F90) with :
(  108)     !              sst_io = sst_io + (nfice - 1) * (tn(:,:,1)+rt0 )
(  109)     !              sss_io = sss_io + (nfice - 1) * sn(:,:,1)
(  110)     !              u_io  = u_io  + (nfice - 1) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  111)     !              v_io  = v_io  + (nfice - 1) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  112)     !    cumulate fields
(  113)     !
(  114)     sst_io(:,:) = sst_io(:,:) + tn(:,:,1) + rt0
(  115)     sss_io(:,:) = sss_io(:,:) + sn(:,:,1)
(  116) 
(  117) 
(  118)     ! vectors at F-point
(  119)     DO jj = 2, jpj
(  120)        DO ji = 2, jpi   ! vector opt.
(  121)           u_io(ji,jj) = u_io(ji,jj) + 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
(  122)           v_io(ji,jj) = v_io(ji,jj) + 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
(  123)        END DO
(  124)     END DO
(  125) 
(  126) !!DB: 2009.09.30
(  127) !    IF( MOD( kt-1, nfice ) == 0 ) THEN
(  128)     IF( MOD( kt-nit000, nfice ) == 0 ) THEN
(  129) 
(  130)        ! The LIM model is going to be call
(  131)        sst_io(:,:) = sst_io(:,:) / FLOAT( nfice ) * tmask(:,:,1)
(  132)        sss_io(:,:) = sss_io(:,:) / FLOAT( nfice )
(  133) 
(  134)        ! stress from ocean U- and V-points to ice U,V point
(  135)        DO jj = 2, jpj
(  136)           DO ji = 2, jpi   ! vector opt.
(  137)              gtaux(ji,jj) = 0.5 * ( taux(ji-1,jj  ) + taux(ji-1,jj-1) )
(  138)              gtauy(ji,jj) = 0.5 * ( tauy(ji  ,jj-1) + tauy(ji-1,jj-1) )
(  139)              u_io  (ji,jj) = u_io(ji,jj) / FLOAT( nfice )
(  140)              v_io  (ji,jj) = v_io(ji,jj) / FLOAT( nfice )
(  141)           END DO
(  142)        END DO
(  143) 
(  144)        ! lateral boundary condition
(  145)        CALL lbc_lnk( gtaux(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  146)        CALL lbc_lnk( gtauy(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  147)        CALL lbc_lnk( u_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 10

(  148)        CALL lbc_lnk( v_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
(  149) 
(  150)        !!gmbug  in the ocean freezing point computed as :
(  151)        !!gm           fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
(  152)        !!gm                                     - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
(  153)        !!gm           !!                        - 7.53e-4 * pressure
(  154)        !!gm
(  155)        !!gm!bug this is much more accurate and efficient computation
(  156)        !!gm       **************************************************
(  157)        !!gm freezing point from unesco:
(  158)        !!gm     real function tf(s,p)
(  159)        !!gm   function to compute the freezing point of seawater
(  160)        !!gm
(  161)        !!gm   reference: unesco tech. papers in the marine science no. 28. 1978
(  162)        !!gm   eighth report jpots
(  163)        !!gm   annex 6 freezing point of seawater f.j. millero pp.29-35.
(  164)        !!gm
(  165)        !!gm  units:
(  166)        !!gm         pressure      p          decibars
(  167)        !!gm         salinity      s          pss-78
(  168)        !!gm         temperature   tf         degrees celsius
(  169)        !!gm         freezing pt.
(  170)        !!gm************************************************************
(  171)        !!gm  checkvalue: tf= -2.588567 deg. c for s=40.0, p=500. decibars
(  172)        !!gm     tf=(-.0575+1.710523e-3*sqrt(abs(s))-2.154996e-4*s)*s-7.53e-4*p
(  173)        !!gm     return
(  174)        !!gm     end
(  175)        !!gm!bug
(  176) 
(  177) 
(  178)        !!gm      DO jj = 1, jpj
(  179)        !!gm         DO ji = 1, jpi
(  180)        !!gm            tfu(ji,jj)  = (  rt0 + ( - 0.0575                              &
(  181)        !!gm               &                     + 1.710523e-3 * SQRT( sss_io(ji,jj) )   &
(  182)        !!gm               &                     - 2.154996e-4 *       sss_io(ji,jj)   ) * sss_io(ji,jj)  ) * tms(ji,jj)
(  183)        !!gm         END DO
(  184)        !!gm      END DO
(  185)        !!gm
(  186)        zsss_io (:,:) = SQRT( sss_io(:,:) ) 
(  187)        zsss2_io(:,:) =  sss_io(:,:) *  sss_io(:,:)
(  188)        zsss3_io(:,:) = zsss_io(:,:) * zsss_io(:,:) * zsss_io(:,:)
(  189) 
(  190)        DO jj = 1, jpj
(  191)           DO ji = 1, jpi
(  192)              tfu(ji,jj)  = ABS ( rt0 - 0.0575       *   sss_io(ji,jj)   &
(  193)                   &                    + 1.710523e-03 * zsss3_io(ji,jj)   &
(  194)                   &                    - 2.154996e-04 * zsss2_io(ji,jj) ) * tms(ji,jj)
(  195)           END DO
(  196)        END DO
(  197) 
(  198) 
(  199) 
(  200)        IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  201)           CALL prt_ctl_info('Ice Forcings ')
(  202)           CALL prt_ctl(tab2d_1=qsr_oce ,clinfo1=' qsr_oce  : ', tab2d_2=qsr_ice , clinfo2=' qsr_ice   : ')
(  203)           CALL prt_ctl(tab2d_1=qnsr_oce,clinfo1=' qnsr_oce : ', tab2d_2=qnsr_ice, clinfo2=' qnsr_ice  : ')
(  204)           CALL prt_ctl(tab2d_1=evap    ,clinfo1=' evap     : ')
(  205)           CALL prt_ctl(tab2d_1=tprecip ,clinfo1=' precip   : ', tab2d_2=sprecip , clinfo2=' Snow      : ')






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 11

(  206)           CALL prt_ctl(tab2d_1=gtaux   ,clinfo1=' u-stress : ', tab2d_2=gtauy   , clinfo2=' v-stress  : ')
(  207)           CALL prt_ctl(tab2d_1=sst_io  ,clinfo1=' sst      : ', tab2d_2=sss_io  , clinfo2=' sss       : ')
(  208)           CALL prt_ctl(tab2d_1=u_io    ,clinfo1=' u_io     : ', tab2d_2=v_io    , clinfo2=' v_io      : ')
(  209)           CALL prt_ctl(tab2d_1=hsnif   ,clinfo1=' hsnif  1 : ', tab2d_2=hicif   , clinfo2=' hicif     : ')
(  210)           CALL prt_ctl(tab2d_1=frld    ,clinfo1=' frld   1 : ', tab2d_2=sist    , clinfo2=' sist      : ')
(  211)        ENDIF
(  212) 
(  213)        ! Ice model call
(  214)        numit = numit + nfice 
(  215) 
(  216)        !                                                           !--------------!
(  217)        CALL lim_dyn                                                ! Ice dynamics !   ( rheology/dynamics )
(  218)        !                                                           !--------------!
(  219)        IF(ln_ctl) THEN
(  220)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  2 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  221)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   2 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  222)        ENDIF
(  223) 
(  224) 
(  225)        !                                                           !---------------!
(  226)        CALL lim_trp                                                ! Ice transport !  ( Advection/diffusion )
(  227)        !                                                           !---------------!
(  228)        IF(ln_ctl) THEN
(  229)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  3 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  230)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   3 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  231)        ENDIF
(  232) 
(  233)        !                                                           !-------------!
(  234)        IF( ln_limdmp ) CALL lim_dmp(kt)                            ! Ice damping !
(  235)        !                                                           !-------------!
(  236) 
(  237)        !                                                           !--------------------!
(  238)        CALL lim_thd                                                ! Ice thermodynamics !
(  239)        !                                                           !--------------------!
(  240)        IF(ln_ctl) THEN
(  241)           CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  4 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
(  242)           CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   4 : ', tab2d_2=sist  , clinfo2=' sist      : ')
(  243)        ENDIF
(  244) 
(  245) 
(  246)        ! Mass and heat fluxes from ice to ocean
(  247)        !                                                           !------------------------------!
(  248)        CALL lim_flx                                                ! Ice/Ocean Mass & Heat fluxes !
(  249)        !                                                           !------------------------------!
(  250) 
(  251)        IF( MOD( numit, ninfo ) == 0 .OR. ntmoy == 1 )  THEN        !-----------------!
(  252)           CALL lim_dia                                             ! Ice Diagnostics !
(  253)        ENDIF                                                       !-----------------!
(  254) 
(  255)        !                                                           !-------------!
(  256)        CALL lim_wri                                                ! Ice outputs !
(  257)        !                                                           !-------------!
(  258) !!DB: 2009.09.30: Moved outside of if-block
(  259) !       IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
(  260) !          !                                                        !------------------!
(  261) !          !!DB: old IOIPSL call
(  262) !          !            CALL lim_rst_write( numit )                              ! Ice restart file !
(  263) !          !                                                        !------------------!






PGF90 (Version     10.5)          09/23/2010  09:50:02      page 12

(  264) !          !!DB
(  265) !          call rst_ice_write(numit)
(  266) !       ENDIF
(  267) 
(  268) 
(  269) 
(  270) # 276
(  276)        ! Re-initialization of forcings
(  277)        qsr_oce (:,:) = 0.e0
(  278)        qsr_ice (:,:) = 0.e0
(  279)        qnsr_oce(:,:) = 0.e0
(  280)        qnsr_ice(:,:) = 0.e0 
(  281)        dqns_ice(:,:) = 0.e0 
(  282)        tprecip (:,:) = 0.e0 
(  283)        sprecip (:,:) = 0.e0
(  284)        fr1_i0  (:,:) = 0.e0
(  285)        fr2_i0  (:,:) = 0.e0
(  286)        evap    (:,:) = 0.e0
(  287) # 291
(  291)        qla_ice (:,:) = 0.e0
(  292)        dqla_ice(:,:) = 0.e0
(  293) 
(  294) # 295
(  295)     ENDIF   !!  IF( MOD( kt-1, nfice ) == 0 ) 
(  296) 
(  297) !!DB: 2009.09.30: Force ice restarts at same kt as restart()
(  298)     if( mod( kt, nstock ) == 0 .OR. kt == nitend ) THEN
(  299)        call rst_ice_write(kt)
(  300)     endif
(  301) 
(  302) 
(  303) 
(  304)   END SUBROUTINE ice_stp
(  305) 
(  306) 
(  307) # 317
(  317)    !!======================================================================
(  318) END MODULE icestp
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (icestp.F90: 131)
