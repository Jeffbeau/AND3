


PGF90 (Version     10.5)          09/23/2010  09:50:03      page 1

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

Filename: istate.F90

(    1) MODULE istate
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  istate  ***
(    4)    !! Ocean state   :  initial state setting
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   istate_init   : initial state setting
(    9)    !!   istate_tem    : analytical profile for initial Temperature
(   10)    !!   istate_sal    : analytical profile for initial Salinity
(   11)    !!   istate_eel    : initial state setting of EEL R5 configuration
(   12)    !!   istate_gyre   : initial state setting of GYRE configuration
(   13)    !!   istate_uvg    : initial velocity in geostropic balance
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and active tracers 
(   17)    USE dom_oce         ! ocean space and time domain 
(   18)    USE daymod          ! 
(   19)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   20)    USE zdf_oce         ! ocean vertical physics
(   21)    USE in_out_manager  ! I/O manager
(   22)    USE phycst          ! physical constants
(   23)    USE wzvmod          ! verctical velocity               (wzv     routine)
(   24)    USE dtatem          ! temperature data                 (dta_tem routine)
(   25)    USE dtasal          ! salinity data                    (dta_sal routine)
(   26)    USE restart         ! ocean restart                   (rst_read routine)
(   27)    USE solisl          ! ???
(   28) 
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Routine accessibility
(   33)    PUBLIC istate_init   ! routine called by step.F90
(   34) 
(   35)    !! * Substitutions
(   36) #  include "domzgr_substitute.h90"
(   37) #  include "vectopt_loop_substitute.h90"
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/istate.F90,v 1.11 2006/04/10 15:46:04 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE istate_init






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !!                   ***  ROUTINE istate_init  ***
(   49)       !! 
(   50)       !! ** Purpose :   Initialization of the dynamics and tracers.
(   51)       !!
(   52)       !! ** Method  :
(   53)       !!
(   54)       !! History :
(   55)       !!   4.0  !  91-03  ()  Original code
(   56)       !!        !  91-11  (G. Madec)
(   57)       !!   9.0  !  03-09  (G. Madec)  F90: Free form, modules, orthogonality
(   58)       !!----------------------------------------------------------------------
(   59)       !! * Local declarations
(   60)       !!----------------------------------------------------------------------
(   61) 
(   62) 
(   63)       ! Initialization to zero
(   64)       ! ----------------------
(   65) 
(   66)       !     before fields       !       now fields        !      after fields       !
(   67)      ub   (:,:,:) = 0.e0   ;   un   (:,:,:) = 0.e0   ;   ua   (:,:,:) = 0.e0
(   68)      vb   (:,:,:) = 0.e0   ;   vn   (:,:,:) = 0.e0   ;   va   (:,:,:) = 0.e0
(   69)      wn   (:,:,:) = 0.e0   ;
(   70)      rotb (:,:,:) = 0.e0   ;   rotn (:,:,:) = 0.e0   ;
(   71)      hdivb(:,:,:) = 0.e0   ;   hdivn(:,:,:) = 0.e0   ;
(   72)      
(   73)      tb   (:,:,:) = 0.e0   ;   tn   (:,:,:) = 0.e0   ;   ta   (:,:,:) = 0.e0
(   74)      sb   (:,:,:) = 0.e0   ;   sn   (:,:,:) = 0.e0   ;   sa   (:,:,:) = 0.e0
(   75)      
(   76)      rhd  (:,:,:) = 0.e0
(   77)      rhop (:,:,:) = 0.e0
(   78)      rn2  (:,:,:) = 0.e0 
(   79)      
(   80) #if defined key_dynspg_rl
(   81)       ! rigid-lid formulation
(   82)       bsfb(:,:) = 0.e0      ! before barotropic stream-function
(   83)       bsfn(:,:) = 0.e0      ! now    barotropic stream-function
(   84)       bsfd(:,:) = 0.e0      ! barotropic stream-function trend
(   85) #endif
(   86)       ! free surface formulation
(   87)       sshb(:,:) = 0.e0      ! before sea-surface height
(   88)       sshn(:,:) = 0.e0      ! now    sea-surface height
(   89) 
(   90) 
(   91)       IF( ln_rstart ) THEN                    ! Restart from a file
(   92)          !                                    ! -------------------
(   93)          neuler = 1                              ! Set time-step indicator at nit000 (leap-frog)
(   94)          CALL rst_read                           ! Read the restart file
(   95)       ELSE
(   96)          !                                    ! Start from rest
(   97)          !                                    ! ---------------
(   98)          neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
(   99)          adatrj = 0._wp
(  100) !!DB: OLD
(  101) !         IF( cp_cfg == 'eel' ) THEN
(  102) !            CALL istate_eel                      ! EEL   configuration : start from pre-defined
(  103) !            !                                    !                       velocity and thermohaline fields
(  104) !         ELSEIF( cp_cfg == 'gyre' ) THEN         






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 3

(  105) !            CALL istate_gyre                     ! GYRE  configuration : start from pre-defined temperature
(  106) !            !                                    !                       and salinity fields 
(  107) !         ELSE
(  108)          !                                       ! Other configurations: Initial temperature and salinity fields
(  109) !#if defined key_dtatem
(  110)          CALL dta_tem( nit000 )                  ! read 3D temperature data
(  111)          tb(:,:,:) = t_dta(:,:,:)                ! use temperature data read
(  112)          tn(:,:,:) = t_dta(:,:,:)
(  113) !#else
(  114) !         IF(lwp) WRITE(numout,*)                 ! analytical temperature profile
(  115) !         IF(lwp) WRITE(numout,*)' Temperature initialization using an analytic profile'
(  116) !         CALL istate_tem
(  117) !#endif
(  118) !#if defined key_dtasal
(  119)          CALL dta_sal( nit000 )                  ! read 3D salinity data
(  120)          sb(:,:,:) = s_dta(:,:,:)                ! use salinity data read
(  121)          sn(:,:,:) = s_dta(:,:,:)
(  122) !#else
(  123) !         ! No salinity data
(  124) !         IF(lwp)WRITE(numout,*)                  ! analytical salinity profile
(  125) !         IF(lwp)WRITE(numout,*)' Salinity initialisation using a constant value'
(  126) !         CALL istate_sal
(  127) !#endif
(  128)          !         ENDIF
(  129) 
(  130)       ENDIF
(  131)       !                                       ! Vertical velocity
(  132)       !                                       ! -----------------
(  133)       CALL wzv( nit000 )                         ! from horizontal divergence
(  134) 
(  135)    END SUBROUTINE istate_init
(  136) 
(  137) 
(  138) 
(  139) 
(  140) 
(  141)    SUBROUTINE istate_uvg
(  142)       !!----------------------------------------------------------------------
(  143)       !!                  ***  ROUTINE istate_uvg  ***
(  144)       !!
(  145)       !! ** Purpose :   Compute the geostrophic velocities from (tn,sn) fields
(  146)       !!
(  147)       !! ** Method  :   Using the hydrostatic hypothesis the now hydrostatic 
(  148)       !!      pressure is computed by integrating the in-situ density from the
(  149)       !!      surface to the bottom.
(  150)       !!                 p=integral [ rau*g dz ]
(  151)       !!
(  152)       !! History :
(  153)       !!   8.1  !  01-09  (M. Levy, M. Ben Jelloul)  Original code
(  154)       !!   8.5  !  02-09  (G. Madec)  F90: Free form
(  155)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  156)       !!----------------------------------------------------------------------
(  157)       !! * Modules used
(  158)       USE eosbn2          ! eq. of state, Brunt Vaisala frequency (eos     routine)
(  159)       USE dynspg          ! surface pressure gradient             (dyn_spg routine)
(  160)       USE divcur          ! hor. divergence & rel. vorticity      (div_cur routine)
(  161)       USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 4

(  163)       !! * Local declarations
(  164)       INTEGER ::   ji, jj, jk        ! dummy loop indices
(  165)       INTEGER ::   indic             ! ???
(  166)       REAL(wp) ::   &
(  167)          zmsv, zphv, zmsu, zphu,  &  ! temporary scalars
(  168)          zalfg
(  169)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(  170)          zprn                        ! workspace
(  171)       !!----------------------------------------------------------------------
(  172) 
(  173)       IF(lwp) WRITE(numout,*) 
(  174)       IF(lwp) WRITE(numout,*) 'istate_uvg : Start from Geostrophy'
(  175)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
(  176) 
(  177)       ! Compute the now hydrostatic pressure
(  178)       ! ------------------------------------
(  179) 
(  180)       zalfg = 0.5 * grav * rau0
(  181)       ! Surface value
(  182)       zprn(:,:,1) = zalfg * fse3w(:,:,1) * ( 1 + rhd(:,:,1) )
(  183) 
(  184)       ! Vertical integration from the surface
(  185)       DO jk = 2, jpkm1
(  186)          zprn(:,:,jk) = zprn(:,:,jk-1)   &
(  187)             &         + zalfg * fse3w(:,:,jk) * ( 2. + rhd(:,:,jk) + rhd(:,:,jk-1) )
(  188)       END DO  
(  189) 
(  190)       ! Compute geostrophic balance
(  191)       ! ---------------------------
(  192) 
(  193)       DO jk = 1, jpkm1
(  194)          DO jj = 2, jpjm1
(  195)             DO ji = fs_2, fs_jpim1   ! vertor opt.
(  196)                zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
(  197)                                + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1.  )
(  198)                zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
(  199)                     + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
(  200)                     + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
(  201)                     + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
(  202)                zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)
(  203) 
(  204)                zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
(  205)                                + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1.  )
(  206)                zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
(  207)                     + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
(  208)                     + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &
(  209)                     + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
(  210)                zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)
(  211) 
(  212)                ! Compute the geostrophic velocities
(  213)                un(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
(  214)                vn(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
(  215)             END DO
(  216)          END DO
(  217)       END DO
(  218) 
(  219)       IF(lwp) WRITE(numout,*) '         we force to zero bottom velocity'
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 5

(  221)       ! Susbtract the bottom velocity (level jpk-1 for flat bottom case)
(  222)       ! to have a zero bottom velocity
(  223) 
(  224)       DO jk = 1, jpkm1
(  225)          un(:,:,jk) = ( un(:,:,jk) - un(:,:,jpkm1) ) * umask(:,:,jk)
(  226)          vn(:,:,jk) = ( vn(:,:,jk) - vn(:,:,jpkm1) ) * vmask(:,:,jk)
(  227)       END DO
(  228) 
(  229)       CALL lbc_lnk( un, 'U', -1. )
(  230)       CALL lbc_lnk( vn, 'V', -1. )
(  231)       
(  232)       ub(:,:,:) = un(:,:,:)
(  233)       vb(:,:,:) = vn(:,:,:)
(  234)       
(  235)       ! WARNING !!!!!
(  236)       ! after initializing u and v, we need to calculate the initial streamfunction bsf.
(  237)       ! Otherwise, only the trend will be computed and the model will blow up (inconsistency).
(  238)       
(  239)       ! to do that, we call dyn_spg with a special trick:
(  240)       ! we fill ua and va with the velocities divided by dt,
(  241)       ! and the streamfunction will be brought to the right
(  242)       ! value assuming the velocities have been set up in
(  243)       ! one time step.
(  244)       ! we then set bsfd to zero (first guess for next step
(  245)       ! is d(psi)/dt = 0.)
(  246) 
(  247)       !  sets up s false trend to calculate the barotropic
(  248)       !  streamfunction.
(  249) 
(  250)       ua(:,:,:) = ub(:,:,:) / rdt
(  251)       va(:,:,:) = vb(:,:,:) / rdt
(  252) 
(  253)       ! calls dyn_spg. we assume euler time step, starting from rest.
(  254)       indic = 0
(  255)       CALL dyn_spg( nit000, indic )       ! surface pressure gradient
(  256) 
(  257)       ! the new velocity is ua*rdt
(  258) 
(  259)       CALL lbc_lnk( ua, 'U', -1. )
(  260)       CALL lbc_lnk( va, 'V', -1. )
(  261) 
(  262)       ub(:,:,:) = ua(:,:,:) * rdt
(  263)       vb(:,:,:) = va(:,:,:) * rdt
(  264)       ua(:,:,:) = 0.e0
(  265)       va(:,:,:) = 0.e0
(  266)       un(:,:,:) = ub(:,:,:)
(  267)       vn(:,:,:) = vb(:,:,:)
(  268)        
(  269) #if defined key_dynspg_rl
(  270)       IF( lk_isl )   bsfb(:,:) = bsfn(:,:)          ! Put bsfb to zero
(  271) #endif
(  272) 
(  273)       ! Compute the divergence and curl
(  274) 
(  275)       CALL div_cur( nit000 )            ! now horizontal divergence and curl
(  276) 
(  277)       hdivb(:,:,:) = hdivn(:,:,:)       ! set the before to the now value
(  278)       rotb (:,:,:) = rotn (:,:,:)       ! set the before to the now value






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 6

(  279) 
(  280)    END SUBROUTINE istate_uvg
(  281) 
(  282)    !!=====================================================================
(  283) END MODULE istate



























































PGF90 (Version     10.5)          09/23/2010  09:50:03      page 7

(    1) # 1 "istate.F90"
(    1) MODULE istate
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  istate  ***
(    4)    !! Ocean state   :  initial state setting
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   istate_init   : initial state setting
(    9)    !!   istate_tem    : analytical profile for initial Temperature
(   10)    !!   istate_sal    : analytical profile for initial Salinity
(   11)    !!   istate_eel    : initial state setting of EEL R5 configuration
(   12)    !!   istate_gyre   : initial state setting of GYRE configuration
(   13)    !!   istate_uvg    : initial velocity in geostropic balance
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and active tracers 
(   17)    USE dom_oce         ! ocean space and time domain 
(   18)    USE daymod          ! 
(   19)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   20)    USE zdf_oce         ! ocean vertical physics
(   21)    USE in_out_manager  ! I/O manager
(   22)    USE phycst          ! physical constants
(   23)    USE wzvmod          ! verctical velocity               (wzv     routine)
(   24)    USE dtatem          ! temperature data                 (dta_tem routine)
(   25)    USE dtasal          ! salinity data                    (dta_sal routine)
(   26)    USE restart         ! ocean restart                   (rst_read routine)
(   27)    USE solisl          ! ???
(   28) 
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Routine accessibility
(   33)    PUBLIC istate_init   ! routine called by step.F90
(   34) 
(   35)    !! * Substitutions
(   36) # 1 "./domzgr_substitute.h90"
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
(   48) # 37 "istate.F90"
(   37) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 8

(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 38 "istate.F90"
(   38) # 38
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/istate.F90,v 1.11 2006/04/10 15:46:04 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE istate_init
(   47)       !!----------------------------------------------------------------------
(   48)       !!                   ***  ROUTINE istate_init  ***
(   49)       !! 
(   50)       !! ** Purpose :   Initialization of the dynamics and tracers.
(   51)       !!
(   52)       !! ** Method  :
(   53)       !!
(   54)       !! History :
(   55)       !!   4.0  !  91-03  ()  Original code
(   56)       !!        !  91-11  (G. Madec)
(   57)       !!   9.0  !  03-09  (G. Madec)  F90: Free form, modules, orthogonality
(   58)       !!----------------------------------------------------------------------
(   59)       !! * Local declarations
(   60)       !!----------------------------------------------------------------------
(   61) 
(   62) 
(   63)       ! Initialization to zero
(   64)       ! ----------------------
(   65) 
(   66)       !     before fields       !       now fields        !      after fields       !
(   67)      ub   (:,:,:) = 0.e0   ;   un   (:,:,:) = 0.e0   ;   ua   (:,:,:) = 0.e0
(   68)      vb   (:,:,:) = 0.e0   ;   vn   (:,:,:) = 0.e0   ;   va   (:,:,:) = 0.e0
(   69)      wn   (:,:,:) = 0.e0   ;
(   70)      rotb (:,:,:) = 0.e0   ;   rotn (:,:,:) = 0.e0   ;
(   71)      hdivb(:,:,:) = 0.e0   ;   hdivn(:,:,:) = 0.e0   ;
(   72)      
(   73)      tb   (:,:,:) = 0.e0   ;   tn   (:,:,:) = 0.e0   ;   ta   (:,:,:) = 0.e0
(   74)      sb   (:,:,:) = 0.e0   ;   sn   (:,:,:) = 0.e0   ;   sa   (:,:,:) = 0.e0
(   75)      
(   76)      rhd  (:,:,:) = 0.e0
(   77)      rhop (:,:,:) = 0.e0
(   78)      rn2  (:,:,:) = 0.e0 
(   79)      
(   80) # 86
(   86)       ! free surface formulation
(   87)       sshb(:,:) = 0.e0      ! before sea-surface height
(   88)       sshn(:,:) = 0.e0      ! now    sea-surface height
(   89) 
(   90) 
(   91)       IF( ln_rstart ) THEN                    ! Restart from a file
(   92)          !                                    ! -------------------






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 9

(   93)          neuler = 1                              ! Set time-step indicator at nit000 (leap-frog)
(   94)          CALL rst_read                           ! Read the restart file
(   95)       ELSE
(   96)          !                                    ! Start from rest
(   97)          !                                    ! ---------------
(   98)          neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
(   99)          adatrj = 0._wp
(  100) !!DB: OLD
(  101) !         IF( cp_cfg == 'eel' ) THEN
(  102) !            CALL istate_eel                      ! EEL   configuration : start from pre-defined
(  103) !            !                                    !                       velocity and thermohaline fields
(  104) !         ELSEIF( cp_cfg == 'gyre' ) THEN         
(  105) !            CALL istate_gyre                     ! GYRE  configuration : start from pre-defined temperature
(  106) !            !                                    !                       and salinity fields 
(  107) !         ELSE
(  108)          !                                       ! Other configurations: Initial temperature and salinity fields
(  109) !#if defined 1
(  110)          CALL dta_tem( nit000 )                  ! read 3D temperature data
(  111)          tb(:,:,:) = t_dta(:,:,:)                ! use temperature data read
(  112)          tn(:,:,:) = t_dta(:,:,:)
(  113) !#else
(  114) !         IF(lwp) WRITE(numout,*)                 ! analytical temperature profile
(  115) !         IF(lwp) WRITE(numout,*)' Temperature initialization using an analytic profile'
(  116) !         CALL istate_tem
(  117) !#endif
(  118) !#if defined 1
(  119)          CALL dta_sal( nit000 )                  ! read 3D salinity data
(  120)          sb(:,:,:) = s_dta(:,:,:)                ! use salinity data read
(  121)          sn(:,:,:) = s_dta(:,:,:)
(  122) !#else
(  123) !         ! No salinity data
(  124) !         IF(lwp)WRITE(numout,*)                  ! analytical salinity profile
(  125) !         IF(lwp)WRITE(numout,*)' Salinity initialisation using a constant value'
(  126) !         CALL istate_sal
(  127) !#endif
(  128)          !         ENDIF
(  129) 
(  130)       ENDIF
(  131)       !                                       ! Vertical velocity
(  132)       !                                       ! -----------------
(  133)       CALL wzv( nit000 )                         ! from horizontal divergence
(  134) 
(  135)    END SUBROUTINE istate_init
(  136) 
(  137) 
(  138) 
(  139) 
(  140) 
(  141)    SUBROUTINE istate_uvg
(  142)       !!----------------------------------------------------------------------
(  143)       !!                  ***  ROUTINE istate_uvg  ***
(  144)       !!
(  145)       !! ** Purpose :   Compute the geostrophic velocities from (tn,sn) fields
(  146)       !!
(  147)       !! ** Method  :   Using the hydrostatic hypothesis the now hydrostatic 
(  148)       !!      pressure is computed by integrating the in-situ density from the
(  149)       !!      surface to the bottom.
(  150)       !!                 p=integral [ rau*g dz ]






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 10

(  151)       !!
(  152)       !! History :
(  153)       !!   8.1  !  01-09  (M. Levy, M. Ben Jelloul)  Original code
(  154)       !!   8.5  !  02-09  (G. Madec)  F90: Free form
(  155)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  156)       !!----------------------------------------------------------------------
(  157)       !! * Modules used
(  158)       USE eosbn2          ! eq. of state, Brunt Vaisala frequency (eos     routine)
(  159)       USE dynspg          ! surface pressure gradient             (dyn_spg routine)
(  160)       USE divcur          ! hor. divergence & rel. vorticity      (div_cur routine)
(  161)       USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(  162) 
(  163)       !! * Local declarations
(  164)       INTEGER ::   ji, jj, jk        ! dummy loop indices
(  165)       INTEGER ::   indic             ! ???
(  166)       REAL(wp) ::   &
(  167)          zmsv, zphv, zmsu, zphu,  &  ! temporary scalars
(  168)          zalfg
(  169)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(  170)          zprn                        ! workspace
(  171)       !!----------------------------------------------------------------------
(  172) 
(  173)       IF(lwp) WRITE(numout,*) 
(  174)       IF(lwp) WRITE(numout,*) 'istate_uvg : Start from Geostrophy'
(  175)       IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
(  176) 
(  177)       ! Compute the now hydrostatic pressure
(  178)       ! ------------------------------------
(  179) 
(  180)       zalfg = 0.5 * grav * rau0
(  181)       ! Surface value
(  182)       zprn(:,:,1) = zalfg * e3w(1) * ( 1 + rhd(:,:,1) )
(  183) 
(  184)       ! Vertical integration from the surface
(  185)       DO jk = 2, jpkm1
(  186)          zprn(:,:,jk) = zprn(:,:,jk-1)   &
(  187)             &         + zalfg * e3w(jk) * ( 2. + rhd(:,:,jk) + rhd(:,:,jk-1) )
(  188)       END DO  
(  189) 
(  190)       ! Compute geostrophic balance
(  191)       ! ---------------------------
(  192) 
(  193)       DO jk = 1, jpkm1
(  194)          DO jj = 2, jpjm1
(  195)             DO ji = 2, jpim1   ! vertor opt.
(  196)                zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
(  197)                                + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1.  )
(  198)                zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
(  199)                     + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
(  200)                     + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
(  201)                     + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
(  202)                zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)
(  203) 
(  204)                zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
(  205)                                + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1.  )
(  206)                zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
(  207)                     + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
(  208)                     + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 11

(  209)                     + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
(  210)                zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)
(  211) 
(  212)                ! Compute the geostrophic velocities
(  213)                un(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
(  214)                vn(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
(  215)             END DO
(  216)          END DO
(  217)       END DO
(  218) 
(  219)       IF(lwp) WRITE(numout,*) '         we force to zero bottom velocity'
(  220) 
(  221)       ! Susbtract the bottom velocity (level jpk-1 for flat bottom case)
(  222)       ! to have a zero bottom velocity
(  223) 
(  224)       DO jk = 1, jpkm1
(  225)          un(:,:,jk) = ( un(:,:,jk) - un(:,:,jpkm1) ) * umask(:,:,jk)
(  226)          vn(:,:,jk) = ( vn(:,:,jk) - vn(:,:,jpkm1) ) * vmask(:,:,jk)
(  227)       END DO
(  228) 
(  229)       CALL lbc_lnk( un, 'U', -1. )
(  230)       CALL lbc_lnk( vn, 'V', -1. )
(  231)       
(  232)       ub(:,:,:) = un(:,:,:)
(  233)       vb(:,:,:) = vn(:,:,:)
(  234)       
(  235)       ! WARNING !!!!!
(  236)       ! after initializing u and v, we need to calculate the initial streamfunction bsf.
(  237)       ! Otherwise, only the trend will be computed and the model will blow up (inconsistency).
(  238)       
(  239)       ! to do that, we call dyn_spg with a special trick:
(  240)       ! we fill ua and va with the velocities divided by dt,
(  241)       ! and the streamfunction will be brought to the right
(  242)       ! value assuming the velocities have been set up in
(  243)       ! one time step.
(  244)       ! we then set bsfd to zero (first guess for next step
(  245)       ! is d(psi)/dt = 0.)
(  246) 
(  247)       !  sets up s false trend to calculate the barotropic
(  248)       !  streamfunction.
(  249) 
(  250)       ua(:,:,:) = ub(:,:,:) / rdt
(  251)       va(:,:,:) = vb(:,:,:) / rdt
(  252) 
(  253)       ! calls dyn_spg. we assume euler time step, starting from rest.
(  254)       indic = 0
(  255)       CALL dyn_spg( nit000, indic )       ! surface pressure gradient
(  256) 
(  257)       ! the new velocity is ua*rdt
(  258) 
(  259)       CALL lbc_lnk( ua, 'U', -1. )
(  260)       CALL lbc_lnk( va, 'V', -1. )
(  261) 
(  262)       ub(:,:,:) = ua(:,:,:) * rdt
(  263)       vb(:,:,:) = va(:,:,:) * rdt
(  264)       ua(:,:,:) = 0.e0
(  265)       va(:,:,:) = 0.e0
(  266)       un(:,:,:) = ub(:,:,:)






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 12

(  267)       vn(:,:,:) = vb(:,:,:)
(  268)        
(  269) 
(  270) # 273
(  273)       ! Compute the divergence and curl
(  274) 
(  275)       CALL div_cur( nit000 )            ! now horizontal divergence and curl
(  276) 
(  277)       hdivb(:,:,:) = hdivn(:,:,:)       ! set the before to the now value
(  278)       rotb (:,:,:) = rotn (:,:,:)       ! set the before to the now value
(  279) 
(  280)    END SUBROUTINE istate_uvg
(  281) 
(  282)    !!=====================================================================
(  283) END MODULE istate
