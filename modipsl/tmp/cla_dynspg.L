


PGF90 (Version     10.5)          09/23/2010  09:49:26      page 1

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

Filename: cla_dynspg.F90

(    1) MODULE cla_dynspg
(    2)    !!======================================================================
(    3)    !!                       ***  cla_dynspg  ***
(    4)    !!======================================================================
(    5)    !!   dyn_spg      : update the momentum trend with the surface pressure
(    6)    !!                  gradient in the free surface constant volume case
(    7)    !!                  with vector optimization
(    8)    !!----------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers 
(   11)    USE dom_oce         ! ocean space and time domain 
(   12)    USE zdf_oce         ! ocean vertical physics
(   13)    USE obc_oce         ! Lateral open boundary condition
(   14)    USE sol_oce         ! solver variables
(   15)    USE phycst          ! physical constants
(   16)    USE ocesbc          ! ocean surface boundary condition (fluxes)
(   17)    USE flxrnf          ! ocean runoffs
(   18)    USE solpcg          ! preconditionned conjugate gradient solver
(   19)    USE solsor          ! Successive Over-relaxation solver
(   20)    USE solfet          ! FETI solver
(   21)    USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
(   22)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE lib_mpp         ! distribued memory computing
(   25)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC dyn_spg_cla   ! routine called by step.F90
(   32) 
(   33)    !! * Substitutions
(   34) #  include "domzgr_substitute.h90"
(   35) #  include "vectopt_loop_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cla_dynspg.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE dyn_spg_cla( kt ) 
(   45)       !!----------------------------------------------------------------------
(   46)       !!              ***  routine dyn_spg_cross_land  ***






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 2

(   47)       !!
(   48)       !! ** Purpose :
(   49)       !!
(   50)       !! ** Method :
(   51)       !!
(   52)       !! ** Action :
(   53)       !!
(   54)       !! History :
(   55)       !!        !         (A. Bozec)  Original code
(   56)       !!   8.5  !  02-11  (A. Bozec)  F90: Free form and module
(   57)       !!---------------------------------------------------------------------
(   58)       !! * Arguments
(   59)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step
(   60)       !! * Local declarations
(   61)       INTEGER  ::   ji, jj, jk                ! dummy loop indices
(   62)       INTEGER  ::   ii0, ii1, ij0, ij1        ! temporary integer
(   63)       REAL(wp) ::    &    
(   64)          zempmed, zempred,   &                ! EMP on Med Sea ans Red Sea
(   65)          zwei,   &                            !              
(   66)          zisw_rs, zurw_rs, zbrw_rs,      &    ! imposed transport Red sea
(   67)          zisw_ms, zurw_ms, zbrw_ms, zmrw_ms   ! imposed transport Med Sea
(   68)       !!----------------------------------------------------------------------
(   69) 
(   70)       ! Different velocities for straits ( Gibraltar, Bab el Mandeb...)
(   71)          
(   72)       ! Control print
(   73)       ! -------------
(   74)       IF( kt == nit000 ) THEN 
(   75)          IF(lwp) WRITE(numout,*)
(   76)          IF(lwp) WRITE(numout,*) 'dynspg_cross_land : cross land advection on surface '
(   77)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~   pressure '
(   78)          IF(lwp) WRITE(numout,*) ' '
(   79)       ENDIF
(   80) 
(   81)       ! EMP on Mediterranean Sea and Red Sea 
(   82)       ! ------------------------------------
(   83)       ! compute the emp in Mediterranean Sea
(   84)       zempmed = 0.e0
(   85)       zwei = 0.e0
(   86)       ij0 =  96   ;   ij1 = 110
(   87)       ii0 = 141   ;   ii1 = 181
(   88)       DO jj = mj0(ij0), mj1(ij1)
(   89)          DO ji = mi0(ii0),mi1(ii1)
(   90)             zwei    = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
(   91)             zempmed = zempmed + emp(ji,jj) * zwei
(   92)          END DO
(   93)       END DO
(   94)       IF( lk_mpp )   CALL mpp_sum( zempmed )      ! sum with other processors value
(   95) 
(   96)       ! minus 2 points in Red Sea and 3 in Atlantic 
(   97)       ij0 =  96   ;   ij1 =  96
(   98)       ii0 = 148   ;   ii1 = 148
(   99)       DO jj = mj0(ij0), mj1(ij1)
(  100)          DO ji = mi0(ii0),mi1(ii1)
(  101)             zempmed = zempmed - emp(ji  ,jj) * tmask(ji  ,jj,1) * e1t(ji  ,jj) * e2t(ji  ,jj)   &
(  102)                &              - emp(ji+1,jj) * tmask(ji+1,jj,1) * e1t(ji+1,jj) * e2t(ji+1,jj)   
(  103)          END DO
(  104)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 3

(  105)       ! we convert in m3
(  106)       zempmed = zempmed * 1.e-3
(  107) 
(  108)       ! compute the emp in Red Sea   
(  109)       zempred = 0.e0
(  110)       zwei = 0.e0
(  111)       ij0 =  87   ;   ij1 =  96
(  112)       ii0 = 148   ;   ii1 = 160
(  113)       DO jj = mj0(ij0), mj1(ij1)
(  114)          DO ji = mi0(ii0),mi1(ii1)
(  115)             zwei      = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
(  116)             zempred   = zempred + emp(ji,jj) * zwei
(  117)          END DO
(  118)       END DO
(  119)       IF( lk_mpp )   CALL mpp_sum( zempred )      ! sum with other processors value
(  120) 
(  121)       ! we convert in m3
(  122)       zempred = zempred * 1.e-3
(  123) 
(  124)       ! New Transport at Bab el Mandeb and Gibraltar
(  125)       ! --------------------------------------------
(  126) 
(  127)       ! imposed transport at Bab el Mandeb
(  128)       zisw_rs = 0.4e6        ! inflow surface water
(  129)       zurw_rs = 0.2e6        ! upper recirculation water
(  130) !!Alex      zbrw_rs = 1.2e6        ! bottom  recirculation water
(  131)       zbrw_rs = 0.5e6        ! bottom  recirculation water
(  132) 
(  133)       ! imposed transport at Gibraltar
(  134)       zisw_ms  = 0.8e6          ! atlantic-mediterranean  water
(  135)       zmrw_ms  = 0.7e6          ! middle recirculation water
(  136)       zurw_ms  = 2.5e6          ! upper  recirculation water 
(  137)       zbrw_ms  = 3.5e6          ! bottom recirculation water 
(  138) 
(  139)       ! Different velocities for straits ( Gibraltar, Bab el Mandeb )
(  140)       ! -------------------------------------------------------------
(  141) 
(  142)       ! Bab el Mandeb
(  143)       ! -------------
(  144)       ! 160,88 north point Bab el Mandeb
(  145)       ij0 =  88   ;   ij1 =  88
(  146)       ii0 = 160   ;   ii1 = 160
(  147)       DO jj = mj0(ij0), mj1(ij1)
(  148)          DO ji = mi0(ii0),mi1(ii1)
(  149)             ua(ji,jj  ,: ) = 0.e0  !  North East Bab el Mandeb 
(  150)          END DO
(  151)       END DO
(  152)       !                              ! surface
(  153)       DO jk = 1,  8                                      
(  154)          DO jj = mj0(ij0), mj1(ij1)
(  155)             DO ji = mi0(ii0),mi1(ii1)
(  156)                ua(ji, jj,jk) = -( ( zisw_rs + zempred ) / 8. ) / ( e2u(ji, jj) * fse3t(ji, jj,jk) )     
(  157)             END DO
(  158)          END DO
(  159)       END DO
(  160)       !                              ! deeper
(  161)       DO jj = mj0(ij0), mj1(ij1)
(  162)          DO ji = mi0(ii0),mi1(ii1)






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 4

(  163)             ua(ji, jj,21) = - zbrw_rs / ( e2u(ji, jj) * fse3t(ji, jj,21) )
(  164)          END DO
(  165)       END DO
(  166) 
(  167)       ! 160,87 south point Bab el Mandeb
(  168)       ij0 =  87   ;   ij1 =  87
(  169)       ii0 = 160   ;   ii1 = 160
(  170)       DO jj = mj0(ij0), mj1(ij1)
(  171)          DO ji = mi0(ii0),mi1(ii1)
(  172)             ua(ji,jj  ,: ) = 0.e0  !  South East Bab el Mandeb 
(  173)          END DO
(  174)       END DO
(  175)       DO jj = mj0(ij0), mj1(ij1)
(  176)          DO ji = mi0(ii0),mi1(ii1)
(  177)             ua(ji, jj,21) =  ( zisw_rs + zbrw_rs ) / ( e2u(ji,jj )*fse3t(ji, jj,21) )      
(  178)          END DO
(  179)       END DO
(  180) 
(  181)       ! Gibraltar
(  182)       ! ---------
(  183) 
(  184)       ! initialisation of velocity at concerned points 
(  185)       ! 139, 101 south point in Gibraltar 
(  186)       ij0 = 101   ;   ij1 = 101
(  187)       ii0 = 139   ;   ii1 = 139
(  188)       DO jj = mj0(ij0), mj1(ij1)
(  189)          DO ji = mi0(ii0),mi1(ii1)
(  190)             ua(ji,jj  ,: ) = 0.e0  !  South West Gibraltar
(  191)             ua(ji,jj+1,: ) = 0.e0  !  North West Gibraltar
(  192)          END DO
(  193)       END DO
(  194)       !                            ! surface
(  195)       DO jk = 1, 14                      
(  196)          DO jj = mj0(ij0), mj1(ij1)
(  197)             DO ji = mi0(ii0),mi1(ii1)
(  198)                ua(ji,jj,jk) =  ( ( zisw_ms + zempmed ) / 14. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
(  199)             END DO
(  200)          END DO
(  201)       END DO
(  202)       !                            ! middle circulation
(  203)       DO jk = 15, 20                      
(  204)          DO jj = mj0(ij0), mj1(ij1)
(  205)             DO ji = mi0(ii0),mi1(ii1)
(  206)                ua(ji,jj,jk) =  ( zmrw_ms / 6. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
(  207)             END DO
(  208)          END DO
(  209)       END DO
(  210)       !                            ! deeper 
(  211)       DO jj = mj0(ij0), mj1(ij1)
(  212)          DO ji = mi0(ii0),mi1(ii1)
(  213)             ua(ji,jj,21) =             zurw_ms   / ( e2u(ji,jj) * fse3t(ji,jj,21) )
(  214)             ua(ji,jj,22) = ( zbrw_ms - zurw_ms ) / ( e2u(ji,jj) * fse3t(ji,jj,22) )
(  215)          END DO
(  216)       END DO
(  217) 
(  218)       ! 139,102 north point in Gibraltar
(  219)       ij0 = 102   ;   ij1 = 102
(  220)       ii0 = 139   ;   ii1 = 139






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 5

(  221)       DO jj = mj0(ij0), mj1(ij1)
(  222)          DO ji = mi0(ii0),mi1(ii1)
(  223)             ua(ji,jj  ,: ) = 0.e0  !  North West Gibraltar
(  224)          END DO
(  225)       END DO
(  226)       DO jk = 15, 20                      
(  227)          DO jj = mj0(ij0), mj1(ij1)
(  228)             DO ji = mi0(ii0),mi1(ii1)
(  229)                ua(ji,jj,jk) = -( zmrw_ms / 6. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
(  230)             END DO
(  231)          END DO
(  232)       END DO
(  233)       !                            ! deeper
(  234)       DO jj = mj0(ij0), mj1(ij1)
(  235)          DO ji = mi0(ii0),mi1(ii1)
(  236)             ua(ji,jj,22) = -( zisw_ms + zbrw_ms ) / ( e2u(ji,jj) * fse3t(ji,jj,22) )
(  237)          END DO
(  238)       END DO
(  239) 
(  240)    END SUBROUTINE dyn_spg_cla
(  241) 
(  242)    !!======================================================================
(  243) END MODULE cla_dynspg









































PGF90 (Version     10.5)          09/23/2010  09:49:26      page 6

(    1) # 1 "cla_dynspg.F90"
(    1) MODULE cla_dynspg
(    2)    !!======================================================================
(    3)    !!                       ***  cla_dynspg  ***
(    4)    !!======================================================================
(    5)    !!   dyn_spg      : update the momentum trend with the surface pressure
(    6)    !!                  gradient in the free surface constant volume case
(    7)    !!                  with vector optimization
(    8)    !!----------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers 
(   11)    USE dom_oce         ! ocean space and time domain 
(   12)    USE zdf_oce         ! ocean vertical physics
(   13)    USE obc_oce         ! Lateral open boundary condition
(   14)    USE sol_oce         ! solver variables
(   15)    USE phycst          ! physical constants
(   16)    USE ocesbc          ! ocean surface boundary condition (fluxes)
(   17)    USE flxrnf          ! ocean runoffs
(   18)    USE solpcg          ! preconditionned conjugate gradient solver
(   19)    USE solsor          ! Successive Over-relaxation solver
(   20)    USE solfet          ! FETI solver
(   21)    USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
(   22)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   23)    USE in_out_manager  ! I/O manager
(   24)    USE lib_mpp         ! distribued memory computing
(   25)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC dyn_spg_cla   ! routine called by step.F90
(   32) 
(   33)    !! * Substitutions
(   34) # 1 "./domzgr_substitute.h90"
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
(   48) # 35 "cla_dynspg.F90"
(   35) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 7

(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 36 "cla_dynspg.F90"
(   36) # 36
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cla_dynspg.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE dyn_spg_cla( kt ) 
(   45)       !!----------------------------------------------------------------------
(   46)       !!              ***  routine dyn_spg_cross_land  ***
(   47)       !!
(   48)       !! ** Purpose :
(   49)       !!
(   50)       !! ** Method :
(   51)       !!
(   52)       !! ** Action :
(   53)       !!
(   54)       !! History :
(   55)       !!        !         (A. Bozec)  Original code
(   56)       !!   8.5  !  02-11  (A. Bozec)  F90: Free form and module
(   57)       !!---------------------------------------------------------------------
(   58)       !! * Arguments
(   59)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step
(   60)       !! * Local declarations
(   61)       INTEGER  ::   ji, jj, jk                ! dummy loop indices
(   62)       INTEGER  ::   ii0, ii1, ij0, ij1        ! temporary integer
(   63)       REAL(wp) ::    &    
(   64)          zempmed, zempred,   &                ! EMP on Med Sea ans Red Sea
(   65)          zwei,   &                            !              
(   66)          zisw_rs, zurw_rs, zbrw_rs,      &    ! imposed transport Red sea
(   67)          zisw_ms, zurw_ms, zbrw_ms, zmrw_ms   ! imposed transport Med Sea
(   68)       !!----------------------------------------------------------------------
(   69) 
(   70)       ! Different velocities for straits ( Gibraltar, Bab el Mandeb...)
(   71)          
(   72)       ! Control print
(   73)       ! -------------
(   74)       IF( kt == nit000 ) THEN 
(   75)          IF(lwp) WRITE(numout,*)
(   76)          IF(lwp) WRITE(numout,*) 'dynspg_cross_land : cross land advection on surface '
(   77)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~   pressure '
(   78)          IF(lwp) WRITE(numout,*) ' '
(   79)       ENDIF
(   80) 
(   81)       ! EMP on Mediterranean Sea and Red Sea 
(   82)       ! ------------------------------------
(   83)       ! compute the emp in Mediterranean Sea
(   84)       zempmed = 0.e0
(   85)       zwei = 0.e0
(   86)       ij0 =  96   ;   ij1 = 110
(   87)       ii0 = 141   ;   ii1 = 181






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 8

(   88)       DO jj = mj0(ij0), mj1(ij1)
(   89)          DO ji = mi0(ii0),mi1(ii1)
(   90)             zwei    = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
(   91)             zempmed = zempmed + emp(ji,jj) * zwei
(   92)          END DO
(   93)       END DO
(   94)       IF( lk_mpp )   CALL mpp_sum( zempmed )      ! sum with other processors value
(   95) 
(   96)       ! minus 2 points in Red Sea and 3 in Atlantic 
(   97)       ij0 =  96   ;   ij1 =  96
(   98)       ii0 = 148   ;   ii1 = 148
(   99)       DO jj = mj0(ij0), mj1(ij1)
(  100)          DO ji = mi0(ii0),mi1(ii1)
(  101)             zempmed = zempmed - emp(ji  ,jj) * tmask(ji  ,jj,1) * e1t(ji  ,jj) * e2t(ji  ,jj)   &
(  102)                &              - emp(ji+1,jj) * tmask(ji+1,jj,1) * e1t(ji+1,jj) * e2t(ji+1,jj)   
(  103)          END DO
(  104)       END DO
(  105)       ! we convert in m3
(  106)       zempmed = zempmed * 1.e-3
(  107) 
(  108)       ! compute the emp in Red Sea   
(  109)       zempred = 0.e0
(  110)       zwei = 0.e0
(  111)       ij0 =  87   ;   ij1 =  96
(  112)       ii0 = 148   ;   ii1 = 160
(  113)       DO jj = mj0(ij0), mj1(ij1)
(  114)          DO ji = mi0(ii0),mi1(ii1)
(  115)             zwei      = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
(  116)             zempred   = zempred + emp(ji,jj) * zwei
(  117)          END DO
(  118)       END DO
(  119)       IF( lk_mpp )   CALL mpp_sum( zempred )      ! sum with other processors value
(  120) 
(  121)       ! we convert in m3
(  122)       zempred = zempred * 1.e-3
(  123) 
(  124)       ! New Transport at Bab el Mandeb and Gibraltar
(  125)       ! --------------------------------------------
(  126) 
(  127)       ! imposed transport at Bab el Mandeb
(  128)       zisw_rs = 0.4e6        ! inflow surface water
(  129)       zurw_rs = 0.2e6        ! upper recirculation water
(  130) !!Alex      zbrw_rs = 1.2e6        ! bottom  recirculation water
(  131)       zbrw_rs = 0.5e6        ! bottom  recirculation water
(  132) 
(  133)       ! imposed transport at Gibraltar
(  134)       zisw_ms  = 0.8e6          ! atlantic-mediterranean  water
(  135)       zmrw_ms  = 0.7e6          ! middle recirculation water
(  136)       zurw_ms  = 2.5e6          ! upper  recirculation water 
(  137)       zbrw_ms  = 3.5e6          ! bottom recirculation water 
(  138) 
(  139)       ! Different velocities for straits ( Gibraltar, Bab el Mandeb )
(  140)       ! -------------------------------------------------------------
(  141) 
(  142)       ! Bab el Mandeb
(  143)       ! -------------
(  144)       ! 160,88 north point Bab el Mandeb
(  145)       ij0 =  88   ;   ij1 =  88






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 9

(  146)       ii0 = 160   ;   ii1 = 160
(  147)       DO jj = mj0(ij0), mj1(ij1)
(  148)          DO ji = mi0(ii0),mi1(ii1)
(  149)             ua(ji,jj  ,: ) = 0.e0  !  North East Bab el Mandeb 
(  150)          END DO
(  151)       END DO
(  152)       !                              ! surface
(  153)       DO jk = 1,  8                                      
(  154)          DO jj = mj0(ij0), mj1(ij1)
(  155)             DO ji = mi0(ii0),mi1(ii1)
(  156)                ua(ji, jj,jk) = -( ( zisw_rs + zempred ) / 8. ) / ( e2u(ji, jj) * e3t(jk) )     
(  157)             END DO
(  158)          END DO
(  159)       END DO
(  160)       !                              ! deeper
(  161)       DO jj = mj0(ij0), mj1(ij1)
(  162)          DO ji = mi0(ii0),mi1(ii1)
(  163)             ua(ji, jj,21) = - zbrw_rs / ( e2u(ji, jj) * e3t(21) )
(  164)          END DO
(  165)       END DO
(  166) 
(  167)       ! 160,87 south point Bab el Mandeb
(  168)       ij0 =  87   ;   ij1 =  87
(  169)       ii0 = 160   ;   ii1 = 160
(  170)       DO jj = mj0(ij0), mj1(ij1)
(  171)          DO ji = mi0(ii0),mi1(ii1)
(  172)             ua(ji,jj  ,: ) = 0.e0  !  South East Bab el Mandeb 
(  173)          END DO
(  174)       END DO
(  175)       DO jj = mj0(ij0), mj1(ij1)
(  176)          DO ji = mi0(ii0),mi1(ii1)
(  177)             ua(ji, jj,21) =  ( zisw_rs + zbrw_rs ) / ( e2u(ji,jj )*e3t(21) )      
(  178)          END DO
(  179)       END DO
(  180) 
(  181)       ! Gibraltar
(  182)       ! ---------
(  183) 
(  184)       ! initialisation of velocity at concerned points 
(  185)       ! 139, 101 south point in Gibraltar 
(  186)       ij0 = 101   ;   ij1 = 101
(  187)       ii0 = 139   ;   ii1 = 139
(  188)       DO jj = mj0(ij0), mj1(ij1)
(  189)          DO ji = mi0(ii0),mi1(ii1)
(  190)             ua(ji,jj  ,: ) = 0.e0  !  South West Gibraltar
(  191)             ua(ji,jj+1,: ) = 0.e0  !  North West Gibraltar
(  192)          END DO
(  193)       END DO
(  194)       !                            ! surface
(  195)       DO jk = 1, 14                      
(  196)          DO jj = mj0(ij0), mj1(ij1)
(  197)             DO ji = mi0(ii0),mi1(ii1)
(  198)                ua(ji,jj,jk) =  ( ( zisw_ms + zempmed ) / 14. ) / ( e2u(ji,jj) * e3t(jk) ) 
(  199)             END DO
(  200)          END DO
(  201)       END DO
(  202)       !                            ! middle circulation
(  203)       DO jk = 15, 20                      






PGF90 (Version     10.5)          09/23/2010  09:49:26      page 10

(  204)          DO jj = mj0(ij0), mj1(ij1)
(  205)             DO ji = mi0(ii0),mi1(ii1)
(  206)                ua(ji,jj,jk) =  ( zmrw_ms / 6. ) / ( e2u(ji,jj) * e3t(jk) ) 
(  207)             END DO
(  208)          END DO
(  209)       END DO
(  210)       !                            ! deeper 
(  211)       DO jj = mj0(ij0), mj1(ij1)
(  212)          DO ji = mi0(ii0),mi1(ii1)
(  213)             ua(ji,jj,21) =             zurw_ms   / ( e2u(ji,jj) * e3t(21) )
(  214)             ua(ji,jj,22) = ( zbrw_ms - zurw_ms ) / ( e2u(ji,jj) * e3t(22) )
(  215)          END DO
(  216)       END DO
(  217) 
(  218)       ! 139,102 north point in Gibraltar
(  219)       ij0 = 102   ;   ij1 = 102
(  220)       ii0 = 139   ;   ii1 = 139
(  221)       DO jj = mj0(ij0), mj1(ij1)
(  222)          DO ji = mi0(ii0),mi1(ii1)
(  223)             ua(ji,jj  ,: ) = 0.e0  !  North West Gibraltar
(  224)          END DO
(  225)       END DO
(  226)       DO jk = 15, 20                      
(  227)          DO jj = mj0(ij0), mj1(ij1)
(  228)             DO ji = mi0(ii0),mi1(ii1)
(  229)                ua(ji,jj,jk) = -( zmrw_ms / 6. ) / ( e2u(ji,jj) * e3t(jk) ) 
(  230)             END DO
(  231)          END DO
(  232)       END DO
(  233)       !                            ! deeper
(  234)       DO jj = mj0(ij0), mj1(ij1)
(  235)          DO ji = mi0(ii0),mi1(ii1)
(  236)             ua(ji,jj,22) = -( zisw_ms + zbrw_ms ) / ( e2u(ji,jj) * e3t(22) )
(  237)          END DO
(  238)       END DO
(  239) 
(  240)    END SUBROUTINE dyn_spg_cla
(  241) 
(  242)    !!======================================================================
(  243) END MODULE cla_dynspg
