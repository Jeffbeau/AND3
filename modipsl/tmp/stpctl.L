


PGF90 (Version     10.5)          09/23/2010  09:50:10      page 1

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

Filename: stpctl.F90

(    1) MODULE stpctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  stpctl  ***
(    4)    !! Ocean run control :  gross check of the ocean time stepping
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   stp_ctl      : Control the run
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE sol_oce         ! ocean space and time domain variables 
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE solisl          ! ???
(   16)    USE diawri          ! ocean output file 
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE lib_mpp         ! distributed memory computing
(   19)    USE dynspg_oce      ! pressure gradient schemes 
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC stp_ctl           ! routine called by step.F90
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE stp_ctl( kt, kindic )
(   31)       !!----------------------------------------------------------------------
(   32)       !!                    ***  ROUTINE stp_ctl  ***
(   33)       !!                     
(   34)       !! ** Purpose :   Control the run
(   35)       !!
(   36)       !! ** Method  : - Save the time step in numstp
(   37)       !!              - Print it each 50 time steps
(   38)       !!              - Print solver statistics in numsol 
(   39)       !!              - Stop the run IF problem for the solver ( indec < 0 )
(   40)       !!
(   41)       !! History :
(   42)       !!        !  91-03  ()
(   43)       !!        !  91-11  (G. Madec)
(   44)       !!        !  92-06  (M. Imbard)
(   45)       !!        !  97-06  (A.M. Treguier)
(   46)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !! * Arguments
(   49)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   50)       INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence
(   51) 
(   52)       !! * local declarations
(   53)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   54)       INTEGER  ::   ii, ij, ik              ! temporary integers
(   55)       REAL(wp) ::   zumax, zsmin            ! temporary scalars
(   56)       INTEGER, DIMENSION(3) ::   ilocu      ! 
(   57)       INTEGER, DIMENSION(2) ::   ilocs      ! 
(   58)       CHARACTER(len=80) :: clname
(   59)       !!----------------------------------------------------------------------
(   60)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   61)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/stpctl.F90,v 1.9 2006/03/21 07:57:56 opalod Exp $ 
(   62)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   63)       !!----------------------------------------------------------------------
(   64) 
(   65)       IF( kt == nit000 .AND. lwp ) THEN
(   66)          WRITE(numout,*)
(   67)          WRITE(numout,*) 'stp_ctl : time-stepping control'
(   68)          WRITE(numout,*) '~~~~~~~'
(   69)          ! open time.step file
(   70)          clname = 'time.step'
(   71)          CALL ctlopn( numstp, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
(   72)       ENDIF
(   73) 
(   74)       ! save the current time step in numstp
(   75)       ! ------------------------------------
(   76) !      IF(lwp) WRITE(numstp,9100) kt
(   77) !      IF(lwp) REWIND(numstp)
(   78) !!DB -- on some machines, the above does not write every timestep
(   79) !!      On drakes, the below made no difference.
(   80)       if(lwp) then
(   81)          WRITE(numstp,9100) kt
(   82)          call flush(numstp)
(   83)          REWIND(numstp)
(   84)       endif
(   85) 
(   86) 9100  FORMAT(1x, i8)
(   87) 
(   88) 
(   89)       ! elliptic solver statistics (if required)
(   90)       ! --------------------------
(   91)       IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
(   92)       ! Solver
(   93)       IF(lwp) WRITE(numsol,9200) kt, niter, res, SQRT(epsr)/eps
(   94) 
(   95)       ! Islands (if exist)
(   96)       IF( lk_isl )   CALL isl_stp_ctl( kt, kindic )
(   97) 
(   98) 
(   99)       ! Output in numwso and numwvo IF kindic<0
(  100)       ! ---------------------------------------
(  101)       !    (i.e. problem for the solver)
(  102)       IF( kindic < 0 ) THEN
(  103)          IF(lwp) THEN
(  104)             WRITE(numout,*) ' stpctl: the elliptic solver DO not converge or explode'






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 3

(  105)             WRITE(numout,*) ' ====== '
(  106)             WRITE(numout,9200) kt, niter, res, sqrt(epsr)/eps
(  107)             WRITE(numout,*)
(  108)             WRITE(numout,*) ' stpctl: output of last fields in numwso'
(  109)             WRITE(numout,*) '                                  numwvo'
(  110)             WRITE(numout,*) ' ======  *******************************'
(  111)          ENDIF
(  112)          CALL dia_wri( kt, kindic )
(  113)       ENDIF
(  114)       ENDIF
(  115) 
(  116) 9200  FORMAT(' it :', i8, ' niter :', i4, ' res :',e20.10,' b :',e20.10)
(  117) 
(  118)       ! Test maximum of velocity (zonal only)
(  119)       ! ------------------------
(  120)       !! zumax = MAXVAL( ABS( un(:,:,:) ) )   ! slower than the following loop on NEC SX5
(  121)       zumax = 0.e0
(  122)       DO jk = 1, jpk
(  123)          DO jj = 1, jpj
(  124)             DO ji = 1, jpi
(  125)                zumax = MAX(zumax,ABS(un(ji,jj,jk)))
(  126)           END DO 
(  127)         END DO 
(  128)       END DO        
(  129)       IF( lk_mpp )   CALL mpp_max( zumax )   ! max over the global domain
(  130) 
(  131)       IF( MOD( kt, nwrite ) == 1 ) THEN
(  132)          IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(U) max: ', zumax
(  133)       ENDIF
(  134)       IF( zumax > 20.) THEN
(  135)          IF( lk_mpp ) THEN
(  136)             CALL mpp_maxloc(ABS(un),umask,zumax,ii,ij,ik)
(  137)          ELSE
(  138)             ilocu = MAXLOC( ABS( un(:,:,:) ) )
(  139)             ii = ilocu(1) + nimpp - 1
(  140)             ij = ilocu(2) + njmpp - 1
(  141)             ik = ilocu(3)
(  142)          ENDIF
(  143)          IF(lwp) THEN
(  144)             WRITE(numout,cform_err)
(  145)             WRITE(numout,*) ' stpctl: the zonal velocity is larger than 20 m/s'
(  146)             WRITE(numout,*) ' ====== '
(  147)             WRITE(numout,9400) kt, zumax, ii, ij, ik
(  148)             WRITE(numout,*)
(  149)             WRITE(numout,*) '          output of last fields in numwso'
(  150)          ENDIF
(  151)          kindic  = -3
(  152) 
(  153)          CALL dia_wri( kt, kindic )
(  154)       ENDIF
(  155) 9400  FORMAT (' kt=',i6,' max abs(U): ',1pg11.4,', i j k: ',3i4)
(  156) 
(  157) 
(  158)       ! Test minimum of salinity
(  159)       ! ------------------------
(  160)       !! zsmin = MINVAL( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )    
(  161)       !                slower than the following loop on NEC SX5
(  162)       zsmin = 100.e0






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 4

(  163)       DO jj = 2, jpjm1
(  164)          DO ji = 1, jpi
(  165)             IF( tmask(ji,jj,1) == 1) zsmin = MIN(zsmin,sn(ji,jj,1))
(  166)          END DO
(  167)       END DO
(  168)       IF( lk_mpp )   CALL mpp_min( zsmin )   ! min over the global domain
(  169) 
(  170)       IF( MOD( kt, nwrite ) == 1 ) THEN
(  171)          IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' SSS min:', zsmin
(  172)       ENDIF
(  173)       IF( zsmin < 0.) THEN 
(  174)          IF (lk_mpp) THEN
(  175)             CALL mpp_minloc ( sn(:,:,1),tmask(:,:,1), zsmin, ii,ij )
(  176)          ELSE
(  177)             ilocs = MINLOC( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )
(  178)             ii = ilocs(1) + nimpp - 1
(  179)             ij = ilocs(2) + njmpp - 1
(  180)          END IF
(  181) 
(  182)          IF(lwp) THEN
(  183)             WRITE(numout,cform_err)
(  184)             WRITE(numout,*) 'stp_ctl : NEGATIVE sea surface salinity'
(  185)             WRITE(numout,*) '======= '
(  186)             WRITE(numout,9500) kt, zsmin, ii, ij
(  187)             WRITE(numout,*)
(  188)             WRITE(numout,*) '          output of last fields in numwso'
(  189)          ENDIF
(  190)          IF( kindic < 0 ) THEN
(  191)             IF(lwp) WRITE(numout,*) ' stpctl diabort done. We wont do it again '
(  192)          ELSE 
(  193)             kindic  = -3
(  194)             CALL dia_wri(kt,kindic)
(  195)          ENDIF
(  196)       ENDIF
(  197) 9500  FORMAT (' kt=',i6,' min SSS: ',1pg11.4,', i j: ',2i4)
(  198) 
(  199)    END SUBROUTINE stp_ctl
(  200) 
(  201)    !!======================================================================
(  202) END MODULE stpctl
























PGF90 (Version     10.5)          09/23/2010  09:50:10      page 5

(    1) # 1 "stpctl.F90"
(    1) MODULE stpctl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  stpctl  ***
(    4)    !! Ocean run control :  gross check of the ocean time stepping
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   stp_ctl      : Control the run
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE sol_oce         ! ocean space and time domain variables 
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE solisl          ! ???
(   16)    USE diawri          ! ocean output file 
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE lib_mpp         ! distributed memory computing
(   19)    USE dynspg_oce      ! pressure gradient schemes 
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC stp_ctl           ! routine called by step.F90
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE stp_ctl( kt, kindic )
(   31)       !!----------------------------------------------------------------------
(   32)       !!                    ***  ROUTINE stp_ctl  ***
(   33)       !!                     
(   34)       !! ** Purpose :   Control the run
(   35)       !!
(   36)       !! ** Method  : - Save the time step in numstp
(   37)       !!              - Print it each 50 time steps
(   38)       !!              - Print solver statistics in numsol 
(   39)       !!              - Stop the run IF problem for the solver ( indec < 0 )
(   40)       !!
(   41)       !! History :
(   42)       !!        !  91-03  ()
(   43)       !!        !  91-11  (G. Madec)
(   44)       !!        !  92-06  (M. Imbard)
(   45)       !!        !  97-06  (A.M. Treguier)
(   46)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   47)       !!----------------------------------------------------------------------
(   48)       !! * Arguments
(   49)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   50)       INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence
(   51) 
(   52)       !! * local declarations
(   53)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   54)       INTEGER  ::   ii, ij, ik              ! temporary integers
(   55)       REAL(wp) ::   zumax, zsmin            ! temporary scalars
(   56)       INTEGER, DIMENSION(3) ::   ilocu      ! 
(   57)       INTEGER, DIMENSION(2) ::   ilocs      ! 






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 6

(   58)       CHARACTER(len=80) :: clname
(   59)       !!----------------------------------------------------------------------
(   60)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   61)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/stpctl.F90,v 1.9 2006/03/21 07:57:56 opalod Exp $ 
(   62)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   63)       !!----------------------------------------------------------------------
(   64) 
(   65)       IF( kt == nit000 .AND. lwp ) THEN
(   66)          WRITE(numout,*)
(   67)          WRITE(numout,*) 'stp_ctl : time-stepping control'
(   68)          WRITE(numout,*) '~~~~~~~'
(   69)          ! open time.step file
(   70)          clname = 'time.step'
(   71)          CALL ctlopn( numstp, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
(   72)       ENDIF
(   73) 
(   74)       ! save the current time step in numstp
(   75)       ! ------------------------------------
(   76) !      IF(lwp) WRITE(numstp,9100) kt
(   77) !      IF(lwp) REWIND(numstp)
(   78) !!DB -- on some machines, the above does not write every timestep
(   79) !!      On drakes, the below made no difference.
(   80)       if(lwp) then
(   81)          WRITE(numstp,9100) kt
(   82)          call flush(numstp)
(   83)          REWIND(numstp)
(   84)       endif
(   85) 
(   86) 9100  FORMAT(1x, i8)
(   87) 
(   88) 
(   89)       ! elliptic solver statistics (if required)
(   90)       ! --------------------------
(   91)       IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
(   92)       ! Solver
(   93)       IF(lwp) WRITE(numsol,9200) kt, niter, res, SQRT(epsr)/eps
(   94) 
(   95)       ! Islands (if exist)
(   96)       IF( lk_isl )   CALL isl_stp_ctl( kt, kindic )
(   97) 
(   98) 
(   99)       ! Output in numwso and numwvo IF kindic<0
(  100)       ! ---------------------------------------
(  101)       !    (i.e. problem for the solver)
(  102)       IF( kindic < 0 ) THEN
(  103)          IF(lwp) THEN
(  104)             WRITE(numout,*) ' stpctl: the elliptic solver DO not converge or explode'
(  105)             WRITE(numout,*) ' ====== '
(  106)             WRITE(numout,9200) kt, niter, res, sqrt(epsr)/eps
(  107)             WRITE(numout,*)
(  108)             WRITE(numout,*) ' stpctl: output of last fields in numwso'
(  109)             WRITE(numout,*) '                                  numwvo'
(  110)             WRITE(numout,*) ' ======  *******************************'
(  111)          ENDIF
(  112)          CALL dia_wri( kt, kindic )
(  113)       ENDIF
(  114)       ENDIF
(  115) 






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 7

(  116) 9200  FORMAT(' it :', i8, ' niter :', i4, ' res :',e20.10,' b :',e20.10)
(  117) 
(  118)       ! Test maximum of velocity (zonal only)
(  119)       ! ------------------------
(  120)       !! zumax = MAXVAL( ABS( un(:,:,:) ) )   ! slower than the following loop on NEC SX5
(  121)       zumax = 0.e0
(  122)       DO jk = 1, jpk
(  123)          DO jj = 1, jpj
(  124)             DO ji = 1, jpi
(  125)                zumax = MAX(zumax,ABS(un(ji,jj,jk)))
(  126)           END DO 
(  127)         END DO 
(  128)       END DO        
(  129)       IF( lk_mpp )   CALL mpp_max( zumax )   ! max over the global domain
(  130) 
(  131)       IF( MOD( kt, nwrite ) == 1 ) THEN
(  132)          IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(U) max: ', zumax
(  133)       ENDIF
(  134)       IF( zumax > 20.) THEN
(  135)          IF( lk_mpp ) THEN
(  136)             CALL mpp_maxloc(ABS(un),umask,zumax,ii,ij,ik)
(  137)          ELSE
(  138)             ilocu = MAXLOC( ABS( un(:,:,:) ) )
(  139)             ii = ilocu(1) + nimpp - 1
(  140)             ij = ilocu(2) + njmpp - 1
(  141)             ik = ilocu(3)
(  142)          ENDIF
(  143)          IF(lwp) THEN
(  144)             WRITE(numout,cform_err)
(  145)             WRITE(numout,*) ' stpctl: the zonal velocity is larger than 20 m/s'
(  146)             WRITE(numout,*) ' ====== '
(  147)             WRITE(numout,9400) kt, zumax, ii, ij, ik
(  148)             WRITE(numout,*)
(  149)             WRITE(numout,*) '          output of last fields in numwso'
(  150)          ENDIF
(  151)          kindic  = -3
(  152) 
(  153)          CALL dia_wri( kt, kindic )
(  154)       ENDIF
(  155) 9400  FORMAT (' kt=',i6,' max abs(U): ',1pg11.4,', i j k: ',3i4)
(  156) 
(  157) 
(  158)       ! Test minimum of salinity
(  159)       ! ------------------------
(  160)       !! zsmin = MINVAL( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )    
(  161)       !                slower than the following loop on NEC SX5
(  162)       zsmin = 100.e0
(  163)       DO jj = 2, jpjm1
(  164)          DO ji = 1, jpi
(  165)             IF( tmask(ji,jj,1) == 1) zsmin = MIN(zsmin,sn(ji,jj,1))
(  166)          END DO
(  167)       END DO
(  168)       IF( lk_mpp )   CALL mpp_min( zsmin )   ! min over the global domain
(  169) 
(  170)       IF( MOD( kt, nwrite ) == 1 ) THEN
(  171)          IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' SSS min:', zsmin
(  172)       ENDIF
(  173)       IF( zsmin < 0.) THEN 






PGF90 (Version     10.5)          09/23/2010  09:50:10      page 8

(  174)          IF (lk_mpp) THEN
(  175)             CALL mpp_minloc ( sn(:,:,1),tmask(:,:,1), zsmin, ii,ij )
(  176)          ELSE
(  177)             ilocs = MINLOC( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )
(  178)             ii = ilocs(1) + nimpp - 1
(  179)             ij = ilocs(2) + njmpp - 1
(  180)          END IF
(  181) 
(  182)          IF(lwp) THEN
(  183)             WRITE(numout,cform_err)
(  184)             WRITE(numout,*) 'stp_ctl : NEGATIVE sea surface salinity'
(  185)             WRITE(numout,*) '======= '
(  186)             WRITE(numout,9500) kt, zsmin, ii, ij
(  187)             WRITE(numout,*)
(  188)             WRITE(numout,*) '          output of last fields in numwso'
(  189)          ENDIF
(  190)          IF( kindic < 0 ) THEN
(  191)             IF(lwp) WRITE(numout,*) ' stpctl diabort done. We wont do it again '
(  192)          ELSE 
(  193)             kindic  = -3
(  194)             CALL dia_wri(kt,kindic)
(  195)          ENDIF
(  196)       ENDIF
(  197) 9500  FORMAT (' kt=',i6,' min SSS: ',1pg11.4,', i j: ',2i4)
(  198) 
(  199)    END SUBROUTINE stp_ctl
(  200) 
(  201)    !!======================================================================
(  202) END MODULE stpctl
