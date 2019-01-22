


PGF90 (Version     10.5)          09/23/2010  09:49:25      page 1

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

Filename: obcvol.F90

(    1) MODULE obcvol
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcvol  ***
(    4)    !! Ocean dynamic :  Volume constraint when OBC and Free surface are used
(    5)    !!=================================================================================
(    6) #if   defined key_obc   &&   ! defined key_dynspg_rl
(    7)    !!---------------------------------------------------------------------------------
(    8)    !!   'key_obc'               and                           open boundary conditions
(    9)    !!   'key_dynspg_flt'                                  constant volume free surface
(   10)    !!---------------------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers 
(   13)    USE dom_oce         ! ocean space and time domain 
(   14)    USE phycst          ! physical constants
(   15)    USE obc_oce         ! ocean open boundary conditions
(   16)    USE lib_mpp         ! for mppsum
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ocesbc          ! ocean surface boundary conditions
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC obc_vol        ! routine called by dynspg_flt
(   25) 
(   26)    !! * Substitutions
(   27) #  include "domzgr_substitute.h90"
(   28) #  include "obc_vectopt_loop_substitute.h90"
(   29)    !!---------------------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcvol.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!---------------------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE obc_vol ( kt )
(   38)       !!------------------------------------------------------------------------------
(   39)       !!                      ***  ROUTINE obcvol  ***
(   40)       !!
(   41)       !! ** Purpose : 
(   42)       !!      This routine is called in dynspg_flt to control 
(   43)       !!      the volume of the system. A correction velocity is calculated
(   44)       !!      to correct the total transport through the OBC. 
(   45)       !!      The total depth used is constant (H0) to be consistent with the 
(   46)       !!      linear free surface coded in OPA 8.2






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 2

(   47)       !!
(   48)       !! ** Method :  
(   49)       !!      The correction velocity (zubtpecor here) is defined calculating
(   50)       !!      the total transport through all open boundaries (trans_obc) minus
(   51)       !!      the cumulate E-P flux (zCflxemp) divided by the total lateral 
(   52)       !!      surface (obcsurftot) of these OBC. 
(   53)       !!
(   54)       !!      zubtpecor = [trans_obc - zCflxemp ]*(1./obcsurftot)
(   55)       !!
(   56)       !!      with zCflxemp => sum of (Evaporation minus Precipitation)
(   57)       !!                       over all the domain in m3/s at each time step.
(   58)       !!
(   59)       !!      zCflxemp < 0 when precipitation dominate
(   60)       !!      zCflxemp > 0 when evaporation dominate
(   61)       !!
(   62)       !!      There are 2 options (user's desiderata): 
(   63)       !!
(   64)       !!         1/ The volume changes according to E-P, this is the default
(   65)       !!            option. In this case the cumulate E-P flux are setting to
(   66)       !!            zero (zCflxemp=0) to calculate the correction velocity. So
(   67)       !!            it will only balance the flux through open boundaries.
(   68)       !!            (set volemp to 0 in tne namelist for this option)
(   69)       !!
(   70)       !!         2/ The volume is constant even with E-P flux. In this case
(   71)       !!            the correction velocity must balance both the flux 
(   72)       !!            through open boundaries and the ones through the free
(   73)       !!            surface. 
(   74)       !!            (set volemp to 1 in tne namelist for this option)
(   75)       !!
(   76)       !! History :
(   77)       !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Original code
(   78)       !!----------------------------------------------------------------------------
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk
(   84)       REAL(wp) ::   zubtpecor
(   85)       REAL(wp) ::   zCflxemp
(   86)       REAL(wp) ::   ztransw, ztranse, ztransn, ztranss, ztranst
(   87)       !!-----------------------------------------------------------------------------
(   88) 
(   89)       IF( kt == nit000 ) THEN 
(   90)          IF(lwp) WRITE(numout,*)'        '
(   91)          IF(lwp) WRITE(numout,*)'obc_vol : Correction of velocities along OBC'
(   92)          IF(lwp) WRITE(numout,*)'~~~~~~~'
(   93)          IF(lwp) WRITE(numout,*)'        '
(   94)       END IF 
(   95) 
(   96)       ! 1. Calculate the cumulate surface Flux zCflxemp (m3/s) over all the domain.
(   97)       ! ---------------------------------------------------------------------------
(   98)  
(   99)       zCflxemp = 0.e0
(  100) 
(  101)       DO jj = 1, jpj
(  102)          DO ji = 1, jpi
(  103)             zCflxemp = zCflxemp + ( (emp(ji,jj)*obctmsk(ji,jj) )/rauw)*e1v(ji,jj)*e2u(ji,jj)
(  104)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 3

(  105)       END DO
(  106)       IF( lk_mpp )   CALL mpp_sum( zCflxemp )   ! sum over the global domain
(  107) 
(  108)       ! 2. Barotropic velocity for each open boundary
(  109)       ! ---------------------------------------------
(  110) 
(  111)       zubtpecor = 0.e0
(  112) 
(  113)       ! ... West open boundary
(  114)       IF( lp_obc_west ) THEN                      ! ... Total transport through the West OBC
(  115)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  116)             DO jk = 1, jpkm1
(  117)                DO jj = 1, jpj
(  118)                   zubtpecor = zubtpecor + ua(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk) * uwmsk(jj,jk)
(  119)                END DO
(  120)             END DO
(  121)          END DO
(  122)       END IF 
(  123) 
(  124)       ! ... East open boundary
(  125)       IF( lp_obc_east ) THEN                      ! ... Total transport through the East OBC
(  126)          DO ji = fs_nie0, fs_nie1 ! Vector opt.
(  127)             DO jk = 1, jpkm1
(  128)                DO jj = 1, jpj
(  129)                   zubtpecor = zubtpecor - ua(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk) * uemsk(jj,jk)
(  130)                END DO
(  131)             END DO
(  132)          END DO
(  133)       END IF 
(  134) 
(  135)       ! ... North open boundary
(  136)       IF( lp_obc_north ) THEN                     ! ... Total transport through the North OBC
(  137)          DO jj = fs_njn0, fs_njn1 ! Vector opt.
(  138)             DO jk = 1, jpkm1
(  139)                DO ji = 1, jpi
(  140)                   zubtpecor = zubtpecor - va(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk) * vnmsk(ji,jk)
(  141)                END DO
(  142)             END DO
(  143)          END DO
(  144)       END IF 
(  145) 
(  146)       ! ... South open boundary
(  147)       IF( lp_obc_south ) THEN                     ! ... Total transport through the South OBC
(  148)          DO jj = fs_njs0, fs_njs1 ! Vector opt.
(  149)             DO jk = 1, jpkm1
(  150)                DO ji = 1, jpi
(  151)                   zubtpecor = zubtpecor + va(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk) * vsmsk(ji,jk)
(  152)                END DO
(  153)             END DO
(  154)          END DO
(  155)       END IF 
(  156) 
(  157)       IF( lk_mpp )   CALL mpp_sum( zubtpecor )   ! sum over the global domain
(  158) 
(  159) 
(  160)       ! 3. The normal velocity correction
(  161)       ! ---------------------------------
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 4

(  163)       zubtpecor = (zubtpecor - zCflxemp*volemp)*(1./obcsurftot)
(  164) 
(  165)       IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  166)          IF(lwp) WRITE(numout,*)'        '
(  167)          IF(lwp) WRITE(numout,*)'obc_vol : time step :', kt
(  168)          IF(lwp) WRITE(numout,*)'~~~~~~~ '
(  169)          IF(lwp) WRITE(numout,*)'        '
(  170)          IF(lwp) WRITE(numout,*)'          cumulate flux EMP :', zCflxemp,' (m3/s)'
(  171)          IF(lwp) WRITE(numout,*)'          total lateral surface of OBC :',obcsurftot,'(m2)'
(  172)          IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor :',zubtpecor,'(m/s)'
(  173)          IF(lwp) WRITE(numout,*)'        '
(  174)       END IF 
(  175) 
(  176)       ! 4. Correction of the total velocity on each open 
(  177)       !    boundary torespect the mass flux conservation
(  178)       ! -------------------------------------------------
(  179) 
(  180)       ztransw = 0.e0
(  181)       ztranse = 0.e0
(  182)       ztransn = 0.e0
(  183)       ztranss = 0.e0
(  184)       ztranst = 0.e0
(  185) 
(  186)       IF( lp_obc_west ) THEN
(  187) 
(  188)          ! ... correction of the west velocity
(  189)          DO ji = fs_niw0, fs_niw1 ! Vector opt.
(  190)             DO jk = 1, jpkm1
(  191)                DO jj = 1, jpj
(  192)                   ua(ji,jj,jk) = ua(ji,jj,jk) - zubtpecor*uwmsk(jj,jk)
(  193)                   ztransw= ztransw + ua(ji,jj,jk)*fse3u(ji,jj,jk)*e2u(ji,jj)*uwmsk(jj,jk)
(  194)                END DO
(  195)             END DO
(  196)          END DO
(  197) 
(  198)          IF( lk_mpp )   CALL mpp_sum( ztransw )   ! sum over the global domain
(  199) 
(  200)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  201)             IF(lwp) WRITE(numout,*)'          West OB transport ztransw :', ztransw,'(m3/s)'
(  202)          END IF 
(  203) 
(  204)       END IF 
(  205) 
(  206)       IF( lp_obc_east ) THEN
(  207) 
(  208)          ! ... correction of the east velocity
(  209)          DO ji = fs_nie0, fs_nie1 ! Vector opt.
(  210)             DO jk = 1, jpkm1
(  211)                DO jj = 1, jpj
(  212)                   ua(ji,jj,jk) = ua(ji,jj,jk) + zubtpecor*uemsk(jj,jk)
(  213)                   ztranse= ztranse + ua(ji,jj,jk)*fse3u(ji,jj,jk)*e2u(ji,jj)*uemsk(jj,jk)
(  214)                END DO
(  215)             END DO
(  216)          END DO
(  217) 
(  218)          IF( lk_mpp )   CALL mpp_sum( ztranse )   ! sum over the global domain
(  219) 
(  220)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 5

(  221)             IF(lwp) WRITE(numout,*)'          East OB transport ztranse :', ztranse,'(m3/s)'
(  222)          END IF 
(  223) 
(  224)       END IF 
(  225) 
(  226)       IF( lp_obc_north ) THEN
(  227) 
(  228)          ! ... correction of the north velocity
(  229)          DO jj = fs_njn0, fs_njn1 ! Vector opt.
(  230)             DO jk = 1, jpkm1
(  231)                DO ji =  1, jpi
(  232)                   va(ji,jj,jk) = va(ji,jj,jk) + zubtpecor*vnmsk(ji,jk)
(  233)                   ztransn= ztransn + va(ji,jj,jk)*fse3v(ji,jj,jk)*e1v(ji,jj)*vnmsk(ji,jk)
(  234)                END DO
(  235)             END DO
(  236)          END DO
(  237)          IF( lk_mpp )   CALL mpp_sum( ztransn )   ! sum over the global domain
(  238) 
(  239)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  240)             IF(lwp) WRITE(numout,*)'          North OB transport ztransn :', ztransn,'(m3/s)'
(  241)          END IF 
(  242) 
(  243)       END IF 
(  244) 
(  245)       IF( lp_obc_south ) THEN
(  246) 
(  247)          ! ... correction of the south velocity
(  248)          DO jj = fs_njs0, fs_njs1 ! Vector opt.
(  249)             DO jk = 1, jpkm1
(  250)                DO ji =  1, jpi
(  251)                   va(ji,jj,jk) = va(ji,jj,jk) - zubtpecor*vsmsk(ji,jk)
(  252)                   ztranss= ztranss + va(ji,jj,jk)*fse3v(ji,jj,jk)*e1v(ji,jj)*vsmsk(ji,jk)
(  253)                END DO
(  254)             END DO
(  255)          END DO
(  256)          IF( lk_mpp )   CALL mpp_sum( ztranss )   ! sum over the global domain
(  257) 
(  258)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  259)             IF(lwp) WRITE(numout,*)'          South OB transport ztranss :', ztranss,'(m3/s)'
(  260)          END IF 
(  261) 
(  262)       END IF 
(  263) 
(  264)       ! 5. Check the cumulate transport through OBC
(  265)       !    once barotropic velocities corrected
(  266)       ! -------------------------------------------
(  267) 
(  268)       ztranst = ztransw - ztranse + ztranss - ztransn
(  269) 
(  270)       IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  271)          IF(lwp) WRITE(numout,*)'        '
(  272)          IF(lwp) WRITE(numout,*)'          Cumulate transport ztranst =', ztranst,'(m3/s)'
(  273)          IF(lwp) WRITE(numout,*)'        '
(  274)       END IF 
(  275) 
(  276)    END SUBROUTINE obc_vol
(  277) 
(  278) #else






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 6

(  279)    !!---------------------------------------------------------------------------------
(  280)    !!  Default option :                                                   Empty module
(  281)    !!---------------------------------------------------------------------------------
(  282) CONTAINS
(  283)    SUBROUTINE obc_vol        ! Empty routine
(  284)    END SUBROUTINE obc_vol
(  285) #endif
(  286) 
(  287)    !!=================================================================================
(  288) END MODULE obcvol






















































PGF90 (Version     10.5)          09/23/2010  09:49:25      page 7

(    1) # 1 "obcvol.F90"
(    1) MODULE obcvol
(    2)    !!=================================================================================
(    3)    !!                       ***  MODULE  obcvol  ***
(    4)    !! Ocean dynamic :  Volume constraint when OBC and Free surface are used
(    5)    !!=================================================================================
(    6) # 7
(    7)    !!---------------------------------------------------------------------------------
(    8)    !!   'key_obc'               and                           open boundary conditions
(    9)    !!   'key_dynspg_flt'                                  constant volume free surface
(   10)    !!---------------------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers 
(   13)    USE dom_oce         ! ocean space and time domain 
(   14)    USE phycst          ! physical constants
(   15)    USE obc_oce         ! ocean open boundary conditions
(   16)    USE lib_mpp         ! for mppsum
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ocesbc          ! ocean surface boundary conditions
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC obc_vol        ! routine called by dynspg_flt
(   25) 
(   26)    !! * Substitutions
(   27) # 1 "./domzgr_substitute.h90"
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
(   48) # 28 "obcvol.F90"
(   28) # 1 "./obc_vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                 ***  obc_vectopt_loop_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_vectopt_loop_substitute.h90,v 1.3 2005/12/12 14:20:26 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 29 "obcvol.F90"
(   29) # 29
(   29)    !!---------------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 8

(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcvol.F90,v 1.5 2005/12/28 09:25:08 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!---------------------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE obc_vol ( kt )
(   38)       !!------------------------------------------------------------------------------
(   39)       !!                      ***  ROUTINE obcvol  ***
(   40)       !!
(   41)       !! ** Purpose : 
(   42)       !!      This routine is called in dynspg_flt to control 
(   43)       !!      the volume of the system. A correction velocity is calculated
(   44)       !!      to correct the total transport through the OBC. 
(   45)       !!      The total depth used is constant (H0) to be consistent with the 
(   46)       !!      linear free surface coded in OPA 8.2
(   47)       !!
(   48)       !! ** Method :  
(   49)       !!      The correction velocity (zubtpecor here) is defined calculating
(   50)       !!      the total transport through all open boundaries (trans_obc) minus
(   51)       !!      the cumulate E-P flux (zCflxemp) divided by the total lateral 
(   52)       !!      surface (obcsurftot) of these OBC. 
(   53)       !!
(   54)       !!      zubtpecor = [trans_obc - zCflxemp ]*(1./obcsurftot)
(   55)       !!
(   56)       !!      with zCflxemp => sum of (Evaporation minus Precipitation)
(   57)       !!                       over all the domain in m3/s at each time step.
(   58)       !!
(   59)       !!      zCflxemp < 0 when precipitation dominate
(   60)       !!      zCflxemp > 0 when evaporation dominate
(   61)       !!
(   62)       !!      There are 2 options (user's desiderata): 
(   63)       !!
(   64)       !!         1/ The volume changes according to E-P, this is the default
(   65)       !!            option. In this case the cumulate E-P flux are setting to
(   66)       !!            zero (zCflxemp=0) to calculate the correction velocity. So
(   67)       !!            it will only balance the flux through open boundaries.
(   68)       !!            (set volemp to 0 in tne namelist for this option)
(   69)       !!
(   70)       !!         2/ The volume is constant even with E-P flux. In this case
(   71)       !!            the correction velocity must balance both the flux 
(   72)       !!            through open boundaries and the ones through the free
(   73)       !!            surface. 
(   74)       !!            (set volemp to 1 in tne namelist for this option)
(   75)       !!
(   76)       !! History :
(   77)       !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Original code
(   78)       !!----------------------------------------------------------------------------
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   81) 
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk
(   84)       REAL(wp) ::   zubtpecor
(   85)       REAL(wp) ::   zCflxemp
(   86)       REAL(wp) ::   ztransw, ztranse, ztransn, ztranss, ztranst
(   87)       !!-----------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 9

(   88) 
(   89)       IF( kt == nit000 ) THEN 
(   90)          IF(lwp) WRITE(numout,*)'        '
(   91)          IF(lwp) WRITE(numout,*)'obc_vol : Correction of velocities along OBC'
(   92)          IF(lwp) WRITE(numout,*)'~~~~~~~'
(   93)          IF(lwp) WRITE(numout,*)'        '
(   94)       END IF 
(   95) 
(   96)       ! 1. Calculate the cumulate surface Flux zCflxemp (m3/s) over all the domain.
(   97)       ! ---------------------------------------------------------------------------
(   98)  
(   99)       zCflxemp = 0.e0
(  100) 
(  101)       DO jj = 1, jpj
(  102)          DO ji = 1, jpi
(  103)             zCflxemp = zCflxemp + ( (emp(ji,jj)*obctmsk(ji,jj) )/rauw)*e1v(ji,jj)*e2u(ji,jj)
(  104)          END DO
(  105)       END DO
(  106)       IF( lk_mpp )   CALL mpp_sum( zCflxemp )   ! sum over the global domain
(  107) 
(  108)       ! 2. Barotropic velocity for each open boundary
(  109)       ! ---------------------------------------------
(  110) 
(  111)       zubtpecor = 0.e0
(  112) 
(  113)       ! ... West open boundary
(  114)       IF( lp_obc_west ) THEN                      ! ... Total transport through the West OBC
(  115)          DO ji = niw0 , niw1 ! Vector opt.
(  116)             DO jk = 1, jpkm1
(  117)                DO jj = 1, jpj
(  118)                   zubtpecor = zubtpecor + ua(ji,jj,jk)*e2u(ji,jj)*e3t(jk) * uwmsk(jj,jk)
(  119)                END DO
(  120)             END DO
(  121)          END DO
(  122)       END IF 
(  123) 
(  124)       ! ... East open boundary
(  125)       IF( lp_obc_east ) THEN                      ! ... Total transport through the East OBC
(  126)          DO ji = nie0 , nie1 ! Vector opt.
(  127)             DO jk = 1, jpkm1
(  128)                DO jj = 1, jpj
(  129)                   zubtpecor = zubtpecor - ua(ji,jj,jk)*e2u(ji,jj)*e3t(jk) * uemsk(jj,jk)
(  130)                END DO
(  131)             END DO
(  132)          END DO
(  133)       END IF 
(  134) 
(  135)       ! ... North open boundary
(  136)       IF( lp_obc_north ) THEN                     ! ... Total transport through the North OBC
(  137)          DO jj = njn0, njn1 ! Vector opt.
(  138)             DO jk = 1, jpkm1
(  139)                DO ji = 1, jpi
(  140)                   zubtpecor = zubtpecor - va(ji,jj,jk)*e1v(ji,jj)*e3t(jk) * vnmsk(ji,jk)
(  141)                END DO
(  142)             END DO
(  143)          END DO
(  144)       END IF 
(  145) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 10

(  146)       ! ... South open boundary
(  147)       IF( lp_obc_south ) THEN                     ! ... Total transport through the South OBC
(  148)          DO jj = njs0, njs1 ! Vector opt.
(  149)             DO jk = 1, jpkm1
(  150)                DO ji = 1, jpi
(  151)                   zubtpecor = zubtpecor + va(ji,jj,jk)*e1v(ji,jj)*e3t(jk) * vsmsk(ji,jk)
(  152)                END DO
(  153)             END DO
(  154)          END DO
(  155)       END IF 
(  156) 
(  157)       IF( lk_mpp )   CALL mpp_sum( zubtpecor )   ! sum over the global domain
(  158) 
(  159) 
(  160)       ! 3. The normal velocity correction
(  161)       ! ---------------------------------
(  162) 
(  163)       zubtpecor = (zubtpecor - zCflxemp*volemp)*(1./obcsurftot)
(  164) 
(  165)       IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  166)          IF(lwp) WRITE(numout,*)'        '
(  167)          IF(lwp) WRITE(numout,*)'obc_vol : time step :', kt
(  168)          IF(lwp) WRITE(numout,*)'~~~~~~~ '
(  169)          IF(lwp) WRITE(numout,*)'        '
(  170)          IF(lwp) WRITE(numout,*)'          cumulate flux EMP :', zCflxemp,' (m3/s)'
(  171)          IF(lwp) WRITE(numout,*)'          total lateral surface of OBC :',obcsurftot,'(m2)'
(  172)          IF(lwp) WRITE(numout,*)'          correction velocity zubtpecor :',zubtpecor,'(m/s)'
(  173)          IF(lwp) WRITE(numout,*)'        '
(  174)       END IF 
(  175) 
(  176)       ! 4. Correction of the total velocity on each open 
(  177)       !    boundary torespect the mass flux conservation
(  178)       ! -------------------------------------------------
(  179) 
(  180)       ztransw = 0.e0
(  181)       ztranse = 0.e0
(  182)       ztransn = 0.e0
(  183)       ztranss = 0.e0
(  184)       ztranst = 0.e0
(  185) 
(  186)       IF( lp_obc_west ) THEN
(  187) 
(  188)          ! ... correction of the west velocity
(  189)          DO ji = niw0 , niw1 ! Vector opt.
(  190)             DO jk = 1, jpkm1
(  191)                DO jj = 1, jpj
(  192)                   ua(ji,jj,jk) = ua(ji,jj,jk) - zubtpecor*uwmsk(jj,jk)
(  193)                   ztransw= ztransw + ua(ji,jj,jk)*e3t(jk)*e2u(ji,jj)*uwmsk(jj,jk)
(  194)                END DO
(  195)             END DO
(  196)          END DO
(  197) 
(  198)          IF( lk_mpp )   CALL mpp_sum( ztransw )   ! sum over the global domain
(  199) 
(  200)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  201)             IF(lwp) WRITE(numout,*)'          West OB transport ztransw :', ztransw,'(m3/s)'
(  202)          END IF 
(  203) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 11

(  204)       END IF 
(  205) 
(  206)       IF( lp_obc_east ) THEN
(  207) 
(  208)          ! ... correction of the east velocity
(  209)          DO ji = nie0 , nie1 ! Vector opt.
(  210)             DO jk = 1, jpkm1
(  211)                DO jj = 1, jpj
(  212)                   ua(ji,jj,jk) = ua(ji,jj,jk) + zubtpecor*uemsk(jj,jk)
(  213)                   ztranse= ztranse + ua(ji,jj,jk)*e3t(jk)*e2u(ji,jj)*uemsk(jj,jk)
(  214)                END DO
(  215)             END DO
(  216)          END DO
(  217) 
(  218)          IF( lk_mpp )   CALL mpp_sum( ztranse )   ! sum over the global domain
(  219) 
(  220)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  221)             IF(lwp) WRITE(numout,*)'          East OB transport ztranse :', ztranse,'(m3/s)'
(  222)          END IF 
(  223) 
(  224)       END IF 
(  225) 
(  226)       IF( lp_obc_north ) THEN
(  227) 
(  228)          ! ... correction of the north velocity
(  229)          DO jj = njn0, njn1 ! Vector opt.
(  230)             DO jk = 1, jpkm1
(  231)                DO ji =  1, jpi
(  232)                   va(ji,jj,jk) = va(ji,jj,jk) + zubtpecor*vnmsk(ji,jk)
(  233)                   ztransn= ztransn + va(ji,jj,jk)*e3t(jk)*e1v(ji,jj)*vnmsk(ji,jk)
(  234)                END DO
(  235)             END DO
(  236)          END DO
(  237)          IF( lk_mpp )   CALL mpp_sum( ztransn )   ! sum over the global domain
(  238) 
(  239)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  240)             IF(lwp) WRITE(numout,*)'          North OB transport ztransn :', ztransn,'(m3/s)'
(  241)          END IF 
(  242) 
(  243)       END IF 
(  244) 
(  245)       IF( lp_obc_south ) THEN
(  246) 
(  247)          ! ... correction of the south velocity
(  248)          DO jj = njs0, njs1 ! Vector opt.
(  249)             DO jk = 1, jpkm1
(  250)                DO ji =  1, jpi
(  251)                   va(ji,jj,jk) = va(ji,jj,jk) - zubtpecor*vsmsk(ji,jk)
(  252)                   ztranss= ztranss + va(ji,jj,jk)*e3t(jk)*e1v(ji,jj)*vsmsk(ji,jk)
(  253)                END DO
(  254)             END DO
(  255)          END DO
(  256)          IF( lk_mpp )   CALL mpp_sum( ztranss )   ! sum over the global domain
(  257) 
(  258)          IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  259)             IF(lwp) WRITE(numout,*)'          South OB transport ztranss :', ztranss,'(m3/s)'
(  260)          END IF 
(  261) 






PGF90 (Version     10.5)          09/23/2010  09:49:25      page 12

(  262)       END IF 
(  263) 
(  264)       ! 5. Check the cumulate transport through OBC
(  265)       !    once barotropic velocities corrected
(  266)       ! -------------------------------------------
(  267) 
(  268)       ztranst = ztransw - ztranse + ztranss - ztransn
(  269) 
(  270)       IF( lwp .AND. MOD( kt, nwrite ) == 0) THEN
(  271)          IF(lwp) WRITE(numout,*)'        '
(  272)          IF(lwp) WRITE(numout,*)'          Cumulate transport ztranst =', ztranst,'(m3/s)'
(  273)          IF(lwp) WRITE(numout,*)'        '
(  274)       END IF 
(  275) 
(  276)    END SUBROUTINE obc_vol
(  277) 
(  278) 
(  279) # 287
(  287)    !!=================================================================================
(  288) END MODULE obcvol
