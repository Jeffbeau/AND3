


PGF90 (Version     10.5)          09/23/2010  09:50:19      page 1

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

Filename: zdfbfr.F90

(    1) MODULE zdfbfr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfbfr  ***
(    4)    !! Ocean physics: Bottom friction
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_bfr      : update momentum Kz at the ocean bottom due to the
(    9)    !!                  type of bottom friction chosen
(   10)    !!   zdf_bfr_init : read in namelist and control the bottom friction
(   11)    !!                  parameters.
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables 
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zdf_bfr    ! called by step.F90
(   26) 
(   27)    !! * Module variables
(   28)    INTEGER ::             & !!! ** bottom friction namelist (nambfr) **
(   29)       nbotfr = 0             ! = 0/1/2/3 type of bottom friction 
(   30)    REAL(wp) ::            & !!! ** bottom friction namelist (nambfr) **
(   31)       bfri1 = 4.0e-4_wp,  &  ! bottom drag coefficient (linear case) 
(   32)       bfri2 = 1.0e-3_wp,  &  ! bottom drag coefficient (non linear case)
(   33)       bfeb2 = 2.5e-3_wp      ! background bottom turbulent kinetic energy  (m2/s2)
(   34) 
(   35)    !! * Substitutions
(   36) #  include "domzgr_substitute.h90"
(   37)    !!----------------------------------------------------------------------
(   38)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfbfr.F90,v 1.5 2005/09/02 15:45:43 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE zdf_bfr( kt )
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 2

(   47)       !!                   ***  ROUTINE zdf_bfr  ***
(   48)       !!                 
(   49)       !! ** Purpose :   Applied the bottom friction through a specification of 
(   50)       !!      Kz at the ocean bottom.
(   51)       !!
(   52)       !! ** Method  :   Update the value of avmu and avmv at the ocean bottom 
(   53)       !!       level following the chosen friction type (no-slip, free-slip, 
(   54)       !!       linear, or quadratic)
(   55)       !!
(   56)       !! History :
(   57)       !!   8.0  !  97-06  (G. Madec, A.-M. Treguier)  Original code
(   58)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   59)       !!----------------------------------------------------------------------
(   60)       !! * Arguments
(   61)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   62) 
(   63)       !! * Local declarations
(   64)       INTEGER ::   &
(   65)          ji, jj,                   &  ! dummy loop indexes
(   66)          ikbu, ikbv,               &  ! temporary integers
(   67)          ikbum1, ikbvm1               !
(   68)       REAL(wp) ::   &
(   69)          zvu, zuv, zecu, zecv         ! temporary scalars
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72) 
(   73)       IF( kt == nit000 )   CALL zdf_bfr_init
(   74) 
(   75) 
(   76)       ! Compute avmu, avmv at the ocean bottom
(   77)       ! --------------------------------------
(   78) 
(   79)       SELECT CASE (nbotfr)
(   80) 
(   81)       CASE( 0 )                 ! no-slip boundary condition
(   82) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(   83)          jj = 1
(   84)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(   85) # else
(   86)          DO jj = 2, jpjm1
(   87)             DO ji = 2, jpim1
(   88) # endif
(   89)                ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(   90)                ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(   91)                ikbum1 = MAX( ikbu-1, 1 )
(   92)                ikbvm1 = MAX( ikbv-1, 1 )
(   93)                avmu(ji,jj,ikbu) = 2. * avmu(ji,jj,ikbum1)
(   94)                avmv(ji,jj,ikbv) = 2. * avmv(ji,jj,ikbvm1)
(   95) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(   96)             END DO
(   97) # endif
(   98)          END DO
(   99) 
(  100)       CASE( 1 )                 ! linear botton friction
(  101) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  102)          jj = 1
(  103)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  104) # else






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 3

(  105)          DO jj = 2, jpjm1
(  106)             DO ji = 2, jpim1
(  107) # endif
(  108)                ikbu = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  109)                ikbv = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  110)                avmu(ji,jj,ikbu) = bfri1 * fse3uw(ji,jj,ikbu)
(  111)                avmv(ji,jj,ikbv) = bfri1 * fse3vw(ji,jj,ikbv)
(  112) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  113)             END DO
(  114) # endif
(  115)          END DO
(  116) 
(  117)       CASE( 2 )                 ! quadratic botton friction
(  118) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  119)          jj = 1
(  120) !CDIR NOVERRCHK
(  121)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  122) # else
(  123) !CDIR NOVERRCHK
(  124)          DO jj = 2, jpjm1
(  125) !CDIR NOVERRCHK
(  126)             DO ji = 2, jpim1
(  127) # endif
(  128)                ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(  129)                ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(  130)                ikbum1 = MAX( ikbu-1, 1 )
(  131)                ikbvm1 = MAX( ikbv-1, 1 )
(  132)                
(  133)                zvu  = 0.25 * (  vn(ji,jj  ,ikbum1) + vn(ji+1,jj  ,ikbum1)     &
(  134)                               + vn(ji,jj-1,ikbum1) + vn(ji+1,jj-1,ikbum1)  )
(  135)                
(  136)                zuv  = 0.25 * (  un(ji,jj  ,ikbvm1) + un(ji-1,jj  ,ikbvm1)     &
(  137)                               + un(ji,jj+1,ikbvm1) + un(ji-1,jj+1,ikbvm1)  )
(  138)                
(  139)                zecu = SQRT(  un(ji,jj,ikbum1) * un(ji,jj,ikbum1) + zvu*zvu + bfeb2  )
(  140)                zecv = SQRT(  vn(ji,jj,ikbvm1) * vn(ji,jj,ikbvm1) + zuv*zuv + bfeb2  )
(  141) 
(  142)                avmu(ji,jj,ikbu) = bfri2 * zecu * fse3uw(ji,jj,ikbu)
(  143)                avmv(ji,jj,ikbv) = bfri2 * zecv * fse3vw(ji,jj,ikbv)
(  144) 
(  145) !!DB: depth dependent bottom friction fudge for tidal run
(  146) !!  bfri2 = bfri2*(1 + F*e(-z_btm/D)^2); take D = 20m???; F=2???
(  147) !               avmu(ji,jj,ikbu) = bfri2 *(1.0+2.*exp(-(gdepw(ikbum1)/20.)**2)) &
(  148) !                                              * zecu * fse3uw(ji,jj,ikbu)
(  149) !               avmv(ji,jj,ikbv) = bfri2 *(1.0+2.*exp(-(gdepw(ikbvm1)/20.)**2)) &
(  150) !                                              * zecv * fse3vw(ji,jj,ikbv)
(  151) 
(  152) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  153)             END DO
(  154) # endif
(  155)          END DO
(  156) 
(  157)       CASE( 3 )                 ! free-slip boundary condition
(  158) # if defined key_vectopt_loop   &&   ! defined key_autotasking
(  159)          jj = 1
(  160)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  161) # else
(  162)          DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 4

(  163)             DO ji = 2, jpim1
(  164) # endif
(  165)                ikbu = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(  166)                ikbv = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(  167)                avmu(ji,jj,ikbu) = 0.e0
(  168)                avmv(ji,jj,ikbv) = 0.e0
(  169) # if ! defined key_vectopt_loop   ||   defined key_autotasking
(  170)             END DO
(  171) # endif
(  172)          END DO
(  173) 
(  174)       END SELECT
(  175) 
(  176)       ! Lateral boundary condition on (avmu,avmv)   (unchanged sign)
(  177)       ! ------------------------------===========
(  178)       CALL lbc_lnk( avmu, 'U', 1. )
(  179)       CALL lbc_lnk( avmv, 'V', 1. )
(  180) 
(  181)       IF(ln_ctl)   THEN
(  182)          CALL prt_ctl(tab3d_1=avmu, clinfo1=' bfr  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
(  183)       ENDIF
(  184) 
(  185)    END SUBROUTINE zdf_bfr
(  186) 
(  187) 
(  188)    SUBROUTINE zdf_bfr_init
(  189)       !!----------------------------------------------------------------------
(  190)       !!                  ***  ROUTINE zdf_bfr_init  ***
(  191)       !!                    
(  192)       !! ** Purpose :   Initialization of the bottom friction
(  193)       !!
(  194)       !! ** Method  :   Read the nammbf namelist and check their consistency
(  195)       !!      called at the first timestep (nit000)
(  196)       !!
(  197)       !! History :
(  198)       !!   9.0  !  02-06  (G. Madec)  Original code
(  199)       !!----------------------------------------------------------------------
(  200)       !! * Local declarations
(  201)       NAMELIST/nambfr/ nbotfr, bfri1, bfri2, bfeb2
(  202)       !!----------------------------------------------------------------------
(  203) 
(  204)       ! Read Namelist nambfr : bottom momentum boundary condition
(  205)       ! --------------------
(  206)       REWIND ( numnam )
(  207)       READ   ( numnam, nambfr )
(  208) 
(  209) 
(  210)       ! Parameter control and print
(  211)       ! ---------------------------
(  212)       IF(lwp) WRITE(numout,*)
(  213)       IF(lwp) WRITE(numout,*) 'zdf_bfr : momentum bottom friction'
(  214)       IF(lwp) WRITE(numout,*) '~~~~~~~'
(  215)       IF(lwp) WRITE(numout,*) '          Namelist nambfr : set bottom friction parameters'
(  216) 
(  217)       SELECT CASE (nbotfr)
(  218) 
(  219)       CASE( 0 )
(  220)          IF(lwp) WRITE(numout,*) '            no-slip '






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 5

(  221) 
(  222)       CASE( 1 )
(  223)          IF(lwp) WRITE(numout,*) '            linear botton friction'
(  224)          IF(lwp) WRITE(numout,*) '            friction coef.   bfri1  = ', bfri1
(  225) 
(  226)       CASE( 2 )
(  227)          IF(lwp) WRITE(numout,*) '            quadratic botton friction'
(  228)          IF(lwp) WRITE(numout,*) '            friction coef.   bfri2  = ', bfri2
(  229)          IF(lwp) WRITE(numout,*) '            background tke   bfeb2  = ', bfeb2
(  230) 
(  231)       CASE( 3 )
(  232)          IF(lwp) WRITE(numout,*) '            free-slip '
(  233) 
(  234)       CASE DEFAULT
(  235)          IF(lwp) WRITE(numout,cform_err)
(  236)          IF(lwp) WRITE(numout,*) '         bad flag value for nbotfr = ', nbotfr
(  237)          nstop = nstop + 1
(  238) 
(  239)       END SELECT
(  240) 
(  241)    END SUBROUTINE zdf_bfr_init
(  242) 
(  243)    !!======================================================================
(  244) END MODULE zdfbfr








































PGF90 (Version     10.5)          09/23/2010  09:50:19      page 6

(    1) # 1 "zdfbfr.F90"
(    1) MODULE zdfbfr
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfbfr  ***
(    4)    !! Ocean physics: Bottom friction
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_bfr      : update momentum Kz at the ocean bottom due to the
(    9)    !!                  type of bottom friction chosen
(   10)    !!   zdf_bfr_init : read in namelist and control the bottom friction
(   11)    !!                  parameters.
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables 
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zdf_bfr    ! called by step.F90
(   26) 
(   27)    !! * Module variables
(   28)    INTEGER ::             & !!! ** bottom friction namelist (nambfr) **
(   29)       nbotfr = 0             ! = 0/1/2/3 type of bottom friction 
(   30)    REAL(wp) ::            & !!! ** bottom friction namelist (nambfr) **
(   31)       bfri1 = 4.0e-4_wp,  &  ! bottom drag coefficient (linear case) 
(   32)       bfri2 = 1.0e-3_wp,  &  ! bottom drag coefficient (non linear case)
(   33)       bfeb2 = 2.5e-3_wp      ! background bottom turbulent kinetic energy  (m2/s2)
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
(   48) # 37 "zdfbfr.F90"
(   37) # 37
(   37)    !!----------------------------------------------------------------------
(   38)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfbfr.F90,v 1.5 2005/09/02 15:45:43 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 7

(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE zdf_bfr( kt )
(   46)       !!----------------------------------------------------------------------
(   47)       !!                   ***  ROUTINE zdf_bfr  ***
(   48)       !!                 
(   49)       !! ** Purpose :   Applied the bottom friction through a specification of 
(   50)       !!      Kz at the ocean bottom.
(   51)       !!
(   52)       !! ** Method  :   Update the value of avmu and avmv at the ocean bottom 
(   53)       !!       level following the chosen friction type (no-slip, free-slip, 
(   54)       !!       linear, or quadratic)
(   55)       !!
(   56)       !! History :
(   57)       !!   8.0  !  97-06  (G. Madec, A.-M. Treguier)  Original code
(   58)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   59)       !!----------------------------------------------------------------------
(   60)       !! * Arguments
(   61)       INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
(   62) 
(   63)       !! * Local declarations
(   64)       INTEGER ::   &
(   65)          ji, jj,                   &  ! dummy loop indexes
(   66)          ikbu, ikbv,               &  ! temporary integers
(   67)          ikbum1, ikbvm1               !
(   68)       REAL(wp) ::   &
(   69)          zvu, zuv, zecu, zecv         ! temporary scalars
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72) 
(   73)       IF( kt == nit000 )   CALL zdf_bfr_init
(   74) 
(   75) 
(   76)       ! Compute avmu, avmv at the ocean bottom
(   77)       ! --------------------------------------
(   78) 
(   79)       SELECT CASE (nbotfr)
(   80) 
(   81)       CASE( 0 )                 ! no-slip boundary condition
(   82) # 86
(   86)          DO jj = 2, jpjm1
(   87)             DO ji = 2, jpim1
(   88) # 89
(   89)                ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(   90)                ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(   91)                ikbum1 = MAX( ikbu-1, 1 )
(   92)                ikbvm1 = MAX( ikbv-1, 1 )
(   93)                avmu(ji,jj,ikbu) = 2. * avmu(ji,jj,ikbum1)
(   94)                avmv(ji,jj,ikbv) = 2. * avmv(ji,jj,ikbvm1)
(   95) # 96
(   96)             END DO
(   97) # 98
(   98)          END DO
(   99) 
(  100)       CASE( 1 )                 ! linear botton friction
(  101) # 105
(  105)          DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 8

(  106)             DO ji = 2, jpim1
(  107) # 108
(  108)                ikbu = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  109)                ikbv = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  110)                avmu(ji,jj,ikbu) = bfri1 * e3w(ikbu)
(  111)                avmv(ji,jj,ikbv) = bfri1 * e3w(ikbv)
(  112) # 113
(  113)             END DO
(  114) # 115
(  115)          END DO
(  116) 
(  117)       CASE( 2 )                 ! quadratic botton friction
(  118) # 123
(  123) !CDIR NOVERRCHK
(  124)          DO jj = 2, jpjm1
(  125) !CDIR NOVERRCHK
(  126)             DO ji = 2, jpim1
(  127) # 128
(  128)                ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(  129)                ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(  130)                ikbum1 = MAX( ikbu-1, 1 )
(  131)                ikbvm1 = MAX( ikbv-1, 1 )
(  132)                
(  133)                zvu  = 0.25 * (  vn(ji,jj  ,ikbum1) + vn(ji+1,jj  ,ikbum1)     &
(  134)                               + vn(ji,jj-1,ikbum1) + vn(ji+1,jj-1,ikbum1)  )
(  135)                
(  136)                zuv  = 0.25 * (  un(ji,jj  ,ikbvm1) + un(ji-1,jj  ,ikbvm1)     &
(  137)                               + un(ji,jj+1,ikbvm1) + un(ji-1,jj+1,ikbvm1)  )
(  138)                
(  139)                zecu = SQRT(  un(ji,jj,ikbum1) * un(ji,jj,ikbum1) + zvu*zvu + bfeb2  )
(  140)                zecv = SQRT(  vn(ji,jj,ikbvm1) * vn(ji,jj,ikbvm1) + zuv*zuv + bfeb2  )
(  141) 
(  142)                avmu(ji,jj,ikbu) = bfri2 * zecu * e3w(ikbu)
(  143)                avmv(ji,jj,ikbv) = bfri2 * zecv * e3w(ikbv)
(  144) 
(  145) !!DB: depth dependent bottom friction fudge for tidal run
(  146) !!  bfri2 = bfri2*(1 + F*e(-z_btm/D)^2); take D = 20m???; F=2???
(  147) !               avmu(ji,jj,ikbu) = bfri2 *(1.0+2.*exp(-(gdepw(ikbum1)/20.)**2)) &
(  148) !                                              * zecu * e3w(ikbu)
(  149) !               avmv(ji,jj,ikbv) = bfri2 *(1.0+2.*exp(-(gdepw(ikbvm1)/20.)**2)) &
(  150) !                                              * zecv * e3w(ikbv)
(  151) 
(  152) # 153
(  153)             END DO
(  154) # 155
(  155)          END DO
(  156) 
(  157)       CASE( 3 )                 ! free-slip boundary condition
(  158) # 162
(  162)          DO jj = 2, jpjm1
(  163)             DO ji = 2, jpim1
(  164) # 165
(  165)                ikbu = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
(  166)                ikbv = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
(  167)                avmu(ji,jj,ikbu) = 0.e0
(  168)                avmv(ji,jj,ikbv) = 0.e0
(  169) # 170
(  170)             END DO






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 9

(  171) # 172
(  172)          END DO
(  173) 
(  174)       END SELECT
(  175) 
(  176)       ! Lateral boundary condition on (avmu,avmv)   (unchanged sign)
(  177)       ! ------------------------------===========
(  178)       CALL lbc_lnk( avmu, 'U', 1. )
(  179)       CALL lbc_lnk( avmv, 'V', 1. )
(  180) 
(  181)       IF(ln_ctl)   THEN
(  182)          CALL prt_ctl(tab3d_1=avmu, clinfo1=' bfr  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
(  183)       ENDIF
(  184) 
(  185)    END SUBROUTINE zdf_bfr
(  186) 
(  187) 
(  188)    SUBROUTINE zdf_bfr_init
(  189)       !!----------------------------------------------------------------------
(  190)       !!                  ***  ROUTINE zdf_bfr_init  ***
(  191)       !!                    
(  192)       !! ** Purpose :   Initialization of the bottom friction
(  193)       !!
(  194)       !! ** Method  :   Read the nammbf namelist and check their consistency
(  195)       !!      called at the first timestep (nit000)
(  196)       !!
(  197)       !! History :
(  198)       !!   9.0  !  02-06  (G. Madec)  Original code
(  199)       !!----------------------------------------------------------------------
(  200)       !! * Local declarations
(  201)       NAMELIST/nambfr/ nbotfr, bfri1, bfri2, bfeb2
(  202)       !!----------------------------------------------------------------------
(  203) 
(  204)       ! Read Namelist nambfr : bottom momentum boundary condition
(  205)       ! --------------------
(  206)       REWIND ( numnam )
(  207)       READ   ( numnam, nambfr )
(  208) 
(  209) 
(  210)       ! Parameter control and print
(  211)       ! ---------------------------
(  212)       IF(lwp) WRITE(numout,*)
(  213)       IF(lwp) WRITE(numout,*) 'zdf_bfr : momentum bottom friction'
(  214)       IF(lwp) WRITE(numout,*) '~~~~~~~'
(  215)       IF(lwp) WRITE(numout,*) '          Namelist nambfr : set bottom friction parameters'
(  216) 
(  217)       SELECT CASE (nbotfr)
(  218) 
(  219)       CASE( 0 )
(  220)          IF(lwp) WRITE(numout,*) '            no-slip '
(  221) 
(  222)       CASE( 1 )
(  223)          IF(lwp) WRITE(numout,*) '            linear botton friction'
(  224)          IF(lwp) WRITE(numout,*) '            friction coef.   bfri1  = ', bfri1
(  225) 
(  226)       CASE( 2 )
(  227)          IF(lwp) WRITE(numout,*) '            quadratic botton friction'
(  228)          IF(lwp) WRITE(numout,*) '            friction coef.   bfri2  = ', bfri2






PGF90 (Version     10.5)          09/23/2010  09:50:19      page 10

(  229)          IF(lwp) WRITE(numout,*) '            background tke   bfeb2  = ', bfeb2
(  230) 
(  231)       CASE( 3 )
(  232)          IF(lwp) WRITE(numout,*) '            free-slip '
(  233) 
(  234)       CASE DEFAULT
(  235)          IF(lwp) WRITE(numout,cform_err)
(  236)          IF(lwp) WRITE(numout,*) '         bad flag value for nbotfr = ', nbotfr
(  237)          nstop = nstop + 1
(  238) 
(  239)       END SELECT
(  240) 
(  241)    END SUBROUTINE zdf_bfr_init
(  242) 
(  243)    !!======================================================================
(  244) END MODULE zdfbfr
