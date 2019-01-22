


PGF90 (Version     10.5)          09/23/2010  09:50:17      page 1

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

Filename: trasbc.F90

(    1) MODULE trasbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trasbc  ***
(    4)    !! Ocean active tracers:  surface boundary condition
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_sbc      : update the tracer trend at ocean surface
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and active tracers
(   12)    USE dom_oce         ! ocean space domain variables
(   13)    USE ocesbc          ! surface thermohaline fluxes
(   14)    USE phycst          ! physical constant
(   15)    USE traqsr          ! solar radiation penetration
(   16)    USE trdmod          ! ocean trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC tra_sbc              ! routine called by step.F90
(   26) 
(   27)    !! * Substitutions
(   28) #  include "domzgr_substitute.h90"
(   29) #  include "vectopt_loop_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trasbc.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE tra_sbc ( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                  ***  ROUTINE tra_sbc  ***
(   41)       !!                   
(   42)       !! ** Purpose :   Compute the tracer surface boundary condition trend of
(   43)       !!      (flux through the interface, concentration/dilution effect)
(   44)       !!      and add it to the general trend of tracer equations.
(   45)       !!
(   46)       !! ** Method :






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 2

(   47)       !!      * flux through the air-sea interface:
(   48)       !!            - temperature : heat flux q (w/m2). If penetrative solar
(   49)       !!         radiation q is only the non solar part of the heat flux, the
(   50)       !!         solar part is added in traqsr.F routine.
(   51)       !!            ta = ta + q /(rau0 rcp e3t)  for k=1
(   52)       !!            - salinity    : no salt flux
(   53)       !!      * concentration/dilution effect:
(   54)       !!            The surface freshwater flux modify the ocean volume
(   55)       !!         and thus the concentration of a tracer and the temperature.
(   56)       !!         First order of the effect of surface freshwater exchange 
(   57)       !!         for salinity, it can be neglected on temperature (especially
(   58)       !!         as the temparature of precipitations and runoffs is usually
(   59)       !!         unknown.
(   60)       !!            - temperature : we assume that the temperature of both
(   61)       !!         precipitations and runoffs is equal to the SST, thus there
(   62)       !!         is no additional flux since in this case, the concentration
(   63)       !!         dilution effect is balanced by the net heat flux associated
(   64)       !!         to the freshwater exchange:
(   65)       !!            (Tp P - Te E) + STT (P-E) = 0 when Tp=Te=SST
(   66)       !!            - salinity    : evaporation, precipitation and runoff
(   67)       !!         water has a zero salinity, thus
(   68)       !!            sa = sa + emp * sn / e3t   for k=1
(   69)       !!         where emp, the surface freshwater budget (evaporation minus
(   70)       !!         precipitation minus runoff) given in kg/m2/s is divided
(   71)       !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
(   72)       !!
(   73)       !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
(   74)       !!                with the tracer surface boundary condition 
(   75)       !!              - save the trend it in ttrd ('key_trdtra')
(   76)       !!
(   77)       !! History :
(   78)       !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
(   79)       !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
(   80)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdta => ua,      & ! use ua as 3D workspace   
(   85)                          ztdsa => va         ! use va as 3D workspace   
(   86) 
(   87)       !! * Arguments
(   88)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER  ::   ji, jj                   ! dummy loop indices
(   92)       REAL(wp) ::   zta, zsa, zsrau, zse3t   ! temporary scalars
(   93)       !!----------------------------------------------------------------------
(   94) 
(   95)       IF( kt == nit000 ) THEN
(   96)          IF(lwp) WRITE(numout,*)
(   97)          IF(lwp) WRITE(numout,*) 'tra_sbc : TRAcer Surface Boundary Condition'
(   98)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   99)       ENDIF
(  100) 
(  101)       ! 0. initialization
(  102)       zsrau = 1. / rauw
(  103) #if ! defined key_s_coord
(  104)       zse3t = 1. / fse3t(1,1,1)






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 3

(  105) #endif
(  106) 
(  107)       ! Save ta and sa trends
(  108)       IF( l_trdtra )   THEN
(  109)          ztdta(:,:,:) = ta(:,:,:) 
(  110)          ztdsa(:,:,:) = sa(:,:,:) 
(  111)       ENDIF
(  112) 
(  113)       IF( .NOT.ln_traqsr )   qsr(:,:) = 0.e0   ! no solar radiation penetration
(  114) 
(  115)       ! 1. Concentration dillution effect on (t,s)
(  116)       DO jj = 2, jpj
(  117)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  118) #if defined key_s_coord
(  119)             zse3t = 1. / fse3t(ji,jj,1)
(  120) #endif
(  121)             ! temperature : heat flux
(  122)             zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t
(  123) 
(  124)             ! salinity :  concent./dilut. effect
(  125)             zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t
(  126)             
(  127)             ! add the trend to the general tracer trend
(  128)             ta(ji,jj,1) = ta(ji,jj,1) + zta
(  129)             sa(ji,jj,1) = sa(ji,jj,1) + zsa
(  130)          END DO
(  131)       END DO
(  132) 
(  133)       ! save the trends for diagnostic
(  134)       ! sea surface boundary condition tracers trends
(  135)       IF( l_trdtra )   THEN
(  136)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  137)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  138)          CALL trd_mod(ztdta, ztdsa, jpttdnsr, 'TRA', kt)
(  139)       ENDIF
(  140)       
(  141)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  142)          CALL prt_ctl(tab3d_1=ta, clinfo1=' sbc  - Ta: ', mask1=tmask, &
(  143)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  144)       ENDIF
(  145) 
(  146)    END SUBROUTINE tra_sbc
(  147) 
(  148)    !!======================================================================
(  149) END MODULE trasbc



















PGF90 (Version     10.5)          09/23/2010  09:50:17      page 4

(    1) # 1 "trasbc.F90"
(    1) MODULE trasbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trasbc  ***
(    4)    !! Ocean active tracers:  surface boundary condition
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_sbc      : update the tracer trend at ocean surface
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and active tracers
(   12)    USE dom_oce         ! ocean space domain variables
(   13)    USE ocesbc          ! surface thermohaline fluxes
(   14)    USE phycst          ! physical constant
(   15)    USE traqsr          ! solar radiation penetration
(   16)    USE trdmod          ! ocean trends 
(   17)    USE trdmod_oce      ! ocean variables trends
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC tra_sbc              ! routine called by step.F90
(   26) 
(   27)    !! * Substitutions
(   28) # 1 "./domzgr_substitute.h90"
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
(   48) # 29 "trasbc.F90"
(   29) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 30 "trasbc.F90"
(   30) # 30






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 5

(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trasbc.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE tra_sbc ( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                  ***  ROUTINE tra_sbc  ***
(   41)       !!                   
(   42)       !! ** Purpose :   Compute the tracer surface boundary condition trend of
(   43)       !!      (flux through the interface, concentration/dilution effect)
(   44)       !!      and add it to the general trend of tracer equations.
(   45)       !!
(   46)       !! ** Method :
(   47)       !!      * flux through the air-sea interface:
(   48)       !!            - temperature : heat flux q (w/m2). If penetrative solar
(   49)       !!         radiation q is only the non solar part of the heat flux, the
(   50)       !!         solar part is added in traqsr.F routine.
(   51)       !!            ta = ta + q /(rau0 rcp e3t)  for k=1
(   52)       !!            - salinity    : no salt flux
(   53)       !!      * concentration/dilution effect:
(   54)       !!            The surface freshwater flux modify the ocean volume
(   55)       !!         and thus the concentration of a tracer and the temperature.
(   56)       !!         First order of the effect of surface freshwater exchange 
(   57)       !!         for salinity, it can be neglected on temperature (especially
(   58)       !!         as the temparature of precipitations and runoffs is usually
(   59)       !!         unknown.
(   60)       !!            - temperature : we assume that the temperature of both
(   61)       !!         precipitations and runoffs is equal to the SST, thus there
(   62)       !!         is no additional flux since in this case, the concentration
(   63)       !!         dilution effect is balanced by the net heat flux associated
(   64)       !!         to the freshwater exchange:
(   65)       !!            (Tp P - Te E) + STT (P-E) = 0 when Tp=Te=SST
(   66)       !!            - salinity    : evaporation, precipitation and runoff
(   67)       !!         water has a zero salinity, thus
(   68)       !!            sa = sa + emp * sn / e3t   for k=1
(   69)       !!         where emp, the surface freshwater budget (evaporation minus
(   70)       !!         precipitation minus runoff) given in kg/m2/s is divided
(   71)       !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
(   72)       !!
(   73)       !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
(   74)       !!                with the tracer surface boundary condition 
(   75)       !!              - save the trend it in ttrd ('key_trdtra')
(   76)       !!
(   77)       !! History :
(   78)       !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
(   79)       !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
(   80)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   81)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   82)       !!----------------------------------------------------------------------
(   83)       !! * Modules used     
(   84)       USE oce, ONLY :    ztdta => ua,      & ! use ua as 3D workspace   
(   85)                          ztdsa => va         ! use va as 3D workspace   
(   86) 
(   87)       !! * Arguments






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 6

(   88)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(   89) 
(   90)       !! * Local declarations
(   91)       INTEGER  ::   ji, jj                   ! dummy loop indices
(   92)       REAL(wp) ::   zta, zsa, zsrau, zse3t   ! temporary scalars
(   93)       !!----------------------------------------------------------------------
(   94) 
(   95)       IF( kt == nit000 ) THEN
(   96)          IF(lwp) WRITE(numout,*)
(   97)          IF(lwp) WRITE(numout,*) 'tra_sbc : TRAcer Surface Boundary Condition'
(   98)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   99)       ENDIF
(  100) 
(  101)       ! 0. initialization
(  102)       zsrau = 1. / rauw
(  103) # 104
(  104)       zse3t = 1. / e3t(1)
(  105) 
(  106) # 107
(  107)       ! Save ta and sa trends
(  108)       IF( l_trdtra )   THEN
(  109)          ztdta(:,:,:) = ta(:,:,:) 
(  110)          ztdsa(:,:,:) = sa(:,:,:) 
(  111)       ENDIF
(  112) 
(  113)       IF( .NOT.ln_traqsr )   qsr(:,:) = 0.e0   ! no solar radiation penetration
(  114) 
(  115)       ! 1. Concentration dillution effect on (t,s)
(  116)       DO jj = 2, jpj
(  117)          DO ji = 2, jpim1   ! vector opt.
(  118) # 121
(  121)             ! temperature : heat flux
(  122)             zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t
(  123) 
(  124)             ! salinity :  concent./dilut. effect
(  125)             zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t
(  126)             
(  127)             ! add the trend to the general tracer trend
(  128)             ta(ji,jj,1) = ta(ji,jj,1) + zta
(  129)             sa(ji,jj,1) = sa(ji,jj,1) + zsa
(  130)          END DO
(  131)       END DO
(  132) 
(  133)       ! save the trends for diagnostic
(  134)       ! sea surface boundary condition tracers trends
(  135)       IF( l_trdtra )   THEN
(  136)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  137)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  138)          CALL trd_mod(ztdta, ztdsa, jpttdnsr, 'TRA', kt)
(  139)       ENDIF
(  140)       
(  141)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  142)          CALL prt_ctl(tab3d_1=ta, clinfo1=' sbc  - Ta: ', mask1=tmask, &
(  143)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  144)       ENDIF
(  145) 
(  146)    END SUBROUTINE tra_sbc
(  147) 






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 7

(  148)    !!======================================================================
(  149) END MODULE trasbc
