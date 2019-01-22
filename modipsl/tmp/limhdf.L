


PGF90 (Version     10.5)          09/23/2010  09:49:54      page 1

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

Filename: limhdf.F90

(    1) MODULE limhdf
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE limhdf   ***
(    4)    !! LIM ice model : horizontal diffusion of sea-ice quantities
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_hdf  : diffusion trend on sea-ice variable
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE dom_oce
(   14)    USE ice_oce         ! ice variables
(   15)    USE in_out_manager
(   16)    USE ice
(   17)    USE lbclnk
(   18)    USE lib_mpp
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_hdf    ! called by lim_tra
(   26) 
(   27)    !! * Module variables
(   28)    LOGICAL  ::   linit = .TRUE.              ! ???
(   29)    REAL(wp) ::   epsi04 = 1e-04              ! constant
(   30)    REAL(wp), DIMENSION(jpi,jpj) ::   zfact   ! ???
(   31) 
(   32)    !! * Substitution 
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limhdf.F90,v 1.8 2005/09/22 13:50:30 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE lim_hdf( ptab )
(   43)       !!-------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE lim_hdf  ***
(   45)       !!
(   46)       !! ** purpose :   Compute and add the diffusive trend on sea-ice






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 2

(   47)       !!      variables
(   48)       !!
(   49)       !! ** method  :   Second order diffusive operator evaluated using a
(   50)       !!      Cranck-Nicholson time Scheme.
(   51)       !!
(   52)       !! ** Action  :    update ptab with the diffusive contribution
(   53)       !!
(   54)       !! History :
(   55)       !!        !  00-01 (LIM) Original code
(   56)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(   57)       !!        !  02-08 (C. Ethe)  F90, free form
(   58)       !!-------------------------------------------------------------------
(   59)       ! * Arguments
(   60)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
(   61)          ptab                 ! Field on which the diffusion is applied  
(   62)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   63)          ptab0                ! ???
(   64) 
(   65)       ! * Local variables
(   66)       INTEGER ::  ji, jj      ! dummy loop indices
(   67)       INTEGER ::  &
(   68)          its, iter            ! temporary integers
(   69)       CHARACTER (len=55) :: charout
(   70)       REAL(wp) ::  &
(   71)          zalfa, zrlxint, zconv, zeps   ! temporary scalars
(   72)       REAL(wp), DIMENSION(jpi,jpj) ::  & 
(   73)          zrlx, zflu, zflv, &  ! temporary workspaces
(   74)          zdiv0, zdiv          !    "           "
(   75)       !!-------------------------------------------------------------------
(   76) 
(   77)       ! Initialisation
(   78)       ! ---------------   
(   79)       ! Time integration parameters
(   80)       zalfa = 0.5       ! =1.0/0.5/0.0 = implicit/Cranck-Nicholson/explicit
(   81)       its   = 100       ! Maximum number of iteration
(   82)       zeps  =  2. * epsi04
(   83) 
(   84)       ! Arrays initialization
(   85)       ptab0 (:, : ) = ptab(:,:)
(   86) !bug  zflu (:,jpj) = 0.e0
(   87) !bug  zflv (:,jpj) = 0.e0
(   88)       zdiv0(:, 1 ) = 0.e0
(   89)       zdiv0(:,jpj) = 0.e0
(   90)       IF( .NOT.lk_vopt_loop ) THEN
(   91)          zflu (jpi,:) = 0.e0   
(   92)          zflv (jpi,:) = 0.e0
(   93)          zdiv0(1,  :) = 0.e0
(   94)          zdiv0(jpi,:) = 0.e0
(   95)       ENDIF
(   96) 
(   97)       ! Metric coefficient (compute at the first call and saved in
(   98)       IF( linit ) THEN
(   99)          DO jj = 2, jpjm1  
(  100)             DO ji = fs_2 , fs_jpim1   ! vector opt.
(  101)                zfact(ji,jj) = ( e2u(ji,jj) + e2u(ji-1,jj  ) + e1v(ji,jj) + e1v(ji,jj-1) ) &
(  102)                   &          / ( e1t(ji,jj) * e2t(ji,jj) )
(  103)             END DO
(  104)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 3

(  105)          linit = .FALSE.
(  106)       ENDIF
(  107) 
(  108) 
(  109)       ! Sub-time step loop
(  110)       zconv = 1.e0
(  111)       iter  = 0
(  112) 
(  113)       !                                                   !===================
(  114)       DO WHILE ( ( zconv > zeps ) .AND. (iter <= its) )   ! Sub-time step loop
(  115)          !                                                !===================
(  116)          ! incrementation of the sub-time step number
(  117)          iter = iter + 1
(  118) 
(  119)          ! diffusive fluxes in U- and V- direction
(  120)          DO jj = 1, jpjm1
(  121)             DO ji = 1 , fs_jpim1   ! vector opt.
(  122)                zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) / e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
(  123)                zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) / e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
(  124)             END DO
(  125)          END DO
(  126) 
(  127)          ! diffusive trend : divergence of the fluxes
(  128)          DO jj= 2, jpjm1
(  129)             DO ji = fs_2 , fs_jpim1   ! vector opt. 
(  130)                zdiv (ji,jj) = (  zflu(ji,jj) - zflu(ji-1,jj  )   &
(  131)                   &            + zflv(ji,jj) - zflv(ji  ,jj-1)  ) / ( e1t (ji,jj) * e2t (ji,jj) )
(  132)             END DO
(  133)          END DO
(  134) 
(  135)          ! save the first evaluation of the diffusive trend in zdiv0
(  136)          IF( iter == 1 )   zdiv0(:,:) = zdiv(:,:)       
(  137) 
(  138)          ! XXXX iterative evaluation?????
(  139)          DO jj = 2, jpjm1
(  140)             DO ji = fs_2 , fs_jpim1   ! vector opt.
(  141)                zrlxint = (   ptab0(ji,jj)    &
(  142)                   &       +  rdt_ice * (           zalfa   * ( zdiv(ji,jj) + zfact(ji,jj) * ptab(ji,jj) )   &
(  143)                   &                      + ( 1.0 - zalfa ) *   zdiv0(ji,jj) )  )                             & 
(  144)                   &    / ( 1.0 + zalfa * rdt_ice * zfact(ji,jj) )
(  145)                zrlx(ji,jj) = ptab(ji,jj) + om * ( zrlxint - ptab(ji,jj) )
(  146)             END DO
(  147)          END DO
(  148) 
(  149)          ! lateral boundary condition on ptab
(  150)          CALL lbc_lnk( zrlx, 'T', 1. )
(  151) 
(  152)          ! convergence test
(  153)          zconv = 0.e0
(  154)          DO jj = 2, jpjm1
(  155)             DO ji = 2, jpim1
(  156)                zconv = MAX( zconv, ABS( zrlx(ji,jj) - ptab(ji,jj) )  )
(  157)             END DO
(  158)          END DO
(  159)          IF( lk_mpp )   CALL mpp_max( zconv )   ! max over the global domain
(  160) 
(  161)          ptab(:,:) = zrlx(:,:)
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 4

(  163)          !                                         !==========================
(  164)       END DO                                       ! end of sub-time step loop
(  165)       !                                            !==========================
(  166) 
(  167)       IF(ln_ctl)   THEN
(  168)          zrlx(:,:) = ptab(:,:) - ptab0(:,:)
(  169)          WRITE(charout,FMT="(' lim_hdf  : zconv =',D23.16, ' iter =',I4,2X)") zconv, iter
(  170)          CALL prt_ctl(tab2d_1=zrlx, clinfo1=charout)
(  171)       ENDIF
(  172) 
(  173)    END SUBROUTINE lim_hdf
(  174) 
(  175) #else
(  176)    !!----------------------------------------------------------------------
(  177)    !!   Default option          Dummy module           NO LIM sea-ice model
(  178)    !!----------------------------------------------------------------------
(  179) CONTAINS
(  180)    SUBROUTINE lim_hdf         ! Empty routine
(  181)    END SUBROUTINE lim_hdf
(  182) #endif
(  183) 
(  184)    !!======================================================================
(  185) END MODULE limhdf









































PGF90 (Version     10.5)          09/23/2010  09:49:54      page 5

(    1) # 1 "limhdf.F90"
(    1) MODULE limhdf
(    2)    !!======================================================================
(    3)    !!                    ***  MODULE limhdf   ***
(    4)    !! LIM ice model : horizontal diffusion of sea-ice quantities
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_hdf  : diffusion trend on sea-ice variable
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE dom_oce
(   14)    USE ice_oce         ! ice variables
(   15)    USE in_out_manager
(   16)    USE ice
(   17)    USE lbclnk
(   18)    USE lib_mpp
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_hdf    ! called by lim_tra
(   26) 
(   27)    !! * Module variables
(   28)    LOGICAL  ::   linit = .TRUE.              ! ???
(   29)    REAL(wp) ::   epsi04 = 1e-04              ! constant
(   30)    REAL(wp), DIMENSION(jpi,jpj) ::   zfact   ! ???
(   31) 
(   32)    !! * Substitution 
(   33) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 34 "limhdf.F90"
(   34) # 34
(   34)    !!----------------------------------------------------------------------
(   35)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limhdf.F90,v 1.8 2005/09/22 13:50:30 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41) 
(   42)    SUBROUTINE lim_hdf( ptab )
(   43)       !!-------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE lim_hdf  ***






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 6

(   45)       !!
(   46)       !! ** purpose :   Compute and add the diffusive trend on sea-ice
(   47)       !!      variables
(   48)       !!
(   49)       !! ** method  :   Second order diffusive operator evaluated using a
(   50)       !!      Cranck-Nicholson time Scheme.
(   51)       !!
(   52)       !! ** Action  :    update ptab with the diffusive contribution
(   53)       !!
(   54)       !! History :
(   55)       !!        !  00-01 (LIM) Original code
(   56)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(   57)       !!        !  02-08 (C. Ethe)  F90, free form
(   58)       !!-------------------------------------------------------------------
(   59)       ! * Arguments
(   60)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
(   61)          ptab                 ! Field on which the diffusion is applied  
(   62)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   63)          ptab0                ! ???
(   64) 
(   65)       ! * Local variables
(   66)       INTEGER ::  ji, jj      ! dummy loop indices
(   67)       INTEGER ::  &
(   68)          its, iter            ! temporary integers
(   69)       CHARACTER (len=55) :: charout
(   70)       REAL(wp) ::  &
(   71)          zalfa, zrlxint, zconv, zeps   ! temporary scalars
(   72)       REAL(wp), DIMENSION(jpi,jpj) ::  & 
(   73)          zrlx, zflu, zflv, &  ! temporary workspaces
(   74)          zdiv0, zdiv          !    "           "
(   75)       !!-------------------------------------------------------------------
(   76) 
(   77)       ! Initialisation
(   78)       ! ---------------   
(   79)       ! Time integration parameters
(   80)       zalfa = 0.5       ! =1.0/0.5/0.0 = implicit/Cranck-Nicholson/explicit
(   81)       its   = 100       ! Maximum number of iteration
(   82)       zeps  =  2. * epsi04
(   83) 
(   84)       ! Arrays initialization
(   85)       ptab0 (:, : ) = ptab(:,:)
(   86) !bug  zflu (:,jpj) = 0.e0
(   87) !bug  zflv (:,jpj) = 0.e0
(   88)       zdiv0(:, 1 ) = 0.e0
(   89)       zdiv0(:,jpj) = 0.e0
(   90)       IF( .NOT.lk_vopt_loop ) THEN
(   91)          zflu (jpi,:) = 0.e0   
(   92)          zflv (jpi,:) = 0.e0
(   93)          zdiv0(1,  :) = 0.e0
(   94)          zdiv0(jpi,:) = 0.e0
(   95)       ENDIF
(   96) 
(   97)       ! Metric coefficient (compute at the first call and saved in
(   98)       IF( linit ) THEN
(   99)          DO jj = 2, jpjm1  
(  100)             DO ji = 2 , jpim1   ! vector opt.
(  101)                zfact(ji,jj) = ( e2u(ji,jj) + e2u(ji-1,jj  ) + e1v(ji,jj) + e1v(ji,jj-1) ) &
(  102)                   &          / ( e1t(ji,jj) * e2t(ji,jj) )






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 7

(  103)             END DO
(  104)          END DO
(  105)          linit = .FALSE.
(  106)       ENDIF
(  107) 
(  108) 
(  109)       ! Sub-time step loop
(  110)       zconv = 1.e0
(  111)       iter  = 0
(  112) 
(  113)       !                                                   !===================
(  114)       DO WHILE ( ( zconv > zeps ) .AND. (iter <= its) )   ! Sub-time step loop
(  115)          !                                                !===================
(  116)          ! incrementation of the sub-time step number
(  117)          iter = iter + 1
(  118) 
(  119)          ! diffusive fluxes in U- and V- direction
(  120)          DO jj = 1, jpjm1
(  121)             DO ji = 1 , jpim1   ! vector opt.
(  122)                zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) / e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
(  123)                zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) / e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
(  124)             END DO
(  125)          END DO
(  126) 
(  127)          ! diffusive trend : divergence of the fluxes
(  128)          DO jj= 2, jpjm1
(  129)             DO ji = 2 , jpim1   ! vector opt. 
(  130)                zdiv (ji,jj) = (  zflu(ji,jj) - zflu(ji-1,jj  )   &
(  131)                   &            + zflv(ji,jj) - zflv(ji  ,jj-1)  ) / ( e1t (ji,jj) * e2t (ji,jj) )
(  132)             END DO
(  133)          END DO
(  134) 
(  135)          ! save the first evaluation of the diffusive trend in zdiv0
(  136)          IF( iter == 1 )   zdiv0(:,:) = zdiv(:,:)       
(  137) 
(  138)          ! XXXX iterative evaluation?????
(  139)          DO jj = 2, jpjm1
(  140)             DO ji = 2 , jpim1   ! vector opt.
(  141)                zrlxint = (   ptab0(ji,jj)    &
(  142)                   &       +  rdt_ice * (           zalfa   * ( zdiv(ji,jj) + zfact(ji,jj) * ptab(ji,jj) )   &
(  143)                   &                      + ( 1.0 - zalfa ) *   zdiv0(ji,jj) )  )                             & 
(  144)                   &    / ( 1.0 + zalfa * rdt_ice * zfact(ji,jj) )
(  145)                zrlx(ji,jj) = ptab(ji,jj) + om * ( zrlxint - ptab(ji,jj) )
(  146)             END DO
(  147)          END DO
(  148) 
(  149)          ! lateral boundary condition on ptab
(  150)          CALL lbc_lnk( zrlx, 'T', 1. )
(  151) 
(  152)          ! convergence test
(  153)          zconv = 0.e0
(  154)          DO jj = 2, jpjm1
(  155)             DO ji = 2, jpim1
(  156)                zconv = MAX( zconv, ABS( zrlx(ji,jj) - ptab(ji,jj) )  )
(  157)             END DO
(  158)          END DO
(  159)          IF( lk_mpp )   CALL mpp_max( zconv )   ! max over the global domain
(  160) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 8

(  161)          ptab(:,:) = zrlx(:,:)
(  162) 
(  163)          !                                         !==========================
(  164)       END DO                                       ! end of sub-time step loop
(  165)       !                                            !==========================
(  166) 
(  167)       IF(ln_ctl)   THEN
(  168)          zrlx(:,:) = ptab(:,:) - ptab0(:,:)
(  169)          WRITE(charout,FMT="(' lim_hdf  : zconv =',D23.16, ' iter =',I4,2X)") zconv, iter
(  170)          CALL prt_ctl(tab2d_1=zrlx, clinfo1=charout)
(  171)       ENDIF
(  172) 
(  173)    END SUBROUTINE lim_hdf
(  174) 
(  175) 
(  176) # 184
(  184)    !!======================================================================
(  185) END MODULE limhdf
