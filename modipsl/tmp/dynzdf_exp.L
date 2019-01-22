


PGF90 (Version     10.5)          09/23/2010  09:49:29      page 1

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

Filename: dynzdf_exp.F90

(    1) MODULE dynzdf_exp
(    2)    !!==============================================================================
(    3)    !!                     ***  MODULE  dynzdf_exp  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_exp  : update the momentum trend with the vertical diffu-
(    9)    !!                  sion using an explicit time-stepping scheme.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE phycst          ! physical constants
(   15)    USE zdf_oce         ! ocean vertical physics
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE taumod          ! surface ocean stress
(   18)    USE trdmod          ! ocean dynamics trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dyn_zdf_exp    ! called by step.F90
(   27) 
(   28)    !! * Substitutions
(   29) #  include "domzgr_substitute.h90"
(   30) #  include "vectopt_loop_substitute.h90"
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_exp.F90,v 1.5 2005/09/02 15:45:24 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE dyn_zdf_exp( kt )
(   40)       !!----------------------------------------------------------------------
(   41)       !!                  ***  ROUTINE dyn_zdf_exp  ***
(   42)       !!                   
(   43)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   44)       !!
(   45)       !! ** Method  :   Explicit forward time stepping with a time splitting
(   46)       !!      technique. The vertical diffusion of momentum is given by:






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 2

(   47)       !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
(   48)       !!      Surface boundary conditions: wind stress input
(   49)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
(   50)       !!      Add this trend to the general trend ua :
(   51)       !!         ua = ua + dz( avmu dz(u) )
(   52)       !!
(   53)       !! ** Action : - Update (ua,va) with the vertical diffusive trend
(   54)       !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
(   55)       !!
(   56)       !! History :
(   57)       !!        !  90-10  (B. Blanke)  Original code
(   58)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   59)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   60)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   61)       !!---------------------------------------------------------------------
(   62)       !! * Modules used     
(   63)       USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
(   64)                          ztdva => sa       ! use sa as 3D workspace   
(   65)       !! * Arguments
(   66)       INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
(   67) 
(   68)       !! * Local declarations
(   69)       INTEGER ::   &
(   70)          ji, jj, jk, jl,                 & ! dummy loop indices
(   71)          ikbu, ikbum1 , ikbv, ikbvm1       ! temporary integers
(   72)       REAL(wp) ::   &
(   73)          zrau0r, zlavmr, z2dt, zua, zva    ! temporary scalars
(   74)       REAL(wp), DIMENSION(jpi,jpk) ::    &
(   75)          zwx, zwy, zwz, zww                ! temporary workspace arrays
(   76)       REAL(wp), DIMENSION(jpi,jpj) ::    &
(   77)          ztsx, ztsy, ztbx, ztby            ! temporary workspace arrays
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80)       IF( kt == nit000 ) THEN
(   81)          IF(lwp) WRITE(numout,*)
(   82)          IF(lwp) WRITE(numout,*) 'dyn_zdf_exp : vertical momentum diffusion explicit operator'
(   83)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   84)       ENDIF
(   85) 
(   86)       ! Local constant initialization
(   87)       ! -----------------------------
(   88)       zrau0r = 1. / rau0                                   ! inverse of the reference density
(   89)       zlavmr = 1. / float( n_zdfexp )                      ! inverse of the number of sub time step
(   90)       z2dt = 2. * rdt                                      ! Leap-frog environnement
(   91)       ztsx(:,:) = 0.e0
(   92)       ztsy(:,:) = 0.e0 
(   93)       ztbx(:,:) = 0.e0
(   94)       ztby(:,:) = 0.e0
(   95) 
(   96)       ! Save ua and va trends
(   97)       IF( l_trddyn )   THEN
(   98)          ztdua(:,:,:) = ua(:,:,:) 
(   99)          ztdva(:,:,:) = va(:,:,:) 
(  100)       ENDIF
(  101) 
(  102)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt    ! Euler time stepping when starting from rest
(  103) 
(  104)       !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 3

(  105)       DO jj = 2, jpjm1                                 !  Vertical slab
(  106)          !                                             ! ===============
(  107) 
(  108)          ! Surface boundary condition
(  109)          DO ji = 2, jpim1
(  110)             zwy(ji,1) = taux(ji,jj) * zrau0r
(  111)             zww(ji,1) = tauy(ji,jj) * zrau0r
(  112)          END DO  
(  113) 
(  114)          ! Initialization of x, z and contingently trends array
(  115)          DO jk = 1, jpk
(  116)             DO ji = 2, jpim1
(  117)                zwx(ji,jk) = ub(ji,jj,jk)
(  118)                zwz(ji,jk) = vb(ji,jj,jk)
(  119)             END DO  
(  120)          END DO  
(  121) 
(  122)          ! Time splitting loop
(  123)          DO jl = 1, n_zdfexp
(  124) 
(  125)             ! First vertical derivative
(  126)             DO jk = 2, jpk
(  127)                DO ji = 2, jpim1
(  128)                   zwy(ji,jk) = avmu(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) / fse3uw(ji,jj,jk) 
(  129)                   zww(ji,jk) = avmv(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) / fse3vw(ji,jj,jk)
(  130)                END DO  
(  131)             END DO  
(  132) 
(  133)             ! Second vertical derivative and trend estimation at kt+l*rdt/n_zdfexp
(  134)             DO jk = 1, jpkm1
(  135)                DO ji = 2, jpim1
(  136)                   zua = zlavmr*( zwy(ji,jk) - zwy(ji,jk+1) ) / fse3u(ji,jj,jk)
(  137)                   zva = zlavmr*( zww(ji,jk) - zww(ji,jk+1) ) / fse3v(ji,jj,jk)
(  138)                   ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  139)                   va(ji,jj,jk) = va(ji,jj,jk) + zva
(  140) 
(  141)                   zwx(ji,jk) = zwx(ji,jk) + z2dt*zua*umask(ji,jj,jk)
(  142)                   zwz(ji,jk) = zwz(ji,jk) + z2dt*zva*vmask(ji,jj,jk)
(  143)                END DO  
(  144)             END DO  
(  145) 
(  146)          END DO  
(  147) 
(  148)          !                                             ! ===============
(  149)       END DO                                           !   End of slab
(  150)       !                                                ! ===============
(  151) 
(  152)       ! save the vertical diffusion trends for diagnostic
(  153)       ! momentum trends
(  154)       IF( l_trddyn )  THEN 
(  155)          ! save the total vertical momentum diffusive trend
(  156)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  157)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  158)  
(  159)          ! subtract and save surface and momentum fluxes
(  160)          !                                                ! ===============
(  161)          DO jj = 2, jpjm1                                 !  Horizontal slab
(  162)             !                                             ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 4

(  163)             DO ji = 2, jpim1
(  164)                ! save the surface momentum fluxes 
(  165)                ztsx(ji,jj) = zwy(ji,1) / fse3u(ji,jj,1)
(  166)                ztsy(ji,jj) = zww(ji,1) / fse3v(ji,jj,1)
(  167)                ! save bottom friction momentum fluxes 
(  168)                ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  169)                ikbum1 = MAX( ikbu-1, 1 )
(  170)                ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  171)                ikbvm1 = MAX( ikbv-1, 1 )
(  172)                ztbx(ji,jj) = avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  173)                                / ( fse3u(ji,jj,ikbum1) * fse3uw(ji,jj,ikbu) )
(  174)                ztby(ji,jj) = avmv(ji,jj,ikbv) * zwz(ji,ikbvm1)   &
(  175)                                / ( fse3v(ji,jj,ikbvm1) * fse3vw(ji,jj,ikbv) )
(  176)                ! subtract surface forcing and bottom friction trend from vertical
(  177)                ! diffusive momentum trend
(  178)                ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  179)                ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  180)                ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  181)                ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  182)             END DO
(  183)             !                                             ! ===============
(  184)          END DO                                           !   End of slab
(  185)          !                                                ! ===============
(  186) 
(  187)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  188)          ztdua(:,:,:) = 0.e0
(  189)          ztdva(:,:,:) = 0.e0
(  190)          ztdua(:,:,1) = ztsx(:,:)
(  191)          ztdva(:,:,1) = ztsy(:,:)
(  192)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  193)          ztdua(:,:,:) = 0.e0
(  194)          ztdva(:,:,:) = 0.e0
(  195)          ztdua(:,:,1) = ztbx(:,:)
(  196)          ztdva(:,:,1) = ztby(:,:)
(  197)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  198)       ENDIF
(  199) 
(  200)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  201)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  202)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  203)       ENDIF
(  204) 
(  205)    END SUBROUTINE dyn_zdf_exp
(  206) 
(  207)    !!==============================================================================
(  208) END MODULE dynzdf_exp


















PGF90 (Version     10.5)          09/23/2010  09:49:29      page 5

(    1) # 1 "dynzdf_exp.F90"
(    1) MODULE dynzdf_exp
(    2)    !!==============================================================================
(    3)    !!                     ***  MODULE  dynzdf_exp  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_exp  : update the momentum trend with the vertical diffu-
(    9)    !!                  sion using an explicit time-stepping scheme.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE phycst          ! physical constants
(   15)    USE zdf_oce         ! ocean vertical physics
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE taumod          ! surface ocean stress
(   18)    USE trdmod          ! ocean dynamics trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dyn_zdf_exp    ! called by step.F90
(   27) 
(   28)    !! * Substitutions
(   29) # 1 "./domzgr_substitute.h90"
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
(   48) # 30 "dynzdf_exp.F90"
(   30) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 31 "dynzdf_exp.F90"






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 6

(   31) # 31
(   31)    !!----------------------------------------------------------------------
(   32)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_exp.F90,v 1.5 2005/09/02 15:45:24 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)    SUBROUTINE dyn_zdf_exp( kt )
(   40)       !!----------------------------------------------------------------------
(   41)       !!                  ***  ROUTINE dyn_zdf_exp  ***
(   42)       !!                   
(   43)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   44)       !!
(   45)       !! ** Method  :   Explicit forward time stepping with a time splitting
(   46)       !!      technique. The vertical diffusion of momentum is given by:
(   47)       !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
(   48)       !!      Surface boundary conditions: wind stress input
(   49)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
(   50)       !!      Add this trend to the general trend ua :
(   51)       !!         ua = ua + dz( avmu dz(u) )
(   52)       !!
(   53)       !! ** Action : - Update (ua,va) with the vertical diffusive trend
(   54)       !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
(   55)       !!
(   56)       !! History :
(   57)       !!        !  90-10  (B. Blanke)  Original code
(   58)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   59)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   60)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   61)       !!---------------------------------------------------------------------
(   62)       !! * Modules used     
(   63)       USE oce, ONLY :    ztdua => ta,    & ! use ta as 3D workspace   
(   64)                          ztdva => sa       ! use sa as 3D workspace   
(   65)       !! * Arguments
(   66)       INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
(   67) 
(   68)       !! * Local declarations
(   69)       INTEGER ::   &
(   70)          ji, jj, jk, jl,                 & ! dummy loop indices
(   71)          ikbu, ikbum1 , ikbv, ikbvm1       ! temporary integers
(   72)       REAL(wp) ::   &
(   73)          zrau0r, zlavmr, z2dt, zua, zva    ! temporary scalars
(   74)       REAL(wp), DIMENSION(jpi,jpk) ::    &
(   75)          zwx, zwy, zwz, zww                ! temporary workspace arrays
(   76)       REAL(wp), DIMENSION(jpi,jpj) ::    &
(   77)          ztsx, ztsy, ztbx, ztby            ! temporary workspace arrays
(   78)       !!----------------------------------------------------------------------
(   79) 
(   80)       IF( kt == nit000 ) THEN
(   81)          IF(lwp) WRITE(numout,*)
(   82)          IF(lwp) WRITE(numout,*) 'dyn_zdf_exp : vertical momentum diffusion explicit operator'
(   83)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   84)       ENDIF
(   85) 
(   86)       ! Local constant initialization
(   87)       ! -----------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 7

(   88)       zrau0r = 1. / rau0                                   ! inverse of the reference density
(   89)       zlavmr = 1. / float( n_zdfexp )                      ! inverse of the number of sub time step
(   90)       z2dt = 2. * rdt                                      ! Leap-frog environnement
(   91)       ztsx(:,:) = 0.e0
(   92)       ztsy(:,:) = 0.e0 
(   93)       ztbx(:,:) = 0.e0
(   94)       ztby(:,:) = 0.e0
(   95) 
(   96)       ! Save ua and va trends
(   97)       IF( l_trddyn )   THEN
(   98)          ztdua(:,:,:) = ua(:,:,:) 
(   99)          ztdva(:,:,:) = va(:,:,:) 
(  100)       ENDIF
(  101) 
(  102)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt    ! Euler time stepping when starting from rest
(  103) 
(  104)       !                                                ! ===============
(  105)       DO jj = 2, jpjm1                                 !  Vertical slab
(  106)          !                                             ! ===============
(  107) 
(  108)          ! Surface boundary condition
(  109)          DO ji = 2, jpim1
(  110)             zwy(ji,1) = taux(ji,jj) * zrau0r
(  111)             zww(ji,1) = tauy(ji,jj) * zrau0r
(  112)          END DO  
(  113) 
(  114)          ! Initialization of x, z and contingently trends array
(  115)          DO jk = 1, jpk
(  116)             DO ji = 2, jpim1
(  117)                zwx(ji,jk) = ub(ji,jj,jk)
(  118)                zwz(ji,jk) = vb(ji,jj,jk)
(  119)             END DO  
(  120)          END DO  
(  121) 
(  122)          ! Time splitting loop
(  123)          DO jl = 1, n_zdfexp
(  124) 
(  125)             ! First vertical derivative
(  126)             DO jk = 2, jpk
(  127)                DO ji = 2, jpim1
(  128)                   zwy(ji,jk) = avmu(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) / e3w(jk) 
(  129)                   zww(ji,jk) = avmv(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) / e3w(jk)
(  130)                END DO  
(  131)             END DO  
(  132) 
(  133)             ! Second vertical derivative and trend estimation at kt+l*rdt/n_zdfexp
(  134)             DO jk = 1, jpkm1
(  135)                DO ji = 2, jpim1
(  136)                   zua = zlavmr*( zwy(ji,jk) - zwy(ji,jk+1) ) / e3t(jk)
(  137)                   zva = zlavmr*( zww(ji,jk) - zww(ji,jk+1) ) / e3t(jk)
(  138)                   ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  139)                   va(ji,jj,jk) = va(ji,jj,jk) + zva
(  140) 
(  141)                   zwx(ji,jk) = zwx(ji,jk) + z2dt*zua*umask(ji,jj,jk)
(  142)                   zwz(ji,jk) = zwz(ji,jk) + z2dt*zva*vmask(ji,jj,jk)
(  143)                END DO  
(  144)             END DO  
(  145) 






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 8

(  146)          END DO  
(  147) 
(  148)          !                                             ! ===============
(  149)       END DO                                           !   End of slab
(  150)       !                                                ! ===============
(  151) 
(  152)       ! save the vertical diffusion trends for diagnostic
(  153)       ! momentum trends
(  154)       IF( l_trddyn )  THEN 
(  155)          ! save the total vertical momentum diffusive trend
(  156)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  157)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  158)  
(  159)          ! subtract and save surface and momentum fluxes
(  160)          !                                                ! ===============
(  161)          DO jj = 2, jpjm1                                 !  Horizontal slab
(  162)             !                                             ! ===============
(  163)             DO ji = 2, jpim1
(  164)                ! save the surface momentum fluxes 
(  165)                ztsx(ji,jj) = zwy(ji,1) / e3t(1)
(  166)                ztsy(ji,jj) = zww(ji,1) / e3t(1)
(  167)                ! save bottom friction momentum fluxes 
(  168)                ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  169)                ikbum1 = MAX( ikbu-1, 1 )
(  170)                ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  171)                ikbvm1 = MAX( ikbv-1, 1 )
(  172)                ztbx(ji,jj) = avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
(  173)                                / ( e3t(ikbum1) * e3w(ikbu) )
(  174)                ztby(ji,jj) = avmv(ji,jj,ikbv) * zwz(ji,ikbvm1)   &
(  175)                                / ( e3t(ikbvm1) * e3w(ikbv) )
(  176)                ! subtract surface forcing and bottom friction trend from vertical
(  177)                ! diffusive momentum trend
(  178)                ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  179)                ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  180)                ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  181)                ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  182)             END DO
(  183)             !                                             ! ===============
(  184)          END DO                                           !   End of slab
(  185)          !                                                ! ===============
(  186) 
(  187)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  188)          ztdua(:,:,:) = 0.e0
(  189)          ztdva(:,:,:) = 0.e0
(  190)          ztdua(:,:,1) = ztsx(:,:)
(  191)          ztdva(:,:,1) = ztsy(:,:)
(  192)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  193)          ztdua(:,:,:) = 0.e0
(  194)          ztdva(:,:,:) = 0.e0
(  195)          ztdua(:,:,1) = ztbx(:,:)
(  196)          ztdva(:,:,1) = ztby(:,:)
(  197)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  198)       ENDIF
(  199) 
(  200)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  201)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  202)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  203)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:29      page 9

(  204) 
(  205)    END SUBROUTINE dyn_zdf_exp
(  206) 
(  207)    !!==============================================================================
(  208) END MODULE dynzdf_exp
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dynzdf_exp.F90: 89)
