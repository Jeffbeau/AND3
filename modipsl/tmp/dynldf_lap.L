


PGF90 (Version     10.5)          09/23/2010  09:49:37      page 1

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

Filename: dynldf_lap.F90

(    1) MODULE dynldf_lap
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynldf_lap  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_ldf_lap  : update the momentum trend with the lateral diffusion
(    9)    !!                  using an iso-level harmonic operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   15)    USE zdf_oce         ! ocean vertical physics
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE trdmod          ! ocean dynamics trends 
(   18)    USE trdmod_oce      ! ocean variables trends
(   19)    USE ldfslp          ! iso-neutral slopes 
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dyn_ldf_lap  ! called by step.F90
(   27) 
(   28)    !! * Substitutions
(   29) #  include "domzgr_substitute.h90"
(   30) #  include "ldfdyn_substitute.h90"
(   31) #  include "vectopt_loop_substitute.h90"
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_lap.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE dyn_ldf_lap( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                     ***  ROUTINE dyn_ldf_lap  ***
(   43)       !!                       
(   44)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   45)       !!      trend and add it to the general trend of tracer equation.
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 2

(   47)       !! ** Method  :   The before horizontal momentum diffusion trend is an
(   48)       !!      harmonic operator (laplacian type) which separates the divergent
(   49)       !!      and rotational parts of the flow.
(   50)       !!      Its horizontal components are computed as follow:
(   51)       !!         difu = 1/e1u di[ahmt hdivb] - 1/(e2u*e3u) dj-1[e3f ahmf rotb]
(   52)       !!         difv = 1/e2v dj[ahmt hdivb] + 1/(e1v*e3v) di-1[e3f ahmf rotb]
(   53)       !!      If 'key_s_coord' key is not activated, the vertical scale factor
(   54)       !!      is simplified in the rotational part of the diffusion.
(   55)       !!      Add this before trend to the general trend (ua,va):
(   56)       !!            (ua,va) = (ua,va) + (diffu,diffv)
(   57)       !!      'key_trddyn' activated: the two components of the horizontal
(   58)       !!                                 diffusion trend are saved.
(   59)       !!
(   60)       !! ** Action : - Update (ua,va) with the before iso-level harmonic 
(   61)       !!               mixing trend.
(   62)       !!             - Save in (ztdua,ztdva) arrays the trends ('key_trddyn')
(   63)       !!
(   64)       !! History :
(   65)       !!        !  90-09 (G. Madec) Original code
(   66)       !!        !  91-11 (G. Madec)
(   67)       !!        !  96-01 (G. Madec) statement function for e3 and ahm
(   68)       !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
(   69)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(   70)       !!----------------------------------------------------------------------
(   71)       !! * Modules used     
(   72)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   73)                          ztdva => sa      ! use sa as 3D workspace   
(   74) 
(   75)       !! * Arguments
(   76)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   77) 
(   78)       !! * Local declarations
(   79)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(   80)       REAL(wp) ::   &
(   81)          zua, zva, ze2u, ze1v             ! temporary scalars
(   82)       !!----------------------------------------------------------------------
(   83) 
(   84)       IF( kt == nit000 ) THEN
(   85)          IF(lwp) WRITE(numout,*)
(   86)          IF(lwp) WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacien) operator'
(   87)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   88)       ENDIF
(   89) 
(   90)       ! Save ua and va trends
(   91)       IF( l_trddyn )   THEN
(   92)          ztdua(:,:,:) = ua(:,:,:) 
(   93)          ztdva(:,:,:) = va(:,:,:) 
(   94)       ENDIF
(   95) 
(   96)       !                                                ! ===============
(   97)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   98)          !                                             ! ===============
(   99)          DO jj = 2, jpjm1
(  100)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  101) #if defined key_s_coord || defined key_partial_steps
(  102)                ze2u = rotb (ji,jj,jk)*fsahmf(ji,jj,jk)*fse3f(ji,jj,jk)
(  103)                ze1v = hdivb(ji,jj,jk)*fsahmt(ji,jj,jk)
(  104)                ! horizontal diffusive trends






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 3

(  105)                zua = - ( ze2u - rotb (ji,jj-1,jk)*fsahmf(ji,jj-1,jk)*fse3f(ji,jj-1,jk) ) / ( e2u(ji,jj) * fse3u(ji,jj,jk) ) 
(  106)                      + ( hdivb(ji+1,jj,jk)*fsahmt(ji+1,jj,jk) - ze1v                   ) / e1u(ji,jj)
(  107) 
(  108)                zva = + ( ze2u - rotb (ji-1,jj,jk)*fsahmf(ji-1,jj,jk)*fse3f(ji-1,jj,jk) ) / ( e1v(ji,jj) * fse3v(ji,jj,jk) ) 
(  109)                      + ( hdivb(ji,jj+1,jk)*fsahmt(ji,jj+1,jk) - ze1v                   ) / e2v(ji,jj)
(  110) #else
(  111)                ! horizontal diffusive trends
(  112)                ze2u = rotb (ji,jj,jk)*fsahmf(ji,jj,jk)
(  113)                ze1v = hdivb(ji,jj,jk)*fsahmt(ji,jj,jk)
(  114)                zua = - (                ze2u                  - rotb (ji,jj-1,jk)*fsahmf(ji,jj-1,jk) ) / e2u(ji,jj)   &
(  115)                      + ( hdivb(ji+1,jj,jk)*fsahmt(ji+1,jj,jk) -                ze1v                  ) / e1u(ji,jj)
(  116) 
(  117)                zva = + (                ze2u                  - rotb (ji-1,jj,jk)*fsahmf(ji-1,jj,jk) ) / e1v(ji,jj)   &
(  118)                      + ( hdivb(ji,jj+1,jk)*fsahmt(ji,jj+1,jk) -                ze1v                  ) / e2v(ji,jj)
(  119) #endif
(  120) 
(  121)                ! add it to the general momentum trends
(  122)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  123)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  124)             END DO
(  125)          END DO
(  126)          !                                             ! ===============
(  127)       END DO                                           !   End of slab
(  128)       !                                                ! ===============
(  129) 
(  130)       ! save the lateral diffusion trends for diagnostic
(  131)       ! momentum trends
(  132)       IF( l_trddyn )   THEN
(  133)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  134)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  135) 
(  136)          CALL trd_mod(ztdua, ztdva, jpdtdldf, 'DYN', kt)
(  137)       ENDIF
(  138) 
(  139)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  140)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
(  141)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  142)       ENDIF
(  143) 
(  144)    END SUBROUTINE dyn_ldf_lap
(  145) 
(  146)    !!======================================================================
(  147) END MODULE dynldf_lap





















PGF90 (Version     10.5)          09/23/2010  09:49:37      page 4

(    1) # 1 "dynldf_lap.F90"
(    1) MODULE dynldf_lap
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynldf_lap  ***
(    4)    !! Ocean dynamics:  lateral viscosity trend
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_ldf_lap  : update the momentum trend with the lateral diffusion
(    9)    !!                  using an iso-level harmonic operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldfdyn_oce      ! ocean dynamics: lateral physics
(   15)    USE zdf_oce         ! ocean vertical physics
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE trdmod          ! ocean dynamics trends 
(   18)    USE trdmod_oce      ! ocean variables trends
(   19)    USE ldfslp          ! iso-neutral slopes 
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dyn_ldf_lap  ! called by step.F90
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
(   48) # 30 "dynldf_lap.F90"
(   30) # 1 "./ldfdyn_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  ldfdyn_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
(    5)    !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_substitute.h90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 19






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 5

(   19) !   ' key_dynldf_c3d' :                  3D coefficient
(   20) # 31 "dynldf_lap.F90"
(   31) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 32 "dynldf_lap.F90"
(   32) # 32
(   32)    !!----------------------------------------------------------------------
(   33)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   34)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_lap.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   35)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE dyn_ldf_lap( kt )
(   41)       !!----------------------------------------------------------------------
(   42)       !!                     ***  ROUTINE dyn_ldf_lap  ***
(   43)       !!                       
(   44)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   45)       !!      trend and add it to the general trend of tracer equation.
(   46)       !!
(   47)       !! ** Method  :   The before horizontal momentum diffusion trend is an
(   48)       !!      harmonic operator (laplacian type) which separates the divergent
(   49)       !!      and rotational parts of the flow.
(   50)       !!      Its horizontal components are computed as follow:
(   51)       !!         difu = 1/e1u di[ahmt hdivb] - 1/(e2u*e3u) dj-1[e3f ahmf rotb]
(   52)       !!         difv = 1/e2v dj[ahmt hdivb] + 1/(e1v*e3v) di-1[e3f ahmf rotb]
(   53)       !!      If 'key_s_coord' key is not activated, the vertical scale factor
(   54)       !!      is simplified in the rotational part of the diffusion.
(   55)       !!      Add this before trend to the general trend (ua,va):
(   56)       !!            (ua,va) = (ua,va) + (diffu,diffv)
(   57)       !!      'key_trddyn' activated: the two components of the horizontal
(   58)       !!                                 diffusion trend are saved.
(   59)       !!
(   60)       !! ** Action : - Update (ua,va) with the before iso-level harmonic 
(   61)       !!               mixing trend.
(   62)       !!             - Save in (ztdua,ztdva) arrays the trends ('key_trddyn')
(   63)       !!
(   64)       !! History :
(   65)       !!        !  90-09 (G. Madec) Original code
(   66)       !!        !  91-11 (G. Madec)
(   67)       !!        !  96-01 (G. Madec) statement function for e3 and ahm
(   68)       !!   8.5  !  02-06 (G. Madec)  F90: Free form and module
(   69)       !!   9.0  !  04-08 (C. Talandier) New trends organization
(   70)       !!----------------------------------------------------------------------
(   71)       !! * Modules used     
(   72)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   73)                          ztdva => sa      ! use sa as 3D workspace   






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 6

(   74) 
(   75)       !! * Arguments
(   76)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   77) 
(   78)       !! * Local declarations
(   79)       INTEGER  ::   ji, jj, jk            ! dummy loop indices
(   80)       REAL(wp) ::   &
(   81)          zua, zva, ze2u, ze1v             ! temporary scalars
(   82)       !!----------------------------------------------------------------------
(   83) 
(   84)       IF( kt == nit000 ) THEN
(   85)          IF(lwp) WRITE(numout,*)
(   86)          IF(lwp) WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacien) operator'
(   87)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   88)       ENDIF
(   89) 
(   90)       ! Save ua and va trends
(   91)       IF( l_trddyn )   THEN
(   92)          ztdua(:,:,:) = ua(:,:,:) 
(   93)          ztdva(:,:,:) = va(:,:,:) 
(   94)       ENDIF
(   95) 
(   96)       !                                                ! ===============
(   97)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   98)          !                                             ! ===============
(   99)          DO jj = 2, jpjm1
(  100)             DO ji = 2, jpim1   ! vector opt.
(  101) # 111
(  111)                ! horizontal diffusive trends
(  112)                ze2u = rotb (ji,jj,jk)*ahm2(ji,jj,jk)
(  113)                ze1v = hdivb(ji,jj,jk)*ahm1(ji,jj,jk)
(  114)                zua = - (                ze2u                  - rotb (ji,jj-1,jk)*ahm2(ji,jj-1,jk) ) / e2u(ji,jj)   &
(  115)                      + ( hdivb(ji+1,jj,jk)*ahm1(ji+1,jj,jk) -                ze1v                  ) / e1u(ji,jj)
(  116) 
(  117)                zva = + (                ze2u                  - rotb (ji-1,jj,jk)*ahm2(ji-1,jj,jk) ) / e1v(ji,jj)   &
(  118)                      + ( hdivb(ji,jj+1,jk)*ahm1(ji,jj+1,jk) -                ze1v                  ) / e2v(ji,jj)
(  119) 
(  120) # 121
(  121)                ! add it to the general momentum trends
(  122)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  123)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  124)             END DO
(  125)          END DO
(  126)          !                                             ! ===============
(  127)       END DO                                           !   End of slab
(  128)       !                                                ! ===============
(  129) 
(  130)       ! save the lateral diffusion trends for diagnostic
(  131)       ! momentum trends
(  132)       IF( l_trddyn )   THEN
(  133)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  134)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  135) 
(  136)          CALL trd_mod(ztdua, ztdva, jpdtdldf, 'DYN', kt)
(  137)       ENDIF
(  138) 
(  139)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  140)          CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 7

(  141)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  142)       ENDIF
(  143) 
(  144)    END SUBROUTINE dyn_ldf_lap
(  145) 
(  146)    !!======================================================================
(  147) END MODULE dynldf_lap
