


PGF90 (Version     10.5)          09/23/2010  09:49:45      page 1

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

Filename: dynzad.F90

(    1) MODULE dynzad
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynzad  ***
(    4)    !! Ocean dynamics : vertical advection trend
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zad      : vertical advection momentum trend
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE trdmod          ! ocean dynamics trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE flxrnf          ! ocean runoffs
(   17)    USE prtctl          ! Print control
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21)    
(   22)    !! * Accessibility
(   23)    PUBLIC dyn_zad                ! routine called by step.F90
(   24) 
(   25)    !! * Substitutions
(   26) #  include "domzgr_substitute.h90"
(   27) #  include "vectopt_loop_substitute.h90"
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzad.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36) #if defined key_autotasking
(   37)    !!----------------------------------------------------------------------
(   38)    !!   'key_autotasking'                              j-k-i loops (j-slab)
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41)    SUBROUTINE dyn_zad( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE dynzad  ***
(   44)       !!
(   45)       !! ** Purpose :   Compute the now vertical momentum advection trend and 
(   46)       !!      add it to the general trend of momentum equation.






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 2

(   47)       !!
(   48)       !! ** Method  :   Use j-slab (j-k-i loops) for auto-tasking
(   49)       !!      The now vertical advection of momentum is given by:
(   50)       !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
(   51)       !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
(   52)       !!      Add this trend to the general trend (ua,va):
(   53)       !!         (ua,va) = (ua,va) + w dz(u,v)
(   54)       !!
(   55)       !! ** Action  : - Update (ua,va) with the vert. momentum advection trends
(   56)       !!              - Save the trends in (utrd,vtrd) ('key_trddyn')
(   57)       !!
(   58)       !! History :
(   59)       !!   6.0  !  91-01  (G. Madec) Original code
(   60)       !!   7.0  !  91-11  (G. Madec)
(   61)       !!   7.5  !  96-01  (G. Madec) statement function for e3
(   62)       !!   8.5  !  02-07  (G. Madec) Free form, F90
(   63)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   64)       !!----------------------------------------------------------------------
(   65)       !! * modules used
(   66)       USE oce, ONLY:   zwuw => ta,   & ! use ta as 3D workspace
(   67)                        zwvw => sa      ! use sa as 3D workspace
(   68) 
(   69)       !! * Arguments
(   70)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step inedx
(   71)       
(   72)       !! * Local declarations
(   73)       INTEGER  ::   ji, jj, jk         ! dummy loop indices
(   74)       REAL(wp) ::   zvn, zua, zva      ! temporary scalars
(   75)       REAL(wp), DIMENSION(jpi) ::   &
(   76)          zww                           ! temporary workspace
(   77)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   78)          ztdua, ztdva                  ! temporary workspace
(   79)       !!----------------------------------------------------------------------
(   80)       
(   81)       IF( kt == nit000 ) THEN
(   82)          IF(lwp) WRITE(numout,*)
(   83)          IF(lwp) WRITE(numout,*) 'dyn_zad : arakawa advection scheme'
(   84)          IF(lwp) WRITE(numout,*) '~~~~~~~   Auto-tasking case, j-slab, no vector opt.'
(   85)       ENDIF
(   86) 
(   87)       ! Save ua and va trends
(   88)       IF( l_trddyn )   THEN
(   89)          ztdua(:,:,:) = ua(:,:,:) 
(   90)          ztdva(:,:,:) = va(:,:,:) 
(   91)       ENDIF
(   92) 
(   93)       !                                                ! ===============
(   94)       DO jj = 2, jpjm1                                 !  Vertical slab
(   95)          !                                             ! ===============
(   96) 
(   97)          ! Vertical momentum advection at level w and u- and v- vertical
(   98)          ! ----------------------------------------------------------------
(   99)          DO jk = 2, jpkm1
(  100)             ! vertical fluxes 
(  101)             DO ji = 2, jpi
(  102)                zww(ji) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
(  103)             END DO
(  104)             ! vertical momentum advection at w-point






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 3

(  105)             DO ji = 2, jpim1
(  106)                zvn = 0.25 * e1t(ji,jj+1) * e2t(ji,jj+1) * wn(ji,jj+1,jk)
(  107)                zwuw(ji,jj,jk) = ( zww(ji+1) + zww(ji) ) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
(  108)                zwvw(ji,jj,jk) = ( zvn       + zww(ji) ) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
(  109)             END DO  
(  110)          END DO   
(  111) 
(  112)          ! Surface and bottom values set to zero
(  113)          DO ji = 2, jpim1
(  114)             zwuw(ji,jj, 1 ) = 0.e0
(  115)             zwvw(ji,jj, 1 ) = 0.e0
(  116)             zwuw(ji,jj,jpk) = 0.e0
(  117)             zwvw(ji,jj,jpk) = 0.e0
(  118)          END DO  
(  119) 
(  120)          ! Vertical momentum advection at u- and v-points
(  121)          ! ----------------------------------------------
(  122)          DO jk = 1, jpkm1
(  123)             DO ji = 2, jpim1
(  124)                ! vertical momentum advective trends
(  125)                zua = - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) )
(  126)                zva = - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) )
(  127)                ! add the trends to the general momentum trends
(  128)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  129)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  130)             END DO  
(  131)          END DO  
(  132)          !                                             ! ===============
(  133)       END DO                                           !   End of slab
(  134)       !                                                ! ===============
(  135) 
(  136)       ! save the vertical advection trends for diagnostic
(  137)       ! momentum trends
(  138)       IF( l_trddyn )   THEN
(  139)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  140)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  141) 
(  142)          CALL trd_mod(ztdua, ztdva, jpdtdzad, 'DYN', kt)
(  143)       ENDIF
(  144) 
(  145)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  146)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask, &
(  147)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  148)       ENDIF
(  149) 
(  150)    END SUBROUTINE dyn_zad
(  151) 
(  152) #else
(  153)    !!----------------------------------------------------------------------
(  154)    !!   Default option                             k-j-i loop (vector opt.)
(  155)    !!----------------------------------------------------------------------
(  156) 
(  157)    SUBROUTINE dyn_zad ( kt )
(  158)       !!----------------------------------------------------------------------
(  159)       !!                  ***  ROUTINE dynzad  ***
(  160)       !! 
(  161)       !! ** Purpose :   Compute the now vertical momentum advection trend and 
(  162)       !!      add it to the general trend of momentum equation.






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 4

(  163)       !!
(  164)       !! ** Method  :   The now vertical advection of momentum is given by:
(  165)       !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
(  166)       !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
(  167)       !!      Add this trend to the general trend (ua,va):
(  168)       !!         (ua,va) = (ua,va) + w dz(u,v)
(  169)       !!
(  170)       !! ** Action  : - Update (ua,va) with the vert. momentum adv. trends
(  171)       !!              - Save the trends in (utrd,vtrd) ('key_trddyn')
(  172)       !!
(  173)       !! History :
(  174)       !!   8.5  !  02-07  (G. Madec)  Original code
(  175)       !!----------------------------------------------------------------------
(  176)       !! * modules used
(  177)       USE oce, ONLY:   zwuw => ta,   & ! use ta as 3D workspace
(  178)                        zwvw => sa      ! use sa as 3D workspace
(  179)       !! * Arguments
(  180)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step inedx
(  181)       
(  182)       !! * Local declarations
(  183)       INTEGER  ::   ji, jj, jk         ! dummy loop indices
(  184)       REAL(wp) ::   zua, zva           ! temporary scalars
(  185)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  186)          zww                           ! temporary  workspace
(  187)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  188)          ztdua, ztdva                  ! temporary workspace
(  189)       !!----------------------------------------------------------------------
(  190)       
(  191)       IF( kt == nit000 ) THEN
(  192)          IF(lwp)WRITE(numout,*)
(  193)          IF(lwp)WRITE(numout,*) 'dyn_zad : arakawa advection scheme'
(  194)          IF(lwp)WRITE(numout,*) '~~~~~~~   vector optimization k-j-i loop'
(  195)       ENDIF
(  196) 
(  197)       ! Save ua and va trends
(  198)       IF( l_trddyn )   THEN
(  199)          ztdua(:,:,:) = ua(:,:,:) 
(  200)          ztdva(:,:,:) = va(:,:,:) 
(  201)       ENDIF
(  202)       
(  203)       ! Vertical momentum advection at level w and u- and v- vertical
(  204)       ! -------------------------------------------------------------
(  205)       DO jk = 2, jpkm1
(  206)          ! vertical fluxes 
(  207)          DO jj = 2, jpj
(  208)             DO ji = fs_2, jpi   ! vector opt.
(  209)                zww(ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
(  210)             END DO
(  211)          END DO
(  212)          ! vertical momentum advection at w-point
(  213)          DO jj = 2, jpjm1
(  214)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  215)                zwuw(ji,jj,jk) = ( zww(ji+1,jj  ) + zww(ji,jj) ) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
(  216)                zwvw(ji,jj,jk) = ( zww(ji  ,jj+1) + zww(ji,jj) ) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
(  217)             END DO  
(  218)          END DO   
(  219)       END DO
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 5

(  221)       ! Surface and bottom values set to zero
(  222)       DO jj = 2, jpjm1
(  223)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  224)             zwuw(ji,jj, 1 ) = 0.e0
(  225)             zwvw(ji,jj, 1 ) = 0.e0
(  226)             zwuw(ji,jj,jpk) = 0.e0
(  227)             zwvw(ji,jj,jpk) = 0.e0
(  228)          END DO  
(  229)       END DO
(  230) 
(  231) 
(  232)       ! Vertical momentum advection at u- and v-points
(  233)       ! ----------------------------------------------
(  234)       DO jk = 1, jpkm1
(  235)          DO jj = 2, jpjm1
(  236)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  237)                ! vertical momentum advective trends
(  238)                zua = - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk) )
(  239)                zva = - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk) )
(  240)                ! add the trends to the general momentum trends
(  241)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  242)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  243)             END DO  
(  244)          END DO  
(  245)       END DO
(  246) 
(  247)       ! save the vertical advection trends for diagnostic
(  248)       ! momentum trends
(  249)       IF( l_trddyn )   THEN
(  250)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  251)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  252) 
(  253)          CALL trd_mod(ztdua, ztdva, jpdtdzad, 'DYN', kt)
(  254)       ENDIF
(  255) 
(  256)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  257)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask, &
(  258)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  259)       ENDIF
(  260) 
(  261)    END SUBROUTINE dyn_zad
(  262) #endif
(  263) 
(  264) !!======================================================================
(  265) END MODULE dynzad



















PGF90 (Version     10.5)          09/23/2010  09:49:45      page 6

(    1) # 1 "dynzad.F90"
(    1) MODULE dynzad
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynzad  ***
(    4)    !! Ocean dynamics : vertical advection trend
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zad      : vertical advection momentum trend
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE trdmod          ! ocean dynamics trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE flxrnf          ! ocean runoffs
(   17)    USE prtctl          ! Print control
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21)    
(   22)    !! * Accessibility
(   23)    PUBLIC dyn_zad                ! routine called by step.F90
(   24) 
(   25)    !! * Substitutions
(   26) # 1 "./domzgr_substitute.h90"
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
(   48) # 27 "dynzad.F90"
(   27) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 28 "dynzad.F90"
(   28) # 28
(   28)    !!----------------------------------------------------------------------
(   29)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 7

(   30)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzad.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   31)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   32)    !!----------------------------------------------------------------------
(   33) 
(   34) CONTAINS
(   35) 
(   36) # 153
(  153)    !!----------------------------------------------------------------------
(  154)    !!   Default option                             k-j-i loop (vector opt.)
(  155)    !!----------------------------------------------------------------------
(  156) 
(  157)    SUBROUTINE dyn_zad ( kt )
(  158)       !!----------------------------------------------------------------------
(  159)       !!                  ***  ROUTINE dynzad  ***
(  160)       !! 
(  161)       !! ** Purpose :   Compute the now vertical momentum advection trend and 
(  162)       !!      add it to the general trend of momentum equation.
(  163)       !!
(  164)       !! ** Method  :   The now vertical advection of momentum is given by:
(  165)       !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
(  166)       !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
(  167)       !!      Add this trend to the general trend (ua,va):
(  168)       !!         (ua,va) = (ua,va) + w dz(u,v)
(  169)       !!
(  170)       !! ** Action  : - Update (ua,va) with the vert. momentum adv. trends
(  171)       !!              - Save the trends in (utrd,vtrd) ('key_trddyn')
(  172)       !!
(  173)       !! History :
(  174)       !!   8.5  !  02-07  (G. Madec)  Original code
(  175)       !!----------------------------------------------------------------------
(  176)       !! * modules used
(  177)       USE oce, ONLY:   zwuw => ta,   & ! use ta as 3D workspace
(  178)                        zwvw => sa      ! use sa as 3D workspace
(  179)       !! * Arguments
(  180)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step inedx
(  181)       
(  182)       !! * Local declarations
(  183)       INTEGER  ::   ji, jj, jk         ! dummy loop indices
(  184)       REAL(wp) ::   zua, zva           ! temporary scalars
(  185)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  186)          zww                           ! temporary  workspace
(  187)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  188)          ztdua, ztdva                  ! temporary workspace
(  189)       !!----------------------------------------------------------------------
(  190)       
(  191)       IF( kt == nit000 ) THEN
(  192)          IF(lwp)WRITE(numout,*)
(  193)          IF(lwp)WRITE(numout,*) 'dyn_zad : arakawa advection scheme'
(  194)          IF(lwp)WRITE(numout,*) '~~~~~~~   vector optimization k-j-i loop'
(  195)       ENDIF
(  196) 
(  197)       ! Save ua and va trends
(  198)       IF( l_trddyn )   THEN
(  199)          ztdua(:,:,:) = ua(:,:,:) 
(  200)          ztdva(:,:,:) = va(:,:,:) 
(  201)       ENDIF
(  202)       
(  203)       ! Vertical momentum advection at level w and u- and v- vertical






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 8

(  204)       ! -------------------------------------------------------------
(  205)       DO jk = 2, jpkm1
(  206)          ! vertical fluxes 
(  207)          DO jj = 2, jpj
(  208)             DO ji = 2, jpi   ! vector opt.
(  209)                zww(ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
(  210)             END DO
(  211)          END DO
(  212)          ! vertical momentum advection at w-point
(  213)          DO jj = 2, jpjm1
(  214)             DO ji = 2, jpim1   ! vector opt.
(  215)                zwuw(ji,jj,jk) = ( zww(ji+1,jj  ) + zww(ji,jj) ) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
(  216)                zwvw(ji,jj,jk) = ( zww(ji  ,jj+1) + zww(ji,jj) ) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
(  217)             END DO  
(  218)          END DO   
(  219)       END DO
(  220) 
(  221)       ! Surface and bottom values set to zero
(  222)       DO jj = 2, jpjm1
(  223)          DO ji = 2, jpim1   ! vector opt.
(  224)             zwuw(ji,jj, 1 ) = 0.e0
(  225)             zwvw(ji,jj, 1 ) = 0.e0
(  226)             zwuw(ji,jj,jpk) = 0.e0
(  227)             zwvw(ji,jj,jpk) = 0.e0
(  228)          END DO  
(  229)       END DO
(  230) 
(  231) 
(  232)       ! Vertical momentum advection at u- and v-points
(  233)       ! ----------------------------------------------
(  234)       DO jk = 1, jpkm1
(  235)          DO jj = 2, jpjm1
(  236)             DO ji = 2, jpim1   ! vector opt.
(  237)                ! vertical momentum advective trends
(  238)                zua = - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * e3t(jk) )
(  239)                zva = - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * e3t(jk) )
(  240)                ! add the trends to the general momentum trends
(  241)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  242)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  243)             END DO  
(  244)          END DO  
(  245)       END DO
(  246) 
(  247)       ! save the vertical advection trends for diagnostic
(  248)       ! momentum trends
(  249)       IF( l_trddyn )   THEN
(  250)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  251)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  252) 
(  253)          CALL trd_mod(ztdua, ztdva, jpdtdzad, 'DYN', kt)
(  254)       ENDIF
(  255) 
(  256)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  257)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask, &
(  258)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  259)       ENDIF
(  260) 
(  261)    END SUBROUTINE dyn_zad






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 9

(  262) 
(  263) # 264
(  264) !!======================================================================
(  265) END MODULE dynzad
