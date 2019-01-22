


PGF90 (Version     10.5)          09/23/2010  09:49:38      page 1

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

Filename: dynnxt.F90

(    1) MODULE dynnxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynnxt  ***
(    4)    !! Ocean dynamics: time stepping
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_nxt      : update the horizontal velocity from the momentum trend
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE obc_oce         ! ocean open boundary conditions
(   15)    USE obcdyn          ! open boundary condition for momentum (obc_dyn routine)
(   16)    USE obcdyn_bt       ! 2D open boundary condition for momentum (obc_dyn_bt routine)
(   17)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   18)    USE dynspg_oce      ! type of surface pressure gradient
(   19)    USE lbclnk          ! lateral boundary condition (or mpp link)
(   20)    USE prtctl          ! Print control
(   21)    USE agrif_opa_update
(   22)    USE agrif_opa_interp
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC dyn_nxt                ! routine called by step.F90
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE dyn_nxt ( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE dyn_nxt  ***
(   36)       !!                   
(   37)       !! ** Purpose :   Compute the after horizontal velocity from the 
(   38)       !!      momentum trend.
(   39)       !!
(   40)       !! ** Method  :   Apply lateral boundary conditions on the trends (ua,va) 
(   41)       !!      through calls to routine lbc_lnk.
(   42)       !!      After velocity is compute using a leap-frog scheme environment:
(   43)       !!         (ua,va) = (ub,vb) + 2 rdt (ua,va)
(   44)       !!      Note that if lk_dynspg_flt=T, the time stepping has already been
(   45)       !!      performed in dynspg module
(   46)       !!      Time filter applied on now horizontal velocity to avoid the






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 2

(   47)       !!      divergence of two consecutive time-steps and swap of dynamics
(   48)       !!      arrays to start the next time step:
(   49)       !!         (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
(   50)       !!         (un,vn) = (ua,va) 
(   51)       !!
(   52)       !! ** Action : - Update ub,vb arrays, the before horizontal velocity
(   53)       !!             - Update un,vn arrays, the now horizontal velocity
(   54)       !!
(   55)       !! History :
(   56)       !!        !  87-02  (P. Andrich, D. L Hostis)  Original code
(   57)       !!        !  90-10  (C. Levy, G. Madec)
(   58)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   59)       !!        !  97-02  (G. Madec & M. Imbard)  opa, release 8.0
(   60)       !!        !  97-04  (A. Weaver)  Euler forward step
(   61)       !!        !  97-06  (G. Madec)  lateral boudary cond., lbc routine
(   62)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   63)       !!        !  02-10  (C. Talandier, A-M. Treguier) Open boundary cond.
(   64)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   65)       !!----------------------------------------------------------------------
(   66)       !! * Arguments
(   67)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   68) 
(   69)       !! * Local declarations
(   70)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(   71)       REAL(wp) ::   z2dt         ! temporary scalar
(   72)       !!----------------------------------------------------------------------
(   73)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   74)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynnxt.F90,v 1.10 2006/03/10 10:55:41 opalod Exp $ 
(   75)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   76)       !!----------------------------------------------------------------------
(   77) 
(   78)       IF( kt == nit000 ) THEN
(   79)          IF(lwp) WRITE(numout,*)
(   80)          IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
(   81)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   82)       ENDIF
(   83) 
(   84)       ! Local constant initialization
(   85)       z2dt = 2. * rdt
(   86)       IF( neuler == 0 .AND. kt == nit000 )  z2dt = rdt
(   87) 
(   88)       ! Lateral boundary conditions on ( ua, va )
(   89)       CALL lbc_lnk( ua, 'U', -1. )
(   90)       CALL lbc_lnk( va, 'V', -1. )
(   91) 
(   92)       !                                                ! ===============
(   93)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   94)          !                                             ! ===============
(   95)          ! Next velocity
(   96)          ! -------------
(   97) #if defined key_dynspg_flt
(   98)          ! Leap-frog time stepping already done in dynspg.F routine
(   99) #else
(  100)          DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  101)             DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  102)                ! Leap-frog time stepping
(  103)                ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
(  104)                va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 3

(  105)             END DO
(  106)          END DO
(  107) # if defined key_obc
(  108)          !                                             ! ===============
(  109)       END DO                                           !   End of slab
(  110)       !                                                ! ===============
(  111)       ! Update (ua,va) along open boundaries (only in the rigid-lid case)
(  112)       CALL obc_dyn( kt )
(  113) 
(  114)       IF ( lk_dynspg_exp .OR. lk_dynspg_ts ) THEN
(  115)          !Flather boundary condition :
(  116)          !        - Update sea surface height on each open boundary
(  117)          !                 sshn (= after ssh) for explicit case
(  118)          !                 sshn_b (= after ssha_b) for time-splitting case
(  119)          !        - Correct the barotropic velocities
(  120)          CALL obc_dyn_bt( kt )
(  121) 
(  122)          !Boundary conditions on sshn ( after ssh)
(  123)          CALL lbc_lnk( sshn, 'T', 1. )
(  124) 
(  125)          IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  126)             CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  127)          ENDIF
(  128) 
(  129)          IF ( ln_vol_cst ) CALL obc_vol( kt )
(  130) 
(  131)       ENDIF
(  132) 
(  133)       !                                                ! ===============
(  134)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  135)          !                                             ! ===============
(  136) # endif
(  137) # if defined key_agrif
(  138)          !                                             ! ===============
(  139)       END DO                                           !   End of slab
(  140)       !                                                ! ===============
(  141)       ! Update (ua,va) along open boundaries (only in the rigid-lid case)
(  142)       CALL Agrif_dyn( kt )
(  143)       !                                                ! ===============
(  144)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  145)          !                                             ! ===============
(  146) # endif
(  147) #endif
(  148)          ! Time filter and swap of dynamics arrays
(  149)          ! ------------------------------------------
(  150)          IF( neuler == 0 .AND. kt == nit000 ) THEN
(  151)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  152)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  153)                   ! Euler (forward) time stepping
(  154)                   ub(ji,jj,jk) = un(ji,jj,jk)
(  155)                   vb(ji,jj,jk) = vn(ji,jj,jk)
(  156)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  157)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  158)                END DO
(  159)             END DO
(  160)          ELSE
(  161)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  162)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 4

(  163)                   ! Leap-frog time stepping
(  164)                   ub(ji,jj,jk) = atfp * ( ub(ji,jj,jk) + ua(ji,jj,jk) ) + atfp1 * un(ji,jj,jk)
(  165)                   vb(ji,jj,jk) = atfp * ( vb(ji,jj,jk) + va(ji,jj,jk) ) + atfp1 * vn(ji,jj,jk)
(  166)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  167)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  168)                END DO
(  169)             END DO
(  170)          ENDIF
(  171)          !                                             ! ===============
(  172)       END DO                                           !   End of slab
(  173)       !                                                ! ===============
(  174) 
(  175)       IF(ln_ctl)   THEN
(  176)          CALL prt_ctl(tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask, &
(  177)             &         tab3d_2=vn, clinfo2=' Vn: ', mask2=vmask)
(  178)       ENDIF
(  179) 
(  180) #if defined key_agrif
(  181)       IF (.NOT.Agrif_Root())    CALL Agrif_Update_Dyn( kt )
(  182) #endif      
(  183) 
(  184)    END SUBROUTINE dyn_nxt
(  185) 
(  186)    !!======================================================================
(  187) END MODULE dynnxt







































PGF90 (Version     10.5)          09/23/2010  09:49:38      page 5

(    1) # 1 "dynnxt.F90"
(    1) MODULE dynnxt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynnxt  ***
(    4)    !! Ocean dynamics: time stepping
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_nxt      : update the horizontal velocity from the momentum trend
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE obc_oce         ! ocean open boundary conditions
(   15)    USE obcdyn          ! open boundary condition for momentum (obc_dyn routine)
(   16)    USE obcdyn_bt       ! 2D open boundary condition for momentum (obc_dyn_bt routine)
(   17)    USE obcvol          ! ocean open boundary condition (obc_vol routines)
(   18)    USE dynspg_oce      ! type of surface pressure gradient
(   19)    USE lbclnk          ! lateral boundary condition (or mpp link)
(   20)    USE prtctl          ! Print control
(   21)    USE agrif_opa_update
(   22)    USE agrif_opa_interp
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Accessibility
(   28)    PUBLIC dyn_nxt                ! routine called by step.F90
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE dyn_nxt ( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE dyn_nxt  ***
(   36)       !!                   
(   37)       !! ** Purpose :   Compute the after horizontal velocity from the 
(   38)       !!      momentum trend.
(   39)       !!
(   40)       !! ** Method  :   Apply lateral boundary conditions on the trends (ua,va) 
(   41)       !!      through calls to routine lbc_lnk.
(   42)       !!      After velocity is compute using a leap-frog scheme environment:
(   43)       !!         (ua,va) = (ub,vb) + 2 rdt (ua,va)
(   44)       !!      Note that if lk_dynspg_flt=T, the time stepping has already been
(   45)       !!      performed in dynspg module
(   46)       !!      Time filter applied on now horizontal velocity to avoid the
(   47)       !!      divergence of two consecutive time-steps and swap of dynamics
(   48)       !!      arrays to start the next time step:
(   49)       !!         (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
(   50)       !!         (un,vn) = (ua,va) 
(   51)       !!
(   52)       !! ** Action : - Update ub,vb arrays, the before horizontal velocity
(   53)       !!             - Update un,vn arrays, the now horizontal velocity
(   54)       !!
(   55)       !! History :
(   56)       !!        !  87-02  (P. Andrich, D. L Hostis)  Original code
(   57)       !!        !  90-10  (C. Levy, G. Madec)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 6

(   58)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   59)       !!        !  97-02  (G. Madec & M. Imbard)  opa, release 8.0
(   60)       !!        !  97-04  (A. Weaver)  Euler forward step
(   61)       !!        !  97-06  (G. Madec)  lateral boudary cond., lbc routine
(   62)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   63)       !!        !  02-10  (C. Talandier, A-M. Treguier) Open boundary cond.
(   64)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   65)       !!----------------------------------------------------------------------
(   66)       !! * Arguments
(   67)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   68) 
(   69)       !! * Local declarations
(   70)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(   71)       REAL(wp) ::   z2dt         ! temporary scalar
(   72)       !!----------------------------------------------------------------------
(   73)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   74)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynnxt.F90,v 1.10 2006/03/10 10:55:41 opalod Exp $ 
(   75)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   76)       !!----------------------------------------------------------------------
(   77) 
(   78)       IF( kt == nit000 ) THEN
(   79)          IF(lwp) WRITE(numout,*)
(   80)          IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
(   81)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   82)       ENDIF
(   83) 
(   84)       ! Local constant initialization
(   85)       z2dt = 2. * rdt
(   86)       IF( neuler == 0 .AND. kt == nit000 )  z2dt = rdt
(   87) 
(   88)       ! Lateral boundary conditions on ( ua, va )
(   89)       CALL lbc_lnk( ua, 'U', -1. )
(   90)       CALL lbc_lnk( va, 'V', -1. )
(   91) 
(   92)       !                                                ! ===============
(   93)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   94)          !                                             ! ===============
(   95)          ! Next velocity
(   96)          ! -------------
(   97) # 100
(  100)          DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  101)             DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  102)                ! Leap-frog time stepping
(  103)                ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
(  104)                va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
(  105)             END DO
(  106)          END DO
(  107) # 108
(  108)          !                                             ! ===============
(  109)       END DO                                           !   End of slab
(  110)       !                                                ! ===============
(  111)       ! Update (ua,va) along open boundaries (only in the rigid-lid case)
(  112)       CALL obc_dyn( kt )
(  113) 
(  114)       IF ( lk_dynspg_exp .OR. lk_dynspg_ts ) THEN
(  115)          !Flather boundary condition :
(  116)          !        - Update sea surface height on each open boundary
(  117)          !                 sshn (= after ssh) for explicit case






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 7

(  118)          !                 sshn_b (= after ssha_b) for time-splitting case
(  119)          !        - Correct the barotropic velocities
(  120)          CALL obc_dyn_bt( kt )
(  121) 
(  122)          !Boundary conditions on sshn ( after ssh)
(  123)          CALL lbc_lnk( sshn, 'T', 1. )
(  124) 
(  125)          IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  126)             CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
(  127)          ENDIF
(  128) 
(  129)          IF ( ln_vol_cst ) CALL obc_vol( kt )
(  130) 
(  131)       ENDIF
(  132) 
(  133)       !                                                ! ===============
(  134)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  135)          !                                             ! ===============
(  136) # 148
(  148)          ! Time filter and swap of dynamics arrays
(  149)          ! ------------------------------------------
(  150)          IF( neuler == 0 .AND. kt == nit000 ) THEN
(  151)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  152)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  153)                   ! Euler (forward) time stepping
(  154)                   ub(ji,jj,jk) = un(ji,jj,jk)
(  155)                   vb(ji,jj,jk) = vn(ji,jj,jk)
(  156)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  157)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  158)                END DO
(  159)             END DO
(  160)          ELSE
(  161)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  162)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  163)                   ! Leap-frog time stepping
(  164)                   ub(ji,jj,jk) = atfp * ( ub(ji,jj,jk) + ua(ji,jj,jk) ) + atfp1 * un(ji,jj,jk)
(  165)                   vb(ji,jj,jk) = atfp * ( vb(ji,jj,jk) + va(ji,jj,jk) ) + atfp1 * vn(ji,jj,jk)
(  166)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  167)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  168)                END DO
(  169)             END DO
(  170)          ENDIF
(  171)          !                                             ! ===============
(  172)       END DO                                           !   End of slab
(  173)       !                                                ! ===============
(  174) 
(  175)       IF(ln_ctl)   THEN
(  176)          CALL prt_ctl(tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask, &
(  177)             &         tab3d_2=vn, clinfo2=' Vn: ', mask2=vmask)
(  178)       ENDIF
(  179) 
(  180) 
(  181) # 184
(  184)    END SUBROUTINE dyn_nxt
(  185) 
(  186)    !!======================================================================
(  187) END MODULE dynnxt
