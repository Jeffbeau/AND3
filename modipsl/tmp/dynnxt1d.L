


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

Filename: dynnxt1d.F90

(    1) MODULE dynnxt1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynnxt1d  ***
(    4)    !! Ocean dynamics: time stepping in 1D configuration
(    5)    !!======================================================================
(    6) #if defined key_cfg_1d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'               1D Configuration
(    9)    !!----------------------------------------------------------------------  
(   10)    !!----------------------------------------------------------------------
(   11)    !!   dyn_nxt_1d   : update the horizontal velocity from the momentum trend
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lbclnk          ! lateral boundary condition (or mpp link)
(   18)    USE prtctl          ! Print control
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC dyn_nxt_1d                ! routine called by step.F90
(   25)    !!----------------------------------------------------------------------
(   26) 
(   27) CONTAINS
(   28) 
(   29)    SUBROUTINE dyn_nxt_1d ( kt )
(   30)       !!----------------------------------------------------------------------
(   31)       !!                  ***  ROUTINE dyn_nxt_1d  ***
(   32)       !!                   
(   33)       !! ** Purpose :   Compute the after horizontal velocity from the 
(   34)       !!      momentum trend.
(   35)       !!
(   36)       !! ** Method  :   Apply lateral boundary conditions on the trends (ua,va) 
(   37)       !!      through calls to routine lbc_lnk.
(   38)       !!      After velocity is compute using a leap-frog scheme environment:
(   39)       !!         (ua,va) = (ub,vb) + 2 rdt (ua,va)
(   40)       !!      Time filter applied on now horizontal velocity to avoid the
(   41)       !!      divergence of two consecutive time-steps and swap of dynamics
(   42)       !!      arrays to start the next time step:
(   43)       !!         (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
(   44)       !!         (un,vn) = (ua,va) 
(   45)       !!
(   46)       !! ** Action : - Update ub,vb arrays, the before horizontal velocity






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 2

(   47)       !!             - Update un,vn arrays, the now horizontal velocity
(   48)       !!
(   49)       !! History :
(   50)       !!        !  87-02  (P. Andrich, D. L Hostis)  Original code
(   51)       !!        !  90-10  (C. Levy, G. Madec)
(   52)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   53)       !!        !  97-02  (G. Madec & M. Imbard)  opa, release 8.0
(   54)       !!        !  97-04  (A. Weaver)  Euler forward step
(   55)       !!        !  97-06  (G. Madec)  lateral boudary cond., lbc routine
(   56)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   57)       !!        !  04-10  (C. Ethe) 1D configuration
(   58)       !!----------------------------------------------------------------------
(   59)       !! * Arguments
(   60)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   61) 
(   62)       !! * Local declarations
(   63)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(   64)       REAL(wp) ::   z2dt         ! temporary scalar
(   65)       !!----------------------------------------------------------------------
(   66)       !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   67)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/dynnxt1d.F90,v 1.3 2005/10/03 09:20:35 opalod Exp $ 
(   68)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   69)       !!----------------------------------------------------------------------
(   70) 
(   71)       IF( kt == nit000 ) THEN
(   72)          IF(lwp) WRITE(numout,*)
(   73)          IF(lwp) WRITE(numout,*) 'dyn_nxt_1d : time stepping on 1D configuation'
(   74)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   75)       ENDIF
(   76) 
(   77)       ! Local constant initialization
(   78)       z2dt = 2. * rdt
(   79)       IF( neuler == 0 .AND. kt == nit000 )  z2dt = rdt
(   80) 
(   81)       ! Lateral boundary conditions on ( ua, va )
(   82)       CALL lbc_lnk( ua, 'U', -1. )
(   83)       CALL lbc_lnk( va, 'V', -1. )
(   84) 
(   85)       !                                                ! ===============
(   86)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   87)          !                                             ! ===============
(   88)          ! Next velocity
(   89)          ! -------------
(   90)          DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(   91)             DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(   92)                ! Leap-frog time stepping
(   93)                ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
(   94)                va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
(   95)             END DO
(   96)          END DO
(   97)          !                                             ! ===============
(   98)       END DO                                           !   End of slab
(   99)       !                                                ! ===============
(  100)  
(  101)      !                                                ! ===============
(  102)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  103)          !                                             ! ===============
(  104)          ! Time filter and swap of dynamics arrays






PGF90 (Version     10.5)          09/23/2010  09:49:37      page 3

(  105)          ! ------------------------------------------
(  106)          IF( neuler == 0 .AND. kt == nit000 ) THEN
(  107)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  108)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  109)                   ! Euler (forward) time stepping
(  110)                   ub(ji,jj,jk) = un(ji,jj,jk)
(  111)                   vb(ji,jj,jk) = vn(ji,jj,jk)
(  112)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  113)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  114)                END DO
(  115)             END DO
(  116)          ELSE
(  117)             DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
(  118)                DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
(  119)                   ! Leap-frog time stepping
(  120)                   ub(ji,jj,jk) = atfp * ( ub(ji,jj,jk) + ua(ji,jj,jk) ) + atfp1 * un(ji,jj,jk)
(  121)                   vb(ji,jj,jk) = atfp * ( vb(ji,jj,jk) + va(ji,jj,jk) ) + atfp1 * vn(ji,jj,jk)
(  122)                   un(ji,jj,jk) = ua(ji,jj,jk)
(  123)                   vn(ji,jj,jk) = va(ji,jj,jk)
(  124)                END DO
(  125)             END DO
(  126)          ENDIF
(  127)          !                                             ! ===============
(  128)       END DO                                           !   End of slab
(  129)       !                                                ! ===============
(  130) 
(  131)       IF(ln_ctl)   THEN
(  132)          CALL prt_ctl(tab3d_1=un, clinfo1=' nxt_1d  - Un: ', mask1=umask, &
(  133)             &         tab3d_2=vn, clinfo2=' Vn: ', mask2=vmask)
(  134)       ENDIF
(  135) 
(  136) !     IF(l_ctl)   WRITE(numout,*) ' nxt  - Un: ', SUM(un(2:nictl,2:njctl,1:jpkm1)*umask(2:nictl,2:njctl,1:jpkm1)), &
(  137) !     &                                  ' Vn: ', SUM(vn(2:nictl,2:njctl,1:jpkm1)*vmask(2:nictl,2:njctl,1:jpkm1))
(  138) 
(  139)    END SUBROUTINE dyn_nxt_1d
(  140) #else
(  141)    !!----------------------------------------------------------------------
(  142)    !!   Default key                                     NO 1D Config
(  143)    !!----------------------------------------------------------------------
(  144)    USE in_out_manager
(  145) CONTAINS
(  146)    SUBROUTINE dyn_nxt_1d ( kt )
(  147)       if(lwp) WRITE(numout,*) 'dyn_nxt_1d: You should not have seen this print! error?', kt
(  148)    END SUBROUTINE dyn_nxt_1d
(  149) #endif
(  150)    !!======================================================================
(  151) END MODULE dynnxt1d

















PGF90 (Version     10.5)          09/23/2010  09:49:37      page 4

(    1) # 1 "dynnxt1d.F90"
(    1) MODULE dynnxt1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynnxt1d  ***
(    4)    !! Ocean dynamics: time stepping in 1D configuration
(    5)    !!======================================================================
(    6) # 141
(  141)    !!----------------------------------------------------------------------
(  142)    !!   Default key                                     NO 1D Config
(  143)    !!----------------------------------------------------------------------
(  144)    USE in_out_manager
(  145) CONTAINS
(  146)    SUBROUTINE dyn_nxt_1d ( kt )
(  147)       if(lwp) WRITE(numout,*) 'dyn_nxt_1d: You should not have seen this print! error?', kt
(  148)    END SUBROUTINE dyn_nxt_1d
(  149) # 150
(  150)    !!======================================================================
(  151) END MODULE dynnxt1d
