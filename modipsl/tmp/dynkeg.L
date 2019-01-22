


PGF90 (Version     10.5)          09/23/2010  09:49:36      page 1

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

Filename: dynkeg.F90

(    1) MODULE dynkeg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynkeg  ***
(    4)    !! Ocean dynamics:  kinetic energy gradient trend
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_keg      : update the momentum trend with the horizontal tke
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE trdmod          ! ocean dynamics trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE prtctl          ! Print control
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC dyn_keg                ! routine called by step.F90
(   23)    
(   24)    !! * Substitutions
(   25) #  include "vectopt_loop_substitute.h90"
(   26)    !!---------------------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynkeg.F90,v 1.7 2005/09/02 15:45:23 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!---------------------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE dyn_keg( kt )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                  ***  ROUTINE dyn_keg  ***
(   37)       !!
(   38)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(   39)       !!      gradient of the horizontal kinetic energy and add it to the 
(   40)       !!      general momentum trend.
(   41)       !!
(   42)       !! ** Method  :   Compute the now horizontal kinetic energy:
(   43)       !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
(   44)       !!      Take its horizontal gradient and add it to the general momentum
(   45)       !!      trend (ua,va).
(   46)       !!         ua = ua - 1/e1u di[ zhke ]






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 2

(   47)       !!         va = va - 1/e2v dj[ zhke ]
(   48)       !!
(   49)       !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
(   50)       !!             - Save the trends in (utrd,vtrd) ('key_trddyn')
(   51)       !!
(   52)       !! History :
(   53)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(   54)       !!   7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
(   55)       !!   9.0  !  02-07  (G. Madec)  F90: Free form and module
(   56)       !!    "   !  04-08  (C. Talandier) New trends organization
(   57)       !!----------------------------------------------------------------------
(   58)       !! * Modules used     
(   59)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   60)                          ztdva => sa      ! use sa as 3D workspace   
(   61) 
(   62)       !! * Arguments
(   63)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(   64) 
(   65)       !! * Local declarations
(   66)       INTEGER  ::   ji, jj, jk          ! dummy loop indices
(   67)       REAL(wp) ::   zua, zva, zu, zv    ! temporary scalars
(   68)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   69)          zhke                           ! temporary workspace
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72)       IF( kt == nit000 ) THEN
(   73)          IF(lwp) WRITE(numout,*)
(   74)          IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend'
(   75)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   76)       ENDIF
(   77) 
(   78)       ! Save ua and va trends
(   79)       IF( l_trddyn )   THEN
(   80)          ztdua(:,:,:) = ua(:,:,:) 
(   81)          ztdva(:,:,:) = va(:,:,:) 
(   82)       ENDIF
(   83)       
(   84)       !                                                ! ===============
(   85)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   86)          !                                             ! ===============
(   87)          ! Horizontal kinetic energy at T-point
(   88)          DO jj = 2, jpj
(   89)             DO ji = fs_2, jpi   ! vector opt.
(   90)                zv = 0.25 * (  vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
(   91)                             + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)  )
(   92)                zu = 0.25 * (  un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
(   93)                             + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)  )
(   94)                zhke(ji,jj,jk) = zv + zu
(   95)             END DO  
(   96)          END DO  
(   97)          
(   98)          ! Horizontal gradient of Horizontal kinetic energy
(   99)          DO jj = 2, jpjm1
(  100)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  101)                ! gradient of kinetic energy
(  102)                zua = -( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
(  103)                zva = -( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
(  104)                ! add to the general momentum trends






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 3

(  105)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  106)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  107)             END DO 
(  108)          END DO
(  109)          !                                             ! ===============
(  110)       END DO                                           !   End of slab
(  111)       !                                                ! ===============
(  112) 
(  113)       ! save the Kinetic Energy trends for diagnostic
(  114)       ! momentum trends
(  115)       IF( l_trddyn )   THEN
(  116)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  117)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  118) 
(  119)          CALL trd_mod(ztdua, ztdva, jpdtdkeg, 'DYN', kt)
(  120)       ENDIF
(  121) 
(  122)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  123)          CALL prt_ctl(tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask, &
(  124)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  125)       ENDIF
(  126) 
(  127)    END SUBROUTINE dyn_keg
(  128) 
(  129)    !!======================================================================
(  130) END MODULE dynkeg






































PGF90 (Version     10.5)          09/23/2010  09:49:36      page 4

(    1) # 1 "dynkeg.F90"
(    1) MODULE dynkeg
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  dynkeg  ***
(    4)    !! Ocean dynamics:  kinetic energy gradient trend
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_keg      : update the momentum trend with the horizontal tke
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers
(   12)    USE dom_oce         ! ocean space and time domain
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE trdmod          ! ocean dynamics trends 
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE prtctl          ! Print control
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC dyn_keg                ! routine called by step.F90
(   23)    
(   24)    !! * Substitutions
(   25) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 26 "dynkeg.F90"
(   26) # 26
(   26)    !!---------------------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynkeg.F90,v 1.7 2005/09/02 15:45:23 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!---------------------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE dyn_keg( kt )
(   35)       !!----------------------------------------------------------------------
(   36)       !!                  ***  ROUTINE dyn_keg  ***
(   37)       !!
(   38)       !! ** Purpose :   Compute the now momentum trend due to the horizontal
(   39)       !!      gradient of the horizontal kinetic energy and add it to the 
(   40)       !!      general momentum trend.
(   41)       !!
(   42)       !! ** Method  :   Compute the now horizontal kinetic energy:
(   43)       !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
(   44)       !!      Take its horizontal gradient and add it to the general momentum






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 5

(   45)       !!      trend (ua,va).
(   46)       !!         ua = ua - 1/e1u di[ zhke ]
(   47)       !!         va = va - 1/e2v dj[ zhke ]
(   48)       !!
(   49)       !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
(   50)       !!             - Save the trends in (utrd,vtrd) ('key_trddyn')
(   51)       !!
(   52)       !! History :
(   53)       !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
(   54)       !!   7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
(   55)       !!   9.0  !  02-07  (G. Madec)  F90: Free form and module
(   56)       !!    "   !  04-08  (C. Talandier) New trends organization
(   57)       !!----------------------------------------------------------------------
(   58)       !! * Modules used     
(   59)       USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
(   60)                          ztdva => sa      ! use sa as 3D workspace   
(   61) 
(   62)       !! * Arguments
(   63)       INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
(   64) 
(   65)       !! * Local declarations
(   66)       INTEGER  ::   ji, jj, jk          ! dummy loop indices
(   67)       REAL(wp) ::   zua, zva, zu, zv    ! temporary scalars
(   68)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   69)          zhke                           ! temporary workspace
(   70)       !!----------------------------------------------------------------------
(   71) 
(   72)       IF( kt == nit000 ) THEN
(   73)          IF(lwp) WRITE(numout,*)
(   74)          IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend'
(   75)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   76)       ENDIF
(   77) 
(   78)       ! Save ua and va trends
(   79)       IF( l_trddyn )   THEN
(   80)          ztdua(:,:,:) = ua(:,:,:) 
(   81)          ztdva(:,:,:) = va(:,:,:) 
(   82)       ENDIF
(   83)       
(   84)       !                                                ! ===============
(   85)       DO jk = 1, jpkm1                                 ! Horizontal slab
(   86)          !                                             ! ===============
(   87)          ! Horizontal kinetic energy at T-point
(   88)          DO jj = 2, jpj
(   89)             DO ji = 2, jpi   ! vector opt.
(   90)                zv = 0.25 * (  vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
(   91)                             + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)  )
(   92)                zu = 0.25 * (  un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
(   93)                             + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)  )
(   94)                zhke(ji,jj,jk) = zv + zu
(   95)             END DO  
(   96)          END DO  
(   97)          
(   98)          ! Horizontal gradient of Horizontal kinetic energy
(   99)          DO jj = 2, jpjm1
(  100)             DO ji = 2, jpim1   ! vector opt.
(  101)                ! gradient of kinetic energy
(  102)                zua = -( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:36      page 6

(  103)                zva = -( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
(  104)                ! add to the general momentum trends
(  105)                ua(ji,jj,jk) = ua(ji,jj,jk) + zua
(  106)                va(ji,jj,jk) = va(ji,jj,jk) + zva
(  107)             END DO 
(  108)          END DO
(  109)          !                                             ! ===============
(  110)       END DO                                           !   End of slab
(  111)       !                                                ! ===============
(  112) 
(  113)       ! save the Kinetic Energy trends for diagnostic
(  114)       ! momentum trends
(  115)       IF( l_trddyn )   THEN
(  116)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  117)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  118) 
(  119)          CALL trd_mod(ztdua, ztdva, jpdtdkeg, 'DYN', kt)
(  120)       ENDIF
(  121) 
(  122)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  123)          CALL prt_ctl(tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask, &
(  124)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  125)       ENDIF
(  126) 
(  127)    END SUBROUTINE dyn_keg
(  128) 
(  129)    !!======================================================================
(  130) END MODULE dynkeg
