


PGF90 (Version     10.5)          09/23/2010  09:50:20      page 1

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

Filename: zdfevd.F90

(    1) MODULE zdfevd
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfevd  ***
(    4)    !! Ocean physics: parameterization of convection through an enhancement
(    5)    !!                of vertical eddy mixing coefficient
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   zdf_evd      : update momentum and tracer Kz at the location of
(   10)    !!                  statically unstable portion of the water column
(   11)    !!                  (called if ln_zdfevd=T)
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE zdfkpp          ! KPP vertical mixing
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zdf_evd      ! called by step.F90
(   26) 
(   27)    !! * Substitutions
(   28) #  include "domzgr_substitute.h90"
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfevd.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE zdf_evd( kt )
(   38)       !!----------------------------------------------------------------------
(   39)       !!                  ***  ROUTINE zdf_evd  ***
(   40)       !!                   
(   41)       !! ** Purpose :   Local increased the vertical eddy viscosity and diffu-
(   42)       !!      sivity coefficients when a static instability is encountered.
(   43)       !!
(   44)       !! ** Method  :   avt, and the 4 neighbouring avmu, avmv coefficients
(   45)       !!      are set to avevd (namelist parameter) if the water column is 
(   46)       !!      statically unstable (i.e. if rn2 < -1.e-12 )






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 2

(   47)       !!
(   48)       !! ** Action  :   Update avt, avmu, avmv in statically instable cases
(   49)       !!                and avt_evd which is avt due to convection
(   50)       !! References :
(   51)       !!      Lazar, A., these de l'universite Paris VI, France, 1997
(   52)       !! History :
(   53)       !!   7.0  !  97-06  (G. Madec, A. Lazar)  Original code
(   54)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   55)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   56)       !!----------------------------------------------------------------------
(   57)       !! * Arguments
(   58)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step
(   59) 
(   60)       !! * Local declarations
(   61)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   62)       !!----------------------------------------------------------------------
(   63) 
(   64)       IF( kt == nit000 ) THEN
(   65)          IF(lwp) WRITE(numout,*)
(   66)          IF(lwp) WRITE(numout,*) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
(   67)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   68)          IF(lwp) WRITE(numout,*)
(   69)       ENDIF
(   70) 
(   71)       ! Initialisation of avt_evd (vertical diffusion due to convection) to avt and avmu_evd to avmu
(   72)       avt_evd  (:,:,:) = avt(:,:,:) 
(   73)       avmu_evd (:,:,:) = avmu(:,:,:) 
(   74) 
(   75)       SELECT CASE ( nevdm )
(   76)  
(   77)       CASE ( 1 )           ! enhance vertical eddy viscosity and diffusivity (if rn2<-1.e-12)
(   78)          !                                                ! ===============
(   79)          DO jk = 1, jpkm1                                 ! Horizontal slab
(   80)             !                                             ! ===============
(   81) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(   82) !!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL!
(   83)             jj = 1                     ! big loop forced
(   84)             DO ji = jpi+2, jpij   
(   85) #   if defined key_zdfkpp
(   86) !! no implicit mixing in the boundary layer with KPP
(   87)                IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) THEN
(   88) #   else
(   89)                IF( rn2(ji,jj,jk) <= -1.e-12 ) THEN
(   90) #   endif
(   91)                   avt (ji  ,jj  ,jk) = avevd * tmask(ji  ,jj  ,jk)
(   92)                   avmu(ji  ,jj  ,jk) = avevd * umask(ji  ,jj  ,jk)
(   93)                   avmu(ji-1,jj  ,jk) = avevd * umask(ji-1,jj  ,jk)
(   94)                   avmv(ji  ,jj  ,jk) = avevd * vmask(ji  ,jj  ,jk)
(   95)                   avmv(ji  ,jj-1,jk) = avevd * vmask(ji  ,jj-1,jk)
(   96)                ENDIF
(   97)             END DO
(   98) #   else
(   99)             DO jj = 2, jpj             ! no vector opt.
(  100)                DO ji = 2, jpi
(  101) #   if defined key_zdfkpp
(  102) !! no implicit mixing in the boundary layer with KPP
(  103)                IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) THEN
(  104) #   else






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 3

(  105)                IF( rn2(ji,jj,jk) <= -1.e-12 ) THEN
(  106) #   endif
(  107)                      avt (ji  ,jj  ,jk) = avevd * tmask(ji  ,jj  ,jk)
(  108)                      avmu(ji  ,jj  ,jk) = avevd * umask(ji  ,jj  ,jk)
(  109)                      avmu(ji-1,jj  ,jk) = avevd * umask(ji-1,jj  ,jk)
(  110)                      avmv(ji  ,jj  ,jk) = avevd * vmask(ji  ,jj  ,jk)
(  111)                      avmv(ji  ,jj-1,jk) = avevd * vmask(ji  ,jj-1,jk)
(  112)                   ENDIF
(  113)                END DO
(  114)             END DO
(  115) #   endif
(  116)             !                                             ! ===============
(  117)          END DO                                           !   End of slab
(  118)          !                                                ! ===============
(  119) 
(  120)          ! Lateral boundary conditions on ( avt, avmu, avmv )   (unchanged sign)
(  121)          ! -------------------------------===================
(  122)          CALL lbc_lnk( avt , 'W', 1. )
(  123)          CALL lbc_lnk( avmu, 'U', 1. )
(  124)          CALL lbc_lnk( avmv, 'V', 1. )
(  125) 
(  126)       CASE DEFAULT         ! enhance vertical eddy diffusivity only (if rn2<-1.e-12) 
(  127)          !                                                ! ===============
(  128)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  129)             !                                             ! ===============
(  130) !!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL! 
(  131) #   if defined key_vectopt_loop   &&   ! defined key_autotasking
(  132)             jj = 1                     ! big loop forced
(  133)             DO ji = 1, jpij   
(  134) #   if defined key_zdfkpp
(  135) !! no implicit mixing in the boundary layer with KPP
(  136)                IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) &              
(  137)                   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
(  138) #   else
(  139)                IF( rn2(ji,jj,jk) <= -1.e-12 )   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
(  140) #   endif
(  141)             END DO
(  142) #   else
(  143)             DO jj = 1, jpj             ! loop over the whole domain (no lbc_lnk call)
(  144)                DO ji = 1, jpi
(  145) #   if defined key_zdfkpp
(  146) !! no implicit mixing in the boundary layer with KPP
(  147)                IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) &          
(  148)                   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
(  149) #   else
(  150)                   IF( rn2(ji,jj,jk) <= -1.e-12 )   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
(  151) #   endif
(  152)                END DO
(  153)             END DO
(  154) #   endif
(  155)             !                                             ! ===============
(  156)          END DO                                           !   End of slab
(  157)          !                                                ! ===============
(  158)       END SELECT 
(  159) 
(  160)       ! update of avt_evd and avmu_evd
(  161)       avt_evd  (:,:,:) = avt (:,:,:)  - avt_evd  (:,:,:) 
(  162)       avmu_evd (:,:,:) = avmu(:,:,:)  - avmu_evd (:,:,:) 






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 4

(  163) 
(  164)    END SUBROUTINE zdf_evd
(  165) 
(  166)    !!======================================================================
(  167) END MODULE zdfevd



























































PGF90 (Version     10.5)          09/23/2010  09:50:20      page 5

(    1) # 1 "zdfevd.F90"
(    1) MODULE zdfevd
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfevd  ***
(    4)    !! Ocean physics: parameterization of convection through an enhancement
(    5)    !!                of vertical eddy mixing coefficient
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   zdf_evd      : update momentum and tracer Kz at the location of
(   10)    !!                  statically unstable portion of the water column
(   11)    !!                  (called if ln_zdfevd=T)
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE zdfkpp          ! KPP vertical mixing
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zdf_evd      ! called by step.F90
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
(   48) # 29 "zdfevd.F90"
(   29) # 29
(   29)    !!----------------------------------------------------------------------
(   30)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   31)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfevd.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
(   32)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE zdf_evd( kt )
(   38)       !!----------------------------------------------------------------------
(   39)       !!                  ***  ROUTINE zdf_evd  ***
(   40)       !!                   
(   41)       !! ** Purpose :   Local increased the vertical eddy viscosity and diffu-






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 6

(   42)       !!      sivity coefficients when a static instability is encountered.
(   43)       !!
(   44)       !! ** Method  :   avt, and the 4 neighbouring avmu, avmv coefficients
(   45)       !!      are set to avevd (namelist parameter) if the water column is 
(   46)       !!      statically unstable (i.e. if rn2 < -1.e-12 )
(   47)       !!
(   48)       !! ** Action  :   Update avt, avmu, avmv in statically instable cases
(   49)       !!                and avt_evd which is avt due to convection
(   50)       !! References :
(   51)       !!      Lazar, A., these de l'universite Paris VI, France, 1997
(   52)       !! History :
(   53)       !!   7.0  !  97-06  (G. Madec, A. Lazar)  Original code
(   54)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   55)       !!   9.0  !  05-06  (C. Ethe) KPP parameterization
(   56)       !!----------------------------------------------------------------------
(   57)       !! * Arguments
(   58)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step
(   59) 
(   60)       !! * Local declarations
(   61)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   62)       !!----------------------------------------------------------------------
(   63) 
(   64)       IF( kt == nit000 ) THEN
(   65)          IF(lwp) WRITE(numout,*)
(   66)          IF(lwp) WRITE(numout,*) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
(   67)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   68)          IF(lwp) WRITE(numout,*)
(   69)       ENDIF
(   70) 
(   71)       ! Initialisation of avt_evd (vertical diffusion due to convection) to avt and avmu_evd to avmu
(   72)       avt_evd  (:,:,:) = avt(:,:,:) 
(   73)       avmu_evd (:,:,:) = avmu(:,:,:) 
(   74) 
(   75)       SELECT CASE ( nevdm )
(   76)  
(   77)       CASE ( 1 )           ! enhance vertical eddy viscosity and diffusivity (if rn2<-1.e-12)
(   78)          !                                                ! ===============
(   79)          DO jk = 1, jpkm1                                 ! Horizontal slab
(   80)             !                                             ! ===============
(   81) # 99
(   99)             DO jj = 2, jpj             ! no vector opt.
(  100)                DO ji = 2, jpi
(  101) # 105
(  105)                IF( rn2(ji,jj,jk) <= -1.e-12 ) THEN
(  106) # 107
(  107)                      avt (ji  ,jj  ,jk) = avevd * tmask(ji  ,jj  ,jk)
(  108)                      avmu(ji  ,jj  ,jk) = avevd * umask(ji  ,jj  ,jk)
(  109)                      avmu(ji-1,jj  ,jk) = avevd * umask(ji-1,jj  ,jk)
(  110)                      avmv(ji  ,jj  ,jk) = avevd * vmask(ji  ,jj  ,jk)
(  111)                      avmv(ji  ,jj-1,jk) = avevd * vmask(ji  ,jj-1,jk)
(  112)                   ENDIF
(  113)                END DO
(  114)             END DO
(  115) # 116
(  116)             !                                             ! ===============
(  117)          END DO                                           !   End of slab
(  118)          !                                                ! ===============
(  119) 






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 7

(  120)          ! Lateral boundary conditions on ( avt, avmu, avmv )   (unchanged sign)
(  121)          ! -------------------------------===================
(  122)          CALL lbc_lnk( avt , 'W', 1. )
(  123)          CALL lbc_lnk( avmu, 'U', 1. )
(  124)          CALL lbc_lnk( avmv, 'V', 1. )
(  125) 
(  126)       CASE DEFAULT         ! enhance vertical eddy diffusivity only (if rn2<-1.e-12) 
(  127)          !                                                ! ===============
(  128)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  129)             !                                             ! ===============
(  130) !!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL! 
(  131) # 143
(  143)             DO jj = 1, jpj             ! loop over the whole domain (no lbc_lnk call)
(  144)                DO ji = 1, jpi
(  145) # 150
(  150)                   IF( rn2(ji,jj,jk) <= -1.e-12 )   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
(  151) # 152
(  152)                END DO
(  153)             END DO
(  154) # 155
(  155)             !                                             ! ===============
(  156)          END DO                                           !   End of slab
(  157)          !                                                ! ===============
(  158)       END SELECT 
(  159) 
(  160)       ! update of avt_evd and avmu_evd
(  161)       avt_evd  (:,:,:) = avt (:,:,:)  - avt_evd  (:,:,:) 
(  162)       avmu_evd (:,:,:) = avmu(:,:,:)  - avmu_evd (:,:,:) 
(  163) 
(  164)    END SUBROUTINE zdf_evd
(  165) 
(  166)    !!======================================================================
(  167) END MODULE zdfevd
