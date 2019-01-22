


PGF90 (Version     10.5)          09/23/2010  09:50:03      page 1

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

Filename: wzvmod.F90

(    1) MODULE wzvmod
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  wzvmod  ***
(    4)    !! Ocean diagnostic variable : vertical velocity
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   wzv        : Compute the vertical velocity
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE prtctl          ! Print control
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC wzv       ! routine called by step.F90 and inidtr.F90
(   21) 
(   22)    !! * Substitutions
(   23) #  include "domzgr_substitute.h90"
(   24)    !!----------------------------------------------------------------------
(   25) 
(   26) CONTAINS
(   27) 
(   28) #if defined key_autotasking
(   29)    !!----------------------------------------------------------------------
(   30)    !!   'key_autotasking'                               j-k-i loop (j-slab)
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33)    SUBROUTINE wzv( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                    ***  ROUTINE wzv  ***
(   36)       !!                     
(   37)       !! ** Purpose :   Compute the now vertical velocity after the array swap
(   38)       !!
(   39)       !! ** Method  :   Using the incompressibility hypothesis, the vertical
(   40)       !!     velocity is computed by integrating the horizontal divergence 
(   41)       !!     from the bottom to the surface.
(   42)       !!       The boundary conditions are w=0 at the bottom (no flux) and,
(   43)       !!     in regid-lid case, w=0 at the sea surface.
(   44)       !!
(   45)       !! ** action  :    wn array : the now vertical velocity
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 2

(   47)       !! History :
(   48)       !!   5.0  !  90-10  (C. Levy, G. Madec)  Original code
(   49)       !!   7.0  !  96-01  (G. Madec)  Statement function for e3
(   50)       !!   8.5  !  02-07  (G. Madec)  Free form, F90
(   51)       !!----------------------------------------------------------------------
(   52)       !! * Arguments
(   53)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   54) 
(   55)       !! * Local declarations
(   56)       INTEGER ::   jj, jk      ! dummy loop indices
(   57)       !!----------------------------------------------------------------------
(   58)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   59)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/wzvmod.F90,v 1.4 2005/09/02 15:45:24 opalod Exp $ 
(   60)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   61)       !!----------------------------------------------------------------------
(   62) 
(   63)       IF( kt == nit000 ) THEN
(   64)          IF(lwp) WRITE(numout,*)
(   65)          IF(lwp) WRITE(numout,*) 'wzv     : vertical velocity from continuity eq.'
(   66)          IF(lwp) WRITE(numout,*) '~~~~~~~   auto-tasking case : j-k-i loop '
(   67) 
(   68)          ! bottom boundary condition: w=0 (set once for all)
(   69)          wn(:,:,jpk) = 0.e0
(   70)       ENDIF
(   71) 
(   72)       !                                                ! ===============
(   73)       DO jj = 1, jpj                                   !  Vertical slab
(   74)          !                                             ! ===============
(   75)          ! Computation from the bottom
(   76)          DO jk = jpkm1, 1, -1
(   77)             wn(:,jj,jk) = wn(:,jj,jk+1) - fse3t(:,jj,jk) * hdivn(:,jj,jk)
(   78)          END DO
(   79)          !                                             ! ===============
(   80)       END DO                                           !   End of slab
(   81)       !                                                ! ===============
(   82) 
(   83)       IF(ln_ctl)   CALL prt_ctl(tab3d_1=wn, clinfo1=' w**2 -   : ', mask1=wn)
(   84) 
(   85)    END SUBROUTINE wzv
(   86) 
(   87) #else
(   88)    !!----------------------------------------------------------------------
(   89)    !!   Default option                                           k-j-i loop
(   90)    !!----------------------------------------------------------------------
(   91) 
(   92)    SUBROUTINE wzv( kt )
(   93)       !!----------------------------------------------------------------------
(   94)       !!                    ***  ROUTINE wzv  ***
(   95)       !!
(   96)       !! ** Purpose :   Compute the now vertical velocity after the array swap
(   97)       !!
(   98)       !! ** Method  :   Using the incompressibility hypothesis, the vertical
(   99)       !!      velocity is computed by integrating the horizontal divergence 
(  100)       !!      from the bottom to the surface.
(  101)       !!        The boundary conditions are w=0 at the bottom (no flux) and,
(  102)       !!      in regid-lid case, w=0 at the sea surface.
(  103)       !!
(  104)       !! ** action  :   wn array : the now vertical velocity






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 3

(  105)       !!
(  106)       !! History :
(  107)       !!   9.0  !  02-07  (G. Madec)  Vector optimization
(  108)       !!----------------------------------------------------------------------
(  109)       !! * Arguments
(  110)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  111) 
(  112)       !! * Local declarations
(  113)       INTEGER ::   jk          ! dummy loop indices
(  114)       !!----------------------------------------------------------------------
(  115)       !!  OPA 8.5, LODYC-IPSL (2002)
(  116)       !!----------------------------------------------------------------------
(  117) 
(  118)       IF( kt == nit000 ) THEN
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*) 'wzv     : vertical velocity from continuity eq.'
(  121)          IF(lwp) WRITE(numout,*) '~~~~~~~ ' 
(  122) 
(  123)          ! bottom boundary condition: w=0 (set once for all)
(  124)          wn(:,:,jpk) = 0.e0
(  125)       ENDIF
(  126) 
(  127)       ! Computation from the bottom
(  128)       DO jk = jpkm1, 1, -1
(  129)          wn(:,:,jk) = wn(:,:,jk+1) - fse3t(:,:,jk) * hdivn(:,:,jk)
(  130)       END DO
(  131) 
(  132)       IF(ln_ctl)   CALL prt_ctl(tab3d_1=wn, clinfo1=' w**2 -   : ', mask1=wn)
(  133) 
(  134)    END SUBROUTINE wzv
(  135) #endif
(  136) 
(  137)    !!======================================================================
(  138) END MODULE wzvmod






























PGF90 (Version     10.5)          09/23/2010  09:50:03      page 4

(    1) # 1 "wzvmod.F90"
(    1) MODULE wzvmod
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  wzvmod  ***
(    4)    !! Ocean diagnostic variable : vertical velocity
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   wzv        : Compute the vertical velocity
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE in_out_manager  ! I/O manager
(   14)    USE prtctl          ! Print control
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC wzv       ! routine called by step.F90 and inidtr.F90
(   21) 
(   22)    !! * Substitutions
(   23) # 1 "./domzgr_substitute.h90"
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
(   48) # 24 "wzvmod.F90"
(   24) # 24
(   24)    !!----------------------------------------------------------------------
(   25) 
(   26) CONTAINS
(   27) 
(   28) # 88
(   88)    !!----------------------------------------------------------------------
(   89)    !!   Default option                                           k-j-i loop
(   90)    !!----------------------------------------------------------------------
(   91) 
(   92)    SUBROUTINE wzv( kt )
(   93)       !!----------------------------------------------------------------------
(   94)       !!                    ***  ROUTINE wzv  ***
(   95)       !!
(   96)       !! ** Purpose :   Compute the now vertical velocity after the array swap
(   97)       !!
(   98)       !! ** Method  :   Using the incompressibility hypothesis, the vertical
(   99)       !!      velocity is computed by integrating the horizontal divergence 
(  100)       !!      from the bottom to the surface.






PGF90 (Version     10.5)          09/23/2010  09:50:03      page 5

(  101)       !!        The boundary conditions are w=0 at the bottom (no flux) and,
(  102)       !!      in regid-lid case, w=0 at the sea surface.
(  103)       !!
(  104)       !! ** action  :   wn array : the now vertical velocity
(  105)       !!
(  106)       !! History :
(  107)       !!   9.0  !  02-07  (G. Madec)  Vector optimization
(  108)       !!----------------------------------------------------------------------
(  109)       !! * Arguments
(  110)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(  111) 
(  112)       !! * Local declarations
(  113)       INTEGER ::   jk          ! dummy loop indices
(  114)       !!----------------------------------------------------------------------
(  115)       !!  OPA 8.5, LODYC-IPSL (2002)
(  116)       !!----------------------------------------------------------------------
(  117) 
(  118)       IF( kt == nit000 ) THEN
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*) 'wzv     : vertical velocity from continuity eq.'
(  121)          IF(lwp) WRITE(numout,*) '~~~~~~~ ' 
(  122) 
(  123)          ! bottom boundary condition: w=0 (set once for all)
(  124)          wn(:,:,jpk) = 0.e0
(  125)       ENDIF
(  126) 
(  127)       ! Computation from the bottom
(  128)       DO jk = jpkm1, 1, -1
(  129)          wn(:,:,jk) = wn(:,:,jk+1) - e3t(jk) * hdivn(:,:,jk)
(  130)       END DO
(  131) 
(  132)       IF(ln_ctl)   CALL prt_ctl(tab3d_1=wn, clinfo1=' w**2 -   : ', mask1=wn)
(  133) 
(  134)    END SUBROUTINE wzv
(  135) 
(  136) # 137
(  137)    !!======================================================================
(  138) END MODULE wzvmod
