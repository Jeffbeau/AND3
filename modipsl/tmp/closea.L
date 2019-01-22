


PGF90 (Version     10.5)          09/23/2010  09:49:27      page 1

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

Filename: closea.F90

(    1) !!DB: 2009.08.31 -- eliminated most code
(    2) MODULE closea
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  closea  ***
(    5)    !! Closed Seas  : 
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dom_clo    : modification of the ocean domain for closed seas cases
(   10)    !!   flx_clo    : Special handling of closed seas
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE ocesbc          ! ocean surface boundary conditions (fluxes)
(   17)    USE flxrnf          ! runoffs
(   18)    USE lib_mpp         ! distributed memory computing library
(   19)    USE lbclnk          ! ???
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC dom_clo      ! routine called by dom_init
(   26)    PUBLIC flx_clo      ! routine called by step
(   27) 
(   28)    !! * Share module variables
(   29)    INTEGER, PUBLIC, PARAMETER ::   &  !:
(   30) !      jpncs   = 4               !: number of closed sea
(   31)       jpncs   = 0               !: number of closed sea
(   32)    INTEGER, PUBLIC ::          & !!: namclo : closed seas and lakes
(   33)       nclosea =  0                !: = 0 no closed sea or lake
(   34)       !                           !  = 1 closed sea or lake in the domain
(   35)    INTEGER, PUBLIC, DIMENSION (jpncs) ::   &  !:
(   36)       ncstt,           &  !: Type of closed sea
(   37)       ncsi1, ncsj1,    &  !: closed sea limits                                                                 
(   38)       ncsi2, ncsj2,    &  !: 
(   39)       ncsnr               !: number of point where run-off pours
(   40)    INTEGER, PUBLIC, DIMENSION (jpncs,4) ::   &
(   41)       ncsir, ncsjr        !: Location of run-off
(   42) 
(   43)    !! * Module variable
(   44)    REAL(wp), DIMENSION (jpncs+1) ::   &
(   45)       surf               ! closed sea surface
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 2

(   47)    !! * Substitutions
(   48) #  include "vectopt_loop_substitute.h90"
(   49)    !!----------------------------------------------------------------------
(   50)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/closea.F90,v 1.5 2005/09/22 14:23:57 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55) CONTAINS
(   56) 
(   57)    SUBROUTINE dom_clo
(   58)       !!---------------------------------------------------------------------
(   59)       !!                  ***  ROUTINE dom_clo  ***
(   60)       !!        
(   61)       !! ** Purpose :   Closed sea domain initialization
(   62)       !!
(   63)       !! ** Method  :   if a closed sea is located only in a model grid point
(   64)       !!      just the thermodynamic processes are applied.
(   65)       !!
(   66)       !! ** Action :   ncsi1(), ncsj1() : south-west closed sea limits (i,j)
(   67)       !!               ncsi2(), ncsj2() : north-east Closed sea limits (i,j)
(   68)       !!               ncsir(), ncsjr() : Location of runoff
(   69)       !!               ncsnr            : number of point where run-off pours
(   70)       !!               ncstt            : Type of closed sea
(   71)       !!                                  =0 spread over the world ocean
(   72)       !!                                  =2 put at location runoff
(   73)       !!
(   74)       !! History :
(   75)       !!        !  01-04  (E. Durand)  Original code
(   76)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Local variables
(   79)       INTEGER ::   jc            ! dummy loop indices
(   80)       !!----------------------------------------------------------------------
(   81)       
(   82) 
(   83)       IF(lwp) WRITE(numout,*)
(   84)       IF(lwp) WRITE(numout,*)'DB: dom_clo : closed seas ---> blank routine'
(   85)       IF(lwp) WRITE(numout,*)'~~~~~~~'
(   86) 
(   87)    END SUBROUTINE dom_clo
(   88) 
(   89) 
(   90)    SUBROUTINE flx_clo( kt )
(   91)       !!---------------------------------------------------------------------
(   92)       !!                  ***  ROUTINE flx_clo  ***
(   93)       !!                    
(   94)       !! ** Purpose :   Special handling of closed seas
(   95)       !!
(   96)       !! ** Method  :   Water flux is forced to zero over closed sea
(   97)       !!      Excess is shared between remaining ocean, or
(   98)       !!      put as run-off in open ocean.
(   99)       !!
(  100)       !! ** Action :
(  101)       !!
(  102)       !! History :
(  103)       !!   8.2  !  00-05  (O. Marti)  Original code
(  104)       !!   8.5  !  02-07  (G. Madec)  Free form, F90






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 3

(  105)       !!----------------------------------------------------------------------
(  106)       !! * Arguments
(  107)       INTEGER, INTENT (in) :: kt
(  108) 
(  109)       !! * Local declarations
(  110)       REAL(wp), DIMENSION (jpncs) :: zemp
(  111)       INTEGER  :: ji, jj, jc, jn
(  112)       REAL(wp) :: zze2
(  113)       !!----------------------------------------------------------------------
(  114) 
(  115)       ! 1 - Initialisation
(  116)       ! ------------------
(  117) 
(  118)       IF( kt == nit000 ) THEN 
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*)'DB: flx_clo : closed seas --- blank routine'
(  121)          IF(lwp) WRITE(numout,*)'~~~~~~~'
(  122) 
(  123)          ! Total surface of ocean
(  124) !         surf(jpncs+1) = SUM( e1t(:,:) * e2t(:,:) * tmask_i(:,:) )
(  125)          surf(jpncs+1) = 0.e0
(  126)       endif
(  127) 
(  128) 
(  129)    END SUBROUTINE flx_clo
(  130) 
(  131)    !!======================================================================
(  132) END MODULE closea




































PGF90 (Version     10.5)          09/23/2010  09:49:27      page 4

(    1) # 1 "closea.F90"
(    1) !!DB: 2009.08.31 -- eliminated most code
(    2) MODULE closea
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  closea  ***
(    5)    !! Closed Seas  : 
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   dom_clo    : modification of the ocean domain for closed seas cases
(   10)    !!   flx_clo    : Special handling of closed seas
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE ocesbc          ! ocean surface boundary conditions (fluxes)
(   17)    USE flxrnf          ! runoffs
(   18)    USE lib_mpp         ! distributed memory computing library
(   19)    USE lbclnk          ! ???
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC dom_clo      ! routine called by dom_init
(   26)    PUBLIC flx_clo      ! routine called by step
(   27) 
(   28)    !! * Share module variables
(   29)    INTEGER, PUBLIC, PARAMETER ::   &  !:
(   30) !      jpncs   = 4               !: number of closed sea
(   31)       jpncs   = 0               !: number of closed sea
(   32)    INTEGER, PUBLIC ::          & !!: namclo : closed seas and lakes
(   33)       nclosea =  0                !: = 0 no closed sea or lake
(   34)       !                           !  = 1 closed sea or lake in the domain
(   35)    INTEGER, PUBLIC, DIMENSION (jpncs) ::   &  !:
(   36)       ncstt,           &  !: Type of closed sea
(   37)       ncsi1, ncsj1,    &  !: closed sea limits                                                                 
(   38)       ncsi2, ncsj2,    &  !: 
(   39)       ncsnr               !: number of point where run-off pours
PGF90-W-0435-Array declared with zero size (closea.F90: 35)
(   40)    INTEGER, PUBLIC, DIMENSION (jpncs,4) ::   &
(   41)       ncsir, ncsjr        !: Location of run-off
(   42) 
(   43)    !! * Module variable
PGF90-W-0435-Array declared with zero size (closea.F90: 40)
(   44)    REAL(wp), DIMENSION (jpncs+1) ::   &
(   45)       surf               ! closed sea surface
(   46) 
(   47)    !! * Substitutions
(   48) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 5

(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 49 "closea.F90"
(   49) # 49
(   49)    !!----------------------------------------------------------------------
(   50)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   51)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/closea.F90,v 1.5 2005/09/22 14:23:57 opalod Exp $ 
(   52)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   53)    !!----------------------------------------------------------------------
(   54) 
(   55) CONTAINS
(   56) 
(   57)    SUBROUTINE dom_clo
(   58)       !!---------------------------------------------------------------------
(   59)       !!                  ***  ROUTINE dom_clo  ***
(   60)       !!        
(   61)       !! ** Purpose :   Closed sea domain initialization
(   62)       !!
(   63)       !! ** Method  :   if a closed sea is located only in a model grid point
(   64)       !!      just the thermodynamic processes are applied.
(   65)       !!
(   66)       !! ** Action :   ncsi1(), ncsj1() : south-west closed sea limits (i,j)
(   67)       !!               ncsi2(), ncsj2() : north-east Closed sea limits (i,j)
(   68)       !!               ncsir(), ncsjr() : Location of runoff
(   69)       !!               ncsnr            : number of point where run-off pours
(   70)       !!               ncstt            : Type of closed sea
(   71)       !!                                  =0 spread over the world ocean
(   72)       !!                                  =2 put at location runoff
(   73)       !!
(   74)       !! History :
(   75)       !!        !  01-04  (E. Durand)  Original code
(   76)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Local variables
(   79)       INTEGER ::   jc            ! dummy loop indices
(   80)       !!----------------------------------------------------------------------
(   81)       
(   82) 
(   83)       IF(lwp) WRITE(numout,*)
(   84)       IF(lwp) WRITE(numout,*)'DB: dom_clo : closed seas ---> blank routine'
(   85)       IF(lwp) WRITE(numout,*)'~~~~~~~'
(   86) 
(   87)    END SUBROUTINE dom_clo
(   88) 
(   89) 
(   90)    SUBROUTINE flx_clo( kt )
(   91)       !!---------------------------------------------------------------------
(   92)       !!                  ***  ROUTINE flx_clo  ***
(   93)       !!                    
(   94)       !! ** Purpose :   Special handling of closed seas
(   95)       !!
(   96)       !! ** Method  :   Water flux is forced to zero over closed sea
(   97)       !!      Excess is shared between remaining ocean, or
(   98)       !!      put as run-off in open ocean.
(   99)       !!
(  100)       !! ** Action :






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 6

(  101)       !!
(  102)       !! History :
(  103)       !!   8.2  !  00-05  (O. Marti)  Original code
(  104)       !!   8.5  !  02-07  (G. Madec)  Free form, F90
(  105)       !!----------------------------------------------------------------------
(  106)       !! * Arguments
(  107)       INTEGER, INTENT (in) :: kt
(  108) 
(  109)       !! * Local declarations
(  110)       REAL(wp), DIMENSION (jpncs) :: zemp
PGF90-W-0435-Array declared with zero size (closea.F90: 110)
(  111)       INTEGER  :: ji, jj, jc, jn
(  112)       REAL(wp) :: zze2
(  113)       !!----------------------------------------------------------------------
(  114) 
(  115)       ! 1 - Initialisation
(  116)       ! ------------------
(  117) 
(  118)       IF( kt == nit000 ) THEN 
(  119)          IF(lwp) WRITE(numout,*)
(  120)          IF(lwp) WRITE(numout,*)'DB: flx_clo : closed seas --- blank routine'
(  121)          IF(lwp) WRITE(numout,*)'~~~~~~~'
(  122) 
(  123)          ! Total surface of ocean
(  124) !         surf(jpncs+1) = SUM( e1t(:,:) * e2t(:,:) * tmask_i(:,:) )
(  125)          surf(jpncs+1) = 0.e0
(  126)       endif
(  127) 
(  128) 
(  129)    END SUBROUTINE flx_clo
(  130) 
(  131)    !!======================================================================
(  132) END MODULE closea
