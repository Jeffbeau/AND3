


PGF90 (Version     10.5)          09/23/2010  09:49:32      page 1

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

Filename: ini1d.F90

(    1) MODULE ini1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ini1D  ***
(    4)    !! Ocean state   :  1D initialization
(    5)    !!=====================================================================
(    6) #if defined key_cfg_1d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'               1D Configuration
(    9)    !!----------------------------------------------------------------------
(   10)    !!   init_1d   : initial mask
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE dom_oce         ! ocean space and time domain 
(   14)    USE phycst
(   15)    USE in_out_manager
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Share Module variables
(   21)    LOGICAL, PUBLIC, PARAMETER ::  lk_cfg_1d = .TRUE.       !: 1D flag
(   22) 
(   23)    !! * Routine accessibility
(   24)    PUBLIC init_1d   ! routine called by OPA.F90
(   25) 
(   26)    !!----------------------------------------------------------------------
(   27)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   28)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/ini1d.F90,v 1.2 2005/09/02 15:33:59 opalod Exp $ 
(   29)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) CONTAINS
(   33) 
(   34)    SUBROUTINE init_1d
(   35)       !!----------------------------------------------------------------------
(   36)       !!                   ***  ROUTINE init_1d  ***
(   37)       !! 
(   38)       !! ** Purpose :   Re-Initialization of masks on 1D configuration
(   39)       !!
(   40)       !! ** Method  :
(   41)       !!
(   42)       !! History :
(   43)       !!   9.0  !  04-09  (C. Ethe) 1D configuration
(   44)       !!----------------------------------------------------------------------
(   45)       !! * Local declarations
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 2

(   47) 
(   48)       IF(lwp) WRITE(numout,*)
(   49)       IF(lwp) WRITE(numout,*) 'init_1d :  masks on 1D configuration'
(   50)       IF(lwp) WRITE(numout,*) '~~~~~~~~~'
(   51) 
(   52)       ! set umask and vmask equal tmask in 1D configuration
(   53)       umask(:,:,:) = tmask(:,:,:)
(   54)       vmask(:,:,:) = tmask(:,:,:)     
(   55) 
(   56)    END SUBROUTINE init_1d
(   57) 
(   58) #else
(   59)    !!----------------------------------------------------------------------
(   60)    !!   Default key                                     NO 1D Config
(   61)    !!----------------------------------------------------------------------
(   62)    LOGICAL , PUBLIC, PARAMETER ::   lk_cfg_1d = .FALSE.    !: internal damping flag
(   63) CONTAINS
(   64)    SUBROUTINE init_1d       ! Empty routine
(   65) 
(   66)    END SUBROUTINE init_1d
(   67) #endif
(   68) 
(   69)    !!=====================================================================
(   70) END MODULE ini1d








































PGF90 (Version     10.5)          09/23/2010  09:49:32      page 3

(    1) # 1 "ini1d.F90"
(    1) MODULE ini1d
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  ini1D  ***
(    4)    !! Ocean state   :  1D initialization
(    5)    !!=====================================================================
(    6) # 59
(   59)    !!----------------------------------------------------------------------
(   60)    !!   Default key                                     NO 1D Config
(   61)    !!----------------------------------------------------------------------
(   62)    LOGICAL , PUBLIC, PARAMETER ::   lk_cfg_1d = .FALSE.    !: internal damping flag
(   63) CONTAINS
(   64)    SUBROUTINE init_1d       ! Empty routine
(   65) 
(   66)    END SUBROUTINE init_1d
(   67) 
(   68) # 69
(   69)    !!=====================================================================
(   70) END MODULE ini1d
