


PGF90 (Version     10.5)          09/23/2010  09:49:02      page 1

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

Filename: dynspg_oce.F90

(    1) MODULE dynspg_oce
(    2)    !!----------------------------------------------------------------------
(    3)    !!                       ***  MODULE dynspg_oce  ***
(    4)    !!       
(    5)    !! ** Purpose :   Define in memory all the ocean space domain variables
(    6)    !!----------------------------------------------------------------------
(    7)    !! Modules used
(    8)    USE par_oce          ! ocean parameters
(    9) 
(   10)    IMPLICIT NONE
(   11)    PUBLIC           
(   12)    !!----------------------------------------------------------------------
(   13)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_oce.F90,v 1.1 2005/12/28 09:25:06 opalod Exp $ 
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   16)    !!----------------------------------------------------------------------
(   17) 
(   18)    !! Surface pressure gradient logicals
(   19)    !! ----------------------------------
(   20) #if   defined key_dynspg_exp   ||  defined key_esopa
(   21)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .TRUE.  !: Explicit free surface flag
(   22) #else
(   23)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .FALSE. !: Explicit free surface flag
(   24) #endif
(   25) #if   defined key_dynspg_ts   ||  defined key_esopa
(   26)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .TRUE.  !: Free surface with time splitting flag
(   27) #else
(   28)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .FALSE. !: Free surface with time splitting flag
(   29) #endif
(   30) #if   defined key_dynspg_flt  ||  defined key_esopa
(   31)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .TRUE.  !: Filtered free surface cst volume flag
(   32) #else
(   33)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .FALSE. !: Filtered free surface cst volume flag
(   34) #endif
(   35) #if   defined key_dynspg_rl
(   36)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl  = .TRUE.  !: Rigid-lid flag
(   37) #else
(   38)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl  = .FALSE. !: Rigid-lid flag
(   39) #endif
(   40) 
(   41) #if   defined key_dynspg_ts   ||  defined key_esopa
(   42)    !! Time splitting variables
(   43)    !! ------------------------
(   44)       REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables averaged over the barotropic loop
(   45)          sshn_b, sshb_b,               &  ! sea surface heigth (now, before)
(   46)          un_b  , vn_b                     ! vertically integrated horizontal velocities (now)






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 2

(   47)       REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables of the explicit barotropic loop
(   48)          sshn_e, ssha_e,               &  ! sea surface heigth (now,after)
(   49)          ua_e  , va_e                     ! vertically integrated horizontal velocities (after)
(   50) #endif
(   51) 
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) END MODULE dynspg_oce
























































PGF90 (Version     10.5)          09/23/2010  09:49:02      page 3

(    1) # 1 "dynspg_oce.F90"
(    1) MODULE dynspg_oce
(    2)    !!----------------------------------------------------------------------
(    3)    !!                       ***  MODULE dynspg_oce  ***
(    4)    !!       
(    5)    !! ** Purpose :   Define in memory all the ocean space domain variables
(    6)    !!----------------------------------------------------------------------
(    7)    !! Modules used
(    8)    USE par_oce          ! ocean parameters
(    9) 
(   10)    IMPLICIT NONE
(   11)    PUBLIC           
(   12)    !!----------------------------------------------------------------------
(   13)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   14)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_oce.F90,v 1.1 2005/12/28 09:25:06 opalod Exp $ 
(   15)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   16)    !!----------------------------------------------------------------------
(   17) 
(   18)    !! Surface pressure gradient logicals
(   19)    !! ----------------------------------
(   20) # 23
(   23)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .FALSE. !: Explicit free surface flag
(   24) # 26
(   26)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .TRUE.  !: Free surface with time splitting flag
(   27) # 33
(   33)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .FALSE. !: Filtered free surface cst volume flag
(   34) # 38
(   38)    LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl  = .FALSE. !: Rigid-lid flag
(   39) 
(   40) # 42
(   42)    !! Time splitting variables
(   43)    !! ------------------------
(   44)       REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables averaged over the barotropic loop
(   45)          sshn_b, sshb_b,               &  ! sea surface heigth (now, before)
(   46)          un_b  , vn_b                     ! vertically integrated horizontal velocities (now)
(   47)       REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables of the explicit barotropic loop
(   48)          sshn_e, ssha_e,               &  ! sea surface heigth (now,after)
(   49)          ua_e  , va_e                     ! vertically integrated horizontal velocities (after)
(   50) 
(   51) # 52
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) END MODULE dynspg_oce
