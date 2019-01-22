


PGF90 (Version     10.5)          09/23/2010  09:49:35      page 1

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

Filename: dom_ice.F90

(    1) MODULE dom_ice
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dom_ice  ***
(    4)    !! LIM Sea Ice :   Domain  variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   2.0  !  03-08  (C. Ethe)  Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/dom_ice.F90,v 1.5 2006/03/21 08:42:22 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_ice
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Share module variables
(   20)    LOGICAL, PUBLIC ::       &  !:
(   21)       l_jeq     = .TRUE. ,  &  !: Equator inside the domain flag
(   22)       ln_limini = .FALSE.,  &  !: Ice initialization state
(   23)       ln_limdmp = .FALSE.      !: Ice damping
(   24) 
(   25)    INTEGER, PUBLIC ::   &  !:
(   26)       njeq , njeqm1        !: j-index of the equator if it is inside the domain
(   27)       !                    !  (otherwise = jpj+10 (SH) or -10 (SH) )
(   28) 
(   29)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   30)       fs2cor ,          &  !: coriolis factor
(   31)       fcor   ,          &  !: coriolis coefficient
(   32)       covrai ,          &  !: sine of geographic latitude
(   33)       area   ,          &  !: surface of grid cell 
(   34)       tms    , tmu         !: temperature and velocity points masks
(   35) 
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2) ::   &  !:
(   37)       wght   ,          &  !: weight of the 4 neighbours to compute averages
(   38)       akappa ,          &  !: first group of metric coefficients
(   39)       bkappa               !: third group of metric coefficients
(   40) 
(   41)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2,2,2) ::   &  !:
(   42)       alambd               !: second group of metric coefficients
(   43) 
(   44)    !!======================================================================
(   45) END MODULE dom_ice







PGF90 (Version     10.5)          09/23/2010  09:49:35      page 2

(    1) # 1 "dom_ice.F90"
(    1) MODULE dom_ice
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE  dom_ice  ***
(    4)    !! LIM Sea Ice :   Domain  variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   2.0  !  03-08  (C. Ethe)  Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/dom_ice.F90,v 1.5 2006/03/21 08:42:22 opalod Exp $
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE par_ice
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Share module variables
(   20)    LOGICAL, PUBLIC ::       &  !:
(   21)       l_jeq     = .TRUE. ,  &  !: Equator inside the domain flag
(   22)       ln_limini = .FALSE.,  &  !: Ice initialization state
(   23)       ln_limdmp = .FALSE.      !: Ice damping
(   24) 
(   25)    INTEGER, PUBLIC ::   &  !:
(   26)       njeq , njeqm1        !: j-index of the equator if it is inside the domain
(   27)       !                    !  (otherwise = jpj+10 (SH) or -10 (SH) )
(   28) 
(   29)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   30)       fs2cor ,          &  !: coriolis factor
(   31)       fcor   ,          &  !: coriolis coefficient
(   32)       covrai ,          &  !: sine of geographic latitude
(   33)       area   ,          &  !: surface of grid cell 
(   34)       tms    , tmu         !: temperature and velocity points masks
(   35) 
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2) ::   &  !:
(   37)       wght   ,          &  !: weight of the 4 neighbours to compute averages
(   38)       akappa ,          &  !: first group of metric coefficients
(   39)       bkappa               !: third group of metric coefficients
(   40) 
(   41)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2,2,2) ::   &  !:
(   42)       alambd               !: second group of metric coefficients
(   43) 
(   44)    !!======================================================================
(   45) END MODULE dom_ice
