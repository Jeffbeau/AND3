


PGF90 (Version     10.5)          09/23/2010  09:48:59      page 1

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

Filename: par_ice.F90

(    1) MODULE par_ice
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE par_ice   ***
(    4)    !! Sea-Ice model : definition of the parameters
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/par_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_oce
(   13) 
(   14)    IMPLICIT NONE
(   15)    PUBLIC               ! allows par_oce and par_kind to be known in ice modules
(   16) 
(   17)    INTEGER, PUBLIC, PARAMETER ::   &  !:
(   18)       jpkmax =  1    ,      &  !: ???
(   19)       jpsmax =  2              !: ???
(   20) 
(   21)    INTEGER, PUBLIC, PARAMETER ::   &  !: 
(   22)       jplayers   = 2 ,           &  !: number of vertical ice layers
(   23)       jplayersp1 = jplayers + 1     !: ???
(   24) 
(   25)    !!======================================================================
(   26) END MODULE par_ice


























PGF90 (Version     10.5)          09/23/2010  09:48:59      page 2

(    1) # 1 "par_ice.F90"
(    1) MODULE par_ice
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE par_ice   ***
(    4)    !! Sea-Ice model : definition of the parameters
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  LIM 2.0, UCL-LOCEAN-IPSL (2005)
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/par_ice.F90,v 1.4 2005/03/27 18:34:42 opalod Exp $
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_oce
(   13) 
(   14)    IMPLICIT NONE
(   15)    PUBLIC               ! allows par_oce and par_kind to be known in ice modules
(   16) 
(   17)    INTEGER, PUBLIC, PARAMETER ::   &  !:
(   18)       jpkmax =  1    ,      &  !: ???
(   19)       jpsmax =  2              !: ???
(   20) 
(   21)    INTEGER, PUBLIC, PARAMETER ::   &  !: 
(   22)       jplayers   = 2 ,           &  !: number of vertical ice layers
(   23)       jplayersp1 = jplayers + 1     !: ???
(   24) 
(   25)    !!======================================================================
(   26) END MODULE par_ice
