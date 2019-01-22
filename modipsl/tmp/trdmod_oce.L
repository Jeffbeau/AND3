


PGF90 (Version     10.5)          09/23/2010  09:49:08      page 1

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

Filename: trdmod_oce.F90

(    1) MODULE trdmod_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdmod_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmod_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE trdicp_oce              ! ocean momentum/tracers bassin properties trends variables
(   13)    USE trdmld_oce              ! ocean active mixed layer tracers trends variables
(   14)    USE trdvor_oce              ! ocean vorticity trends variables
(   15) 
(   16)    !! Control parameters
(   17)    !!----------------------------------------------------------------------
(   18)    LOGICAL, PUBLIC ::   l_trdtra = .FALSE.    !: tracers  trend flag
(   19)    LOGICAL, PUBLIC ::   l_trddyn = .FALSE.    !: momentum trend flag
(   20) 
(   21)   !!======================================================================
(   22) END MODULE trdmod_oce






























PGF90 (Version     10.5)          09/23/2010  09:49:08      page 2

(    1) # 1 "trdmod_oce.F90"
(    1) MODULE trdmod_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdmod_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmod_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE trdicp_oce              ! ocean momentum/tracers bassin properties trends variables
(   13)    USE trdmld_oce              ! ocean active mixed layer tracers trends variables
(   14)    USE trdvor_oce              ! ocean vorticity trends variables
(   15) 
(   16)    !! Control parameters
(   17)    !!----------------------------------------------------------------------
(   18)    LOGICAL, PUBLIC ::   l_trdtra = .FALSE.    !: tracers  trend flag
(   19)    LOGICAL, PUBLIC ::   l_trddyn = .FALSE.    !: momentum trend flag
(   20) 
(   21)   !!======================================================================
(   22) END MODULE trdmod_oce
