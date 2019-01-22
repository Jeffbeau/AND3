


PGF90 (Version     10.5)          09/23/2010  09:49:14      page 1

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

Filename: par_trc_trp.F90

(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE par_trc_trp
(    3)    !!======================================================================
(    4)    !!                        ***  par_trc_trp  ***
(    5)    !! passive tracers :   set the number of passive tracers
(    6)    !!======================================================================
(    7)    !! History :
(    8)    !!   9.0  !  04-03  (C. Ethe)  Orignal
(    9)    !!----------------------------------------------------------------------
(   10)    !!  TOP 1.0,  LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc_trp.F90,v 1.7 2006/04/10 15:40:28 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14) #if defined key_passivetrc
(   15)    !!---------------------------------------------------------------------
(   16)    !!   'key_passivetrc'   :                               Passive tracer
(   17)    !!---------------------------------------------------------------------
(   18) 
(   19)    IMPLICIT NONE
(   20)    PUBLIC
(   21)    
(   22)    !! jptra   : number of passive tracers
(   23)    !! jpdia2d : additional 2d output
(   24)    !! jpdia3d : additional 3d output
(   25) 
(   26) #if defined key_BGCM_01
(   27)    !!---------------------------------------------------------------------
(   28)    !!   'key_BGCM_01     model in development
(   29)    !!---------------------------------------------------------------------
(   30)    INTEGER, PUBLIC, PARAMETER :: jptra   = 3  !! at this time
(   31) #if defined key_trc_diaadd    
(   32)    INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
(   33)    INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
(   34) #endif
(   35) 
(   36) #else
(   37)    !!---------------------------------------------------------------------
(   38)    !!   'default'   :          temperature and salinity as passive tracers
(   39)    !!---------------------------------------------------------------------
(   40)    INTEGER, PUBLIC, PARAMETER :: jptra   = 2
(   41) #if defined key_trc_diaadd
(   42)    INTEGER, PUBLIC, PARAMETER :: jpdia2d = 1
(   43)    INTEGER, PUBLIC, PARAMETER :: jpdia3d = 1
(   44) #endif
(   45) #endif
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:14      page 2

(   47) #else
(   48)    !!======================================================================
(   49)    !!  Empty module : No passive tracer 
(   50)    !!======================================================================
(   51) #endif
(   52) 
(   53) END MODULE par_trc_trp

























































PGF90 (Version     10.5)          09/23/2010  09:49:14      page 3

(    1) # 1 "par_trc_trp.F90"
(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE par_trc_trp
(    3)    !!======================================================================
(    4)    !!                        ***  par_trc_trp  ***
(    5)    !! passive tracers :   set the number of passive tracers
(    6)    !!======================================================================
(    7)    !! History :
(    8)    !!   9.0  !  04-03  (C. Ethe)  Orignal
(    9)    !!----------------------------------------------------------------------
(   10)    !!  TOP 1.0,  LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc_trp.F90,v 1.7 2006/04/10 15:40:28 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14) # 48
(   48)    !!======================================================================
(   49)    !!  Empty module : No passive tracer 
(   50)    !!======================================================================
(   51) 
(   52) # 53
(   53) END MODULE par_trc_trp
