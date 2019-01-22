


PGF90 (Version     10.5)          09/23/2010  09:49:15      page 1

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

Filename: par_trc.F90

(    1) MODULE par_trc
(    2)    !!======================================================================
(    3)    !!                        ***  par_trc  ***
(    4)    !! passive tracers :   set the passive tracers parameters
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.2  !  96-01  (M. Levy)  Original code
(    8)    !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
(    9)    !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
(   10)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  TOP 1.0,  LOCEAN-IPSL (2005)
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc.F90,v 1.4 2005/09/12 09:04:53 opalod Exp $
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17) #if defined key_passivetrc
(   18) 
(   19)    USE par_trc_trp
(   20) 
(   21)    IMPLICIT NONE
(   22)    PUBLIC
(   23) 
(   24) 
(   25) #if defined key_trc_diatrd
(   26) 
(   27) !! number of dynamical trends
(   28) #  if defined key_trc_ldfeiv
(   29) !! we keep 3 more trends for eddy induced flux (gent velocity)
(   30) #    if defined key_trcdmp
(   31)    INTEGER , PARAMETER :: jpdiatrc = 10
(   32) #    else
(   33)    INTEGER , PARAMETER :: jpdiatrc = 9
(   34) #    endif
(   35) #  else
(   36) #    if defined key_trcdmp
(   37)    INTEGER , PARAMETER :: jpdiatrc = 7
(   38) #    else
(   39)    INTEGER , PARAMETER :: jpdiatrc = 6
(   40) #    endif
(   41) #  endif
(   42) # endif
(   43) 
(   44) #else
(   45)    !!======================================================================
(   46)    !!  Empty module : No passive tracer 






PGF90 (Version     10.5)          09/23/2010  09:49:15      page 2

(   47)    !!======================================================================
(   48) #endif
(   49) 
(   50) END MODULE par_trc




























































PGF90 (Version     10.5)          09/23/2010  09:49:15      page 3

(    1) # 1 "par_trc.F90"
(    1) MODULE par_trc
(    2)    !!======================================================================
(    3)    !!                        ***  par_trc  ***
(    4)    !! passive tracers :   set the passive tracers parameters
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.2  !  96-01  (M. Levy)  Original code
(    8)    !!        !  99-07  (M. Levy)  for LOBSTER1 or NPZD model
(    9)    !!        !  00-04  (O. Aumont, M.A. Foujols)  HAMOCC3 and P3ZD
(   10)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(   11)    !!----------------------------------------------------------------------
(   12)    !!  TOP 1.0,  LOCEAN-IPSL (2005)
(   13)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/par_trc.F90,v 1.4 2005/09/12 09:04:53 opalod Exp $
(   14)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   15)    !!----------------------------------------------------------------------
(   16)    !! * Modules used
(   17) # 45
(   45)    !!======================================================================
(   46)    !!  Empty module : No passive tracer 
(   47)    !!======================================================================
(   48) 
(   49) # 50
(   50) END MODULE par_trc
