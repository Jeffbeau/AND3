


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

Filename: sms.F90

(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE sms
(    3)    !!======================================================================
(    4)    !!                        ***  sms  ***
(    5)    !! passive tracers :   set the passive tracers variables
(    6)    !!======================================================================
(    7)    !! History :
(    8)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(    9)    !!----------------------------------------------------------------------
(   10)    !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/sms.F90,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15) #if defined key_passivetrc
(   16) 
(   17)    USE par_oce
(   18)    USE par_trc
(   19) !!DB
(   20) !   USE par_sms
(   21)    IMPLICIT NONE
(   22) 
(   23) #else
(   24)    !!======================================================================
(   25)    !!  Empty module : No passive tracer 
(   26)    !!======================================================================
(   27) #endif
(   28) 
(   29) END MODULE sms























PGF90 (Version     10.5)          09/23/2010  09:49:15      page 2

(    1) # 1 "sms.F90"
(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE sms
(    3)    !!======================================================================
(    4)    !!                        ***  sms  ***
(    5)    !! passive tracers :   set the passive tracers variables
(    6)    !!======================================================================
(    7)    !! History :
(    8)    !!   9.0  !  04-03  (C. Ethe)  Free form and module
(    9)    !!----------------------------------------------------------------------
(   10)    !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(   11)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/sms.F90,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
(   12)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15) # 24
(   24)    !!======================================================================
(   25)    !!  Empty module : No passive tracer 
(   26)    !!======================================================================
(   27) 
(   28) # 29
(   29) END MODULE sms
