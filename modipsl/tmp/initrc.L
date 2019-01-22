


PGF90 (Version     10.5)          09/23/2010  09:50:02      page 1

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

Filename: initrc.F90

(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE initrc
(    3)    !!================================================
(    4)    !!
(    5)    !!                       *** MODULE initrc ***
(    6)    !! Initialisation the tracer model
(    7)    !!================================================
(    8)                                                                                                                             
(    9) #if defined key_passivetrc
(   10) 
(   11) #if defined key_BGCM_01
(   12) !!lib_bgcm_01 contains ini_trc() and associated routines
(   13)   USE lib_bgcm_01
(   14) 
(   15) #else
(   16)    !!======================================================================
(   17)    !!  Empty module : No passive tracer
(   18)    !!======================================================================
(   19) CONTAINS
(   20)    SUBROUTINE ini_trc
(   21)       
(   22)    END SUBROUTINE ini_trc
(   23) #endif
(   24) 
(   25) #endif
(   26) 
(   27) END MODULE initrc 
(   28) 
(   29) 























PGF90 (Version     10.5)          09/23/2010  09:50:02      page 2

(    1) # 1 "initrc.F90"
(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) MODULE initrc
(    3)    !!================================================
(    4)    !!
(    5)    !!                       *** MODULE initrc ***
(    6)    !! Initialisation the tracer model
(    7)    !!================================================
(    8)                                                                                                                             
(    9) 
(   10) # 27
(   27) END MODULE initrc 
(   28) 
(   29) 
