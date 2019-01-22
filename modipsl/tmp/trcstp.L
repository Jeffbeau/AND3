


PGF90 (Version     10.5)          09/23/2010  09:50:18      page 1

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

Filename: trcstp.F90

(    1) !!DB 2009.08 
(    2) !!Eliminate non-BGCM code options
(    3) MODULE trcstp
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE trcstp  ***
(    6)    !! Time-stepping    : time loop of opa for passive tracer
(    7)    !!======================================================================
(    8) 
(    9) 
(   10) #if defined key_passivetrc
(   11) 
(   12) !!DB
(   13) #if defined key_BGCM_01
(   14) 
(   15)   USE bgcm_01_model
(   16) 
(   17) #else
(   18)    !!----------------------------------------------------------------------
(   19)    !!   Default key                                     NO passive tracers
(   20)    !!----------------------------------------------------------------------
(   21) CONTAINS
(   22)    SUBROUTINE trc_stp( kt )        ! Empty routine
(   23)       INTEGER, INTENT(in) :: kt
(   24) !      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
(   25)    END SUBROUTINE trc_stp
(   26) #endif
(   27) #endif
(   28)    !!======================================================================
(   29) END MODULE trcstp























PGF90 (Version     10.5)          09/23/2010  09:50:18      page 2

(    1) # 1 "trcstp.F90"
(    1) !!DB 2009.08 
(    2) !!Eliminate non-BGCM code options
(    3) MODULE trcstp
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE trcstp  ***
(    6)    !! Time-stepping    : time loop of opa for passive tracer
(    7)    !!======================================================================
(    8) 
(    9) 
(   10) # 28
(   28)    !!======================================================================
(   29) END MODULE trcstp
