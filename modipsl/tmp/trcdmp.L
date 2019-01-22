


PGF90 (Version     10.5)          09/23/2010  09:49:19      page 1

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

Filename: trcdmp.F90

(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) !!NB: Code is broken if tracer damping is desired 
(    3) !!     ====> must modify this code, or recover older version
(    4) MODULE trcdmp
(    5)    !!======================================================================
(    6)    !!                       ***  MODULE  trcdmp  ***
(    7)    !! Ocean physics: internal restoring trend on passive tracers
(    8)    !!======================================================================
(    9) #if  defined key_passivetrc && defined key_trcdmp 
(   10)   USE in_out_manager
(   11) !!DB
(   12)    LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
(   13) CONTAINS
(   14)    SUBROUTINE trc_dmp( kt )        ! Empty routine
(   15)       INTEGER, INTENT(in) :: kt
(   16) 
(   17)       if(lwp) WRITE(numout,*) 'DBG -- trc_dmp: key_trcdmp inactivated ===> code will not work', kt
(   18) 
(   19)    END SUBROUTINE trc_dmp
(   20) 
(   21) #else
(   22)    !!----------------------------------------------------------------------
(   23)    !!   Default key                                     NO internal damping
(   24)    !!----------------------------------------------------------------------
(   25)    LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
(   26) CONTAINS
(   27)    SUBROUTINE trc_dmp( kt )        ! Empty routine
(   28)       INTEGER, INTENT(in) :: kt
(   29) !      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
(   30)    END SUBROUTINE trc_dmp
(   31) #endif
(   32) 
(   33)    !!======================================================================
(   34) END MODULE trcdmp


















PGF90 (Version     10.5)          09/23/2010  09:49:19      page 2

(    1) # 1 "trcdmp.F90"
(    1) !!DB 2009.08.24 -- eliminate non BGCM code options 
(    2) !!NB: Code is broken if tracer damping is desired 
(    3) !!     ====> must modify this code, or recover older version
(    4) MODULE trcdmp
(    5)    !!======================================================================
(    6)    !!                       ***  MODULE  trcdmp  ***
(    7)    !! Ocean physics: internal restoring trend on passive tracers
(    8)    !!======================================================================
(    9) # 22
(   22)    !!----------------------------------------------------------------------
(   23)    !!   Default key                                     NO internal damping
(   24)    !!----------------------------------------------------------------------
(   25)    LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
(   26) CONTAINS
(   27)    SUBROUTINE trc_dmp( kt )        ! Empty routine
(   28)       INTEGER, INTENT(in) :: kt
(   29) !      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
(   30)    END SUBROUTINE trc_dmp
(   31) 
(   32) # 33
(   33)    !!======================================================================
(   34) END MODULE trcdmp
