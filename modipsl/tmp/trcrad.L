


PGF90 (Version     10.5)          09/23/2010  09:49:22      page 1

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

Filename: trcrad.F90

(    1) !!DB: 2009.09.09 -- deleted old DBG code 
(    2) !!and disable the routine. If you want it then dig up an older version
(    3) MODULE trcrad
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  trcrad  ***
(    6)    !! Ocean passive tracers:  correction of negative concentrations
(    7)    !!======================================================================
(    8) #if defined key_passivetrc
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_rad    : correction of negative concentrations
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce_trc             ! ocean dynamics and tracers variables
(   14)    USE trc                 ! ocean passive tracers variables
(   15)    USE lib_mpp
(   16)    USE prtctl_trc          ! Print control for debbuging
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC trc_rad        ! routine called by trcstp.F90
(   23)    !! * Substitutions
(   24) #  include "passivetrc_substitute.h90"
(   25)    !!----------------------------------------------------------------------
(   26)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcrad.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) CONTAINS
(   31) 
(   32)    SUBROUTINE trc_rad( kt )
(   33)       !!----------------------------------------------------------------------
(   34)       !!                  ***  ROUTINE trc_rad  ***
(   35)       !!
(   36)       !! ** Purpose : "crappy" routine to correct artificial negative
(   37)       !!      concentrations due to isopycnal scheme
(   38)       !!
(   39)       !! ** Method  : Set negative concentrations to zero
(   40)       !!              compute the corresponding mass added to the tracers
(   41)       !!              and remove it when possible 
(   42)       !!
(   43)       !! History :
(   44)       !!   8.2  !  01-01  (O. Aumont & E. Kestenare)  Original code
(   45)       !!   9.0  !  04-03  (C. Ethe)  free form F90
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 2

(   47)       !! * Arguments
(   48)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   49)       
(   50)       !! * Local declarations
(   51)       INTEGER ::  ji, jj, jk, jn             ! dummy loop indices
(   52)       !!----------------------------------------------------------------------
(   53) 
(   54)       IF( kt == nittrc000 ) THEN
(   55)          IF(lwp) WRITE(numout,*)
(   56)          IF(lwp) WRITE(numout,*) 'DB: trc_rad : ROUTINE DISABLED -- find an old one if you need it'
(   57)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   58)       ENDIF
(   59)       
(   60)    END SUBROUTINE trc_rad
(   61) 
(   62) #else
(   63)    !!----------------------------------------------------------------------
(   64)    !!   Dummy module :                      NO passive tracer
(   65)    !!----------------------------------------------------------------------
(   66) CONTAINS
(   67)    SUBROUTINE trc_rad (kt )              ! Empty routine
(   68)       INTEGER, INTENT(in) :: kt
(   69) !      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
(   70)    END SUBROUTINE trc_rad
(   71) #endif
(   72)    
(   73)    !!======================================================================
(   74) END MODULE trcrad




































PGF90 (Version     10.5)          09/23/2010  09:49:22      page 3

(    1) # 1 "trcrad.F90"
(    1) !!DB: 2009.09.09 -- deleted old DBG code 
(    2) !!and disable the routine. If you want it then dig up an older version
(    3) MODULE trcrad
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  trcrad  ***
(    6)    !! Ocean passive tracers:  correction of negative concentrations
(    7)    !!======================================================================
(    8) # 63
(   63)    !!----------------------------------------------------------------------
(   64)    !!   Dummy module :                      NO passive tracer
(   65)    !!----------------------------------------------------------------------
(   66) CONTAINS
(   67)    SUBROUTINE trc_rad (kt )              ! Empty routine
(   68)       INTEGER, INTENT(in) :: kt
(   69) !      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
(   70)    END SUBROUTINE trc_rad
(   71) # 72
(   72)    
(   73)    !!======================================================================
(   74) END MODULE trcrad
