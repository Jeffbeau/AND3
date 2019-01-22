


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

Filename: trcsbc.F90

(    1) MODULE trcsbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcsbc  ***
(    4)    !! Ocean passive tracers:  surface boundary condition
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_sbc      : update the tracer trend at ocean surface
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce_trc             ! ocean dynamics and active tracers variables
(   12)    USE trc                 ! ocean  passive tracers variables
(   13)    USE prtctl_trc          ! Print control for debbuging
(   14) 
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC trc_sbc              ! routine called by step.F90
(   21) 
(   22)    !! * Substitutions
(   23) #  include "passivetrc_substitute.h90"
(   24)    !!----------------------------------------------------------------------
(   25)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   26)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcsbc.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
(   27)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   28)    !!----------------------------------------------------------------------
(   29) 
(   30) CONTAINS
(   31) 
(   32)    SUBROUTINE trc_sbc ( kt )
(   33)       !!----------------------------------------------------------------------
(   34)       !!                  ***  ROUTINE trc_sbc  ***
(   35)       !!                   
(   36)       !! ** Purpose :   Compute the tracer surface boundary condition trend of
(   37)       !!      (concentration/dilution effect) and add it to the general 
(   38)       !!       trend of tracer equations.
(   39)       !!
(   40)       !! ** Method :
(   41)       !!      * concentration/dilution effect:
(   42)       !!            The surface freshwater flux modify the ocean volume
(   43)       !!         and thus the concentration of a tracer as :
(   44)       !!            tra = tra + emp * trn / e3t   for k=1
(   45)       !!         where emp, the surface freshwater budget (evaporation minus
(   46)       !!         precipitation minus runoff) given in kg/m2/s is divided






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 2

(   47)       !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
(   48)       !!
(   49)       !! ** Action  : - Update the 1st level of tra with the trend associated
(   50)       !!                with the tracer surface boundary condition 
(   51)       !!
(   52)       !! History :
(   53)       !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
(   54)       !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
(   55)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   56)       !!   9.0  !  04-03  (C. Ethe)  adapted for passive tracers
(   57)       !!----------------------------------------------------------------------
(   58)       !! * Arguments
(   59)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(   60) 
(   61)       !! * Local declarations
(   62)       INTEGER  ::   ji, jj, jn           ! dummy loop indices
(   63)       REAL(wp) ::   ztra, zsrau, zse3t   ! temporary scalars
(   64)       CHARACTER (len=22) :: charout
(   65)       !!----------------------------------------------------------------------
(   66) 
(   67)       IF( kt == nittrc000 ) THEN
(   68)          IF(lwp) WRITE(numout,*)
(   69)          IF(lwp) WRITE(numout,*) 'DB: trc_sbc : Potentially dangerous Passive tracer sfce BC routine -- disabled'
(   70)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(   71)       ENDIF
(   72) 
(   73) !!DB 2009.09.09
(   74)       return 
(   75) 
(   76) 
(   77) 
(   78)       ! 0. initialization
(   79)       zsrau = 1. / rauw
(   80) #if ! defined key_s_coord
(   81)       zse3t = 1. / fse3t(1,1,1)
(   82) #endif
(   83) 
(   84)       DO jn = 1, jptra
(   85)          ! 1. Concentration dillution effect on tra
(   86)          DO jj = 2, jpj
(   87)             DO ji = fs_2, fs_jpim1   ! vector opt.
(   88) #if defined key_s_coord
(   89)                zse3t = 1. / fse3t(ji,jj,1)
(   90) #endif
(   91)                ! concent./dilut. effect
(   92)                ztra = emps(ji,jj) * zsrau * trn(ji,jj,1,jn) * zse3t * tmask(ji,jj,1)
(   93)                
(   94)                ! add the trend to the general tracer trend
(   95)                tra(ji,jj,1,jn) = tra(ji,jj,1,jn) + ztra
(   96)             END DO
(   97)          END DO
(   98)          
(   99)       END DO
(  100) 
(  101)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  102)          WRITE(charout, FMT="('sbc')")
(  103)          CALL prt_ctl_trc_info(charout)
(  104)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 3

(  105)       ENDIF
(  106) 
(  107)    END SUBROUTINE trc_sbc
(  108) 
(  109) #else
(  110)    !!----------------------------------------------------------------------
(  111)    !!   Dummy module :                      NO passive tracer
(  112)    !!----------------------------------------------------------------------
(  113) CONTAINS
(  114)    SUBROUTINE trc_sbc (kt)              ! Empty routine
(  115)       INTEGER, INTENT(in) :: kt
(  116) !      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
(  117)    END SUBROUTINE trc_sbc
(  118) #endif
(  119)    
(  120)    !!======================================================================
(  121) END MODULE trcsbc















































PGF90 (Version     10.5)          09/23/2010  09:49:22      page 4

(    1) # 1 "trcsbc.F90"
(    1) MODULE trcsbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcsbc  ***
(    4)    !! Ocean passive tracers:  surface boundary condition
(    5)    !!==============================================================================
(    6) # 110
(  110)    !!----------------------------------------------------------------------
(  111)    !!   Dummy module :                      NO passive tracer
(  112)    !!----------------------------------------------------------------------
(  113) CONTAINS
(  114)    SUBROUTINE trc_sbc (kt)              ! Empty routine
(  115)       INTEGER, INTENT(in) :: kt
(  116) !      WRITE(*,*) 'trc_sbc: You should not have seen this print! error?', kt
(  117)    END SUBROUTINE trc_sbc
(  118) # 119
(  119)    
(  120)    !!======================================================================
(  121) END MODULE trcsbc
