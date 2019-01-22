


PGF90 (Version     10.5)          09/23/2010  09:49:27      page 1

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

Filename: diafwb.F90

(    1) !!DB: 2009.08.31 -- Removed ORCA-code
(    2) MODULE diafwb
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  diafwb  ***
(    5)    !! Ocean diagnostics: freshwater budget
(    6)    !!======================================================================
(    7) !#if ( defined key_orca_r2 || defined  key_orca_r4 ) && ! defined key_dynspg_rl && ! defined key_coupled
(    8) 
(    9)    LOGICAL, PUBLIC, PARAMETER ::   lk_diafwb = .FALSE.    !: fresh water budget flag
(   10) CONTAINS
(   11)    SUBROUTINE dia_fwb( kt )        ! Empty routine
(   12)      USE in_out_manager
(   13)       
(   14)      if(lwp .AND. kt==nit000) write(numout,*) 'DB: dia_fwb: disabled '
(   15) 
(   16)    END SUBROUTINE dia_fwb
(   17) END MODULE diafwb



































PGF90 (Version     10.5)          09/23/2010  09:49:27      page 2

(    1) # 1 "diafwb.F90"
(    1) !!DB: 2009.08.31 -- Removed ORCA-code
(    2) MODULE diafwb
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE  diafwb  ***
(    5)    !! Ocean diagnostics: freshwater budget
(    6)    !!======================================================================
(    7) !#if ( defined key_orca_r2 || defined  key_orca_r4 ) && ! defined key_dynspg_rl && ! defined key_coupled
(    8) 
(    9)    LOGICAL, PUBLIC, PARAMETER ::   lk_diafwb = .FALSE.    !: fresh water budget flag
(   10) CONTAINS
(   11)    SUBROUTINE dia_fwb( kt )        ! Empty routine
(   12)      USE in_out_manager
(   13)       
(   14)      if(lwp .AND. kt==nit000) write(numout,*) 'DB: dia_fwb: disabled '
(   15) 
(   16)    END SUBROUTINE dia_fwb
(   17) END MODULE diafwb
