


PGF90 (Version     10.5)          09/23/2010  09:49:55      page 1

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

Filename: limtab.F90

(    1) MODULE limtab
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limtab   ***
(    4)    !!              transform 1D (2D) array to a 2D (1D) table
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tab_2d_1d  : 2-D to 1-D
(    9)    !!   tab_1d_2d  : 1-D to 2-D
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_kind
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    !! * Routine accessibility
(   18)    PUBLIC tab_2d_1d  ! called by lim_ther
(   19)    PUBLIC tab_1d_2d  ! called by lim_ther
(   20) 
(   21)    !!----------------------------------------------------------------------
(   22)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtab.F90,v 1.2 2005/03/27 18:34:42 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26) CONTAINS
(   27) 
(   28)    SUBROUTINE tab_2d_1d ( ndim1d, tab1d, tab2d, ndim2d_x, ndim2d_y, tab_ind )
(   29) 
(   30)       INTEGER, INTENT(in) :: &
(   31)          ndim1d, ndim2d_x, ndim2d_y
(   32) 
(   33)       REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT(in) ::  &
(   34)          tab2d
(   35) 
(   36)       INTEGER, DIMENSION ( ndim1d), INTENT ( in) :: &
(   37)          tab_ind
(   38) 
(   39)       REAL(wp), DIMENSION(ndim1d), INTENT ( out) ::  & 
(   40)          tab1d
(   41) 
(   42)       INTEGER ::  &
(   43)          jn , jid, jjd
(   44)         
(   45)       DO jn = 1, ndim1d
(   46)          jid        = MOD( tab_ind(jn) - 1, ndim2d_x ) + 1






PGF90 (Version     10.5)          09/23/2010  09:49:55      page 2

(   47)          jjd        = ( tab_ind(jn) - 1 ) / ndim2d_x + 1
(   48)          tab1d( jn) = tab2d( jid, jjd)
(   49)       END DO 
(   50) 
(   51)    END SUBROUTINE tab_2d_1d
(   52) 
(   53) 
(   54)    SUBROUTINE tab_1d_2d ( ndim1d, tab2d, tab_ind, tab1d, ndim2d_x, ndim2d_y )
(   55) 
(   56)       INTEGER, INTENT ( in) :: &
(   57)          ndim1d, ndim2d_x, ndim2d_y
(   58) 
(   59)       INTEGER, DIMENSION (ndim1d) , INTENT (in) :: &
(   60)          tab_ind
(   61) 
(   62)       REAL(wp), DIMENSION(ndim1d), INTENT (in) ::  &
(   63)          tab1d  
(   64) 
(   65)       REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT ( out) :: &
(   66)          tab2d
(   67) 
(   68)       INTEGER :: &
(   69)          jn, jid, jjd
(   70) 
(   71)       DO jn = 1, ndim1d
(   72)          jid             = MOD( tab_ind(jn) - 1, ndim2d_x) + 1
(   73)          jjd             =    ( tab_ind(jn) - 1 ) / ndim2d_x  + 1
(   74)          tab2d(jid, jjd) = tab1d( jn)
(   75)       END DO
(   76) 
(   77)    END SUBROUTINE tab_1d_2d
(   78) 
(   79) #endif
(   80) END MODULE limtab






























PGF90 (Version     10.5)          09/23/2010  09:49:55      page 3

(    1) # 1 "limtab.F90"
(    1) MODULE limtab
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limtab   ***
(    4)    !!              transform 1D (2D) array to a 2D (1D) table
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tab_2d_1d  : 2-D to 1-D
(    9)    !!   tab_1d_2d  : 1-D to 2-D
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE par_kind
(   13) 
(   14)    IMPLICIT NONE
(   15)    PRIVATE
(   16) 
(   17)    !! * Routine accessibility
(   18)    PUBLIC tab_2d_1d  ! called by lim_ther
(   19)    PUBLIC tab_1d_2d  ! called by lim_ther
(   20) 
(   21)    !!----------------------------------------------------------------------
(   22)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   23)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtab.F90,v 1.2 2005/03/27 18:34:42 opalod Exp $ 
(   24)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   25)    !!----------------------------------------------------------------------
(   26) CONTAINS
(   27) 
(   28)    SUBROUTINE tab_2d_1d ( ndim1d, tab1d, tab2d, ndim2d_x, ndim2d_y, tab_ind )
(   29) 
(   30)       INTEGER, INTENT(in) :: &
(   31)          ndim1d, ndim2d_x, ndim2d_y
(   32) 
(   33)       REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT(in) ::  &
(   34)          tab2d
(   35) 
(   36)       INTEGER, DIMENSION ( ndim1d), INTENT ( in) :: &
(   37)          tab_ind
(   38) 
(   39)       REAL(wp), DIMENSION(ndim1d), INTENT ( out) ::  & 
(   40)          tab1d
(   41) 
(   42)       INTEGER ::  &
(   43)          jn , jid, jjd
(   44)         
(   45)       DO jn = 1, ndim1d
(   46)          jid        = MOD( tab_ind(jn) - 1, ndim2d_x ) + 1
(   47)          jjd        = ( tab_ind(jn) - 1 ) / ndim2d_x + 1
(   48)          tab1d( jn) = tab2d( jid, jjd)
(   49)       END DO 
(   50) 
(   51)    END SUBROUTINE tab_2d_1d
(   52) 
(   53) 
(   54)    SUBROUTINE tab_1d_2d ( ndim1d, tab2d, tab_ind, tab1d, ndim2d_x, ndim2d_y )
(   55) 
(   56)       INTEGER, INTENT ( in) :: &
(   57)          ndim1d, ndim2d_x, ndim2d_y






PGF90 (Version     10.5)          09/23/2010  09:49:55      page 4

(   58) 
(   59)       INTEGER, DIMENSION (ndim1d) , INTENT (in) :: &
(   60)          tab_ind
(   61) 
(   62)       REAL(wp), DIMENSION(ndim1d), INTENT (in) ::  &
(   63)          tab1d  
(   64) 
(   65)       REAL(wp), DIMENSION (ndim2d_x, ndim2d_y), INTENT ( out) :: &
(   66)          tab2d
(   67) 
(   68)       INTEGER :: &
(   69)          jn, jid, jjd
(   70) 
(   71)       DO jn = 1, ndim1d
(   72)          jid             = MOD( tab_ind(jn) - 1, ndim2d_x) + 1
(   73)          jjd             =    ( tab_ind(jn) - 1 ) / ndim2d_x  + 1
(   74)          tab2d(jid, jjd) = tab1d( jn)
(   75)       END DO
(   76) 
(   77)    END SUBROUTINE tab_1d_2d
(   78) 
(   79) # 80
(   80) END MODULE limtab
