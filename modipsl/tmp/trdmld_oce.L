


PGF90 (Version     10.5)          09/23/2010  09:49:08      page 1

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

Filename: trdmld_oce.F90

(    1) MODULE trdmld_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdmld_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmld_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce         ! ocean parameters
(   14) 
(   15)    IMPLICIT NONE
(   16)    PUBLIC
(   17) 
(   18)    INTEGER, PARAMETER ::            &  !: mixed layer trends index
(   19)       jpmldxad = 1,   &  !: zonal advection
(   20)       jpmldyad = 2,   &  !: meridionnal advection
(   21)       jpmldzad = 3,   &  !: vertical advection
(   22)       jpmldldf = 4,   &  !: lateral diffusion (horiz. component+Beckman)
(   23)       jpmldfor = 5,   &  !: forcing 
(   24)       jpmldevd = 6,   &  !: entrainment due to vertical diffusion (TKE)
(   25)       jpmldzdf = 7,   &  !: explicit vertical part if isopycnal diffusion
(   26)       jpmldxei = 8,   &  !: eddy induced zonal advection
(   27)       jpmldyei = 9,   &  !: eddy induced meridional advection
(   28)       jpmldzei =10       !: eddy induced vertical advection
(   29) 
(   30) #if   defined  key_trdmld   ||   defined key_esopa
(   31)    !!----------------------------------------------------------------------
(   32)    !!   'key_trdmld'                         mixed layer trends diagnostics
(   33)    !!----------------------------------------------------------------------
(   34) 
(   35)    !! Trends diagnostics parameters
(   36)    !!---------------------------------------------------------------------
(   37)    INTEGER, PARAMETER ::            &  !:
(   38) # if defined key_traldf_eiv
(   39)       jpltrd = 10,  &  !: number of mixed-layer trends arrays
(   40)       jpktrd = jpk     !: max level for mixed-layer trends diag.
(   41) # else
(   42)       jpltrd = 7,   &  !: number of mixed-layer trends arrays
(   43)       jpktrd = jpk     !: max level for mixed-layer trends diag.
(   44) # endif
(   45) 
(   46) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:08      page 2

(   47)   !!======================================================================
(   48) END MODULE trdmld_oce






























































PGF90 (Version     10.5)          09/23/2010  09:49:08      page 3

(    1) # 1 "trdmld_oce.F90"
(    1) MODULE trdmld_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdmld_oce  ***
(    4)    !! Ocean trends :   set tracer and momentum trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdmld_oce.F90,v 1.2 2005/03/27 18:35:23 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce         ! ocean parameters
(   14) 
(   15)    IMPLICIT NONE
(   16)    PUBLIC
(   17) 
(   18)    INTEGER, PARAMETER ::            &  !: mixed layer trends index
(   19)       jpmldxad = 1,   &  !: zonal advection
(   20)       jpmldyad = 2,   &  !: meridionnal advection
(   21)       jpmldzad = 3,   &  !: vertical advection
(   22)       jpmldldf = 4,   &  !: lateral diffusion (horiz. component+Beckman)
(   23)       jpmldfor = 5,   &  !: forcing 
(   24)       jpmldevd = 6,   &  !: entrainment due to vertical diffusion (TKE)
(   25)       jpmldzdf = 7,   &  !: explicit vertical part if isopycnal diffusion
(   26)       jpmldxei = 8,   &  !: eddy induced zonal advection
(   27)       jpmldyei = 9,   &  !: eddy induced meridional advection
(   28)       jpmldzei =10       !: eddy induced vertical advection
(   29) 
(   30) # 47
(   47)   !!======================================================================
(   48) END MODULE trdmld_oce
