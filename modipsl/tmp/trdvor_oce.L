


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

Filename: trdvor_oce.F90

(    1) MODULE trdvor_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdvor_oce  ***
(    4)    !! Ocean trends :   set vorticity trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdvor_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce      ! ocean parameters
(   14) 
(   15)    IMPLICIT NONE
(   16)    PUBLIC
(   17) 
(   18)    INTEGER,PARAMETER :: jplvor = 11     ! Number of vorticity trend terms
(   19) 
(   20)    INTEGER, PARAMETER ::            &  !: vorticity trends index
(   21)       jpvorprg = 1,   &  !: Pressure Gradient Trend
(   22)       jpvorkeg = 2,   &  !: KE Gradient Trend
(   23)       jpvorrvo = 3,   &  !: Relative Vorticity Trend
(   24)       jpvorpvo = 4,   &  !: Planetary Vorticity Term Trend
(   25)       jpvorldf = 5,   &  !: Horizontal Diffusion Trend
(   26)       jpvorzad = 6,   &  !: Vertical Advection Trend
(   27)       jpvorzdf = 7,   &  !: Vertical Diffusion Trend
(   28)       jpvorspg = 8,   &  !: Surface Pressure Grad. Trend
(   29)       jpvorbev = 9,   &  !: Beta V
(   30)       jpvorswf =10,   &  !: wind stress forcing term
(   31)       jpvorbfr =11       !: bottom friction term
(   32) 
(   33)   !!======================================================================
(   34) END MODULE trdvor_oce


















PGF90 (Version     10.5)          09/23/2010  09:49:08      page 2

(    1) # 1 "trdvor_oce.F90"
(    1) MODULE trdvor_oce
(    2)    !!======================================================================
(    3)    !!                   ***  MODULE trdvor_oce  ***
(    4)    !! Ocean trends :   set vorticity trend variables
(    5)    !!======================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    8)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRD/trdvor_oce.F90,v 1.2 2005/03/27 18:35:24 opalod Exp $ 
(    9)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce      ! ocean parameters
(   14) 
(   15)    IMPLICIT NONE
(   16)    PUBLIC
(   17) 
(   18)    INTEGER,PARAMETER :: jplvor = 11     ! Number of vorticity trend terms
(   19) 
(   20)    INTEGER, PARAMETER ::            &  !: vorticity trends index
(   21)       jpvorprg = 1,   &  !: Pressure Gradient Trend
(   22)       jpvorkeg = 2,   &  !: KE Gradient Trend
(   23)       jpvorrvo = 3,   &  !: Relative Vorticity Trend
(   24)       jpvorpvo = 4,   &  !: Planetary Vorticity Term Trend
(   25)       jpvorldf = 5,   &  !: Horizontal Diffusion Trend
(   26)       jpvorzad = 6,   &  !: Vertical Advection Trend
(   27)       jpvorzdf = 7,   &  !: Vertical Diffusion Trend
(   28)       jpvorspg = 8,   &  !: Surface Pressure Grad. Trend
(   29)       jpvorbev = 9,   &  !: Beta V
(   30)       jpvorswf =10,   &  !: wind stress forcing term
(   31)       jpvorbfr =11       !: bottom friction term
(   32) 
(   33)   !!======================================================================
(   34) END MODULE trdvor_oce
