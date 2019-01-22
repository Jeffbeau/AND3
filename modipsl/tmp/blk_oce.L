


PGF90 (Version     10.5)          09/23/2010  09:49:01      page 1

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

Filename: blk_oce.F90

(    1) MODULE blk_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  blk_oce  ***
(    4)    !! Bulk   :  bulk parameter and  variables defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/blk_oce.F90,v 1.4 2005/09/22 10:58:15 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(   14)    !!----------------------------------------------------------------------
(   15)    !! ' key_flx_bulk_monthly or defined key_flx_bulk_daily             bulk
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce          ! ocean parameters
(   19) 
(   20)    IMPLICIT NONE
(   21) 
(   22)    LOGICAL, PUBLIC ::   l_bulk = .TRUE.   !: 
(   23)    
(   24)    !!----------------------------------------------------------------------
(   25)    !! bulk common variables
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   29)       watm     ,      &  !: precipitation
(   30)       tatm     ,      &  !: atmospheric temperature
(   31)       hatm     ,      &  !: relative humidity
(   32)       vatm     ,      &  !: wind speed
(   33)       catm               !: percent of cloud cover
(   34) 
(   35)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
(   36)       gsst               !: SST mean on nfbulk ocean time step
(   37) 
(   38)    REAL(wp) ::        &
(   39)       yearday  ,      &  !: number of days per year
(   40)       rdtbs2             !: bulk time step divide by 2
(   41) #else
(   42)    !!----------------------------------------------------------------------
(   43)    !!   Default option                                         Empty module
(   44)    !!----------------------------------------------------------------------
(   45)    LOGICAL, PUBLIC ::   l_bulk = .FALSE.  !:
(   46) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:01      page 2

(   47)    
(   48)    INTEGER ::         & !!: namdom : space/time domain (namlist)
(   49)       nfbulk =  5        !: bulk computation frequency 
(   50)    !!----------------------------------------------------------------------
(   51) END MODULE blk_oce



























































PGF90 (Version     10.5)          09/23/2010  09:49:01      page 3

(    1) # 1 "blk_oce.F90"
(    1) MODULE blk_oce
(    2)    !!======================================================================
(    3)    !!                 ***  MODULE  blk_oce  ***
(    4)    !! Bulk   :  bulk parameter and  variables defined in memory 
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   10)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/blk_oce.F90,v 1.4 2005/09/22 10:58:15 opalod Exp $ 
(   11)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   12)    !!----------------------------------------------------------------------
(   13) # 14
(   14)    !!----------------------------------------------------------------------
(   15)    !! ' key_flx_bulk_monthly or defined key_flx_bulk_daily             bulk
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE par_oce          ! ocean parameters
(   19) 
(   20)    IMPLICIT NONE
(   21) 
(   22)    LOGICAL, PUBLIC ::   l_bulk = .TRUE.   !: 
(   23)    
(   24)    !!----------------------------------------------------------------------
(   25)    !! bulk common variables
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   29)       watm     ,      &  !: precipitation
(   30)       tatm     ,      &  !: atmospheric temperature
(   31)       hatm     ,      &  !: relative humidity
(   32)       vatm     ,      &  !: wind speed
(   33)       catm               !: percent of cloud cover
(   34) 
(   35)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
(   36)       gsst               !: SST mean on nfbulk ocean time step
(   37) 
(   38)    REAL(wp) ::        &
(   39)       yearday  ,      &  !: number of days per year
(   40)       rdtbs2             !: bulk time step divide by 2
(   41) # 47
(   47)    
(   48)    INTEGER ::         & !!: namdom : space/time domain (namlist)
(   49)       nfbulk =  5        !: bulk computation frequency 
(   50)    !!----------------------------------------------------------------------
(   51) END MODULE blk_oce
