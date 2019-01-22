


PGF90 (Version     10.5)          09/23/2010  09:49:15      page 1

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

Filename: trc_oce.F90

(    1) MODULE trc_oce
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  trc_oce  ***
(    4)    !! Ocean passive tracer  :  share SMS/Ocean variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9) #if defined key_passivetrc && defined key_trc_pisces
(   10)    !!----------------------------------------------------------------------
(   11)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/trc_oce.F90,v 1.3 2005/03/27 18:34:49 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15)    USE sms , ONLY :  &
(   16)       etot3    =>   etot3   !!:  Biological fluxes for light
(   17)    !! Shared module variables
(   18)    LOGICAL, PUBLIC, PARAMETER ::   lk_qsr_sms = .TRUE. 
(   19) #else
(   20)    !!----------------------------------------------------------------------
(   21)    !! Default option                         No Biological fluxes for light          
(   22)    !!----------------------------------------------------------------------
(   23)    USE par_oce
(   24)    LOGICAL, PUBLIC, PARAMETER ::   lk_qsr_sms = .FALSE. 
(   25)    REAL(wp), PUBLIC , DIMENSION (jpi,jpj,jpk) :: &
(   26)       etot3
(   27) #endif
(   28) 
(   29) END MODULE trc_oce























PGF90 (Version     10.5)          09/23/2010  09:49:15      page 2

(    1) # 1 "trc_oce.F90"
(    1) MODULE trc_oce
(    2)    !!======================================================================
(    3)    !!                      ***  MODULE  trc_oce  ***
(    4)    !! Ocean passive tracer  :  share SMS/Ocean variables
(    5)    !!======================================================================
(    6)    !! History :
(    7)    !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
(    8)    !!----------------------------------------------------------------------
(    9) # 20
(   20)    !!----------------------------------------------------------------------
(   21)    !! Default option                         No Biological fluxes for light          
(   22)    !!----------------------------------------------------------------------
(   23)    USE par_oce
(   24)    LOGICAL, PUBLIC, PARAMETER ::   lk_qsr_sms = .FALSE. 
(   25)    REAL(wp), PUBLIC , DIMENSION (jpi,jpj,jpk) :: &
(   26)       etot3
(   27) 
(   28) # 29
(   29) END MODULE trc_oce
