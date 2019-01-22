


PGF90 (Version     10.5)          09/23/2010  09:49:49      page 1

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

Filename: floats.F90

(    1) !!DB 2008.05.16
(    2) !!Modified routine for more specific DB particle tracking needs
(    3) !!For older code, see OLD_CODE/ or look elsewhere for an older version
(    4) !!Modifications: 
(    5) !!(1) Make jpnfl = #-floats = something that is determined from the input
(    6) !!file "init_float". This means that arrays are dynamically allocated,
(    7) !!allowing also for zero floats (if init_float does not exist),
(    8) !!and that significant code changes occur, starting in flo_oce ...
(    9) !!(2) Note the re-definition of ln_flork4:
(   10) !!if TRUE (see namelist file) then  CALL flo_RDM(kt) which is a DB-written
(   11) !!Random Displacement Model for particle tracking. Note that at this time
(   12) !!this model is written for constant z-level particle tracking. Also,
(   13) !!if ln_flork4 is FALSE then existing Blanke routine is called.
(   14) !!(3) DB has eliminated most of the original float restart stuff, and the
(   15) !!argo float routine. 
(   16) !!Also, flo4rk.F90 has not been re-coded (it did not work anyways)
(   17) 
(   18) MODULE floats
(   19)    !!======================================================================
(   20)    !!                       ***  MODULE  floats  ***
(   21)    !! Ocean floats : floats
(   22)    !!======================================================================
(   23) #if   defined key_floats   ||   defined key_esopa
(   24)    !!----------------------------------------------------------------------
(   25)    !!   'key_floats'                                     float trajectories
(   26)    !!----------------------------------------------------------------------
(   27)    !!   flo_stp   : float trajectories computation
(   28)    !!   flo_init  : initialization of float trajectories computation
(   29)    !!----------------------------------------------------------------------
(   30)    !! * Modules used
(   31)    USE flo_oce         ! floats variables
(   32)    USE lib_mpp         ! distributed memory computing
(   33)    USE flodom          ! initialisation Module 
(   34)    USE flowri          ! float output                     (flo_wri routine)
(   35)    USE flo4rk          ! Trajectories, Runge Kutta scheme (flo_4rk routine)
(   36)    USE floblk          ! Trajectories, Blanke scheme      (flo_blk routine)
(   37) 
(   38)    IMPLICIT NONE
(   39)    PRIVATE  
(   40) 
(   41)    !! * Routine accessibility
(   42)    PUBLIC flo_stp    ! routine called by step.F90
(   43)    !!----------------------------------------------------------------------
(   44)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/floats.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 2

(   47)    !!----------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE flo_stp( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                   ***  ROUTINE flo_stp  ***
(   54)       !!                    
(   55)       !! ** Purpose :   Compute the geographical position (lat., long., depth)
(   56)       !!      of each float at each time step with one of the algorithm.
(   57)       !! 
(   58)       !! ** Method  :   The position of a float is computed with Bruno Blanke 
(   59)       !!        algorithm by default and with a 4th order Runge-Kutta scheme
(   60)       !!        if ln_flork4 =T
(   61)       !!      
(   62)       !! History :
(   63)       !!   8.5  !  02-06  (A. Bozec, G. Madec )  F90: Free form and module
(   64)       !!----------------------------------------------------------------------
(   65)       !! * arguments
(   66)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(   67)       !!----------------------------------------------------------------------
(   68) 
(   69)       IF( kt == nit000 ) THEN
(   70)          IF(lwp) WRITE(numout,*)
(   71)          IF(lwp) WRITE(numout,*) 'flo_stp : call floats routine '
(   72)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   73) 
(   74)          CALL flo_init           ! read the namelist of floats             
(   75) 
(   76)          CALL flo_dom            ! compute/read initial position of floats
(   77) 
(   78)          ! Initialisation of wb for computation of floats trajectories at the first time step
(   79)          wb(:,:,:) = wn(:,:,:)
(   80)       ENDIF
(   81) 
(   82)       if(jpnfl == 0) return 
(   83) 
(   84) 
(   85) !!DB 2008.03.19
(   86)       if( ln_flork4 ) THEN
(   87)          CALL flo_RDM( kt )       ! Trajectories using DB RDM 
(   88)       else
(   89)          CALL flo_blk( kt )        ! Trajectories using Blanke' algorithme
(   90)       endif
(   91) 
(   92)       IF( lk_mpp )   CALL mppsync   ! synchronization of all the processor
(   93) 
(   94) 
(   95)       ! Writing and restart      
(   96)       
(   97)       ! trajectories file 
(   98)       IF( kt == nit000 .OR. MOD( kt, nwritefl ) == 0 )   CALL flo_wri( kt )
(   99)       ! restart file 
(  100)       IF( kt == nitend .OR. MOD( kt, nstockfl ) == 0 )   CALL flo_wri( kt )
(  101) 
(  102)       ! Save the old vertical velocity field
(  103)       wb(:,:,:) = wn(:,:,:)
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:49      page 3

(  105)    END SUBROUTINE flo_stp
(  106) 
(  107) 
(  108)    SUBROUTINE flo_init
(  109)       !!----------------------------------------------------------------
(  110)       !!                 ***  ROUTINE flo_init  ***
(  111)       !!                   
(  112)       !! ** Purpose :   Read the namelist of floats
(  113)       !!      
(  114)       !! History :
(  115)       !!   8.0  !         (CLIPPER)   original Code
(  116)       !!   8.5  !  02-06  (A. Bozec)  F90, Free form and module
(  117)       !!----------------------------------------------------------------------
(  118)       !! * Modules used
(  119) !!DB: 2009.09.02
(  120) !      USE ioipsl
(  121) 
(  122)       !! * Local declarations
(  123)       NAMELIST/namflo/ ln_rstflo, nwritefl, nstockfl, ln_argo, ln_flork4
(  124)       !!---------------------------------------------------------------------
(  125)       ! Namelist namflo : floats
(  126)       
(  127)       ! default values
(  128)       ln_rstflo  = .FALSE.
(  129)       nwritefl  = 150
(  130)       nstockfl  = 450
(  131)       
(  132)       ! lecture of namflo
(  133)       REWIND( numnam )
(  134)       READ  ( numnam, namflo )
(  135) 
(  136)       IF(lwp) THEN
(  137)          WRITE(numout,*) ' '
(  138)          WRITE(numout,*) '         Namelist floats :'
(  139)          WRITE(numout,*) '            restart                          ln_rstflo = ', ln_rstflo
(  140)          WRITE(numout,*) '            frequency of float output file   nwritefl  = ', nwritefl
(  141)          WRITE(numout,*) '            frequency of float restart file  nstockfl  = ', nstockfl
(  142)          WRITE(numout,*) ' '
(  143)       ENDIF
(  144) 
(  145)    END SUBROUTINE flo_init
(  146) 
(  147) #  else
(  148)    !!----------------------------------------------------------------------
(  149)    !!   Default option :                                       Empty module
(  150)    !!----------------------------------------------------------------------
(  151)    USE in_out_manager
(  152) CONTAINS
(  153)    SUBROUTINE flo_stp( kt )          ! Empty routine
(  154)       if(lwp) WRITE(numout,*) 'flo_stp: You should not have seen this print! error?', kt
(  155)    END SUBROUTINE flo_stp
(  156) #endif
(  157) 
(  158)    !!======================================================================
(  159)  END MODULE floats









PGF90 (Version     10.5)          09/23/2010  09:49:49      page 4

(    1) # 1 "floats.F90"
(    1) !!DB 2008.05.16
(    2) !!Modified routine for more specific DB particle tracking needs
(    3) !!For older code, see OLD_CODE/ or look elsewhere for an older version
(    4) !!Modifications: 
(    5) !!(1) Make jpnfl = #-floats = something that is determined from the input
(    6) !!file "init_float". This means that arrays are dynamically allocated,
(    7) !!allowing also for zero floats (if init_float does not exist),
(    8) !!and that significant code changes occur, starting in flo_oce ...
(    9) !!(2) Note the re-definition of ln_flork4:
(   10) !!if TRUE (see namelist file) then  CALL flo_RDM(kt) which is a DB-written
(   11) !!Random Displacement Model for particle tracking. Note that at this time
(   12) !!this model is written for constant z-level particle tracking. Also,
(   13) !!if ln_flork4 is FALSE then existing Blanke routine is called.
(   14) !!(3) DB has eliminated most of the original float restart stuff, and the
(   15) !!argo float routine. 
(   16) !!Also, flo4rk.F90 has not been re-coded (it did not work anyways)
(   17) 
(   18) MODULE floats
(   19)    !!======================================================================
(   20)    !!                       ***  MODULE  floats  ***
(   21)    !! Ocean floats : floats
(   22)    !!======================================================================
(   23) # 148
(  148)    !!----------------------------------------------------------------------
(  149)    !!   Default option :                                       Empty module
(  150)    !!----------------------------------------------------------------------
(  151)    USE in_out_manager
(  152) CONTAINS
(  153)    SUBROUTINE flo_stp( kt )          ! Empty routine
(  154)       if(lwp) WRITE(numout,*) 'flo_stp: You should not have seen this print! error?', kt
(  155)    END SUBROUTINE flo_stp
(  156) 
(  157) # 158
(  158)    !!======================================================================
(  159)  END MODULE floats
