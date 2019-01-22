


PGF90 (Version     10.5)          09/23/2010  09:49:02      page 1

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

Filename: ocfzpt.F90

(    1) MODULE ocfzpt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ocfzpt  ***
(    4)    !! Ocean active tracers:  freezing point computation and freezing area
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !! * Modules used
(    9)    USE oce             ! ocean dynamics and tracers 
(   10)    USE dom_oce         ! ocean space and time domain 
(   11) 
(   12)    IMPLICIT NONE
(   13)    PRIVATE
(   14) 
(   15)    !! * Routine accessibility
(   16)    PUBLIC oc_fz_pt        ! called by opa.F90 and step.F90
(   17) 
(   18)    !! * Shared module variables   
(   19)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   20)       freeze, freezn,  &  !: after and now ice mask (0 or 1)
(   21)       fzptb, fzptn        !: before and now freezing point
(   22)    !!----------------------------------------------------------------------
(   23)    !! OPA 9.0 , LOCEAN-IPSL (2005) 
(   24)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ocfzpt.F90,v 1.3 2005/03/27 18:34:48 opalod Exp $ 
(   25)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE oc_fz_pt
(   31)       !!----------------------------------------------------------------------
(   32)       !!                     ***  ROUTINE oc_fz_pt  ***
(   33)       !!
(   34)       !! ** Purpose : - Calculate ocean surface freezing temperature
(   35)       !!              - Calculate related boolean for now ice distribution 
(   36)       !!
(   37)       !! ** Method  :   Caution, freezing point only for jackett & McDougall eos
(   38)       !!
(   39)       !! ** Action  : - fzptn  : now freezing temperature at ocean surface
(   40)       !!              - fzptb  : before freezing temperature at ocean surface
(   41)       !!              - freezn : boolean indicating freezing conditions at the 
(   42)       !!                ocean surface at time step "now"
(   43)       !!
(   44)       !! History :
(   45)       !!        !  94-08  (M.-A. Filiberti)  Original code
(   46)       !!   8.5  !  02-08  (G. Madec)  F90: Free form






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 2

(   47)       !!----------------------------------------------------------------------
(   48)       !! * Local declarations
(   49)       INTEGER ::  ji, jj      ! dummy loop indices
(   50)       !!----------------------------------------------------------------------      
(   51) 
(   52) !CDIR NOVERRCHK
(   53)       DO jj = 1, jpj
(   54) !CDIR NOVERRCHK
(   55)          DO ji = 1, jpi
(   56)             fzptb(ji,jj) = fzptn(ji,jj)        ! swap freezing point array
(   57) 
(   58)             !                                  ! ocean local surface freezing temperature
(   59)             !                                  ! sn >= 0 : this is verified in stpctl at each time step
(   60)             fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
(   61)                                       - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
(   62)             !!                        - 7.53e-4 * pressure
(   63) 
(   64)             !                                  ! Define boolean related to freezing conditions
(   65)             freezn(ji,jj) = tmask(ji,jj,1) * MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) )  )
(   66)          END DO
(   67)       END DO
(   68) 
(   69)    END SUBROUTINE oc_fz_pt
(   70) 
(   71)    !!======================================================================
(   72) END MODULE ocfzpt






































PGF90 (Version     10.5)          09/23/2010  09:49:02      page 3

(    1) # 1 "ocfzpt.F90"
(    1) MODULE ocfzpt
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  ocfzpt  ***
(    4)    !! Ocean active tracers:  freezing point computation and freezing area
(    5)    !!======================================================================
(    6)    
(    7)    !!----------------------------------------------------------------------
(    8)    !! * Modules used
(    9)    USE oce             ! ocean dynamics and tracers 
(   10)    USE dom_oce         ! ocean space and time domain 
(   11) 
(   12)    IMPLICIT NONE
(   13)    PRIVATE
(   14) 
(   15)    !! * Routine accessibility
(   16)    PUBLIC oc_fz_pt        ! called by opa.F90 and step.F90
(   17) 
(   18)    !! * Shared module variables   
(   19)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   20)       freeze, freezn,  &  !: after and now ice mask (0 or 1)
(   21)       fzptb, fzptn        !: before and now freezing point
(   22)    !!----------------------------------------------------------------------
(   23)    !! OPA 9.0 , LOCEAN-IPSL (2005) 
(   24)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ocfzpt.F90,v 1.3 2005/03/27 18:34:48 opalod Exp $ 
(   25)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   26)    !!----------------------------------------------------------------------
(   27) 
(   28) CONTAINS
(   29) 
(   30)    SUBROUTINE oc_fz_pt
(   31)       !!----------------------------------------------------------------------
(   32)       !!                     ***  ROUTINE oc_fz_pt  ***
(   33)       !!
(   34)       !! ** Purpose : - Calculate ocean surface freezing temperature
(   35)       !!              - Calculate related boolean for now ice distribution 
(   36)       !!
(   37)       !! ** Method  :   Caution, freezing point only for jackett & McDougall eos
(   38)       !!
(   39)       !! ** Action  : - fzptn  : now freezing temperature at ocean surface
(   40)       !!              - fzptb  : before freezing temperature at ocean surface
(   41)       !!              - freezn : boolean indicating freezing conditions at the 
(   42)       !!                ocean surface at time step "now"
(   43)       !!
(   44)       !! History :
(   45)       !!        !  94-08  (M.-A. Filiberti)  Original code
(   46)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   47)       !!----------------------------------------------------------------------
(   48)       !! * Local declarations
(   49)       INTEGER ::  ji, jj      ! dummy loop indices
(   50)       !!----------------------------------------------------------------------      
(   51) 
(   52) !CDIR NOVERRCHK
(   53)       DO jj = 1, jpj
(   54) !CDIR NOVERRCHK
(   55)          DO ji = 1, jpi
(   56)             fzptb(ji,jj) = fzptn(ji,jj)        ! swap freezing point array
(   57) 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 4

(   58)             !                                  ! ocean local surface freezing temperature
(   59)             !                                  ! sn >= 0 : this is verified in stpctl at each time step
(   60)             fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
(   61)                                       - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
(   62)             !!                        - 7.53e-4 * pressure
(   63) 
(   64)             !                                  ! Define boolean related to freezing conditions
(   65)             freezn(ji,jj) = tmask(ji,jj,1) * MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) )  )
(   66)          END DO
(   67)       END DO
(   68) 
(   69)    END SUBROUTINE oc_fz_pt
(   70) 
(   71)    !!======================================================================
(   72) END MODULE ocfzpt
