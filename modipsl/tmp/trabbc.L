


PGF90 (Version     10.5)          09/23/2010  09:50:13      page 1

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

Filename: trabbc.F90

(    1) MODULE trabbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trabbc  ***
(    4)    !! Ocean active tracers:  bottom boundary condition
(    5)    !!==============================================================================
(    6) #if   defined key_trabbc   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_trabbc'                                  geothermal heat flux
(    9)    !!----------------------------------------------------------------------
(   10)    !!   tra_bbc      : update the tracer trend at ocean bottom 
(   11)    !!   tra_bbc_init : initialization of geothermal heat flux trend
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and active tracers
(   15)    USE dom_oce         ! ocean space and time domain
(   16)    USE phycst          ! physical constants
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE prtctl          ! Print control
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Accessibility
(   24)    PUBLIC tra_bbc          ! routine called by step.F90
(   25) 
(   26)    !! to be transfert in the namelist ???!   
(   27)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbc = .TRUE.   !: bbc flag
(   28) 
(   29)    !! * Module variables
(   30)    INTEGER ::                       & !!! ** bbc namelist (nambbc) **
(   31)       ngeo_flux = 1                    ! Geothermal flux (0:no flux, 1:constant flux,
(   32)       !                                !                  2:read in file )
(   33)    REAL(wp) ::                      & !!! ** bbc namlist **
(   34)       ngeo_flux_const = 86.4e-3        ! Constant value of geothermal heat flux
(   35) 
(   36)    INTEGER, DIMENSION(jpi,jpj) ::   &
(   37)       nbotlevt                         ! ocean bottom level index at T-pt
(   38)    REAL(wp), DIMENSION(jpi,jpj) ::  &
(   39)       qgh_trd                          ! geothermal heating trend
(   40)  
(   41)    !! * Substitutions
(   42) #  include "domzgr_substitute.h90"
(   43)    !!----------------------------------------------------------------------
(   44)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   45)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trabbc.F90,v 1.8 2006/03/10 10:55:44 opalod Exp $ 
(   46)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 2

(   47)    !!----------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE tra_bbc( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                  ***  ROUTINE tra_bbc  ***
(   54)       !!
(   55)       !! ** Purpose :   Compute the bottom boundary contition on temperature 
(   56)       !!      associated with geothermal heating and add it to the general
(   57)       !!      trend of temperature equations.
(   58)       !!
(   59)       !! ** Method  :   The geothermal heat flux set to its constant value of 
(   60)       !!       86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
(   61)       !!       The temperature trend associated to this heat flux through the
(   62)       !!       ocean bottom can be computed once and is added to the temperature
(   63)       !!       trend juste above the bottom at each time step:
(   64)       !!            ta = ta + Qsf / (rau0 rcp e3T) for k= mbathy -1
(   65)       !!       Where Qsf is the geothermal heat flux.
(   66)       !!
(   67)       !! ** Action  : - update the temperature trends (ta) with the trend of
(   68)       !!                the ocean bottom boundary condition
(   69)       !!
(   70)       !! References :
(   71)       !!      Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
(   72)       !!
(   73)       !! History :
(   74)       !!   8.1  !  99-10  (G. Madec)  original code
(   75)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(   76)       !!----------------------------------------------------------------------
(   77)       !! * Arguments
(   78)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
(   79) 
(   80)       !! * Local declarations
(   81) #if defined key_vectopt_loop   &&   ! defined key_autotasking
(   82)       INTEGER ::   ji                  ! dummy loop indices
(   83) #else
(   84)       INTEGER ::   ji, jj              ! dummy loop indices
(   85) #endif
(   86)       !!----------------------------------------------------------------------
(   87) 
(   88)       ! 0. Initialization
(   89)       IF( kt == nit000 )   CALL tra_bbc_init
(   90) 
(   91)       ! 1. Add the geothermal heat flux trend on temperature
(   92) 
(   93)       SELECT CASE ( ngeo_flux )
(   94) 
(   95)       CASE ( 1:2 )                !  geothermal heat flux
(   96) 
(   97) #if defined key_vectopt_loop   &&   ! defined key_autotasking
(   98)          DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(   99)             ta(ji,1,nbotlevt(ji,1)) = ta(ji,1,nbotlevt(ji,1)) + qgh_trd(ji,1)
(  100)          END DO
(  101) #else
(  102)          DO jj = 2, jpjm1
(  103)             DO ji = 2, jpim1
(  104)                ta(ji,jj,nbotlevt(ji,jj)) = ta(ji,jj,nbotlevt(ji,jj)) + qgh_trd(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 3

(  105)             END DO
(  106)          END DO
(  107) #endif
(  108) 
(  109)       IF(ln_ctl) THEN
(  110)          CALL prt_ctl(tab3d_1=ta, clinfo1=' bbc  - Ta: ', mask1=tmask, clinfo3='tra-ta')
(  111)       ENDIF
(  112) 
(  113)       END SELECT
(  114) 
(  115)    END SUBROUTINE tra_bbc
(  116) 
(  117) 
(  118)    SUBROUTINE tra_bbc_init
(  119)       !!----------------------------------------------------------------------
(  120)       !!                  ***  ROUTINE tra_bbc_init  ***
(  121)       !!
(  122)       !! ** Purpose :   Compute once for all the trend associated with geo-
(  123)       !!      thermal heating that will be applied at each time step at the
(  124)       !!      bottom ocean level
(  125)       !!
(  126)       !! ** Method  :   Read the nambbc namelist and check the parameters.
(  127)       !!      called at the first time step (nit000)
(  128)       !!
(  129)       !! ** Input   : - Namlist nambbc
(  130)       !!              - NetCDF file  : geothermal_heating.nc ( if necessary )
(  131)       !!
(  132)       !! ** Action  : - compute the heat geothermal trend qgh_trd
(  133)       !!              - compute the bottom ocean level nbotlevt
(  134)       !!
(  135)       !! history :
(  136)       !!  8.5  ! 02-11 (A. Bozec) original code
(  137)       !!----------------------------------------------------------------------
(  138)       !! * Modules used
(  139)       USE ioipsl
(  140) 
(  141)       !! * local declarations
(  142)       CHARACTER (len=32) ::   clname
(  143)       INTEGER  ::   ji, jj              ! dummy loop indices
(  144)       INTEGER  ::   inum = 11           ! temporary logical unit
(  145)       INTEGER  ::   itime               ! temporary integers
(  146)       REAL(wp) ::   zdate0, zdt         ! temporary scalars
(  147)       REAL(wp), DIMENSION(1) :: zdept   ! temporary workspace
(  148)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  149)          zlamt, zphit, zdta   ! temporary workspace
(  150) 
(  151)       NAMELIST/nambbc/ngeo_flux, ngeo_flux_const 
(  152)       !!----------------------------------------------------------------------
(  153) 
(  154)       ! Read Namelist nambbc : bottom momentum boundary condition
(  155)       REWIND ( numnam )
(  156)       READ   ( numnam, nambbc )
(  157) 
(  158)       ! Control print
(  159)       IF(lwp) WRITE(numout,*)
(  160)       IF(lwp) WRITE(numout,*) 'tra_bbc : tempearture Bottom Boundary Condition (bbc)'
(  161)       IF(lwp) WRITE(numout,*) '~~~~~~~   Geothermal heatflux'
(  162)       IF(lwp) WRITE(numout,*) '          Namelist nambbc : set bbc parameters'






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 4

(  163)       IF(lwp) WRITE(numout,*)
(  164)       IF(lwp) WRITE(numout,*) '             Geothermal flux           ngeo_flux       = ', ngeo_flux
(  165)       IF(lwp) WRITE(numout,*) '             Constant geothermal flux  ngeo_flux_const = ', ngeo_flux_const
(  166)       IF(lwp) WRITE(numout,*)
(  167) 
(  168)       ! level of the ocean bottom at T-point
(  169) 
(  170)       DO jj = 1, jpj
(  171)          DO ji = 1, jpi
(  172)             nbotlevt(ji,jj) = MAX( mbathy(ji,jj)-1, 1 )
(  173)          END DO
(  174)       END DO
(  175) 
(  176)       ! initialization of geothermal heat flux
(  177) 
(  178)       SELECT CASE ( ngeo_flux )
(  179) 
(  180)       CASE ( 0 )                ! no geothermal heat flux
(  181)          IF(lwp) WRITE(numout,*)
(  182)          IF(lwp) WRITE(numout,*) '             *** no geothermal heat flux'
(  183) 
(  184)       CASE ( 1 )                ! constant flux
(  185)          IF(lwp) WRITE(numout,*) '             *** constant heat flux  =   ', ngeo_flux_const
(  186)          qgh_trd(:,:) = ngeo_flux_const
(  187) 
(  188)       CASE ( 2 )                ! variable geothermal heat flux
(  189)          ! read the geothermal fluxes in mW/m2
(  190)          clname = 'geothermal_heating'
(  191) #if defined key_agrif
(  192)       if ( .NOT. Agrif_Root() ) then
(  193)          clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
(  194)       endif
(  195) #endif 	 
(  196)          itime = 1
(  197)          zlamt(:,:) = 0.
(  198)          zphit(:,:) = 0.
(  199)          IF(lwp) WRITE(numout,*) '             *** variable geothermal heat flux read in ', clname, ' file'
(  200)          CALL restini( clname, jpidta, jpjdta, zlamt, zphit, 1, zdept , 'NONE',   &
(  201)             &          itime, zdate0, zdt, inum, domain_id=nidom )
(  202)          CALL restget( inum, 'heatflow', jpidta, jpjdta, 1, itime, .FALSE., zdta )
(  203)          DO jj = 1, nlcj
(  204)             DO ji = 1, nlci
(  205)               qgh_trd(ji,jj) = zdta(mig(ji),mjg(jj))
(  206)             END DO
(  207)          END DO
(  208) 
(  209)          CALL restclo( inum )
(  210)          qgh_trd(:,:) = qgh_trd(:,:) * 1.e-3 ! conversion in W/m2
(  211) 
(  212)       CASE DEFAULT
(  213)          IF(lwp) WRITE(numout,cform_err)
(  214)          IF(lwp) WRITE(numout,*) '     bad flag value for ngeo_flux = ', ngeo_flux
(  215)          nstop = nstop + 1
(  216) 
(  217)       END SELECT
(  218) 
(  219)       ! geothermal heat flux trend
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 5

(  221)       SELECT CASE ( ngeo_flux )
(  222) 
(  223)       CASE ( 1:2 )                !  geothermal heat flux
(  224) 
(  225) #if defined key_vectopt_loop   &&   ! defined key_autotasking
(  226)          DO ji = 1, jpij   ! vector opt. (forced unrolling)
(  227)             qgh_trd(ji,1) = ro0cpr * qgh_trd(ji,1) / fse3t(ji,1,nbotlevt(ji,1) )
(  228)          END DO
(  229) #else
(  230)          DO jj = 1, jpj
(  231)             DO ji = 1, jpi
(  232)                qgh_trd(ji,jj) = ro0cpr * qgh_trd(ji,jj) / fse3t(ji,jj,nbotlevt(ji,jj))
(  233)             END DO
(  234)          END DO
(  235) #endif
(  236) 
(  237)       END SELECT
(  238) 
(  239)    END SUBROUTINE tra_bbc_init
(  240) 
(  241) #else
(  242)    !!----------------------------------------------------------------------
(  243)    !!   Default option                                         Empty module
(  244)    !!----------------------------------------------------------------------
(  245)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbc = .FALSE.  !: bbc flag
(  246) CONTAINS
(  247)    SUBROUTINE tra_bbc( kt )           ! Empty routine
(  248) !      WRITE(*,*) 'tra_bbc: You should not have seen this print! error?', kt
(  249)    END SUBROUTINE tra_bbc
(  250) #endif
(  251) 
(  252)    !!======================================================================
(  253) END MODULE trabbc































PGF90 (Version     10.5)          09/23/2010  09:50:13      page 6

(    1) # 1 "trabbc.F90"
(    1) MODULE trabbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trabbc  ***
(    4)    !! Ocean active tracers:  bottom boundary condition
(    5)    !!==============================================================================
(    6) # 242
(  242)    !!----------------------------------------------------------------------
(  243)    !!   Default option                                         Empty module
(  244)    !!----------------------------------------------------------------------
(  245)    LOGICAL, PUBLIC, PARAMETER ::   lk_trabbc = .FALSE.  !: bbc flag
(  246) CONTAINS
(  247)    SUBROUTINE tra_bbc( kt )           ! Empty routine
(  248) !      WRITE(*,*) 'tra_bbc: You should not have seen this print! error?', kt
(  249)    END SUBROUTINE tra_bbc
(  250) 
(  251) # 252
(  252)    !!======================================================================
(  253) END MODULE trabbc
