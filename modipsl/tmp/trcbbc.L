


PGF90 (Version     10.5)          09/23/2010  09:50:23      page 1

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

Filename: trcbbc.F90

(    1) MODULE trcbbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcbbc  ***
(    4)    !! Ocean passive tracers:  bottom boundary condition
(    5)    !!==============================================================================
(    6) #if   defined key_passivetrc && defined key_trcbbc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_trcbbc'                                  geothermal heat flux
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_bbc      : update the tracer trend at ocean bottom 
(   11)    !!   trc_bbc_init : initialization of geothermal heat flux trend
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce_trc             ! ocean dynamics and active tracers variables
(   15)    USE trc                 ! ocean passive tracers variables
(   16)    USE prtctl_trc          ! Print control for debbuging
(   17)  
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC trc_bbc          ! routine called by trcstp.F90
(   23) 
(   24)    !! to be transfert in the namelist ???!   
(   25)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbc = .TRUE.   !: bbc flag
(   26) 
(   27)    !! * Module variables
(   28)    INTEGER ::                       & !!! ** bbc namelist (nambbc) **
(   29)       ngeo_trc_flux = 1                    ! Geothermal flux (0:no flux, 1:constant flux,
(   30)       !                                !                  2:read in file )
(   31)    REAL(wp) ::                      & !!! ** bbc namlist **
(   32)       ngeo_trc_flux_const = 86.4e-3        ! Constant value of geothermal heat flux
(   33) 
(   34)    INTEGER, DIMENSION(jpi,jpj) ::   &
(   35)       nbotlevt                         ! ocean bottom level index at T-pt
(   36)    REAL(wp), DIMENSION(jpi,jpj) ::  &
(   37)       qgh_trd                          ! geothermal heating trend
(   38)  
(   39)    !! * Substitutions
(   40) #  include "passivetrc_substitute.h90"
(   41)    !!----------------------------------------------------------------------
(   42)    !!  TOP 1.0 , LOCEAN-IPSL (2005) 
(   43)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcbbc.F90,v 1.9 2005/12/12 14:18:10 opalod Exp $ 
(   44)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   45)    !!----------------------------------------------------------------------
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 2

(   47) CONTAINS
(   48) 
(   49)    SUBROUTINE trc_bbc( kt )
(   50)       !!----------------------------------------------------------------------
(   51)       !!                  ***  ROUTINE trc_bbc  ***
(   52)       !!
(   53)       !! ** Purpose :   Compute the bottom boundary contition on passive tracer 
(   54)       !!      associated with geothermal heating and add it to the general
(   55)       !!      trend of tracers equations.
(   56)       !!
(   57)       !! ** Method  :   The geothermal heat flux set to its constant value of 
(   58)       !!       86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
(   59)       !!       The temperature trend associated to this heat flux through the
(   60)       !!       ocean bottom can be computed once and is added to the temperature
(   61)       !!       trend juste above the bottom at each time step:
(   62)       !!            tra = tra + Qsf / (rau0 rcp e3T) for k= mbathy -1
(   63)       !!       Where Qsf is the geothermal heat flux.
(   64)       !!
(   65)       !! ** Action  : - update the temperature trends tra with the trend of
(   66)       !!                the ocean bottom boundary condition
(   67)       !!
(   68)       !! References :
(   69)       !!      Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
(   70)       !!
(   71)       !! History :
(   72)       !!   8.1  !  99-10  (G. Madec)  original code
(   73)       !!   8.5  !  02-08  (G. Madec)  free form + modules
(   74)       !!   9.0  !  04-03  (C. Ethe)  adpated for passive tracers
(   75)       !!----------------------------------------------------------------------
(   76)       !! * Arguments
(   77)       INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
(   78) 
(   79)       !! * Local declarations
(   80) #if defined key_vectopt_loop   &&   ! defined key_autotasking
(   81)       INTEGER ::   ji, jn                  ! dummy loop indices
(   82) #else
(   83)       INTEGER ::   ji, jj, jn              ! dummy loop indices
(   84) #endif
(   85)       REAL(wp) ::   ztra                ! temporary scalar
(   86)       CHARACTER (len=22) :: charout
(   87)       !!----------------------------------------------------------------------
(   88) 
(   89)       ! 0. Initialization
(   90)       IF( kt == nittrc000 )   CALL trc_bbc_init
(   91) 
(   92)       ! 1. Add the geothermal heat flux trend on temperature
(   93) 
(   94)       SELECT CASE ( ngeo_trc_flux )
(   95) 
(   96)       CASE ( 1:2 )                !  geothermal heat flux
(   97) 
(   98)          DO jn = 1, jptra
(   99) #if defined key_vectopt_loop   &&   ! defined key_autotasking
(  100)             DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
(  101)                tra(ji,1,nbotlevt(ji,1),jn) = tra(ji,1,nbotlevt(ji,1),jn) + qgh_trd(ji,1)
(  102)             END DO
(  103) #else
(  104)             DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 3

(  105)                DO ji = 2, jpim1
(  106)                   tra(ji,jj,nbotlevt(ji,jj),jn) = tra(ji,jj,nbotlevt(ji,jj),jn) + qgh_trd(ji,jj)
(  107)                END DO
(  108)             END DO
(  109) #endif
(  110)          END DO
(  111) 
(  112)          IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  113)             WRITE(charout, FMT="('bbc')")
(  114)             CALL prt_ctl_trc_info(charout)
(  115)             CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  116)          ENDIF
(  117)       END SELECT
(  118) 
(  119)    END SUBROUTINE trc_bbc
(  120) 
(  121) 
(  122)    SUBROUTINE trc_bbc_init
(  123)       !!----------------------------------------------------------------------
(  124)       !!                  ***  ROUTINE trc_bbc_init  ***
(  125)       !!
(  126)       !! ** Purpose :   Compute once for all the trend associated with geo-
(  127)       !!      thermal heating that will be applied at each time step at the
(  128)       !!      bottom ocean level
(  129)       !!
(  130)       !! ** Method  :   Read the namtrabbc namelist and check the parameters.
(  131)       !!      called at the first time step (nittrc000)
(  132)       !!
(  133)       !! ** Input   : - Namlist namtrcbbc
(  134)       !!              - NetCDF file  : passivetrc_geothermal_heating.nc 
(  135)       !!                               ( if necessary )
(  136)       !!
(  137)       !! ** Action  : - compute the heat geothermal trend qgh_trd
(  138)       !!              - compute the bottom ocean level nbotlevt
(  139)       !!
(  140)       !! history :
(  141)       !!  8.5  ! 02-11 (A. Bozec) original code
(  142)       !!----------------------------------------------------------------------
(  143)       !! * Modules used
(  144)       USE ioipsl
(  145) 
(  146)       !! * local declarations
(  147)       CHARACTER (len=32) ::   clname
(  148)       INTEGER  ::   ji, jj              ! dummy loop indices
(  149)       INTEGER  ::   inum = 11           ! temporary logical unit
(  150)       INTEGER  ::   itime               ! temporary integers
(  151)       REAL(wp) ::   zdate0, zdt         ! temporary scalars
(  152)       REAL(wp), DIMENSION(1) :: zdept   ! temporary workspace
(  153)       REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
(  154)          zlamt, zphit, zdta   ! temporary workspace
(  155) 
(  156)       NAMELIST/namtrcbbc/ngeo_trc_flux, ngeo_trc_flux_const 
(  157)       !!----------------------------------------------------------------------
(  158) 
(  159)       ! Read Namelist nambbc : bottom momentum boundary condition
(  160)       REWIND ( numnamtra )
(  161)       READ   ( numnamtra, namtrcbbc )
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 4

(  163)       ! Control print
(  164)       IF(lwp) WRITE(numout,*)
(  165)       IF(lwp) WRITE(numout,*) 'trc_bbc : Passive tracers Bottom Boundary Condition (bbc)'
(  166)       IF(lwp) WRITE(numout,*) '~~~~~~~   Geothermal heatflux'
(  167)       IF(lwp) WRITE(numout,*) '          Namelist namtrcbbc : set bbc parameters'
(  168)       IF(lwp) WRITE(numout,*)
(  169)       IF(lwp) WRITE(numout,*) '             Geothermal flux           ngeo_trc_flux       = ', ngeo_trc_flux
(  170)       IF(lwp) WRITE(numout,*) '             Constant geothermal flux  ngeo_trc_flux_const = ', ngeo_trc_flux_const
(  171)       IF(lwp) WRITE(numout,*)
(  172) 
(  173)       ! level of the ocean bottom at T-point
(  174) 
(  175)       DO jj = 1, jpj
(  176)          DO ji = 1, jpi
(  177)             nbotlevt(ji,jj) = MAX( mbathy(ji,jj)-1, 1 )
(  178)          END DO
(  179)       END DO
(  180) 
(  181)       ! initialization of geothermal heat flux
(  182) 
(  183)       SELECT CASE ( ngeo_trc_flux )
(  184) 
(  185)       CASE ( 0 )                ! no geothermal heat flux
(  186)          IF(lwp) WRITE(numout,*)
(  187)          IF(lwp) WRITE(numout,*) '             *** no geothermal heat flux'
(  188) 
(  189)       CASE ( 1 )                ! constant flux
(  190)          IF(lwp) WRITE(numout,*) '             *** constant heat flux  =   ', ngeo_trc_flux_const
(  191)          qgh_trd(:,:) = ngeo_trc_flux_const
(  192) 
(  193)       CASE ( 2 )                ! variable geothermal heat flux
(  194)          ! read the geothermal fluxes in mW/m2
(  195)          clname = 'passivetrc_geothermal_heating'
(  196)          itime = 1
(  197)          zlamt(:,:) = 0.
(  198)          zphit(:,:) = 0.
(  199)          IF(lwp) WRITE(numout,*) '             *** variable geothermal heat flux read in ', clname, ' file'
(  200)          CALL restini( clname, jpidta, jpjdta, zlamt, zphit, 1, zdept , clname,   &
(  201)                        itime, zdate0, zdt, inum , domain_id=nidom )
(  202)          CALL restget( inum, 'heatflow', jpidta, jpjdta, 1, 0, .FALSE., zdta )
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
(  214)          IF(lwp) WRITE(numout,*) '     bad flag value for ngeo_trc_flux = ', ngeo_trc_flux
(  215)          nstop = nstop + 1
(  216) 
(  217)       END SELECT
(  218) 
(  219)       ! geothermal heat flux trend
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:50:23      page 5

(  221)       SELECT CASE ( ngeo_trc_flux )
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
(  239)    END SUBROUTINE trc_bbc_init
(  240) 
(  241) #else
(  242)    !!----------------------------------------------------------------------
(  243)    !!   Default option                                         Empty module
(  244)    !!----------------------------------------------------------------------
(  245)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbc = .FALSE.  !: bbc flag
(  246) CONTAINS
(  247)    SUBROUTINE trc_bbc( kt )           ! Empty routine
(  248)       INTEGER, INTENT(in) :: kt
(  249) !      WRITE(*,*) 'trc_bbc: You should not have seen this print! error?', kt
(  250)    END SUBROUTINE trc_bbc
(  251) #endif
(  252) 
(  253)    !!======================================================================
(  254) END MODULE trcbbc






























PGF90 (Version     10.5)          09/23/2010  09:50:23      page 6

(    1) # 1 "trcbbc.F90"
(    1) MODULE trcbbc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcbbc  ***
(    4)    !! Ocean passive tracers:  bottom boundary condition
(    5)    !!==============================================================================
(    6) # 242
(  242)    !!----------------------------------------------------------------------
(  243)    !!   Default option                                         Empty module
(  244)    !!----------------------------------------------------------------------
(  245)    LOGICAL, PUBLIC, PARAMETER ::   lk_trcbbc = .FALSE.  !: bbc flag
(  246) CONTAINS
(  247)    SUBROUTINE trc_bbc( kt )           ! Empty routine
(  248)       INTEGER, INTENT(in) :: kt
(  249) !      WRITE(*,*) 'trc_bbc: You should not have seen this print! error?', kt
(  250)    END SUBROUTINE trc_bbc
(  251) 
(  252) # 253
(  253)    !!======================================================================
(  254) END MODULE trcbbc
