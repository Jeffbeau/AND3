


PGF90 (Version     10.5)          09/23/2010  09:49:11      page 1

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

Filename: dtasst.F90

(    1) !!DB: 2009.09.02 -- disabled thie key
(    2) !!sends message and increments nstop flag if this key is on
(    3) MODULE dtasst
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  dtasst  ***
(    6)    !! Data : Sea Surface Temperature (SST)
(    7)    
(    8)    !!      BUG initialisation  nyearsst !!!!!!bug
(    9)    
(   10)    !!======================================================================
(   11)    
(   12)    !!----------------------------------------------------------------------
(   13)    !!   dta_sst      : Reynolds sst data
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE ocfzpt          ! ???
(   20)    USE daymod          ! calendar
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Shared routine
(   26)    PUBLIC dta_sst
(   27) 
(   28)    !! * Shared module variables
(   29) #if defined key_dtasst
(   30)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasst = .TRUE.   !: sst data flag
(   31) #else
(   32)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasst = .FALSE.  !: sst data flag
(   33) #endif
(   34)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   35)       sst             !: surface temperature
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2) ::   &  !:
(   37)       rclice          !: climatological ice index (0/1) (2 months)
(   38) 
(   39) CONTAINS
(   40) 
(   41) #if defined key_dtasst
(   42)    !!----------------------------------------------------------------------
(   43)    !!   'key_dtasst'                                               SST data
(   44)    !!----------------------------------------------------------------------
(   45) 
(   46)    SUBROUTINE dta_sst( kt )






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 2

(   47)       
(   48)       !! * Arguments
(   49)       INTEGER ::   kt
(   50) 
(   51)       IF( kt == nit000 ) THEN
(   52)          IF(lwp) WRITE(numout,*)
(   53)          IF(lwp) WRITE(numout,*) 'DB: Defunct key dtasst defined ----> stop'
(   54)          IF(lwp) WRITE(numout,*) '    If you want this, find old code'
(   55)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   56)          nstop = nstop + 1
(   57)       ENDIF
(   58)    END SUBROUTINE dta_sst
(   59) 
(   60) #else
(   61)    !!----------------------------------------------------------------------
(   62)    !!   Default option :                                        NO SST data
(   63)    !!----------------------------------------------------------------------
(   64) 
(   65)    SUBROUTINE dta_sst( kt )
(   66)       !!---------------------------------------------------------------------
(   67)       !!                  ***  ROUTINE dta_sst  ***
(   68)       !!                    
(   69)       !! ** Purpose :   sea surface temperature data and update it
(   70)       !!     at each time step.   ???
(   71)       !!
(   72)       !! ** Method  : - sst   = tn
(   73)       !!              - rclice = 1. IF tn =< ztgel
(   74)       !!
(   75)       !! History :
(   76)       !!        !  91-03  ()  Original code
(   77)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   78)       !!----------------------------------------------------------------------
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt      ! ocean timestep
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER :: ji, jj
(   84)       !!---------------------------------------------------------------------
(   85)       
(   86)       IF( kt == nit000 ) THEN
(   87)          IF(lwp) WRITE(numout,*)
(   88)          IF(lwp) WRITE(numout,*) 'dta_sst : No SST data'
(   89)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   90)       ENDIF
(   91)       
(   92)       ! 1. Update at each time step
(   93)       ! ---------------------------
(   94) 
(   95)       sst   (:,:)   = tn   (:,:,1)
(   96)       rclice(:,:,1) = tmask(:,:,1)
(   97)       DO jj = 1, jpj
(   98)          DO ji = 1, jpi
(   99)             IF( tn(ji,jj,1) >= fzptn(ji,jj) ) rclice(ji,jj,1) = 0.e0
(  100)          END DO
(  101)       END DO
(  102)       rclice(:,:,2) = rclice(:,:,1)
(  103)       
(  104)    END SUBROUTINE dta_sst






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 3

(  105) #endif
(  106) 
(  107)    !!======================================================================
(  108) END MODULE dtasst




























































PGF90 (Version     10.5)          09/23/2010  09:49:11      page 4

(    1) # 1 "dtasst.F90"
(    1) !!DB: 2009.09.02 -- disabled thie key
(    2) !!sends message and increments nstop flag if this key is on
(    3) MODULE dtasst
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  dtasst  ***
(    6)    !! Data : Sea Surface Temperature (SST)
(    7)    
(    8)    !!      BUG initialisation  nyearsst !!!!!!bug
(    9)    
(   10)    !!======================================================================
(   11)    
(   12)    !!----------------------------------------------------------------------
(   13)    !!   dta_sst      : Reynolds sst data
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE ocfzpt          ! ???
(   20)    USE daymod          ! calendar
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Shared routine
(   26)    PUBLIC dta_sst
(   27) 
(   28)    !! * Shared module variables
(   29) # 32
(   32)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasst = .FALSE.  !: sst data flag
(   33) # 34
(   34)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   35)       sst             !: surface temperature
(   36)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2) ::   &  !:
(   37)       rclice          !: climatological ice index (0/1) (2 months)
(   38) 
(   39) CONTAINS
(   40) 
(   41) # 61
(   61)    !!----------------------------------------------------------------------
(   62)    !!   Default option :                                        NO SST data
(   63)    !!----------------------------------------------------------------------
(   64) 
(   65)    SUBROUTINE dta_sst( kt )
(   66)       !!---------------------------------------------------------------------
(   67)       !!                  ***  ROUTINE dta_sst  ***
(   68)       !!                    
(   69)       !! ** Purpose :   sea surface temperature data and update it
(   70)       !!     at each time step.   ???
(   71)       !!
(   72)       !! ** Method  : - sst   = tn
(   73)       !!              - rclice = 1. IF tn =< ztgel
(   74)       !!
(   75)       !! History :
(   76)       !!        !  91-03  ()  Original code
(   77)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   78)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 5

(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt      ! ocean timestep
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER :: ji, jj
(   84)       !!---------------------------------------------------------------------
(   85)       
(   86)       IF( kt == nit000 ) THEN
(   87)          IF(lwp) WRITE(numout,*)
(   88)          IF(lwp) WRITE(numout,*) 'dta_sst : No SST data'
(   89)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   90)       ENDIF
(   91)       
(   92)       ! 1. Update at each time step
(   93)       ! ---------------------------
(   94) 
(   95)       sst   (:,:)   = tn   (:,:,1)
(   96)       rclice(:,:,1) = tmask(:,:,1)
(   97)       DO jj = 1, jpj
(   98)          DO ji = 1, jpi
(   99)             IF( tn(ji,jj,1) >= fzptn(ji,jj) ) rclice(ji,jj,1) = 0.e0
(  100)          END DO
(  101)       END DO
(  102)       rclice(:,:,2) = rclice(:,:,1)
(  103)       
(  104)    END SUBROUTINE dta_sst
(  105) 
(  106) # 107
(  107)    !!======================================================================
(  108) END MODULE dtasst
