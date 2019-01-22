


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

Filename: dtasss.F90

(    1) !!DB: 2009.09.02 -- disabled thie key
(    2) !!sends message and increments nstop flag if this key is on
(    3) MODULE dtasss
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  dtasss  ***
(    6)    !! Data : Sea Surface Salinity (SSS)
(    7)    !!======================================================================
(    8)    
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dta_sss      : sss data
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE ocfzpt          ! ???
(   17)    USE daymod          ! calendar
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Shared routine
(   23)    PUBLIC dta_sss
(   24) 
(   25)    !! * Shared module variables
(   26) #if defined key_dtasss
(   27)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasss = .TRUE.   !: sss data flag
(   28) #else
(   29)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasss = .FALSE.  !: sss data flag
(   30) #endif
(   31) !   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   32) !      sss             !: surface salinity
(   33) 
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , IPSL-LODYC  (2005)
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40) #if defined key_dtasss
(   41)    !!----------------------------------------------------------------------
(   42)    !!   'key_dtasss'                                               SSS data
(   43)    !!----------------------------------------------------------------------
(   44) 
(   45)    SUBROUTINE dta_sss( kt )
(   46)       !!---------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 2

(   47)       !!                  ***  ROUTINE dta_sss  ***
(   48)       !!                    
(   49)       !! ** Purpose :   Read surface salinity data 
(   50)       !!
(   51)       !! ** Method  : - Read a specific sss.
(   52)       !!
(   53)       !! ** Action  : - sss 
(   54)       !!
(   55)       !! History :
(   56)       !!        !  90-03  (O. Marti and Ph Dandin)  Original code
(   57)       !!        !  92-07  (M. Imbard)
(   58)       !!        !  96-11  (E. Guilyardi)  Daily AGCM input files
(   59)       !!        !  00-04  (M. Imbard)  NetCDF FORMAT
(   60)       !!        !  00-10  (J.-P. Boulanger)  passage ORCA a TDH
(   61)       !!        !  01-10  (A. Lazar)  Reynolds default
(   62)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   63)       !!        !  02-11  (C. Levy)  MPP/MPI NetCDF read
(   64)       !!        !  05-03  (M. Levy) adapt SST to SSS
(   65)       !!----------------------------------------------------------------------
(   66)       
(   67)       !! * Arguments
(   68)       INTEGER ::   kt
(   69) 
(   70)       IF( kt == nit000 ) THEN
(   71)          IF(lwp) WRITE(numout,*)
(   72)          IF(lwp) WRITE(numout,*) 'DB: Defunct key dtasss defined ----> stop'
(   73)          IF(lwp) WRITE(numout,*) '    If you want this, find old code'
(   74)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(   75)          nstop = nstop + 1
(   76)       ENDIF
(   77) 
(   78)    END SUBROUTINE dta_sss
(   79) 
(   80) #else
(   81)    !!----------------------------------------------------------------------
(   82)    !!   Default option :                                        NO SSS data
(   83)    !!----------------------------------------------------------------------
(   84) 
(   85)    SUBROUTINE dta_sss( kt )
(   86)       !!---------------------------------------------------------------------
(   87)       !!                  ***  ROUTINE dta_sss  ***
(   88)       !!                    
(   89)       !! ** Purpose :   sea surface salinity data and update it
(   90)       !!     at each time step.   ???
(   91)       !!
(   92)       !! ** Method  : - sss  
(   93)       !!
(   94)       !! History :
(   95)       !!        !  91-03  ()  Original code
(   96)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   97)       !!----------------------------------------------------------------------
(   98)       !! * Arguments
(   99)       INTEGER, INTENT( in ) ::   kt      ! ocean timestep
(  100)       !!---------------------------------------------------------------------
(  101)       
(  102)       IF( kt == nit000 ) THEN
(  103)          IF(lwp) WRITE(numout,*)
(  104)          IF(lwp) WRITE(numout,*) 'dta_sss : No SSS data'






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 3

(  105)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  106)       ENDIF
(  107) 
(  108)       
(  109)    END SUBROUTINE dta_sss
(  110) #endif
(  111) 
(  112)    !!======================================================================
(  113) END MODULE dtasss























































PGF90 (Version     10.5)          09/23/2010  09:49:11      page 4

(    1) # 1 "dtasss.F90"
(    1) !!DB: 2009.09.02 -- disabled thie key
(    2) !!sends message and increments nstop flag if this key is on
(    3) MODULE dtasss
(    4)    !!======================================================================
(    5)    !!                       ***  MODULE  dtasss  ***
(    6)    !! Data : Sea Surface Salinity (SSS)
(    7)    !!======================================================================
(    8)    
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dta_sss      : sss data
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE ocfzpt          ! ???
(   17)    USE daymod          ! calendar
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Shared routine
(   23)    PUBLIC dta_sss
(   24) 
(   25)    !! * Shared module variables
(   26) # 29
(   29)    LOGICAL , PUBLIC, PARAMETER ::   lk_dtasss = .FALSE.  !: sss data flag
(   30) # 31
(   31) !   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
(   32) !      sss             !: surface salinity
(   33) 
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , IPSL-LODYC  (2005)
(   36)    !!----------------------------------------------------------------------
(   37) 
(   38) CONTAINS
(   39) 
(   40) # 81
(   81)    !!----------------------------------------------------------------------
(   82)    !!   Default option :                                        NO SSS data
(   83)    !!----------------------------------------------------------------------
(   84) 
(   85)    SUBROUTINE dta_sss( kt )
(   86)       !!---------------------------------------------------------------------
(   87)       !!                  ***  ROUTINE dta_sss  ***
(   88)       !!                    
(   89)       !! ** Purpose :   sea surface salinity data and update it
(   90)       !!     at each time step.   ???
(   91)       !!
(   92)       !! ** Method  : - sss  
(   93)       !!
(   94)       !! History :
(   95)       !!        !  91-03  ()  Original code
(   96)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   97)       !!----------------------------------------------------------------------
(   98)       !! * Arguments
(   99)       INTEGER, INTENT( in ) ::   kt      ! ocean timestep






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 5

(  100)       !!---------------------------------------------------------------------
(  101)       
(  102)       IF( kt == nit000 ) THEN
(  103)          IF(lwp) WRITE(numout,*)
(  104)          IF(lwp) WRITE(numout,*) 'dta_sss : No SSS data'
(  105)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  106)       ENDIF
(  107) 
(  108)       
(  109)    END SUBROUTINE dta_sss
(  110) 
(  111) # 112
(  112)    !!======================================================================
(  113) END MODULE dtasss
