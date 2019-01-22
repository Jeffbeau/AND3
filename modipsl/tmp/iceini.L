


PGF90 (Version     10.5)          09/23/2010  09:49:53      page 1

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

Filename: iceini.F90

(    1) MODULE iceini
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE iceini   ***
(    4)    !!   Sea-ice model : LIM Sea ice model Initialization
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ice_init       : sea-ice model initialization
(   11)    !!----------------------------------------------------------------------
(   12)    USE dom_oce
(   13)    USE in_out_manager
(   14)    USE ice_oce         ! ice variables
(   15)    USE flx_oce
(   16)    USE phycst          ! Define parameters for the routines
(   17)    USE ocfzpt
(   18)    USE ice
(   19)    USE limmsh
(   20)    USE limistate
(   21)    USE limrst
(   22)    USE ini1d           ! initialization of the 1D configuration
(   23) !!DB
(   24)    USE restart
(   25) 
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC ice_init                 ! called by opa.F90
(   32) 
(   33)    !! * Share Module variables
(   34)    LOGICAL , PUBLIC  ::   & !!! ** init namelist (namicerun) **
(   35)       ln_limdyn   = .TRUE.   !: flag for ice dynamics (T) or not (F)
(   36)    INTEGER , PUBLIC  ::   &  !:
(   37)       nstart ,            &  !: iteration number of the begining of the run 
(   38)       nlast  ,            &  !: iteration number of the end of the run 
(   39)       nitrun ,            &  !: number of iteration
(   40)       numit                  !: iteration number
(   41)    REAL(wp), PUBLIC  ::   &  !:
(   42)       hsndif = 0.e0 ,     &  !: computation of temp. in snow (0) or not (9999)
(   43)       hicdif = 0.e0 ,     &  !: computation of temp. in ice (0) or not (9999)
(   44)       tpstot                 !: time of the run in seconds
(   45)    REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:
(   46)       acrit  = (/ 1.e-06 , 1.e-06 /)    !: minimum fraction for leads in 






PGF90 (Version     10.5)          09/23/2010  09:49:53      page 2

(   47)       !                                   !  north and south hemisphere
(   48)    !!----------------------------------------------------------------------
(   49)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   50)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/iceini.F90,v 1.6 2006/03/10 10:35:42 opalod Exp $ 
(   51)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) CONTAINS
(   55) 
(   56)    SUBROUTINE ice_init
(   57)       !!----------------------------------------------------------------------
(   58)       !!                  ***  ROUTINE ice_init  ***
(   59)       !!
(   60)       !! ** purpose :   
(   61)       !!
(   62)       !! History :
(   63)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and modules
(   64)       !!----------------------------------------------------------------------
(   65)        CHARACTER(len=80) :: namelist_icename
(   66)        
(   67)       ! Open the namelist file 
(   68)       namelist_icename = 'namelist_ice'
(   69)            
(   70)       CALL ctlopn(numnam_ice,namelist_icename,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
(   71)                      1,numout,.FALSE.,1)      
(   72) 
(   73)       CALL ice_run                    !  read in namelist some run parameters
(   74)                  
(   75)       ! Louvain la Neuve Ice model
(   76)       IF( nacc == 1 ) THEN
(   77)           dtsd2   = nfice * rdtmin * 0.5
(   78)           rdt_ice = nfice * rdtmin
(   79)       ELSE
(   80)           dtsd2   = nfice * rdt * 0.5
(   81)           rdt_ice = nfice * rdt
(   82)       ENDIF
(   83) 
(   84)       CALL lim_msh                    ! ice mesh initialization
(   85)      
(   86)       ! Initial sea-ice state
(   87)       IF( .NOT.ln_rstart ) THEN
(   88)          numit = 0
(   89)          CALL lim_istate              ! start from rest: sea-ice deduced from sst
(   90)       ELSE
(   91) !!DB
(   92)          call rst_ice_read(numit)
(   93)       ENDIF
(   94)       
(   95)       tn_ice(:,:) = sist(:,:)         ! initialisation of ice temperature   
(   96)       freeze(:,:) = 1.0 - frld(:,:)   ! initialisation of sea/ice cover    
(   97) # if defined key_coupled
(   98)       alb_ice(:,:) = albege(:,:)      ! sea-ice albedo
(   99) # endif
(  100)       
(  101)       nstart = numit  + nfice      
(  102)       nitrun = nitend - nit000 + 1 
(  103)       nlast  = numit  + nitrun 
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:53      page 3

(  105)       IF( nstock == 0  )  nstock = nlast + 1
(  106) 
(  107)    END SUBROUTINE ice_init
(  108) 
(  109) 
(  110)    SUBROUTINE ice_run
(  111)       !!-------------------------------------------------------------------
(  112)       !!                  ***  ROUTINE ice_run ***
(  113)       !!                 
(  114)       !! ** Purpose :   Definition some run parameter for ice model
(  115)       !!
(  116)       !! ** Method  :   Read the namicerun namelist and check the parameter 
(  117)       !!       values called at the first timestep (nit000)
(  118)       !!
(  119)       !! ** input   :   Namelist namicerun
(  120)       !!
(  121)       !! history :
(  122)       !!   2.0  !  03-08 (C. Ethe)  Original code
(  123)       !!-------------------------------------------------------------------
(  124) 
(  125)       NAMELIST/namicerun/ ln_limdyn, acrit, hsndif, hicdif
(  126)       !!-------------------------------------------------------------------
(  127) 
(  128)       !                                           ! Read Namelist namicerun 
(  129)       REWIND ( numnam_ice )
(  130)       READ   ( numnam_ice , namicerun )
(  131) 
(  132)       IF( lk_cfg_1d  )  ln_limdyn = .FALSE.       ! No ice transport in 1D configuration
(  133) 
(  134)       IF(lwp) THEN
(  135)          WRITE(numout,*)
(  136)          WRITE(numout,*) 'ice_run : ice share parameters for dynamics/advection/thermo of sea-ice'
(  137)          WRITE(numout,*) ' ~~~~~~'
(  138)          WRITE(numout,*) '   switch for ice dynamics (1) or not (0)      ln_limdyn   = ', ln_limdyn
(  139)          WRITE(numout,*) '   minimum fraction for leads in the NH (SH)  acrit(1/2)   = ', acrit(:)
(  140)          WRITE(numout,*) '   computation of temp. in snow (=0) or not (=9999) hsndif = ', hsndif
(  141)          WRITE(numout,*) '   computation of temp. in ice  (=0) or not (=9999) hicdif = ', hicdif
(  142)       ENDIF
(  143)    END SUBROUTINE ice_run
(  144) 
(  145) #else
(  146)    !!----------------------------------------------------------------------
(  147)    !!   Default option :        Empty module           NO LIM sea-ice model
(  148)    !!----------------------------------------------------------------------
(  149) CONTAINS
(  150)    SUBROUTINE ice_init        ! Empty routine
(  151)    END SUBROUTINE ice_init
(  152) #endif
(  153) 
(  154)    !!======================================================================
(  155) END MODULE iceini













PGF90 (Version     10.5)          09/23/2010  09:49:53      page 4

(    1) # 1 "iceini.F90"
(    1) MODULE iceini
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE iceini   ***
(    4)    !!   Sea-ice model : LIM Sea ice model Initialization
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   ice_init       : sea-ice model initialization
(   11)    !!----------------------------------------------------------------------
(   12)    USE dom_oce
(   13)    USE in_out_manager
(   14)    USE ice_oce         ! ice variables
(   15)    USE flx_oce
(   16)    USE phycst          ! Define parameters for the routines
(   17)    USE ocfzpt
(   18)    USE ice
(   19)    USE limmsh
(   20)    USE limistate
(   21)    USE limrst
(   22)    USE ini1d           ! initialization of the 1D configuration
(   23) !!DB
(   24)    USE restart
(   25) 
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC ice_init                 ! called by opa.F90
(   32) 
(   33)    !! * Share Module variables
(   34)    LOGICAL , PUBLIC  ::   & !!! ** init namelist (namicerun) **
(   35)       ln_limdyn   = .TRUE.   !: flag for ice dynamics (T) or not (F)
(   36)    INTEGER , PUBLIC  ::   &  !:
(   37)       nstart ,            &  !: iteration number of the begining of the run 
(   38)       nlast  ,            &  !: iteration number of the end of the run 
(   39)       nitrun ,            &  !: number of iteration
(   40)       numit                  !: iteration number
(   41)    REAL(wp), PUBLIC  ::   &  !:
(   42)       hsndif = 0.e0 ,     &  !: computation of temp. in snow (0) or not (9999)
(   43)       hicdif = 0.e0 ,     &  !: computation of temp. in ice (0) or not (9999)
(   44)       tpstot                 !: time of the run in seconds
(   45)    REAL(wp), PUBLIC, DIMENSION(2)  ::  &  !:
(   46)       acrit  = (/ 1.e-06 , 1.e-06 /)    !: minimum fraction for leads in 
(   47)       !                                   !  north and south hemisphere
(   48)    !!----------------------------------------------------------------------
(   49)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   50)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/iceini.F90,v 1.6 2006/03/10 10:35:42 opalod Exp $ 
(   51)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) CONTAINS
(   55) 
(   56)    SUBROUTINE ice_init
(   57)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:53      page 5

(   58)       !!                  ***  ROUTINE ice_init  ***
(   59)       !!
(   60)       !! ** purpose :   
(   61)       !!
(   62)       !! History :
(   63)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and modules
(   64)       !!----------------------------------------------------------------------
(   65)        CHARACTER(len=80) :: namelist_icename
(   66)        
(   67)       ! Open the namelist file 
(   68)       namelist_icename = 'namelist_ice'
(   69)            
(   70)       CALL ctlopn(numnam_ice,namelist_icename,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
(   71)                      1,numout,.FALSE.,1)      
(   72) 
(   73)       CALL ice_run                    !  read in namelist some run parameters
(   74)                  
(   75)       ! Louvain la Neuve Ice model
(   76)       IF( nacc == 1 ) THEN
(   77)           dtsd2   = nfice * rdtmin * 0.5
(   78)           rdt_ice = nfice * rdtmin
(   79)       ELSE
(   80)           dtsd2   = nfice * rdt * 0.5
(   81)           rdt_ice = nfice * rdt
(   82)       ENDIF
(   83) 
(   84)       CALL lim_msh                    ! ice mesh initialization
(   85)      
(   86)       ! Initial sea-ice state
(   87)       IF( .NOT.ln_rstart ) THEN
(   88)          numit = 0
(   89)          CALL lim_istate              ! start from rest: sea-ice deduced from sst
(   90)       ELSE
(   91) !!DB
(   92)          call rst_ice_read(numit)
(   93)       ENDIF
(   94)       
(   95)       tn_ice(:,:) = sist(:,:)         ! initialisation of ice temperature   
(   96)       freeze(:,:) = 1.0 - frld(:,:)   ! initialisation of sea/ice cover    
(   97) # 100
(  100)       
(  101)       nstart = numit  + nfice      
(  102)       nitrun = nitend - nit000 + 1 
(  103)       nlast  = numit  + nitrun 
(  104) 
(  105)       IF( nstock == 0  )  nstock = nlast + 1
(  106) 
(  107)    END SUBROUTINE ice_init
(  108) 
(  109) 
(  110)    SUBROUTINE ice_run
(  111)       !!-------------------------------------------------------------------
(  112)       !!                  ***  ROUTINE ice_run ***
(  113)       !!                 
(  114)       !! ** Purpose :   Definition some run parameter for ice model
(  115)       !!
(  116)       !! ** Method  :   Read the namicerun namelist and check the parameter 
(  117)       !!       values called at the first timestep (nit000)






PGF90 (Version     10.5)          09/23/2010  09:49:53      page 6

(  118)       !!
(  119)       !! ** input   :   Namelist namicerun
(  120)       !!
(  121)       !! history :
(  122)       !!   2.0  !  03-08 (C. Ethe)  Original code
(  123)       !!-------------------------------------------------------------------
(  124) 
(  125)       NAMELIST/namicerun/ ln_limdyn, acrit, hsndif, hicdif
(  126)       !!-------------------------------------------------------------------
(  127) 
(  128)       !                                           ! Read Namelist namicerun 
(  129)       REWIND ( numnam_ice )
(  130)       READ   ( numnam_ice , namicerun )
(  131) 
(  132)       IF( lk_cfg_1d  )  ln_limdyn = .FALSE.       ! No ice transport in 1D configuration
(  133) 
(  134)       IF(lwp) THEN
(  135)          WRITE(numout,*)
(  136)          WRITE(numout,*) 'ice_run : ice share parameters for dynamics/advection/thermo of sea-ice'
(  137)          WRITE(numout,*) ' ~~~~~~'
(  138)          WRITE(numout,*) '   switch for ice dynamics (1) or not (0)      ln_limdyn   = ', ln_limdyn
(  139)          WRITE(numout,*) '   minimum fraction for leads in the NH (SH)  acrit(1/2)   = ', acrit(:)
(  140)          WRITE(numout,*) '   computation of temp. in snow (=0) or not (=9999) hsndif = ', hsndif
(  141)          WRITE(numout,*) '   computation of temp. in ice  (=0) or not (=9999) hicdif = ', hicdif
(  142)       ENDIF
(  143)    END SUBROUTINE ice_run
(  144) 
(  145) 
(  146) # 154
(  154)    !!======================================================================
(  155) END MODULE iceini
