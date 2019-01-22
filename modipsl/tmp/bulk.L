


PGF90 (Version     10.5)          09/23/2010  09:49:12      page 1

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

Filename: bulk.F90

(    1) MODULE bulk
(    2)    !!======================================================================
(    3)    !!                           ***  bulk  ***
(    4)    !!======================================================================
(    5) #if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
(    6)    !!----------------------------------------------------------------------
(    7)    !!   'key_flx_bulk_monthly'                        monthly bulk formulea
(    8)    !!   'key_flx_bulk_daily'                          daily bulk formulea
(    9)    !!----------------------------------------------------------------------
(   10)    !!   bulk          : computation of fluxes using bulk formulation
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used   
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE ice_oce         ! bulk variable  
(   16)    USE ocfzpt          ! ocean freezing point
(   17)    USE flxblk          ! bulk formulae
(   18)    USE blk_oce         ! bulk variable 
(   19)    USE flx_oce
(   20)    USE taumod
(   21)    USE phycst          ! physical constants
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC blk        ! called by flx.F90   
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/bulk.F90,v 1.9 2005/09/22 10:58:15 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE blk( kt )
(   38)       !!---------------------------------------------------------------------
(   39)       !!                    ***  ROUTINE blk  ***
(   40)       !!        
(   41)       !! ** Purpose :   provide the heat fluxes on ice and ocean 
(   42)       !!                using bulk formulation
(   43)       !!
(   44)       !! History :
(   45)       !!   9.0  !  03-11  (C. Ethe and G. Madec)  F90: Free form and MODULE 
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 2

(   47)       !! * arguments
(   48)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(   49) 
(   50)       !! * Local declarations    
(   51)       REAL(wp), DIMENSION(jpi,jpj) ::   zsst 
(   52) # if ! defined key_ice_lim
(   53)       INTEGER  ::   ji, jj         ! dummy loop indices  
(   54)       REAL(wp) ::   ztgel, zicopa
(   55) # endif
(   56)       !!---------------------------------------------------------------------
(   57) 
(   58)      ! Initialisation
(   59)      IF( kt == nit000) THEN
(   60)       ! computation of rdtbs2
(   61)         IF( nacc == 1 ) THEN
(   62)            rdtbs2 = nfbulk * rdtmin * 0.5
(   63)         ELSE
(   64)            rdtbs2 = nfbulk * rdt * 0.5
(   65)         ENDIF
(   66)         IF ( .NOT.ln_rstart ) THEN
(   67)            gsst(:,:) =  ( nfbulk - 1 ) * ( tn(:,:,1) + rt0 )
(   68)         ENDIF
(   69)      ENDIF
(   70) 
(   71) # if ! defined key_ice_lim
(   72)       ! opa model ice freeze()      
(   73)       DO jj = 1, jpj
(   74)          DO ji = 1, jpi
(   75)             ztgel  = fzptn(ji,jj)
(   76)             zicopa = tmask(ji,jj,1)
(   77)             IF( tn(ji,jj,1) >= ztgel ) zicopa = 0.
(   78)             freeze(ji,jj) = zicopa
(   79)          END DO
(   80)       END DO
(   81) # endif
(   82) 
(   83)       gsst(:,:) = gsst(:,:) + tn(:,:,1) + rt0  
(   84) 
(   85)       !  Computation of the fluxes       
(   86)       IF( MOD( kt - 1 , nfbulk ) == 0 ) THEN
(   87) 
(   88)          zsst(:,:) = gsst(:,:) / REAL( nfbulk ) * tmask(:,:,1)
(   89)          CALL flx_blk( zsst )    
(   90)          gsst(:,:) = 0.    
(   91) 
(   92) # if ! defined key_ice_lim
(   93)          IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(   94)             CALL prt_ctl_info(' Forcings ')
(   95)             CALL prt_ctl(tab2d_1=qsr_oce , clinfo1=' qsr_oce   : ', mask1=tmask, ovlap=1)
(   96)             CALL prt_ctl(tab2d_1=qsr_ice , clinfo1=' qsr_ice   : ', mask1=tmask, ovlap=1)
(   97)             CALL prt_ctl(tab2d_1=qnsr_oce, clinfo1=' qnsr_oce  : ', mask1=tmask, ovlap=1)
(   98)             CALL prt_ctl(tab2d_1=qnsr_ice, clinfo1=' qnsr_ice  : ', mask1=tmask, ovlap=1)
(   99)             CALL prt_ctl(tab2d_1=evap    , clinfo1=' evap      : ', mask1=tmask, ovlap=1)
(  100)             CALL prt_ctl(tab2d_1=tprecip , clinfo1=' precip    : ', mask1=tmask, ovlap=1)
(  101)             CALL prt_ctl(tab2d_1=sprecip , clinfo1=' Snow      : ', mask1=tmask, ovlap=1)
(  102)             CALL prt_ctl(tab2d_1=taux    , clinfo1=' u-stress  : ', mask1=tmask, ovlap=1)
(  103)             CALL prt_ctl(tab2d_1=tauy    , clinfo1=' v-stress  : ', mask1=tmask, ovlap=1)
(  104)             CALL prt_ctl(tab2d_1=zsst    , clinfo1=' sst       : ', mask1=tmask, ovlap=1)






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 3

(  105)          ENDIF
(  106) # endif   
(  107)       ENDIF
(  108)  
(  109)    END SUBROUTINE blk
(  110) 
(  111) #else
(  112)    !!----------------------------------------------------------------------
(  113)    !!   Dummy module :                                     NO bulk formulea
(  114)    !!----------------------------------------------------------------------
(  115) CONTAINS
(  116)    SUBROUTINE blk( kt )          ! Dummy routine
(  117)       if(lwp) WRITE(numout,*) 'blk: You should not see this print! error? ', kt
(  118)    END SUBROUTINE blk
(  119) #endif
(  120)  
(  121)    !!======================================================================
(  122) END MODULE bulk














































PGF90 (Version     10.5)          09/23/2010  09:49:12      page 4

(    1) # 1 "bulk.F90"
(    1) MODULE bulk
(    2)    !!======================================================================
(    3)    !!                           ***  bulk  ***
(    4)    !!======================================================================
(    5) # 6
(    6)    !!----------------------------------------------------------------------
(    7)    !!   'key_flx_bulk_monthly'                        monthly bulk formulea
(    8)    !!   'key_flx_bulk_daily'                          daily bulk formulea
(    9)    !!----------------------------------------------------------------------
(   10)    !!   bulk          : computation of fluxes using bulk formulation
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used   
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE ice_oce         ! bulk variable  
(   16)    USE ocfzpt          ! ocean freezing point
(   17)    USE flxblk          ! bulk formulae
(   18)    USE blk_oce         ! bulk variable 
(   19)    USE flx_oce
(   20)    USE taumod
(   21)    USE phycst          ! physical constants
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE prtctl          ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC blk        ! called by flx.F90   
(   30)    !!----------------------------------------------------------------------
(   31)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/bulk.F90,v 1.9 2005/09/22 10:58:15 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) CONTAINS
(   36) 
(   37)    SUBROUTINE blk( kt )
(   38)       !!---------------------------------------------------------------------
(   39)       !!                    ***  ROUTINE blk  ***
(   40)       !!        
(   41)       !! ** Purpose :   provide the heat fluxes on ice and ocean 
(   42)       !!                using bulk formulation
(   43)       !!
(   44)       !! History :
(   45)       !!   9.0  !  03-11  (C. Ethe and G. Madec)  F90: Free form and MODULE 
(   46)       !!----------------------------------------------------------------------
(   47)       !! * arguments
(   48)       INTEGER, INTENT( in  ) ::   kt   ! ocean time step
(   49) 
(   50)       !! * Local declarations    
(   51)       REAL(wp), DIMENSION(jpi,jpj) ::   zsst 
(   52) # 56
(   56)       !!---------------------------------------------------------------------
(   57) 
(   58)      ! Initialisation
(   59)      IF( kt == nit000) THEN
(   60)       ! computation of rdtbs2






PGF90 (Version     10.5)          09/23/2010  09:49:12      page 5

(   61)         IF( nacc == 1 ) THEN
(   62)            rdtbs2 = nfbulk * rdtmin * 0.5
(   63)         ELSE
(   64)            rdtbs2 = nfbulk * rdt * 0.5
(   65)         ENDIF
(   66)         IF ( .NOT.ln_rstart ) THEN
(   67)            gsst(:,:) =  ( nfbulk - 1 ) * ( tn(:,:,1) + rt0 )
(   68)         ENDIF
(   69)      ENDIF
(   70) 
(   71) 
(   72) # 83
(   83)       gsst(:,:) = gsst(:,:) + tn(:,:,1) + rt0  
(   84) 
(   85)       !  Computation of the fluxes       
(   86)       IF( MOD( kt - 1 , nfbulk ) == 0 ) THEN
(   87) 
(   88)          zsst(:,:) = gsst(:,:) / REAL( nfbulk ) * tmask(:,:,1)
(   89)          CALL flx_blk( zsst )    
(   90)          gsst(:,:) = 0.    
(   91) 
(   92) # 107
(  107)       ENDIF
(  108)  
(  109)    END SUBROUTINE blk
(  110) 
(  111) # 120
(  120)  
(  121)    !!======================================================================
(  122) END MODULE bulk
