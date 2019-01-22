


PGF90 (Version     10.5)          09/23/2010  09:49:58      page 1

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

Filename: limflx.F90

(    1) MODULE limflx
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limflx   ***
(    4)    !!           computation of the flux at the sea ice/ocean interface
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_flx  : flux at the ice / ocean interface
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce          ! ocean parameters
(   14)    USE phycst           ! physical constants
(   15)    USE ocfzpt           ! surface ocean freezing point
(   16)    USE ice_oce          ! sea-ice variable
(   17)    USE flx_oce          ! sea-ice/ocean forcings variables
(   18)    USE ice              ! LIM sea-ice variables
(   19)    USE flxblk           ! bulk formulea
(   20)    USE lbclnk           ! ocean lateral boundary condition
(   21)    USE in_out_manager   ! I/O manager
(   22)    USE albedo           ! albedo parameters
(   23)    USE prtctl           ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC lim_flx       ! called by lim_step
(   30) 
(   31)    !! * Module variables
(   32)    REAL(wp)  ::           &  ! constant values
(   33)       epsi16 = 1.e-16  ,  &
(   34)       rzero  = 0.e0    ,  &
(   35)       rone   = 1.e0
(   36) 
(   37)    !! * Substitutions
(   38) #  include "vectopt_loop_substitute.h90"
(   39)    !!----------------------------------------------------------------------
(   40)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   41)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limflx.F90,v 1.8 2005/12/21 10:46:27 opalod Exp $ 
(   42)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   43)    !!----------------------------------------------------------------------
(   44) CONTAINS
(   45) 
(   46)    SUBROUTINE lim_flx






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 2

(   47)       !!-------------------------------------------------------------------
(   48)       !!                ***  ROUTINE lim_flx ***
(   49)       !!  
(   50)       !! ** Purpose : Computes the mass and heat fluxes to the ocean
(   51)       !!         
(   52)       !! ** Action  : - Initialisation of some variables
(   53)       !!              - comput. of the fluxes at the sea ice/ocean interface
(   54)       !!     
(   55)       !! ** Outputs : - fsolar  : solar heat flux at sea ice/ocean interface
(   56)       !!              - fnsolar : non solar heat flux 
(   57)       !!              - fsalt   : salt flux at sea ice/ocean interface
(   58)       !!              - fmass   : freshwater flux at sea ice/ocean interface
(   59)       !!
(   60)       !! ** References :
(   61)       !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
(   62)       !!         original    : 00-01 (LIM)
(   63)       !!         addition    : 02-07 (C. Ethe, G. Madec)
(   64)       !!---------------------------------------------------------------------
(   65)       !! * Local variables
(   66)       INTEGER ::   ji, jj         ! dummy loop indices
(   67) 
(   68)       INTEGER ::   &
(   69)          ifvt, i1mfr, idfr ,   &  ! some switches
(   70)          iflt, ial, iadv, ifral, ifrdv
(   71)       
(   72)       REAL(wp) ::   &
(   73)          zinda  ,              &  ! switch for testing the values of ice concentration
(   74)          z1mthcm                  ! 1 - thcm
(   75) !!         zfcm1  ,              &  ! solar  heat fluxes
(   76) !!         zfcm2  ,              &  !  non solar heat fluxes
(   77) #if defined key_lim_fdd   
(   78)       REAL(wp) ::   &
(   79)          zfons,                &  ! salt exchanges at the ice/ocean interface
(   80)          zpme                     ! freshwater exchanges at the ice/ocean interface
(   81) #else
(   82)       REAL(wp) ::   &
(   83)          zprs  , zfons,        &  ! salt exchanges at the ice/ocean interface
(   84)          zpmess                   ! freshwater exchanges at the ice/ocean interface
(   85) #endif
(   86)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   87)          zfcm1  ,              &  ! solar  heat fluxes
(   88)          zfcm2                    !  non solar heat fluxes      
(   89) #if defined key_coupled    
(   90)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   91)          zalb  ,               &  ! albedo of ice under overcast sky
(   92)          zalcn ,               &  ! albedo of ocean under overcast sky
(   93)          zalbp ,               &  ! albedo of ice under clear sky
(   94)          zaldum                   ! albedo of ocean under clear sky
(   95) #endif
(   96)       !!---------------------------------------------------------------------
(   97)      
(   98)       !---------------------------------!
(   99)       !      Sea ice/ocean interface    !
(  100)       !---------------------------------!
(  101)        
(  102)        
(  103)       !      heat flux at the ocean surface
(  104)       !-------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 3

(  105)        
(  106)       DO jj = 1, jpj
(  107)          DO ji = 1, jpi
(  108)             zinda   = 1.0 - MAX( rzero , SIGN( rone , - ( 1.0 - pfrld(ji,jj) ) ) )
(  109)             ifvt    = zinda  *  MAX( rzero , SIGN( rone, -phicif(ji,jj) ) )
(  110)             i1mfr   = 1.0 - MAX( rzero , SIGN( rone ,  - ( 1.0 - frld(ji,jj) ) ) )
(  111)             idfr    = 1.0 - MAX( rzero , SIGN( rone , frld(ji,jj) - pfrld(ji,jj) ) )
(  112)             iflt    = zinda  * (1 - i1mfr) * (1 - ifvt )
(  113)             ial     = ifvt   * i1mfr + ( 1 - ifvt ) * idfr
(  114)             iadv    = ( 1  - i1mfr ) * zinda
(  115)             ifral   = ( 1  - i1mfr * ( 1 - ial ) )   
(  116)             ifrdv   = ( 1  - ifral * ( 1 - ial ) ) * iadv 
(  117)             z1mthcm =   1. - thcm(ji,jj)       
(  118)             !   computation the solar flux at ocean surface
(  119)             zfcm1(ji,jj)   = pfrld(ji,jj) * qsr_oce(ji,jj)  + ( 1. - pfrld(ji,jj) ) * fstric(ji,jj)
(  120)             !  computation the non solar heat flux at ocean surface
(  121)             zfcm2(ji,jj) =  - z1mthcm * zfcm1(ji,jj)   &
(  122)                &           + iflt    * ( fscmbq(ji,jj) + ffltbif(ji,jj) )                            &
(  123)                &           + ifral   * ( ial * qcmif(ji,jj) + (1 - ial) * qldif(ji,jj) ) / rdt_ice   &
(  124)                &           + ifrdv   * ( qfvbq(ji,jj) + qdtcn(ji,jj) ) / rdt_ice
(  125) 
(  126)             fsbbq(ji,jj) = ( 1.0 - ( ifvt + iflt ) ) * fscmbq(ji,jj)     ! ???
(  127)             
(  128)             fsolar (ji,jj) = zfcm1(ji,jj)                                       ! solar heat flux 
(  129) 
(  130)             fnsolar(ji,jj) = zfcm2(ji,jj) - fdtcn(ji,jj)                        ! non solar heat flux
(  131)          END DO
(  132)       END DO
(  133)   
(  134)        
(  135)       !      mass flux at the ocean surface
(  136)       !-------------------------------------------------------
(  137)        
(  138)       DO jj = 1, jpj
(  139)          DO ji = 1, jpi
(  140) #if defined key_lim_fdd
(  141)             !  case of realistic freshwater flux (Tartinville et al., 2001)
(  142)             
(  143)             !  computing freshwater exchanges at the ice/ocean interface
(  144)             zpme = - evap(ji,jj) * frld(ji,jj)            &   !  evaporation over oceanic fraction
(  145)                &   + tprecip(ji,jj)                            &   !  total precipitation
(  146)                &   - sprecip(ji,jj) * ( 1. - pfrld(ji,jj) )  &   !  remov. snow precip over ice
(  147)                &   - rdmsnif(ji,jj) / rdt_ice                   !  freshwaterflux due to snow melting 
(  148)             
(  149)             !  computing salt exchanges at the ice/ocean interface
(  150)             zfons =  ( soce - sice ) * ( rdmicif(ji,jj) / rdt_ice ) 
(  151)             
(  152)             !  converting the salt flux from ice to a freshwater flux from ocean
(  153)             fsalt(ji,jj) = zfons / ( sss_io(ji,jj) + epsi16 )
(  154)             
(  155)             !  freshwater masses
(  156)             fmass(ji,jj) = - zpme 
(  157) #else
(  158)             !  case of freshwater flux equivalent as salt flux
(  159)             !  dilution effect due to evaporation and precipitation
(  160)             zprs  = ( tprecip(ji,jj) - sprecip(ji,jj) * ( 1. - pfrld(ji,jj) ) ) * soce  
(  161)             !  freshwater flux
(  162)             zfons = rdmicif(ji,jj) * ( soce - sice )  &  !  fwf : ice formation and melting






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 4

(  163)                &   -  dmgwi(ji,jj) * sice             &  !  fwf : salt flx needed to bring the fresh snow to sea/ice salinit
(  164)                &   + rdmsnif(ji,jj) * soce               !  fwf to ocean due to snow melting 
(  165)             !  salt exchanges at the ice/ocean interface
(  166)             zpmess         =  zprs - zfons / rdt_ice - evap(ji,jj) * soce * frld(ji,jj)
(  167)             fsalt(ji,jj) =  - zpmess
(  168) #endif
(  169)          END DO
(  170)       END DO
(  171) 
(  172) 
(  173)       !-------------------------------------------------------------------!
(  174)       !  computation of others transmitting variables from ice to ocean   !
(  175)       !------------------------------------------ ------------------------!
(  176) 
(  177)       !-----------------------------------------------!
(  178)       !   Storing the transmitted variables           !
(  179)       !-----------------------------------------------!
(  180) 
(  181)       ftaux (:,:) = - tio_u(:,:) * rau0   ! taux ( ice: N/m2/rau0, ocean: N/m2 )
(  182)       ftauy (:,:) = - tio_v(:,:) * rau0   ! tauy ( ice: N/m2/rau0, ocean: N/m2 )                
(  183)       freeze(:,:) = 1.0 - frld(:,:)       ! Sea ice cover            
(  184)       tn_ice(:,:) = sist(:,:)             ! Ice surface temperature                      
(  185) 
(  186) #if defined key_coupled            
(  187)       zalb  (:,:) = 0.e0
(  188)       zalcn (:,:) = 0.e0
(  189)       zalbp (:,:) = 0.e0
(  190)       zaldum(:,:) = 0.e0
(  191) 
(  192)       !------------------------------------------------!
(  193)       !  2) Computation of snow/ice and ocean albedo   !
(  194)       !------------------------------------------------!
(  195)       CALL flx_blk_albedo( zalb, zalcn, zalbp, zaldum )
(  196) 
(  197)       alb_ice(:,:) =  0.5 * zalbp(:,:) + 0.5 * zalb (:,:)   ! Ice albedo                       
(  198) #endif
(  199) 
(  200)       IF(ln_ctl) THEN
(  201)          CALL prt_ctl(tab2d_1=fsolar, clinfo1=' lim_flx: fsolar : ', tab2d_2=fnsolar, clinfo2=' fnsolar : ')
(  202)          CALL prt_ctl(tab2d_1=fmass , clinfo1=' lim_flx: fmass  : ', tab2d_2=fsalt  , clinfo2=' fsalt   : ')
(  203)          CALL prt_ctl(tab2d_1=ftaux , clinfo1=' lim_flx: ftaux  : ', tab2d_2=ftauy  , clinfo2=' ftauy   : ')
(  204)          CALL prt_ctl(tab2d_1=freeze, clinfo1=' lim_flx: freeze : ', tab2d_2=tn_ice , clinfo2=' tn_ice  : ')
(  205)       ENDIF 
(  206)    
(  207)     END SUBROUTINE lim_flx
(  208) 
(  209) #else
(  210)    !!----------------------------------------------------------------------
(  211)    !!   Default option :        Dummy module           NO LIM sea-ice model
(  212)    !!----------------------------------------------------------------------
(  213) CONTAINS
(  214)    SUBROUTINE lim_flx         ! Dummy routine
(  215)    END SUBROUTINE lim_flx
(  216) #endif 
(  217) 
(  218)    !!======================================================================
(  219) END MODULE limflx







PGF90 (Version     10.5)          09/23/2010  09:49:58      page 5

(    1) # 1 "limflx.F90"
(    1) MODULE limflx
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limflx   ***
(    4)    !!           computation of the flux at the sea ice/ocean interface
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_flx  : flux at the ice / ocean interface
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce          ! ocean parameters
(   14)    USE phycst           ! physical constants
(   15)    USE ocfzpt           ! surface ocean freezing point
(   16)    USE ice_oce          ! sea-ice variable
(   17)    USE flx_oce          ! sea-ice/ocean forcings variables
(   18)    USE ice              ! LIM sea-ice variables
(   19)    USE flxblk           ! bulk formulea
(   20)    USE lbclnk           ! ocean lateral boundary condition
(   21)    USE in_out_manager   ! I/O manager
(   22)    USE albedo           ! albedo parameters
(   23)    USE prtctl           ! Print control
(   24) 
(   25)    IMPLICIT NONE
(   26)    PRIVATE
(   27) 
(   28)    !! * Routine accessibility
(   29)    PUBLIC lim_flx       ! called by lim_step
(   30) 
(   31)    !! * Module variables
(   32)    REAL(wp)  ::           &  ! constant values
(   33)       epsi16 = 1.e-16  ,  &
(   34)       rzero  = 0.e0    ,  &
(   35)       rone   = 1.e0
(   36) 
(   37)    !! * Substitutions
(   38) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 39 "limflx.F90"
(   39) # 39
(   39)    !!----------------------------------------------------------------------
(   40)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   41)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limflx.F90,v 1.8 2005/12/21 10:46:27 opalod Exp $ 
(   42)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   43)    !!----------------------------------------------------------------------
(   44) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 6

(   45) 
(   46)    SUBROUTINE lim_flx
(   47)       !!-------------------------------------------------------------------
(   48)       !!                ***  ROUTINE lim_flx ***
(   49)       !!  
(   50)       !! ** Purpose : Computes the mass and heat fluxes to the ocean
(   51)       !!         
(   52)       !! ** Action  : - Initialisation of some variables
(   53)       !!              - comput. of the fluxes at the sea ice/ocean interface
(   54)       !!     
(   55)       !! ** Outputs : - fsolar  : solar heat flux at sea ice/ocean interface
(   56)       !!              - fnsolar : non solar heat flux 
(   57)       !!              - fsalt   : salt flux at sea ice/ocean interface
(   58)       !!              - fmass   : freshwater flux at sea ice/ocean interface
(   59)       !!
(   60)       !! ** References :
(   61)       !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
(   62)       !!         original    : 00-01 (LIM)
(   63)       !!         addition    : 02-07 (C. Ethe, G. Madec)
(   64)       !!---------------------------------------------------------------------
(   65)       !! * Local variables
(   66)       INTEGER ::   ji, jj         ! dummy loop indices
(   67) 
(   68)       INTEGER ::   &
(   69)          ifvt, i1mfr, idfr ,   &  ! some switches
(   70)          iflt, ial, iadv, ifral, ifrdv
(   71)       
(   72)       REAL(wp) ::   &
(   73)          zinda  ,              &  ! switch for testing the values of ice concentration
(   74)          z1mthcm                  ! 1 - thcm
(   75) !!         zfcm1  ,              &  ! solar  heat fluxes
(   76) !!         zfcm2  ,              &  !  non solar heat fluxes
(   77) # 78
(   78)       REAL(wp) ::   &
(   79)          zfons,                &  ! salt exchanges at the ice/ocean interface
(   80)          zpme                     ! freshwater exchanges at the ice/ocean interface
(   81) # 86
(   86)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   87)          zfcm1  ,              &  ! solar  heat fluxes
(   88)          zfcm2                    !  non solar heat fluxes      
(   89) # 96
(   96)       !!---------------------------------------------------------------------
(   97)      
(   98)       !---------------------------------!
(   99)       !      Sea ice/ocean interface    !
(  100)       !---------------------------------!
(  101)        
(  102)        
(  103)       !      heat flux at the ocean surface
(  104)       !-------------------------------------------------------
(  105)        
(  106)       DO jj = 1, jpj
(  107)          DO ji = 1, jpi
(  108)             zinda   = 1.0 - MAX( rzero , SIGN( rone , - ( 1.0 - pfrld(ji,jj) ) ) )
(  109)             ifvt    = zinda  *  MAX( rzero , SIGN( rone, -phicif(ji,jj) ) )
(  110)             i1mfr   = 1.0 - MAX( rzero , SIGN( rone ,  - ( 1.0 - frld(ji,jj) ) ) )
(  111)             idfr    = 1.0 - MAX( rzero , SIGN( rone , frld(ji,jj) - pfrld(ji,jj) ) )
(  112)             iflt    = zinda  * (1 - i1mfr) * (1 - ifvt )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 7

(  113)             ial     = ifvt   * i1mfr + ( 1 - ifvt ) * idfr
(  114)             iadv    = ( 1  - i1mfr ) * zinda
(  115)             ifral   = ( 1  - i1mfr * ( 1 - ial ) )   
(  116)             ifrdv   = ( 1  - ifral * ( 1 - ial ) ) * iadv 
(  117)             z1mthcm =   1. - thcm(ji,jj)       
(  118)             !   computation the solar flux at ocean surface
(  119)             zfcm1(ji,jj)   = pfrld(ji,jj) * qsr_oce(ji,jj)  + ( 1. - pfrld(ji,jj) ) * fstric(ji,jj)
(  120)             !  computation the non solar heat flux at ocean surface
(  121)             zfcm2(ji,jj) =  - z1mthcm * zfcm1(ji,jj)   &
(  122)                &           + iflt    * ( fscmbq(ji,jj) + ffltbif(ji,jj) )                            &
(  123)                &           + ifral   * ( ial * qcmif(ji,jj) + (1 - ial) * qldif(ji,jj) ) / rdt_ice   &
(  124)                &           + ifrdv   * ( qfvbq(ji,jj) + qdtcn(ji,jj) ) / rdt_ice
(  125) 
(  126)             fsbbq(ji,jj) = ( 1.0 - ( ifvt + iflt ) ) * fscmbq(ji,jj)     ! ???
(  127)             
(  128)             fsolar (ji,jj) = zfcm1(ji,jj)                                       ! solar heat flux 
(  129) 
(  130)             fnsolar(ji,jj) = zfcm2(ji,jj) - fdtcn(ji,jj)                        ! non solar heat flux
(  131)          END DO
(  132)       END DO
(  133)   
(  134)        
(  135)       !      mass flux at the ocean surface
(  136)       !-------------------------------------------------------
(  137)        
(  138)       DO jj = 1, jpj
(  139)          DO ji = 1, jpi
(  140) # 141
(  141)             !  case of realistic freshwater flux (Tartinville et al., 2001)
(  142)             
(  143)             !  computing freshwater exchanges at the ice/ocean interface
(  144)             zpme = - evap(ji,jj) * frld(ji,jj)            &   !  evaporation over oceanic fraction
(  145)                &   + tprecip(ji,jj)                            &   !  total precipitation
(  146)                &   - sprecip(ji,jj) * ( 1. - pfrld(ji,jj) )  &   !  remov. snow precip over ice
(  147)                &   - rdmsnif(ji,jj) / rdt_ice                   !  freshwaterflux due to snow melting 
(  148)             
(  149)             !  computing salt exchanges at the ice/ocean interface
(  150)             zfons =  ( soce - sice ) * ( rdmicif(ji,jj) / rdt_ice ) 
(  151)             
(  152)             !  converting the salt flux from ice to a freshwater flux from ocean
(  153)             fsalt(ji,jj) = zfons / ( sss_io(ji,jj) + epsi16 )
(  154)             
(  155)             !  freshwater masses
(  156)             fmass(ji,jj) = - zpme 
(  157) # 169
(  169)          END DO
(  170)       END DO
(  171) 
(  172) 
(  173)       !-------------------------------------------------------------------!
(  174)       !  computation of others transmitting variables from ice to ocean   !
(  175)       !------------------------------------------ ------------------------!
(  176) 
(  177)       !-----------------------------------------------!
(  178)       !   Storing the transmitted variables           !
(  179)       !-----------------------------------------------!
(  180) 
(  181)       ftaux (:,:) = - tio_u(:,:) * rau0   ! taux ( ice: N/m2/rau0, ocean: N/m2 )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 8

(  182)       ftauy (:,:) = - tio_v(:,:) * rau0   ! tauy ( ice: N/m2/rau0, ocean: N/m2 )                
(  183)       freeze(:,:) = 1.0 - frld(:,:)       ! Sea ice cover            
(  184)       tn_ice(:,:) = sist(:,:)             ! Ice surface temperature                      
(  185) 
(  186) 
(  187) # 200
(  200)       IF(ln_ctl) THEN
(  201)          CALL prt_ctl(tab2d_1=fsolar, clinfo1=' lim_flx: fsolar : ', tab2d_2=fnsolar, clinfo2=' fnsolar : ')
(  202)          CALL prt_ctl(tab2d_1=fmass , clinfo1=' lim_flx: fmass  : ', tab2d_2=fsalt  , clinfo2=' fsalt   : ')
(  203)          CALL prt_ctl(tab2d_1=ftaux , clinfo1=' lim_flx: ftaux  : ', tab2d_2=ftauy  , clinfo2=' ftauy   : ')
(  204)          CALL prt_ctl(tab2d_1=freeze, clinfo1=' lim_flx: freeze : ', tab2d_2=tn_ice , clinfo2=' tn_ice  : ')
(  205)       ENDIF 
(  206)    
(  207)     END SUBROUTINE lim_flx
(  208) 
(  209) 
(  210) # 218
(  218)    !!======================================================================
(  219) END MODULE limflx
