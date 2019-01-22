


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

Filename: albedo.F90

(    1) MODULE albedo
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  albedo  ***
(    4)    !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
(    5)    !!=====================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!   flx_blk_albedo : albedo for ocean and ice (clear and overcast skies)
(    8)    !!----------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE cpl_oce         ! ???
(   13)    USE phycst          ! physical constants
(   14)    USE daymod
(   15)    USE blk_oce         ! bulk variables
(   16)    USE flx_oce         ! forcings variables
(   17)    USE ocfzpt          ! ???
(   18)    USE in_out_manager
(   19)    USE lbclnk
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC flx_blk_albedo ! routine called by limflx.F90 in coupled
(   26)                          ! and in flxblk.F90 in forced
(   27)    !! * Module variables
(   28)    INTEGER  ::             &  !: nameos : ocean physical parameters
(   29)       albd_init = 0           !: control flag for initialization
(   30) 
(   31)    REAL(wp)  ::            &  ! constant values
(   32)       zzero   = 0.e0    ,  &
(   33)       zone    = 1.0
(   34) 
(   35)    !! * constants for albedo computation (flx_blk_albedo)
(   36)    REAL(wp) ::   &
(   37)       c1     = 0.05  ,     &   ! constants values
(   38)       c2     = 0.10  ,     &
(   39)       albice = 0.50  ,     &   !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
(   40)       cgren  = 0.06  ,     &   !  correction of the snow or ice albedo to take into account
(   41)                                !  effects of cloudiness (Grenfell & Perovich, 1984)
(   42)       alphd  = 0.80  ,     &   !  coefficients for linear interpolation used to compute
(   43)       alphdi = 0.72  ,     &   !  albedo between two extremes values (Pyane, 1972)
(   44)       alphc  = 0.65  ,     &
(   45)       zmue   = 0.40            !  cosine of local solar altitude
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 2

(   47)    !!----------------------------------------------------------------------
(   48)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/albedo.F90,v 1.3 2005/03/27 18:35:12 opalod Exp $ 
(   50)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   51)    !!----------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55) #if defined key_ice_lim
(   56)    !!----------------------------------------------------------------------
(   57)    !!   'key_ice_lim'                                         LIM ice model
(   58)    !!----------------------------------------------------------------------
(   59) 
(   60)    SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
(   61)       !!----------------------------------------------------------------------
(   62)       !!               ***  ROUTINE flx_blk_albedo  ***
(   63)       !!          
(   64)       !! ** Purpose :   Computation of the albedo of the snow/ice system 
(   65)       !!      as well as the ocean one
(   66)       !!       
(   67)       !! ** Method  : - Computation of the albedo of snow or ice (choose the 
(   68)       !!      rignt one by a large number of tests
(   69)       !!              - Computation of the albedo of the ocean
(   70)       !!
(   71)       !! References :
(   72)       !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
(   73)       !!
(   74)       !! History :
(   75)       !!  8.0   !  01-04  (LIM 1.0)
(   76)       !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Modules used
(   79)       USE ice                   ! ???
(   80) 
(   81)       !! * Arguments
(   82)       REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
(   83)          palb         ,     &    !  albedo of ice under overcast sky
(   84)          palcn        ,     &    !  albedo of ocean under overcast sky
(   85)          palbp        ,     &    !  albedo of ice under clear sky 
(   86)          palcnp                  !  albedo of ocean under clear sky
(   87) 
(   88)       !! * Local variables
(   89)       INTEGER ::    &
(   90)          ji, jj                   ! dummy loop indices
(   91)       REAL(wp) ::   & 
(   92)          zmue14         ,     &   !  zmue**1.4
(   93)          zalbpsnm       ,     &   !  albedo of ice under clear sky when snow is melting
(   94)          zalbpsnf       ,     &   !  albedo of ice under clear sky when snow is freezing
(   95)          zalbpsn        ,     &   !  albedo of snow/ice system when ice is coverd by snow
(   96)          zalbpic        ,     &   !  albedo of snow/ice system when ice is free of snow
(   97)          zithsn         ,     &   !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
(   98)          zitmlsn        ,     &   !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
(   99)          zihsc1         ,     &   !  = 1 hsn <= c1 ; = 0 hsn > c1
(  100)          zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2
(  101)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(  102)          zalbfz         ,     &   !  ( = alphdi for freezing ice ; = albice for melting ice )
(  103)          zficeth                  !  function of ice thickness
(  104)       LOGICAL , DIMENSION(jpi,jpj) ::  &






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 3

(  105)          llmask
(  106)       !!---------------------------------------------------------------------
(  107)       
(  108)       ! initialization 
(  109)       IF( albd_init == 0 )   CALL albedo_init
(  110) 
(  111)       !-------------------------                                                             
(  112)       !  Computation of  zficeth
(  113)       !-------------------------- 
(  114)       
(  115)       llmask = (hsnif == 0.e0) .AND. ( sist >= rt0_ice )
(  116)       WHERE ( llmask )   !  ice free of snow and melts
(  117)          zalbfz = albice
(  118)       ELSEWHERE                   
(  119)          zalbfz = alphdi
(  120)       END WHERE
(  121)       
(  122)       DO jj = 1, jpj
(  123)          DO ji = 1, jpi
(  124)             IF( hicif(ji,jj) > 1.5 ) THEN
(  125)                zficeth(ji,jj) = zalbfz(ji,jj)
(  126)             ELSEIF( hicif(ji,jj) > 1.0  .AND. hicif(ji,jj) <= 1.5 ) THEN
(  127)                zficeth(ji,jj) = 0.472 + 2.0 * ( zalbfz(ji,jj) - 0.472 ) * ( hicif(ji,jj) - 1.0 )
(  128)             ELSEIF( hicif(ji,jj) > 0.05 .AND. hicif(ji,jj) <= 1.0 ) THEN
(  129)                zficeth(ji,jj) = 0.2467 + 0.7049 * hicif(ji,jj)                                &
(  130)                   &                    - 0.8608 * hicif(ji,jj) * hicif(ji,jj)                 &
(  131)                   &                    + 0.3812 * hicif(ji,jj) * hicif(ji,jj) * hicif (ji,jj)
(  132)             ELSE
(  133)                zficeth(ji,jj) = 0.1 + 3.6 * hicif(ji,jj) 
(  134)             ENDIF
(  135)          END DO
(  136)       END DO
(  137)       
(  138)       !----------------------------------------------- 
(  139)       !    Computation of the snow/ice albedo system 
(  140)       !-------------------------- ---------------------
(  141)       
(  142)       !    Albedo of snow-ice for clear sky.
(  143)       !-----------------------------------------------    
(  144)       DO jj = 1, jpj
(  145)          DO ji = 1, jpi
(  146)             !  Case of ice covered by snow.             
(  147)             
(  148)             !  melting snow        
(  149)             zihsc1       = 1.0 - MAX ( zzero , SIGN ( zone , - ( hsnif(ji,jj) - c1 ) ) )
(  150)             zalbpsnm     = ( 1.0 - zihsc1 ) * ( zficeth(ji,jj) + hsnif(ji,jj) * ( alphd - zficeth(ji,jj) ) / c1 ) &
(  151)                &                 + zihsc1   * alphd  
(  152)             !  freezing snow                
(  153)             zihsc2       = MAX ( zzero , SIGN ( zone , hsnif(ji,jj) - c2 ) )
(  154)             zalbpsnf     = ( 1.0 - zihsc2 ) * ( albice + hsnif(ji,jj) * ( alphc - albice ) / c2 )                 &
(  155)                &                 + zihsc2   * alphc 
(  156)             
(  157)             zitmlsn      =  MAX ( zzero , SIGN ( zone , sist(ji,jj) - rt0_snow ) )   
(  158)             zalbpsn      =  zitmlsn * zalbpsnf + ( 1.0 - zitmlsn ) * zalbpsnm 
(  159)             
(  160)             !  Case of ice free of snow.
(  161)             zalbpic      = zficeth(ji,jj) 
(  162)             






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 4

(  163)             ! albedo of the system   
(  164)             zithsn       = 1.0 - MAX ( zzero , SIGN ( zone , - hsnif(ji,jj) ) )
(  165)             palbp(ji,jj) =  zithsn * zalbpsn + ( 1.0 - zithsn ) *  zalbpic
(  166)          END DO
(  167)       END DO
(  168)       
(  169)       !    Albedo of snow-ice for overcast sky.
(  170)       !----------------------------------------------  
(  171)       palb(:,:)   = palbp(:,:) + cgren                                           
(  172)       
(  173)       !--------------------------------------------
(  174)       !    Computation of the albedo of the ocean 
(  175)       !-------------------------- -----------------                                                          
(  176)       
(  177)       !  Parameterization of Briegled and Ramanathan, 1982 
(  178)       zmue14      = zmue**1.4                                       
(  179)       palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )                
(  180)       
(  181)       !  Parameterization of Kondratyev, 1969 and Payne, 1972
(  182)       palcn(:,:)  = 0.06                                                 
(  183)       
(  184)    END SUBROUTINE flx_blk_albedo
(  185) 
(  186) # else
(  187)    !!----------------------------------------------------------------------
(  188)    !!   Default option :                                   NO sea-ice model
(  189)    !!----------------------------------------------------------------------
(  190) 
(  191)    SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
(  192)       !!----------------------------------------------------------------------
(  193)       !!               ***  ROUTINE flx_blk_albedo  ***
(  194)       !! 
(  195)       !! ** Purpose :   Computation of the albedo of the snow/ice system
(  196)       !!      as well as the ocean one
(  197)       !!
(  198)       !! ** Method  :   Computation of the albedo of snow or ice (choose the
(  199)       !!      wright one by a large number of tests Computation of the albedo
(  200)       !!      of the ocean
(  201)       !!
(  202)       !! History :
(  203)       !!  8.0   !  01-04  (LIM 1.0)
(  204)       !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
(  205)       !!----------------------------------------------------------------------
(  206)       !! * Arguments
(  207)       REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
(  208)          palb         ,     &    !  albedo of ice under overcast sky
(  209)          palcn        ,     &    !  albedo of ocean under overcast sky
(  210)          palbp        ,     &    !  albedo of ice under clear sky
(  211)          palcnp                  !  albedo of ocean under clear sky
(  212) 
(  213)       REAL(wp) ::   &
(  214)          zmue14                 !  zmue**1.4
(  215)       !!----------------------------------------------------------------------
(  216) 
(  217)       !--------------------------------------------
(  218)       !    Computation of the albedo of the ocean
(  219)       !-------------------------- -----------------
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 5

(  221)       !  Parameterization of Briegled and Ramanathan, 1982
(  222)       zmue14      = zmue**1.4
(  223)       palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )
(  224) 
(  225)       !  Parameterization of Kondratyev, 1969 and Payne, 1972
(  226)       palcn(:,:)  = 0.06
(  227) 
(  228)       palb (:,:)  = palcn(:,:)
(  229)       palbp(:,:)  = palcnp(:,:)
(  230) 
(  231)    END SUBROUTINE flx_blk_albedo
(  232) 
(  233) #endif
(  234) 
(  235)    SUBROUTINE albedo_init
(  236)       !!----------------------------------------------------------------------
(  237)       !!                 ***  ROUTINE albedo_init  ***
(  238)       !!
(  239)       !! ** Purpose :   initializations for the albedo parameters
(  240)       !!
(  241)       !! ** Method  :   Read the namelist namalb
(  242)       !!
(  243)       !! ** Action  :  
(  244)       !!
(  245)       !!
(  246)       !! History :
(  247)       !!   9.0  !  04-11  (C. Talandier)  Original code
(  248)       !!----------------------------------------------------------------------
(  249)       NAMELIST/namalb/ cgren, albice, alphd, alphdi, alphc
(  250)       !!----------------------------------------------------------------------
(  251)       !!  OPA 9.0, LODYC-IPSL (2004)
(  252)       !!----------------------------------------------------------------------
(  253) 
(  254)       ! set the initialization flag to 1
(  255)       albd_init = 1           ! indicate that the initialization has been done
(  256) 
(  257)       ! Read Namelist namalb : albedo parameters
(  258)       REWIND( numnam )
(  259)       READ  ( numnam, namalb )
(  260) 
(  261)       ! Control print
(  262)       IF(lwp) THEN
(  263)          WRITE(numout,*)
(  264)          WRITE(numout,*) 'albedo_init : albedo '
(  265)          WRITE(numout,*) '~~~~~~~~~~~'
(  266)          WRITE(numout,*) '          Namelist namalb : set albedo parameters'
(  267)          WRITE(numout,*)
(  268)          WRITE(numout,*) '             correction of the snow or ice albedo to take into account cgren = ', cgren
(  269)          WRITE(numout,*) '             albedo of melting ice in the arctic and antarctic        albice = ', albice
(  270)          WRITE(numout,*) '             coefficients for linear                                   alphd = ', alphd
(  271)          WRITE(numout,*) '             interpolation used to compute albedo                     alphdi = ', alphdi
(  272)          WRITE(numout,*) '             between two extremes values (Pyane, 1972)                 alphc = ', alphc
(  273)          WRITE(numout,*)
(  274)       ENDIF
(  275) 
(  276)    END SUBROUTINE albedo_init
(  277)    !!======================================================================
(  278) END MODULE albedo






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 6

(    1) # 1 "albedo.F90"
(    1) MODULE albedo
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  albedo  ***
(    4)    !! Ocean forcing:  bulk thermohaline forcing of the ocean (or ice)
(    5)    !!=====================================================================
(    6)    !!----------------------------------------------------------------------
(    7)    !!   flx_blk_albedo : albedo for ocean and ice (clear and overcast skies)
(    8)    !!----------------------------------------------------------------------
(    9)    !! * Modules used
(   10)    USE oce             ! ocean dynamics and tracers
(   11)    USE dom_oce         ! ocean space and time domain
(   12)    USE cpl_oce         ! ???
(   13)    USE phycst          ! physical constants
(   14)    USE daymod
(   15)    USE blk_oce         ! bulk variables
(   16)    USE flx_oce         ! forcings variables
(   17)    USE ocfzpt          ! ???
(   18)    USE in_out_manager
(   19)    USE lbclnk
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Accessibility
(   25)    PUBLIC flx_blk_albedo ! routine called by limflx.F90 in coupled
(   26)                          ! and in flxblk.F90 in forced
(   27)    !! * Module variables
(   28)    INTEGER  ::             &  !: nameos : ocean physical parameters
(   29)       albd_init = 0           !: control flag for initialization
(   30) 
(   31)    REAL(wp)  ::            &  ! constant values
(   32)       zzero   = 0.e0    ,  &
(   33)       zone    = 1.0
(   34) 
(   35)    !! * constants for albedo computation (flx_blk_albedo)
(   36)    REAL(wp) ::   &
(   37)       c1     = 0.05  ,     &   ! constants values
(   38)       c2     = 0.10  ,     &
(   39)       albice = 0.50  ,     &   !  albedo of melting ice in the arctic and antarctic (Shine & Hendersson-Sellers)
(   40)       cgren  = 0.06  ,     &   !  correction of the snow or ice albedo to take into account
(   41)                                !  effects of cloudiness (Grenfell & Perovich, 1984)
(   42)       alphd  = 0.80  ,     &   !  coefficients for linear interpolation used to compute
(   43)       alphdi = 0.72  ,     &   !  albedo between two extremes values (Pyane, 1972)
(   44)       alphc  = 0.65  ,     &
(   45)       zmue   = 0.40            !  cosine of local solar altitude
(   46) 
(   47)    !!----------------------------------------------------------------------
(   48)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   49)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/albedo.F90,v 1.3 2005/03/27 18:35:12 opalod Exp $ 
(   50)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   51)    !!----------------------------------------------------------------------
(   52) 
(   53) CONTAINS
(   54) 
(   55) # 56
(   56)    !!----------------------------------------------------------------------
(   57)    !!   'key_ice_lim'                                         LIM ice model






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 7

(   58)    !!----------------------------------------------------------------------
(   59) 
(   60)    SUBROUTINE flx_blk_albedo( palb , palcn , palbp , palcnp )
(   61)       !!----------------------------------------------------------------------
(   62)       !!               ***  ROUTINE flx_blk_albedo  ***
(   63)       !!          
(   64)       !! ** Purpose :   Computation of the albedo of the snow/ice system 
(   65)       !!      as well as the ocean one
(   66)       !!       
(   67)       !! ** Method  : - Computation of the albedo of snow or ice (choose the 
(   68)       !!      rignt one by a large number of tests
(   69)       !!              - Computation of the albedo of the ocean
(   70)       !!
(   71)       !! References :
(   72)       !!      Shine and Hendersson-Sellers 1985, JGR, 90(D1), 2243-2250.
(   73)       !!
(   74)       !! History :
(   75)       !!  8.0   !  01-04  (LIM 1.0)
(   76)       !!  8.5   !  03-07  (C. Ethe, G. Madec)  Optimization (old name:shine)
(   77)       !!----------------------------------------------------------------------
(   78)       !! * Modules used
(   79)       USE ice                   ! ???
(   80) 
(   81)       !! * Arguments
(   82)       REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::  &
(   83)          palb         ,     &    !  albedo of ice under overcast sky
(   84)          palcn        ,     &    !  albedo of ocean under overcast sky
(   85)          palbp        ,     &    !  albedo of ice under clear sky 
(   86)          palcnp                  !  albedo of ocean under clear sky
(   87) 
(   88)       !! * Local variables
(   89)       INTEGER ::    &
(   90)          ji, jj                   ! dummy loop indices
(   91)       REAL(wp) ::   & 
(   92)          zmue14         ,     &   !  zmue**1.4
(   93)          zalbpsnm       ,     &   !  albedo of ice under clear sky when snow is melting
(   94)          zalbpsnf       ,     &   !  albedo of ice under clear sky when snow is freezing
(   95)          zalbpsn        ,     &   !  albedo of snow/ice system when ice is coverd by snow
(   96)          zalbpic        ,     &   !  albedo of snow/ice system when ice is free of snow
(   97)          zithsn         ,     &   !  = 1 for hsn >= 0 ( ice is cov. by snow ) ; = 0 otherwise (ice is free of snow)
(   98)          zitmlsn        ,     &   !  = 1 freezinz snow (sist >=rt0_snow) ; = 0 melting snow (sist<rt0_snow)
(   99)          zihsc1         ,     &   !  = 1 hsn <= c1 ; = 0 hsn > c1
(  100)          zihsc2                   !  = 1 hsn >= c2 ; = 0 hsn < c2
(  101)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(  102)          zalbfz         ,     &   !  ( = alphdi for freezing ice ; = albice for melting ice )
(  103)          zficeth                  !  function of ice thickness
(  104)       LOGICAL , DIMENSION(jpi,jpj) ::  &
(  105)          llmask
(  106)       !!---------------------------------------------------------------------
(  107)       
(  108)       ! initialization 
(  109)       IF( albd_init == 0 )   CALL albedo_init
(  110) 
(  111)       !-------------------------                                                             
(  112)       !  Computation of  zficeth
(  113)       !-------------------------- 
(  114)       
(  115)       llmask = (hsnif == 0.e0) .AND. ( sist >= rt0_ice )






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 8

(  116)       WHERE ( llmask )   !  ice free of snow and melts
(  117)          zalbfz = albice
(  118)       ELSEWHERE                   
(  119)          zalbfz = alphdi
(  120)       END WHERE
(  121)       
(  122)       DO jj = 1, jpj
(  123)          DO ji = 1, jpi
(  124)             IF( hicif(ji,jj) > 1.5 ) THEN
(  125)                zficeth(ji,jj) = zalbfz(ji,jj)
(  126)             ELSEIF( hicif(ji,jj) > 1.0  .AND. hicif(ji,jj) <= 1.5 ) THEN
(  127)                zficeth(ji,jj) = 0.472 + 2.0 * ( zalbfz(ji,jj) - 0.472 ) * ( hicif(ji,jj) - 1.0 )
(  128)             ELSEIF( hicif(ji,jj) > 0.05 .AND. hicif(ji,jj) <= 1.0 ) THEN
(  129)                zficeth(ji,jj) = 0.2467 + 0.7049 * hicif(ji,jj)                                &
(  130)                   &                    - 0.8608 * hicif(ji,jj) * hicif(ji,jj)                 &
(  131)                   &                    + 0.3812 * hicif(ji,jj) * hicif(ji,jj) * hicif (ji,jj)
(  132)             ELSE
(  133)                zficeth(ji,jj) = 0.1 + 3.6 * hicif(ji,jj) 
(  134)             ENDIF
(  135)          END DO
(  136)       END DO
(  137)       
(  138)       !----------------------------------------------- 
(  139)       !    Computation of the snow/ice albedo system 
(  140)       !-------------------------- ---------------------
(  141)       
(  142)       !    Albedo of snow-ice for clear sky.
(  143)       !-----------------------------------------------    
(  144)       DO jj = 1, jpj
(  145)          DO ji = 1, jpi
(  146)             !  Case of ice covered by snow.             
(  147)             
(  148)             !  melting snow        
(  149)             zihsc1       = 1.0 - MAX ( zzero , SIGN ( zone , - ( hsnif(ji,jj) - c1 ) ) )
(  150)             zalbpsnm     = ( 1.0 - zihsc1 ) * ( zficeth(ji,jj) + hsnif(ji,jj) * ( alphd - zficeth(ji,jj) ) / c1 ) &
(  151)                &                 + zihsc1   * alphd  
(  152)             !  freezing snow                
(  153)             zihsc2       = MAX ( zzero , SIGN ( zone , hsnif(ji,jj) - c2 ) )
(  154)             zalbpsnf     = ( 1.0 - zihsc2 ) * ( albice + hsnif(ji,jj) * ( alphc - albice ) / c2 )                 &
(  155)                &                 + zihsc2   * alphc 
(  156)             
(  157)             zitmlsn      =  MAX ( zzero , SIGN ( zone , sist(ji,jj) - rt0_snow ) )   
(  158)             zalbpsn      =  zitmlsn * zalbpsnf + ( 1.0 - zitmlsn ) * zalbpsnm 
(  159)             
(  160)             !  Case of ice free of snow.
(  161)             zalbpic      = zficeth(ji,jj) 
(  162)             
(  163)             ! albedo of the system   
(  164)             zithsn       = 1.0 - MAX ( zzero , SIGN ( zone , - hsnif(ji,jj) ) )
(  165)             palbp(ji,jj) =  zithsn * zalbpsn + ( 1.0 - zithsn ) *  zalbpic
(  166)          END DO
(  167)       END DO
(  168)       
(  169)       !    Albedo of snow-ice for overcast sky.
(  170)       !----------------------------------------------  
(  171)       palb(:,:)   = palbp(:,:) + cgren                                           
(  172)       
(  173)       !--------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 9

(  174)       !    Computation of the albedo of the ocean 
(  175)       !-------------------------- -----------------                                                          
(  176)       
(  177)       !  Parameterization of Briegled and Ramanathan, 1982 
(  178)       zmue14      = zmue**1.4                                       
(  179)       palcnp(:,:) = 0.05 / ( 1.1 * zmue14 + 0.15 )                
(  180)       
(  181)       !  Parameterization of Kondratyev, 1969 and Payne, 1972
(  182)       palcn(:,:)  = 0.06                                                 
(  183)       
(  184)    END SUBROUTINE flx_blk_albedo
(  185) 
(  186) 
(  187) # 235
(  235)    SUBROUTINE albedo_init
(  236)       !!----------------------------------------------------------------------
(  237)       !!                 ***  ROUTINE albedo_init  ***
(  238)       !!
(  239)       !! ** Purpose :   initializations for the albedo parameters
(  240)       !!
(  241)       !! ** Method  :   Read the namelist namalb
(  242)       !!
(  243)       !! ** Action  :  
(  244)       !!
(  245)       !!
(  246)       !! History :
(  247)       !!   9.0  !  04-11  (C. Talandier)  Original code
(  248)       !!----------------------------------------------------------------------
(  249)       NAMELIST/namalb/ cgren, albice, alphd, alphdi, alphc
(  250)       !!----------------------------------------------------------------------
(  251)       !!  OPA 9.0, LODYC-IPSL (2004)
(  252)       !!----------------------------------------------------------------------
(  253) 
(  254)       ! set the initialization flag to 1
(  255)       albd_init = 1           ! indicate that the initialization has been done
(  256) 
(  257)       ! Read Namelist namalb : albedo parameters
(  258)       REWIND( numnam )
(  259)       READ  ( numnam, namalb )
(  260) 
(  261)       ! Control print
(  262)       IF(lwp) THEN
(  263)          WRITE(numout,*)
(  264)          WRITE(numout,*) 'albedo_init : albedo '
(  265)          WRITE(numout,*) '~~~~~~~~~~~'
(  266)          WRITE(numout,*) '          Namelist namalb : set albedo parameters'
(  267)          WRITE(numout,*)
(  268)          WRITE(numout,*) '             correction of the snow or ice albedo to take into account cgren = ', cgren
(  269)          WRITE(numout,*) '             albedo of melting ice in the arctic and antarctic        albice = ', albice
(  270)          WRITE(numout,*) '             coefficients for linear                                   alphd = ', alphd
(  271)          WRITE(numout,*) '             interpolation used to compute albedo                     alphdi = ', alphdi
(  272)          WRITE(numout,*) '             between two extremes values (Pyane, 1972)                 alphc = ', alphc
(  273)          WRITE(numout,*)
(  274)       ENDIF
(  275) 
(  276)    END SUBROUTINE albedo_init
(  277)    !!======================================================================
(  278) END MODULE albedo



