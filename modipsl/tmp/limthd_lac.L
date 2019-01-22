


PGF90 (Version     10.5)          09/23/2010  09:49:56      page 1

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

Filename: limthd_lac.F90

(    1) MODULE limthd_lac
(    2) #if defined key_ice_lim
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE limthd_lac   ***
(    5)    !!                lateral thermodynamic growth of the ice 
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   lim_lat_acr    : lateral accretion of ice
(   10)    !! * Modules used
(   11)    USE par_oce          ! ocean parameters
(   12)    USE phycst
(   13)    USE ice_oce         ! ice variables
(   14)    USE thd_ice
(   15)    USE iceini
(   16)    USE limistate 
(   17)      
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC lim_thd_lac     ! called by lim_thd
(   23) 
(   24)    !! * Module variables
(   25)    REAL(wp)  ::           &  ! constant values
(   26)       epsi20 = 1.e-20  ,  &
(   27)       epsi13 = 1.e-13  ,  &
(   28)       zzero  = 0.e0    ,  &
(   29)       zone   = 1.e0
(   30)    !!----------------------------------------------------------------------
(   31)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd_lac.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) CONTAINS
(   36)     
(   37)    SUBROUTINE lim_thd_lac( kideb, kiut )
(   38)       !!-------------------------------------------------------------------
(   39)       !!               ***   ROUTINE lim_thd_lac  ***
(   40)       !!  
(   41)       !! ** Purpose : Computation of the evolution of the ice thickness and 
(   42)       !!      concentration as a function of the heat balance in the leads.
(   43)       !!      It is only used for lateral accretion
(   44)       !!       
(   45)       !! ** Method  : Ice is formed in the open water when ocean lose heat
(   46)       !!      (heat budget of open water Bl is negative) .






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 2

(   47)       !!      Computation of the increase of 1-A (ice concentration) fol-
(   48)       !!      lowing the law :
(   49)       !!      (dA/dt)acc = F[ (1-A)/(1-a) ] * [ Bl / (Li*h0) ]
(   50)       !!       where - h0 is the thickness of ice created in the lead
(   51)       !!             - a is a minimum fraction for leads
(   52)       !!             - F is a monotonic non-increasing function defined as:
(   53)       !!                  F(X)=( 1 - X**exld )**(1.0/exld)
(   54)       !!             - exld is the exponent closure rate (=2 default val.)
(   55)       !! 
(   56)       !! ** Action : - Adjustment of snow and ice thicknesses and heat
(   57)       !!                content in brine pockets
(   58)       !!             - Updating ice internal temperature
(   59)       !!             - Computation of variation of ice volume and mass
(   60)       !!             - Computation of frldb after lateral accretion and 
(   61)       !!               update h_snow_1d, h_ice_1d and tbif_1d(:,:)      
(   62)       !! 
(   63)       !! ** References :
(   64)       !!      M. Maqueda, 1995, PhD Thesis, Univesidad Complutense Madrid
(   65)       !!      Fichefet T. and M. Maqueda 1997, J. Geo. Res., 102(C6), 
(   66)       !!                                                12609 -12646   
(   67)       !! History :
(   68)       !!   1.0  !  01-04 (LIM)  original code
(   69)       !!   2.0  !  02-08 (C. Ethe, G. Madec)  F90, mpp
(   70)       !!-------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER , INTENT(IN)::  &
(   73)          kideb          ,   &  ! start point on which the the computation is applied
(   74)          kiut                  ! end point on which the the computation is applied
(   75) 
(   76)       !! * Local variables
(   77)       INTEGER ::            &
(   78)          ji             ,   &  !  dummy loop indices
(   79)          iicefr         ,   &  !  1 = existing ice ; 0 = no ice
(   80)          iiceform       ,   &  !  1 = ice formed   ; 0 = no ice formed
(   81)          ihemis                !  dummy indice
(   82)       REAL(wp), DIMENSION(jpij) :: &
(   83)          zqbgow           ,  &  !  heat budget of the open water (negative)
(   84)          zfrl_old         ,  &  !  previous sea/ice fraction
(   85)          zhice_old        ,  &  !  previous ice thickness
(   86)          zhice0           ,  &  !  thickness of newly formed ice in leads
(   87)          zfrlmin          ,  &  !  minimum fraction for leads
(   88)          zdhicbot               !  part of thickness of newly formed ice in leads which 
(   89)                                 !  has been already used in transport for example
(   90)       REAL(wp)  ::  &
(   91)          zhemis           ,  &  !  hemisphere (0 = North, 1 = South)
(   92)          zhicenew         ,  &  !  new ice thickness
(   93)          zholds2          ,  &  !  ratio of previous ice thickness and 2 
(   94)          zhnews2          ,  &  !  ratio of new ice thickness and 2 
(   95)          zfrlnew          ,  &  !  new sea/ice fraction
(   96)          zfrld            ,  &  !  ratio of sea/ice fraction and minimum fraction for leads
(   97)          zfrrate          ,  &  !  leads-closure rate
(   98)          zdfrl                  !  sea-ice fraction increment
(   99)       REAL(wp)  ::  &
(  100)          zdh1 , zdh2 , zdh3 , zdh4, zdh5   , &   ! tempory scalars
(  101)          ztint , zta1 , zta2 , zta3 , zta4 , &
(  102)          zah, zalpha , zbeta
(  103)       !!---------------------------------------------------------------------      
(  104)                    






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 3

(  105)       !--------------------------------------------------------------
(  106)       !   Computation of the heat budget of the open water (negative)
(  107)       !--------------------------------------------------------------
(  108)       
(  109)       DO ji = kideb , kiut      
(  110)          zqbgow(ji) = qldif_1d(ji) - qcmif_1d(ji)
(  111)       END DO
(  112)       
(  113)       !-----------------------------------------------------------------
(  114)       !   Taking the appropriate values for the corresponding hemisphere
(  115)       !-----------------------------------------------------------------
(  116)       DO ji = kideb , kiut
(  117)          zhemis       = MAX( zzero , SIGN( zone , frld_1d(ji) - 2.0 ) ) 
(  118)          ihemis       = INT( 1 + zhemis )
(  119)          zhice0  (ji) = hiccrit( ihemis ) 
(  120)          zfrlmin (ji) = acrit  ( ihemis )   
(  121)          frld_1d (ji) = frld_1d(ji) - 2.0 * zhemis
(  122)          zfrl_old(ji) = frld_1d(ji)
(  123)       END DO
(  124)       
(  125)       !-------------------------------------------------------------------
(  126)       !     Lateral Accretion (modification of the fraction of open water)
(  127)       !     The ice formed in the leads has always a thickness zhice0, but
(  128)       !     only a fraction zfrrate of the ice formed contributes to the 
(  129)       !     increase of the ice fraction. The remaining part (1-zfrrate)
(  130)       !     is rather assumed to lead to an increase in the thickness of the
(  131)       !     pre-existing ice (transport for example). 
(  132)       !     Morales Maqueda, 1995 - Fichefet and Morales Maqueda, 1997
(  133)       !---------------------------------------------------------------------
(  134)       
(  135)       !CDIR NOVERRCHK
(  136)       DO ji = kideb , kiut
(  137)          iicefr       = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  138)          !---computation of the leads-closure rate
(  139)          zfrld        = MIN( zone , ( 1.0 - frld_1d(ji) ) / ( 1.0 - zfrlmin(ji) ) )
(  140)          zfrrate      = ( 1.0 - zfrld**exld )**( 1.0 / exld )
(  141)          !--computation of the sea-ice fraction increment and the new fraction 
(  142)          zdfrl        = ( zfrrate / zhice0(ji) )  * ( zqbgow(ji) / xlic )
(  143)          zfrlnew      = zfrl_old(ji) + zdfrl
(  144)          !--update the sea-ice fraction 
(  145)          frld_1d   (ji) = MAX( zfrlnew , zfrlmin(ji) )
(  146)          !--computation of the remaining part of ice thickness which has been already used
(  147)          zdhicbot(ji) =  ( frld_1d(ji) - zfrlnew ) * zhice0(ji) / ( 1.0 - zfrlmin(ji) ) & 
(  148)                       -  (  ( 1.0 - zfrrate ) / ( 1.0 - frld_1d(ji) ) )  * ( zqbgow(ji) / xlic ) 
(  149)       END DO
(  150)  
(  151)       !----------------------------------------------------------------------------------------
(  152)       !      Ajustement of snow and ice thicknesses and updating the total heat stored in brine pockets  
(  153)       !      The thickness of newly formed ice is averaged with that of the pre-existing
(  154)       !         (1-Anew) * hinew = (1-Aold) * hiold + ((1-Anew)-(1-Aold)) * h0
(  155)       !      Snow is distributed over the new ice-covered area 
(  156)       !         (1-Anew) * hsnew = (1-Aold) * hsold            
(  157)       !--------------------------------------------------------------------------------------------
(  158)       
(  159)       DO ji = kideb , kiut
(  160)          iicefr       = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  161)          zhice_old(ji) = h_ice_1d(ji)
(  162)          zhicenew      = iicefr * zhice_old(ji) + ( 1 - iicefr ) * zhice0(ji)






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 4

(  163)          zalpha        = ( 1. - zfrl_old(ji) ) / ( 1.- frld_1d(ji) )
(  164)          h_snow_1d(ji) = zalpha * h_snow_1d(ji)
(  165)          h_ice_1d (ji) = zalpha * zhicenew + ( 1.0 - zalpha ) * zhice0(ji)
(  166)          qstbif_1d(ji) = zalpha * qstbif_1d(ji) 
(  167)       END DO
(  168)       
(  169)       !-------------------------------------------------------
(  170)       !   Ajustement of ice internal temperatures
(  171)       !-------------------------------------------------------
(  172)       
(  173)       DO ji = kideb , kiut
(  174)          iicefr      = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  175)          iiceform    = 1 - MAX( 0 ,INT( SIGN( 1.5 * zone , zhice0(ji) - h_ice_1d(ji) ) ) )
(  176)          zholds2     = zhice_old(ji)/ 2.
(  177)          zhnews2     = h_ice_1d(ji) / 2.
(  178)          zdh1        = MAX( zzero ,  zhice_old(ji)   - zhnews2 )
(  179)          zdh2        = MAX( zzero , -zhice_old(ji)   + zhnews2 )
(  180)          zdh3        = MAX( zzero ,  h_ice_1d(ji) - zholds2 )
(  181)          zdh4        = MAX( zzero , -h_ice_1d(ji) + zholds2 )
(  182)          zdh5        = MAX( zzero , zhice0(ji)      - zholds2 )
(  183)          ztint       =       iiceform   * (  ( zholds2 - zdh3 ) * tbif_1d(ji,3) + zdh4 * tbif_1d(ji,2) )      &
(  184)             &                           / MAX( epsi20 , h_ice_1d(ji) - zhice0(ji) )                           &
(  185)             &                 + ( 1 - iiceform ) * tfu_1d(ji)
(  186)          zta1        = iicefr * ( 1.  - zfrl_old(ji) ) * tbif_1d(ji,2) 
(  187)          zta2        = iicefr * ( 1.  - zfrl_old(ji) ) * tbif_1d(ji,3)
(  188)          zta3        = iicefr * ( 1.  - zfrl_old(ji) ) * ztint
(  189)          zta4        = ( zfrl_old(ji) - frld_1d   (ji) ) * tfu_1d(ji)
(  190)          zah         = ( 1. - frld_1d(ji) ) * zhnews2 
(  191) 
(  192)          tbif_1d(ji,2) = (  MIN( zhnews2 , zholds2 )                                              * zta1   &
(  193)             &          + ( 1 - iiceform ) * ( zholds2 - zdh1 )                                    * zta2   &
(  194)             &          + ( iiceform * ( zhnews2 - zhice0(ji) + zdh5 ) + ( 1 - iiceform ) * zdh2 ) * zta3   & 
(  195)             &          + MIN ( zhnews2 , zhice0(ji) )                                             * zta4   &
(  196)             &          ) / zah
(  197)          
(  198)          tbif_1d(ji,3) =     (  iiceform * ( zhnews2 - zdh3 )                                          * zta1  &
(  199)             &              + ( iiceform * zdh3 + ( 1 - iiceform ) * zdh1 )                             * zta2  &
(  200)             &              + ( iiceform * ( zhnews2 - zdh5 ) + ( 1 - iiceform ) * ( zhnews2 - zdh1 ) ) * zta3  & 
(  201)             &              + ( iiceform * zdh5 + ( 1 - iiceform ) * zhnews2 )                          * zta4  &
(  202)             &            ) / zah
(  203)          !---removing the remaining part of ice formed which has been already used
(  204)          zbeta         = h_ice_1d(ji) / ( h_ice_1d(ji) + zdhicbot(ji) )
(  205)          h_ice_1d(ji)  = h_ice_1d(ji) + zdhicbot(ji)
(  206)          tbif_1d (ji,2)= zbeta * tbif_1d(ji,2) + ( 1.0 - zbeta ) * tbif_1d(ji,3)
(  207)          tbif_1d (ji,3)= ( 2. * zbeta - 1.0 ) * tbif_1d(ji,3) + ( 2. * zdhicbot(ji) / h_ice_1d(ji) ) * tfu_1d(ji)
(  208)          
(  209)       END DO
(  210)       
(  211)       !-------------------------------------------------------------
(  212)       !    Computation of variation of ice volume and ice mass 
(  213)       !           Vold = (1-Aold) * hiold ; Vnew = (1-Anew) * hinew
(  214)       !           dV = Vnew - Vold
(  215)       !-------------------------------------------------------------
(  216)       
(  217)       DO ji = kideb , kiut
(  218)          dvlbq_1d  (ji) = ( 1. - frld_1d(ji) ) * h_ice_1d(ji) - ( 1. - zfrl_old(ji) ) * zhice_old(ji)
(  219)          rdmicif_1d(ji) = rdmicif_1d(ji) + rhoic * dvlbq_1d(ji)
(  220)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 5

(  221)       
(  222)    END SUBROUTINE lim_thd_lac
(  223) #else
(  224)    !!======================================================================
(  225)    !!                       ***  MODULE limthd_lac   ***
(  226)    !!                           no sea ice model
(  227)    !!======================================================================
(  228) CONTAINS
(  229)    SUBROUTINE lim_thd_lac           ! Empty routine
(  230)    END SUBROUTINE lim_thd_lac
(  231) #endif
(  232) END MODULE limthd_lac




















































PGF90 (Version     10.5)          09/23/2010  09:49:56      page 6

(    1) # 1 "limthd_lac.F90"
(    1) MODULE limthd_lac
(    2) # 3
(    3)    !!======================================================================
(    4)    !!                       ***  MODULE limthd_lac   ***
(    5)    !!                lateral thermodynamic growth of the ice 
(    6)    !!======================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   lim_lat_acr    : lateral accretion of ice
(   10)    !! * Modules used
(   11)    USE par_oce          ! ocean parameters
(   12)    USE phycst
(   13)    USE ice_oce         ! ice variables
(   14)    USE thd_ice
(   15)    USE iceini
(   16)    USE limistate 
(   17)      
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC lim_thd_lac     ! called by lim_thd
(   23) 
(   24)    !! * Module variables
(   25)    REAL(wp)  ::           &  ! constant values
(   26)       epsi20 = 1.e-20  ,  &
(   27)       epsi13 = 1.e-13  ,  &
(   28)       zzero  = 0.e0    ,  &
(   29)       zone   = 1.e0
(   30)    !!----------------------------------------------------------------------
(   31)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd_lac.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) CONTAINS
(   36)     
(   37)    SUBROUTINE lim_thd_lac( kideb, kiut )
(   38)       !!-------------------------------------------------------------------
(   39)       !!               ***   ROUTINE lim_thd_lac  ***
(   40)       !!  
(   41)       !! ** Purpose : Computation of the evolution of the ice thickness and 
(   42)       !!      concentration as a function of the heat balance in the leads.
(   43)       !!      It is only used for lateral accretion
(   44)       !!       
(   45)       !! ** Method  : Ice is formed in the open water when ocean lose heat
(   46)       !!      (heat budget of open water Bl is negative) .
(   47)       !!      Computation of the increase of 1-A (ice concentration) fol-
(   48)       !!      lowing the law :
(   49)       !!      (dA/dt)acc = F[ (1-A)/(1-a) ] * [ Bl / (Li*h0) ]
(   50)       !!       where - h0 is the thickness of ice created in the lead
(   51)       !!             - a is a minimum fraction for leads
(   52)       !!             - F is a monotonic non-increasing function defined as:
(   53)       !!                  F(X)=( 1 - X**exld )**(1.0/exld)
(   54)       !!             - exld is the exponent closure rate (=2 default val.)
(   55)       !! 
(   56)       !! ** Action : - Adjustment of snow and ice thicknesses and heat
(   57)       !!                content in brine pockets






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 7

(   58)       !!             - Updating ice internal temperature
(   59)       !!             - Computation of variation of ice volume and mass
(   60)       !!             - Computation of frldb after lateral accretion and 
(   61)       !!               update h_snow_1d, h_ice_1d and tbif_1d(:,:)      
(   62)       !! 
(   63)       !! ** References :
(   64)       !!      M. Maqueda, 1995, PhD Thesis, Univesidad Complutense Madrid
(   65)       !!      Fichefet T. and M. Maqueda 1997, J. Geo. Res., 102(C6), 
(   66)       !!                                                12609 -12646   
(   67)       !! History :
(   68)       !!   1.0  !  01-04 (LIM)  original code
(   69)       !!   2.0  !  02-08 (C. Ethe, G. Madec)  F90, mpp
(   70)       !!-------------------------------------------------------------------
(   71)       !! * Arguments
(   72)       INTEGER , INTENT(IN)::  &
(   73)          kideb          ,   &  ! start point on which the the computation is applied
(   74)          kiut                  ! end point on which the the computation is applied
(   75) 
(   76)       !! * Local variables
(   77)       INTEGER ::            &
(   78)          ji             ,   &  !  dummy loop indices
(   79)          iicefr         ,   &  !  1 = existing ice ; 0 = no ice
(   80)          iiceform       ,   &  !  1 = ice formed   ; 0 = no ice formed
(   81)          ihemis                !  dummy indice
(   82)       REAL(wp), DIMENSION(jpij) :: &
(   83)          zqbgow           ,  &  !  heat budget of the open water (negative)
(   84)          zfrl_old         ,  &  !  previous sea/ice fraction
(   85)          zhice_old        ,  &  !  previous ice thickness
(   86)          zhice0           ,  &  !  thickness of newly formed ice in leads
(   87)          zfrlmin          ,  &  !  minimum fraction for leads
(   88)          zdhicbot               !  part of thickness of newly formed ice in leads which 
(   89)                                 !  has been already used in transport for example
(   90)       REAL(wp)  ::  &
(   91)          zhemis           ,  &  !  hemisphere (0 = North, 1 = South)
(   92)          zhicenew         ,  &  !  new ice thickness
(   93)          zholds2          ,  &  !  ratio of previous ice thickness and 2 
(   94)          zhnews2          ,  &  !  ratio of new ice thickness and 2 
(   95)          zfrlnew          ,  &  !  new sea/ice fraction
(   96)          zfrld            ,  &  !  ratio of sea/ice fraction and minimum fraction for leads
(   97)          zfrrate          ,  &  !  leads-closure rate
(   98)          zdfrl                  !  sea-ice fraction increment
(   99)       REAL(wp)  ::  &
(  100)          zdh1 , zdh2 , zdh3 , zdh4, zdh5   , &   ! tempory scalars
(  101)          ztint , zta1 , zta2 , zta3 , zta4 , &
(  102)          zah, zalpha , zbeta
(  103)       !!---------------------------------------------------------------------      
(  104)                    
(  105)       !--------------------------------------------------------------
(  106)       !   Computation of the heat budget of the open water (negative)
(  107)       !--------------------------------------------------------------
(  108)       
(  109)       DO ji = kideb , kiut      
(  110)          zqbgow(ji) = qldif_1d(ji) - qcmif_1d(ji)
(  111)       END DO
(  112)       
(  113)       !-----------------------------------------------------------------
(  114)       !   Taking the appropriate values for the corresponding hemisphere
(  115)       !-----------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 8

(  116)       DO ji = kideb , kiut
(  117)          zhemis       = MAX( zzero , SIGN( zone , frld_1d(ji) - 2.0 ) ) 
(  118)          ihemis       = INT( 1 + zhemis )
(  119)          zhice0  (ji) = hiccrit( ihemis ) 
(  120)          zfrlmin (ji) = acrit  ( ihemis )   
(  121)          frld_1d (ji) = frld_1d(ji) - 2.0 * zhemis
(  122)          zfrl_old(ji) = frld_1d(ji)
(  123)       END DO
(  124)       
(  125)       !-------------------------------------------------------------------
(  126)       !     Lateral Accretion (modification of the fraction of open water)
(  127)       !     The ice formed in the leads has always a thickness zhice0, but
(  128)       !     only a fraction zfrrate of the ice formed contributes to the 
(  129)       !     increase of the ice fraction. The remaining part (1-zfrrate)
(  130)       !     is rather assumed to lead to an increase in the thickness of the
(  131)       !     pre-existing ice (transport for example). 
(  132)       !     Morales Maqueda, 1995 - Fichefet and Morales Maqueda, 1997
(  133)       !---------------------------------------------------------------------
(  134)       
(  135)       !CDIR NOVERRCHK
(  136)       DO ji = kideb , kiut
(  137)          iicefr       = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  138)          !---computation of the leads-closure rate
(  139)          zfrld        = MIN( zone , ( 1.0 - frld_1d(ji) ) / ( 1.0 - zfrlmin(ji) ) )
(  140)          zfrrate      = ( 1.0 - zfrld**exld )**( 1.0 / exld )
(  141)          !--computation of the sea-ice fraction increment and the new fraction 
(  142)          zdfrl        = ( zfrrate / zhice0(ji) )  * ( zqbgow(ji) / xlic )
(  143)          zfrlnew      = zfrl_old(ji) + zdfrl
(  144)          !--update the sea-ice fraction 
(  145)          frld_1d   (ji) = MAX( zfrlnew , zfrlmin(ji) )
(  146)          !--computation of the remaining part of ice thickness which has been already used
(  147)          zdhicbot(ji) =  ( frld_1d(ji) - zfrlnew ) * zhice0(ji) / ( 1.0 - zfrlmin(ji) ) & 
(  148)                       -  (  ( 1.0 - zfrrate ) / ( 1.0 - frld_1d(ji) ) )  * ( zqbgow(ji) / xlic ) 
(  149)       END DO
(  150)  
(  151)       !----------------------------------------------------------------------------------------
(  152)       !      Ajustement of snow and ice thicknesses and updating the total heat stored in brine pockets  
(  153)       !      The thickness of newly formed ice is averaged with that of the pre-existing
(  154)       !         (1-Anew) * hinew = (1-Aold) * hiold + ((1-Anew)-(1-Aold)) * h0
(  155)       !      Snow is distributed over the new ice-covered area 
(  156)       !         (1-Anew) * hsnew = (1-Aold) * hsold            
(  157)       !--------------------------------------------------------------------------------------------
(  158)       
(  159)       DO ji = kideb , kiut
(  160)          iicefr       = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  161)          zhice_old(ji) = h_ice_1d(ji)
(  162)          zhicenew      = iicefr * zhice_old(ji) + ( 1 - iicefr ) * zhice0(ji)
(  163)          zalpha        = ( 1. - zfrl_old(ji) ) / ( 1.- frld_1d(ji) )
(  164)          h_snow_1d(ji) = zalpha * h_snow_1d(ji)
(  165)          h_ice_1d (ji) = zalpha * zhicenew + ( 1.0 - zalpha ) * zhice0(ji)
(  166)          qstbif_1d(ji) = zalpha * qstbif_1d(ji) 
(  167)       END DO
(  168)       
(  169)       !-------------------------------------------------------
(  170)       !   Ajustement of ice internal temperatures
(  171)       !-------------------------------------------------------
(  172)       
(  173)       DO ji = kideb , kiut






PGF90 (Version     10.5)          09/23/2010  09:49:56      page 9

(  174)          iicefr      = 1 - MAX( 0, INT( SIGN( 1.5 * zone , zfrl_old(ji) - 1.0 + epsi13 ) ) )
(  175)          iiceform    = 1 - MAX( 0 ,INT( SIGN( 1.5 * zone , zhice0(ji) - h_ice_1d(ji) ) ) )
(  176)          zholds2     = zhice_old(ji)/ 2.
(  177)          zhnews2     = h_ice_1d(ji) / 2.
(  178)          zdh1        = MAX( zzero ,  zhice_old(ji)   - zhnews2 )
(  179)          zdh2        = MAX( zzero , -zhice_old(ji)   + zhnews2 )
(  180)          zdh3        = MAX( zzero ,  h_ice_1d(ji) - zholds2 )
(  181)          zdh4        = MAX( zzero , -h_ice_1d(ji) + zholds2 )
(  182)          zdh5        = MAX( zzero , zhice0(ji)      - zholds2 )
(  183)          ztint       =       iiceform   * (  ( zholds2 - zdh3 ) * tbif_1d(ji,3) + zdh4 * tbif_1d(ji,2) )      &
(  184)             &                           / MAX( epsi20 , h_ice_1d(ji) - zhice0(ji) )                           &
(  185)             &                 + ( 1 - iiceform ) * tfu_1d(ji)
(  186)          zta1        = iicefr * ( 1.  - zfrl_old(ji) ) * tbif_1d(ji,2) 
(  187)          zta2        = iicefr * ( 1.  - zfrl_old(ji) ) * tbif_1d(ji,3)
(  188)          zta3        = iicefr * ( 1.  - zfrl_old(ji) ) * ztint
(  189)          zta4        = ( zfrl_old(ji) - frld_1d   (ji) ) * tfu_1d(ji)
(  190)          zah         = ( 1. - frld_1d(ji) ) * zhnews2 
(  191) 
(  192)          tbif_1d(ji,2) = (  MIN( zhnews2 , zholds2 )                                              * zta1   &
(  193)             &          + ( 1 - iiceform ) * ( zholds2 - zdh1 )                                    * zta2   &
(  194)             &          + ( iiceform * ( zhnews2 - zhice0(ji) + zdh5 ) + ( 1 - iiceform ) * zdh2 ) * zta3   & 
(  195)             &          + MIN ( zhnews2 , zhice0(ji) )                                             * zta4   &
(  196)             &          ) / zah
(  197)          
(  198)          tbif_1d(ji,3) =     (  iiceform * ( zhnews2 - zdh3 )                                          * zta1  &
(  199)             &              + ( iiceform * zdh3 + ( 1 - iiceform ) * zdh1 )                             * zta2  &
(  200)             &              + ( iiceform * ( zhnews2 - zdh5 ) + ( 1 - iiceform ) * ( zhnews2 - zdh1 ) ) * zta3  & 
(  201)             &              + ( iiceform * zdh5 + ( 1 - iiceform ) * zhnews2 )                          * zta4  &
(  202)             &            ) / zah
(  203)          !---removing the remaining part of ice formed which has been already used
(  204)          zbeta         = h_ice_1d(ji) / ( h_ice_1d(ji) + zdhicbot(ji) )
(  205)          h_ice_1d(ji)  = h_ice_1d(ji) + zdhicbot(ji)
(  206)          tbif_1d (ji,2)= zbeta * tbif_1d(ji,2) + ( 1.0 - zbeta ) * tbif_1d(ji,3)
(  207)          tbif_1d (ji,3)= ( 2. * zbeta - 1.0 ) * tbif_1d(ji,3) + ( 2. * zdhicbot(ji) / h_ice_1d(ji) ) * tfu_1d(ji)
(  208)          
(  209)       END DO
(  210)       
(  211)       !-------------------------------------------------------------
(  212)       !    Computation of variation of ice volume and ice mass 
(  213)       !           Vold = (1-Aold) * hiold ; Vnew = (1-Anew) * hinew
(  214)       !           dV = Vnew - Vold
(  215)       !-------------------------------------------------------------
(  216)       
(  217)       DO ji = kideb , kiut
(  218)          dvlbq_1d  (ji) = ( 1. - frld_1d(ji) ) * h_ice_1d(ji) - ( 1. - zfrl_old(ji) ) * zhice_old(ji)
(  219)          rdmicif_1d(ji) = rdmicif_1d(ji) + rhoic * dvlbq_1d(ji)
(  220)       END DO
(  221)       
(  222)    END SUBROUTINE lim_thd_lac
(  223) # 232
(  232) END MODULE limthd_lac
