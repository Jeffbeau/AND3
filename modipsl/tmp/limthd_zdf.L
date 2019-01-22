


PGF90 (Version     10.5)          09/23/2010  09:49:57      page 1

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

Filename: limthd_zdf.F90

(    1) MODULE limthd_zdf
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limthd_zdf ***
(    4)    !!                thermodynamic growth and decay of the ice 
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_thd_zdf  : vertical accr./abl. and lateral ablation of sea ice
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce          ! ocean parameters
(   14)    USE phycst           ! ???
(   15)    USE ice_oce          ! ice variables
(   16)    USE thd_ice
(   17)    USE iceini
(   18)    USE limistate
(   19)    USE in_out_manager
(   20)       
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_thd_zdf        ! called by lim_thd
(   26) 
(   27)    !! * Module variables
(   28)    REAL(wp)  ::           &  ! constant values
(   29)       epsi20 = 1.e-20  ,  &
(   30)       epsi13 = 1.e-13  ,  &
(   31)       zzero  = 0.e0    ,  &
(   32)       zone   = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd_zdf.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE lim_thd_zdf( kideb , kiut )
(   41)       !!------------------------------------------------------------------
(   42)       !!                ***  ROUTINE lim_thd_zdf  ***
(   43)       !!              
(   44)       !! ** Purpose : This routine determines the time evolution of snow 
(   45)       !!      and sea-ice thicknesses, concentration and heat content 
(   46)       !!      due to the vertical and lateral thermodynamic accretion-






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 2

(   47)       !!      ablation processes. One only treats the case of lat. abl.
(   48)       !!      For lateral accretion, see routine lim_lat_accr 
(   49)       !! 
(   50)       !! ** Method  : The representation of vertical growth and decay of 
(   51)       !!      the sea-ice model is based upon the diffusion of heat 
(   52)       !!      through the external and internal boundaries of a 
(   53)       !!      three-layer system (two layers of ice and one layer and 
(   54)       !!      one layer of snow, if present, on top of the ice).
(   55)       !! 
(   56)       !! ** Action  : - Calculation of some intermediates variables
(   57)       !!              - Calculation of surface temperature
(   58)       !!              - Calculation of available heat for surface ablation
(   59)       !!              - Calculation of the changes in internal temperature 
(   60)       !!                of the three-layer system, due to vertical diffusion 
(   61)       !!                processes
(   62)       !!              - Performs surface ablation and bottom accretion-ablation
(   63)       !!              - Performs snow-ice formation
(   64)       !!              - Performs lateral ablation
(   65)       !!
(   66)       !! References :
(   67)       !!   Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646   
(   68)       !!   Fichefet T. and M. Maqueda 1999, Clim. Dyn, 15(4), 251-268  
(   69)       !!
(   70)       !! History :
(   71)       !!   original    : 01-04 (LIM)
(   72)       !!   addition    : 02-08 (C. Ethe, G. Madec)
(   73)       !!------------------------------------------------------------------
(   74)       !! * Arguments
(   75)       INTEGER , INTENT (in) ::  &
(   76)          kideb ,  &  ! Start point on which the  the computation is applied
(   77)          kiut        ! End point on which the  the computation is applied
(   78) 
(   79)       !! * Local variables
(   80)       INTEGER ::   ji       ! dummy loop indices
(   81) 
(   82)       REAL(wp) , DIMENSION(jpij,2) ::  &
(   83)          zqcmlt        ! energy due to surface( /1 ) and bottom melting( /2 )
(   84) 
(   85)       REAL(wp), DIMENSION(jpij) ::  &
(   86)          ztsmlt      &    ! snow/ice surface melting temperature
(   87)          ,ztbif      &    ! int. temp. at the mid-point of the 1st layer of the snow/ice sys. 
(   88)          ,zksn       &    ! effective conductivity of snow
(   89)          ,zkic       &    ! effective conductivity of ice
(   90)          ,zksndh     &    ! thermal cond. at the mid-point of the 1st layer of the snow/ice sys. 
(   91)          , zfcsu     &    ! conductive heat flux at the surface of the snow/ice system 
(   92)          , zfcsudt   &    ! = zfcsu * dt
(   93)          , zi0       &    ! frac. of the net SW rad. which is not absorbed at the surface
(   94)          , z1mi0     &    ! fraction of the net SW radiation absorbed at the surface
(   95)          , zqmax     &    ! maximum energy stored in brine pockets
(   96)          , zrcpdt    &    ! h_su*rho_su*cp_su/dt(h_su being the thick. of surf. layer)
(   97)          , zts_old   &    ! previous surface temperature
(   98)          , zidsn , z1midsn , zidsnic ! tempory variables
(   99) 
(  100)       REAL(wp), DIMENSION(jpij) :: &
(  101)           zfnet       &  ! net heat flux at the top surface( incl. conductive heat flux)
(  102)           , zsprecip  &    ! snow accumulation
(  103)           , zhsnw_old &    ! previous snow thickness
(  104)           , zdhictop  &    ! change in ice thickness due to top surf ablation/accretion






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 3

(  105)           , zdhicbot  &    ! change in ice thickness due to bottom surf abl/acc
(  106)           , zqsup     &    ! energy transmitted to ocean (coming from top surf abl/acc)
(  107)           , zqocea    &    ! energy transmitted to ocean (coming from bottom sur abl/acc)
(  108)           , zfrl_old  &    ! previous sea/ice fraction
(  109)           , zfrld_1d    &    ! new sea/ice fraction
(  110)           , zep            ! internal temperature of the 2nd layer of the snow/ice system
(  111) 
(  112)        REAL(wp), DIMENSION(3) :: & 
(  113)           zplediag  &    ! principle diagonal, subdiag. and supdiag. of the 
(  114)           , zsubdiag  &    ! tri-diagonal matrix coming from the computation
(  115)           , zsupdiag  &    ! of the temperatures inside the snow-ice system
(  116)           , zsmbr          ! second member
(  117) 
(  118)        REAL(wp) :: & 
(  119)           zhsu     &     ! thickness of surface layer
(  120)           , zhe      &     ! effective thickness for compu. of equ. thermal conductivity
(  121)           , zheshth  &     ! = zhe / thth
(  122)           , zghe     &     ! correction factor of the thermal conductivity
(  123)           , zumsb    &     ! parameter for numerical method to solve heat-diffusion eq.
(  124)           , zkhsn    &     ! conductivity at the snow layer
(  125)           , zkhic    &     ! conductivity at the ice layers
(  126)           , zkint    &     ! equivalent conductivity at the snow-ice interface
(  127)           , zkhsnint &     ! = zkint*dt / (hsn*rhosn*cpsn)  
(  128)           , zkhicint &     ! = 2*zkint*dt / (hic*rhoic*cpic)
(  129)           , zpiv1 , zpiv2  &       ! tempory scalars used to solve the tri-diagonal system
(  130)           , zb2 , zd2 , zb3 , zd3 &
(  131)           , ztint          ! equivalent temperature at the snow-ice interface
(  132) 
(  133)        REAL(wp) :: & 
(  134)           zexp      &     ! exponential function of the ice thickness
(  135)           , zfsab     &     ! part of solar radiation stored in brine pockets
(  136)           , zfts      &     ! value of energy balance function when the temp. equal surf. temp.
(  137)           , zdfts     &     ! value of derivative of ztfs when the temp. equal surf. temp.
(  138)           , zdts      &     ! surface temperature increment
(  139)           , zqsnw_mlt &     ! energy needed to melt snow
(  140)           , zdhsmlt   &     ! change in snow thickness due to melt
(  141)           , zhsn      &     ! snow thickness (previous+accumulation-melt)
(  142)           , zqsn_mlt_rem &  ! remaining heat coming from snow melting
(  143)           , zqice_top_mlt & ! energy used to melt ice at top surface
(  144)           , zdhssub      &  ! change in snow thick. due to sublimation or evaporation
(  145)           , zdhisub      &  ! change in ice thick. due to sublimation or evaporation    
(  146)           , zdhsn        &  ! snow ice thickness increment
(  147)           , zdtsn        &  ! snow internal temp. increment
(  148)           , zdtic        &  ! ice internal temp. increment
(  149)           , zqnes          ! conductive energy due to ice melting in the first ice layer
(  150) 
(  151)        REAL(wp) :: & 
(  152)           ztbot     &      ! temperature at the bottom surface
(  153)           , zfcbot    &      ! conductive heat flux at bottom surface
(  154)           , zqice_bot &      ! energy used for bottom melting/growing
(  155)           , zqice_bot_mlt &  ! energy used for bottom melting
(  156)           , zqstbif_bot  &  ! part of energy stored in brine pockets used for bottom melting
(  157)           , zqstbif_old  &  ! tempory var. for zqstbif_bot
(  158)           , zdhicmlt      &  ! change in ice thickness due to bottom melting
(  159)           , zdhicm        &  ! change in ice thickness var. 
(  160)           , zdhsnm        &  ! change in snow thickness var. 
(  161)           , zhsnfi        &  ! snow thickness var. 
(  162)           , zc1, zpc1, zc2, zpc2, zp1, zp2 & ! tempory variables






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 4

(  163)           , ztb2, ztb3
(  164) 
(  165)        REAL(wp) :: & 
(  166)           zdrmh         &   ! change in snow/ice thick. after snow-ice formation
(  167)           , zhicnew       &   ! new ice thickness
(  168)           , zhsnnew       &   ! new snow thickness
(  169)           , zquot , ztneq &   ! tempory temp. variables
(  170)           , zqice, zqicetot & ! total heat inside the snow/ice system
(  171)           , zdfrl         &   ! change in ice concentration
(  172)           , zdvsnvol      &   ! change in snow volume
(  173)           , zdrfrl1, zdrfrl2 &  ! tempory scalars
(  174)           , zihsn, zidhb, zihic, zihe, zihq, ziexp, ziqf, zihnf, zibmlt, ziqr, zihgnew, zind
(  175)        !!----------------------------------------------------------------------
(  176) 
(  177)        !-----------------------------------------------------------------------
(  178)        !  1. Boundaries conditions for snow/ice system internal temperature
(  179)        !       - If tbif_1d(ji,1) > rt0_snow, tbif_1d(ji,1) = rt0_snow 
(  180)        !       - If tbif_1d(ji,2/3) > rt0_ice, tbif_1d(ji,2/3) = rt0_ice 
(  181)        !     Computation of energies due to surface and bottom melting 
(  182)        !-----------------------------------------------------------------------
(  183)        
(  184)        
(  185)        DO ji = kideb , kiut
(  186)           zihsn = MAX( zzero , SIGN( zone , hsndif - h_snow_1d(ji) ) )
(  187)           zihic = MAX( zzero , SIGN( zone , hicdif - h_ice_1d(ji) ) )
(  188)           !--computation of energy due to surface melting
(  189)           zqcmlt(ji,1) = ( MAX ( zzero ,  &
(  190)              &                   rcpsn * h_snow_1d(ji) * ( tbif_1d(ji,1) - rt0_snow ) ) ) * ( 1.0 - zihsn )
(  191)           !--computation of energy due to bottom melting
(  192)           zqcmlt(ji,2) = ( MAX( zzero , &
(  193)              &                  rcpic * ( tbif_1d(ji,2) - rt0_ice ) * ( h_ice_1d(ji) / 2. ) ) &
(  194)              &           + MAX( zzero , &
(  195)              &                  rcpic * ( tbif_1d(ji,3) - rt0_ice ) * ( h_ice_1d(ji) / 2. ) ) &
(  196)              &           ) * ( 1.0 - zihic  )
(  197)           !--limitation of  snow/ice system internal temperature
(  198)           tbif_1d(ji,1)   = MIN( rt0_snow, tbif_1d(ji,1) )
(  199)           tbif_1d(ji,2)   = MIN( rt0_ice , tbif_1d(ji,2) )
(  200)           tbif_1d(ji,3)   = MIN( rt0_ice , tbif_1d(ji,3) )
(  201)        END DO
(  202) 
(  203) 
(  204)        !-------------------------------------------
(  205)        !  2. Calculate some intermediate variables.  
(  206)        !-------------------------------------------
(  207)        
(  208)        ! initialisation of the thickness of surface layer
(  209)        zhsu = hnzst  
(  210) 
(  211)        DO ji = kideb , kiut
(  212)           zind   = MAX( zzero , SIGN( zone , zhsu - h_snow_1d(ji) ) )
(  213)           zihsn  = MAX( zzero , SIGN( zone , hsndif - h_snow_1d(ji) ) )
(  214)           zihsn  = MAX( zihsn , zind )
(  215)           zihic  = MAX( zzero , sign( zone , hicdif - h_ice_1d(ji) ) )
(  216)           !     2.1. Computation of surface melting temperature
(  217)           !----------------------------------------------------
(  218)           zind  = MAX( zzero , SIGN( zone , -h_snow_1d(ji) ) )
(  219)           ztsmlt(ji) = ( 1.0 - zind ) * rt0_snow + zind * rt0_ice
(  220)           !






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 5

(  221)           !     2.2. Effective conductivity of snow and ice
(  222)           !-----------------------------------------------
(  223) 
(  224)           !---computation of the correction factor on the thermal conductivity
(  225)           !-- (Morales Maqueda, 1995 ; Fichefet and Morales Maqueda, 1997)
(  226)           zhe      =  ( rcdsn / ( rcdsn + rcdic ) ) * h_ice_1d(ji)   &
(  227)              &     + ( rcdic / ( rcdsn + rcdic ) ) * h_snow_1d(ji) 
(  228)           zihe     = MAX( zzero , SIGN( zone , 2.0 * zhe - thth ) )
(  229)           zheshth  = zhe / thth
(  230)           zghe     = ( 1.0 - zihe ) * zheshth * ( 2.0 - zheshth )   &
(  231)              &     +         zihe   * 0.5 * ( 1.5 + LOG( 2.0 * zheshth ) )
(  232) #if defined key_lim_cp3
(  233)           zghe = 1.0
(  234) #endif 
(  235) 
(  236)           !---effective conductivities 
(  237)           zksn(ji)  = zghe * rcdsn  
(  238)           zkic(ji)  = zghe * rcdic
(  239) 
(  240)           !
(  241)           !     2.3. Computation of the conductive heat flux from the snow/ice 
(  242)           !          system interior toward the top surface
(  243)           !------------------------------------------------------------------
(  244) 
(  245)           !---Thermal conductivity at the mid-point of the first snow/ice system layer
(  246)           zksndh(ji) =   ( ( 1.0 - zihsn ) * 2.0 * zksn(ji) + zihsn * 4.0 * zkic(ji) )   &
(  247)              &         / ( ( 1.0 - zihsn ) *  h_snow_1d(ji)                              &
(  248)              &           +        zihsn   *  ( ( 1.0 + 3.0 * zihic ) * h_ice_1d(ji)      &
(  249)              &           + 4.0 * zkic(ji)/zksn(ji) * h_snow_1d(ji) ) )
(  250) 
(  251)           !---internal temperature at the mid-point of the first snow/ice system layer
(  252)           ztbif(ji)  = ( 1.0 - zihsn ) * tbif_1d(ji,1)                       &
(  253)              &       +         zihsn   * ( ( 1.0 - zihic ) * tbif_1d(ji,2)   &
(  254)              &       +         zihic   * tfu_1d(ji)   )
(  255)           !---conductive heat flux 
(  256)           zfcsu(ji) = zksndh(ji) * ( ztbif(ji) - sist_1d(ji) )
(  257) 
(  258)        END DO
(  259) 
(  260)        !--------------------------------------------------------------------
(  261)        !  3. Calculate : 
(  262)        !     - fstbif_1d, part of solar radiation absorbing inside the ice
(  263)        !       assuming an exponential absorption (Grenfell and Maykut, 1977)
(  264)        !     - zqmax,  maximum energy stored in brine pockets
(  265)        !     - qstbif_1d, total energy stored in brine pockets (updating)
(  266)        !-------------------------------------------------------------------
(  267) 
(  268) 
(  269)        DO ji = kideb , kiut
(  270)           zihsn  = MAX( zzero , SIGN (zone , -h_snow_1d(ji) ) )
(  271)           zihic  = MAX( zzero , 1.0 - ( h_ice_1d(ji) / zhsu ) )     
(  272)           zind   = MAX( zzero , SIGN (zone , hicdif - h_ice_1d(ji) ) )
(  273)           !--Computation of the fraction of the net shortwave radiation which
(  274)           !--penetrates inside the ice cover ( See Forcat)
(  275)           zi0(ji)  = zihsn * ( fr1_i0_1d(ji) + zihic * fr2_i0_1d(ji) )
(  276)           zexp     = MIN( zone , EXP( -1.5 * ( h_ice_1d(ji) - zhsu ) ) )
(  277)           fstbif_1d(ji) = zi0(ji) * qsr_ice_1d(ji) * zexp
(  278)           !--Computation of maximum energy stored in brine pockets zqmax and update






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 6

(  279)           !--the total energy stored in brine pockets, if less than zqmax
(  280)           zqmax(ji) = MAX( zzero , 0.5 * xlic * ( h_ice_1d(ji) - hicmin ) )
(  281)           zfsab   = zi0(ji) * qsr_ice_1d(ji) * ( 1.0 - zexp )
(  282)           zihq    = ( 1.0 - zind ) * MAX(zzero, SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) ) &
(  283)              &    +         zind   * zone
(  284)           qstbif_1d(ji) = ( qstbif_1d(ji) + ( 1.0 - zihq ) * zfsab * rdt_ice ) * swiqst
(  285)           !--fraction of shortwave radiation absorbed at surface
(  286)           ziexp = zihq * zexp + ( 1.0 - zihq ) * ( swiqst + ( 1.0 - swiqst ) * zexp )
(  287)           z1mi0(ji) = 1.0 - zi0(ji) * ziexp
(  288)        END DO
(  289) 
(  290) 
(  291)        !--------------------------------------------------------------------------------
(  292)        !  4. Computation of the surface temperature : determined by considering the 
(  293)        !     budget of a thin layer of thick. zhsu at the top surface (H. Grenier, 1995) 
(  294)        !     and based on a surface energy balance : 
(  295)        !     hsu * rcp * dT/dt = Fsr + Fnsr(T) + Fcs(T),
(  296)        !     where - Fsr is the net absorbed solar radiation, 
(  297)        !           - Fnsr is the total non solar radiation (incoming and outgoing long-wave,
(  298)        !             sensible and latent heat fluxes)
(  299)        !           - Fcs the conductive heat flux at the top of surface
(  300)        !------------------------------------------------------------------------------
(  301) 
(  302)        !     4.1. Computation of intermediate values
(  303)        !---------------------------------------------
(  304)        DO ji = kideb, kiut
(  305)           zrcpdt(ji) = ( rcpsn * MIN( h_snow_1d(ji) , zhsu )    &
(  306)              &       + rcpic * MAX( zhsu - h_snow_1d(ji) , zzero ) ) / rdt_ice
(  307)           zts_old(ji) =  sist_1d(ji)
(  308)        END DO
(  309) 
(  310)        !     4.2. Computation of surface temperature by expanding the eq. of energy balance
(  311)        !          with Ts = Tp + DT. One obtain , F(Tp) + DT * DF(Tp) = 0
(  312)        !          where  - F(Tp) = Fsr + Fnsr(Tp) + Fcs(Tp) 
(  313)        !                 - DF(Tp)= (dFnsr(Tp)/dT) + (dFcs(Tp)/dT) - hsu*rcp/dt
(  314)        !---------------------------------------------------------------------------------
(  315) 
(  316)        DO ji = kideb, kiut
(  317)           !---computation of the derivative of energy balance function 
(  318) #if defined key_coupled
(  319) #   if defined key_lim_cp2
(  320)           zdfts   =   zksndh(ji)   & ! contribution of the conductive heat flux
(  321)              &      + zrcpdt(ji)   & ! contribution of hsu * rcp / dt
(  322)              &      - dqns_ice_1d(ji)      ! contribution of the total non solar radiation 
(  323) #   else
(  324)           zdfts   =   zksndh(ji)   & ! contribution of the conductive heat flux
(  325)              &      + zrcpdt(ji)    ! contribution of hsu * rcp / dt
(  326) #   endif
(  327) 
(  328) #else
(  329)           zdfts    =  zksndh(ji)   & ! contribution of the conductive heat flux
(  330)              &      + zrcpdt(ji)   & ! contribution of hsu * rcp / dt
(  331)              &      - dqns_ice_1d (ji)     ! contribution of the total non solar radiation 
(  332) #endif
(  333)           !---computation of the energy balance function 
(  334)           zfts    = - z1mi0 (ji) * qsr_ice_1d(ji)   & ! net absorbed solar radiation
(  335)              &      - qnsr_ice_1d(ji)                & ! total non solar radiation
(  336)              &      - zfcsu (ji)                  ! conductive heat flux from the surface






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 7

(  337)           !---computation of surface temperature increment  
(  338)           zdts    = -zfts / zdfts
(  339) #if defined key_lim_cp3
(  340)           zdts = zdts / 3.0
(  341) #endif
(  342)           !---computation of the new surface temperature 
(  343)           sist_1d(ji) = sist_1d(ji) + zdts
(  344) 
(  345)        END DO
(  346) 
(  347)        !----------------------------------------------------------------------------
(  348)        !  5. Boundary condition at the top surface
(  349)        !--    IF Tsb < Tmelt, Fnet = Fcs (the net heat flux equal the conductive heat flux)
(  350)        !      Otherwise Tsb = Tmelt and Qnet(Tmelt) > 0 
(  351)        !      Fnet(Tmelt) is therefore the net surface flux needed for melting
(  352)        !----------------------------------------------------------------------------
(  353)        
(  354)        
(  355)        !     5.1.  Limitation of surface temperature and update total non solar fluxes,
(  356)        !          latent heat flux and conductive flux at the top surface 
(  357)        !----------------------------------------------------------------------  
(  358)                      
(  359)        DO ji = kideb, kiut
(  360)           sist_1d(ji) = MIN( ztsmlt(ji) , sist_1d(ji) )
(  361) #if ! defined key_coupled
(  362)           qnsr_ice_1d(ji) = qnsr_ice_1d(ji) + dqns_ice_1d(ji) * ( sist_1d(ji) - zts_old(ji) )
(  363)           qla_ice_1d (ji) = qla_ice_1d (ji) + dqla_ice_1d(ji) * ( sist_1d(ji) - zts_old(ji) )
(  364) #endif
(  365)           zfcsu(ji)  = zksndh(ji) * ( ztbif(ji) - sist_1d(ji) )
(  366)        END DO
(  367) 
(  368)      
(  369) 
(  370)        !     5.2. Calculate available heat for surface ablation. 
(  371)        !---------------------------------------------------------------------
(  372) 
(  373)        DO ji = kideb, kiut
(  374)           zfnet(ji) = qnsr_ice_1d(ji) + z1mi0(ji) * qsr_ice_1d(ji) + zfcsu(ji)          
(  375)           zfnet(ji) = MAX( zzero , zfnet(ji) )
(  376)           zfnet(ji) = zfnet(ji) * MAX( zzero , SIGN( zone , sist_1d(ji) - ztsmlt(ji) ) )
(  377)        END DO
(  378) 
(  379)        !-------------------------------------------------------------------------
(  380)        !  6. Calculate changes in internal temperature due to vertical diffusion   
(  381)        !     processes. The evolution of this temperature is governed by the one-
(  382)        !     dimensionnal heat-diffusion equation. 
(  383)        !     Given the temperature tbif(1/2/3), at time m we solve a set
(  384)        !     of finite difference equations to obtain new tempe. Each tempe is coupled
(  385)        !     to the temp. immediatly above and below by heat conduction terms. Thus 
(  386)        !     we have a set of equations of the form A * T = B, where A is a tridiagonal
(  387)        !     matrix, T a vector whose components are the unknown new temp.
(  388)        !-------------------------------------------------------------------------
(  389)        
(  390)        !--parameter for the numerical methode use to solve the heat-diffusion equation
(  391)        !- implicit, explicit or Crank-Nicholson
(  392)        zumsb = 1.0 - sbeta  
(  393)        DO ji = kideb, kiut
(  394)           zidsn(ji)   = MAX ( zzero, SIGN( zone, hsndif - h_snow_1d(ji) ) ) 






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 8

(  395)           z1midsn(ji) = 1.0 - zidsn(ji)
(  396)           zihic       = MAX ( zzero, SIGN( zone, hicdif - h_ice_1d(ji) ) ) 
(  397)           zidsnic(ji) = zidsn(ji) *  zihic 
(  398)           zfcsudt(ji) = zfcsu(ji) * rdt_ice 
(  399)        END DO 
(  400)    
(  401)        DO ji = kideb, kiut
(  402) 
(  403)           !     6.1 Calculate intermediate variables.
(  404)           !----------------------------------------
(  405) 
(  406)           !--conductivity at the snow surface
(  407)           zkhsn = 2.0 * zksn(ji) * rdt_ice / rcpsn
(  408)           !--conductivity at the ice surface
(  409)           zkhic = 4.0 * zkic(ji) * rdt_ice / MAX( h_ice_1d(ji) * h_ice_1d(ji) * rcpic , epsi20 )
(  410)           !--conductivity at the snow/ice interface 
(  411)           zkint = 4.0 * zksn(ji) * zkic(ji)  &
(  412)              &        / ( zksn(ji) * h_ice_1d(ji) + 2.0 * zkic(ji) * h_snow_1d(ji) * z1midsn(ji)) 
(  413)           zkhsnint = zkint * rdt_ice / rcpsn
(  414)           zkhicint = zkint * 2.0 * rdt_ice / MAX( h_ice_1d(ji) * rcpic , epsi20 )
(  415)           
(  416)           !     6.2. Fulfill the linear system matrix.
(  417)           !-----------------------------------------
(  418) !$$$          zplediag(1) = 1 + sbeta * z1midsn(ji) * ( zkhsn + zkhsnint )       
(  419)           zplediag(1) =   zidsn(ji) + z1midsn(ji) * h_snow_1d(ji)   &
(  420)              &          + sbeta * z1midsn(ji) * zkhsnint 
(  421)           zplediag(2) = 1 + sbeta * ( z1midsn(ji) * zkhicint + zkhic ) 
(  422)           zplediag(3) = 1 + 3.0 * sbeta * zkhic   
(  423) 
(  424)           zsubdiag(1) =  0.e0              
(  425)           zsubdiag(2) = -1.e0 * z1midsn(ji) * sbeta * zkhicint
(  426)           zsubdiag(3) = -1.e0 * sbeta * zkhic 
(  427) 
(  428)           zsupdiag(1) = -1.e0 * z1midsn(ji) * sbeta * zkhsnint 
(  429)           zsupdiag(2) = zsubdiag(3)
(  430)           zsupdiag(3) =  0.e0
(  431)           
(  432)           !     6.3. Fulfill the idependent term vector.
(  433)           !-------------------------------------------
(  434)           
(  435) !$$$          zsmbr(1) = zidsn(ji) * sist_1d(ji) + z1midsn(ji) *   &
(  436) !$$$             &         ( tbif_1d(ji,1) + zkhsn * sist_1d(ji)
(  437) !$$$             &         - zumsb * ( zkhsn * tbif_1d(ji,1)
(  438) !$$$             &                   + zkhsnint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) ) )
(  439)           zsmbr(1) = zidsn(ji) * sist_1d(ji) + z1midsn(ji) *    &
(  440)              &       ( h_snow_1d(ji) * tbif_1d(ji,1) - ( zfcsudt(ji) / rcpsn )  &
(  441)              &       - zumsb * zkhsnint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) )
(  442) 
(  443)           zsmbr(2) =  tbif_1d(ji,2)  &
(  444)              &      - zidsn(ji) * ( 1.0 - zidsnic(ji) ) &
(  445)              &        * ( zfcsudt(ji) / MAX( h_ice_1d(ji) * rcpic , epsi20 ) ) &
(  446)              &      + zumsb * ( zkhicint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) &
(  447)              &                   - zkhic * ( tbif_1d(ji,2) - tbif_1d(ji,3) )  )
(  448) 
(  449)           zsmbr(3) =  tbif_1d(ji,3)  &
(  450)              &      + zkhic * ( 2.0 * tfu_1d(ji) &
(  451)              &                + zumsb * ( tbif_1d(ji,2) - 3.0 * tbif_1d(ji,3) ) )
(  452)           






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 9

(  453)           !     6.4. Solve the system (Gauss elimination method).
(  454)           !----------------------------------------------------
(  455)           
(  456)           zpiv1 = zsubdiag(2) / zplediag(1) 
(  457)           zb2   = zplediag(2) - zpiv1 * zsupdiag(1)
(  458)           zd2   = zsmbr(2) - zpiv1 * zsmbr(1)
(  459) 
(  460)           zpiv2 = zsubdiag(3) / zb2
(  461)           zb3   = zplediag(3) - zpiv2 * zsupdiag(2)
(  462)           zd3   = zsmbr(3) - zpiv2 * zd2
(  463) 
(  464)           tbif_1d(ji,3) = zd3 / zb3
(  465)           tbif_1d(ji,2) = ( zd2 - zsupdiag(2) * tbif_1d(ji,3) ) / zb2
(  466)           tbif_1d(ji,1) = ( zsmbr(1) - zsupdiag(1) * tbif_1d(ji,2) ) / zplediag(1)            
(  467) 
(  468)           !--- taking into account the particular case of  zidsnic(ji) = 1
(  469)           ztint =  (  zkic(ji) * h_snow_1d(ji) * tfu_1d (ji)    &
(  470)              &      + zksn(ji) * h_ice_1d(ji) * sist_1d(ji) )   &
(  471)              &   / ( zkic(ji) * h_snow_1d(ji) + zksn(ji) * h_ice_1d(ji) ) 
(  472) 
(  473)           tbif_1d(ji,1) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,1)   &
(  474)              &                + zidsnic(ji)   * ( ztint + sist_1d(ji) ) / 2.0
(  475)           tbif_1d(ji,2) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,2)   &
(  476)              &                + zidsnic(ji)   * ( 3.0 * ztint + tfu_1d(ji) ) / 4.0
(  477)           tbif_1d(ji,3) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,3)   &
(  478)              &                + zidsnic(ji)   * ( ztint + 3.0 * tfu_1d(ji) ) / 4.0     
(  479)        END DO
(  480)  
(  481)        !----------------------------------------------------------------------
(  482)        !  9. Take into account surface ablation and bottom accretion-ablation.|
(  483)        !----------------------------------------------------------------------
(  484)        
(  485)        !---Snow accumulation in one thermodynamic time step
(  486)        zsprecip(kideb:kiut) = sprecip_1d(kideb:kiut) * rdt_ice / rhosn
(  487) 
(  488) 
(  489)        DO ji = kideb, kiut
(  490)           
(  491)           !      9.1. Surface ablation and update of snow thickness and qstbif_1d
(  492)           !--------------------------------------------------------------------
(  493)           
(  494)           !--------------------------------------------------------------------------
(  495)           !--      Melting snow processes :
(  496)           !--      Melt at the upper surface is computed from the difference between 
(  497)           !--      the net heat flux (including the conductive heat flux) at the upper 
(  498)           !--      surface and the pre-existing energy due to surface melting
(  499)           !------------------------------------------------------------------------------
(  500)           
(  501)           !-- store the snow thickness
(  502)           zhsnw_old(ji) =  h_snow_1d(ji)
(  503)           !--computation of the energy needed to melt snow
(  504)           zqsnw_mlt  = zfnet(ji) * rdt_ice - zqcmlt(ji,1)
(  505)           !--change in snow thickness due to melt
(  506)           zdhsmlt = - zqsnw_mlt / xlsn
(  507)           
(  508)           !-- compute new snow thickness, taking into account the part of snow accumulation
(  509)           !   (as snow precipitation) and the part of snow lost due to melt
(  510)           zhsn =  h_snow_1d(ji) + zsprecip(ji) + zdhsmlt






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 10

(  511)           h_snow_1d(ji) = MAX( zzero , zhsn )
(  512)           !-- compute the volume of snow lost after surface melting and the associated mass
(  513)           dvsbq_1d(ji) =  ( 1.0 - frld_1d(ji) ) * ( h_snow_1d(ji) - zhsnw_old(ji) - zsprecip(ji) )
(  514)           dvsbq_1d(ji) =  MIN( zzero , dvsbq_1d(ji) )
(  515)           rdmsnif_1d(ji) =  rhosn * dvsbq_1d(ji)
(  516)           !-- If the snow is completely melted the remaining heat is used to melt ice
(  517)           zqsn_mlt_rem  = MAX( zzero , -zhsn ) * xlsn
(  518)           zqice_top_mlt = zqsn_mlt_rem 
(  519)           zqstbif_old   = qstbif_1d(ji)
(  520) 
(  521)           !--------------------------------------------------------------------------
(  522)           !--      Melting ice processes at the top surface :
(  523)           !--      The energy used to melt ice, zqice_top_mlt, is taken from the energy
(  524)           !--      stored in brine pockets qstbif_1d and the remaining energy coming 
(  525)           !--      from the melting snow process zqsn_mlt_rem.
(  526)           !--      If qstbif_1d > zqsn_mlt_rem then, one uses only a zqsn_mlt_rem part
(  527)           !--      of qstbif_1d to melt ice,
(  528)           !--         zqice_top_mlt = zqice_top_mlt + zqsn_mlt_rem
(  529)           !--         qstbif_1d = qstbif_1d - zqsn_mlt_rem
(  530)           !--      Otherwise one uses all qstbif_1d to melt ice
(  531)           !--         zqice_top_mlt = zqice_top_mlt + qstbif_1d
(  532)           !--         qstbif_1d = 0
(  533)           !------------------------------------------------------
(  534)           
(  535)           ziqf =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqsn_mlt_rem  ) )
(  536)           zqice_top_mlt =         ziqf   * ( zqice_top_mlt + zqsn_mlt_rem )   &
(  537)              &          + ( 1.0 - ziqf ) * ( zqice_top_mlt + qstbif_1d(ji)  )
(  538) 
(  539)           qstbif_1d(ji) =         ziqf   * ( qstbif_1d(ji) - zqsn_mlt_rem )   &
(  540)              &        + ( 1.0 - ziqf ) * ( qstbif_1d(ji) - qstbif_1d(ji)  )
(  541) 
(  542)           !--    The contribution of the energy stored in brine pockets qstbif_1d to melt
(  543)           !--    ice is taking into account only when qstbif_1d is less than zqmax. 
(  544)           !--    Otherwise, only the remaining energy coming from the melting snow 
(  545)           !--    process is used 
(  546)           zihq =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) )
(  547) 
(  548)           zqice_top_mlt =         zihq   * zqice_top_mlt   &
(  549)              &          + ( 1.0 - zihq ) * zqsn_mlt_rem
(  550) 
(  551)           qstbif_1d(ji) =         zihq   * qstbif_1d(ji)   &
(  552)              &        + ( 1.0 - zihq ) * zqstbif_old
(  553) 
(  554)           !--change in ice thickness due to melt at the top surface
(  555)           zdhictop(ji) = -zqice_top_mlt / xlic
(  556)           !--compute the volume formed after surface melting
(  557)           dvsbq_1d(ji) =  zdhictop(ji) * ( 1.0 - frld_1d(ji) )
(  558) 
(  559)           !-------------------------------------------------------------------------
(  560)           !--      A small variation at the surface also occurs because of sublimation
(  561)           !--      associated with the latent flux. If qla_ice_1d is negative, snow condensates at 
(  562)           !        the surface. Otherwise, snow evaporates
(  563)           !-----------------------------------------------------------------------
(  564)           !----change in snow and ice thicknesses due to sublimation or evaporation
(  565)           zdhssub  = parsub * ( qla_ice_1d(ji) / ( rhosn * xsn ) ) * rdt_ice 
(  566)           zhsn     = h_snow_1d(ji) - zdhssub
(  567)           zdhisub  = MAX( zzero , -zhsn ) * rhosn/rhoic
(  568)           zdhictop(ji) =  zdhictop(ji) - zdhisub






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 11

(  569)           h_snow_1d(ji)  =  MAX( zzero , zhsn )
(  570)           !-------------------------------------------------
(  571)           !--  Update Internal temperature and qstbif_1d.
(  572)           !-------------------------------------------
(  573)           zihsn  =  MAX( zzero , SIGN( zone, -h_snow_1d(ji) ) )
(  574)           tbif_1d(ji,1) = ( 1.0 - zihsn ) * tbif_1d(ji,1) + zihsn   * tfu_1d(ji)
(  575)           !--change in snow internal temperature if snow has increased
(  576)           zihnf = MAX( zzero , SIGN( zone , h_snow_1d(ji) - zhsnw_old(ji) ) )
(  577)           zdhsn = 1.0 - zhsnw_old(ji) / MAX( h_snow_1d(ji) , epsi20 )
(  578)           zdtsn = zdhsn * ( sist_1d(ji) - tbif_1d(ji,1) )
(  579)           tbif_1d(ji,1) = tbif_1d(ji,1) + z1midsn(ji) * zihnf * zdtsn
(  580)           !--energy created due to ice melting in the first ice layer
(  581)           zqnes  = ( rt0_ice - tbif_1d(ji,2) ) * rcpic * ( h_ice_1d(ji) / 2. )
(  582)           !--change in first ice layer internal temperature
(  583)           ziqr  = MAX( zzero , SIGN( zone , qstbif_1d(ji) - zqnes ) )
(  584)           zdtic = qstbif_1d(ji) / ( rcpic * ( h_ice_1d(ji) / 2. ) )
(  585)           tbif_1d(ji,2) =  ziqr * rt0_ice + ( 1 - ziqr ) * ( tbif_1d(ji,2) + zdtic )
(  586)           !--update qstbif_1d
(  587)           qstbif_1d(ji) = ziqr * ( qstbif_1d(ji) - zqnes ) * swiqst
(  588) 
(  589) 
(  590)           !--      9.2. Calculate bottom accretion-ablation and update qstbif_1d.
(  591)           !             Growth and melting at bottom ice surface are governed by  
(  592)           !                 -xlic * Dh = (Fcb - Fbot ) * Dt 
(  593)           !             where Fbot is the net downward heat flux from ice to the ocean
(  594)           !            and Fcb is the conductive heat flux at the bottom surface
(  595)           !---------------------------------------------------------------------------
(  596)           ztbot = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,3) + zidsnic(ji) * sist_1d(ji)
(  597)           !---computes conductive heat flux at bottom surface 
(  598)           zfcbot =  4.0 * zkic(ji) * ( tfu_1d(ji) - ztbot )   &
(  599)              &   / ( h_ice_1d(ji) + zidsnic(ji) * ( 3. * h_ice_1d(ji) &
(  600)              &   + 4.0 * zkic(ji)/zksn(ji) * h_snow_1d(ji) ) )
(  601)           !---computation of net energy needed for bottom melting/growing
(  602)           zqice_bot = ( zfcbot - ( fbif_1d(ji) + qlbbq_1d(ji) ) ) * rdt_ice
(  603)           zqstbif_bot = qstbif_1d(ji)
(  604)           !---switch to know if bottom surface melts ( = 1 ) or grows ( = 0 )occurs
(  605)           zibmlt = MAX( zzero , SIGN( zone , -zqice_bot ) )
(  606)           !--particular case of melting (in the same way as the top surface)
(  607)           zqice_bot_mlt = zqice_bot 
(  608)           zqstbif_old = zqstbif_bot
(  609) 
(  610)           ziqf =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) + zqice_bot_mlt  ) )
(  611)           zqice_bot_mlt =         ziqf   * ( zqice_bot_mlt + zqice_bot_mlt ) &
(  612)              &          + ( 1.0 - ziqf ) * ( zqice_bot_mlt + qstbif_1d(ji)  )          
(  613)           qstbif_1d(ji)   =         ziqf   * ( qstbif_1d(ji) + zqice_bot_mlt ) &
(  614)              &          + ( 1.0 - ziqf ) * ( qstbif_1d(ji) - qstbif_1d(ji)  )
(  615)           !--    The contribution of the energy stored in brine pockets qstbif_1d to melt
(  616)           !--    ice is taking into account only when qstbif_1d is less than zqmax. 
(  617)           zihq =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) )
(  618)           zqice_bot_mlt =         zihq   * zqice_bot_mlt   &
(  619)              &          + ( 1.0 - zihq ) * zqice_bot
(  620)           qstbif_1d(ji)   =         zihq   * qstbif_1d(ji)   &
(  621)              &             + ( 1.0 - zihq ) * zqstbif_old
(  622) 
(  623)           !---treatment of the case of melting/growing
(  624)           zqice_bot   =         zibmlt   * ( zqice_bot_mlt - zqcmlt(ji,2) )   &
(  625)              &        + ( 1.0 - zibmlt ) * ( zqice_bot - zqcmlt(ji,2)  )
(  626)           qstbif_1d(ji) =         zibmlt   * qstbif_1d(ji)   &






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 12

(  627)              &           + ( 1.0 - zibmlt ) * zqstbif_bot
(  628) 
(  629)           !--computes change in ice thickness due to melt or growth
(  630)           zdhicbot(ji) = zqice_bot / xlic
(  631)           !--limitation of bottom melting if so : hmelt maximum melting at bottom
(  632)           zdhicmlt  = MAX( hmelt , zdhicbot(ji) ) 
(  633)           !--energy after bottom melting/growing
(  634)           zqsup(ji) = ( 1.0 - frld_1d(ji) ) * xlic * ( zdhicmlt - zdhicbot(ji) )
(  635)           !-- compute the new thickness and the newly formed volume after bottom melting/growing
(  636)           zdhicbot(ji)  = zdhicmlt
(  637)           dvbbq_1d(ji) = ( 1.0 - frld_1d(ji) ) * zdhicbot(ji)
(  638) 
(  639) 
(  640)           !        9.3.  Updating ice thickness after top surface ablation 
(  641)           !              and bottom surface accretion/ablation
(  642)           !---------------------------------------------------------------
(  643)           zhicnew  = h_ice_1d(ji) + zdhictop(ji) + zdhicbot(ji)
(  644) 
(  645)           !
(  646)           !        9.4. Case of total ablation (ice is gone but snow may be left)
(  647)           !-------------------------------------------------------------------
(  648)           zhsn  = h_snow_1d(ji)
(  649)           zihgnew = 1.0 - MAX( zzero , SIGN( zone , -zhicnew ) )
(  650)           zihsn   = MAX( zzero , SIGN( zone , -zhsn ) )
(  651)           !---convert 
(  652)           zdhicm  = ( 1.0 - zihgnew ) * ( zhicnew - qstbif_1d(ji) / xlic )
(  653)           zdhsnm  = ( 1.0 - zihsn ) * zdhicm * rhoic / rhosn
(  654)           !---updating new ice thickness and computing the newly formed ice mass
(  655)           zhicnew   =  zihgnew * zhicnew
(  656)           rdmicif_1d(ji) =  rdmicif_1d(ji) + ( 1.0 - frld_1d(ji) ) * ( zhicnew - h_ice_1d(ji) ) * rhoic
(  657)           !---updating new snow thickness and computing the newly formed snow mass
(  658)           zhsnfi   = zhsn + zdhsnm
(  659)           h_snow_1d(ji) = MAX( zzero , zhsnfi )
(  660)           rdmsnif_1d(ji) =  rdmsnif_1d(ji) + ( 1.0 - frld_1d(ji) ) * ( h_snow_1d(ji) - zhsn ) * rhosn
(  661)           !--remaining energy in case of total ablation
(  662)           zqocea(ji) = - ( zihsn * xlic * zdhicm + xlsn * ( zhsnfi - h_snow_1d(ji) ) ) * ( 1.0 - frld_1d(ji) )
(  663)           qstbif_1d(ji) = zihgnew * qstbif_1d(ji)
(  664) 
(  665)           !
(  666)           !        9.5. Update internal temperature and ice thickness.
(  667)           !-------------------------------------------------------
(  668)           !
(  669)           sist_1d(ji) = zihgnew * sist_1d(ji) + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  670)           zidhb  = MAX( zzero , SIGN( zone , - zdhicbot(ji) ) )
(  671)           zc1    = - zhicnew * 0.5
(  672)           zpc1   = MIN( 0.5 * zone , - h_ice_1d(ji) * 0.5 - zdhictop(ji) )
(  673)           zc2    = - zhicnew
(  674)           zpc2   =  zidhb * zc2 + ( 1.0 - zidhb ) * ( - h_ice_1d(ji) - zdhictop(ji) )
(  675)           zp1    =  MAX( zpc1 , zc1 )
(  676)           zp2    =  MAX( zpc2 , zc1 )
(  677)           zep(ji) =  tbif_1d(ji,2)
(  678)           ztb2  = 2.0 * (         - zp1   * tbif_1d(ji,2)  &
(  679)              &  + ( zp1 - zp2 ) * tbif_1d(ji,3)  &
(  680)              &  + ( zp2 - zc1 ) * tfu_1d(ji) ) / MAX( zhicnew , epsi20 ) 
(  681)           tbif_1d(ji,2) = zihgnew * ztb2 + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  682)           !---
(  683)           zp1  =  MIN( zpc1 , zc1 )
(  684)           zp2  =  MIN( zpc2 , zc1 )






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 13

(  685)           zp1  =  MAX( zc2  , zp1 )
(  686)           ztb3 =  2.0 * (   ( 1.0 - zidhb ) * (  ( zc1 - zp2 ) * tbif_1d(ji,3)  &
(  687)              &                                 + ( zp2 - zc2 ) * tfu_1d(ji) )   &
(  688)              &               +      zidhb   * (  ( zc1 - zp1 ) * zep(ji)      &
(  689)              &                                 + ( zp1 - zc2 ) * tbif_1d(ji,3))  ) / MAX( zhicnew , epsi20 )
(  690)           tbif_1d(ji,3) =  zihgnew * ztb3 + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  691)           h_ice_1d(ji)  =  zhicnew
(  692)        END DO
(  693) 
(  694) 
(  695)        !----------------------------------------------------------------------------
(  696)        !  10. Surface accretion. 
(  697)        !      The change of ice thickness after snow/ice formation is such that
(  698)        !      the interface between snow and ice is located at the same height
(  699)        !      as the ocean surface. It is given by (Fichefet and Morales Maqueda 1999)
(  700)        !          D(h_ice) = (- D(hsn)/alph) =  [rhosn*hsn - (rau0 - rhoic)*hic]
(  701)        !                                     / [alph*rhosn+rau0 - rhoic]
(  702)        !----------------------------------------------------------------------------
(  703)        !
(  704)        DO ji = kideb , kiut
(  705) 
(  706)           !--  Computation of the change of ice thickness after snow-ice formation
(  707)           zdrmh =  ( rhosn * h_snow_1d(ji) + ( rhoic - rau0 ) * h_ice_1d(ji) )  &
(  708)              &  / ( alphs * rhosn + rau0 - rhoic )
(  709)           zdrmh = MAX( zzero , zdrmh )
(  710) 
(  711)           !--New ice and snow thicknesses Fichefet and Morales Maqueda (1999)
(  712)           zhicnew  = MAX( h_ice_1d(ji) , h_ice_1d(ji) + zdrmh )
(  713)           zhsnnew  = MIN( h_snow_1d(ji) , h_snow_1d(ji) - alphs * zdrmh )
(  714)           !---Compute new ice temperatures. snow temperature remains unchanged
(  715)           !   Lepparanta (1983):
(  716)           zihic = 1.0 - MAX( zzero , SIGN( zone , -zhicnew ) )
(  717)           zquot  = ( 1.0 - zihic ) &
(  718)              &   +         zihic * MIN( zone , h_ice_1d(ji) / MAX( zhicnew , epsi20 ) ) 
(  719)           ztneq  =         alphs * cnscg * tbif_1d(ji,1)    &
(  720)              &   + ( 1.0 - alphs * ( rhosn/rhoic ) ) * tfu_1d(ji)
(  721)           zep(ji) = tbif_1d(ji,2)
(  722)           tbif_1d(ji,2) = ztneq - zquot * zquot * ( ztneq - tbif_1d(ji,2) )
(  723)           tbif_1d(ji,3) = 2.0 * ztneq &
(  724)              &        + zquot * ( tbif_1d(ji,3) + zep(ji) - 2.0 * ztneq ) - tbif_1d(ji,2)
(  725) 
(  726)           !---  Lepparanta (1983) (latent heat released during white ice formation
(  727)           !     goes to the ocean -for lateral ablation-)
(  728)           qldif_1d(ji)  = qldif_1d(ji) + zdrmh * ( 1.0 - alphs * ( rhosn/rhoic ) ) * xlic * ( 1.0 - frld_1d(ji) )
(  729)           !--   Changes in ice volume and ice mass Lepparanta (1983):
(  730)           dvnbq_1d(ji) = ( 1.0 - frld_1d(ji) ) * ( zhicnew - h_ice_1d(ji) )
(  731)           dmgwi_1d(ji) = dmgwi_1d(ji) + ( 1.0 -frld_1d(ji) ) * ( h_snow_1d(ji) - zhsnnew ) * rhosn
(  732)           !  case of natural freshwater flux
(  733) #if defined key_lim_fdd  
(  734)           rdmicif_1d(ji) = rdmicif_1d(ji) + ( 1.0 - frld_1d(ji) )   * ( zhicnew - h_ice_1d(ji) )  * rhoic
(  735)           rdmsnif_1d(ji) = rdmsnif_1d(ji) + ( 1.0 - frld_1d(ji) )   * ( zhsnnew - h_snow_1d(ji) ) * rhosn
(  736) 
(  737) #else
(  738)           rdmicif_1d(ji) = rdmicif_1d(ji) + ( 1.0 - frld_1d(ji) ) * (  ( zhicnew - h_ice_1d(ji) ) * rhoic   &
(  739)              &                                                    + ( zhsnnew - h_snow_1d(ji) ) * rhosn )
(  740) #endif
(  741) 
(  742)           !---  Actualize new snow and ice thickness.






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 14

(  743) 
(  744)           h_snow_1d(ji)  = zhsnnew
(  745)           h_ice_1d(ji)  = zhicnew
(  746) 
(  747)        END DO
(  748) 
(  749)        !----------------------------------------------------
(  750)        !  11. Lateral ablation (Changes in sea/ice fraction) 
(  751)        !----------------------------------------------------
(  752)        DO ji = kideb , kiut
(  753)           zfrl_old(ji)  = frld_1d(ji)
(  754)           zihic   = 1.0 - MAX( zzero , SIGN( zone , -h_ice_1d(ji) ) )
(  755)           zihsn   = 1.0 - MAX( zzero , SIGN( zone , -h_snow_1d(ji) ) )
(  756)           !--In the case of total ablation (all the ice ice has melted) frld = 1
(  757)           frld_1d(ji)  = ( 1.0 - zihic ) + zihic * zfrl_old(ji)
(  758)           !--Part of solar radiation absorbing inside the ice and going
(  759)           !--through the ocean
(  760)           fscbq_1d(ji) = ( 1.0 - zfrl_old(ji) ) * ( 1.0 - thcm_1d(ji) ) * fstbif_1d(ji)
(  761)           !--Total remaining energy after bottom melting/growing 
(  762)           qfvbq_1d(ji) = zqsup(ji) + ( 1.0 - zihic ) * zqocea(ji)
(  763)           !--Updating of total heat from the ocean
(  764)           qldif_1d(ji)  = qldif_1d(ji) + qfvbq_1d(ji) + ( 1.0 - zihic ) * fscbq_1d(ji) * rdt_ice
(  765)           !--Computation of total heat inside the snow/ice system
(  766)           zqice  = h_snow_1d(ji) * xlsn + h_ice_1d(ji) * xlic
(  767)           zqicetot  = ( 1.0 - frld_1d(ji) ) * zqice
(  768)           !--The concentration of ice is reduced (frld increases) if the heat 
(  769)           !--exchange between ice and ocean is positive
(  770)           ziqf = MAX( zzero , SIGN( zone ,  zqicetot - qldif_1d(ji) ) )
(  771)           zdfrl = qldif_1d(ji) / MAX( epsi20 , zqice ) 
(  772)           frld_1d(ji)  = ( 1.0 - ziqf )    &
(  773)              &       +         ziqf * ( frld_1d(ji) + MAX( zzero , zdfrl ) ) 
(  774)           fltbif_1d(ji) = ( ( 1.0 - zfrl_old(ji) ) * qstbif_1d(ji) - zqicetot  ) / rdt_ice
(  775)           !--  Opening of leads: Hakkinen & Mellor, 1992.
(  776)           zdfrl = - ( zdhictop(ji) + zdhicbot(ji) ) * hakspl * ( 1.0 - zfrl_old(ji) ) &
(  777)              &  / MAX( epsi13 , h_ice_1d(ji) + h_snow_1d(ji) * rhosn/rhoic ) 
(  778)           zfrld_1d(ji) =  frld_1d(ji) + MAX( zzero , zdfrl )
(  779)           !--Limitation of sea-ice fraction <= 1
(  780)           zfrld_1d(ji) = ziqf * MIN( 0.99 * zone , zfrld_1d(ji) ) + ( 1 - ziqf )
(  781)           !---Update surface and internal temperature and snow/ice thicknesses.
(  782)           sist_1d(ji)   = sist_1d(ji)   + ( 1.0 - ziqf ) * ( tfu_1d(ji) - sist_1d(ji)   )
(  783)           tbif_1d(ji,1) = tbif_1d(ji,1) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,1) )
(  784)           tbif_1d(ji,2) = tbif_1d(ji,2) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,2) )
(  785)           tbif_1d(ji,3) = tbif_1d(ji,3) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,3) )
(  786)           !--variation of ice volume and ice mass 
(  787)           dvlbq_1d(ji)   = zihic * ( zfrl_old(ji) - frld_1d(ji) ) * h_ice_1d(ji)
(  788)           rdmicif_1d(ji) = rdmicif_1d(ji) + dvlbq_1d(ji) * rhoic
(  789)           !--variation of snow volume and snow mass 
(  790)           zdvsnvol    = zihsn * ( zfrl_old(ji) - frld_1d(ji) ) * h_snow_1d(ji)
(  791)           rdmsnif_1d(ji) = rdmsnif_1d(ji) + zdvsnvol * rhosn
(  792)           h_snow_1d(ji)  = ziqf * h_snow_1d(ji)
(  793) 
(  794)           zdrfrl1 = ziqf * ( 1.0 -  frld_1d(ji) ) / MAX( epsi20 , 1.0 - zfrld_1d(ji) )
(  795)           zdrfrl2 = ziqf * ( 1.0 - zfrl_old(ji) ) / MAX( epsi20 , 1.0 - zfrld_1d(ji) )
(  796) 
(  797)           h_snow_1d (ji) = zdrfrl1 * h_snow_1d(ji)
(  798)           h_ice_1d  (ji) = zdrfrl1 * h_ice_1d(ji)
(  799)           qstbif_1d(ji) = zdrfrl2 * qstbif_1d(ji)
(  800)           frld_1d(ji)    = zfrld_1d(ji)






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 15

(  801) 
(  802)        END DO
(  803)        
(  804)     END SUBROUTINE lim_thd_zdf
(  805) #else
(  806)    !!======================================================================
(  807)    !!                       ***  MODULE limthd_zdf   ***
(  808)    !!                              no sea ice model  
(  809)    !!======================================================================
(  810) CONTAINS
(  811)    SUBROUTINE lim_thd_zdf          ! Empty routine
(  812)    END SUBROUTINE lim_thd_zdf
(  813) #endif
(  814)  END MODULE limthd_zdf


















































PGF90 (Version     10.5)          09/23/2010  09:49:57      page 16

(    1) # 1 "limthd_zdf.F90"
(    1) MODULE limthd_zdf
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limthd_zdf ***
(    4)    !!                thermodynamic growth and decay of the ice 
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_thd_zdf  : vertical accr./abl. and lateral ablation of sea ice
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE par_oce          ! ocean parameters
(   14)    USE phycst           ! ???
(   15)    USE ice_oce          ! ice variables
(   16)    USE thd_ice
(   17)    USE iceini
(   18)    USE limistate
(   19)    USE in_out_manager
(   20)       
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_thd_zdf        ! called by lim_thd
(   26) 
(   27)    !! * Module variables
(   28)    REAL(wp)  ::           &  ! constant values
(   29)       epsi20 = 1.e-20  ,  &
(   30)       epsi13 = 1.e-13  ,  &
(   31)       zzero  = 0.e0    ,  &
(   32)       zone   = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd_zdf.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) CONTAINS
(   39) 
(   40)    SUBROUTINE lim_thd_zdf( kideb , kiut )
(   41)       !!------------------------------------------------------------------
(   42)       !!                ***  ROUTINE lim_thd_zdf  ***
(   43)       !!              
(   44)       !! ** Purpose : This routine determines the time evolution of snow 
(   45)       !!      and sea-ice thicknesses, concentration and heat content 
(   46)       !!      due to the vertical and lateral thermodynamic accretion-
(   47)       !!      ablation processes. One only treats the case of lat. abl.
(   48)       !!      For lateral accretion, see routine lim_lat_accr 
(   49)       !! 
(   50)       !! ** Method  : The representation of vertical growth and decay of 
(   51)       !!      the sea-ice model is based upon the diffusion of heat 
(   52)       !!      through the external and internal boundaries of a 
(   53)       !!      three-layer system (two layers of ice and one layer and 
(   54)       !!      one layer of snow, if present, on top of the ice).
(   55)       !! 
(   56)       !! ** Action  : - Calculation of some intermediates variables
(   57)       !!              - Calculation of surface temperature






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 17

(   58)       !!              - Calculation of available heat for surface ablation
(   59)       !!              - Calculation of the changes in internal temperature 
(   60)       !!                of the three-layer system, due to vertical diffusion 
(   61)       !!                processes
(   62)       !!              - Performs surface ablation and bottom accretion-ablation
(   63)       !!              - Performs snow-ice formation
(   64)       !!              - Performs lateral ablation
(   65)       !!
(   66)       !! References :
(   67)       !!   Fichefet T. and M. Maqueda 1997, J. Geophys. Res., 102(C6), 12609-12646   
(   68)       !!   Fichefet T. and M. Maqueda 1999, Clim. Dyn, 15(4), 251-268  
(   69)       !!
(   70)       !! History :
(   71)       !!   original    : 01-04 (LIM)
(   72)       !!   addition    : 02-08 (C. Ethe, G. Madec)
(   73)       !!------------------------------------------------------------------
(   74)       !! * Arguments
(   75)       INTEGER , INTENT (in) ::  &
(   76)          kideb ,  &  ! Start point on which the  the computation is applied
(   77)          kiut        ! End point on which the  the computation is applied
(   78) 
(   79)       !! * Local variables
(   80)       INTEGER ::   ji       ! dummy loop indices
(   81) 
(   82)       REAL(wp) , DIMENSION(jpij,2) ::  &
(   83)          zqcmlt        ! energy due to surface( /1 ) and bottom melting( /2 )
(   84) 
(   85)       REAL(wp), DIMENSION(jpij) ::  &
(   86)          ztsmlt      &    ! snow/ice surface melting temperature
(   87)          ,ztbif      &    ! int. temp. at the mid-point of the 1st layer of the snow/ice sys. 
(   88)          ,zksn       &    ! effective conductivity of snow
(   89)          ,zkic       &    ! effective conductivity of ice
(   90)          ,zksndh     &    ! thermal cond. at the mid-point of the 1st layer of the snow/ice sys. 
(   91)          , zfcsu     &    ! conductive heat flux at the surface of the snow/ice system 
(   92)          , zfcsudt   &    ! = zfcsu * dt
(   93)          , zi0       &    ! frac. of the net SW rad. which is not absorbed at the surface
(   94)          , z1mi0     &    ! fraction of the net SW radiation absorbed at the surface
(   95)          , zqmax     &    ! maximum energy stored in brine pockets
(   96)          , zrcpdt    &    ! h_su*rho_su*cp_su/dt(h_su being the thick. of surf. layer)
(   97)          , zts_old   &    ! previous surface temperature
(   98)          , zidsn , z1midsn , zidsnic ! tempory variables
(   99) 
(  100)       REAL(wp), DIMENSION(jpij) :: &
(  101)           zfnet       &  ! net heat flux at the top surface( incl. conductive heat flux)
(  102)           , zsprecip  &    ! snow accumulation
(  103)           , zhsnw_old &    ! previous snow thickness
(  104)           , zdhictop  &    ! change in ice thickness due to top surf ablation/accretion
(  105)           , zdhicbot  &    ! change in ice thickness due to bottom surf abl/acc
(  106)           , zqsup     &    ! energy transmitted to ocean (coming from top surf abl/acc)
(  107)           , zqocea    &    ! energy transmitted to ocean (coming from bottom sur abl/acc)
(  108)           , zfrl_old  &    ! previous sea/ice fraction
(  109)           , zfrld_1d    &    ! new sea/ice fraction
(  110)           , zep            ! internal temperature of the 2nd layer of the snow/ice system
(  111) 
(  112)        REAL(wp), DIMENSION(3) :: & 
(  113)           zplediag  &    ! principle diagonal, subdiag. and supdiag. of the 
(  114)           , zsubdiag  &    ! tri-diagonal matrix coming from the computation
(  115)           , zsupdiag  &    ! of the temperatures inside the snow-ice system






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 18

(  116)           , zsmbr          ! second member
(  117) 
(  118)        REAL(wp) :: & 
(  119)           zhsu     &     ! thickness of surface layer
(  120)           , zhe      &     ! effective thickness for compu. of equ. thermal conductivity
(  121)           , zheshth  &     ! = zhe / thth
(  122)           , zghe     &     ! correction factor of the thermal conductivity
(  123)           , zumsb    &     ! parameter for numerical method to solve heat-diffusion eq.
(  124)           , zkhsn    &     ! conductivity at the snow layer
(  125)           , zkhic    &     ! conductivity at the ice layers
(  126)           , zkint    &     ! equivalent conductivity at the snow-ice interface
(  127)           , zkhsnint &     ! = zkint*dt / (hsn*rhosn*cpsn)  
(  128)           , zkhicint &     ! = 2*zkint*dt / (hic*rhoic*cpic)
(  129)           , zpiv1 , zpiv2  &       ! tempory scalars used to solve the tri-diagonal system
(  130)           , zb2 , zd2 , zb3 , zd3 &
(  131)           , ztint          ! equivalent temperature at the snow-ice interface
(  132) 
(  133)        REAL(wp) :: & 
(  134)           zexp      &     ! exponential function of the ice thickness
(  135)           , zfsab     &     ! part of solar radiation stored in brine pockets
(  136)           , zfts      &     ! value of energy balance function when the temp. equal surf. temp.
(  137)           , zdfts     &     ! value of derivative of ztfs when the temp. equal surf. temp.
(  138)           , zdts      &     ! surface temperature increment
(  139)           , zqsnw_mlt &     ! energy needed to melt snow
(  140)           , zdhsmlt   &     ! change in snow thickness due to melt
(  141)           , zhsn      &     ! snow thickness (previous+accumulation-melt)
(  142)           , zqsn_mlt_rem &  ! remaining heat coming from snow melting
(  143)           , zqice_top_mlt & ! energy used to melt ice at top surface
(  144)           , zdhssub      &  ! change in snow thick. due to sublimation or evaporation
(  145)           , zdhisub      &  ! change in ice thick. due to sublimation or evaporation    
(  146)           , zdhsn        &  ! snow ice thickness increment
(  147)           , zdtsn        &  ! snow internal temp. increment
(  148)           , zdtic        &  ! ice internal temp. increment
(  149)           , zqnes          ! conductive energy due to ice melting in the first ice layer
(  150) 
(  151)        REAL(wp) :: & 
(  152)           ztbot     &      ! temperature at the bottom surface
(  153)           , zfcbot    &      ! conductive heat flux at bottom surface
(  154)           , zqice_bot &      ! energy used for bottom melting/growing
(  155)           , zqice_bot_mlt &  ! energy used for bottom melting
(  156)           , zqstbif_bot  &  ! part of energy stored in brine pockets used for bottom melting
(  157)           , zqstbif_old  &  ! tempory var. for zqstbif_bot
(  158)           , zdhicmlt      &  ! change in ice thickness due to bottom melting
(  159)           , zdhicm        &  ! change in ice thickness var. 
(  160)           , zdhsnm        &  ! change in snow thickness var. 
(  161)           , zhsnfi        &  ! snow thickness var. 
(  162)           , zc1, zpc1, zc2, zpc2, zp1, zp2 & ! tempory variables
(  163)           , ztb2, ztb3
(  164) 
(  165)        REAL(wp) :: & 
(  166)           zdrmh         &   ! change in snow/ice thick. after snow-ice formation
(  167)           , zhicnew       &   ! new ice thickness
(  168)           , zhsnnew       &   ! new snow thickness
(  169)           , zquot , ztneq &   ! tempory temp. variables
(  170)           , zqice, zqicetot & ! total heat inside the snow/ice system
(  171)           , zdfrl         &   ! change in ice concentration
(  172)           , zdvsnvol      &   ! change in snow volume
(  173)           , zdrfrl1, zdrfrl2 &  ! tempory scalars






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 19

(  174)           , zihsn, zidhb, zihic, zihe, zihq, ziexp, ziqf, zihnf, zibmlt, ziqr, zihgnew, zind
(  175)        !!----------------------------------------------------------------------
(  176) 
(  177)        !-----------------------------------------------------------------------
(  178)        !  1. Boundaries conditions for snow/ice system internal temperature
(  179)        !       - If tbif_1d(ji,1) > rt0_snow, tbif_1d(ji,1) = rt0_snow 
(  180)        !       - If tbif_1d(ji,2/3) > rt0_ice, tbif_1d(ji,2/3) = rt0_ice 
(  181)        !     Computation of energies due to surface and bottom melting 
(  182)        !-----------------------------------------------------------------------
(  183)        
(  184)        
(  185)        DO ji = kideb , kiut
(  186)           zihsn = MAX( zzero , SIGN( zone , hsndif - h_snow_1d(ji) ) )
(  187)           zihic = MAX( zzero , SIGN( zone , hicdif - h_ice_1d(ji) ) )
(  188)           !--computation of energy due to surface melting
(  189)           zqcmlt(ji,1) = ( MAX ( zzero ,  &
(  190)              &                   rcpsn * h_snow_1d(ji) * ( tbif_1d(ji,1) - rt0_snow ) ) ) * ( 1.0 - zihsn )
(  191)           !--computation of energy due to bottom melting
(  192)           zqcmlt(ji,2) = ( MAX( zzero , &
(  193)              &                  rcpic * ( tbif_1d(ji,2) - rt0_ice ) * ( h_ice_1d(ji) / 2. ) ) &
(  194)              &           + MAX( zzero , &
(  195)              &                  rcpic * ( tbif_1d(ji,3) - rt0_ice ) * ( h_ice_1d(ji) / 2. ) ) &
(  196)              &           ) * ( 1.0 - zihic  )
(  197)           !--limitation of  snow/ice system internal temperature
(  198)           tbif_1d(ji,1)   = MIN( rt0_snow, tbif_1d(ji,1) )
(  199)           tbif_1d(ji,2)   = MIN( rt0_ice , tbif_1d(ji,2) )
(  200)           tbif_1d(ji,3)   = MIN( rt0_ice , tbif_1d(ji,3) )
(  201)        END DO
(  202) 
(  203) 
(  204)        !-------------------------------------------
(  205)        !  2. Calculate some intermediate variables.  
(  206)        !-------------------------------------------
(  207)        
(  208)        ! initialisation of the thickness of surface layer
(  209)        zhsu = hnzst  
(  210) 
(  211)        DO ji = kideb , kiut
(  212)           zind   = MAX( zzero , SIGN( zone , zhsu - h_snow_1d(ji) ) )
(  213)           zihsn  = MAX( zzero , SIGN( zone , hsndif - h_snow_1d(ji) ) )
(  214)           zihsn  = MAX( zihsn , zind )
(  215)           zihic  = MAX( zzero , sign( zone , hicdif - h_ice_1d(ji) ) )
(  216)           !     2.1. Computation of surface melting temperature
(  217)           !----------------------------------------------------
(  218)           zind  = MAX( zzero , SIGN( zone , -h_snow_1d(ji) ) )
(  219)           ztsmlt(ji) = ( 1.0 - zind ) * rt0_snow + zind * rt0_ice
(  220)           !
(  221)           !     2.2. Effective conductivity of snow and ice
(  222)           !-----------------------------------------------
(  223) 
(  224)           !---computation of the correction factor on the thermal conductivity
(  225)           !-- (Morales Maqueda, 1995 ; Fichefet and Morales Maqueda, 1997)
(  226)           zhe      =  ( rcdsn / ( rcdsn + rcdic ) ) * h_ice_1d(ji)   &
(  227)              &     + ( rcdic / ( rcdsn + rcdic ) ) * h_snow_1d(ji) 
(  228)           zihe     = MAX( zzero , SIGN( zone , 2.0 * zhe - thth ) )
(  229)           zheshth  = zhe / thth
(  230)           zghe     = ( 1.0 - zihe ) * zheshth * ( 2.0 - zheshth )   &
(  231)              &     +         zihe   * 0.5 * ( 1.5 + LOG( 2.0 * zheshth ) )






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 20

(  232) 
(  233) # 236
(  236)           !---effective conductivities 
(  237)           zksn(ji)  = zghe * rcdsn  
(  238)           zkic(ji)  = zghe * rcdic
(  239) 
(  240)           !
(  241)           !     2.3. Computation of the conductive heat flux from the snow/ice 
(  242)           !          system interior toward the top surface
(  243)           !------------------------------------------------------------------
(  244) 
(  245)           !---Thermal conductivity at the mid-point of the first snow/ice system layer
(  246)           zksndh(ji) =   ( ( 1.0 - zihsn ) * 2.0 * zksn(ji) + zihsn * 4.0 * zkic(ji) )   &
(  247)              &         / ( ( 1.0 - zihsn ) *  h_snow_1d(ji)                              &
(  248)              &           +        zihsn   *  ( ( 1.0 + 3.0 * zihic ) * h_ice_1d(ji)      &
(  249)              &           + 4.0 * zkic(ji)/zksn(ji) * h_snow_1d(ji) ) )
(  250) 
(  251)           !---internal temperature at the mid-point of the first snow/ice system layer
(  252)           ztbif(ji)  = ( 1.0 - zihsn ) * tbif_1d(ji,1)                       &
(  253)              &       +         zihsn   * ( ( 1.0 - zihic ) * tbif_1d(ji,2)   &
(  254)              &       +         zihic   * tfu_1d(ji)   )
(  255)           !---conductive heat flux 
(  256)           zfcsu(ji) = zksndh(ji) * ( ztbif(ji) - sist_1d(ji) )
(  257) 
(  258)        END DO
(  259) 
(  260)        !--------------------------------------------------------------------
(  261)        !  3. Calculate : 
(  262)        !     - fstbif_1d, part of solar radiation absorbing inside the ice
(  263)        !       assuming an exponential absorption (Grenfell and Maykut, 1977)
(  264)        !     - zqmax,  maximum energy stored in brine pockets
(  265)        !     - qstbif_1d, total energy stored in brine pockets (updating)
(  266)        !-------------------------------------------------------------------
(  267) 
(  268) 
(  269)        DO ji = kideb , kiut
(  270)           zihsn  = MAX( zzero , SIGN (zone , -h_snow_1d(ji) ) )
(  271)           zihic  = MAX( zzero , 1.0 - ( h_ice_1d(ji) / zhsu ) )     
(  272)           zind   = MAX( zzero , SIGN (zone , hicdif - h_ice_1d(ji) ) )
(  273)           !--Computation of the fraction of the net shortwave radiation which
(  274)           !--penetrates inside the ice cover ( See Forcat)
(  275)           zi0(ji)  = zihsn * ( fr1_i0_1d(ji) + zihic * fr2_i0_1d(ji) )
(  276)           zexp     = MIN( zone , EXP( -1.5 * ( h_ice_1d(ji) - zhsu ) ) )
(  277)           fstbif_1d(ji) = zi0(ji) * qsr_ice_1d(ji) * zexp
(  278)           !--Computation of maximum energy stored in brine pockets zqmax and update
(  279)           !--the total energy stored in brine pockets, if less than zqmax
(  280)           zqmax(ji) = MAX( zzero , 0.5 * xlic * ( h_ice_1d(ji) - hicmin ) )
(  281)           zfsab   = zi0(ji) * qsr_ice_1d(ji) * ( 1.0 - zexp )
(  282)           zihq    = ( 1.0 - zind ) * MAX(zzero, SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) ) &
(  283)              &    +         zind   * zone
(  284)           qstbif_1d(ji) = ( qstbif_1d(ji) + ( 1.0 - zihq ) * zfsab * rdt_ice ) * swiqst
(  285)           !--fraction of shortwave radiation absorbed at surface
(  286)           ziexp = zihq * zexp + ( 1.0 - zihq ) * ( swiqst + ( 1.0 - swiqst ) * zexp )
(  287)           z1mi0(ji) = 1.0 - zi0(ji) * ziexp
(  288)        END DO
(  289) 
(  290) 
(  291)        !--------------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 21

(  292)        !  4. Computation of the surface temperature : determined by considering the 
(  293)        !     budget of a thin layer of thick. zhsu at the top surface (H. Grenier, 1995) 
(  294)        !     and based on a surface energy balance : 
(  295)        !     hsu * rcp * dT/dt = Fsr + Fnsr(T) + Fcs(T),
(  296)        !     where - Fsr is the net absorbed solar radiation, 
(  297)        !           - Fnsr is the total non solar radiation (incoming and outgoing long-wave,
(  298)        !             sensible and latent heat fluxes)
(  299)        !           - Fcs the conductive heat flux at the top of surface
(  300)        !------------------------------------------------------------------------------
(  301) 
(  302)        !     4.1. Computation of intermediate values
(  303)        !---------------------------------------------
(  304)        DO ji = kideb, kiut
(  305)           zrcpdt(ji) = ( rcpsn * MIN( h_snow_1d(ji) , zhsu )    &
(  306)              &       + rcpic * MAX( zhsu - h_snow_1d(ji) , zzero ) ) / rdt_ice
(  307)           zts_old(ji) =  sist_1d(ji)
(  308)        END DO
(  309) 
(  310)        !     4.2. Computation of surface temperature by expanding the eq. of energy balance
(  311)        !          with Ts = Tp + DT. One obtain , F(Tp) + DT * DF(Tp) = 0
(  312)        !          where  - F(Tp) = Fsr + Fnsr(Tp) + Fcs(Tp) 
(  313)        !                 - DF(Tp)= (dFnsr(Tp)/dT) + (dFcs(Tp)/dT) - hsu*rcp/dt
(  314)        !---------------------------------------------------------------------------------
(  315) 
(  316)        DO ji = kideb, kiut
(  317)           !---computation of the derivative of energy balance function 
(  318) # 329
(  329)           zdfts    =  zksndh(ji)   & ! contribution of the conductive heat flux
(  330)              &      + zrcpdt(ji)   & ! contribution of hsu * rcp / dt
(  331)              &      - dqns_ice_1d (ji)     ! contribution of the total non solar radiation 
(  332) # 333
(  333)           !---computation of the energy balance function 
(  334)           zfts    = - z1mi0 (ji) * qsr_ice_1d(ji)   & ! net absorbed solar radiation
(  335)              &      - qnsr_ice_1d(ji)                & ! total non solar radiation
(  336)              &      - zfcsu (ji)                  ! conductive heat flux from the surface
(  337)           !---computation of surface temperature increment  
(  338)           zdts    = -zfts / zdfts
(  339) # 342
(  342)           !---computation of the new surface temperature 
(  343)           sist_1d(ji) = sist_1d(ji) + zdts
(  344) 
(  345)        END DO
(  346) 
(  347)        !----------------------------------------------------------------------------
(  348)        !  5. Boundary condition at the top surface
(  349)        !--    IF Tsb < Tmelt, Fnet = Fcs (the net heat flux equal the conductive heat flux)
(  350)        !      Otherwise Tsb = Tmelt and Qnet(Tmelt) > 0 
(  351)        !      Fnet(Tmelt) is therefore the net surface flux needed for melting
(  352)        !----------------------------------------------------------------------------
(  353)        
(  354)        
(  355)        !     5.1.  Limitation of surface temperature and update total non solar fluxes,
(  356)        !          latent heat flux and conductive flux at the top surface 
(  357)        !----------------------------------------------------------------------  
(  358)                      
(  359)        DO ji = kideb, kiut
(  360)           sist_1d(ji) = MIN( ztsmlt(ji) , sist_1d(ji) )
(  361) # 362






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 22

(  362)           qnsr_ice_1d(ji) = qnsr_ice_1d(ji) + dqns_ice_1d(ji) * ( sist_1d(ji) - zts_old(ji) )
(  363)           qla_ice_1d (ji) = qla_ice_1d (ji) + dqla_ice_1d(ji) * ( sist_1d(ji) - zts_old(ji) )
(  364) # 365
(  365)           zfcsu(ji)  = zksndh(ji) * ( ztbif(ji) - sist_1d(ji) )
(  366)        END DO
(  367) 
(  368)      
(  369) 
(  370)        !     5.2. Calculate available heat for surface ablation. 
(  371)        !---------------------------------------------------------------------
(  372) 
(  373)        DO ji = kideb, kiut
(  374)           zfnet(ji) = qnsr_ice_1d(ji) + z1mi0(ji) * qsr_ice_1d(ji) + zfcsu(ji)          
(  375)           zfnet(ji) = MAX( zzero , zfnet(ji) )
(  376)           zfnet(ji) = zfnet(ji) * MAX( zzero , SIGN( zone , sist_1d(ji) - ztsmlt(ji) ) )
(  377)        END DO
(  378) 
(  379)        !-------------------------------------------------------------------------
(  380)        !  6. Calculate changes in internal temperature due to vertical diffusion   
(  381)        !     processes. The evolution of this temperature is governed by the one-
(  382)        !     dimensionnal heat-diffusion equation. 
(  383)        !     Given the temperature tbif(1/2/3), at time m we solve a set
(  384)        !     of finite difference equations to obtain new tempe. Each tempe is coupled
(  385)        !     to the temp. immediatly above and below by heat conduction terms. Thus 
(  386)        !     we have a set of equations of the form A * T = B, where A is a tridiagonal
(  387)        !     matrix, T a vector whose components are the unknown new temp.
(  388)        !-------------------------------------------------------------------------
(  389)        
(  390)        !--parameter for the numerical methode use to solve the heat-diffusion equation
(  391)        !- implicit, explicit or Crank-Nicholson
(  392)        zumsb = 1.0 - sbeta  
(  393)        DO ji = kideb, kiut
(  394)           zidsn(ji)   = MAX ( zzero, SIGN( zone, hsndif - h_snow_1d(ji) ) ) 
(  395)           z1midsn(ji) = 1.0 - zidsn(ji)
(  396)           zihic       = MAX ( zzero, SIGN( zone, hicdif - h_ice_1d(ji) ) ) 
(  397)           zidsnic(ji) = zidsn(ji) *  zihic 
(  398)           zfcsudt(ji) = zfcsu(ji) * rdt_ice 
(  399)        END DO 
(  400)    
(  401)        DO ji = kideb, kiut
(  402) 
(  403)           !     6.1 Calculate intermediate variables.
(  404)           !----------------------------------------
(  405) 
(  406)           !--conductivity at the snow surface
(  407)           zkhsn = 2.0 * zksn(ji) * rdt_ice / rcpsn
(  408)           !--conductivity at the ice surface
(  409)           zkhic = 4.0 * zkic(ji) * rdt_ice / MAX( h_ice_1d(ji) * h_ice_1d(ji) * rcpic , epsi20 )
(  410)           !--conductivity at the snow/ice interface 
(  411)           zkint = 4.0 * zksn(ji) * zkic(ji)  &
(  412)              &        / ( zksn(ji) * h_ice_1d(ji) + 2.0 * zkic(ji) * h_snow_1d(ji) * z1midsn(ji)) 
(  413)           zkhsnint = zkint * rdt_ice / rcpsn
(  414)           zkhicint = zkint * 2.0 * rdt_ice / MAX( h_ice_1d(ji) * rcpic , epsi20 )
(  415)           
(  416)           !     6.2. Fulfill the linear system matrix.
(  417)           !-----------------------------------------
(  418) !$$$          zplediag(1) = 1 + sbeta * z1midsn(ji) * ( zkhsn + zkhsnint )       
(  419)           zplediag(1) =   zidsn(ji) + z1midsn(ji) * h_snow_1d(ji)   &






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 23

(  420)              &          + sbeta * z1midsn(ji) * zkhsnint 
(  421)           zplediag(2) = 1 + sbeta * ( z1midsn(ji) * zkhicint + zkhic ) 
(  422)           zplediag(3) = 1 + 3.0 * sbeta * zkhic   
(  423) 
(  424)           zsubdiag(1) =  0.e0              
(  425)           zsubdiag(2) = -1.e0 * z1midsn(ji) * sbeta * zkhicint
(  426)           zsubdiag(3) = -1.e0 * sbeta * zkhic 
(  427) 
(  428)           zsupdiag(1) = -1.e0 * z1midsn(ji) * sbeta * zkhsnint 
(  429)           zsupdiag(2) = zsubdiag(3)
(  430)           zsupdiag(3) =  0.e0
(  431)           
(  432)           !     6.3. Fulfill the idependent term vector.
(  433)           !-------------------------------------------
(  434)           
(  435) !$$$          zsmbr(1) = zidsn(ji) * sist_1d(ji) + z1midsn(ji) *   &
(  436) !$$$             &         ( tbif_1d(ji,1) + zkhsn * sist_1d(ji)
(  437) !$$$             &         - zumsb * ( zkhsn * tbif_1d(ji,1)
(  438) !$$$             &                   + zkhsnint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) ) )
(  439)           zsmbr(1) = zidsn(ji) * sist_1d(ji) + z1midsn(ji) *    &
(  440)              &       ( h_snow_1d(ji) * tbif_1d(ji,1) - ( zfcsudt(ji) / rcpsn )  &
(  441)              &       - zumsb * zkhsnint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) )
(  442) 
(  443)           zsmbr(2) =  tbif_1d(ji,2)  &
(  444)              &      - zidsn(ji) * ( 1.0 - zidsnic(ji) ) &
(  445)              &        * ( zfcsudt(ji) / MAX( h_ice_1d(ji) * rcpic , epsi20 ) ) &
(  446)              &      + zumsb * ( zkhicint * ( tbif_1d(ji,1) - tbif_1d(ji,2) ) &
(  447)              &                   - zkhic * ( tbif_1d(ji,2) - tbif_1d(ji,3) )  )
(  448) 
(  449)           zsmbr(3) =  tbif_1d(ji,3)  &
(  450)              &      + zkhic * ( 2.0 * tfu_1d(ji) &
(  451)              &                + zumsb * ( tbif_1d(ji,2) - 3.0 * tbif_1d(ji,3) ) )
(  452)           
(  453)           !     6.4. Solve the system (Gauss elimination method).
(  454)           !----------------------------------------------------
(  455)           
(  456)           zpiv1 = zsubdiag(2) / zplediag(1) 
(  457)           zb2   = zplediag(2) - zpiv1 * zsupdiag(1)
(  458)           zd2   = zsmbr(2) - zpiv1 * zsmbr(1)
(  459) 
(  460)           zpiv2 = zsubdiag(3) / zb2
(  461)           zb3   = zplediag(3) - zpiv2 * zsupdiag(2)
(  462)           zd3   = zsmbr(3) - zpiv2 * zd2
(  463) 
(  464)           tbif_1d(ji,3) = zd3 / zb3
(  465)           tbif_1d(ji,2) = ( zd2 - zsupdiag(2) * tbif_1d(ji,3) ) / zb2
(  466)           tbif_1d(ji,1) = ( zsmbr(1) - zsupdiag(1) * tbif_1d(ji,2) ) / zplediag(1)            
(  467) 
(  468)           !--- taking into account the particular case of  zidsnic(ji) = 1
(  469)           ztint =  (  zkic(ji) * h_snow_1d(ji) * tfu_1d (ji)    &
(  470)              &      + zksn(ji) * h_ice_1d(ji) * sist_1d(ji) )   &
(  471)              &   / ( zkic(ji) * h_snow_1d(ji) + zksn(ji) * h_ice_1d(ji) ) 
(  472) 
(  473)           tbif_1d(ji,1) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,1)   &
(  474)              &                + zidsnic(ji)   * ( ztint + sist_1d(ji) ) / 2.0
(  475)           tbif_1d(ji,2) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,2)   &
(  476)              &                + zidsnic(ji)   * ( 3.0 * ztint + tfu_1d(ji) ) / 4.0
(  477)           tbif_1d(ji,3) = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,3)   &






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 24

(  478)              &                + zidsnic(ji)   * ( ztint + 3.0 * tfu_1d(ji) ) / 4.0     
(  479)        END DO
(  480)  
(  481)        !----------------------------------------------------------------------
(  482)        !  9. Take into account surface ablation and bottom accretion-ablation.|
(  483)        !----------------------------------------------------------------------
(  484)        
(  485)        !---Snow accumulation in one thermodynamic time step
(  486)        zsprecip(kideb:kiut) = sprecip_1d(kideb:kiut) * rdt_ice / rhosn
(  487) 
(  488) 
(  489)        DO ji = kideb, kiut
(  490)           
(  491)           !      9.1. Surface ablation and update of snow thickness and qstbif_1d
(  492)           !--------------------------------------------------------------------
(  493)           
(  494)           !--------------------------------------------------------------------------
(  495)           !--      Melting snow processes :
(  496)           !--      Melt at the upper surface is computed from the difference between 
(  497)           !--      the net heat flux (including the conductive heat flux) at the upper 
(  498)           !--      surface and the pre-existing energy due to surface melting
(  499)           !------------------------------------------------------------------------------
(  500)           
(  501)           !-- store the snow thickness
(  502)           zhsnw_old(ji) =  h_snow_1d(ji)
(  503)           !--computation of the energy needed to melt snow
(  504)           zqsnw_mlt  = zfnet(ji) * rdt_ice - zqcmlt(ji,1)
(  505)           !--change in snow thickness due to melt
(  506)           zdhsmlt = - zqsnw_mlt / xlsn
(  507)           
(  508)           !-- compute new snow thickness, taking into account the part of snow accumulation
(  509)           !   (as snow precipitation) and the part of snow lost due to melt
(  510)           zhsn =  h_snow_1d(ji) + zsprecip(ji) + zdhsmlt
(  511)           h_snow_1d(ji) = MAX( zzero , zhsn )
(  512)           !-- compute the volume of snow lost after surface melting and the associated mass
(  513)           dvsbq_1d(ji) =  ( 1.0 - frld_1d(ji) ) * ( h_snow_1d(ji) - zhsnw_old(ji) - zsprecip(ji) )
(  514)           dvsbq_1d(ji) =  MIN( zzero , dvsbq_1d(ji) )
(  515)           rdmsnif_1d(ji) =  rhosn * dvsbq_1d(ji)
(  516)           !-- If the snow is completely melted the remaining heat is used to melt ice
(  517)           zqsn_mlt_rem  = MAX( zzero , -zhsn ) * xlsn
(  518)           zqice_top_mlt = zqsn_mlt_rem 
(  519)           zqstbif_old   = qstbif_1d(ji)
(  520) 
(  521)           !--------------------------------------------------------------------------
(  522)           !--      Melting ice processes at the top surface :
(  523)           !--      The energy used to melt ice, zqice_top_mlt, is taken from the energy
(  524)           !--      stored in brine pockets qstbif_1d and the remaining energy coming 
(  525)           !--      from the melting snow process zqsn_mlt_rem.
(  526)           !--      If qstbif_1d > zqsn_mlt_rem then, one uses only a zqsn_mlt_rem part
(  527)           !--      of qstbif_1d to melt ice,
(  528)           !--         zqice_top_mlt = zqice_top_mlt + zqsn_mlt_rem
(  529)           !--         qstbif_1d = qstbif_1d - zqsn_mlt_rem
(  530)           !--      Otherwise one uses all qstbif_1d to melt ice
(  531)           !--         zqice_top_mlt = zqice_top_mlt + qstbif_1d
(  532)           !--         qstbif_1d = 0
(  533)           !------------------------------------------------------
(  534)           
(  535)           ziqf =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqsn_mlt_rem  ) )






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 25

(  536)           zqice_top_mlt =         ziqf   * ( zqice_top_mlt + zqsn_mlt_rem )   &
(  537)              &          + ( 1.0 - ziqf ) * ( zqice_top_mlt + qstbif_1d(ji)  )
(  538) 
(  539)           qstbif_1d(ji) =         ziqf   * ( qstbif_1d(ji) - zqsn_mlt_rem )   &
(  540)              &        + ( 1.0 - ziqf ) * ( qstbif_1d(ji) - qstbif_1d(ji)  )
(  541) 
(  542)           !--    The contribution of the energy stored in brine pockets qstbif_1d to melt
(  543)           !--    ice is taking into account only when qstbif_1d is less than zqmax. 
(  544)           !--    Otherwise, only the remaining energy coming from the melting snow 
(  545)           !--    process is used 
(  546)           zihq =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) )
(  547) 
(  548)           zqice_top_mlt =         zihq   * zqice_top_mlt   &
(  549)              &          + ( 1.0 - zihq ) * zqsn_mlt_rem
(  550) 
(  551)           qstbif_1d(ji) =         zihq   * qstbif_1d(ji)   &
(  552)              &        + ( 1.0 - zihq ) * zqstbif_old
(  553) 
(  554)           !--change in ice thickness due to melt at the top surface
(  555)           zdhictop(ji) = -zqice_top_mlt / xlic
(  556)           !--compute the volume formed after surface melting
(  557)           dvsbq_1d(ji) =  zdhictop(ji) * ( 1.0 - frld_1d(ji) )
(  558) 
(  559)           !-------------------------------------------------------------------------
(  560)           !--      A small variation at the surface also occurs because of sublimation
(  561)           !--      associated with the latent flux. If qla_ice_1d is negative, snow condensates at 
(  562)           !        the surface. Otherwise, snow evaporates
(  563)           !-----------------------------------------------------------------------
(  564)           !----change in snow and ice thicknesses due to sublimation or evaporation
(  565)           zdhssub  = parsub * ( qla_ice_1d(ji) / ( rhosn * xsn ) ) * rdt_ice 
(  566)           zhsn     = h_snow_1d(ji) - zdhssub
(  567)           zdhisub  = MAX( zzero , -zhsn ) * rhosn/rhoic
(  568)           zdhictop(ji) =  zdhictop(ji) - zdhisub
(  569)           h_snow_1d(ji)  =  MAX( zzero , zhsn )
(  570)           !-------------------------------------------------
(  571)           !--  Update Internal temperature and qstbif_1d.
(  572)           !-------------------------------------------
(  573)           zihsn  =  MAX( zzero , SIGN( zone, -h_snow_1d(ji) ) )
(  574)           tbif_1d(ji,1) = ( 1.0 - zihsn ) * tbif_1d(ji,1) + zihsn   * tfu_1d(ji)
(  575)           !--change in snow internal temperature if snow has increased
(  576)           zihnf = MAX( zzero , SIGN( zone , h_snow_1d(ji) - zhsnw_old(ji) ) )
(  577)           zdhsn = 1.0 - zhsnw_old(ji) / MAX( h_snow_1d(ji) , epsi20 )
(  578)           zdtsn = zdhsn * ( sist_1d(ji) - tbif_1d(ji,1) )
(  579)           tbif_1d(ji,1) = tbif_1d(ji,1) + z1midsn(ji) * zihnf * zdtsn
(  580)           !--energy created due to ice melting in the first ice layer
(  581)           zqnes  = ( rt0_ice - tbif_1d(ji,2) ) * rcpic * ( h_ice_1d(ji) / 2. )
(  582)           !--change in first ice layer internal temperature
(  583)           ziqr  = MAX( zzero , SIGN( zone , qstbif_1d(ji) - zqnes ) )
(  584)           zdtic = qstbif_1d(ji) / ( rcpic * ( h_ice_1d(ji) / 2. ) )
(  585)           tbif_1d(ji,2) =  ziqr * rt0_ice + ( 1 - ziqr ) * ( tbif_1d(ji,2) + zdtic )
(  586)           !--update qstbif_1d
(  587)           qstbif_1d(ji) = ziqr * ( qstbif_1d(ji) - zqnes ) * swiqst
(  588) 
(  589) 
(  590)           !--      9.2. Calculate bottom accretion-ablation and update qstbif_1d.
(  591)           !             Growth and melting at bottom ice surface are governed by  
(  592)           !                 -xlic * Dh = (Fcb - Fbot ) * Dt 
(  593)           !             where Fbot is the net downward heat flux from ice to the ocean






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 26

(  594)           !            and Fcb is the conductive heat flux at the bottom surface
(  595)           !---------------------------------------------------------------------------
(  596)           ztbot = ( 1.0 - zidsnic(ji) ) * tbif_1d(ji,3) + zidsnic(ji) * sist_1d(ji)
(  597)           !---computes conductive heat flux at bottom surface 
(  598)           zfcbot =  4.0 * zkic(ji) * ( tfu_1d(ji) - ztbot )   &
(  599)              &   / ( h_ice_1d(ji) + zidsnic(ji) * ( 3. * h_ice_1d(ji) &
(  600)              &   + 4.0 * zkic(ji)/zksn(ji) * h_snow_1d(ji) ) )
(  601)           !---computation of net energy needed for bottom melting/growing
(  602)           zqice_bot = ( zfcbot - ( fbif_1d(ji) + qlbbq_1d(ji) ) ) * rdt_ice
(  603)           zqstbif_bot = qstbif_1d(ji)
(  604)           !---switch to know if bottom surface melts ( = 1 ) or grows ( = 0 )occurs
(  605)           zibmlt = MAX( zzero , SIGN( zone , -zqice_bot ) )
(  606)           !--particular case of melting (in the same way as the top surface)
(  607)           zqice_bot_mlt = zqice_bot 
(  608)           zqstbif_old = zqstbif_bot
(  609) 
(  610)           ziqf =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) + zqice_bot_mlt  ) )
(  611)           zqice_bot_mlt =         ziqf   * ( zqice_bot_mlt + zqice_bot_mlt ) &
(  612)              &          + ( 1.0 - ziqf ) * ( zqice_bot_mlt + qstbif_1d(ji)  )          
(  613)           qstbif_1d(ji)   =         ziqf   * ( qstbif_1d(ji) + zqice_bot_mlt ) &
(  614)              &          + ( 1.0 - ziqf ) * ( qstbif_1d(ji) - qstbif_1d(ji)  )
(  615)           !--    The contribution of the energy stored in brine pockets qstbif_1d to melt
(  616)           !--    ice is taking into account only when qstbif_1d is less than zqmax. 
(  617)           zihq =  MAX ( zzero , SIGN( zone , qstbif_1d(ji) - zqmax(ji) ) )
(  618)           zqice_bot_mlt =         zihq   * zqice_bot_mlt   &
(  619)              &          + ( 1.0 - zihq ) * zqice_bot
(  620)           qstbif_1d(ji)   =         zihq   * qstbif_1d(ji)   &
(  621)              &             + ( 1.0 - zihq ) * zqstbif_old
(  622) 
(  623)           !---treatment of the case of melting/growing
(  624)           zqice_bot   =         zibmlt   * ( zqice_bot_mlt - zqcmlt(ji,2) )   &
(  625)              &        + ( 1.0 - zibmlt ) * ( zqice_bot - zqcmlt(ji,2)  )
(  626)           qstbif_1d(ji) =         zibmlt   * qstbif_1d(ji)   &
(  627)              &           + ( 1.0 - zibmlt ) * zqstbif_bot
(  628) 
(  629)           !--computes change in ice thickness due to melt or growth
(  630)           zdhicbot(ji) = zqice_bot / xlic
(  631)           !--limitation of bottom melting if so : hmelt maximum melting at bottom
(  632)           zdhicmlt  = MAX( hmelt , zdhicbot(ji) ) 
(  633)           !--energy after bottom melting/growing
(  634)           zqsup(ji) = ( 1.0 - frld_1d(ji) ) * xlic * ( zdhicmlt - zdhicbot(ji) )
(  635)           !-- compute the new thickness and the newly formed volume after bottom melting/growing
(  636)           zdhicbot(ji)  = zdhicmlt
(  637)           dvbbq_1d(ji) = ( 1.0 - frld_1d(ji) ) * zdhicbot(ji)
(  638) 
(  639) 
(  640)           !        9.3.  Updating ice thickness after top surface ablation 
(  641)           !              and bottom surface accretion/ablation
(  642)           !---------------------------------------------------------------
(  643)           zhicnew  = h_ice_1d(ji) + zdhictop(ji) + zdhicbot(ji)
(  644) 
(  645)           !
(  646)           !        9.4. Case of total ablation (ice is gone but snow may be left)
(  647)           !-------------------------------------------------------------------
(  648)           zhsn  = h_snow_1d(ji)
(  649)           zihgnew = 1.0 - MAX( zzero , SIGN( zone , -zhicnew ) )
(  650)           zihsn   = MAX( zzero , SIGN( zone , -zhsn ) )
(  651)           !---convert 






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 27

(  652)           zdhicm  = ( 1.0 - zihgnew ) * ( zhicnew - qstbif_1d(ji) / xlic )
(  653)           zdhsnm  = ( 1.0 - zihsn ) * zdhicm * rhoic / rhosn
(  654)           !---updating new ice thickness and computing the newly formed ice mass
(  655)           zhicnew   =  zihgnew * zhicnew
(  656)           rdmicif_1d(ji) =  rdmicif_1d(ji) + ( 1.0 - frld_1d(ji) ) * ( zhicnew - h_ice_1d(ji) ) * rhoic
(  657)           !---updating new snow thickness and computing the newly formed snow mass
(  658)           zhsnfi   = zhsn + zdhsnm
(  659)           h_snow_1d(ji) = MAX( zzero , zhsnfi )
(  660)           rdmsnif_1d(ji) =  rdmsnif_1d(ji) + ( 1.0 - frld_1d(ji) ) * ( h_snow_1d(ji) - zhsn ) * rhosn
(  661)           !--remaining energy in case of total ablation
(  662)           zqocea(ji) = - ( zihsn * xlic * zdhicm + xlsn * ( zhsnfi - h_snow_1d(ji) ) ) * ( 1.0 - frld_1d(ji) )
(  663)           qstbif_1d(ji) = zihgnew * qstbif_1d(ji)
(  664) 
(  665)           !
(  666)           !        9.5. Update internal temperature and ice thickness.
(  667)           !-------------------------------------------------------
(  668)           !
(  669)           sist_1d(ji) = zihgnew * sist_1d(ji) + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  670)           zidhb  = MAX( zzero , SIGN( zone , - zdhicbot(ji) ) )
(  671)           zc1    = - zhicnew * 0.5
(  672)           zpc1   = MIN( 0.5 * zone , - h_ice_1d(ji) * 0.5 - zdhictop(ji) )
(  673)           zc2    = - zhicnew
(  674)           zpc2   =  zidhb * zc2 + ( 1.0 - zidhb ) * ( - h_ice_1d(ji) - zdhictop(ji) )
(  675)           zp1    =  MAX( zpc1 , zc1 )
(  676)           zp2    =  MAX( zpc2 , zc1 )
(  677)           zep(ji) =  tbif_1d(ji,2)
(  678)           ztb2  = 2.0 * (         - zp1   * tbif_1d(ji,2)  &
(  679)              &  + ( zp1 - zp2 ) * tbif_1d(ji,3)  &
(  680)              &  + ( zp2 - zc1 ) * tfu_1d(ji) ) / MAX( zhicnew , epsi20 ) 
(  681)           tbif_1d(ji,2) = zihgnew * ztb2 + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  682)           !---
(  683)           zp1  =  MIN( zpc1 , zc1 )
(  684)           zp2  =  MIN( zpc2 , zc1 )
(  685)           zp1  =  MAX( zc2  , zp1 )
(  686)           ztb3 =  2.0 * (   ( 1.0 - zidhb ) * (  ( zc1 - zp2 ) * tbif_1d(ji,3)  &
(  687)              &                                 + ( zp2 - zc2 ) * tfu_1d(ji) )   &
(  688)              &               +      zidhb   * (  ( zc1 - zp1 ) * zep(ji)      &
(  689)              &                                 + ( zp1 - zc2 ) * tbif_1d(ji,3))  ) / MAX( zhicnew , epsi20 )
(  690)           tbif_1d(ji,3) =  zihgnew * ztb3 + ( 1.0 - zihgnew ) * tfu_1d(ji)
(  691)           h_ice_1d(ji)  =  zhicnew
(  692)        END DO
(  693) 
(  694) 
(  695)        !----------------------------------------------------------------------------
(  696)        !  10. Surface accretion. 
(  697)        !      The change of ice thickness after snow/ice formation is such that
(  698)        !      the interface between snow and ice is located at the same height
(  699)        !      as the ocean surface. It is given by (Fichefet and Morales Maqueda 1999)
(  700)        !          D(h_ice) = (- D(hsn)/alph) =  [rhosn*hsn - (rau0 - rhoic)*hic]
(  701)        !                                     / [alph*rhosn+rau0 - rhoic]
(  702)        !----------------------------------------------------------------------------
(  703)        !
(  704)        DO ji = kideb , kiut
(  705) 
(  706)           !--  Computation of the change of ice thickness after snow-ice formation
(  707)           zdrmh =  ( rhosn * h_snow_1d(ji) + ( rhoic - rau0 ) * h_ice_1d(ji) )  &
(  708)              &  / ( alphs * rhosn + rau0 - rhoic )
(  709)           zdrmh = MAX( zzero , zdrmh )






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 28

(  710) 
(  711)           !--New ice and snow thicknesses Fichefet and Morales Maqueda (1999)
(  712)           zhicnew  = MAX( h_ice_1d(ji) , h_ice_1d(ji) + zdrmh )
(  713)           zhsnnew  = MIN( h_snow_1d(ji) , h_snow_1d(ji) - alphs * zdrmh )
(  714)           !---Compute new ice temperatures. snow temperature remains unchanged
(  715)           !   Lepparanta (1983):
(  716)           zihic = 1.0 - MAX( zzero , SIGN( zone , -zhicnew ) )
(  717)           zquot  = ( 1.0 - zihic ) &
(  718)              &   +         zihic * MIN( zone , h_ice_1d(ji) / MAX( zhicnew , epsi20 ) ) 
(  719)           ztneq  =         alphs * cnscg * tbif_1d(ji,1)    &
(  720)              &   + ( 1.0 - alphs * ( rhosn/rhoic ) ) * tfu_1d(ji)
(  721)           zep(ji) = tbif_1d(ji,2)
(  722)           tbif_1d(ji,2) = ztneq - zquot * zquot * ( ztneq - tbif_1d(ji,2) )
(  723)           tbif_1d(ji,3) = 2.0 * ztneq &
(  724)              &        + zquot * ( tbif_1d(ji,3) + zep(ji) - 2.0 * ztneq ) - tbif_1d(ji,2)
(  725) 
(  726)           !---  Lepparanta (1983) (latent heat released during white ice formation
(  727)           !     goes to the ocean -for lateral ablation-)
(  728)           qldif_1d(ji)  = qldif_1d(ji) + zdrmh * ( 1.0 - alphs * ( rhosn/rhoic ) ) * xlic * ( 1.0 - frld_1d(ji) )
(  729)           !--   Changes in ice volume and ice mass Lepparanta (1983):
(  730)           dvnbq_1d(ji) = ( 1.0 - frld_1d(ji) ) * ( zhicnew - h_ice_1d(ji) )
(  731)           dmgwi_1d(ji) = dmgwi_1d(ji) + ( 1.0 -frld_1d(ji) ) * ( h_snow_1d(ji) - zhsnnew ) * rhosn
(  732)           !  case of natural freshwater flux
(  733) # 734
(  734)           rdmicif_1d(ji) = rdmicif_1d(ji) + ( 1.0 - frld_1d(ji) )   * ( zhicnew - h_ice_1d(ji) )  * rhoic
(  735)           rdmsnif_1d(ji) = rdmsnif_1d(ji) + ( 1.0 - frld_1d(ji) )   * ( zhsnnew - h_snow_1d(ji) ) * rhosn
(  736) 
(  737) 
(  738) # 742
(  742)           !---  Actualize new snow and ice thickness.
(  743) 
(  744)           h_snow_1d(ji)  = zhsnnew
(  745)           h_ice_1d(ji)  = zhicnew
(  746) 
(  747)        END DO
(  748) 
(  749)        !----------------------------------------------------
(  750)        !  11. Lateral ablation (Changes in sea/ice fraction) 
(  751)        !----------------------------------------------------
(  752)        DO ji = kideb , kiut
(  753)           zfrl_old(ji)  = frld_1d(ji)
(  754)           zihic   = 1.0 - MAX( zzero , SIGN( zone , -h_ice_1d(ji) ) )
(  755)           zihsn   = 1.0 - MAX( zzero , SIGN( zone , -h_snow_1d(ji) ) )
(  756)           !--In the case of total ablation (all the ice ice has melted) frld = 1
(  757)           frld_1d(ji)  = ( 1.0 - zihic ) + zihic * zfrl_old(ji)
(  758)           !--Part of solar radiation absorbing inside the ice and going
(  759)           !--through the ocean
(  760)           fscbq_1d(ji) = ( 1.0 - zfrl_old(ji) ) * ( 1.0 - thcm_1d(ji) ) * fstbif_1d(ji)
(  761)           !--Total remaining energy after bottom melting/growing 
(  762)           qfvbq_1d(ji) = zqsup(ji) + ( 1.0 - zihic ) * zqocea(ji)
(  763)           !--Updating of total heat from the ocean
(  764)           qldif_1d(ji)  = qldif_1d(ji) + qfvbq_1d(ji) + ( 1.0 - zihic ) * fscbq_1d(ji) * rdt_ice
(  765)           !--Computation of total heat inside the snow/ice system
(  766)           zqice  = h_snow_1d(ji) * xlsn + h_ice_1d(ji) * xlic
(  767)           zqicetot  = ( 1.0 - frld_1d(ji) ) * zqice
(  768)           !--The concentration of ice is reduced (frld increases) if the heat 
(  769)           !--exchange between ice and ocean is positive
(  770)           ziqf = MAX( zzero , SIGN( zone ,  zqicetot - qldif_1d(ji) ) )






PGF90 (Version     10.5)          09/23/2010  09:49:57      page 29

(  771)           zdfrl = qldif_1d(ji) / MAX( epsi20 , zqice ) 
(  772)           frld_1d(ji)  = ( 1.0 - ziqf )    &
(  773)              &       +         ziqf * ( frld_1d(ji) + MAX( zzero , zdfrl ) ) 
(  774)           fltbif_1d(ji) = ( ( 1.0 - zfrl_old(ji) ) * qstbif_1d(ji) - zqicetot  ) / rdt_ice
(  775)           !--  Opening of leads: Hakkinen & Mellor, 1992.
(  776)           zdfrl = - ( zdhictop(ji) + zdhicbot(ji) ) * hakspl * ( 1.0 - zfrl_old(ji) ) &
(  777)              &  / MAX( epsi13 , h_ice_1d(ji) + h_snow_1d(ji) * rhosn/rhoic ) 
(  778)           zfrld_1d(ji) =  frld_1d(ji) + MAX( zzero , zdfrl )
(  779)           !--Limitation of sea-ice fraction <= 1
(  780)           zfrld_1d(ji) = ziqf * MIN( 0.99 * zone , zfrld_1d(ji) ) + ( 1 - ziqf )
(  781)           !---Update surface and internal temperature and snow/ice thicknesses.
(  782)           sist_1d(ji)   = sist_1d(ji)   + ( 1.0 - ziqf ) * ( tfu_1d(ji) - sist_1d(ji)   )
(  783)           tbif_1d(ji,1) = tbif_1d(ji,1) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,1) )
(  784)           tbif_1d(ji,2) = tbif_1d(ji,2) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,2) )
(  785)           tbif_1d(ji,3) = tbif_1d(ji,3) + ( 1.0 - ziqf ) * ( tfu_1d(ji) - tbif_1d(ji,3) )
(  786)           !--variation of ice volume and ice mass 
(  787)           dvlbq_1d(ji)   = zihic * ( zfrl_old(ji) - frld_1d(ji) ) * h_ice_1d(ji)
(  788)           rdmicif_1d(ji) = rdmicif_1d(ji) + dvlbq_1d(ji) * rhoic
(  789)           !--variation of snow volume and snow mass 
(  790)           zdvsnvol    = zihsn * ( zfrl_old(ji) - frld_1d(ji) ) * h_snow_1d(ji)
(  791)           rdmsnif_1d(ji) = rdmsnif_1d(ji) + zdvsnvol * rhosn
(  792)           h_snow_1d(ji)  = ziqf * h_snow_1d(ji)
(  793) 
(  794)           zdrfrl1 = ziqf * ( 1.0 -  frld_1d(ji) ) / MAX( epsi20 , 1.0 - zfrld_1d(ji) )
(  795)           zdrfrl2 = ziqf * ( 1.0 - zfrl_old(ji) ) / MAX( epsi20 , 1.0 - zfrld_1d(ji) )
(  796) 
(  797)           h_snow_1d (ji) = zdrfrl1 * h_snow_1d(ji)
(  798)           h_ice_1d  (ji) = zdrfrl1 * h_ice_1d(ji)
(  799)           qstbif_1d(ji) = zdrfrl2 * qstbif_1d(ji)
(  800)           frld_1d(ji)    = zfrld_1d(ji)
(  801) 
(  802)        END DO
(  803)        
(  804)     END SUBROUTINE lim_thd_zdf
(  805) # 814
(  814)  END MODULE limthd_zdf
