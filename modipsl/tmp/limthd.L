


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

Filename: limthd.F90

(    1) MODULE limthd
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE limthd   ***
(    4)    !!              LIM thermo ice model : ice thermodynamic
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_thd      : thermodynamic of sea ice
(   11)    !!   lim_thd_init : initialisation of sea-ice thermodynamic
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst          ! physical constants
(   15)    USE dom_oce         ! ocean space and time domain variables
(   16)    USE lbclnk
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ice             ! LIM sea-ice variables
(   19)    USE ice_oce         ! sea-ice/ocean variables
(   20)    USE flx_oce         ! sea-ice/ocean forcings variables 
(   21)    USE thd_ice         ! LIM thermodynamic sea-ice variables
(   22)    USE dom_ice         ! LIM sea-ice domain
(   23)    USE iceini
(   24)    USE limthd_zdf
(   25)    USE limthd_lac
(   26)    USE limtab
(   27)    USE prtctl          ! Print control
(   28)       
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Routine accessibility
(   33)    PUBLIC lim_thd       ! called by lim_step
(   34) 
(   35)    !! * Module variables
(   36)    REAL(wp)  ::            &  ! constant values
(   37)       epsi20 = 1.e-20   ,  &
(   38)       epsi16 = 1.e-16   ,  &
(   39)       epsi04 = 1.e-04   ,  &
(   40)       zzero  = 0.e0     ,  &
(   41)       zone   = 1.e0
(   42) 
(   43)    !! * Substitutions
(   44) #  include "domzgr_substitute.h90"
(   45) #  include "vectopt_loop_substitute.h90"
(   46)    !!-------- -------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 2

(   47)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   48)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd.F90,v 1.8 2006/03/21 08:44:32 opalod Exp $ 
(   49)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)    SUBROUTINE lim_thd
(   55)       !!-------------------------------------------------------------------
(   56)       !!                ***  ROUTINE lim_thd  ***       
(   57)       !!  
(   58)       !! ** Purpose : This routine manages the ice thermodynamic.
(   59)       !!         
(   60)       !! ** Action : - Initialisation of some variables
(   61)       !!             - Some preliminary computation (oceanic heat flux
(   62)       !!               at the ice base, snow acc.,heat budget of the leads)
(   63)       !!             - selection of the icy points and put them in an array
(   64)       !!             - call lim_vert_ther for vert ice thermodynamic
(   65)       !!             - back to the geographic grid
(   66)       !!             - selection of points for lateral accretion
(   67)       !!             - call lim_lat_acc  for the ice accretion
(   68)       !!             - back to the geographic grid
(   69)       !!
(   70)       !! ** References :
(   71)       !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
(   72)       !!
(   73)       !! History :
(   74)       !!   1.0  !  00-01 (LIM)
(   75)       !!   2.0  !  02-07 (C. Ethe, G. Madec) F90
(   76)       !!---------------------------------------------------------------------
(   77)       !! * Local variables
(   78)       INTEGER  ::   ji, jj,    &   ! dummy loop indices
(   79)          nbpb  ,               &   ! nb of icy pts for thermo. cal.
(   80)          nbpac                     ! nb of pts for lateral accretion 
(   81)       CHARACTER (len=22) :: charout
(   82)       REAL(wp) ::  &
(   83)          zfric_umin = 5e-03 ,  &   ! lower bound for the friction velocity
(   84)          zfric_umax = 2e-02        ! upper bound for the friction velocity
(   85)       REAL(wp) ::   &
(   86)          zinda              ,  &   ! switch for test. the val. of concen.
(   87)          zindb, zindg       ,  &   ! switches for test. the val of arg
(   88)          za , zh, zthsnice  ,  &
(   89)          zfric_u            ,  &   ! friction velocity 
(   90)          zfnsol             ,  &   ! total non solar heat
(   91)          zfontn             ,  &   ! heat flux from snow thickness
(   92)          zfntlat, zpareff          ! test. the val. of lead heat budget
(   93)       REAL(wp), DIMENSION(jpi,jpj) :: &
(   94)          zhicifp            ,  &   ! ice thickness for outputs
(   95)          zqlbsbq                   ! link with lead energy budget qldif
(   96)       REAL(wp), DIMENSION(jpi,jpj,2) :: &
(   97)          zmsk                      ! working array
(   98)       !!-------------------------------------------------------------------
(   99) 
(  100)       IF( numit == nstart  )   CALL lim_thd_init  ! Initialization (first time-step only)
(  101)    
(  102)       !-------------------------------------------!
(  103)       !   Initilization of diagnostic variables   !
(  104)       !-------------------------------------------!






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 3

(  105)       
(  106) !i est-ce utile?  oui au moins en partie
(  107)       rdvosif(:,:) = 0.e0   ! variation of ice volume at surface
(  108)       rdvobif(:,:) = 0.e0   ! variation of ice volume at bottom
(  109)       fdvolif(:,:) = 0.e0   ! total variation of ice volume
(  110)       rdvonif(:,:) = 0.e0   ! lateral variation of ice volume
(  111)       fstric (:,:) = 0.e0   ! part of solar radiation absorbing inside the ice
(  112)       fscmbq (:,:) = 0.e0   ! linked with fstric
(  113)       ffltbif(:,:) = 0.e0   ! linked with fstric
(  114)       qfvbq  (:,:) = 0.e0   ! linked with fstric
(  115)       rdmsnif(:,:) = 0.e0   ! variation of snow mass per unit area
(  116)       rdmicif(:,:) = 0.e0   ! variation of ice mass per unit area
(  117)       hicifp (:,:) = 0.e0   ! daily thermodynamic ice production. 
(  118)       zmsk (:,:,:) = 0.e0
(  119) 
(  120)       DO jj = 1, jpj
(  121)          DO ji = 1, jpi
(  122)             hsnif(ji,jj)  = hsnif(ji,jj) *  MAX( zzero, SIGN( zone , hsnif(ji,jj) - epsi04 ) )
(  123)          END DO
(  124)       END DO
(  125) 
(  126)       IF(ln_ctl)   CALL prt_ctl(tab2d_1=hsnif     , clinfo1=' lim_thd: hsnif   : ')
(  127)       
(  128)       !-----------------------------------!
(  129)       !   Treatment of particular cases   !
(  130)       !-----------------------------------!
(  131)       
(  132)       DO jj = 1, jpj
(  133)          DO ji = 1, jpi
(  134)             !  snow is transformed into ice if the original ice cover disappears.
(  135)             zindg         = tms(ji,jj) *  MAX( zzero , SIGN( zone , -hicif(ji,jj) ) )
(  136)             hicif(ji,jj)  = hicif(ji,jj) + zindg * rhosn * hsnif(ji,jj) / rau0
(  137)             hsnif(ji,jj)  = ( zone - zindg ) * hsnif(ji,jj) + zindg * hicif(ji,jj) * ( rau0 - rhoic ) / rhosn
(  138)             dmgwi(ji,jj)  = zindg * (1.0 - frld(ji,jj)) * rhoic * hicif(ji,jj)   ! snow/ice mass
(  139)             
(  140)             !  the lead fraction, frld, must be little than or equal to amax (ice ridging).
(  141)             zthsnice      = hsnif(ji,jj) + hicif(ji,jj)
(  142)             zindb         = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
(  143)             za            = zindb * MIN( zone, ( 1.0 - frld(ji,jj) ) * uscomi )
(  144)             hsnif (ji,jj) = hsnif(ji,jj)  * za
(  145)             hicif (ji,jj) = hicif(ji,jj)  * za
(  146)             qstoif(ji,jj) = qstoif(ji,jj) * za
(  147)             frld  (ji,jj) = 1.0 - zindb * ( 1.0 - frld(ji,jj) ) / MAX( za , epsi20 )
(  148)             
(  149)             !  the in situ ice thickness, hicif, must be equal to or greater than hiclim.
(  150)             zh            = MAX( zone , zindb * hiclim  / MAX( hicif(ji,jj) , epsi20 ) )
(  151)             hsnif (ji,jj) = hsnif(ji,jj)  * zh
(  152)             hicif (ji,jj) = hicif(ji,jj)  * zh
(  153)             qstoif(ji,jj) = qstoif(ji,jj) * zh
(  154)             frld  (ji,jj) = ( frld(ji,jj) + ( zh - 1.0 ) ) / zh
(  155)          END DO
(  156)       END DO
(  157) 
(  158)       IF(ln_ctl) THEN
(  159)          CALL prt_ctl(tab2d_1=hicif  , clinfo1=' lim_thd: hicif   : ')
(  160)          CALL prt_ctl(tab2d_1=hsnif  , clinfo1=' lim_thd: hsnif   : ')
(  161)          CALL prt_ctl(tab2d_1=dmgwi  , clinfo1=' lim_thd: dmgwi   : ')
(  162)          CALL prt_ctl(tab2d_1=qstoif , clinfo1=' lim_thd: qstoif  : ')






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 4

(  163)          CALL prt_ctl(tab2d_1=frld   , clinfo1=' lim_thd: frld    : ')
(  164)       ENDIF
(  165) 
(  166)       
(  167)       !-------------------------------!
(  168)       !   Thermodynamics of sea ice   !
(  169)       !-------------------------------!
(  170)       
(  171)       !      Partial computation of forcing for the thermodynamic sea ice model.
(  172)       !--------------------------------------------------------------------------
(  173) 
(  174)       !CDIR NOVERRCHK
(  175)       DO jj = 1, jpj
(  176)          !CDIR NOVERRCHK
(  177)          DO ji = 1, jpi
(  178)             zthsnice       = hsnif(ji,jj) + hicif(ji,jj)
(  179)             zindb          = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
(  180)             pfrld(ji,jj)   = frld(ji,jj)
(  181)             zinda          = 1.0 - MAX( zzero , SIGN( zone , - ( 1.0 - pfrld(ji,jj) ) ) )
(  182)             
(  183)             !  solar irradiance transmission at the mixed layer bottom and used in the lead heat budget
(  184)             thcm(ji,jj)    = 0.e0 
(  185)             
(  186)             !  net downward heat flux from the ice to the ocean, expressed as a function of ocean 
(  187)             !  temperature and turbulent mixing (McPhee, 1992)
(  188)             zfric_u        = MAX ( MIN( SQRT( ust2s(ji,jj) ) , zfric_umax ) , zfric_umin )  ! friction velocity
(  189)             fdtcn(ji,jj)  = zindb * rau0 * rcp * 0.006  * zfric_u * ( sst_io(ji,jj) - tfu(ji,jj) ) 
(  190)             qdtcn(ji,jj)  = zindb * fdtcn(ji,jj) * frld(ji,jj) * rdt_ice
(  191)                         
(  192)             !  partial computation of the lead energy budget (qldif)
(  193)             zfontn         = ( sprecip(ji,jj) / rhosn ) * xlsn  !   energy for melting
(  194)             zfnsol         = qnsr_oce(ji,jj)  !  total non solar flux
(  195)             qldif(ji,jj)   = tms(ji,jj) * ( qsr_oce(ji,jj) * ( 1.0 - thcm(ji,jj) )   &
(  196)                &                               + zfnsol + fdtcn(ji,jj) - zfontn     &
(  197)                &                               + ( 1.0 - zindb ) * fsbbq(ji,jj) )   &
(  198)                &                               * frld(ji,jj) * rdt_ice    
(  199)             !  parlat : percentage of energy used for lateral ablation (0.0) 
(  200)             zfntlat        = 1.0 - MAX( zzero , SIGN( zone ,  - qldif(ji,jj) ) )
(  201)             zpareff        = 1.0 + ( parlat - 1.0 ) * zinda * zfntlat
(  202)             zqlbsbq(ji,jj) = qldif(ji,jj) * ( 1.0 - zpareff ) / MAX( (1.0 - frld(ji,jj)) * rdt_ice , epsi16 )
(  203)             qldif  (ji,jj) = zpareff *  qldif(ji,jj)
(  204)             qdtcn  (ji,jj) = zpareff * qdtcn(ji,jj)
(  205)             
(  206)             !  energy needed to bring ocean surface layer until its freezing
(  207)             qcmif  (ji,jj) =  rau0 * rcp * fse3t(ji,jj,1) * ( tfu(ji,jj) - sst_io(ji,jj) ) * ( 1 - zinda )
(  208)             
(  209)             !  calculate oceanic heat flux.
(  210)             fbif   (ji,jj) = zindb * (  fsbbq(ji,jj) / MAX( (1.0 - frld(ji,jj)) , epsi20 ) + fdtcn(ji,jj) )
(  211)             
(  212)             ! computation of the daily thermodynamic ice production (only needed for output)
(  213)             zhicifp(ji,jj) = hicif(ji,jj) * ( 1.0 - frld(ji,jj) )
(  214)          END DO
(  215)       END DO
(  216)       
(  217)       
(  218)       !         Select icy points and fulfill arrays for the vectorial grid.
(  219)       !----------------------------------------------------------------------
(  220)       nbpb = 0






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 5

(  221)       DO jj = 1, jpj
(  222)          DO ji = 1, jpi
(  223)             IF ( frld(ji,jj) < 1.0 ) THEN     
(  224)                nbpb      = nbpb + 1
(  225)                npb(nbpb) = (jj - 1) * jpi + ji
(  226)             ENDIF
(  227)          END DO
(  228)       END DO
(  229) 
(  230)       IF(ln_ctl) THEN
(  231)          CALL prt_ctl(tab2d_1=pfrld, clinfo1=' lim_thd: pfrld   : ', tab2d_2=thcm   , clinfo2='  thcm    : ')
(  232)          CALL prt_ctl(tab2d_1=fdtcn, clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn  , clinfo2='  qdtcn   : ')
(  233)          CALL prt_ctl(tab2d_1=qldif, clinfo1=' lim_thd: qldif   : ', tab2d_2=zqlbsbq, clinfo2='  zqlbsbq : ')
(  234)          CALL prt_ctl(tab2d_1=qcmif, clinfo1=' lim_thd: qcmif   : ', tab2d_2=fbif   , clinfo2='  fbif    : ')
(  235)          zmsk(:,:,1) = tms(:,:)
(  236)          CALL prt_ctl(tab2d_1=qcmif  , clinfo1=' lim_thd: qcmif   : ', mask1=zmsk)
(  237)          CALL prt_ctl(tab2d_1=zhicifp, clinfo1=' lim_thd: zhicifp : ')
(  238)          WRITE(charout, FMT="('lim_thd: nbpb = ',I4)") nbpb
(  239)          CALL prt_ctl_info(charout)
(  240)       ENDIF
(  241)       
(  242)       
(  243)       ! If there is no ice, do nothing. Otherwise, compute Top and Bottom accretion/ablation 
(  244)       !------------------------------------------------------------------------------------ 
(  245) 
(  246)       IF ( nbpb > 0) THEN
(  247)          
(  248)          !  put the variable in a 1-D array for thermodynamics process
(  249)          CALL tab_2d_1d( nbpb, frld_1d    (1:nbpb)     , frld       , jpi, jpj, npb(1:nbpb) )
(  250)          CALL tab_2d_1d( nbpb, h_ice_1d   (1:nbpb)     , hicif      , jpi, jpj, npb(1:nbpb) )
(  251)          CALL tab_2d_1d( nbpb, h_snow_1d  (1:nbpb)     , hsnif      , jpi, jpj, npb(1:nbpb) )
(  252)          CALL tab_2d_1d( nbpb, sist_1d    (1:nbpb)     , sist       , jpi, jpj, npb(1:nbpb) )
(  253)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 1 ), tbif(:,:,1), jpi, jpj, npb(1:nbpb) )
(  254)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 2 ), tbif(:,:,2), jpi, jpj, npb(1:nbpb) )
(  255)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 3 ), tbif(:,:,3), jpi, jpj, npb(1:nbpb) )
(  256)          CALL tab_2d_1d( nbpb, qsr_ice_1d (1:nbpb)     , qsr_ice    , jpi, jpj, npb(1:nbpb) )
(  257)          CALL tab_2d_1d( nbpb, fr1_i0_1d  (1:nbpb)     , fr1_i0     , jpi, jpj, npb(1:nbpb) )
(  258)          CALL tab_2d_1d( nbpb, fr2_i0_1d  (1:nbpb)     , fr2_i0     , jpi, jpj, npb(1:nbpb) )
(  259)          CALL tab_2d_1d( nbpb, qnsr_ice_1d(1:nbpb)     , qnsr_ice   , jpi, jpj, npb(1:nbpb) )
(  260) #if ! defined key_coupled
(  261)          CALL tab_2d_1d( nbpb, qla_ice_1d (1:nbpb)     , qla_ice    , jpi, jpj, npb(1:nbpb) )
(  262)          CALL tab_2d_1d( nbpb, dqla_ice_1d(1:nbpb)     , dqla_ice   , jpi, jpj, npb(1:nbpb) )
(  263) #endif
(  264)          CALL tab_2d_1d( nbpb, dqns_ice_1d(1:nbpb)     , dqns_ice   , jpi, jpj, npb(1:nbpb) )
(  265)          CALL tab_2d_1d( nbpb, tfu_1d     (1:nbpb)     , tfu        , jpi, jpj, npb(1:nbpb) )
(  266)          CALL tab_2d_1d( nbpb, sprecip_1d (1:nbpb)     , sprecip    , jpi, jpj, npb(1:nbpb) ) 
(  267)          CALL tab_2d_1d( nbpb, fbif_1d    (1:nbpb)     , fbif       , jpi, jpj, npb(1:nbpb) )
(  268)          CALL tab_2d_1d( nbpb, thcm_1d    (1:nbpb)     , thcm       , jpi, jpj, npb(1:nbpb) )
(  269)          CALL tab_2d_1d( nbpb, qldif_1d   (1:nbpb)     , qldif      , jpi, jpj, npb(1:nbpb) )
(  270)          CALL tab_2d_1d( nbpb, qstbif_1d  (1:nbpb)     , qstoif     , jpi, jpj, npb(1:nbpb) )
(  271)          CALL tab_2d_1d( nbpb, rdmicif_1d (1:nbpb)     , rdmicif    , jpi, jpj, npb(1:nbpb) )
(  272)          CALL tab_2d_1d( nbpb, dmgwi_1d   (1:nbpb)     , dmgwi      , jpi, jpj, npb(1:nbpb) )
(  273)          CALL tab_2d_1d( nbpb, qlbbq_1d   (1:nbpb)     , zqlbsbq    , jpi, jpj, npb(1:nbpb) )
(  274)  
(  275)          CALL lim_thd_zdf( 1, nbpb )       !  compute ice growth
(  276)          
(  277)          !  back to the geographic grid.
(  278)          CALL tab_1d_2d( nbpb, frld       , npb, frld_1d   (1:nbpb)     , jpi, jpj )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 6

(  279)          CALL tab_1d_2d( nbpb, hicif      , npb, h_ice_1d  (1:nbpb)     , jpi, jpj )
(  280)          CALL tab_1d_2d( nbpb, hsnif      , npb, h_snow_1d (1:nbpb)     , jpi, jpj )
(  281)          CALL tab_1d_2d( nbpb, sist       , npb, sist_1d   (1:nbpb)     , jpi, jpj )
(  282)          CALL tab_1d_2d( nbpb, tbif(:,:,1), npb, tbif_1d   (1:nbpb , 1 ), jpi, jpj )   
(  283)          CALL tab_1d_2d( nbpb, tbif(:,:,2), npb, tbif_1d   (1:nbpb , 2 ), jpi, jpj )   
(  284)          CALL tab_1d_2d( nbpb, tbif(:,:,3), npb, tbif_1d   (1:nbpb , 3 ), jpi, jpj )   
(  285)          CALL tab_1d_2d( nbpb, fscmbq     , npb, fscbq_1d  (1:nbpb)     , jpi, jpj )
(  286)          CALL tab_1d_2d( nbpb, ffltbif    , npb, fltbif_1d (1:nbpb)     , jpi, jpj )
(  287)          CALL tab_1d_2d( nbpb, fstric     , npb, fstbif_1d (1:nbpb)     , jpi, jpj )
(  288)          CALL tab_1d_2d( nbpb, qldif      , npb, qldif_1d  (1:nbpb)     , jpi, jpj )
(  289)          CALL tab_1d_2d( nbpb, qfvbq      , npb, qfvbq_1d  (1:nbpb)     , jpi, jpj )
(  290)          CALL tab_1d_2d( nbpb, qstoif     , npb, qstbif_1d (1:nbpb)     , jpi, jpj )
(  291)          CALL tab_1d_2d( nbpb, rdmicif    , npb, rdmicif_1d(1:nbpb)     , jpi, jpj )
(  292)          CALL tab_1d_2d( nbpb, dmgwi      , npb, dmgwi_1d  (1:nbpb)     , jpi, jpj )
(  293)          CALL tab_1d_2d( nbpb, rdmsnif    , npb, rdmsnif_1d(1:nbpb)     , jpi, jpj )
(  294)          CALL tab_1d_2d( nbpb, rdvosif    , npb, dvsbq_1d  (1:nbpb)     , jpi, jpj )
(  295)          CALL tab_1d_2d( nbpb, rdvobif    , npb, dvbbq_1d  (1:nbpb)     , jpi, jpj )
(  296)          CALL tab_1d_2d( nbpb, fdvolif    , npb, dvlbq_1d  (1:nbpb)     , jpi, jpj )
(  297)          CALL tab_1d_2d( nbpb, rdvonif    , npb, dvnbq_1d  (1:nbpb)     , jpi, jpj ) 
(  298) 
(  299)  
(  300)       ENDIF
(  301) 
(  302)       
(  303)       !      Up-date sea ice thickness.
(  304)       !---------------------------------
(  305)       DO jj = 1, jpj
(  306)          DO ji = 1, jpi
(  307)             phicif(ji,jj) = hicif(ji,jj)  
(  308)             hicif(ji,jj)  = hicif(ji,jj) *  ( zone -  MAX( zzero, SIGN( zone, - ( 1.0 - frld(ji,jj) ) ) ) )
(  309)          END DO
(  310)       END DO
(  311) 
(  312)       
(  313)       !      Tricky trick : add 2 to frld in the Southern Hemisphere.
(  314)       !----------------------------------------------------------
(  315)       IF( fcor(1,1) < 0.e0 ) THEN
(  316)          DO jj = 1, njeqm1
(  317)             DO ji = 1, jpi
(  318)                frld(ji,jj) = frld(ji,jj) + 2.0
(  319)             END DO
(  320)          END DO
(  321)       ENDIF
(  322)       
(  323)       
(  324)       !     Select points for lateral accretion (this occurs when heat exchange
(  325)       !     between ice and ocean is negative; ocean losing heat) 
(  326)       !-----------------------------------------------------------------
(  327)       nbpac = 0
(  328)       DO jj = 1, jpj
(  329)          DO ji = 1, jpi
(  330) !i yes!     IF ( ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
(  331)             IF ( tms(ji,jj) * ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
(  332)                nbpac = nbpac + 1
(  333)                npac( nbpac ) = (jj - 1) * jpi + ji
(  334)             ENDIF
(  335)          END DO
(  336)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 7

(  337)       
(  338)       IF(ln_ctl) THEN
(  339)          CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=hicif, clinfo2=' hicif : ')
(  340)          WRITE(charout, FMT="('lim_thd: nbpac = ',I4)") nbpac
(  341)          CALL prt_ctl_info(charout)
(  342)       ENDIF
(  343) 
(  344)       
(  345)       !
(  346)       !     If ocean gains heat do nothing ; otherwise, one performs lateral accretion
(  347)       !--------------------------------------------------------------------------------
(  348) 
(  349)       IF( nbpac > 0 ) THEN
(  350)          
(  351)          !...Put the variable in a 1-D array for lateral accretion
(  352)          CALL tab_2d_1d( nbpac, frld_1d   (1:nbpac)     , frld       , jpi, jpj, npac(1:nbpac) )
(  353)          CALL tab_2d_1d( nbpac, h_snow_1d (1:nbpac)     , hsnif      , jpi, jpj, npac(1:nbpac) )
(  354)          CALL tab_2d_1d( nbpac, h_ice_1d  (1:nbpac)     , hicif      , jpi, jpj, npac(1:nbpac) )
(  355)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 1 ), tbif(:,:,1), jpi, jpj, npac(1:nbpac) )   
(  356)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 2 ), tbif(:,:,2), jpi, jpj, npac(1:nbpac) )   
(  357)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 3 ), tbif(:,:,3), jpi, jpj, npac(1:nbpac) )   
(  358)          CALL tab_2d_1d( nbpac, qldif_1d  (1:nbpac)     , qldif      , jpi, jpj, npac(1:nbpac) )
(  359)          CALL tab_2d_1d( nbpac, qcmif_1d  (1:nbpac)     , qcmif      , jpi, jpj, npac(1:nbpac) )
(  360)          CALL tab_2d_1d( nbpac, qstbif_1d (1:nbpac)     , qstoif     , jpi, jpj, npac(1:nbpac) )
(  361)          CALL tab_2d_1d( nbpac, rdmicif_1d(1:nbpac)     , rdmicif    , jpi, jpj, npac(1:nbpac) )
(  362)          CALL tab_2d_1d( nbpac, dvlbq_1d  (1:nbpac)     , fdvolif    , jpi, jpj, npac(1:nbpac) )
(  363)          CALL tab_2d_1d( nbpac, tfu_1d    (1:nbpac)     , tfu        , jpi, jpj, npac(1:nbpac) )
(  364)         
(  365)          !  call lateral accretion routine.
(  366)          CALL lim_thd_lac( 1 , nbpac )
(  367)          
(  368)          !   back to the geographic grid
(  369)          CALL tab_1d_2d( nbpac, frld       , npac(1:nbpac), frld_1d   (1:nbpac)     , jpi, jpj )
(  370)          CALL tab_1d_2d( nbpac, hsnif      , npac(1:nbpac), h_snow_1d (1:nbpac)     , jpi, jpj )
(  371)          CALL tab_1d_2d( nbpac, hicif      , npac(1:nbpac), h_ice_1d  (1:nbpac)     , jpi, jpj )
(  372)          CALL tab_1d_2d( nbpac, tbif(:,:,1), npac(1:nbpac), tbif_1d   (1:nbpac , 1 ), jpi, jpj )
(  373)          CALL tab_1d_2d( nbpac, tbif(:,:,2), npac(1:nbpac), tbif_1d   (1:nbpac , 2 ), jpi, jpj )
(  374)          CALL tab_1d_2d( nbpac, tbif(:,:,3), npac(1:nbpac), tbif_1d   (1:nbpac , 3 ), jpi, jpj )
(  375)          CALL tab_1d_2d( nbpac, qstoif     , npac(1:nbpac), qstbif_1d (1:nbpac)     , jpi, jpj )
(  376)          CALL tab_1d_2d( nbpac, rdmicif    , npac(1:nbpac), rdmicif_1d(1:nbpac)     , jpi, jpj )
(  377)          CALL tab_1d_2d( nbpac, fdvolif    , npac(1:nbpac), dvlbq_1d  (1:nbpac)     , jpi, jpj )
(  378)         
(  379)       ENDIF
(  380)        
(  381)        
(  382)       !      Recover frld values between 0 and 1 in the Southern Hemisphere (tricky trick)
(  383)       !      Update daily thermodynamic ice production.    
(  384)       !------------------------------------------------------------------------------
(  385)        
(  386)       DO jj = 1, jpj
(  387)          DO ji = 1, jpi
(  388)             frld  (ji,jj) = MIN( frld(ji,jj), ABS( frld(ji,jj) - 2.0 ) )
(  389)             hicifp(ji,jj) =  hicif(ji,jj) * ( 1.0 - frld(ji,jj) ) - zhicifp(ji,jj) + hicifp(ji,jj)
(  390)          END DO
(  391)       END DO
(  392) 
(  393)       IF(ln_ctl) THEN
(  394)          CALL prt_ctl_info(' lim_thd  end  ')






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 8

(  395)          CALL prt_ctl(tab2d_1=hicif , clinfo1=' lim_thd: hicif   : ', tab2d_2=hsnif , clinfo2=' hsnif  : ')
(  396)          CALL prt_ctl(tab2d_1=frld  , clinfo1=' lim_thd: frld    : ', tab2d_2=hicifp, clinfo2=' hicifp : ')
(  397)          CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=pfrld , clinfo2=' pfrld  : ')
(  398)          CALL prt_ctl(tab2d_1=sist  , clinfo1=' lim_thd: sist    : ')
(  399)          CALL prt_ctl(tab2d_1=tbif(:,:,1), clinfo1=' lim_thd: tbif 1  : ')
(  400)          CALL prt_ctl(tab2d_1=tbif(:,:,2), clinfo1=' lim_thd: tbif 2  : ')
(  401)          CALL prt_ctl(tab2d_1=tbif(:,:,3), clinfo1=' lim_thd: tbif 3  : ')
(  402)          CALL prt_ctl(tab2d_1=fdtcn , clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn , clinfo2=' qdtcn  : ')
(  403)          CALL prt_ctl(tab2d_1=qstoif, clinfo1=' lim_thd: qstoif  : ', tab2d_2=fsbbq , clinfo2=' fsbbq  : ')
(  404)       ENDIF
(  405) 
(  406)     END SUBROUTINE lim_thd
(  407) 
(  408) 
(  409)     SUBROUTINE lim_thd_init
(  410)       !!-------------------------------------------------------------------
(  411)       !!                   ***  ROUTINE lim_thd_init *** 
(  412)       !!                 
(  413)       !! ** Purpose :   Physical constants and parameters linked to the ice 
(  414)       !!      thermodynamics
(  415)       !!
(  416)       !! ** Method  :   Read the namicethd namelist and check the ice-thermo
(  417)       !!       parameter values called at the first timestep (nit000)
(  418)       !!
(  419)       !! ** input   :   Namelist namicether
(  420)       !!
(  421)       !! history :
(  422)       !!  8.5  ! 03-08 (C. Ethe) original code
(  423)       !!-------------------------------------------------------------------
(  424)       NAMELIST/namicethd/ hmelt , hiccrit, hicmin, hiclim, amax  ,        &
(  425)          &                swiqst, sbeta  , parlat, hakspl, hibspl, exld,  &
(  426)          &                hakdif, hnzst  , thth  , parsub, alphs
(  427)       !!-------------------------------------------------------------------
(  428)       
(  429) 
(  430)       ! Define the initial parameters
(  431)       ! -------------------------
(  432)       REWIND( numnam_ice )
(  433)       READ  ( numnam_ice , namicethd )
(  434)       IF(lwp) THEN
(  435)          WRITE(numout,*)
(  436)          WRITE(numout,*)'lim_thd_init : ice parameters for ice thermodynamic computation '
(  437)          WRITE(numout,*)'~~~~~~~~~~~~'
(  438)          WRITE(numout,*)'       maximum melting at the bottom                           hmelt        = ', hmelt
(  439)          WRITE(numout,*)'       ice thick. for lateral accretion in NH (SH)             hiccrit(1/2) = ', hiccrit
(  440)          WRITE(numout,*)'       ice thick. corr. to max. energy stored in brine pocket  hicmin       = ', hicmin  
(  441)          WRITE(numout,*)'       minimum ice thickness                                   hiclim       = ', hiclim 
(  442)          WRITE(numout,*)'       maximum lead fraction                                   amax         = ', amax 
(  443)          WRITE(numout,*)'       energy stored in brine pocket (=1) or not (=0)	  swiqst       = ', swiqst 
(  444)          WRITE(numout,*)'       numerical carac. of the scheme for diffusion in ice '
(  445)          WRITE(numout,*)'       Cranck-Nicholson (=0.5), implicit (=1), explicit (=0)   sbeta        = ', sbeta
(  446)          WRITE(numout,*)'       percentage of energy used for lateral ablation          parlat       = ', parlat
(  447)          WRITE(numout,*)'       slope of distr. for Hakkinen-Mellor lateral melting     hakspl       = ', hakspl  
(  448)          WRITE(numout,*)'       slope of distribution for Hibler lateral melting        hibspl       = ', hibspl
(  449)          WRITE(numout,*)'       exponent for leads-closure rate                         exld         = ', exld
(  450)          WRITE(numout,*)'       coefficient for diffusions of ice and snow              hakdif       = ', hakdif
(  451)          WRITE(numout,*)'       threshold thick. for comp. of eq. thermal conductivity  zhth         = ', thth 
(  452)          WRITE(numout,*)'       thickness of the surf. layer in temp. computation       hnzst        = ', hnzst






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 9

(  453)          WRITE(numout,*)'       switch for snow sublimation  (=1) or not (=0)           parsub       = ', parsub  
(  454)          WRITE(numout,*)'       coefficient for snow density when snow ice formation    alphs        = ', alphs
(  455)       ENDIF
(  456)             
(  457)       uscomi = 1.0 / ( 1.0 - amax )   ! inverse of minimum lead fraction
(  458)       rcdsn = hakdif * rcdsn 
(  459)       rcdic = hakdif * rcdic
(  460)       
(  461)       IF ( ( hsndif > 100.e0 ) .OR. ( hicdif > 100.e0 ) ) THEN
(  462)          cnscg = 0.e0
(  463)       ELSE
(  464)          cnscg = rcpsn / rcpic   ! ratio  rcpsn/rcpic
(  465)       ENDIF
(  466)  
(  467)    END SUBROUTINE lim_thd_init
(  468) 
(  469) #else
(  470)    !!----------------------------------------------------------------------
(  471)    !!   Default option          Dummy module           NO LIM sea-ice model
(  472)    !!----------------------------------------------------------------------
(  473) CONTAINS
(  474)    SUBROUTINE lim_thd         ! Dummy routine
(  475)    END SUBROUTINE lim_thd
(  476) #endif
(  477) 
(  478)    !!======================================================================
(  479) END MODULE limthd





































PGF90 (Version     10.5)          09/23/2010  09:49:58      page 10

(    1) # 1 "limthd.F90"
(    1) MODULE limthd
(    2)    !!======================================================================
(    3)    !!                  ***  MODULE limthd   ***
(    4)    !!              LIM thermo ice model : ice thermodynamic
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_thd      : thermodynamic of sea ice
(   11)    !!   lim_thd_init : initialisation of sea-ice thermodynamic
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst          ! physical constants
(   15)    USE dom_oce         ! ocean space and time domain variables
(   16)    USE lbclnk
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ice             ! LIM sea-ice variables
(   19)    USE ice_oce         ! sea-ice/ocean variables
(   20)    USE flx_oce         ! sea-ice/ocean forcings variables 
(   21)    USE thd_ice         ! LIM thermodynamic sea-ice variables
(   22)    USE dom_ice         ! LIM sea-ice domain
(   23)    USE iceini
(   24)    USE limthd_zdf
(   25)    USE limthd_lac
(   26)    USE limtab
(   27)    USE prtctl          ! Print control
(   28)       
(   29)    IMPLICIT NONE
(   30)    PRIVATE
(   31) 
(   32)    !! * Routine accessibility
(   33)    PUBLIC lim_thd       ! called by lim_step
(   34) 
(   35)    !! * Module variables
(   36)    REAL(wp)  ::            &  ! constant values
(   37)       epsi20 = 1.e-20   ,  &
(   38)       epsi16 = 1.e-16   ,  &
(   39)       epsi04 = 1.e-04   ,  &
(   40)       zzero  = 0.e0     ,  &
(   41)       zone   = 1.e0
(   42) 
(   43)    !! * Substitutions
(   44) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 11

(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 45 "limthd.F90"
(   45) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 46 "limthd.F90"
(   46) # 46
(   46)    !!-------- -------------------------------------------------------------
(   47)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   48)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limthd.F90,v 1.8 2006/03/21 08:44:32 opalod Exp $ 
(   49)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   50)    !!----------------------------------------------------------------------
(   51) 
(   52) CONTAINS
(   53) 
(   54)    SUBROUTINE lim_thd
(   55)       !!-------------------------------------------------------------------
(   56)       !!                ***  ROUTINE lim_thd  ***       
(   57)       !!  
(   58)       !! ** Purpose : This routine manages the ice thermodynamic.
(   59)       !!         
(   60)       !! ** Action : - Initialisation of some variables
(   61)       !!             - Some preliminary computation (oceanic heat flux
(   62)       !!               at the ice base, snow acc.,heat budget of the leads)
(   63)       !!             - selection of the icy points and put them in an array
(   64)       !!             - call lim_vert_ther for vert ice thermodynamic
(   65)       !!             - back to the geographic grid
(   66)       !!             - selection of points for lateral accretion
(   67)       !!             - call lim_lat_acc  for the ice accretion
(   68)       !!             - back to the geographic grid
(   69)       !!
(   70)       !! ** References :
(   71)       !!       H. Goosse et al. 1996, Bul. Soc. Roy. Sc. Liege, 65, 87-90
(   72)       !!
(   73)       !! History :
(   74)       !!   1.0  !  00-01 (LIM)
(   75)       !!   2.0  !  02-07 (C. Ethe, G. Madec) F90
(   76)       !!---------------------------------------------------------------------
(   77)       !! * Local variables
(   78)       INTEGER  ::   ji, jj,    &   ! dummy loop indices
(   79)          nbpb  ,               &   ! nb of icy pts for thermo. cal.
(   80)          nbpac                     ! nb of pts for lateral accretion 
(   81)       CHARACTER (len=22) :: charout
(   82)       REAL(wp) ::  &
(   83)          zfric_umin = 5e-03 ,  &   ! lower bound for the friction velocity
(   84)          zfric_umax = 2e-02        ! upper bound for the friction velocity
(   85)       REAL(wp) ::   &
(   86)          zinda              ,  &   ! switch for test. the val. of concen.
(   87)          zindb, zindg       ,  &   ! switches for test. the val of arg






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 12

(   88)          za , zh, zthsnice  ,  &
(   89)          zfric_u            ,  &   ! friction velocity 
(   90)          zfnsol             ,  &   ! total non solar heat
(   91)          zfontn             ,  &   ! heat flux from snow thickness
(   92)          zfntlat, zpareff          ! test. the val. of lead heat budget
(   93)       REAL(wp), DIMENSION(jpi,jpj) :: &
(   94)          zhicifp            ,  &   ! ice thickness for outputs
(   95)          zqlbsbq                   ! link with lead energy budget qldif
(   96)       REAL(wp), DIMENSION(jpi,jpj,2) :: &
(   97)          zmsk                      ! working array
(   98)       !!-------------------------------------------------------------------
(   99) 
(  100)       IF( numit == nstart  )   CALL lim_thd_init  ! Initialization (first time-step only)
(  101)    
(  102)       !-------------------------------------------!
(  103)       !   Initilization of diagnostic variables   !
(  104)       !-------------------------------------------!
(  105)       
(  106) !i est-ce utile?  oui au moins en partie
(  107)       rdvosif(:,:) = 0.e0   ! variation of ice volume at surface
(  108)       rdvobif(:,:) = 0.e0   ! variation of ice volume at bottom
(  109)       fdvolif(:,:) = 0.e0   ! total variation of ice volume
(  110)       rdvonif(:,:) = 0.e0   ! lateral variation of ice volume
(  111)       fstric (:,:) = 0.e0   ! part of solar radiation absorbing inside the ice
(  112)       fscmbq (:,:) = 0.e0   ! linked with fstric
(  113)       ffltbif(:,:) = 0.e0   ! linked with fstric
(  114)       qfvbq  (:,:) = 0.e0   ! linked with fstric
(  115)       rdmsnif(:,:) = 0.e0   ! variation of snow mass per unit area
(  116)       rdmicif(:,:) = 0.e0   ! variation of ice mass per unit area
(  117)       hicifp (:,:) = 0.e0   ! daily thermodynamic ice production. 
(  118)       zmsk (:,:,:) = 0.e0
(  119) 
(  120)       DO jj = 1, jpj
(  121)          DO ji = 1, jpi
(  122)             hsnif(ji,jj)  = hsnif(ji,jj) *  MAX( zzero, SIGN( zone , hsnif(ji,jj) - epsi04 ) )
(  123)          END DO
(  124)       END DO
(  125) 
(  126)       IF(ln_ctl)   CALL prt_ctl(tab2d_1=hsnif     , clinfo1=' lim_thd: hsnif   : ')
(  127)       
(  128)       !-----------------------------------!
(  129)       !   Treatment of particular cases   !
(  130)       !-----------------------------------!
(  131)       
(  132)       DO jj = 1, jpj
(  133)          DO ji = 1, jpi
(  134)             !  snow is transformed into ice if the original ice cover disappears.
(  135)             zindg         = tms(ji,jj) *  MAX( zzero , SIGN( zone , -hicif(ji,jj) ) )
(  136)             hicif(ji,jj)  = hicif(ji,jj) + zindg * rhosn * hsnif(ji,jj) / rau0
(  137)             hsnif(ji,jj)  = ( zone - zindg ) * hsnif(ji,jj) + zindg * hicif(ji,jj) * ( rau0 - rhoic ) / rhosn
(  138)             dmgwi(ji,jj)  = zindg * (1.0 - frld(ji,jj)) * rhoic * hicif(ji,jj)   ! snow/ice mass
(  139)             
(  140)             !  the lead fraction, frld, must be little than or equal to amax (ice ridging).
(  141)             zthsnice      = hsnif(ji,jj) + hicif(ji,jj)
(  142)             zindb         = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
(  143)             za            = zindb * MIN( zone, ( 1.0 - frld(ji,jj) ) * uscomi )
(  144)             hsnif (ji,jj) = hsnif(ji,jj)  * za
(  145)             hicif (ji,jj) = hicif(ji,jj)  * za






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 13

(  146)             qstoif(ji,jj) = qstoif(ji,jj) * za
(  147)             frld  (ji,jj) = 1.0 - zindb * ( 1.0 - frld(ji,jj) ) / MAX( za , epsi20 )
(  148)             
(  149)             !  the in situ ice thickness, hicif, must be equal to or greater than hiclim.
(  150)             zh            = MAX( zone , zindb * hiclim  / MAX( hicif(ji,jj) , epsi20 ) )
(  151)             hsnif (ji,jj) = hsnif(ji,jj)  * zh
(  152)             hicif (ji,jj) = hicif(ji,jj)  * zh
(  153)             qstoif(ji,jj) = qstoif(ji,jj) * zh
(  154)             frld  (ji,jj) = ( frld(ji,jj) + ( zh - 1.0 ) ) / zh
(  155)          END DO
(  156)       END DO
(  157) 
(  158)       IF(ln_ctl) THEN
(  159)          CALL prt_ctl(tab2d_1=hicif  , clinfo1=' lim_thd: hicif   : ')
(  160)          CALL prt_ctl(tab2d_1=hsnif  , clinfo1=' lim_thd: hsnif   : ')
(  161)          CALL prt_ctl(tab2d_1=dmgwi  , clinfo1=' lim_thd: dmgwi   : ')
(  162)          CALL prt_ctl(tab2d_1=qstoif , clinfo1=' lim_thd: qstoif  : ')
(  163)          CALL prt_ctl(tab2d_1=frld   , clinfo1=' lim_thd: frld    : ')
(  164)       ENDIF
(  165) 
(  166)       
(  167)       !-------------------------------!
(  168)       !   Thermodynamics of sea ice   !
(  169)       !-------------------------------!
(  170)       
(  171)       !      Partial computation of forcing for the thermodynamic sea ice model.
(  172)       !--------------------------------------------------------------------------
(  173) 
(  174)       !CDIR NOVERRCHK
(  175)       DO jj = 1, jpj
(  176)          !CDIR NOVERRCHK
(  177)          DO ji = 1, jpi
(  178)             zthsnice       = hsnif(ji,jj) + hicif(ji,jj)
(  179)             zindb          = tms(ji,jj) * ( 1.0 - MAX( zzero , SIGN( zone , - zthsnice ) ) ) 
(  180)             pfrld(ji,jj)   = frld(ji,jj)
(  181)             zinda          = 1.0 - MAX( zzero , SIGN( zone , - ( 1.0 - pfrld(ji,jj) ) ) )
(  182)             
(  183)             !  solar irradiance transmission at the mixed layer bottom and used in the lead heat budget
(  184)             thcm(ji,jj)    = 0.e0 
(  185)             
(  186)             !  net downward heat flux from the ice to the ocean, expressed as a function of ocean 
(  187)             !  temperature and turbulent mixing (McPhee, 1992)
(  188)             zfric_u        = MAX ( MIN( SQRT( ust2s(ji,jj) ) , zfric_umax ) , zfric_umin )  ! friction velocity
(  189)             fdtcn(ji,jj)  = zindb * rau0 * rcp * 0.006  * zfric_u * ( sst_io(ji,jj) - tfu(ji,jj) ) 
(  190)             qdtcn(ji,jj)  = zindb * fdtcn(ji,jj) * frld(ji,jj) * rdt_ice
(  191)                         
(  192)             !  partial computation of the lead energy budget (qldif)
(  193)             zfontn         = ( sprecip(ji,jj) / rhosn ) * xlsn  !   energy for melting
(  194)             zfnsol         = qnsr_oce(ji,jj)  !  total non solar flux
(  195)             qldif(ji,jj)   = tms(ji,jj) * ( qsr_oce(ji,jj) * ( 1.0 - thcm(ji,jj) )   &
(  196)                &                               + zfnsol + fdtcn(ji,jj) - zfontn     &
(  197)                &                               + ( 1.0 - zindb ) * fsbbq(ji,jj) )   &
(  198)                &                               * frld(ji,jj) * rdt_ice    
(  199)             !  parlat : percentage of energy used for lateral ablation (0.0) 
(  200)             zfntlat        = 1.0 - MAX( zzero , SIGN( zone ,  - qldif(ji,jj) ) )
(  201)             zpareff        = 1.0 + ( parlat - 1.0 ) * zinda * zfntlat
(  202)             zqlbsbq(ji,jj) = qldif(ji,jj) * ( 1.0 - zpareff ) / MAX( (1.0 - frld(ji,jj)) * rdt_ice , epsi16 )
(  203)             qldif  (ji,jj) = zpareff *  qldif(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 14

(  204)             qdtcn  (ji,jj) = zpareff * qdtcn(ji,jj)
(  205)             
(  206)             !  energy needed to bring ocean surface layer until its freezing
(  207)             qcmif  (ji,jj) =  rau0 * rcp * e3t(1) * ( tfu(ji,jj) - sst_io(ji,jj) ) * ( 1 - zinda )
(  208)             
(  209)             !  calculate oceanic heat flux.
(  210)             fbif   (ji,jj) = zindb * (  fsbbq(ji,jj) / MAX( (1.0 - frld(ji,jj)) , epsi20 ) + fdtcn(ji,jj) )
(  211)             
(  212)             ! computation of the daily thermodynamic ice production (only needed for output)
(  213)             zhicifp(ji,jj) = hicif(ji,jj) * ( 1.0 - frld(ji,jj) )
(  214)          END DO
(  215)       END DO
(  216)       
(  217)       
(  218)       !         Select icy points and fulfill arrays for the vectorial grid.
(  219)       !----------------------------------------------------------------------
(  220)       nbpb = 0
(  221)       DO jj = 1, jpj
(  222)          DO ji = 1, jpi
(  223)             IF ( frld(ji,jj) < 1.0 ) THEN     
(  224)                nbpb      = nbpb + 1
(  225)                npb(nbpb) = (jj - 1) * jpi + ji
(  226)             ENDIF
(  227)          END DO
(  228)       END DO
(  229) 
(  230)       IF(ln_ctl) THEN
(  231)          CALL prt_ctl(tab2d_1=pfrld, clinfo1=' lim_thd: pfrld   : ', tab2d_2=thcm   , clinfo2='  thcm    : ')
(  232)          CALL prt_ctl(tab2d_1=fdtcn, clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn  , clinfo2='  qdtcn   : ')
(  233)          CALL prt_ctl(tab2d_1=qldif, clinfo1=' lim_thd: qldif   : ', tab2d_2=zqlbsbq, clinfo2='  zqlbsbq : ')
(  234)          CALL prt_ctl(tab2d_1=qcmif, clinfo1=' lim_thd: qcmif   : ', tab2d_2=fbif   , clinfo2='  fbif    : ')
(  235)          zmsk(:,:,1) = tms(:,:)
(  236)          CALL prt_ctl(tab2d_1=qcmif  , clinfo1=' lim_thd: qcmif   : ', mask1=zmsk)
(  237)          CALL prt_ctl(tab2d_1=zhicifp, clinfo1=' lim_thd: zhicifp : ')
(  238)          WRITE(charout, FMT="('lim_thd: nbpb = ',I4)") nbpb
(  239)          CALL prt_ctl_info(charout)
(  240)       ENDIF
(  241)       
(  242)       
(  243)       ! If there is no ice, do nothing. Otherwise, compute Top and Bottom accretion/ablation 
(  244)       !------------------------------------------------------------------------------------ 
(  245) 
(  246)       IF ( nbpb > 0) THEN
(  247)          
(  248)          !  put the variable in a 1-D array for thermodynamics process
(  249)          CALL tab_2d_1d( nbpb, frld_1d    (1:nbpb)     , frld       , jpi, jpj, npb(1:nbpb) )
(  250)          CALL tab_2d_1d( nbpb, h_ice_1d   (1:nbpb)     , hicif      , jpi, jpj, npb(1:nbpb) )
(  251)          CALL tab_2d_1d( nbpb, h_snow_1d  (1:nbpb)     , hsnif      , jpi, jpj, npb(1:nbpb) )
(  252)          CALL tab_2d_1d( nbpb, sist_1d    (1:nbpb)     , sist       , jpi, jpj, npb(1:nbpb) )
(  253)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 1 ), tbif(:,:,1), jpi, jpj, npb(1:nbpb) )
(  254)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 2 ), tbif(:,:,2), jpi, jpj, npb(1:nbpb) )
(  255)          CALL tab_2d_1d( nbpb, tbif_1d    (1:nbpb , 3 ), tbif(:,:,3), jpi, jpj, npb(1:nbpb) )
(  256)          CALL tab_2d_1d( nbpb, qsr_ice_1d (1:nbpb)     , qsr_ice    , jpi, jpj, npb(1:nbpb) )
(  257)          CALL tab_2d_1d( nbpb, fr1_i0_1d  (1:nbpb)     , fr1_i0     , jpi, jpj, npb(1:nbpb) )
(  258)          CALL tab_2d_1d( nbpb, fr2_i0_1d  (1:nbpb)     , fr2_i0     , jpi, jpj, npb(1:nbpb) )
(  259)          CALL tab_2d_1d( nbpb, qnsr_ice_1d(1:nbpb)     , qnsr_ice   , jpi, jpj, npb(1:nbpb) )
(  260) # 261
(  261)          CALL tab_2d_1d( nbpb, qla_ice_1d (1:nbpb)     , qla_ice    , jpi, jpj, npb(1:nbpb) )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 15

(  262)          CALL tab_2d_1d( nbpb, dqla_ice_1d(1:nbpb)     , dqla_ice   , jpi, jpj, npb(1:nbpb) )
(  263) # 264
(  264)          CALL tab_2d_1d( nbpb, dqns_ice_1d(1:nbpb)     , dqns_ice   , jpi, jpj, npb(1:nbpb) )
(  265)          CALL tab_2d_1d( nbpb, tfu_1d     (1:nbpb)     , tfu        , jpi, jpj, npb(1:nbpb) )
(  266)          CALL tab_2d_1d( nbpb, sprecip_1d (1:nbpb)     , sprecip    , jpi, jpj, npb(1:nbpb) ) 
(  267)          CALL tab_2d_1d( nbpb, fbif_1d    (1:nbpb)     , fbif       , jpi, jpj, npb(1:nbpb) )
(  268)          CALL tab_2d_1d( nbpb, thcm_1d    (1:nbpb)     , thcm       , jpi, jpj, npb(1:nbpb) )
(  269)          CALL tab_2d_1d( nbpb, qldif_1d   (1:nbpb)     , qldif      , jpi, jpj, npb(1:nbpb) )
(  270)          CALL tab_2d_1d( nbpb, qstbif_1d  (1:nbpb)     , qstoif     , jpi, jpj, npb(1:nbpb) )
(  271)          CALL tab_2d_1d( nbpb, rdmicif_1d (1:nbpb)     , rdmicif    , jpi, jpj, npb(1:nbpb) )
(  272)          CALL tab_2d_1d( nbpb, dmgwi_1d   (1:nbpb)     , dmgwi      , jpi, jpj, npb(1:nbpb) )
(  273)          CALL tab_2d_1d( nbpb, qlbbq_1d   (1:nbpb)     , zqlbsbq    , jpi, jpj, npb(1:nbpb) )
(  274)  
(  275)          CALL lim_thd_zdf( 1, nbpb )       !  compute ice growth
(  276)          
(  277)          !  back to the geographic grid.
(  278)          CALL tab_1d_2d( nbpb, frld       , npb, frld_1d   (1:nbpb)     , jpi, jpj )
(  279)          CALL tab_1d_2d( nbpb, hicif      , npb, h_ice_1d  (1:nbpb)     , jpi, jpj )
(  280)          CALL tab_1d_2d( nbpb, hsnif      , npb, h_snow_1d (1:nbpb)     , jpi, jpj )
(  281)          CALL tab_1d_2d( nbpb, sist       , npb, sist_1d   (1:nbpb)     , jpi, jpj )
(  282)          CALL tab_1d_2d( nbpb, tbif(:,:,1), npb, tbif_1d   (1:nbpb , 1 ), jpi, jpj )   
(  283)          CALL tab_1d_2d( nbpb, tbif(:,:,2), npb, tbif_1d   (1:nbpb , 2 ), jpi, jpj )   
(  284)          CALL tab_1d_2d( nbpb, tbif(:,:,3), npb, tbif_1d   (1:nbpb , 3 ), jpi, jpj )   
(  285)          CALL tab_1d_2d( nbpb, fscmbq     , npb, fscbq_1d  (1:nbpb)     , jpi, jpj )
(  286)          CALL tab_1d_2d( nbpb, ffltbif    , npb, fltbif_1d (1:nbpb)     , jpi, jpj )
(  287)          CALL tab_1d_2d( nbpb, fstric     , npb, fstbif_1d (1:nbpb)     , jpi, jpj )
(  288)          CALL tab_1d_2d( nbpb, qldif      , npb, qldif_1d  (1:nbpb)     , jpi, jpj )
(  289)          CALL tab_1d_2d( nbpb, qfvbq      , npb, qfvbq_1d  (1:nbpb)     , jpi, jpj )
(  290)          CALL tab_1d_2d( nbpb, qstoif     , npb, qstbif_1d (1:nbpb)     , jpi, jpj )
(  291)          CALL tab_1d_2d( nbpb, rdmicif    , npb, rdmicif_1d(1:nbpb)     , jpi, jpj )
(  292)          CALL tab_1d_2d( nbpb, dmgwi      , npb, dmgwi_1d  (1:nbpb)     , jpi, jpj )
(  293)          CALL tab_1d_2d( nbpb, rdmsnif    , npb, rdmsnif_1d(1:nbpb)     , jpi, jpj )
(  294)          CALL tab_1d_2d( nbpb, rdvosif    , npb, dvsbq_1d  (1:nbpb)     , jpi, jpj )
(  295)          CALL tab_1d_2d( nbpb, rdvobif    , npb, dvbbq_1d  (1:nbpb)     , jpi, jpj )
(  296)          CALL tab_1d_2d( nbpb, fdvolif    , npb, dvlbq_1d  (1:nbpb)     , jpi, jpj )
(  297)          CALL tab_1d_2d( nbpb, rdvonif    , npb, dvnbq_1d  (1:nbpb)     , jpi, jpj ) 
(  298) 
(  299)  
(  300)       ENDIF
(  301) 
(  302)       
(  303)       !      Up-date sea ice thickness.
(  304)       !---------------------------------
(  305)       DO jj = 1, jpj
(  306)          DO ji = 1, jpi
(  307)             phicif(ji,jj) = hicif(ji,jj)  
(  308)             hicif(ji,jj)  = hicif(ji,jj) *  ( zone -  MAX( zzero, SIGN( zone, - ( 1.0 - frld(ji,jj) ) ) ) )
(  309)          END DO
(  310)       END DO
(  311) 
(  312)       
(  313)       !      Tricky trick : add 2 to frld in the Southern Hemisphere.
(  314)       !----------------------------------------------------------
(  315)       IF( fcor(1,1) < 0.e0 ) THEN
(  316)          DO jj = 1, njeqm1
(  317)             DO ji = 1, jpi
(  318)                frld(ji,jj) = frld(ji,jj) + 2.0
(  319)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 16

(  320)          END DO
(  321)       ENDIF
(  322)       
(  323)       
(  324)       !     Select points for lateral accretion (this occurs when heat exchange
(  325)       !     between ice and ocean is negative; ocean losing heat) 
(  326)       !-----------------------------------------------------------------
(  327)       nbpac = 0
(  328)       DO jj = 1, jpj
(  329)          DO ji = 1, jpi
(  330) !i yes!     IF ( ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
(  331)             IF ( tms(ji,jj) * ( qcmif(ji,jj) - qldif(ji,jj) ) > 0.e0 ) THEN
(  332)                nbpac = nbpac + 1
(  333)                npac( nbpac ) = (jj - 1) * jpi + ji
(  334)             ENDIF
(  335)          END DO
(  336)       END DO
(  337)       
(  338)       IF(ln_ctl) THEN
(  339)          CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=hicif, clinfo2=' hicif : ')
(  340)          WRITE(charout, FMT="('lim_thd: nbpac = ',I4)") nbpac
(  341)          CALL prt_ctl_info(charout)
(  342)       ENDIF
(  343) 
(  344)       
(  345)       !
(  346)       !     If ocean gains heat do nothing ; otherwise, one performs lateral accretion
(  347)       !--------------------------------------------------------------------------------
(  348) 
(  349)       IF( nbpac > 0 ) THEN
(  350)          
(  351)          !...Put the variable in a 1-D array for lateral accretion
(  352)          CALL tab_2d_1d( nbpac, frld_1d   (1:nbpac)     , frld       , jpi, jpj, npac(1:nbpac) )
(  353)          CALL tab_2d_1d( nbpac, h_snow_1d (1:nbpac)     , hsnif      , jpi, jpj, npac(1:nbpac) )
(  354)          CALL tab_2d_1d( nbpac, h_ice_1d  (1:nbpac)     , hicif      , jpi, jpj, npac(1:nbpac) )
(  355)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 1 ), tbif(:,:,1), jpi, jpj, npac(1:nbpac) )   
(  356)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 2 ), tbif(:,:,2), jpi, jpj, npac(1:nbpac) )   
(  357)          CALL tab_2d_1d( nbpac, tbif_1d   (1:nbpac , 3 ), tbif(:,:,3), jpi, jpj, npac(1:nbpac) )   
(  358)          CALL tab_2d_1d( nbpac, qldif_1d  (1:nbpac)     , qldif      , jpi, jpj, npac(1:nbpac) )
(  359)          CALL tab_2d_1d( nbpac, qcmif_1d  (1:nbpac)     , qcmif      , jpi, jpj, npac(1:nbpac) )
(  360)          CALL tab_2d_1d( nbpac, qstbif_1d (1:nbpac)     , qstoif     , jpi, jpj, npac(1:nbpac) )
(  361)          CALL tab_2d_1d( nbpac, rdmicif_1d(1:nbpac)     , rdmicif    , jpi, jpj, npac(1:nbpac) )
(  362)          CALL tab_2d_1d( nbpac, dvlbq_1d  (1:nbpac)     , fdvolif    , jpi, jpj, npac(1:nbpac) )
(  363)          CALL tab_2d_1d( nbpac, tfu_1d    (1:nbpac)     , tfu        , jpi, jpj, npac(1:nbpac) )
(  364)         
(  365)          !  call lateral accretion routine.
(  366)          CALL lim_thd_lac( 1 , nbpac )
(  367)          
(  368)          !   back to the geographic grid
(  369)          CALL tab_1d_2d( nbpac, frld       , npac(1:nbpac), frld_1d   (1:nbpac)     , jpi, jpj )
(  370)          CALL tab_1d_2d( nbpac, hsnif      , npac(1:nbpac), h_snow_1d (1:nbpac)     , jpi, jpj )
(  371)          CALL tab_1d_2d( nbpac, hicif      , npac(1:nbpac), h_ice_1d  (1:nbpac)     , jpi, jpj )
(  372)          CALL tab_1d_2d( nbpac, tbif(:,:,1), npac(1:nbpac), tbif_1d   (1:nbpac , 1 ), jpi, jpj )
(  373)          CALL tab_1d_2d( nbpac, tbif(:,:,2), npac(1:nbpac), tbif_1d   (1:nbpac , 2 ), jpi, jpj )
(  374)          CALL tab_1d_2d( nbpac, tbif(:,:,3), npac(1:nbpac), tbif_1d   (1:nbpac , 3 ), jpi, jpj )
(  375)          CALL tab_1d_2d( nbpac, qstoif     , npac(1:nbpac), qstbif_1d (1:nbpac)     , jpi, jpj )
(  376)          CALL tab_1d_2d( nbpac, rdmicif    , npac(1:nbpac), rdmicif_1d(1:nbpac)     , jpi, jpj )
(  377)          CALL tab_1d_2d( nbpac, fdvolif    , npac(1:nbpac), dvlbq_1d  (1:nbpac)     , jpi, jpj )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 17

(  378)         
(  379)       ENDIF
(  380)        
(  381)        
(  382)       !      Recover frld values between 0 and 1 in the Southern Hemisphere (tricky trick)
(  383)       !      Update daily thermodynamic ice production.    
(  384)       !------------------------------------------------------------------------------
(  385)        
(  386)       DO jj = 1, jpj
(  387)          DO ji = 1, jpi
(  388)             frld  (ji,jj) = MIN( frld(ji,jj), ABS( frld(ji,jj) - 2.0 ) )
(  389)             hicifp(ji,jj) =  hicif(ji,jj) * ( 1.0 - frld(ji,jj) ) - zhicifp(ji,jj) + hicifp(ji,jj)
(  390)          END DO
(  391)       END DO
(  392) 
(  393)       IF(ln_ctl) THEN
(  394)          CALL prt_ctl_info(' lim_thd  end  ')
(  395)          CALL prt_ctl(tab2d_1=hicif , clinfo1=' lim_thd: hicif   : ', tab2d_2=hsnif , clinfo2=' hsnif  : ')
(  396)          CALL prt_ctl(tab2d_1=frld  , clinfo1=' lim_thd: frld    : ', tab2d_2=hicifp, clinfo2=' hicifp : ')
(  397)          CALL prt_ctl(tab2d_1=phicif, clinfo1=' lim_thd: phicif  : ', tab2d_2=pfrld , clinfo2=' pfrld  : ')
(  398)          CALL prt_ctl(tab2d_1=sist  , clinfo1=' lim_thd: sist    : ')
(  399)          CALL prt_ctl(tab2d_1=tbif(:,:,1), clinfo1=' lim_thd: tbif 1  : ')
(  400)          CALL prt_ctl(tab2d_1=tbif(:,:,2), clinfo1=' lim_thd: tbif 2  : ')
(  401)          CALL prt_ctl(tab2d_1=tbif(:,:,3), clinfo1=' lim_thd: tbif 3  : ')
(  402)          CALL prt_ctl(tab2d_1=fdtcn , clinfo1=' lim_thd: fdtcn   : ', tab2d_2=qdtcn , clinfo2=' qdtcn  : ')
(  403)          CALL prt_ctl(tab2d_1=qstoif, clinfo1=' lim_thd: qstoif  : ', tab2d_2=fsbbq , clinfo2=' fsbbq  : ')
(  404)       ENDIF
(  405) 
(  406)     END SUBROUTINE lim_thd
(  407) 
(  408) 
(  409)     SUBROUTINE lim_thd_init
(  410)       !!-------------------------------------------------------------------
(  411)       !!                   ***  ROUTINE lim_thd_init *** 
(  412)       !!                 
(  413)       !! ** Purpose :   Physical constants and parameters linked to the ice 
(  414)       !!      thermodynamics
(  415)       !!
(  416)       !! ** Method  :   Read the namicethd namelist and check the ice-thermo
(  417)       !!       parameter values called at the first timestep (nit000)
(  418)       !!
(  419)       !! ** input   :   Namelist namicether
(  420)       !!
(  421)       !! history :
(  422)       !!  8.5  ! 03-08 (C. Ethe) original code
(  423)       !!-------------------------------------------------------------------
(  424)       NAMELIST/namicethd/ hmelt , hiccrit, hicmin, hiclim, amax  ,        &
(  425)          &                swiqst, sbeta  , parlat, hakspl, hibspl, exld,  &
(  426)          &                hakdif, hnzst  , thth  , parsub, alphs
(  427)       !!-------------------------------------------------------------------
(  428)       
(  429) 
(  430)       ! Define the initial parameters
(  431)       ! -------------------------
(  432)       REWIND( numnam_ice )
(  433)       READ  ( numnam_ice , namicethd )
(  434)       IF(lwp) THEN
(  435)          WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 18

(  436)          WRITE(numout,*)'lim_thd_init : ice parameters for ice thermodynamic computation '
(  437)          WRITE(numout,*)'~~~~~~~~~~~~'
(  438)          WRITE(numout,*)'       maximum melting at the bottom                           hmelt        = ', hmelt
(  439)          WRITE(numout,*)'       ice thick. for lateral accretion in NH (SH)             hiccrit(1/2) = ', hiccrit
(  440)          WRITE(numout,*)'       ice thick. corr. to max. energy stored in brine pocket  hicmin       = ', hicmin  
(  441)          WRITE(numout,*)'       minimum ice thickness                                   hiclim       = ', hiclim 
(  442)          WRITE(numout,*)'       maximum lead fraction                                   amax         = ', amax 
(  443)          WRITE(numout,*)'       energy stored in brine pocket (=1) or not (=0)	  swiqst       = ', swiqst 
(  444)          WRITE(numout,*)'       numerical carac. of the scheme for diffusion in ice '
(  445)          WRITE(numout,*)'       Cranck-Nicholson (=0.5), implicit (=1), explicit (=0)   sbeta        = ', sbeta
(  446)          WRITE(numout,*)'       percentage of energy used for lateral ablation          parlat       = ', parlat
(  447)          WRITE(numout,*)'       slope of distr. for Hakkinen-Mellor lateral melting     hakspl       = ', hakspl  
(  448)          WRITE(numout,*)'       slope of distribution for Hibler lateral melting        hibspl       = ', hibspl
(  449)          WRITE(numout,*)'       exponent for leads-closure rate                         exld         = ', exld
(  450)          WRITE(numout,*)'       coefficient for diffusions of ice and snow              hakdif       = ', hakdif
(  451)          WRITE(numout,*)'       threshold thick. for comp. of eq. thermal conductivity  zhth         = ', thth 
(  452)          WRITE(numout,*)'       thickness of the surf. layer in temp. computation       hnzst        = ', hnzst
(  453)          WRITE(numout,*)'       switch for snow sublimation  (=1) or not (=0)           parsub       = ', parsub  
(  454)          WRITE(numout,*)'       coefficient for snow density when snow ice formation    alphs        = ', alphs
(  455)       ENDIF
(  456)             
(  457)       uscomi = 1.0 / ( 1.0 - amax )   ! inverse of minimum lead fraction
(  458)       rcdsn = hakdif * rcdsn 
(  459)       rcdic = hakdif * rcdic
(  460)       
(  461)       IF ( ( hsndif > 100.e0 ) .OR. ( hicdif > 100.e0 ) ) THEN
(  462)          cnscg = 0.e0
(  463)       ELSE
(  464)          cnscg = rcpsn / rcpic   ! ratio  rcpsn/rcpic
(  465)       ENDIF
(  466)  
(  467)    END SUBROUTINE lim_thd_init
(  468) 
(  469) 
(  470) # 478
(  478)    !!======================================================================
(  479) END MODULE limthd
