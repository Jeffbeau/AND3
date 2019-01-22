


PGF90 (Version     10.5)          09/23/2010  09:50:22      page 1

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

Filename: step.F90

(    1) MODULE step
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE step  ***
(    4)    !! Time-stepping    : manager of the ocean, tracer and ice time stepping
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   stp            : OPA system time-stepping
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE ldftra_oce
(   15)    USE ldfdyn_oce
(   16)    USE cpl_oce         ! coupled ocean-atmosphere variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk
(   19) 
(   20)    USE daymod          ! calendar                         (day     routine)
(   21) 
(   22)    USE dtatem          ! ocean temperature data           (dta_tem routine)
(   23)    USE dtasal          ! ocean salinity    data           (dta_sal routine)
(   24)    USE dtasst          ! ocean sea surface temperature    (dta_sst routine)
(   25)    USE dtasss          ! ocean sea surface salinity       (dta_sss routine)
(   26)    USE taumod          ! surface stress                   (tau     routine)
(   27)    USE flxmod          ! thermohaline fluxes              (flx     routine)
(   28)    USE ocesbc          ! thermohaline fluxes              (oce_sbc routine)
(   29)    USE flxrnf          ! runoffs                          (flx_rnf routine)
(   30)    USE flxfwb          ! freshwater budget correction     (flx_fwb routine)
(   31)    USE closea          ! closed sea freshwater budget     (flx_clo routine)
(   32)    USE ocfzpt          ! surface ocean freezing point    (oc_fz_pt routine)
(   33) 
(   34)    USE trcstp          ! passive tracer time-stepping      (trc_stp routine)
(   35) 
(   36)    USE dynhpg          ! hydrostatic pressure grad.       (dyn_hpg routine)
(   37)    USE dynhpg_atsk     ! hydrostatic pressure grad.  (dyn_hpg_atsk routine)
(   38)    USE dynspg_oce      ! surface pressure gradient        (dyn_spg routine)
(   39)    USE dynspg          ! surface pressure gradient        (dyn_spg routine)
(   40)    USE dynkeg          ! kinetic energy gradient          (dyn_keg routine)
(   41)    USE dynvor          ! vorticity term              (dyn_vor_... routines)
(   42)    USE dynzad          ! vertical advection               (dyn_adv routine)
(   43)    USE dynldf_bilapg   ! lateral mixing            (dyn_ldf_bilapg routine)
(   44)    USE dynldf_bilap    ! lateral mixing             (dyn_ldf_bilap routine)
(   45)    USE dynldf_iso      ! lateral mixing               (dyn_ldf_iso routine)
(   46)    USE dynldf_lap      ! lateral mixing               (dyn_ldf_lap routine)






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 2

(   47)    USE dynzdf_imp      ! vertical diffusion: implicit     (dyn_zdf routine)
(   48)    USE dynzdf_imp_atsk ! vertical diffusion: implicit     (dyn_zdf routine)
(   49)    USE dynzdf_iso      ! vertical diffusion: isopycnal    (dyn_zdf routine)
(   50)    USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp routine)
(   51)    USE dynnxt          ! time-stepping                    (dyn_nxt routine)
(   52) 
(   53)    USE trabbc          ! bottom boundary condition        (tra_bbc routine)
(   54)    USE trabbl          ! bottom boundary layer            (tra_bbl routine)
(   55)    USE tradmp          ! internal damping                 (tra_dmp routine)
(   56)    USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
(   57)    USE traldf_bilap    ! lateral mixing             (tra_ldf_bilap routine)
(   58)    USE traldf_iso      ! lateral mixing               (tra_ldf_iso routine)
(   59)    USE traldf_iso_zps  ! lateral mixing           (tra_ldf_iso_zps routine)
(   60)    USE traldf_lap      ! lateral mixing               (tra_ldf_lap routine)
(   61)    USE traqsr          ! solar radiation penetration      (tra_qsr routine)
(   62)    USE tranpc          ! non-penetrative convection       (tra_npc routine)
(   63)    USE tranxt          ! time-stepping                    (tra_nxt routine)
(   64)    USE traadv_ctl      ! advection scheme control     (tra_adv_ctl routine)
(   65)    USE traadv_cen2     ! 2nd order centered scheme   (tra_adv_cen2 routine)
(   66)    USE traadv_tvd      ! TVD scheme                (tra_adv_tvd    routine)
(   67)    USE traadv_muscl    ! MUSCL scheme              (tra_adv_muscl  routine)
(   68)    USE traadv_muscl2   ! MUSCL2 scheme             (tra_adv_muscl2 routine)
(   69) !   USE cla             ! cross land advection             (tra_cla routine)
(   70)    USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp routine)
(   71)    USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp routine)
(   72)    USE trazdf_iso      ! vertical diffusion           (tra_zdf_exp routine)
(   73)    USE trazdf_iso_vopt ! vertical diffusion           (tra_zdf_exp routine)
(   74)    USE trasbc          ! surface boundary condition       (tra_sbc routine)
(   75) 
(   76)    USE eosbn2          ! equation of state                (eos_bn2 routine)
(   77) 
(   78)    USE obc_par         ! open boundary condition variables
(   79)    USE obcdta          ! open boundary condition data     (obc_dta routine)
(   80)    USE obcrst          ! open boundary cond. restart      (obc_rst routine)
(   81)    USE obcrad          ! open boundary cond. radiation    (obc_rad routine)
(   82)    USE obcspg          ! open boundary cond  spg          (obc_spg routine)
(   83) 
(   84)    USE divcur          ! hor. divergence and curl      (div & cur routines)
(   85) !   USE cla_div         ! cross land: hor. divergence      (div_cla routine)
(   86)    USE wzvmod          ! vertical velocity                (wzv     routine)
(   87) 
(   88)    USE ldfslp          ! iso-neutral slopes               (ldf_slp routine)
(   89)    USE ldfeiv          ! eddy induced velocity coef.      (ldf_eiv routine)
(   90)    USE ldfdyn          ! 
(   91)    USE ldftra          
(   92) 
(   93) 
(   94)    USE zdfbfr          ! bottom friction                  (zdf_bfr routine)
(   95)    USE zdftke          ! TKE vertical mixing              (zdf_tke routine)
(   96)    USE zdfkpp          ! KPP vertical mixing              (zdf_kpp routine)
(   97)    USE zdfddm          ! double diffusion mixing          (zdf_ddm routine)
(   98)    USE zdfevd          ! enhanced vertical diffusion      (zdf_evd routine)
(   99)    USE zdfric          ! Richardson vertical mixing       (zdf_ric routine)
(  100)    USE zdfmxl          ! Mixed-layer depth                (zdf_mxl routine)
(  101) 
(  102)    USE zpshde          ! partial step: hor. derivative     (zps_hde routine)
(  103)    USE ice_oce         ! sea-ice variable
(  104)    USE icestp          ! sea-ice time-stepping             (ice_stp routine)






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 3

(  105) 
(  106)    USE diawri          ! Standard run outputs             (dia_wri routine)
(  107)    USE trdicp          ! Ocean momentum/tracers trends    (trd_wri routine)
(  108)    USE trdmld          ! mixed-layer trends               (trd_mld routine)
(  109)    USE trdvor          ! vorticity budget                 (trd_vor routine)
(  110) !   USE diagap          ! hor. mean model-data gap         (dia_gap routine)
(  111)    USE diahdy          ! dynamic height                   (dia_hdy routine)
(  112)    USE diaptr          ! poleward transports              (dia_ptr routine)
(  113)    USE diahth          ! thermocline depth                (dia_hth routine)
(  114)    USE diafwb          ! freshwater budget                (dia_fwb routine)
(  115)    USE diaspr          ! suface pressure (rigid-lid)      (dia_spr routine)
(  116)    USE flo_oce         ! floats variables
(  117)    USE floats          ! floats computation               (flo_stp routine)
(  118) 
(  119)    USE stpctl          ! time stepping control            (stp_ctl routine)
(  120)    USE restart         ! ocean restart                    (rst_wri routine)
(  121)    USE cpl             ! exchanges in coupled mode        (cpl_stp routine)
(  122)    USE prtctl          ! Print control                    (prt_ctl routine)
(  123) 
(  124) !!DB 2009.08.06
(  125)    USE sopa_mc
(  126) 
(  127) 
(  128) #if defined key_agrif
(  129)    USE agrif_opa_sponge ! Momemtum and tracers sponges
(  130) #endif
(  131) 
(  132) 
(  133)    IMPLICIT NONE
(  134)    PRIVATE
(  135) 
(  136)    !! * Routine accessibility
(  137)    PUBLIC stp            ! called by opa.F90
(  138) 
(  139)    !! * Substitutions
(  140) #  include "domzgr_substitute.h90"
(  141) #  include "zdfddm_substitute.h90"
(  142)    !!----------------------------------------------------------------------
(  143)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(  144)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/step.F90,v 1.25 2006/04/26 09:26:48 opalod Exp $ 
(  145)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(  146)    !!----------------------------------------------------------------------
(  147) 
(  148) CONTAINS
(  149) 
(  150)    SUBROUTINE stp( &
(  151) #if !defined key_agrif
(  152)    kstp &
(  153) #endif   
(  154)    )      !!----------------------------------------------------------------------
(  155)       !!                     ***  ROUTINE stp  ***
(  156)       !!                      
(  157)       !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
(  158)       !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
(  159)       !!              - Tme stepping  of TRC (passive tracer eqs.)
(  160)       !! 
(  161)       !! ** Method  : -1- Update forcings and data  
(  162)       !!              -2- Update ocean physics 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 4

(  163)       !!              -3- Compute the t and s trends 
(  164)       !!              -4- Update t and s 
(  165)       !!              -5- Compute the momentum trends
(  166)       !!              -6- Update the horizontal velocity
(  167)       !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
(  168)       !!              -8- Outputs and diagnostics
(  169)       !!
(  170)       !! History :
(  171)       !!        !  91-03  ()  Original code
(  172)       !!        !  91-11  (G. Madec)
(  173)       !!        !  92-06  (M. Imbard)  add a first output record
(  174)       !!        !  96-04  (G. Madec)  introduction of dynspg
(  175)       !!        !  96-04  (M.A. Foujols)  introduction of passive tracer
(  176)       !!   8.0  !  97-06  (G. Madec)  new architecture of call
(  177)       !!   8.2  !  97-06  (G. Madec, M. Imbard, G. Roullet)  free surface
(  178)       !!   8.2  !  99-02  (G. Madec, N. Grima)  hpg implicit
(  179)       !!   8.2  !  00-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
(  180)       !!   9.0  !  02-06  (G. Madec)  free form, suppress macro-tasking
(  181)       !!    "   !  04-08  (C. Talandier) New trends organization
(  182)       !!    "   !  05-01  (C. Ethe) Add the KPP closure scheme
(  183)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  184)       !!----------------------------------------------------------------------
(  185)       !! * Arguments
(  186)       INTEGER &
(  187) #if !defined key_agrif   
(  188)       , INTENT( in ) &
(  189) #endif      
(  190)       ::   kstp   ! ocean time-step index
(  191) 
(  192)       !! * local declarations
(  193)       INTEGER ::   indic    ! error indicator if < 0
(  194)       !! ---------------------------------------------------------------------
(  195)       INTEGER :: ji, jj, jk 
(  196) 
(  197) 
(  198) 
(  199) #if defined key_agrif
(  200)       kstp = nit000 + Agrif_Nb_Step()
(  201) #endif   
(  202)       indic = 1                    ! reset to no error condition
(  203)       adatrj = adatrj + rdt/86400._wp
(  204) 
(  205) !!DB 2008.04.16 -- hardwired ramp function currently used in obcdta, tau_forced_*
(  206)      ramp=tanh(kstp*rdt/(2.0*86400.0))
(  207) 
(  208) !!DB: perpetual forcing on ====> perpetual_forcing /= 0 ====> do NOT call day()
(  209)       if(perpetual_forcing == 0) then
(  210)          CALL day( kstp )             ! Calendar
(  211)       endif
(  212) 
(  213) !!DBG: 2009.06.12
(  214) !      call OB_LIMITER(kstp,1)
(  215) 
(  216) 
(  217) 
(  218)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  219)       ! Update data, open boundaries and Forcings
(  220)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 5

(  221) 
(  222)       IF( lk_dtatem  )   CALL dta_tem( kstp )         ! update 3D temperature data
(  223) 
(  224)       IF( lk_dtasal  )   CALL dta_sal( kstp )         ! Salinity data
(  225) 
(  226)       IF( lk_dtasst  )   CALL dta_sst( kstp )         ! Sea Surface Temperature data
(  227) 
(  228)       IF( lk_dtasss  )   CALL dta_sss( kstp )         ! Sea Surface salinity data
(  229) 
(  230)       IF( lk_obc     )   CALL obc_dta( kstp )         ! update dynamic and tracer data at open boundaries
(  231) 
(  232)       IF( lk_obc     )   CALL obc_rad( kstp )         ! compute phase velocities at open boundaries
(  233)       
(  234)       CALL tau( kstp )             ! wind stress
(  235)       
(  236)       CALL flx_rnf( kstp )         ! runoff data
(  237)       
(  238)       CALL flx( kstp )             ! heat and freshwater fluxes
(  239) 
(  240)       IF( lk_ice_lim )   CALL ice_stp( kstp )         ! sea-ice model (Update stress & fluxes)
(  241)       
(  242)       CALL oce_sbc( kstp )         ! ocean surface boundaries
(  243)       
(  244)       IF( ln_fwb     )   CALL flx_fwb( kstp )         ! freshwater budget
(  245)       
(  246)       IF( nclosea == 1 ) CALL flx_clo( kstp )         ! closed sea in the domain (update freshwater fluxes)
(  247)       
(  248)       IF( kstp == nit000 ) THEN 
(  249)          IF( ninist == 1 ) THEN                       ! Output the initial state and forcings
(  250)             CALL dia_wri_state( 'output.init' )
(  251)          ENDIF
(  252)       ENDIF
(  253)       
(  254)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  255)          CALL prt_ctl(tab2d_1=emp    , clinfo1=' emp  -   : ', mask1=tmask, ovlap=1)
(  256)          CALL prt_ctl(tab2d_1=emps   , clinfo1=' emps -   : ', mask1=tmask, ovlap=1)
(  257)          CALL prt_ctl(tab2d_1=qt     , clinfo1=' qt   -   : ', mask1=tmask, ovlap=1)
(  258)          CALL prt_ctl(tab2d_1=qsr    , clinfo1=' qsr  -   : ', mask1=tmask, ovlap=1)
(  259)          CALL prt_ctl(tab2d_1=runoff , clinfo1=' runoff   : ', mask1=tmask, ovlap=1)
(  260)          CALL prt_ctl(tab3d_1=tmask  , clinfo1=' tmask    : ', mask1=tmask, ovlap=1, kdim=jpk)
(  261)          CALL prt_ctl(tab3d_1=tn     , clinfo1=' sst  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  262)          CALL prt_ctl(tab3d_1=sn     , clinfo1=' sss  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  263)          CALL prt_ctl(tab2d_1=taux   , clinfo1=' tau  - x : ', tab2d_2=tauy, clinfo2='      - y : ', ovlap=1)
(  264)       ENDIF
(  265) 
(  266) 
(  267)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  268)       ! Ocean physics update
(  269)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  270)       !-----------------------------------------------------------------------
(  271)       !  VERTICAL PHYSICS
(  272)       !-----------------------------------------------------------------------
(  273)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  274)       !-----------------------------------------------------------------------
(  275)       
(  276)       CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
(  277)       
(  278)       !                                                     ! Vertical eddy viscosity and diffusivity coefficients






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 6

(  279)       IF( lk_zdfric )   CALL zdf_ric( kstp )                       ! Richardson number dependent Kz
(  280)       IF( lk_zdftke )   CALL zdf_tke( kstp )                       ! TKE closure scheme for Kz
(  281)       IF( lk_zdfkpp )   CALL zdf_kpp( kstp )                       ! KPP closure scheme for Kz
(  282)       IF( lk_zdfcst )   avt (:,:,:) = avt0 * tmask(:,:,:)          ! Constant Kz (reset avt to the background value)
(  283) 
(  284)       IF( ln_zdfevd )   CALL zdf_evd( kstp )                 ! enhanced vertical eddy diffusivity
(  285) 
(  286)       IF( lk_zdfddm .AND. .NOT. lk_zdfkpp)   &
(  287)            &              CALL zdf_ddm( kstp )                 ! double diffusive mixing
(  288) 
(  289)       CALL zdf_bfr( kstp )                 ! bottom friction
(  290)       
(  291)       CALL zdf_mxl( kstp )                 ! mixed layer depth
(  292) 
(  293) 
(  294)       !-----------------------------------------------------------------------
(  295)       !  LATERAL PHYSICS
(  296)       !-----------------------------------------------------------------------
(  297)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  298)       !-----------------------------------------------------------------------
(  299) 
(  300)       IF( lk_ldfslp     )   CALL ldf_slp( kstp, rhd, rn2 )       ! before slope of the lateral mixing
(  301) 
(  302) #if defined key_traldf_c2d
(  303)       IF( lk_traldf_eiv )   CALL ldf_eiv( kstp )                 ! eddy induced velocity coefficient
(  304) #endif
(  305) 
(  306) 
(  307) #if defined key_passivetrc
(  308)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  309)       ! Passive Tracer Model
(  310)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  311)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  312)       !-----------------------------------------------------------------------
(  313) 
(  314)       CALL trc_stp( kstp, indic )            ! time-stepping
(  315) 
(  316) #endif
(  317) 
(  318) 
(  319)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  320)       ! Active tracers
(  321)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  322)       ! N.B. ua, va arrays are used as workspace in this section
(  323)       !-----------------------------------------------------------------------
(  324)       
(  325)       ta(:,:,:) = 0.e0               ! set tracer trends to zero
(  326)       sa(:,:,:) = 0.e0
(  327)       
(  328)       CALL tra_sbc( kstp )           ! surface boundary condition
(  329) 
(  330)       IF( ln_traqsr        )   CALL tra_qsr( kstp )           ! penetrative solar radiation qsr
(  331) 
(  332)       IF( lk_trabbc        )   CALL tra_bbc( kstp )           ! bottom heat flux
(  333)       
(  334)       IF( lk_trabbl_dif    )   CALL tra_bbl_dif( kstp )           ! diffusive bottom boundary layer scheme
(  335)       IF( lk_trabbl_adv    )   CALL tra_bbl_adv( kstp )           ! advective (and/or diffusive) bottom boundary layer schem
(  336)       






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 7

(  337)       IF( lk_tradmp        )   CALL tra_dmp( kstp )           ! internal damping trends
(  338) 
(  339)       !                                                       ! horizontal & vertical advection
(  340)       IF( kstp == nit000   )   CALL tra_adv_ctl                    ! chose/control the scheme used
(  341)       IF( ln_traadv_cen2   )   CALL tra_adv_cen2  ( kstp )         ! 2nd order centered scheme
(  342)       IF( ln_traadv_tvd    )   CALL tra_adv_tvd   ( kstp )         ! TVD scheme
(  343)       IF( ln_traadv_muscl  )   CALL tra_adv_muscl ( kstp )         ! MUSCL scheme
(  344)       IF( ln_traadv_muscl2 )   CALL tra_adv_muscl2( kstp )         ! MUSCL2 scheme
(  345) 
(  346) !!DB: orca-related
(  347) !      IF( n_cla == 1       )   CALL tra_cla( kstp )           ! Cross Land Advection (Update Hor. advection)
(  348) 
(  349)       !                                                       ! lateral mixing 
(  350)       IF( l_traldf_lap     )   CALL tra_ldf_lap    ( kstp )           ! iso-level laplacian
(  351)       IF( l_traldf_bilap   )   CALL tra_ldf_bilap  ( kstp )           ! iso-level bilaplacian 
(  352)       IF( l_traldf_bilapg  )   CALL tra_ldf_bilapg ( kstp )           ! s-coord. horizontal bilaplacian
(  353)       IF( l_traldf_iso     )   CALL tra_ldf_iso    ( kstp )           ! iso-neutral/geopot. laplacian 
(  354)       IF( l_traldf_iso_zps )   CALL tra_ldf_iso_zps( kstp )           ! partial step iso-neutral/geopot. laplacian
(  355) 
(  356) #if defined key_agrif
(  357)       IF (.NOT. Agrif_Root())  CALL Agrif_Sponge_tra( kstp )          ! tracers sponge
(  358) #endif
(  359)       !                                                       ! vertical diffusion
(  360)       IF( l_trazdf_exp     )   CALL tra_zdf_exp     ( kstp )          ! explicit time stepping (time splitting scheme)
(  361)       IF( l_trazdf_imp     )   CALL tra_zdf_imp     ( kstp )          ! implicit time stepping (euler backward)
(  362)       IF( l_trazdf_iso     )   CALL tra_zdf_iso     ( kstp )          ! isopycnal
(  363)       IF( l_trazdf_iso_vo  )   CALL tra_zdf_iso_vopt( kstp )          ! vector opt. isopycnal
(  364)       
(  365)       CALL tra_nxt( kstp )           ! tracer fields at next time step
(  366)       
(  367)       IF( ln_zdfnpc        )   CALL tra_npc( kstp )           ! update the new (t,s) fields by non
(  368)       !                                                       ! penetrative convective adjustment
(  369)       
(  370)       IF( ln_dynhpg_imp    ) THEN                             ! semi-implicit hpg 
(  371)          CALL eos( ta, sa, rhd, rhop )   ! Time-filtered in situ density used in dynhpg module
(  372)          IF( lk_zps    )          CALL zps_hde( kstp, ta, sa, rhd,  & ! Partial steps: time filtered hor. gradient 
(  373)               &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
(  374)               &                                        gtv, gsv, grv )  
(  375)       ELSE                                                    ! centered hpg (default case)
(  376)          CALL eos( tb, sb, rhd, rhop )       ! now (swap=before) in situ density for dynhpg module
(  377)          IF( lk_zps    )          CALL zps_hde( kstp, tb, sb, rhd,  & ! Partial steps: now horizontal gradient
(  378)               &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
(  379)               &                                        gtv, gsv, grv )  
(  380)       ENDIF
(  381)       
(  382)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  383)       ! Dynamics
(  384)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  385)       ! N.B. ta, sa arrays are used as workspace in this section 
(  386)       !-----------------------------------------------------------------------
(  387) 
(  388)       ua(:,:,:) = 0.e0               ! set dynamics trends to zero
(  389)       va(:,:,:) = 0.e0
(  390) 
(  391)       CALL dyn_keg( kstp )           ! horizontal gradient of kinetic energy
(  392) 
(  393)       !                                                       ! vorticity term including Coriolis
(  394)       IF( kstp == nit000   )   CALL dyn_vor_ctl                      ! chose/control the scheme used






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 8

(  395)       IF( ln_dynvor_ens    )   CALL dyn_vor_enstrophy( kstp )        ! enstrophy conserving scheme
(  396)       IF( ln_dynvor_ene    )   CALL dyn_vor_energy   ( kstp )        ! energy conserving scheme
(  397)       IF( ln_dynvor_mix    )   CALL dyn_vor_mixed    ( kstp )        ! mixed energy/enstrophy conserving scheme
(  398)       IF( ln_dynvor_een    )   CALL dyn_vor_ene_ens  ( kstp )        ! combined energy/enstrophy conserving scheme
(  399)       
(  400)       !                                                       ! lateral mixing 
(  401)       IF( l_dynldf_lap     )   CALL dyn_ldf_lap    ( kstp )          ! iso-level laplacian
(  402)       IF( l_dynldf_bilap   )   CALL dyn_ldf_bilap  ( kstp )          ! iso-level bilaplacian 
(  403)       IF( l_dynldf_bilapg  )   CALL dyn_ldf_bilapg ( kstp )          ! s-coord. horizontal bilaplacian
(  404)       IF( l_dynldf_iso     )   CALL dyn_ldf_iso    ( kstp )          ! iso-neutral laplacian 
(  405)       
(  406) #if defined key_agrif
(  407)       IF (.NOT. Agrif_Root())  CALL Agrif_Sponge_dyn( kstp )         ! momemtum sponge
(  408) #endif
(  409)       !                                                       ! horizontal gradient of Hydrostatic pressure 
(  410)       IF ( lk_jki ) THEN
(  411)          CALL dyn_hpg_atsk( kstp )             ! autotask case (j-k-i loop)
(  412)       ELSE
(  413)          CALL dyn_hpg     ( kstp )             ! default case  (k-j-i loop)
(  414)       ENDIF
(  415)       
(  416)       CALL dyn_zad    ( kstp )       ! vertical advection       
(  417)       
(  418)       !                                                       ! vertical diffusion
(  419)       IF( l_dynzdf_exp     )   CALL dyn_zdf_exp    ( kstp )          ! explicit time stepping (time splitting scheme)
(  420)       IF( l_dynzdf_imp     )   CALL dyn_zdf_imp    ( kstp )          ! implicit time stepping (euler backward)
(  421)       IF( l_dynzdf_imp_tsk )   CALL dyn_zdf_imp_tsk( kstp )          ! autotask implicit time stepping (euler backward)
(  422)       IF( l_dynzdf_iso     )   CALL dyn_zdf_iso    ( kstp )          ! iso-neutral case
(  423)       
(  424)       IF( lk_dynspg_rl ) THEN 
(  425)          IF( lk_obc    )       CALL obc_spg( kstp )           ! surface pressure gradient at open boundaries
(  426)       ENDIF
(  427)       indic=0
(  428)       !i bug lbc sur emp
(  429)       CALL lbc_lnk( emp, 'T', 1. )
(  430)       !i
(  431)       CALL dyn_spg( kstp, indic )    ! surface pressure gradient
(  432)       
(  433)       CALL dyn_nxt( kstp )           ! velocity at next time step 
(  434)       
(  435)       
(  436)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  437)       ! Computation of diagnostic variables
(  438)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  439)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  440)       !-----------------------------------------------------------------------
(  441)       
(  442)       CALL oc_fz_pt                        ! ocean surface freezing temperature
(  443)       
(  444)       CALL div_cur( kstp )                 ! Horizontal divergence & Relative vorticity
(  445) 
(  446) !!DB: orca-related      
(  447) !      IF( n_cla == 1 ) CALL div_cla( kstp )                 ! Cross Land Advection (Update Hor. divergence)
(  448)       
(  449)       CALL wzv( kstp )                     ! Vertical velocity
(  450)       
(  451) !!DB 2008.04.07 -- the below ultimately calls the appropriate routine to update the
(  452) !! relevant diffusivity coeffs (see ldfdyn.F90). The prob is that it performs numerous






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 9

(  453) !! calcs every timestep (including reading namelist) that are not necessary
(  454) !! =========> modification ... DONE ...
(  455) 
(  456) !      !ZW
(  457) !      CALL ldf_dyn_init
(  458) !      CALL ldf_tra_init
(  459) !      !ZW
(  460) 
(  461) !!DB 04.09 -- modify above so that init routines are not called every dt
(  462) !!Note that the smag keys are separated to the end of this code fragment
(  463) !!Also, there is a potential problem if there is a mistake in the keys
(  464) !!so that (e.g.) no or more-than-one tracer key is defined. I do not check
(  465) !!for this possibility
(  466) !!Momentum diffusivity updates
(  467) #if defined key_dynldf_c3d
(  468)       CALL ldf_dyn_c3d   ! ahm = 3D coef. = F( longitude, latitude, depth )
(  469) #elif defined key_dynldf_c2d
(  470)       CALL ldf_dyn_c2d   ! ahm = 1D coef. = F( longitude, latitude )
(  471) #elif defined key_dynldf_c1d
(  472)       CALL ldf_dyn_c1d   ! ahm = 1D coef. = F( depth )
(  473) #else
(  474)       !do nothing unless smag is on -- see below
(  475) #endif
(  476) !!tracer diffusivity updates
(  477) #if defined key_traldf_c3d
(  478)       CALL ldf_tra_c3d           ! aht = 3D coef. = F( longitude, latitude, depth )
(  479) #elif defined key_traldf_c2d
(  480)       CALL ldf_tra_c2d           ! aht = 2D coef. = F( longitude, latitude )
(  481) #elif defined key_traldf_c1d
(  482)       CALL ldf_tra_c1d           ! aht = 1D coef. = F( depth )
(  483) #else
(  484)       !do nothing unless smag is on -- see below
(  485) #endif
(  486) 
(  487) #if defined key_dynldf_smag || defined key_traldf_smag
(  488)       call ldf_smag( kstp )
(  489) #endif
(  490) 
(  491)       
(  492)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  493)       ! Control, diagnostics and outputs
(  494)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  495)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  496)       !-----------------------------------------------------------------------
(  497)       
(  498)       !                                            ! Time loop: control and print
(  499)       CALL stp_ctl( kstp, indic )
(  500)       IF ( indic < 0 )   nstop = nstop + 1
(  501) 
(  502) 
(  503)       IF ( nstop == 0 ) THEN
(  504)          !                                         ! Diagnostics:
(  505)          IF( lk_floats  )   CALL flo_stp( kstp )                 ! drifting Floats
(  506)          IF( lk_trddyn  )   CALL trd_dwr( kstp )                 ! trends: dynamics 
(  507)          IF( lk_trdtra  )   CALL trd_twr( kstp )                 ! trends: active tracers
(  508)          IF( lk_trdmld  )   CALL trd_mld( kstp )                 ! trends: Mixed-layer 
(  509)          IF( lk_trdvor  )   CALL trd_vor( kstp )                 ! trends: vorticity budget
(  510)          IF( lk_diaspr  )   CALL dia_spr( kstp )                 ! Surface pressure diagnostics






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 10

(  511)          IF( lk_diahth  )   CALL dia_hth( kstp )                 ! Thermocline depth (20 degres isotherm depth)
(  512) !         IF( lk_diagap  )   CALL dia_gap( kstp )                 ! basin averaged diagnostics
(  513)          IF( lk_diahdy  )   CALL dia_hdy( kstp )                 ! dynamical heigh diagnostics
(  514)          IF( lk_diafwb  )   CALL dia_fwb( kstp )                 ! Fresh water budget diagnostics
(  515)          IF( ln_diaptr  )   CALL dia_ptr( kstp )                 ! Poleward TRansports diagnostics
(  516)          
(  517)          !                                         ! save and outputs
(  518)          CALL rst_write  ( kstp )             ! ocean model: restart file output
(  519) !!DB
(  520) !         IF( lk_obc     )   CALL obc_rst_wri( kstp )             ! ocean model: open boundary restart file output
(  521)          
(  522) !!DB 2009.08.06 -- These are the *_grid_*.nc files which we never look at
(  523) !! (Replaced by M2 aves for U,V and *_aveTSUV.nc (see below) for other variables
(  524) !         CALL dia_wri    ( kstp, indic )      ! ocean model: outputs
(  525) 
(  526) !!DB
(  527)          if(M2_ave > 0) call output_special(kstp, indic)      ! special M2 time-averaged fields
(  528) !!DB
(  529)          if(ioutput_ave /= 0) call output_aveTSUV(kstp, indic)      ! special time-averaged fields
(  530) 
(  531)          
(  532)       ENDIF
(  533) 
(  534) !DB Output some BoF vels -- NB: code moved to sopa_mc module
(  535) !! Model diagnostics (somewhat) specific to sopa MC domain 
(  536)       call sopa_mc_diagnostics(kstp)
(  537) 
(  538) 
(  539) 
(  540) 
(  541)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  542)       ! Coupled mode
(  543)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  544) 
(  545)       IF( lk_cpl    )   CALL cpl_stp( kstp )                 ! coupled mode : field exchanges
(  546) 
(  547)    END SUBROUTINE stp
(  548) 
(  549)    !!======================================================================
(  550) 
(  551) 
(  552) END MODULE step






















PGF90 (Version     10.5)          09/23/2010  09:50:22      page 11

(    1) # 1 "step.F90"
(    1) MODULE step
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE step  ***
(    4)    !! Time-stepping    : manager of the ocean, tracer and ice time stepping
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   stp            : OPA system time-stepping
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables 
(   13)    USE zdf_oce         ! ocean vertical physics variables
(   14)    USE ldftra_oce
(   15)    USE ldfdyn_oce
(   16)    USE cpl_oce         ! coupled ocean-atmosphere variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk
(   19) 
(   20)    USE daymod          ! calendar                         (day     routine)
(   21) 
(   22)    USE dtatem          ! ocean temperature data           (dta_tem routine)
(   23)    USE dtasal          ! ocean salinity    data           (dta_sal routine)
(   24)    USE dtasst          ! ocean sea surface temperature    (dta_sst routine)
(   25)    USE dtasss          ! ocean sea surface salinity       (dta_sss routine)
(   26)    USE taumod          ! surface stress                   (tau     routine)
(   27)    USE flxmod          ! thermohaline fluxes              (flx     routine)
(   28)    USE ocesbc          ! thermohaline fluxes              (oce_sbc routine)
(   29)    USE flxrnf          ! runoffs                          (flx_rnf routine)
(   30)    USE flxfwb          ! freshwater budget correction     (flx_fwb routine)
(   31)    USE closea          ! closed sea freshwater budget     (flx_clo routine)
(   32)    USE ocfzpt          ! surface ocean freezing point    (oc_fz_pt routine)
(   33) 
(   34)    USE trcstp          ! passive tracer time-stepping      (trc_stp routine)
(   35) 
(   36)    USE dynhpg          ! hydrostatic pressure grad.       (dyn_hpg routine)
(   37)    USE dynhpg_atsk     ! hydrostatic pressure grad.  (dyn_hpg_atsk routine)
(   38)    USE dynspg_oce      ! surface pressure gradient        (dyn_spg routine)
(   39)    USE dynspg          ! surface pressure gradient        (dyn_spg routine)
(   40)    USE dynkeg          ! kinetic energy gradient          (dyn_keg routine)
(   41)    USE dynvor          ! vorticity term              (dyn_vor_... routines)
(   42)    USE dynzad          ! vertical advection               (dyn_adv routine)
(   43)    USE dynldf_bilapg   ! lateral mixing            (dyn_ldf_bilapg routine)
(   44)    USE dynldf_bilap    ! lateral mixing             (dyn_ldf_bilap routine)
(   45)    USE dynldf_iso      ! lateral mixing               (dyn_ldf_iso routine)
(   46)    USE dynldf_lap      ! lateral mixing               (dyn_ldf_lap routine)
(   47)    USE dynzdf_imp      ! vertical diffusion: implicit     (dyn_zdf routine)
(   48)    USE dynzdf_imp_atsk ! vertical diffusion: implicit     (dyn_zdf routine)
(   49)    USE dynzdf_iso      ! vertical diffusion: isopycnal    (dyn_zdf routine)
(   50)    USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp routine)
(   51)    USE dynnxt          ! time-stepping                    (dyn_nxt routine)
(   52) 
(   53)    USE trabbc          ! bottom boundary condition        (tra_bbc routine)
(   54)    USE trabbl          ! bottom boundary layer            (tra_bbl routine)
(   55)    USE tradmp          ! internal damping                 (tra_dmp routine)
(   56)    USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
(   57)    USE traldf_bilap    ! lateral mixing             (tra_ldf_bilap routine)






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 12

(   58)    USE traldf_iso      ! lateral mixing               (tra_ldf_iso routine)
(   59)    USE traldf_iso_zps  ! lateral mixing           (tra_ldf_iso_zps routine)
(   60)    USE traldf_lap      ! lateral mixing               (tra_ldf_lap routine)
(   61)    USE traqsr          ! solar radiation penetration      (tra_qsr routine)
(   62)    USE tranpc          ! non-penetrative convection       (tra_npc routine)
(   63)    USE tranxt          ! time-stepping                    (tra_nxt routine)
(   64)    USE traadv_ctl      ! advection scheme control     (tra_adv_ctl routine)
(   65)    USE traadv_cen2     ! 2nd order centered scheme   (tra_adv_cen2 routine)
(   66)    USE traadv_tvd      ! TVD scheme                (tra_adv_tvd    routine)
(   67)    USE traadv_muscl    ! MUSCL scheme              (tra_adv_muscl  routine)
(   68)    USE traadv_muscl2   ! MUSCL2 scheme             (tra_adv_muscl2 routine)
(   69) !   USE cla             ! cross land advection             (tra_cla routine)
(   70)    USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp routine)
(   71)    USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp routine)
(   72)    USE trazdf_iso      ! vertical diffusion           (tra_zdf_exp routine)
(   73)    USE trazdf_iso_vopt ! vertical diffusion           (tra_zdf_exp routine)
(   74)    USE trasbc          ! surface boundary condition       (tra_sbc routine)
(   75) 
(   76)    USE eosbn2          ! equation of state                (eos_bn2 routine)
(   77) 
(   78)    USE obc_par         ! open boundary condition variables
(   79)    USE obcdta          ! open boundary condition data     (obc_dta routine)
(   80)    USE obcrst          ! open boundary cond. restart      (obc_rst routine)
(   81)    USE obcrad          ! open boundary cond. radiation    (obc_rad routine)
(   82)    USE obcspg          ! open boundary cond  spg          (obc_spg routine)
(   83) 
(   84)    USE divcur          ! hor. divergence and curl      (div & cur routines)
(   85) !   USE cla_div         ! cross land: hor. divergence      (div_cla routine)
(   86)    USE wzvmod          ! vertical velocity                (wzv     routine)
(   87) 
(   88)    USE ldfslp          ! iso-neutral slopes               (ldf_slp routine)
(   89)    USE ldfeiv          ! eddy induced velocity coef.      (ldf_eiv routine)
(   90)    USE ldfdyn          ! 
(   91)    USE ldftra          
(   92) 
(   93) 
(   94)    USE zdfbfr          ! bottom friction                  (zdf_bfr routine)
(   95)    USE zdftke          ! TKE vertical mixing              (zdf_tke routine)
(   96)    USE zdfkpp          ! KPP vertical mixing              (zdf_kpp routine)
(   97)    USE zdfddm          ! double diffusion mixing          (zdf_ddm routine)
(   98)    USE zdfevd          ! enhanced vertical diffusion      (zdf_evd routine)
(   99)    USE zdfric          ! Richardson vertical mixing       (zdf_ric routine)
(  100)    USE zdfmxl          ! Mixed-layer depth                (zdf_mxl routine)
(  101) 
(  102)    USE zpshde          ! partial step: hor. derivative     (zps_hde routine)
(  103)    USE ice_oce         ! sea-ice variable
(  104)    USE icestp          ! sea-ice time-stepping             (ice_stp routine)
(  105) 
(  106)    USE diawri          ! Standard run outputs             (dia_wri routine)
(  107)    USE trdicp          ! Ocean momentum/tracers trends    (trd_wri routine)
(  108)    USE trdmld          ! mixed-layer trends               (trd_mld routine)
(  109)    USE trdvor          ! vorticity budget                 (trd_vor routine)
(  110) !   USE diagap          ! hor. mean model-data gap         (dia_gap routine)
(  111)    USE diahdy          ! dynamic height                   (dia_hdy routine)
(  112)    USE diaptr          ! poleward transports              (dia_ptr routine)
(  113)    USE diahth          ! thermocline depth                (dia_hth routine)
(  114)    USE diafwb          ! freshwater budget                (dia_fwb routine)
(  115)    USE diaspr          ! suface pressure (rigid-lid)      (dia_spr routine)






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 13

(  116)    USE flo_oce         ! floats variables
(  117)    USE floats          ! floats computation               (flo_stp routine)
(  118) 
(  119)    USE stpctl          ! time stepping control            (stp_ctl routine)
(  120)    USE restart         ! ocean restart                    (rst_wri routine)
(  121)    USE cpl             ! exchanges in coupled mode        (cpl_stp routine)
(  122)    USE prtctl          ! Print control                    (prt_ctl routine)
(  123) 
(  124) !!DB 2009.08.06
(  125)    USE sopa_mc
(  126) 
(  127) 
(  128) 
(  129) 
(  130) # 133
(  133)    IMPLICIT NONE
(  134)    PRIVATE
(  135) 
(  136)    !! * Routine accessibility
(  137)    PUBLIC stp            ! called by opa.F90
(  138) 
(  139)    !! * Substitutions
(  140) # 1 "./domzgr_substitute.h90"
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
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 141 "step.F90"
(  141) # 1 "./zdfddm_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** zdfddm_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm_substitute.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 16
(   16) !   Defautl option :                     avs = avt
(   17) # 142 "step.F90"
(  142) # 142
(  142)    !!----------------------------------------------------------------------
(  143)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(  144)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/step.F90,v 1.25 2006/04/26 09:26:48 opalod Exp $ 
(  145)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 14

(  146)    !!----------------------------------------------------------------------
(  147) 
(  148) CONTAINS
(  149) 
(  150)    SUBROUTINE stp( &
(  152)    kstp &
(  154)    )      !!----------------------------------------------------------------------
(  155)       !!                     ***  ROUTINE stp  ***
(  156)       !!                      
(  157)       !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
(  158)       !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
(  159)       !!              - Tme stepping  of TRC (passive tracer eqs.)
(  160)       !! 
(  161)       !! ** Method  : -1- Update forcings and data  
(  162)       !!              -2- Update ocean physics 
(  163)       !!              -3- Compute the t and s trends 
(  164)       !!              -4- Update t and s 
(  165)       !!              -5- Compute the momentum trends
(  166)       !!              -6- Update the horizontal velocity
(  167)       !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
(  168)       !!              -8- Outputs and diagnostics
(  169)       !!
(  170)       !! History :
(  171)       !!        !  91-03  ()  Original code
(  172)       !!        !  91-11  (G. Madec)
(  173)       !!        !  92-06  (M. Imbard)  add a first output record
(  174)       !!        !  96-04  (G. Madec)  introduction of dynspg
(  175)       !!        !  96-04  (M.A. Foujols)  introduction of passive tracer
(  176)       !!   8.0  !  97-06  (G. Madec)  new architecture of call
(  177)       !!   8.2  !  97-06  (G. Madec, M. Imbard, G. Roullet)  free surface
(  178)       !!   8.2  !  99-02  (G. Madec, N. Grima)  hpg implicit
(  179)       !!   8.2  !  00-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
(  180)       !!   9.0  !  02-06  (G. Madec)  free form, suppress macro-tasking
(  181)       !!    "   !  04-08  (C. Talandier) New trends organization
(  182)       !!    "   !  05-01  (C. Ethe) Add the KPP closure scheme
(  183)       !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
(  184)       !!----------------------------------------------------------------------
(  185)       !! * Arguments
(  186)       INTEGER &
(  188)       , INTENT( in ) &
(  190)       ::   kstp   ! ocean time-step index
(  191) 
(  192)       !! * local declarations
(  193)       INTEGER ::   indic    ! error indicator if < 0
(  194)       !! ---------------------------------------------------------------------
(  195)       INTEGER :: ji, jj, jk 
(  196) 
(  197) 
(  198) 
(  199) # 202
(  202)       indic = 1                    ! reset to no error condition
(  203)       adatrj = adatrj + rdt/86400._wp
(  204) 
(  205) !!DB 2008.04.16 -- hardwired ramp function currently used in obcdta, tau_forced_*
(  206)      ramp=tanh(kstp*rdt/(2.0*86400.0))
(  207) 
(  208) !!DB: perpetual forcing on ====> perpetual_forcing /= 0 ====> do NOT call day()
(  209)       if(perpetual_forcing == 0) then






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 15

(  210)          CALL day( kstp )             ! Calendar
(  211)       endif
(  212) 
(  213) !!DBG: 2009.06.12
(  214) !      call OB_LIMITER(kstp,1)
(  215) 
(  216) 
(  217) 
(  218)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  219)       ! Update data, open boundaries and Forcings
(  220)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  221) 
(  222)       IF( lk_dtatem  )   CALL dta_tem( kstp )         ! update 3D temperature data
(  223) 
(  224)       IF( lk_dtasal  )   CALL dta_sal( kstp )         ! Salinity data
(  225) 
(  226)       IF( lk_dtasst  )   CALL dta_sst( kstp )         ! Sea Surface Temperature data
(  227) 
(  228)       IF( lk_dtasss  )   CALL dta_sss( kstp )         ! Sea Surface salinity data
(  229) 
(  230)       IF( lk_obc     )   CALL obc_dta( kstp )         ! update dynamic and tracer data at open boundaries
(  231) 
(  232)       IF( lk_obc     )   CALL obc_rad( kstp )         ! compute phase velocities at open boundaries
(  233)       
(  234)       CALL tau( kstp )             ! wind stress
(  235)       
(  236)       CALL flx_rnf( kstp )         ! runoff data
(  237)       
(  238)       CALL flx( kstp )             ! heat and freshwater fluxes
(  239) 
(  240)       IF( lk_ice_lim )   CALL ice_stp( kstp )         ! sea-ice model (Update stress & fluxes)
(  241)       
(  242)       CALL oce_sbc( kstp )         ! ocean surface boundaries
(  243)       
(  244)       IF( ln_fwb     )   CALL flx_fwb( kstp )         ! freshwater budget
(  245)       
(  246)       IF( nclosea == 1 ) CALL flx_clo( kstp )         ! closed sea in the domain (update freshwater fluxes)
(  247)       
(  248)       IF( kstp == nit000 ) THEN 
(  249)          IF( ninist == 1 ) THEN                       ! Output the initial state and forcings
(  250)             CALL dia_wri_state( 'output.init' )
(  251)          ENDIF
(  252)       ENDIF
(  253)       
(  254)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  255)          CALL prt_ctl(tab2d_1=emp    , clinfo1=' emp  -   : ', mask1=tmask, ovlap=1)
(  256)          CALL prt_ctl(tab2d_1=emps   , clinfo1=' emps -   : ', mask1=tmask, ovlap=1)
(  257)          CALL prt_ctl(tab2d_1=qt     , clinfo1=' qt   -   : ', mask1=tmask, ovlap=1)
(  258)          CALL prt_ctl(tab2d_1=qsr    , clinfo1=' qsr  -   : ', mask1=tmask, ovlap=1)
(  259)          CALL prt_ctl(tab2d_1=runoff , clinfo1=' runoff   : ', mask1=tmask, ovlap=1)
(  260)          CALL prt_ctl(tab3d_1=tmask  , clinfo1=' tmask    : ', mask1=tmask, ovlap=1, kdim=jpk)
(  261)          CALL prt_ctl(tab3d_1=tn     , clinfo1=' sst  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  262)          CALL prt_ctl(tab3d_1=sn     , clinfo1=' sss  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  263)          CALL prt_ctl(tab2d_1=taux   , clinfo1=' tau  - x : ', tab2d_2=tauy, clinfo2='      - y : ', ovlap=1)
(  264)       ENDIF
(  265) 
(  266) 
(  267)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 16

(  268)       ! Ocean physics update
(  269)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  270)       !-----------------------------------------------------------------------
(  271)       !  VERTICAL PHYSICS
(  272)       !-----------------------------------------------------------------------
(  273)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  274)       !-----------------------------------------------------------------------
(  275)       
(  276)       CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
(  277)       
(  278)       !                                                     ! Vertical eddy viscosity and diffusivity coefficients
(  279)       IF( lk_zdfric )   CALL zdf_ric( kstp )                       ! Richardson number dependent Kz
(  280)       IF( lk_zdftke )   CALL zdf_tke( kstp )                       ! TKE closure scheme for Kz
(  281)       IF( lk_zdfkpp )   CALL zdf_kpp( kstp )                       ! KPP closure scheme for Kz
(  282)       IF( lk_zdfcst )   avt (:,:,:) = avt0 * tmask(:,:,:)          ! Constant Kz (reset avt to the background value)
(  283) 
(  284)       IF( ln_zdfevd )   CALL zdf_evd( kstp )                 ! enhanced vertical eddy diffusivity
(  285) 
(  286)       IF( lk_zdfddm .AND. .NOT. lk_zdfkpp)   &
(  287)            &              CALL zdf_ddm( kstp )                 ! double diffusive mixing
(  288) 
(  289)       CALL zdf_bfr( kstp )                 ! bottom friction
(  290)       
(  291)       CALL zdf_mxl( kstp )                 ! mixed layer depth
(  292) 
(  293) 
(  294)       !-----------------------------------------------------------------------
(  295)       !  LATERAL PHYSICS
(  296)       !-----------------------------------------------------------------------
(  297)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  298)       !-----------------------------------------------------------------------
(  299) 
(  300)       IF( lk_ldfslp     )   CALL ldf_slp( kstp, rhd, rn2 )       ! before slope of the lateral mixing
(  301) 
(  302) 
(  303) 
(  304) 
(  305) 
(  306) # 319
(  319)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  320)       ! Active tracers
(  321)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  322)       ! N.B. ua, va arrays are used as workspace in this section
(  323)       !-----------------------------------------------------------------------
(  324)       
(  325)       ta(:,:,:) = 0.e0               ! set tracer trends to zero
(  326)       sa(:,:,:) = 0.e0
(  327)       
(  328)       CALL tra_sbc( kstp )           ! surface boundary condition
(  329) 
(  330)       IF( ln_traqsr        )   CALL tra_qsr( kstp )           ! penetrative solar radiation qsr
(  331) 
(  332)       IF( lk_trabbc        )   CALL tra_bbc( kstp )           ! bottom heat flux
(  333)       
(  334)       IF( lk_trabbl_dif    )   CALL tra_bbl_dif( kstp )           ! diffusive bottom boundary layer scheme
(  335)       IF( lk_trabbl_adv    )   CALL tra_bbl_adv( kstp )           ! advective (and/or diffusive) bottom boundary layer scheme
(  336)       
(  337)       IF( lk_tradmp        )   CALL tra_dmp( kstp )           ! internal damping trends






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 17

(  338) 
(  339)       !                                                       ! horizontal & vertical advection
(  340)       IF( kstp == nit000   )   CALL tra_adv_ctl                    ! chose/control the scheme used
(  341)       IF( ln_traadv_cen2   )   CALL tra_adv_cen2  ( kstp )         ! 2nd order centered scheme
(  342)       IF( ln_traadv_tvd    )   CALL tra_adv_tvd   ( kstp )         ! TVD scheme
(  343)       IF( ln_traadv_muscl  )   CALL tra_adv_muscl ( kstp )         ! MUSCL scheme
(  344)       IF( ln_traadv_muscl2 )   CALL tra_adv_muscl2( kstp )         ! MUSCL2 scheme
(  345) 
(  346) !!DB: orca-related
(  347) !      IF( n_cla == 1       )   CALL tra_cla( kstp )           ! Cross Land Advection (Update Hor. advection)
(  348) 
(  349)       !                                                       ! lateral mixing 
(  350)       IF( l_traldf_lap     )   CALL tra_ldf_lap    ( kstp )           ! iso-level laplacian
(  351)       IF( l_traldf_bilap   )   CALL tra_ldf_bilap  ( kstp )           ! iso-level bilaplacian 
(  352)       IF( l_traldf_bilapg  )   CALL tra_ldf_bilapg ( kstp )           ! s-coord. horizontal bilaplacian
(  353)       IF( l_traldf_iso     )   CALL tra_ldf_iso    ( kstp )           ! iso-neutral/geopot. laplacian 
(  354)       IF( l_traldf_iso_zps )   CALL tra_ldf_iso_zps( kstp )           ! partial step iso-neutral/geopot. laplacian
(  355) 
(  356) # 359
(  359)       !                                                       ! vertical diffusion
(  360)       IF( l_trazdf_exp     )   CALL tra_zdf_exp     ( kstp )          ! explicit time stepping (time splitting scheme)
(  361)       IF( l_trazdf_imp     )   CALL tra_zdf_imp     ( kstp )          ! implicit time stepping (euler backward)
(  362)       IF( l_trazdf_iso     )   CALL tra_zdf_iso     ( kstp )          ! isopycnal
(  363)       IF( l_trazdf_iso_vo  )   CALL tra_zdf_iso_vopt( kstp )          ! vector opt. isopycnal
(  364)       
(  365)       CALL tra_nxt( kstp )           ! tracer fields at next time step
(  366)       
(  367)       IF( ln_zdfnpc        )   CALL tra_npc( kstp )           ! update the new (t,s) fields by non
(  368)       !                                                       ! penetrative convective adjustment
(  369)       
(  370)       IF( ln_dynhpg_imp    ) THEN                             ! semi-implicit hpg 
(  371)          CALL eos( ta, sa, rhd, rhop )   ! Time-filtered in situ density used in dynhpg module
(  372)          IF( lk_zps    )          CALL zps_hde( kstp, ta, sa, rhd,  & ! Partial steps: time filtered hor. gradient 
(  373)               &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
(  374)               &                                        gtv, gsv, grv )  
(  375)       ELSE                                                    ! centered hpg (default case)
(  376)          CALL eos( tb, sb, rhd, rhop )       ! now (swap=before) in situ density for dynhpg module
(  377)          IF( lk_zps    )          CALL zps_hde( kstp, tb, sb, rhd,  & ! Partial steps: now horizontal gradient
(  378)               &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
(  379)               &                                        gtv, gsv, grv )  
(  380)       ENDIF
(  381)       
(  382)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  383)       ! Dynamics
(  384)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  385)       ! N.B. ta, sa arrays are used as workspace in this section 
(  386)       !-----------------------------------------------------------------------
(  387) 
(  388)       ua(:,:,:) = 0.e0               ! set dynamics trends to zero
(  389)       va(:,:,:) = 0.e0
(  390) 
(  391)       CALL dyn_keg( kstp )           ! horizontal gradient of kinetic energy
(  392) 
(  393)       !                                                       ! vorticity term including Coriolis
(  394)       IF( kstp == nit000   )   CALL dyn_vor_ctl                      ! chose/control the scheme used
(  395)       IF( ln_dynvor_ens    )   CALL dyn_vor_enstrophy( kstp )        ! enstrophy conserving scheme
(  396)       IF( ln_dynvor_ene    )   CALL dyn_vor_energy   ( kstp )        ! energy conserving scheme
(  397)       IF( ln_dynvor_mix    )   CALL dyn_vor_mixed    ( kstp )        ! mixed energy/enstrophy conserving scheme






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 18

(  398)       IF( ln_dynvor_een    )   CALL dyn_vor_ene_ens  ( kstp )        ! combined energy/enstrophy conserving scheme
(  399)       
(  400)       !                                                       ! lateral mixing 
(  401)       IF( l_dynldf_lap     )   CALL dyn_ldf_lap    ( kstp )          ! iso-level laplacian
(  402)       IF( l_dynldf_bilap   )   CALL dyn_ldf_bilap  ( kstp )          ! iso-level bilaplacian 
(  403)       IF( l_dynldf_bilapg  )   CALL dyn_ldf_bilapg ( kstp )          ! s-coord. horizontal bilaplacian
(  404)       IF( l_dynldf_iso     )   CALL dyn_ldf_iso    ( kstp )          ! iso-neutral laplacian 
(  405)       
(  406) # 409
(  409)       !                                                       ! horizontal gradient of Hydrostatic pressure 
(  410)       IF ( lk_jki ) THEN
(  411)          CALL dyn_hpg_atsk( kstp )             ! autotask case (j-k-i loop)
(  412)       ELSE
(  413)          CALL dyn_hpg     ( kstp )             ! default case  (k-j-i loop)
(  414)       ENDIF
(  415)       
(  416)       CALL dyn_zad    ( kstp )       ! vertical advection       
(  417)       
(  418)       !                                                       ! vertical diffusion
(  419)       IF( l_dynzdf_exp     )   CALL dyn_zdf_exp    ( kstp )          ! explicit time stepping (time splitting scheme)
(  420)       IF( l_dynzdf_imp     )   CALL dyn_zdf_imp    ( kstp )          ! implicit time stepping (euler backward)
(  421)       IF( l_dynzdf_imp_tsk )   CALL dyn_zdf_imp_tsk( kstp )          ! autotask implicit time stepping (euler backward)
(  422)       IF( l_dynzdf_iso     )   CALL dyn_zdf_iso    ( kstp )          ! iso-neutral case
(  423)       
(  424)       IF( lk_dynspg_rl ) THEN 
(  425)          IF( lk_obc    )       CALL obc_spg( kstp )           ! surface pressure gradient at open boundaries
(  426)       ENDIF
(  427)       indic=0
(  428)       !i bug lbc sur emp
(  429)       CALL lbc_lnk( emp, 'T', 1. )
(  430)       !i
(  431)       CALL dyn_spg( kstp, indic )    ! surface pressure gradient
(  432)       
(  433)       CALL dyn_nxt( kstp )           ! velocity at next time step 
(  434)       
(  435)       
(  436)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  437)       ! Computation of diagnostic variables
(  438)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  439)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  440)       !-----------------------------------------------------------------------
(  441)       
(  442)       CALL oc_fz_pt                        ! ocean surface freezing temperature
(  443)       
(  444)       CALL div_cur( kstp )                 ! Horizontal divergence & Relative vorticity
(  445) 
(  446) !!DB: orca-related      
(  447) !      IF( n_cla == 1 ) CALL div_cla( kstp )                 ! Cross Land Advection (Update Hor. divergence)
(  448)       
(  449)       CALL wzv( kstp )                     ! Vertical velocity
(  450)       
(  451) !!DB 2008.04.07 -- the below ultimately calls the appropriate routine to update the
(  452) !! relevant diffusivity coeffs (see ldfdyn.F90). The prob is that it performs numerous
(  453) !! calcs every timestep (including reading namelist) that are not necessary
(  454) !! =========> modification ... DONE ...
(  455) 
(  456) !      !ZW
(  457) !      CALL ldf_dyn_init






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 19

(  458) !      CALL ldf_tra_init
(  459) !      !ZW
(  460) 
(  461) !!DB 04.09 -- modify above so that init routines are not called every dt
(  462) !!Note that the smag keys are separated to the end of this code fragment
(  463) !!Also, there is a potential problem if there is a mistake in the keys
(  464) !!so that (e.g.) no or more-than-one tracer key is defined. I do not check
(  465) !!for this possibility
(  466) !!Momentum diffusivity updates
(  467) # 474
(  474)       !do nothing unless smag is on -- see below
(  475) # 476
(  476) !!tracer diffusivity updates
(  477) # 484
(  484)       !do nothing unless smag is on -- see below
(  485) 
(  486) # 488
(  488)       call ldf_smag( kstp )
(  489) 
(  490) # 491
(  491)       
(  492)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  493)       ! Control, diagnostics and outputs
(  494)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  495)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  496)       !-----------------------------------------------------------------------
(  497)       
(  498)       !                                            ! Time loop: control and print
(  499)       CALL stp_ctl( kstp, indic )
(  500)       IF ( indic < 0 )   nstop = nstop + 1
(  501) 
(  502) 
(  503)       IF ( nstop == 0 ) THEN
(  504)          !                                         ! Diagnostics:
(  505)          IF( lk_floats  )   CALL flo_stp( kstp )                 ! drifting Floats
(  506)          IF( lk_trddyn  )   CALL trd_dwr( kstp )                 ! trends: dynamics 
(  507)          IF( lk_trdtra  )   CALL trd_twr( kstp )                 ! trends: active tracers
(  508)          IF( lk_trdmld  )   CALL trd_mld( kstp )                 ! trends: Mixed-layer 
(  509)          IF( lk_trdvor  )   CALL trd_vor( kstp )                 ! trends: vorticity budget
(  510)          IF( lk_diaspr  )   CALL dia_spr( kstp )                 ! Surface pressure diagnostics
(  511)          IF( lk_diahth  )   CALL dia_hth( kstp )                 ! Thermocline depth (20 degres isotherm depth)
(  512) !         IF( lk_diagap  )   CALL dia_gap( kstp )                 ! basin averaged diagnostics
(  513)          IF( lk_diahdy  )   CALL dia_hdy( kstp )                 ! dynamical heigh diagnostics
(  514)          IF( lk_diafwb  )   CALL dia_fwb( kstp )                 ! Fresh water budget diagnostics
(  515)          IF( ln_diaptr  )   CALL dia_ptr( kstp )                 ! Poleward TRansports diagnostics
(  516)          
(  517)          !                                         ! save and outputs
(  518)          CALL rst_write  ( kstp )             ! ocean model: restart file output
(  519) !!DB
(  520) !         IF( lk_obc     )   CALL obc_rst_wri( kstp )             ! ocean model: open boundary restart file output
(  521)          
(  522) !!DB 2009.08.06 -- These are the *_grid_*.nc files which we never look at
(  523) !! (Replaced by M2 aves for U,V and *_aveTSUV.nc (see below) for other variables
(  524) !         CALL dia_wri    ( kstp, indic )      ! ocean model: outputs
(  525) 
(  526) !!DB
(  527)          if(M2_ave > 0) call output_special(kstp, indic)      ! special M2 time-averaged fields
(  528) !!DB






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 20

(  529)          if(ioutput_ave /= 0) call output_aveTSUV(kstp, indic)      ! special time-averaged fields
(  530) 
(  531)          
(  532)       ENDIF
(  533) 
(  534) !DB Output some BoF vels -- NB: code moved to sopa_mc module
(  535) !! Model diagnostics (somewhat) specific to sopa MC domain 
(  536)       call sopa_mc_diagnostics(kstp)
(  537) 
(  538) 
(  539) 
(  540) 
(  541)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  542)       ! Coupled mode
(  543)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  544) 
(  545)       IF( lk_cpl    )   CALL cpl_stp( kstp )                 ! coupled mode : field exchanges
(  546) 
(  547)    END SUBROUTINE stp
(  548) 
(  549)    !!======================================================================
(  550) 
(  551) 
(  552) END MODULE step
