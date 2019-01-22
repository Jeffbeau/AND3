


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

Filename: step1d.F90

(    1) MODULE step1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE step1D  ***
(    4)    !! Time-stepping    : manager of the ocean, tracer and ice time stepping
(    5)    !!======================================================================
(    6) #if defined key_cfg_1d
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_cfg_1d'               1D Configuration
(    9)    !!----------------------------------------------------------------------  
(   10)    !!----------------------------------------------------------------------
(   11)    !!   stp_1d           : OPA system time-stepping on 1 direction
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables 
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE ldftra_oce
(   18)    USE ldfdyn_oce
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE lbclnk
(   21) 
(   22)    USE daymod          ! calendar                         (day     routine)
(   23) 
(   24)    USE dtatem          ! ocean temperature data           (dta_tem routine)
(   25)    USE dtasal          ! ocean salinity    data           (dta_sal routine)
(   26)    USE dtasst          ! ocean sea surface temerature     (dta_sst routine)
(   27)    USE taumod          ! surface stress                   (tau     routine)
(   28)    USE flxmod          ! thermohaline fluxes              (flx     routine)
(   29)    USE ocesbc          ! thermohaline fluxes              (oce_sbc routine)
(   30)    USE flxrnf          ! runoffs                          (flx_rnf routine)
(   31)    USE flxfwb          ! freshwater budget correction     (flx_fwb routine)
(   32)    USE ocfzpt          ! surface ocean freezing point    (oc_fz_pt routine)
(   33) 
(   34)    USE trcstp          ! passive tracer time-stepping     (trc_stp routine)
(   35) 
(   36)    USE dynzdf_imp      ! vertical diffusion: implicit     (dyn_zdf routine)
(   37)    USE dynzdf_imp_atsk ! vertical diffusion: implicit     (dyn_zdf routine)
(   38)    USE dynzdf_iso      ! vertical diffusion: isopycnal    (dyn_zdf routine)
(   39)    USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp routine)
(   40)  
(   41) 
(   42)    USE traqsr          ! solar radiation penetration      (tra_qsr routine)
(   43)    USE tranxt          ! time-stepping                    (tra_nxt routine)
(   44)    USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp routine)
(   45)    USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp routine)
(   46)    USE trazdf_iso      ! vertical diffusion           (tra_zdf_exp routine)






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 2

(   47)    USE trazdf_iso_vopt ! vertical diffusion           (tra_zdf_exp routine)
(   48)    USE trasbc          ! surface boundary condition       (tra_sbc routine)
(   49) 
(   50)    USE eosbn2
(   51) 
(   52)    USE zdfbfr          ! bottom friction                  (zdf_bfr routine)
(   53)    USE zdftke          ! TKE vertical mixing              (zdf_tke routine)
(   54)    USE zdfkpp          ! KPP vertical mixing              (zdf_kpp routine)
(   55)    USE zdfddm          ! double diffusion mixing          (zdf_ddm routine)
(   56)    USE zdfevd          ! enhanced vertical diffusion      (zdf_evd routine)
(   57)    USE zdfric          ! Richardson vertical mixing       (zdf_ric routine)
(   58)    USE zdfmxl          ! Mixed-layer depth                (zdf_mxl routine)
(   59) 
(   60)    USE dyncor1d
(   61)    USE dynnxt1d
(   62)    USE diawri1d        ! Standard run outputs             (dia_wri_1d routine)
(   63) 
(   64)    USE ice_oce         ! sea-ice variable
(   65)    USE icestp1d        ! sea-ice time-stepping             (ice_stp routine)
(   66) 
(   67)    USE diawri          ! Standard run outputs             (dia_wri_state routine)
(   68) 
(   69) 
(   70)    USE stpctl          ! time stepping control            (stp_ctl routine)
(   71)    USE restart         ! ocean restart                    (rst_wri routine)
(   72)    USE prtctl          ! Print control                    (prt_ctl routine)
(   73)    IMPLICIT NONE
(   74)    PRIVATE
(   75) 
(   76)    !! * Routine accessibility
(   77)    PUBLIC stp_1d            ! called by opa.F90
(   78) 
(   79)    !! * Substitutions
(   80) #  include "domzgr_substitute.h90"
(   81) #  include "zdfddm_substitute.h90"
(   82)    !!----------------------------------------------------------------------
(   83)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   84)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/step1d.F90,v 1.5 2006/04/11 13:52:28 opalod Exp $ 
(   85)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   86)    !!----------------------------------------------------------------------
(   87) 
(   88) CONTAINS
(   89) 
(   90)    SUBROUTINE stp_1d( kstp )
(   91)       !!----------------------------------------------------------------------
(   92)       !!                     ***  ROUTINE stp1D  ***
(   93)       !!                      
(   94)       !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
(   95)       !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
(   96)       !! 
(   97)       !! ** Method  : -1- Update forcings and data  
(   98)       !!              -2- Update ocean physics 
(   99)       !!              -3- Compute the t and s trends 
(  100)       !!              -4- Update t and s 
(  101)       !!              -5- Compute the momentum trends
(  102)       !!              -6- Update the horizontal velocity
(  103)       !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
(  104)       !!              -8- Outputs and diagnostics






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 3

(  105)       !!
(  106)       !! History :
(  107)       !!        !  91-03  ()  Original code
(  108)       !!        !  91-11  (G. Madec)
(  109)       !!        !  92-06  (M. Imbard)  add a first output record
(  110)       !!        !  96-04  (G. Madec)  introduction of dynspg
(  111)       !!        !  96-04  (M.A. Foujols)  introduction of passive tracer
(  112)       !!   8.0  !  97-06  (G. Madec)  new architecture of call
(  113)       !!   8.2  !  97-06  (G. Madec, M. Imbard, G. Roullet)  free surface
(  114)       !!   8.2  !  99-02  (G. Madec, N. Grima)  hpg implicit
(  115)       !!   8.2  !  00-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
(  116)       !!   9.0  !  02-06  (G. Madec)  free form, suppress macro-tasking
(  117)       !!        !  04-10  (C. Ethe) 1D configuration
(  118)       !!----------------------------------------------------------------------
(  119)       !! * Arguments
(  120)       INTEGER, INTENT( in ) ::   kstp   ! ocean time-step index
(  121) 
(  122)       !! * local declarations
(  123)       INTEGER ::   indic    ! error indicator if < 0
(  124) !!      INTEGER ::   ii0, ii1, ij0, ij1   ! temporary integers
(  125)       !! ---------------------------------------------------------------------
(  126) 
(  127)       indic = 1                    ! reset to no error condition
(  128)       adatrj = adatrj + rdt/86400._wp
(  129) 
(  130)       CALL day( kstp )             ! Calendar
(  131) 
(  132)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  133)       ! Update data, open boundaries and Forcings
(  134)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  135) 
(  136)       IF( lk_dtatem  )   CALL dta_tem( kstp )         ! update 3D temperature data
(  137) 
(  138)       IF( lk_dtasal  )   CALL dta_sal( kstp )         ! Salinity data
(  139) 
(  140)       IF( lk_dtasst  )   CALL dta_sst( kstp )         ! Sea Surface Temperature data
(  141) 
(  142)                          CALL tau( kstp )             ! wind stress
(  143) 
(  144)                          CALL flx_rnf( kstp )         ! runoff data
(  145) 
(  146)                          CALL flx( kstp )             ! heat and freshwater fluxes
(  147) 
(  148)       IF( lk_ice_lim )  THEN 
(  149)                         CALL ice_stp_1d( kstp )      ! sea-ice model (Update stress & fluxes)
(  150)       ELSE
(  151)                         CALL oce_sbc( kstp )         ! ocean surface boudaries
(  152)       ENDIF
(  153) 
(  154)       IF( ln_fwb     )   CALL flx_fwb( kstp )         ! freshwater budget
(  155) 
(  156) 
(  157)       IF( kstp == nit000 ) THEN 
(  158)          IF( ninist == 1 ) THEN                       ! Output the initial state and forcings
(  159)             CALL dia_wri_state( 'output.init' )
(  160)          ENDIF
(  161)       ENDIF
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 4

(  163)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  164)          CALL prt_ctl(tab2d_1=emp    , clinfo1=' emp  -   : ', mask1=tmask, ovlap=1)
(  165)          CALL prt_ctl(tab2d_1=emps   , clinfo1=' emps -   : ', mask1=tmask, ovlap=1)
(  166)          CALL prt_ctl(tab2d_1=qt     , clinfo1=' qt   -   : ', mask1=tmask, ovlap=1)
(  167)          CALL prt_ctl(tab2d_1=qsr    , clinfo1=' qsr  -   : ', mask1=tmask, ovlap=1)
(  168)          CALL prt_ctl(tab2d_1=runoff , clinfo1=' runoff   : ', mask1=tmask, ovlap=1)
(  169)          CALL prt_ctl(tab3d_1=tmask  , clinfo1=' tmask    : ', mask1=tmask, ovlap=1, kdim=jpk)
(  170)          CALL prt_ctl(tab3d_1=tn     , clinfo1=' sst  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  171)          CALL prt_ctl(tab3d_1=sn     , clinfo1=' sss  -   : ', mask1=tmask, ovlap=1, kdim=1)
(  172)          CALL prt_ctl(tab2d_1=taux   , clinfo1=' tau  - x : ', tab2d_2=tauy, clinfo2='      - y : ', ovlap=1)
(  173)       ENDIF
(  174) 
(  175) 
(  176) 
(  177)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  178)       ! Ocean physics update
(  179)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  180)       !-----------------------------------------------------------------------
(  181)       !  VERTICAL PHYSICS
(  182)       !-----------------------------------------------------------------------
(  183)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  184)       !-----------------------------------------------------------------------
(  185) 
(  186)                        CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
(  187)       
(  188)       !                                                     ! Vertical eddy viscosity and diffusivity coefficients
(  189)       IF( lk_zdfric )   CALL zdf_ric( kstp )                       ! Richardson number dependent Kz
(  190)       IF( lk_zdftke )   CALL zdf_tke( kstp )                       ! TKE closure scheme for Kz
(  191)       IF( lk_zdfkpp )   CALL zdf_kpp( kstp )                       ! KPP scheme for Kz
(  192)       IF( lk_zdfcst )   avt (:,:,:) = avt0 * tmask(:,:,:)          ! Constant Kz (reset avt to the background value)
(  193) 
(  194) 
(  195)       IF( ln_zdfevd )   CALL zdf_evd( kstp )                 ! enhanced vertical eddy diffusivity
(  196) 
(  197)       IF( lk_zdfddm .AND. .NOT. lk_zdfkpp)   &
(  198)          &              CALL zdf_ddm( kstp )                 ! double diffusive mixing
(  199) 
(  200)                         CALL zdf_bfr( kstp )                 ! bottom friction
(  201) 
(  202)                         CALL zdf_mxl( kstp )                 ! mixed layer depth
(  203) 
(  204) 
(  205) #if defined key_passivetrc
(  206)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  207)       ! Passive Tracer Model
(  208)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  209)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  210)       !-----------------------------------------------------------------------
(  211) 
(  212)                                CALL trc_stp( kstp, indic )            ! time-stepping
(  213) 
(  214) #endif
(  215) 
(  216) 
(  217)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  218)       ! Active tracers
(  219)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  220)       ! N.B. ua, va arrays are used as workspace in this section






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 5

(  221)       !-----------------------------------------------------------------------
(  222) 
(  223)                                ta(:,:,:) = 0.e0               ! set tracer trends to zero
(  224)                                sa(:,:,:) = 0.e0
(  225) 
(  226)                                CALL tra_sbc( kstp )           ! surface boundary condition
(  227) 
(  228)       IF( ln_traqsr        )   CALL tra_qsr( kstp )           ! penetrative solar radiation qsr
(  229) 
(  230)       !                                                       ! vertical diffusion
(  231)       IF( l_trazdf_exp     )   CALL tra_zdf_exp     ( kstp )          ! explicit time stepping (time splitting scheme)
(  232)       IF( l_trazdf_imp     )   CALL tra_zdf_imp     ( kstp )          ! implicit time stepping (euler backward)
(  233)       IF( l_trazdf_iso     )   CALL tra_zdf_iso     ( kstp )          ! isopycnal
(  234)       IF( l_trazdf_iso_vo  )   CALL tra_zdf_iso_vopt( kstp )          ! vector opt. isopycnal
(  235) 
(  236)                                CALL tra_nxt( kstp )           ! tracer fields at next time step
(  237) 
(  238)                                CALL eos( tb, sb, rhd, rhop )       ! now (swap=before) in situ density for dynhpg module
(  239) 
(  240)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  241)       ! Dynamics
(  242)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  243)       ! N.B. ta, sa arrays are used as workspace in this section 
(  244)       !-----------------------------------------------------------------------
(  245) 
(  246)                                ua(:,:,:) = 0.e0               ! set dynamics trends to zero
(  247)                                va(:,:,:) = 0.e0
(  248)   
(  249)                                CALL dyn_cor_1d     ( kstp )
(  250)       !                                                       ! vertical diffusion
(  251)       IF( l_dynzdf_exp     )   CALL dyn_zdf_exp    ( kstp )         ! explicit time stepping (time splitting scheme)
(  252)       IF( l_dynzdf_imp     )   CALL dyn_zdf_imp    ( kstp )         ! implicit time stepping (euler backward)
(  253)       IF( l_dynzdf_imp_tsk )   CALL dyn_zdf_imp_tsk( kstp )         ! autotask implicit time stepping (euler backward)
(  254)       IF( l_dynzdf_iso     )   CALL dyn_zdf_iso    ( kstp )         ! iso-neutral case
(  255) 
(  256) !i bug lbc sur emp
(  257)       CALL lbc_lnk( emp, 'T', 1. )
(  258) !i
(  259)  
(  260)                                 CALL dyn_nxt_1d( kstp )          ! lateral velocity at next time step 
(  261) 
(  262) 
(  263)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  264)       ! Computation of diagnostic variables
(  265)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  266)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  267)       !-----------------------------------------------------------------------
(  268) 
(  269)                                CALL oc_fz_pt                    ! ocean surface freezing temperature
(  270) 
(  271) 
(  272)       !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
(  273)       ! Control, diagnostics and outputs
(  274)       !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
(  275)       ! N.B. ua, va, ta, sa arrays are used as workspace in this section
(  276)       !-----------------------------------------------------------------------
(  277) 
(  278)       !                                            ! Time loop: control and print






PGF90 (Version     10.5)          09/23/2010  09:50:22      page 6

(  279)                        CALL stp_ctl( kstp, indic )
(  280)                        IF ( indic < 0 )   nstop = nstop + 1
(  281) 
(  282)       IF ( nstop == 0 ) THEN
(  283)          !                                         ! Diagnostics:
(  284)          !                                         ! save and outputs
(  285)                            CALL rst_write  ( kstp )              ! ocean model: restart file output
(  286)                            CALL dia_wri_1d ( kstp, indic )       ! ocean model: outputs
(  287) 
(  288)       ENDIF
(  289) 
(  290) 
(  291)    END SUBROUTINE stp_1d
(  292) #else
(  293)    !!----------------------------------------------------------------------
(  294)    !!   Default key                                     NO 1D Config
(  295)    !!----------------------------------------------------------------------
(  296) CONTAINS
(  297)    SUBROUTINE stp_1d ( kt )
(  298) !      WRITE(*,*) 'stp_1d: You should not have seen this print! error?', kt
(  299)    END SUBROUTINE stp_1d
(  300) #endif
(  301)    !!======================================================================
(  302) END MODULE step1d








































PGF90 (Version     10.5)          09/23/2010  09:50:22      page 7

(    1) # 1 "step1d.F90"
(    1) MODULE step1d
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE step1D  ***
(    4)    !! Time-stepping    : manager of the ocean, tracer and ice time stepping
(    5)    !!======================================================================
(    6) # 293
(  293)    !!----------------------------------------------------------------------
(  294)    !!   Default key                                     NO 1D Config
(  295)    !!----------------------------------------------------------------------
(  296) CONTAINS
(  297)    SUBROUTINE stp_1d ( kt )
(  298) !      WRITE(*,*) 'stp_1d: You should not have seen this print! error?', kt
(  299)    END SUBROUTINE stp_1d
(  300) # 301
(  301)    !!======================================================================
(  302) END MODULE step1d
