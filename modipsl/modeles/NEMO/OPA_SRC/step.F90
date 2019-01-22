MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   stp            : OPA system time-stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE ldftra_oce
   USE ldfdyn_oce
   USE cpl_oce         ! coupled ocean-atmosphere variables
   USE in_out_manager  ! I/O manager
   USE lbclnk

   USE daymod          ! calendar                         (day     routine)

   USE dtatem          ! ocean temperature data           (dta_tem routine)
   USE dtasal          ! ocean salinity    data           (dta_sal routine)
   USE dtasst          ! ocean sea surface temperature    (dta_sst routine)
   USE dtasss          ! ocean sea surface salinity       (dta_sss routine)
   USE taumod          ! surface stress                   (tau     routine)
   USE flxmod          ! thermohaline fluxes              (flx     routine)
   USE ocesbc          ! thermohaline fluxes              (oce_sbc routine)
   USE flxrnf          ! runoffs                          (flx_rnf routine)
   USE flxfwb          ! freshwater budget correction     (flx_fwb routine)
   USE closea          ! closed sea freshwater budget     (flx_clo routine)
   USE ocfzpt          ! surface ocean freezing point    (oc_fz_pt routine)

   USE trcstp          ! passive tracer time-stepping      (trc_stp routine)

   USE dynhpg          ! hydrostatic pressure grad.       (dyn_hpg routine)
   USE dynhpg_atsk     ! hydrostatic pressure grad.  (dyn_hpg_atsk routine)
   USE dynspg_oce      ! surface pressure gradient        (dyn_spg routine)
   USE dynspg          ! surface pressure gradient        (dyn_spg routine)
   USE dynkeg          ! kinetic energy gradient          (dyn_keg routine)
   USE dynvor          ! vorticity term              (dyn_vor_... routines)
   USE dynzad          ! vertical advection               (dyn_adv routine)
   USE dynldf_bilapg   ! lateral mixing            (dyn_ldf_bilapg routine)
   USE dynldf_bilap    ! lateral mixing             (dyn_ldf_bilap routine)
   USE dynldf_iso      ! lateral mixing               (dyn_ldf_iso routine)
   USE dynldf_lap      ! lateral mixing               (dyn_ldf_lap routine)
   USE dynzdf_imp      ! vertical diffusion: implicit     (dyn_zdf routine)
   USE dynzdf_imp_atsk ! vertical diffusion: implicit     (dyn_zdf routine)
   USE dynzdf_iso      ! vertical diffusion: isopycnal    (dyn_zdf routine)
   USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp routine)
   USE dynnxt          ! time-stepping                    (dyn_nxt routine)

   USE trabbc          ! bottom boundary condition        (tra_bbc routine)
   USE trabbl          ! bottom boundary layer            (tra_bbl routine)
   USE tradmp          ! internal damping                 (tra_dmp routine)
   USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
   USE traldf_bilap    ! lateral mixing             (tra_ldf_bilap routine)
   USE traldf_iso      ! lateral mixing               (tra_ldf_iso routine)
   USE traldf_iso_zps  ! lateral mixing           (tra_ldf_iso_zps routine)
   USE traldf_lap      ! lateral mixing               (tra_ldf_lap routine)
   USE traqsr          ! solar radiation penetration      (tra_qsr routine)
   USE tranpc          ! non-penetrative convection       (tra_npc routine)
   USE tranxt          ! time-stepping                    (tra_nxt routine)
   USE traadv_ctl      ! advection scheme control     (tra_adv_ctl routine)
   USE traadv_cen2     ! 2nd order centered scheme   (tra_adv_cen2 routine)
   USE traadv_tvd      ! TVD scheme                (tra_adv_tvd    routine)
   USE traadv_muscl    ! MUSCL scheme              (tra_adv_muscl  routine)
   USE traadv_muscl2   ! MUSCL2 scheme             (tra_adv_muscl2 routine)
!   USE cla             ! cross land advection             (tra_cla routine)
   USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp routine)
   USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp routine)
   USE trazdf_iso      ! vertical diffusion           (tra_zdf_exp routine)
   USE trazdf_iso_vopt ! vertical diffusion           (tra_zdf_exp routine)
   USE trasbc          ! surface boundary condition       (tra_sbc routine)

   USE eosbn2          ! equation of state                (eos_bn2 routine)

   USE obc_par         ! open boundary condition variables
   USE obcdta          ! open boundary condition data     (obc_dta routine)
   USE obcrst          ! open boundary cond. restart      (obc_rst routine)
   USE obcrad          ! open boundary cond. radiation    (obc_rad routine)
   USE obcspg          ! open boundary cond  spg          (obc_spg routine)

   USE divcur          ! hor. divergence and curl      (div & cur routines)
!   USE cla_div         ! cross land: hor. divergence      (div_cla routine)
   USE wzvmod          ! vertical velocity                (wzv     routine)

   USE ldfslp          ! iso-neutral slopes               (ldf_slp routine)
   USE ldfeiv          ! eddy induced velocity coef.      (ldf_eiv routine)
   USE ldfdyn          ! 
   USE ldftra          


   USE zdfbfr          ! bottom friction                  (zdf_bfr routine)
   USE zdftke          ! TKE vertical mixing              (zdf_tke routine)
   USE zdfkpp          ! KPP vertical mixing              (zdf_kpp routine)
   USE zdfddm          ! double diffusion mixing          (zdf_ddm routine)
   USE zdfevd          ! enhanced vertical diffusion      (zdf_evd routine)
   USE zdfric          ! Richardson vertical mixing       (zdf_ric routine)
   USE zdfmxl          ! Mixed-layer depth                (zdf_mxl routine)

   USE zpshde          ! partial step: hor. derivative     (zps_hde routine)
   USE ice_oce         ! sea-ice variable
   USE icestp          ! sea-ice time-stepping             (ice_stp routine)

   USE diawri          ! Standard run outputs             (dia_wri routine)
   USE trdicp          ! Ocean momentum/tracers trends    (trd_wri routine)
   USE trdmld          ! mixed-layer trends               (trd_mld routine)
   USE trdvor          ! vorticity budget                 (trd_vor routine)
!   USE diagap          ! hor. mean model-data gap         (dia_gap routine)
   USE diahdy          ! dynamic height                   (dia_hdy routine)
   USE diaptr          ! poleward transports              (dia_ptr routine)
   USE diahth          ! thermocline depth                (dia_hth routine)
   USE diafwb          ! freshwater budget                (dia_fwb routine)
   USE diaspr          ! suface pressure (rigid-lid)      (dia_spr routine)
   USE flo_oce         ! floats variables
   USE floats          ! floats computation               (flo_stp routine)

   USE stpctl          ! time stepping control            (stp_ctl routine)
   USE restart         ! ocean restart                    (rst_wri routine)
   USE cpl             ! exchanges in coupled mode        (cpl_stp routine)
   USE prtctl          ! Print control                    (prt_ctl routine)

!!DB 2009.08.06
   USE sopa_mc


#if defined key_agrif
   USE agrif_opa_sponge ! Momemtum and tracers sponges
#endif
#ifdef key_RIVER_INPUT
!!DB
   USE rivers
#endif



   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC stp            ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/step.F90,v 1.25 2006/04/26 09:26:48 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp( &
#if !defined key_agrif
   kstp &
#endif   
   )      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!                      
      !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
      !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
      !!              - Tme stepping  of TRC (passive tracer eqs.)
      !! 
      !! ** Method  : -1- Update forcings and data  
      !!              -2- Update ocean physics 
      !!              -3- Compute the t and s trends 
      !!              -4- Update t and s 
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  add a first output record
      !!        !  96-04  (G. Madec)  introduction of dynspg
      !!        !  96-04  (M.A. Foujols)  introduction of passive tracer
      !!   8.0  !  97-06  (G. Madec)  new architecture of call
      !!   8.2  !  97-06  (G. Madec, M. Imbard, G. Roullet)  free surface
      !!   8.2  !  99-02  (G. Madec, N. Grima)  hpg implicit
      !!   8.2  !  00-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
      !!   9.0  !  02-06  (G. Madec)  free form, suppress macro-tasking
      !!    "   !  04-08  (C. Talandier) New trends organization
      !!    "   !  05-01  (C. Ethe) Add the KPP closure scheme
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER &
#if !defined key_agrif   
      , INTENT( in ) &
#endif      
      ::   kstp   ! ocean time-step index

      !! * local declarations
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------
      INTEGER :: ji, jj, jk 



#if defined key_agrif
      kstp = nit000 + Agrif_Nb_Step()
#endif   
      indic = 1                    ! reset to no error condition
      adatrj = adatrj + rdt/86400._wp

!!DB 2008.04.16 -- hardwired ramp function currently used in obcdta, tau_forced_*
     ramp=tanh(kstp*rdt/(2.0*86400.0))

!!DB: perpetual forcing on ====> perpetual_forcing /= 0 ====> do NOT call day()
      if(perpetual_forcing == 0) then
         CALL day( kstp )             ! Calendar
      endif

!!DBG: 2009.06.12
!      call OB_LIMITER(kstp,1)



      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries and Forcings
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( lk_dtatem  )   CALL dta_tem( kstp )         ! update 3D temperature data

      IF( lk_dtasal  )   CALL dta_sal( kstp )         ! Salinity data

      IF( lk_dtasst  )   CALL dta_sst( kstp )         ! Sea Surface Temperature data

      IF( lk_dtasss  )   CALL dta_sss( kstp )         ! Sea Surface salinity data

      IF( lk_obc     )   CALL obc_dta( kstp )         ! update dynamic and tracer data at open boundaries

      IF( lk_obc     )   CALL obc_rad( kstp )         ! compute phase velocities at open boundaries
      
      CALL tau( kstp )             ! wind stress
      
      CALL flx_rnf( kstp )         ! runoff data
      
      CALL flx( kstp )             ! heat and freshwater fluxes

      IF( lk_ice_lim )   CALL ice_stp( kstp )         ! sea-ice model (Update stress & fluxes)
      
      CALL oce_sbc( kstp )         ! ocean surface boundaries
      
      IF( ln_fwb     )   CALL flx_fwb( kstp )         ! freshwater budget
      
      IF( nclosea == 1 ) CALL flx_clo( kstp )         ! closed sea in the domain (update freshwater fluxes)
      
      IF( kstp == nit000 ) THEN 
         IF( ninist == 1 ) THEN                       ! Output the initial state and forcings
            CALL dia_wri_state( 'output.init' )
         ENDIF
      ENDIF
      
      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab2d_1=emp    , clinfo1=' emp  -   : ', mask1=tmask, ovlap=1)
         CALL prt_ctl(tab2d_1=emps   , clinfo1=' emps -   : ', mask1=tmask, ovlap=1)
         CALL prt_ctl(tab2d_1=qt     , clinfo1=' qt   -   : ', mask1=tmask, ovlap=1)
         CALL prt_ctl(tab2d_1=qsr    , clinfo1=' qsr  -   : ', mask1=tmask, ovlap=1)
         CALL prt_ctl(tab2d_1=runoff , clinfo1=' runoff   : ', mask1=tmask, ovlap=1)
         CALL prt_ctl(tab3d_1=tmask  , clinfo1=' tmask    : ', mask1=tmask, ovlap=1, kdim=jpk)
         CALL prt_ctl(tab3d_1=tn     , clinfo1=' sst  -   : ', mask1=tmask, ovlap=1, kdim=1)
         CALL prt_ctl(tab3d_1=sn     , clinfo1=' sss  -   : ', mask1=tmask, ovlap=1, kdim=1)
         CALL prt_ctl(tab2d_1=taux   , clinfo1=' tau  - x : ', tab2d_2=tauy, clinfo2='      - y : ', ovlap=1)
      ENDIF


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Ocean physics update
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !-----------------------------------------------------------------------
      !  VERTICAL PHYSICS
      !-----------------------------------------------------------------------
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------
      
      CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency
      
      !                                                     ! Vertical eddy viscosity and diffusivity coefficients
      IF( lk_zdfric )   CALL zdf_ric( kstp )                       ! Richardson number dependent Kz
      IF( lk_zdftke )   CALL zdf_tke( kstp )                       ! TKE closure scheme for Kz
      IF( lk_zdfkpp )   CALL zdf_kpp( kstp )                       ! KPP closure scheme for Kz
      IF( lk_zdfcst )   avt (:,:,:) = avt0 * tmask(:,:,:)          ! Constant Kz (reset avt to the background value)

      IF( ln_zdfevd )   CALL zdf_evd( kstp )                 ! enhanced vertical eddy diffusivity

      IF( lk_zdfddm .AND. .NOT. lk_zdfkpp)   &
           &              CALL zdf_ddm( kstp )                 ! double diffusive mixing

      CALL zdf_bfr( kstp )                 ! bottom friction
      
      CALL zdf_mxl( kstp )                 ! mixed layer depth


      !-----------------------------------------------------------------------
      !  LATERAL PHYSICS
      !-----------------------------------------------------------------------
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------

      IF( lk_ldfslp     )   CALL ldf_slp( kstp, rhd, rn2 )       ! before slope of the lateral mixing

#if defined key_traldf_c2d
      IF( lk_traldf_eiv )   CALL ldf_eiv( kstp )                 ! eddy induced velocity coefficient
#endif


#if defined key_passivetrc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------

      CALL trc_stp( kstp, indic )            ! time-stepping

#endif


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Active tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va arrays are used as workspace in this section
      !-----------------------------------------------------------------------
      
      ta(:,:,:) = 0.e0               ! set tracer trends to zero
      sa(:,:,:) = 0.e0
      
      CALL tra_sbc( kstp )           ! surface boundary condition

      IF( ln_traqsr        )   CALL tra_qsr( kstp )           ! penetrative solar radiation qsr

      IF( lk_trabbc        )   CALL tra_bbc( kstp )           ! bottom heat flux
      
      IF( lk_trabbl_dif    )   CALL tra_bbl_dif( kstp )           ! diffusive bottom boundary layer scheme
      IF( lk_trabbl_adv    )   CALL tra_bbl_adv( kstp )           ! advective (and/or diffusive) bottom boundary layer scheme
      
      IF( lk_tradmp        )   CALL tra_dmp( kstp )           ! internal damping trends

      !                                                       ! horizontal & vertical advection
      IF( kstp == nit000   )   CALL tra_adv_ctl                    ! chose/control the scheme used
#ifdef key_RIVER_INPUT
      call riv_tra(kstp)  !JC: Put T/S at river points
      call riv_dyn(kstp)  !JC:  assign un/vn at River points
#endif

!-------------------------------------------------------------------------------------------
      IF( ln_traadv_cen2   )   CALL tra_adv_cen2  ( kstp )         ! 2nd order centered scheme
      IF( ln_traadv_tvd    )   CALL tra_adv_tvd   ( kstp )         ! TVD scheme
      IF( ln_traadv_muscl  )   CALL tra_adv_muscl ( kstp )         ! MUSCL scheme
      IF( ln_traadv_muscl2 )   CALL tra_adv_muscl2( kstp )         ! MUSCL2 scheme

!!DB: orca-related
!      IF( n_cla == 1       )   CALL tra_cla( kstp )           ! Cross Land Advection (Update Hor. advection)

      !                                                       ! lateral mixing 
      IF( l_traldf_lap     )   CALL tra_ldf_lap    ( kstp )           ! iso-level laplacian
      IF( l_traldf_bilap   )   CALL tra_ldf_bilap  ( kstp )           ! iso-level bilaplacian 
      IF( l_traldf_bilapg  )   CALL tra_ldf_bilapg ( kstp )           ! s-coord. horizontal bilaplacian
      IF( l_traldf_iso     )   CALL tra_ldf_iso    ( kstp )           ! iso-neutral/geopot. laplacian 
      IF( l_traldf_iso_zps )   CALL tra_ldf_iso_zps( kstp )           ! partial step iso-neutral/geopot. laplacian

#if defined key_agrif
      IF (.NOT. Agrif_Root())  CALL Agrif_Sponge_tra( kstp )          ! tracers sponge
#endif
      !                                                       ! vertical diffusion
      IF( l_trazdf_exp     )   CALL tra_zdf_exp     ( kstp )          ! explicit time stepping (time splitting scheme)
      IF( l_trazdf_imp     )   CALL tra_zdf_imp     ( kstp )          ! implicit time stepping (euler backward)
      IF( l_trazdf_iso     )   CALL tra_zdf_iso     ( kstp )          ! isopycnal
      IF( l_trazdf_iso_vo  )   CALL tra_zdf_iso_vopt( kstp )          ! vector opt. isopycnal
      
      CALL tra_nxt( kstp )           ! tracer fields at next time step
      
      IF( ln_zdfnpc        )   CALL tra_npc( kstp )           ! update the new (t,s) fields by non
      !                                                       ! penetrative convective adjustment
      
      IF( ln_dynhpg_imp    ) THEN                             ! semi-implicit hpg 
         CALL eos( ta, sa, rhd, rhop )   ! Time-filtered in situ density used in dynhpg module
         IF( lk_zps    )          CALL zps_hde( kstp, ta, sa, rhd,  & ! Partial steps: time filtered hor. gradient 
              &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
              &                                        gtv, gsv, grv )  
      ELSE                                                    ! centered hpg (default case)
         CALL eos( tb, sb, rhd, rhop )       ! now (swap=before) in situ density for dynhpg module
         IF( lk_zps    )          CALL zps_hde( kstp, tb, sb, rhd,  & ! Partial steps: now horizontal gradient
              &                                        gtu, gsu, gru, & ! of t, s, rd at the bottom ocean level
              &                                        gtv, gsv, grv )  
      ENDIF
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Dynamics
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ta, sa arrays are used as workspace in this section 
      !-----------------------------------------------------------------------

      ua(:,:,:) = 0.e0               ! set dynamics trends to zero
      va(:,:,:) = 0.e0

      CALL dyn_keg( kstp )           ! horizontal gradient of kinetic energy

      !                                                       ! vorticity term including Coriolis
      IF( kstp == nit000   )   CALL dyn_vor_ctl                      ! chose/control the scheme used
      IF( ln_dynvor_ens    )   CALL dyn_vor_enstrophy( kstp )        ! enstrophy conserving scheme
      IF( ln_dynvor_ene    )   CALL dyn_vor_energy   ( kstp )        ! energy conserving scheme
      IF( ln_dynvor_mix    )   CALL dyn_vor_mixed    ( kstp )        ! mixed energy/enstrophy conserving scheme
      IF( ln_dynvor_een    )   CALL dyn_vor_ene_ens  ( kstp )        ! combined energy/enstrophy conserving scheme
      
      !                                                       ! lateral mixing 
      IF( l_dynldf_lap     )   CALL dyn_ldf_lap    ( kstp )          ! iso-level laplacian
      IF( l_dynldf_bilap   )   CALL dyn_ldf_bilap  ( kstp )          ! iso-level bilaplacian 
      IF( l_dynldf_bilapg  )   CALL dyn_ldf_bilapg ( kstp )          ! s-coord. horizontal bilaplacian
      IF( l_dynldf_iso     )   CALL dyn_ldf_iso    ( kstp )          ! iso-neutral laplacian 
      
#if defined key_agrif
      IF (.NOT. Agrif_Root())  CALL Agrif_Sponge_dyn( kstp )         ! momemtum sponge
#endif
      !                                                       ! horizontal gradient of Hydrostatic pressure 
      IF ( lk_jki ) THEN
         CALL dyn_hpg_atsk( kstp )             ! autotask case (j-k-i loop)
      ELSE
         CALL dyn_hpg     ( kstp )             ! default case  (k-j-i loop)
      ENDIF
      
      CALL dyn_zad    ( kstp )       ! vertical advection       
      
      !                                                       ! vertical diffusion
      IF( l_dynzdf_exp     )   CALL dyn_zdf_exp    ( kstp )          ! explicit time stepping (time splitting scheme)
      IF( l_dynzdf_imp     )   CALL dyn_zdf_imp    ( kstp )          ! implicit time stepping (euler backward)
      IF( l_dynzdf_imp_tsk )   CALL dyn_zdf_imp_tsk( kstp )          ! autotask implicit time stepping (euler backward)
      IF( l_dynzdf_iso     )   CALL dyn_zdf_iso    ( kstp )          ! iso-neutral case
      
      IF( lk_dynspg_rl ) THEN 
         IF( lk_obc    )       CALL obc_spg( kstp )           ! surface pressure gradient at open boundaries
      ENDIF
      indic=0
      !i bug lbc sur emp
      CALL lbc_lnk( emp, 'T', 1. )
      !i
      CALL dyn_spg( kstp, indic )    ! surface pressure gradient
      
      CALL dyn_nxt( kstp )           ! velocity at next time step 

      
      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Computation of diagnostic variables
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------
      
      CALL oc_fz_pt                        ! ocean surface freezing temperature
      
      CALL div_cur( kstp )                 ! Horizontal divergence & Relative vorticity

!!DB: orca-related      
!      IF( n_cla == 1 ) CALL div_cla( kstp )                 ! Cross Land Advection (Update Hor. divergence)
      
      CALL wzv( kstp )                     ! Vertical velocity
      
!!DB 2008.04.07 -- the below ultimately calls the appropriate routine to update the
!! relevant diffusivity coeffs (see ldfdyn.F90). The prob is that it performs numerous
!! calcs every timestep (including reading namelist) that are not necessary
!! =========> modification ... DONE ...

!      !ZW
!      CALL ldf_dyn_init
!      CALL ldf_tra_init
!      !ZW

!!DB 04.09 -- modify above so that init routines are not called every dt
!!Note that the smag keys are separated to the end of this code fragment
!!Also, there is a potential problem if there is a mistake in the keys
!!so that (e.g.) no or more-than-one tracer key is defined. I do not check
!!for this possibility
!!Momentum diffusivity updates
#if defined key_dynldf_c3d
      CALL ldf_dyn_c3d   ! ahm = 3D coef. = F( longitude, latitude, depth )
#elif defined key_dynldf_c2d
      CALL ldf_dyn_c2d   ! ahm = 1D coef. = F( longitude, latitude )
#elif defined key_dynldf_c1d
      CALL ldf_dyn_c1d   ! ahm = 1D coef. = F( depth )
#else
      !do nothing unless smag is on -- see below
#endif
!!tracer diffusivity updates
#if defined key_traldf_c3d
      CALL ldf_tra_c3d           ! aht = 3D coef. = F( longitude, latitude, depth )
#elif defined key_traldf_c2d
      CALL ldf_tra_c2d           ! aht = 2D coef. = F( longitude, latitude )
#elif defined key_traldf_c1d
      CALL ldf_tra_c1d           ! aht = 1D coef. = F( depth )
#else
      !do nothing unless smag is on -- see below
#endif

#if defined key_dynldf_smag || defined key_traldf_smag
      call ldf_smag( kstp )
#endif

      
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control, diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------
      
      !                                            ! Time loop: control and print
      CALL stp_ctl( kstp, indic )
      IF ( indic < 0 )   nstop = nstop + 1


      IF ( nstop == 0 ) THEN
         !                                         ! Diagnostics:
         IF( lk_floats  )   CALL flo_stp( kstp )                 ! drifting Floats
         IF( lk_trddyn  )   CALL trd_dwr( kstp )                 ! trends: dynamics 
         IF( lk_trdtra  )   CALL trd_twr( kstp )                 ! trends: active tracers
         IF( lk_trdmld  )   CALL trd_mld( kstp )                 ! trends: Mixed-layer 
         IF( lk_trdvor  )   CALL trd_vor( kstp )                 ! trends: vorticity budget
         IF( lk_diaspr  )   CALL dia_spr( kstp )                 ! Surface pressure diagnostics
         IF( lk_diahth  )   CALL dia_hth( kstp )                 ! Thermocline depth (20 degres isotherm depth)
!         IF( lk_diagap  )   CALL dia_gap( kstp )                 ! basin averaged diagnostics
         IF( lk_diahdy  )   CALL dia_hdy( kstp )                 ! dynamical heigh diagnostics
         IF( lk_diafwb  )   CALL dia_fwb( kstp )                 ! Fresh water budget diagnostics
         IF( ln_diaptr  )   CALL dia_ptr( kstp )                 ! Poleward TRansports diagnostics
         
         !                                         ! save and outputs
         CALL rst_write  ( kstp )             ! ocean model: restart file output
!!DB
!         IF( lk_obc     )   CALL obc_rst_wri( kstp )             ! ocean model: open boundary restart file output
         
!!DB 2009.08.06 -- These are the *_grid_*.nc files which we never look at
!! (Replaced by M2 aves for U,V and *_aveTSUV.nc (see below) for other variables
!         CALL dia_wri    ( kstp, indic )      ! ocean model: outputs

!!DB
         if(M2_ave > 0) call output_special(kstp, indic)      ! special M2 time-averaged fields
!!DB
         if(ioutput_ave /= 0) call output_aveTSUV(kstp, indic)      ! special time-averaged fields

         
      ENDIF

!DB Output some BoF vels -- NB: code moved to sopa_mc module
!! Model diagnostics (somewhat) specific to sopa MC domain 
      !call sopa_mc_diagnostics(kstp)




      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Coupled mode
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( lk_cpl    )   CALL cpl_stp( kstp )                 ! coupled mode : field exchanges

   END SUBROUTINE stp

   !!======================================================================


END MODULE step
