MODULE step1d
   !!======================================================================
   !!                       ***  MODULE step1D  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!======================================================================
#if defined key_cfg_1d
   !!----------------------------------------------------------------------
   !!   'key_cfg_1d'               1D Configuration
   !!----------------------------------------------------------------------  
   !!----------------------------------------------------------------------
   !!   stp_1d           : OPA system time-stepping on 1 direction
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE ldftra_oce
   USE ldfdyn_oce
   USE in_out_manager  ! I/O manager
   USE lbclnk

   USE daymod          ! calendar                         (day     routine)

   USE dtatem          ! ocean temperature data           (dta_tem routine)
   USE dtasal          ! ocean salinity    data           (dta_sal routine)
   USE dtasst          ! ocean sea surface temerature     (dta_sst routine)
   USE taumod          ! surface stress                   (tau     routine)
   USE flxmod          ! thermohaline fluxes              (flx     routine)
   USE ocesbc          ! thermohaline fluxes              (oce_sbc routine)
   USE flxrnf          ! runoffs                          (flx_rnf routine)
   USE flxfwb          ! freshwater budget correction     (flx_fwb routine)
   USE ocfzpt          ! surface ocean freezing point    (oc_fz_pt routine)

   USE trcstp          ! passive tracer time-stepping     (trc_stp routine)

   USE dynzdf_imp      ! vertical diffusion: implicit     (dyn_zdf routine)
   USE dynzdf_imp_atsk ! vertical diffusion: implicit     (dyn_zdf routine)
   USE dynzdf_iso      ! vertical diffusion: isopycnal    (dyn_zdf routine)
   USE dynzdf_exp      ! vertical diffusion: explicit (dyn_zdf_exp routine)
 

   USE traqsr          ! solar radiation penetration      (tra_qsr routine)
   USE tranxt          ! time-stepping                    (tra_nxt routine)
   USE trazdf_exp      ! vertical diffusion: explicit (tra_zdf_exp routine)
   USE trazdf_imp      ! vertical diffusion: implicit (tra_zdf_imp routine)
   USE trazdf_iso      ! vertical diffusion           (tra_zdf_exp routine)
   USE trazdf_iso_vopt ! vertical diffusion           (tra_zdf_exp routine)
   USE trasbc          ! surface boundary condition       (tra_sbc routine)

   USE eosbn2

   USE zdfbfr          ! bottom friction                  (zdf_bfr routine)
   USE zdftke          ! TKE vertical mixing              (zdf_tke routine)
   USE zdfkpp          ! KPP vertical mixing              (zdf_kpp routine)
   USE zdfddm          ! double diffusion mixing          (zdf_ddm routine)
   USE zdfevd          ! enhanced vertical diffusion      (zdf_evd routine)
   USE zdfric          ! Richardson vertical mixing       (zdf_ric routine)
   USE zdfmxl          ! Mixed-layer depth                (zdf_mxl routine)

   USE dyncor1d
   USE dynnxt1d
   USE diawri1d        ! Standard run outputs             (dia_wri_1d routine)

   USE ice_oce         ! sea-ice variable
   USE icestp1d        ! sea-ice time-stepping             (ice_stp routine)

   USE diawri          ! Standard run outputs             (dia_wri_state routine)


   USE stpctl          ! time stepping control            (stp_ctl routine)
   USE restart         ! ocean restart                    (rst_wri routine)
   USE prtctl          ! Print control                    (prt_ctl routine)
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC stp_1d            ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/step1d.F90,v 1.5 2006/04/11 13:52:28 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp_1d( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp1D  ***
      !!                      
      !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
      !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
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
      !!        !  04-10  (C. Ethe) 1D configuration
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kstp   ! ocean time-step index

      !! * local declarations
      INTEGER ::   indic    ! error indicator if < 0
!!      INTEGER ::   ii0, ii1, ij0, ij1   ! temporary integers
      !! ---------------------------------------------------------------------

      indic = 1                    ! reset to no error condition
      adatrj = adatrj + rdt/86400._wp

      CALL day( kstp )             ! Calendar

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Update data, open boundaries and Forcings
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      IF( lk_dtatem  )   CALL dta_tem( kstp )         ! update 3D temperature data

      IF( lk_dtasal  )   CALL dta_sal( kstp )         ! Salinity data

      IF( lk_dtasst  )   CALL dta_sst( kstp )         ! Sea Surface Temperature data

                         CALL tau( kstp )             ! wind stress

                         CALL flx_rnf( kstp )         ! runoff data

                         CALL flx( kstp )             ! heat and freshwater fluxes

      IF( lk_ice_lim )  THEN 
                        CALL ice_stp_1d( kstp )      ! sea-ice model (Update stress & fluxes)
      ELSE
                        CALL oce_sbc( kstp )         ! ocean surface boudaries
      ENDIF

      IF( ln_fwb     )   CALL flx_fwb( kstp )         ! freshwater budget


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
      IF( lk_zdfkpp )   CALL zdf_kpp( kstp )                       ! KPP scheme for Kz
      IF( lk_zdfcst )   avt (:,:,:) = avt0 * tmask(:,:,:)          ! Constant Kz (reset avt to the background value)


      IF( ln_zdfevd )   CALL zdf_evd( kstp )                 ! enhanced vertical eddy diffusivity

      IF( lk_zdfddm .AND. .NOT. lk_zdfkpp)   &
         &              CALL zdf_ddm( kstp )                 ! double diffusive mixing

                        CALL zdf_bfr( kstp )                 ! bottom friction

                        CALL zdf_mxl( kstp )                 ! mixed layer depth


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

      !                                                       ! vertical diffusion
      IF( l_trazdf_exp     )   CALL tra_zdf_exp     ( kstp )          ! explicit time stepping (time splitting scheme)
      IF( l_trazdf_imp     )   CALL tra_zdf_imp     ( kstp )          ! implicit time stepping (euler backward)
      IF( l_trazdf_iso     )   CALL tra_zdf_iso     ( kstp )          ! isopycnal
      IF( l_trazdf_iso_vo  )   CALL tra_zdf_iso_vopt( kstp )          ! vector opt. isopycnal

                               CALL tra_nxt( kstp )           ! tracer fields at next time step

                               CALL eos( tb, sb, rhd, rhop )       ! now (swap=before) in situ density for dynhpg module

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Dynamics
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ta, sa arrays are used as workspace in this section 
      !-----------------------------------------------------------------------

                               ua(:,:,:) = 0.e0               ! set dynamics trends to zero
                               va(:,:,:) = 0.e0
  
                               CALL dyn_cor_1d     ( kstp )
      !                                                       ! vertical diffusion
      IF( l_dynzdf_exp     )   CALL dyn_zdf_exp    ( kstp )         ! explicit time stepping (time splitting scheme)
      IF( l_dynzdf_imp     )   CALL dyn_zdf_imp    ( kstp )         ! implicit time stepping (euler backward)
      IF( l_dynzdf_imp_tsk )   CALL dyn_zdf_imp_tsk( kstp )         ! autotask implicit time stepping (euler backward)
      IF( l_dynzdf_iso     )   CALL dyn_zdf_iso    ( kstp )         ! iso-neutral case

!i bug lbc sur emp
      CALL lbc_lnk( emp, 'T', 1. )
!i
 
                                CALL dyn_nxt_1d( kstp )          ! lateral velocity at next time step 


      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Computation of diagnostic variables
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------

                               CALL oc_fz_pt                    ! ocean surface freezing temperature


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
         !                                         ! save and outputs
                           CALL rst_write  ( kstp )              ! ocean model: restart file output
                           CALL dia_wri_1d ( kstp, indic )       ! ocean model: outputs

      ENDIF


   END SUBROUTINE stp_1d
#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Config
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE stp_1d ( kt )
!      WRITE(*,*) 'stp_1d: You should not have seen this print! error?', kt
   END SUBROUTINE stp_1d
#endif
   !!======================================================================
END MODULE step1d
