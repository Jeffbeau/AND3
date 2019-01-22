MODULE trctrp
   !!======================================================================
   !!                       ***  MODULE trctrp  ***
   !! Ocean Physics    : manage the passive tracer transport
   !!======================================================================

!!DB: Do not think I need this anymore so hack-it-out with a bogus key
#if defined key_OLD_MODULE

!#if defined key_passivetrc
!#if defined key_passivetrc && defined key_BGCM_01
   !!----------------------------------------------------------------------
   !!   trc_trp        : passive tracer transport
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and active tracers variables
   USE trc             ! ocean passive tracers variables 

!   USE trctrp_lec      ! passive tracers transport

   USE trcbbl          ! bottom boundary layer               (trc_bbl routine)
   USE trcdmp          ! internal damping                    (trc_dmp routine)

   USE trcldf_bilapg   ! lateral mixing               (trc_ldf_bilapg routine)
   USE trcldf_bilap    ! lateral mixing                (trc_ldf_bilap routine)
   USE trcldf_iso      ! lateral mixing                  (trc_ldf_iso routine)
   USE trcldf_iso_zps  ! lateral mixing              (trc_ldf_iso_zps routine)
   USE trcldf_lap      ! lateral mixing                  (trc_ldf_lap routine)
 
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcrad          ! positivity                          (trc_rad routine)

   USE trcadv_cen2     ! 2nd order centered advection   (trc_adv_cen2 routine)
   USE trcadv_muscl    ! MUSCL advection               (trc_adv_muscl routine)
   USE trcadv_muscl2   ! MUSCL2 advection             (trc_adv_muscl2 routine)
   USE trcadv_tvd      ! TVD advection                   (trc_adv_tvd routine)
   USE trcadv_smolar   ! SMOLAR advection             (trc_adv_smolar routine)

   USE trczdf_exp      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_imp      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_iso      ! vertical diffusion              (trc_zdf_exp routine)
   USE trczdf_iso_vopt ! vertical diffusion              (trc_zdf_exp routine)
   USE trcsbc          ! surface boundary condition          (trc_sbc routine)

   USE zpshde_trc      ! partial step: hor. derivative   (zps_hde_trc routine)

#ifdef key_BGCM_01
   USE lib_bgcm_01
   USE in_out_manager
#endif



   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_trp            ! called by trc_stp

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp.F90,v 1.10 2006/04/11 13:49:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE trc_trp  ***
      !!                      
      !! ** Purpose : Management of passive tracers transport
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe)  Original
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      !! ---------------------------------------------------------------------

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !-----------------------------------------------------------------------

      CALL trc_sbc( kt )            ! surface boundary condition

# if defined key_trcbbc
      IF(lwp) WRITE(numout,cform_err)
      IF(lwp) WRITE(numout,*) ' Bottom heat flux not yet implemented'
      IF(lwp) WRITE(numout,*) ' With passive tracers. '
      IF(lwp) WRITE(numout,*) ' Check trc_trp routine'
      nstop = nstop + 1
# endif 
      !                                                      ! bottom boundary condition
      IF( lk_trcbbl_dif    )   CALL trc_bbl_dif( kt )                ! diffusive bottom boundary layer scheme
      IF( lk_trcbbl_adv    )   CALL trc_bbl_adv( kt )                ! advective (and/or diffusive) bottom boundary layer scheme

      IF( lk_trcdmp        )   CALL trc_dmp( kt )            ! internal damping trends

!!DB 2008.11.19 -- Enclose changes in ifdef just in case
#ifdef key_BGCM_01
!!DBG: 10.03 -- force this call as for some reason the tvd scheme --> nan
!!All of the below work, but cen2 seems to be slightly faster
!!DBG: no advection
!      if(lwp .AND. kt==nit000)write(numout2,*)'DBG -- BGCM: NO ADVECTION'
      CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme -- works

!      CALL trc_adv_muscl ( kt )             ! MUSCL scheme -- works
!      CALL trc_adv_muscl2 ( kt )             ! MUSCL scheme -- works
#else
!!DB: NB same as above but other possibilities are more apparent
      !                                                      ! horizontal & vertical advection
!      IF( ln_trcadv_cen2   )   CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme
!      IF( ln_trcadv_muscl  )   CALL trc_adv_muscl ( kt )             ! MUSCL scheme
!      IF( ln_trcadv_muscl2 )   CALL trc_adv_muscl2( kt )             ! MUSCL2 scheme
!      IF( ln_trcadv_tvd    )   CALL trc_adv_tvd   ( kt )             ! TVD scheme
!      IF( ln_trcadv_smolar )   CALL trc_adv_smolar( kt )             ! SMOLARKIEWICZ scheme

!!DBG: 10.03 -- force this call as for some reason the tvd scheme --> nan
!!All of the below work, but cen2 seems to be slightly faster
      CALL trc_adv_cen2  ( kt )             ! 2nd order centered scheme -- works
!      CALL trc_adv_muscl ( kt )             ! MUSCL scheme -- works
!      CALL trc_adv_muscl2 ( kt )             ! MUSCL scheme -- works
#endif


      IF( n_cla == 1   )   THEN 
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          Cross Land Advection not yet implemented'
         IF(lwp) WRITE(numout,*) '          With Passive tracers. n_cla = ', n_cla
         IF(lwp) WRITE(numout,*) '          Check trc_trp routine'
         nstop = nstop + 1
      ENDIF

!!DBG
!      if(lwp .AND. kt==nit000)write(numout2,*)'DBG -- BGCM: NO HORIZONTAL DIFFUSION'
      !                                                      ! lateral mixing 
      IF( l_trcldf_bilapg  )   CALL trc_ldf_bilapg ( kt )            ! s-coord. horizontal bilaplacian
      IF( l_trcldf_bilap   )   CALL trc_ldf_bilap  ( kt )            ! iso-level bilaplacian 
      IF( l_trcldf_iso     )   CALL trc_ldf_iso    ( kt )            ! iso-neutral laplacian 
      IF( l_trcldf_iso_zps )   CALL trc_ldf_iso_zps( kt )            ! partial step iso-neutral laplacian
      IF( l_trcldf_lap     )   CALL trc_ldf_lap    ( kt )            ! iso-level laplacian

#ifdef key_BGCM_01
!!DB -- force implicit scheme
      call trc_zdf_imp( kt )                ! implicit time stepping (euler backward)
#else
      !                                                      ! vertical diffusion
      IF( l_trczdf_exp     )   CALL trc_zdf_exp( kt )                ! explicit time stepping (time splitting scheme)
      IF( l_trczdf_imp     )   CALL trc_zdf_imp( kt )                ! implicit time stepping (euler backward)
      IF( l_trczdf_iso     )   CALL trc_zdf_iso( kt )                ! isopycnal
      IF( l_trczdf_iso_vo  )   CALL trc_zdf_iso_vopt( kt )           ! vector opt. isopycnal
#endif

      CALL trc_nxt( kt )            ! tracer fields at next time step
 
!!DB: 2009.09.09 routine disabled
!      CALL trc_rad( kt )            ! Correct artificial negative concentrations for isopycnal scheme
      !                                                      

      IF( lk_zps .AND. .NOT. lk_trccfg_1d ) &
         &                     CALL zps_hde_trc( kt, trb, gtru, gtrv )  ! Partial steps: now horizontal gradient
      !                                                                 ! of passive tracers at the bottom ocean level


    END SUBROUTINE trc_trp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_trp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_trp
#endif
   
   !!======================================================================
END MODULE trctrp
