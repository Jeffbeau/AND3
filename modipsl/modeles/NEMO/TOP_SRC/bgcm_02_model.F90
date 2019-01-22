!!2010.08
!!Skeleton for the IML model = key_BGCM_02
!!See bgcm_01 programs for other code fragments of interest

!!DB: 
#define is_3D  !! ====> call adv/diff routines

#if defined key_passivetrc && defined key_BGCM_02
MODULE bgcm_02_model
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================


   !!----------------------------------------------------------------------
   !!   trc_stp      : passive tracer system time-stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc              ! ocean passive tracers variables 
!!DBG see ca. line-1
!   USE trctrp           ! passive tracers transport
   USE trcadv_cen2     ! 2nd order centered advection   (trc_adv_cen2 routine)
!   USE trcadv_tvd     ! 2nd order centered advection   (trc_adv_cen2 routine)
   USE trcldf_lap      ! lateral mixing                  (trc_ldf_lap routine)
!   USE trcldf_iso      ! 
   USE trczdf_imp      ! vertical diffusion              (trc_zdf_exp routine)
!   USE trczdf_iso      ! 
   USE trcnxt          ! time-stepping                       (trc_nxt routine)
   USE trcdmp          ! internal damping                 (trc_dmp routine)   ! NL#7

   USE lib_bgcm_02
   USE lifemaker3D
   USE lib_ncdf

!   USE sopa_mc, ONLY : divN

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_stp           ! called by step

!!DB
   PUBLIC bgcm_sms


 CONTAINS

   SUBROUTINE trc_stp( kt, kindic )
      !!-------------------------------------------------------------------
      !!                     ***  ROUTINE trc_stp  ***
      !!                      
      !! ** Purpose : Time loop of opa for passive tracer
      !! 
      !! ** Method  : 
      !!              Compute the passive tracers trends 
      !!              Update the passive tracers
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe)  Original
      !!-------------------------------------------------------------------

      USE lbclnk

      !! * Arguments
      INTEGER, INTENT( in ) ::  kt  ! ocean time-step index
      INTEGER, INTENT( in ) ::  kindic
      CHARACTER (len=25) :: charout

!!DB
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zfield, transx_no3,transy_no3,transz_no3

      integer :: ji,jj,jk,jn, rec_num, status
      CHARACTER (len=80) :: fname
      REAL(wp) :: tmp1
!      LOGICAL :: is_1D = .false. 

!!DB write into to file (HARDWIRED):
      if(lwp .AND. kt==nit000) then
         write(numout2,*)'----------------------'
         write(numout2,*)'BGCM 02 MODULE '
         write(numout2,*)jptra,' Tracers'
         write(numout2,*)'----------------------'
      endif

   !   if(kt==nit000) then
   !      transx_no3 = 0 ;transy_no3 = 0 ;transz_no3 = 0
   !   endif 

      ! this ROUTINE is called only every ndttrc time step
      IF( MOD( kt , ndttrc ) /= 0 ) RETURN

!!This is the routine that calls the BGCM Source-Minus-Sink code
      CALL bgcm_sms( kt )   !lifemaker
!      if (nproc.eq.0) print*, 'after life trn= ', trn(3,3,1,1:8)
!      if (nproc.eq.0) print*, 'after life trb= ', trb(3,3,1,1:8)

      IF( lk_trcdmp        )   CALL trc_dmp( kt )            ! internal damping trends  ! NL#7

!      IF( ln_trcadv_cen2   )   CALL trc_adv_cen2  ( kt )         ! 2nd order centered scheme
!      IF( ln_trcadv_tvd    )   CALL trc_adv_tvd   ( kt )         ! TVD scheme
!      IF( ln_trcadv_muscl  )   CALL trc_adv_muscl ( kt )         ! MUSCL scheme
!      IF( l_trcldf_lap     ) CALL trc_ldf_lap    ( kt )           ! iso-level laplacian
!      IF( l_trcldf_bilap   )   CALL trc_ldf_bilap  ( kt )           ! iso-level bilaplacian
!      IF( l_trcldf_bilapg  )   CALL trc_ldf_bilapg ( kt )           ! s-coord. horizontal bilaplacian
!      IF( l_trcldf_iso     )
!      CALL trc_ldf_iso    ( kt )           ! iso-neutral/geopot. laplacian
!      IF( l_trczdf_exp     )   CALL trc_zdf_exp     ( kt )          ! explicit time stepping (time splitting scheme)
!      IF( l_trczdf_imp     )   CALL trc_zdf_imp     ( kt )          ! implicit time stepping (euler backward)
!      IF( l_trczdf_iso     )
!      CALL trc_zdf_iso     ( kt )          ! isopycnal
!      IF( l_trczdf_iso_vo  )   CALL trc_zdf_iso_vopt( kt )          ! vector opt. isopycnal
!------------------------------------

#if defined is_3D
      CALL trc_adv_cen2 ( kt )            ! 2nd order centered scheme --fastest
!      CALL trc_adv_tvd ( kt )            ! 2nd order centered scheme --fastest
      CALL trc_ldf_lap  ( kt )            ! iso-level laplacian
#endif
      CALL trc_zdf_imp  ( kt )            ! implicit time stepping (euler backward)
      CALL trc_nxt      ( kt )            ! tracer fields at next time step

!      if (nproc.eq.0) print*, 'after tnxt trn= ', trn(3,3,1,1:8)
!      if (nproc.eq.0) print*, 'after tnxt trb= ', trb(3,3,1,1:8)
!!DB: Open Boundary code 
!!This is where special OBC-related operations could also be done
#ifdef key_obc

#endif
!print*, kt, nrsttr, nitend
!print*, mod(kt, nrsttr)
!print*, 'before call to restart'

!!DB: write restart -- not yet written
!!Restart file will be v.similar to regular output file 
!!the last frame of which can be used for the time-being. 
!      CALL trc_wri( kt )            ! outputs
      if( mod( kt, nrsttr ) == 0 .OR. kt == nitend ) THEN     !NL#2
        print*, 'In call to restart'
        call rst_bgcm_write(kt)     !NL#2
      endif     !NL#2

!print*, 'after call to restart'
      ! NL#3 : sum the flux of the nitrate for this step (mmol m-2)
!      transx_no3= transx_no3 + trn(:,:,:,3)*un*rdt
!      transy_no3= transy_no3 + trn(:,:,:,3)*vn*rdt
!      transz_no3= transz_no3 + trn(:,:,:,3)*wn*rdt

!!DB: output model variables
      if(mod( (kt-nit000+1), nwritetrc ) == 0 ) then
!!DB: write to ncdf file
         rec_num = (kt-nit000+1)/nwritetrc

         if(lwp) write(numout2,*)'PASSIVE TRACER: writing rec_num ',rec_num
         CALL ncdf_write(BGCM_fname, 'time_counter', REAL(kt * rdt), rec_num, status)

!!DB: 2008.10.17: Both of the below work 
!!DB: Note that to force a write of a 2D,3D,4D variable directly to a record number
!!    use a -ve value for rec_num

         CALL ncdf_write(BGCM_fname, 'trn', trn, jptra, -rec_num, status)
!         CALL ncdf_write(BGCM_fname, 'transx_no3', transx_no3, -rec_num, status)  ! NL#3
!         CALL ncdf_write(BGCM_fname, 'transy_no3', transy_no3, -rec_num, status)  ! NL#3
!         CALL ncdf_write(BGCM_fname, 'transz_no3', transz_no3, -rec_num, status)  ! NL#3
!         transx_no3 = 0 ;transy_no3 = 0 ;transz_no3 = 0      ! NL#3  reset because it accumulate for the 180 time step

!         CALL ncdf_write(BGCM_fname, 'parz', parz, -rec_num, status)
!         CALL ncdf_write(BGCM_fname, 'par_surf', (PAR_surf/24)*60.386473, -rec_num, status)
         CALL ncdf_write(BGCM_fname, 'par_surf', PAR_surf_2h, -rec_num, status)
!         CALL ncdf_write(BGCM_fname, 'kcdom', kcdom, -rec_num, status)
!         CALL ncdf_write(BGCM_fname, 'ze01',  ze01 , -rec_num, status)  ! NL#4

#if defined (key_carbon)
         CALL ncdf_write(BGCM_fname, 'pH', pHW3D, -rec_num, status)
         CALL ncdf_write(BGCM_fname, 'satca', SatCa3D, -rec_num, status)
         CALL ncdf_write(BGCM_fname, 'satar', SatAr3D, -rec_num, status)
         CALL ncdf_write(BGCM_fname, 'flxCO2',  flx_CO2 , -rec_num, status)
         flx_CO2 = 0.  !reset 
#endif
         CALL ncdf_write(BGCM_fname, 'ndastp',REAL(ndastp), rec_num, status)
         CALL ncdf_write(BGCM_fname, 'model_time_step',REAL(kt), rec_num, status)
         CALL ncdf_write(BGCM_fname, 'model_time',model_time, rec_num, status)
      endif

!!DB: tra is done in trc_trp ---> trc_nxt. 
!!To be safe do all of them (as done in p4zprg)
!!NB: Likely do not have to do this
      DO jn=1 , jptra
        CALL lbc_lnk(trn(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(trb(:,:,:,jn), 'T', 1. )
        CALL lbc_lnk(tra(:,:,:,jn), 'T', 1. )
      END DO

   END SUBROUTINE trc_stp

!!DB: BGCM_02 SMS model
!!    Placeholder for future use
   SUBROUTINE BGCM_sms( kt )
      !!===========================================================================================

      !! * Arguments
      !! -----------
     INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      

!!DB: The below must go here, i.e. after last local declaration and before first assignment statement
!!    The file contains the assignment of variable names as pointers to trn() locations
!AD: defined inside lifemaker3D module#   include "bgcm_02_pointers.h90"

     !! this ROUTINE is called only every ndttrc time step
     !! --------------------------------------------------
!     IF ( MOD(kt,ndttrc) /= 0) RETURN
     
      call lifemaker(kt)
     
   END SUBROUTINE BGCM_sms


END MODULE bgcm_02_model

#endif
