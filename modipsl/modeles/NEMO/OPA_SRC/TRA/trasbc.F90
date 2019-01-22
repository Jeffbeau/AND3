MODULE trasbc
   !!==============================================================================
   !!                       ***  MODULE  trasbc  ***
   !! Ocean active tracers:  surface boundary condition
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space domain variables
   USE ocesbc          ! surface thermohaline fluxes
   USE phycst          ! physical constant
   USE traqsr          ! solar radiation penetration
   USE trdmod          ! ocean trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_sbc              ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trasbc.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_sbc ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sbc  ***
      !!                   
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (flux through the interface, concentration/dilution effect)
      !!      and add it to the general trend of tracer equations.
      !!
      !! ** Method :
      !!      * flux through the air-sea interface:
      !!            - temperature : heat flux q (w/m2). If penetrative solar
      !!         radiation q is only the non solar part of the heat flux, the
      !!         solar part is added in traqsr.F routine.
      !!            ta = ta + q /(rau0 rcp e3t)  for k=1
      !!            - salinity    : no salt flux
      !!      * concentration/dilution effect:
      !!            The surface freshwater flux modify the ocean volume
      !!         and thus the concentration of a tracer and the temperature.
      !!         First order of the effect of surface freshwater exchange 
      !!         for salinity, it can be neglected on temperature (especially
      !!         as the temparature of precipitations and runoffs is usually
      !!         unknown.
      !!            - temperature : we assume that the temperature of both
      !!         precipitations and runoffs is equal to the SST, thus there
      !!         is no additional flux since in this case, the concentration
      !!         dilution effect is balanced by the net heat flux associated
      !!         to the freshwater exchange:
      !!            (Tp P - Te E) + STT (P-E) = 0 when Tp=Te=SST
      !!            - salinity    : evaporation, precipitation and runoff
      !!         water has a zero salinity, thus
      !!            sa = sa + emp * sn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation minus runoff) given in kg/m2/s is divided
      !!         by 1000 kg/m3 (density of plain water) to obtain m/s.
      !!
      !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
      !!                with the tracer surface boundary condition 
      !!              - save the trend it in ttrd ('key_trdtra')
      !!
      !! History :
      !!   8.2  !  98-10  (G. Madec, G. Roullet, M. Imbard)  Original code
      !!   8.2  !  01-02  (D. Ludicone)  sea ice and free surface
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,      & ! use ua as 3D workspace   
                         ztdsa => va         ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj                   ! dummy loop indices
      REAL(wp) ::   zta, zsa, zsrau, zse3t   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_sbc : TRAcer Surface Boundary Condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! 0. initialization
      zsrau = 1. / rauw
#if ! defined key_s_coord
      zse3t = 1. / fse3t(1,1,1)
#endif

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      IF( .NOT.ln_traqsr )   qsr(:,:) = 0.e0   ! no solar radiation penetration

      ! 1. Concentration dillution effect on (t,s)
      DO jj = 2, jpj
         DO ji = fs_2, fs_jpim1   ! vector opt.
#if defined key_s_coord
            zse3t = 1. / fse3t(ji,jj,1)
#endif
            ! temperature : heat flux
            zta = ro0cpr * ( qt(ji,jj) - qsr(ji,jj) ) * zse3t

            ! salinity :  concent./dilut. effect
            zsa = emps(ji,jj) * zsrau * sn(ji,jj,1) * zse3t
            
            ! add the trend to the general tracer trend
            ta(ji,jj,1) = ta(ji,jj,1) + zta
            sa(ji,jj,1) = sa(ji,jj,1) + zsa
         END DO
      END DO

      ! save the trends for diagnostic
      ! sea surface boundary condition tracers trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
         CALL trd_mod(ztdta, ztdsa, jpttdnsr, 'TRA', kt)
      ENDIF
      
      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' sbc  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

   END SUBROUTINE tra_sbc

   !!======================================================================
END MODULE trasbc
