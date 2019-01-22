MODULE dynspg
   !!======================================================================
   !!                       ***  MODULE  dynspg  ***
   !! Ocean dynamics:  surface pressure gradient control
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dyn_spg     : update the dynamics trend with the lateral diffusion
   !!   dyn_spg_ctl : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE obc_oce        ! ocean open boundary conditions
   USE dynspg_oce     ! surface pressure gradient variables
   USE dynspg_exp     ! surface pressure gradient     (dyn_spg_exp routine)
   USE dynspg_ts      ! surface pressure gradient     (dyn_spg_ts  routine)
   USE dynspg_flt     ! surface pressure gradient     (dyn_spg_flt routine)
   USE dynspg_rl      ! surface pressure gradient     (dyn_spg_rl  routine)
   USE dynspg_exp_jki ! surface pressure gradient (dyn_spg_exp_jki routine)
   USE dynspg_ts_jki  ! surface pressure gradient (dyn_spg_ts_jki  routine)
   USE dynspg_flt_jki ! surface pressure gradient (dyn_spg_flt_jki routine)
   USE trdmod         ! ocean dynamics trends
   USE trdmod_oce     ! ocean variables trends
   USE prtctl         ! Print control                     (prt_ctl routine)
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg         ! routine called by step module

   !! * module variables
   INTEGER ::                        &
      nspg = 0                         ! type of surface pressure gradient scheme
      !                                ! defined from lk_dynspg_... 

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg.F90,v 1.4 2005/12/29 10:51:26 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg  ***
      !!
      !! ** Purpose :   compute the lateral ocean dynamics physics.
      !!
      !! History :
      !!   9.0  !  05-12  (C. Talandier, G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in  ) ::   kt     ! ocean time-step index
      INTEGER, INTENT( out ) ::   kindic ! solver flag

      !! * local declarations
      REAL(wp) ::   z2dt                      ! temporary scalar
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztrdu, ztrdv                         ! 3D temporary workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 )   CALL dyn_spg_ctl      ! initialisation & control of options

      IF( l_trddyn )   THEN                      ! temporary save of ta and sa trends
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF

      SELECT CASE ( nspg )                       ! compute surf. pressure gradient trend and add it to the general trend
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL dyn_spg_exp    ( kt )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg0 - Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL dyn_spg_ts     ( kt )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg1 - Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL dyn_spg_flt  ( kt, kindic )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg2 - Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL dyn_spg_exp_jki( kt )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg10- Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL dyn_spg_ts_jki ( kt )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg12- Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         CALL dyn_spg_flt_jki( kt, kindic )
         IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg13- Ua: ', mask1=umask, &
               &                    tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      CASE (  0 )                                       ! explicit
         CALL dyn_spg_exp    ( kt )
      CASE (  1 )                                       ! time-splitting
         CALL dyn_spg_ts     ( kt )
      CASE (  2 )                                       ! filtered
         CALL dyn_spg_flt    ( kt, kindic )
      CASE (  3 )                                       ! rigid lid
         CALL dyn_spg_rl     ( kt, kindic )

      CASE ( 10 )                                       ! explicit with j-k-i loop
         CALL dyn_spg_exp_jki( kt )
      CASE ( 11 )                                       ! time-splitting with j-k-i loop
         CALL dyn_spg_ts_jki ( kt )
      CASE ( 12 )                                       ! filtered with j-k-i loop
         CALL dyn_spg_flt_jki( kt, kindic )
      END SELECT

      !                                          ! save the horizontal diffusive trends for further diagnostics
      IF( l_trddyn )   THEN
         SELECT CASE ( nspg )
         CASE ( 0, 1, 3, 10, 11 )
            ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
            ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CASE( 2, 12 )
            z2dt = 2. * rdt
            IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt
            ztrdu(:,:,:) = ( ua(:,:,:) - ub(:,:,:) ) / z2dt - ztrdu(:,:,:)
            ztrdv(:,:,:) = ( va(:,:,:) - vb(:,:,:) ) / z2dt - ztrdv(:,:,:)
         END SELECT
         CALL trd_mod( ztrdu, ztrdv, jpdtdspg, 'DYN', kt )
      ENDIF

      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' spg  - Ua: ', mask1=umask, &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )

   END SUBROUTINE dyn_spg


   SUBROUTINE dyn_spg_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between cpp options for 
      !!      surface pressure gradient schemes
      !!
      !! History :
      !!   9.0  !  05-10  (V. Garnier)  Original code : spg re-organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ioptio
      !!----------------------------------------------------------------------

      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_spg_ctl : choice of the surface pressure gradient scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '     Explicit free surface                  lk_dynspg_exp = ', lk_dynspg_exp
         WRITE(numout,*) '     Free surface with time splitting       lk_dynspg_ts  = ', lk_dynspg_ts
         WRITE(numout,*) '     Filtered free surface cst volume       lk_dynspg_flt = ', lk_dynspg_flt
         WRITE(numout,*) '     Rigid-lid case                         lk_dynspg_rl  = ', lk_dynspg_rl
      ENDIF

      ! Control of surface pressure gradient scheme options
      ! ---------------------------------------------------
      ioptio = 0
      IF(lk_dynspg_exp)   ioptio = ioptio + 1
      IF(lk_dynspg_ts )   ioptio = ioptio + 1
      IF(lk_dynspg_flt)   ioptio = ioptio + 1
      IF(lk_dynspg_rl )   ioptio = ioptio + 1

      IF( ( ioptio > 1 .AND. .NOT. lk_esopa ) .OR. ioptio == 0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) ' Choose only one surface pressure gradient scheme with a key cpp'
         nstop = nstop + 1
      ENDIF

      IF( lk_esopa     )   nspg = -1
      IF( lk_dynspg_exp)   nspg =  0
      IF( lk_dynspg_ts )   nspg =  1
      IF( lk_dynspg_flt)   nspg =  2
      IF( lk_dynspg_rl )   nspg =  3
      IF( lk_jki       )   nspg =  nspg + 10
      IF( nspg == 13   )   nspg =  3

      IF( lk_esopa     )   nspg = -1

     IF(lwp) THEN
         WRITE(numout,*)
         IF( nspg == -1 )   WRITE(numout,*) '     ESOPA test All scheme used except rigid-lid'
         IF( nspg ==  0 )   WRITE(numout,*) '     explicit free surface'
         IF( nspg ==  1 )   WRITE(numout,*) '     free surface with time splitting scheme'
         IF( nspg ==  2 )   WRITE(numout,*) '     filtered free surface'
         IF( nspg ==  3 )   WRITE(numout,*) '     rigid-lid'
         IF( nspg == 10 )   WRITE(numout,*) '     explicit free surface with j-k-i loop'
         IF( nspg == 11 )   WRITE(numout,*) '     time splitting free surface with j-k-i loop'
         IF( nspg == 12 )   WRITE(numout,*) '     filtered free surface with j-k-i loop'
      ENDIF

      ! Control of timestep choice
      ! --------------------------
      IF( lk_dynspg_ts ) THEN
         IF( MOD( rdt , rdtbt ) /= 0. ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' The barotropic timestep must be an integer divisor of the baroclinic timestep'
            nstop = nstop + 1
         ENDIF
      ENDIF

#if key_obc
      ! Conservation of ocean volume (key_dynspg_flt)
      ! ---------------------------------------------
      IF( lk_dynspg_flt ) ln_vol_cst = .true.

      ! Application of Flather's algorithm at open boundaries
      ! -----------------------------------------------------
      IF( lk_dynspg_flt ) ln_obc_fla = .false.
      IF( lk_dynspg_exp ) ln_obc_fla = .true.
      IF( lk_dynspg_ts  ) ln_obc_fla = .true.
#endif

   END SUBROUTINE dyn_spg_ctl

  !!======================================================================
END MODULE dynspg
