MODULE cla_dynspg
   !!======================================================================
   !!                       ***  cla_dynspg  ***
   !!======================================================================
   !!   dyn_spg      : update the momentum trend with the surface pressure
   !!                  gradient in the free surface constant volume case
   !!                  with vector optimization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE zdf_oce         ! ocean vertical physics
   USE obc_oce         ! Lateral open boundary condition
   USE sol_oce         ! solver variables
   USE phycst          ! physical constants
   USE ocesbc          ! ocean surface boundary condition (fluxes)
   USE flxrnf          ! ocean runoffs
   USE solpcg          ! preconditionned conjugate gradient solver
   USE solsor          ! Successive Over-relaxation solver
   USE solfet          ! FETI solver
   USE obcdyn          ! ocean open boundary condition (obc_dyn routines)
   USE obcvol          ! ocean open boundary condition (obc_vol routines)
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_cla   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/cla_dynspg.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_cla( kt ) 
      !!----------------------------------------------------------------------
      !!              ***  routine dyn_spg_cross_land  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method :
      !!
      !! ** Action :
      !!
      !! History :
      !!        !         (A. Bozec)  Original code
      !!   8.5  !  02-11  (A. Bozec)  F90: Free form and module
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step
      !! * Local declarations
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1        ! temporary integer
      REAL(wp) ::    &    
         zempmed, zempred,   &                ! EMP on Med Sea ans Red Sea
         zwei,   &                            !              
         zisw_rs, zurw_rs, zbrw_rs,      &    ! imposed transport Red sea
         zisw_ms, zurw_ms, zbrw_ms, zmrw_ms   ! imposed transport Med Sea
      !!----------------------------------------------------------------------

      ! Different velocities for straits ( Gibraltar, Bab el Mandeb...)
         
      ! Control print
      ! -------------
      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dynspg_cross_land : cross land advection on surface '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~   pressure '
         IF(lwp) WRITE(numout,*) ' '
      ENDIF

      ! EMP on Mediterranean Sea and Red Sea 
      ! ------------------------------------
      ! compute the emp in Mediterranean Sea
      zempmed = 0.e0
      zwei = 0.e0
      ij0 =  96   ;   ij1 = 110
      ii0 = 141   ;   ii1 = 181
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            zwei    = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
            zempmed = zempmed + emp(ji,jj) * zwei
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( zempmed )      ! sum with other processors value

      ! minus 2 points in Red Sea and 3 in Atlantic 
      ij0 =  96   ;   ij1 =  96
      ii0 = 148   ;   ii1 = 148
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            zempmed = zempmed - emp(ji  ,jj) * tmask(ji  ,jj,1) * e1t(ji  ,jj) * e2t(ji  ,jj)   &
               &              - emp(ji+1,jj) * tmask(ji+1,jj,1) * e1t(ji+1,jj) * e2t(ji+1,jj)   
         END DO
      END DO
      ! we convert in m3
      zempmed = zempmed * 1.e-3

      ! compute the emp in Red Sea   
      zempred = 0.e0
      zwei = 0.e0
      ij0 =  87   ;   ij1 =  96
      ii0 = 148   ;   ii1 = 160
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            zwei      = tmask(ji,jj,1) * e1t(ji,jj) * e2t(ji,jj)
            zempred   = zempred + emp(ji,jj) * zwei
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum( zempred )      ! sum with other processors value

      ! we convert in m3
      zempred = zempred * 1.e-3

      ! New Transport at Bab el Mandeb and Gibraltar
      ! --------------------------------------------

      ! imposed transport at Bab el Mandeb
      zisw_rs = 0.4e6        ! inflow surface water
      zurw_rs = 0.2e6        ! upper recirculation water
!!Alex      zbrw_rs = 1.2e6        ! bottom  recirculation water
      zbrw_rs = 0.5e6        ! bottom  recirculation water

      ! imposed transport at Gibraltar
      zisw_ms  = 0.8e6          ! atlantic-mediterranean  water
      zmrw_ms  = 0.7e6          ! middle recirculation water
      zurw_ms  = 2.5e6          ! upper  recirculation water 
      zbrw_ms  = 3.5e6          ! bottom recirculation water 

      ! Different velocities for straits ( Gibraltar, Bab el Mandeb )
      ! -------------------------------------------------------------

      ! Bab el Mandeb
      ! -------------
      ! 160,88 north point Bab el Mandeb
      ij0 =  88   ;   ij1 =  88
      ii0 = 160   ;   ii1 = 160
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj  ,: ) = 0.e0  !  North East Bab el Mandeb 
         END DO
      END DO
      !                              ! surface
      DO jk = 1,  8                                      
         DO jj = mj0(ij0), mj1(ij1)
            DO ji = mi0(ii0),mi1(ii1)
               ua(ji, jj,jk) = -( ( zisw_rs + zempred ) / 8. ) / ( e2u(ji, jj) * fse3t(ji, jj,jk) )     
            END DO
         END DO
      END DO
      !                              ! deeper
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji, jj,21) = - zbrw_rs / ( e2u(ji, jj) * fse3t(ji, jj,21) )
         END DO
      END DO

      ! 160,87 south point Bab el Mandeb
      ij0 =  87   ;   ij1 =  87
      ii0 = 160   ;   ii1 = 160
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj  ,: ) = 0.e0  !  South East Bab el Mandeb 
         END DO
      END DO
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji, jj,21) =  ( zisw_rs + zbrw_rs ) / ( e2u(ji,jj )*fse3t(ji, jj,21) )      
         END DO
      END DO

      ! Gibraltar
      ! ---------

      ! initialisation of velocity at concerned points 
      ! 139, 101 south point in Gibraltar 
      ij0 = 101   ;   ij1 = 101
      ii0 = 139   ;   ii1 = 139
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj  ,: ) = 0.e0  !  South West Gibraltar
            ua(ji,jj+1,: ) = 0.e0  !  North West Gibraltar
         END DO
      END DO
      !                            ! surface
      DO jk = 1, 14                      
         DO jj = mj0(ij0), mj1(ij1)
            DO ji = mi0(ii0),mi1(ii1)
               ua(ji,jj,jk) =  ( ( zisw_ms + zempmed ) / 14. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
            END DO
         END DO
      END DO
      !                            ! middle circulation
      DO jk = 15, 20                      
         DO jj = mj0(ij0), mj1(ij1)
            DO ji = mi0(ii0),mi1(ii1)
               ua(ji,jj,jk) =  ( zmrw_ms / 6. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
            END DO
         END DO
      END DO
      !                            ! deeper 
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj,21) =             zurw_ms   / ( e2u(ji,jj) * fse3t(ji,jj,21) )
            ua(ji,jj,22) = ( zbrw_ms - zurw_ms ) / ( e2u(ji,jj) * fse3t(ji,jj,22) )
         END DO
      END DO

      ! 139,102 north point in Gibraltar
      ij0 = 102   ;   ij1 = 102
      ii0 = 139   ;   ii1 = 139
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj  ,: ) = 0.e0  !  North West Gibraltar
         END DO
      END DO
      DO jk = 15, 20                      
         DO jj = mj0(ij0), mj1(ij1)
            DO ji = mi0(ii0),mi1(ii1)
               ua(ji,jj,jk) = -( zmrw_ms / 6. ) / ( e2u(ji,jj) * fse3t(ji,jj,jk) ) 
            END DO
         END DO
      END DO
      !                            ! deeper
      DO jj = mj0(ij0), mj1(ij1)
         DO ji = mi0(ii0),mi1(ii1)
            ua(ji,jj,22) = -( zisw_ms + zbrw_ms ) / ( e2u(ji,jj) * fse3t(ji,jj,22) )
         END DO
      END DO

   END SUBROUTINE dyn_spg_cla

   !!======================================================================
END MODULE cla_dynspg
