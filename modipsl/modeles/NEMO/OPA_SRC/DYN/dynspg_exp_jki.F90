MODULE dynspg_exp_jki
   !!======================================================================
   !!                   ***  MODULE  dynspg_exp_jki  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
#if ( defined key_dynspg_exp && defined key_autotasking )   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_exp'                              explicit free surface
   !!   'key_autotasking'                                        j-k-i loop
   !!----------------------------------------------------------------------
   !!   dyn_spg_exp_jki  : update the momentum trend with the surface 
   !!                      pressure gradient in the free surface constant  
   !!                      volume case with vector optimization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE ocesbc          ! ocean surface boundary condition
   USE obc_oce         ! Lateral open boundary condition
   USE obc_par         ! open boundary condition parameters
   USE obcdta          ! open boundary condition data     (obc_dta_bt routine)
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_exp_jki  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_exp_jki.F90,v 1.1 2005/12/29 10:44:29 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_exp_jki( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_exp_jki  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!      gradient in case of explicit free surface formulation and 
      !!      add it to the general trend of momentum equation. Compute
      !!      the free surface.
      !!
      !! ** Method  :   Explicit free surface formulation. The surface pressure
      !!      gradient is given by:
      !!         spgu = 1/rau0 d/dx(ps) =  g/e1u di( sshn )
      !!         spgv = 1/rau0 d/dy(ps) =  g/e2v dj( sshn )
      !!      -1- Compute the now surface pressure gradient
      !!      -2- Add it to the general trend
      !!      -3- Compute the horizontal divergence of velocities
      !!      - the now divergence is given by :
      !!         zhdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      - integrate the horizontal divergence from the bottom 
      !!         to the surface
      !!      - apply lateral boundary conditions on zhdivn
      !!      -4- Estimate the after sea surface elevation from the kinematic
      !!         surface boundary condition:
      !!         zssha = sshb - 2 rdt ( zhdiv + emp )
      !!      - Time filter applied on now sea surface elevation to avoid
      !!         the divergence of two consecutive time-steps and swap of free
      !!         surface arrays to start the next time step:
      !!         sshb = sshn + atfp * [ sshb + zssha - 2 sshn ]
      !!         sshn = zssha
      !!      - apply lateral boundary conditions on sshn
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References :
      !!
      !! History :
      !!   9.0  !  05-11  (V. Garnier, G. Madec, L. Bessieres) Original code
      !!
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in )  ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   z2dt, zraur, zssha       ! temporary scalars 
      REAL(wp), DIMENSION(jpi,jpj)    ::  &  ! temporary arrays
         &         zhdiv
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_exp_jki : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   (explicit free surface, j-k-i loop)'

         ! set to zero free surface specific arrays
         spgu(:,:) = 0.e0                     ! surface pressure gradient (i-direction)
         spgv(:,:) = 0.e0                     ! surface pressure gradient (j-direction)
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      ! read or estimate sea surface height and vertically integrated velocities
      IF( lk_obc )   CALL obc_dta_bt( kt, 0 )
      z2dt = 2. * rdt                                       ! time step: leap-frog
      IF( neuler == 0 .AND. kt == nit000 ) z2dt = rdt       ! time step: Euler if restart from rest
      zraur = 1.e0 / rauw

!CDIR PARALLEL DO
!$OMP PARALLEL DO
      !                                                ! =============== !
      DO jj = 2, jpjm1                                 !  Vertical slab  !
         !                                             ! =============== !

         ! Surface pressure gradient (now)
         DO ji = 2, jpim1
            spgu(ji,jj) = - grav * ( sshn(ji+1,jj) - sshn(ji,jj) ) / e1u(ji,jj)
            spgv(ji,jj) = - grav * ( sshn(ji,jj+1) - sshn(ji,jj) ) / e2v(ji,jj)
         END DO

         ! Add the surface pressure trend to the general trend
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ua(ji,jj,jk) + spgu(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) + spgv(ji,jj)
            END DO
         END DO
     
         !  Vertical integration of horizontal divergence of velocities
         zhdiv(:,jj) = 0.e0
         DO jk = jpkm1, 1, -1
            DO ji = 2, jpim1
               zhdiv(ji,jj) = zhdiv(ji,jj) + (  e2u(ji  ,jj  ) * fse3u(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)      &
                  &                           - e2u(ji-1,jj  ) * fse3u(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)      &
                  &                           + e1v(ji  ,jj  ) * fse3v(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)      &
                  &                           - e1v(ji  ,jj-1) * fse3v(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)  )   &
                  &                        / ( e1t(ji,jj) * e2t(ji,jj) )
            END DO
         END DO

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  ) THEN
            IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
         ENDIF
         IF( lp_obc_west  ) THEN
            IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
         ENDIF
         IF( lp_obc_north ) THEN
            IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
         ENDIF
         IF( lp_obc_south ) THEN
            IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
         ENDIF
#endif

         ! Sea surface elevation time stepping
         IF( neuler == 0 .AND. kt == nit000 ) THEN       ! Euler (forward) time stepping, no time filter
            DO ji = 2, jpim1
               ! after free surface elevation
               zssha = sshb(ji,jj) - rdt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
               ! swap of arrays
               sshb(ji,jj) = sshn(ji,jj)
               sshn(ji,jj) = zssha
            END DO
         ELSE                                            ! Leap-frog time stepping and time filter
            DO ji = 2, jpim1
               ! after free surface elevation
               zssha = sshb(ji,jj) - z2dt * ( zraur * emp(ji,jj) + zhdiv(ji,jj) ) * tmask(ji,jj,1)
               ! time filter and array swap
               sshb(ji,jj) = atfp * ( sshb(ji,jj) + zssha ) + atfp1 * sshn(ji,jj)
               sshn(ji,jj) = zssha
            END DO
         ENDIF
         !                                             ! =============== !
      END DO                                           !     end slab    !
      !                                                ! =============== !

      ! Boundary conditions on sshn
      IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )
 
      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
      ENDIF
      
   END SUBROUTINE dyn_spg_exp_jki

#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart explicit free surface 
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_spg_exp_jki( kt )       ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_spg_exp_jki: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_jki
#endif
   
   !!======================================================================
END MODULE dynspg_exp_jki
