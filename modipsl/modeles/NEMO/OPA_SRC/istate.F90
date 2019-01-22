MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_tem    : analytical profile for initial Temperature
   !!   istate_sal    : analytical profile for initial Salinity
   !!   istate_eel    : initial state setting of EEL R5 configuration
   !!   istate_gyre   : initial state setting of GYRE configuration
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE daymod          ! 
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE wzvmod          ! verctical velocity               (wzv     routine)
   USE dtatem          ! temperature data                 (dta_tem routine)
   USE dtasal          ! salinity data                    (dta_sal routine)
   USE restart         ! ocean restart                   (rst_read routine)
   USE solisl          ! ???

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC istate_init   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/istate.F90,v 1.11 2006/04/10 15:46:04 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE istate_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers.
      !!
      !! ** Method  :
      !!
      !! History :
      !!   4.0  !  91-03  ()  Original code
      !!        !  91-11  (G. Madec)
      !!   9.0  !  03-09  (G. Madec)  F90: Free form, modules, orthogonality
      !!----------------------------------------------------------------------
      !! * Local declarations
      !!----------------------------------------------------------------------


      ! Initialization to zero
      ! ----------------------

      !     before fields       !       now fields        !      after fields       !
     ub   (:,:,:) = 0.e0   ;   un   (:,:,:) = 0.e0   ;   ua   (:,:,:) = 0.e0
     vb   (:,:,:) = 0.e0   ;   vn   (:,:,:) = 0.e0   ;   va   (:,:,:) = 0.e0
     wn   (:,:,:) = 0.e0   ;
     rotb (:,:,:) = 0.e0   ;   rotn (:,:,:) = 0.e0   ;
     hdivb(:,:,:) = 0.e0   ;   hdivn(:,:,:) = 0.e0   ;
     
     tb   (:,:,:) = 0.e0   ;   tn   (:,:,:) = 0.e0   ;   ta   (:,:,:) = 0.e0
     sb   (:,:,:) = 0.e0   ;   sn   (:,:,:) = 0.e0   ;   sa   (:,:,:) = 0.e0
     
     rhd  (:,:,:) = 0.e0
     rhop (:,:,:) = 0.e0
     rn2  (:,:,:) = 0.e0 
     
#if defined key_dynspg_rl
      ! rigid-lid formulation
      bsfb(:,:) = 0.e0      ! before barotropic stream-function
      bsfn(:,:) = 0.e0      ! now    barotropic stream-function
      bsfd(:,:) = 0.e0      ! barotropic stream-function trend
#endif
      ! free surface formulation
      sshb(:,:) = 0.e0      ! before sea-surface height
      sshn(:,:) = 0.e0      ! now    sea-surface height


      IF( ln_rstart ) THEN                    ! Restart from a file
         !                                    ! -------------------
         neuler = 1                              ! Set time-step indicator at nit000 (leap-frog)
         CALL rst_read                           ! Read the restart file
      ELSE
         !                                    ! Start from rest
         !                                    ! ---------------
         neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
         adatrj = 0._wp
!!DB: OLD
!         IF( cp_cfg == 'eel' ) THEN
!            CALL istate_eel                      ! EEL   configuration : start from pre-defined
!            !                                    !                       velocity and thermohaline fields
!         ELSEIF( cp_cfg == 'gyre' ) THEN         
!            CALL istate_gyre                     ! GYRE  configuration : start from pre-defined temperature
!            !                                    !                       and salinity fields 
!         ELSE
         !                                       ! Other configurations: Initial temperature and salinity fields
!#if defined key_dtatem
         CALL dta_tem( nit000 )                  ! read 3D temperature data
         tb(:,:,:) = t_dta(:,:,:)                ! use temperature data read
         tn(:,:,:) = t_dta(:,:,:)
!#else
!         IF(lwp) WRITE(numout,*)                 ! analytical temperature profile
!         IF(lwp) WRITE(numout,*)' Temperature initialization using an analytic profile'
!         CALL istate_tem
!#endif
!#if defined key_dtasal
         CALL dta_sal( nit000 )                  ! read 3D salinity data
         sb(:,:,:) = s_dta(:,:,:)                ! use salinity data read
         sn(:,:,:) = s_dta(:,:,:)
!#else
!         ! No salinity data
!         IF(lwp)WRITE(numout,*)                  ! analytical salinity profile
!         IF(lwp)WRITE(numout,*)' Salinity initialisation using a constant value'
!         CALL istate_sal
!#endif
         !         ENDIF

      ENDIF
      !                                       ! Vertical velocity
      !                                       ! -----------------
      CALL wzv( nit000 )                         ! from horizontal divergence

   END SUBROUTINE istate_init





   SUBROUTINE istate_uvg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE istate_uvg  ***
      !!
      !! ** Purpose :   Compute the geostrophic velocities from (tn,sn) fields
      !!
      !! ** Method  :   Using the hydrostatic hypothesis the now hydrostatic 
      !!      pressure is computed by integrating the in-situ density from the
      !!      surface to the bottom.
      !!                 p=integral [ rau*g dz ]
      !!
      !! History :
      !!   8.1  !  01-09  (M. Levy, M. Ben Jelloul)  Original code
      !!   8.5  !  02-09  (G. Madec)  F90: Free form
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE eosbn2          ! eq. of state, Brunt Vaisala frequency (eos     routine)
      USE dynspg          ! surface pressure gradient             (dyn_spg routine)
      USE divcur          ! hor. divergence & rel. vorticity      (div_cur routine)
      USE lbclnk          ! ocean lateral boundary condition (or mpp link)

      !! * Local declarations
      INTEGER ::   ji, jj, jk        ! dummy loop indices
      INTEGER ::   indic             ! ???
      REAL(wp) ::   &
         zmsv, zphv, zmsu, zphu,  &  ! temporary scalars
         zalfg
      REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
         zprn                        ! workspace
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*) 
      IF(lwp) WRITE(numout,*) 'istate_uvg : Start from Geostrophy'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      ! Compute the now hydrostatic pressure
      ! ------------------------------------

      zalfg = 0.5 * grav * rau0
      ! Surface value
      zprn(:,:,1) = zalfg * fse3w(:,:,1) * ( 1 + rhd(:,:,1) )

      ! Vertical integration from the surface
      DO jk = 2, jpkm1
         zprn(:,:,jk) = zprn(:,:,jk-1)   &
            &         + zalfg * fse3w(:,:,jk) * ( 2. + rhd(:,:,jk) + rhd(:,:,jk-1) )
      END DO  

      ! Compute geostrophic balance
      ! ---------------------------

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vertor opt.
               zmsv = 1. / MAX(  umask(ji-1,jj+1,jk) + umask(ji  ,jj+1,jk)   &
                               + umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) , 1.  )
               zphv = ( zprn(ji  ,jj+1,jk) - zprn(ji-1,jj+1,jk) ) * umask(ji-1,jj+1,jk) / e1u(ji-1,jj+1)   &
                    + ( zprn(ji+1,jj+1,jk) - zprn(ji  ,jj+1,jk) ) * umask(ji  ,jj+1,jk) / e1u(ji  ,jj+1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji-1,jj  ,jk) ) * umask(ji-1,jj  ,jk) / e1u(ji-1,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji  ,jj  ,jk) ) * umask(ji  ,jj  ,jk) / e1u(ji  ,jj  )
               zphv = 1. / rau0 * zphv * zmsv * vmask(ji,jj,jk)

               zmsu = 1. / MAX(  vmask(ji+1,jj  ,jk) + vmask(ji  ,jj  ,jk)   &
                               + vmask(ji+1,jj-1,jk) + vmask(ji  ,jj-1,jk) , 1.  )
               zphu = ( zprn(ji+1,jj+1,jk) - zprn(ji+1,jj  ,jk) ) * vmask(ji+1,jj  ,jk) / e2v(ji+1,jj  )   &
                    + ( zprn(ji  ,jj+1,jk) - zprn(ji  ,jj  ,jk) ) * vmask(ji  ,jj  ,jk) / e2v(ji  ,jj  )   &
                    + ( zprn(ji+1,jj  ,jk) - zprn(ji+1,jj-1,jk) ) * vmask(ji+1,jj-1,jk) / e2v(ji+1,jj-1)   &
                    + ( zprn(ji  ,jj  ,jk) - zprn(ji  ,jj-1,jk) ) * vmask(ji  ,jj-1,jk) / e2v(ji  ,jj-1)
               zphu = 1. / rau0 * zphu * zmsu * umask(ji,jj,jk)

               ! Compute the geostrophic velocities
               un(ji,jj,jk) = -2. * zphu / ( ff(ji,jj) + ff(ji  ,jj-1) )
               vn(ji,jj,jk) =  2. * zphv / ( ff(ji,jj) + ff(ji-1,jj  ) )
            END DO
         END DO
      END DO

      IF(lwp) WRITE(numout,*) '         we force to zero bottom velocity'

      ! Susbtract the bottom velocity (level jpk-1 for flat bottom case)
      ! to have a zero bottom velocity

      DO jk = 1, jpkm1
         un(:,:,jk) = ( un(:,:,jk) - un(:,:,jpkm1) ) * umask(:,:,jk)
         vn(:,:,jk) = ( vn(:,:,jk) - vn(:,:,jpkm1) ) * vmask(:,:,jk)
      END DO

      CALL lbc_lnk( un, 'U', -1. )
      CALL lbc_lnk( vn, 'V', -1. )
      
      ub(:,:,:) = un(:,:,:)
      vb(:,:,:) = vn(:,:,:)
      
      ! WARNING !!!!!
      ! after initializing u and v, we need to calculate the initial streamfunction bsf.
      ! Otherwise, only the trend will be computed and the model will blow up (inconsistency).
      
      ! to do that, we call dyn_spg with a special trick:
      ! we fill ua and va with the velocities divided by dt,
      ! and the streamfunction will be brought to the right
      ! value assuming the velocities have been set up in
      ! one time step.
      ! we then set bsfd to zero (first guess for next step
      ! is d(psi)/dt = 0.)

      !  sets up s false trend to calculate the barotropic
      !  streamfunction.

      ua(:,:,:) = ub(:,:,:) / rdt
      va(:,:,:) = vb(:,:,:) / rdt

      ! calls dyn_spg. we assume euler time step, starting from rest.
      indic = 0
      CALL dyn_spg( nit000, indic )       ! surface pressure gradient

      ! the new velocity is ua*rdt

      CALL lbc_lnk( ua, 'U', -1. )
      CALL lbc_lnk( va, 'V', -1. )

      ub(:,:,:) = ua(:,:,:) * rdt
      vb(:,:,:) = va(:,:,:) * rdt
      ua(:,:,:) = 0.e0
      va(:,:,:) = 0.e0
      un(:,:,:) = ub(:,:,:)
      vn(:,:,:) = vb(:,:,:)
       
#if defined key_dynspg_rl
      IF( lk_isl )   bsfb(:,:) = bsfn(:,:)          ! Put bsfb to zero
#endif

      ! Compute the divergence and curl

      CALL div_cur( nit000 )            ! now horizontal divergence and curl

      hdivb(:,:,:) = hdivn(:,:,:)       ! set the before to the now value
      rotb (:,:,:) = rotn (:,:,:)       ! set the before to the now value

   END SUBROUTINE istate_uvg

   !!=====================================================================
END MODULE istate
