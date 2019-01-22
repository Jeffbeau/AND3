MODULE stpctl
   !!==============================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean space and time domain variables 
   USE in_out_manager  ! I/O manager
   USE solisl          ! ???
   USE diawri          ! ocean output file 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE dynspg_oce      ! pressure gradient schemes 

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC stp_ctl           ! routine called by step.F90
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp_ctl( kt, kindic )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol 
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !! History :
      !!        !  91-03  ()
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  97-06  (A.M. Treguier)
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence

      !! * local declarations
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ii, ij, ik              ! temporary integers
      REAL(wp) ::   zumax, zsmin            ! temporary scalars
      INTEGER, DIMENSION(3) ::   ilocu      ! 
      INTEGER, DIMENSION(2) ::   ilocs      ! 
      CHARACTER(len=80) :: clname
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/stpctl.F90,v 1.9 2006/03/21 07:57:56 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!----------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         ! open time.step file
         clname = 'time.step'
         CALL ctlopn( numstp, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
      ENDIF

      ! save the current time step in numstp
      ! ------------------------------------
!      IF(lwp) WRITE(numstp,9100) kt
!      IF(lwp) REWIND(numstp)
!!DB -- on some machines, the above does not write every timestep
!!      On drakes, the below made no difference.
      if(lwp) then
         WRITE(numstp,9100) kt
         call flush(numstp)
         REWIND(numstp)
      endif

9100  FORMAT(1x, i8)


      ! elliptic solver statistics (if required)
      ! --------------------------
      IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
      ! Solver
      IF(lwp) WRITE(numsol,9200) kt, niter, res, SQRT(epsr)/eps

      ! Islands (if exist)
      IF( lk_isl )   CALL isl_stp_ctl( kt, kindic )


      ! Output in numwso and numwvo IF kindic<0
      ! ---------------------------------------
      !    (i.e. problem for the solver)
      IF( kindic < 0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' stpctl: the elliptic solver DO not converge or explode'
            WRITE(numout,*) ' ====== '
            WRITE(numout,9200) kt, niter, res, sqrt(epsr)/eps
            WRITE(numout,*)
            WRITE(numout,*) ' stpctl: output of last fields in numwso'
            WRITE(numout,*) '                                  numwvo'
            WRITE(numout,*) ' ======  *******************************'
         ENDIF
         CALL dia_wri( kt, kindic )
      ENDIF
      ENDIF

9200  FORMAT(' it :', i8, ' niter :', i4, ' res :',e20.10,' b :',e20.10)

      ! Test maximum of velocity (zonal only)
      ! ------------------------
      !! zumax = MAXVAL( ABS( un(:,:,:) ) )   ! slower than the following loop on NEC SX5
      zumax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zumax = MAX(zumax,ABS(un(ji,jj,jk)))
          END DO 
        END DO 
      END DO        
      IF( lk_mpp )   CALL mpp_max( zumax )   ! max over the global domain

      IF( MOD( kt, nwrite ) == 1 ) THEN
         IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(U) max: ', zumax
      ENDIF
      IF( zumax > 20.) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(un),umask,zumax,ii,ij,ik)
         ELSE
            ilocu = MAXLOC( ABS( un(:,:,:) ) )
            ii = ilocu(1) + nimpp - 1
            ij = ilocu(2) + njmpp - 1
            ik = ilocu(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl: the zonal velocity is larger than 20 m/s'
            WRITE(numout,*) ' ====== '
            WRITE(numout,9400) kt, zumax, ii, ij, ik
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         kindic  = -3

         CALL dia_wri( kt, kindic )
      ENDIF
9400  FORMAT (' kt=',i6,' max abs(U): ',1pg11.4,', i j k: ',3i4)


      ! Test minimum of salinity
      ! ------------------------
      !! zsmin = MINVAL( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )    
      !                slower than the following loop on NEC SX5
      zsmin = 100.e0
      DO jj = 2, jpjm1
         DO ji = 1, jpi
            IF( tmask(ji,jj,1) == 1) zsmin = MIN(zsmin,sn(ji,jj,1))
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_min( zsmin )   ! min over the global domain

      IF( MOD( kt, nwrite ) == 1 ) THEN
         IF(lwp) WRITE(numout,*) ' ==>> time-step= ',kt,' SSS min:', zsmin
      ENDIF
      IF( zsmin < 0.) THEN 
         IF (lk_mpp) THEN
            CALL mpp_minloc ( sn(:,:,1),tmask(:,:,1), zsmin, ii,ij )
         ELSE
            ilocs = MINLOC( sn(:,:,1), mask = tmask(:,:,1) == 1.e0 )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
         END IF

         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) 'stp_ctl : NEGATIVE sea surface salinity'
            WRITE(numout,*) '======= '
            WRITE(numout,9500) kt, zsmin, ii, ij
            WRITE(numout,*)
            WRITE(numout,*) '          output of last fields in numwso'
         ENDIF
         IF( kindic < 0 ) THEN
            IF(lwp) WRITE(numout,*) ' stpctl diabort done. We wont do it again '
         ELSE 
            kindic  = -3
            CALL dia_wri(kt,kindic)
         ENDIF
      ENDIF
9500  FORMAT (' kt=',i6,' min SSS: ',1pg11.4,', i j: ',2i4)

   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
