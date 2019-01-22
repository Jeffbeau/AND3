MODULE zdfric
   !!======================================================================
   !!                       ***  MODULE  zdfric  ***
   !! Ocean physics:  vertical mixing coefficient compute from the local
   !!                 Richardson number dependent formulation
   !!======================================================================
#if defined key_zdfric   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_zdfric'                                             Kz = f(Ri)
   !!----------------------------------------------------------------------
   !!   zdf_ric      : update momentum and tracer Kz from the Richardson
   !!                  number computation
   !!   zdf_ric_init : initialization, namelist read, & parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics
!  USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary condition (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_ric   ! called by step.F90

   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .TRUE.    !: Richardson vertical mixing flag

   !! * Module variables
   INTEGER ::               & !!! namric   richardson number dependent Kz
      nric  = 2                ! coefficient of the parameterization
   REAL(wp) ::              & !!! namric   richardson number dependent Kz
      avmri = 100.e-4_wp ,  &  ! maximum value of the vertical eddy viscosity
      alp   =   5._wp          ! coefficient of the parameterization
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
      tmric                    ! coef. for the horizontal mean at t-point

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfric.F90,v 1.3 2005/03/27 18:35:26 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_ric( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdfric  ***
      !!                    
      !! ** Purpose :   Compute the before eddy viscosity and diffusivity as
      !!      a function of the local richardson number.
      !!
      !! ** Method  :   Local richardson number dependent formulation of the 
      !!      vertical eddy viscosity and diffusivity coefficients. the eddy
      !!      coefficients are given by:
      !!              avm = avm0 + avmb
      !!              avt = avm0 / (1 + alp*ri)
      !!      with    ri  = N^2 / dz(u)**2
      !!                  = e3w**2 * rn2/[ mi( dk(ub) )+mj( dk(vb) ) ]
      !!              avm0= avmri / (1 + alp*ri)**nric
      !!      Where ri is the before local Richardson number, avmri the maximum
      !!      value reaches by the vertical eddy coefficients, avmb and avtb
      !!      the background (or minimum) values of these coefficients for
      !!      momemtum and tracers, and alp, nric are adjustable parameters.
      !!      typical values used are : avm0=1.e-2 m2/s, avmb=1.e-6 m2/s
      !!      avtb=1.e-7 m2/s, alp=5. and nric=2.
      !!      this formulation needs ri>=0 : ri is set to zero if dz(rau)<0.
      !!      a numerical threshold is impose on the vertical shear (1.e-20)
      !!        N.B. the mask are required for implicit scheme, and surface
      !!      and bottom value already set in inimix.F
      !!
      !! References :
      !!      pacanowski & philander 1981, j. phys. oceanogr., 1441-1451.
      !! History :
      !!        !  87-09  (P. Andrich)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  96-01  (G. Madec)  complet rewriting of multitasking
      !!                                  suppression of common work arrays
      !!        !  97-06 (G. Madec)  complete rewriting of zdfmix
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      REAL(wp) ::   &
         zcoef, zdku, zdkv, zri, z05alp     ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   zwx ! temporary workspace

      IF( kt == nit000  ) CALL zdf_ric_init            ! Initialization (first time-step only)

      !                                                ! ===============
      DO jk = 2, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Richardson number (put in zwx(ji,jj))
         ! -----------------
         ! minimum value set to zero
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zcoef = 0.5 / fse3w(ji,jj,jk)
               ! shear of horizontal velocity
               zdku = zcoef * (  ub(ji-1,jj,jk-1) + ub(ji,jj,jk-1)   &
                                -ub(ji-1,jj,jk  ) - ub(ji,jj,jk  )  )
               zdkv = zcoef * (  vb(ji,jj-1,jk-1) + vb(ji,jj,jk-1)   &
                                -vb(ji,jj-1,jk  ) - vb(ji,jj,jk  )  )
               ! richardson number (minimum value set to zero)
               zri = rn2(ji,jj,jk) / ( zdku*zdku + zdkv*zdkv + 1.e-20 )
               zwx(ji,jj) = MAX( zri, 0.e0 )
            END DO
         END DO

         ! Boundary condition on zwx   (sign unchanged)
         CALL lbc_lnk( zwx, 'W', 1. )


         ! Vertical eddy viscosity and diffusivity coefficients
         ! -------------------------------------------------------
         ! Eddy viscosity coefficients
         z05alp = 0.5 * alp
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               avmu(ji,jj,jk) = umask(ji,jj,jk)   &
                              * avmri / ( 1. + z05alp*( zwx(ji+1,jj)+zwx(ji,jj) ) )**nric
               avmv(ji,jj,jk) = vmask(ji,jj,jk)   &
                              * avmri / ( 1. + z05alp*( zwx(ji,jj+1)+zwx(ji,jj) ) )**nric
            END DO
         END DO

         ! Eddy diffusivity coefficients
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               avt(ji,jj,jk) = tmric(ji,jj,jk) / ( 1. + alp * zwx(ji,jj) )   &
                             * (  avmu(ji,jj,jk) + avmu(ji-1, jj ,jk)        &
                                + avmv(ji,jj,jk) + avmv( ji ,jj-1,jk)  )     &
                             + avtb(jk) * tmask(ji,jj,jk)
            END DO
         END DO

         ! Add the background coefficient on eddy viscosity
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               avmu(ji,jj,jk) = avmu(ji,jj,jk) + avmb(jk) * umask(ji,jj,jk)
               avmv(ji,jj,jk) = avmv(ji,jj,jk) + avmb(jk) * vmask(ji,jj,jk)
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! Boundary conditions on (avt,avmu,avmv)   (unchanged sign)
      ! -----------------------===============
      CALL lbc_lnk( avt , 'W', 1. )
      CALL lbc_lnk( avmu, 'U', 1. )
      CALL lbc_lnk( avmv, 'V', 1. )

   END SUBROUTINE zdf_ric


   SUBROUTINE zdf_ric_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdfbfr_init  ***
      !!                    
      !! ** Purpose :   Initialization of the vertical eddy diffusivity and
      !!      viscosity coef. for the Richardson number dependent formulation.
      !!
      !! ** Method  :   Read the namric namelist and check the parameter values
      !!
      !! ** input   :   Namelist namric
      !!
      !! ** Action  :   increase by 1 the nstop flag is setting problem encounter
      !!
      !! history :
      !!  8.5  !  02-06  (G. Madec)  original code
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER :: ji, jj, jk        ! dummy loop indices

      NAMELIST/namric/ avmri, alp, nric
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      ! Read Namelist namric : richardson number dependent Kz
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namric )


      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zdf_ric : Ri depend vertical mixing scheme'
      IF(lwp) WRITE(numout,*) '======='
      IF(lwp) WRITE(numout,*) '          Namelist namric : set Kz(Ri) parameters'

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '             maximum vertical viscosity     avmri  = ', avmri
         WRITE(numout,*) '             coefficient                    alp    = ', alp
         WRITE(numout,*) '             coefficient                    nric   = ', nric
         WRITE(numout,*)
      ENDIF


      ! Work arrays for Ri number formulation
      ! -------------------------------------

      ! background eddy viscosity and diffusivity profiles
      avmb(:) = avm0
      avtb(:) = avt0

      ! background profile of avm (fit the theoretical/observational
      !     profile shown by Krauss (1990) and avt
!!!   avtb(:) = 1.e-5 + 2.8e-8 * gdepw(:) ! m2/s

      ! Increase the background in the surface layers
      avmb(1) = 10.  * avmb(1)      ;      avtb(1) = 10.  * avtb(1)
      avmb(2) = 10.  * avmb(2)      ;      avtb(2) = 10.  * avtb(2)
      avmb(3) =  5.  * avmb(3)      ;      avtb(3) =  5.  * avtb(3)
      avmb(4) =  2.5 * avmb(4)      ;      avtb(4) =  2.5 * avtb(4)

      ! weighting mean array tmric for 4 T-points which accounts for coastal boundary conditions.
      DO jk = 1, jpk
         DO jj = 2, jpj
            DO ji = 2, jpi
               tmric(ji,jj,jk) =  tmask(ji,jj,jk)                                  &
                               / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
                                         + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
            END DO
         END DO
      END DO

      tmric(:,1,:) = 0.e0

      ! Initialization of vertical eddy coef. to the background value
      DO jk = 1, jpk
         avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
         avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
         avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
      END DO

   END SUBROUTINE zdf_ric_init

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :              NO Richardson dependent vertical mixing
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .FALSE.   !: Richardson mixing flag
CONTAINS
   SUBROUTINE zdf_ric( kt )        ! Dummy routine
!      WRITE(*,*) 'zdf_ric: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ric
#endif

   !!======================================================================
END MODULE zdfric
