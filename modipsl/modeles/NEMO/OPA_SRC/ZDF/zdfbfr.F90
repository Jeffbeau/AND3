MODULE zdfbfr
   !!======================================================================
   !!                       ***  MODULE  zdfbfr  ***
   !! Ocean physics: Bottom friction
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   zdf_bfr      : update momentum Kz at the ocean bottom due to the
   !!                  type of bottom friction chosen
   !!   zdf_bfr_init : read in namelist and control the bottom friction
   !!                  parameters.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_bfr    ! called by step.F90

   !! * Module variables
   INTEGER ::             & !!! ** bottom friction namelist (nambfr) **
      nbotfr = 0             ! = 0/1/2/3 type of bottom friction 
   REAL(wp) ::            & !!! ** bottom friction namelist (nambfr) **
      bfri1 = 4.0e-4_wp,  &  ! bottom drag coefficient (linear case) 
      bfri2 = 1.0e-3_wp,  &  ! bottom drag coefficient (non linear case)
      bfeb2 = 2.5e-3_wp      ! background bottom turbulent kinetic energy  (m2/s2)

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfbfr.F90,v 1.5 2005/09/02 15:45:43 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zdf_bfr( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_bfr  ***
      !!                 
      !! ** Purpose :   Applied the bottom friction through a specification of 
      !!      Kz at the ocean bottom.
      !!
      !! ** Method  :   Update the value of avmu and avmv at the ocean bottom 
      !!       level following the chosen friction type (no-slip, free-slip, 
      !!       linear, or quadratic)
      !!
      !! History :
      !!   8.0  !  97-06  (G. Madec, A.-M. Treguier)  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index

      !! * Local declarations
      INTEGER ::   &
         ji, jj,                   &  ! dummy loop indexes
         ikbu, ikbv,               &  ! temporary integers
         ikbum1, ikbvm1               !
      REAL(wp) ::   &
         zvu, zuv, zecu, zecv         ! temporary scalars
      !!----------------------------------------------------------------------


      IF( kt == nit000 )   CALL zdf_bfr_init


      ! Compute avmu, avmv at the ocean bottom
      ! --------------------------------------

      SELECT CASE (nbotfr)

      CASE( 0 )                 ! no-slip boundary condition
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
# endif
               ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
               ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
               ikbum1 = MAX( ikbu-1, 1 )
               ikbvm1 = MAX( ikbv-1, 1 )
               avmu(ji,jj,ikbu) = 2. * avmu(ji,jj,ikbum1)
               avmv(ji,jj,ikbv) = 2. * avmv(ji,jj,ikbvm1)
# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO

      CASE( 1 )                 ! linear botton friction
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
# endif
               ikbu = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
               ikbv = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
               avmu(ji,jj,ikbu) = bfri1 * fse3uw(ji,jj,ikbu)
               avmv(ji,jj,ikbv) = bfri1 * fse3vw(ji,jj,ikbv)
# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO

      CASE( 2 )                 ! quadratic botton friction
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
!CDIR NOVERRCHK
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = 2, jpim1
# endif
               ikbu   = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
               ikbv   = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
               ikbum1 = MAX( ikbu-1, 1 )
               ikbvm1 = MAX( ikbv-1, 1 )
               
               zvu  = 0.25 * (  vn(ji,jj  ,ikbum1) + vn(ji+1,jj  ,ikbum1)     &
                              + vn(ji,jj-1,ikbum1) + vn(ji+1,jj-1,ikbum1)  )
               
               zuv  = 0.25 * (  un(ji,jj  ,ikbvm1) + un(ji-1,jj  ,ikbvm1)     &
                              + un(ji,jj+1,ikbvm1) + un(ji-1,jj+1,ikbvm1)  )
               
               zecu = SQRT(  un(ji,jj,ikbum1) * un(ji,jj,ikbum1) + zvu*zvu + bfeb2  )
               zecv = SQRT(  vn(ji,jj,ikbvm1) * vn(ji,jj,ikbvm1) + zuv*zuv + bfeb2  )

               avmu(ji,jj,ikbu) = bfri2 * zecu * fse3uw(ji,jj,ikbu)
               avmv(ji,jj,ikbv) = bfri2 * zecv * fse3vw(ji,jj,ikbv)

!!DB: depth dependent bottom friction fudge for tidal run
!!  bfri2 = bfri2*(1 + F*e(-z_btm/D)^2); take D = 20m???; F=2???
!               avmu(ji,jj,ikbu) = bfri2 *(1.0+2.*exp(-(gdepw(ikbum1)/20.)**2)) &
!                                              * zecu * fse3uw(ji,jj,ikbu)
!               avmv(ji,jj,ikbv) = bfri2 *(1.0+2.*exp(-(gdepw(ikbvm1)/20.)**2)) &
!                                              * zecv * fse3vw(ji,jj,ikbv)

# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO

      CASE( 3 )                 ! free-slip boundary condition
# if defined key_vectopt_loop   &&   ! defined key_autotasking
         jj = 1
         DO ji = jpi+2, jpij-jpi-1   ! vector opt. (forced unrolling)
# else
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
# endif
               ikbu = MIN( mbathy(ji+1,jj  ), mbathy(ji,jj) )
               ikbv = MIN( mbathy(ji  ,jj+1), mbathy(ji,jj) )
               avmu(ji,jj,ikbu) = 0.e0
               avmv(ji,jj,ikbv) = 0.e0
# if ! defined key_vectopt_loop   ||   defined key_autotasking
            END DO
# endif
         END DO

      END SELECT

      ! Lateral boundary condition on (avmu,avmv)   (unchanged sign)
      ! ------------------------------===========
      CALL lbc_lnk( avmu, 'U', 1. )
      CALL lbc_lnk( avmv, 'V', 1. )

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab3d_1=avmu, clinfo1=' bfr  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
      ENDIF

   END SUBROUTINE zdf_bfr


   SUBROUTINE zdf_bfr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_bfr_init  ***
      !!                    
      !! ** Purpose :   Initialization of the bottom friction
      !!
      !! ** Method  :   Read the nammbf namelist and check their consistency
      !!      called at the first timestep (nit000)
      !!
      !! History :
      !!   9.0  !  02-06  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      NAMELIST/nambfr/ nbotfr, bfri1, bfri2, bfeb2
      !!----------------------------------------------------------------------

      ! Read Namelist nambfr : bottom momentum boundary condition
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, nambfr )


      ! Parameter control and print
      ! ---------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zdf_bfr : momentum bottom friction'
      IF(lwp) WRITE(numout,*) '~~~~~~~'
      IF(lwp) WRITE(numout,*) '          Namelist nambfr : set bottom friction parameters'

      SELECT CASE (nbotfr)

      CASE( 0 )
         IF(lwp) WRITE(numout,*) '            no-slip '

      CASE( 1 )
         IF(lwp) WRITE(numout,*) '            linear botton friction'
         IF(lwp) WRITE(numout,*) '            friction coef.   bfri1  = ', bfri1

      CASE( 2 )
         IF(lwp) WRITE(numout,*) '            quadratic botton friction'
         IF(lwp) WRITE(numout,*) '            friction coef.   bfri2  = ', bfri2
         IF(lwp) WRITE(numout,*) '            background tke   bfeb2  = ', bfeb2

      CASE( 3 )
         IF(lwp) WRITE(numout,*) '            free-slip '

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '         bad flag value for nbotfr = ', nbotfr
         nstop = nstop + 1

      END SELECT

   END SUBROUTINE zdf_bfr_init

   !!======================================================================
END MODULE zdfbfr
