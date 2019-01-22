MODULE obccli
   !!===================================================================================
   !!                       ***  MODULE  obccli  ***
   !! Ocean dynamics:   Baroclinic componant of velocities on each open boundary
   !!===================================================================================
#if defined key_obc && defined key_dynspg_rl
   !!-----------------------------------------------------------------------------------
   !!   'key_obc'               and 
   !!   'key_dynspg_rl'
   !!-----------------------------------------------------------------------------------
   !!   obc_cli_dyn : Compute the baroclinic componant after the radiation phase
   !!   obc_cli_dta : Compute the baroclinic componant for the climatological velocities
   !!-----------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_cli    ! routine called in obcdyn.F90 and obcdta.F90 (rigid lid case)

   INTERFACE obc_cli
     MODULE PROCEDURE obc_cli_dyn, obc_cli_dta
   END INTERFACE

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!-----------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obccli.F90,v 1.3 2005/12/28 09:25:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!-----------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_cli_dyn( obvel, velcli, obd, obf, obtyp, obl)
      !!--------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_cli_dyn  ***
      !!                   
      !! ** Purpose :   Compute the baroclinic velocities at the open boundaries.
      !!
      !! ** Method  :
      !!      - Compute de barotropic velocity along the considered Open Boundary 
      !!        and substract it to the total velocity to have baroclinic velotity.
      !!      - obtyp must be set to | 0 when traiting an East or West OB 
      !!                             | 1 when traiting a North or South OB.
      !!      - obl is the lenght of the OB (jpi or jpj) 
      !!
      !! History :
      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!--------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   & ! OB localization:jpieob or jpiwob for East or West 
         obd, obf,               & !                 jpjnob or jpjsob for North or South 
         obl,                    & ! Lenght of the Open Boundary
         obtyp                     ! Type of Open Boundary: zonal or Meridional 
      REAL(wp), DIMENSION(:,:), INTENT( out) ::   &
         velcli                    ! Baroclinic velocity calculated
      REAL(wp), DIMENSION(:,:,:), INTENT( in ) ::   &
         obvel                     ! ua or va velocities from obcdyn.F90 routine

      !! * Local declarations
      INTEGER ::   &   
         ji, jj, jk, jle, jol         ! loop indices  
      REAL(wp) ::   zcbl              ! Temporary Baroclinic velocity 
      REAL(wp), DIMENSION(obl) ::   & 
         zvelbtpe,                  & ! Barotropic velocity 
         zhinv                        ! Invert of the local depth 1/H
      REAL(wp), DIMENSION(obl,jpk) ::   &
         zmskob,                    & ! Velocity mask
         zvel                         ! 2D Local velocity on OB
# if defined key_partial_steps
      REAL(wp), DIMENSION(obl,jpk) ::   &
         ze3ob                        ! Vertical scale factor
# else
      REAL(wp), DIMENSION(jpk) ::   &
         ze3ob                        ! Vertical scale factor
# endif
      !!--------------------------------------------------------------------------------

      ! 0. Array initialization
      ! -----------------------

      zhinv(:) = 0.e0
      zmskob(:,:) = 0.e0
      zvel(:,:) = 0.e0
# if defined key_partial_steps
      ze3ob(:,:) = 0.e0
# else
      ze3ob(:) = 0.e0
# endif

      IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
         DO ji = obd, obf
            zhinv(:) = hur(ji,:)
            zmskob(:,:) = umask(ji,:,:)
            zvel(:,:) = obvel(ji,:,:)
# if defined key_partial_steps
            ze3ob(:,:) = fse3u(ji,:,:)
# else
            ze3ob(:) = fse3u(:,:,:)
# endif
         END DO
      ELSE                             ! Zonal Open Boundary ( North or South OB )
         DO jj = obd, obf
            zhinv(:) = hvr(:,jj)
            zmskob(:,:) = vmask(:,jj,:)
            zvel(:,:) = obvel(:,jj,:)
# if defined key_partial_steps
            ze3ob(:,:) = fse3v(:,jj,:)
# else
            ze3ob(:) = fse3v(:,:,:)
# endif
         END DO
      END IF

      zvelbtpe(:) = 0.e0

      ! 1. vertical sum
      ! ----------------
# if defined key_vectopt_loop
!CDIR NOLOOPCHG
# endif
      DO jol = obd, obf ! Vector opt.
         DO jk = 1, jpkm1
            DO jle = 1, obl
               zvelbtpe(jle) = zvelbtpe(jle) + zvel(jle,jk)*zmskob(jle,jk) &
# if defined key_partial_steps
                                           * ze3ob(jol,jle,jk)
# else
                                           * ze3ob(jk)
# endif
            END DO
         END DO
      END DO

      ! 2. divide by the depth
      ! -----------------------
      DO jle = 1, obl
         zvelbtpe(jle) = zvelbtpe(jle) * zhinv(jle) * zmskob(jle,1) 
      END DO

      ! 3. substract zvelbtpe to the total velocity
      !    and save the baroclinic velocity in velcli()
      ! ------------------------------------------------
      DO jk = 1, jpkm1
         DO jle = 1, obl
            zcbl = zvel(jle,jk) - zvelbtpe(jle)*zmskob(jle,jk)
            velcli(jle,jk) = zcbl * zmskob(jle,jk)
         END DO
      END DO

   END SUBROUTINE obc_cli_dyn


   SUBROUTINE obc_cli_dta( obvel, velcli, obd, obf, obtyp, obl, mpp )
      !!--------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_cli_dta  ***
      !!                   
      !! ** Purpose :
      !!      Compute the baroclinic velocities for the climatological velocities.
      !!
      !! ** Method  :
      !!      - Compute de barotropic velocity along the considered Open Boundary 
      !!        and substract it to the total velocity to have baroclinic velotity.
      !!      - obtyp must be set to | 0 when traiting an East or West OB 
      !!                             | 1 when traiting a North or South OB.
      !!      - obl is the lenght of the OB (jpi or jpj) 
      !!
      !! History :
      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
      !!--------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   & ! OB localization: jpieob or jpiwob for East or West 
         obd, obf,               & !                  jpjnob or jpjsob for North or South 
         obl,                    & ! Lenght of the Open Boundary
         mpp,                    & ! MPP index
         obtyp                     ! Type of Open Boundary: zonal or Meridional 
      REAL(wp), INTENT( out), DIMENSION(:,:) ::   &
         velcli                    ! Baroclinic velocity calculated
      REAL(wp), INTENT( inout ), DIMENSION(:,:,:) ::   &
         obvel                     ! uXdta or vXdta climatological velocities from 
                                   ! obcdta.F90 routine

      !! * Local declarations
      INTEGER ::   &
         ji, jj, jk, jle, jol, ij     ! loop indices  
      REAL(wp), DIMENSION(obl) ::   & 
         zvelbtpe,                  & ! Barotropic velocity 
         zhinv                        ! Invert of the local depth 1/H
      REAL(wp), DIMENSION(obl,jpk) ::   &
         zmskob                       ! Velocity mask
# if defined key_partial_steps
      REAL(wp), DIMENSION(obl,jpk) ::   &
         ze3ob                        ! Vertical scale factor
# else
      REAL(wp), DIMENSION(jpk) ::   &
         ze3ob                        ! Vertical scale factor
# endif
      !!--------------------------------------------------------------------------------

      ! 0. Array initialization
      ! -----------------------

      zhinv(:) = 0.e0
      zmskob(:,:) = 0.e0
# if defined key_partial_steps
      ze3ob(:,:) = 0.e0
# else
      ze3ob(:) = 0.e0
# endif

      IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
         DO ji = obd, obf
            zhinv(:) = hur(ji,:)
            zmskob(:,:) = umask(ji,:,:)
# if defined key_partial_steps
            ze3ob(:,:) = fse3u(ji,:,:)
# else
            ze3ob(:) = fse3u(:,:,:)
# endif
         END DO
      ELSE                             ! Zonal Open Boundary ( North or South OB )
         DO jj = obd, obf
            zhinv(:) = hvr(:,jj)
            zmskob(:,:) = vmask(:,jj,:)
# if defined key_partial_steps
            ze3ob(:,:) = fse3v(:,jj,:)
# else
            ze3ob(:) = fse3v(:,:,:)
# endif
         END DO
      END IF

      zvelbtpe(:) = 0.e0

      ! 1. vertical sum
      ! ----------------
# if defined key_vectopt_loop
!CDIR NOLOOPCHG
# endif
      DO jol = obd, obf ! Vector opt.
         DO jk = 1, jpkm1
            DO jle = 1, obl
               ij = jle -1 + mpp
               zvelbtpe(jle) = zvelbtpe(jle) + obvel(ij,jk,1)*zmskob(jle,jk) &
# if defined key_partial_steps
                                           * ze3ob(jol,jle,jk)
# else
                                           * ze3ob(jk)
# endif
            END DO
         END DO
      END DO

      ! 2. divide by the depth
      ! -----------------------
      DO jle = 1, obl
         zvelbtpe(jle) = zvelbtpe(jle) * zhinv(jle) * zmskob(jle,1) 
      END DO 

      ! 3. substract zvelbtpe to the total velocity
      !    and save the baroclinic velocity in velcli()
      ! ------------------------------------------------
      DO jk = 1, jpkm1
         DO jle = 1, obl
            ij = jle -1 + mpp
            obvel(ij,jk,1) = obvel(ij,jk,1) - zvelbtpe(jle)*zmskob(jle,jk)
            velcli(jle,jk) = obvel(ij,jk,1) * zmskob(jle,jk)
         END DO
      END DO

   END SUBROUTINE obc_cli_dta

#else
   !!----------------------------------------------------------------------------------
   !!   Default options :                                                  Empty module
   !!----------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_cli_dyn       ! Empty routine
   END SUBROUTINE obc_cli_dyn
   SUBROUTINE obc_cli_dta       ! Empty routine
   END SUBROUTINE obc_cli_dta
#endif

   !!==================================================================================
END MODULE obccli
