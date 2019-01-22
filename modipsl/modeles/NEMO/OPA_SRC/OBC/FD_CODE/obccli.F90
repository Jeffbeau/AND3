MODULE obccli
   !!======================================================================
   !!                       ***  MODULE  obccli  ***
   !! Ocean dynamics:   Baroclinic velocities on each open boundary
   !!======================================================================
   !! History :
   !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
   !!   9.0  !  06-04 (R.Benshila, G. Madec)  zco, zps, sco coordinate
   !!----------------------------------------------------------------------
#if defined key_obc && defined key_dynspg_rl
   !!----------------------------------------------------------------------
   !!   'key_obc' and 'key_dynspg_rl' open boundary condition and rigid-lid
   !!----------------------------------------------------------------------
   !!   obc_cli_dyn : baroclinic componant after the radiation phase
   !!   obc_cli_dta : baroclinic componant for the climatological velocities
   !!----------------------------------------------------------------------
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
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obccli.F90,v 1.4 2006/05/10 17:38:46 opalod Exp $ 
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_cli_dyn( obvel, velcli, obd, obf, obtyp, obl)
      !!----------------------------------------------------------------------
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
      !!-----------------------------------------------------------------------
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
      REAL(wp), DIMENSION(obl,jpk) ::   &
         ze3ob                        ! Vertical scale factor
      !!--------------------------------------------------------------------------------

      ! 0. Array initialization
      ! -----------------------

      zhinv (:)   = 0.e0
      zmskob(:,:) = 0.e0
      zvel  (:,:) = 0.e0
      ze3ob (:,:) = 0.e0

      IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
         DO ji = obd, obf
            zhinv (:)   = hur  (ji,:)
            zmskob(:,:) = umask(ji,:,:)
            zvel  (:,:) = obvel(ji,:,:)
            ze3ob (:,:) = fse3u(ji,:,:)
         END DO
      ELSE                             ! Zonal Open Boundary ( North or South OB )
         DO jj = obd, obf
            zhinv (:)   = hvr  (:,jj)
            zmskob(:,:) = vmask(:,jj,:)
            zvel  (:,:) = obvel(:,jj,:)
            ze3ob (:,:) = fse3v(:,jj,:)
         END DO
      END IF


      ! 1. vertical sum
      ! ----------------
      zvelbtpe(1) = 0.e0
!CDIR NOLOOPCHG
      DO jk = 1, jpkm1
         DO jle = 1, obl
            zvelbtpe(jle) = zvelbtpe(jle) + zvel(jle,jk) * zmskob(jle,jk) * ze3ob(jle,jk)
         END DO
      END DO

      ! 2. divide by the depth
      ! -----------------------
      zvelbtpe(:) = zvelbtpe(:) * zhinv(:) * zmskob(:,1) 

      ! 3. substract zvelbtpe to the total velocity
      !    and save the baroclinic velocity in velcli()
      ! ------------------------------------------------
      DO jk = 1, jpkm1
         velcli(:,jk) = ( zvel(:,jk) - zvelbtpe(:) ) * zmskob(:,jk)
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
      REAL(wp), DIMENSION(obl,jpk) ::   &
         ze3ob                        ! Vertical scale factor
      !!--------------------------------------------------------------------------------

      ! 0. Array initialization
      ! -----------------------

      zhinv (:)   = 0.e0
      zmskob(:,:) = 0.e0
      ze3ob (:,:) = 0.e0

      IF( obtyp == 0 ) THEN            ! Meridional Open Boundary ( East or West OB )
         DO ji = obd, obf
            zhinv (:)   = hur  (ji,:)
            zmskob(:,:) = umask(ji,:,:)
            ze3ob (:,:) = fse3u(ji,:,:)
         END DO
      ELSE                             ! Zonal Open Boundary ( North or South OB )
         DO jj = obd, obf
            zhinv (:)   = hvr  (:,jj)
            zmskob(:,:) = vmask(:,jj,:)
            ze3ob (:,:) = fse3v(:,jj,:)
         END DO
      END IF

      ! 1. vertical sum
      ! ----------------
      zvelbtpe(1) = 0.e0
!CDIR NOLOOPCHG
      DO jk = 1, jpkm1
         DO jle = 1, obl
            ij = jle -1 + mpp
            zvelbtpe(jle) = zvelbtpe(jle) + obvel(ij,jk,1)*zmskob(jle,jk) * ze3ob(jle,jk)
         END DO
      END DO

      ! 2. divide by the depth
      ! -----------------------
         zvelbtpe(:) = zvelbtpe(:) * zhinv(:) * zmskob(:,1) 

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
