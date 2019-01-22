MODULE obcfla
#if defined key_obc && defined key_dynspg_ts
   !!=================================================================================
   !!                       ***  MODULE  obcfla  ***
   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
   !!=================================================================================

   !!---------------------------------------------------------------------------------
   !!   obc_fla_ts        : call the subroutine for each open boundary
   !!   obc_fla_ts_east   : Flather on the east  open boundary velocities & ssh
   !!   obc_fla_ts_west   : Flather on the west  open boundary velocities & ssh
   !!   obc_fla_ts_north  : Flather on the north open boundary velocities & ssh
   !!   obc_fla_ts_south  : Flather on the south open boundary velocities & ssh
   !!----------------------------------------------------------------------------------

   !!----------------------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE dynspg_oce      ! surface pressure gradient variables
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE obcdta          ! ocean open boundary conditions: climatology

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_fla_ts  ! routine called in dynspg_ts (free surface time splitting case)

   !!---------------------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcfla.F90,v 1.3 2006/06/19 13:05:47 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_fla_ts
      !!------------------------------------------------------------------------------
      !!                      SUBROUTINE obc_fla_ts
      !!                     **********************
      !! ** Purpose :
      !!      Apply Flather's algorithm at open boundaries for the time-splitting
      !!      free surface case (barotropic variables)
      !!
      !!      This routine is called in dynspg_ts.F90 routine 
      !!
      !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
      !!      and/or lp_obc_south allow the user to determine which boundary is an
      !!      open one (must be done in the obc_par.F90 file).
      !!
      !! ** Reference :
      !!         Flather, R. A., 1976, Mem. Soc. R. Sci. Liege, Ser. 6, 10, 141-164
      !!
      !! History :
      !!   9.0  !  05-12  (V. Garnier) original
      !!------------------------------------------------------------------------------

      IF( lp_obc_east  )   CALL obc_fla_ts_east 
      IF( lp_obc_west  )   CALL obc_fla_ts_west 
      IF( lp_obc_north )   CALL obc_fla_ts_north
      IF( lp_obc_south )   CALL obc_fla_ts_south

   END SUBROUTINE obc_fla_ts


   SUBROUTINE obc_fla_ts_east
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_fla_ts_east  ***
      !!
      !! ** Purpose :
      !!      Apply Flather's algorithm on east OBC velocities ua, va
      !!      Fix sea surface height (sshn_e) on east open boundary
      !!
      !!  History :
      !!   9.0  !  05-12  (V. Garnier) original
      !!------------------------------------------------------------------------------
      !! * Local declaration
! FD modifications: start
      INTEGER ::   ji, jj ! dummy loop indices
      !!------------------------------------------------------------------------------

      DO ji = nie0, nie1
         DO jj = 1, jpj
            ua_e(ji,jj) = (  ubtfoe(jj) + sqrt( grav*hu(ji,jj) )	   &
               &	    * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
               &	    - sshfoe(jj) )  ) * uemsk(jj,1)
! DW Note that sshfoe_b is a velocity -- not a sea level!!
! DW It is averaged in time over the BT loop.
            sshfoe_b(ji,jj) = sshfoe_b(ji,jj) + sqrt( grav*hur(ji,jj) )     &
               &             * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
               &                 - sshfoe(jj) ) * uemsk(jj,1)
         END DO
      END DO
      DO ji = nie0p1, nie1p1
         DO jj = 1, jpj
! FD modifications: end
            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - temsk(jj,1) ) &
               &            + temsk(jj,1) * sshfoe(jj)
            va_e(ji,jj) = vbtfoe(jj) * uemsk(jj,1)
         END DO
      END DO

   END SUBROUTINE obc_fla_ts_east


   SUBROUTINE obc_fla_ts_west
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_fla_ts_west  ***
      !! 
      !! ** Purpose :
      !!      Apply Flather's algorithm on west OBC velocities ua, va
      !!      Fix sea surface height (sshn_e) on west open boundary
      !!
      !!  History :
      !!   9.0  !  05-12  (V. Garnier) original
      !!------------------------------------------------------------------------------
      !! * Local declaration
! FD modifications: start
      INTEGER ::   ji, jj, jk ! dummy loop indices
      !!------------------------------------------------------------------------------

      DO ji = niw0, niw1
         DO jj = 1, jpj
            ua_e(ji,jj) = ( ubtfow(jj) - sqrt( grav * hu(ji,jj) )	   &
               &	    * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
               &		- sshfow(jj) ) ) * uwmsk(jj,1)
            va_e(ji,jj) = vbtfow(jj) * uwmsk(jj,1)
            sshfow_b(ji,jj) = sshfow_b(ji,jj) - sqrt( grav * hur(ji,jj) )     &
                              * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5   &
                                 - sshfow(jj) ) * uwmsk(jj,1)
! FD modifications: end
            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - twmsk(jj,1) ) &
               &            + twmsk(jj,1)*sshfow(jj)
         END DO
      END DO

   END SUBROUTINE obc_fla_ts_west

   SUBROUTINE obc_fla_ts_north
      !!------------------------------------------------------------------------------
      !!                     SUBROUTINE obc_fla_ts_north
      !!                    *************************
      !! ** Purpose :
      !!      Apply Flather's algorithm on north OBC velocities ua, va
      !!      Fix sea surface height (sshn_e) on north open boundary
      !!
      !!  History :
      !!   9.0  !  05-12  (V. Garnier) original
      !!------------------------------------------------------------------------------
      !! * Local declaration
! FD modifications: start
      INTEGER ::   ji, jj ! dummy loop indices
      !!------------------------------------------------------------------------------

      DO jj = njn0, njn1
         DO ji = 1, jpi
            va_e(ji,jj) = ( vbtfon(ji) + sqrt( grav * hv(ji,jj) )	    &
               &	    * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5   &
               &		- sshfon(ji) ) ) * vnmsk(ji,1)
            sshfon_b(ji,jj) = sshfon_b(ji,jj) + sqrt( grav * hvr(ji,jj) )  &
               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
               &                  - sshfon(ji) ) * vnmsk(ji,1)
         END DO
      END DO
      DO jj = njn0p1, njn1p1
         DO ji = 1, jpi
            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - tnmsk(ji,1) ) &
               &            + sshfon(ji) * tnmsk(ji,1)
            ua_e(ji,jj) = ubtfon(ji) * vnmsk(ji,1)
         END DO
      END DO
! FD modifications: end

   END SUBROUTINE obc_fla_ts_north

   SUBROUTINE obc_fla_ts_south
      !!------------------------------------------------------------------------------
      !!                     SUBROUTINE obc_fla_ts_south
      !!                    *************************
      !! ** Purpose :
      !!      Apply Flather's algorithm on south OBC velocities ua, va
      !!      Fix sea surface height (sshn_e) on south open boundary
      !!
      !!  History :
      !!   9.0  !  05-12  (V. Garnier) original
      !!------------------------------------------------------------------------------
      !! * Local declaration
! FD modifications: start
      INTEGER ::   ji, jj ! dummy loop indices

      !!------------------------------------------------------------------------------

      DO jj = njs0, njs1
         DO ji = 1, jpi
            va_e(ji,jj) = ( vbtfos(ji) - sqrt( grav * hv(ji,jj) )            &
                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
                  &                - sshfos(ji) ) ) * vsmsk(ji,1)
            ua_e(ji,jj) = ubtfos(ji) * vsmsk(ji,1)
            sshfos_b(ji,jj) = sshfos_b(ji,jj) - sqrt( grav * hvr(ji,jj) )      &
               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
               &                  - sshfos(ji) ) * vsmsk(ji,1)
            ssha_e(ji,jj) = ssha_e(ji,jj) * (1. - tsmsk(ji,1) ) &
               &            + tsmsk(ji,1) * sshfos(ji)
         END DO
      END DO
! FD modifications: end

   END SUBROUTINE obc_fla_ts_south
#else
   !!=================================================================================
   !!                       ***  MODULE  obcfla  ***
   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
   !!=================================================================================
CONTAINS

   SUBROUTINE obc_fla_ts
      WRITE(*,*) 'obc_fla_ts: You should not have seen this print! error?'
   END SUBROUTINE obc_fla_ts
#endif

END MODULE obcfla
