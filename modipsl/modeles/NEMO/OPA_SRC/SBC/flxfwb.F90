MODULE flxfwb
   !!======================================================================
   !!                       ***  MODULE  flxfwb  ***
   !! Ocean fluxes   : domain averaged freshwater budget
   !!======================================================================
#if ! defined key_dynspg_rl
   !!----------------------------------------------------------------------
   !!   NOT 'key_dynspg_rl'                        Free surface formulation
   !!----------------------------------------------------------------------
   !!   flx_fwb      : freshwater budget for global ocean configurations
   !!                  in free surface and forced mode
   !!   flx_fwb_init : freshwater budget for global ocean configurations
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE cpl_oce         ! coupled atmosphere/ocean
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE flxrnf          ! ocean runoffs
   USE ocesbc          ! ocean surface boudaries conditions
   USE blk_oce
   USE flxblk          ! bulk formulea
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC flx_fwb      ! routine called by step
   PUBLIC flx_fwb_init ! routine called by opa

   !! * Share module variables
   LOGICAL, PUBLIC ::    & !!: * namelist *
      ln_fwb = .TRUE.       !: Flag to activate the fwb computation
   REAL(wp), PUBLIC ::   &  !:
      a_fwb_b  ,         &  !: annual domain averaged freshwater budget
      a_fwb                 !: for 2 year before (_b) and before year.

   !! * Module variables
   REAL(wp) ::   &
      a_emp   ,  & ! domain averaged evaporation minus precipitation
      a_precip,  & ! domain averaged precipitation
      a_rnf   ,  & ! domain averaged runoff
      a_sshb  ,  & ! domain averaged sea surface heigh at nit000
      a_sshend,  & ! domain averaged sea surface heigh at nitend
      a_sal000,  & ! domain averaged ocean salinity at nit000 (before time step)
      a_salend,  & ! domain averaged ocean salinity at nitend (now time step)
      a_aminus,  & ! 
      a_aplus
      REAL(wp), DIMENSION(jpi,jpj) ::  &
         e1e2_i    ! area of the interior domain (e1t*e2t*tmask_i)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/flxfwb.F90,v 1.6 2005/12/21 10:46:40 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE flx_fwb( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flx_fwb  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!	
      !! History :
      !!   8.2  !  01-02  (E. Durand)  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk        ! dummy loop indices
      INTEGER  ::   inum = 11         ! temporary logical unit
      INTEGER  ::   ikty              ! 
      REAL(wp) ::   &
         zarea, zvol, zwei,       &  ! temporary scalars
         zsm0, zempnew                !    "         "
      !!----------------------------------------------------------------------

      ! Mean global salinity
      zsm0 = 34.72654

      ! To compute emp mean value mean emp

      IF( kt == nit000 ) THEN
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'flx_fwb : FreshWater Budget correction'
            WRITE(numout,*) '~~~~~~~'
         ENDIF

         a_emp    = 0.e0
         a_precip = 0.e0
         a_rnf    = 0.e0
         a_sshb   = 0.e0   ! averaged sea surface heigh at nit000 
         a_sal000 = 0.e0   ! averaged ocean salinity at nit000
         a_aminus = 0.e0
         a_aplus  = 0.e0
         
         e1e2_i(:,:) = e1t(:,:) * e2t(:,:) * tmask_i(:,:)
         a_sshb = SUM( e1e2_i(:,:) * sshn(:,:) )
         IF( lk_mpp )   CALL  mpp_sum( a_sshb   )   ! sum over all the global domain
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
                  a_sal000 = a_sal000 + ( sb(ji,jj,jk) - zsm0 ) * zwei
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL  mpp_sum( a_sal000 )   ! sum over the global domain

      ENDIF
      
      ! cumulate surface freshwater budget at each time-step
      ! --------------------------------------====----------
      a_emp    = SUM( e1e2_i(:,:) * emp   (:,:) )
      IF( lk_mpp )   CALL  mpp_sum( a_emp    )   ! sum over the global domain
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
      a_precip = SUM( e1e2_i(:,:) * watm  (:,:) )
      IF( lk_mpp )   CALL  mpp_sum( a_precip )   ! sum over the global domain
#endif
      a_rnf    = SUM( e1e2_i(:,:) * runoff(:,:) )
      IF( lk_mpp )   CALL  mpp_sum( a_rnf    )   ! sum over the global domain

      IF( aminus /= 0.e0 ) a_aminus = a_aminus + ( MIN( aplus, aminus ) / aminus )
      IF( aplus  /= 0.e0 ) a_aplus  = a_aplus  + ( MIN( aplus, aminus ) / aplus  )


      ! Update empold if new year start
      ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
      IF( MOD( kt, ikty ) == 0 ) THEN
         zarea    = SUM( e1e2_i(:,:)            )
         IF( lk_mpp )   CALL  mpp_sum( zarea    )   ! sum over the global domain
         a_fwb_b = a_fwb
         a_fwb   = SUM( e1e2_i(:,:) * sshn(:,:) )
         IF( lk_mpp )   CALL  mpp_sum( a_fwb    )   ! sum over the global domain

         a_fwb   = a_fwb * 1.e+3 / ( zarea * 86400. * 365. )    ! convert in Kg/m3/s = mm/s
         !                                                      !!bug 365d year 
         empold =  a_fwb                 ! current year freshwater budget correction
         !                               ! estimate from the previous year budget
      ENDIF


      IF( kt == nitend ) THEN
         zvol  = 0.e0
         zempnew = 0.e0
         ! Mean sea level at nitend
         a_sshend = SUM( e1e2_i(:,:) * sshn(:,:) )
         zarea    = SUM( e1e2_i(:,:)             )
         
         a_salend = 0.e0
         DO jk = 1, jpkm1   
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwei  = e1e2_i(ji,jj) * fse3t(ji,jj,jk) * tmask(ji,jj,jk)
                  a_salend = a_salend + ( sn(ji,jj,jk) - zsm0 ) * zwei
                  zvol  = zvol  + zwei
               END DO
            END DO
         END DO
         IF( lk_mpp ) THEN
            CALL  mpp_sum( zarea    )   ! sums over all the global domain
            CALL  mpp_sum( a_sshend )        
            CALL  mpp_sum( a_salend )        
         ENDIF

         a_aminus = a_aminus / ( nitend - nit000 + 1 )
         a_aplus  = a_aplus  / ( nitend - nit000 + 1 )

         ! Conversion in m3
         a_emp    = a_emp    * rdttra(1) * 1.e-3 
         a_precip = a_precip * rdttra(1) * 1.e-3 / rday
         a_rnf    = a_rnf    * rdttra(1) * 1.e-3
         
      ENDIF


      ! Ecriture des diagnostiques 
      ! --------------------------

      IF( kt == nitend ) THEN

         OPEN( inum, FILE='EMPave.dat' )
         WRITE(inum, "(24X,I8,2ES24.16)" ) nyear, a_fwb_b, a_fwb
         WRITE(inum,*)
         WRITE(inum,*)    'Net freshwater budget '
         WRITE(inum,9010) '  emp    = ', a_emp   , ' m3 =', a_emp   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
         WRITE(inum,9010) '  precip = ', a_precip, ' m3 =', a_precip/((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
         WRITE(inum,9010) '  a_rnf  = ', a_rnf   , ' m3 =', a_rnf   /((nitend-nit000+1)*rdttra(1)) * 1.e-6,' Sv'
         WRITE(inum,*)
         WRITE(inum,9010) '  zarea =',zarea
         WRITE(inum,9010) '  zvol  =',zvol
         WRITE(inum,*)
         WRITE(inum,*)    'Mean sea level : '
         WRITE(inum,9010) '  at nit000 = ',a_sshb           ,' m3 '
         WRITE(inum,9010) '  at nitend = ',a_sshend         ,' m3 '
         WRITE(inum,9010) '  diff      = ',(a_sshend-a_sshb),' m3 =', (a_sshend-a_sshb)/((nitend-nit000+1)*rdt) * 1.e-6,' Sv'
         WRITE(inum,9020) '  mean sea level elevation    =', a_sshend/zarea,' m'
         WRITE(inum,*)
         WRITE(inum,*)    'Anomaly of salinity content : '
         WRITE(inum,9010) '  at nit000 = ', a_sal000           ,' psu.m3 '
         WRITE(inum,9010) '  at nitend = ', a_salend           ,' psu.m3 '
         WRITE(inum,9010) '  diff      = ', a_salend - a_sal000,' psu.m3'
         WRITE(inum,*)
         WRITE(inum,*)    'Mean salinity : '
         WRITE(inum,9020) '  at nit000 = ',  a_sal000/zvol+zsm0     ,' psu '
         WRITE(inum,9020) '  at nitend = ',  a_salend/zvol+zsm0     ,' psu '
         WRITE(inum,9020) '  diff      = ', (a_salend-a_sal000)/zvol,' psu'
         WRITE(inum,9020) '  S-SLevitus= ',  a_salend          /zvol,' psu'
         WRITE(inum,*)
         WRITE(inum,*)    'Coeff : '
         WRITE(inum,9030) '  Alpha+   =  ', a_aplus
         WRITE(inum,9030) '  Alpha-   =  ', a_aminus
         WRITE(inum,*)
      ENDIF

 9006 FORMAT(1X,A,ES24.16)
 9010 FORMAT(1X,A,ES12.5,A,F10.5,A)
 9020 FORMAT(1X,A,F10.5,A)
 9030 FORMAT(1X,A,F8.2,A)

   END SUBROUTINE flx_fwb


   SUBROUTINE flx_fwb_init
      !!---------------------------------------------------------------------
      !!                ***  ROUTINE flx_fwb_init  ***
      !!
      !! ** Purpose :
      !!
      !! ** Method  :
      !!  
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      LOGICAL ::   llbon
      CHARACTER (len=32) ::   &
         clname = 'EMPave_old.dat'
      INTEGER ::   inum = 11         ! temporary logical unit
      INTEGER ::   iyear

      NAMELIST/namfwb/ ln_fwb 
      !!----------------------------------------------------------------------
         
      ! Read Namelist namfwb : freshWater Budget correction
      ! --------------------
      REWIND( numnam )
      READ  ( numnam, namfwb )            
      
      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'flx_fwb_init : FreshWater Budget correction'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '               Namelist namfwb : set fwb parameters'
         WRITE(numout,*) '                  use or not fwb correction      ln_fwb   = ', ln_fwb
      ENDIF
      ! Option consistency
#if defined key_dynspg_rl
      IF(lwp) WRITE '               Rigid-lid option, fwb correction is useless, but valid'
#endif
      IF( lk_cpl ) THEN
         IF(lwp) WRITE(numout,*) '               Coupled option, fwb correction is a flux correction ! '
         IF(lwp) WRITE(numout,*) '               ln_fwb = .FALSE. is recommanded'
      ENDIF

      !                                        ! ==============================
      IF( ln_fwb ) THEN                        !  Freshwater budget correction 
         !                                     ! ==============================
         ! Read the corrective factor on precipitations (empold)
         INQUIRE( FILE=clname, EXIST=llbon )
         IF( llbon ) THEN
            OPEN ( inum, FILE=clname)
            READ ( inum, "(24X,I8,2ES24.16)" ) iyear, a_fwb_b, a_fwb
            CLOSE( inum )
            empold = a_fwb                  ! current year freshwater budget correction
            !                               ! estimate from the previous year budget
            IF(lwp)WRITE(numout,*)
            IF(lwp)WRITE(numout,*)'flx_fwb_init : year = ',iyear  , ' freshwater budget correction = ', empold
            IF(lwp)WRITE(numout,*)'               year = ',iyear-1, ' freshwater budget read       = ', a_fwb
            IF(lwp)WRITE(numout,*)'               year = ',iyear-2, ' freshwater budget read       = ', a_fwb_b
         ELSE
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*)'flx_fwb_init : unable to read the file', clname
            nstop = nstop + 1
         ENDIF
         !                                    ! ==============================
      ELSE                                    !      NO  budget correction 
         !                                    ! ==============================
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'flx_fwb_init : NO freshwater budget correction'
         empold = 0.e0
      ENDIF

   END SUBROUTINE flx_fwb_init

#else
   !!----------------------------------------------------------------------
   !!   Default case :                       
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL, PUBLIC ::   ln_fwb = .FALSE.   !: no fwb forced
CONTAINS
   SUBROUTINE flx_fwb( kt )                ! dummy routine
      if(lwp) WRITE(numout,*) 'flx_fwb: You should not have seen this print! error?', kt
   END SUBROUTINE flx_fwb
   SUBROUTINE flx_fwb_init
   END SUBROUTINE flx_fwb_init
#endif
   !!======================================================================
END MODULE flxfwb
