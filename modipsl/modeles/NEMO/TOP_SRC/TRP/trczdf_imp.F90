MODULE trczdf_imp
   !!==============================================================================
   !!                    ***  MODULE  trczdf_imp  ***
   !! Ocean passive tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_zdf_imp  : update the tracer trend with the vertical diffusion
   !!                  using an implicit time-stepping.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
!!DB
#ifdef key_BGCM_01
!   USE bgcm_01_initrc
   USE lib_bgcm_01
#else
   USE trctrp_lec      ! passive tracers transport
#endif
   USE prtctl_trc

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_zdf_imp          ! routine called by step.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      rdttrc                     ! vertical profile of 2 x tracer time-step

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_imp.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_zdf_imp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_imp  ***
      !!
      !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
      !!      using an implicit time stepping and add it to the general trend
      !!      of the tracer equations.
      !!
      !! ** Method  :   The vertical diffusion of tracers tra is given by:
      !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(tra) )
      !!      It is thus evaluated using a backward time scheme
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend tra :
      !!          tra = tra + dz( avt dz(t) )
      !!         (tra = tra + dz( avs dz(t) ) if lk_zdfddmtrc=T)
      !!
      !! ** Action  : - Update tra with the before vertical diffusion trend
      !!              - save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  correction on tracer trend loops
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
      !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
      !!        !                     avt multiple correction
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
      INTEGER ::   ikst, ikenm2, ikstp1
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwd, zws, zwi,          &  ! ???
         zwx, zwy, zwt              ! ???
      REAL(wp) ::  ztra      ! temporary scalars

      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   &
         ztrd
      CHARACTER (len=22) :: charout
      !!---------------------------------------------------------------------

      IF(lwp .AND. kt == nittrc000 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_zdf_implicit : vertical tracer mixing'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
         ! time step = 2 rdttra with Arakawa or TVD advection scheme
         IF( neuler == 0 .AND. kt == nittrc000 ) THEN
            rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
         ELSEIF( kt <= nittrc000 + 1 ) THEN
            rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
         ENDIF
      ELSE
         rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
      ENDIF

      DO jn = 1 , jptra

	 ! Initialisation     
	 zwd( 1 ,:,:) = 0.e0     ;     zwd(jpi,:,:) = 0.e0
	 zws( 1 ,:,:) = 0.e0     ;     zws(jpi,:,:) = 0.e0
	 zwi( 1 ,:,:) = 0.e0     ;     zwi(jpi,:,:) = 0.e0
	 zwt( 1 ,:,:) = 0.e0     ;     zwt(jpi,:,:) = 0.e0     
         zwt(  :,:,1) = 0.e0     ;     zwt(  :,:,jpk) = 0.e0
         !                                          
         ! 0. Matrix construction
         ! ----------------------

         ! Diagonal, inferior, superior
         ! (including the bottom boundary condition via avs masked
         DO jk = 1, jpkm1                                                     
            DO jj = 2, jpjm1                                      
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zwi(ji,jj,jk) = - rdttrc(jk) * fstravs(ji,jj,jk  ) /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
                  zws(ji,jj,jk) = - rdttrc(jk) * fstravs(ji,jj,jk+1) /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
                  zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
               END DO
            END DO
         END DO

         ! Surface boudary conditions
         DO jj = 2, jpjm1        
            DO ji = fs_2, fs_jpim1
               zwi(ji,jj,1) = 0.e0
               zwd(ji,jj,1) = 1. - zws(ji,jj,1)
            END DO
         END DO
         
         ! Second member construction
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1     
               DO ji = fs_2, fs_jpim1
                  zwy(ji,jj,jk) = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)
               END DO
            END DO
         END DO
         
 
	! Matrix inversion from the first level
	ikst = 1

#   include "zdf.matrixsolver.vopt.h90"        
 
         
#if defined key_trc_diatrd
         ! Compute and save the vertical diffusive of tracers trends
#  if defined key_trc_ldfiso
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztra = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn) + trtrd(ji,jj,jk,ikeep(jn),6)
               END DO
            END DO
         END DO
#  else
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztra = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn)
               END DO
            END DO
         END DO
#  endif
#endif  
         ! Save the masked passive tracer after in tra
         ! (c a u t i o n:  tracer not its trend, Leap-frog scheme done
         !                  it will not be done in tranxt)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1
                  tra(ji,jj,jk,jn) = zwx(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

         IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
            ztrd(:,:,:,:) = 0.
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     ztrd(ji,jj,jk,jn) = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
                  END DO
               END DO
            END DO
         ENDIF

      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - imp')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=ztrd, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_zdf_imp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf_imp (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_zdf_imp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf_imp
#endif
   
!!==============================================================================
END MODULE trczdf_imp
