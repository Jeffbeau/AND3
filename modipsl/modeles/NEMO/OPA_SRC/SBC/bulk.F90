MODULE bulk
   !!======================================================================
   !!                           ***  bulk  ***
   !!======================================================================
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   !!----------------------------------------------------------------------
   !!   'key_flx_bulk_monthly'                        monthly bulk formulea
   !!   'key_flx_bulk_daily'                          daily bulk formulea
   !!----------------------------------------------------------------------
   !!   bulk          : computation of fluxes using bulk formulation
   !!----------------------------------------------------------------------
   !! * Modules used   
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ice_oce         ! bulk variable  
   USE ocfzpt          ! ocean freezing point
   USE flxblk          ! bulk formulae
   USE blk_oce         ! bulk variable 
   USE flx_oce
   USE taumod
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC blk        ! called by flx.F90   
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/bulk.F90,v 1.9 2005/09/22 10:58:15 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE blk( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE blk  ***
      !!        
      !! ** Purpose :   provide the heat fluxes on ice and ocean 
      !!                using bulk formulation
      !!
      !! History :
      !!   9.0  !  03-11  (C. Ethe and G. Madec)  F90: Free form and MODULE 
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * Local declarations    
      REAL(wp), DIMENSION(jpi,jpj) ::   zsst 
# if ! defined key_ice_lim
      INTEGER  ::   ji, jj         ! dummy loop indices  
      REAL(wp) ::   ztgel, zicopa
# endif
      !!---------------------------------------------------------------------

     ! Initialisation
     IF( kt == nit000) THEN
      ! computation of rdtbs2
        IF( nacc == 1 ) THEN
           rdtbs2 = nfbulk * rdtmin * 0.5
        ELSE
           rdtbs2 = nfbulk * rdt * 0.5
        ENDIF
        IF ( .NOT.ln_rstart ) THEN
           gsst(:,:) =  ( nfbulk - 1 ) * ( tn(:,:,1) + rt0 )
        ENDIF
     ENDIF

# if ! defined key_ice_lim
      ! opa model ice freeze()      
      DO jj = 1, jpj
         DO ji = 1, jpi
            ztgel  = fzptn(ji,jj)
            zicopa = tmask(ji,jj,1)
            IF( tn(ji,jj,1) >= ztgel ) zicopa = 0.
            freeze(ji,jj) = zicopa
         END DO
      END DO
# endif

      gsst(:,:) = gsst(:,:) + tn(:,:,1) + rt0  

      !  Computation of the fluxes       
      IF( MOD( kt - 1 , nfbulk ) == 0 ) THEN

         zsst(:,:) = gsst(:,:) / REAL( nfbulk ) * tmask(:,:,1)
         CALL flx_blk( zsst )    
         gsst(:,:) = 0.    

# if ! defined key_ice_lim
         IF(ln_ctl) THEN         ! print mean trends (used for debugging)
            CALL prt_ctl_info(' Forcings ')
            CALL prt_ctl(tab2d_1=qsr_oce , clinfo1=' qsr_oce   : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=qsr_ice , clinfo1=' qsr_ice   : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=qnsr_oce, clinfo1=' qnsr_oce  : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=qnsr_ice, clinfo1=' qnsr_ice  : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=evap    , clinfo1=' evap      : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=tprecip , clinfo1=' precip    : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=sprecip , clinfo1=' Snow      : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=taux    , clinfo1=' u-stress  : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=tauy    , clinfo1=' v-stress  : ', mask1=tmask, ovlap=1)
            CALL prt_ctl(tab2d_1=zsst    , clinfo1=' sst       : ', mask1=tmask, ovlap=1)
         ENDIF
# endif   
      ENDIF
 
   END SUBROUTINE blk

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                                     NO bulk formulea
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE blk( kt )          ! Dummy routine
      if(lwp) WRITE(numout,*) 'blk: You should not see this print! error? ', kt
   END SUBROUTINE blk
#endif
 
   !!======================================================================
END MODULE bulk
