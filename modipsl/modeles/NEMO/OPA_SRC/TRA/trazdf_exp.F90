MODULE trazdf_exp
   !!==============================================================================
   !!                    ***  MODULE  trazdf_exp  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend using
   !!                        an explicit time-stepping (time spllitting scheme)
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_zdf_exp  : update the tracer trend with the vertical diffusion
   !!                  using an explicit time stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE zdf_oce         ! ocean vertical physics
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_zdf_exp          ! routine called by step.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      r2dt                     ! vertical profile of 2 x tracer time-step

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_exp.F90,v 1.4 2005/09/02 15:45:34 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_zdf_exp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
      !!      using an explicit time stepping and add it to the general trend 
      !!      of the tracer equations.
      !!
      !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
      !!         difft = dz( avt dz(tb) ) = 1/e3t dk+1( avt/e3w dk(tb) )
      !!      It is evaluated with an Euler scheme, using a time splitting
      !!      technique.
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!          ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm= T)
      !!
      !! ** Action : - Update (ta,sa) with the before vertical diffusion trend
      !!             - Save the trends  in (ztdta,ztdsa) ('key_trdtra')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  correction on tracer trend loops
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!---------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
                         ztdsa => va          ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jl             ! dummy loop indices
      REAL(wp) ::   &
         zlavmr,                            & ! temporary scalars
         zave3r, ze3tr,                     & !    "         "
         zta, zsa                             !    "         " 
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwx, zwy, zwz, zww
      !!---------------------------------------------------------------------


      ! 0. Local constant initialization
      ! --------------------------------
      ! time step = 2 rdttra 
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
      ELSEIF( kt <= nit000 + 1) THEN
         r2dt(:) = 2. * rdttra(:)          ! leapfrog
      ENDIF
      zlavmr = 1. / float( n_zdfexp )

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 1. Initializations
         ! ------------------

         ! Surface & bottom boundary conditions: no flux
         DO ji = 2, jpim1
            zwy(ji, 1 ) = 0.e0
            zwy(ji,jpk) = 0.e0
            zww(ji, 1 ) = 0.e0
            zww(ji,jpk) = 0.e0
         END DO

         ! zwx and zwz arrays set to before tracer values
         DO jk = 1, jpk
            DO ji = 2, jpim1
               zwx(ji,jk) = tb(ji,jj,jk)
               zwz(ji,jk) = sb(ji,jj,jk)
            END DO
         END DO


         ! 2. Time splitting loop
         ! ----------------------

         DO jl = 1, n_zdfexp

            ! first vertical derivative
            IF( lk_zdfddm ) THEN       ! double diffusion: avs /= avt
               DO jk = 2, jpk
                  DO ji = 2, jpim1
                     zave3r = 1.e0 / fse3w(ji,jj,jk) 
                     zwy(ji,jk) =   avt(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) * zave3r
                     zww(ji,jk) = fsavs(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) * zave3r
                  END DO
               END DO
            ELSE                      ! default : avs = avt
               DO jk = 2, jpk
                  DO ji = 2, jpim1
                     zave3r = avt(ji,jj,jk) / fse3w(ji,jj,jk)
                     zwy(ji,jk) = zave3r *(zwx(ji,jk-1) - zwx(ji,jk) )
                     zww(ji,jk) = zave3r *(zwz(ji,jk-1) - zwz(ji,jk) )
                  END DO
               END DO
            ENDIF

            ! trend estimation at kt+l*2*rdt/n_zdfexp
            DO jk = 1, jpkm1
               DO ji = 2, jpim1
                  ze3tr = zlavmr / fse3t(ji,jj,jk)
                  ! 2nd vertical derivative
                  zta = ( zwy(ji,jk) - zwy(ji,jk+1) ) * ze3tr
                  zsa = ( zww(ji,jk) - zww(ji,jk+1) ) * ze3tr
                  ! update the tracer trends
                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
                  ! update tracer fields at kt+l*2*rdt/n_zdfexp
                  zwx(ji,jk) = zwx(ji,jk) + r2dt(jk) * zta * tmask(ji,jj,jk)
                  zwz(ji,jk) = zwz(ji,jk) + r2dt(jk) * zsa * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      ! save the trends for diagnostic
      ! vertical diffusive tracers trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
         CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

   END SUBROUTINE tra_zdf_exp

   !!==============================================================================
END MODULE trazdf_exp
