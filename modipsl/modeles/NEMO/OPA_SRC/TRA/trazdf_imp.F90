MODULE trazdf_imp
   !!==============================================================================
   !!                    ***  MODULE  trazdf_imp  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_zdf_imp  : update the tracer trend with the vertical diffusion
   !!                  using an implicit time-stepping.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE zdf_oce         ! ocean vertical physics
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE zdfkpp          ! KPP parameterisation
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_zdf_imp          ! routine called by step.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      r2dt                     ! vertical profile of 2 x tracer time-step

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "zdfddm_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_imp.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_zdf_imp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_imp  ***
      !!
      !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
      !!      using an implicit time stepping and add it to the general trend
      !!      of the tracer equations.
      !!
      !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
      !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(ta) )
      !!      It is thus evaluated using a backward time scheme
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!          ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T)
      !!
      !! ** Action  : - Update (ta,sa) with the before vertical diffusion trend
      !!              - save the trends in (ttrd,strd) ('key_trdtra')
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
      !!   9.0  !  05-01  (C. Ethe )  non-local flux in KPP vertical mixing scheme
      !!---------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
                         ztdsa => va          ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk                 ! dummy loop indices
      INTEGER ::   ikst, ikenm2, ikstp1
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwd, zws, zwi,          &  ! ???
         zwx, zwy, zwz, zwt         ! ???
      !!---------------------------------------------------------------------


      ! 0. Local constant initialization
      ! --------------------------------

      ! time step = 2 rdttra ex
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
      ELSEIF( kt <= nit000 + 1) THEN
         r2dt(:) = 2. * rdttra(:)          ! leapfrog
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 0. Matrix construction
         ! ----------------------

         ! Diagonal, inferior, superior (including the bottom boundary condition via avt masked)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwi(ji,jk) = - r2dt(jk) * avt(ji,jj,jk  )   &
                                       / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jk) = - r2dt(jk) * avt(ji,jj,jk+1)   &
                                       / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
            END DO
         END DO

         ! Surface boudary conditions
         DO ji = 2, jpim1
            zwi(ji,1) = 0.e0
            zwd(ji,1) = 1. - zws(ji,1)
         END DO

         ! 1. Vertical diffusion on temperature
         ! -------------------------===========

         ! Second member construction
#if defined key_zdfkpp
         ! add non-local temperature flux ( in convective case only)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1  
               zwy(ji,jk) = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)  &
                  &  - r2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
                  &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
            END DO
         END DO
#else
         DO jk = 1, jpkm1
            DO ji = 2, jpim1             
               zwy(ji,jk) = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)
            END DO
         END DO
#endif

         ! Matrix inversion from the first level
         ikst = 1

#   include "zdf.matrixsolver.h90"

         ! Save the masked temperature after in ta
         ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
         !                  it will not be done in tranxt)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
            END DO
         END DO


         ! 2. Vertical diffusion on salinity
         ! -------------------------========

#if defined key_zdfddm
         ! Rebuild the Matrix as avt /= avs

         ! Diagonal, inferior, superior
         ! (including the bottom boundary condition via avs masked)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zwi(ji,jk) = - r2dt(jk) * fsavs(ji,jj,jk  )   &
                  /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jk) = - r2dt(jk) * fsavs(ji,jj,jk+1)   &
                  /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
            END DO
         END DO

         ! Surface boudary conditions
         DO ji = 2, jpim1
            zwi(ji,1) = 0.e0
            zwd(ji,1) = 1. - zws(ji,1)
         END DO
#endif
         ! Second member construction
#if defined key_zdfkpp
         ! add non-local salinity flux ( in convective case only)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1  
               zwy(ji,jk) = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)  &
                  &  - r2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
                  &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
            END DO
         END DO
#else
         DO jk = 1, jpkm1
            DO ji = 2, jpim1             
               zwy(ji,jk) = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)
            END DO
         END DO
#endif
 
         ! Matrix inversion from the first level
         ikst = 1

#   include "zdf.matrixsolver.h90"


         ! Save the masked salinity after in sa
         ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
         !                  it will not be done in tranxt)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the trends for diagnostic
      ! Compute and save the vertical diffusive temperature & salinity trends
      IF( l_trdtra )   THEN
         ! compute the vertical diffusive trends in substracting the previous 
         ! trends ztdta()/ztdsa() to the new one computed (dT/dt or dS/dt) 
         ! with the new temperature/salinity ta/sa
         DO jk = 1, jpkm1
            ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) )   & ! new trend
                &           - ztdta(:,:,jk)                                ! old trend
            ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) )   & ! new trend
                &           - ztdsa(:,:,jk)                                ! old trend
         END DO

         CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

   END SUBROUTINE tra_zdf_imp

   !!==============================================================================
END MODULE trazdf_imp
