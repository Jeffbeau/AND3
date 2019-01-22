MODULE dynzdf_iso
   !!==============================================================================
   !!                       ***  MODULE  dynzdf_iso  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!==============================================================================
#if   defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                          rotation of the mixing tensor
   !!----------------------------------------------------------------------
   !!   dyn_zdf_iso  : update the momentum trend with the vertical diffusion
   !!                  (vertical mixing + vertical component of lateral
   !!                  mixing) (rotated lateral operator case)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE taumod          ! surface ocean stress
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_zdf_iso    ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_zdf_iso( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_iso  ***
      !!                   
      !! ** Purpose :
      !!         Compute the vertical momentum trend due to both vertical and
      !!      lateral mixing (only for second order lateral operator, for
      !!      fourth order it is already computed and add to the general trend
      !!      in dynldf.F) and the surface forcing, and add it to the general
      !!      trend of the momentum equations.
      !!
      !! ** Method :
      !!         The vertical component of the lateral diffusive trends is
      !!      provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      First part: vertical trends associated with the lateral mixing
      !!      ==========  (excluding the vertical flux proportional to dk[U] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zfuw =-ahm {  e2t*mi(wslpi) di[ mi(mk(ub)) ]
      !!                     + e1t*mj(wslpj) dj[ mj(mk(ub)) ]  }
      !!      update and save in zavt the vertical eddy viscosity coefficient:
      !!         avmu = avmu + mi(wslpi)^2 + mj(wslj)^2
      !!      take the horizontal divergence of the fluxes:
      !!         diffu = 1/(e1u*e2u*e3u) dk[ zfuw ]
      !!      Add this trend to the general trend (ta,sa):
      !!         ua = ua + difft
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[U]
      !!                    associated with the lateral mixing, through the
      !!                    update of avmu)
      !!      The vertical diffusion of momentum is given by:
      !!             diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
      !!      using a backward (implicit) time stepping.
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !!      'key_trddyn' defined: trend saved for further diagnostics.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
      !!               mixing trend.
      !!             - Save the trends in (ztdua,ztdva) ('key_trddyn')
      !!
      !! History :
      !!        !  90-10  (B. Blanke)  Original code
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier)  New trends organization
      !!---------------------------------------------------------------------
      !! * Modules used
      USE ldfslp    , ONLY : wslpi, wslpj
      USE ldftra_oce, ONLY : aht0
      USE oce, ONLY :    ztdua => ta,        & ! use ta as 3D workspace   
                         ztdva => sa           ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt            ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk                  ! dummy loop indices
      INTEGER ::   &
         ikst, ikenm2, ikstp1,               & ! temporary integers
         ikbu, ikbum1 , ikbv, ikbvm1           !    "        "      
      REAL(wp) ::   &
         zrau0r, z2dt,                       & ! temporary scalars
         z2dtf, zcoef, zzws
      REAL(wp) ::   &
         zcoef0, zcoef3, zcoef4, zbu, zbv, zmkt, zmkf,   &
         zuav, zvav, zuwslpi, zuwslpj, zvwslpi, zvwslpj
      REAL(wp), DIMENSION(jpi,jpk) ::        &
         zwx, zwy, zwz,                      & ! workspace arrays
         zwd, zws, zwi, zwt,                 & !    "        "
         zfuw, zdiu, zdju, zdj1u,            & !    "        "
         zfvw, zdiv, zdjv, zdj1v
      REAL(wp), DIMENSION(jpi,jpj) ::        &
         ztsx, ztsy, ztbx, ztby                ! temporary workspace arrays
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_iso : vertical momentum diffusion isopycnal operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      
      ! inverse of the reference density
      zrau0r = 1. / rau0
      ! Leap-frog environnement
      z2dt = 2. * rdt
      ! workspace arrays
      ztsx(:,:)   = 0.e0
      ztsy(:,:)   = 0.e0 
      ztbx(:,:)   = 0.e0
      ztby(:,:)   = 0.e0
      ! Euler time stepping when starting from rest
      IF ( neuler == 0 .AND. kt == nit000 ) z2dt = rdt

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============

  
         ! I. vertical trends associated with the lateral mixing
         ! =====================================================
         !  (excluding the vertical flux proportional to dk[t]


         ! I.1 horizontal momentum gradient
         ! --------------------------------

         DO jk = 1, jpk
            DO ji = 2, jpi
               ! i-gradient of u at jj
               zdiu (ji,jk) = tmask(ji,jj  ,jk) * ( ub(ji,jj  ,jk) - ub(ji-1,jj  ,jk) )
               ! j-gradient of u and v at jj
               zdju (ji,jk) = fmask(ji,jj  ,jk) * ( ub(ji,jj+1,jk) - ub(ji  ,jj  ,jk) )
               zdjv (ji,jk) = tmask(ji,jj  ,jk) * ( vb(ji,jj  ,jk) - vb(ji  ,jj-1,jk) )
               ! j-gradient of u and v at jj+1
               zdj1u(ji,jk) = fmask(ji,jj-1,jk) * ( ub(ji,jj  ,jk) - ub(ji  ,jj-1,jk) )
               zdj1v(ji,jk) = tmask(ji,jj+1,jk) * ( vb(ji,jj+1,jk) - vb(ji  ,jj  ,jk) )
            END DO
         END DO
         DO jk = 1, jpk
            DO ji = 1, jpim1
               ! i-gradient of v at jj
               zdiv (ji,jk) = fmask(ji,jj  ,jk) * ( vb(ji+1,jj,jk) - vb(ji  ,jj  ,jk) )
            END DO
         END DO


         ! I.2 Vertical fluxes
         ! -------------------

         ! Surface and bottom vertical fluxes set to zero
         DO ji = 1, jpi
            zfuw(ji, 1 ) = 0.e0
            zfvw(ji, 1 ) = 0.e0
            zfuw(ji,jpk) = 0.e0
            zfvw(ji,jpk) = 0.e0
         END DO

         ! interior (2=<jk=<jpk-1) on U field
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               zcoef0= 0.5 * aht0 * umask(ji,jj,jk)

               zuwslpi = zcoef0 * ( wslpi(ji+1,jj,jk) + wslpi(ji,jj,jk) )
               zuwslpj = zcoef0 * ( wslpj(ji+1,jj,jk) + wslpj(ji,jj,jk) )

               zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji+1,jj,jk-1)   &
                             + tmask(ji,jj,jk  )+tmask(ji+1,jj,jk  ), 1. )
               zmkf = 1./MAX(  fmask(ji,jj-1,jk-1)+fmask(ji,jj,jk-1)   &
                             + fmask(ji,jj-1,jk  )+fmask(ji,jj,jk  ), 1. )

               zcoef3 = - e2u(ji,jj) * zmkt * zuwslpi
               zcoef4 = - e1u(ji,jj) * zmkf * zuwslpj
               ! vertical flux on u field
               zfuw(ji,jk) = zcoef3 * ( zdiu (ji,jk-1) + zdiu (ji+1,jk-1)     &
                                       +zdiu (ji,jk  ) + zdiu (ji+1,jk  ) )   &
                           + zcoef4 * ( zdj1u(ji,jk-1) + zdju (ji  ,jk-1)     &
                                       +zdj1u(ji,jk  ) + zdju (ji  ,jk  ) )
               ! update avmu (add isopycnal vertical coefficient to avmu)
               avmu(ji,jj,jk) = avmu(ji,jj,jk) + ( zuwslpi * zuwslpi          &
                                                 + zuwslpj * zuwslpj ) / aht0
            END DO
         END DO

         ! interior (2=<jk=<jpk-1) on V field
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               zcoef0= 0.5 * aht0 * vmask(ji,jj,jk)

               zvwslpi = zcoef0 * ( wslpi(ji,jj+1,jk) + wslpi(ji,jj,jk) )
               zvwslpj = zcoef0 * ( wslpj(ji,jj+1,jk) + wslpj(ji,jj,jk) )

               zmkf = 1./MAX(  fmask(ji-1,jj,jk-1)+fmask(ji,jj,jk-1)   &
                             + fmask(ji-1,jj,jk  )+fmask(ji,jj,jk  ), 1. )
               zmkt = 1./MAX(  tmask(ji,jj,jk-1)+tmask(ji,jj+1,jk-1)   &
                             + tmask(ji,jj,jk  )+tmask(ji,jj+1,jk  ), 1. )

               zcoef3 = - e2v(ji,jj) * zmkf * zvwslpi
               zcoef4 = - e1v(ji,jj) * zmkt * zvwslpj
               ! vertical flux on v field
               zfvw(ji,jk) = zcoef3 * ( zdiv (ji,jk-1) + zdiv (ji-1,jk-1)     &
                                       +zdiv (ji,jk  ) + zdiv (ji-1,jk  ) )   &
                           + zcoef4 * ( zdjv (ji,jk-1) + zdj1v(ji  ,jk-1)     &
                                       +zdjv (ji,jk  ) + zdj1v(ji  ,jk  ) )
               ! update avmv (add isopycnal vertical coefficient to avmv)
               avmv(ji,jj,jk) = avmv(ji,jj,jk) + ( zvwslpi * zvwslpi          &
                                                 + zvwslpj * zvwslpj ) / aht0
            END DO
         END DO


         ! I.3 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------

         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ! volume elements
               zbu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
               zbv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
               ! part of the k-component of isopycnal momentum diffusive trends
               zuav = ( zfuw(ji,jk) - zfuw(ji,jk+1) ) / zbu
               zvav = ( zfvw(ji,jk) - zfvw(ji,jk+1) ) / zbv
               ! add the trends to the general trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zuav
               va(ji,jj,jk) = va(ji,jj,jk) + zvav
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      IF( l_trddyn )   THEN
         ! save these trends in addition to the lateral diffusion one for diagnostics
         uldftrd(:,:,:) = uldftrd(:,:,:) + ua(:,:,:) - ztdua(:,:,:)
         vldftrd(:,:,:) = vldftrd(:,:,:) + va(:,:,:) - ztdva(:,:,:)

         ! save new trends ua and va
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 1. Vertical diffusion on u
         ! ---------------------------

         ! 1.0 Matrix and second member construction
         ! bottom boundary condition: only zws must be masked as avmu can take
         ! non zero value at the ocean bottom depending on the bottom friction
         ! used (see zdfmix.F)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zcoef = - z2dt / fse3u(ji,jj,jk)
               zwi(ji,jk) = zcoef * avmu(ji,jj,jk  ) / fse3uw(ji,jj,jk  )
               zzws       = zcoef * avmu(ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
               zws(ji,jk) = zzws * umask(ji,jj,jk+1)
               zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
               zwy(ji,jk) = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)
            END DO  
         END DO  

         ! 1.1 Surface boudary conditions
         DO ji = 2, jpim1
            z2dtf = z2dt / ( fse3u(ji,jj,1)*rau0 )
            zwi(ji,1) = 0.
            zwd(ji,1) = 1. - zws(ji,1)
            zwy(ji,1) = zwy(ji,1) + z2dtf * taux(ji,jj)
         END DO  

         ! 1.2 Matrix inversion starting from the first level
         ikst = 1
#include "zdf.matrixsolver.h90"

         ! 1.3 Normalization to obtain the general momentum trend ua
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ( zwx(ji,jk) - ub(ji,jj,jk) ) / z2dt
            END DO  
         END DO  

         ! 1.4 diagnose surface and bottom momentum fluxes
         DO ji = 2, jpim1
            ! save the surface forcing momentum fluxes
            ztsx(ji,jj) = taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 )
            ! save bottom friction momentum fluxes
            ikbu   = min( mbathy(ji+1,jj), mbathy(ji,jj) )
            ikbum1 = max( ikbu-1, 1 )
            ztbx(ji,jj) = - avmu(ji,jj,ikbu) * zwx(ji,ikbum1)   &
                            / ( fse3u(ji,jj,ikbum1)*fse3uw(ji,jj,ikbu) )
            ! subtract surface forcing and bottom friction trend from vertical
            ! diffusive momentum trend
            ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
            ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
         END DO

         ! 2. Vertical diffusion on v
         ! ---------------------------

         ! 2.0 Matrix and second member construction
         ! bottom boundary condition: only zws must be masked as avmv can take
         ! non zero value at the ocean bottom depending on the bottom friction
         ! used (see zdfmix.F)
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               zcoef = -z2dt/fse3v(ji,jj,jk)
               zwi(ji,jk) = zcoef * avmv(ji,jj,jk  ) / fse3vw(ji,jj,jk  )
               zzws       = zcoef * avmv(ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
               zws(ji,jk) =  zzws * vmask(ji,jj,jk+1)
               zwd(ji,jk) = 1. - zwi(ji,jk) - zzws
               zwy(ji,jk) = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)
            END DO  
         END DO  

         ! 2.1 Surface boudary conditions
         DO ji = 2, jpim1
            z2dtf = z2dt / ( fse3v(ji,jj,1)*rau0 )
            zwi(ji,1) = 0.e0
            zwd(ji,1) = 1. - zws(ji,1)
            zwy(ji,1) = zwy(ji,1) + z2dtf * tauy(ji,jj)
         END DO  

         ! 2.2 Matrix inversion starting from the first level
         ikst = 1
#include "zdf.matrixsolver.h90"

         ! 2.3 Normalization to obtain the general momentum trend va
         DO jk = 1, jpkm1
            DO ji = 2, jpim1
               va(ji,jj,jk) = ( zwx(ji,jk) - vb(ji,jj,jk) ) / z2dt
            END DO    
         END DO   

         ! 2.4 diagnose surface and bottom momentum fluxes
         DO ji = 2, jpim1
            ! save the surface forcing momentum fluxes
            ztsy(ji,jj) = tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 )
            ! save bottom friction momentum fluxes
            ikbv   = min( mbathy(ji,jj+1), mbathy(ji,jj) )
            ikbvm1 = max( ikbv-1, 1 )
            ztby(ji,jj) = - avmv(ji,jj,ikbv) * zwx(ji,ikbvm1)   &
                            / ( fse3v(ji,jj,ikbvm1)*fse3vw(ji,jj,ikbv) )
            ! subtract surface forcing and bottom friction trend from vertical
            ! diffusive momentum trend
            ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
            ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the vertical diffusive trends for diagnostic
      ! momentum trends
      IF( l_trddyn )  THEN 
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(uldftrd, vldftrd, jpdtdldf, 'DYN', kt)
         CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
         ztdua(:,:,:) = 0.e0
         ztdva(:,:,:) = 0.e0
         ztdua(:,:,1) = ztsx(:,:)
         ztdva(:,:,1) = ztsy(:,:)
         CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
         ztdua(:,:,:) = 0.e0
         ztdva(:,:,:) = 0.e0
         ztdua(:,:,1) = ztbx(:,:)
         ztdva(:,:,1) = ztby(:,:)
         CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_zdf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module                       NO rotation of the mixing tensor
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_zdf_iso( kt )                        ! Dummy routine
      if(lwp) WRITE(numout,*) 'dyn_zdf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_zdf_iso
#endif

   !!==============================================================================
END MODULE dynzdf_iso
