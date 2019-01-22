MODULE dynldf_iso
   !!======================================================================
   !!                     ***  MODULE  dynldf_iso  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
#if defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                      slopes of the direction of mixing
   !!----------------------------------------------------------------------
   !!   dyn_ldf_iso  : update the momentum trend with the horizontal part
   !!                  of the lateral diffusion using isopycnal or horizon-
   !!                  tal s-coordinate laplacian operator.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE ldfslp          ! iso-neutral slopes 
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_ldf_iso           ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldfdyn_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynldf_iso.F90,v 1.6 2005/09/02 15:45:24 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_ldf_iso( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_iso  ***
      !!                       
      !! ** Purpose :   Compute the before trend of the horizontal part of the
      !!      lateral momentum diffusion and add it to the general trend of
      !!      momentum equation.
      !!
      !! ** Method :
      !!        The horizontal component of the lateral diffusive trends on
      !!      momentum is provided by a 2nd order operator rotated along neu-
      !!      tral or geopotential surfaces (in s-coordinates).
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp or inildf.
      !!      Here, u and v components are considered as 2 independent scalar
      !!      fields. Therefore, the property of splitting divergent and rota-
      !!      tional part of the flow of the standard, z-coordinate laplacian
      !!      momentum diffusion is lost.
      !!      horizontal fluxes associated with the rotated lateral mixing:
      !!      u-component:
      !!         ziut = ( ahmt + ahmb0 ) e2t * e3t / e1t  di[ ub ]
      !!               -      ahmt       e2t * mi-1(uslp) dk[ mi(mk(ub)) ]
      !!         zjuf = ( ahmf + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
      !!               -      ahmf       e1f * mi(vslp)   dk[ mj(mk(ub)) ]
      !!      v-component:
      !!         zivf = ( ahmf + ahmb0 ) e2t * e3t / e1t  di[ vb ]
      !!               -      ahmf       e2t * mj(uslp)   dk[ mi(mk(vb)) ]
      !!         zjvt = ( ahmt + ahmb0 ) e1f * e3f / e2f  dj[ ub ]
      !!               -      ahmt       e1f * mj-1(vslp) dk[ mj(mk(vb)) ]
      !!      take the horizontal divergence of the fluxes:
      !!         diffu = 1/(e1u*e2u*e3u) {  di  [ ziut ] + dj-1[ zjuf ]  }
      !!         diffv = 1/(e1v*e2v*e3v) {  di-1[ zivf ] + dj  [ zjvt ]  }
      !!      Add this trend to the general trend (ua,va):
      !!         ua = ua + diffu
      !!      'key_trddyn' defined: the trends are saved for diagnostics.
      !!
      !! ** Action :
      !!        Update (ua,va) arrays with the before geopotential biharmonic
      !!      mixing trend.
      !!        Save in (uldftrd,vldftrd) arrays the trends if 'key_trddyn' defined
      !!
      !! History :
      !!   8.0  !  97-07  (G. Madec)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
                         ztdva => sa      ! use sa as 3D workspace   
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   &
         zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
         zmskt, zmskf, zbu, zbv,       &
         zuah, zvah
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         ziut, zjuf, zjvt, zivf,        & ! temporary workspace
         zdku, zdk1u, zdkv, zdk1v
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_iso : iso-neutral laplacian diffusive operator or '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate horizontal diffusive operator'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Vertical u- and v-shears at level jk and jk+1
         ! ---------------------------------------------
         ! surface boundary condition: zdku(jk=1)=zdku(jk=2)
         !                             zdkv(jk=1)=zdkv(jk=2)

         zdk1u(:,:) = ( ub(:,:,jk) -ub(:,:,jk+1) ) * umask(:,:,jk+1)
         zdk1v(:,:) = ( vb(:,:,jk) -vb(:,:,jk+1) ) * vmask(:,:,jk+1)

         IF( jk == 1 ) THEN
            zdku(:,:) = zdk1u(:,:)
            zdkv(:,:) = zdk1v(:,:)
         ELSE
            zdku(:,:) = ( ub(:,:,jk-1) - ub(:,:,jk) ) * umask(:,:,jk)
            zdkv(:,:) = ( vb(:,:,jk-1) - vb(:,:,jk) ) * vmask(:,:,jk)
         ENDIF

         !                               -----f-----
         ! Horizontal fluxes on U             |  
         ! --------------------===        t   u   t
         !                                    |  
         ! i-flux at t-point             -----f-----

         DO jj = 2, jpjm1
            DO ji = fs_2, jpi   ! vector opt.
               zabe1 = ( fsahmt(ji,jj,jk) + ahmb0 )   &
#if defined key_partial_steps
                     * e2t(ji,jj) * MIN( fse3u(ji,jj,jk), fse3u(ji-1, jj,jk) ) / e1t(ji,jj)
#else
                     * e2t(ji,jj) * fse3t(ji,jj,jk) / e1t(ji,jj)
#endif

               zmskt = 1./MAX(  umask(ji-1,jj,jk  )+umask(ji,jj,jk+1)   &
                              + umask(ji-1,jj,jk+1)+umask(ji,jj,jk  ), 1. )

               zcof1 = - aht0 * e2t(ji,jj) * zmskt   &
                     * 0.5  * ( uslp(ji-1,jj,jk) + uslp(ji,jj,jk) )

               ziut(ji,jj) = tmask(ji,jj,jk) *   &
                           (  zabe1 * ( ub(ji,jj,jk) - ub(ji-1,jj,jk) )   &
                            + zcof1 * ( zdku (ji,jj) + zdk1u(ji-1,jj)     &
                                       +zdk1u(ji,jj) + zdku (ji-1,jj) )  )
            END DO
         END DO

         ! j-flux at f-point
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zabe2 = ( fsahmf(ji,jj,jk) + ahmb0 )   &
                     * e1f(ji,jj) * fse3f(ji,jj,jk) / e2f(ji,jj)

               zmskf = 1./MAX(  umask(ji,jj+1,jk  )+umask(ji,jj,jk+1)   &
                              + umask(ji,jj+1,jk+1)+umask(ji,jj,jk  ), 1. )

               zcof2 = - aht0 * e1f(ji,jj) * zmskf   &
                     * 0.5  * ( vslp(ji+1,jj,jk) + vslp(ji,jj,jk) )

               zjuf(ji,jj) = fmask(ji,jj,jk) *   &
                           (  zabe2 * ( ub(ji,jj+1,jk) - ub(ji,jj,jk) )   &
                            + zcof2 * ( zdku (ji,jj+1) + zdk1u(ji,jj)     &
                                       +zdk1u(ji,jj+1) + zdku (ji,jj) )  )
            END DO
         END DO

         !                                |   t   |
         ! Horizontal fluxes on V         |       |
         ! --------------------===        f---v---f
         !                                |       |
         ! i-flux at f-point              |   t   |

         DO jj = 2, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               zabe1 = ( fsahmf(ji,jj,jk) + ahmb0 )   &
                     * e2f(ji,jj) * fse3f(ji,jj,jk) / e1f(ji,jj)

               zmskf = 1./MAX(  vmask(ji+1,jj,jk  )+vmask(ji,jj,jk+1)   &
                              + vmask(ji+1,jj,jk+1)+vmask(ji,jj,jk  ), 1. )

               zcof1 = - aht0 * e2f(ji,jj) * zmskf   &
                     * 0.5 * ( uslp(ji,jj+1,jk) + uslp(ji,jj,jk) )

               zivf(ji,jj) = fmask(ji,jj,jk) *   &
                           (  zabe1 * ( vb(ji+1,jj,jk) - vb(ji,jj,jk) )   &
                            + zcof1 * ( zdkv (ji,jj) + zdk1v(ji+1,jj)     &
                                       +zdk1v(ji,jj) + zdkv (ji+1,jj) )  )
            END DO
         END DO

         ! j-flux at t-point
         DO jj = 2, jpj
            DO ji = 1, fs_jpim1   ! vector opt.
               zabe2 = ( fsahmt(ji,jj,jk) + ahmb0 )   &
#if defined key_partial_steps
                     * e1t(ji,jj) * MIN( fse3v(ji,jj,jk), fse3v(ji, jj-1, jk) ) / e2t(ji,jj)
#else
                     * e1t(ji,jj) * fse3t(ji,jj,jk) / e2t(ji,jj)
#endif

               zmskt = 1./MAX(  vmask(ji,jj-1,jk  )+vmask(ji,jj,jk+1)   &
                              + vmask(ji,jj-1,jk+1)+vmask(ji,jj,jk  ), 1. )

               zcof2 = - aht0 * e1t(ji,jj) * zmskt   &
                     * 0.5 * ( vslp(ji,jj-1,jk) + vslp(ji,jj,jk) )

               zjvt(ji,jj) = tmask(ji,jj,jk) *   &
                           (  zabe2 * ( vb(ji,jj,jk) - vb(ji,jj-1,jk) )   &
                            + zcof2 * ( zdkv (ji,jj-1) + zdk1v(ji,jj)     &
                                       +zdk1v(ji,jj-1) + zdkv (ji,jj) )  )
            END DO
         END DO


         ! Second derivative (divergence) and add to the general trend
         ! -----------------------------------------------------------

         DO jj = 2, jpjm1
            DO ji = 2, jpim1          !! Question vectop possible??? !!bug
               ! volume elements
               zbu = e1u(ji,jj) * e2u(ji,jj) * fse3u(ji,jj,jk)
               zbv = e1v(ji,jj) * e2v(ji,jj) * fse3v(ji,jj,jk)
               ! horizontal component of isopycnal momentum diffusive trends
               zuah =( ziut (ji+1,jj) - ziut (ji,jj  ) +   &
                       zjuf (ji  ,jj) - zjuf (ji,jj-1)  ) / zbu
               zvah =( zivf (ji,jj  ) - zivf (ji-1,jj) +   &
                       zjvt (ji,jj+1) - zjvt (ji,jj  )  ) / zbv
               ! add the trends to the general trends
               ua (ji,jj,jk) = ua (ji,jj,jk) + zuah
               va (ji,jj,jk) = va (ji,jj,jk) + zvah
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the lateral diffusion trends for diagnostic
      ! momentum trends will be saved in dynzdf_iso.F90
      IF( l_trddyn )   THEN
         uldftrd(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         vldftrd(:,:,:) = va(:,:,:) - ztdva(:,:,:)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' ldf  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_ldf_iso

# else
   !!----------------------------------------------------------------------
   !!   Dummy module                           NO rotation of mixing tensor
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_ldf_iso( kt )               ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_ldf_iso: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_ldf_iso
#endif

   !!======================================================================
END MODULE dynldf_iso
