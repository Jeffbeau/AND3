MODULE trazdf_iso_vopt
   !!==============================================================================
   !!                 ***  MODULE  trazdf_iso_vopt  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!==============================================================================
#if defined key_ldfslp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   tra_zdf_iso_vopt : Update the tracer trend with the vertical part of 
   !!                  the isopycnal or geopotential s-coord. operator and
   !!                  the vertical diffusion. vector optimization, use
   !!                  k-j-i loops.
   !!   tra_zdf_iso  :
   !!   tra_zdf_zdf  :
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE ldfslp          ! iso-neutral slopes 
   USE zdfddm          ! ocean vertical physics: double diffusion
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE zdfkpp          ! KPP parameterisation
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_zdf_iso_vopt   !  routine called by step.F90

   !! * Module variables
   REAL(wp), DIMENSION(jpk) ::  &
      r2dt                          ! vertical profile of 2 x time-step
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
      tavg, savg,                     & ! workspace arrays
      tdta, tdsa                        ! workspace arrays

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "ldftra_substitute.h90"
#  include "ldfeiv_substitute.h90"
#  include "zdfddm_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_iso_vopt.F90,v 1.10 2005/09/22 10:29:06 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE tra_zdf_iso_vopt( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_iso_vopt  ***
      !!
      !! ** Purpose :
      !! ** Method  :
      !! ** Action  :
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!   9.0  !  05-06  (C. Ethe) KPP parameterization
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      IF( kt == nit000 ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'tra_zdf_iso_vopt : vertical mixing computation'
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~~   is  iso-neutral diffusion : implicit vertical time stepping'
#if defined key_diaeiv 
         w_eiv(:,:,:) = 0.e0
#endif
      ENDIF

      ! initialization step
      tavg(:,:,:) = 0.e0
      savg(:,:,:) = 0.e0

      ! I. vertical extra-diagonal part of the rotated tensor
      ! -----------------------------------------------------

      CALL tra_zdf_iso( kt )

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 1- Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF

      ! II. vertical diffusion (including the vertical diagonal part of the rotated tensor)
      ! ----------------------

      CALL tra_zdf_zdf( kt )

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf 2- Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask)
      ENDIF

   END SUBROUTINE tra_zdf_iso_vopt


   SUBROUTINE tra_zdf_zdf( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_zdf_zdf  ***
      !!                    
      !! ** Purpose :   Compute the trend due to the vertical tracer diffusion
      !!     including the vertical component of lateral mixing (only for 2nd
      !!     order operator, for fourth order it is already computed and add
      !!     to the general trend in traldf.F) and add it to the general trend
      !!     of the tracer equations.
      !!
      !! ** Method  :   The vertical component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geo-
      !!      potential surfaces to which an eddy induced advection can be 
      !!      added. It is computed using before fields (forward in time) and 
      !!      isopycnal or geopotential slopes computed in routine ldfslp.
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers (t & s) is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=ta).
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!         ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in zdfmxl.F90).
      !!         avt = zavt
      !!         (avs = zavs if lk_zdfddm=T )
      !!
      !!      'key_trdtra' defined: trend saved for futher diagnostics.
      !!
      !!      macro-tasked on vertical slab (jj-loop)
      !!
      !! ** Action  : - Update (ta,sa) with before vertical diffusion trend
      !!              - Save the trend in (ztdta,ztdsa) ('key_trdtra')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11 (G. Madec)
      !!        !  92-06 (M. Imbard) correction on tracer trend loops
      !!        !  96-01 (G. Madec) statement function for e3
      !!        !  97-05 (G. Madec) vertical component of isopycnal
      !!        !  97-07 (G. Madec) geopotential diffusion in s-coord
      !!        !  00-08 (G. Madec) double diffusive mixing
      !!   8.5  !  02-08 (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08 (C. Talandier) New trends organization
      !!   9.0  !  05-06 (C. Ethe )  non-local flux in KPP vertical mixing scheme
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce    , ONLY :   zwd   => ua,  &  ! ua, va used as
                            zws   => va      ! workspace
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   &
         zavi, zrhs                          ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwi, zwt, zavsi                     ! temporary workspace arrays



      ! I. Local constant initialization
      ! --------------------------------
      ! ... time step = 2 rdttra ex
      IF( neuler == 0 .AND. kt == nit000 ) THEN
         r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
      ELSEIF( kt <= nit000 + 1) THEN
         r2dt(:) = 2. * rdttra(:)          ! leapfrog
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         tdta(:,:,:) = ta(:,:,:) 
         tdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      zwd  (1,:, : ) = 0.e0     ;     zwd  (jpi,:,:) = 0.e0
      zws  (1,:, : ) = 0.e0     ;     zws  (jpi,:,:) = 0.e0
      zwi  (1,:, : ) = 0.e0     ;     zwi  (jpi,:,:) = 0.e0
      zwt  (1,:, : ) = 0.e0     ;     zwt  (jpi,:,:) = 0.e0
      zavsi(1,:, : ) = 0.e0     ;     zavsi(jpi,:,:) = 0.e0
      zwt  (:,:,jpk) = 0.e0     ;     zwt  ( : ,:,1) = 0.e0
      zavsi(:,:,jpk) = 0.e0     ;     zavsi( : ,:,1) = 0.e0

      ! II. Vertical trend associated with the vertical physics
      ! =======================================================
      !     (including the vertical flux proportional to dk[t] associated
      !      with the lateral mixing, through the avt update)
      !     dk[ avt dk[ (t,s) ] ] diffusive trends


      ! II.0 Matrix construction
      ! ------------------------

      ! update and save of avt (and avs if double diffusive mixing)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zavi = fsahtw(ji,jj,jk) * (                 &   ! vertical mixing coef. due to lateral mixing
                    wslpi(ji,jj,jk) * wslpi(ji,jj,jk)      &
                  + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
               zwt(ji,jj,jk) = avt(ji,jj,jk) + zavi            ! zwt=avt+zavi (total vertical mixing coef. on temperature)
#if defined key_zdfddm
               zavsi(ji,jj,jk) = fsavs(ji,jj,jk) + zavi        ! dd mixing: zavsi = total vertical mixing coef. on salinity
#endif
            END DO
         END DO
      END DO

      ! Diagonal, inferior, superior  (including the bottom boundary condition via avt masked)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwi(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jj,jk) = - r2dt(jk) * zwt(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Surface boudary conditions
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0.e0
            zwd(ji,jj,1) = 1. - zws(ji,jj,1)
         END DO
      END DO


      ! II.1. Vertical diffusion on t
      ! ---------------------------

      !! Matrix inversion from the first level
      !!----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The second member is in 2d array zwy
      !   The solution is in 2d array zwx
      !   The 3d arry zwt is a work space array
      !   zwy is used and then used as a work space array : its value is modified!

      ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            zwt(ji,jj,1) = zwd(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
#if defined key_zdfkpp
         ! add non-local temperature flux ( in convective case only)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            ta(ji,jj,1) = tb(ji,jj,1) + r2dt(1) * ta(ji,jj,1) &
                  &  - r2dt(1) * ( ghats(ji,jj,1) * avt(ji,jj,1) - ghats(ji,jj,2) * avt(ji,jj,2) ) &
                  &               * wt0(ji,jj) / fse3t(ji,jj,2) 
         END DO
      END DO

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ! zrhs=right hand side 
               zrhs = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)  &
                  &  - r2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
                  &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
               ta(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *ta(ji,jj,jk-1)
            END DO
         END DO
      END DO
#else
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            ta(ji,jj,1) = tb(ji,jj,1) + r2dt(1) * ta(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zrhs = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)   ! zrhs=right hand side 
               ta(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *ta(ji,jj,jk-1)
            END DO
         END DO
      END DO
#endif

      ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
      ! Save the masked temperature after in ta
      ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            ta(ji,jj,jpkm1) = ta(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ta(ji,jj,jk) = ( ta(ji,jj,jk) - zws(ji,jj,jk) * ta(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      ! II.2 Vertical diffusion on salinity
      ! -----------------------------------

#if defined key_zdfddm
      ! Rebuild the Matrix as avt /= avs

      ! Diagonal, inferior, superior  (including the bottom boundary condition via avs masked)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zwi(ji,jj,jk) = - r2dt(jk) * zavsi(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
               zws(ji,jj,jk) = - r2dt(jk) * zavsi(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
               zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Surface boudary conditions
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zwi(ji,jj,1) = 0.e0
            zwd(ji,jj,1) = 1. - zws(ji,jj,1)
         END DO
      END DO
#endif


      !! Matrix inversion from the first level
      !!----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular
      !   matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The second member is in 2d array zwy
      !   The solution is in 2d array zwx
      !   The 3d arry zwt is a work space array
      !   zwy is used and then used as a work space array : its value is modified!

      ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            zwt(ji,jj,1) = zwd(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
            END DO
         END DO
      END DO

      ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1

#if defined key_zdfkpp
         ! add non-local temperature flux ( in convective case only)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            sa(ji,jj,1) = sb(ji,jj,1) + r2dt(1) * sa(ji,jj,1) &
                  &  - r2dt(1) * ( ghats(ji,jj,1) * fsavs(ji,jj,1) - ghats(ji,jj,2) * fsavs(ji,jj,2) ) &
                  &               * ws0(ji,jj) / fse3t(ji,jj,2) 
         END DO
      END DO

      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               ! zrhs=right hand side 
               zrhs = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)  &
                  &  - r2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
                  &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
               sa(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *sa(ji,jj,jk-1)
            END DO
         END DO
      END DO
#else
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            sa(ji,jj,1) = sb(ji,jj,1) + r2dt(1) * sa(ji,jj,1)
         END DO
      END DO
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               zrhs = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)   ! zrhs=right hand side
               sa(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) *sa(ji,jj,jk-1)
            END DO
         END DO
      END DO
#endif

      ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
      ! Save the masked temperature after in ta
      ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            sa(ji,jj,jpkm1) = sa(ji,jj,jpkm1) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
         END DO
      END DO
      DO jk = jpk-2, 1, -1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1
               sa(ji,jj,jk) = ( sa(ji,jj,jk) - zws(ji,jj,jk) * sa(ji,jj,jk+1) ) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
            END DO
         END DO
      END DO


      ! save the trends for diagnostic
      ! compute the vertical diffusive trends in substracting the previous 
      ! trends tdta()/tdsa() to the new one computed via dT/dt or dS/dt 
      ! i.e. with the new temperature/salinity ta/sa computed above
      IF( l_trdtra )   THEN
         IF( ln_traldf_iso)   THEN
            DO jk = 1, jpkm1
               tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk) + tavg(:,:,jk) 
               tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk) + savg(:,:,jk) 
            END DO
         ELSE
            DO jk = 1, jpkm1
               tdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) ) - tdta(:,:,jk)                             
               tdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) ) - tdsa(:,:,jk)                             
            END DO
         ENDIF

         CALL trd_mod(tdta, tdsa, jpttdzdf, 'TRA', kt)
      ENDIF

   END SUBROUTINE tra_zdf_zdf


   SUBROUTINE tra_zdf_iso( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_iso  ***
      !!
      !! ** Purpose :
      !!     Compute the trend due to the vertical tracer diffusion inclu-
      !!     ding the vertical component of lateral mixing (only for second
      !!     order operator, for fourth order it is already computed and
      !!     add to the general trend in traldf.F) and add it to the general
      !!     trend of the tracer equations.
      !!
      !! ** Method :
      !!         The vertical component of the lateral diffusive trends is
      !!      provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      First part: vertical trends associated with the lateral mixing
      !!      ==========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      save avt coef. resulting from vertical physics alone in zavt:
      !!         zavt = avt
      !!      update and save in zavt the vertical eddy viscosity coefficient:
      !!         avt = avt + wslpi^2+wslj^2
      !!      add vertical Eddy Induced advective fluxes (lk_traldf_eiv=T):
      !!         zftw = zftw + { di[aht e2u mi(wslpi)]
      !!                    +dj[aht e1v mj(wslpj)] } mk(tb)
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !! ** Action :
      !!         Update (ta,sa) arrays with the before vertical diffusion trend
      !!         Save in (ztdta,ztdsa) arrays the trends if 'key_trdtra' defined
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard) correction on tracer trend loops
      !!        !  96-01  (G. Madec) statement function for e3
      !!        !  97-05  (G. Madec) vertical component of isopycnal
      !!        !  97-07  (G. Madec) geopotential diffusion in s-coord
      !!        !  00-08  (G. Madec) double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!---------------------------------------------------------------------
      !! * Modules used
      USE oce    , ONLY :   zwx => ua,  &  ! use ua, va as
                            zwy => va      ! workspace arrays

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk       ! dummy loop indices
#if defined key_partial_steps
      INTEGER ::   iku, ikv
#endif
      REAL(wp) ::   &
         zcoef0, zcoef3,         &  !    "         "
         zcoef4,                 &  !    "         "
         zbtr, zmku, zmkv,       &  !    "         "
#if defined key_traldf_eiv
         zcoeg3,                 &  !    "         "
         zuwki, zvwki,           &  !    "         "
         zuwk, zvwk,             &  !    "         "
#endif
         ztav, zsav
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zwz, zwt, ztfw, zsfw       ! temporary workspace arrays


      ! 0. Local constant initialization
      ! --------------------------------
      zwx (1,:,:) = 0.e0     ;     zwx (jpi,:,:) = 0.e0
      zwy (1,:,:) = 0.e0     ;     zwy (jpi,:,:) = 0.e0
      zwz (1,:,:) = 0.e0     ;     zwz (jpi,:,:) = 0.e0
      zwt (1,:,:) = 0.e0     ;     zwt (jpi,:,:) = 0.e0
      ztfw(1,:,:) = 0.e0     ;     ztfw(jpi,:,:) = 0.e0
      zsfw(1,:,:) = 0.e0     ;     zsfw(jpi,:,:) = 0.e0

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         tdta(:,:,:) = ta(:,:,:) 
         tdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      ! I. Vertical trends associated with lateral mixing
      ! -------------------------------------------------
      !    (excluding the vertical flux proportional to dk[t] )


      ! I.1 horizontal tracer gradient
      ! ------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               ! i-gradient of T and S at jj
               zwx (ji,jj,jk) = ( tb(ji+1,jj,jk)-tb(ji,jj,jk) ) * umask(ji,jj,jk)
               zwy (ji,jj,jk) = ( sb(ji+1,jj,jk)-sb(ji,jj,jk) ) * umask(ji,jj,jk)
               ! j-gradient of T and S at jj
               zwz (ji,jj,jk) = ( tb(ji,jj+1,jk)-tb(ji,jj,jk) ) * vmask(ji,jj,jk)
               zwt (ji,jj,jk) = ( sb(ji,jj+1,jk)-sb(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
#  if defined key_partial_steps
      ! partial steps correction at the bottom ocean level 
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vector opt.
            ! last ocean level
            iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
            ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
            ! i-gradient of T and S
            zwx (ji,jj,iku) = gtu(ji,jj)
            zwy (ji,jj,iku) = gsu(ji,jj)
            ! j-gradient of T and S
            zwz (ji,jj,ikv) = gtv(ji,jj) 
            zwt (ji,jj,ikv) = gsv(ji,jj) 
         END DO
         END DO
#endif


      ! I.2 Vertical fluxes
      ! -------------------

      ! Surface and bottom vertical fluxes set to zero
      ztfw(:,:, 1 ) = 0.e0
      zsfw(:,:, 1 ) = 0.e0
      ztfw(:,:,jpk) = 0.e0
      zsfw(:,:,jpk) = 0.e0

      ! interior (2=<jk=<jpk-1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)

               zmku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
                              + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )

               zmkv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
                              + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )

               zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
               zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)

               ztfw(ji,jj,jk) = zcoef3 * (   zwx(ji  ,jj  ,jk-1) + zwx(ji-1,jj  ,jk)      &
                                           + zwx(ji-1,jj  ,jk-1) + zwx(ji  ,jj  ,jk)  )   &
                              + zcoef4 * (   zwz(ji  ,jj  ,jk-1) + zwz(ji  ,jj-1,jk)      &
                                           + zwz(ji  ,jj-1,jk-1) + zwz(ji  ,jj  ,jk)  )

               zsfw(ji,jj,jk) = zcoef3 * (   zwy(ji  ,jj  ,jk-1) + zwy(ji-1,jj  ,jk)      &
                                           + zwy(ji-1,jj  ,jk-1) + zwy(ji  ,jj  ,jk)  )   &
                              + zcoef4 * (   zwt(ji  ,jj  ,jk-1) + zwt(ji  ,jj-1,jk)      &
                                           + zwt(ji  ,jj-1,jk-1) + zwt(ji  ,jj  ,jk)  )
            END DO
         END DO
      END DO

#if defined key_traldf_eiv
      !                              ! ---------------------------------------!
      !                              ! Eddy induced vertical advective fluxes !
      !                              ! ---------------------------------------!
         zwx(:,:, 1 ) = 0.e0
         zwy(:,:, 1 ) = 0.e0
         zwx(:,:,jpk) = 0.e0
         zwy(:,:,jpk) = 0.e0

         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
#   if defined key_traldf_c2d || defined key_traldf_c3d
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * fsaeiu(ji-1,jj,jk) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * fsaeiu(ji  ,jj,jk) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * fsaeiv(ji,jj-1,jk) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * fsaeiv(ji,jj  ,jk) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
#   else
                  zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
                     &  * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
                     &  * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
                     &  * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
                     &  * e1v(ji  ,jj) * vmask(ji  ,jj,jk)

                  zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
                     &            * ( zuwk - zuwki + zvwk - zvwki )
#   endif
                  zwx(ji,jj,jk) = + zcoeg3 * ( tb(ji,jj,jk) + tb(ji,jj,jk-1) )
                  zwy(ji,jj,jk) = + zcoeg3 * ( sb(ji,jj,jk) + sb(ji,jj,jk-1) )

                  ztfw(ji,jj,jk) = ztfw(ji,jj,jk) + zwx(ji,jj,jk)
                  zsfw(ji,jj,jk) = zsfw(ji,jj,jk) + zwy(ji,jj,jk)
#   if defined key_diaeiv
                  w_eiv(ji,jj,jk) = -2. * zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
#   endif
               END DO
            END DO
         END DO
#endif

      ! I.5 Divergence of vertical fluxes added to the general tracer trend
      ! -------------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
               ztav = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
               zsav = (  zsfw(ji,jj,jk) - zsfw(ji,jj,jk+1)  ) * zbtr
               ta(ji,jj,jk) = ta(ji,jj,jk) + ztav
               sa(ji,jj,jk) = sa(ji,jj,jk) + zsav
            END DO
         END DO
      END DO

      ! save the trends for diagnostic
      !  WARNING jpttddoe is used here for vertical Gent velocity trend not for damping !!!
      IF( l_trdtra )   THEN
#   if defined key_traldf_eiv
         ! Compute the vertical Gent velocity trend
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
                  tavg(ji,jj,jk) = ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) ) * zbtr
                  savg(ji,jj,jk) = ( zwy(ji,jj,jk) - zwy(ji,jj,jk+1) ) * zbtr
               END DO
            END DO
         END DO

         CALL trd_mod(tavg, savg, jpttddoe, 'TRA', kt)
#   endif
         ! Recompute the divergence of vertical fluxes ztav & zsav trends 
         ! computed at step 1.5 above in making the difference between the new 
         ! trend ta()/sa() and the previous one tdta()/tdsa() and substract 
         ! the vertical Gent velocity trend tavg()/savg() (zero if not used)
         tavg(:,:,:) = ta(:,:,:) - tdta(:,:,:) - tavg(:,:,:) 
         savg(:,:,:) = sa(:,:,:) - tdsa(:,:,:) - savg(:,:,:) 
      ENDIF

   END SUBROUTINE tra_zdf_iso

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :             NO rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE tra_zdf_iso_vopt( kt )              ! empty routine
!      WRITE(*,*) 'tra_zdf_iso_vopt: You should not have seen this print! error?', kt
   END SUBROUTINE tra_zdf_iso_vopt
#endif

   !!==============================================================================
END MODULE trazdf_iso_vopt
