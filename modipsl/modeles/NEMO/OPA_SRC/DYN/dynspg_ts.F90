MODULE dynspg_ts
   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
#if ( defined key_dynspg_ts && ! defined key_autotasking ) ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_ts'     free surface cst volume with time splitting
   !!   NOT 'key_autotasking'                      k-j-i loop (vector opt.)
   !!----------------------------------------------------------------------
   !!   dyn_spg_ts  : compute surface pressure gradient trend using a time-
   !!                 splitting scheme and add to the general trend 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ocesbc          ! ocean surface boundary condition
   USE obcdta          ! open boundary condition data     
   USE obcfla          ! Flather open boundary condition  
   USE dynvor          ! vorticity term
   USE obc_oce         ! Lateral open boundary condition
   USE obc_par         ! open boundary condition parameters
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE dynspg_oce      ! surface pressure gradient variables
   USE in_out_manager  ! I/O manager
#ifdef key_RIVER_INPUT
!!DB
   USE rivers
#endif


   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_ts  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_ts.F90,v 1.6 2006/01/03 15:04:14 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_ts( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts  ***
      !!
      !! ** Purpose :   Compute the now trend due to the surface pressure
      !!      gradient in case of free surface formulation with time-splitting.
      !!      Add it to the general trend of momentum equation.
      !!      Compute the free surface.
      !!
      !! ** Method  :   Free surface formulation with time-splitting
      !!      -1- Save the vertically integrated trend. This general trend is
      !!          held constant over the barotropic integration.
      !!          The Coriolis force is removed from the general trend as the
      !!          surface gradient and the Coriolis force are updated within
      !!          the barotropic integration.
      !!      -2- Barotropic loop : updates of sea surface height (ssha_e) and 
      !!          barotropic transports (ua_e and va_e) through barotropic 
      !!          momentum and continuity integration. Barotropic former 
      !!          variables are time averaging over the full barotropic cycle
      !!          (= 2 * baroclinic time step) and saved in zsshX_b, zuX_b 
      !!          and zvX_b (X specifying after, now or before).
      !!      -3- Update of sea surface height from time averaged barotropic 
      !!          variables.
      !!        - apply lateral boundary conditions on sshn.
      !!      -4- The new general trend becomes :
      !!          ua = ua - sum_k(ua)/H + ( zua_b - sum_k(ub) )/H
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References :
      !!   Griffies et al., (2003): A technical guide to MOM4. NOAA/GFDL
      !!
      !! History :
      !!   9.0  !  04-12  (L. Bessieres, G. Madec)  Original code
      !!        !  05-11  (V. Garnier, G. Madec)  optimization
      !!---------------------------------------------------------------------


      !! * Arguments
     INTEGER, INTENT( in )  ::   kt           ! ocean time-step index

      !! * Local declarations
      INTEGER  ::  ji, jj, jk, jit             ! dummy loop indices
      INTEGER  ::  icycle                      ! temporary scalar
      REAL(wp) ::                           &
         zraur, zcoef, z2dt_e, z2dt_b, zfac25,   &  ! temporary scalars
         zfact1, zspgu, zcubt, zx1, zy1,    &  !     "        "
         zfact2, zspgv, zcvbt, zx2, zy2        !     "        "
      REAL(wp), DIMENSION(jpi,jpj) ::       &
         zcu, zcv, zwx, zwy, zhdiv,         &  ! temporary arrays
         zua, zva, zub, zvb,                &  !     "        "
         zssha_b, zua_b, zva_b,             &  !     "        "
         zsshb_e, zub_e, zvb_e,             &  !     "        "
         zun_e, zvn_e                          !     "        "
      REAL(wp), DIMENSION(jpi,jpj),SAVE ::  &
         ztnw, ztne, ztsw, ztse
      !!----------------------------------------------------------------------

      ! Arrays initialization
      ! ---------------------
      zua_b(:,:) = 0.e0   ;   zub_e(:,:) = 0.e0   ;   zun_e(:,:) = 0.e0
      zva_b(:,:) = 0.e0   ;   zvb_e(:,:) = 0.e0   ;   zvn_e(:,:) = 0.e0
      zhdiv(:,:) = 0.e0


      IF( kt == nit000 ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*) ' Number of sub cycle in 1 time-step (2 rdt) : icycle = ', FLOOR( 2*rdt/rdtbt )

         IF( .NOT. ln_rstart ) THEN
            ! initialize barotropic specific arrays
            sshb_b(:,:) = sshb(:,:)
            sshn_b(:,:) = sshn(:,:)
            un_b(:,:)   = 0.e0
            vn_b(:,:)   = 0.e0
            ! vertical sum
            IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
               DO jk = 1, jpkm1
                  DO ji = 1, jpij
                     un_b(ji,1) = un_b(ji,1) + fse3u(ji,1,jk) * un(ji,1,jk)
                     vn_b(ji,1) = vn_b(ji,1) + fse3v(ji,1,jk) * vn(ji,1,jk)
                  END DO
               END DO
            ELSE                             ! No  vector opt.
               DO jk = 1, jpkm1
                  un_b(:,:) = un_b(:,:) + fse3u(:,:,jk) * un(:,:,jk)
                  vn_b(:,:) = vn_b(:,:) + fse3v(:,:,jk) * vn(:,:,jk)
               END DO
            ENDIF
         ENDIF
         ssha_e(:,:) = sshn(:,:)
         ua_e(:,:)   = un_b(:,:)
         va_e(:,:)   = vn_b(:,:)

         IF( ln_dynvor_een ) THEN
            ztne(1,:) = 0.e0   ;   ztnw(1,:) = 0.e0   ;   ztse(1,:) = 0.e0   ;   ztsw(1,:) = 0.e0
            DO jj = 2, jpj
               DO ji = fs_2, jpi   ! vector opt.
                  ztne(ji,jj) = ( ff(ji-1,jj  ) + ff(ji  ,jj  ) + ff(ji  ,jj-1) ) / 3.
                  ztnw(ji,jj) = ( ff(ji-1,jj-1) + ff(ji-1,jj  ) + ff(ji  ,jj  ) ) / 3.
                  ztse(ji,jj) = ( ff(ji  ,jj  ) + ff(ji  ,jj-1) + ff(ji-1,jj-1) ) / 3.
                  ztsw(ji,jj) = ( ff(ji  ,jj-1) + ff(ji-1,jj-1) + ff(ji-1,jj  ) ) / 3.
               END DO
            END DO
         ENDIF

      ENDIF      !!kt==nit000
    
      ! Local constant initialization
      ! --------------------------------
      z2dt_b = 2.0 * rdt                                    ! baroclinic time step
      IF ( neuler == 0 .AND. kt == nit000 ) z2dt_b = rdt
      zfact1 = 0.5 * 0.25                                   ! coefficient for vorticity estimates
      zfact2 = 0.5 * 0.5
      zraur  = 1. / rauw                                    ! 1 / volumic mass of pure water
      
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------

      ! Vertically integrated quantities
      ! --------------------------------
      zua(:,:) = 0.e0
      zva(:,:) = 0.e0
      zub(:,:) = 0.e0
      zvb(:,:) = 0.e0
      zwx(:,:) = 0.e0
      zwy(:,:) = 0.e0

      ! vertical sum
      IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
         DO jk = 1, jpkm1
            DO ji = 1, jpij
               !                                                           ! Vertically integrated momentum trends
               zua(ji,1) = zua(ji,1) + fse3u(ji,1,jk) * umask(ji,1,jk) * ua(ji,1,jk)
               zva(ji,1) = zva(ji,1) + fse3v(ji,1,jk) * vmask(ji,1,jk) * va(ji,1,jk)
               !                                                           ! Vertically integrated transports (before)
               zub(ji,1) = zub(ji,1) + fse3u(ji,1,jk) * ub(ji,1,jk)
               zvb(ji,1) = zvb(ji,1) + fse3v(ji,1,jk) * vb(ji,1,jk)
               !                                                           ! Planetary vorticity transport fluxes (now)
               zwx(ji,1) = zwx(ji,1) + e2u(ji,1) * fse3u(ji,1,jk) * un(ji,1,jk)
               zwy(ji,1) = zwy(ji,1) + e1v(ji,1) * fse3v(ji,1,jk) * vn(ji,1,jk)
            END DO
         END DO
      ELSE                             ! No  vector opt.
         DO jk = 1, jpkm1
            !                                                           ! Vertically integrated momentum trends
            zua(:,:) = zua(:,:) + fse3u(:,:,jk) * umask(:,:,jk) * ua(:,:,jk)
            zva(:,:) = zva(:,:) + fse3v(:,:,jk) * vmask(:,:,jk) * va(:,:,jk)
            !                                                           ! Vertically integrated transports (before)
            zub(:,:) = zub(:,:) + fse3u(:,:,jk) * ub(:,:,jk)
            zvb(:,:) = zvb(:,:) + fse3v(:,:,jk) * vb(:,:,jk)
            !                                                           ! Planetary vorticity (now)
            zwx(:,:) = zwx(:,:) + e2u(:,:) * fse3u(:,:,jk) * un(:,:,jk)
            zwy(:,:) = zwy(:,:) + e1v(:,:) * fse3v(:,:,jk) * vn(:,:,jk)
         END DO
      ENDIF

      IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               ! energy conserving formulation for planetary vorticity term
               zcu(ji,jj) = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
               zcv(ji,jj) =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
            END DO
         END DO

      ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
                              + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
                              + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               zcu(ji,jj)  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
               zcv(ji,jj)  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
            END DO
         END DO

      ELSEIF ( ln_dynvor_een ) THEN                    ! enstrophy and energy conserving scheme
         zfac25 = 0.25
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zcu(ji,jj) = + zfac25 / e1u(ji,jj)   &
                  &       * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                  &          + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zcv(ji,jj) = - zfac25 / e2v(ji,jj)   &
                  &       * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                  &          + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
            END DO
         END DO

      ENDIF


      ! Remove barotropic trend from general momentum trend
      ! ---------------------------------------------------
      DO jk = 1 , jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - zua(ji,jj) * hur(ji,jj) !JC:zua/zva aren't updated at river points
               va(ji,jj,jk) = va(ji,jj,jk) - zva(ji,jj) * hvr(ji,jj) !but ua/vb will be in dynnxt
            END DO
         END DO
      END DO

      ! Remove coriolis term from barotropic trend
      ! ------------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1
            zua(ji,jj) = zua(ji,jj) - zcu(ji,jj)
            zva(ji,jj) = zva(ji,jj) - zcv(ji,jj)
         END DO
      END DO

      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations with time splitting
      ! -----------------------------------------------------------------------

      ! Initialisations
      !----------------
      ! Number of iteration of the barotropic loop
      icycle = FLOOR( z2dt_b / rdtbt )

      ! variables for the barotropic equations
      zsshb_e(:,:) = sshn_b(:,:)       ! (barotropic) sea surface height (before and now)
      sshn_e (:,:) = sshn_b(:,:)
      zub_e  (:,:) = un_b  (:,:)       ! barotropic transports issued from the barotropic equations (before and now)
      zvb_e  (:,:) = vn_b  (:,:)
      zun_e  (:,:) = un_b  (:,:)
      zvn_e  (:,:) = vn_b  (:,:)
      zssha_b(:,:) = sshn  (:,:)       ! time averaged variables over all sub-timesteps
      zua_b  (:,:) = un_b  (:,:)   
      zva_b  (:,:) = vn_b  (:,:)

      ! set ssh corrections to 0
      ! ssh corrections are applied to normal velocities (Flather's algorithm) and averaged over the barotropic loop
#if defined key_obc
      IF( lp_obc_east  )   sshfoe_b(:,:) = 0.e0
      IF( lp_obc_west  )   sshfow_b(:,:) = 0.e0
      IF( lp_obc_south )   sshfos_b(:,:) = 0.e0
      IF( lp_obc_north )   sshfon_b(:,:) = 0.e0
#endif

      ! Barotropic integration over 2 baroclinic time steps
      ! ---------------------------------------------------

      !                                                    ! ==================== !
      DO jit = 1, icycle                                   !  sub-time-step loop  !
         !                                                 ! ==================== !

         z2dt_e = 2. * rdtbt
         IF ( jit == 1 )   z2dt_e = rdtbt

         ! Time interpolation of open boundary condition data
         IF( lk_obc )   CALL obc_dta_bt( kt, jit )

         ! Horizontal divergence of barotropic transports
         !--------------------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zhdiv(ji,jj) = ( e2u(ji  ,jj  ) * zun_e(ji  ,jj)              &
                  &            -e2u(ji-1,jj  ) * zun_e(ji-1,jj)              &
                  &            +e1v(ji  ,jj  ) * zvn_e(ji  ,jj)              &
                  &            -e1v(ji  ,jj-1) * zvn_e(ji  ,jj-1) )          &
                  &           / (e1t(ji,jj)*e2t(ji,jj))
            END DO
         END DO

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0      ! east
         IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0      ! west
         IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0      ! north
         IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0      ! south
#endif

         ! Sea surface height from the barotropic system
         !----------------------------------------------
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ssha_e(ji,jj) = ( zsshb_e(ji,jj) - z2dt_e *  ( zraur * emp(ji,jj)  &
            &            +  zhdiv(ji,jj) ) ) * tmask(ji,jj,1)
            END DO
         END DO

         ! evolution of the barotropic transport ( following the vorticity scheme used)
         ! ----------------------------------------------------------------------------
         zwx(:,:) = e2u(:,:) * zun_e(:,:)
         zwy(:,:) = e1v(:,:) * zvn_e(:,:)

         IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! surface pressure gradient
                  zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
                  zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
                  ! energy conserving formulation for planetary vorticity term
                  zy1 = ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
                  zy2 = ( zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
                  zx1 = ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
                  zx2 = ( zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
                  zcubt = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
                  zcvbt =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
                  ! after transports
                  ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
                  va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
               END DO
            END DO

         ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! surface pressure gradient
                  zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
                  zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
                  ! enstrophy conserving formulation for planetary vorticity term
                  zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
                                 + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
                  zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
                                 + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
                  zcubt  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
                  zcvbt  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
                  ! after transports
                  ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
                  va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
               END DO
            END DO

         ELSEIF ( ln_dynvor_een ) THEN                    ! energy and enstrophy conserving scheme
            zfac25 = 0.25
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ! surface pressure gradient
                  zspgu = -grav * ( sshn_e(ji+1,jj) - sshn_e(ji,jj) ) * hu(ji,jj) / e1u(ji,jj)
                  zspgv = -grav * ( sshn_e(ji,jj+1) - sshn_e(ji,jj) ) * hv(ji,jj) / e2v(ji,jj)
                  ! energy/enstrophy conserving formulation for planetary vorticity term
                  zcubt = + zfac25 / e1u(ji,jj) * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                     &                             + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
                  zcvbt = - zfac25 / e2v(ji,jj) * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                     &                             + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
                  ! after transports
                  ua_e(ji,jj) = ( zub_e(ji,jj) + z2dt_e * ( zcubt + zspgu + zua(ji,jj) ) ) * umask(ji,jj,1)
                  va_e(ji,jj) = ( zvb_e(ji,jj) + z2dt_e * ( zcvbt + zspgv + zva(ji,jj) ) ) * vmask(ji,jj,1)
               END DO
            END DO

         ENDIF

!!DB: add St. Lawrence River runoff to barotropic velocity
#ifdef key_RIVER_INPUT
         call riv_ts(ua_e,va_e)     !JC: General case for many rivers
#endif


         ! Flather's boundary condition for the barotropic loop :
         !         - Update sea surface height on each open boundary
         !         - Correct the barotropic transports
         IF( lk_obc )   CALL obc_fla_ts


         ! ... Boundary conditions on ua_e, va_e, ssha_e
         CALL lbc_lnk( ua_e  , 'U', -1. )
         CALL lbc_lnk( va_e  , 'V', -1. )
         CALL lbc_lnk( ssha_e, 'T',  1. )

         ! temporal sum
         !-------------
         zssha_b(:,:) = zssha_b(:,:) + ssha_e(:,:)
         zua_b  (:,:) = zua_b  (:,:) + ua_e  (:,:)
         zva_b  (:,:) = zva_b  (:,:) + va_e  (:,:) 

         ! Time filter and swap of dynamics arrays
         ! ---------------------------------------
         IF( neuler == 0 .AND. kt == nit000 ) THEN   ! Euler (forward) time stepping
            zsshb_e(:,:) = sshn_e(:,:)
            zub_e  (:,:) = zun_e (:,:)
            zvb_e  (:,:) = zvn_e (:,:)
            sshn_e (:,:) = ssha_e(:,:)
            zun_e  (:,:) = ua_e  (:,:)
            zvn_e  (:,:) = va_e  (:,:)
         ELSE                                        ! Asselin filtering
            zsshb_e(:,:) = atfp * ( zsshb_e(:,:) + ssha_e(:,:) ) + atfp1 * sshn_e(:,:)
            zub_e  (:,:) = atfp * ( zub_e  (:,:) + ua_e  (:,:) ) + atfp1 * zun_e  (:,:)
            zvb_e  (:,:) = atfp * ( zvb_e  (:,:) + va_e  (:,:) ) + atfp1 * zvn_e  (:,:)
            sshn_e (:,:) = ssha_e(:,:)
            zun_e  (:,:) = ua_e  (:,:)
            zvn_e  (:,:) = va_e  (:,:)
         ENDIF

         !                                               ! ==================== !
      END DO       !   END OF  do jit = 1, icycle        !  sub-time-step loop  !
      !                                                  ! ==================== !


      ! Time average of after barotropic variables
      zcoef =  1.e0 / (  FLOAT( icycle +1 )  )
      zssha_b(:,:) = zcoef * zssha_b(:,:) 
      zua_b  (:,:) = zcoef *  zua_b (:,:) 
      zva_b  (:,:) = zcoef *  zva_b (:,:) 
#if defined key_obc
         IF( lp_obc_east  )   sshfoe_b(:,:) = zcoef * sshfoe_b(:,:)
         IF( lp_obc_west  )   sshfow_b(:,:) = zcoef * sshfow_b(:,:)
         IF( lp_obc_north )   sshfon_b(:,:) = zcoef * sshfon_b(:,:)
         IF( lp_obc_south )   sshfos_b(:,:) = zcoef * sshfos_b(:,:)
#endif
     

      ! ---------------------------------------------------------------------------
      ! Phase 3 : Update sea surface height from time averaged barotropic variables
      ! ---------------------------------------------------------------------------

 
      ! Horizontal divergence of time averaged barotropic transports
      !-------------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            zhdiv(ji,jj) = ( e2u(ji,jj) * un_b(ji,jj) - e2u(ji-1,jj  ) * un_b(ji-1,jj  )     &
           &                +e1v(ji,jj) * vn_b(ji,jj) - e1v(ji  ,jj-1) * vn_b(ji  ,jj-1) )   &
           &             / ( e1t(ji,jj) * e2t(ji,jj) )
         END DO
      END DO

#if defined key_obc
      ! open boundaries (div must be zero behind the open boundary)
      !  mpp remark: The zeroing of hdiv can probably be extended to 1->jpi/jpj for the correct row/column
      IF( lp_obc_east  )   zhdiv(nie0p1:nie1p1,nje0  :nje1)   = 0.e0    ! east
      IF( lp_obc_west  )   zhdiv(niw0  :niw1  ,njw0  :njw1)   = 0.e0    ! west
      IF( lp_obc_north )   zhdiv(nin0  :nin1  ,njn0p1:njn1p1) = 0.e0    ! north
      IF( lp_obc_south )   zhdiv(nis0  :nis1  ,njs0  :njs1)   = 0.e0    ! south
#endif

      ! sea surface height
      !-------------------
      sshb(:,:) = sshn(:,:)
      sshn(:,:) = (  sshb_b(:,:) - z2dt_b * ( zraur * emp(:,:) + zhdiv(:,:) )  ) * tmask(:,:,1)

      ! ... Boundary conditions on sshn
      IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )


      ! -----------------------------------------------------------------------------
      ! Phase 4. Coupling between general trend and barotropic estimates - (2nd step)
      ! -----------------------------------------------------------------------------

      ! Swap on time averaged barotropic variables
      ! ------------------------------------------
      sshb_b(:,:) = sshn_b (:,:)
      sshn_b(:,:) = zssha_b(:,:)
      un_b  (:,:) = zua_b  (:,:) 
      vn_b  (:,:) = zva_b  (:,:) 
   
      ! add time averaged barotropic coriolis and surface pressure gradient
      ! terms to the general momentum trend
      ! --------------------------------------------------------------------

         ! Sea surface height from the barotropic system
      DO jk=1,jpkm1
         ua(:,:,jk) = ua(:,:,jk) + hur(:,:) * ( zua_b(:,:) - zub(:,:) ) / z2dt_b
         va(:,:,jk) = va(:,:,jk) + hvr(:,:) * ( zva_b(:,:) - zvb(:,:) ) / z2dt_b
      END DO



      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
      ENDIF
      
   END SUBROUTINE dyn_spg_ts
#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart free surface cst volume
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_spg_ts( kt )       ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_spg_ts: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_ts
#endif
   
   !!======================================================================
END MODULE dynspg_ts
