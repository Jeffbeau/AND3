MODULE dynspg_ts_jki
   !!======================================================================
   !!                   ***  MODULE  dynspg_ts_jki  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
#if ( defined key_dynspg_ts && defined key_autotasking )   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_dynspg_ts'                    free surface with time splitting
   !!   'key_autotasking'                          j-k-i loop (vector opt.)
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

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_ts_jki  ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LODYC-IPSL  (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_ts_jki.F90,v 1.2 2006/01/03 15:04:14 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_spg_ts_jki( kt )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_jki  ***
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
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts_jki : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   free surface with time splitting with j-k-i loop'
         IF(lwp) WRITE(numout,*) ' Number of sub cycle in 1 time-step (2 rdt) : icycle = ', FLOOR( 2*rdt/rdtbt )

         IF( .NOT. ln_rstart ) THEN
            ! initialize barotropic specific arrays
            sshb_b(:,:) = sshb(:,:)
            sshn_b(:,:) = sshn(:,:)
            un_b(:,:)   = 0.e0
            vn_b(:,:)   = 0.e0
            ! vertical sum
            DO jk = 1, jpkm1
               un_b(:,:) = un_b(:,:) + fse3u(:,:,jk) * un(:,:,jk)
               vn_b(:,:) = vn_b(:,:) + fse3v(:,:,jk) * vn(:,:,jk)
            END DO
         ENDIF
         ssha_e(:,:) = sshn(:,:)
         ua_e  (:,:) = un_b(:,:)
         va_e  (:,:) = vn_b(:,:)

         IF( ln_dynvor_een ) THEN
            ztne(1,:) = 0.e0   ;   ztnw(1,:) = 0.e0   ;   ztse(1,:) = 0.e0   ;   ztsw(1,:) = 0.e0
            DO jj = 2, jpj
               DO ji = 2, jpi
                  ztne(ji,jj) = ( ff(ji-1,jj  ) + ff(ji  ,jj  ) + ff(ji  ,jj-1) ) / 3.
                  ztnw(ji,jj) = ( ff(ji-1,jj-1) + ff(ji-1,jj  ) + ff(ji  ,jj  ) ) / 3.
                  ztse(ji,jj) = ( ff(ji  ,jj  ) + ff(ji  ,jj-1) + ff(ji-1,jj-1) ) / 3.
                  ztsw(ji,jj) = ( ff(ji  ,jj-1) + ff(ji-1,jj-1) + ff(ji-1,jj  ) ) / 3.
               END DO
            END DO
         ENDIF

      ENDIF
    
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

      DO jj = 1, jpj

         ! variables for the barotropic equations
         zsshb_e(:,jj) = sshn_b(:,jj)       ! (barotropic) sea surface height (before and now)
         sshn_e (:,jj) = sshn_b(:,jj)
         zub_e  (:,jj) = un_b  (:,jj)       ! barotropic transports issued from the barotropic equations (before and now)
         zvb_e  (:,jj) = vn_b  (:,jj)
         zun_e  (:,jj) = un_b  (:,jj)
         zvn_e  (:,jj) = vn_b  (:,jj)
         zssha_b(:,jj) = sshn  (:,jj)        ! time averaged variables over all sub-timesteps
         zua_b  (:,jj) = un_b  (:,jj)   
         zva_b  (:,jj) = vn_b  (:,jj)

         ! Vertically integrated quantities
         ! --------------------------------
         zua(:,jj) = 0.e0
         zva(:,jj) = 0.e0
         zub(:,jj) = 0.e0
         zvb(:,jj) = 0.e0
         zwx(:,jj) = 0.e0
         zwy(:,jj) = 0.e0

         ! vertical sum
         DO jk = 1, jpkm1
            !                                                           ! Vertically integrated momentum trends
            zua(:,jj) = zua(:,jj) + fse3u(:,jj,jk) * umask(:,jj,jk) * ua(:,jj,jk)
            zva(:,jj) = zva(:,jj) + fse3v(:,jj,jk) * vmask(:,jj,jk) * va(:,jj,jk)
            !                                                           ! Vertically integrated transports (before)
            zub(:,jj) = zub(:,jj) + fse3u(:,jj,jk) * ub(:,jj,jk)
            zvb(:,jj) = zvb(:,jj) + fse3v(:,jj,jk) * vb(:,jj,jk)
            !                                                           ! Planetary vorticity (now)
            zwx(:,jj) = zwx(:,jj) + e2u(:,jj) * fse3u(:,jj,jk) * un(:,jj,jk)
            zwy(:,jj) = zwy(:,jj) + e1v(:,jj) * fse3v(:,jj,jk) * vn(:,jj,jk)
         END DO

      END DO

      DO jj = 2, jpjm1

         IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
            DO ji = 2, jpim1
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               ! energy conserving formulation for planetary vorticity term
               zcu(ji,jj) = zfact2 * ( ff(ji  ,jj-1) * zy1 + ff(ji,jj) * zy2 )
               zcv(ji,jj) =-zfact2 * ( ff(ji-1,jj  ) * zx1 + ff(ji,jj) * zx2 )
            END DO

         ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
            DO ji = 2, jpim1
               zy1 = zfact1 * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
                              + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 =-zfact1 * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
                              + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               zcu(ji,jj)  = zy1 * ( ff(ji  ,jj-1) + ff(ji,jj) )
               zcv(ji,jj)  = zx1 * ( ff(ji-1,jj  ) + ff(ji,jj) )
            END DO

         ELSEIF ( ln_dynvor_een ) THEN                    ! enstrophy and energy conserving scheme
         zfac25 = 0.25
            DO ji = 2, jpim1
               zcu(ji,jj) = + zfac25 / e1u(ji,jj)   &
                  &       * (  ztne(ji,jj  ) * zwy(ji  ,jj  ) + ztnw(ji+1,jj) * zwy(ji+1,jj  )   &
                  &          + ztse(ji,jj  ) * zwy(ji  ,jj-1) + ztsw(ji+1,jj) * zwy(ji+1,jj-1) )
               zcv(ji,jj) = - zfac25 / e2v(ji,jj)   &
                  &       * (  ztsw(ji,jj+1) * zwx(ji-1,jj+1) + ztse(ji,jj+1) * zwx(ji  ,jj+1)   &
                  &          + ztnw(ji,jj  ) * zwx(ji-1,jj  ) + ztne(ji,jj  ) * zwx(ji  ,jj  ) )
            END DO

         ENDIF


         ! Remove barotropic trend from general momentum trend
         DO jk = 1 , jpkm1
            DO ji = 2, jpim1
               ua(ji,jj,jk) = ua(ji,jj,jk) - zua(ji,jj) * hur(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - zva(ji,jj) * hvr(ji,jj)
            END DO
         END DO

         ! Remove coriolis term from barotropic trend
         ! ------------------------------------------
         DO ji = 2, jpim1
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

         DO jj = 2, jpjm1

         ! Horizontal divergence of barotropic transports
         !--------------------------------------------------
            DO ji = 2, jpim1
               zhdiv(ji,jj) = ( e2u(ji  ,jj  ) * zun_e(ji  ,jj)              &
                  &            -e2u(ji-1,jj  ) * zun_e(ji-1,jj)              &
                  &            +e1v(ji  ,jj  ) * zvn_e(ji  ,jj)              &
                  &            -e1v(ji  ,jj-1) * zvn_e(ji  ,jj-1) )          &
                  &           / (e1t(ji,jj)*e2t(ji,jj))
            END DO

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  ) THEN
            IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
         ENDIF
         IF( lp_obc_west  ) THEN
            IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
         ENDIF
         IF( lp_obc_north ) THEN
            IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
         ENDIF
         IF( lp_obc_south ) THEN
            IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
         ENDIF
#endif

         ! Sea surface height from the barotropic system
         !----------------------------------------------
            DO ji = 2, jpim1
               ssha_e(ji,jj) = ( zsshb_e(ji,jj) - z2dt_e *  ( zraur * emp(ji,jj)  &
                  &            +  zhdiv(ji,jj) ) ) * tmask(ji,jj,1)
            END DO

         END DO

         ! evolution of the barotropic transport ( following the vorticity scheme used)
         ! ----------------------------------------------------------------------------
         zwx(:,:) = e2u(:,:) * zun_e(:,:)
         zwy(:,:) = e1v(:,:) * zvn_e(:,:)

         DO jj = 2, jpjm1

            IF( ln_dynvor_ene .OR. ln_dynvor_mix ) THEN      ! energy conserving or mixed scheme
               DO ji = 2, jpim1
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

            ELSEIF ( ln_dynvor_ens ) THEN                    ! enstrophy conserving scheme
               DO ji = 2, jpim1
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

            ELSEIF ( ln_dynvor_een ) THEN                    ! energy and enstrophy conserving scheme
               zfac25 = 0.25
               DO ji = 2, jpim1
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
            ENDIF

         END DO


         ! Flather's boundary condition for the barotropic loop :
         !         - Update sea surface height on each open boundary
         !         - Correct the barotropic transports
         IF( lk_obc )   CALL obc_fla_ts

         ! ... Boundary conditions on ua_e, va_e, ssha_e
         CALL lbc_lnk( ua_e  , 'U', -1. )
         CALL lbc_lnk( va_e  , 'V', -1. )
         CALL lbc_lnk( ssha_e, 'T',  1. )

         DO jj = 1, jpj

            ! temporal sum
            !-------------
            zssha_b(:,jj) = zssha_b(:,jj) + ssha_e(:,jj)
            zua_b  (:,jj) = zua_b  (:,jj) + ua_e  (:,jj)
            zva_b  (:,jj) = zva_b  (:,jj) + va_e  (:,jj) 

            ! Time filter and swap of dynamics arrays
            ! ---------------------------------------
            IF( neuler == 0 .AND. kt == nit000 ) THEN   ! Euler (forward) time stepping
               zsshb_e(:,jj) = sshn_e(:,jj)
               zub_e  (:,jj) = zun_e (:,jj)
               zvb_e  (:,jj) = zvn_e (:,jj)
               sshn_e (:,jj) = ssha_e(:,jj)
               zun_e  (:,jj) = ua_e  (:,jj)
               zvn_e  (:,jj) = va_e  (:,jj)
            ELSE                                        ! Asselin filtering
               zsshb_e(:,jj) = atfp * ( zsshb_e(:,jj) + ssha_e(:,jj) ) + atfp1 * sshn_e (:,jj)
               zub_e  (:,jj) = atfp * ( zub_e  (:,jj) + ua_e  (:,jj) ) + atfp1 * zun_e  (:,jj)
               zvb_e  (:,jj) = atfp * ( zvb_e  (:,jj) + va_e  (:,jj) ) + atfp1 * zvn_e  (:,jj)
               sshn_e (:,jj) = ssha_e(:,jj)
               zun_e  (:,jj) = ua_e  (:,jj)
               zvn_e  (:,jj) = va_e  (:,jj)
            ENDIF

         END DO

         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !


      ! Time average of after barotropic variables
      zcoef =  1.e0 / (  FLOAT( icycle +1 )  )
      zssha_b(:,:) = zcoef * zssha_b(:,:) 
      zua_b  (:,:) = zcoef * zua_b  (:,:) 
      zva_b  (:,:) = zcoef * zva_b  (:,:) 
#if defined key_obc
         IF( lp_obc_east  )   sshfoe_b(:,:) = zcoef * sshfoe_b(:,:)
         IF( lp_obc_west  )   sshfow_b(:,:) = zcoef * sshfow_b(:,:)
         IF( lp_obc_north )   sshfon_b(:,:) = zcoef * sshfon_b(:,:)
         IF( lp_obc_south )   sshfos_b(:,:) = zcoef * sshfos_b(:,:)
#endif
     

      ! ---------------------------------------------------------------------------
      ! Phase 3 : Update sea surface height from time averaged barotropic variables
      ! ---------------------------------------------------------------------------

      sshb(:,:) = sshn(:,:)
 
      ! Horizontal divergence of time averaged barotropic transports
      !-------------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zhdiv(ji,jj) = ( e2u(ji,jj) * un_b(ji,jj) - e2u(ji-1,jj  ) * un_b(ji-1,jj  )     &
           &                +e1v(ji,jj) * vn_b(ji,jj) - e1v(ji  ,jj-1) * vn_b(ji  ,jj-1) )   &
           &             / ( e1t(ji,jj) * e2t(ji,jj) )
         END DO

#if defined key_obc
         ! open boundaries (div must be zero behind the open boundary)
         !  mpp remark: The zeroing of zhdiv can probably be extended to 1->jpi/jpj for the correct row/column
         IF( lp_obc_east  ) THEN
            IF( nje0   <= jj .AND. jj <= nje1   )   zhdiv(nie0p1:nie1p1,jj) = 0.e0      ! east
         ENDIF
         IF( lp_obc_west  ) THEN
            IF( njw0   <= jj .AND. jj <= njw1   )   zhdiv(niw0  :niw1  ,jj) = 0.e0      ! west
         ENDIF
         IF( lp_obc_north ) THEN
            IF( njn0p1 <= jj .AND. jj <= njn1p1 )   zhdiv(nin0  :nin1  ,jj) = 0.e0      ! north
         ENDIF
         IF( lp_obc_south ) THEN
            IF( njs0   <= jj .AND. jj <= njs1   )   zhdiv(nis0  :nis1  ,jj) = 0.e0      ! south
         ENDIF
#endif

         ! sea surface height
         !-------------------
         DO ji = 2, jpim1
            sshn(ji,jj) = (  sshb_b(ji,jj) - z2dt_b * ( zraur * emp(ji,jj) + zhdiv(ji,jj) )  ) * tmask(ji,jj,1)
         END DO

      END DO

      ! ... Boundary conditions on sshn
      IF( .NOT. lk_obc ) CALL lbc_lnk( sshn, 'T', 1. )


      ! -----------------------------------------------------------------------------
      ! Phase 4. Coupling between general trend and barotropic estimates - (2nd step)
      ! -----------------------------------------------------------------------------

      DO jj = 1, jpj

         ! Swap on time averaged barotropic variables
         ! ------------------------------------------
         sshb_b(:,jj) = sshn_b (:,jj)
         sshn_b(:,jj) = zssha_b(:,jj)
         un_b  (:,jj) = zua_b  (:,jj) 
         vn_b  (:,jj) = zva_b  (:,jj) 
   
         ! add time averaged barotropic coriolis and surface pressure gradient
         ! terms to the general momentum trend
         ! --------------------------------------------------------------------
         DO jk = 1, jpkm1
            ua(:,jj,jk) = ua(:,jj,jk) + hur(:,jj) * ( zua_b(:,jj) - zub(:,jj) ) / z2dt_b
            va(:,jj,jk) = va(:,jj,jk) + hvr(:,jj) * ( zva_b(:,jj) - zvb(:,jj) ) / z2dt_b
         END DO
      END DO

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
      ENDIF
      
   END SUBROUTINE dyn_spg_ts_jki
#else
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart free surface cst volume
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_spg_ts_jki( kt )       ! Empty routine
      if(lwp) WRITE(numout,*) 'dyn_spg_ts_jki: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_ts_jki
#endif
   
   !!======================================================================
END MODULE dynspg_ts_jki
