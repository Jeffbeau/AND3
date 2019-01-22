MODULE obcspg 
   !!======================================================================
   !!                       ***  MODULE  obcspg  ***
   !! Open Boundaries  :   Radiation of barotropic stream function on each
   !!                      open boundary
   !!======================================================================
#if   defined key_obc   &&   defined key_dynspg_rl
   !!----------------------------------------------------------------------
   !!   'key_obc'    and                            Open Boundary Condition
   !!   'key_dynspg_rl'                                 Rigid-Lid
   !!----------------------------------------------------------------------
   !!   obc_spg       : call the subroutine for each open boundary
   !!   obc_spg_east  : radiation of the east open boundary streamfunction
   !!   obc_spg_west  : radiation of the west open boundary streamfunction
   !!   obc_spg_north : radiation of the north open boundary streamfunction
   !!   obc_spg_south : radiation of the south open boundary streamfunction
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE phycst          ! physical constants
   USE obc_oce         ! ocean open boundary conditions
   USE lib_mpp         ! for mppobc
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC obc_spg     ! routine called in step.F90 (rigid lid case)

   !! * Module variables
   INTEGER ::   ji, jj, jk, jnic   ! dummy loop indices

   INTEGER ::      & ! ... boundary space indices 
      nib   = 1,   & ! nib   = boundary point
      nibm  = 2,   & ! nibm  = 1st interior point
      nibm2 = 3,   & ! nibm2 = 2nd interior point
      !              ! ... boundary time indices 
      nit   = 1,   & ! nit    = now
      nitm  = 2,   & ! nitm   = before
      nitm2 = 3      ! nitm2  = before-before

   REAL(wp) ::   rtaue  , rtauw  , rtaun  , rtaus  ,   &  !
                 rtauein, rtauwin, rtaunin, rtausin       !

   !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcspg.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_spg ( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE obc_spg  ***
      !!
      !! **  Purpose :
      !!       Compute now barotropic stream function at the open boundaries.
      !!       (lp_obc_east, and/or lp_obc_west, and/or lp_obc_north, and/or lp_obc_south).
      !!       Deduce the correct bsf trend on the open boundaries and isolated 
      !!       coastlines previous to the call to the barotropic solver.
      !!
      !! ** Method :
      !!      In case of open boundaries, there can be a net barotropic flow
      !!      through the boundaries, hence the potential on the coastlines
      !!      on each side of the OBC is different.
      !!      This routine:
      !!           1. compute the contribution of the isolated coastlines to the
      !!              rhs of the barotropic equation
      !!           2. compute the contribution of the OBC to the rhs of the
      !!              barotropic equation using a radiation equation as explained
      !!              in the OBC routine.
      !!
      !! Reference : 
      !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
      !! History :
      !!        ! 95-03 (J.-M. Molines) Original, SPEM
      !!        ! 97-07 (G. Madec, J.-M. Molines) additions
      !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt
      !!----------------------------------------------------------------------

      IF( kt == nit000 .OR. ln_rstart ) THEN      ! Initialization
         ! ... Boundary restoring coefficient
         rtaue = 2. * rdt / rdpeob
         rtauw = 2. * rdt / rdpwob
         rtaun = 2. * rdt / rdpnob
         rtaus = 2. * rdt / rdpsob
         ! ... Boundary restoring coefficient for inflow ( all boundaries)
         rtauein = 2. * rdt / rdpein 
         rtauwin = 2. * rdt / rdpwin
         rtaunin = 2. * rdt / rdpnin
         rtausin = 2. * rdt / rdpsin 
      ENDIF

      ! right hand side of the barotropic elliptic equation
      ! ---------------------------------------------------

      ! Isolated coastline contribution to the RHS of the barotropic Eq.
      gcbob(:,:) = 0.e0
      DO jnic = 1, nbobc-1
         gcbob(:,:) = gcbob(:,:) + gcfobc(:,:,jnic) * gcbic(jnic)
      END DO

      IF( lp_obc_east  )   CALL obc_spg_east ( kt )    ! East open boundary

      IF( lp_obc_west  )   CALL obc_spg_west ( kt )    ! West open boundary

      IF( lp_obc_north )   CALL obc_spg_north( kt )    ! North open boundary

      IF( lp_obc_south )   CALL obc_spg_south( kt )    ! South open boundary

      IF( lk_mpp )   CALL lbc_lnk( gcbob, 'G', 1. )
 
   END SUBROUTINE obc_spg


   SUBROUTINE obc_spg_east ( kt )
      !!------------------------------------------------------------------------------
      !!                ***  SUBROUTINE obc_spg_east  ***
      !!                 
      !! ** Purpose :   Apply the radiation algorithm on east OBC stream function.
      !!      If lfbceast=T , there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER ::   ij
      REAL(wp) ::   z2dtr, ztau, zin
      REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbceast is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart) .OR. lfbceast ) THEN

         ! 1.1 Fixed barotropic stream function
         ! ------------------------------------
         DO jj = nje0m1, nje1 
            ij = jj -1 + njmpp
            bsfeob(jj)=bfoe(ij)
         END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbceast is .FALSE.
      ! -----------------------------------------------------

         ! 2.1. Barotropic stream function radiation
         ! ----------------------------------------
         !
         !          nibm2      nibm      nib
         !            |   nibm  |   nib   |///
         !            |    |    |    |    |///
         !  jj-line --f----v----f----v----f---
         !            |         |         |///
         !            |    |    |    |    |///
         !         jpieob-2   jpieob-1    jpieob
         !                 |         |        
         !              jpieob-1    jpieob      
         !
         ! ... radiative conditions plus restoring term toward climatology
         ! ... Initialize bsfeob to clim in any case, at the first step
         !     to ensure proper values at the ends of the open line.
         ! ... Take care that the j indices starts at nje0 (jpjed) and finish 
         !     at nje1 (jpjef) to be sure that jpjefm1 and jpjef are set OK.
         DO ji = fs_nie0-1, fs_nie1-1 ! Vector opt.
            DO jj = nje0p1, nje1m2 
               ij = jj -1 + njmpp
         ! ... 2* gradi(bsf) (v-point i=nibm, time mean)
               z2dx = ( bebnd(jj,nibm ,nit) + bebnd(jj,nibm ,nitm2) - 2.*bebnd(jj,nibm2,nitm) ) &
                      / e1v(ji,jj)
         ! ... 2* gradj(bsf) (f-point i=nibm, time nitm)
               z2dy = ( bebnd(jj+1,nibm,nitm) - bebnd(jj-1,nibm,nitm) ) / e2f(ji,jj)
         ! ... square of the norm of grad(bsf)
               z4nor2 = z2dx * z2dx + z2dy * z2dy
         ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
               zdt = bebnd(jj,nibm,nitm2) - bebnd(jj,nibm,nit)
         ! ... i-phase speed ratio (bounded by 1) and MASKED!
               IF( z4nor2 == 0 ) THEN
                  IF(lwp) WRITE(numout,*)' PB dans obc_spg_east au pt ',jj,' : z4nor=0'
                  z4nor2 = 0.001
               ENDIF
               z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
               z05cx = z05cx / e1v(ji+1,jj)
               z05cx = min( z05cx, 1. )
         ! ... z05cx < 0, inflow  zin=0, ztau=1  
         !          => 0, outflow zin=1, ztau=rtaue
               zin = sign( 1., z05cx )
               zin = 0.5*( zin + abs(zin) )
         ! ... Modification JM:  We maintain a restoring term toward
         !                   bsfb even in case of inflow 
         ! But restoring is stronger when in flow (10 days) (ztau in set in obcspg.F)
               ztau = (1.-zin ) * rtauein + zin * rtaue
               z05cx = z05cx * zin
         ! ... update bsfn with radiative or climatological bsf (not mask!)
               bsfeob(jj) = ( ( 1. - z05cx - ztau ) * bebnd(jj,nib ,nitm) + 2.*z05cx  &
                               * bebnd(jj,nibm,nit) + ztau * bfoe (ij) )              &
                            / (1. + z05cx)
            END DO
         END DO

      ENDIF
      IF( lk_mpp )   CALL mppobc(bsfeob,jpjed,jpjef,jpieob-1,1,2,jpj)


      ! 3. right hand side of the barotropic elliptic equation
      ! ------------------------------------------------------
 
      IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
         z2dtr = 1.0 / rdt
      ELSE
         z2dtr = 0.5 / rdt
      ENDIF
      DO ji = fs_nie0-1, fs_nie1-1 ! Vector opt.
         DO jj = nje0m1, nje1 
            gcbob(ji,jj) = gcbob(ji,jj) - hvr(ji+1,jj) * e2v(ji+1,jj) / e1v(ji+1,jj) &
                           * ( bsfeob(jj) - bsfb(ji+1,jj) ) * z2dtr * bmask(ji,jj) 
         END DO
      END DO

   END SUBROUTINE obc_spg_east

   SUBROUTINE obc_spg_west ( kt )
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_spg_west  ***
      !!                    
      !! ** Purpose :
      !!      Apply the radiation algorithm on west OBC stream function.
      !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER ::   ij

      REAL(wp) ::   z2dtr, ztau, zin
      REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy

      !!------------------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcwest is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcwest ) THEN

         ! 1.1 Fixed barotropic stream function
         ! ------------------------------------
         DO jj = njw0m1, njw1
            ij = jj -1 + njmpp
            bsfwob(jj)=bfow(ij)
         END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
      ! -----------------------------------------------------

         ! 2.1. Barotropic stream function radiation
         ! ----------------------------------------
         !
         !         nib       nibm     nibm2
         !       ///|   nib   |   nibm  |
         !       ///|    |    |    |    |
         !       ---f----v----f----v----f-- jj-line
         !       ///|         |         |
         !       ///|    |    |    |    |
         !        jpiwob    jpiwob+1    jpiwob+2
         !               |         |        
         !             jpiwob+1    jpiwob+2     
         !
         ! ... radiative conditions plus restoring term toward climatology
         ! ... Initialize bsfwob to clim in any case, at the first step
         !     to ensure proper values at the ends of the open line.
         ! ... Take care that the j indices starts at njw0 (jpjwd) and finish
         ! ... at njw1 (jpjwf) to be sure that jpjwfm1 and jpjwf are set OK.
         DO ji = fs_niw0+1, fs_niw1+1 ! Vector opt.
            DO jj = njw0p1, njw1m2
               ij = jj -1 + njmpp
         ! ... 2* gradi(bsf) (v-point i=nibm, time mean)
               z2dx = ( -  bwbnd(jj,nibm ,nit ) - bwbnd(jj,nibm ,nitm2) + 2.*bwbnd(jj,nibm2,nitm ) ) &
                      / e1v(ji+1,jj)
         ! ... 2* gradj(bsf) (f-point i=nibm, time nitm)
               z2dy = ( bwbnd(jj+1,nibm,nitm) - bwbnd(jj-1,nibm,nitm) ) / e2f(ji,jj)
         ! ... square of the norm of grad(bsf)
               z4nor2 = z2dx * z2dx + z2dy * z2dy
         ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
               zdt = bwbnd(jj,nibm,nitm2) - bwbnd(jj,nibm,nit)
         ! ... i-phase speed ratio (bounded by 1) and MASKED!
               IF( z4nor2 == 0 ) THEN
                  IF(lwp) WRITE(numout,*)' PB dans obc_spg_west au pt ',jj,' : z4nor =0'
                  z4nor2=0.0001
               ENDIF
               z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
               z05cx = z05cx / e1v(ji,jj)
               z05cx = max( z05cx, -1. )
         ! ... z05cx => 0, inflow  zin=0, ztau=1  
         !           <  0, outflow zin=1, ztau=rtauw
               zin = sign( 1., -1. * z05cx )
               zin = 0.5*( zin + abs(zin) )
               ztau = (1.-zin )*rtauwin + zin * rtauw
               z05cx = z05cx * zin
         !  ... update bsfn with radiative or climatological bsf (not mask!)
               bsfwob(jj) = ( ( 1. + z05cx - ztau ) * bwbnd(jj,nib ,nitm) - 2.*z05cx &
                               * bwbnd(jj,nibm,nit) + ztau * bfow (ij) )             &
                            / (1. - z05cx)
            END DO
         END DO

      ENDIF
      IF( lk_mpp )   CALL mppobc(bsfwob,jpjwd,jpjwf,jpiwob+1,1,2,jpj) 


      ! 3. right hand side of the barotropic elliptic equation
      ! -------------------------------------------------------

      IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
         z2dtr = 1.0 / rdt
      ELSE
         z2dtr = 0.5 / rdt
      ENDIF
      DO ji = fs_niw0+1, fs_niw1+1 ! Vector opt.
         DO jj = njw0m1, njw1
            gcbob(ji,jj) = gcbob(ji,jj) - hvr(ji,jj) * e2v(ji,jj) / e1v(ji,jj)     &
                           * ( bsfwob(jj) - bsfb(ji-1,jj) ) * z2dtr * bmask(ji,jj)
         END DO
      END DO

   END SUBROUTINE obc_spg_west

   SUBROUTINE obc_spg_north ( kt )
      !!------------------------------------------------------------------------------
      !!                 ***  SUBROUTINE obc_spg_north  ***
      !! 
      !! ** Purpose :   Apply the radiation algorithm on north OBC stream function.
      !!      If lfbcnorth=T, there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt

      !! * Local declarations
      INTEGER ::   ii
      REAL(wp) ::   z2dtr, ztau, zin
      REAL(wp) ::   z05cx, zdt, z4nor2, z2dx, z2dy
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcnorth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------

      IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth ) THEN

         ! 1.1 Fixed barotropic stream function
         ! ------------------------------------
         DO ji = nin0m1, nin1
            ii = ji -1 + nimpp
            bsfnob(ji)=bfon(ii)
         END DO

      ELSE      

      ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
      ! -----------------------------------------------------

         ! 2.1. Barotropic stream function radiation
         ! -----------------------------------------
         !
         !           ji-row
         !             |
         !        ////////////
         !        ////////////
         !   nib  -----f------  jpjnob
         !             |    
         !      nib--  u   ---- jpjnob
         !             |        
         !  nibm  -----f-----   jpjnob-1
         !             |        
         !     nibm--  u   ---- jpjnob-1
         !             |        
         !  nibm2 -----f-----   jpjnob-2
         !             |
         !
         ! ... radiative conditions plus restoring term toward climatology
         ! ... z05cx is always the cross boundary phase velocity
         ! ... Initialize bsfnob to clim in any case, at the first step
         !     to ensure proper values at the ends of the open line.
         ! ... Take care that the i indices starts at nin0 (jpind) and finish
         ! ... at nin1 (jpinf) to be sure that jpinfm1 and jpinf are set OK.
         DO jj = fs_njn0-1, fs_njn1-1 ! Vector opt.
            DO ji = nin0p1, nin1m2
               ii = ji -1 + nimpp
         ! ... 2* gradj(bsf) (u-point i=nibm, time mean)
               z2dx = ( bnbnd(ji,nibm ,nit) + bnbnd(ji,nibm ,nitm2) - 2.*bnbnd(ji,nibm2,nitm) ) &
                      / e2u(ji,jj)
         ! ... 2* gradi(bsf) (f-point i=nibm, time nitm)
               z2dy = ( bnbnd(ji+1,nibm,nitm) - bnbnd(ji-1,nibm,nitm) ) / e1f(ji,jj)
         ! ... square of the norm of grad(bsf)
               z4nor2 = z2dx * z2dx + z2dy * z2dy
         ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
               zdt = bnbnd(ji,nibm,nitm2) - bnbnd(ji,nibm,nit)
         ! ... j-phase speed ratio (bounded by 1) and MASKED!
               IF( z4nor2 == 0 ) THEN
                  IF(lwp) WRITE(numout,*)' PB dans obc_spg_north au pt',ji,' : z4nor =0'
               ENDIF
               z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
               z05cx = z05cx / e2u(ji,jj+1)
               z05cx = min( z05cx, 1. )
         ! ... z05cx < 0, inflow  zin=0, ztau=1 
         !          => 0, outflow zin=1, ztau=rtaun
               zin = sign( 1., z05cx )
               zin = 0.5*( zin + abs(zin) )
               ztau = (1.-zin ) * rtaunin + zin * rtaun
               z05cx = z05cx * zin
         ! ... update bsfn with radiative or climatological bsf (not mask!)
               bsfnob(ji) = ( ( 1. - z05cx - ztau ) * bnbnd(ji,nib ,nitm) + 2.*z05cx  &
                               * bnbnd(ji,nibm,nit) + ztau * bfon (ii) )              &
                            / (1. + z05cx)
            END DO
         END DO

      ENDIF
      IF( lk_mpp )   CALL mppobc(bsfnob,jpind,jpinf,jpjnob-1,1,1,jpi)


      ! 3. right hand side of the barotropic elliptic equation
      !-------------------------------------------------------

      IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
         z2dtr = 1.0 / rdt
      ELSE
         z2dtr = 0.5 / rdt
      ENDIF
      DO jj = fs_njn0-1, fs_njn1-1 ! Vector opt.
         DO ji = nin0m1, nin1
            gcbob(ji,jj) = gcbob(ji,jj) - hur(ji,jj+1) *  e1u(ji,jj+1) / e2u(ji,jj+1) &
                           * ( bsfnob(ji) - bsfb(ji,jj+1) ) * z2dtr * bmask(ji,jj)
         END DO
      END DO

   END SUBROUTINE obc_spg_north


   SUBROUTINE obc_spg_south ( kt )
      !!------------------------------------------------------------------------------
      !!                  ***  SUBROUTINE obc_spg_south  ***
      !!                
      !! ** Purpose :   Apply the radiation algorithm on south OBC stream function.
      !!      If lfbcsouth=T, there is no radiation but only fixed OBC
      !!
      !!  History :
      !!         ! 95-03 (J.-M. Molines) Original from SPEM
      !!         ! 97-07 (G. Madec, J.-M. Molines) additions
      !!         ! 97-12 (M. Imbard) Mpp adaptation
      !!         ! 00-06 (J.-M. Molines) 
      !!    8.5  ! 02-10 (C. Talandier, A-M Treguier) F90
      !!------------------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in) ::   kt  ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ii             ! temporary integers
      REAL(wp) ::   &
         z2dtr, ztau, zin   ,   &  ! temporary scalars
         z05cx, zdt , z4nor2,   &  !    "         "
         z2dx , z2dy               !    "         "
      !!------------------------------------------------------------------------------

      ! 1. First three time steps and more if lfbcsouth is .TRUE.
      !    In that case open boundary conditions are FIXED.
      ! --------------------------------------------------------
      
      IF( ( kt < nit000 + 3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth ) THEN

         ! 1.1 Fixed barotropic stream function
         ! ------------------------------------
         DO ji = nis0m1, nis1 
            ii = ji -1 + nimpp 
            bsfsob(ji)=bfos(ii)
         END DO

      ELSE

      ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
      ! -----------------------------------------------------

         ! 2.1. Barotropic stream function radiation
         ! -----------------------------------------
         !
         !           ji-row
         !             |
         ! nibm2  -----f------  jpjsob + 2
         !             |    
         !   nibm  --  u  ----- jpjsob + 2
         !             |        
         !  nibm  -----f-----   jpjsob + 1
         !             |        
         !    nib  --  u  ----- jpjsob + 1
         !             |        
         !    nib -----f-----   jpjsob
         !        ///////////     
         !        ///////////
         !
         ! ... radiative conditions plus restoring term toward climatology
         ! ... z05cx is always the cross boundary phase velocity
         ! ... Initialize bsfsob to clim in any case, at the first step
         !     to ensure proper values at the ends of the open line.
         ! ... Take care that the i indices starts at nis0 (jpisd) and finish
         ! ... at nis1 (jpisf) to be sure that jpisfm1 and jpisf are set OK.
         DO jj = fs_njs0+1, fs_njs1+1 ! Vector opt.
            DO ji = nis0p1, nis1m2
               ii = ji -1 + nimpp
         ! ... 2* gradj(bsf) (u-point i=nibm, time mean)
               z2dx = ( - bsbnd(ji,nibm ,nit) - bsbnd(ji,nibm ,nitm2) + 2.*bsbnd(ji,nibm2,nitm) ) &
                      / e2u(ji,jj+1)
         ! ... 2* gradi(bsf) (f-point i=nibm, time nitm)
               z2dy = ( bsbnd(ji+1,nibm,nitm) - bsbnd(ji-1,nibm,nitm) ) / e1f(ji,jj)
         ! ... square of the norm of grad(bsf)
               z4nor2 = z2dx * z2dx + z2dy * z2dy
         ! ... minus time derivative (leap-frog) at nibm, without / 2 dt
               zdt = bsbnd(ji,nibm,nitm2) - bsbnd(ji,nibm,nit)
         ! ... j-phase speed ratio (bounded by -1) and MASKED!
               IF( z4nor2 == 0 ) THEN
                  IF(lwp) WRITE(numout,*)' PB dans obc_spg_south au pt ',ji,' : z4nor =0'
               ENDIF
               z05cx = zdt * z2dx / z4nor2 * bmask(ji,jj)
               z05cx = z05cx / e2u(ji,jj)
               z05cx = max( z05cx, -1. )
         ! ... z05cx => 0, inflow  zin=0, ztau=1
         !           <  0, outflow zin=1, ztau=rtaus
               zin = sign( 1., -1. * z05cx )
               zin = 0.5*( zin + abs(zin) )
               ztau = (1.-zin ) *rtausin  + zin * rtaus
               z05cx = z05cx * zin
         ! ... update bsfn with radiative or climatological bsf (not mask!)
               bsfsob(ji) = ( ( 1. + z05cx - ztau ) * bsbnd(ji,nib ,nitm) - 2.*z05cx  & 
                               * bsbnd(ji,nibm,nit) + ztau * bfos (ii) )              &
                            / (1. - z05cx)
            END DO
         END DO

      ENDIF
      IF( lk_mpp )   CALL mppobc(bsfsob,jpisd,jpisf,jpjsob+1,1,1,jpi)

 
      ! 3. right hand side of the barotropic elliptic equation
      ! -------------------------------------------------------

      IF( ( neuler == 0 ) .AND. ( kt == nit000 ) ) THEN
         z2dtr = 1.0 / rdt
      ELSE
         z2dtr = 0.5 / rdt
      ENDIF
      DO jj = fs_njs0+1, fs_njs1+1 ! Vector opt.
         DO ji = nis0m1, nis1 
            gcbob(ji,jj) = gcbob(ji,jj) - hur(ji,jj) * e1u(ji,jj) / e2u(ji,jj) &
                           * ( bsfsob(ji) - bsfb(ji,jj-1) ) * z2dtr * bmask(ji,jj)
         END DO
      END DO

   END SUBROUTINE obc_spg_south

#else
   !!----------------------------------------------------------------------
   !!   Default case :                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_spg( kt )        ! Empty routine
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'obc_spg: You should not have seen this print! error?', kt
   END SUBROUTINE obc_spg
#endif

   !!======================================================================
END MODULE obcspg
