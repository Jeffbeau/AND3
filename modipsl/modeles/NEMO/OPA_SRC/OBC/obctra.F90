MODULE obctra
  !!=================================================================================
  !!                       ***  MODULE  obctra  ***
  !! Ocean tracers:   Radiation of tracers on each open boundary
  !!=================================================================================
#if defined key_obc
  !!---------------------------------------------------------------------------------
  !!   'key_obc'      :                                      Open Boundary Conditions
  !!---------------------------------------------------------------------------------
  !!   obc_tra        : call the subroutine for each open boundary
  !!   obc_tra_east   : radiation of the east open boundary tracers
  !!   obc_tra_west   : radiation of the west open boundary tracers
  !!   obc_tra_north  : radiation of the north open boundary tracers
  !!   obc_tra_south  : radiation of the south open boundary tracers
  !!----------------------------------------------------------------------------------
  !! * Modules used
  USE oce             ! ocean dynamics and tracers variables
  USE dom_oce         ! ocean space and time domain variables
  USE phycst          ! physical constants
  USE obc_oce         ! ocean open boundary conditions
  USE lib_mpp         ! ???
  USE lbclnk          ! ???
  USE in_out_manager  ! I/O manager

  IMPLICIT NONE
  PRIVATE

  !! * Accessibility
  PUBLIC obc_tra     ! routine called in tranxt.F90

  !! * Module variables
  INTEGER ::      & ! ... boundary space indices
     nib   = 1,   & ! nib   = boundary point
     nibm  = 2,   & ! nibm  = 1st interior point
     nibm2 = 3,   & ! nibm2 = 2nd interior point
                    ! ... boundary time indices
     nit   = 1,   & ! nit    = now
     nitm  = 2,   & ! nitm   = before
     nitm2 = 3      ! nitm2  = before-before

  REAL(wp) ::     &
     rtaue  , rtauw  , rtaun  , rtaus  ,  &  ! Boundary restoring coefficient
     rtauein, rtauwin, rtaunin, rtausin      ! Boundary restoring coefficient for inflow

  !! * Substitutions
#  include "obc_vectopt_loop_substitute.h90"
  !!---------------------------------------------------------------------------------
  !!   OPA 9.0 , LOCEAN-IPSL (2005)
  !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obctra.F90,v 1.4 2005/03/27 18:35:10 opalod Exp $
  !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE obc_tra( kt )
     !!-------------------------------------------------------------------------------
     !!                 ***  SUBROUTINE obc_tra  ***
     !!
     !! ** Purpose :   Compute tracer fields (t,s) along the open boundaries.
     !!      This routine is called by the tranxt.F routine and updates ta,sa
     !!      which are the actual temperature and salinity fields.
     !!        The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
     !!      and/or lp_obc_south allow the user to determine which boundary is an
     !!      open one (must be done in the param_obc.h90 file).
     !!
     !! Reference :
     !!   Marchesiello P., 1995, these de l'universite J. Fourier, Grenoble, France.
     !!
     !!  History :
     !!        !  95-03 (J.-M. Molines) Original, SPEM
     !!        !  97-07 (G. Madec, J.-M. Molines) addition
     !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) F90
     !!----------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt
     !!----------------------------------------------------------------------

     ! 0. Local constant initialization

     IF( kt == nit000 .OR. ln_rstart) THEN
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
     END IF

     IF( lp_obc_east  )   CALL obc_tra_east ( kt )    ! East open boundary

     IF( lp_obc_west  )   CALL obc_tra_west ( kt )    ! West open boundary

     IF( lp_obc_north )   CALL obc_tra_north( kt )    ! North open boundary

     IF( lp_obc_south )   CALL obc_tra_south( kt )    ! South open boundary

     IF( lk_mpp ) THEN                  !!bug ???
        IF( kt >= nit000+3 .AND. ln_rstart ) THEN
           CALL lbc_lnk( tb, 'T', 1. )
           CALL lbc_lnk( sb, 'T', 1. )
        END IF
        CALL lbc_lnk( ta, 'T', 1. )
        CALL lbc_lnk( sa, 'T', 1. )
     ENDIF

  END SUBROUTINE obc_tra


  SUBROUTINE obc_tra_east ( kt )
     !!------------------------------------------------------------------------------
     !!                ***  SUBROUTINE obc_tra_east  ***
     !!
     !! ** Purpose :
     !!      Apply the radiation algorithm on east OBC tracers ta, sa using the
     !!      phase velocities calculated in obc_rad_east subroutine in obcrad.F90 module
     !!      If the logical lfbceast is .TRUE., there is no radiation but only fixed OBC
     !!
     !!  History :
     !!         ! 95-03 (J.-M. Molines) Original from SPEM
     !!         ! 97-07 (G. Madec, J.-M. Molines) additions
     !!         ! 97-12 (M. Imbard) Mpp adaptation
     !!         ! 00-06 (J.-M. Molines)
     !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
     !!------------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt

     !! * Local declaration
     INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
     REAL(wp) ::   z05cx, ztau, zin
!!DB
     INTEGER :: jji, jjj, jjk

     !!------------------------------------------------------------------------------

     ! 1. First three time steps and more if lfbceast is .TRUE.
     !    In that case open boundary conditions are FIXED.
     ! --------------------------------------------------------

     IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbceast ) THEN
        DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                 ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-temsk(jj,jk)) + &
                               temsk(jj,jk) * tfoe(jj,jk,1)
                 sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-temsk(jj,jk)) + &
                               temsk(jj,jk) * sfoe(jj,jk,1)
              END DO
           END DO
        END DO

!!DB DBG TS sponge; 5 cells; timescale = 1day(?)
!!
        zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
!        zin = 2. * rdt / (20.0*rday)
!        do jji = nie0, nie1   !!isolate processor 
!DBG -- skip 
        do jji = 1, 0 
           do ji = nie0-4, nie0-1 !!5 cells - 1
!restore based on position -- but do not use
              ztau = zin * exp(float( -(nie0-ji)/4 )) 
              do jk = 1, jpkm1
                 do jj = 1, jpj
                    ta(ji,jj,jk) = ta(ji,jj,jk) + temsk(jj,jk)*&
                         zin*(tfoe(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
                    sa(ji,jj,jk) = sa(ji,jj,jk) + temsk(jj,jk)*&
                         zin*(sfoe(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
                 enddo
              enddo
           enddo
        enddo







     ELSE

     ! 2. Beyond the fourth time step if lfbceast is .FALSE.
     ! -----------------------------------------------------

        ! Temperature and salinity radiation
        ! ----------------------------------
        !
        !            nibm2      nibm      nib
        !              |   nibm  |   nib///|///
        !              |    |    |    |////|///
        !  jj   line --v----f----v----f----v---
        !              |    |    |    |////|///
        !                   |         |///   //
        !  jj   line   T    u    T    u/// T //
        !                   |         |///   //
        !              |    |    |    |////|///
        !  jj-1 line --v----f----v----f----v---
        !              |    |    |    |////|///
        !                jpieob-1    jpieob / ///
        !              |         |         |
        !           jpieob-1    jpieob     jpieob+1
        !
        ! ... radiative conditions + relaxation toward a climatology
        !     the phase velocity is taken as the phase velocity of the tangen-
        !     tial velocity (here vn), which have been saved in (u_cxebnd,v_cxebnd)
        ! ... (jpjedp1, jpjefm1), jpieob+1
        DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 2, jpjm1
        ! ... i-phase speed ratio (from averaged of v_cxebnd)
                 z05cx = ( 0.5 * ( v_cxebnd(jj,jk) + v_cxebnd(jj-1,jk) ) ) / e1t(ji-1,jj)
                 z05cx = min( z05cx, 1. )
        ! ... z05cx=< 0, inflow  zin=0, ztau=1
        !           > 0, outflow zin=1, ztau=rtaue
                 zin = sign( 1., z05cx )
                 zin = 0.5*( zin + abs(zin) )
        ! ... for inflow rtauein is used for relaxation coefficient else rtaue
                 ztau = (1.-zin ) * rtauein  + zin * rtaue
                 z05cx = z05cx * zin
        ! ... update ( ta, sa ) with radiative or climatological (t, s)
                 ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
                                temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
                                * tebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
!byoung                                 * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk) ) &
                                * tebnd(ji,jk,nibm,nit ) + ztau * tfoe (ji,jk,1) ) &
                                / (1. + z05cx)
                 sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - temsk(jj,jk)) +           &
                                temsk(jj,jk) * ( ( 1. - z05cx - ztau )         &
                                * sebnd(jj,jk,nib ,nitm) + 2.*z05cx              &
!byoung                                * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk) ) &
                                * sebnd(jj,jk,nibm,nit ) + ztau * sfoe (jj,jk,1) ) &
                                / (1. + z05cx)
              END DO
           END DO
        END DO

     END IF

  END SUBROUTINE obc_tra_east


  SUBROUTINE obc_tra_west ( kt )
     !!------------------------------------------------------------------------------
     !!                 ***  SUBROUTINE obc_tra_west  ***
     !!
     !! ** Purpose :
     !!      Apply the radiation algorithm on west OBC tracers ta, sa using the
     !!      phase velocities calculated in obc_rad_west subroutine in obcrad.F90 module
     !!      If the logical lfbcwest is .TRUE., there is no radiation but only fixed OBC
     !!
     !!  History :
     !!         ! 95-03 (J.-M. Molines) Original from SPEM
     !!         ! 97-07 (G. Madec, J.-M. Molines) additions
     !!         ! 97-12 (M. Imbard) Mpp adaptation
     !!         ! 00-06 (J.-M. Molines)
     !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
     !!------------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt

     !! * Local declaration
     INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
     REAL(wp) ::   z05cx, ztau, zin
!!DB
     INTEGER :: jji, jjj, jjk

     !!------------------------------------------------------------------------------

     ! 1. First three time steps and more if lfbcwest is .TRUE.
     !    In that case open boundary conditions are FIXED.
     ! --------------------------------------------------------

     IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcwest ) THEN

        DO ji = fs_niw0, fs_niw1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                 ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
                               twmsk(jj,jk) * tfow(jj,jk,1)
                 sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-twmsk(jj,jk)) + &
                               twmsk(jj,jk) * sfow(jj,jk,1)
              END DO
           END DO
        END DO
!!DB DBG TS sponge; 5 cells; timescale = 1day(?)
!!
        zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
!        zin = 2. * rdt / (20.0*rday)
!        do jji = niw0, niw1   !!isolate processor 
        do jji = 1, 0    !!isolate processor 
           do ji = niw0+1, niw0+4 !!5 cells - 1
!!restore based on position -- but do not use
              ztau = zin * exp(float( -(ji-niw0)/4 )) 
              do jk = 1, jpkm1
                 do jj = 1, jpj
                    ta(ji,jj,jk) = ta(ji,jj,jk) + twmsk(jj,jk)*&
                         zin*(tfow(jj,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
                    sa(ji,jj,jk) = sa(ji,jj,jk) + twmsk(jj,jk)*&
                         zin*(sfow(jj,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
                 enddo
              enddo
           enddo
        enddo

     ELSE

     ! 2. Beyond the fourth time step if lfbcwest is .FALSE.
     ! -----------------------------------------------------

        ! Temperature and salinity radiation
        ! ----------------------------------
        !
        !          nib       nibm     nibm2
        !     nib///|   nibm  |  nibm2  |
        !   ///|////|    |    |    |    |
        !   ---v----f----v----f----v----f-- jj   line
        !   ///|////|    |    |    |    |
        !   //   ///|         |         |
        !   // T ///u    T    u    T    u   jj   line
        !   //   ///|         |         |
        !   ///|////|    |    |    |    |
        !   ---v----f----v----f----v----f-- jj-1 line
        !   ///|////|    |    |    |    |
        !         jpiwob    jpiwob+1    jpiwob+2
        !      |         |         |
        !    jpiwob    jpiwob+1   jpiwob+2
        !
        ! ... radiative conditions + relaxation toward a climatology
        ! ... the phase velocity is taken as the phase velocity of the tangen-
        ! ... tial velocity (here vn), which have been saved in (v_cxwbnd)
        DO ji = fs_niw0, fs_niw1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 2, jpjm1
        ! ... i-phase speed ratio (from averaged of v_cxwbnd)
                 z05cx = (  0.5 * ( v_cxwbnd(jj,jk) + v_cxwbnd(jj-1,jk) ) ) / e1t(ji+1,jj)
                 z05cx = max( z05cx, -1. )
        ! ... z05cx > 0, inflow  zin=0, ztau=1
        !           < 0, outflow zin=1, ztau=rtauw
                 zin = sign( 1., -1.* z05cx )
                 zin = 0.5*( zin + abs(zin) )
                 ztau = (1.-zin )*rtauwin + zin * rtauw
                 z05cx = z05cx * zin
        ! ... update (ta,sa) with radiative or climatological (t, s)
                 ta(ji,jj,jk) = ta(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
                                twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
                                * twbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
!byoung                                * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk) ) &
                                * twbnd(jj,jk,nibm,nit ) + ztau * tfow (jj,jk,1) ) &
                                / (1. - z05cx)
                 sa(ji,jj,jk) = sa(ji,jj,jk) * (1. - twmsk(jj,jk)) +           &
                                twmsk(jj,jk) * ( ( 1. + z05cx - ztau )         &
                                * swbnd(jj,jk,nib ,nitm) - 2.*z05cx              &
!byoung                                * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk) ) &
                                * swbnd(jj,jk,nibm,nit ) + ztau * sfow (jj,jk,1) ) &
                                / (1. - z05cx)
              END DO
           END DO
        END DO

     END IF

  END SUBROUTINE obc_tra_west


  SUBROUTINE obc_tra_north ( kt )
     !!------------------------------------------------------------------------------
     !!                 ***  SUBROUTINE obc_tra_north  ***
     !!
     !! ** Purpose :
     !!      Apply the radiation algorithm on north OBC tracers ta, sa using the
     !!      phase velocities calculated in obc_rad_north subroutine in obcrad.F90 module
     !!      If the logical lfbcnorth is .TRUE., there is no radiation but only fixed OBC
     !!
     !!  History :
     !!         ! 95-03 (J.-M. Molines) Original from SPEM
     !!         ! 97-07 (G. Madec, J.-M. Molines) additions
     !!         ! 97-12 (M. Imbard) Mpp adaptation
     !!         ! 00-06 (J.-M. Molines)
     !!    8.5  ! 02-10 (C. Talandier, A-M. Treguier) F90
     !!------------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt

     !! * Local declaration
!sujie      INTEGER ::   ji, jj, jk      ! dummy loop indices
     INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
     REAL(wp) ::   z05cx, ztau, zin
     !!------------------------------------------------------------------------------

     ! 1. First three time steps and more if lfbcnorth is .TRUE.
     !    In that case open boundary conditions are FIXED.
     ! --------------------------------------------------------

     IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcnorth ) THEN

!sujie
        DO jj = fs_njn0+1, fs_njn1+1  ! Vector opt.
           DO jk = 1, jpkm1
              DO ji = 1, jpi
!                  ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
!                                tnmsk(ji,jk) * tfon(ji,jk)
!                  sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
!                                tnmsk(ji,jk) * sfon(ji,jk)
             if(sfon(ji,jk,1).gt.0.) then
                 ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
                               tnmsk(ji,jk) * tfon(ji,jk,1)
                 sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) + &
                               tnmsk(ji,jk) * sfon(ji,jk,1)
             endif
!             do j1=1,4
!             if(sfon(ji,jk,j1+1).gt.0.) then
!             ta(ji,jj-j1,jk)= ta(ji,jj-j1,jk)  - &
!                          (ta(ji,jj-j1,jk)  - tfon(ji,jk,j1+1)) &
!                         * rdt/(86400.* exp(float(j1+1))/15.0) &
!                         * tnmsk5(ji,jk,j1+1)
!             sa(ji,jj-j1,jk)= sa(ji,jj-j1,jk)  - &
!                          (sa(ji,jj-j1,jk)  - sfon(ji,jk,j1+1)) &
!                         * rdt/(86400.* exp(float(j1+1))/15.0) &
!                         * tnmsk5(ji,jk,j1+1)
!             endif
!             end do

              END DO
           END DO
        END DO

!ylu
        j1=0
        DO jj = mj0(jpjglo-2), mj1(jpjglo-5),-1
           j1=j1+1
           DO jk = 1, jpkm1
              DO ji = 2, jpi-1
            if(sfon(ji,jk,j1+1).gt.0.) then
            ta(ji,jj,jk)= ta(ji,jj,jk)  - &
                         (ta(ji,jj,jk)  - tfon(ji,jk,j1+1)) &
                        * rdt/(86400.* exp(float(j1+1))/15.0) &
                        * tnmsk5(ji,jk,j1+1)
            sa(ji,jj,jk)= sa(ji,jj,jk)  - &
                         (sa(ji,jj,jk)  - sfon(ji,jk,j1+1)) &
                        * rdt/(86400.* exp(float(j1+1))/15.0) &
                        * tnmsk5(ji,jk,j1+1)
            endif
              END DO
           END DO
        END DO

     ELSE

     ! 2. Beyond the fourth time step if lfbcnorth is .FALSE.
     ! -------------------------------------------------------

        ! Temperature and salinity radiation
        ! ----------------------------------
        !
        !           ji-1   ji   ji   ji +1
        !             |
        !    nib //// u // T // u // T //   jpjnob + 1
        !        /////|//////////////////
        !    nib  ----f----v----f----v---   jpjnob
        !             |         |
        !      nibm-- u -- T -- u -- T --   jpjnob
        !             |         |
        !   nibm  ----f----v----f----v---  jpjnob-1
        !             |         |
        !     nibm2-- u -- T -- T -- T --  jpjnob-1
        !             |         |
        !   nibm2 ----f----v----f----v---  jpjnob-2
        !             |         |
        !
        ! ... radiative conditions + relaxation toward a climatology
        ! ... the phase velocity is taken as the normal phase velocity of the tangen-
        ! ... tial velocity (here un), which has been saved in (u_cynbnd)
        ! ... jpjnob+1,(jpindp1, jpinfm1)
        DO jj = fs_njn0+1, fs_njn1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO ji = 2, jpim1
        ! ... j-phase speed ratio (from averaged of vtnbnd)
        !        (bounded by 1)
                 z05cx = ( 0.5 * ( u_cynbnd(ji,jk) + u_cynbnd(ji-1,jk) ) ) / e2t(ji,jj-1)
                 z05cx = min( z05cx, 1. )
        ! ... z05cx=< 0, inflow  zin=0, ztau=1
        !           > 0, outflow zin=1, ztau=rtaun
                 zin = sign( 1., z05cx )
                 zin = 0.5*( zin + abs(zin) )
        ! ... for inflow rtaunin is used for relaxation coefficient else rtaun
                 ztau = (1.-zin ) * rtaunin + zin * rtaun
                 z05cx = z05cx * zin
        ! ... update (ta,sa) with radiative or climatological (t, s)
                 ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
                                tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
                                * tnbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
!sujie                                 * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk) ) &
                                * tnbnd(ji,jk,nibm,nit ) + ztau * tfon (ji,jk,1) ) &
                                / (1. + z05cx)
                 sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tnmsk(ji,jk)) +             &
                                tnmsk(ji,jk) * ( ( 1. - z05cx - ztau )         &
                                * snbnd(ji,jk,nib ,nitm) + 2.*z05cx              &
!sujie                                 * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk) ) &
                                * snbnd(ji,jk,nibm,nit ) + ztau * sfon (ji,jk,1) ) &
                                / (1. + z05cx)
              END DO
           END DO
        END DO

     END IF

  END SUBROUTINE obc_tra_north


  SUBROUTINE obc_tra_south ( kt )
     !!------------------------------------------------------------------------------
     !!                ***  SUBROUTINE obc_tra_south  ***
     !!
     !! ** Purpose :
     !!      Apply the radiation algorithm on south OBC tracers ta, sa using the
     !!      phase velocities calculated in obc_rad_south subroutine in obcrad.F90 module
     !!      If the logical lfbcsouth is .TRUE., there is no radiation but only fixed OBC
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

     !! * Local declaration
     INTEGER ::   ji, jj, jk, j1      ! dummy loop indices
     REAL(wp) ::   z05cx, ztau, zin
!!DB
     INTEGER :: jji, jjj, jjk

     !!------------------------------------------------------------------------------

     ! 1. First three time steps and more if lfbcsouth is .TRUE.
     !    In that case open boundary conditions are FIXED.
     ! --------------------------------------------------------

     IF( ( kt < nit000+3 .AND. .NOT.ln_rstart ) .OR. lfbcsouth ) THEN


        DO jj = fs_njs0, fs_njs1  ! Vector opt.
           DO jk = 1, jpkm1
              DO ji = 1, jpi
                 ta(ji,jj,jk)= ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
                               tsmsk(ji,jk) * tfos(ji,jk,1)
                 sa(ji,jj,jk)= sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) + &
                               tsmsk(ji,jk) * sfos(ji,jk,1)
              END DO
           END DO
        END DO

!!DB DBG velocity sponge; 5 cells; timescale = 1day(?)
!!
        zin = 2. * rdt / (1.0*rday)   !!blows up < 30dt
!        zin = 2. * rdt / (20.0*rday)
!        do jji = njs0, njs1   !!isolate processor 
        do jji = 1, 0   !!isolate processor 
           do jj = njs0+1, njs0+4 !!5 cells - 1
!!restore based on position -- but do not use
              ztau = zin * exp(float( -(jj-njs0)/4 )) 
!!DBG
!              zin = ztau
              do jk = 1, jpkm1
                 do ji = 1, jpi
                    ta(ji,jj,jk) = ta(ji,jj,jk) + tsmsk(ji,jk)*&
                         zin*(tfos(ji,jk,1)-ta(ji,jj,jk)) * tmask(ji,jj,jk)
                    sa(ji,jj,jk) = sa(ji,jj,jk) + tsmsk(ji,jk)*&
                         zin*(sfos(ji,jk,1)-sa(ji,jj,jk)) * tmask(ji,jj,jk)
                 enddo
              enddo
           enddo
        enddo



     ELSE

     ! 2. Beyond the fourth time step if lfbcsouth is .FALSE.
     ! -------------------------------------------------------

        ! Temperature and salinity radiation
        ! ----------------------------------
        !
        !           ji-1   ji   ji   ji +1
        !             |         |
        !   nibm2 ----f----v----f----v---   jpjsob+2
        !             |         |
        !   nibm2 --  u -- T -- u -- T --   jpjsob+2
        !             |         |
        !   nibm  ----f----v----f----v---   jpjsob+1
        !             |         |
        !    nibm --  u -- T -- T -- T --   jpjsob+1
        !             |         |
        !   nib  -----f----v----f----v---   jpjsob
        !       //////|/////////|////////
        !    nib //// u // T // u // T //   jpjsob
        !
        !... radiative conditions + relaxation toward a climatology
        !... the phase velocity is taken as the phase velocity of the tangen-
        !... tial velocity (here un), which has been saved in (u_cysbnd)
        !... jpjsob,(jpisdp1, jpisfm1)
        DO jj = fs_njs0, fs_njs1  ! Vector opt.
           DO jk = 1, jpkm1
              DO ji = 2, jpim1
        !... j-phase speed ratio (from averaged of u_cysbnd)
        !       (bounded by 1)
                 z05cx = ( 0.5 * ( u_cysbnd(ji,jk) + u_cysbnd(ji-1,jk) ) ) / e2t(ji,jj+1)
                 z05cx = max( z05cx, -1. )
        !... z05cx > 0, inflow  zin=0, ztau=1
        !          < 0, outflow zin=1, ztau=rtaus
                 zin = sign( 1., -1.* z05cx )
                 zin = 0.5*( zin + abs(zin) )
                 ztau = (1.-zin ) + zin * rtaus
                 z05cx = z05cx * zin
        !... update (ta,sa) with radiative or climatological (t, s)
                 ta(ji,jj,jk) = ta(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
                                tsmsk(ji,jk) * ( ( 1. + z05cx - ztau )         &
                                * tsbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
!sujie                                 * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk) ) &
                                * tsbnd(ji,jk,nibm,nit ) + ztau * tfos (ji,jk,1) ) &
                                / (1. - z05cx)
                 sa(ji,jj,jk) = sa(ji,jj,jk) * (1.-tsmsk(ji,jk)) +             &
                                tsmsk(ji,jk) * (  ( 1. + z05cx - ztau )        &
                                * ssbnd(ji,jk,nib ,nitm) - 2.*z05cx              &
!sujie                                 * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk) ) &
                                * ssbnd(ji,jk,nibm,nit ) + ztau * sfos (ji,jk,1) ) &
                                / (1. - z05cx)
              END DO
           END DO
        END DO

     END IF

  END SUBROUTINE obc_tra_south

#else
  !!---------------------------------------------------------------------------------
  !!   Default option                                                    Empty module
  !!---------------------------------------------------------------------------------
CONTAINS
  SUBROUTINE obc_tra      ! Empty routine
  END SUBROUTINE obc_tra
#endif

  !!=================================================================================
END MODULE obctra
