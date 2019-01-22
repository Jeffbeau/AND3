MODULE eosbn2
   !!==============================================================================
   !!                       ***  MODULE  eosbn2  ***
   !! Ocean diagnostic variable : equation of state - in situ and potential density
   !!                                               - Brunt-Vaisala frequency 
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   eos            : generic interface of the equation of state
   !!   eos_insitu     : Compute the in situ density
   !!   eos_insitu_pot : Compute the insitu and surface referenced potential
   !!                    volumic mass
   !!   eos_insitu_2d  : Compute the in situ density for 2d fields
   !!   eos_bn2        : Compute the Brunt-Vaisala frequency
   !!   eos_init       : set eos parameters (namelist)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE zdfddm          ! vertical physics: double diffusion
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Interface 
   INTERFACE eos
      MODULE PROCEDURE eos_insitu, eos_insitu_pot, eos_insitu_2d
   END INTERFACE 
   INTERFACE bn2
      MODULE PROCEDURE eos_bn2
   END INTERFACE 

   !! * Routine accessibility
   PUBLIC eos        ! called by step.F90, inidtr.F90, tranpc.F90 and intgrd.F90
   PUBLIC bn2        ! called by step.F90

   !! * Share module variables
   INTEGER , PUBLIC ::   &  !: nameos : ocean physical parameters
      neos      = 0,     &  !: = 0/1/2 type of eq. of state and Brunt-Vaisala frequ.
      neos_init = 0         !: control flag for initialization

   REAL(wp), PUBLIC ::   &  !: nameos : ocean physical parameters
      ralpha = 2.0e-4,   &  !: thermal expension coeff. (linear equation of state)
      rbeta  = 7.7e-4       !: saline  expension coeff. (linear equation of state)
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/eosbn2.F90,v 1.9 2005/09/02 15:45:19 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE eos_insitu ( ptem, psal, prd )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE eos_insitu  ***
      !! 
      !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
      !!       potential temperature and salinity using an equation of state
      !!       defined through the namelist parameter neos.
      !!
      !! ** Method  :   3 cases:
      !!      neos = 0 : Jackett and McDougall (1994) equation of state.
      !!         the in situ density is computed directly as a function of
      !!         potential temperature relative to the surface (the opa t
      !!         variable), salt and pressure (assuming no pressure variation
      !!         along geopotential surfaces, i.e. the pressure p in decibars
      !!         is approximated by the depth in meters.
      !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
      !!         with pressure                      p        decibars
      !!              potential temperature         t        deg celsius
      !!              salinity                      s        psu
      !!              reference volumic mass        rau0     kg/m**3
      !!              in situ volumic mass          rho      kg/m**3
      !!              in situ density anomalie      prd      no units
      !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
      !!          t = 40 deg celcius, s=40 psu
      !!      neos = 1 : linear equation of state function of temperature only
      !!              prd(t) = 0.0285 - ralpha * t
      !!      neos = 2 : linear equation of state function of temperature and
      !!               salinity
      !!              prd(t,s) = rbeta * s - ralpha * tn - 1.
      !!      Note that no boundary condition problem occurs in this routine
      !!      as (ptem,psal) are defined over the whole domain.
      !!
      !! ** Action  :   compute prd , the in situ density (no units)
      !!
      !! References :
      !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
      !!
      !! History :
      !!        !  89-03 (o. Marti)  Original code
      !!        ! 94-08 (G. Madec)
      !!        !  96-01 (G. Madec) statement function for e3
      !!        !  97-07 (G. Madec) introduction of neos, OPA8.1
      !!        !  97-07 (G. Madec) density instead of volumic mass
      !!        !  99-02 (G. Madec, N. Grima) semi-implicit pressure gradient
      !!        !  01-09 (M. Ben Jelloul) bugfix   
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         ptem,                 &  ! potential temperature
         psal                     ! salinity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
         prd                      ! potential density (surface referenced)

      !! * Local declarations
      INTEGER ::  ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   &          
         zt , zs , zh , zsr,   &  ! temporary scalars
         zr1, zr2, zr3, zr4,   &  !    "         " 
         zrhop, ze, zbw, zb,   &  !    "         "
         zd , zc , zaw, za ,   &  !    "         "
         zb1, za1, zkw, zk0       !    "         "
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zws                      ! temporary workspace
      !!----------------------------------------------------------------------


      ! initialization (in not already done)
      IF( neos_init == 0 ) CALL eos_init


      SELECT CASE ( neos )

      CASE ( 0 )               ! Jackett and McDougall (1994) formulation

!CDIR NOVERRCHK
         zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  ! depth
                  zh = fsdept(ji,jj,jk)
                  ! square root salinity
                  zsr= zws(ji,jj,jk)
                  ! compute volumic mass pure water at atm pressure
                  zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                     -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
                  ! seawater volumic mass atm pressure
                  zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
                     -4.0899e-3 ) *zt+0.824493
                  zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
                  zr4= 4.8314e-4

                  ! potential volumic mass (reference to the surface)
                  zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

                  ! add the compression terms
                  ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
                  zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
                  zb = zbw + ze * zs

                  zd = -2.042967e-2
                  zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
                  zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
                  za = ( zd*zsr + zc ) *zs + zaw

                  zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
                  za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
                  zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
                  zk0= ( zb1*zsr + za1 )*zs + zkw

                  ! masked in situ density anomaly
                  prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
                     - rau0 ) / rau0 * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 1 )               ! Linear formulation function of temperature only

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  !   ... density and potential volumic mass
                  prd(ji,jj,jk) = ( 0.0285 - ralpha * zt ) * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  !   ... density and potential volumic mass
                  prd(ji,jj,jk) = (   rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab3d_1=prd, clinfo1=' eos  : ', ovlap=1, kdim=jpk)
      ENDIF

   END SUBROUTINE eos_insitu


   SUBROUTINE eos_insitu_pot ( ptem, psal, prd, prhop)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_insitu_pot  ***
      !!           
      !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
      !!      potential volumic mass (Kg/m3) from potential temperature and
      !!      salinity fields using an equation of state defined through the 
      !!     namelist parameter neos.
      !!
      !! ** Method  :
      !!      neos = 0 : Jackett and McDougall (1994) equation of state.
      !!         the in situ density is computed directly as a function of
      !!         potential temperature relative to the surface (the opa t
      !!         variable), salt and pressure (assuming no pressure variation
      !!         along geopotential surfaces, i.e. the pressure p in decibars
      !!         is approximated by the depth in meters.
      !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
      !!              rhop(t,s)  = rho(t,s,0)
      !!         with pressure                      p        decibars
      !!              potential temperature         t        deg celsius
      !!              salinity                      s        psu
      !!              reference volumic mass        rau0     kg/m**3
      !!              in situ volumic mass          rho      kg/m**3
      !!              in situ density anomalie      prd      no units
      !!
      !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
      !!          t = 40 deg celcius, s=40 psu
      !!
      !!      neos = 1 : linear equation of state function of temperature only
      !!              prd(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
      !!              rhop(t,s)  = rho(t,s)
      !!
      !!      neos = 2 : linear equation of state function of temperature and
      !!               salinity
      !!              prd(t,s) = ( rho(t,s) - rau0 ) / rau0 
      !!                       = rbeta * s - ralpha * tn - 1.
      !!              rhop(t,s)  = rho(t,s)
      !!      Note that no boundary condition problem occurs in this routine
      !!      as (tn,sn) or (ta,sa) are defined over the whole domain.
      !!
      !! ** Action  : - prd  , the in situ density (no units)
      !!              - prhop, the potential volumic mass (Kg/m3)
      !!
      !! References :
      !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
      !!      Brown, J. A. and K. A. Campana. Mon. Weather Rev., 1978
      !!
      !! History :
      !!   4.0  !  89-03  (O. Marti)
      !!        !  94-08  (G. Madec)
      !!        !  96-01  (G. Madec) statement function for e3
      !!        !  97-07  (G. Madec) introduction of neos, OPA8.1
      !!        !  97-07  (G. Madec) density instead of volumic mass
      !!        !  99-02  (G. Madec, N. Grima) semi-implicit pressure gradient
      !!        !  01-09  (M. Ben Jelloul) bugfix   
      !!   9.0  !  03-08  (G. Madec)  F90, free form
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         ptem,   &  ! potential temperature
         psal       ! salinity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
         prd,    &  ! potential density (surface referenced)
         prhop      ! potential density (surface referenced)

      !! * Local declarations
      INTEGER  ::  ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   &             ! temporary scalars
         zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
         zb, zd, zc, zaw, za, zb1, za1, zkw, zk0
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zws
      !!----------------------------------------------------------------------

      ! initialization (in not already done)
      IF( neos_init == 0 ) CALL eos_init


      SELECT CASE ( neos )

      CASE ( 0 )               ! Jackett and McDougall (1994) formulation

!CDIR NOVERRCHK
         zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  ! depth
                  zh = fsdept(ji,jj,jk)
                  ! square root salinity
!!Edmee           zsr= SQRT( ABS( zs ) )
                  zsr= zws(ji,jj,jk)
                  ! compute volumic mass pure water at atm pressure
                  zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                     -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
                  ! seawater volumic mass atm pressure
                  zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
                     -4.0899e-3 ) *zt+0.824493
                  zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
                  zr4= 4.8314e-4

                  ! potential volumic mass (reference to the surface)
                  zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

                  ! save potential volumic mass
                  prhop(ji,jj,jk) = zrhop * tmask(ji,jj,jk)

                  ! add the compression terms
                  ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
                  zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
                  zb = zbw + ze * zs

                  zd = -2.042967e-2
                  zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
                  zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
                  za = ( zd*zsr + zc ) *zs + zaw

                  zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
                  za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
                  zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
                  zk0= ( zb1*zsr + za1 )*zs + zkw

                  ! masked in situ density anomaly
                  prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
                     - rau0 ) / rau0 * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 1 )               ! Linear formulation function of temperature only

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  !   ... density and potential volumic mass
                  prd  (ji,jj,jk) = ( 0.0285 - ralpha * zt )        * tmask(ji,jj,jk)
                  prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 ) * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zt = ptem(ji,jj,jk)
                  zs = psal(ji,jj,jk)
                  !   ... density and potential volumic mass
                  prd  (ji,jj,jk) = ( rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
                  prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 )    * tmask(ji,jj,jk)
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab3d_1=prd, clinfo1=' eos-p: ', tab3d_2=prhop, clinfo2=' pot : ', ovlap=1, kdim=jpk)
      ENDIF

   END SUBROUTINE eos_insitu_pot

   SUBROUTINE eos_insitu_2d ( ptem, psal, pdep, prd )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_insitu_2d  ***
      !!
      !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
      !!      potential temperature and salinity using an equation of state
      !!      defined through the namelist parameter neos. * 2D field case
      !!
      !! ** Method :
      !!      neos = 0 : Jackett and McDougall (1994) equation of state.
      !!         the in situ density is computed directly as a function of
      !!         potential temperature relative to the surface (the opa t
      !!         variable), salt and pressure (assuming no pressure variation
      !!         along geopotential surfaces, i.e. the pressure p in decibars
      !!         is approximated by the depth in meters.
      !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
      !!         with pressure                      p        decibars
      !!              potential temperature         t        deg celsius
      !!              salinity                      s        psu
      !!              reference volumic mass        rau0     kg/m**3
      !!              in situ volumic mass          rho      kg/m**3
      !!              in situ density anomalie      prd      no units
      !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
      !!          t = 40 deg celcius, s=40 psu
      !!      neos = 1 : linear equation of state function of temperature only
      !!              prd(t) = 0.0285 - ralpha * t
      !!      neos = 2 : linear equation of state function of temperature and
      !!               salinity
      !!              prd(t,s) = rbeta * s - ralpha * tn - 1.
      !!      Note that no boundary condition problem occurs in this routine
      !!      as (ptem,psal) are defined over the whole domain.
      !!
      !! ** Action  : - prd , the in situ density (no units)
      !!
      !! References :
      !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
      !!
      !! History :
      !!   8.5  !  02-11  (G. Madec, A. Bozec)  partial step
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT( in ) ::   &
         ptem,                           &  ! potential temperature
         psal,                           &  ! salinity
         pdep                               ! depth
      REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   &
         prd                                ! potential density (surface referenced)

      !! * Local declarations
      INTEGER ::  ji, jj                    ! dummy loop indices
      REAL(wp) ::   &             ! temporary scalars
         zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
         zb, zd, zc, zaw, za, zb1, za1, zkw, zk0,               &
         zmask
      REAL(wp), DIMENSION(jpi,jpj) :: zws
      !!----------------------------------------------------------------------


      ! initialization (in not already done)
      IF( neos_init == 0 ) CALL eos_init

      prd(:,:) = 0.e0

      SELECT CASE ( neos )

      CASE ( 0 )               ! Jackett and McDougall (1994) formulation

!CDIR NOVERRCHK
         DO jj = 1, jpjm1
!CDIR NOVERRCHK
#if defined key_autotasking
            DO ji = 1, jpim1
#else
            DO ji = 1, fs_jpim1   ! vector opt.
#endif
               zws(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
            END DO
         END DO

         !                                                ! ===============
         DO jj = 1, jpjm1                                 ! Horizontal slab
            !                                             ! ===============
#if defined key_autotasking
            DO ji = 1, jpim1
#else
            DO ji = 1, fs_jpim1   ! vector opt.
#endif

               zmask = tmask(ji,jj,1)      ! land/sea bottom mask = surf. mask

               zt = ptem (ji,jj)            ! interpolated T
               zs = psal (ji,jj)            ! interpolated S
               zsr= zws(ji,jj)            ! square root of interpolated S
               zh = pdep(ji,jj)            ! depth at the partial step level

               ! compute volumic mass pure water at atm pressure
               zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
                         -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
               ! seawater volumic mass atm pressure
               zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
                         -4.0899e-3 ) *zt+0.824493
               zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
               zr4= 4.8314e-4

               ! potential volumic mass (reference to the surface)
               zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

               ! add the compression terms
               ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
               zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
               zb = zbw + ze * zs

               zd = -2.042967e-2
               zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
               zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt -4.721788
               za = ( zd*zsr + zc ) *zs + zaw

               zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
               za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
               zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt   &
                         +2098.925 ) *zt+190925.6
               zk0= ( zb1*zsr + za1 )*zs + zkw

               ! masked in situ density anomaly
               prd(ji,jj) = ( zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rau0 )   &
                          / rau0 * zmask
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 1 )               ! Linear formulation function of temperature only

         !                                                ! ===============
         DO jj = 1, jpjm1                                 ! Horizontal slab
            !                                             ! ===============
#if defined key_autotasking
            DO ji = 1, jpim1
#else
            DO ji = 1, fs_jpim1   ! vector opt.
#endif
               prd(ji,jj) = ( 0.0285 - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1)
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         !                                                ! ===============
         DO jj = 1, jpjm1                                 ! Horizontal slab
            !                                             ! ===============
#if defined key_autotasking
            DO ji = 1, jpim1
#else
            DO ji = 1, fs_jpim1   ! vector opt.
#endif
               prd(ji,jj) = ( rbeta * psal(ji,jj) - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1) 
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

      IF(ln_ctl)   CALL prt_ctl(tab2d_1=prd, clinfo1=' eos2d: ')

   END SUBROUTINE eos_insitu_2d


   SUBROUTINE eos_bn2( ptem, psal, pn2 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE eos_bn2  ***
      !!
      !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the time-
      !!      step of the input arguments
      !!      
      !! ** Method :
      !!       * neos = 0  : UNESCO sea water properties
      !!         The brunt-vaisala frequency is computed using the polynomial
      !!      polynomial expression of McDougall (1987):
      !!            N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
      !!      If lk_zdfddm=T, the heat/salt buoyancy flux ratio Rrau is
      !!      computed and used in zdfddm module :
      !!              Rrau = alpha/beta * ( dk[ t ] / dk[ s ] )
      !!       * neos = 1  : linear equation of state (temperature only)
      !!            N^2 = grav * ralpha * dk[ t ]/e3w
      !!       * neos = 2  : linear equation of state (temperature & salinity)
      !!            N^2 = grav * (ralpha * dk[ t ] - rbeta * dk[ s ] ) / e3w
      !!      The use of potential density to compute N^2 introduces e r r o r
      !!      in the sign of N^2 at great depths. We recommand the use of 
      !!      neos = 0, except for academical studies.
      !!        Macro-tasked on horizontal slab (jk-loop)
      !!      N.B. N^2 is set to zero at the first level (JK=1) in inidtr
      !!      and is never used at this level.
      !!
      !! ** Action  : - pn2 : the brunt-vaisala frequency
      !!
      !! References :
      !!      McDougall, T. J., J. Phys. Oceanogr., 17, 1950-1964, 1987.
      !!
      !! History :
      !!   6.0  !  94-07  (G. Madec, M. Imbard)  Original code
      !!   8.0  !  97-07  (G. Madec) introduction of statement functions
      !!   8.5  !  02-07  (G. Madec) Free form, F90
      !!   8.5  !  02-08  (G. Madec) introduction of arguments
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
         ptem,                           &  ! potential temperature
         psal                               ! salinity
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
         pn2                               ! Brunt-Vaisala frequency

      !! * Local declarations
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   &
         zgde3w, zt, zs, zh,  &  ! temporary scalars 
         zalbet, zbeta           !    "         "
#if defined key_zdfddm
      REAL(wp) ::   zds          ! temporary scalars
#endif
      !!----------------------------------------------------------------------
      !!  OPA8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      ! pn2 : first and last levels
      ! ---------------------------
      ! bn^2=0. at jk=1 and jpk set in inidtr.F : no computation


      ! pn2 : interior points only (2=< jk =< jpkm1 )
      ! -------------------------- 

      SELECT CASE ( neos )

      CASE ( 0 )               ! Jackett and McDougall (1994) formulation

         !                                                ! ===============
         DO jk = 2, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zgde3w = grav / fse3w(ji,jj,jk)
                  zt = 0.5 * ( ptem(ji,jj,jk) + ptem(ji,jj,jk-1) )          ! potential temperature at w-point
                  zs = 0.5 * ( psal(ji,jj,jk) + psal(ji,jj,jk-1) ) - 35.0   ! salinity anomaly (s-35) at w-point
                  zh = fsdepw(ji,jj,jk)                                     ! depth in meters  at w-point

                  zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
                     &                               - 0.203814e-03 ) * zt   &
                     &                               + 0.170907e-01 ) * zt   &
                     &   + 0.665157e-01                                      &
                     &   +     ( - 0.678662e-05 * zs                         &
                     &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
                     &   +   ( ( - 0.302285e-13 * zh                         &
                     &           - 0.251520e-11 * zs                         &
                     &           + 0.512857e-12 * zt * zt           ) * zh   &
                     &           - 0.164759e-06 * zs                         &
                     &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
                     &                               + 0.380374e-04 ) * zh

                  zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
                     &                            - 0.301985e-05 ) * zt      &
                     &   + 0.785567e-03                                      &
                     &   + (     0.515032e-08 * zs                           &
                     &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
                     &   +(  (   0.121551e-17 * zh                           &
                     &         - 0.602281e-15 * zs                           &
                     &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
                     &                             + 0.408195e-10   * zs     &
                     &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
                     &                             - 0.121555e-07 ) * zh

                  pn2(ji,jj,jk) = zgde3w * zbeta * tmask(ji,jj,jk)           &   ! N^2 
                     &          * ( zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &
                     &                     - ( psal(ji,jj,jk-1) - psal(ji,jj,jk) ) )
#if defined key_zdfddm
                  !                                                         !!bug **** caution a traiter zds=dk[S]= 0 !!!!
                  zds = ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )                    ! Rrau = (alpha / beta) (dk[t] / dk[s])
                  IF ( ABS( zds) <= 1.e-20 ) zds = 1.e-20
                  rrau(ji,jj,jk) = zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) ) / zds
#endif
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 1 )               ! Linear formulation function of temperature only

         !                                                ! ===============
         DO jk = 2, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zgde3w = grav / fse3w(ji,jj,jk) * tmask(ji,jj,jk)
                  pn2(ji,jj,jk) = zgde3w * ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============


      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         !                                                ! ===============
         DO jk = 2, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zgde3w = grav / fse3w(ji,jj,jk) * tmask(ji,jj,jk)
                  pn2(ji,jj,jk) = zgde3w * (  ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &
                     &                      - rbeta  * ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )  )
               END DO
            END DO
#if defined key_zdfddm
            !                                                ! Rrau = (alpha / beta) (dk[t] / dk[s])
            zalbet = ralpha / rbeta
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zds = ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )
                  IF ( ABS( zds ) <= 1.e-20 ) zds = 1.e-20
                  rrau(ji,jj,jk) = zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) ) / zds
               END DO
            END DO
#endif
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab3d_1=pn2, clinfo1=' bn2  : ', ovlap=1, kdim=jpk)
#if defined key_zdfddm
         CALL prt_ctl(tab3d_1=rrau, clinfo1=' rrau : ', ovlap=1, kdim=jpk)
#endif
      ENDIF

   END SUBROUTINE eos_bn2


   SUBROUTINE eos_init
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE eos_init  ***
      !!
      !! ** Purpose :   initializations for the equation of state
      !!
      !! ** Method  :   Read the namelist nameos
      !!
      !! ** Action  :   blahblah....
      !!
      !! History :
      !!   8.5  !  02-10  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      NAMELIST/nameos/ neos, ralpha, rbeta
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      ! set the initialization flag to 1
      neos_init = 1           ! indicate that the initialization has been done

      ! namelist nameos : ocean physical parameters

      ! Read Namelist nameos : equation of state
      REWIND( numnam )
      READ  ( numnam, nameos )

      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'eos_init : equation of state'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '          Namelist nameos : set eos parameters'
         WRITE(numout,*)
         WRITE(numout,*) '             flag for eq. of state and N^2  neos   = ', neos
         WRITE(numout,*) '             thermal exp. coef. (linear)    ralpha = ', ralpha
         WRITE(numout,*) '             saline  exp. coef. (linear)    rbeta  = ', rbeta
         WRITE(numout,*)
      ENDIF

      SELECT CASE ( neos )

      CASE ( 0 )               ! Jackett and McDougall (1994) formulation

         IF(lwp) WRITE(numout,*) '          use of Jackett & McDougall (1994) equation of state and'
         IF(lwp) WRITE(numout,*) '                 McDougall (1987) Brunt-Vaisala frequency'

      CASE ( 1 )               ! Linear formulation function of temperature only

         IF(lwp) WRITE(numout,*) '          use of linear eos rho(T) = rau0 * ( 1.0285 - ralpha * T )'
         IF( lk_zdfddm ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          double diffusive mixing parameterization requires',   &
                                             ' that T and S are used as state variables'
            nstop = nstop + 1
         ENDIF

      CASE ( 2 )               ! Linear formulation function of temperature and salinity

         IF(lwp) WRITE(numout,*) '          use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'

      CASE DEFAULT

         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
         nstop = nstop + 1

      END SELECT

   END SUBROUTINE eos_init

   !!======================================================================
END MODULE eosbn2  
