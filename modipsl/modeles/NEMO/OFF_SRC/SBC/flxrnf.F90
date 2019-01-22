MODULE flxrnf
   !!======================================================================
   !!                       ***  MODULE  flxrnf  ***
   !! Ocean forcing:  runoff
   !!=====================================================================
#if defined key_orca_r05
   !!----------------------------------------------------------------------
   !!   'key_orca_r05'                               ORCA R05 configuration
   !!----------------------------------------------------------------------
#  include "flxrnf_ORCA_R05.h90"
#else
   !!----------------------------------------------------------------------
   !!   Default option                                     Standard runoffs
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   flx_rnf      : monthly runoff read in a NetCDF file
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar
   USE ioipsl          ! NetCDF IPSL library

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC flx_rnf          ! routine call in step module

   !! * Shared module variables
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      runoff,           &  !: monthly runoff (kg/m2/s)
      upsadv,           &  !: mixed adv scheme in straits vicinity (hori.)
      upsrnfh              !: mixed adv scheme in runoffs vicinity (hori.)
   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &  !:
      upsrnfz              !: mixed adv scheme in runoffs vicinity (vert.)
   INTEGER, PUBLIC ::   &  !:
      nrunoff =  0 ,    &  !: runoff option (namelist)
      nrnf1, nrnf2         !: first and second record used

   !! * Module variable
   REAL(wp), DIMENSION(jpi,jpj,2) ::   &  !:
      rnfdta               !: monthly runoff data array (kg/m2/s)

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/SBC/flxrnf.F90,v 1.2 2005/11/16 16:14:39 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flx_rnf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE flx_rnf  ***
      !!       
      !! ** Purpose :   Introduce a climatological run off forcing
      !!
      !! ** Method :
      !!      Initialze each mouth of river with a monthly climatology 
      !!      provided from different data.
      !!     C a u t i o n : upward water flux, runoff is negative
      !!                     set at the last loop of the routine
      !!
      !! ** Action :
      !!
      !! References : 
      !!       J. D. Milliman and R. H. Meade, 1983 : world-wide delivery
      !!          of river sediment to the oceans, journal of geology vol 91
      !!          pp 1-21.
      !!       G. L. Russell and J. R. Miller, 1990 : global river runoff
      !!          calculated from a global atmospheric general circulation
      !!          model, journal of hydrology, 117(1990), pp 241-254.
      !!       F. Van Der Leeden, Troise F. L., Todd D. K. : the water
      !!          encyclopedia, second edition, lewis publishers.
      !!       J. W. Weatherly, J. E. Walsh : The effects of precipitation
      !!          and river runoff in a coupled ice-ocean model of Arctic
      !!          Climate dynamics 1996 12:785,798
      !!       Jacobs et al. 1992. J. Glaciol. 38 (130) 375-387.
      !!
      !! History :
      !!        !  94-10  (G.Madec, M. Pontaud, M. Imbard)  Original code
      !!        !  97-03  (G.Madec)  time dependent version
      !!        !  98-06  (J.M. Molines)  exact computation of zxy 
      !!                         for months that are not 30 days
      !!        !  98-07  (M. Imbard)  ORCA and mpp option
      !!        !  99-08  (J.P. Boulanger H.L.Ayina)  New rivers and 
      !!                         values given in m3/s 
      !!        !  00-04  (G. Madec, K. Roberts) add antarctica ice discharge.
      !!        !  00-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt       ! ocean time step

      !! * Local declarations
# if ! defined key_coupled
      INTEGER  ::   ji, jj                 ! dummy loop indices
      INTEGER ::   &
         i15 , imois , iman,            &  ! temporary integers
         idbd, idmeom                      !    "          "
      REAL(wp) ::   zxy
# endif
      CHARACTER (len=32) ::   &
         clname = 'runoff_1m_nomask'       ! monthly runoff filename
      INTEGER, PARAMETER :: jpmois = 12
      INTEGER  ::   ipi, ipj, ipk          ! temporary integers
      INTEGER  ::   ii0, ii1, ij0, ij1     !    "          "
      INTEGER, DIMENSION(jpmois) ::     &
         istep                             ! temporary workspace
      REAL(wp) ::   zdate0, zdt            ! temporary scalars
      REAL(wp), DIMENSION(jpk) ::       &
         zlev                              ! temporary workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zlon, zlat,                    &  ! temporary workspace
         zcoefr                            ! coeff of advection link to runoff
      !!----------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN

         SELECT CASE ( nrunoff )

         CASE ( 0 )
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'flx_rnf : No runoff in this simulation (nrunoff=0)'
            IF(lwp) WRITE(numout,*) '~~~~~~~'
            
         CASE ( 1 )
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff (nrunoff=1)'
            IF(lwp) WRITE(numout,*) '~~~~~~~'

         CASE ( 2 )
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'flx_rnf : monthly runoff with upsteam advection'
            IF(lwp) WRITE(numout,*) '~~~~~~~   in the vicinity of river mouths (nrunoff=2)'

         CASE DEFAULT
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' Error nrunoff = ', nrunoff, ' /= 0, 1 or 2'
            nstop = nstop + 1

         END SELECT

         ! Set runoffs and upstream coeff to zero
         runoff (:,:) = 0.e0
         upsrnfh(:,:) = 0.e0
         upsrnfz(:)   = 0.e0 
         upsadv (:,:) = 0.e0

      ENDIF


      ! 1. Initialization
      ! -----------------

      IF( nrunoff == 1 .OR. nrunoff == 2 ) THEN
# if ! defined key_coupled

         ! year, month, day

         iman  = jpmois
         i15   = nday / 16
         imois = nmonth + i15 - 1
         IF( imois == 0 )   imois = iman
         ! Number of days in the month
         IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
            idbd = nbiss(imois)
         ELSEIF( nleapy > 1 ) THEN
            idbd = nleapy
         ELSE
            idbd = nobis(imois)
         ENDIF
         ! Number of days between imois, 15 and the end of month
         idmeom = idbd - 15
# endif
         ipi = jpiglo
         ipj = jpjglo
         ipk = jpk
         zdt = rdt
         
         ! Open file

         IF( kt == nit000 ) THEN
            CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj,    &
               &           .false., ipi, ipj, ipk, zlon,        &
               &           zlat, zlev, jpmois, istep, zdate0,   &
               &           zdt, numrnf )
            !   Title, dimensions and tests
# if ! defined key_coupled
            IF( iman /= jpmois ) THEN
               IF(lwp) WRITE(numout,*)
               IF(lwp) WRITE(numout,*) 'problem with time coordinates'
               IF(lwp) WRITE(numout,*) ' iman ', iman, ' jpmois ', jpmois
               nstop = nstop + 1
            ENDIF
            IF(lwp) WRITE(numout,*) iman, istep, zdate0, rdt, numrnf
            IF(lwp) WRITE(numout,*) 'numrnf=', numrnf
            IF(lwp) WRITE(numout,*) 'jpmois=', jpmois
            IF(lwp) WRITE(numout,*) 'zdt=', zdt
# endif
            IF(ipi /= jpidta .AND. ipj /= jpjdta .AND. ipk /= 1) THEN
               IF(lwp)WRITE(numout,*) ' '
               IF(lwp)WRITE(numout,*) 'problem with dimensions'
               IF(lwp)WRITE(numout,*) ' ipi ', ipi, ' jpidta ', jpidta
               IF(lwp)WRITE(numout,*) ' ipj ', ipj, ' jpjdta ', jpjdta
               IF(lwp)WRITE(numout,*) ' ipk ', ipk, ' =? 1'
               nstop = nstop + 1
            ENDIF
            IF(lwp)WRITE(numout,*) 'ipi=', ipi, ' ipj=', ipj, ' ipk=', ipk
         ENDIF
         
# if ! defined key_coupled

         ! 2. Read monthly file of runoff
         ! ------------------------------

         IF( kt == nit000 .OR. imois /= nrnf1 ) THEN

            ! Calendar computation for interpolation
            !     nrnf1 number of the first array record used in the simulation
            !     nrnf2 number of the last  array record

            nrnf1 = imois
            nrnf2 = nrnf1 + 1
            nrnf1 = MOD( nrnf1, iman )
            IF( nrnf1 == 0 ) nrnf1 = iman
            nrnf2 = MOD( nrnf2, iman )
            IF( nrnf2 == 0 ) nrnf2 = iman
            
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) ' runoff monthly field'
               WRITE(numout,*) ' --------------------'
               WRITE(numout,*) ' NetCDF format'
               WRITE(numout,*)
               WRITE(numout,*) 'first array record used nrnf1 ',nrnf1
               WRITE(numout,*) 'last  array record used nrnf2 ',nrnf2
               WRITE(numout,*)
            ENDIF
            
            ! Read monthly runoff data in kg/m2/s
!ibug
            IF( kt == nit000 )   rnfdta(:,:,:) = 0.e0
!ibug
            CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
               &        , nrnf1, nrnf1, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,1) )
            CALL flinget( numrnf, 'sorunoff', jpidta, jpjdta, 1, jpmois   &
               &        , nrnf2, nrnf2, mig(1), nlci, mjg(1), nlcj, rnfdta(1:nlci,1:nlcj,2) )

            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) ' read runoff field ok'
            IF(lwp) WRITE(numout,*)

         ENDIF

         ! Linear interpolation and conversion in upward water flux
         ! C a u t i o n : runoff is negative and in kg/m2/s 

         zxy = FLOAT( nday + idmeom - idbd * i15 ) / idbd

         runoff(:,:) = -( ( 1.e0 - zxy ) * rnfdta(:,:,1) + zxy * rnfdta(:,:,2) )

         ! Runoff reduction
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( gphit(ji,jj) > 40 .AND. gphit(ji,jj) < 65 )   runoff(ji,jj) = 0.85 * runoff(ji,jj)
            END DO
         END DO
         
# endif

      ENDIF


      ! 3. Mixed advection scheme 
      ! -------------------------

      IF( nrunoff == 2 .AND. kt == nit000 ) THEN

         ! Upstream and centered scheme in the vicinity of river mouths

         !  Creates the array coef that contains the coefficient to affect to
         !  the upstream scheme. advection scheme will be:
         !  coefr * upstream + (1- coefr) centered
         !  coefr must be between 0 and 1.
!ibug
         zcoefr(:,:) = 0.e0
!ibug

         CALL flinget( numrnf, 'socoefr', jpidta, jpjdta, 1, jpmois, nrnf1,   &
            &          nrnf1, mig(1), nlci, mjg(1), nlcj, zcoefr(1:nlci,1:nlcj) )

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' read coefr for advection ok'
         IF(lwp) WRITE(numout,*)
         
         upsrnfh(:,:) = zcoefr(:,:)
         upsrnfz(:)   = 0.e0
         upsrnfz(1)   = 1.0
         upsrnfz(2)   = 1.0
         upsrnfz(3)   = 0.5
         upsrnfz(4)   = 0.25
         upsrnfz(5)   = 0.125
         
         IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
            ! ORCA_R2 configuration : upstream scheme in the Sound Strait
            ij0 = 116   ;   ij1 = 116
            ii0 = 144   ;   ii1 = 144   ;   upsrnfh( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
            ii0 = 145   ;   ii1 = 147   ;   upsrnfh( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
            ii0 = 148   ;   ii1 = 148   ;   upsrnfh( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
         ENDIF

      ENDIF

      ! Upstream and centered scheme in the vicinity of some straits

      IF( kt == nit000 ) THEN 

         IF( cp_cfg == "orca" ) THEN

            SELECT CASE ( jp_cfg )
            !                                        ! =======================
            CASE ( 4 )                               !  ORCA_R4 configuration 
               !                                     ! =======================

               !                                          ! Gibraltar Strait
               ii0 =  70   ;   ii1 =  71
               ij0 =  52   ;   ij1 =  53   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         
               !                                     ! =======================
            CASE ( 2 )                               !  ORCA_R2 configuration 
               !                                     ! =======================

               !                                          ! Gibraltar Strait
               ij0 = 102   ;   ij1 = 102
               ii0 = 138   ;   ii1 = 138   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.20
               ii0 = 139   ;   ii1 = 139   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40
               ii0 = 140   ;   ii1 = 140   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
               ij0 = 101   ;   ij1 = 102
               ii0 = 141   ;   ii1 = 141   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50

               !                                          ! Bab el Mandeb Strait
               ij0 =  87   ;   ij1 =  88
               ii0 = 164   ;   ii1 = 164   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.10
               ij0 =  88   ;   ij1 =  88
               ii0 = 163   ;   ii1 = 163   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
               ii0 = 162   ;   ii1 = 162   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40
               ii0 = 160   ;   ii1 = 161   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
               ij0 =  89   ;   ij1 =  89
               ii0 = 158   ;   ii1 = 160   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25
               ij0 =  90   ;   ij1 =  90
               ii0 = 160   ;   ii1 = 160   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.25

               !                                          ! Sound Strait
               ij0 = 116   ;   ij1 = 116
               ii0 = 145   ;   ii1 = 147   ;   upsadv( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.50
         
            END SELECT 

         ENDIF

      ENDIF
     
      ! 4. Closing all files
      ! --------------------

      IF( kt == nitend .AND. nrunoff >= 1 )   CALL flinclo( numrnf )

   END SUBROUTINE flx_rnf

#endif
   !!======================================================================
END MODULE flxrnf
