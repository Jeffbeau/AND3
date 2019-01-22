!!DB 2008.09.30 -- Eliminated IOIPSL in hgr_read
!!   2009.08.31 -- Eliminated GYRE (jphgr_msh=5)
MODULE domhgr
   !!==============================================================================
   !!                       ***  MODULE domhgr   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_hgr        : initialize the horizontal mesh 
   !!   hgr_read       : read "coordinate" NetCDF file 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   !! * Module variables
   REAL(wp) ::   glam0, gphi0           ! variables corresponding to parameters
      !                                 ! ppglam0 ppgphi0 set in par_oce

   !! * Routine accessibility
   PUBLIC dom_hgr        ! called by domain.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domhgr.F90,v 1.15 2006/04/10 15:46:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dom_hgr
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_hgr  ***
      !!
      !! ** Purpose :   Compute the geographical position (in degre) of the 
      !!      model grid-points,  the horizontal scale factors (in meters) and 
      !!      the Coriolis factor (in s-1).
      !!
      !! ** Method  :   The geographical position of the model grid-points is
      !!      defined from analytical functions, fslam and fsphi, the deriva-
      !!      tives of which gives the horizontal scale factors e1,e2.
      !!      Defining two function fslam and fsphi and their derivatives in 
      !!      the two horizontal directions (fse1 and fse2), the model grid-
      !!      point position and scale factors are given by:
      !!         t-point:
      !!      glamt(i,j) = fslam(i    ,j    )   e1t(i,j) = fse1(i    ,j    )
      !!      gphit(i,j) = fsphi(i    ,j    )   e2t(i,j) = fse2(i    ,j    )
      !!         u-point:
      !!      glamu(i,j) = fslam(i+1/2,j    )   e1u(i,j) = fse1(i+1/2,j    )
      !!      gphiu(i,j) = fsphi(i+1/2,j    )   e2u(i,j) = fse2(i+1/2,j    )
      !!         v-point:
      !!      glamv(i,j) = fslam(i    ,j+1/2)   e1v(i,j) = fse1(i    ,j+1/2)
      !!      gphiv(i,j) = fsphi(i    ,j+1/2)   e2v(i,j) = fse2(i    ,j+1/2)
      !!            f-point:
      !!      glamf(i,j) = fslam(i+1/2,j+1/2)   e1f(i,j) = fse1(i+1/2,j+1/2)
      !!      gphif(i,j) = fsphi(i+1/2,j+1/2)   e2f(i,j) = fse2(i+1/2,j+1/2)
      !!      Where fse1 and fse2 are defined by:
      !!         fse1(i,j) = ra * rad * SQRT( (cos(phi) di(fslam))**2
      !!                                     +          di(fsphi) **2 )(i,j)
      !!         fse2(i,j) = ra * rad * SQRT( (cos(phi) dj(fslam))**2
      !!                                     +          dj(fsphi) **2 )(i,j)
      !!
      !!        The coriolis factor is given at z-point by:
      !!                     ff = 2.*omega*sin(gphif)      (in s-1)
      !!
      !!        This routine is given as an example, it must be modified
      !!      following the user s desiderata. nevertheless, the output as
      !!      well as the way to compute the model grid-point position and
      !!      horizontal scale factors must be respected in order to insure
      !!      second order accuracy schemes.
      !!
      !! N.B. If the domain is periodic, verify that scale factors are also
      !!      periodic, and the coriolis term again.
      !!
      !! ** Action  : - define  glamt, glamu, glamv, glamf: longitude of t-, 
      !!                u-, v- and f-points (in degre)
      !!              - define  gphit, gphiu, gphiv, gphit: latitude  of t-,
      !!               u-, v-  and f-points (in degre)
      !!        define e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f: horizontal
      !!      scale factors (in meters) at t-, u-, v-, and f-points.
      !!        define ff: coriolis factor at f-point
      !!
      !! References :
      !!      Marti, Madec and Delecluse, 1992, j. geophys. res., in press.
      !!
      !! History :
      !!        !  88-03  (G. Madec)
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  96-01  (G. Madec)  terrain following coordinates
      !!        !  97-02  (G. Madec)  print mesh informations
      !!        !  01-09  (M. Levy)  eel config: grid in km, beta-plane
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module, namelist
      !!   9.0  !  04-01  (A.M. Treguier, J.M. Molines) Case 4 (Mercator mesh)
      !!                  use of parameters in par_CONFIG-Rxx.h90, not in namelist
      !!        !  04-05  (A. Koch-Larrouy) Add Gyre configuration 
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER  ::   ji, jj              ! dummy loop indices
      INTEGER  ::   ii0, ii1, ij0, ij1  ! temporary integers
      INTEGER  ::   ijeq                ! index of equator T point (used in case 4)
      REAL(wp) ::   &
         zti, zui, zvi, zfi,         &  ! temporary scalars
         ztj, zuj, zvj, zfj,         &  !
         zphi0, zbeta, znorme,       &  !
         zarg, zf0, zminff, zmaxff
      REAL(wp) ::   &
         zlam1, zcos_alpha, zim1 , zjm1 , ze1, ze1deg,   &
         zphi1, zsin_alpha, zim05, zjm05
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_hgr : define the horizontal mesh from ithe following par_oce parameters '
         WRITE(numout,*) '~~~~~~~      type of horizontal mesh           jphgr_msh = ', jphgr_msh
         WRITE(numout,*) '             position of the first row and     ppglam0  = ', ppglam0
         WRITE(numout,*) '             column grid-point (degrees)       ppgphi0  = ', ppgphi0
         WRITE(numout,*) '             zonal      grid-spacing (degrees) ppe1_deg = ', ppe1_deg
         WRITE(numout,*) '             meridional grid-spacing (degrees) ppe2_deg = ', ppe2_deg
         WRITE(numout,*) '             zonal      grid-spacing (meters)  ppe1_m   = ', ppe1_m  
         WRITE(numout,*) '             meridional grid-spacing (meters)  ppe2_m   = ', ppe2_m  
      ENDIF


      SELECT CASE( jphgr_msh )   ! type of horizontal mesh

      CASE ( 0 )                     !  curvilinear coordinate on the sphere read in coordinate.nc file

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          curvilinear coordinate on the sphere read in "coordinate" file'

         CALL hgr_read           ! Defaultl option  :   NetCDF file

         !                                                ! =====================
!!DB: delete ORCA
!         IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN    ! ORCA R2 configuration

         !                                                ! ======================
!         IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05 configuration
         ! N.B. :  General case, lat and long function of both i and j indices:
         !     e1t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdila( zti, ztj ) )**2   &
         !                                  + (                           fsdiph( zti, ztj ) )**2  )
         !     e1u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdila( zui, zuj ) )**2   &
         !                                  + (                           fsdiph( zui, zuj ) )**2  )
         !     e1v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdila( zvi, zvj ) )**2   &
         !                                  + (                           fsdiph( zvi, zvj ) )**2  )
         !     e1f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdila( zfi, zfj ) )**2   &
         !                                  + (                           fsdiph( zfi, zfj ) )**2  )
         !
         !     e2t(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphit(ji,jj) ) * fsdjla( zti, ztj ) )**2   &
         !                                  + (                           fsdjph( zti, ztj ) )**2  )
         !     e2u(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiu(ji,jj) ) * fsdjla( zui, zuj ) )**2   &
         !                                  + (                           fsdjph( zui, zuj ) )**2  )
         !     e2v(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphiv(ji,jj) ) * fsdjla( zvi, zvj ) )**2   &
         !                                  + (                           fsdjph( zvi, zvj ) )**2  )
         !     e2f(ji,jj) = ra * rad * SQRT(  ( cos( rad*gphif(ji,jj) ) * fsdjla( zfi, zfj ) )**2   &
         !                                  + (                           fsdjph( zfi, zfj ) )**2  )


      CASE ( 1 )                     ! geographical mesh on the sphere with regular grid-spacing

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere with regular grid-spacing'
         IF(lwp) WRITE(numout,*) '          given by ppe1_deg and ppe2_deg' 

         DO jj = 1, jpj
            DO ji = 1, jpi
               zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - 1 + njmpp - 1 )
               zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - 1 + njmpp - 1 )
               zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
               zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - 1 + njmpp - 1 ) + 0.5
         ! Longitude
               glamt(ji,jj) = ppglam0 + ppe1_deg * zti
               glamu(ji,jj) = ppglam0 + ppe1_deg * zui
               glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
               glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
         ! Latitude
               gphit(ji,jj) = ppgphi0 + ppe2_deg * ztj
               gphiu(ji,jj) = ppgphi0 + ppe2_deg * zuj
               gphiv(ji,jj) = ppgphi0 + ppe2_deg * zvj
               gphif(ji,jj) = ppgphi0 + ppe2_deg * zfj
         ! e1
               e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
               e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
               e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
               e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
         ! e2
               e2t(ji,jj) = ra * rad * ppe2_deg
               e2u(ji,jj) = ra * rad * ppe2_deg
               e2v(ji,jj) = ra * rad * ppe2_deg
               e2f(ji,jj) = ra * rad * ppe2_deg
            END DO
         END DO


      CASE ( 2:3 )                   ! f- or beta-plane with regular grid-spacing

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          f- or beta-plane with regular grid-spacing'
         IF(lwp) WRITE(numout,*) '          given by ppe1_m and ppe2_m' 

         ! Position coordinates (in kilometers)
         !                          ==========
         glam0 = 0.e0
         gphi0 = - ppe2_m * 1.e-3
         
#if defined key_agrif && defined key_eel_r6
         IF (.Not.Agrif_Root()) THEN
           glam0  = Agrif_Parent(glam0) + (Agrif_ix())*Agrif_Parent(ppe1_m) * 1.e-3
           gphi0  = Agrif_Parent(gphi0) + (Agrif_iy())*Agrif_Parent(ppe2_m) * 1.e-3
           ppe1_m = Agrif_Parent(ppe1_m)/Agrif_Rhox()
           ppe2_m = Agrif_Parent(ppe2_m)/Agrif_Rhoy()          
         ENDIF
#endif         
         DO jj = 1, jpj
            DO ji = 1, jpi
               glamt(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 )       )
               glamu(ji,jj) = glam0 + ppe1_m * 1.e-3 * ( FLOAT( ji - 1 + nimpp - 1 ) + 0.5 )
               glamv(ji,jj) = glamt(ji,jj)
               glamf(ji,jj) = glamu(ji,jj)
   
               gphit(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 )       )
               gphiu(ji,jj) = gphit(ji,jj)
               gphiv(ji,jj) = gphi0 + ppe2_m * 1.e-3 * ( FLOAT( jj - 1 + njmpp - 1 ) + 0.5 )
               gphif(ji,jj) = gphiv(ji,jj)
            END DO
         END DO

         ! Horizontal scale factors (in meters)
         !                              ======
         e1t(:,:) = ppe1_m      ;      e2t(:,:) = ppe2_m
         e1u(:,:) = ppe1_m      ;      e2u(:,:) = ppe2_m
         e1v(:,:) = ppe1_m      ;      e2v(:,:) = ppe2_m
         e1f(:,:) = ppe1_m      ;      e2f(:,:) = ppe2_m

      CASE ( 4 )                     ! geographical mesh on the sphere, isotropic MERCATOR type

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          geographical mesh on the sphere, MERCATOR type'
         IF(lwp) WRITE(numout,*) '          longitudinal/latitudinal spacing given by ppe1_deg'
         IF ( ppgphi0 == -90 ) THEN
                IF(lwp) WRITE(numout,*) ' Mercator grid cannot start at south pole !!!! '
                IF(lwp) WRITE(numout,*) ' We stop '
                STOP
         ENDIF

         !  Find index corresponding to the equator, given the grid spacing e1_deg
         !  and the (approximate) southern latitude ppgphi0.
         !  This way we ensure that the equator is at a "T / U" point, when in the domain.
         !  The formula should work even if the equator is outside the domain.
         zarg = rpi / 4. - rpi / 180. * ppgphi0 / 2.
         ijeq = ABS( 180./rpi * LOG( COS( zarg ) / SIN( zarg ) ) / ppe1_deg )
         IF(  ppgphi0 > 0 )  ijeq = -ijeq

         IF(lwp) WRITE(numout,*) '          Index of the equator on the MERCATOR grid:', ijeq

         DO jj = 1, jpj
            DO ji = 1, jpi
               zti = FLOAT( ji - 1 + nimpp - 1 )         ;   ztj = FLOAT( jj - ijeq + njmpp - 1 )
               zui = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zuj = FLOAT( jj - ijeq + njmpp - 1 )
               zvi = FLOAT( ji - 1 + nimpp - 1 )         ;   zvj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
               zfi = FLOAT( ji - 1 + nimpp - 1 ) + 0.5   ;   zfj = FLOAT( jj - ijeq + njmpp - 1 ) + 0.5
         ! Longitude
               glamt(ji,jj) = ppglam0 + ppe1_deg * zti
               glamu(ji,jj) = ppglam0 + ppe1_deg * zui
               glamv(ji,jj) = ppglam0 + ppe1_deg * zvi
               glamf(ji,jj) = ppglam0 + ppe1_deg * zfi
         ! Latitude
               gphit(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* ztj ) )
               gphiu(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zuj ) )
               gphiv(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zvj ) )
               gphif(ji,jj) = 1./rad * ASIN ( TANH( ppe1_deg *rad* zfj ) )
         ! e1
               e1t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
               e1u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
               e1v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
               e1f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
         ! e2
               e2t(ji,jj) = ra * rad * COS( rad * gphit(ji,jj) ) * ppe1_deg
               e2u(ji,jj) = ra * rad * COS( rad * gphiu(ji,jj) ) * ppe1_deg
               e2v(ji,jj) = ra * rad * COS( rad * gphiv(ji,jj) ) * ppe1_deg
               e2f(ji,jj) = ra * rad * COS( rad * gphif(ji,jj) ) * ppe1_deg
            END DO
         END DO

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for jphgr_msh = ', jphgr_msh
         nstop = nstop + 1

      END SELECT


      ! Control printing : Grid informations (if not restart)
      ! ----------------

      IF(lwp .AND. .NOT.ln_rstart ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          longitude and e1 scale factors'
         WRITE(numout,*) '          ------------------------------'
         WRITE(numout,9300) ( ji, glamt(ji,1), glamu(ji,1),   &
            glamv(ji,1), glamf(ji,1),   &
            e1t(ji,1), e1u(ji,1),   &
            e1v(ji,1), e1f(ji,1), ji = 1, jpi,10)
9300     FORMAT( 1x, i4, f8.2,1x, f8.2,1x, f8.2,1x, f8.2, 1x,    &
            f19.10, 1x, f19.10, 1x, f19.10, 1x, f19.10 )
         
         WRITE(numout,*)
         WRITE(numout,*) '          latitude and e2 scale factors'
         WRITE(numout,*) '          -----------------------------'
         WRITE(numout,9300) ( jj, gphit(1,jj), gphiu(1,jj),   &
            &                     gphiv(1,jj), gphif(1,jj),   &
            &                     e2t  (1,jj), e2u  (1,jj),   &
            &                     e2v  (1,jj), e2f  (1,jj), jj = 1, jpj, 10 )
      ENDIF

      
      IF( nprint == 1 .AND. lwp ) THEN
         WRITE(numout,*) '          e1u e2u '
         CALL prihre( e1u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         CALL prihre( e2u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         WRITE(numout,*) '          e1v e2v  '
         CALL prihre( e1v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         CALL prihre( e2v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         WRITE(numout,*) '          e1f e2f  '
         CALL prihre( e1f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         CALL prihre( e2f,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
      ENDIF


      ! ================= !
      !  Coriolis factor  !
      ! ================= !

      SELECT CASE( jphgr_msh )   ! type of horizontal mesh

      CASE ( 0, 1, 4 )               ! mesh on the sphere

         ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 

      CASE ( 2 )                     ! f-plane at ppgphi0 

         ff(:,:) = 2. * omega * SIN( rad * ppgphi0 )

         IF(lwp) WRITE(numout,*) '          f-plane: Coriolis parameter = constant = ', ff(1,1)

      CASE ( 3 )                     ! beta-plane

         zbeta   = 2. * omega * COS( rad * ppgphi0 ) / ra                       ! beta at latitude ppgphi0
         zphi0   = ppgphi0 - FLOAT( jpjglo/2) * ppe2_m / ( ra * rad )           ! latitude of the first row F-points
         
#if defined key_agrif && defined key_eel_r6
         IF (.Not.Agrif_Root()) THEN
           zphi0 = ppgphi0 - FLOAT( Agrif_Parent(jpjglo)/2)*Agrif_Parent(ppe2_m) / (ra * rad)
         ENDIF
#endif         
         zf0     = 2. * omega * SIN( rad * zphi0 )                              ! compute f0 1st point south

         ff(:,:) = ( zf0  + zbeta * gphif(:,:) * 1.e+3 )                        ! f = f0 +beta* y ( y=0 at south)
         
         IF(lwp) WRITE(numout,*) 
         IF(lwp) WRITE(numout,*) ' Beta-plane: Beta parameter = constant = ', ff(nldi,nldj)
         IF(lwp) WRITE(numout,*) ' Coriolis parameter varies from ', ff(nldi,nldj),' to ', ff(nldi,nlej)
         IF( lk_mpp ) THEN 
            zminff=ff(nldi,nldj)
            zmaxff=ff(nldi,nlej)
            CALL mpp_min( zminff )   ! min over the global domain
            CALL mpp_max( zmaxff )   ! max over the global domain
            IF(lwp) WRITE(numout,*) ' Coriolis parameter varies globally from ', zminff,' to ', zmaxff
         END IF


      END SELECT


      ! Control of domain for symetrical condition
      ! ------------------------------------------
      ! The equator line must be the latitude coordinate axe

      IF( nperio == 2 ) THEN
         znorme = SQRT( SUM( gphiu(:,2) * gphiu(:,2) ) ) / FLOAT( jpi )
         IF( znorme > 1.e-13 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' ===>>>> : symmetrical condition: rerun with good equator line'
            nstop = nstop + 1
         ENDIF
      ENDIF

   END SUBROUTINE dom_hgr


   SUBROUTINE hgr_read
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE hgr_read  ***
      !!
      !! ** Purpose :   Read a coordinate file in NetCDF format 
      !!
      !! ** Method  :   The mesh file has been defined trough a analytical 
      !!      or semi-analytical method. It is read in a NetCDF file. 
      !!     
      !! References :
      !!      Marti, Madec and Delecluse, 1992, JGR, 97, 12,763-12,766.
      !!      Madec, Imbard, 1996, Clim. Dyn., 12, 381-388.
      !!
      !! History :
      !!        !         (O. Marti)  Original code
      !!        !  91-03  (G. Madec)
      !!        !  92-07  (M. Imbard)
      !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
      !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
      !!   8.5  !  02-06  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf


      !! * Local declarations
      LOGICAL ::   llog = .FALSE.
      CHARACTER(len=21) ::   clname
      INTEGER  ::   ji, jj              ! dummy loop indices
      INTEGER  ::   inum                ! temporary logical unit
      INTEGER  ::   ilev, itime         ! temporary integers
      REAL(wp) ::   zdt, zdate0         ! temporary scalars
      REAL(wp) ::   zdept(1)            ! temporary workspace
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
         zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
!!DB
      INTEGER :: status
      REAL(wp), DIMENSION(jpidta,jpjdta,1) ::   data_3d

      !!----------------------------------------------------------------------

!!DB
      clname = 'coordinates.nc'


#if defined key_agrif
      if ( .NOT. Agrif_Root() ) then
         clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
      endif
#endif         


      ! 1. Read of the grid coordinates and scale factors
      ! -------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'hgr_read : read the horizontal coordinates'
         WRITE(numout,*) '~~~~~~~~~~~      jpiglo = ', jpiglo, ' jpjglo = ', jpjglo, ' jpk = ', jpk
      ENDIF

      ! read the file
      itime = 0
      ilev = 1   
      zlamt(:,:) = 0.e0
      zphit(:,:) = 0.e0

!!DB: replace IOIPSL with lib_ncdf
      call ncdf_read_global(clname,'glamt',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            glamt(ji,jj) =  zdta(mig(ji),mjg(jj))
         END DO
      END DO
      call ncdf_read_global(clname,'glamu',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            glamu(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'glamv',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            glamv(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'glamf',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            glamf(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'gphit',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            gphit(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
!!DB: Initialize coriolis arrays along OBC
!!Uses f at T-pt
!!REM:   ff(:,:) = 2. * omega * SIN( rad * gphif(:,:) ) 
      do ji = 1, jpidta
         ff_s(ji) = 2. * omega * sin(rad*zdta(ji,1))
         ff_n(ji) = 2. * omega * sin(rad*zdta(ji,jpjdta))  !!for future use
      enddo
      do jj = 1, jpjdta
         ff_e(jj) = 2. * omega * sin(rad*zdta(jpidta,jj))
         ff_w(jj) =  2. * omega * sin(rad*zdta(1,jj))
      enddo

      call ncdf_read_global(clname,'gphiu',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            gphiu(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'gphiv',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            gphiv(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'gphif',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            gphif(ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO

      call ncdf_read_global(clname,'e1t',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e1t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'e1u',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e1u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'e1v',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e1v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'e1f',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e1f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'e2t',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e2t  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
      call ncdf_read_global(clname,'e2u',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e2u  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
!!DB: 2007.12.11 -- Initialize boundary array
      do ji = 1, jpidta
         e2u_s(ji) = zdta(ji,1)
         e2u_n(ji) = zdta(ji,jpjdta)
      enddo

      call ncdf_read_global(clname,'e2v',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e2v  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO
!!DB: 2007.12.11 -- Initialize boundary arrays
!!   Not sure at this time exactly what i I should use, but 2 should be OK
      do jj = 1, jpjdta
         e2v_e(jj) = zdta(jpidta,jj)
         e2v_w(jj) = zdta(2,jj)
      enddo

      call ncdf_read_global(clname,'e2f',data_3d,-1,status)
      zdta(:,:) = data_3d(:,:,1)
      DO jj = 1, nlcj
         DO ji = 1, nlci
            e2f  (ji,jj) = zdta(mig(ji),mjg(jj))                    
         END DO
      END DO


      ! set extra rows add in mpp to none zero values
      DO jj = nlcj+1, jpj
         DO ji = 1, nlci
            glamt(ji,jj) = glamt(ji,1)   ;   gphit(ji,jj) = gphit(ji,1)
            glamu(ji,jj) = glamu(ji,1)   ;   gphiu(ji,jj) = gphiu(ji,1)
            glamv(ji,jj) = glamv(ji,1)   ;   gphiv(ji,jj) = gphiv(ji,1)
            glamf(ji,jj) = glamf(ji,1)   ;   gphif(ji,jj) = gphif(ji,1)
            e1t  (ji,jj) = e1t  (ji,1)   ;   e2t  (ji,jj) = e2t  (ji,1)
            e1u  (ji,jj) = e1u  (ji,1)   ;   e2u  (ji,jj) = e2u  (ji,1)
            e1v  (ji,jj) = e1v  (ji,1)   ;   e2v  (ji,jj) = e2v  (ji,1)
            e1f  (ji,jj) = e1f  (ji,1)   ;   e2f  (ji,jj) = e2f  (ji,1)
         END DO
      END DO

      ! set extra columns add in mpp to none zero values
      DO ji = nlci+1, jpi
         glamt(ji,:) = glamt(1,:)   ;   gphit(ji,:) = gphit(1,:)
         glamu(ji,:) = glamu(1,:)   ;   gphiu(ji,:) = gphiu(1,:)
         glamv(ji,:) = glamv(1,:)   ;   gphiv(ji,:) = gphiv(1,:)
         glamf(ji,:) = glamf(1,:)   ;   gphif(ji,:) = gphif(1,:)
         e1t  (ji,:) = e1t  (1,:)   ;   e2t  (ji,:) = e2t  (1,:)
         e1u  (ji,:) = e1u  (1,:)   ;   e2u  (ji,:) = e2u  (1,:)
         e1v  (ji,:) = e1v  (1,:)   ;   e2v  (ji,:) = e2v  (1,:)
         e1f  (ji,:) = e1f  (1,:)   ;   e2f  (ji,:) = e2f  (1,:)
      END DO

   END SUBROUTINE hgr_read

   !!======================================================================
END MODULE domhgr
