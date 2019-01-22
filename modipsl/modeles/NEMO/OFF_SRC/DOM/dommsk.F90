MODULE dommsk
   !!==============================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_msk        : compute land/ocean mask
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dom_msk        ! routine called by inidom.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/DOM/dommsk.F90,v 1.2 2005/11/16 16:12:12 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dom_msk
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) and baro-
      !!      tropic stream function  points (b).
      !!        Set mbathy to the number of non-zero w-levels of a water column
      !!      (if island in the domain (lk_isl=T), this is done latter in
      !!      routine solver_init)
      !!
      !! ** Method  :   The ocean/land mask is computed from the basin bathy-
      !!      metry in level (mbathy) which is defined or read in dommba.
      !!      mbathy equals 0 over continental T-point, -n over the nth 
      !!      island T-point, and the number of ocean level over the ocean.
      !!
      !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
      !!      t-point : 0. IF mbathy( ji ,jj) =< 0
      !!                1. IF mbathy( ji ,jj) >= jk
      !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
      !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
      !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
      !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
      !!                and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
      !!      b-point : the same definition as for f-point of the first ocean
      !!                level (surface level) but with 0 along coastlines.
      !!
      !!        The lateral friction is set through the value of fmask along
      !!      the coast and topography. This value is defined by shlat, a
      !!      namelist parameter:
      !!         shlat = 0, free slip  (no shear along the coast)
      !!         shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < shlat, strong slip        | in the lateral boundary layer
      !!
      !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
      !!      are defined with the proper value at lateral domain boundaries,
      !!      but bmask. indeed, bmask defined the domain over which the
      !!      barotropic stream function is computed. this domain cannot
      !!      contain identical columns because the matrix associated with
      !!      the barotropic stream function equation is then no more inverti-
      !!      ble. therefore bmask is set to 0 along lateral domain boundaries
      !!      even IF nperio is not zero.
      !!
      !!      In case of open boundaries (lk_obc=T):
      !!        - tmask is set to 1 on the points to be computed bay the open
      !!          boundaries routines.
      !!        - bmask is  set to 0 on the open boundaries.
      !!
      !!      Set mbathy to the number of non-zero w-levels of a water column
      !!                  mbathy = min( mbathy, 1 ) + 1
      !!      (note that the minimum value of mbathy is 2).
      !!
      !! ** Action :
      !!                     tmask    : land/ocean mask at t-point (=0. or 1.)
      !!                     umask    : land/ocean mask at u-point (=0. or 1.)
      !!                     vmask    : land/ocean mask at v-point (=0. or 1.)
      !!                     fmask    : land/ocean mask at f-point (=0. or 1.)
      !!                          =shlat along lateral boundaries
      !!                     bmask    : land/ocean mask at barotropic stream
      !!                          function point (=0. or 1.) and set to
      !!                          0 along lateral boundaries
      !!                   mbathy   : number of non-zero w-levels 
      !!
      !! History :
      !!        !  87-07  (G. Madec)  Original code
      !!        !  91-12  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
      !!        !  96-01  (G. Madec)  suppression of common work arrays
      !!        !  96-05  (G. Madec)  mask computed from tmask and sup-
      !!                 pression of the double computation of bmask
      !!        !  97-02  (G. Madec)  mesh information put in domhgr.F
      !!        !  97-07  (G. Madec)  modification of mbathy and fmask
      !!        !  98-05  (G. Roullet)  free surface
      !!        !  00-03  (G. Madec)  no slip accurate
      !!        !  01-09  (J.-M. Molines)  Open boundaries
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! *Local declarations
      INTEGER  ::   ji, jk     ! dummy loop indices
      INTEGER  ::   iif, iil, ijf, ijl
      INTEGER, DIMENSION(jpi,jpj) ::  imsk

      !!---------------------------------------------------------------------
      


      ! Interior domain mask (used for global sum)
      ! --------------------

      tmask_i(:,:) = tmask(:,:,1)
      iif = jpreci                         ! ???
      iil = nlci - jpreci + 1
      ijf = jprecj                         ! ???
      ijl = nlcj - jprecj + 1

      tmask_i( 1 :iif,   :   ) = 0.e0      ! first columns
      tmask_i(iil:jpi,   :   ) = 0.e0      ! last  columns (including mpp extra columns)
      tmask_i(   :   , 1 :ijf) = 0.e0      ! first rows
      tmask_i(   :   ,ijl:jpj) = 0.e0      ! last  rows (including mpp extra rows)


      ! north fold mask
      tpol(1:jpiglo) = 1.e0 
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0.e0
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0.e0
      ENDIF

      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot: only half of the nlcj-1 row
         if (mjg(ijl-1) == jpjglo-1) then
         DO ji = iif+1, iil-1
            tmask_i(ji,ijl-1) = tmask_i(ji,ijl-1) * tpol(mig(ji))
         END DO
         endif
      ENDIF 

      ! Control print
      ! -------------
      IF( nprint == 1 .AND. lwp ) THEN
         imsk(:,:) = INT( tmask_i(:,:) )
         WRITE(numout,*) ' tmask_i : '
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                           1, jpj, 1, 1, numout)
         WRITE (numout,*)
         WRITE (numout,*) ' dommsk: tmask for each level'
         WRITE (numout,*) ' ----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( tmask(:,:,jk) )

            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
      ENDIF

   END SUBROUTINE dom_msk

END MODULE dommsk
