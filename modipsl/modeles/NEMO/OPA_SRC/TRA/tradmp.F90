MODULE tradmp
   !!======================================================================
   !!                       ***  MODULE  tradmp  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!======================================================================
#if   defined key_tradmp   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   key_tradmp                                         internal damping
   !!----------------------------------------------------------------------
   !!   tra_dmp      : update the tracer trend with the internal damping
   !!   tra_dmp_init : initialization, namlist read, parameters control
   !!   dtacof_zoom  : restoring coefficient for zoom domain
   !!   dtacof       : restoring coefficient for global domain
   !!   cofdis       : compute the distance to the coastline
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! Define parameters for the routines
   USE dtatem          ! temperature data
   USE dtasal          ! salinity data
   USE zdfmxl          ! mixed layer depth
   USE lib_mpp         ! distribued memory computing
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_dmp   ! routine called by step.F90

   !! * Shared module variables
   LOGICAL , PUBLIC &
#if ! defined key_agrif
   , PARAMETER  &
#endif
   ::   lk_tradmp = .TRUE.    !: internal damping flag

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &
      strdmp,              &  ! damping salinity trend (psu/s)
      resto                   ! restoring coeff. on T and S (s-1)

   !! * Module variables
   INTEGER  ::             & !!! * newtonian damping namelist (mandmp) *
      ndmp   =   -1 ,      &  ! = 0/-1/'latitude' for damping over T and S
      ndmpf  =    2 ,      &  ! = 1 create a damping.coeff NetCDF file 
      nmldmp =    0           ! = 0/1/2 flag for damping in the mixed layer
   REAL(wp) ::             & !!!  * newtonian damping namelist *
      sdmp   =   50.,      &  ! surface time scale for internal damping (days)
      bdmp   =  360.,      &  ! bottom time scale for internal damping (days)
      hdmp   =  800.          ! depth of transition between sdmp and bdmp (meters)

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tradmp.F90,v 1.15 2006/04/19 14:43:17 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp  ***
      !!                  
      !! ** Purpose :   Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed 
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - update the tracer trends (ta,sa) with the newtonian 
      !!                damping trends.
      !!              - save the trends in (ttrd,strd) ('key_trdtra')
      !!
      !! History :
      !!   7.0  !         (G. Madec)  Original code
      !!        !  96-01  (G. Madec) 
      !!        !  97-05  (G. Madec)  macro-tasked on jk-slab
      !!   8.5  !  02-08  (G. Madec)  free form + modules
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
                         ztdsa => va      ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   ztest, zta, zsa       ! temporary scalars
      !!----------------------------------------------------------------------

      ! 0. Initialization (first time-step only)
      !    --------------
      IF( kt == nit000 ) CALL tra_dmp_init

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = sa(:,:,:) 
      ENDIF

      ! 1. Newtonian damping trends on tracer fields
      ! --------------------------------------------
      !    compute the newtonian damping trends depending on nmldmp

      SELECT CASE ( nmldmp )

      CASE( 0 )                ! newtonian damping throughout the water column
!!DBG: 2009.06.24 -- HARDWIRE strong sfce S restoring
!         DO jk = 1, 1
!            DO jj = 2, jpjm1
!               DO ji = fs_2, fs_jpim1   ! vector opt.
!                  zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
!!                  zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
!                  zsa = (1./(3.*rday)) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
!                  ! add the trends to the general tracer trends
!                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
!                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
!                  ! save the salinity trend (used in flx to close the salt budget)
!                  strdmp(ji,jj,jk) = zsa
!               END DO
!            END DO
!         END DO
!         DO jk = 2, jpkm1

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
                  zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
                  ! add the trends to the general tracer trends
                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp(ji,jj,jk) = zsa
               END DO
            END DO
         END DO

      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  ztest = avt(ji,jj,jk) - 5.e-4
                  IF( ztest < 0. ) THEN
                     zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
                     zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
                  ELSE
                     zta = 0.e0
                     zsa = 0.e0
                  ENDIF
                  ! add the trends to the general tracer trends
                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp(ji,jj,jk) = zsa
               END DO
            END DO
         END DO

      CASE ( 2 )                ! no damping in the mixed layer 
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     zta = resto(ji,jj,jk) * ( t_dta(ji,jj,jk) - tb(ji,jj,jk) )
                     zsa = resto(ji,jj,jk) * ( s_dta(ji,jj,jk) - sb(ji,jj,jk) )
                  ELSE
                     zta = 0.e0
                     zsa = 0.e0
                  ENDIF
                  ! add the trends to the general tracer trends
                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta
                  sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp(ji,jj,jk) = zsa
               END DO
            END DO
         END DO

      END SELECT

      ! save the trends for diagnostic
      ! damping salinity trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
         ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
         CALL trd_mod(ztdta, ztdsa, jpttddoe, 'TRA', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' dmp  - Ta: ', mask1=tmask, &
            &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
      ENDIF


   END SUBROUTINE tra_dmp


   SUBROUTINE tra_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the nammbf namelist and check the parameters
      !!      called by tra_dmp at the first timestep (nit000)
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      NAMELIST/namtdp/ ndmp, ndmpf, nmldmp, sdmp, bdmp, hdmp
      !!----------------------------------------------------------------------

      ! Read Namelist namtdp : temperature and salinity damping term
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namtdp )
      IF( lzoom )   nmldmp = 0           ! restoring to climatology at closed north or south boundaries

      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dmp : T and S newtonian damping'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist namtdp : set damping parameter'
         WRITE(numout,*)
         WRITE(numout,*) '             T and S damping option         ndmp   = ', ndmp
         WRITE(numout,*) '             create a damping.coeff file    ndmpf  = ', ndmpf
         WRITE(numout,*) '             mixed layer damping option     nmldmp = ', nmldmp, '(zoom: forced to 0)'
         WRITE(numout,*) '             surface time scale (days)      sdmp   = ', sdmp
         WRITE(numout,*) '             bottom time scale (days)       bdmp   = ', bdmp
         WRITE(numout,*) '             depth of transition (meters)   hdmp   = ', hdmp
         WRITE(numout,*)
      ENDIF

      SELECT CASE ( ndmp )

!byoung------------------------------------------------------------------------
      CASE ( 0 )               ! SS: damping in the Scotian Shelf only
         IF(lwp) WRITE(numout,*) '          tracer damping in the Scotian Shelf sea only'
!----------------------------------------------------------------------------	 
      CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
         IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'

      CASE ( 1:90 )             ! Damping poleward of 'ndmp' degrees
         IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmp, ' degrees'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for ndmp = ', ndmp
         nstop = nstop + 1

      END SELECT


      SELECT CASE ( nmldmp )

      CASE ( 0 )                ! newtonian damping throughout the water column
         IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'

      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'

      CASE ( 2 )                ! no damping in the mixed layer 
         IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for nmldmp = ', nmldmp
         nstop = nstop + 1

      END SELECT

      IF( .NOT.lk_dtasal .OR. .NOT.lk_dtatem ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          no temperature and/or salinity data '
         IF(lwp) WRITE(numout,*) '          define key_dtatem and key_dtasal'
         nstop = nstop + 1
      ENDIF


      strdmp(:,:,:) = 0.e0       ! internal damping salinity trend (used in ocesbc)

      ! Damping coefficients initialization
      ! -----------------------------------

      IF( lzoom ) THEN
         CALL dtacof_zoom
      ELSE
         CALL dtacof
      ENDIF

   END SUBROUTINE tra_dmp_init


   SUBROUTINE dtacof_zoom
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dtacof_zoom  ***
      !!
      !! ** Purpose :   Compute the damping coefficient for zoom domain
      !!
      !! ** Method  : - set along closed boundary due to zoom a damping over
      !!      6 points with a max time scale of 5 days.
      !!              - ORCA arctic/antarctic zoom: set the damping along
      !!      south/north boundary over a latitude strip.
      !!
      !! ** Action  : - resto, the damping coeff. for T and S
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn      ! dummy loop indices
      REAL(wp) ::   &
         zlat, zlat0, zlat1, zlat2     ! temporary scalar
      REAL(wp), DIMENSION(6)  ::   &
         zfact                         ! temporary workspace
      !!----------------------------------------------------------------------

      zfact(1) =  1.
      zfact(2) =  1. 
      zfact(3) = 11./12.
      zfact(4) =  8./12.
      zfact(5) =  4./12.
      zfact(6) =  1./12.
      zfact(:) = zfact(:) / ( 5. * rday )    ! 5 days max restoring time scale

      resto(:,:,:) = 0.e0

      ! damping along the forced closed boundary over 6 grid-points
      DO jn = 1, 6
         IF( lzoom_w )   resto( mi0(jn+jpizoom):mi1(jn+jpizoom), : , : ) = zfact(jn) ! west  closed
         IF( lzoom_s )   resto( : , mj0(jn+jpjzoom):mj1(jn+jpjzoom), : ) = zfact(jn) ! south closed 
         IF( lzoom_e )   resto( mi0(jpiglo+jpizoom-1-jn):mi1(jpiglo+jpizoom-1-jn) , : , : ) &
                       &              = zfact(jn)                                 ! east  closed 
         IF( lzoom_n )   resto( : , mj0(jpjglo+jpjzoom-1-jn):mj1(jpjglo+jpjzoom-1-jn) , : ) &
                       &              = zfact(jn)                                 ! north closed
      END DO


      IF( lzoom_arct .AND. lzoom_anta ) THEN

         ! ====================================================
         !  ORCA configuration : arctic zoom or antarctic zoom
         ! ====================================================

         IF(lwp) WRITE(numout,*)
         IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA    Arctic zoom'
         IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              dtacof_zoom : ORCA Antarctic zoom'
         IF(lwp) WRITE(numout,*)

         ! ... Initialization : 
         !     zlat0 : latitude strip where resto decreases
         !     zlat1 : resto = 1 before zlat1
         !     zlat2 : resto decreases from 1 to 0 between zlat1 and zlat2
         resto(:,:,:) = 0.e0
         zlat0 = 10.
         zlat1 = 30.
         zlat2 = zlat1 + zlat0

         ! ... Compute arrays resto ; value for internal damping : 5 days
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zlat = ABS( gphit(ji,jj) )
                  IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
                     resto(ji,jj,jk) = 0.5 * ( 1./(5.*rday) ) *   &
                        ( 1. - cos(rpi*(zlat2-zlat)/zlat0) ) 
                  ELSE IF ( zlat < zlat1 ) THEN
                     resto(ji,jj,jk) = 1./(5.*rday)
                  ENDIF
               END DO
            END DO
         END DO

      ENDIF

      ! ... Mask resto array
      resto(:,:,:) = resto(:,:,:) * tmask(:,:,:)

   END SUBROUTINE dtacof_zoom

   SUBROUTINE dtacof
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dtacof  ***
      !!
      !! ** Purpose :   Compute the damping coefficient
      !!
      !! ** Method  :   Arrays defining the damping are computed for each grid
      !!      point for temperature and salinity (resto)
      !!      Damping depends on distance to coast, depth and latitude
      !!
      !! ** Action  : - resto, the damping coeff. for T and S
      !!
      !! History :
      !!   5.0  !  91-03  (O. Marti, G. Madec)  Original code
      !!        !  92-06  (M. Imbard)  doctor norme
      !!        !  96-01  (G. Madec) statement function for e3
      !!        !  98-07  (M. Imbard, G. Madec) ORCA version
      !!        !  00-08  (G. Madec, D. Ludicone) 
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Local declarations
      INTEGER ::   ji, jj, jk, je      ! dummy loop indices
      INTEGER, PARAMETER ::   jpmois=1
      INTEGER ::   ipi, ipj, ipk       ! temporary integers
      INTEGER ::   ii0, ii1, ij0, ij1  !    "          "
      INTEGER ::   &
         idmp,     &  ! logical unit for file restoring damping term
         icot         ! logical unit for file distance to the coast
      INTEGER :: itime, istep(jpmois), ie
      LOGICAL :: llbon
      CHARACTER (len=32) ::  clname, clname2, clname3
      REAL(wp) ::   &
         zdate0, zinfl, zlon,         & ! temporary scalars
         zlat, zlat0, zlat1, zlat2,   & !    "         "
         zsdmp, zbdmp                   !    "         "
      REAL(wp), DIMENSION(jpk) ::   &
         zdept, zhfac
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zmrs, zlamt, zphit
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zdct
!!DB
      INTEGER :: k_hdmp

      !!----------------------------------------------------------------------

      ! ====================================
      !  ORCA configuration : global domain
      ! ====================================

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '              dtacof : Global domain of ORCA'
      IF(lwp) WRITE(numout,*) '              ------------------------------'

      ! ... Initialization : 
      !   zdct()      : distant to the coastline
      !   resto()     : array of restoring coeff. on T and S

      zdct (:,:,:) = 0.e0
      resto(:,:,:) = 0.e0

      IF ( ndmp > 0 ) THEN

         !    ------------------------------------
         !     Damping poleward of 'ndmp' degrees
         !    ------------------------------------

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              Damping poleward of ', ndmp,' deg.'
         IF(lwp) WRITE(numout,*)

         ! ... Distance to coast (zdct)

         !   ... Test the existance of distance-to-coast file
         itime = jpmois
         ipi = jpiglo
         ipj = jpjglo
         ipk = jpk
         clname = 'dist.coast'
         DO je = 1,32
            IF( clname(je:je) == ' ' ) go to 140
         END DO
140      CONTINUE
         ie = je
         clname2 = clname(1:ie-1)//".nc"
         inquire( FILE = clname2, EXIST = llbon )

         IF ( llbon ) THEN

            !   ... Read file distance to coast if possible
            CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj, .false.,   &
               ipi, ipj, ipk, zlamt, zphit, zdept, jpmois,   &
               istep, zdate0, rdt, icot )
            CALL flinget( icot, 'Tcoast', jpidta, jpjdta, jpk,    &
               jpmois, 1, 1, mig(1), nlci, mjg(1), nlcj, zdct(1:nlci,1:nlcj,1:jpk) )
            CALL flinclo( icot )
            IF(lwp)WRITE(numout,*) '    ** : File dist.coast.nc read'

         ELSE

            !   ... Compute and save the distance-to-coast array (output in zdct)
            CALL cofdis ( zdct )

         ENDIF

         ! ... Compute arrays resto 
         !      zinfl : distance of influence for damping term
         !      zlat0 : latitude strip where resto decreases
         !      zlat1 : resto = 0 between -zlat1 and zlat1
         !      zlat2 : resto increases from 0 to 1 between |zlat1| and |zlat2|
         !          and resto = 1 between |zlat2| and 90 deg.
         zinfl = 1000.e3
         zlat0 = 10
         zlat1 = ndmp
         zlat2 = zlat1 + zlat0

         DO jj = 1, jpj
            DO ji = 1, jpi
               zlat = ABS( gphit(ji,jj) )
               IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
                  resto(ji,jj,1) = 0.5 * ( 1. - cos(rpi*(zlat-zlat1)/zlat0 ) )
               ELSEIF ( zlat > zlat2 ) THEN
                  resto(ji,jj,1) = 1.
               ENDIF
            END DO
         END DO

         !   ... North Indian ocean (20N/30N x 45E/100E) : resto=0
         IF ( ndmp == 20 ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zlat = gphit(ji,jj)
                  zlon = MOD( glamt(ji,jj), 360. )
                  IF ( zlat1 < zlat .AND. zlat < zlat2 .AND.   &
                     45.  < zlon .AND. zlon < 100. ) THEN
                     resto(ji,jj,1) = 0.
                  ENDIF
               END DO
            END DO
         ENDIF

         zsdmp = 1./(sdmp * rday)
         zbdmp = 1./(bdmp * rday)
         DO jk = 2, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zdct(ji,jj,jk) = MIN( zinfl, zdct(ji,jj,jk) )

                  !   ... Decrease the value in the vicinity of the coast
                  resto(ji,jj,jk) = resto(ji,jj,1)*0.5   &
                     &            * ( 1. - COS( rpi*zdct(ji,jj,jk)/zinfl) )

                  !   ... Vertical variation from zsdmp (sea surface) to zbdmp (bottom)
                  resto(ji,jj,jk) = resto(ji,jj,jk)   &
                     &            * ( zbdmp + (zsdmp-zbdmp)*EXP(-fsdept(ji,jj,jk)/hdmp) )
               END DO
            END DO
         END DO

      ENDIF


!!DB: deleted ORCA
!      IF( cp_cfg == "orca" .AND. ( ndmp > 0 .OR. ndmp == -1 ) ) THEN
!!DB -----------------------------------------------------------------------
      if( cp_cfg == "SS" .AND. ndmp == 0 ) THEN
      
         resto(:,:,:) = 0.e0

!!DB 2008.06.30 
!Determine depth at which hdmp applies
!ifdef defaults to old system if key_flux_bulk_* not defined
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
         k_hdmp = jpk
         do jk = 1, jpk
            if(hdmp > gdept(jk)) k_hdmp = jk
         enddo
!!DBG
         if(lwp) write(numout2,*)'DBG: hdmp, k_hdmp: ', hdmp, k_hdmp 
#else
         k_hdmp = 0
#endif
         
         DO jk = 1, k_hdmp
            DO jj = 1, jpj
               DO ji = 1, jpi
                  resto(ji,jj,jk) = 1./(sdmp * rday)   
               END DO
            END DO
         END DO
         DO jk = k_hdmp+1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  resto(ji,jj,jk) = 1./(bdmp * rday)   
               END DO
            END DO
         END DO
         
         resto(:,:, : ) = resto(:,:,:) * tmask(:,:,:)
!!DB
         resto(:,:, 1 ) = 1./(sdmp *rday) 
         
         !-------------------------------------------------------------------------------         	 
      ELSE
         !    ------------
         !     No damping
         !    ------------
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) 'Choose a correct value of ndmp or DO NOT defined key_tradmp'
         nstop = nstop + 1
      ENDIF

      !    ----------------------------
      !     Create Print damping array
      !    ----------------------------

      ! ndmpf   : = 1 create a damping.coeff NetCDF file

      IF( ndmpf == 1 ) THEN
         IF(lwp) WRITE(numout,*) '              create damping.coeff.nc file'
         itime   = 0
         clname3 = 'damping.coeff'
         CALL ymds2ju( 0     , 1     , 1      , 0.e0 , zdate0 )
         CALL restini( 'NONE', jpi   , jpj    , glamt, gphit,    &
                       jpk   , gdept , clname3, itime, zdate0,   &
                       rdt   , idmp, domain_id=nidom )
         CALL restput( idmp, 'Resto', jpi, jpj, jpk,   &
                       0   , resto  )
         CALL restclo( idmp )
      ENDIF

   END SUBROUTINE dtacof


   SUBROUTINE cofdis ( pdct )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cofdis  ***
      !!
      !! ** Purpose :   Compute the distance between ocean T-points and the
      !!      ocean model coastlines. Save the distance in a NetCDF file.
      !!
      !! ** Method  :   For each model level, the distance-to-coast is 
      !!      computed as follows : 
      !!       - The coastline is defined as the serie of U-,V-,F-points
      !!      that are at the ocean-land bound.
      !!       - For each ocean T-point, the distance-to-coast is then 
      !!      computed as the smallest distance (on the sphere) between the 
      !!      T-point and all the coastline points.
      !!       - For land T-points, the distance-to-coast is set to zero.
      !!      C A U T I O N : Computation not yet implemented in mpp case.
      !!
      !! ** Action  : - pdct, distance to the coastline (argument)
      !!              - NetCDF file 'dist.coast.nc' 
      !!        
      !! History :
      !!   7.0  !  01-02  (M. Imbard)  Original code
      !!   8.1  !  01-02  (G. Madec, E. Durand)
      !!   8.5  !  02-08  (G. Madec, E. Durand)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
         pdct                     ! distance to the coastline

      !! * local declarations
      INTEGER :: ji, jj, jk, jl      ! dummy loop indices
      INTEGER :: iju, ijt            ! temporary integers
      INTEGER :: icoast, itime
      INTEGER ::   &
         icot         ! logical unit for file distance to the coast
      LOGICAL, DIMENSION(jpi,jpj) ::   &
         llcotu, llcotv, llcotf   ! ???
      CHARACTER (len=32) ::   clname
      REAL(wp) ::   zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zxt, zyt, zzt,                 &  ! cartesian coordinates for T-points
         zmask                             
      REAL(wp), DIMENSION(3*jpi*jpj) ::   &
         zxc, zyc, zzc, zdis      ! temporary workspace
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cofdis : compute the distance to coastline'
      IF(lwp) WRITE(numout,*) '~~~~~~'
      IF(lwp) WRITE(numout,*)
      IF( lk_mpp ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '         Computation not yet implemented with key_mpp_...'
         IF(lwp) WRITE(numout,*) '         Rerun the code on another computer or '
         IF(lwp) WRITE(numout,*) '         create the "dist.coast.nc" file using IDL'
         nstop = nstop + 1
      ENDIF

      pdct(:,:,:) = 0.e0
      zxt(:,:) = cos( rad * gphit(:,:) ) * cos( rad * glamt(:,:) )
      zyt(:,:) = cos( rad * gphit(:,:) ) * sin( rad * glamt(:,:) )
      zzt(:,:) = sin( rad * gphit(:,:) )


      ! 1. Loop on vertical levels
      ! --------------------------
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Define the coastline points (U, V and F)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zmask(ji,jj) =  ( tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk) &
                   &           + tmask(ji,jj  ,jk) + tmask(ji+1,jj  ,jk) )
               llcotu(ji,jj) = ( tmask(ji,jj,  jk) + tmask(ji+1,jj  ,jk) == 1. ) 
               llcotv(ji,jj) = ( tmask(ji,jj  ,jk) + tmask(ji  ,jj+1,jk) == 1. ) 
               llcotf(ji,jj) = ( zmask(ji,jj) > 0. ) .AND. ( zmask(ji,jj) < 4. )
            END DO
         END DO

         ! Lateral boundaries conditions
         llcotu(:, 1 ) = umask(:,  2  ,jk) == 1
         llcotu(:,jpj) = umask(:,jpjm1,jk) == 1
         llcotv(:, 1 ) = vmask(:,  2  ,jk) == 1
         llcotv(:,jpj) = vmask(:,jpjm1,jk) == 1
         llcotf(:, 1 ) = fmask(:,  2  ,jk) == 1
         llcotf(:,jpj) = fmask(:,jpjm1,jk) == 1

         IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
            llcotu( 1 ,:) = llcotu(jpim1,:)
            llcotu(jpi,:) = llcotu(  2  ,:)
            llcotv( 1 ,:) = llcotv(jpim1,:)
            llcotv(jpi,:) = llcotv(  2  ,:)
            llcotf( 1 ,:) = llcotf(jpim1,:)
            llcotf(jpi,:) = llcotf(  2  ,:)
         ELSE
            llcotu( 1 ,:) = umask(  2  ,:,jk) == 1
            llcotu(jpi,:) = umask(jpim1,:,jk) == 1
            llcotv( 1 ,:) = vmask(  2  ,:,jk) == 1
            llcotv(jpi,:) = vmask(jpim1,:,jk) == 1
            llcotf( 1 ,:) = fmask(  2  ,:,jk) == 1
            llcotf(jpi,:) = fmask(jpim1,:,jk) == 1
         ENDIF
         IF( nperio == 3 .OR. nperio == 4 ) THEN
            DO ji = 1, jpim1
               iju = jpi - ji + 1
               llcotu(ji,jpj  ) = llcotu(iju,jpj-2)
               llcotf(ji,jpj-1) = llcotf(iju,jpj-2)
               llcotf(ji,jpj  ) = llcotf(iju,jpj-3)
            END DO
            DO ji = jpi/2, jpi-1
               iju = jpi - ji + 1
               llcotu(ji,jpjm1) = llcotu(iju,jpjm1)
            END DO
            DO ji = 2, jpi
               ijt = jpi - ji + 2
               llcotv(ji,jpj-1) = llcotv(ijt,jpj-2)
               llcotv(ji,jpj  ) = llcotv(ijt,jpj-3)
            END DO
         ENDIF
         IF( nperio == 5 .OR. nperio == 6 ) THEN
            DO ji = 1, jpim1
               iju = jpi - ji
               llcotu(ji,jpj  ) = llcotu(iju,jpj-1)
               llcotf(ji,jpj  ) = llcotf(iju,jpj-2)
            END DO
            DO ji = jpi/2, jpi-1
               iju = jpi - ji
               llcotf(ji,jpjm1) = llcotf(iju,jpjm1)
            END DO
            DO ji = 1, jpi
               ijt = jpi - ji + 1
               llcotv(ji,jpj  ) = llcotv(ijt,jpj-1)
            END DO
            DO ji = jpi/2+1, jpi
               ijt = jpi - ji + 1
               llcotv(ji,jpjm1) = llcotv(ijt,jpjm1)
            END DO
         ENDIF

         ! Compute cartesian coordinates of coastline points
         ! and the number of coastline points

         icoast = 0
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( llcotf(ji,jj) ) THEN
                  icoast = icoast + 1
                  zxc(icoast) = COS( rad*gphif(ji,jj) ) * COS( rad*glamf(ji,jj) )
                  zyc(icoast) = COS( rad*gphif(ji,jj) ) * SIN( rad*glamf(ji,jj) )
                  zzc(icoast) = SIN( rad*gphif(ji,jj) )
               ENDIF
               IF( llcotu(ji,jj) ) THEN
                  icoast = icoast+1
                  zxc(icoast) = COS( rad*gphiu(ji,jj) ) * COS( rad*glamu(ji,jj) )
                  zyc(icoast) = COS( rad*gphiu(ji,jj) ) * SIN( rad*glamu(ji,jj) )
                  zzc(icoast) = SIN( rad*gphiu(ji,jj) )
               ENDIF
               IF( llcotv(ji,jj) ) THEN
                  icoast = icoast+1
                  zxc(icoast) = COS( rad*gphiv(ji,jj) ) * COS( rad*glamv(ji,jj) )
                  zyc(icoast) = COS( rad*gphiv(ji,jj) ) * SIN( rad*glamv(ji,jj) )
                  zzc(icoast) = SIN( rad*gphiv(ji,jj) )
               ENDIF
            END DO
         END DO

         ! Distance for the T-points

         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tmask(ji,jj,jk) == 0. ) THEN
                  pdct(ji,jj,jk) = 0.
               ELSE
                  DO jl = 1, icoast
                     zdis(jl) = ( zxt(ji,jj) - zxc(jl) )**2   &
                              + ( zyt(ji,jj) - zyc(jl) )**2   &
                              + ( zzt(ji,jj) - zzc(jl) )**2
                  END DO
                  pdct(ji,jj,jk) = ra * SQRT( MINVAL( zdis(1:icoast) ) )
               ENDIF
            END DO
         END DO
         !                                                ! ===============
      END DO                                              !   End of slab
      !                                                   ! ===============


      ! 2. Create the  distance to the coast file in NetCDF format
      ! ----------------------------------------------------------    
      clname = 'dist.coast'
      itime = 0
      CALL ymds2ju( 0     , 1     , 1     , 0.e0 , zdate0 )
      CALL restini( 'NONE', jpi   , jpj   , glamt, gphit ,   &
                    jpk   , gdept , clname, itime, zdate0,   &
                    rdt   , icot                         )
      CALL restput( icot, 'Tcoast', jpi, jpj, jpk, 0, pdct )
      CALL restclo( icot )

   END SUBROUTINE cofdis

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO internal damping
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_tradmp = .FALSE.    !: internal damping flag
CONTAINS
   SUBROUTINE tra_dmp( kt )        ! Empty routine
!      WRITE(*,*) 'tra_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE tra_dmp
#endif

   !!======================================================================
END MODULE tradmp
