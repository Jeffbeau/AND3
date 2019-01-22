MODULE trcdmp
   !!======================================================================
   !!                       ***  MODULE  trcdmp  ***
   !! Ocean physics: internal restoring trend on passive tracers
   !!======================================================================
#if  defined key_passivetrc && defined key_trcdmp 
   !!----------------------------------------------------------------------
   !!   key_trcdmp                                         internal damping
   !!----------------------------------------------------------------------
   !!   trc_dmp      : update the tracer trend with the internal damping
   !!   trc_dmp_init : initialization, namlist read, parameters control
   !!   trccof_zoom  : restoring coefficient for zoom domain
   !!   trccof       : restoring coefficient for global domain
   !!   cofdis       : compute the distance to the coastline
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE trctrp_lec      ! passive tracers transport
   USE trcdta
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_dmp   ! routine called by step.F90

   !! * Shared module variables
   LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .TRUE.    !: internal damping flag

   REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   &
      restotr         ! restoring coeff. on tracers (s-1)

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcdmp.F90,v 1.10 2006/04/10 15:38:54 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_dmp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trc_dmp  ***
      !!                  
      !! ** Purpose :   Compute the passive tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards trdta computed 
      !!      and add to the general tracer trends:
      !!                     trn = tra + restotr * (trdta - trb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmptr=0) or in area of weak vertical mixing (nlmdmptr=1) or
      !!      below the well mixed layer (nlmdmptr=2)
      !!
      !! ** Action  : - update the tracer trends tra with the newtonian 
      !!                damping trends.
      !!              - save the trends in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   7.0  !         (G. Madec)  Original code
      !!        !  96-01  (G. Madec) 
      !!        !  97-05  (H. Loukos)  adapted for passive tracers
      !!   8.5  !  02-08  (G. Madec )  free form + modules
      !!   9.0  !  04-03  (C. Ethe)    free form + modules
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk, jn     ! dummy loop indices
      REAL(wp) ::   ztest, ztra, zdt   ! temporary scalars
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      ! 0. Initialization (first time-step only)
      !    --------------
      IF( kt == nittrc000 ) CALL trc_dmp_init

      ! 1. Newtonian damping trends on tracer fields
      ! --------------------------------------------
      !    compute the newtonian damping trends depending on nmldmptr

!!!      zdt  = rdt * FLOAT( ndttrc )

      ! Initialize the input fields for newtonian damping
      CALL dta_trc( kt )

      DO jn = 1, jptra

         IF( lutini(jn) ) THEN

            SELECT CASE ( nmldmptr )

            CASE( 0 )                ! newtonian damping throughout the water column

               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
                        ztra = restotr(ji,jj,jk,jn) * ( trdta(ji,jj,jk,jn) - trb(ji,jj,jk,jn) )
                        ! add the trends to the general tracer trends
!!                        trn(ji,jj,jk,jn) = trn(ji,jj,jk,jn) + ztra * zdt
                        tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
#    if defined key_trc_diatrd
                        ! save the trends for diagnostics
                        IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),jpdiatrc) = ztra
#    endif
                     END DO
                  END DO
               END DO

            CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
                        ztest = avt(ji,jj,jk) - 5.e-4
                        IF( ztest < 0. ) THEN
                           ztra = restotr(ji,jj,jk,jn) * ( trdta(ji,jj,jk,jn) - trb(ji,jj,jk,jn) )
                        ELSE
                           ztra = 0.e0
                        ENDIF
                        ! add the trends to the general tracer trends
!!                        trn(ji,jj,jk,jn) = trn(ji,jj,jk,jn) + ztra * zdt
                        tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra 
#    if defined key_trc_diatrd
                        ! save the trends for diagnostics
                        IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),jpdiatrc) = ztra
#    endif
                     END DO
                  END DO
               END DO

            CASE ( 2 )                ! no damping in the mixed layer 
               DO jk = 1, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = fs_2, fs_jpim1   ! vector opt.
                        IF( fsdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                           ztra = restotr(ji,jj,jk,jn) * ( trdta(ji,jj,jk,jn) - trb(ji,jj,jk,jn) )
                        ELSE
                           ztra = 0.e0
                        ENDIF
                        ! add the trends to the general tracer trends
!!                        trn(ji,jj,jk,jn) = trn(ji,jj,jk,jn) + ztra * zdt
                        tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
#    if defined key_trc_diatrd
                        ! save the trends for diagnostics
                        IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),jpdiatrc) = ztra
#    endif
                     END DO
                  END DO
               END DO
               
            END SELECT

         ENDIF

      END DO

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('dmp')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF
  
      trb(:,:,:,:) = trn(:,:,:,:)
   
   END SUBROUTINE trc_dmp


   SUBROUTINE trc_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the newtonian damping 
      !!
      !! ** Method  :   read the nammbf namelist and check the parameters
      !!      called by trc_dmp at the first timestep (nit000)
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------

      SELECT CASE ( ndmptr )

      CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
         IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'

      CASE ( 1:90 )             ! Damping poleward of 'ndmptr' degrees
         IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmptr, ' degrees'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for ndmptr = ', ndmptr
         nstop = nstop + 1

      END SELECT


      SELECT CASE ( nmldmptr )

      CASE ( 0 )                ! newtonian damping throughout the water column
         IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'

      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'

      CASE ( 2 )                ! no damping in the mixed layer 
         IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for nmldmptr = ', nmldmptr
         nstop = nstop + 1

      END SELECT


      ! 3. Damping coefficients initialization
      ! --------------------------------------

         IF( lzoom ) THEN
            CALL trccof_zoom
         ELSE
            CALL trccof
         ENDIF
 
   END SUBROUTINE trc_dmp_init


   SUBROUTINE trccof_zoom
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trccof_zoom  ***
      !!
      !! ** Purpose :   Compute the damping coefficient for zoom domain
      !!
      !! ** Method  : - set along closed boundary due to zoom a damping over
      !!      6 points with a max time scale of 5 days.
      !!              - ORCA arctic/antarctic zoom: set the damping along
      !!      south/north boundary over a latitude strip.
      !!
      !! ** Action  : - restotr, the damping coeff. passive tracers
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!   9.0  !  04-03  (C. Ethe)   adapted for passive tracers
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn       ! dummy loop indices
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

      restotr(:,:,:,:) = 0.e0

      ! damping along the forced closed boundary over 6 grid-points
      DO jn = 1, 6
         IF( lzoom_w )   restotr( mi0(jn+jpizoom):mi1(jn+jpizoom), : , : , : ) = zfact(jn) ! west  closed
         IF( lzoom_s )   restotr( : , mj0(jn+jpjzoom):mj1(jn+jpjzoom), : , : ) = zfact(jn) ! south closed 
         IF( lzoom_e )   restotr( mi0(jpiglo+jpizoom-1-jn):mi1(jpiglo+jpizoom-1-jn) , : , : , : ) &
                       &              = zfact(jn)                                 ! east  closed 
         IF( lzoom_n )   restotr( : , mj0(jpjglo+jpjzoom-1-jn):mj1(jpjglo+jpjzoom-1-jn) , : , : ) &
                       &              = zfact(jn)                                 ! north closed
      END DO


      IF( lzoom_arct .AND. lzoom_anta ) THEN

         ! ====================================================
         !  ORCA configuration : arctic zoom or antarctic zoom
         ! ====================================================

         IF(lwp) WRITE(numout,*)
         IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              trccof_zoom : ORCA    Arctic zoom'
         IF(lwp .AND. lzoom_arct ) WRITE(numout,*) '              trccof_zoom : ORCA Antarctic zoom'
         IF(lwp) WRITE(numout,*)

         ! ... Initialization : 
         !     zlat0 : latitude strip where resto decreases
         !     zlat1 : resto = 1 before zlat1
         !     zlat2 : resto decreases from 1 to 0 between zlat1 and zlat2
         restotr(:,:,:,:) = 0.e0
         zlat0 = 10.
         zlat1 = 30.
         zlat2 = zlat1 + zlat0

         ! ... Compute arrays resto ; value for internal damping : 5 days
         DO jn = 1, jptra
            DO jk = 2, jpkm1 
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zlat = ABS( gphit(ji,jj) )
                     IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
                        restotr(ji,jj,jk,jn) = 0.5 * ( 1./(5.*rday) ) *   &
                           ( 1. - COS(rpi*(zlat2-zlat)/zlat0) ) 
                     ELSE IF ( zlat < zlat1 ) THEN
                        restotr(ji,jj,jk,jn) = 1./(5.*rday)
                     ENDIF
                  END DO
               END DO
            END DO
         END DO

      ENDIF

      ! ... Mask resto array
        DO jn = 1, jptra
           restotr(:,:,:,jn) = restotr(:,:,:,jn) * tmask(:,:,:)
        END DO


   END SUBROUTINE trccof_zoom

   SUBROUTINE trccof
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trccof  ***
      !!
      !! ** Purpose :   Compute the damping coefficient
      !!
      !! ** Method  :   Arrays defining the damping are computed for each grid
      !!      point passive tracers (restotr)
      !!      Damping depends on distance to coast, depth and latitude
      !!
      !! ** Action  : - restotr, the damping coeff. for passive tracers
      !!
      !! History :
      !!   5.0  !  91-03  (O. Marti, G. Madec)  Original code
      !!        !  92-06  (M. Imbard)  doctor norme
      !!        !  96-01  (G. Madec) statement function for e3
      !!        !  98-07  (M. Imbard, G. Madec) ORCA version
      !!        !  00-08  (G. Madec, D. Ludicone) 
      !!   8.2  !  04-03  (H. Loukos) adapted for passive tracers
      !!        !  04-02  (O. Aumont, C. Ethe) rewritten for debuging and update
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Local declarations
      INTEGER ::  ji, jj, jk, je, jn     ! dummy loop indices
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
      !!----------------------------------------------------------------------

      ! ====================================
      !  ORCA configuration : global domain
      ! ====================================

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '              trccof : Global domain of ORCA'
      IF(lwp) WRITE(numout,*) '              ------------------------------'


      ! ... Initialization : 
      !   zdct()      : distant to the coastline
      !   resto()     : array of restoring coeff. 
      
      zdct (:,:,:) = 0.e0
      restotr(:,:,:,:) = 0.e0


      IF ( ndmptr > 0 ) THEN

         !    ------------------------------------
         !     Damping poleward of 'ndmptr' degrees
         !    ------------------------------------

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              Damping poleward of ', ndmptr,' deg.'
         IF(lwp) WRITE(numout,*)

         ! ... Distance to coast (zdct)

         !   ... Test the existance of distance-to-coast file
         itime = jpmois
         ipi = jpiglo
         ipj = jpjglo
         ipk = jpk
         clname = 'dist.coast.trc'
         DO je = 1,32
            IF( clname(je:je) == ' ' ) go to 140
         END DO
140      CONTINUE
         ie = je
         clname2 = clname(1:ie-1)//".nc"
         INQUIRE( FILE = clname2, EXIST = llbon )

         IF ( llbon ) THEN

            !   ... Read file distance to coast if possible
            CALL flinopen( clname, mig(1), nlci, mjg(1), nlcj, .FALSE.,   &
               ipi, ipj, ipk, zlamt, zphit, zdept, jpmois,   &
               istep, zdate0, rdt, icot )
            CALL flinget( icot, 'Tcoast', jpidta, jpjdta, jpk,    &
               jpmois, 1, 1, mig(1), nlci, mjg(1), nlcj, zdct(1:nlci,1:nlcj,1:jpk) )
            CALL flinclo( icot )
            IF(lwp)WRITE(numout,*) '    ** : File trc.dist.coast.nc read'

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
         zlat1 = ndmptr
         zlat2 = zlat1 + zlat0

         DO jn = 1, jptra
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zlat = ABS( gphit(ji,jj) )
                  IF ( zlat1 <= zlat .AND. zlat <= zlat2 ) THEN
                     restotr(ji,jj,1,jn) = 0.5 * ( 1. - COS(rpi*(zlat-zlat1)/zlat0 ) )
                  ELSEIF ( zlat > zlat2 ) THEN
                     restotr(ji,jj,1,jn) = 1.
                  ENDIF
               END DO
            END DO
         END DO

         !   ... North Indian ocean (20N/30N x 45E/100E) : resto=0
         IF ( ndmptr == 20 ) THEN
            DO jn = 1, jptra
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zlat = gphit(ji,jj)
                     zlon = MOD( glamt(ji,jj), 360. )
                     IF ( zlat1 < zlat .AND. zlat < zlat2 .AND.   &
                        45.  < zlon .AND. zlon < 100. ) THEN
                        restotr(ji,jj,1,jn) = 0.
                     ENDIF
                  END DO
               END DO
            END DO
         ENDIF

         zsdmp = 1./(sdmptr * rday)
         zbdmp = 1./(bdmptr * rday)
         DO jn = 1, jptra
            DO jk = 2, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zdct(ji,jj,jk) = MIN( zinfl, zdct(ji,jj,jk) )

                     !   ... Decrease the value in the vicinity of the coast
                     restotr(ji,jj,jk,jn) = restotr(ji,jj,1,jn)*0.5   &
                        &                 * ( 1. - COS( rpi*zdct(ji,jj,jk)/zinfl) )

                     !   ... Vertical variation from zsdmp (sea surface) to zbdmp (bottom)
                     restotr(ji,jj,jk,jn) = restotr(ji,jj,jk,jn)   &
                        &                 * ( zbdmp + (zsdmp-zbdmp)*EXP(-fsdept(ji,jj,jk)/hdmptr) )
                  END DO
               END DO
            END DO
         END DO

      ENDIF


      IF( cp_cfg == "orca" .AND. ( ndmptr > 0 .OR. ndmptr == -1 ) ) THEN

         !                                         ! =========================
         !                                         !  Med and Red Sea damping
         !                                         ! =========================
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) '              ORCA configuration: Damping in Med and Red Seas'


         zmrs(:,:) = 0.e0                             ! damping term on the Med or Red Sea

         SELECT CASE ( jp_cfg )
            !                                           ! =======================
         CASE ( 4 )                                     !  ORCA_R4 configuration 
            !                                           ! =======================

            ! Mediterranean Sea
            ij0 =  50   ;   ij1 =  56
            ii0 =  81   ;   ii1 =  91   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ij0 =  50   ;   ij1 =  55
            ii0 =  75   ;   ii1 =  80   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ij0 =  52   ;   ij1 =  53
            ii0 =  70   ;   ii1 =  74   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Smooth transition from 0 at surface to 1./rday at the 18th level in Med and Red Sea
            DO jk = 1, 17
               zhfac (jk) = 0.5*( 1.- COS( rpi*(jk-1)/16. ) ) / rday
            END DO
            DO jk = 18, jpkm1
               zhfac (jk) = 1./rday
            END DO

            !                                        ! =======================
         CASE ( 2 )                                  !  ORCA_R2 configuration 
            !                                        ! =======================

            ! Mediterranean Sea
            ij0 =  96   ;   ij1 = 110
            ii0 = 157   ;   ii1 = 181   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ij0 = 100   ;   ij1 = 110
            ii0 = 144   ;   ii1 = 156   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ij0 = 100   ;   ij1 = 103
            ii0 = 139   ;   ii1 = 143   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Decrease before Gibraltar Strait
            ij0 = 101   ;   ij1 = 102
            ii0 = 139   ;   ii1 = 141   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.e0
            ii0 = 142   ;   ii1 = 142   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0 / 90.e0
            ii0 = 143   ;   ii1 = 143   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40e0
            ii0 = 144   ;   ii1 = 144   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.75e0
            ! Red Sea
            ij0 =  87   ;   ij1 =  96
            ii0 = 147   ;   ii1 = 163   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Decrease before Bab el Mandeb Strait
            ij0 =  91   ;   ij1 =  91
            ii0 = 153   ;   ii1 = 160   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.80e0
            ij0 =  90   ;   ij1 =  90
            ii0 = 153   ;   ii1 = 160   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40e0
            ij0 =  89   ;   ij1 =  89
            ii0 = 158   ;   ii1 = 160   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0 / 90.e0
            ij0 =  88   ;   ij1 =  88
            ii0 = 160   ;   ii1 = 163   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.e0
            ! Smooth transition from 0 at surface to 1./rday at the 18th level in Med and Red Sea
            DO jk = 1, 17
               zhfac (jk) = 0.5*( 1.- COS( rpi*(jk-1)/16. ) ) / rday
            END DO
            DO jk = 18, jpkm1
               zhfac (jk) = 1./rday
            END DO

            !                                        ! =======================
         CASE ( 05 )                                 !  ORCA_R05 configuration
            !                                        ! =======================

            ! Mediterranean Sea
            ii0 = 568   ;   ii1 = 574 
            ij0 = 324   ;   ij1 = 333   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ii0 = 575   ;   ii1 = 658
            ij0 = 314   ;   ij1 = 366   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Black Sea (remaining part
            ii0 = 641   ;   ii1 = 651
            ij0 = 367   ;   ij1 = 372   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Decrease before Gibraltar Strait
            ii0 = 324   ;   ii1 = 333
            ij0 = 565   ;   ij1 = 565   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0 / 90.e0
            ij0 = 566   ;   ij1 = 566   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.40
            ij0 = 567   ;   ij1 = 567   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 0.75
            ! Red Sea
            ii0 = 641   ;   ii1 = 665
            ij0 = 270   ;   ij1 = 310   ;   zmrs( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) ) = 1.e0
            ! Decrease before Bab el Mandeb Strait
            ii0 = 666   ;   ii1 = 675
            ij0 = 270   ;   ij1 = 290   
            DO ji = mi0(ii0), mi1(ii1)
               zmrs( ji , mj0(ij0):mj1(ij1) ) = 0.1 * ABS( FLOAT(ji - mi1(ii1)) )
            END DO
            zsdmp = 1./(sdmptr * rday)
            zbdmp = 1./(bdmptr * rday)
            DO jk = 1, jpk
               zhfac (jk) = ( zbdmp + (zsdmp-zbdmp) * EXP(-fsdept(1,1,jk)/hdmptr) )
            END DO

            !                                       ! ========================
         CASE ( 025 )                               !  ORCA_R025 configuration 
            !                                       ! ========================
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*)' Not yet implemented in ORCA_R025'
            nstop = nstop + 1

         END SELECT

         DO jn = 1, jptra
            DO jk = 1, jpkm1
               restotr(:,:,jk,jn) = zmrs(:,:) * zhfac(jk) + ( 1. - zmrs(:,:) ) * restotr(:,:,jk,jn)
            END DO

            ! Mask resto array and set to 0 first and last levels
            restotr(:,:, : ,jn) = restotr(:,:,:,jn) * tmask(:,:,:)
            restotr(:,:, 1 ,jn) = 0.e0
            restotr(:,:,jpk,jn) = 0.e0
         END DO

      ELSE
         !    ------------
         !     No damping
         !    ------------
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) 'Choose a correct value of ndmptr or DO NOT defined key_trcdmp'
         nstop = nstop + 1
      ENDIF

        !    ----------------------------
         !     Create Print damping array
         !    ----------------------------
         
         ! ndmpftr   : = 1 create a damping.coeff NetCDF file

      IF( ndmpftr == 1 ) THEN
         DO jn = 1, jptra
            IF(lwp) WRITE(numout,*) '  create damping.coeff.nc file  ',jn
            itime   = 0
            clname3 = 'damping.coeff'//ctrcnm(jn)
            CALL ymds2ju( 0     , 1     , 1      , 0.e0 , zdate0 )
            CALL restini( 'NONE', jpi   , jpj    , glamt, gphit,    &
           &              jpk   , gdept , clname3, itime, zdate0,   &
           &              rdt   , idmp  , domain_id=nidom)
            CALL restput( idmp, 'Resto', jpi, jpj, jpk, 0 , restotr(:,:,:,jn)  )
            CALL restclo( idmp )
         END DO
      ENDIF


   END SUBROUTINE trccof


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
      !!              - NetCDF file 'trc.dist.coast.nc' 
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
      clname = 'trc.dist.coast'
      itime = 0
      CALL ymds2ju( 0     , 1     , 1     , 0.e0 , zdate0 )
      CALL restini( 'NONE', jpi   , jpj   , glamt, gphit ,   &
                    jpk   , gdept , clname, itime, zdate0,   &
                    rdt   , icot , domain_id=nidom         )
      CALL restput( icot, 'Tcoast', jpi, jpj, jpk, 0, pdct )
      CALL restclo( icot )

   END SUBROUTINE cofdis

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO internal damping
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
CONTAINS
   SUBROUTINE trc_dmp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dmp
#endif

   !!======================================================================
END MODULE trcdmp
