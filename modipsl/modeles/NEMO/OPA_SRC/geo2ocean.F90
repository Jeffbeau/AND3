MODULE geo2ocean
   !!======================================================================
   !!                     ***  MODULE  geo2ocean  ***
   !! Ocean mesh    :  ???
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   repcmo      : 
   !!   angle       :
   !!   geo2oce     :
   !!   repere      :   old routine suppress it ???
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! mesh and scale factors
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE

   !! * Accessibility
   PRIVATE
   PUBLIC repcmo   ! routine called by ???.F90
   PUBLIC geo2oce  ! routine called by ???.F90
   PUBLIC repere   ! routine called by ???.F90

   !! * Module variables
   REAL(wp), DIMENSION(jpi,jpj) ::   &
      gsinu , gcosu ,   &  ! matrix element for change grid u (repcmo.F)
      gsinv , gcosv ,   &  ! matrix element for change grid v (repcmo.F)
      gsinus, gcosin      ! matrix element for change grid (repere.F)

  !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/geo2ocean.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE repcmo ( pxu1, pyu1, pxv1, pyv1,   &
                       px2 , py2 , kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE repcmo  ***
      !!
      !! ** Purpose :   Change vector componantes from a geographic grid to a
      !!      stretched coordinates grid.
      !!
      !! ** Method  :   Initialization of arrays at the first call.
      !!
      !! ** Action  : - px2 : first componante (defined at u point)
      !!              - py2 : second componante (defined at v point)
      !!
      !! History :
      !!   7.0  !  07-96  (O. Marti)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Arguments 
      INTEGER,  INTENT( in ) ::   &
         kt                ! ocean time-step
      REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   & 
         pxu1, pyu1,     & ! geographic vector componantes at u-point
         pxv1, pyv1        ! geographic vector componantes at v-point
      REAL(wp), INTENT( out ), DIMENSION(jpi,jpj) ::   & 
         px2,            & ! i-componante (defined at u-point)
         py2               ! j-componante (defined at v-point)
      !!----------------------------------------------------------------------


      ! Initialization of gsin* and gcos* at first call
      ! -----------------------------------------------

      IF( kt <= nit000 + 1 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'repcmo : use the geographic to stretched'
         IF(lwp) WRITE(numout,*) ' ~~~~~   coordinate transformation'

         CALL angle       ! initialization of the transformation
      ENDIF
      
      ! Change from geographic to stretched coordinate
      ! ----------------------------------------------
      
      px2(:,:) = pxu1(:,:) * gcosu(:,:) + pyu1(:,:) * gsinu(:,:)
      py2(:,:) = pyv1(:,:) * gcosv(:,:) - pxv1(:,:) * gsinv(:,:)   
      
   END SUBROUTINE repcmo


   SUBROUTINE angle
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE angle  ***
      !! 
      !! ** Purpose :   Compute angles between model grid lines and the 
      !!      direction of the North
      !!
      !! ** Method  :
      !!
      !! ** Action  :   Compute (gsinu, gcosu, gsinv, gcosv) arrays: sinus and 
      !!      cosinus of the angle between the north-south axe and the 
      !!      j-direction at u and v-points
      !!
      !! History :
      !!   7.0  !  96-07  (O. Marti)  Original code
      !!   8.0  !  98-06  (G. Madec)
      !!   8.5  !  98-06  (G. Madec)  Free form, F90 + opt.
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::   ji, jj      ! dummy loop indices

      REAL(wp) ::   &
         zlam, zphi,             &  ! temporary scalars
         zlan, zphh,             &  !    "         "
         zxnpu, zxnpv , znnpu,   &  !    "         "
         zynpu, zynpv , znnpv,   &  !    "         "
         zxffu, zmnpfu, zxffv,   &  !    "         "
         zyffu, zmnpfv, zyffv       !    "         "
      !!----------------------------------------------------------------------

      ! ============================= !
      ! Compute the cosinus and sinus !
      ! ============================= !
      ! (computation done on the north stereographic polar plane)

      DO jj = 2, jpj
!CDIR NOVERRCHK
         DO ji = fs_2, jpi   ! vector opt.

            ! north pole direction & modulous (at u-point)
            zlam = glamu(ji,jj)
            zphi = gphiu(ji,jj)
            zxnpu = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpu = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu

            ! north pole direction & modulous (at v-point)
            zlam = glamv(ji,jj)
            zphi = gphiv(ji,jj)
            zxnpv = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpv = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv

            ! j-direction: f-point segment direction (u-point)
            zlam = glamf(ji,jj  )
            zphi = gphif(ji,jj  )
            zlan = glamf(ji,jj-1)
            zphh = gphif(ji,jj-1)
            zxffu =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffu =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zmnpfu = SQRT ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            zmnpfu = MAX( zmnpfu, 1.e-14 )

            ! i-direction: f-point segment direction (v-point)
            zlam = glamf(ji  ,jj)
            zphi = gphif(ji  ,jj)
            zlan = glamf(ji-1,jj)
            zphh = gphif(ji-1,jj)
            zxffv =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffv =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zmnpfv = SQRT ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            zmnpfv = MAX( zmnpfv, 1.e-14 )

            ! cosinus and sinus using scalar and vectorial products
            gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
            gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu

            ! cosinus and sinus using scalar and vectorial products
            ! (caution, rotation of 90 degres)
            gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
            gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv

         END DO
      END DO

      ! =============== !
      ! Geographic mesh !
      ! =============== !

      DO jj = 2, jpj
         DO ji = fs_2, jpi   ! vector opt.
            IF( ABS( glamf(ji,jj) - glamf(ji,jj-1) ) < 1.e-8 ) THEN
               gsinu(ji,jj) = 0.
               gcosu(ji,jj) = 1.
            ENDIF
            IF( ABS( gphif(ji,jj) - gphif(ji-1,jj) ) < 1.e-8 ) THEN
               gsinv(ji,jj) = 0.
               gcosv(ji,jj) = 1.
            ENDIF
         END DO
      END DO

      ! =========================== !
      ! Lateral boundary conditions !
      ! =========================== !

      ! lateral boundary cond.: U-, V-pts, sgn
      CALL lbc_lnk ( gsinu, 'U', -1. )   ;   CALL lbc_lnk( gsinv, 'V', -1. )
      CALL lbc_lnk ( gcosu, 'U', -1. )   ;   CALL lbc_lnk( gcosv, 'V', -1. )

   END SUBROUTINE angle


   SUBROUTINE geo2oce ( pxx , pyy , pzz, cgrid,     &
                        plon, plat, pte, ptn  , ptv )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE geo2oce  ***
      !!      
      !! ** Purpose :
      !!
      !! ** Method  :   Change wind stress from geocentric to east/north
      !!
      !! History :
      !!        !         (O. Marti)  Original code
      !!        !  91-03  (G. Madec)
      !!        !  92-07  (M. Imbard)
      !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
      !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
      !!   8.5  !  02-06  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Local declarations
      REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
         pxx, pyy, pzz
      CHARACTER (len=1), INTENT( in) ::   &
         cgrid
      REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
         plon, plat
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::    &
         pte, ptn, ptv
      REAL(wp), PARAMETER :: rpi = 3.141592653E0
      REAL(wp), PARAMETER :: rad = rpi / 180.e0

      !! * Local variables
      INTEGER ::   ig     !

      !! * Local save
      REAL(wp), SAVE, DIMENSION(jpi,jpj,4) ::   &
         zsinlon, zcoslon,   &
         zsinlat, zcoslat
      LOGICAL, SAVE, DIMENSION (4) ::   &
         linit = .FALSE.
      !!----------------------------------------------------------------------

      SELECT CASE( cgrid)

         CASE ( 't' ) ;; ig = 1
         CASE ( 'u' ) ;; ig = 2
         CASE ( 'v' ) ;; ig = 3
         CASE ( 'f' ) ;; ig = 4

         CASE default
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) 'geo2oce : bad grid argument : ', cgrid
            nstop = nstop + 1
       END SELECT
      
      IF( .NOT. linit(ig) ) THEN 
         zsinlon (:,:,ig) = SIN (rad * plon)
         zcoslon (:,:,ig) = COS (rad * plon)
         zsinlat (:,:,ig) = SIN (rad * plat)
         zcoslat (:,:,ig) = COS (rad * plat)
         linit (ig) = .TRUE.
      ENDIF
      
      pte = - zsinlon (:,:,ig) * pxx + zcoslon (:,:,ig) * pyy
      ptn = - zcoslon (:,:,ig) * zsinlat (:,:,ig) * pxx    &
            - zsinlon (:,:,ig) * zsinlat (:,:,ig) * pyy    &
            + zcoslat (:,:,ig) * pzz
      ptv =   zcoslon (:,:,ig) * zcoslat (:,:,ig) * pxx    &
            + zsinlon (:,:,ig) * zcoslat (:,:,ig) * pyy    &
            + zsinlat (:,:,ig) * pzz

   END SUBROUTINE geo2oce


   SUBROUTINE repere ( px1, py1, px2, py2, kchoix )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE repere  ***
      !!        
      !! ** Purpose :   Change vector componantes between a geopgraphic grid 
      !!      and a stretched coordinates grid.
      !!
      !! ** Method  :   initialization of arrays at the first call.
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  89-03  (O. Marti)  original code
      !!        !  92-02  (M. Imbard)
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  98-05  (B. Blanke)
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT( in   ), DIMENSION(jpi,jpj) ::   &
         px1, py1          ! two horizontal components to be rotated
      REAL(wp), INTENT( out  ), DIMENSION(jpi,jpj) ::   &
         px2, py2          ! the two horizontal components in the model repere
      INTEGER, INTENT( inout ) ::   &
         kchoix   ! type of transformation
                  ! = 1 change from geographic to model grid.
                  ! =-1 change from model to geographic grid
                  ! = 0 same as the previous call
      !! * Local declarations
      INTEGER, SAVE :: nmem

      INTEGER ::   ji, jj                    ! dummy loop indices

      REAL(wp) :: zxx, zcof1, zcof2,    &
         ze1t, ze2t
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zlamdu, zphiu,   &
         zlamdv, zphiv
      !!----------------------------------------------------------------------


      ! 0. Initialization of gsinus and gcosin IF first call
      ! ----------------------------------------------------
      
      ! 0.1 control argument
      
      IF( kchoix == 0 ) THEN
         IF( nmem == 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
            IF(lwp) WRITE(numout,*) ' for the first call , you must indicate '
            IF(lwp) WRITE(numout,*) ' the direction of change '
            IF(lwp) WRITE(numout,*) ' kchoix = 1 geo       --> stretched '
            IF(lwp) WRITE(numout,*) ' kchoix =-1 stretched --> geo '
            nstop = nstop + 1
         ELSE
            kchoix = nmem
         ENDIF
      ELSEIF( kchoix == 1 .OR. kchoix == -1 ) THEN
         nmem = kchoix
      ELSE
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
         IF(lwp) WRITE(numout,*) ' kchoix must be equal to -1, 0 or 1 '
         nstop = nstop + 1
      ENDIF

      ! 0.2 Initialization

      zxx = gsinus(jpi/2,jpj/2)**2+gcosin(jpi/2,jpj/2)**2
      IF( zxx <= 0.1 ) THEN
         IF(lwp) WRITE(numout,*) 'repere : initialization '
         DO jj = 1, jpj
            DO ji = 2, jpi
               zlamdu(ji,jj) = glamu(ji,jj) - glamu(ji-1,jj)
               zlamdu(ji,jj) = ASIN( SIN( rad*zlamdu(ji,jj) ) )/rad
               zphiu (ji,jj) = gphiu(ji,jj) - gphiu(ji-1,jj)
            END DO
         END DO
         DO jj = 2, jpj
            zlamdv(:,jj) = glamv(:,jj)-glamv(:,jj-1)
            zlamdv(:,jj) = ASIN(SIN(rad*zlamdv(:,jj)))/rad
            zphiv (:,jj) = gphiv(:,jj)-gphiv(:,jj-1)
         END DO
         
         ! 0.3 Boudary conditions and periodicity
         
         IF( nperio == 1 .OR.nperio == 4.OR.nperio == 6 ) THEN
            DO jj = 1, jpj
               zlamdu(1,jj) = zlamdu(jpim1,jj)
               zphiu (1,jj) = zphiu (jpim1,jj)
            END DO
         ELSE
            DO jj = 1, jpj
               zlamdu(1,jj) = zlamdu(2,jj)
               zphiu (1,jj) = zphiu (2,jj)
            END DO
         ENDIF
         
         DO ji = 1, jpi
            zlamdv(ji,1) = zlamdv(ji,2)
            zphiv (ji,1) = zphiv (ji,2)
         END DO
         
         IF( nperio == 2 ) THEN
            DO ji = 1, jpi
               zphiv (ji,1) = zphiv (ji,3)
            END DO
         ENDIF
         
         ! 0.4 gsinus gcosin
         
!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               zcof1 = rad * ra * COS( rad * gphit(ji,jj) )
               zcof2 = rad * ra
               zlamdu(ji,jj) = zlamdu(ji,jj) * zcof1
               zlamdv(ji,jj) = zlamdv(ji,jj) * zcof1
               zphiu (ji,jj) = zphiu (ji,jj) * zcof2
               zphiv (ji,jj) = zphiv (ji,jj) * zcof2
            END DO
         END DO

!CDIR NOVERRCHK
         DO jj = 1, jpj
!CDIR NOVERRCHK
            DO ji = 1, jpi
               ze1t = SQRT( zlamdu(ji,jj)*zlamdu(ji,jj) + zphiu(ji,jj)*zphiu(ji,jj) )
               ze2t = SQRT( zlamdv(ji,jj)*zlamdv(ji,jj) + zphiv(ji,jj)*zphiv(ji,jj) )
               gsinus(ji,jj) = 0.5*( zphiu(ji,jj)/ze1t - zlamdv(ji,jj)/ze2t )
               gcosin(ji,jj) = 0.5*( zphiv(ji,jj)/ze2t + zlamdu(ji,jj)/ze1t )
            END DO
         END DO
         CALL lbc_lnk( gsinus, 'T', -1. )
         CALL lbc_lnk( gcosin, 'T', -1. )
         
      ENDIF


      ! 1. Change from geographic to tretched
      ! -------------------------------------
      
      IF( kchoix == 1 ) THEN
          px2(:,:) =  px1(:,:) * gcosin(:,:) + py1(:,:) * gsinus(:,:)
          py2(:,:) = -px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
      ENDIF
      

      ! 2. Change from tretched to geographic
      ! -------------------------------------
      
      IF( kchoix == -1 ) THEN
          px2(:,:) =  px1(:,:) * gcosin(:,:) - py1(:,:) * gsinus(:,:)
          py2(:,:) =  px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
      ENDIF
      
   END SUBROUTINE repere

  !!======================================================================
END MODULE geo2ocean
