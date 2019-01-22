MODULE domrea
   !!======================================================================
   !!                       ***  MODULE domrea  ***
   !! Ocean initialization : read the ocean domain meshmask file(s)
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dom_rea        : read mesh and mask file(s)
   !!                    nmsh = 1  :   mesh_mask file
   !!                         = 2  :   mesh and mask file
   !!                         = 3  :   mesh_hgr, mesh_zgr and mask
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dom_rea        ! routine called by inidom.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/DOM/domrea.F90,v 1.5 2006/04/11 13:48:18 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_fdir
   !!----------------------------------------------------------------------
   !!   'key_fdir' :                                     direct access file
   !!----------------------------------------------------------------------
#  include "domrea_fdir.h90"

#elif ( defined key_mpp_mpi || defined key_mpp_shmem ) && defined key_dimgout
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'     OR
   !!   'key_mpp_shmem'
   !!   'key_dimgout' :         each processor makes its own direct access file 
   !!                      use build_nc_meshmask off line to retrieve 
   !!                      a ioipsl compliant meshmask file
   !!----------------------------------------------------------------------
#  include "domrea_dimg.h90"


#else
   !!----------------------------------------------------------------------
   !!   Default option :                                        NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE dom_rea
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_rea  ***
      !!                   
      !! ** Purpose :  Read the NetCDF file(s) which contain(s) all the
      !!      ocean domain informations (mesh and mask arrays). This (these)
      !!      file(s) is (are) used for visualisation (SAXO software) and
      !!      diagnostic computation.
      !!
      !! ** Method  :   Read in a file all the arrays generated in routines
      !!      domhgr, domzgr, and dommsk. Note: the file contain depends on
      !!      the vertical coord. used (z-coord, partial steps, s-coord)
      !!                    nmsh = 1  :   'mesh_mask.nc' file
      !!                         = 2  :   'mesh.nc' and mask.nc' files
      !!                         = 3  :   'mesh_hgr.nc', 'mesh_zgr.nc' and
      !!                                  'mask.nc' files
      !!      For huge size domain, use option 2 or 3 depending on your 
      !!      vertical coordinate.
      !!
      !! ** input file : 
      !!      meshmask.nc  : domain size, horizontal grid-point position,
      !!                     masks, depth and vertical scale factors
      !!
      !! History :
      !!        !  97-02  (G. Madec)  Original code
      !!        !  99-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
      !!   9.0  !  02-08  (G. Madec)  F90 and several file
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Local declarations
      LOGICAL ::   llog
      INTEGER  ::   ji, jj, jk, ik
      INTEGER  ::                & !!! * temprary units for :
         inum0 ,                 &  ! 'mesh_mask.nc' file
         inum1 ,                 &  ! 'mesh.nc'      file
         inum2 ,                 &  ! 'mask.nc'      file
         inum3 ,                 &  ! 'mesh_hgr.nc'  file
         inum4                      ! 'mesh_zgr.nc'  file
      INTEGER  ::   itime           !  output from restini ???
      REAL(wp) ::   zdate0, zdt
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
         zta, zlamt, zphit       ! dummy array for bathymetry 
      REAL(wp) , DIMENSION(jpidta,jpjdta,jpk) :: &
         zt3a      ! dummy array for bathymetry 
      REAL(wp), DIMENSION(jpi,jpj) :: &
         zprt = 0.

      CHARACTER (len=21) ::      &
         clnam0 = 'mesh_mask',   &  ! filename (mesh and mask informations)
         clnam1 = 'mesh'     ,   &  ! filename (mesh informations)
         clnam2 = 'mask'     ,   &  ! filename (mask informations)
         clnam3 = 'mesh_hgr' ,   &  ! filename (horizontal mesh informations)
         clnam4 = 'mesh_zgr'        ! filename (vertical   mesh informations)
      !!----------------------------------------------------------------------

       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,*) 'dom_rea : read NetCDF mesh and mask information file(s)'
       IF(lwp) WRITE(numout,*) '~~~~~~~'

      llog  = .FALSE.
      zlamt(:,:) = 0.e0
      zphit(:,:) = 0.e0

      CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )    ! calendar initialization

!       note that mbathy has been modified in dommsk or in solver.
!       it is the number of non-zero "w" levels in the water, and the minimum 
!       value (on land) is 2. We define zprt as the number of "T" points in the ocean 
!       at any location, and zero on land. 
!

      SELECT CASE (nmsh)
         !                                     ! ============================
         CASE ( 1 )                            !  create 'mesh_mask.nc' file
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          one file in "mesh_mask.nc" '
            CALL restini( clnam0, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mesh_mask.nc' file
            &             jpk   , gdept , trim(clnam0)        ,  &   ! in unit inum0
            &             itime , zdate0, zdt   , inum0, domain_id=nidom )
            inum2 = inum0                                            ! put all the informations
            inum3 = inum0                                            ! in unit inum0
            inum4 = inum0

            !                                  ! ============================
         CASE ( 2 )                            !  create 'mesh.nc' and 
            !                                  !         'mask.nc' files
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          two files in "mesh.nc" and "mask.nc" '
            CALL restini( clnam1, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mesh.nc' file 
            &             jpk   , gdept , trim(clnam1)        ,  &   ! in unit inum1 
            &             itime , zdate0, zdt   , inum1, domain_id=nidom )
            CALL restini( clnam2, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mask.nc' file 
            &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2 
            &             itime , zdate0, zdt   , inum2, domain_id=nidom )
            inum3 = inum1                                            ! put mesh informations 
            inum4 = inum1                                            ! in unit inum1 

            !                                  ! ============================
         CASE ( 3 )                            !  create 'mesh_hgr.nc'
            !                                  !         'mesh_zgr.nc' and
            !                                  !         'mask.nc'     files
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          three files in "mesh_hgr.nc" , mesh_zgr.nc" and "mask.nc" '
            CALL restini( clnam3, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mesh_hgr.nc' file
            &             jpk   , gdept , trim(clnam3)        ,  &   ! in unit inum3
            &             itime , zdate0, zdt   , inum3, domain_id=nidom )
            CALL restini( clnam4, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mesh_zgr.nc' file
            &             jpk   , gdept , trim(clnam4)        ,  &   ! in unit inum4
            &             itime , zdate0, zdt   , inum4, domain_id=nidom )
            CALL restini( clnam2, jpidta   , jpjdta   , zlamt, zphit,  &   ! create 'mask.nc' file
            &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2
            &             itime , zdate0, zdt   , inum2, domain_id=nidom )

         END SELECT

         !                                                         ! masks (inum2) 
         CALL restget( inum2, 'tmask', jpidta, jpjdta, jpk, 0, llog, zt3a ) 
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               tmask(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum2, 'umask', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               umask(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum2, 'vmask', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               vmask(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum2, 'fmask', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               fmask(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO

#if defined key_cfg_1d
      IF(lwp) WRITE(numout,*) '**********  1D configuration : set umask and vmask equal tmask ********'
      IF(lwp) WRITE(numout,*) '**********                                                     ********'
      ! set umask and vmask equal tmask in 1D configuration
      umask(:,:,:) = tmask(:,:,:)
      vmask(:,:,:) = tmask(:,:,:)
#endif

#if defined key_off_degrad
         CALL restget( inum2, 'facvolt', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               facvol(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
#endif

         !                                                         ! horizontal mesh (inum3)
         CALL restget( inum3, 'glamt', jpidta, jpjdta, 1, 0, llog, zta )     !    ! latitude
           DO jj = 1, nlcj
             DO ji = 1, nlci
               glamt(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'glamu', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               glamu(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'glamv', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               glamv(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'glamf', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               glamf(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         CALL restget( inum3, 'gphit', jpidta, jpjdta, 1, 0, llog, zta )     !    ! longitude
           DO jj = 1, nlcj
             DO ji = 1, nlci
               gphit(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'gphiu', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               gphiu(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'gphiv', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               gphiv(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'gphif', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               gphif(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         CALL restget( inum3, 'e1t', jpidta, jpjdta, 1, 0, llog, zta )         !    ! e1 scale factors
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e1t(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'e1u', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e1u(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'e1v', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e1v(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'e2t', jpidta, jpjdta, 1, 0, llog, zta )         !    ! e2 scale factors
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e2t(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'e2u', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e2u(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'e2v', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e2v(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum3, 'ff', jpidta, jpjdta, 1, 0, llog, zta )           !    ! coriolis factor
           DO jj = 1, nlcj
             DO ji = 1, nlci
               ff(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         CALL restget( inum4, 'mbathy', jpidta, jpjdta, 1, 0, llog, zta )
! Bathymetry
           DO jj = 1, nlcj
             DO ji = 1, nlci
               zprt(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         mbathy(:,:)=zprt(:,:)*tmask(:,:,1)+1

# if defined key_s_coord
         !                                                         ! s-coordinate
         CALL restget( inum4, 'hbatt', jpidta, jpjdta, 1, 0, llog, zta )      !    ! depth
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hbatt(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum4, 'hbatu', jpidta, jpjdta, 1, 0, llog, zta ) 
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hbatu(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum4, 'hbatv', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hbatv(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum4, 'hbatf', jpidta, jpjdta, 1, 0, llog, zta )
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hbatf(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         CALL restget( inum4, 'gsigt', 1, 1, jpk, 0, llog, gsigt )        !    ! scaling coef.
         CALL restget( inum4, 'gsigw', 1, 1, jpk, 0, llog, gsigw )  
         CALL restget( inum4, 'gsi3w', 1, 1, jpk, 0, llog, gsi3w )
         CALL restget( inum4, 'esigt', 1, 1, jpk, 0, llog, esigt )
         CALL restget( inum4, 'esigw', 1, 1, jpk, 0, llog, esigw )

# elif defined key_partial_steps
         !                                                          ! z-coordinate with partial steps
         CALL restget( inum4, 'hdept' , jpidta, jpjdta, 1, 0, llog, zta  )    !    ! depth
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hdept(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO
         CALL restget( inum4, 'hdepw' , jpidta, jpjdta, 1, 0, llog, zta  ) 
           DO jj = 1, nlcj
             DO ji = 1, nlci
               hdepw(ji,jj) = zta(mig(ji),mjg(jj))
             END DO
           END DO

         CALL restget( inum4, 'e3t_ps', jpidta, jpjdta, jpk, 0, llog, zt3a )  !    ! scale factors
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e3t_ps(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum4, 'e3u_ps', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e3u_ps(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum4, 'e3v_ps', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e3v_ps(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO
         CALL restget( inum4, 'e3w_ps', jpidta, jpjdta, jpk, 0, llog, zt3a )
         DO jk = 1, jpk
           DO jj = 1, nlcj
             DO ji = 1, nlci
               e3w_ps(ji,jj,jk) = zt3a(mig(ji),mjg(jj),jk)
             END DO
           END DO
         END DO

         CALL restget( inum4, 'gdept' , 1, 1, jpk, 0, llog, gdept )       !    ! reference z-coord.
         CALL restget( inum4, 'gdepw' , 1, 1, jpk, 0, llog, gdepw )
         CALL restget( inum4, 'e3t'   , 1, 1, jpk, 0, llog, e3t   )
         CALL restget( inum4, 'e3w'   , 1, 1, jpk, 0, llog, e3w   )

         DO jk=1,jpk
            gdept_ps(:,:,jk) = gdept(jk)
            gdepw_ps(:,:,jk) = gdepw(jk)
         END DO

         DO jj = 1, jpj
            DO ji = 1, jpi
               ik = mbathy(ji,jj) - 1
               ! ocean point only 
               IF( ik > 0 ) THEN
                  ! max ocean level case
                  gdepw_ps(ji,jj,ik+1) = hdepw(ji,jj)
                  gdept_ps(ji,jj,ik  ) = hdept(ji,jj)
                  gdept_ps(ji,jj,ik+1) = gdept_ps(ji,jj,ik) + e3t_ps(ji,jj,ik)
               ENDIF
            END DO
         END DO
         

# else
         !                                                          ! z-coordinate 
         CALL restget( inum4, 'gdept', 1, 1, jpk, 0, llog, gdept )        !    ! depth
         CALL restget( inum4, 'gdepw', 1, 1, jpk, 0, llog, gdepw )
         CALL restget( inum4, 'e3t'  , 1, 1, jpk, 0, llog, e3t   )        !    ! scale factors
         CALL restget( inum4, 'e3w'  , 1, 1, jpk, 0, llog, e3w   )
# endif

      ! Control printing : Grid informations (if not restart)
      ! ----------------

      IF(lwp .AND. .NOT.ln_rstart ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          longitude and e1 scale factors'
         WRITE(numout,*) '          ------------------------------'
         WRITE(numout,9300) ( ji, glamt(ji,1), glamu(ji,1),   &
            glamv(ji,1), glamf(ji,1),   &
            e1t(ji,1), e1u(ji,1),   &
            e1v(ji,1), ji = 1, jpi,10)
9300     FORMAT( 1x, i4, f8.2,1x, f8.2,1x, f8.2,1x, f8.2, 1x,    &
            f19.10, 1x, f19.10, 1x, f19.10 )

         WRITE(numout,*)
         WRITE(numout,*) '          latitude and e2 scale factors'
         WRITE(numout,*) '          -----------------------------'
         WRITE(numout,9300) ( jj, gphit(1,jj), gphiu(1,jj),   &
            &                     gphiv(1,jj), gphif(1,jj),   &
            &                     e2t  (1,jj), e2u  (1,jj),   &
            &                     e2v  (1,jj), jj = 1, jpj, 10 )
      ENDIF


      IF( nprint == 1 .AND. lwp ) THEN
         WRITE(numout,*) '          e1u e2u '
         CALL prihre( e1u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         CALL prihre( e2u,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         WRITE(numout,*) '          e1v e2v  '
         CALL prihre( e1v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
         CALL prihre( e2v,jpi,jpj,jpi-5,jpi,1,jpj-5,jpj,1,0.,numout )
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '              Reference z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level   gdept    gdepw     e3t      e3w  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, gdept(jk), gdepw(jk), e3t(jk), e3w(jk), jk = 1, jpk )
      ENDIF

      DO jk = 1, jpk
         IF( e3w(jk) <= 0. .OR. e3t(jk) <= 0. ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' e3w or e3t =< 0 '
            nstop = nstop + 1
         ENDIF
         IF( gdepw(jk) < 0. .OR. gdept(jk) < 0.) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' gdepw or gdept < 0 '
            nstop = nstop + 1
         ENDIF
      END DO

         !                                     ! ============================
         !                                     !        close the files 
         !                                     ! ============================
         SELECT CASE ( nmsh )
            CASE ( 1 )                
               CALL restclo( inum0 )
            CASE ( 2 )
               CALL restclo( inum1 )
               CALL restclo( inum2 )
            CASE ( 3 )
               CALL restclo( inum2 )
               CALL restclo( inum3 )
               CALL restclo( inum4 )
         END SELECT

   END SUBROUTINE dom_rea

#endif

   !!======================================================================
END MODULE domrea
