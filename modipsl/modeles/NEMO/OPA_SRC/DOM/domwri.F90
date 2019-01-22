!!DB -- 2009.09.04 -- key_diadimg eliminated
MODULE domwri
   !!======================================================================
   !!                       ***  MODULE domwri  ***
   !! Ocean initialization : write the ocean domain mesh ask file(s)
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dom_wri        : create and write mesh and mask file(s)
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
   PUBLIC dom_wri        ! routine called by inidom.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domwri.F90,v 1.9 2006/03/10 10:55:38 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS


   !!----------------------------------------------------------------------
   !!   Default option :                                        NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE dom_wri
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_wri  ***
      !!                   
      !! ** Purpose :   Create the NetCDF file(s) which contain(s) all the
      !!      ocean domain informations (mesh and mask arrays). This (these)
      !!      file(s) is (are) used for visualisation (SAXO software) and
      !!      diagnostic computation.
      !!
      !! ** Method  :   Write in a file all the arrays generated in routines
      !!      domhgr, domzgr, and dommsk. Note: the file contain depends on
      !!      the vertical coord. used (z-coord, partial steps, s-coord)
      !!                    nmsh = 1  :   'mesh_mask.nc' file
      !!                         = 2  :   'mesh.nc' and mask.nc' files
      !!                         = 3  :   'mesh_hgr.nc', 'mesh_zgr.nc' and
      !!                                  'mask.nc' files
      !!      For huge size domain, use option 2 or 3 depending on your 
      !!      vertical coordinate.
      !!
      !! ** output file : 
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
      INTEGER  ::                & !!! * temprary units for :
         inum0 ,                 &  ! 'mesh_mask.nc' file
         inum1 ,                 &  ! 'mesh.nc'      file
         inum2 ,                 &  ! 'mask.nc'      file
         inum3 ,                 &  ! 'mesh_hgr.nc'  file
         inum4                      ! 'mesh_zgr.nc'  file
      INTEGER  ::   itime           !  output from restini ???
      REAL(wp) ::   zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         zprt                       ! temporary array for bathymetry 

      CHARACTER (len=21) ::      &
         clnam0  ,   &  ! filename (mesh and mask informations)
         clnam1  ,   &  ! filename (mesh informations)
         clnam2  ,   &  ! filename (mask informations)
         clnam3  ,   &  ! filename (horizontal mesh informations)
         clnam4         ! filename (vertical   mesh informations)
      !!----------------------------------------------------------------------

       IF(lwp) WRITE(numout,*)
       IF(lwp) WRITE(numout,*) 'dom_wri : create NetCDF mesh and mask information file(s)'
       IF(lwp) WRITE(numout,*) '~~~~~~~'

         clnam0 = 'mesh_mask'  ! filename (mesh and mask informations)
         clnam1 = 'mesh'       ! filename (mesh informations)
         clnam2 = 'mask'       ! filename (mask informations)
         clnam3 = 'mesh_hgr'   ! filename (horizontal mesh informations)
         clnam4 = 'mesh_zgr'   ! filename (vertical   mesh informations)

#if defined key_agrif
      if ( .NOT. Agrif_Root() ) then
        clnam0 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam0)
        clnam1 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam1)
        clnam2 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam2)
        clnam3 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam3)
        clnam4 = TRIM(Agrif_CFixed())//'_'//TRIM(clnam4)
      endif
#endif

      CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )    ! calendar initialization

!       note that mbathy has been modified in dommsk or in solver.
!       it is the number of non-zero "w" levels in the water, and the minimum 
!       value (on land) is 2. We define zprt as the number of "T" points in the ocean 
!       at any location, and zero on land. 
!
      zprt = tmask(:,:,1)*(mbathy-1)

      SELECT CASE (nmsh)
         !                                     ! ============================
         CASE ( 1 )                            !  create 'mesh_mask.nc' file
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          one file in "mesh_mask.nc" '
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_mask.nc' file
            &             jpk   , gdept , trim(clnam0)        ,  &   ! in unit inum0
            &             itime , zdate0, rdt   , inum0 , domain_id=nidom )
            inum2 = inum0                                            ! put all the informations
            inum3 = inum0                                            ! in unit inum0
            inum4 = inum0

            !                                  ! ============================
         CASE ( 2 )                            !  create 'mesh.nc' and 
            !                                  !         'mask.nc' files
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          two files in "mesh.nc" and "mask.nc" '
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh.nc' file 
            &             jpk   , gdept , trim(clnam1)        ,  &   ! in unit inum1 
            &             itime , zdate0, rdt   , inum1, domain_id=nidom )
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file 
            &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2 
            &             itime , zdate0, rdt   , inum2, domain_id=nidom )
            inum3 = inum1                                            ! put mesh informations 
            inum4 = inum1                                            ! in unit inum1 

            !                                  ! ============================
         CASE ( 3 )                            !  create 'mesh_hgr.nc'
            !                                  !         'mesh_zgr.nc' and
            !                                  !         'mask.nc'     files
            !                                  ! ============================

            IF(lwp) WRITE(numout,*) '          three files in "mesh_hgr.nc" , mesh_zgr.nc" and "mask.nc" '
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_hgr.nc' file
            &             jpk   , gdept , trim(clnam3)        ,  &   ! in unit inum3
            &             itime , zdate0, rdt   , inum3, domain_id=nidom )
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mesh_zgr.nc' file
            &             jpk   , gdept , trim(clnam4)        ,  &   ! in unit inum4
            &             itime , zdate0, rdt   , inum4, domain_id=nidom )
            CALL restini( 'NONE', jpi   , jpj   , glamt, gphit,  &   ! create 'mask.nc' file
            &             jpk   , gdept , trim(clnam2)        ,  &   ! in unit inum2
            &             itime , zdate0, rdt   , inum2, domain_id=nidom )

         END SELECT

         !                                                         ! masks (inum2) 
         CALL restput( inum2, 'tmask', jpi, jpj, jpk, 0, tmask ) 
         CALL restput( inum2, 'umask', jpi, jpj, jpk, 0, umask )
         CALL restput( inum2, 'vmask', jpi, jpj, jpk, 0, vmask )
         CALL restput( inum2, 'fmask', jpi, jpj, jpk, 0, fmask )

         !                                                         ! horizontal mesh (inum3)
         CALL restput( inum3, 'glamt', jpi, jpj, 1, 0, glamt )     !    ! latitude
         CALL restput( inum3, 'glamu', jpi, jpj, 1, 0, glamu )
         CALL restput( inum3, 'glamv', jpi, jpj, 1, 0, glamv )
         CALL restput( inum3, 'glamf', jpi, jpj, 1, 0, glamf )

         CALL restput( inum3, 'gphit', jpi, jpj, 1, 0, gphit )     !    ! longitude
         CALL restput( inum3, 'gphiu', jpi, jpj, 1, 0, gphiu )
         CALL restput( inum3, 'gphiv', jpi, jpj, 1, 0, gphiv )
         CALL restput( inum3, 'gphif', jpi, jpj, 1, 0, gphif )

         CALL restput( inum3, 'e1t', jpi, jpj, 1, 0, e1t )         !    ! e1 scale factors
         CALL restput( inum3, 'e1u', jpi, jpj, 1, 0, e1u )
         CALL restput( inum3, 'e1v', jpi, jpj, 1, 0, e1v )
         CALL restput( inum3, 'e1f', jpi, jpj, 1, 0, e1f )

         CALL restput( inum3, 'e2t', jpi, jpj, 1, 0, e2t )         !    ! e2 scale factors
         CALL restput( inum3, 'e2u', jpi, jpj, 1, 0, e2u )
         CALL restput( inum3, 'e2v', jpi, jpj, 1, 0, e2v )
         CALL restput( inum3, 'e2f', jpi, jpj, 1, 0, e2f )

         CALL restput( inum3, 'ff', jpi, jpj, 1, 0, ff )           !    ! coriolis factor

         CALL restput( inum4, 'mbathy', jpi, jpj, 1, 0, zprt )

# if defined key_s_coord
         !                                                         ! s-coordinate
         CALL restput( inum4, 'hbatt', jpi, jpj, 1, 0, hbatt )      !    ! depth
         CALL restput( inum4, 'hbatu', jpi, jpj, 1, 0, hbatu ) 
         CALL restput( inum4, 'hbatv', jpi, jpj, 1, 0, hbatv )
         CALL restput( inum4, 'hbatf', jpi, jpj, 1, 0, hbatf )

         CALL restput( inum4, 'gsigt', 1, 1, jpk, 0, gsigt )        !    ! scaling coef.
         CALL restput( inum4, 'gsigw', 1, 1, jpk, 0, gsigw )  
         CALL restput( inum4, 'gsi3w', 1, 1, jpk, 0, gsi3w )
         CALL restput( inum4, 'esigt', 1, 1, jpk, 0, esigt )
         CALL restput( inum4, 'esigw', 1, 1, jpk, 0, esigw )

# elif defined key_partial_steps
         !                                                          ! z-coordinate with partial steps
         CALL restput( inum4, 'hdept' , jpi, jpj, 1, 0, hdept  )    !    ! depth
         CALL restput( inum4, 'hdepw' , jpi, jpj, 1, 0, hdepw  ) 

         CALL restput( inum4, 'e3t_ps', jpi, jpj, jpk, 0, e3t_ps )  !    ! scale factors
         CALL restput( inum4, 'e3u_ps', jpi, jpj, jpk, 0, e3u_ps )
         CALL restput( inum4, 'e3v_ps', jpi, jpj, jpk, 0, e3v_ps )
         CALL restput( inum4, 'e3w_ps', jpi, jpj, jpk, 0, e3w_ps )

         CALL restput( inum4, 'gdept' , 1, 1, jpk, 0, gdept )       !    ! reference z-coord.
         CALL restput( inum4, 'gdepw' , 1, 1, jpk, 0, gdepw )
         CALL restput( inum4, 'e3t'   , 1, 1, jpk, 0, e3t   )
         CALL restput( inum4, 'e3w'   , 1, 1, jpk, 0, e3w   )

# else
         !                                                          ! z-coordinate 
         CALL restput( inum4, 'gdept', 1, 1, jpk, 0, gdept )        !    ! depth
         CALL restput( inum4, 'gdepw', 1, 1, jpk, 0, gdepw )
         CALL restput( inum4, 'e3t'  , 1, 1, jpk, 0, e3t   )        !    ! scale factors
         CALL restput( inum4, 'e3w'  , 1, 1, jpk, 0, e3w   )
# endif

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

   END SUBROUTINE dom_wri


   !!======================================================================
END MODULE domwri
